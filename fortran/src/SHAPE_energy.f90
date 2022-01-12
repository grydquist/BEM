MODULE SHAPEMOD
USE HARMMOD
IMPLICIT NONE

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! The actual cell shape object
TYPE cellType

!   Current cell
    INTEGER :: id

!   Material Properties
    REAL(KIND = 8) :: mu, lam, B, C, Eb, c0, Ca, int_pres

!   Harmonics info
    INTEGER :: p, q, ftot
    TYPE(YType), POINTER :: Y, Yf, Ys
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:)

!   Geometric info
    REAL(KIND = 8), ALLOCATABLE :: J(:,:), x(:,:,:), xf(:,:,:), k(:,:)
    REAL(KIND = 8) :: V0, h
!   Derivatives
    REAL(KIND = 8), ALLOCATABLE :: dxt(:,:,:), dxp(:,:,:), dxp2(:,:,:), &
    dxt2(:,:,:), dxtp(:,:,:), dxp3(:,:,:), dxt3(:,:,:), dxt2p(:,:,:), &
    dxtp2(:,:,:), dxp4(:,:,:), dxtp3(:,:,:), dxt2p2(:,:,:), dxt3p(:,:,:), dxt4(:,:,:)

!   Reference state variables
    REAL(KIND = 8), ALLOCATABLE :: kR(:,:), kdR(:,:,:), kd2R(:,:,:), &
    c1R(:,:,:), c2R(:,:,:), c1tR(:,:,:), c1pR(:,:,:), c2tR(:,:,:), &
    c2pR(:,:,:), JR(:,:)

!   Force variables
    REAL(KIND = 8), ALLOCATABLE :: fab(:,:,:), ff(:,:,:), fc(:,:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:), nkmn(:,:), fmn2(:,:), nkt(:,:)

!   Cell velocity constants
    COMPLEX(KIND = 8), ALLOCATABLE :: umn(:,:)    

!   Has this object been initialized?
    LOGICAL :: init = .false.
    LOGICAL :: writ = .false.

!   Name of output file
    CHARACTER(:), ALLOCATABLE :: fileout
    
    CONTAINS
    PROCEDURE :: Derivs  => Derivscell
    PROCEDURE :: Stress  => Stresscell
    PROCEDURE :: Fluid   => Fluidcell
    PROCEDURE :: Relax   => RelaxCell
    PROCEDURE :: Dealias => Dealiascell
    PROCEDURE :: Vol     => Volcell
    PROCEDURE :: SA      => SAcell
    PROCEDURE :: Intg    => Intgcell
END TYPE cellType

! -------------------------------------------------------------------------!
! Contains all the miscellaneous info about the problem
TYPE probType
    INTEGER :: cts, NT, dtinc, NCell, Nmat, NmatT
    REAL(KIND = 8) :: dt, t
!   Velocity gradient & its file location
    REAL(KIND = 8) :: dU(3,3)
    CHARACTER(len = 15) gradfile

!   Normalized Legendres/exp's at GPs
    COMPLEX(KIND =8), ALLOCATABLE :: es(:,:), esf(:,:)
    REAL(KIND = 8), ALLOCATABLE :: cPmn(:,:,:), cPmnf(:,:,:), xsf(:)

!   Harmonics info
    INTEGER :: p, q, ftot
    TYPE(YType) :: Y, Yf, Ys
    REAL(KIND = 8) :: thtc, k, h

!   Pointer to the cells
    TYPE(cellType), POINTER :: cell(:)

!   MPI stuff
    TYPE(cmType), POINTER :: cm
    INTEGER :: PCells(2)

!   Are we continuing from somewhere?
    LOGICAL :: cont = .false.

    CONTAINS
    PROCEDURE :: Update  => UpdateProb
    PROCEDURE :: Write   => WriteProb
    PROCEDURE :: Output  => OutputProb
    PROCEDURE :: Continue=> ContinueProb
END TYPE probType

! -------------------------------------------------------------------------!
INTERFACE cellType
    PROCEDURE :: newcell
END INTERFACE cellType

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

! Constructs cell object, takes in order
FUNCTION newcell(filein, reduce, prob) RESULT(cell)
    CHARACTER(len = *), INTENT(IN) :: filein
    LOGICAL, INTENT(IN) :: reduce
    TYPE(probType), TARGET, INTENT(INOUT) :: prob
    TYPE(cellType), ALLOCATABLE, TARGET :: cell(:)
    CHARACTER(len = 3) :: restart, cont
    CHARACTER(len = 30) :: restfile, icC, contfile, cfile2
    CHARACTER(:), ALLOCATABLE :: fileout
    REAL(KIND = 8), ALLOCATABLE :: cPt(:,:), ths(:,:), phs(:,:), thts(:), phis(:), xs(:), wg(:)
    REAL(KIND = 8) lam, Ca, C, Eb, c0, dphi, int_pres
    INTEGER :: nt, np, ntf, npf, fali, p, m, ind, n, it, im2, im, ic, stat

!   General problem parameters
    prob%NT = READ_GRINT_INT(filein, 'Max_time_steps')
    prob%NCell = READ_GRINT_INT(filein, 'Number_cells')
    prob%dt = READ_GRINT_DOUB(filein, 'Time_step')
    prob%dtinc = READ_GRINT_INT(filein, 'Time_inc')
    prob%cts = 0
    prob%t = 0D0
!   Gradient file location
    prob%gradfile = READ_GRINT_CHAR(filein, 'Gradient_file')

!   Check if there are more processors than cells (can't handle)
    IF(prob%cm%np() .gt. prob%NCell) THEN
        print *, 'ERROR: More processors than cells'
        STOP
    ENDIF

!   Material properties from input
    lam = READ_GRINT_DOUB(filein, 'Viscosity_Ratio')
    Ca = READ_GRINT_DOUB(filein, 'Capillary')
    C = READ_GRINT_DOUB(filein, 'Dilatation_Ratio')
    Eb = READ_GRINT_DOUB(filein, 'Bending_Modulus')
    c0 = READ_GRINT_DOUB(filein, 'Spont_Curvature')
    int_pres = READ_GRINT_DOUB(filein, 'Internal_Pressure')

!   Write location
    fileout = TRIM(READ_GRINT_CHAR(filein, 'Output'))

!   Coarse and fine grids    
    p = READ_GRINT_INT(filein, 'Harmonic_order')
    prob%p = p
    fali   = READ_GRINT_INT(filein, 'Refinement_factor')
    prob%q = p*fali

!   Restart location
    restart = READ_GRINT_CHAR(filein, 'Restart')
!   Choose file to restart from (assuming that this is after deflation)
    restfile = READ_GRINT_CHAR(filein, 'Restart_Loc')
!   Should we continue from a  previous file?
    cont = READ_GRINT_CHAR(filein, 'Continue')

!   Note that the file MUST be in the restart directory!
    restfile = 'restart/'//trim(restfile)

!   Make harmonics(order, # of derivs, if we'll rotate or not)
    prob%Y = YType(p, 1, .true.)
    prob%Yf = YType(p*fali, 4, .true., p)
    nt  = prob%Y%nt
    np  = prob%Y%np
    ntf = prob%Yf%nt
    npf = prob%Yf%np

!   Harmonics for the singular integration, slightly trickier
!   Construct the singular integration grid, with patch. Essentially 2 grids at once,
!   a fine patch with a sinh transform and a coarse one.
    IF(prob%NCell.gt.1) THEN
        ALLOCATE(thts(nt + ntf), phis(np + npf), &
                 ths (nt + ntf, np + npf), phs(nt + ntf, np + npf))
!       Cutoff theta, can affect accuracy. Depends on spacing, but want consistent. Just hardcode for now
        prob%thtc = pi/6D0 !!!

!       Coarser part away from near-singularity
        ALLOCATE(xs(nt), wg(nt))
        CALL lgwt(nt, xs, wg)
        thts(1:nt) = xs*(pi - prob%thtc)/2D0 + (pi + prob%thtc)/2D0
        DEALLOCATE(xs, wg)

!       Finer part near near-singularity
        ALLOCATE(xs(ntf), wg(ntf))
        CALL lgwt(ntf, xs, wg)
!       The integration rule is based on the separation distance. However, this changes and
!       we don't want to recalculate the grid every time. Instead, just choose a distance
!       where it's accurate (approx spacing/10 here)
        prob%h = SQRT(PI)/10D0/nt
!       Sinh normalizing constant (to get -1, 1 range)
        prob%k = -0.5D0*LOG(prob%thtc/prob%h + sqrt((prob%thtc/prob%h)**2D0 + 1D0));
        thts(nt + 1:nt + ntf) = prob%h*SINH(prob%k*(xs - 1D0))
        prob%xsf = xs
        DEALLOCATE(xs, wg)

!       Calculate phis and construct mesh grid
        phis(1) = 0D0
        dphi = prob%Y%dphi
        DO ic = 1,(np + npf)
            IF(ic .gt. 1) phis(ic) = phis(ic - 1) + dphi

            IF(ic .eq. np + 1) THEN
                dphi = prob%Yf%dphi
                phis(ic) = 0D0
            ENDIF
            phs(:,ic) = phis(ic)
            ths(:,ic) = thts
        ENDDO

!       Create a bare harmonics object evaluated at this grid
        prob%Ys = YType(ths, phs, p)
    ENDIF

    ALLOCATE(prob%es(2*(p-1)+1, np), prob%cPmn(p, 2*(p-1)+1, nt), cPt(nt, p*(p+1)/2))
    ALLOCATE(prob%esf(2*(p-1)+1, npf + np), prob%cPmnf(p, 2*(p-1)+1, ntf + nt))

!   Matrix size
    prob%Nmat = 3*(prob%Y%p)*(prob%Y%p)
    prob%NmatT= prob%Nmat*prob%NCell

    ALLOCATE(cell(prob%NCell))

    DO ic = 1, prob%NCell
        cell(ic)%id = ic
        cell(ic)%mu = 1D0
        cell(ic)%lam = lam

!       To fit dimesnionless parameters, we set a0 = 1, flow timescale = 1, mu = 1
!       and fit the rest from there
!       Ca = (shear rate)*mu*(cell radius)/(B/2)
        cell(ic)%Ca = Ca
        cell(ic)%B = 2D0/Ca
!       A note on the dilation modulus: many papers use (B*C) in the spot where
!       I use C, so I just make that correction here
        cell(ic)%C = C*cell(ic)%B
!       Similar situation for beinding modulus: the input is the non-dim
!       parameter E*b = Eb/(a_0^2*(B/2))
        cell(ic)%Eb = Eb*2D0*cell(ic)%B
        cell(ic)%c0 = c0
        write(icC, "(I0.1)") ic
        cell(ic)%fileout = TRIM(fileout//icC)
!       Internal (osmotic) pressure
        cell(ic)%int_pres = int_pres

!       Coarse and fine grids/harms
        cell(ic)%p = p
        cell(ic)%q = p*fali
        cell(ic)%Y => prob%Y
        cell(ic)%Yf => prob%Yf
        cell(ic)%Ys => prob%Ys

!       Allocating everything, starting with derivatives
        ALLOCATE(cell(ic)%dxt(3,ntf,npf), cell(ic)%dxp(3,ntf,npf), cell(ic)%dxp2(3,ntf,npf), &
                cell(ic)%dxt2(3,ntf,npf), cell(ic)%dxtp(3,ntf,npf), cell(ic)%dxp3(3,ntf,npf), &
                cell(ic)%dxt3(3,ntf,npf), cell(ic)%dxt2p(3,ntf,npf), cell(ic)%dxtp2(3,ntf,npf), &
                cell(ic)%dxp4(3,ntf,npf), cell(ic)%dxtp3(3,ntf,npf), cell(ic)%dxt2p2(3,ntf,npf), &
                cell(ic)%dxt3p(3,ntf,npf), cell(ic)%dxt4(3,ntf,npf), &
                cell(ic)%J(ntf,npf), cell(ic)%x(3,nt,np), cell(ic)%xf(3,ntf,npf), &
                cell(ic)%xmn(3,(p+1)*(p+1)), cell(ic)%umn(3,(p+1)*(p+1)))
!       Reference items
        ALLOCATE(cell(ic)%kR(ntf,npf), cell(ic)%kdR(2,ntf,npf), cell(ic)%kd2R(3,ntf,npf), &
                cell(ic)%c1R(3,ntf,npf), cell(ic)%c2R(3,ntf,npf), cell(ic)%c1tR(3,ntf,npf), &
                cell(ic)%c1pR(3,ntf,npf), cell(ic)%c2tR(3,ntf,npf), cell(ic)%c2pR(3,ntf,npf))

!       Force items
        ALLOCATE(cell(ic)%fab(3,ntf,npf), cell(ic)%ff(3,ntf,npf), cell(ic)%fc(3,nt,np), &
                 cell(ic)%fmn(3,(cell(ic)%q+1)*(cell(ic)%q+1)), cell(ic)%nkmn(3,(cell(ic)%q+1)*(cell(ic)%q+1)), &
                 cell(ic)%fmn2(3,(cell(ic)%q+1)*(cell(ic)%q+1)), cell(ic)%nkt(3,(cell(ic)%q+1)*(cell(ic)%q+1)))

        cell(ic)%ff = 0D0
        cell(ic)%fmn = 0D0
        cell(ic)%umn = 0D0

!       For a = 1, V = 4.18904795321178, SA = 16.8447913187040, sphere 6.50088174342271

!       First, we need to get the reference shape for the shear stress
        ! cell(ic)%xmn = RBCcoeff(cell(ic)%Y)
        ! cell(ic)%xmn = Cubecoeff(cell(ic)%Y)
        ! cell(ic)%xmn = Bactcoeff(cell(ic)%Y, 1D0)
        cell(ic)%xmn = Spherecoeff(cell(ic)%Y, .76D0) !1D0)!!!! Reduced volume .997, .98, .95-> .9,.76,.65

!       Initial surfce derivatives/reference state
        CALL cell(ic)%Derivs()
        CALL cell(ic)%Stress()

!       Now scale to get an equivalent radius of 1
        IF(.not. reduce) THEN
            cell(ic)%xmn = cell(ic)%xmn/(3D0*cell(ic)%vol()/(4D0*pi))**(1D0/3D0)
!       If we want to do deflation, comment the above and scale to get right SA
        ELSE
            cell(ic)%xmn = cell(ic)%xmn*(16.8447913187040D0/cell(ic)%SA())**(1D0/2D0) 
        ENDIF

!       Initialize derivs again... I do this because I need the derivs to get the volume
        CALL cell(ic)%Derivs()
        CALL cell(ic)%Stress()

!       Now with a reference shape made, should we choose an intermediate point
!       to restart from?
        IF (restart .eq. "Yes") THEN
!           This gives us the right surface area with a blown up volume
!           The 16.8 is for equiv radius of 1
            cell(ic)%xmn = cell(ic)%xmn*sqrt(16.8447913187040D0/cell(ic)%SA())
!           Get referennce state
            CALL cell(ic)%Derivs()
            CALL cell(ic)%Stress()

!           And just set the constants equal, scaling for equiv radius
!           This assumes the equiv radius in input is 1.
            cell(ic)%xmn = Readcoeff(restfile, cell(ic)%Y%p)
        ENDIF

!       Should we continue from a spot where we left off? !! Not working right now,
        !! doesn't return the same values for some reason
        IF(cont .eq. "Yes") prob%cont = .true.

        !! Test                                        !!!!!
        ! SELECT CASE(ic)
        ! CASE(1)
        !     cell(ic)%xmn(3,1) =  0D0!  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = -5D0*SQRT(pi)*3D0/5D0
        ! CASE(2)
        !     cell(ic)%xmn(3,1) = 0D0!  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = 5D0*SQRT(pi)*3D0/5D0
        ! CASE(3)
        !     cell(ic)%xmn(3,1) =  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = -5D0*SQRT(pi)*3D0/2D0
        ! CASE(4)
        !     cell(ic)%xmn(3,1) =  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  0D0! -5D0*SQRT(pi)*3D0/3D0
        ! CASE(5)
        !     cell(ic)%xmn(3,1) =  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  5D0*SQRT(pi)*3D0/2D0
        ! CASE(6)
        !     cell(ic)%xmn(3,1) = -5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = -5D0*SQRT(pi)*3D0/2D0
        ! CASE(7)
        !     cell(ic)%xmn(3,1) = -5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  0D0! -5D0*SQRT(pi)*3D0/3D0
        ! CASE(8)
        !     cell(ic)%xmn(3,1) = -5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  5D0*SQRT(pi)*3D0/2D0
        ! END SELECT

!       Initial positions
        cell(ic)%x(1,:,:) = cell(ic)%Y%backward(cell(ic)%xmn(1,:))
        cell(ic)%x(2,:,:) = cell(ic)%Y%backward(cell(ic)%xmn(2,:))
        cell(ic)%x(3,:,:) = cell(ic)%Y%backward(cell(ic)%xmn(3,:))

!       Characteristic grid spacing
        cell(ic)%h = SQRT(cell(ic)%SA()/nt**2D0)
    ENDDO
    
!   Exponential calculation part (coarse and fine)
    DO m = -(p-1),(p-1)
        ind = m + p
        prob%es(ind,:) = EXP(ii*DBLE(m)*prob%Y%phi)*prob%Y%dphi
    ENDDO

!   Legendre polynomial calculation part (coarse and fine)
    cPt = Alegendre(p-1,COS(prob%Y%tht))
    it = 0
    DO n = 0, p-1
        ind = n+1
        im = 0
        DO m = -(p-1),p-1
            im = im + 1
            IF(ABS(m) .gt. n) THEN
                prob%cPmn(ind,im,:) = 0D0
            ELSEIF(m .le. 0) THEN
                it = it + 1
                IF(m.eq.-n) im2 = it
                prob%cPmn(ind,im,:) = (-1D0)**m*cPt(:, im2 + abs(m))
            ELSE
                prob%cPmn(ind,im,:) = (-1D0)**m*prob%cPmn(ind, im - 2*m, :)
            ENDIF
        ENDDO
    ENDDO

    IF(prob%NCell .gt. 1) THEN
!   Fine part, essentially done on 2 grids
!   Technically the grid goes up to order q, but we only calculate up to p
    DO m = -(p-1),(p-1)
        ind = m + p
        prob%esf(ind,:) = EXP(ii*DBLE(m)*prob%Ys%ph(1,:))
    ENDDO

!   Manage the dphi
    DO ic = 1,np + npf
        IF(ic .le. np) THEN
            prob%esf(:,ic) = prob%esf(:,ic)*prob%Y%dphi
        ELSE
            prob%esf(:,ic) = prob%esf(:,ic)*prob%Yf%dphi
        ENDIF
    ENDDO
    DEALLOCATE(cPt)
!   Technically the grid goes up to order q, but we only calculate up to p
    ALLOCATE(cPt(nt + ntf, p*(p+1)/2))
    cPt = Alegendre(p-1,COS(prob%Ys%th(:,1)))
    it = 0
    DO n = 0, p-1
        ind = n+1
        im = 0
        DO m = -(p-1),p-1
            im = im + 1
            IF(ABS(m) .gt. n) THEN
                prob%cPmnf(ind,im,:) = 0D0
            ELSEIF(m .le. 0) THEN
                it = it + 1
                IF(m.eq.-n) im2 = it
                prob%cPmnf(ind,im,:) = (-1D0)**m*cPt(:, im2 + abs(m))
            ELSE
                prob%cPmnf(ind,im,:) = (-1D0)**m*prob%cPmnf(ind, im - 2*m, :)
            ENDIF
        ENDDO
    ENDDO
    ENDIF

!   Which cells will the current processor handle? Integer div. rounds down
    n = MOD(prob%NCell, prob%cm%np())
    m = prob%NCell/prob%cm%np()

!   Get the top and bottom cell indices for this processor,
!   giving the remainder to the first processors
    IF (prob%cm%id() .lt. n) THEN
        prob%PCells(1) = (prob%cm%id()    )*(m + 1) + 1
        prob%PCells(2) = (prob%cm%id() + 1)*(m + 1)
    ELSE
        prob%PCells(1) = n*(m + 1) + (prob%cm%id()     - n)*m + 1
        prob%PCells(2) = n*(m + 1) + (prob%cm%id() + 1 - n)*m
    ENDIF

    cell%init = .true.
END FUNCTION newcell

! -------------------------------------------------------------------------!
! Writes xmn to a text file. Could be better
SUBROUTINE WriteProb(prob)
    TYPE(cellType), POINTER :: cell
    CLASS(probType), INTENT(IN) :: prob
    CHARACTER (LEN = 25) ctsst, datdir, filename
    INTEGER ic
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:)

    IF(prob%cm%slv()) RETURN
!   Don't write if it's not a timestep to write
    IF(.not.((prob%cts .eq. 1) .or. &
         (MOD(prob%cts,prob%dtinc)) .eq. 0)) RETURN

    DO ic = 1,prob%NCell
    cell => prob%cell(ic)
    IF(NOT(ALLOCATED(fmn))) THEN
        ALLOCATE(fmn(3,(cell%q+1)*(cell%q+1)))
    ENDIF

!   Formatting pain
    write(ctsst, "(I0.5)") prob%cts
    datdir = TRIM('dat/'//cell%fileout//'/')
    filename = TRIM('x_'//ctsst)

!   Write position - first half - real, x,y,z groups, second half imag
    IF(.not. cell%writ) THEN
        CALL MAKEDIRQQ(datdir)
    ENDIF
    OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    WRITE(88,*) REAL(cell%xmn)
    WRITE(88,*) AIMAG(cell%xmn)
    CLOSE(88)

!   Write velocity
    filename = TRIM('u_'//ctsst)
    OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    WRITE(88,*) REAL(cell%umn)!REAL(cell%fmn(:,1:((cell%p+1)*(cell%p+1)))) !
    WRITE(88,*) AIMAG(cell%umn)!AIMAG(cell%fmn(:,1:((cell%p+1)*(cell%p+1))))!
    CLOSE(88)

!   Write force
    ! fmn(1,:) = cell%Yf%forward(cell%ff(1,:,:), cell%q)
    ! fmn(2,:) = cell%Yf%forward(cell%ff(2,:,:), cell%q)
    ! fmn(3,:) = cell%Yf%forward(cell%ff(3,:,:), cell%q)

    ! filename = TRIM('f_'//ctsst)
    ! OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    ! WRITE(88,*) REAL(fmn(:,1:((cell%p+1)*(cell%p+1))))
    ! WRITE(88,*) AIMAG(fmn(:,1:((cell%p+1)*(cell%p+1))))
    ! CLOSE(88)

!   Another file that just has all the material constants for the simulation
    IF(.not. cell%writ) THEN
        filename = 'Params'
        OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
        WRITE(88,*) "p"
        WRITE(88,*) cell%p
        WRITE(88,*) "dt"
        WRITE(88,*) prob%dt
        WRITE(88,*) "dt_inc"
        WRITE(88,*) prob%dtinc
        WRITE(88,*) "Ed"
        WRITE(88,*) cell%C/cell%B
        WRITE(88,*) "Eb"
        WRITE(88,*) cell%Eb/2D0/cell%B
        WRITE(88,*) "Ca"
        WRITE(88,*) cell%Ca
        WRITE(88,*) "c0"
        WRITE(88,*) cell%c0
        WRITE(88,*) "lambda"
        WRITE(88,*) cell%lam
        CLOSE(88)
    ENDIF
    
!   Not optimal, but just a file for the max timestep
    filename = 'maxdt'
    OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    WRITE(88,*) prob%cts
    CLOSE(88)
    
    IF(.not. cell%writ) cell%writ = .true.
    ENDDO

END SUBROUTINE WriteProb

! -------------------------------------------------------------------------!
! Updates the values of the derivatives on the surface of the cell on fine grid.
SUBROUTINE Derivscell(cell)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell

    TYPE(nmType), POINTER :: nm
    TYPE(YType), POINTER :: Y
    INTEGER :: ih, it, im, n, m, p
    COMPLEX(KIND = 8) :: f1, f2, f3
    COMPLEX(KIND = 8), ALLOCATABLE :: vcur(:,:), td1(:,:), td2(:,:), td3(:,:), td4(:,:)

    Y => cell%Yf
    p = cell%Y%p
    ALLOCATE(vcur(Y%nt, Y%np), td1(Y%nt, Y%np), td2(Y%nt, Y%np), &
    td3(Y%nt, Y%np), td4(Y%nt, Y%np))

!   Initialize variables
    cell%xf     = 0D0
    cell%dxt    = 0D0
    cell%dxp    = 0D0
    cell%dxt2   = 0D0
    cell%dxp2   = 0D0
    cell%dxtp   = 0D0
    cell%dxt3   = 0D0
    cell%dxp3   = 0D0
    cell%dxt2p  = 0D0
    cell%dxtp2  = 0D0
    cell%dxp4   = 0D0
    cell%dxtp3  = 0D0
    cell%dxt2p2 = 0D0
    cell%dxt3p  = 0D0
    cell%dxt4   = 0D0
    
    ih = 0
    it = 0
!   Loop over harmonics, but only up to coarse grid order, since the rest are 0
    DO n = 0, p
        im = 0
        ih = ih + 1

!       Values at current order        
        nm =>Y%nm(ih)
        DO m = -n,n ! Could exploit symmetry here... but it isn't a big deal since it takes so little time
            it = it + 1
            im = im + 1

!           Current harrmonic coeff
            f1 = cell%xmn(1,it)
            f2 = cell%xmn(2,it)
            f3 = cell%xmn(3,it)

!           Values at current order and degree
            vcur = nm%v(im,:,:)
!           Theta derivatives
            td1 = nm%dY(1,im,:,:)
            td2 = nm%dY(2,im,:,:)
            td3 = nm%dY(3,im,:,:)
            td4 = nm%dY(4,im,:,:)

!           Distance from zero
            cell%xf(1,:,:) = cell%xf(1,:,:) + REAL(f1*vcur)
            cell%xf(2,:,:) = cell%xf(2,:,:) + REAL(f2*vcur)
            cell%xf(3,:,:) = cell%xf(3,:,:) + REAL(f3*vcur)

!           Theta 1
            cell%dxt(1,:,:) = cell%dxt(1,:,:) + REAL(f1*td1)
            cell%dxt(2,:,:) = cell%dxt(2,:,:) + REAL(f2*td1)
            cell%dxt(3,:,:) = cell%dxt(3,:,:) + REAL(f3*td1)

!           Phi 1
            cell%dxp(1,:,:) = cell%dxp(1,:,:) + REAL(f1*ii*m*vcur)
            cell%dxp(2,:,:) = cell%dxp(2,:,:) + REAL(f2*ii*m*vcur)
            cell%dxp(3,:,:) = cell%dxp(3,:,:) + REAL(f3*ii*m*vcur)

!           Theta 2 
            cell%dxt2(1,:,:) = cell%dxt2(1,:,:) + REAL(f1*td2)
            cell%dxt2(2,:,:) = cell%dxt2(2,:,:) + REAL(f2*td2)
            cell%dxt2(3,:,:) = cell%dxt2(3,:,:) + REAL(f3*td2)
            
!           Phi 2
            cell%dxp2(1,:,:) = cell%dxp2(1,:,:) + REAL(-f1*m*m*vcur)
            cell%dxp2(2,:,:) = cell%dxp2(2,:,:) + REAL(-f2*m*m*vcur)
            cell%dxp2(3,:,:) = cell%dxp2(3,:,:) + REAL(-f3*m*m*vcur)
            
!           Theta 1 phi 1
            cell%dxtp(1,:,:) = cell%dxtp(1,:,:) + REAL(f1*ii*m*td1)
            cell%dxtp(2,:,:) = cell%dxtp(2,:,:) + REAL(f2*ii*m*td1)
            cell%dxtp(3,:,:) = cell%dxtp(3,:,:) + REAL(f3*ii*m*td1)

!           Theta 3
            cell%dxt3(1,:,:) = cell%dxt3(1,:,:) + REAL(f1*td3)
            cell%dxt3(2,:,:) = cell%dxt3(2,:,:) + REAL(f2*td3)
            cell%dxt3(3,:,:) = cell%dxt3(3,:,:) + REAL(f3*td3)
            
!           Theta 2 Phi 1
            cell%dxt2p(1,:,:) = cell%dxt2p(1,:,:) + REAL(f1*ii*m*td2)
            cell%dxt2p(2,:,:) = cell%dxt2p(2,:,:) + REAL(f2*ii*m*td2)
            cell%dxt2p(3,:,:) = cell%dxt2p(3,:,:) + REAL(f3*ii*m*td2)
            
!           Theta 1 Phi 2
            cell%dxtp2(1,:,:) = cell%dxtp2(1,:,:) + REAL(-f1*m*m*td1)
            cell%dxtp2(2,:,:) = cell%dxtp2(2,:,:) + REAL(-f2*m*m*td1)
            cell%dxtp2(3,:,:) = cell%dxtp2(3,:,:) + REAL(-f3*m*m*td1)
            
!           Phi 3
            cell%dxp3(1,:,:) = cell%dxp3(1,:,:) + REAL(-f1*ii*m*m*m*vcur)
            cell%dxp3(2,:,:) = cell%dxp3(2,:,:) + REAL(-f2*ii*m*m*m*vcur)
            cell%dxp3(3,:,:) = cell%dxp3(3,:,:) + REAL(-f3*ii*m*m*m*vcur)
            
!           Theta 4
            cell%dxt4(1,:,:) = cell%dxt4(1,:,:) + REAL(f1*td4)
            cell%dxt4(2,:,:) = cell%dxt4(2,:,:) + REAL(f2*td4)
            cell%dxt4(3,:,:) = cell%dxt4(3,:,:) + REAL(f3*td4)

!           Theta 3 Phi 1
            cell%dxt3p(1,:,:) = cell%dxt3p(1,:,:) + REAL(f1*ii*m*td3)
            cell%dxt3p(2,:,:) = cell%dxt3p(2,:,:) + REAL(f2*ii*m*td3)
            cell%dxt3p(3,:,:) = cell%dxt3p(3,:,:) + REAL(f3*ii*m*td3)

!           Theta 2 Phi 2
            cell%dxt2p2(1,:,:) = cell%dxt2p2(1,:,:) + REAL(-f1*m*m*td2)
            cell%dxt2p2(2,:,:) = cell%dxt2p2(2,:,:) + REAL(-f2*m*m*td2)
            cell%dxt2p2(3,:,:) = cell%dxt2p2(3,:,:) + REAL(-f3*m*m*td2)

!           Theta 1 Phi 3
            cell%dxtp3(1,:,:) = cell%dxtp3(1,:,:) + REAL(-f1*ii*m*m*m*td1)
            cell%dxtp3(2,:,:) = cell%dxtp3(2,:,:) + REAL(-f2*ii*m*m*m*td1)
            cell%dxtp3(3,:,:) = cell%dxtp3(3,:,:) + REAL(-f3*ii*m*m*m*td1)
            
!           Phi 4
            cell%dxp4(1,:,:) = cell%dxp4(1,:,:) + REAL(f1*m*m*m*m*vcur)
            cell%dxp4(2,:,:) = cell%dxp4(2,:,:) + REAL(f2*m*m*m*m*vcur)
            cell%dxp4(3,:,:) = cell%dxp4(3,:,:) + REAL(f3*m*m*m*m*vcur)
        ENDDO
    ENDDO

END SUBROUTINE Derivscell

! -------------------------------------------------------------------------!
! Gets the force jump at the locations on the cell, de-aliased!
SUBROUTINE Stresscell(cell)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell

    INTEGER :: i, j, q, rc
    REAL(KIND = 8) :: nk(3), E, F, G, L, M, N, D, k, dnt(3), dnp(3), &
                      dnt2(3), dnp2(3), dntp(3), gv(2,2), gn(2,2), &
                      dgt(2,2), dgp(2,2), dgt2(2,2), dgp2(2,2), dgtp(2,2), &
                      c1(3), c2(3), c1t(3), c1p(3), c2t(3), c2p(3)
    REAL(KIND = 8) :: Et, Ep, Et2, Ep2, Etp, Gt, Gp, Gt2, Gp2, Gtp, &
                      Ft, Fp, Ft2, Fp2, Ftp, Lt, Lp, Lt2, Lp2, Ltp, &
                      Nt, Np, Nt2, Np2, Ntp, Mt, Mp, Mt2, Mp2, Mtp, &
                      Jt, Jp, Jt2, Jp2, Jtp, Dt, Dp, Dt2, Dp2, Dtp, &
                      kt, kp, kt2, kp2, ktp
    REAL(KIND = 8) :: ct2_p(3), ct_tp(3), ctp_p(3), ct_p2(3), &
                      ct3_p(3), ct2_tp(3), ct_t2p(3), &
                      ct2p_p(3), ct2_p2(3), ct_tp2(3), &
                      ctp2_p(3), ctp_p2(3), ct_p3(3)
    REAL(KIND = 8) :: c111,  c112,  c122,  c222,  c221,  c211, &
                      c111t, c112t, c122t, c222t, c221t, &
                      c111p, c112p,        c222p, c221p, c211p
    REAL(KIND = 8) :: Fd(3,3), Prj(3,3), V2(3,3), &
                      lams(3), es(2), ev(3,3), ev1(3), ev2(3), I1, I2, &
                      tau(3,3), tauab(2,2), B, C, Eb, wrk(102), tV2(3,3)
    REAL(KIND = 8) :: Fdt(3,3), Fdp(3,3), V2t(3,3), V2p(3,3), lam1dV2(3,3), &
                      lam2dV2(3,3), es1t, es1p, es2t, es2p, I1t, I1p, I2t, I2p, &
                      Prjt(3,3), Prjp(3,3), taut(3,3), taup(3,3), dtauab(2,2), &
                      bv(2,2), bm(2,2), cvt, cvp, &
                      fb, LBk, kG, c0, fbt(3)

!!  Integrating strain energy
    REAL(KIND = 8), ALLOCATABLE :: Ws(:,:), Wb(:,:)
    ALLOCATE(Ws(cell%Yf%nt, cell%Yf%np), Wb(cell%Yf%nt, cell%Yf%np))
    IF(.not.ALLOCATED(cell%JR)) ALLOCATE(cell%JR(cell%Yf%nt,cell%Yf%np))
    
    B = cell%B
    C = cell%C
    Eb = cell%Eb
    c0 = cell%c0
    
!   Big loop over all points in grid
    DO i = 1, cell%Yf%nt
        inner:DO j = 1,cell%Yf%np
!           Normal vector (inward)
            nk = cross(cell%dxt(:,i,j), cell%dxp(:,i,j))
            nk = -nk/sqrt(nk(1)*nk(1) + nk(2)*nk(2) + nk(3)*nk(3))

!           Fundamental forms/area element
            E = dot(cell%dxt(:,i,j),cell%dxt(:,i,j))
            F = dot(cell%dxt(:,i,j),cell%dxp(:,i,j))
            G = dot(cell%dxp(:,i,j),cell%dxp(:,i,j))
            cell%J(i,j) = sqrt(E*G - F*F)

!           Second fundamental forms/curvature
            L = dot(cell%dxt2(:,i,j),-nk)
            M = dot(cell%dxtp(:,i,j),-nk)
            N = dot(cell%dxp2(:,i,j),-nk)
!           Curvature numerator
            D = E*N - 2D0*F*M + G*L
            k = 0.5D0*D/(cell%J(i,j)*cell%J(i,j))
!           Gaussian curvature
            kG = (L*N - M*M)/(cell%J(i,j)*cell%J(i,j))

!           We need several cross products in the below equations
!           so I do them all here
            ct2_p = cross(cell%dxt2(:,i,j),cell%dxp(:,i,j))
            ct_tp = cross(cell%dxt (:,i,j),cell%dxtp(:,i,j))
            ctp_p = cross(cell%dxtp(:,i,j),cell%dxp(:,i,j))
            ct_p2 = cross(cell%dxt (:,i,j),cell%dxp2(:,i,j))
            
            ct3_p = cross(cell%dxt3(:,i,j),cell%dxp(:,i,j))
            ct2_tp = cross(cell%dxt2(:,i,j),cell%dxtp(:,i,j))
            ct_t2p = cross(cell%dxt(:,i,j),cell%dxt2p(:,i,j))
            
            ct2p_p = cross(cell%dxt2p(:,i,j),cell%dxp(:,i,j))
            ct2_p2 = cross(cell%dxt2(:,i,j),cell%dxp2(:,i,j))
            ct_tp2 = cross(cell%dxt(:,i,j),cell%dxtp2(:,i,j))
            
            ctp2_p = cross(cell%dxtp2(:,i,j),cell%dxp(:,i,j))
            ctp_p2 = cross(cell%dxtp(:,i,j),cell%dxp2(:,i,j))
            ct_p3 = cross( cell%dxt(:,i,j),cell%dxp3(:,i,j))

!           Partials of the area element
            Jt = dot(-nk,ct2_p + ct_tp);
            Jp = dot(-nk,ctp_p + ct_p2);

!           Outward normal vector partials            
            dnt = (1/cell%J(i,j))*(ct2_p + ct_tp - Jt*(-nk))
            dnp = (1/cell%J(i,j))*(ctp_p + ct_p2 - Jp*(-nk))
            
!           Second Jacobian derivs
            Jt2 = dot(dnt,ct2_p + ct_tp)  &
                + dot(-nk, ct3_p + 2D0*ct2_tp + ct_t2p)
            Jtp = dot(dnp,ct2_p + ct_tp) &
                + dot(-nk, ct2p_p +  ct2_p2 + ct_tp2)
            Jp2 = dot(dnp,ctp_p + ct_p2) &
                + dot(-nk, ctp2_p + 2D0*ctp_p2 + ct_p3)

!           Second normal derivatives (the fact that nk is inward
!           is already taken into account)
            dnt2 = (1D0/cell%J(i,j))*(-2D0*Jt*dnt &
                 + nk*Jt2 + ct3_p + 2D0*ct2_tp + ct_t2p);
            dnp2 = (1D0/cell%J(i,j))*(-2D0*Jp*dnp &
                 + nk*Jp2 + ct_p3 + 2D0*ctp_p2 + ctp2_p);
            dntp = (1D0/cell%J(i,j))*(-Jp*dnt - Jt*dnp &
                 + nk*Jtp + ct2p_p + ct2_p2 + ct_tp2);

!           Partials of the fundamental forms
            Et  = 2D0*dot(cell%dxt2(:,i,j),cell%dxt(:,i,j))
            Ep  = 2D0*dot(cell%dxtp(:,i,j),cell%dxt(:,i,j))
            Et2 = 2D0*(dot(cell%dxt3(:,i,j),cell%dxt(:,i,j))  &
                + dot(cell%dxt2(:,i,j),cell%dxt2(:,i,j)))
            Ep2 = 2D0*(dot(cell%dxtp2(:,i,j),cell%dxt(:,i,j)) &
                + dot(cell%dxtp(:,i,j),cell%dxtp(:,i,j)))
            Etp = 2D0*(dot(cell%dxt2p(:,i,j),cell%dxt(:,i,j)) &
                + dot(cell%dxt2(:,i,j),cell%dxtp(:,i,j)))
            
            Gt  = 2D0*dot(cell%dxtp(:,i,j),cell%dxp(:,i,j))
            Gp  = 2D0*dot(cell%dxp (:,i,j),cell%dxp2(:,i,j))
            Gt2 = 2D0*(dot(cell%dxt2p(:,i,j),cell%dxp(:,i,j)) &
                + dot(cell%dxtp(:,i,j),cell%dxtp(:,i,j)))
            Gp2 = 2D0*(dot(cell%dxp3(:,i,j),cell%dxp(:,i,j))  &
                + dot(cell%dxp2(:,i,j),cell%dxp2(:,i,j)))
            Gtp = 2D0*(dot(cell%dxtp2(:,i,j),cell%dxp(:,i,j)) &
                + dot(cell%dxtp(:,i,j),cell%dxp2(:,i,j)))
            
            Ft  = dot(cell%dxt2(:,i,j),cell%dxp(:,i,j))      &
                + dot(cell%dxt(:,i,j),cell%dxtp(:,i,j))
            Fp  = dot(cell%dxtp(:,i,j),cell%dxp(:,i,j))      &
                + dot(cell%dxt(:,i,j),cell%dxp2(:,i,j))
            Ft2 = dot(cell%dxt3(:,i,j),cell%dxp(:,i,j))      &
                + 2D0*dot(cell%dxt2(:,i,j),cell%dxtp(:,i,j)) &
                + dot(cell%dxt(:,i,j),cell%dxt2p(:,i,j))
            Fp2 = dot(cell%dxtp2(:,i,j),cell%dxp(:,i,j))     &
                + 2D0*dot(cell%dxtp(:,i,j),cell%dxp2(:,i,j)) & 
                + dot(cell%dxt(:,i,j),cell%dxp3(:,i,j))
            Ftp = dot(cell%dxt2p(:,i,j),cell%dxp(:,i,j))     &
                + dot(cell%dxtp(:,i,j),cell%dxtp(:,i,j))     &
                + dot(cell%dxt2(:,i,j),cell%dxp2(:,i,j))     &
                + dot(cell%dxt(:,i,j),cell%dxtp2(:,i,j))
            
            Lt  = dot(cell%dxt3(:,i,j),-nk)      &
                + dot(cell%dxt2(:,i,j),dnt)
            Lp  = dot(cell%dxt2p(:,i,j),-nk)     &
                + dot(cell%dxt2(:,i,j),dnp)
            Lt2 = dot(cell%dxt4(:,i,j),-nk)      &
                + 2D0*dot(cell%dxt3(:,i,j),dnt)  & 
                + dot(cell%dxt2(:,i,j),dnt2)
            Lp2 = dot(cell%dxt2p2(:,i,j),-nk)    &
                + 2D0*dot(cell%dxt2p(:,i,j),dnp) &
                + dot(cell%dxt2(:,i,j),dnp2)
            Ltp = dot(cell%dxt3p(:,i,j),-nk)     &
                + dot(cell%dxt2p(:,i,j),dnt)     &
                + dot(cell%dxt3(:,i,j),dnp) + dot(cell%dxt2(:,i,j),dntp)
            
            Nt  = dot(cell%dxtp2(:,i,j),-nk) + dot(cell%dxp2(:,i,j),dnt)
            Np  = dot(cell%dxp3(:,i,j),-nk)  + dot(cell%dxp2(:,i,j),dnp)
            Nt2 = dot(cell%dxt2p2(:,i,j),-nk)+ 2D0*dot(cell%dxtp2(:,i,j),dnt) &
                + dot(cell%dxp2(:,i,j),dnt2)
            Np2 = dot(cell%dxp4(:,i,j),-nk)  + 2D0*dot(cell%dxp3(:,i,j),dnp)  &
                + dot(cell%dxp2(:,i,j),dnp2)
            Ntp = dot(cell%dxtp3(:,i,j),-nk) + dot(cell%dxtp2(:,i,j),dnp)     & 
                + dot(cell%dxp3(:,i,j),dnt) + dot(cell%dxp2(:,i,j),dntp)
            
            Mt  = dot(cell%dxt2p(:,i,j),-nk) + dot(cell%dxtp(:,i,j),dnt)
            Mp  = dot(cell%dxtp2(:,i,j),-nk) + dot(cell%dxtp(:,i,j),dnp)
            Mt2 = dot(cell%dxt3p(:,i,j),-nk) + 2D0*dot(cell%dxt2p(:,i,j),dnt) &
                + dot(cell%dxtp(:,i,j),dnt2)
            Mp2 = dot(cell%dxtp3(:,i,j),-nk) + 2D0*dot(cell%dxtp2(:,i,j),dnp) &
                + dot(cell%dxtp(:,i,j),dnp2)
            Mtp = dot(cell%dxt2p2(:,i,j),-nk)+ dot(cell%dxt2p(:,i,j),dnp)     &
                + dot(cell%dxtp2(:,i,j),dnt) + dot(cell%dxtp(:,i,j),dntp)

!           Curvature derivatives, with numerator derivatives
            Dt = Et*N + E*Nt - 2D0*Ft*M - 2D0*F*Mt + Gt*L + G*Lt
            Dp = Ep*N + E*Np - 2D0*Fp*M - 2D0*F*Mp + Gp*L + G*Lp
            
            Dt2 = Et2*N + 2D0*Et*Nt + E*Nt2 - 2D0*Ft2*M - 4D0*Ft*Mt - 2D0*F*Mt2     &
                + Gt2*L + 2D0*Gt*Lt + G*Lt2
            Dp2 = Ep2*N + 2D0*Ep*Np + E*Np2 - 2D0*Fp2*M - 4D0*Fp*Mp - 2D0*F*Mp2     &
                + Gp2*L + 2D0*Gp*Lp + G*Lp2
            Dtp = Etp*N + Et*Np + Ep*Nt + E*Ntp - 2D0*Ftp*M - 2D0*Ft*Mp - 2D0*Fp*Mt &
                - 2D0*F*Mtp + Gtp*L + Gt*Lp + Gp*Lt + G*Ltp
            
            kt = 0.5D0*Dt/(cell%J(i,j)*cell%J(i,j)) - 2D0*k*Jt/cell%J(i,j)
            kp = 0.5D0*Dp/(cell%J(i,j)*cell%J(i,j)) - 2D0*k*Jp/cell%J(i,j)
                
            kt2 = 0.5D0*Dt2/(cell%J(i,j)*cell%J(i,j)) - 2D0*Jt*(cell%J(i,j)**(-3D0))*Dt &
                + 3D0*(cell%J(i,j)**-(4D0))*Jt*Jt*D - (cell%J(i,j)**(-3D0))*Jt2*D
            kp2 = 0.5D0*Dp2/(cell%J(i,j)*cell%J(i,j)) - 2D0*Jp*(cell%J(i,j)**-3D0)*Dp &
                + 3D0*(cell%J(i,j)**(-4D0))*Jp*Jp*D - (cell%J(i,j)**(-3D0))*Jp2*D
            ktp = 0.5D0*Dtp/(cell%J(i,j)*cell%J(i,j)) - Jt*(cell%J(i,j)**(-3D0))*Dp &
                - Jp*(cell%J(i,j)**(-3D0))*Dt +3D0*(cell%J(i,j)**(-4D0))*Jt*Jp*D - (cell%J(i,j)**(-3D0))*Jtp*D

!           Metric tensor, covariant followed by contravariant (inverse)
            gv = RESHAPE(((/E,F,F,G/)), (/2,2/))
            gn = RESHAPE(((/G,-F,-F,E/)), (/2,2/))/(cell%J(i,j)*cell%J(i,j))

!           Partials of contravariant metric tensor (via properties of inverse matrices)
            dgt = -MATMUL(MATMUL(gn,RESHAPE(((/Et ,Ft ,Ft ,Gt /)), (/2,2/))), gn)
            dgp = -MATMUL(MATMUL(gn,RESHAPE(((/Ep ,Fp ,Fp ,Gp /)), (/2,2/))), gn)
            dgt2= -MATMUL(MATMUL(gn,RESHAPE(((/Et2,Ft2,Ft2,Gt2/)), (/2,2/))), gn) &
                 - 2D0*MATMUL(MATMUL(dgt,RESHAPE(((/Et ,Ft ,Ft ,Gt /)), (/2,2/))), gn)
            dgp2= -MATMUL(MATMUL(gn,RESHAPE(((/Ep2,Fp2,Fp2,Gp2/)), (/2,2/))), gn) &
                 - 2D0*MATMUL(MATMUL(dgp,RESHAPE(((/Ep ,Fp ,Fp ,Gp /)), (/2,2/))), gn)
            dgtp= -MATMUL(MATMUL(gn,RESHAPE(((/Etp,Ftp,Ftp,Gtp/)), (/2,2/))), gn) &
                 - 2D0*MATMUL(MATMUL(dgt,RESHAPE(((/Ep ,Fp ,Fp ,Gp /)), (/2,2/))), gn)

!           Contravariant basis vectors (raising index)
            c1 = cell%dxt(:,i,j)*gn(1,1) + cell%dxp(:,i,j)*gn(1,2)
            c2 = cell%dxt(:,i,j)*gn(2,1) + cell%dxp(:,i,j)*gn(2,2)

!           Derivative of contravariant basis (chain rule across raise index operation)
            c1t = cell%dxt2(:,i,j)*gn(1,1) + cell%dxtp(:,i,j)*gn(2,1) &
                + cell%dxt(:,i,j)*dgt(1,1) + cell%dxp(:,i,j)*dgt(2,1)
            c1p = cell%dxtp(:,i,j)*gn(1,1) + cell%dxp2(:,i,j)*gn(2,1) &
                + cell%dxt(:,i,j)*dgp(1,1) + cell%dxp(:,i,j)*dgp(2,1)
            c2t = cell%dxt2(:,i,j)*gn(1,2) + cell%dxtp(:,i,j)*gn(2,2) &
                + cell%dxt(:,i,j)*dgt(1,2) + cell%dxp(:,i,j)*dgt(2,2)
            c2p = cell%dxtp(:,i,j)*gn(1,2) + cell%dxp2(:,i,j)*gn(2,2) &
                + cell%dxt(:,i,j)*dgp(1,2) + cell%dxp(:,i,j)*dgp(2,2)

!           Laplace-Beltrami of mean curvature. Three chain rules with 4 things to sum each
!           First: deriv of J
            LBk = 1D0/cell%J(i,j)*(Jt*gn(1,1)*kt + Jt*gn(1,2)*kp + Jp*gn(2,1)*kt + Jp*gn(2,2)*kp) &
                + dgt(1,1)*kt + dgt(1,2)*kp + dgp(2,1)*kt + dgp(2,2)*kp &
                + gn(1,1)*kt2 + gn(1,2)*ktp + gn(2,1)*ktp + gn(2,2)*kp2

!           We can finally get the reference shapes. Cycle loop if we're getting them
            IF(.not. cell%init) THEN
                cell%c1R(:,i,j) = c1
                cell%c1R(:,i,j) = c1;
                cell%c2R(:,i,j) = c2;
                cell%kR(i,j) = k;
                cell%kdR(1,i,j) = kt;
                cell%kdR(2,i,j) = kp;
                cell%kd2R(1,i,j)= kt2;
                cell%kd2R(2,i,j)= kp2;
                cell%kd2R(3,i,j)= ktp;
                
                cell%c1tR(:,i,j) = c1t;
                cell%c1pR(:,i,j) = c1p;
                cell%c2tR(:,i,j) = c2t;
                cell%c2pR(:,i,j) = c2p;
                
                cell%JR(i,j) = cell%J(i,j)

                CYCLE inner
            ENDIF

!           Christoffel symbols
            c111 = dot(cell%dxt2(:,i,j),c1)
            c112 = dot(cell%dxtp(:,i,j),c1)
            c122 = dot(cell%dxp2(:,i,j),c1)
            c222 = dot(cell%dxp2(:,i,j),c2)
            c221 = dot(cell%dxtp(:,i,j),c2)
            c211 = dot(cell%dxt2(:,i,j),c2)
            
!           Christoffel derivatives: need 5 each
            c111t = dot(cell%dxt3(:,i,j) ,c1) + dot(cell%dxt2(:,i,j),c1t)
            c112t = dot(cell%dxt2p(:,i,j),c1) + dot(cell%dxtp(:,i,j),c1t)
            c122t = dot(cell%dxtp2(:,i,j),c1) + dot(cell%dxp2(:,i,j),c1t)
            c221t = dot(cell%dxt2p(:,i,j),c2) + dot(cell%dxtp(:,i,j),c2t)
            c222t = dot(cell%dxtp2(:,i,j),c2) + dot(cell%dxp2(:,i,j),c2t)
            
            c111p = dot(cell%dxt2p(:,i,j),c1) + dot(cell%dxt2(:,i,j),c1p)
            c112p = dot(cell%dxtp2(:,i,j),c1) + dot(cell%dxtp(:,i,j),c1p)
            c211p = dot(cell%dxt2p(:,i,j),c2) + dot(cell%dxt2(:,i,j),c2p)
            c221p = dot(cell%dxtp2(:,i,j),c2) + dot(cell%dxtp(:,i,j),c2p)
            c222p = dot(cell%dxp3(:,i,j) ,c2) + dot(cell%dxp2(:,i,j),c2p)

!           Surface projection operator
            Prj = 0D0
            FORALL(q = 1:3) Prj(q,q) = 1D0
            Prj = Prj - outer(nk,nk)

!           Now finally some mechanics
!           New Bending: Helfrich
!           Normal component (commented part assumes homogeneous spontaneous curvature)
            fb = Eb*(2D0*LBk + (2D0*k + c0)*(2D0*k*k - 2D0*kG  - c0*k)) ! + Eb*cell%LBc0(i,j)
            Wb(i,j) = Eb/2D0*(2D0*k-c0)*(2D0*k-c0)

!           Tangential component (in Cartesian). Just surface gradient of curvature.
!           Start w/ just gradient, Prjt is just temp storage
!           Most people ignore, really doesn't do  much
            fbt(1) = 4D0*k*(kt*c1(1) + kp*c2(1))
            fbt(2) = 4D0*k*(kt*c1(2) + kp*c2(2))
            fbt(3) = 4D0*k*(kt*c1(3) + kp*c2(3))
            Prjt = TRANSPOSE(Prj)
            fbt = INNER3_33(fbt, Prjt)
            
!           In plane tension and derivatives
!           Deformation gradient tensor                  
            Fd = outer(cell%dxt(:,i,j), cell%c1R(:,i,j)) + outer(cell%dxp(:,i,j), cell%c2R(:,i,j))

!           Left Cauchy-Green
            V2 = MATMUL(Fd, TRANSPOSE(Fd))

!           Principal strains/directions
            tV2 = V2
            CALL dsyev('v', 'u', 3, V2, 3, lams, wrk, 102, rc)
            ev = V2
            V2 = tV2
            
            es(1) = sqrt(lams(3))
            es(2) = sqrt(lams(2))

!           Eigenvectors I actually need
            ev1 = ev(:,3)
            ev2 = ev(:,2)

!           Strain invariants
            I1 = es(1)*es(1) + es(2)*es(2) - 2D0
            I2 = es(1)*es(1)*es(2)*es(2) - 1D0

            Ws(i,j) = 0.5D0*B*(I1*I1 + 2D0*I2*I2 - 2D0*I2 + C*I2*I2/B)

!           In plane tension (Cartesian)
            tau = 0.5D0*(B/(es(1)*es(2))*(I1 + 1)*V2 &
                + es(1)*es(2)*(C*I2 - B)*Prj)
                
!           Tau in surface coordinates
            tauab(1,1) = DOT(INNER3_33(c1,tau), c1)
            tauab(1,2) = DOT(INNER3_33(c1,tau), c2)
            tauab(2,1) = DOT(INNER3_33(c2,tau), c1)
            tauab(2,2) = DOT(INNER3_33(c2,tau), c2)

!           We need the derivatives of tau with respect to phi/theta. It isn't really written
!           in terms of these, though, so we need several partials in order to do chain rule
            Fdt = outer(cell%dxt2(:,i,j), cell%c1R(:,i,j)) + outer(cell%dxtp(:,i,j), cell%c2R(:,i,j)) &
                + outer(cell%dxt(:,i,j), cell%c1tR(:,i,j)) + outer(cell%dxp(:,i,j), cell%c2tR(:,i,j))
            Fdp = outer(cell%dxtp(:,i,j), cell%c1R(:,i,j)) + outer(cell%dxp2(:,i,j), cell%c2R(:,i,j)) &
                + outer(cell%dxt(:,i,j), cell%c1pR(:,i,j)) + outer(cell%dxp(:,i,j), cell%c2pR(:,i,j))

            V2t = MATMUL(Fdt, TRANSPOSE(Fd)) + MATMUL(Fd, TRANSPOSE(Fdt))
            V2p = MATMUL(Fdp, TRANSPOSE(Fd)) + MATMUL(Fd, TRANSPOSE(Fdp))

!           Eigenvalues w.r.t. their matrix
            lam1dV2 = outer(ev1, ev1)
            lam2dV2 = outer(ev2, ev2)

!           Principal strains to angles.
            es1t = 0.5D0/es(1)*SUM(lam1dV2*V2t)
            es1p = 0.5D0/es(1)*SUM(lam1dV2*V2p)
            es2t = 0.5D0/es(2)*SUM(lam2dV2*V2t)
            es2p = 0.5D0/es(2)*SUM(lam2dV2*V2p)

!           Strain invariant derivatives
            I1t = 2D0*es(1)*es1t + 2D0*es(2)*es2t
            I1p = 2D0*es(1)*es1p + 2D0*es(2)*es2p
            I2t = (2D0*es(1)*es(2)*es(2))*es1t + (es(1)*es(1)*2D0*es(2))*es2t
            I2p = (2D0*es(1)*es(2)*es(2))*es1p + (es(1)*es(1)*2D0*es(2))*es2p

!           Projection operator deriv (double neg in 2nd part)
            Prjt = outer(-dnt, -nk) + outer(nk, dnt)
            Prjp = outer(-dnp, -nk) + outer(nk, dnp)

!           And finally the big one: derivatives of tau, separated via chain rule
!           starting with es(1), es(2)
            taut = (-B/(2D0*es(1)*es(1)*es(2))*(I1 + 1D0)*V2 + 0.5D0*es(2)*(C*I2 - B)*Prj)*es1t & 
                 + (-B/(2D0*es(1)*es(2)*es(2))*(I1 + 1D0)*V2 + 0.5D0*es(1)*(C*I2 - B)*Prj)*es2t &
                 + 0.5D0*B/(es(1)*es(2))*V2*I1t & ! Strain invariants
                 + 0.5D0*es(1)*es(2)*C*Prj*I2t &
                 + 0.5D0*B/(es(1)*es(2))*(I1 + 1D0)*V2t & ! Tensors
                 + 0.5D0*es(1)*es(2)*(I2 - B)*Prjt

            taup = (-B/(2D0*es(1)*es(1)*es(2))*(I1 + 1D0)*V2 + 0.5D0*es(2)*(C*I2 - B)*Prj)*es1p & 
                 + (-B/(2D0*es(1)*es(2)*es(2))*(I1 + 1D0)*V2 + 0.5D0*es(1)*(C*I2 - B)*Prj)*es2p &
                 + 0.5D0*B/(es(1)*es(2))*V2*I1p & ! Strain invariants
                 + 0.5D0*es(1)*es(2)*C*Prj*I2p &
                 + 0.5D0*B/(es(1)*es(2))*(I1 + 1D0)*V2p & ! Tensors
                 + 0.5D0*es(1)*es(2)*(I2 - B)*Prjp
                 
!           Put into local surface coordinates, keeping only needed derivs
            dtauab(1,1) = DOT(INNER3_33(c1t,tau), c1) + DOT(INNER3_33(c1,taut), c1) + DOT(INNER3_33(c1,tau), c1t)
            dtauab(1,2) = DOT(INNER3_33(c1t,tau), c2) + DOT(INNER3_33(c1,taut), c2) + DOT(INNER3_33(c1,tau), c2t)
            dtauab(2,1) = DOT(INNER3_33(c2p,tau), c1) + DOT(INNER3_33(c2,taup), c1) + DOT(INNER3_33(c2,tau), c1p)
            dtauab(2,2) = DOT(INNER3_33(c2p,tau), c2) + DOT(INNER3_33(c2,taup), c2) + DOT(INNER3_33(c2,tau), c2p)
                 
!           Covariant curvature tensor
            bv(1,:) = (/L,M/)
            bv(2,:) = (/M,N/)

!           Mixed curvature tensor (raise the column)
            bm(1,1) = bv(1,1)*gn(1,1) + bv(1,2)*gn(2,1)
            bm(1,2) = bv(1,1)*gn(1,2) + bv(1,2)*gn(2,2)
            bm(2,1) = bv(2,1)*gn(1,1) + bv(2,2)*gn(2,1)
            bm(2,2) = bv(2,1)*gn(1,2) + bv(2,2)*gn(2,2)

!           Covariant divergences of tau, in theta and phi
            cvt = dtauab(1,1) + dtauab(2,1) + tauab(1,1)*(2D0*c111 + c221) &
                + tauab(2,1)*(2D0*c112 + c222) + tauab(1,2)*c112 + tauab(2,2)*c122
            cvp = dtauab(1,2) + dtauab(2,2) + tauab(1,2)*(c111 + 2D0*c221) &
                + tauab(2,2)*(c112 + 2D0*c222) + tauab(1,1)*c211 + tauab(2,1)*c221


!           Forces in surface coordinates
            cell%fab(1, i, j) = -cvt !   + bm(1,1)*q1 + bm(2,1)*q2
            cell%fab(2, i, j) = -cvp !   + bm(1,2)*q1 + bm(2,2)*q2
            cell%fab(3, i, j) =     - tauab(1,1)*bv(1,1) - tauab(1,2)*bv(1,2) &
                                    - tauab(2,1)*bv(2,1) - tauab(2,2)*bv(2,2) &
                                    + fb - cell%int_pres
                                    !  - cq

!           Forces in Cartesian, on fine grid
            cell%ff(:, i, j) = cell%fab(1, i, j)*cell%dxt(:,i,j) &
                              + cell%fab(2, i, j)*cell%dxp(:,i,j) &
                              + cell%fab(3, i, j)*(-nk) &
                              + 0D0 !fbt
            
!           This is pretty lazy, but I need to save the normal vector and not the
!           surface forces, when it was opposite at one time, so I just swap them
            cell%fab(:,i,j) = -nk

        ENDDO inner
    ENDDO

    ! cvt = 0D0

!     DO i = 1,cell%Yf%nt
!         DO j = 1,cell%Yf%np
! !           Integrate via Gauss quad
!             cvt = cvt + Ws(i,j)*cell%Yf%wg(i)*cell%JR(i,j)*cell%Yf%dphi/sin(cell%Yf%tht(i))
!         ENDDO
!     ENDDO

!     print *, cvt, cell%intg(Wb), cvt + cell%intg(Wb), cell%intg(Ws)

!   Now we need to filter for anti-aliasing. Do the transform of the force into spectral
!   space, cut the highest modes, and transform back into physical space
    IF(cell%init) THEN
        cell%fmn(1,:) = cell%Yf%forward(cell%ff(1,:,:)*cell%J/SIN(cell%Yf%th), cell%q)
        cell%fmn(2,:) = cell%Yf%forward(cell%ff(2,:,:)*cell%J/SIN(cell%Yf%th), cell%q)
        cell%fmn(3,:) = cell%Yf%forward(cell%ff(3,:,:)*cell%J/SIN(cell%Yf%th), cell%q)

!       Normal vector for volume correction
        cell%nkmn(1,:) = cell%Yf%forward(cell%fab(1,:,:), cell%q) 
        cell%nkmn(2,:) = cell%Yf%forward(cell%fab(2,:,:), cell%q)
        cell%nkmn(3,:) = cell%Yf%forward(cell%fab(3,:,:), cell%q)

!       Normal and area for fluid
        cell%nkt(1,:) = cell%Yf%forward(cell%fab(1,:,:)*cell%J/SIN(cell%Yf%th), cell%q) 
        cell%nkt(2,:) = cell%Yf%forward(cell%fab(2,:,:)*cell%J/SIN(cell%Yf%th), cell%q)
        cell%nkt(3,:) = cell%Yf%forward(cell%fab(3,:,:)*cell%J/SIN(cell%Yf%th), cell%q)
    ENDIF
END SUBROUTINE Stresscell

! -------------------------------------------------------------------------!
! Knowing the force jump get the velocity on the surface of the cell via
! the fluid problem
SUBROUTINE Fluidcell(cell, prob, A2, b2, celli)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell
    TYPE(probType), INTENT(IN), TARGET :: prob
    COMPLEX(KIND = 8), INTENT(OUT), ALLOCATABLE :: A2(:,:), b2(:)
    TYPE(cellType), INTENT(IN), POINTER, OPTIONAL :: celli

    INTEGER :: ip, ic, i, j, i2, j2, n, m, it, im, row, col, im2, n2, m2, &
               colm, ind, im3, nt ,np, indi = 1, indj = 1
    LOGICAL sing
    COMPLEX(KIND = 8) :: At(3,3), bt(3), tmpsum(3,3)
    REAL(KIND = 8) :: Uc(3), xcr(3), Utmp(3,3), r(3), eye(3,3), minr, dphi
    COMPLEX(KIND = 8), ALLOCATABLE :: b(:), Ci(:,:,:,:), Ei(:,:,:,:), & 
                                      Dr(:,:,:), Fi(:,:,:,:), Ai(:,:,:,:,:,:), &
                                      fmnR(:,:), xmnR(:,:), nmnR(:,:)
    REAL(KIND = 8), ALLOCATABLE :: frot(:,:,:), xcg(:,:,:), nJt(:,:,:), &
                                   ft(:),ft2(:,:), Bi(:,:,:,:), wgi(:)
    COMPLEX(KIND = 8), POINTER :: vcurn(:,:), es(:,:)
    REAL(KIND = 8), POINTER :: cPmn(:,:, :)
    TYPE(YType), POINTER :: Y
    TYPE(nmType), POINTER :: nm
    INTEGER(KIND = 8) :: tic, toc, rate

    Y => cell%Y

!   Allocate things
    ALLOCATE(A2(prob%Nmat, prob%Nmat), &
             b(3*Y%nt*Y%np), &
             b2(prob%Nmat), &
             frot(3, Y%nt, Y%np), &
             xcg(3, Y%nt, Y%np), &
             nJt(3, Y%nt, Y%np), &
             ft(3), ft2(3,3), &
             Bi(3,3,Y%nt,Y%np), &
             Ci(3,3, 2*(Y%p-1)+1, Y%nt), &
             Ei(3,3, 2*(Y%p-1)+1, Y%p), &
             Fi(3,3, 2*(Y%p-1)+1, Y%nt), &
             Ai(3,3, 2*(Y%p-1)+1, Y%p, Y%nt, Y%np), &
             Dr(3,3,Y%p*Y%p),  &
             es(2*(Y%p-1)+1, Y%np), &
             cPmn(Y%p, 2*(Y%p-1)+1, Y%nt), &
             fmnR(3, (cell%p+1)*(cell%p+1)), &
             nmnR(3, (cell%p+1)*(cell%p+1)), &
             xmnR(3, (cell%p+1)*(cell%p+1)), &
             wgi(Y%nt))

    eye = 0D0
    Ai = 0D0
    FORALL(j = 1:3) eye(j,j) = 1D0
!   We need to do two integral loops: the first calculates the integral 
!   that is the convolution of the Stokeslets and stresslets. We need those
!   values at the integration points for the next step. Next we do the
!   Galerkin integral over a sphere, using the inner integrals calculated
!   in the previous step

!   Doing the Galerkin projection of the single/double layer integrals could require
!   O(p^5) operations if we did it all at once (it's a p^4 size matrix), the two
!   integrals require p^2 sized sums each. However, some of these sums have components
!   that are independent of each other, so we can do them 1 at a time and store them
!   in temporary arrays, leading down to O(p^5). The 1st big matrix is the components of
!   the spherical harmonic functions evaluated at the Gauss points (rows=>GP,
!   cols=>Sph fns).

!!! Net force calculation - should be zero with no external force
    ! Uc = 0D0
    ! DO i = 1,cell%Yf%nt
    !     DO j = 1,cell%Yf%np
    !         Uc = Uc + cell%Yf%wg(i)*cell%J(i,j)*cell%Yf%dphi*cell%ff(:,i,j)
    !     ENDDO
    ! ENDDO
    ! print *, Uc
    ! Uc = 0D0

    ip = 0
    ic = 0
    CALL SYSTEM_CLOCK(tic,rate)
!   First loops: singular integral points
    DO i = 1,Y%nt
        DO j = 1,Y%np
!           Bookkeeping
            ic = ic + 1
            row = 3*(ic - 1)  + 1

!           Location of north pole in unrotated frame (just current integration point)
            xcr = cell%x(:,i,j)

            sing = .false.
!           If the integration and target surfaces are different, check minimum spacing
            IF(PRESENT(celli)) THEN
                minr = celli%h + 1D0
                DO i2 = 1,celli%Yf%nt
                    DO j2 = 1,celli%Yf%np
                        r = celli%xf(:,i2,j2) - xcr
                        minr = MIN(norm2(r), minr)
!                       Save indices of min spacing
                        IF(minr .eq. norm2(r)) THEN
                            indi = i2
                            indj = j2
                        ENDIF
                    ENDDO
                ENDDO

!               If min spacing is small, we need to do near-singular integration
                IF(minr .lt. celli%h) THEN
                    sing = .true.

!                   Need to integrate on finer grid
                    nt = cell%Yf%nt + cell%Y%nt
                    np = cell%Yf%np + cell%Y%np

!                   Deallocate integ. quants
                    DEALLOCATE(frot, xcg, nJt, Bi, Ci, wgi)
                    ALLOCATE( &
                    frot(3, nt, np), &
                    xcg(3,  nt, np), &
                    nJt(3,  nt, np), &
                    Bi(3,3, nt, np), &
                    Ci(3,3, 2*(Y%p-1)+1, nt), &
                    wgi(nt))

                    es   => prob%esf
                    cPmn => prob%cPmnf
                    dphi = celli%Y%dphi

!                   Manage the additional prefactors stemming from the integrals
                    wgi(1:Y%nt)  = celli%Y%wg*(pi - prob%thtc)/2D0
                    wgi(Y%nt + 1:Y%nt + celli%Yf%nt)  = celli%Yf%wg*prob%h*-prob%k*COSH(prob%k*prob%xsf - prob%k)

!                   We don't do cosine transformation to cluster points near near-singularity, mult sine back in
                    wgi = wgi*SIN(celli%Ys%th(:,1))

!                   The below formulation is slightly inefficient. To remain general, I want to just have a single
!                   grid. However, the singular integral is calculated on 2 grids, one fine and one coarse.
!                   I put this in one grid and there is some overlap, so that there are points that aren't used.
!                   When calculating the integrals, I just cycle past these points.

!                   Rotate about nearest point to projected singularity
                    xmnR(1,:) = celli%Yf%rotate(celli%xmn(1,:), indi, indj, -celli%Yf%phi(indj))
                    xmnR(2,:) = celli%Yf%rotate(celli%xmn(2,:), indi, indj, -celli%Yf%phi(indj))
                    xmnR(3,:) = celli%Yf%rotate(celli%xmn(3,:), indi, indj, -celli%Yf%phi(indj))

                    xcg(1,:,:) = celli%Ys%backward(xmnR(1,:))
                    xcg(2,:,:) = celli%Ys%backward(xmnR(2,:))
                    xcg(3,:,:) = celli%Ys%backward(xmnR(3,:))

!                   Forces on rotated grid
                    fmnR(1,:) = celli%Yf%rotate(celli%fmn(1,1:(celli%p+1)*(celli%p+1)), indi, indj, -celli%Yf%phi(indj))
                    fmnR(2,:) = celli%Yf%rotate(celli%fmn(2,1:(celli%p+1)*(celli%p+1)), indi, indj, -celli%Yf%phi(indj))
                    fmnR(3,:) = celli%Yf%rotate(celli%fmn(3,1:(celli%p+1)*(celli%p+1)), indi, indj, -celli%Yf%phi(indj))

                    frot(1,:,:) = celli%Ys%backward(fmnR(1,:), celli%p)
                    frot(2,:,:) = celli%Ys%backward(fmnR(2,:), celli%p)
                    frot(3,:,:) = celli%Ys%backward(fmnR(3,:), celli%p)

!                   Rotate the normal vector total constants
                    nmnR(1,:) = celli%Yf%rotate(celli%nkt(1,1:(celli%p+1)*(celli%p+1)), indi, indj, -celli%Yf%phi(indj))
                    nmnR(2,:) = celli%Yf%rotate(celli%nkt(2,1:(celli%p+1)*(celli%p+1)), indi, indj, -celli%Yf%phi(indj))
                    nmnR(3,:) = celli%Yf%rotate(celli%nkt(3,1:(celli%p+1)*(celli%p+1)), indi, indj, -celli%Yf%phi(indj))

                    nJt(1,:,:) = celli%Ys%backward(nmnR(1,:), celli%p)
                    nJt(2,:,:) = celli%Ys%backward(nmnR(2,:), celli%p)
                    nJt(3,:,:) = celli%Ys%backward(nmnR(3,:), celli%p)
!               Well-separated, normal grid/integration
                ELSE
                    sing = .false.

!                   We can use the coarse grid
                    nt = Y%nt
                    np = Y%np

!                   Deallocate integ. quants
                    DEALLOCATE(frot, xcg, nJt, Bi, Ci, wgi)
                    ALLOCATE( &
                    frot(3, nt, np), &
                    xcg(3,  nt, np), &
                    nJt(3,  nt, np), &
                    Bi(3,3, nt, np), &
                    Ci(3,3, 2*(Y%p-1)+1, nt), &
                    wgi(nt))
                    
                    es   => prob%es
                    cPmn => prob%cPmn
                    dphi = celli%Y%dphi
                    wgi  = celli%Y%wg

                    xcg(1,:,:) = celli%x(1,:,:)
                    xcg(2,:,:) = celli%x(2,:,:)
                    xcg(3,:,:) = celli%x(3,:,:)

                    fmnR(1,:) = celli%fmn(1,1:(cell%p+1)*(cell%p+1))
                    fmnR(2,:) = celli%fmn(2,1:(cell%p+1)*(cell%p+1))
                    fmnR(3,:) = celli%fmn(3,1:(cell%p+1)*(cell%p+1))

                    frot(1,:,:) = Y%backward(fmnR(1,:), cell%p)
                    frot(2,:,:) = Y%backward(fmnR(2,:), cell%p)
                    frot(3,:,:) = Y%backward(fmnR(3,:), cell%p)
    
                    nmnR(1,:) = celli%nkt(1,1:(cell%p+1)*(cell%p+1))
                    nmnR(2,:) = celli%nkt(2,1:(cell%p+1)*(cell%p+1))
                    nmnR(3,:) = celli%nkt(3,1:(cell%p+1)*(cell%p+1))

                    nJt(1,:,:) = Y%backward(nmnR(1,:), cell%p)
                    nJt(2,:,:) = Y%backward(nmnR(2,:), cell%p)
                    nJt(3,:,:) = Y%backward(nmnR(3,:), cell%p)
                ENDIF

!           Fully singular integration on same cell: Get rotated constants
            ELSE
                sing = .false.

!               Velocity at integration point
                Utmp = TRANSPOSE(prob%dU)
                Uc = INNER3_33(cell%x(:,i,j), Utmp)

!               For integration (changes if need a finer grid)
                dphi = Y%dphi

!               Exponential  part
                es   => prob%es
!               Legendre polynomial part
                cPmn => prob%cPmn

                nt = Y%nt
                np = Y%np

                xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
                xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
                xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

!               Rotated integration points in unrotated frame
                xcg(1,:,:) = Y%backward(xmnR(1,:))
                xcg(2,:,:) = Y%backward(xmnR(2,:))
                xcg(3,:,:) = Y%backward(xmnR(3,:))

!               Forces on rotated grid
                fmnR(1,:) = Y%rotate(cell%fmn(1,1:(cell%p+1)*(cell%p+1)), i, j, -Y%phi(j))
                fmnR(2,:) = Y%rotate(cell%fmn(2,1:(cell%p+1)*(cell%p+1)), i, j, -Y%phi(j))
                fmnR(3,:) = Y%rotate(cell%fmn(3,1:(cell%p+1)*(cell%p+1)), i, j, -Y%phi(j))

                frot(1,:,:) = Y%backward(fmnR(1,:), cell%p)
                frot(2,:,:) = Y%backward(fmnR(2,:), cell%p)
                frot(3,:,:) = Y%backward(fmnR(3,:), cell%p)

!               Rotate the normal vector total constants
                nmnR(1,:) = Y%rotate(cell%nkt(1,1:(cell%p+1)*(cell%p+1)), i, j, -Y%phi(j))
                nmnR(2,:) = Y%rotate(cell%nkt(2,1:(cell%p+1)*(cell%p+1)), i, j, -Y%phi(j))
                nmnR(3,:) = Y%rotate(cell%nkt(3,1:(cell%p+1)*(cell%p+1)), i, j, -Y%phi(j))

                nJt(1,:,:) = Y%backward(nmnR(1,:), cell%p)
                nJt(2,:,:) = Y%backward(nmnR(2,:), cell%p)
                nJt(3,:,:) = Y%backward(nmnR(3,:), cell%p)

                wgi = Y%ws
            ENDIF

!           First matrix: integral components at the i,jth - i2,j2 grid.
            Bi = 0D0
            bt = 0D0

            DO j2 = 1,np
                DO i2 = 1,nt
!                   Need an exception for near-singular int for if we're in mis-matched grids
                    IF(sing .and. ((i2.gt.Y%nt .and. j2.le.Y%np) .or. (i2.le.Y%nt .and. j2.gt.Y%np))) CYCLE
                    IF(sing .and. j2 .eq. Y%np + 1) dphi = celli%Yf%dphi

                    r = xcg(:,i2,j2) - xcr
!                   Matrix of integration grid about i,j-th point 
                    Bi(1:3,1:3,i2,j2) = Tij(r, nJt(:,i2,j2))
                    
!                   RHS vector
                    ft = frot(:,i2,j2)
                    ft2 = Gij(r, eye)
                    bt = bt + INNER3_33(ft,ft2)*wgi(i2)*dphi
                ENDDO
            ENDDO

            b(row:row+2) = bt/(1D0 + cell%lam)
            IF(.not. PRESENT(celli)) b(row:row+2) = b(row:row+2) &
                                                  - Uc*8D0*pi/(1D0 + cell%lam)

!           Next intermediate matrices: over phi's and theta's
            im2 = 0
            DO m2 = -(Y%p-1), Y%p-1
                im2 = im2 + 1
                DO i2 = 1, nt
                    tmpsum = 0D0
                    DO j2 = 1, np
                        IF(sing .and. ((i2.gt.Y%nt .and. j2.le.Y%np) .or. (i2.le.Y%nt .and. j2.gt.Y%np))) CYCLE
                        tmpsum = tmpsum + Bi(1:3, 1:3, i2, j2)*es(im2, j2) ! dphis incorporated into es
                    ENDDO
                    Ci(1:3,1:3,im2,i2) = tmpsum
                ENDDO
            ENDDO

            DO n = 0, Y%p-1
                ind = n+1
                im2 = Y%p-1
                DO m2 = 0,(Y%p-1)
                    im2 = im2+1
                    colm = im2 - 2*m2
                    tmpsum = 0D0
                    DO i2 = 1,nt
                        tmpsum = tmpsum + Ci(1:3,1:3, im2, i2)*cPmn(ind,im2,i2)*wgi(i2)
                    ENDDO
                    Ei(1:3,1:3, im2, ind) = tmpsum
!                   Symmetry
                    IF(m2.gt.0) THEN
                        Ei(1:3,1:3, colm, ind) = CONJG(tmpsum)*(-1D0)**m2
                    ENDIF
                ENDDO
            ENDDO

!           Last loop to bring it all together and get the row
            it = 0
            Dr = 0D0
            DO n = 0,Y%p-1
                ind = n + 1
                im = 0
                DO m = -n,n
                    im = im + 1
                    it = it + 1
                    tmpsum = 0D0
                    im3 = 0
!                   If this is a self-integral, we need to rotate.
!                   If not, and not near-singular, just do normal SpHarms
                    IF(PRESENT(celli) .and. .not. sing) THEN
                        im2 = im + Y%p - n - 1
                        tmpsum = Ei(1:3,1:3, im2, ind)
                    ELSE
!                       The near-singular integral has a rotation, need to account for                        
                        IF(sing) THEN
                        DO m2 = -n,n
                            im2 = Y%p + im3 - n
                            im3 = im3 + 1
                            tmpsum = tmpsum &
                                + Ei(1:3,1:3, im2, ind) &
                                * prob%Yf%rot(indi,indj,ind)%dmms(im,im3) &
                                * EXP(ii*(m-m2)*prob%Yf%phi(indj))
                        ENDDO
                        ELSE
                        DO m2 = -n,n
                            im2 = Y%p + im3 - n
                            im3 = im3 + 1
                            tmpsum = tmpsum &
                                + Ei(1:3,1:3, im2, ind) &
                                * Y%rot(i,j,ind)%dmms(im,im3) &
                                * EXP(ii*(m-m2)*Y%phi(j))
                        ENDDO
                        ENDIF
                    ENDIF
                    Dr(1:3,1:3,it) = tmpsum
                ENDDO
            ENDDO

!           Now let's put this in the matrix
            it = 0
            ind = 0
            DO n = 0,Y%p-1
                ind = ind+1
                nm  => Y%nm(n+1)
                im = 0
                DO m = -n,n
                    im = im + 1
                    im2 = m + (Y%p)
                    vcurn => nm%v(im,:,:)
                    it = it+1
                    col = 3*(it-1) + 1
                    Ai(1:3,1:3,im2, ind, i, j) = Dr(1:3,1:3, it)*(1D0-cell%lam)/(1D0+cell%lam)
                    IF(.not. PRESENT(celli)) Ai(1:3,1:3, im2, ind, i, j) = &
                                             Ai(1:3,1:3, im2, ind, i, j) - vcurn(i,j)*4D0*pi*eye
                ENDDO
            ENDDO

        ENDDO
    ENDDO

!   Second loop is over just the normal grid, redo pre-allocated parts for this grid
    es   => prob%es
    cPmn => prob%cPmn

!   Second integral: The outer loops go over the order and degree of the previous integrals
    it = 0
    DO n = 0,Y%p - 1
        nm => cell%Y%nm(n+1)
        ! it = it + n
        DO m = -n,n
            im = m + Y%p
            it = it + 1
            col = 3*it - 2
            ! colm= col - 2*3*m

!           First loop: m2 (Galerkin order), theta, sum phi, 
            DO m2 = -(Y%p-1), Y%p-1
                im2 = m2 + Y%p
                DO i2 = 1,Y%nt
                    tmpsum = 0D0
                    DO j2 = 1,Y%np
                        tmpsum = tmpsum + Ai(1:3,1:3, im, n+1, i2, j2)*CONJG(es(im2,j2))
                    ENDDO
                    Fi(1:3,1:3,im2,i2) = tmpsum
                ENDDO
            ENDDO

!           Second loop: n2, m2, sum theta
            im2 = 0
            DO n2 = 0, Y%p-1
                DO m2 = -n2, n2
                    im3 = m2 + Y%p
                    im2 = im2 + 1
                    row = 3*im2 - 2
                    ! rowm= row - 2*3*m2
                    At = 0D0
                    ! At2 = 0D0
                    DO i2 = 1,Y%nt
                        At  = At  + Fi(1:3,1:3,im3, i2)*cPmn(n2+1,im3,i2)*cell%Y%wg(i2)
                        ! At2 = At2 + CONJG(Fi(1:3,1:3,im3, i2))*cPmn(n2+1,im3,i2)*cell%Y%wg(i2)*(-1D0)**m2
                    ENDDO
                    A2(row:row+2,col:col+2) = At

!                   Can't get symmetry working as before, because I loop over cols, then rows now...
!                   Exploit symmetry (This and the calculation of At2 are a little overkill,
!                   but it's finnicky to get the right if statements so I'll just leave them)
                    ! A2(rowm:rowm+2, col :col +2) = At2
                    ! A2(row :row +2, colm:colm+2) = (-1D0)**(m - m2)*CONJG(At2)
                    ! A2(rowm:rowm+2, colm:colm+2) = (-1D0)**(m + m2)*CONJG(At)
                ENDDO
            ENDDO

            ic = 0
!           Loop over integration points to calc integral
            bt = 0D0
            vcurn => nm%v(m + n + 1,:,:)
            DO i =1,Y%nt
                DO j = 1,Y%np
                ic = ic+1

!               Intg. b (essentially forward transform of RHS!)
                bt = bt + b(3*ic-2:3*ic)*CONJG(vcurn(i,j)) &
                    *cell%Y%wg(i)*cell%Y%dphi
                ENDDO
            ENDDO
            b2(col:col+2) = bt

        ENDDO
    ENDDO
    CALL SYSTEM_CLOCK(toc)
    ! print *, REAL(toc-tic)/REAL(rate)
END SUBROUTINE Fluidcell

! -------------------------------------------------------------------------!
! Problem advancement
SUBROUTINE UpdateProb(prob, ord, reduce)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell, celli
    INTEGER, INTENT(IN) :: ord
    LOGICAL, INTENT(IN) :: reduce
    REAL(KIND = 8) :: zm
    COMPLEX(KIND = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:), ut(:), wrk(:), &
                                      A2(:,:), b2(:), A(:,:), b(:)
    INTEGER :: ic, ic2, i, row, col, iter, info
    REAL(KIND = 8), ALLOCATABLE :: rwrk(:), utr(:), uti(:)
    COMPLEX, ALLOCATABLE :: swrk(:)
    COMPLEX(KIND = 8), POINTER :: xmn(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER(KIND = 8) tic, toc, rate
    
    ALLOCATE(ut(prob%NmatT), &
             uti(prob%NmatT), &
             utr(prob%NmatT), &
             IPIV(prob%NmatT), wrk(prob%NmatT), &
             swrk(prob%NmatT*(prob%NmatT+1)), &
             rwrk(prob%NmatT), &
             A(prob%NmatT, prob%NmatT), &
             b(prob%NmatT))
    
    CALL SYSTEM_CLOCK(tic,rate)

    A = 0D0
    b = 0D0
    b2 = 0D0
    A2 = 0D0
!   Construct the matrix
!   First loop: velocity surface
    DO ic = prob%PCells(1), prob%PCells(2)
        row = (ic -1)*prob%Nmat + 1
        cell => prob%cell(ic)

!       Second loop: integral surface
        DO ic2 = 1,prob%NCell
            celli => prob%cell(ic2)
            
!           Get geometric info about new state, get stress in new state if first go round
            IF(ic.eq.prob%PCells(1)) THEN
                CALL celli%derivs()
                CALL celli%stress()
!               Characteristic grid spacing
                celli%h = SQRT(celli%SA()/celli%Y%nt**2D0)
            ENDIF

            col = (ic2-1)*prob%Nmat + 1
!           Get velocity sub-matrix for cell-cell combo (Same or diff cell)
            IF(ic .eq. ic2) THEN
                CALL cell%fluid(prob, A2, b2)
            ELSE
                CALL cell%fluid(prob, A2, b2, celli)
            ENDIF

!           Put sub matrix into big matrix
            A(row:row + prob%Nmat - 1, col:col + prob%Nmat - 1) = A2
!           Sum over all the integrals
            b(row:row + prob%Nmat - 1) = b(row:row + prob%Nmat - 1) + b2
        ENDDO
    ENDDO

    DO i = 1,prob%NmatT
        A(:,i) = prob%cm%reduce(REAL(A(:,i))) + prob%cm%reduce(AIMAG(A(:,i)))*ii
    ENDDO
    b = prob%cm%reduce(REAL(b)) + prob%cm%reduce(AIMAG(b))*ii

!   Invert big matrix to get a list of all the vel constants of all cells
    ut = 0D0
    IF(prob%cm%mas()) THEN
        CALL zcgesv(prob%NmatT, 1, A, prob%NmatT, IPIV, b, prob%NmatT, &
                    ut, prob%NmatT, wrk, swrk, rwrk, iter, info)
    ENDIF

!   Broadcast to all processors
    IF(prob%cm%np() .gt. 1) THEN
        utr = REAL(ut)
        uti = AIMAG(ut)
        CALL prob%cm%bcast(utr)
        CALL prob%cm%bcast(uti)
        ut = utr + uti*ii
    ENDIF

!   Advance in time now
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)

!       Reconstruct individual vels
        row = (ic-1)*prob%Nmat + 1
        cell%umn = 0D0
        cell%umn(1,1:prob%Nmat/3) = ut((/(i, i=row    , row + prob%Nmat-2, 3)/))
        cell%umn(2,1:prob%Nmat/3) = ut((/(i, i=row + 1, row + prob%Nmat-1, 3)/))
        cell%umn(3,1:prob%Nmat/3) = ut((/(i, i=row + 2, row + prob%Nmat  , 3)/))

!       Volume correction: small, inward normal velocity based off current volume/SA/time step
!       Removed for reduce, because that keeps things at constant volume
        if(.not. reduce) THEN
            zm = -(cell%Vol() - cell%V0)/(cell%SA()*prob%dt)
            cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
        ENDIF

        umnt = cell%umn
        xmnt = cell%xmn

!       Volume reduction (add small inward normal vel every timestep)
        IF(reduce) THEN
            IF(cell%Vol().gt. 4.22  ) umnt = umnt - 0.10D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
            IF(cell%Vol().lt. 4.185 ) umnt = umnt + 0.01D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
            IF(cell%Vol().gt. 4.1894) umnt = umnt - 0.01D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
        ENDIF
    
!       Second part for midpoint
        IF(ord .eq. 1) THEN
            ! CALL cell%derivs()
            ! CALL cell%stress()
            ! CALL cell%fluid(prob, A2, b2)
            ! zm = -(cell%Vol() - cell%V0)/(3D0*cell%SA()*prob%dt)
            ! cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
            ! umnt = 0.5D0*(umnt + cell%umn)
            ! cell%xmn = xmnt + umnt*prob%dt
            ! cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
            ! cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
            ! cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
        ELSEIF(ord .ne. 1) THEN
            print*, "ERROR: time advancement of order >1 not supported"
            stop
        ENDIF
    
!       Update position and current time step
        cell%xmn = xmnt + umnt*prob%dt


! !       Simple periodic
!         IF(REAL(cell%xmn(1,1)).lt.-12D0*sqrt(pi)) THEN
!             cell%xmn(1,1) = cell%xmn(1,1) + 24D0*sqrt(PI)
!         ELSEIF(REAL(cell%xmn(1,1)).gt.12D0*sqrt(pi)) THEN
!             cell%xmn(1,1) = cell%xmn(1,1) - 24D0*sqrt(PI)
!         ENDIF

        cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
    ENDDO

    CALL SYSTEM_CLOCK(toc)
    !print *, REAL(toc-tic)/REAL(rate)

    prob%cts = prob%cts + 1
    prob%t = prob%t + prob%dt

!   Check if there's any funny business
    IF(isNaN(MAXVAL(ABS(cell%umn))) .or. MAXVAL(ABS(cell%umn)) .gt. HUGE(zm)) THEN
            print *, 'ERROR: inftys or NaNs'
            STOP
    ENDIF
END SUBROUTINE UpdateProb

! -------------------------------------------------------------------------!
! Continue from where we left off at? 
!!! Should make this more robust (read Params, but at the top, etc.)
SUBROUTINE ContinueProb(prob)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell
    CHARACTER (LEN = 25) ctsst, contfile, cfile2
    INTEGER ic, endt, stat, p, i, jmp
    REAL(KIND = 8), ALLOCATABLE :: xmnraw(:,:)
    REAL(KIND = 8) :: xmnrawind
    LOGICAL :: exist_in

    DO ic = 1,prob%NCell
        cell => prob%cell(ic)
        INQUIRE(FILE = TRIM('dat/'//TRIM(cell%fileout)//'/maxdt'), EXIST = exist_in)
        IF(.not. exist_in) THEN
            print *, "ERROR: Files not found for continue. Have you run previous cases with the same number of cells?"
            stop
        ENDIF

!       Get the timestep of where we left off
        contfile = TRIM('dat/'//TRIM(cell%fileout))
        cfile2 = TRIM('dat/'//TRIM(cell%fileout)//'/maxdt')
        OPEN(unit = 13, file = TRIM(cfile2), action = 'read')
        READ(13, '(I16)', iostat = stat) endt
        CLOSE(13)
        prob%cts = endt

!       Open the file of where we left off and get the shape
        write(cfile2, "(I0.5)") endt
        cfile2 = TRIM('dat/'//TRIM(cell%fileout)//'/x_'//TRIM(cfile2))

!       Read in the files
        p = (cell%p + 1)*(cell%p + 1)*2
        ALLOCATE(xmnraw(3,p))
        OPEN(unit = 13, file = cfile2, action = 'read')
        DO i = 1,p
            READ(13, *, iostat = stat) xmnraw(:,i)
        ENDDO
        CLOSE(13)

!       Text file format: all real, then imag
        p = cell%p
        jmp = (p+1)*(p+1)
    
        cell%xmn = 0D0

        cell%xmn(1,1:(p+1)*(p+1)) = xmnraw(1,1:(p+1)*(p+1))
        cell%xmn(2,1:(p+1)*(p+1)) = xmnraw(2,1:(p+1)*(p+1))
        cell%xmn(3,1:(p+1)*(p+1)) = xmnraw(3,1:(p+1)*(p+1))

        cell%xmn(1,1:(p+1)*(p+1)) = cell%xmn(1,1:(p+1)*(p+1)) + xmnraw(1, jmp+1: 2*jmp)*ii
        cell%xmn(2,1:(p+1)*(p+1)) = cell%xmn(2,1:(p+1)*(p+1)) + xmnraw(2, jmp+1: 2*jmp)*ii
        cell%xmn(3,1:(p+1)*(p+1)) = cell%xmn(3,1:(p+1)*(p+1)) + xmnraw(3, jmp+1: 2*jmp)*ii
        
    ENDDO
    prob%t = prob%cts*prob%dt
END SUBROUTINE ContinueProb

! -------------------------------------------------------------------------!
! Runs until initial cell is relaxed
SUBROUTINE RelaxCell(cell, prob, tol)
    CLASS(cellType), INTENT(INOUT) :: cell
    TYPE(probType), INTENT(IN) :: prob
    REAL(KIND = 8), INTENT(IN) :: tol
    COMPLEX(KIND = 8), ALLOCATABLE :: A2(:,:), b2(:), ut(:), wrk(:)
    INTEGER iter, info, i
    REAL(KIND = 8), ALLOCATABLE :: rwrk(:)
    COMPLEX, ALLOCATABLE :: swrk(:)
    INTEGER, ALLOCATABLE :: IPIV(:)

    ALLOCATE(ut(prob%Nmat), &
             IPIV(prob%Nmat), wrk(prob%Nmat), &
             swrk(prob%Nmat*(prob%Nmat+1)), &
             rwrk(prob%Nmat))

    cell%umn(1,1) = 1/cell%Ca
    DO WHILE(MAXVAL(ABS(cell%umn))*cell%Ca .gt. tol)
            CALL cell%derivs()
            CALL cell%stress() 
            CALL cell%fluid(prob, A2, b2)
            CALL zcgesv(prob%Nmat, 1, A2, prob%Nmat, IPIV, b2, prob%Nmat, ut, prob%Nmat, wrk, swrk, rwrk, iter, info)
            cell%umn = 0D0
            cell%umn(1,1:prob%Nmat/3) = ut((/(i, i=1,prob%Nmat-2, 3)/))
            cell%umn(2,1:prob%Nmat/3) = ut((/(i, i=2,prob%Nmat-1, 3)/))
            cell%umn(3,1:prob%Nmat/3) = ut((/(i, i=3,prob%Nmat  , 3)/))
            cell%umn(:,1) = 0D0
            cell%xmn = cell%xmn + cell%umn*prob%dt
            cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
            cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:)) 
            cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
            write(*,'(F8.6)') MAXVAL(ABS(cell%umn))*cell%Ca
    ENDDO
END SUBROUTINE RelaxCell
! -------------------------------------------------------------------------!
! Some output
SUBROUTINE OutputProb(prob)
    CLASS(probType), INTENT(IN) :: prob
    INTEGER ic
    IF(prob%cm%slv()) RETURN

    DO ic = 1, prob%NCell
        write(*,'(I5,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') & 
        prob%cts, prob%t, 1D0*2D0*MAXVAL(ABS(prob%cell(ic)%ff))/prob%cell(ic)%B, &
        MAXVAL(ABS(prob%cell(ic)%umn)), prob%cell(ic)%vol(), prob%cell(ic)%SA()
    ENDDO
END SUBROUTINE OutputProb

! -------------------------------------------------------------------------!
! Takes an input and cell and de-alises it
FUNCTION DealiasCell(cell, f) RESULT(fc)
    CLASS(cellType), INTENT(IN) :: cell
    REAL(KIND = 8), INTENT(IN) :: f(:,:)
    REAL(KIND = 8), ALLOCATABLE :: fc(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:)

    ALLOCATE(fc(cell%Y%nt, cell%Y%np))
!   Check to make sure f is the right size
    IF(SIZE(f) .ne. cell%Yf%nt*cell%Yf%np) THEN
        print *, 'ERROR: trying to de-alias something of wrong size'
        STOP
    ENDIF

    fmn = cell%Yf%forward(f, cell%q)
    fc = cell%Yf%backward(fmn, cell%p)
END FUNCTION

! -------------------------------------------------------------------------!
! Calculate the surface area of a cell
FUNCTION SAcell(cell) RESULT(SA)
    CLASS(cellType), INTENT(IN) :: cell
    REAL(KIND = 8) :: SA
    INTEGER :: i, j
    SA = 0D0

!   Basically the integrand is just the infinitesimal area element
    DO i = 1,cell%Yf%nt
        DO j = 1,cell%Yf%np
!           Integrate via Gauss quad
            SA = SA + cell%Yf%wg(i)*cell%J(i,j)*cell%Yf%dphi/sin(cell%Yf%tht(i))
        ENDDO
    ENDDO
END FUNCTION SAcell

! -------------------------------------------------------------------------!
! Calculate volume of the cell
! https://math.stackexchange.com/questions/709566/compute-the-volume-bounded-by-a-parametric-surface
FUNCTION Volcell(cell) RESULT(V)
    CLASS(cellType), INTENT(IN), TARGET :: cell
    REAL(KIND = 8):: V, intgd
    REAL(KIND = 8), POINTER :: x(:,:,:), xt(:,:,:), xp(:,:,:)
    INTEGER :: i, j

!   Integral we need is int_0^2pi int_0^pi r(theta)^3/3 sin(theta) d theta d phi.
!   Calculate r at int points
    V = 0D0

    x  => cell%xf
    xt => cell%dxt
    xp => cell%dxp

    DO i = 1,cell%Yf%nt
        DO j = 1,cell%Yf%np
!           Integrand
            intgd = x(3,i,j)*(xt(1,i,j)*xp(2,i,j) - xp(1,i,j)*xt(2,i,j))

!           Integrate via Gauss quad
            V = V + intgd* &
            cell%Yf%wg(i)*cell%Yf%dphi/sin(cell%Yf%tht(i))
        ENDDO
    ENDDO
END FUNCTION Volcell

! -------------------------------------------------------------------------!
! Integrates quantity over the surface of the cell (fine grid)
FUNCTION Intgcell(cell, x) RESULT(intg)
    CLASS(cellType), INTENT(IN), TARGET :: cell
    REAL(KIND = 8), INTENT(IN) :: x(:,:)
    TYPE(YType), POINTER :: Y
    REAL(KIND = 8) :: intg
    INTEGER :: i, j
    Y => cell%Yf
    intg = 0D0

!   Basically the integrand is just the infinitesimal area element
    DO i = 1,Y%nt
        DO j = 1,Y%np
!           Integrate via Gauss quad
            intg = intg + x(i,j)*Y%wg(i)*cell%J(i,j)*Y%dphi/sin(Y%tht(i))
        ENDDO
    ENDDO
END FUNCTION Intgcell

! -------------------------------------------------------------------------!
! Functions to calculate the kernels
FUNCTION Gij(r,eye) RESULT(A)
    REAL(KIND = 8) r(3), A(3,3), eye(3,3), mri
    mri = 1/(sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
    A(1,1) = r(1)*r(1)
    A(2,2) = r(2)*r(2)
    A(3,3) = r(3)*r(3)
    A(1,2) = r(1)*r(2)
    A(1,3) = r(1)*r(3)
    A(2,3) = r(2)*r(3)
    A(3,2) = A(2,3)
    A(3,1) = A(1,3)
    A(2,1) = A(1,2)
    A = A*mri*mri*mri + eye*mri
END FUNCTION Gij

FUNCTION Tij(r, n) RESULT(A)
    REAL(KIND = 8) r(3), A(3,3), n(3), mri
    mri = 1/(sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
    A(1,1) = r(1)*r(1)
    A(2,2) = r(2)*r(2)
    A(3,3) = r(3)*r(3)
    A(1,2) = r(1)*r(2)
    A(1,3) = r(1)*r(3)
    A(2,3) = r(2)*r(3)
    A(3,2) = A(2,3)
    A(3,1) = A(1,3)
    A(2,1) = A(1,2)
    A = -6D0*A*(mri*mri*mri*mri*mri)*(r(1)*n(1) + r(2)*n(2) + r(3)*n(3))
END FUNCTION Tij
! -------------------------------------------------------------------------!
! Periodic functions to calculate the kernels
! Arguments are (in order) distance vector, wavenumber vector, lattice vectors,
!   inverse lattice vectors, Ewald parameter, cutoff point (for loop), dij
! FUNCTION PGij(r, k, lv, ilv, xi, cut, eye) RESULT(A)
!     REAL(KIND = 8) r(3), k(3), lv(3,3), ilv(3,3), xi, &
!                    eye(3,3), A(3,3), rcur(3), kcur(3), V
!     INTEGER cut, i, j, k

! !   Calculate volume

! !   We do the real and Fourier sums in the same loops
!     DO i = -cut, cut
!         DO j = -cut,cut
!             DO k = -cut cut
! !               Skip current box(???)

! !               Real part (get current vector first)
!                 rcur = r + i*lv(:,1) + j*lv(:,2) + k*lv(:,3)

! !               Fourier part
!                 kcur = k + i*ilv(:,1) + j*ilv(:,2) + k*ilv(:,3)


!             ENDDO
!         ENDDO
!     ENDDO



!     mri = 1/(sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
!     A(1,1) = r(1)*r(1)
!     A(2,2) = r(2)*r(2)
!     A(3,3) = r(3)*r(3)
!     A(1,2) = r(1)*r(2)
!     A(1,3) = r(1)*r(3)
!     A(2,3) = r(2)*r(3)
!     A(3,2) = A(2,3)
!     A(3,1) = A(1,3)
!     A(2,1) = A(1,2)
!     A = A*mri*mri*mri + eye*mri
! END FUNCTION PGij

END MODULE SHAPEMOD