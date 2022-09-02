MODULE PROBMOD
USE SHAPEMOD
USE CMMOD
USE PERIODICMOD
IMPLICIT NONE

!==============================================================================!
!              The purpose of this module is perform calculations              !
!            that solve the overall math problem, including updating           !
!               the timestep, performing matrix inversions, and                !
!           facilitating a loop to calculate all the cells' responses          !
!==============================================================================!

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! Contains all the miscellaneous info about the problem
TYPE probType
    INTEGER :: cts, NT, dtinc, NCell, nts
    REAL(KIND = 8) :: t, kdt, kfr
    REAL(KIND = 8), ALLOCATABLE :: G(:,:,:)

!   The location of the shared container
    TYPE(sharedType), POINTER :: info

!   Pointer to the cells
    TYPE(cellType), POINTER :: cell(:)
    
!   Name of output file
    CHARACTER(:), ALLOCATABLE :: fileout

!   MPI stuff
    TYPE(cmType), POINTER :: cm
    INTEGER :: PCells(2)

    CONTAINS
    PROCEDURE :: Update  => UpdateProb
    PROCEDURE :: Write   => WriteProb
    PROCEDURE :: Output  => OutputProb
    PROCEDURE :: Continue=> ContinueProb
    PROCEDURE :: StrnRt  => StrnRtProb
END TYPE probType

! -------------------------------------------------------------------------!
INTERFACE probType
    PROCEDURE :: newprob
END INTERFACE probType

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

! Constructs the problem information
FUNCTION newprob(filein, reduce, cm, info) RESULT(prob)
    CHARACTER(len = *), INTENT(IN) :: filein
    LOGICAL, INTENT(IN) :: reduce
    TYPE(cmType), INTENT(IN), TARGET :: cm
    TYPE(sharedType), TARGET :: info
    TYPE(probType), TARGET :: prob
    TYPE(cellType) :: celltmp
    CHARACTER(:), ALLOCATABLE :: fileout, cont, gfile
    CHARACTER(len = 30) :: icC
    INTEGER :: m, n, ic, pthline
    REAL(KIND = 8) :: Gfac
    REAL(KIND = 8), ALLOCATABLE :: Gtmp(:,:,:,:)
    LOGICAL :: fl = .false.
    TYPE(cellprops) props

    prob%cm => cm
    prob%info => info

!   General problem parameters
    CALL READ_MFS(prob%NT, filein, 'Max_time_steps')
    CALL READ_MFS(prob%NCell, filein, 'Number_cells')
    CALL READ_MFS(prob%dtinc, filein, 'Time_inc')
    CALL READ_MFS(gfile, filein, 'Gradient_file')
    CALL READ_MFS(fileout, filein, 'Output')
    CALL READ_MFS(cont, filein, 'Continue')
    CALL READ_MFS(prob%kdt, filein, 'Shear_dim_time')
    CALL READ_MFS(prob%kfr, filein, 'Gradient_timestep')
    prob%kfr = prob%kfr/prob%kdt ! Fraction of ts between velgrads
    CALL READ_MFS(pthline, filein, 'Path_line')
    prob%cts = 0
    prob%t = 0D0

!   The rest of the parameters specific to the cell go here
    CALL READ_MFS(props%lam, filein, 'Viscosity_Ratio')
    CALL READ_MFS(props%Ca, filein, 'Capillary')
    CALL READ_MFS(props%C, filein, 'Dilatation_Ratio')
    CALL READ_MFS(props%Eb, filein, 'Bending_Modulus')
    CALL READ_MFS(props%c0, filein, 'Spont_Curvature')
    CALL READ_MFS(props%int_pres, filein, 'Internal_Pressure')

!   Initialize info now
    IF(prob%cm%mas()) print *, 'Done with init file!'
    CALL info%init()

!   Check if there are more processors than cells (can't handle)
    IF(prob%cm%np() .gt. prob%NCell) THEN
        print *, 'ERROR: More processors than cells'
        STOP
    ENDIF
    
    ALLOCATE(prob%cell(prob%NCell))

!   Make a master cell of sorts
    celltmp = cellType(filein, reduce, info, props)

!   Loop and make individual cells
    DO ic = 1, prob%NCell
        prob%cell(ic) = celltmp
        prob%cell(ic)%id = ic
        write(icC, "(I0.1)") ic
        prob%cell(ic)%fileout = TRIM('cell_'//icC)
    ENDDO
    prob%fileout=TRIM(fileout)
    
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

!   Velocity gradient information
    IF(prob%cm%mas()) print *, 'Reading in velocity gradient...'
    OPEN(1,FILE = gfile, ACCESS = 'stream', ACTION = 'read')
    READ(1) prob%nts, ic
!   Read into a temporary array so we don't hold onto this big one
    ALLOCATE(Gtmp(prob%nts,3,3,ic))
    READ(1) Gtmp
    CLOSE(1)

!   How many timesteps from G do we actually need?
    Gfac = prob%nts*prob%kfr/(prob%NT*prob%info%dt) ! Gfac: Ratio of total time in velgrad to total time requested
!   If we have fewer velocity gradient time steps than requested, set down to gradient
    IF(Gfac .le. 1) THEN
        prob%NT = FLOOR(prob%NT*Gfac)
        print *, "Warning: Fewer gradient time steps than requested total time steps"
        print *, "Setting total timesteps to "
        print *, prob%NT
        fl = .true.
    ELSE
        prob%nts = CEILING(prob%nts/Gfac)
    ENDIF
    
!   To prevent only having 2 timesteps
    IF(prob%nts .eq. 1) THEN
        ALLOCATE(prob%G(3, 3, 3))
        prob%G = Gtmp(1:3, :, :, pthline)
    ELSE
!       When we need to use the whole grad file
        IF(.not. fl) THEN
            ALLOCATE(prob%G(prob%nts+1, 3, 3))
            prob%G = Gtmp(1:prob%nts+1, :, :, pthline)
        ELSE
            ALLOCATE(prob%G(prob%nts, 3, 3))
            prob%G = Gtmp(1:prob%nts, :, :, pthline)
        ENDIF
    ENDIF

!   Inital support point array
    IF(prob%info%periodic) THEN
        CALL SuppPoints(prob%info)
    ENDIF

!   Relax cell to get it into an equilibrium state, get volumes
!!  ============================
    print *, 'Getting initial shape -'
    print *, 'Max velocity coefficient w.r.t. membrane time (want less than 0.005*Ca):'
    Gfac = (((4D0*PI/3D0)/((0.95D0**(3D0))*(PI/6D0)))**(1D0/3D0))/2D0
    DO ic = 1, prob%NCell
            ! CALL prob%cell(ic)%relax(0.005D0)!!!
            CALL prob%cell(ic)%derivs()
            CALL prob%cell(ic)%stress()
            prob%cell(ic)%V0 = prob%cell(ic)%Vol()
            !! Test                                        !!!!!
            SELECT CASE(MOD(ic, 4))!9) )
            CASE(1)
                prob%cell(ic)%xmn(1,1) = (Gfac)/(ispi*0.5D0)!0.50D0/(ispi*0.5D0)!
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!0.50D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (Gfac)/(ispi*0.5D0)!0.50D0/(ispi*0.5D0)!
            CASE(3)
                prob%cell(ic)%xmn(1,1) = (Gfac)/(ispi*0.5D0) !+ info%bvl/(ispi*0.5D0)!3D0/(ispi*0.5D0)! 
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!0.5D0/(ispi*0.5D0)!
            CASE(2)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)!
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!0.5D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!1.5D0/(ispi*0.5D0)!
            CASE(0)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!0.5D0/(ispi*0.5D0)!
                prob%cell(ic)%xmn(2,1) = (-1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (Gfac)/(ispi*0.5D0)!1.5D0/(ispi*0.5D0)!
            CASE(5)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!3.5D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!2.75D0/(ispi*0.5D0)
            CASE(6)
                prob%cell(ic)%xmn(1,1) = (5D0*Gfac)/(ispi*0.5D0)!4.5D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!4.25D0/(ispi*0.5D0)
            CASE(7)
                prob%cell(ic)%xmn(1,1) = (Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!0.1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
            CASE(8)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!0.1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
            CASE(4)
                prob%cell(ic)%xmn(1,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!0.1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
            END SELECT
            prob%cell(ic)%x(1,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(1,:))
            prob%cell(ic)%x(2,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(2,:))
            prob%cell(ic)%x(3,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(3,:))
    ENDDO
!    ============================
    
!   Should we continue from a spot where we left off? !! Not working right now,
    !! doesn't return the same values for some reason
    IF(cont .eq. "Yes") CALL prob%Continue()
    
!   Write initial configuration
    CALL prob%write()

END FUNCTION newprob

! -------------------------------------------------------------------------!
! Writes xmn to a text file. Could be better
SUBROUTINE WriteProb(prob)
    TYPE(cellType), POINTER :: cell
    CLASS(probType), INTENT(IN) :: prob
    CHARACTER (LEN = 45) ctsst, datdir, filename
    INTEGER ic
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:)

    IF(prob%cm%slv()) RETURN
!   Don't write if it's not a timestep to write
    IF(.not.((prob%cts .eq. 1) .or. &
         (MOD(prob%cts,prob%dtinc)) .eq. 0)) RETURN

!   Write top-level directory containing all lower level directoriesx
    IF(.not. prob%cell(1)%writ) THEN
        datdir = TRIM('dat/'//prob%fileout//'/')
        CALL MAKEDIRQQ(datdir)
    ENDIF

    DO ic = 1,prob%NCell
    cell => prob%cell(ic)
    IF(NOT(ALLOCATED(fmn))) THEN
        ALLOCATE(fmn(3,(cell%info%q+1)*(cell%info%q+1)))
    ENDIF

!   Formatting pain
    write(ctsst, "(I0.5)") prob%cts
    datdir = TRIM('dat/'//prob%fileout//'/'//cell%fileout//'/')
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
        WRITE(88,*) cell%info%p
        WRITE(88,*) "dt"
        WRITE(88,*) prob%info%dt
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
! Problem advancement
SUBROUTINE UpdateProb(prob, ord, reduce)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell, celli
    INTEGER, INTENT(IN) :: ord
    LOGICAL, INTENT(IN) :: reduce
    TYPE(YType), POINTER :: Y
    TYPE(sharedType), POINTER :: info
    REAL(KIND = 8) :: zm
    COMPLEX(KIND = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:), ut(:), wrk(:), &
                                      A2(:,:), b2(:), A(:,:), b(:), &
                                      fmnR(:,:), xmnR(:,:), &
                                      fv1(:), fv2(:), fv3(:), fv(:,:), &
                                      b2f(:), b2r(:), YVt(:), um(:,:,:,:,:), &
                                      A2f(:,:,:,:,:,:,:,:), A2r(:,:,:,:,:,:)
    INTEGER :: ic, ic2, i, row, col, iter, p_info, m, n, im, im2, it
    REAL(KIND = 8), ALLOCATABLE :: rwrk(:), utr(:), uti(:), nJt(:,:,:), xv(:,:)
    COMPLEX, ALLOCATABLE :: swrk(:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER(KIND = 8) tic, toc, rate
    TYPE(nmType), ALLOCATABLE :: YV(:)

    Y => prob%info%Y
    info => prob%info
    ALLOCATE(ut(info%NmatT), &
             uti(info%NmatT), &
             utr(info%NmatT), &
             IPIV(info%NmatT), wrk(info%NmatT), &
             swrk(info%NmatT*(info%NmatT+1)), &
             rwrk(info%NmatT), &
             A(info%NmatT, info%NmatT), &
             b(info%NmatT), &
             fv(3,Y%nt*Y%np*prob%NCell), &
             YV(Y%p), &
             nJt(3, Y%nt, Y%np), &
             um(3,3,Y%nt,Y%np, prob%NCell), &
             fmnR(3, (info%p+1)*(info%p+1)), &
             xmnR(3, (info%p+1)*(info%p+1)), &
             A2r(3,3, 2*(Y%p-1)+1, Y%p, Y%nt, Y%np), &
             b2r(3*info%NmatT))

!   If periodic, make variables for the Fourier part
    IF(info%periodic) THEN
        ALLOCATE( &
        A2f(3,3, 2*(Y%p-1)+1, Y%p, Y%nt, Y%np, prob%Ncell, prob%Ncell), & !!! Note that in parallel, many entries will be 0 here
        b2f(3*Y%nt*Y%np*prob%NCell))
        A2f = 0D0
    ENDIF


!   Do interpolation to get current grad tensor, then normalize by kolm time
    info%dU = VelInterp(prob%G,prob%t,prob%nts,prob%kfr)*prob%kdt
!! ============================
!   Hardcoded shear
    ! info%dU = 0D0 !!!
    ! info%dU(1,3) = 1D0
    ! prob%dU(1,1) =  1D0
    ! prob%dU(3,3) = -1D0
!! ============================

!   First, calculate geometric quantities and stresses in all of the cells
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)

!       For periodic, check if cells are outside primary cell, and put them in if so
        IF(info%periodic) THEN
            CALL cell%inPrim()
        ENDIF

        CALL cell%derivs()
        CALL cell%stress()

!       Construct a large vector containing all the forces/locations in
!       preparation for the fast Ewald summation
        IF(info%periodic) THEN
!           Forces and locations
            fmnR = cell%fmn(:,1:(info%p+1)*(info%p+1))
            xmnR = cell%xmn(:,1:(info%p+1)*(info%p+1))
!           Append vectors
            CALL EwaldPreCalc(fv1, Y, xv, fmn=fmnR(1,:), xmn=xmnR)
            CALL EwaldPreCalc(fv2, Y, fmn=fmnR(2,:))
            CALL EwaldPreCalc(fv3, Y, fmn=fmnR(3,:))
        ENDIF
    ENDDO

!   Now we need to calculate the SpH velocity coefficients of the cells.
!   For periodic, we calculate LHS matrix and RHS vector in two parts:
!   Long range (via FFTs) and short range (directly, truncated at low r).
!   We calculate these two matrices individually at integration points,
!   add them, and then perform the Galerkin transfo on these summed matrices.

    CALL SYSTEM_CLOCK(tic,rate)
    IF(info%periodic) THEN
!       Calculate RHS vector
        fv(1,:) = fv1 
        fv(2,:) = fv2
        fv(3,:) = fv3
        CALL EwaldG(info=info, x0=xv, f=fv, u1=b2f, full=.true.)

!       Start with the long range part, as this calculates all points for one integral surface
!       The "forces" here are the spherical harmonics multiplied by the normal vector
        DO ic = prob%PCells(1), prob%PCells(2)
            cell => prob%cell(ic)

!           We need to make a vector for all of the spherical harmonics,
!           multiplied by normal vector, for this integration surface
            fmnR(1,:) = cell%nkt(1,1:(info%p+1)*(info%p+1))
            fmnR(2,:) = cell%nkt(2,1:(info%p+1)*(info%p+1))
            fmnR(3,:) = cell%nkt(3,1:(info%p+1)*(info%p+1))

            nJt(1,:,:) = Y%backward(fmnR(1,:), info%p)
            nJt(2,:,:) = Y%backward(fmnR(2,:), info%p)
            nJt(3,:,:) = Y%backward(fmnR(3,:), info%p)

            DO n = 0, Y%p-1
                IF(.not.(ALLOCATED(YV(n+1)%v))) ALLOCATE(YV(n+1)%v(n+1, 3, Y%nt*Y%np))
                im = 0
                DO m = -n,n
                    im = im + 1
                    im2 = m + (Y%p)
                    IF(m .gt. 0) THEN ! Symmetry
                        A2f(:,:,im2,n+1,:,:,:, ic) = CONJG(A2f(:,:,im2 - 2*m,n+1,:,:,:, ic))*(-1D0)**m
                        CYCLE
                    ENDIF

                    CALL EwaldPreCalc(YVt, Y, fin=Y%nm(n+1)%v(im,:,:)*nJt(1,:,:))
                    YV(n+1)%v(im,1,:) = YVt
                    DEALLOCATE(YVt)

                    CALL EwaldPreCalc(YVt, Y, fin=Y%nm(n+1)%v(im,:,:)*nJt(2,:,:))
                    YV(n+1)%v(im,2,:) = YVt
                    DEALLOCATE(YVt)

                    CALL EwaldPreCalc(YVt, Y, fin=Y%nm(n+1)%v(im,:,:)*nJt(3,:,:))
                    YV(n+1)%v(im,3,:) = YVt
                    DEALLOCATE(YVt)

!                   Now do actual calcs for all points at this int surface/nm combo
                    CALL EwaldT(info=info, x0=xv, f=YV(n+1)%v(im,:,:), full=.true., um=um, strt=ic)
                    A2f(:,:, im2, n+1,:,:,:, ic) = um
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    CALL SYSTEM_CLOCK(toc)
    ! print *, REAL(toc-tic)/REAL(rate)

    A = 0D0
    b = 0D0
    b2 = 0D0
    A2 = 0D0

    CALL SYSTEM_CLOCK(tic,rate)
!   Now the real space part
!   First loop: velocity surface
    DO ic = 1,prob%NCell
        row = (ic - 1)*info%Nmat + 1
        cell => prob%cell(ic)

!       Second loop: integral surface
        DO ic2 = prob%PCells(1), prob%PCells(2)
            col = (ic2-1)*info%Nmat + 1
            celli => prob%cell(ic2)
            
!           Get velocity sub-matrix for cell-cell combo (Same or diff cell)
            IF(ic .eq. ic2) THEN
                CALL cell%fluid(Ao = A2r, bo = b2r)
            ELSE
                CALL cell%fluid(Ao = A2r, bo = b2r, celli = celli)
            ENDIF

!           Add long range interactions if periodic
            IF(info%periodic) THEN
                A2r = A2r + A2f(:,:,:,:,:,:,ic,ic2)*(1D0-cell%lam)/(1D0+cell%lam)
!               This is b/c b2f is the vector at the velocity points integrated across all surfaces,
!               so we only add in once
                IF(ic .eq. ic2) &
                    b2r = b2r + b2f( (ic-1)*(3*Y%nt*Y%np)+1 : (ic)*(3*Y%nt*Y%np))/(1D0 + cell%lam)
            ENDIF

!           Perform the Galerkin
            CALL info%Gal(A2r, b2r, A2, b2)
!           Put sub matrix into big matrix
            A(row:row + info%Nmat - 1, col:col + info%Nmat - 1) = A2
!           Sum over all the integrals
            b(row:row + info%Nmat - 1) = b(row:row + info%Nmat - 1) + b2
        ENDDO
    ENDDO
    CALL SYSTEM_CLOCK(toc)
    ! print *, REAL(toc-tic)/REAL(rate)

    DO i = 1,info%NmatT
        A(:,i) = prob%cm%reduce(REAL(A(:,i))) + prob%cm%reduce(AIMAG(A(:,i)))*ii
    ENDDO
    b = prob%cm%reduce(REAL(b)) + prob%cm%reduce(AIMAG(b))*ii

!   Invert big matrix to get a list of all the vel constants of all cells
    ut = 0D0
    IF(prob%cm%mas()) THEN
        CALL zcgesv(info%NmatT, 1, A, info%NmatT, IPIV, b, info%NmatT, &
                    ut, info%NmatT, wrk, swrk, rwrk, iter, p_info)
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
        row = (ic-1)*info%Nmat + 1
        cell%umn = 0D0
        cell%umn(1,1:info%Nmat/3) = ut((/(i, i=row    , row + info%Nmat-2, 3)/))
        cell%umn(2,1:info%Nmat/3) = ut((/(i, i=row + 1, row + info%Nmat-1, 3)/))
        cell%umn(3,1:info%Nmat/3) = ut((/(i, i=row + 2, row + info%Nmat  , 3)/))

!       Volume correction: small, inward normal velocity based off current volume/SA/time step
!       Removed for reduce, because that keeps things at constant volume
        if(.not. reduce) THEN
            zm = -(cell%Vol() - cell%V0)/(cell%SA()*info%dt)
            cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
        ENDIF

        umnt = cell%umn
        xmnt = cell%xmn

!       Volume reduction (add small inward normal vel every timestep)
        IF(reduce) THEN
            IF(cell%Vol().gt. 4.22  ) umnt = umnt - 0.10D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
            IF(cell%Vol().lt. 4.185 ) umnt = umnt + 0.01D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
            IF(cell%Vol().gt. 4.1894) umnt = umnt - 0.01D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
        ENDIF
    
!       Second part for midpoint
        IF(ord .eq. 1) THEN
            ! CALL cell%derivs()
            ! CALL cell%stress()
            ! CALL cell%fluid(prob, A2, b2)
            ! zm = -(cell%Vol() - cell%V0)/(3D0*cell%SA()*info%dt)
            ! cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
            ! umnt = 0.5D0*(umnt + cell%umn)
            ! cell%xmn = xmnt + umnt*info%dt
            ! cell%x(1,:,:) = Y%backward(cell%xmn(1,:))
            ! cell%x(2,:,:) = Y%backward(cell%xmn(2,:))
            ! cell%x(3,:,:) = Y%backward(cell%xmn(3,:))
        ELSEIF(ord .ne. 1) THEN
            print*, "ERROR: time advancement of order >1 not supported"
            stop
        ENDIF
    
!       Update position and current time step
        cell%xmn = xmnt + umnt*info%dt

        cell%x(1,:,:) = Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = Y%backward(cell%xmn(3,:))
    ENDDO

!   Advance periodic basis vectors
    IF(info%periodic) THEN
        CALL info%bvAdvance()
!       Update array containing support points
        CALL SuppPoints(info)
    ENDIF

    prob%cts = prob%cts + 1
    prob%t = prob%t + info%dt

!   Check if there's any funny business
    IF(ANY(isNaN((ABS(cell%umn)))) .or. ANY((ABS(cell%umn)) .gt. HUGE(zm))) THEN
            print *, 'ERROR: inftys or NaNs'
            STOP
    ENDIF

!   Write and display some output
    CALL prob%write()
    CALL prob%output()
END SUBROUTINE UpdateProb

! -------------------------------------------------------------------------!
! Finds equivalent strain rate
FUNCTION StrnRtProb(prob) RESULT(Sig)
    CLASS(probType), TARGET :: prob
    TYPE(sharedType), POINTER :: info
    TYPE(YType), POINTER :: Y
    TYPE(cellType), POINTER :: cell
    REAL(KIND = 8) :: Sig(3,3), tmp(3,3), lam, tau
    REAL(KIND = 8), ALLOCATABLE :: fMxPuMn(:,:,:,:), u(:,:,:), n(:,:,:)
    INTEGER :: i, j, ic, nt, np, d1, d2

    info => prob%info
    Y => info%Yf
    nt = Y%nt
    np = Y%np
    
    tau = DOT(CROSS(info%bv(:,1), info%bv(:,2)), info%bv(:,3))

    tmp = 0D0

    ALLOCATE(fMxPuMn(3,3,nt, np), u(3,nt, np), n(3,nt, np))

!   One paper says to use x - x_c,
!   but I shouldn't need to as long as net force in cell=0

    DO ic = 1, prob%NCell
        cell => prob%cell(ic)
        lam = cell%lam
        u(1,:,:) = Y%backward(cell%umn(1,:), info%P)
        u(2,:,:) = Y%backward(cell%umn(2,:), info%P)
        u(3,:,:) = Y%backward(cell%umn(3,:), info%P)

        n(1,:,:) = Y%backward(cell%nkmn(1,:))
        n(2,:,:) = Y%backward(cell%nkmn(2,:))
        n(3,:,:) = Y%backward(cell%nkmn(3,:))
!       Get velocity and normalss at integration points
        DO i = 1,nt
            DO j = 1,np
                fMxPuMn(:,:,i,j) = OUTER(cell%ff(:,i,j), cell%xf(:,i,j))  & !! Need Ca here? don't seem to need to.
                                 + (lam - 1D0)*(OUTER(u(:,i,j), n(:,i,j)) &     !! Run with sphere array, compare to Zinchenko
                                 + OUTER(n(:,i,j), u(:,i,j)))
            ENDDO
        ENDDO
        DO d1 = 1,3
        DO d2 = 1,3
            tmp(d1, d2) = tmp(d1, d2) + cell%intg(fMxPuMn(d1,d2,:,:))
        ENDDO
        ENDDO
        ! print *, u(1,(nt+1)/2,1)! - cell%x(3,2,1)
        ! print *, ' '
        ! print *, cell%intg(fMxPuMn(1,3,:,:)), u(1,1,1), SUM(u(1,:,:))/nt/np
    ENDDO
    ! stop

    Sig = tmp/tau + info%dU + TRANSPOSE(info%dU)

END FUNCTION StrnRtProb

! -------------------------------------------------------------------------!
! Continue from where we left off at? 
!!! Should make this more robust (read Params, but at the top, etc.)
SUBROUTINE ContinueProb(prob)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell
    CHARACTER (LEN = 25) contfile, cfile2
    INTEGER ic, endt, stat, p, i, jmp
    REAL(KIND = 8), ALLOCATABLE :: xmnraw(:,:)
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
        p = (cell%info%p + 1)*(cell%info%p + 1)*2
        ALLOCATE(xmnraw(3,p))
        OPEN(unit = 13, file = cfile2, action = 'read')
        DO i = 1,p
            READ(13, *, iostat = stat) xmnraw(:,i)
        ENDDO
        CLOSE(13)

!       Text file format: all real, then imag
        p = cell%info%p
        jmp = (p+1)*(p+1)
    
        cell%xmn = 0D0

        cell%xmn(1,1:(p+1)*(p+1)) = xmnraw(1,1:(p+1)*(p+1))
        cell%xmn(2,1:(p+1)*(p+1)) = xmnraw(2,1:(p+1)*(p+1))
        cell%xmn(3,1:(p+1)*(p+1)) = xmnraw(3,1:(p+1)*(p+1))

        cell%xmn(1,1:(p+1)*(p+1)) = cell%xmn(1,1:(p+1)*(p+1)) + xmnraw(1, jmp+1: 2*jmp)*ii
        cell%xmn(2,1:(p+1)*(p+1)) = cell%xmn(2,1:(p+1)*(p+1)) + xmnraw(2, jmp+1: 2*jmp)*ii
        cell%xmn(3,1:(p+1)*(p+1)) = cell%xmn(3,1:(p+1)*(p+1)) + xmnraw(3, jmp+1: 2*jmp)*ii
        
    ENDDO
    prob%t = prob%cts*prob%info%dt
END SUBROUTINE ContinueProb

! -------------------------------------------------------------------------!
! Some output
SUBROUTINE OutputProb(prob)
    CLASS(probType), INTENT(IN) :: prob
    INTEGER ic
    IF(prob%cm%slv()) RETURN

!   What the output means
    IF(MOD(prob%cts,50).eq.1) THEN
        PRINT *, "  CTS  Time    Cell  Max F     Max umn  Vol     SA"
        PRINT *, "----------------------------------------------------"
    ENDIF

!   The output
    DO ic = 1, prob%NCell
        write(*,'(I5,X,F8.4,X,I5, X, F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') & 
        prob%cts, prob%t, ic, MAXVAL(ABS(prob%cell(ic)%ff))*prob%cell(ic)%Ca, &
        MAXVAL(ABS(prob%cell(ic)%umn)), prob%cell(ic)%vol(), prob%cell(ic)%SA()
    ENDDO
END SUBROUTINE OutputProb

END MODULE PROBMOD