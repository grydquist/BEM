MODULE SHAPEMOD
USE SHAREDMOD
IMPLICIT NONE

!==============================================================================!
!              The purpose of this module is perform calculations              !
!            on the actual cell, including calculations of internal            !
!               stress and the forcing of the fluid on the cell                !
!==============================================================================!

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! The actual cell shape object
TYPE cellType

!   Current cell
    INTEGER :: id

!   Material Properties
    REAL(KIND = 8) :: mu, lam, B, C, Eb, c0, Ca, int_pres

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
    c2pR(:,:,:)

!   Force variables
    REAL(KIND = 8), ALLOCATABLE :: fab(:,:,:), ff(:,:,:), fc(:,:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:), nkmn(:,:), fmn2(:,:), &
                                      nkt(:,:), Jtmn(:)

!   Cell velocity constants
    COMPLEX(KIND = 8), ALLOCATABLE :: umn(:,:), xmn(:,:)
    REAL(KIND = 8), ALLOCATABLE :: u(:,:,:)

!   Shared info pointer
    TYPE(sharedType), POINTER :: info

!   Has this object been initialized?
    LOGICAL :: init = .false.
    LOGICAL :: writ = .false.

!   Logicals to tell which modes to use in reconstruction
    LOGICAL, ALLOCATABLE :: xrecon(:,:), urecon(:,:), nJrecon(:,:), frecon(:,:)

!   Tabulates which indices to rotate about for singular intg.
!   so I don't have to find it every GMRES iteratio
    INTEGER, ALLOCATABLE :: sing_inds(:,:,:,:), cut_box(:,:,:,:,:,:)

!   Name of output file
    CHARACTER(:), ALLOCATABLE :: fileout
    
    CONTAINS
    PROCEDURE :: Derivs  => Derivscell
    PROCEDURE :: Stress  => Stresscell
    PROCEDURE :: Relax   => RelaxCell
    PROCEDURE :: Dealias => Dealiascell
    PROCEDURE :: Vol     => Volcell
    PROCEDURE :: SA      => SAcell
    PROCEDURE :: Intg    => Intgcell
    PROCEDURE :: InPrim  => InPrimcell
    PROCEDURE :: Reparam => ReparamCell
    PROCEDURE :: LHS_real=> LHS_realCell
    PROCEDURE :: RHS_real=> RHS_realCell
END TYPE cellType

! -------------------------------------------------------------------------!
INTERFACE cellType
    PROCEDURE :: newcell
END INTERFACE cellType

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

! Constructs cell object, takes in order
FUNCTION newcell(filein, reduce, info, props) RESULT(cell)
    CHARACTER(len = *), INTENT(IN) :: filein
    LOGICAL, INTENT(IN) :: reduce
    TYPE(sharedType), INTENT(INOUT), TARGET :: info
    TYPE(cellType) :: cell
    CHARACTER(:), ALLOCATABLE :: restart, restfile
    REAL(KIND = 8) lam, Ca, C, Eb, c0, int_pres
    INTEGER :: nt, np, ntf, npf
    TYPE(cellprops), INTENT(IN) :: props

    cell%info => info

!   Material properties from input
    lam      = props%lam
    Ca       = props%Ca
    C        = props%C
    Eb       = props%Eb
    c0       = props%c0
    int_pres = props%int_pres

!   Restart location
    CALL READ_MFS(restart, filein, 'Restart')
!   Choose file to restart from (assuming that this is after deflation)
    CALL READ_MFS(restfile, filein, 'Restart_Loc')

!   Note that the file MUST be in the restart directory!
    restfile = 'restart/'//trim(restfile)

    nt  = info%Y%nt
    np  = info%Y%np
    ntf = info%Yf%nt
    npf = info%Yf%np

    cell%mu = 1D0
    cell%lam = lam

!   To fit dimesnionless parameters, we set a0 = 1, flow timescale = 1, mu = 1
!   and fit the rest from there
!   Ca = (shear rate)*mu*(cell radius)/(B/2)
    cell%Ca = Ca
    cell%B = 2D0/Ca
!   A note on the dilation modulus: many papers use (B*C) in the spot where
!   I use C, so I just make that correction here
    cell%C = C*cell%B
!   Similar situation for beinding modulus: the input is the non-dim
!   parameter E*b = Eb/(a_0^2*(B/2))
    cell%Eb = Eb*2D0*cell%B
    cell%c0 = c0
!   Internal (osmotic) pressure
    cell%int_pres = int_pres

!   Allocating everything, starting with derivatives
    ALLOCATE(cell%dxt(3,ntf,npf), cell%dxp(3,ntf,npf), cell%dxp2(3,ntf,npf), &
            cell%dxt2(3,ntf,npf), cell%dxtp(3,ntf,npf), cell%dxp3(3,ntf,npf), &
            cell%dxt3(3,ntf,npf), cell%dxt2p(3,ntf,npf), cell%dxtp2(3,ntf,npf), &
            cell%dxp4(3,ntf,npf), cell%dxtp3(3,ntf,npf), cell%dxt2p2(3,ntf,npf), &
            cell%dxt3p(3,ntf,npf), cell%dxt4(3,ntf,npf), &
            cell%J(ntf,npf), cell%x(3,nt,np), cell%xf(3,ntf,npf), cell%u(3,nt,np), &
            cell%xmn(3,(info%p+1)*(info%p+1)), cell%umn(3,(info%p + 1)*(info%p + 1)))

!   Reference items
    ALLOCATE(cell%kR(ntf,npf), cell%kdR(2,ntf,npf), cell%kd2R(3,ntf,npf), &
            cell%c1R(3,ntf,npf), cell%c2R(3,ntf,npf), cell%c1tR(3,ntf,npf), &
            cell%c1pR(3,ntf,npf), cell%c2tR(3,ntf,npf), cell%c2pR(3,ntf,npf))

!   Force items
    ALLOCATE(cell%fab(3,ntf,npf), cell%ff(3,ntf,npf), cell%fc(3,nt,np), &
                cell%fmn(3,(info%q+1)*(info%q+1)), cell%nkmn(3,(info%q+1)*(info%q+1)), &
                cell%fmn2(3,(info%q+1)*(info%q+1)), cell%nkt(3,(info%q+1)*(info%q+1)), &
                cell%Jtmn((info%q+1)*(info%q+1)))

!   Reconstruction items
    ALLOCATE(cell%xrecon(3, info%p + 1), cell%urecon (3, info%p + 1), &
             cell%frecon(3, info%p + 1), cell%nJrecon(3, info%p + 1))
    cell%xrecon  = .true.
    cell%urecon  = .true.
    cell%frecon  = .true.
    cell%nJrecon = .true.
    
    cell%ff = 0D0
    cell%fmn = 0D0
    cell%umn = 0D0

!   For a = 1, V = 4.18904795321178, SA = 16.8447913187040, sphere 6.50088174342271

!   First, we need to get the reference shape for the shear stress
    ! cell%xmn = RBCcoeff(cell%Y)
    ! cell%xmn = Cubecoeff(cell%Y)
    ! cell%xmn = Bactcoeff(cell%Y, 1D0)
    cell%xmn = Spherecoeff(info%Y, .76D0) !1D0)!!!! Reduced volume .997, .98, .95-> .9,.76,.65

!   Initial surfce derivatives/reference state
    CALL cell%Derivs()
    CALL cell%Stress()

!   Now scale to get an equivalent radius of 1
    IF(.not. reduce) THEN
        cell%xmn = cell%xmn/(3D0*cell%vol()/(4D0*pi))**(1D0/3D0)
!   If we want to do deflation, comment the above and scale to get right SA
    ELSE
        cell%xmn = cell%xmn*(16.8447913187040D0/cell%SA())**(1D0/2D0) 
    ENDIF

!   Initialize derivs again... I do this because I need the derivs to get the volume
    CALL cell%Derivs()
    CALL cell%Stress()

!   Now with a reference shape made, should we choose an intermediate point
!   to restart from?
    IF (restart .eq. "Yes") THEN
!       This gives us the right surface area with a blown up volume
!       The 16.8 is for equiv radius of 1
        cell%xmn = cell%xmn*sqrt(16.8447913187040D0/cell%SA())
!       Get referennce state
        CALL cell%Derivs()
        CALL cell%Stress()

!       And just set the constants equal, scaling for equiv radius
!       This assumes the equiv radius in input is 1.
        cell%xmn = Readcoeff(restfile, info%Y%p)
    ENDIF

!   Initial positions
    cell%x(1,:,:) = info%Y%backward(cell%xmn(1,:))
    cell%x(2,:,:) = info%Y%backward(cell%xmn(2,:))
    cell%x(3,:,:) = info%Y%backward(cell%xmn(3,:))

!   Characteristic grid spacing
    cell%h = SQRT(cell%SA()/nt**2D0)

    cell%init = .true.
END FUNCTION newcell

! -------------------------------------------------------------------------!
! Updates the values of the derivatives on the surface of the cell on fine grid.
SUBROUTINE Derivscell(cell, m_in, n_in, cm)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell
    TYPE(cmType), INTENT(IN), POINTER, OPTIONAL :: cm
    TYPE(nmType), POINTER :: nm
    TYPE(YType), POINTER :: Y
    INTEGER, INTENT(IN), OPTIONAL :: m_in(2), n_in(2)
    INTEGER :: ih, it, im, n, m, p, n_low, n_up, m_low, m_up, ml_cur, mu_cur
    COMPLEX(KIND = 8) :: f1, f2, f3
    COMPLEX(KIND = 8), ALLOCATABLE :: vcur(:,:), td1(:,:), td2(:,:), td3(:,:), td4(:,:)

    Y => cell%info%Yf
    p = cell%info%Y%p
    ALLOCATE(vcur(Y%nt, Y%np), td1(Y%nt, Y%np), td2(Y%nt, Y%np), &
    td3(Y%nt, Y%np), td4(Y%nt, Y%np))

    IF(PRESENT(m_in)) THEN
        IF(.not. PRESENT(cm)) THEN
            PRINT *, 'ERROR: must pass communicator into derivs cell' 
            PRINT *, 'if passing parallel indices'
            stop
        ENDIF
        IF(.not. PRESENT(n_in)) THEN
            PRINT *, 'ERROR: must pass both m and n into derivs cell' 
            PRINT *, 'if passing parallel indices'
            stop
        ENDIF
        n_low = n_in(1)
        n_up  = n_in(2)
        m_low = m_in(1)
        m_up  = m_in(2)
    ELSE
        n_low = 0
        n_up  = p
        m_low = 0
        m_up  = p
    ENDIF

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
    
!   Characteristic grid spacing
    cell%h = SQRT(cell%SA()/cell%info%Y%nt**2D0)
    
    it = n_low*n_low + m_low + n_low
!   Loop over harmonics, but only up to coarse grid order, since the rest are 0
    DO n = n_low, n_up
        ih = n + 1

!       Values at current order        
        nm =>Y%nm(ih)

        IF(n .eq. n_low) THEN
            ml_cur = m_low
            im = ml_cur + n
        ELSE
            ml_cur = -n
            im = 0
        ENDIF

        IF(n .eq. n_up)  THEN
            mu_cur = m_up
        ELSE
            mu_cur =  n
        ENDIF

        DO m = ml_cur, mu_cur ! Could exploit symmetry here... but it isn't a big deal since it takes so little time
            
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

!   Reduce above
    IF(PRESENT(cm)) THEN
        cell%xf     = cm%reduce(cell%xf)
        cell%dxt    = cm%reduce(cell%dxt)
        cell%dxp    = cm%reduce(cell%dxp)
        cell%dxt2   = cm%reduce(cell%dxt2)
        cell%dxp2   = cm%reduce(cell%dxp2)
        cell%dxtp   = cm%reduce(cell%dxtp)
        cell%dxt3   = cm%reduce(cell%dxt3)
        cell%dxt2p  = cm%reduce(cell%dxt2p)
        cell%dxtp2  = cm%reduce(cell%dxtp2)
        cell%dxp3   = cm%reduce(cell%dxp3)
        cell%dxt4   = cm%reduce(cell%dxt4)
        cell%dxt3p  = cm%reduce(cell%dxt3p)
        cell%dxt2p2 = cm%reduce(cell%dxt2p2)
        cell%dxtp3  = cm%reduce(cell%dxtp3)
        cell%dxp4   = cm%reduce(cell%dxp4)
    ENDIF

    nm => NULL()
    Y => NULL()
END SUBROUTINE Derivscell

! -------------------------------------------------------------------------!
! Gets the force jump at the locations on the cell, de-aliased!
SUBROUTINE Stresscell(cell, thti, cm)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell
    TYPE(YType), POINTER :: Y
    TYPE(cmType), INTENT(IN), POINTER, OPTIONAL :: cm
    INTEGER :: i, j, q, rc, it1, it2
    INTEGER, INTENT(IN), OPTIONAL :: thti(2)
    REAL(KIND = 8), ALLOCATABLE :: nks(:,:,:)
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
    
    B = cell%B
    C = cell%C
    Eb = cell%Eb
    c0 = cell%c0

    Y => cell%info%Yf
    ALLOCATE(nks(3, Y%nt, Y%np))
    
    IF(PRESENT(thti)) THEN
        IF(.not. PRESENT(cm)) THEN
            PRINT *, 'ERROR: must pass communicator into stress cell' 
            PRINT *, 'if passing parallel indices'
            stop
        ENDIF
        it1 = thti(1)
        it2 = thti(2)
    ELSE
        it1 = 1
        it2 = Y%nt
    ENDIF

    nks = 0D0
    cell%ff = 0D0
!   Big loop over all points in grid
    DO i = it1, it2
        inner:DO j = 1, Y%np
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
            
            nks(:,i,j) = -nk
            
        ENDDO inner
    ENDDO

    IF(PRESENT(thti)) THEN ! No clue why the generic isn't working here
        nks = cm%REDUCER3(nks)
        cell%ff = cm%REDUCER3(cell%ff)
    ENDIF

!   Now we need to filter for anti-aliasing. Do the transform of the force into spectral
!   space, cut the highest modes, and transform back into physical space
    IF(cell%init) THEN
        cell%fmn(1,:) = Y%forward(cell%ff(1,:,:)*cell%J/SIN(Y%th), cell%info%q)
        cell%fmn(2,:) = Y%forward(cell%ff(2,:,:)*cell%J/SIN(Y%th), cell%info%q)
        cell%fmn(3,:) = Y%forward(cell%ff(3,:,:)*cell%J/SIN(Y%th), cell%info%q)

!       Normal vector for volume correction
        cell%nkmn(1,:) = Y%forward(nks(1,:,:), cell%info%q) 
        cell%nkmn(2,:) = Y%forward(nks(2,:,:), cell%info%q)
        cell%nkmn(3,:) = Y%forward(nks(3,:,:), cell%info%q)

!       Normal and area for fluid
        cell%nkt(1,:) = Y%forward(nks(1,:,:)*cell%J/SIN(Y%th), cell%info%q) 
        cell%nkt(2,:) = Y%forward(nks(2,:,:)*cell%J/SIN(Y%th), cell%info%q)
        cell%nkt(3,:) = Y%forward(nks(3,:,:)*cell%J/SIN(Y%th), cell%info%q)

!       Area/sin(theta)
        cell%Jtmn = Y%forward(cell%J/SIN(Y%th), cell%info%q)
    ENDIF
END SUBROUTINE Stresscell

! -------------------------------------------------------------------------!
! Runs until initial cell is relaxed
SUBROUTINE RelaxCell(cell, tol)
    CLASS(cellType), INTENT(INOUT) :: cell
    TYPE(sharedType), POINTER :: info
    REAL(KIND = 8), INTENT(IN) :: tol
    COMPLEX(KIND = 8), ALLOCATABLE :: A2(:,:), b2(:), ut(:), wrk(:)
    INTEGER iter, p_info, i
    REAL(KIND = 8), ALLOCATABLE :: rwrk(:)
    COMPLEX, ALLOCATABLE :: swrk(:)
    INTEGER, ALLOCATABLE :: IPIV(:)

    print *, 'WARNING: Relax cell not supported right now for GMRES'
    RETURN

    info => cell%info

    ALLOCATE(ut(info%Nmat), &
             IPIV(info%Nmat), wrk(info%Nmat), &
             swrk(info%Nmat*(info%Nmat+1)), &
             rwrk(info%Nmat))

    cell%umn(1,1) = 1/cell%Ca
    DO WHILE(MAXVAL(ABS(cell%umn))*cell%Ca .gt. tol)
            CALL cell%derivs()
            CALL cell%stress() 
            ! CALL cell%fluid(A2, b2, .false.)
            CALL zcgesv(info%Nmat, 1, A2, info%Nmat, IPIV, b2, info%Nmat, ut, info%Nmat, wrk, swrk, rwrk, iter, p_info)
            cell%umn = 0D0
            cell%umn(1,1:info%Nmat/3) = ut((/(i, i=1,info%Nmat-2, 3)/))
            cell%umn(2,1:info%Nmat/3) = ut((/(i, i=2,info%Nmat-1, 3)/))
            cell%umn(3,1:info%Nmat/3) = ut((/(i, i=3,info%Nmat  , 3)/))
            cell%umn(:,1) = 0D0
            cell%xmn = cell%xmn + cell%umn*info%dt
            cell%x(1,:,:) = info%Y%backward(cell%xmn(1,:))
            cell%x(2,:,:) = info%Y%backward(cell%xmn(2,:)) 
            cell%x(3,:,:) = info%Y%backward(cell%xmn(3,:))
            write(*,'(F8.6)') MAXVAL(ABS(cell%umn))*cell%Ca
    ENDDO
END SUBROUTINE RelaxCell

! -------------------------------------------------------------------------!
! Takes an input and cell and de-alises it
FUNCTION DealiasCell(cell, f) RESULT(fc)
    CLASS(cellType), INTENT(IN) :: cell
    REAL(KIND = 8), INTENT(IN) :: f(:,:)
    REAL(KIND = 8), ALLOCATABLE :: fc(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:)

    ALLOCATE(fc(cell%info%Y%nt, cell%info%Y%np))
!   Check to make sure f is the right size
    IF(SIZE(f) .ne. cell%info%Yf%nt*cell%info%Yf%np) THEN
        print *, 'ERROR: trying to de-alias something of wrong size'
        STOP
    ENDIF

    fmn = cell%info%Yf%forward(f, cell%info%q)
    fc = cell%info%Yf%backward(fmn, cell%info%p)
END FUNCTION

! -------------------------------------------------------------------------!
! Calculate the surface area of a cell
FUNCTION SAcell(cell) RESULT(SA)
    CLASS(cellType), INTENT(IN) :: cell
    REAL(KIND = 8) :: SA
    INTEGER :: i, j
    SA = 0D0

!   Basically the integrand is just the infinitesimal area element
    DO i = 1,cell%info%Yf%nt
        DO j = 1,cell%info%Yf%np
!           Integrate via Gauss quad
            SA = SA + cell%info%Yf%wg(i)*cell%J(i,j)*cell%info%Yf%dphi/sin(cell%info%Yf%tht(i))
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

    DO i = 1,cell%info%Yf%nt
        DO j = 1,cell%info%Yf%np
!           Integrand
            intgd = x(3,i,j)*(xt(1,i,j)*xp(2,i,j) - xp(1,i,j)*xt(2,i,j))

!           Integrate via Gauss quad
            V = V + intgd* &
            cell%info%Yf%wg(i)*cell%info%Yf%dphi/sin(cell%info%Yf%tht(i))
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
    Y => cell%info%Yf
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
! Checks if the cell is in the primary cell, and gets the image that is if not
SUBROUTINE InPrimcell(cell)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell
    REAL(KIND = 8) :: xc(3), bv1(3), bv2(3), bv3(3)
    INTEGER :: move(3)
    TYPE(sharedType), POINTER:: info

    info => cell%info

!   Get center of cell
    xc = REAL(cell%xmn(:,1)*0.5D0*ispi)

!   Get number of boxes to move based on partial coordinates
    move = FLOOR(bvPartial(info%bv, xc))
    bv1 = info%bv(:,1)
    bv2 = info%bv(:,2)
    bv3 = info%bv(:,3)

!   Now move them back the specified number of cells, accounting for spherical hamronics
    cell%xmn(:,1) = cell%xmn(:,1) &
                  - (REAL(move(1))*bv1)/(0.5D0*ispi) &
                  - (REAL(move(2))*bv2)/(0.5D0*ispi) &
                  - (REAL(move(3))*bv3)/(0.5D0*ispi)

    cell%x(1,:,:) = info%Y%backward(cell%xmn(1,:))
    cell%x(2,:,:) = info%Y%backward(cell%xmn(2,:))
    cell%x(3,:,:) = info%Y%backward(cell%xmn(3,:))

END SUBROUTINE InPrimcell

! -------------------------------------------------------------------------!
! Calculates LHS given a vec, e.g. vel, assumed on intg grid
! Composed of two components: double layer integral & essentially identity mat
SUBROUTINE LHS_realCell(cell, v_input, v_input_mn, v, celli, &
                        periodic_in, dbg, itt1, itt2)
    CLASS(cellType), INTENT(IN), TARGET :: cell
    CLASS(cellType), INTENT(IN), TARGET, OPTIONAL :: celli
    TYPE(sharedType), POINTER:: info
    TYPE(YType), POINTER :: Y, Ys, Yf
    REAL(KIND = 8), ALLOCATABLE, INTENT(OUT) :: v(:)
    REAL(KIND = 8), ALLOCATABLE, INTENT(IN ), OPTIONAL :: v_input(:,:,:)
    REAL(KIND = 8), ALLOCATABLE :: xcg(:,:,:), nJt(:,:,:), urot(:,:,:), &
                                   ft(:), ft2(:,:), wgi(:), v_in(:,:,:), nJuR(:,:,:)
    REAL(KIND = 8) :: xcr(3), r(3), dphi, rcur(3), u_pt(3), v_tmp(3), minr, &
                      rt(3), rho, t, nn, tht_rot, phi_rot
    COMPLEX(KIND = 8), ALLOCATABLE, INTENT(IN), OPTIONAL :: v_input_mn(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: xmnR(:,:), nmnR(:,:), umnR(:,:), v_in_mn(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: periodic_in, dbg
    LOGICAL, ALLOCATABLE :: vrecon(:, :), xrecon(:,:), nJrecon(:,:)
    LOGICAL ::  sing, periodic, dbgflg=.false.
    INTEGER, INTENT(IN), OPTIONAL :: itt1, itt2
    INTEGER :: i, j, nt, np, i2, j2, ic, row, iter, &
               iper, jper, kper, indi, indj, it1, it2, cid
    INTEGER(KIND = 8) :: tic, toc, rate

    info => cell%info
    Y    => info%Y

!   ====== Check optional inputs ======
    IF(PRESENT(dbg)) dbgflg = dbg

    IF(PRESENT(periodic_in)) THEN
        periodic = periodic_in
    ELSE
        periodic = info%periodic
    ENDIF

    IF(PRESENT(v_input)) THEN
        v_in = v_input
        ALLOCATE(v_in_mn(3, (info%p+1)*(info%p+1)))
!       First, calculate the spherical harmonic coefficients of the input vector
!       (this is typically either a residual or a velocity vector)
        v_in_mn(1,:) = Y%forward(v_in(1,:,:))
        v_in_mn(2,:) = Y%forward(v_in(2,:,:))
        v_in_mn(3,:) = Y%forward(v_in(3,:,:))
    ELSEIF(PRESENT(v_input_mn)) THEN
        v_in_mn = v_input_mn
        ALLOCATE(v_in(3, Y%nt, Y%np))
!       Instead calculate points on the grid
        v_in(1,:,:) = Y%backward(v_in_mn(1,:))
        v_in(2,:,:) = Y%backward(v_in_mn(2,:))
        v_in(3,:,:) = Y%backward(v_in_mn(3,:))
    ELSEIF(.not. (PRESENT(v_input_mn) .or. PRESENT(v_input))) THEN
        print *, 'Need at least one input into LHS_real'
        stop
    ELSE
        v_in_mn = v_input_mn
        v_in = v_input
    ENDIF

    IF(PRESENT(itt1)) THEN
        it1 = itt1
    ELSE
        it1 = 1
    ENDIF

    IF(PRESENT(itt2)) THEN
        it2 = itt2
    ELSE
        it2 = Y%nt
    ENDIF

    ALLOCATE(&
    xcg(3, Y%nt, Y%np), &
    nJt(3, Y%nt, Y%np), &
    urot(3, Y%nt, Y%np), &
    ft(3), ft2(3,3), &
    nmnR(3, (info%p+1)*(info%p+1)), &
    xmnR(3, (info%p+1)*(info%p+1)), &
    umnR(3, (info%p+1)*(info%p+1)), &
    wgi(Y%nt), &
    v(3*Y%nt*Y%np), &
    vrecon(3, Y%p + 1), &
    xrecon(3, Y%p + 1), &
    nJrecon(3, Y%p + 1))

    vrecon(1,:) = Y%reconConsts(v_in_mn(1,:))
    vrecon(2,:) = Y%reconConsts(v_in_mn(2,:))
    vrecon(3,:) = Y%reconConsts(v_in_mn(3,:))

    IF(PRESENT(celli)) THEN
        nJrecon =celli%nJrecon
        xrecon = celli%xrecon
        Yf => info%Yf
        cid = celli%id
    ELSE
        nJrecon =cell%nJrecon
        xrecon = cell%xrecon
        cid = cell%id
    ENDIF

!   Want normals on unrotated grid if I might do integration on unrotated grid
    IF(PRESENT(celli) .or. periodic) THEN
        ALLOCATE(nJur(3,Y%nt, Y%np))
    
        IF(PRESENT(celli)) THEN
            nmnR(1,:) = celli%nkt(1,1:(info%p+1)*(info%p+1))
            nmnR(2,:) = celli%nkt(2,1:(info%p+1)*(info%p+1))
            nmnR(3,:) = celli%nkt(3,1:(info%p+1)*(info%p+1))
        ELSE
            nmnR(1,:) = cell%nkt(1,1:(info%p+1)*(info%p+1))
            nmnR(2,:) = cell%nkt(2,1:(info%p+1)*(info%p+1))
            nmnR(3,:) = cell%nkt(3,1:(info%p+1)*(info%p+1))
        ENDIF
        nJur(1,:,:) = Y%backward(nmnR(1,:), info%p, nJrecon(1,:))
        nJur(2,:,:) = Y%backward(nmnR(2,:), info%p, nJrecon(2,:))
        nJur(3,:,:) = Y%backward(nmnR(3,:), info%p, nJrecon(3,:))
    ENDIF

    v  = 0D0
    ic = (it1 - 1)*Y%np
    iter = 0
!   Loop over all of the integration points now and calculate in the standard way
    DO i = it1, it2
        DO j = 1,Y%np
            xcr = cell%x(:, i, j)
            v_tmp = 0D0
            sing = .false.
            ic = ic + 1
            row = 3*(ic - 1)  + 1

!           Conditinal statements: integration method varies depending on problem
!           parameters (periodic,. mult cells, etc.)
!           Construct grid appropriate to problem
            IF(PRESENT(celli)) THEN

!               If min spacing is small, we need to do near-singular integration
                IF(cell%sing_inds(1, i, j, celli%id) .gt. 0) THEN
                    indi = cell%sing_inds(1, i, j, celli%id)
                    indj = cell%sing_inds(2, i, j, celli%id)
                    tht_rot = Yf%tht(indi)
                    phi_rot = Yf%phi(indj)
                    
                    sing = .true.
                    Ys => info%Ys

!                   Need to integrate on finer grid
                    nt = Ys%nt
                    np = Ys%np

!                   Deallocate integ. quants
                    DEALLOCATE(xcg, nJt, urot, wgi)
                    ALLOCATE( &
                    xcg(3,  nt, np), &
                    urot(3, nt, np), &
                    nJt(3,  nt, np), &
                    wgi(nt))

                    dphi = Ys%dphi

!                   Manage the additional prefactors stemming from the integrals
                    wgi = Ys%wg*info%h*(-info%k)*COSH(info%k*info%xsf - info%k)

!                   We don't do cosine transformation to cluster points near near-singularity, mult sine back in
                    wgi = wgi*SIN(Ys%th(:,1))

!                   Rotate about nearest point to projected singularity
                    xmnR(1,:) = Yf%rotate(celli%xmn(1,:), indi, indj, -Yf%phi(indj))
                    xmnR(2,:) = Yf%rotate(celli%xmn(2,:), indi, indj, -Yf%phi(indj))
                    xmnR(3,:) = Yf%rotate(celli%xmn(3,:), indi, indj, -Yf%phi(indj))

!                   These backward trnsforms are quite expensive
                    xcg(1,:,:) = Ys%backward(xmnR(1,:), info%p, xrecon(1,:))
                    xcg(2,:,:) = Ys%backward(xmnR(2,:), info%p, xrecon(2,:))
                    xcg(3,:,:) = Ys%backward(xmnR(3,:), info%p, xrecon(3,:))

!                   Rotate the normal vector total constants
                    nmnR(1,:) = Yf%rotate(celli%nkt(1,1:(info%p+1)*(info%p+1)), indi, indj, -Yf%phi(indj))
                    nmnR(2,:) = Yf%rotate(celli%nkt(2,1:(info%p+1)*(info%p+1)), indi, indj, -Yf%phi(indj))
                    nmnR(3,:) = Yf%rotate(celli%nkt(3,1:(info%p+1)*(info%p+1)), indi, indj, -Yf%phi(indj))

                    nJt(1,:,:) = Ys%backward(nmnR(1,:), info%p, nJrecon(1,:))
                    nJt(2,:,:) = Ys%backward(nmnR(2,:), info%p, nJrecon(2,:))
                    nJt(3,:,:) = Ys%backward(nmnR(3,:), info%p, nJrecon(3,:))

!                   Velocities on rotated grid
                    umnR(1,:) = Yf%rotate(v_in_mn(1,:), indi, indj, -Yf%phi(indj))
                    umnR(2,:) = Yf%rotate(v_in_mn(2,:), indi, indj, -Yf%phi(indj))
                    umnR(3,:) = Yf%rotate(v_in_mn(3,:), indi, indj, -Yf%phi(indj))

                    urot(1,:,:) = Ys%backward(umnR(1,:), info%p, vrecon(1,:))
                    urot(2,:,:) = Ys%backward(umnR(2,:), info%p, vrecon(2,:))
                    urot(3,:,:) = Ys%backward(umnR(3,:), info%p, vrecon(3,:))

!               Well-separated, normal grid/integration
                ELSE
                    sing = .false.

!                   We can use the coarse grid
                    nt = Y%nt
                    np = Y%np

!                   Deallocate integ. quants
                    DEALLOCATE(urot, xcg, nJt, wgi)
                    ALLOCATE( &
                    xcg(3,  nt, np), &
                    urot(3,  nt, np), &
                    nJt(3,  nt, np), &
                    wgi(nt))
                    
                    dphi = Y%dphi
                    wgi  = Y%wg

                    xcg = celli%x
                    urot = v_in
                    nJt = nJur
                ENDIF

            ELSE
!               For periodic, the singular integration doesn't work and we must use near-sing
                IF(periodic) THEN                    
                    sing = .true.
                    Ys => info%Ys
                    tht_rot = Y%tht(i)
                    phi_rot = Y%phi(j)

!                   Need to integrate on finer grid
                    nt = Ys%nt
                    np = Ys%np

!                   Deallocate integ. quants
                    DEALLOCATE(xcg, nJt, urot, wgi)
                    ALLOCATE( &
                    xcg(3,  nt, np), &
                    urot(3, nt, np), &
                    nJt(3,  nt, np), &
                    wgi(nt))

                    dphi = Ys%dphi

!                   Manage the additional prefactors stemming from the integrals
                    wgi = Ys%wg*info%h*(-info%k)*COSH(info%k*info%xsf - info%k)

!                   We don't do cosine transformation to cluster points near near-singularity, mult sine back in
                    wgi = wgi*SIN(Ys%th(:,1))

!                   Rotate about nearest point to projected singularity
                    xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
                    xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
                    xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

!                   These backward trnsforms are quite expensive
                    xcg(1,:,:) = Ys%backward(xmnR(1,:), info%p, xrecon(1,:))
                    xcg(2,:,:) = Ys%backward(xmnR(2,:), info%p, xrecon(2,:))
                    xcg(3,:,:) = Ys%backward(xmnR(3,:), info%p, xrecon(3,:))

!                   Rotate the normal vector total constants
                    nmnR(1,:) = Y%rotate(cell%nkt(1,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    nmnR(2,:) = Y%rotate(cell%nkt(2,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    nmnR(3,:) = Y%rotate(cell%nkt(3,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))

                    nJt(1,:,:) = Ys%backward(nmnR(1,:), info%p, nJrecon(1,:))
                    nJt(2,:,:) = Ys%backward(nmnR(2,:), info%p, nJrecon(2,:))
                    nJt(3,:,:) = Ys%backward(nmnR(3,:), info%p, nJrecon(3,:))

!                   Velocities on rotated grid
                    umnR(1,:) = Y%rotate(v_in_mn(1,:), i, j, -Y%phi(j))
                    umnR(2,:) = Y%rotate(v_in_mn(2,:), i, j, -Y%phi(j))
                    umnR(3,:) = Y%rotate(v_in_mn(3,:), i, j, -Y%phi(j))

                    urot(1,:,:) = Ys%backward(umnR(1,:), info%p, vrecon(1,:))
                    urot(2,:,:) = Ys%backward(umnR(2,:), info%p, vrecon(2,:))
                    urot(3,:,:) = Ys%backward(umnR(3,:), info%p, vrecon(3,:))
!               Fully singular integration on same cell: Get rotated constants
                ELSE
                    sing = .false.

!                   For integration (changes if need a finer grid)
                    dphi = Y%dphi

!                   Number intg points to loop over for this case
                    nt = Y%nt
                    np = Y%np

                    DEALLOCATE(urot, xcg, nJt, wgi)
                    ALLOCATE( &
                    xcg(3,  nt, np), &
                    urot(3,  nt, np), &
                    nJt(3,  nt, np), &
                    wgi(nt))

                    xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
                    xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
                    xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

!                   Rotated integration points in unrotated frame
                    xcg(1,:,:) = Y%backward(xmnR(1,:), info%p, xrecon(1,:))
                    xcg(2,:,:) = Y%backward(xmnR(2,:), info%p, xrecon(2,:))
                    xcg(3,:,:) = Y%backward(xmnR(3,:), info%p, xrecon(3,:))

!                   Rotate the normal vector total constants
                    nmnR(1,:) = Y%rotate(cell%nkt(1,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    nmnR(2,:) = Y%rotate(cell%nkt(2,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    nmnR(3,:) = Y%rotate(cell%nkt(3,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))

                    nJt(1,:,:) = Y%backward(nmnR(1,:), info%p, nJrecon(1,:))
                    nJt(2,:,:) = Y%backward(nmnR(2,:), info%p, nJrecon(2,:))
                    nJt(3,:,:) = Y%backward(nmnR(3,:), info%p, nJrecon(3,:))

!                   Velocities on rotated grid
                    umnR(1,:) = Y%rotate(v_in_mn(1,:), i, j, -Y%phi(j))
                    umnR(2,:) = Y%rotate(v_in_mn(2,:), i, j, -Y%phi(j))
                    umnR(3,:) = Y%rotate(v_in_mn(3,:), i, j, -Y%phi(j))

                    urot(1,:,:) = Y%backward(umnR(1,:), info%p, vrecon(1,:))
                    urot(2,:,:) = Y%backward(umnR(2,:), info%p, vrecon(2,:))
                    urot(3,:,:) = Y%backward(umnR(3,:), info%p, vrecon(3,:))

                    wgi = Y%ws
                ENDIF
            ENDIF

!           Now perform the actual integrations
            DO j2 = 1,np
                INNER:DO i2 = 1,nt

!                   Use periodic Greens (or not)
                    IF(periodic) THEN
!                       Three cases need consideration here: singular (integ on same cell), near-singular, non-singular
!                       Each has unique ways that these integrals need to be calc'd

!                       We checked if this point is in the cutoff range in RHS,
!                       cycle if not
                        IF(cell%cut_box(1, i, j, i2, j2, cid) .eq. 2) CYCLE INNER

!                       Bring to nearest distance
                        rt = xcg(:,i2,j2) - xcr

                        IF(PRESENT(celli))&
                        rt = rt &
                           - cell%cut_box(1, i, j, i2, j2, cid)*info%bv(:,1) &
                           - cell%cut_box(2, i, j, i2, j2, cid)*info%bv(:,2) &
                           - cell%cut_box(3, i, j, i2, j2, cid)*info%bv(:,3)

                        ! IF(NORM2(rt)*info%xi.gt.350) THEN
                        IF(NORM2(rt)*info%xi.gt.3.5) THEN
                            CYCLE INNER
                        ELSE
!                           This intg is taking place on the patch,
                            ft2 = PTij(rt, 0, info%bv, nJt(:,i2,j2), info%eye, info%xi)
                            ! ft2 = PTij(rt, 4, info%bv, nJt(:,i2,j2), info%eye, .3545d0, fourier = .true.)
                        ENDIF
                    ELSE
                        r = xcg(:,i2,j2) - xcr
                        ft2 = Tij(r, nJt(:,i2,j2))
                    ENDIF
                    IF(sing) ft2 = ft2*info%n(i2)
                    ft = urot(:,i2,j2)

                    v_tmp = v_tmp + INNER3_33(ft,ft2)*dphi*wgi(i2)
                ENDDO INNER
            ENDDO

!           There is a linear periodic part in the stresslet that enforces 0 mean
!           flow through the box. It is non-singular, but it requires integration
!           over the whole cell, so I perform this integration here
!           This is also applicable to the (near-) singular intg, where
!           the partitioned, nonsingular fn must be intgd as well
            IF(periodic .or. sing) THEN
                DEALLOCATE(xcg, nJt, urot, wgi)
                ALLOCATE( &
                xcg(3,  Y%nt, Y%np), &
                urot(3, Y%nt, Y%np), &
                nJt(3,  Y%nt, Y%np), &
                wgi(Y%nt))

                IF(     PRESENT(celli)) xcg = celli%x
                IF(.not.PRESENT(celli)) xcg = cell%x
                nJt = nJur
                urot = v_in
                dphi = Y%dphi
                wgi  = Y%wg

                DO i2 = 1, Y%nt
                    DO j2 = 1, Y%np
                        ft2 = 0D0
                        r = xcg(:,i2,j2) - xcr
                        IF(periodic) ft2 = -8D0*PI/info%tau*OUTER(nJt(:,i2,j2), r)

!                       Partition fn calculation
                        IF(sing)THEN
!                           Smallest angle along unit sphere
                            rho = ACOS(COS(Y%tht(i2))*COS(tht_rot) &
                                +      SIN(Y%tht(i2))*SIN(tht_rot) &
                                *      COS(ABS(Y%phi(j2) -phi_rot)))
                            t = rho/info%thtc
                            nn = EXP(2D0*EXP(-1D0/t)/(t - 1D0))

                            IF(t .gt. 1) nn = 0D0
                            IF(periodic) THEN
!                               Deal with limiting value
!                               1 box to cover if celli present and near point in diff box
                                IF(NORM2(r) .gt. EPSILON(r)*1D1) &
                                ft2 = ft2 + PTij(r, 1, info%bv, nJt(:,i2,j2), info%eye, info%xi)*(1D0 - nn)
                            ELSE
                                IF(NORM2(r) .gt. EPSILON(r)*1D1) &
                                ft2 = Tij(r, nJt(:,i2,j2))*(1D0 - nn)
                            ENDIF
                        ENDIF
                        ft = urot(:,i2,j2)
                        v_tmp = v_tmp + INNER3_33(ft,ft2)*dphi*wgi(i2)
                    ENDDO
                ENDDO
            ENDIF

!           Now the integrations are done: perform final operations and add in
            v_tmp = v_tmp*(1D0 - cell%lam)/(4D0*pi*(1D0 + cell%lam))
            IF(PRESENT(celli)) THEN
                u_pt = 0D0
            ELSE
                u_pt = v_in(:, i, j)
            ENDIF
            v(row:row + 2) = u_pt - v_tmp
        ENDDO
    ENDDO

END SUBROUTINE LHS_realCell

! -------------------------------------------------------------------------!
! Calculates RHS given input vec (can be mn consts), e.g. force, on intg grid
! Composed of two components: Single layer integral & background vel
SUBROUTINE RHS_realCell(cell, v_input, v_input_mn, v, celli, periodic_in, &
                        itt1, itt2)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell
    CLASS(cellType), INTENT(IN), TARGET, OPTIONAL :: celli
    TYPE(YType), POINTER :: Y, Ys, Yf
    TYPE(sharedType), POINTER:: info
    REAL(KIND = 8), ALLOCATABLE, INTENT(OUT):: v(:)
    REAL(KIND = 8), ALLOCATABLE, INTENT(IN), OPTIONAL :: v_input(:,:,:)
    REAL(KIND = 8), ALLOCATABLE :: xcg(:,:,:), frot(:,:,:), &
                                   ft(:), ft2(:,:), wgi(:), &
                                   v_in(:,:,:), tht_t(:)
    REAL(KIND = 8) :: xcr(3), r(3), dphi, rn, rcur(3), Uc(3), v_tmp(3), &
                      Utmp(3,3), minr, rt(3), t, rho, nn, tht_rot, phi_rot
    COMPLEX(KIND = 8), ALLOCATABLE, INTENT(IN), OPTIONAL :: v_input_mn(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fmnR(:,:), xmnR(:,:), v_in_mn(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: periodic_in
    LOGICAL sing, periodic
    LOGICAL, ALLOCATABLE :: vrecon(:, :), xrecon(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: itt1, itt2
    INTEGER :: i, j, nt, np, i2, j2, ic, row, &
               iper, jper, kper, indi, indj, it1, it2, cid
    INTEGER(KIND = 8) :: tic, toc, rate

    info => cell%info
    Y    => info%Y

!   ====== Check optional inputs ======
    IF(PRESENT(periodic_in)) THEN
        periodic = periodic_in
    ELSE
        periodic = info%periodic
    ENDIF

    IF(PRESENT(v_input)) THEN
        v_in = v_input
        ALLOCATE(v_in_mn(3, (info%p+1)*(info%p+1)))

!       First, calculate the spherical harmonic coefficients of the input vector
!       (this is typically either a residual or a velocity vector)
        v_in_mn(1,:) = Y%forward(v_in(1,:,:))
        v_in_mn(2,:) = Y%forward(v_in(2,:,:))
        v_in_mn(3,:) = Y%forward(v_in(3,:,:))
    ELSEIF(PRESENT(v_input_mn)) THEN
        v_in_mn = v_input_mn
        ALLOCATE(v_in(3, Y%nt, Y%np))
!       Instead calculate points on the grid
        v_in(1,:,:) = Y%backward(v_in_mn(1,:))
        v_in(2,:,:) = Y%backward(v_in_mn(2,:))
        v_in(3,:,:) = Y%backward(v_in_mn(3,:))
    ELSEIF(.not. (PRESENT(v_input_mn) .or. PRESENT(v_input))) THEN
        print *, 'Need at least one input into LHS_real'
        stop
    ELSE
        v_in_mn = v_input_mn
        v_in = v_input
    ENDIF

    IF(PRESENT(itt1)) THEN
        it1 = itt1
    ELSE
        it1 = 1
    ENDIF

    IF(PRESENT(itt2)) THEN
        it2 = itt2
    ELSE
        it2 = Y%nt
    ENDIF

    ALLOCATE(&
    xcg(3, Y%nt, Y%np), &
    frot(3, Y%nt, Y%np), &
    ft(3), ft2(3,3), &
    fmnR(3, (info%p+1)*(info%p+1)), &
    xmnR(3, (info%p+1)*(info%p+1)), &
    wgi(Y%nt), &
    tht_t(Y%nt), &
    v(3*Y%nt*Y%np), &
    vrecon(3, Y%p + 1), &
    xrecon(3, Y%p + 1))

    vrecon(1,:) = Y%reconConsts(v_in_mn(1,:))
    vrecon(2,:) = Y%reconConsts(v_in_mn(2,:))
    vrecon(3,:) = Y%reconConsts(v_in_mn(3,:))
    
    IF(PRESENT(celli)) THEN
        xrecon = celli%xrecon
        Yf => info%Yf
        cid = celli%id
    ELSE
        xrecon = cell%xrecon
        cid = cell%id
    ENDIF

    Utmp = TRANSPOSE(info%dU)

    v = 0D0
    ic = (it1 - 1)*Y%np
!   Loop over all of the integration points now and calculate in the standard way
    DO i = it1, it2
        DO j = 1,Y%np
            xcr = cell%x(:, i, j)
            v_tmp = 0D0
            ic = ic + 1
            row = 3*(ic - 1)  + 1

!           Conditinal statements: integration method varies depending on problem
!           parameters (periodic,. mult cells, etc.)
            
            IF(PRESENT(celli)) THEN

!               Check minimum spacing between eval point and integration cell
                minr = celli%h + 1D0
                DO i2 = 1, Yf%nt
                    DO j2 = 1, Yf%np
                        r = celli%xf(:,i2,j2) - xcr
                        minr = MIN(NORM2(r), minr)
!                       We need to check all of the periodic images of the cell as well
!                       Just check immediately surrounding boxes the point is close to
                        IF(periodic) THEN
                            DO iper = -1,1
                                DO jper = -1,1
                                    DO kper = -1,1
                                        rcur = r &
                                             - iper*info%bv(:,1) &
                                             - jper*info%bv(:,2) &
                                             - kper*info%bv(:,3) 
                                        minr = MIN(NORM2(rcur), minr)
                                        IF(minr .eq. NORM2(rcur)) THEN
                                            indi = i2
                                            indj = j2
                                            cell%sing_inds(1, i, j, celli%id) = indi
                                            cell%sing_inds(2, i, j, celli%id) = indj
                                        ENDIF
                                    ENDDO
                                ENDDO
                            ENDDO
                        ENDIF
!                       Save indices of min spacing, as well as for GMRES its
                        IF(minr .eq. NORM2(r)) THEN
                            indi = i2
                            indj = j2
                            cell%sing_inds(1, i, j, celli%id) = indi
                            cell%sing_inds(2, i, j, celli%id) = indj
                        ENDIF
                    ENDDO
                ENDDO

                IF(minr .lt. celli%h) THEN
                    sing = .true.
                    Ys => info%Ys
                    tht_rot = Yf%tht(indi)
                    phi_rot = Yf%phi(indj)

!                   Need to integrate on finer grid
                    nt = Ys%nt
                    np = Ys%np

!                   Deallocate integ. quants
                    DEALLOCATE(frot, xcg, wgi, tht_t)
                    ALLOCATE( &
                    frot(3, nt, np), &
                    xcg(3,  nt, np), &
                    wgi(nt), &
                    tht_t(nt))

                    dphi = Ys%dphi

!                   Manage the additional prefactors stemming from the integrals
                    wgi = Ys%wg*info%h*(-info%k)*COSH(info%k*info%xsf - info%k)

!                   We don't do cosine transformation to cluster points near near-singularity, mult sine back in
                    wgi = wgi*SIN(info%Ys%th(:,1))
                    tht_t = info%Ys%th(:,1)

!                   Rotate about nearest point to projected singularity
                    xmnR(1,:) = Yf%rotate(celli%xmn(1,:), indi, indj, -Yf%phi(indj))
                    xmnR(2,:) = Yf%rotate(celli%xmn(2,:), indi, indj, -Yf%phi(indj))
                    xmnR(3,:) = Yf%rotate(celli%xmn(3,:), indi, indj, -Yf%phi(indj))

!                   These backward trnsforms are quite expensive
                    xcg(1,:,:) = Ys%backward(xmnR(1,:), info%p, xrecon(1,:))
                    xcg(2,:,:) = Ys%backward(xmnR(2,:), info%p, xrecon(2,:))
                    xcg(3,:,:) = Ys%backward(xmnR(3,:), info%p, xrecon(3,:))

!                   Forces on rotated grid
                    fmnR(1,:) = Yf%rotate(v_in_mn(1,1:(info%p+1)*(info%p+1)), indi, indj, -Yf%phi(indj))
                    fmnR(2,:) = Yf%rotate(v_in_mn(2,1:(info%p+1)*(info%p+1)), indi, indj, -Yf%phi(indj))
                    fmnR(3,:) = Yf%rotate(v_in_mn(3,1:(info%p+1)*(info%p+1)), indi, indj, -Yf%phi(indj))

                    frot(1,:,:) = Ys%backward(fmnR(1,:), info%p)
                    frot(2,:,:) = Ys%backward(fmnR(2,:), info%p)
                    frot(3,:,:) = Ys%backward(fmnR(3,:), info%p)

!               Well-separated, normal grid/integration
                ELSE
                    sing = .false.
                    cell%sing_inds(1, i, j, celli%id) = -1

!                   We can use the coarse grid
                    nt = Y%nt
                    np = Y%np

!                   Deallocate integ. quants
                    DEALLOCATE(frot, xcg, wgi, tht_t)
                    ALLOCATE( &
                    frot(3, nt, np), &
                    xcg(3,  nt, np), &
                    wgi(nt), &
                    tht_t(nt))
                    
                    dphi = Y%dphi
                    wgi  = Y%wg
                    tht_t= Y%tht

                    xcg = celli%x
                    frot = v_in
                ENDIF
            ELSE
                IF(periodic) THEN
                    sing = .true.
                    Ys => info%Ys
                    tht_rot = Y%tht(i)
                    phi_rot = Y%phi(j)

!                   Need to integrate on finer grid
                    nt = Ys%nt
                    np = Ys%np

!                   Deallocate integ. quants
                    DEALLOCATE(frot, xcg, wgi, tht_t)
                    ALLOCATE( &
                    frot(3, nt, np), &
                    xcg(3,  nt, np), &
                    wgi(nt), &
                    tht_t(nt))

                    dphi = Ys%dphi

!                   Manage the additional prefactors stemming from the integrals
                    wgi = Ys%wg*info%h*(-info%k)*COSH(info%k*info%xsf - info%k)

!                   We don't do cosine transformation to cluster points near near-singularity, mult sine back in
                    wgi = wgi*SIN(Ys%th(:,1))
                    tht_t = Ys%th(:,1)

!                   Rotate about nearest point to projected singularity
                    xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
                    xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
                    xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

!                   These backward trnsforms are quite expensive
                    xcg(1,:,:) = Ys%backward(xmnR(1,:), info%p, xrecon(1,:))
                    xcg(2,:,:) = Ys%backward(xmnR(2,:), info%p, xrecon(2,:))
                    xcg(3,:,:) = Ys%backward(xmnR(3,:), info%p, xrecon(3,:))

!                   Forces on rotated grid
                    fmnR(1,:) = Y%rotate(v_in_mn(1,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    fmnR(2,:) = Y%rotate(v_in_mn(2,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    fmnR(3,:) = Y%rotate(v_in_mn(3,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))

                    frot(1,:,:) = Ys%backward(fmnR(1,:), info%p)
                    frot(2,:,:) = Ys%backward(fmnR(2,:), info%p)
                    frot(3,:,:) = Ys%backward(fmnR(3,:), info%p)

!               Fully singular integration on same cell: Get rotated constants
                ELSE
                    sing = .false.

!                   For integration (changes if need a finer grid)
                    dphi = Y%dphi

!                   Number intg points to loop over for this case
                    nt = Y%nt
                    np = Y%np
                    DEALLOCATE(frot, xcg, wgi)
                    ALLOCATE( &
                    frot(3, nt, np), &
                    xcg(3,  nt, np), &
                    wgi(nt))

                    xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
                    xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
                    xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

!                   Rotated integration points in unrotated frame
                    xcg(1,:,:) = Y%backward(xmnR(1,:), info%p, xrecon(1,:))
                    xcg(2,:,:) = Y%backward(xmnR(2,:), info%p, xrecon(2,:))
                    xcg(3,:,:) = Y%backward(xmnR(3,:), info%p, xrecon(3,:))

!                   Forces on rotated grid
                    fmnR(1,:) = Y%rotate(v_in_mn(1,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    fmnR(2,:) = Y%rotate(v_in_mn(2,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
                    fmnR(3,:) = Y%rotate(v_in_mn(3,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))

                    frot(1,:,:) = Y%backward(fmnR(1,:), info%p, vrecon(1,:))
                    frot(2,:,:) = Y%backward(fmnR(2,:), info%p, vrecon(2,:))
                    frot(3,:,:) = Y%backward(fmnR(3,:), info%p, vrecon(3,:))

                    wgi = Y%ws
                ENDIF
            ENDIF

            Uc = INNER3_33(xcr, Utmp)
!           Now perform the actual integrations
            DO j2 = 1,np
                INNER:DO i2 = 1,nt
!                   Need an exception for near-singular int for if we're in mis-matched grids
                    r = xcg(:,i2,j2) - xcr
                    rn = NORM2(r)

!                   Use periodic Greens (or not)
                    IF(periodic) THEN
!                       Three cases need consideration here: singular (integ on same cell), near-singular, non-singular
!                       Each has unique ways that these integrals need to be calc'd

                        rt = r
                        IF(PRESENT(celli)) THEN
!                           Check over surrounding boxes to see if there's a nearby point
                            minr = rn + 1D0
                            DO iper = -1,1
                                DO jper = -1,1
                                    DO kper = -1,1
                                        rcur = r &
                                             - iper*info%bv(:,1) &
                                             - jper*info%bv(:,2) &
                                             - kper*info%bv(:,3)
                                        IF(NORM2(rcur) .le. minr) THEN
                                            rt = rcur
                                            minr = NORM2(rcur)
                                            cell%cut_box(:, i ,j ,i2 ,j2 , cid) = (/iper,jper,kper/)
                                        ENDIF
                                    ENDDO
                                ENDDO
                            ENDDO
                        ENDIF
                            
                        IF(NORM2(rt)*info%xi.gt.3.5) THEN
                        ! IF(NORM2(rt)*info%xi.gt.350) THEN
                            cell%cut_box(:, i ,j ,i2 ,j2 , cid) = 2
                            CYCLE INNER
                        ELSE
                            ft2 = PGij(rt, 0, info%bv, info%eye, info%xi)
                            ! ft2 = PGij(rt, 4, info%bv, info%eye, .3454d0, fourier=.true.)
                        ENDIF

!                       We need to check if the periodic images give short range cell-cell interactions.
!                       This needs to be added separately, because it uses the normal Stokeslet, not the periodic.
!                       Go to each of the surrounding boxes, and check if the image point is within the cutoff distance
!                       If it is, add it directly to b with non-periodic Green's function
!                   !!!!This is incorrect right now: Needs to be multiplied by area
                        !!!!!
                        IF(info%CellCell .and. PRESENT(celli)) THEN
                            IF(sing) THEN
                                v_tmp = v_tmp + PeriodicCellCell(info, r)*wgi(i2)*dphi/SIN(tht_t(i2))!*NORM2(nJt(:,i2,j2))
                            ELSE
                                v_tmp = v_tmp + PeriodicCellCell(info, r)*wgi(i2)*dphi!*NORM2(nJt(:,i2,j2))
                            ENDIF
                        ENDIF
                    ELSE
                        ft2 = Gij(r, info%eye)
!                       Add in Morse potential if cells are close enough
                        IF(PRESENT(celli) .and. (rn .lt. 3D0*info%r0) .and. info%CellCell) THEN 
                            ft = ft + Morse(r, rn, info%D, info%r0, info%Beta)
                            ft = ft +    LJ(r, rn, info%epsi, info%r0)
                            IF(.not.sing) ft = ft/SIN(tht_t(i2))
                        ENDIF
                    ENDIF
                    ft  = frot(:,i2,j2)
                    IF(sing) ft2 = ft2*info%n(i2)

                    v_tmp = v_tmp + INNER3_33(ft,ft2)*dphi*wgi(i2)
                ENDDO INNER
            ENDDO

!           If this is (near-) singular and we're doing partition of unity,
!           need to integrate that and add in
            IF(sing) THEN
                DEALLOCATE(xcg, frot, wgi)
                ALLOCATE( &
                xcg(3,  Y%nt, Y%np), &
                frot(3, Y%nt, Y%np), &
                wgi(Y%nt))

                IF(     PRESENT(celli)) xcg = celli%x
                IF(.not.PRESENT(celli)) xcg = cell%x
                frot = v_in
                dphi = Y%dphi
                wgi  = Y%wg

                DO i2 = 1, Y%nt
                    DO j2 = 1, Y%np
                        ft2 = 0D0
                        r = xcg(:,i2,j2) - xcr

                        rho = ACOS(COS(Y%tht(i2))*COS(tht_rot) &
                            +      SIN(Y%tht(i2))*SIN(tht_rot) &
                            *      COS(ABS(Y%phi(j2) -phi_rot)))
                        t = rho/info%thtc
                        nn=EXP(2D0*EXP(-1D0/t)/(t - 1D0))
                        IF(t .gt. 1) nn = 0D0

!                       Happens a lot that r = 0, causing Gij to be 1/0
!                       but with nn this tends to 0.
                        IF(periodic) THEN
!                           1 box to cover if celli present and near point in diff box
                            IF(NORM2(r) .gt. EPSILON(r)*1D1) &
                            ft2 = PGij(r, 1, info%bv, info%eye, info%xi)*(1D0 - nn)
                        ELSE
                            IF(NORM2(r) .gt. EPSILON(r)*1D1) &
                            ft2 = Gij(r, info%eye)*(1D0 - nn)
                        ENDIF
                        ft  = frot(:,i2,j2)

                        v_tmp = v_tmp + INNER3_33(ft,ft2)*dphi*wgi(i2)
                    ENDDO
                ENDDO
            ENDIF

            v_tmp = -v_tmp/(4D0*pi*(1D0 + cell%lam))
            IF(.not. PRESENT(celli)) v_tmp = Uc*2D0/(1D0 + cell%lam) + v_tmp
            v(row:row + 2) = v_tmp
        ENDDO
    ENDDO

END SUBROUTINE RHS_realCell

! -------------------------------------------------------------------------!
! Commented to improve compilation speeds...
! Does LHS for Spherical harmonic input (usefel for comparison w/ old way)
! SUBROUTINE LHS_SpHCell(cell, v, n, m)
!     CLASS(cellType), INTENT(IN), TARGET :: cell
!     COMPLEX(KIND = 8), ALLOCATABLE :: v_in(:,:,:)
!     COMPLEX(KIND = 8), INTENT(OUT), ALLOCATABLE :: v(:,:)
!     INTEGER, INTENT(IN) :: n, m
!     TYPE(sharedType), POINTER:: info
!     INTEGER :: i, j, nt, np, i2, j2, ic, row, iter, iper, jper, kper, indi, indj, tot
!     LOGICAL sing, periodic
!     REAL(KIND = 8) :: xcr(3), r(3), dphi, rn, rcur(3), minr
!     COMPLEX(KIND = 8), ALLOCATABLE :: xmnR(:,:), nmnR(:,:), v_in_mn(:,:), umnR(:,:), &
!                                  urot(:,:,:), ft(:)
!     REAL(KIND = 8), ALLOCATABLE :: xcg(:,:,:), nJt(:,:,:), &
!                                    wgi(:), ft2(:,:)
!     COMPLEX(KIND = 8) ::  v_tmp(3), u_pt(3)
!     TYPE(YType), POINTER :: Y, Yfi
!     INTEGER(KIND = 8) :: tic, toc, rate

!     info => cell%info
!     Y    => info%Y

!     ALLOCATE(&
!     xcg(3, Y%nt, Y%np), &
!     nJt(3, Y%nt, Y%np), &
!     urot(3, Y%nt, Y%np), &
!     ft(3), ft2(3,3), &
!     nmnR(3, (info%p+1)*(info%p+1)), &
!     xmnR(3, (info%p+1)*(info%p+1)), &
!     umnR(3, (info%p+1)*(info%p+1)), &
!     v_in_mn(3, (info%p+1)*(info%p+1)), &
!     wgi(Y%nt), &
!     v(3, 3*Y%nt*Y%np), &
!     v_in(3, Y%nt, Y%np))

!     v  = 0D0
!     DO tot = 1,3

!     v_in = 0D0
!     v_in(tot,:,:) = Y%nm(n+1)%v(m + n + 1, :,:)
!     v_in_mn = 0D0
!     v_in_mn(tot, n*n + 1 + n + m) = 1D0

!     ic = 0
!     iter = 0
! !   Loop over all of the integration points now and calculate in the standard way
!     DO i = 1, Y%nt
!         DO j = 1,Y%np
!             xcr = cell%x(:, i, j)
!             v_tmp = 0D0
!             sing = .false.
!             ic = ic + 1
!             row = 3*(ic - 1)  + 1

! !           Conditinal statements: integration method varies depending on problem
! !           parameters (periodic,. mult cells, etc.)
! !           Construct grid appropriate to problem
! !           For periodic, the singular integration doesn't work and we must use near-sing
!             IF(info%periodic) THEN                    
!                 sing = .true.
!                 Yfi => cell%info%Yf

! !               Need to integrate on finer grid
!                 nt = Yfi%nt + Y%nt
!                 np = Yfi%np + Y%np

! !               Deallocate integ. quants
!                 DEALLOCATE(xcg, nJt, urot, wgi)
!                 ALLOCATE( &
!                 xcg(3,  nt, np), &
!                 urot(3, nt, np), &
!                 nJt(3,  nt, np), &
!                 wgi(nt))

!                 dphi = info%Y%dphi

! !               Manage the additional prefactors stemming from the integrals
!                 wgi(1:Y%nt)  = info%Y%wg*(pi - info%thtc)/2D0
!                 wgi(Y%nt + 1:Y%nt + Yfi%nt)  = Yfi%wg*info%h*(-info%k)*COSH(info%k*info%xsf - info%k)

! !               We don't do cosine transformation to cluster points near near-singularity, mult sine back in
!                 wgi = wgi*SIN(info%Ys%th(:,1))

! !               See above note on inefficient formulation

! !               Rotate about nearest point to projected singularity
!                 xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
!                 xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
!                 xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

! !               These backward trnsforms are quite expensive
!                 xcg(1,:,:) = info%Ys%backward(xmnR(1,:), Y%nt, Y%np, info%p, cell%xrecon(1,:))
!                 xcg(2,:,:) = info%Ys%backward(xmnR(2,:), Y%nt, Y%np, info%p, cell%xrecon(2,:))
!                 xcg(3,:,:) = info%Ys%backward(xmnR(3,:), Y%nt, Y%np, info%p, cell%xrecon(3,:))

! !               Rotate the normal vector total constants
!                 nmnR(1,:) = Y%rotate(cell%nkt(1,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
!                 nmnR(2,:) = Y%rotate(cell%nkt(2,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
!                 nmnR(3,:) = Y%rotate(cell%nkt(3,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))

!                 nJt(1,:,:) = info%Ys%backward(nmnR(1,:), Y%nt, Y%np, info%p, cell%nJrecon(1,:))
!                 nJt(2,:,:) = info%Ys%backward(nmnR(2,:), Y%nt, Y%np, info%p, cell%nJrecon(2,:))
!                 nJt(3,:,:) = info%Ys%backward(nmnR(3,:), Y%nt, Y%np, info%p, cell%nJrecon(3,:))

! !               Velocities on rotated grid

! !               Rotate the normal vector total constants
!                 umnR(1,:) = rotatecmpY(Y,v_in_mn(1,:), Y%phi(j), -Y%tht(i), -Y%phi(j))
!                 umnR(2,:) = rotatecmpY(Y,v_in_mn(2,:), Y%phi(j), -Y%tht(i), -Y%phi(j))
!                 umnR(3,:) = rotatecmpY(Y,v_in_mn(3,:), Y%phi(j), -Y%tht(i), -Y%phi(j))

!                 urot(1,:,:) = backwardfastcmpY(info%Ys, umnR(1,:), Y%nt, Y%np, info%p)
!                 urot(2,:,:) = backwardfastcmpY(info%Ys, umnR(2,:), Y%nt, Y%np, info%p)
!                 urot(3,:,:) = backwardfastcmpY(info%Ys, umnR(3,:), Y%nt, Y%np, info%p)
! !           Fully singular integration on same cell: Get rotated constants
!             ELSE
!                 sing = .false.

! !               For integration (changes if need a finer grid)
!                 dphi = Y%dphi

! !               Number intg points to loop over for this case
!                 nt = Y%nt
!                 np = Y%np

!                 DEALLOCATE(urot, xcg, nJt, wgi)
!                 ALLOCATE( &
!                 xcg(3,  nt, np), &
!                 urot(3,  nt, np), &
!                 nJt(3,  nt, np), &
!                 wgi(nt))

!                 xmnR(1,:) = Y%rotate(cell%xmn(1,:), i, j, -Y%phi(j))
!                 xmnR(2,:) = Y%rotate(cell%xmn(2,:), i, j, -Y%phi(j))
!                 xmnR(3,:) = Y%rotate(cell%xmn(3,:), i, j, -Y%phi(j))

! !               Rotated integration points in unrotated frame
!                 xcg(1,:,:) = Y%backward(xmnR(1,:), info%p, cell%xrecon(1,:))
!                 xcg(2,:,:) = Y%backward(xmnR(2,:), info%p, cell%xrecon(2,:))
!                 xcg(3,:,:) = Y%backward(xmnR(3,:), info%p, cell%xrecon(3,:))

! !               Rotate the normal vector total constants
!                 nmnR(1,:) = Y%rotate(cell%nkt(1,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
!                 nmnR(2,:) = Y%rotate(cell%nkt(2,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))
!                 nmnR(3,:) = Y%rotate(cell%nkt(3,1:(info%p+1)*(info%p+1)), i, j, -Y%phi(j))

!                 nJt(1,:,:) = Y%backward(nmnR(1,:), info%p, cell%nJrecon(1,:))
!                 nJt(2,:,:) = Y%backward(nmnR(2,:), info%p, cell%nJrecon(2,:))
!                 nJt(3,:,:) = Y%backward(nmnR(3,:), info%p, cell%nJrecon(3,:))
                    
! !               Velocities on rotated grid
!                 umnR(1,:) = rotatecmpY(Y,v_in_mn(1,:), Y%phi(j), -Y%tht(i), -Y%phi(j))
!                 umnR(2,:) = rotatecmpY(Y,v_in_mn(2,:), Y%phi(j), -Y%tht(i), -Y%phi(j))
!                 umnR(3,:) = rotatecmpY(Y,v_in_mn(3,:), Y%phi(j), -Y%tht(i), -Y%phi(j))

!                 urot(1,:,:) = backwardcmpY(Y, umnR(1,:), info%p)
!                 urot(2,:,:) = backwardcmpY(Y, umnR(2,:), info%p)
!                 urot(3,:,:) = backwardcmpY(Y, umnR(3,:), info%p)

!                 wgi = Y%ws
!             ENDIF

! !           Now perform the actual integrations
!             DO j2 = 1,np
!                 DO i2 = 1,nt

! !                   Need an exception for near-singular int for if we're in mis-matched grids
!                     IF(sing .and. ((i2.gt.Y%nt .and. j2.le.Y%np) .or. (i2.le.Y%nt .and. j2.gt.Y%np))) CYCLE
!                     IF(sing .and. j2 .eq. Y%np + 1) dphi = cell%info%Yf%dphi
!                     r = xcg(:,i2,j2) - xcr

! !                   Use periodic Greens (or not)
!                     IF(info%periodic) THEN
! !                       Three cases need consideration here: singular (integ on same cell), near-singular, non-singular
! !                       Each has unique ways that these integrals need to be calc'd

!                         ! IF(NORM2(r)*info%xi.gt.3.5) THEN
!                         IF(NORM2(r)*info%xi.gt.350) THEN
!                             ft2 = -8D0*PI/info%tau*OUTER(r, nJt(:,i2,j2))*wgi(i2)
!                         ELSE
!                             ! ft2 = PTij(r, 1, info%bv, nJt(:,i2,j2), info%eye, info%xi)*wgi(i2)
!                             ft2 = PTij(r, 4, info%bv, nJt(:,i2,j2), info%eye, .3545d0, fourier = .true.)*wgi(i2)
!                         ENDIF
!                     ELSE
!                         ft2 = Tij(r, nJt(:,i2,j2))*wgi(i2)
!                     ENDIF
!                     ft2 = TRANSPOSE(ft2)

!                     ft = urot(:,i2,j2)
                    
!                     ft(1) = urot(1,i2,j2)*ft2(1,1) + urot(2,i2,j2)*ft2(2,1) + urot(3,i2,j2)*ft2(3,1)
!                     ft(2) = urot(1,i2,j2)*ft2(1,2) + urot(2,i2,j2)*ft2(2,2) + urot(3,i2,j2)*ft2(3,2)
!                     ft(3) = urot(1,i2,j2)*ft2(1,3) + urot(2,i2,j2)*ft2(2,3) + urot(3,i2,j2)*ft2(3,3)
!                     v_tmp = v_tmp + ft*dphi
!                 ENDDO
!             ENDDO
!             v_tmp = v_tmp*(1D0 - cell%lam)/(4D0*pi*(1D0 + cell%lam))
!             u_pt = v_in(:, i, j)
!             v(tot, row:row + 2) = u_pt - v_tmp
!         ENDDO
!     ENDDO
!     ENDDO

! END SUBROUTINE LHS_SpHCell

! -------------------------------------------------------------------------!
! Gets very close range interactions for a cell and periodic images of other cells
FUNCTION PeriodicCellCell(info, r) RESULT(fG)
    TYPE(sharedType), POINTER:: info
    REAL(KIND = 8) :: r(3), fG(3), rcur(3), f(3), rn, G(3,3)
    INTEGER :: i, j, k

    f = 0D0
!   Loop over all boxes
!   Check distance
!   If distance below cutoff, calculate morse and LJ
!   Multiply it with the Stokeslet of r and return
!   Then that should be the only box that has a close enough length, so we can leave fn
    DO i = -1,1
        DO j = -1,1
            DO k = -1,1
                rcur = r + i*info%bv(:,1) + j*info%bv(:,2) + k*info%bv(:,3)
                rn = NORM2(rcur)
                IF(rn .lt. 3D0*info%r0) THEN
                    f =  Morse(rcur, rn, info%D, info%r0, info%Beta)
                    f = f + LJ(rcur, rn, info%epsi, info%r0)
                    G = Gij(rcur, info%eye)
                    fG = INNER3_33(f,G)
                    RETURN
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    fG = 0D0
    RETURN

END FUNCTION PeriodicCellCell

! -------------------------------------------------------------------------!
! Reparametrize the surface if the surface grid gets distorted
!   Note that this is usually not neccessary for RBCs, because the in-plane
!   shear resistance is a built-in mechanism for preventing distortions
!   As a result, this function is somewhat underdeveloped since I only used
!   it or a single validation case
!   (e.g., we could measure grid to know when to reparam)
!   Combo of Sorgentone and Rahimian
SUBROUTINE ReparamCell(cell)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell
    TYPE(YType), POINTER :: Y, Yc
    TYPE(nmType), POINTER :: nm
    REAL(KIND = 8) :: dtau, N1, N2, nkc(3), g(3,3), delE(3)
    COMPLEX(KIND = 8), ALLOCATABLE :: xmntmp(:)
    INTEGER :: ncut, no, n, m, it, i, j, t, ih, im

    Y => cell%info%Yf
    Yc=> cell%info%Y
    ALLOCATE(xmntmp((Y%p+1)*(Y%p+1)))
    
!   Pseudo-time step (chosen somewhat arbitrarily)
    dtau = 0.25D0

!   Get the cutoff mode (see Sorgentone).
    DO no = 1, Yc%p
        N2 = 0
        it = no*no
        DO n = no, Yc%p
            DO m = -n, n
                it = it + 1
                N2 = N2 + NORM2(ABS(cell%xmn(:,it)))
            ENDDO
        ENDDO
        IF(no .eq. 1) N1 = N2
        IF(N2/N1 .lt. 0.2) THEN
            ncut = no
            EXIT
        ENDIF
    ENDDO

!   Psuedo time-step loop
    DO t = 1, 15
    DO i = 1,Y%nt
        DO j = 1,Y%np
!           Normal vector at this point
            nkc = CROSS(cell%dxt(:,i,j), cell%dxp(:,i,j))
            nkc = -nkc/NORM2(nkc)
            ! print *, cell%xf(:,i,j)
            ! print *, NKc
            ! stop

            delE = 0D0
            DO n = ncut, Yc%p
                ih = n+1
                it = n*n
                im = 0
                nm => Y%nm(ih)

!               Calc optimization gradient
                DO m = -n, n
                    im = im + 1
                    it = it + 1
                    delE = delE + REAL(cell%xmn(:,it)*nm%v(im,i,j))
                ENDDO
            ENDDO
            g = (cell%info%eye - OUTER(nkc,nkc))
            nkc = INNER3_33(delE, g)

!           Step this x-point forward
            cell%xf(:,i,j) = cell%xf(:,i,j) - dtau*nkc

        ENDDO
    ENDDO
!   Need to get new shape quantities associated with the new grid to step again,
!   starting with shape coefficients
!       A bit overkill, but, again, not worth optimizing
    xmntmp = Y%forward(cell%xf(1,:,:))
    cell%xmn(1,:) = xmntmp(1:(Yc%p+1)*(Yc%p+1))
    xmntmp = Y%forward(cell%xf(2,:,:))
    cell%xmn(2,:) = xmntmp(1:(Yc%p+1)*(Yc%p+1))
    xmntmp = Y%forward(cell%xf(3,:,:))
    cell%xmn(3,:) = xmntmp(1:(Yc%p+1)*(Yc%p+1))
    CALL cell%derivs()

    ENDDO

!   Make sure to get the coarse grid too!
    cell%x(1,:,:) = cell%info%Y%backward(cell%xmn(1,:))
    cell%x(2,:,:) = cell%info%Y%backward(cell%xmn(2,:))
    cell%x(3,:,:) = cell%info%Y%backward(cell%xmn(3,:))

    RETURN
END SUBROUTINE ReparamCell

END MODULE SHAPEMOD