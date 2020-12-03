MODULE SHAPEMOD
USE HARMMOD
IMPLICIT NONE

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! The actual cell shape object
TYPE cellType

!   Time
    INTEGER :: cts, NT
    REAL(KIND = 8) :: dt

!   Material Properties
    REAL(KIND = 8) :: mu, lam, B, C, Eb

!   Harmonics info
    INTEGER :: p, q, ftot
    TYPE(YType) :: Y, Yf
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:), xmnR(:,:)

!   Geometric info
    REAL(KIND = 8), ALLOCATABLE :: J(:,:), x(:,:,:), xf(:,:,:), k(:,:)
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
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:)

!   Velocity gradient
    REAL(KIND = 8) :: dU(3,3)

!   Cell velocity constants
    COMPLEX(KIND = 8), ALLOCATABLE :: umn(:,:)    

!   Has this object been initialized?
    LOGICAL :: init = .false.

!   Name of output file
    CHARACTER(len = 15) fileout
    
    CONTAINS
    PROCEDURE :: Write   => Writecell
    PROCEDURE :: Derivs  => Derivscell
    PROCEDURE :: Stress  => Stresscell
    PROCEDURE :: Fluid   => Fluidcell
    PROCEDURE :: Dealias => Dealiascell
END TYPE cellType

INTERFACE cellType
    PROCEDURE :: newcell
END INTERFACE cellType

CONTAINS
!=============================================================================!
!================================= ROUTIUNES =================================!
!=============================================================================!

! Constructs cell object, takes in order and alias amount
FUNCTION newcell(filein) RESULT(cell)
    CHARACTER(len = *), INTENT(IN) :: filein
    TYPE(cellType) :: cell

    INTEGER :: nt, np, ntf, npf, fali,p

!   Material properties from input
    cell%mu = READ_GRINT_DOUB(filein, 'Viscosity')
    cell%lam = READ_GRINT_DOUB(filein, 'Viscosity_ratio')
    ! cell%B = 12.4D0
    ! cell%C = 200D0
    cell%B = READ_GRINT_DOUB(filein, 'Shear_Modulus')
    cell%C = READ_GRINT_DOUB(filein, 'Dilatation_Modulus')
    cell%Eb = READ_GRINT_DOUB(filein, 'Bending_Modulus')
    cell%cts = 0
    cell%NT = READ_GRINT_INT(filein, 'Max_time_steps')
    cell%dt = READ_GRINT_DOUB(filein, 'Time_step')
    cell%fileout = READ_GRINT_CHAR(filein, 'Output')

!   Coarse and fine grids    
    p = READ_GRINT_INT(filein, 'Harmonic_order')
    cell%p = p
    fali   = READ_GRINT_INT(filein, 'Refinement_factor')
    cell%q = cell%p*fali

    cell%Y = YType(p, 1)
    cell%Yf = YType(p*fali, 4)
    print *, 'Harmonics made'

    nt  = cell%Y%nt
    np  = cell%Y%np
    ntf = cell%Yf%nt
    npf = cell%Yf%np

!   Velocity gradient (from input eventually)
    cell%dU = 0
    cell%dU(1,3) = 1

!   Allocating everything, starting with derivatives
    ALLOCATE(cell%dxt(3,ntf,npf), cell%dxp(3,ntf,npf), cell%dxp2(3,ntf,npf), &
             cell%dxt2(3,ntf,npf), cell%dxtp(3,ntf,npf), cell%dxp3(3,ntf,npf), &
             cell%dxt3(3,ntf,npf), cell%dxt2p(3,ntf,npf), cell%dxtp2(3,ntf,npf), &
             cell%dxp4(3,ntf,npf), cell%dxtp3(3,ntf,npf), cell%dxt2p2(3,ntf,npf), &
             cell%dxt3p(3,ntf,npf), cell%dxt4(3,ntf,npf), &
             cell%J(ntf,npf), cell%x(3,nt,np), cell%xf(3,ntf,npf), &
             cell%xmn(3,(p+1)*(p+1)), cell%xmnR(3,(p+1)*(p+1)), cell%umn(3,(p+1)*(p+1)))
!   Reference items
    ALLOCATE(cell%kR(ntf,npf), cell%kdR(2,ntf,npf), cell%kd2R(3,ntf,npf), &
             cell%c1R(3,ntf,npf), cell%c2R(3,ntf,npf), cell%c1tR(3,ntf,npf), &
             cell%c1pR(3,ntf,npf), cell%c2tR(3,ntf,npf), cell%c2pR(3,ntf,npf))

!   Force items
    ALLOCATE(cell%fab(3,ntf,npf), cell%ff(3,ntf,npf), cell%fc(3,nt,np), &
             cell%fmn(3,(cell%q+1)*(cell%q+1)))
    cell%ff = 0D0  

!   Get the harmonic constants to start out
    cell%xmn = RBCcoeff(cell%Y)
    
!   Initial positions
    cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
    cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
    cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))

!   Initial surfce derivatives
    CALL cell%Derivs()

!   Initial forces. Forces aren't needed and should be zero, but this serves
!   the important function of getting reference state variables
    CALL cell%Stress()

!   Write initial configuration
    CALL cell%write()    
    
    cell%init = .true.
END FUNCTION


! -------------------------------------------------------------------------!
! Writes xmn to a text file. Could be better
SUBROUTINE Writecell(cell)
    CLASS(cellType), INTENT(IN) ::cell
    IF(cell%init) THEN
        OPEN (UNIT = 88, STATUS = "old", POSITION = "append", FILE = "xmn.txt")
    ELSE
        OPEN (UNIT = 88, FILE = "xmn.txt")
    ENDIF
    WRITE(88,*) REAL(cell%xmn)
    WRITE(88,*) AIMAG(cell%xmn)
    WRITE(88,*) cell%cts
    CLOSE(88)

END SUBROUTINE Writecell


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
        DO m = -n,n
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

!   De-alias all derivs
!   First order
!     CALL cell%dealias(cell%dxt(1,:,:))
!     CALL cell%dealias(cell%dxt(2,:,:))
!     CALL cell%dealias(cell%dxt(3,:,:))
    
!     CALL cell%dealias(cell%dxp(1,:,:))
!     CALL cell%dealias(cell%dxp(2,:,:))
!     CALL cell%dealias(cell%dxp(3,:,:))
    
! !   Second order
!     CALL cell%dealias(cell%dxt2(1,:,:))
!     CALL cell%dealias(cell%dxt2(2,:,:))
!     CALL cell%dealias(cell%dxt2(3,:,:))
    
!     CALL cell%dealias(cell%dxtp(1,:,:))
!     CALL cell%dealias(cell%dxtp(2,:,:))
!     CALL cell%dealias(cell%dxtp(3,:,:))

!     CALL cell%dealias(cell%dxp2(1,:,:))
!     CALL cell%dealias(cell%dxp2(2,:,:))
!     CALL cell%dealias(cell%dxp2(3,:,:))

! !   Third order
!     CALL cell%dealias(cell%dxt3(1,:,:))
!     CALL cell%dealias(cell%dxt3(2,:,:))
!     CALL cell%dealias(cell%dxt3(3,:,:))
    
!     CALL cell%dealias(cell%dxt2p(1,:,:))
!     CALL cell%dealias(cell%dxt2p(2,:,:))
!     CALL cell%dealias(cell%dxt2p(3,:,:))

!     CALL cell%dealias(cell%dxtp2(1,:,:))
!     CALL cell%dealias(cell%dxtp2(2,:,:))
!     CALL cell%dealias(cell%dxtp2(3,:,:))

!     CALL cell%dealias(cell%dxp3(1,:,:))
!     CALL cell%dealias(cell%dxp3(2,:,:))
!     CALL cell%dealias(cell%dxp3(3,:,:))

! !   Fourth order
!     CALL cell%dealias(cell%dxt4(1,:,:))
!     CALL cell%dealias(cell%dxt4(2,:,:))
!     CALL cell%dealias(cell%dxt4(3,:,:))
    
!     CALL cell%dealias(cell%dxt3p(1,:,:))
!     CALL cell%dealias(cell%dxt3p(2,:,:))
!     CALL cell%dealias(cell%dxt3p(3,:,:))

!     CALL cell%dealias(cell%dxt2p2(1,:,:))
!     CALL cell%dealias(cell%dxt2p2(2,:,:))
!     CALL cell%dealias(cell%dxt2p2(3,:,:))

!     CALL cell%dealias(cell%dxtp3(1,:,:))
!     CALL cell%dealias(cell%dxtp3(2,:,:))
!     CALL cell%dealias(cell%dxtp3(3,:,:))

!     CALL cell%dealias(cell%dxp4(1,:,:))
!     CALL cell%dealias(cell%dxp4(2,:,:))
!     CALL cell%dealias(cell%dxp4(3,:,:))
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
    REAL(KIND = 8) :: mab(2,2), mabt(2,2), mabp(2,2), dmab2(2,2), &
                      q1, q2, dq(2), Fd(3,3), Prj(3,3), V2(3,3), &
                      lams(3), es(2), ev(3,3), ev1(3), ev2(3), I1, I2, &
                      tau(3,3), tauab(2,2), B, C, Eb, wrk(102), tV2(3,3)
    REAL(KIND = 8) :: Fdt(3,3), Fdp(3,3), V2t(3,3), V2p(3,3), lam1dV2(3,3), &
                      lam2dV2(3,3), es1t, es1p, es2t, es2p, I1t, I1p, I2t, I2p, &
                      Prjt(3,3), Prjp(3,3), taut(3,3), taup(3,3), dtauab(2,2), &
                      bv(2,2), bm(2,2), cvt, cvp, cq, normev1, normev2
                      


                      COMPLEX(KIND = 8), allocatable :: tmp(:), tmp2(:,:)
    
    B = cell%B
    C = cell%C
    Eb = cell%Eb
    
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

!           Now finally some mechanics, starting with bending moments/derivs
            mab = (k - cell%kR(i,j))*gn
            mabt = (kt - cell%kdR(1,i,j))*gn + (k - cell%kR(i,j))*dgt
            mabp = (kp - cell%kdR(2,i,j))*gn + (k - cell%kR(i,j))*dgp
            
            dmab2(1,1) = (kt2 - cell%kd2R(1,i,j))*gn(1,1) + 2D0*(kt - cell%kdR(1,i,j))*dgt(1,1) &
                       + (k - cell%kR(i,j))*dgt2(1,1)
            dmab2(1,2) = (ktp - cell%kd2R(3,i,j))*gn(1,2) + (kt - cell%kdR(1,i,j))*dgp(1,2) &
                       + (kp  - cell%kdR(2,i,j))*dgt(1,2) + (k  - cell%kR(i,j))*dgtp(1,2)
            dmab2(2,1) = (ktp - cell%kd2R(3,i,j))*gn(2,1) + (kt - cell%kdR(1,i,j))*dgp(2,1) &
                       + (kp  - cell%kdR(2,i,j))*dgt(2,1) + (k  - cell%kR(i,j))*dgtp(2,1)
            dmab2(2,2) = (kp2 - cell%kd2R(2,i,j))*gn(2,2) + 2D0*(kp - cell%kdR(2,i,j))*dgp(2,2) &
                       + (k - cell%kR(i,j))*dgp2(2,2)

!           Transverse shears and partials
            q1 = Eb*(mabt(1,1) + mabp(2,1) + mab(1,1)*(2D0*c111 + c221) &
               + mab(2,1)*(2D0*c112 + c222) + mab(1,2)*c112 + mab(2,2)*c122)
            q2 = Eb*(mabt(1,2) + mabp(2,2) + mab(1,2)*(c111 + 2D0*c221) &
               + mab(2,2)*(c112 + 2D0*c222) + mab(1,1)*c211 + mab(2,1)*c221)
            
!           Partials of q
            dq(1) = Eb*(dmab2(1,1) + dmab2(2,1) &
                  + mabt(1,1)*(2D0*c111 + c221) + mab(1,1)*(2D0*c111t + c221t) &
                  + mabt(2,1)*(2D0*c112 + c222) + mab(2,1)*(2D0*c112t + c222t) &
                  + mabt(1,2)*c112 + mab(1,2)*c112t + mabt(2,2)*c122 + mab(2,2)*c122t)
            
            dq(2) = Eb*(dmab2(1,2) + dmab2(2,2) &
                  + mabp(1,2)*(c111 + 2D0*c221) + mab(1,2)*(c111p + 2D0*c221p) &
                  + mabp(2,2)*(c112 + 2D0*c222) + mab(2,2)*(c112p + 2D0*c222p) &
                  + mabp(1,1)*c211 + mab(1,1)*c211p + mabp(2,1)*c221 + mab(2,1)*c221p)

!           In plane tension and derivatives
!           Deformation gradient tensor                  
            Fd = outer(cell%dxt(:,i,j), cell%c1R(:,i,j)) + outer(cell%dxp(:,i,j), cell%c2R(:,i,j))

!           Surface projection operator
            Prj = 0D0
            FORALL(q = 1:3) Prj(q,q) = 1D0
            Prj = Prj - outer(nk,nk)

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

! !           Old way of doing it
!             CALL EIG3(V2, lams, ev)
!             es(1) = sqrt(lams(1))
!             es(2) = sqrt(lams(2))

! !           Eigenvectors I actually need
!             ev1 = ev(:,1)/norm2(ev(:,1))
!             ev2 = ev(:,2)/norm2(ev(:,3))

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
            dtauab(1,1) = DOT(INNER3_33(c1,taut), c1)
            dtauab(1,2) = DOT(INNER3_33(c1,taut), c2)
            dtauab(2,1) = DOT(INNER3_33(c2,taup), c1)
            dtauab(2,2) = DOT(INNER3_33(c2,taup), c2)

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

!           Covariant divergence of q 
            cq = dq(1) + dq(2) + c111*q1 + c112*q2 + c221*q1 + c222*q2

!           Forces in surface coordinates
            cell%fab(1, i, j) = bm(1,1)*q1 + bm(2,1)*q2 - cvt
            cell%fab(2, i, j) = bm(1,2)*q1 + bm(2,2)*q2 - cvp
            cell%fab(3, i, j) = -cq - tauab(1,1)*bv(1,1) - tauab(1,2)*bv(1,2) &
                                    - tauab(2,1)*bv(2,1) - tauab(2,2)*bv(2,2)

!           Forces in Cartesian, on fine grid
            cell%ff(:, i, j) = cell%fab(1, i, j)*cell%dxt(:,i,j) &
                              + cell%fab(2, i, j)*cell%dxp(:,i,j) &
                              + cell%fab(3, i, j)*(-nk)
        ENDDO inner
    ENDDO

!!!! Either my thetder algorithm or harmonics algorithm appears to be slightly unstable.
    ! At higher order harmonics, they start to differ by about 1e-14 from matlab, which
    ! gets exacerbated by the eigenvalue algorithm. Perhaps the problem is in the legendre
    ! algorithm?

!   Now we need to filter for anti-aliasing. Do the transform of the force into spectral
!   space, cut the highest modes, and transform back into physical space
    IF(cell%init) THEN
        cell%fmn(1,:) = cell%Yf%forward(cell%ff(1,:,:), cell%q)
        cell%fmn(2,:) = cell%Yf%forward(cell%ff(2,:,:), cell%q)
        cell%fmn(3,:) = cell%Yf%forward(cell%ff(3,:,:), cell%q)

        cell%fc(1,:,:) = cell%Y%backward(cell%fmn(1,:), cell%p)
        cell%fc(2,:,:) = cell%Y%backward(cell%fmn(2,:), cell%p)
        cell%fc(3,:,:) = cell%Y%backward(cell%fmn(3,:), cell%p)

        tmp = cell%fmn(1,:)
        ! tmp2 = cell%ff(1,:,:)
        ! print *, tmp !!!!!!! precision dips: from cvt=>tauab,dtauab, I2t,es1t, FROM C=100!
    ENDIF
END SUBROUTINE Stresscell

! -------------------------------------------------------------------------!
! Knowing the force jump get the velocity on the surface of the cell via
! the fluid problem
SUBROUTINE Fluidcell(cell)
    CLASS(cellType), INTENT(INOUT), TARGET :: cell

    INTEGER :: ip, ic, i, j, i2, j2, n, m, ih, it, im, row, col, im2, n2, m2
    COMPLEX(KIND = 8) :: A(1,1), &!A(3*cell%Y%nt*cell%Y%np, 3*(cell%Y%p + 1)*(cell%Y%p + 1)), &
                         b(3*cell%Y%nt*cell%Y%np), At(3,3), bt(3), &
                         A2(1,1), &!A2(3*(cell%Y%p + 1)*(cell%Y%p + 1), 3*(cell%Y%p + 1)*(cell%Y%p + 1)), &
                         b2(3*(cell%Y%p + 1)*(cell%Y%p + 1)), td1, &
                         vcur, v(3,3)
    REAL(KIND = 8) :: dxtg(3), dxpg(3), Uc(3), t1(3,3), t2(3,3), t3(3,3), Tx(3,3), &
                      xcr(3), gp(3), Utmp(3,3), &
                      thet(cell%Y%nt, cell%Y%np), phit(cell%Y%nt, cell%Y%np), nkg(3), &
                      Jg(cell%Y%nt, cell%Y%np), frot(3, cell%Y%nt, cell%Y%np), &
                      r(3), vT(3,3,cell%Y%nt, cell%Y%np), vG(3,3,cell%Y%nt, cell%Y%np), &
                      eye(3,3) , xcg(3, cell%Y%nt, cell%Y%np), Jgf(cell%Yf%nt, cell%Yf%np), &
                      Vtt(3,3), Vtg(3,3)
    COMPLEX(KIND = 8), POINTER :: vcurn(:,:), vcurt(:,:)
    TYPE(YType), POINTER :: Y, Yf
    TYPE(Ytype), TARGET :: Yt
    TYPE(nmType), POINTER :: nm, nmt
    
            COMPLEX(KIND = 8), allocatable :: tmp(:)

    Y => cell%Y
    Yf=> cell%Yf
    eye = 0D0
    FORALL(j = 1:3) eye(j,j) = 1D0
!   We need to do two integral loops: the first calculates the integral 
!   that is the convolution of the Stokeslets and stresslets. We need those
!   values at the integration points for the next step. Next we do the
!   Galerkin integral over a sphere, using the inner integrals calculated
!   in the previous step

    ip = 0
    ic = 0
    A = 0D0
    b = 0D0
!   First loop: inner integrals at GPs
    DO i = 1, Y%nt
        DO j = 1,Y%np
!           Total Galerkin mode count
            ic = ic + 1

!           Velocity at integration point
            Utmp = TRANSPOSE(cell%dU)
            Uc = INNER3_33(cell%x(:,i,j), Utmp)

!           To integrate these nasty integrals, we need to rotate the 
!           spherical harmonic coefficients so the integration point is at
!           the north pole. Process as follows:

!           Construct rotation matrix
            t1(1,:) = (/COS( Y%phi(j)), -SIN( Y%phi(j)), 0D0/)
            t2(1,:) = (/COS(-Y%tht(i)), 0D0, SIN(-Y%tht(i))/)
            t3(1,:) = (/COS(-Y%phi(j)), -SIN(-Y%phi(j)), 0D0/)

            t1(2,:) = (/SIN( Y%phi(j)), COS( Y%phi(j)), 0D0/)
            t2(2,:) = (/0D0, 1D0, 0D0/)
            t3(2,:) = (/SIN(-Y%phi(j)), COS(-Y%phi(j)), 0D0/)

            t1(3,:) = (/0D0, 0D0, 1D0/)
            t2(3,:) = (/-SIN(-Y%tht(i)), 0D0, COS(-Y%tht(i))/)
            t3(3,:) = (/0D0, 0D0, 1D0/)

            Tx = MATMUL(MATMUL(t1,t2),t3)

!           Get the rotated constants
            cell%xmnR(1,:) = Y%rotate(cell%xmn(1,:), Y%phi(j), -Y%tht(i), -Y%phi(j))
            cell%xmnR(2,:) = Y%rotate(cell%xmn(2,:), Y%phi(j), -Y%tht(i), -Y%phi(j))
            cell%xmnR(3,:) = Y%rotate(cell%xmn(3,:), Y%phi(j), -Y%tht(i), -Y%phi(j))

!           Rotated integration points in unrotated frame
            xcg(1,:,:) = Y%backward(cell%xmnR(1,:))
            xcg(2,:,:) = Y%backward(cell%xmnR(2,:))
            xcg(3,:,:) = Y%backward(cell%xmnR(3,:))

!           Location of north pole in unrotated frame (just current integration point)
            xcr = cell%x(:,i,j)
            
!           Area element, needed for integration, in rotated reference frame
!           Additionally, we want rotated points in nonrotated reference frame
            DO i2 = 1,Y%nt
                DO j2 = 1,Y%np
!                   Gauss point rotated in parameter space
                    gp = (/SIN(Y%tht(i2))*COS(Y%phi(j2)), &
                           SIN(Y%tht(i2))*SIN(Y%phi(j2)), &
                           COS(Y%tht(i2))/)

!                   Rotate this Gauss point to nonrotated parameter space !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check this is right
                    gp = INNER3_33(gp, Tx)
!                   Sometimes precision can be an issue...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Don't need to do this everytime (calcs same spot always)
                    IF(gp(3).gt.1)  gp(3) =  1D0
                    IF(gp(3).lt.-1) gp(3) = -1D0
                    phit(i2, j2) = ATAN2(gp(2), gp(1))
                    thet(i2, j2) = ACOS(gp(3))

!                   We need the basis vectors in the rotated parameter space !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Perhaps theres a way to speed this up
!                   to get the area element
                    it = 0
                    ih = 0
                    dxtg = 0
                    dxpg = 0
                    DO n = 0,Y%p
                        ih = ih + 1
                        nm => Y%nm(ih)
                        im = 0
                        DO m = -n,n
                            im = im + 1
                            it = it + 1
                            
                            td1 = nm%dY(1,im,i2,j2)
                            dxtg(1) = dxtg(1) + REAL(cell%xmnR(1,it)*td1)
                            dxtg(2) = dxtg(2) + REAL(cell%xmnR(2,it)*td1)
                            dxtg(3) = dxtg(3) + REAL(cell%xmnR(3,it)*td1)

                            vcur = nm%v(im,i2,j2)
                            dxpg(1) = dxpg(1) + REAL(cell%xmnR(1,it)*ii*m*vcur)
                            dxpg(2) = dxpg(2) + REAL(cell%xmnR(2,it)*ii*m*vcur)
                            dxpg(3) = dxpg(3) + REAL(cell%xmnR(3,it)*ii*m*vcur)
                        ENDDO
                    ENDDO
!                   Inward normal                    
                    nkg = CROSS(dxtg, dxpg)
                    nkg = -nkg/(sqrt(nkg(1)*nkg(1) + nkg(2)*nkg(2) + nkg(3)*nkg(3)))

!                   Jacobian via fundamental forms
                    Jg(i2,j2) = sqrt(DOT(dxtg,dxtg)*DOT(dxpg,dxpg) &
                              -      DOT(dxtg,dxpg)*DOT(dxtg,dxpg))

!                   Calculate kernels now, since we're already looping over these points
                    r = xcr - xcg(:,i2,j2)
                    vG(:,:,i2,j2) = Gij(r, eye)
                    vT(:,:,i2,j2) = Tij(r, nkg)
                ENDDO
            ENDDO
!           Harmonics of integration points of rotated frame in nonrotated frame
            Yt = Ytype(thet, phit)

!           Forces on rotated grid
            frot(1,:,:) = Yt%backward(cell%fmn(1,:))
            frot(2,:,:) = Yt%backward(cell%fmn(2,:))
            frot(3,:,:) = Yt%backward(cell%fmn(3,:))

!           Bookkeeping
            row = 3*(ic - 1)  + 1
            bt = 0
            ih = 0

!           Loop over harmonics
            DO n = 0, Y%p
                nm  => Y%nm(n+1)
                nmt => Yt%nm(n+1)
                im = 0
                DO m = -n,n
                    ih = ih + 1
                    im = im + 1
                    col = 3*(ih-1) + 1
                    vcurn => nm%v(im,:,:)
                    vcurt => nmt%v(im,:,:)
                    At = 0
                    
!                   Here's where the actul integral is performed,
!                   centered on point i, j
                    DO i2 = 1,Y%nt
                        DO j2 = 1,Y%np
!                           Add in integral parts
                            At = At + vT(:,:,i2,j2)*vcurt(i2,j2)*Jg(i2,j2)&
                               * cell%Y%ws(i2)*cell%Y%dphi
!                           RHS part, only need to do once
                            IF(n .eq. 0) THEN
                                bt = bt + INNER3_33(frot(:,i2,j2),vG(:,:,i2,j2)) &
                                   * Jg(i2,j2)*cell%Y%ws(i2)
                            ENDIF
                        ENDDO
                    ENDDO
!                   Add in the rest of the LHS that isn't an integral
                    At = At*(1D0-cell%lam)/(1D0+cell%lam) - vcurn(i,j)*4D0*pi*eye

!                   LHS at integration point/harmonic combo, put in big matrix
                    !A(row:row+2, col:col+2) = A(row:row+2, col:col+2) + At
                ENDDO
            ENDDO
!           Put RHS into big vector
            b(row:row+2) = b(row:row+2) + bt*cell%Y%dphi/cell%mu/(1+cell%lam) &
                         - Uc*8D0*pi/(1+cell%lam)
        ENDDO
    ENDDO

!   Second integral loop, Galerkin
    A2 = 0D0
    b2 = 0D0
    it = 0

!   Loop over outer product harmonics
    DO n = 0,Y%p
        nm => cell%Y%nm(n+1)
        im = 0
        DO m = -n,n
            im = im + 1
            it = it + 1
            row = 3*it - 2
            bt(:) = 0D0
            vcurn => nm%v(im,:,:)
            im2 = 0
            
!           Loop over inner harmonics (constant value is a column)
            DO n2 = 0,Y%p
                DO m2 = -n2,n2
                    im2 = im2 + 1
                    col = 3*im2 - 2
                    At = 0D0
                    ic = 0
!                   Loop over integration points to calc integral
                    DO i =1,Y%nt
                        DO j = 1,Y%np
                            ic = ic+1
!                           Value of inner integral at IP
                            !v = A(3*ic-2:3*ic, 3*im2-2:3*im2)
                            At = At + v*CONJG(vcurn(i,j))*cell%Y%wg(i)*cell%Y%dphi

!                           Intg. b (essentially forward transform of RHS!)
                            IF(n2 .eq. 0) THEN
                                bt = bt + b(3*ic-2:3*ic)*CONJG(vcurn(i,j)) &
                                   *cell%Y%wg(i)*cell%Y%dphi
                            ENDIF
                        ENDDO
                    ENDDO
                    !A2(row:row+2,col:col+2) = At
                ENDDO
            ENDDO
            b2(row:row+2) = bt
        ENDDO
    ENDDO

!   The highest order coefficients are calculated with large error,
!   so we need to throw them out
!   We could just not even calculate these at all to save some time!!!!!!!!!!!!!!!!!!!!!!!!!!
    b2(3*Y%p*Y%p+1:3*(Y%p+1)*(Y%p+1)) = 0D0

    cell%umn(1,:) = b2((/(i, i=1,3*(Y%p+1)*(Y%p+1)-2, 3)/))/(-4D0*pi)
    cell%umn(2,:) = b2((/(i, i=2,3*(Y%p+1)*(Y%p+1)-1, 3)/))/(-4D0*pi)
    cell%umn(3,:) = b2((/(i, i=3,3*(Y%p+1)*(Y%p+1)  , 3)/))/(-4D0*pi)

    ! tmp = cell%umn(3,:)
    ! print *, tmp

!   Velocity doesn't match after first time step. Surface derivs are ok though. Could be in SpT
END SUBROUTINE Fluidcell


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

END MODULE SHAPEMOD