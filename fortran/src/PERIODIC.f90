MODULE PERIODICMOD
USE SHAREDMOD
IMPLICIT NONE
    
!==============================================================================!
!               The purpose of this module is to the calculations              !
!                specifically related to the fast Ewald sum part               !
!                                                                              !
!==============================================================================!

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!
INTERFACE HtInt
PROCEDURE HtIntG, HtIntT
END INTERFACE HtInt

INTERFACE HtCalc
    PROCEDURE HtG, HtT
END INTERFACE HtCalc

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

! -------------------------------------------------------------------------!
! Some precalculations that only need to be done once per timestep to run HtCell
! Mostly taking cell forces and locations and constructing a single vector
! Only need to do this once for x0, so it's optional.
! Note that this will append fin to f
SUBROUTINE HtPreCalc(fin, f, Ja, Y, x0in, x0)
    COMPLEX(KIND = 8), INTENT(IN) :: fin(:,:)  !!!!!!! Pointers???
    COMPLEX(KIND = 8), ALLOCATABLE, INTENT(INOUT):: f(:)    !!!!!!! Pointers???
    REAL(KIND = 8), INTENT(IN) :: Ja(:,:)
    TYPE(YType), INTENT(IN) :: Y
    REAL(KIND = 8), INTENT(IN), OPTIONAL :: x0in(:,:,:)
    REAL(KIND = 8), INTENT(INOUT), ALLOCATABLE, OPTIONAL:: x0(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fapp(:), ftmp(:)
    REAL(KIND = 8), ALLOCATABLE :: xapp(:,:), xtmp(:,:)
    INTEGER :: shp(2), nt, np, i, j, it, tot

    shp = SHAPE(fin)
    nt = shp(1)
    np = shp(2)
    tot = SIZE(f)
    ALLOCATE(fapp(nt*np), xapp(3, nt*np), ftmp(nt*np + tot), xtmp(3, nt*np + tot))
    it = 0

!   Replace with "effective" point force strength
    DO i = 1,nt
        DO j = 1,np
            it = it + 1
            fapp(it) = fin(i,j)*Y%dphi*Y%wg(i)*Ja(i,j)
        ENDDO
    ENDDO

!   Now some shuffling/appending
    IF(ALLOCATED(f)) ftmp(1:tot) = f
    ftmp(tot + 1: tot + nt*np) = fapp
    IF(ALLOCATED(f)) DEALLOCATE(f)
    ALLOCATE(f(nt*np + tot))
    f = ftmp

!   Append source locations (only needs to be done once per timestep)
!   Same process as above, except it's purely just reshuffling indices
    IF(PRESENT(x0in)) THEN
        it = 0 
        DO i = 1,nt
            DO j = 1,np
                it = it + 1
                xapp(:,it) = x0in(:,i,j)
            ENDDO
        ENDDO
        IF(ALLOCATED(x0)) xtmp(:, 1:tot) = x0
        xtmp(:, tot + 1: tot + nt*np) = xapp
        IF(ALLOCATED(x0)) DEALLOCATE(x0)
        ALLOCATE(x0(3, nt*np + tot))
        x0 = xtmp
    ENDIF

END SUBROUTINE HtPreCalc

! -------------------------------------------------------------------------!
! For a given set of cells with point forces, construct the H_tilde grid
! for the single layer
! See: Spectrally accurate fast summation for periodic Stokes potentials
!       D. Lindbo, AK Tornberg
! Can be used to interpolate back to grid.
SUBROUTINE HtG(Ht, info, x0, f)
    REAL(KIND = 8), INTENT(OUT) :: Ht(:,:,:,:)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    COMPLEX(KIND = 8), INTENT(IN) :: f(:,:)    !!!!!!! Pointers???
    COMPLEX(KIND = 8), ALLOCATABLE ::  H(:,:,:,:), Hh(:,:,:,:), Hht(:,:,:,:), Htmp(:,:,:)
    REAL(KIND = 8) :: bvi_gr(3,3), bv_gr(3,3), rr(3), r2, k3(3), &
                      kn, B(3,3), xcur(3)
    REAL(KIND = 8), ALLOCATABLE :: wrk(:)
    INTEGER :: i, j, k, pts, curijk(3), gp, inds(3), &
               iper, jper, kper, iw, jw, kw, nt, lwrk

    gp = info%gp

    lwrk = gp*2*3
    ALLOCATE(wrk(lwrk))

!   Total number of points
    nt = size(x0)/3

    ALLOCATE(H(3,gp,gp,gp), Hh(3,gp,gp,gp), Hht(3,gp,gp,gp), Htmp(gp,gp,gp))
    H = 0D0

!   Subgrid basis vectors and inverse
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)

!   Now we can start the actual construction of the Ht matrix. Start with the H matrix,
!   which is just the point forces smeared onto grid points. Go to each force point and
!   smear to supporting points
    DO i = 1, nt
!       At a point, find grid index & loop over supporting grid points
        xcur = x0(:,i)
        curijk(1) = FLOOR(bvi_gr(1,1)*xcur(1) + bvi_gr(1,2)*xcur(2) + bvi_gr(1,3)*xcur(3))
        curijk(2) = FLOOR(bvi_gr(2,1)*xcur(1) + bvi_gr(2,2)*xcur(2) + bvi_gr(2,3)*xcur(3))
        curijk(3) = FLOOR(bvi_gr(3,1)*xcur(1) + bvi_gr(3,2)*xcur(2) + bvi_gr(3,3)*xcur(3))
        DO pts = 1,info%suppPoints
            inds = curijk + info%suppmat(:, pts)

!           Manage periodicity
            IF(inds(1) .lt. 1) THEN
                iper = inds(1) + gp
            ELSEIF(inds(1) .gt. gp) THEN
                iper = inds(1) - gp
            ELSE
                iper = inds(1)
            ENDIF
            IF(inds(2) .lt. 1) THEN
                jper = inds(2) + gp
            ELSEIF(inds(2) .gt. gp) THEN
                jper = inds(2) - gp
            ELSE
                jper = inds(2)
            ENDIF
            IF(inds(3) .lt. 1) THEN
                kper = inds(3) + gp
            ELSEIF(inds(3) .gt. gp) THEN
                kper = inds(3) - gp
            ELSE
                kper = inds(3)
            ENDIF

!           Add smeared point forces to grid point
            rr(1) = xcur(1) - (bv_gr(1,1)*inds(1) + bv_gr(1,2)*inds(2) + bv_gr(1,3)*inds(3)) 
            rr(2) = xcur(2) - (bv_gr(2,1)*inds(1) + bv_gr(2,2)*inds(2) + bv_gr(2,3)*inds(3)) 
            rr(3) = xcur(3) - (bv_gr(3,1)*inds(1) + bv_gr(3,2)*inds(2) + bv_gr(3,3)*inds(3)) 

            r2 = rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3)

            H(:,iper,jper,kper) = H(:,iper,jper,kper) + f(:,i)*info%par*EXP(info%parexp*r2)
        ENDDO
    ENDDO

!   FFT
    !!!! Make pointers!!!!!!!
    Htmp = H(1,:,:,:)
    Htmp = FFT3(Htmp, info%WSAVE)
    CALL FFTSHIFT(Htmp)
    Hh(1,:,:,:) = Htmp

    Htmp = H(2,:,:,:)
    Htmp = FFT3(Htmp, info%WSAVE)
    CALL FFTSHIFT(Htmp)
    Hh(2,:,:,:) = Htmp

    Htmp = H(3,:,:,:)
    Htmp = FFT3(Htmp, info%WSAVE)
    CALL FFTSHIFT(Htmp)
    Hh(3,:,:,:) = Htmp

!   Now perform the convolution with the kernel (which is a truncated sum in spectral space)
!   (gp needs to be even)
!   Also assumes the Hh matrix has negative and positive parts
    DO iw = -gp/2, gp/2 - 1
        DO jw = -gp/2, gp/2 - 1
            DO kw = -gp/2, gp/2 - 1
                i = iw + gp/2 + 1
                j = jw + gp/2 + 1
                k = kw + gp/2 + 1

                k3 = info%kv(:,1)*iw + info%kv(:,2)*jw + info%kv(:,3)*kw
                kn = NORM2(k3)

!               Truncate
                IF(kn .eq. 0 .or. kn .gt. 40) CYCLE

!               Amplification factor, gets multiplied onto spectral point forces
                B = BG(k3, info%xi, kn, info%eta, info%eye)
                Hht(1, i, j, k) = B(1,1)*Hh(1,i,j,k) &
                                + B(1,2)*Hh(2,i,j,k) &
                                + B(1,3)*Hh(3,i,j,k)
                Hht(2, i, j, k) = B(2,1)*Hh(1,i,j,k) &
                                + B(2,2)*Hh(2,i,j,k) &
                                + B(2,3)*Hh(3,i,j,k)
                Hht(3, i, j, k) = B(3,1)*Hh(1,i,j,k) &
                                + B(3,2)*Hh(2,i,j,k) &
                                + B(3,3)*Hh(3,i,j,k)

            ENDDO
        ENDDO
    ENDDO
    Ht = 0D0
    ! Double layer... DON'T FORGET TO MULTIPLY ii!!!!!! TAU????

!!!!!!!!!!!! Note: FFTshift only for even number grid. Need to make iFFTshift for odds
    Htmp = Hht(1,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Ht(1,:,:,:) = REAL(Htmp)

    Htmp = Hht(2,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Ht(2,:,:,:) = REAL(Htmp)

    Htmp = Hht(3,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Ht(3,:,:,:) = REAL(Htmp)

!   Return Ht into real space, then we can integrate/interpolate back to the point forces in the
!   routine (although maybe this should be optional prob, and we canm just add to mat at the end???)
!   Could hypothetically just do all this and add it to the big mat????
END SUBROUTINE HtG

! -------------------------------------------------------------------------!
! For a given set of cells with point forces, construct the H_tilde grid
! for the single layer
! See: Spectrally accurate fast summation for periodic Stokes potentials
!       D. Lindbo, AK Tornberg
! Can be used to interpolate back to grid.
! Note that in most cases, f will be the normal vector x a sphHarm
SUBROUTINE HtT(Ht, info, x0, f)
    REAL(KIND = 8), INTENT(OUT) :: Ht(:,:,:,:,:)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    COMPLEX(KIND = 8), INTENT(IN) :: f(:,:)    !!!!!!! Pointers???
    COMPLEX(KIND = 8), ALLOCATABLE ::  H(:,:,:,:), Hh(:,:,:,:), Hht(:,:,:,:,:), Htmp(:,:,:)
    COMPLEX(KIND = 8) :: B(3,3,3)
    REAL(KIND = 8) :: bvi_gr(3,3), bv_gr(3,3), rr(3), r2, k3(3), &
                      kn, xcur(3)
    REAL(KIND = 8), ALLOCATABLE :: wrk(:)
    INTEGER :: i, j, k, pts, curijk(3), gp, inds(3), &
               iper, jper, kper, iw, jw, kw, nt, lwrk
    gp = info%gp

    lwrk = gp*2*3
    ALLOCATE(wrk(lwrk))

!   Total number of points
    nt = size(x0)/3

    ALLOCATE(H(3,gp,gp,gp), Hh(3,gp,gp,gp), Hht(3,3,gp,gp,gp), Htmp(gp,gp,gp))
    H = 0D0

!   Subgrid basis vectors and inverse
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)

!   Now we can start the actual construction of the Ht matrix. Start with the H matrix,
!   which is just the point forces smeared onto grid points. Go to each force point and
!   smear to supporting points
    DO i = 1, nt
!       At a point, find grid index & loop over supporting grid points
        xcur = x0(:,i)
        curijk(1) = FLOOR(bvi_gr(1,1)*xcur(1) + bvi_gr(1,2)*xcur(2) + bvi_gr(1,3)*xcur(3))
        curijk(2) = FLOOR(bvi_gr(2,1)*xcur(1) + bvi_gr(2,2)*xcur(2) + bvi_gr(2,3)*xcur(3))
        curijk(3) = FLOOR(bvi_gr(3,1)*xcur(1) + bvi_gr(3,2)*xcur(2) + bvi_gr(3,3)*xcur(3))
        DO pts = 1,info%suppPoints
            inds = curijk + info%suppmat(:, pts)

!           Manage periodicity
            IF(inds(1) .lt. 1) THEN
                iper = inds(1) + gp
            ELSEIF(inds(1) .gt. gp) THEN
                iper = inds(1) - gp
            ELSE
                iper = inds(1)
            ENDIF
            IF(inds(2) .lt. 1) THEN
                jper = inds(2) + gp
            ELSEIF(inds(2) .gt. gp) THEN
                jper = inds(2) - gp
            ELSE
                jper = inds(2)
            ENDIF
            IF(inds(3) .lt. 1) THEN
                kper = inds(3) + gp
            ELSEIF(inds(3) .gt. gp) THEN
                kper = inds(3) - gp
            ELSE
                kper = inds(3)
            ENDIF

!           Add smeared point forces to grid point
            rr(1) = xcur(1) - (bv_gr(1,1)*inds(1) + bv_gr(1,2)*inds(2) + bv_gr(1,3)*inds(3)) 
            rr(2) = xcur(2) - (bv_gr(2,1)*inds(1) + bv_gr(2,2)*inds(2) + bv_gr(2,3)*inds(3)) 
            rr(3) = xcur(3) - (bv_gr(3,1)*inds(1) + bv_gr(3,2)*inds(2) + bv_gr(3,3)*inds(3)) 

            r2 = rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3)

            H(:,iper,jper,kper) = H(:,iper,jper,kper) + f(:,i)*info%par*EXP(info%parexp*r2)
        ENDDO
    ENDDO

!   FFT
    !!!! Make pointers!!!!!!!
    Htmp = H(1,:,:,:)
    Htmp = FFT3(Htmp, info%WSAVE)
    CALL FFTSHIFT(Htmp)
    Hh(1,:,:,:) = Htmp

    Htmp = H(2,:,:,:)
    Htmp = FFT3(Htmp, info%WSAVE)
    CALL FFTSHIFT(Htmp)
    Hh(2,:,:,:) = Htmp

    Htmp = H(3,:,:,:)
    Htmp = FFT3(Htmp, info%WSAVE)
    CALL FFTSHIFT(Htmp)
    Hh(3,:,:,:) = Htmp

!   Now perform the convolution with the kernel (which is a truncated sum in spectral space)
!   (gp needs to be even)
!   Also assumes the Hh matrix has negative and positive parts
    DO iw = -gp/2, gp/2 - 1
        DO jw = -gp/2, gp/2 - 1
            DO kw = -gp/2, gp/2 - 1
                i = iw + gp/2 + 1
                j = jw + gp/2 + 1
                k = kw + gp/2 + 1

                k3 = info%kv(:,1)*iw + info%kv(:,2)*jw + info%kv(:,3)*kw
                kn = NORM2(k3)

!               Truncate
                IF(kn .eq. 0 .or. kn .gt. 40) CYCLE

!               Amplification factor, gets multiplied onto spectral point forces
                B = BT(k3, info%xi, kn, info%eta)*ii
                Hht(:,:, i, j, k) = B(:,:,1)*Hh(1,i,j,k) &
                                  + B(:,:,2)*Hh(2,i,j,k) &
                                  + B(:,:,3)*Hh(3,i,j,k)

            ENDDO
        ENDDO
    ENDDO
    Ht = 0D0
    ! Double layer... DON'T FORGET TO MULTIPLY ii!!!!!! TAU????

!!!!!!!!!!!! Note: FFTshift only for even number grid. Need to make iFFTshift for odds
    DO i = 1,3
        DO j = 1,3
            Htmp = Hht(i,j,:,:,:)
            CALL FFTSHIFT(Htmp)
            Htmp = iFFT3(Htmp, info%WSAVE)
            Ht(i,j,:,:,:) = REAL(Htmp)
        ENDDO
    ENDDO
END SUBROUTINE HtT

! -------------------------------------------------------------------------!
! Given an Ht grid as above and a set of points, integrates back to the points
! For single layer
SUBROUTINE HtIntG(Ht, info, x0, u)
    REAL(KIND = 8), INTENT(IN) :: Ht(:,:,:,:), x0(:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    REAL(KIND = 8), INTENT(OUT) :: u(:,:)
    REAL(KIND = 8) :: Jbv, xcur(3), rr(3), r2, &
                      bv_gr(3,3), bvi_gr(3,3), h
    INTEGER :: nt, i, curijk(3), inds(3), &
               iper, jper, kper, pts, gp

    nt = size(x0)/3
    gp = info%gp
    Jbv = DET3(info%bv)
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)
    h = 1D0/REAL(gp)
    u = 0D0

!   This is the "interpolate back" step, where we have and send it back to discrete points
    DO i = 1, nt
!       At a point, find grid index & loop over supporting grid points
        xcur = x0(:,i)
        curijk(1) = FLOOR(bvi_gr(1,1)*xcur(1) + bvi_gr(1,2)*xcur(2) + bvi_gr(1,3)*xcur(3))
        curijk(2) = FLOOR(bvi_gr(2,1)*xcur(1) + bvi_gr(2,2)*xcur(2) + bvi_gr(2,3)*xcur(3))
        curijk(3) = FLOOR(bvi_gr(3,1)*xcur(1) + bvi_gr(3,2)*xcur(2) + bvi_gr(3,3)*xcur(3))
        DO pts = 1,info%suppPoints
            inds = curijk + info%suppmat(:, pts)

!           Manage periodicity
            IF(inds(1) .lt. 1) THEN
                iper = inds(1) + gp
            ELSEIF(inds(1) .gt. gp) THEN
                iper = inds(1) - gp
            ELSE
                iper = inds(1)
            ENDIF
            IF(inds(2) .lt. 1) THEN
                jper = inds(2) + gp
            ELSEIF(inds(2) .gt. gp) THEN
                jper = inds(2) - gp
            ELSE
                jper = inds(2)
            ENDIF
            IF(inds(3) .lt. 1) THEN
                kper = inds(3) + gp
            ELSEIF(inds(3) .gt. gp) THEN
                kper = inds(3) - gp
            ELSE
                kper = inds(3)
            ENDIF

!           Add smeared point forces to grid point
            rr(1) = xcur(1) - (bv_gr(1,1)*inds(1) + bv_gr(1,2)*inds(2) + bv_gr(1,3)*inds(3)) 
            rr(2) = xcur(2) - (bv_gr(2,1)*inds(1) + bv_gr(2,2)*inds(2) + bv_gr(2,3)*inds(3)) 
            rr(3) = xcur(3) - (bv_gr(3,1)*inds(1) + bv_gr(3,2)*inds(2) + bv_gr(3,3)*inds(3)) 

            r2 = rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3)

            u(:,i) = u(:,i) + Ht(:, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv
        ENDDO
    ENDDO
END SUBROUTINE HtIntG

! -------------------------------------------------------------------------!
! Given an Ht grid as above and a set of points, integrates back to the points
! For double layer
SUBROUTINE HtIntT(Ht, info, x0, u)
    REAL(KIND = 8), INTENT(IN) :: Ht(:,:,:,:,:), x0(:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    REAL(KIND = 8), INTENT(OUT) :: u(:,:,:)
    REAL(KIND = 8) :: Jbv, xcur(3), rr(3), r2, &
                      bv_gr(3,3), bvi_gr(3,3), h
    INTEGER :: nt, i, curijk(3), inds(3), &
               iper, jper, kper, pts, gp

    nt = size(x0)/3
    gp = info%gp
    Jbv = DET3(info%bv)
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)
    h = 1D0/REAL(gp)
    u = 0D0

!   This is the "interpolate back" step, where we have and send it back to discrete points
    DO i = 1, nt
!       At a point, find grid index & loop over supporting grid points
        xcur = x0(:,i)
        curijk(1) = FLOOR(bvi_gr(1,1)*xcur(1) + bvi_gr(1,2)*xcur(2) + bvi_gr(1,3)*xcur(3))
        curijk(2) = FLOOR(bvi_gr(2,1)*xcur(1) + bvi_gr(2,2)*xcur(2) + bvi_gr(2,3)*xcur(3))
        curijk(3) = FLOOR(bvi_gr(3,1)*xcur(1) + bvi_gr(3,2)*xcur(2) + bvi_gr(3,3)*xcur(3))
        DO pts = 1,info%suppPoints
            inds = curijk + info%suppmat(:, pts)

!           Manage periodicity
            IF(inds(1) .lt. 1) THEN
                iper = inds(1) + gp
            ELSEIF(inds(1) .gt. gp) THEN
                iper = inds(1) - gp
            ELSE
                iper = inds(1)
            ENDIF
            IF(inds(2) .lt. 1) THEN
                jper = inds(2) + gp
            ELSEIF(inds(2) .gt. gp) THEN
                jper = inds(2) - gp
            ELSE
                jper = inds(2)
            ENDIF
            IF(inds(3) .lt. 1) THEN
                kper = inds(3) + gp
            ELSEIF(inds(3) .gt. gp) THEN
                kper = inds(3) - gp
            ELSE
                kper = inds(3)
            ENDIF

!           Add smeared point forces to grid point
            rr(1) = xcur(1) - (bv_gr(1,1)*inds(1) + bv_gr(1,2)*inds(2) + bv_gr(1,3)*inds(3)) 
            rr(2) = xcur(2) - (bv_gr(2,1)*inds(1) + bv_gr(2,2)*inds(2) + bv_gr(2,3)*inds(3)) 
            rr(3) = xcur(3) - (bv_gr(3,1)*inds(1) + bv_gr(3,2)*inds(2) + bv_gr(3,3)*inds(3)) 

            r2 = rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3)

            u(:,:,i) = u(:,:,i) + Ht(:, :, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv
        ENDDO
    ENDDO
END SUBROUTINE HtIntT

! ! -------------------------------------------------------------------------!
! Double Layer amplification factor (all unrolled)
FUNCTION BT(k, xi, kn, eta) RESULT(B)
    REAL(KIND = 8) :: kn, k(3), xi, eta, B(3,3,3), D(3,3,3)
    D(1,1,1) = -2D0/(kn*kn*kn*kn)*k(1)*k(1)*k(1) + 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(2,1,1) = -2D0/(kn*kn*kn*kn)*k(2)*k(1)*k(1) + 1/(kn*kn)*(       k(2)       )
    D(3,1,1) = -2D0/(kn*kn*kn*kn)*k(3)*k(1)*k(1) + 1/(kn*kn)*(       k(3)       )
    D(1,2,1) = -2D0/(kn*kn*kn*kn)*k(1)*k(2)*k(1) + 1/(kn*kn)*(              k(2))
    D(2,2,1) = -2D0/(kn*kn*kn*kn)*k(2)*k(2)*k(1) + 1/(kn*kn)*(k(1)              )
    D(3,2,1) = -2D0/(kn*kn*kn*kn)*k(3)*k(2)*k(1)!+ 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(1,3,1) = -2D0/(kn*kn*kn*kn)*k(1)*k(3)*k(1) + 1/(kn*kn)*(              k(3))
    D(2,3,1) = -2D0/(kn*kn*kn*kn)*k(2)*k(3)*k(1)!+ 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(3,3,1) = -2D0/(kn*kn*kn*kn)*k(3)*k(3)*k(1) + 1/(kn*kn)*(k(1)              )
    D(1,1,2) = -2D0/(kn*kn*kn*kn)*k(1)*k(1)*k(2) + 1/(kn*kn)*(k(2)              )
    D(2,1,2) = -2D0/(kn*kn*kn*kn)*k(2)*k(1)*k(2) + 1/(kn*kn)*(              k(1))
    D(3,1,2) = -2D0/(kn*kn*kn*kn)*k(3)*k(1)*k(2)!+ 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(1,2,2) = -2D0/(kn*kn*kn*kn)*k(1)*k(2)*k(2) + 1/(kn*kn)*(       k(1)       )
    D(2,2,2) = -2D0/(kn*kn*kn*kn)*k(2)*k(2)*k(2) + 1/(kn*kn)*(k(2) + k(2) + k(2))
    D(3,2,2) = -2D0/(kn*kn*kn*kn)*k(3)*k(2)*k(2) + 1/(kn*kn)*(       k(3)       )
    D(1,3,2) = -2D0/(kn*kn*kn*kn)*k(1)*k(3)*k(2)!+ 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(2,3,2) = -2D0/(kn*kn*kn*kn)*k(2)*k(3)*k(2) + 1/(kn*kn)*(              k(3))
    D(3,3,2) = -2D0/(kn*kn*kn*kn)*k(3)*k(3)*k(2) + 1/(kn*kn)*(k(2)              )
    D(1,1,3) = -2D0/(kn*kn*kn*kn)*k(1)*k(1)*k(3) + 1/(kn*kn)*(k(3)              )
    D(2,1,3) = -2D0/(kn*kn*kn*kn)*k(2)*k(1)*k(3)!+ 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(3,1,3) = -2D0/(kn*kn*kn*kn)*k(3)*k(1)*k(3) + 1/(kn*kn)*(              k(1))
    D(1,2,3) = -2D0/(kn*kn*kn*kn)*k(1)*k(2)*k(3)!+ 1/(kn*kn)*(k(1) + k(1) + k(1))
    D(2,2,3) = -2D0/(kn*kn*kn*kn)*k(2)*k(2)*k(3) + 1/(kn*kn)*(k(3)              )
    D(3,2,3) = -2D0/(kn*kn*kn*kn)*k(3)*k(2)*k(3) + 1/(kn*kn)*(              k(2))
    D(1,3,3) = -2D0/(kn*kn*kn*kn)*k(1)*k(3)*k(3) + 1/(kn*kn)*(       k(1)       )
    D(2,3,3) = -2D0/(kn*kn*kn*kn)*k(2)*k(3)*k(3) + 1/(kn*kn)*(       k(2)       )
    D(3,3,3) = -2D0/(kn*kn*kn*kn)*k(3)*k(3)*k(3) + 1/(kn*kn)*(k(3) + k(3) + k(3))
    B = -8D0*pi*(1 + kn*kn*0.25D0/(xi*xi))*D*EXP((eta - 1D0)*kn*kn*0.25D0/(xi*xi))
END FUNCTION BT

! -------------------------------------------------------------------------!
! Single Layer amplification factor
FUNCTION BG(k, xi, kn, eta, eye) RESULT(B)
    REAL(KIND = 8) :: kn, k(3), xi, eta, B(3,3), eye(3,3)
    B = 8D0*pi*(1 + kn*kn*0.25D0/(xi*xi))/(kn*kn*kn*kn)*(kn*kn*eye - OUTER(k,k))
    B = B*EXP((eta - 1D0)*kn*kn*0.25D0/(xi*xi))
END FUNCTION BG

! -------------------------------------------------------------------------!
! Single Layer amplification factor
SUBROUTINE SuppPoints(info)
    TYPE(sharedType), INTENT(INOUT), POINTER :: info
    REAL(KIND = 8) :: gp, bv_gr(3,3), rr(3), r, cut, cutpar
    INTEGER :: i, j, k, it
    INTEGER, ALLOCATABLE :: tmpmat(:,:)

    it = 0
!!!!! Maybe wnt to put cut into info?
    cut = 1D-12
    cutpar = SQRT(-LOG(cut)*info%eta/(info%xi*info%xi*2D0))

    gp = info%gp
    bv_gr = info%bv/REAL(gp)
    ALLOCATE(tmpmat(3,FLOOR(gp*gp*gp)))
    IF(ALLOCATED(info%suppMat)) DEALLOCATE(info%suppMat)

!   Find points that fall within cutoff distance
    DO i = -FLOOR(gp/2),FLOOR(gp/2)
        DO j = -FLOOR(gp/2),FLOOR(gp/2)
            DO k = -FLOOR(gp/2),FLOOR(gp/2)
                rr = bv_gr(:,1)*i + bv_gr(:,2)*j + bv_gr(:,3)*k
                r = NORM2(rr)
                IF(r .lt. cutpar) THEN
                    it = it + 1
                    tmpmat(:,it) = (/i,j,k/)
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    ALLOCATE(info%suppMat(3,it))
    info%suppMat = tmpmat(:,:it)
    info%suppPoints = it
END SUBROUTINE SuppPoints

! -------------------------------------------------------------------------!
! Test to see if things are working
SUBROUTINE PeriodicTest(info)
    TYPE(sharedType), INTENT(INOUT), POINTER :: info
    TYPE(YType), POINTER :: Y
    COMPLEX(KIND = 8) :: f(3,2)
    REAL(KIND = 8) ::x0(3,2), Ht(3,info%gp,info%gp,info%gp), HtT(3,3,info%gp,info%gp,info%gp)
    REAL(KIND = 8), ALLOCATABLE :: x(:,:,:,:), xv(:,:), Jtmp(:,:), x2(:,:), u2(:,:), u3(:,:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fs(:,:,:,:), fv1(:), fv2(:), fv3(:), ft(:,:)
    INTEGER :: i, j, ic, nt, np

    Y => info%Y
    
!   Support points
    CALL SuppPoints(info)
    ALLOCATE(x2(3,1), u2(3,1), u3(3,3,1))
    x2(:,1) = (/0.5D0, 0.75D0, 0.5D0/)

!   Two point forces
    f=0D0
    f(1,1) = 1D0
    f(1,2) = -1D0
    x0 = 0.5D0
    x0(1,1) = 0.25D0
    x0(1,2) = 0.75D0

!   Get the integration matrix
    CALL HtCalc(Ht, info, x0, f)
    CALL HtCalc(HtT, info, x0, f)
    CALL HtInt (Ht, info, x2, u2)
    print *, u2
    print *, ' '
    CALL HtInt(HtT, info, x2, u3)
    print *, u3
    print *, ' '
    print *, HtT(:,:,16,16,16)
    print *, ' '

!   Now let's do a test on two spheres
    nt = Y%nt
    np = Y%np
    ALLOCATE(x(2, 3, nt, np), fs(2, 3, nt, np), Jtmp(nt, np),ft(3, nt*np*2))
    fs = 0D0

    ic = 0
!   Sphere setup
    DO i = 1,nt
        DO j = 1,np
            ic = ic + 1
            x(1,1,i,j) = SIN(Y%tht(i))*COS(Y%phi(j))/4D0 + 0.5D0
            x(1,2,i,j) = SIN(Y%tht(i))*SIN(Y%phi(j))/4D0 + 0.5D0
            x(1,3,i,j) = COS(Y%tht(i))/4D0 + 0.5D0
            
            Jtmp(i,j) = SIN(Y%tht(i))*(0.25D0*0.25D0)

            fs(1,1,i,j) = CMPLX(MODULO(ic,2)*2-1)
        ENDDO
    ENDDO

    fs(2,:,:,:) = -fs(1,:,:,:)*0D0

    x(2,1,:,:) = x(1,1,:,:)
    x(2,3,:,:) = x(1,2,:,:)

!   Now the actual process of running the calculations, starting with putting matrices into vectors
    DO ic = 1,2
        CALL HtPreCalc(fs(ic,1,:,:), fv1, Jtmp, Y, x(ic,:,:,:), xv)
        CALL HtPreCalc(fs(ic,2,:,:), fv2, Jtmp, Y)
        CALL HtPreCalc(fs(ic,3,:,:), fv3, Jtmp, Y)
    ENDDO
    ft(1,:) = fv1
    ft(2,:) = fv2
    ft(3,:) = fv3
    
    CALL HtG(Ht,info, xv, ft)
    CALL HtInt (Ht, info, x2, u2)
    print *, u2

!   Try it with evenly spaced points
    DEALLOCATE(xv, ft)
    ALLOCATE(xv(3,100), ft(3,100))
    ft = 0D0
    ic = 0
    DO i = 1,10
        DO j = 1,10
            ic = ic+1
            xv(1,ic) = (i-1)*0.1D0
            xv(2,ic) = (j-1)*0.1D0
            ft(1,ic) = CMPLX(MODULO(ic,2)*2-1)
        ENDDO
    ENDDO
    xv(3,:) = 0.5D0

    CALL HtCalc(Ht,info, xv, ft)
    CALL HtInt (Ht, info, x2, u2)
    print *, u2


END SUBROUTINE PeriodicTest

END MODULE PERIODICMOD