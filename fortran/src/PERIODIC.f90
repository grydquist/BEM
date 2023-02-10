MODULE PERIODICMOD
USE SHAREDMOD
USE CMMOD
IMPLICIT NONE
    
!==============================================================================!
!               The purpose of this module is to the calculations              !
!                specifically related to the fast Ewald sum part               !
!                                                                              !
!==============================================================================!

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!
INTERFACE Ewaldint
PROCEDURE EwaldintG, EwaldintT
END INTERFACE Ewaldint

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

! -------------------------------------------------------------------------!
! Some precalculations that only need to be done once per timestep to run HtCell
! Mostly taking cell forces and locations and constructing a single vector
! Only need to do this once for x0, so it's optional.
! Note that this will append fin to f
SUBROUTINE EwaldPreCalc(f, Y, x0, fin, x0in, Ja, fmn, xmn)
    COMPLEX(KIND = 8), INTENT(IN), OPTIONAL :: fmn(:), xmn(:,:)
    REAL(KIND = 8), ALLOCATABLE, INTENT(INOUT):: f(:)
    REAL(KIND = 8), INTENT(IN), OPTIONAL :: fin(:,:), Ja(:,:)
    TYPE(YType), INTENT(IN) :: Y
    REAL(KIND = 8), INTENT(IN), OPTIONAL :: x0in(:,:,:)
    REAL(KIND = 8), INTENT(INOUT), ALLOCATABLE, OPTIONAL:: x0(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fapp(:), ftmp(:), fint(:,:)
    REAL(KIND = 8), ALLOCATABLE :: xapp(:,:), xtmp(:,:), x0int(:,:,:)
    INTEGER :: shp(2), nt, np, i, j, it, tot

!   We want to be able to take in both force and location harmonics
!   but this needs some exception handling
    IF(PRESENT(fin) .and. PRESENT(fmn)) THEN
        print *, 'ERROR: For pre-calculation of Ewald sum,', &
                 'cant have both force harmonic constants and scalars'
        STOP
    ELSEIF(PRESENT(fin)) THEN
        shp = SHAPE(fin)
        nt = shp(1)
        np = shp(2)
        ALLOCATE(fint(nt, np))
        fint = fin
    ELSEIF(PRESENT(fmn)) THEN
        nt = Y%nt
        np = Y%np
        ALLOCATE(fint(nt, np))
        fint = Y%backward(fmn, Y%p)
    ELSE
        print *, 'ERROR: For pre-calculation of Ewald sum,', &
                 'must have either force harmonic constants and scalars'
        STOP
    ENDIF

    IF(PRESENT(x0in) .and. PRESENT(xmn)) THEN
        print *, 'ERROR: For pre-calculation of Ewald sum,', &
                 'cant have both location harmonic constants and scalars'
        STOP
    ELSEIF(PRESENT(xmn)) THEN
        ALLOCATE(x0int(3,nt, np))
        x0int(1,:,:) = Y%backward(xmn(1,:), Y%p)
        x0int(2,:,:) = Y%backward(xmn(2,:), Y%p)
        x0int(3,:,:) = Y%backward(xmn(3,:), Y%p)
    ELSEIF(PRESENT(x0in)) THEN
        ALLOCATE(x0int(3,nt, np))
        x0int = x0in
    ENDIF

    IF(ALLOCATED(f)) THEN
        tot = SIZE(f)
    ELSE
        tot = 0
    ENDIF
    ALLOCATE(fapp(nt*np), xapp(3, nt*np), ftmp(nt*np + tot), xtmp(3, nt*np + tot))
    it = 0

!   Replace with "effective" point force strength
    DO i = 1,nt
        DO j = 1,np
            it = it + 1
!           Area element often built into f
            IF(PRESENT(Ja)) THEN
                fapp(it) = fint(i,j)*Y%dphi*Y%wg(i)*Ja(i,j)
            ELSE
                fapp(it) = fint(i,j)*Y%dphi*Y%wg(i)
            ENDIF
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
    IF(PRESENT(x0in) .or. PRESENT(xmn)) THEN
        it = 0 
        DO i = 1,nt
            DO j = 1,np
                it = it + 1
                xapp(:,it) = x0int(:,i,j)
            ENDDO
        ENDDO
        IF(ALLOCATED(x0)) xtmp(:, 1:tot) = x0
        xtmp(:, tot + 1: tot + nt*np) = xapp
        IF(ALLOCATED(x0)) DEALLOCATE(x0)
        ALLOCATE(x0(3, nt*np + tot))
        x0 = xtmp
    ENDIF

END SUBROUTINE EwaldPreCalc

! -------------------------------------------------------------------------!
! For a given set of cells with point forces, construct the H_tilde grid
! for the single layer
! See: Spectrally accurate fast summation for periodic Stokes potentials
!       D. Lindbo, AK Tornberg
! Can be used to interpolate back to grid.
SUBROUTINE EwaldG(Ht, info, x0, f, strt, full, u3, u1)
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: Ht(:,:,:,:), u3(:,:), u1(:)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:)
    REAL(KIND = 8) :: bvi_gr(3,3), bv_gr(3,3), rr(3), r2, k3(3), &
                      kn, B(3,3), xcur(3)
    REAL(KIND = 8), ALLOCATABLE :: wrk(:), Htout(:,:,:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    REAL(KIND = 8), INTENT(IN) :: f(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: full
    INTEGER, INTENT(IN), OPTIONAL :: strt
    COMPLEX(KIND = 8), ALLOCATABLE ::  H(:,:,:,:), Hh(:,:,:,:), Hht(:,:,:,:), Htmp(:,:,:)
    INTEGER :: i, j, k, pts, curijk(3), gp, inds(3), &
               iper, jper, kper, iw, jw, kw, nt, lwrk

    gp = info%gp

    lwrk = gp*2*3
    ALLOCATE(wrk(lwrk))

!   Total number of points
    nt = size(x0)/3

    ALLOCATE(H(3,gp,gp,gp), Hh(3,gp,gp,gp), &
           Hht(3,gp,gp,gp), Htmp(gp,gp,gp), Htout(3,gp,gp,gp))
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

    Hht = 0D0
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
                IF(kn .eq. 0) CYCLE

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

!!!! Note: FFTshift only for even number grid. Need to make iFFTshift for odds
    Htmp = Hht(1,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Htout(1,:,:,:) = REAL(Htmp)

    Htmp = Hht(2,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Htout(2,:,:,:) = REAL(Htmp)

    Htmp = Hht(3,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Htout(3,:,:,:) = REAL(Htmp)

    IF(PRESENT(Ht)) Ht = Htout

!   If we want, we can just do the next step and integrate
    IF(PRESENT(full) .and. full) THEN
        IF(.not. PRESENT(u3) .and. .not.PRESENT(u1)) THEN
            print *, "Warning: no output matrix for full Ewald G, exiting routine"
            RETURN
        ENDIF
        IF(PRESENT(u3) .and. PRESENT(u1)) THEN
            print *, "Warning: too many output matrices for full Ewald G, exiting routine"
            RETURN
        ENDIF
        IF(PRESENT(u1)) THEN
            CALL Ewaldint(Htout, info, x0, u1=u1)
        ENDIF
        IF(PRESENT(u3)) THEN
            CALL Ewaldint(Htout, info, x0, u3=u3)
        ENDIF
    ENDIF
END SUBROUTINE EwaldG

! -------------------------------------------------------------------------!
! For a given set of cells with point forces, construct the H_tilde grid
! for the single layer
! See: Spectrally accurate fast summation for periodic Stokes potentials
!       D. Lindbo, AK Tornberg
! Can be used to interpolate back to grid.
! Note that in most cases, f will be the normal vector x a sphHarm
SUBROUTINE EwaldTold(Ht, info, x0, f, full, um, u3, strt)
    COMPLEX(KIND = 8), INTENT(OUT), OPTIONAL :: Ht(:,:,:,:,:)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    COMPLEX(KIND = 8), INTENT(IN) :: f(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: full
    INTEGER, INTENT(IN), OPTIONAL :: strt
    COMPLEX(KIND = 8), INTENT(OUT), OPTIONAL :: um(:,:,:,:,:), u3(:,:,:)
    COMPLEX(KIND = 8), ALLOCATABLE ::  H(:,:,:,:), Hh(:,:,:,:), &
        Hht(:,:,:,:,:), Htmp(:,:,:)
    COMPLEX(KIND = 8) :: B(3,3,3)
    REAL(KIND = 8) :: bvi_gr(3,3), bv_gr(3,3), rr(3), r2, k3(3), &
                      kn, xcur(3)
    REAL(KIND = 8), ALLOCATABLE :: wrk(:)
    INTEGER :: i, j, k, pts, curijk(3), gp, inds(3), &
               iper, jper, kper, iw, jw, kw, nt, lwrk, strti
    gp = info%gp

    lwrk = gp*2*3
    ALLOCATE(wrk(lwrk))

!   Total number of points
    nt = size(f)/3

    ALLOCATE(H(3,gp,gp,gp), Hh(3,gp,gp,gp), &
           Hht(3,3,gp,gp,gp), Htmp(gp,gp,gp))
    H = 0D0
    IF(PRESENT(strt)) strti=strt

!   Subgrid basis vectors and inverse
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)

!   Now we can start the actual construction of the Ht matrix. Start with the H matrix,
!   which is just the point forces smeared onto grid points. Go to each force point and
!   smear to supporting points
    DO i = 1, nt
!       At a point, find grid index & loop over supporting grid points
!       Since we only do this for one surface, start at the points at this surface
        xcur = x0(:,i + (strti-1)*nt)
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

    Hht = 0D0
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
                IF(kn .eq. 0) CYCLE

!               Amplification factor, gets multiplied onto spectral point forces
                B = -BT(k3, info%xi, kn, info%eta)*ii
                Hht(:,:, i, j, k) = B(:,:,1)*Hh(1,i,j,k) &
                                  + B(:,:,2)*Hh(2,i,j,k) &
                                  + B(:,:,3)*Hh(3,i,j,k)

            ENDDO
        ENDDO
    ENDDO

!!!! Note: FFTshift only for even number grid. Need to make iFFTshift for odds
    DO i = 1,3
        DO j = 1,3
            Htmp = Hht(i,j,:,:,:)
            CALL FFTSHIFT(Htmp)
            Htmp = iFFT3(Htmp, info%WSAVE)
            Hht(i,j,:,:,:) = Htmp
        ENDDO
    ENDDO

    IF(PRESENT(Ht)) Ht = Hht

!   If we want, we can just do the next step and integrate
    IF(PRESENT(full) .and. full) THEN
        IF(.not. PRESENT(u3) .and. .not.PRESENT(um)) THEN
            print *, "Warning: no output matrix for full Ewald G, exiting routine"
            RETURN
        ENDIF
        IF(PRESENT(u3) .and. PRESENT(um)) THEN
            print *, "Warning: too many output matrices for full Ewald G, exiting routine"
            RETURN
        ENDIF
        IF(PRESENT(um)) THEN
            ! CALL Ewaldint(Hht, info, x0, um=um)
        ENDIF
        IF(PRESENT(u3)) THEN
            ! CALL Ewaldint(Hht, info, x0, u3=u3)
        ENDIF
    ENDIF
END SUBROUTINE EwaldTold

! -------------------------------------------------------------------------!
! For a given set of cells with point forces, construct the H_tilde grid
! for the double layer
! See: Spectrally accurate fast summation for periodic Stokes potentials
!       D. Lindbo, AK Tornberg
! Can be used to interpolate back to grid.
SUBROUTINE EwaldT(Ht, info, x0, f1, f3, n, full, u3, u1, strt, cm)
    COMPLEX(KIND = 8), INTENT(OUT), OPTIONAL :: Ht(:,:,:,:)
    COMPLEX(KIND = 8) :: B(3,3,3)
    COMPLEX(KIND = 8), ALLOCATABLE ::  H(:,:,:,:,:), Hh(:,:,:,:,:), &
        Hht(:,:,:,:), Htmp(:,:,:)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:),  n(:,:)
    REAL(KIND = 8), INTENT(IN) , OPTIONAL :: f1(:), f3(:,:)
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: u1(:), u3(:,:,:)
    REAL(KIND = 8), ALLOCATABLE :: f(:,:)
    REAL(KIND = 8) :: bvi_gr(3,3), bv_gr(3,3), rr(3), r2, k3(3), &
                      kn, xcur(3)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    TYPE(cmType), INTENT(IN), POINTER :: cm
    LOGICAL, INTENT(IN), OPTIONAL :: full
    INTEGER, INTENT(IN), OPTIONAL :: strt
    INTEGER :: i, j, k, pts, curijk(3), gp, inds(3), &
               iper, jper, kper, iw, jw, kw, nt, strti, i1, i2
    gp = info%gp

    IF(PRESENT(f1)) THEN
        IF(PRESENT(f3)) THEN
            print *, "ERROR: too many arguments for force for EwaldT"
            stop
        ENDIF
        ALLOCATE(f(3, SIZE(f1)/3))
        j = 1
        DO i = 1, SIZE(f1)/3
            f(1,i) = f1(j)
            f(2,i) = f1(j + 1)
            f(3,i) = f1(j + 2)
            j = j + 3
        ENDDO
    ELSEIF(.not.PRESENT(f3)) THEN
        print *, "ERROR: No forcing argument for EwaldT"
        stop
    ELSE
        f = f3
    ENDIF

!   Total number of points
    nt = SIZE(f)/3

    i1 = (info%PCells(1,1) - 1)*info%Y%nt*info%Y%np + (info%PCells(2,1) - 1)*info%Y%np + 1 ! 1 !
    i2 = (info%PCells(1,2) - 1)*info%Y%nt*info%Y%np + (info%PCells(2,2)    )*info%Y%np ! nt !

    ALLOCATE(H(3,3,gp,gp,gp), Hh(3,3,gp,gp,gp), &
           Hht(3,gp,gp,gp), Htmp(gp,gp,gp))
    H = 0D0
    IF(PRESENT(strt)) strti=strt

!   Subgrid basis vectors and inverse
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)

!   Now we can start the actual construction of the Ht matrix. Start with the H matrix,
!   which is just the point forces smeared onto grid points. Go to each force point and
!   smear to supporting points
    DO i = i1, i2
!       At a point, find grid index & loop over supporting grid points
!       Since we only do this for one surface, start at the points at this surface
        xcur = x0(:,i + (strti-1)*nt)
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

            H(:,:,iper,jper,kper) = H(:,:,iper,jper,kper) &
                                  + OUTER(f(:,i), n(:,i))*info%par*EXP(info%parexp*r2)
        ENDDO
    ENDDO

    H = cm%reduce(H)

!   This is a constant cost (same size FFT regardless of problem size)
!   so perhaps no need to parallelize
    DO i = 1,3
        DO j = 1,3
            Htmp = H(i,j,:,:,:)
            Htmp = FFT3(Htmp, info%WSAVE)
            CALL FFTSHIFT(Htmp)
            Hh(i,j,:,:,:) = Htmp
        ENDDO
    ENDDO

    Hht = 0D0
!   Now perform the convolution with the kernel (which is a truncated sum in spectral space)
!   (gp needs to be even)
!   Also assumes the Hh matrix has negative and positive parts
!   Again, constant cost
    DO iw = -gp/2, gp/2 - 1
        DO jw = -gp/2, gp/2 - 1
            DO kw = -gp/2, gp/2 - 1
                i = iw + gp/2 + 1
                j = jw + gp/2 + 1
                k = kw + gp/2 + 1

                k3 = info%kv(:,1)*iw + info%kv(:,2)*jw + info%kv(:,3)*kw
                kn = NORM2(k3)

!               Truncate
                IF(kn .eq. 0) CYCLE

!               Amplification factor, gets multiplied onto spectral point forces
                B = -BT(k3, info%xi, kn, info%eta)*ii
                Hht(:, i, j, k) = B(1,:,1)*Hh(1,1,i,j,k) &
                                + B(1,:,2)*Hh(1,2,i,j,k) &
                                + B(1,:,3)*Hh(1,3,i,j,k) &
                                + B(2,:,1)*Hh(2,1,i,j,k) &
                                + B(2,:,2)*Hh(2,2,i,j,k) &
                                + B(2,:,3)*Hh(2,3,i,j,k) &
                                + B(3,:,1)*Hh(3,1,i,j,k) &
                                + B(3,:,2)*Hh(3,2,i,j,k) &
                                + B(3,:,3)*Hh(3,3,i,j,k) 

            ENDDO
        ENDDO
    ENDDO

!!!! Note: FFTshift only for even number grid. Need to make iFFTshift for odds

    Htmp = Hht(1,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Hht(1,:,:,:) = Htmp

    Htmp = Hht(2,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Hht(2,:,:,:) = Htmp

    Htmp = Hht(3,:,:,:)
    CALL FFTSHIFT(Htmp)
    Htmp = iFFT3(Htmp, info%WSAVE)
    Hht(3,:,:,:) = Htmp

    IF(PRESENT(Ht)) Ht = Hht

!   If we want, we can just do the next step and integrate
    IF(PRESENT(full) .and. full) THEN
        IF(.not. PRESENT(u3) .and. .not.PRESENT(u1)) THEN
            print *, "Warning: no output matrix for full Ewald G, exiting routine"
            RETURN
        ENDIF
        IF(PRESENT(u3) .and. PRESENT(u1)) THEN
            print *, "Warning: too many output matrices for full Ewald G, exiting routine"
            RETURN
        ENDIF
        IF(PRESENT(u1)) THEN
            CALL Ewaldint(Hht, info, x0, u1=u1, cm=cm)
        ENDIF
        IF(PRESENT(u3)) THEN
            CALL Ewaldint(Hht, info, x0, u3=u3, cm=cm)
        ENDIF
    ENDIF
END SUBROUTINE EwaldT

! -------------------------------------------------------------------------!
! Given an Ht grid as above and a set of points, integrates back to the points
! For single layer
SUBROUTINE EwaldintG(Ht, info, x0, u3, u1)
    REAL(KIND = 8), INTENT(IN) :: Ht(:,:,:,:), x0(:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: u3(:,:), u1(:)
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

!   Sometimes we want the output in different ranks
    IF(PRESENT(u3)) u3 = 0D0
    IF(PRESENT(u1)) u1 = 0D0
    IF(PRESENT(u3) .and. PRESENT(u1)) THEN
        print *, 'ERROR: For integration on periodic Fourier matrix,',  &
                 'choose rank 1 or rank 2 matrices, not both'
        stop
    ENDIF

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

            IF(PRESENT(u3)) u3(:,i) = u3(:,i) &
                           + Ht(:, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv

            IF(PRESENT(u1)) u1(3*(i-1) + 1:3*(i-1) + 3) = u1(3*(i-1) + 1:3*(i-1) + 3) &
                           + Ht(:, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv
        ENDDO
    ENDDO
END SUBROUTINE EwaldintG

! -------------------------------------------------------------------------!
! Given an Ht grid as above and a set of points, integrates back to the points
! For double layer
SUBROUTINE EwaldintT(Ht, info, x0, u3, u1, cm)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:)
    COMPLEX(KIND = 8), OPTIONAL, INTENT(IN) :: Ht(:,:,:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    TYPE(cmType), INTENT(IN), OPTIONAL, POINTER :: cm
    REAL(KIND = 8), OPTIONAL, INTENT(OUT) :: u3(:,:,:), u1(:)
    REAL(KIND = 8) :: Jbv, xcur(3), rr(3), r2, &
                      bv_gr(3,3), bvi_gr(3,3), h
    INTEGER :: nt, np, nc, tpts, i, curijk(3), inds(3), &
               iper, jper, kper, pts, gp, it, ip, ic, i1, i2

    tpts = size(x0)/3
    gp = info%gp
    Jbv = DET3(info%bv)
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)
    h = 1D0/REAL(gp)

!   Get indices for looping parallel (was previous way to parallelize)
!   The master list of points goes cell -> theta -> phi
    i1 = (info%PCells(1,1) - 1)*info%Y%nt*info%Y%np + (info%PCells(2,1) - 1)*info%Y%np + 1 ! 1
    i2 = (info%PCells(1,2) - 1)*info%Y%nt*info%Y%np + (info%PCells(2,2)    )*info%Y%np ! tpts

!   Often, we want an output matrix of size (3,3,nt,np,ncell) as output.
!   Exception checking for this option
    IF(PRESENT(u3)) THEN
        print *, 'ERROR:3D vector not currently supported for Ewald back-integration'
        stop
    ENDIF
    IF(PRESENT(u1)) THEN
        u1 = 0D0
    ENDIF
    IF(PRESENT(u3) .and. PRESENT(u1)) THEN
        print *, 'ERROR: For integration on periodic fourier matrix,',  &
                 'choose rank 3 or rank 5 matrices, not both'
        stop
    ENDIF

!   This is the "interpolate back" step, where we have and send it back to discrete points
    DO i = i1, i2
        
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

            IF(PRESENT(u1)) &
            u1(3*(i-1) + 1:3*(i-1) + 3) = u1(3*(i-1) + 1:3*(i-1) + 3) &
                           + Ht(:, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv

        ENDDO
    ENDDO
    u1 = cm%reduce(u1)
END SUBROUTINE EwaldintT

! -------------------------------------------------------------------------!
! Given an Ht grid as above and a set of points, integrates back to the points
! For double layer
SUBROUTINE EwaldintTold(Ht, info, x0, u3, um)
    REAL(KIND = 8), INTENT(IN) :: x0(:,:)
    COMPLEX(KIND = 8), OPTIONAL, INTENT(IN) :: Ht(:,:,:,:,:)
    TYPE(sharedType), INTENT(IN), POINTER :: info
    COMPLEX(KIND = 8), OPTIONAL, INTENT(OUT) :: u3(:,:,:), um(:,:,:,:,:)
    REAL(KIND = 8) :: Jbv, xcur(3), rr(3), r2, &
                      bv_gr(3,3), bvi_gr(3,3), h
    INTEGER :: nt, np, nc, tpts, i, curijk(3), inds(3), &
               iper, jper, kper, pts, gp, sz(5), it, ip, ic, i1, i2

    tpts = size(x0)/3
    gp = info%gp
    Jbv = DET3(info%bv)
    bv_gr = info%bv/REAL(gp)
    bvi_gr = INVERT33(bv_gr)
    h = 1D0/REAL(gp)

!   Get indices for looping parallel (was previous way to parallelize)
!   The master list of points goes cell -> theta -> phi
    i1 = 1!(info%PCells(1,1) - 1)*info%Y%nt*info%Y%np + (info%PCells(2,1) - 1)*info%Y%np + 1
    i2 = tpts!(info%PCells(1,2) - 1)*info%Y%nt*info%Y%np + (info%PCells(2,2)    )*info%Y%np

!   Often, we want an output matrix of size (3,3,nt,np,ncell) as output.
!   Exception checking for this option
    IF(PRESENT(u3)) u3 = 0D0
    IF(PRESENT(um)) THEN
        um = 0D0
        sz = SHAPE(um)
        nt = sz(3)
        np = sz(4)
        nc = sz(5)
        ic = 1!info%PCells(1,1) ! Again, previous, less efficient way to parallelize
        it = 1!info%PCells(2,1)
        ip = 0
    ENDIF
    IF(PRESENT(u3) .and. PRESENT(um)) THEN
        print *, 'ERROR: For integration on periodic fourier matrix,',  &
                 'choose rank 3 or rank 5 matrices, not both'
        stop
    ENDIF

!   This is the "interpolate back" step, where we have and send it back to discrete points
    DO i = i1, i2
!       Lots of managing of indices if we want matrix-type output
        IF(PRESENT(um)) THEN
            ip = ip + 1
            IF(ip .gt. np) THEN
                ip = 1
                it = it + 1
            ENDIF
            IF(it .gt. nt) THEN
                it = 1
                ic = ic + 1
            ENDIF
        ENDIF
        
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

            IF(PRESENT(u3)) &
            u3(:,:,i) = u3(:,:,i) + Ht(:, :, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv

            IF(PRESENT(um)) &
            um(:,:,it,ip,ic) = um(:,:,it,ip,ic) &
                             + Ht(:, :, iper, jper, kper)*info%par*EXP(r2*info%parexp)*h*h*h*Jbv

        ENDDO
    ENDDO
END SUBROUTINE EwaldintTold

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
    REAL(KIND = 8) :: gp, bv_gr(3,3), rr(3), r, w, P, m
    INTEGER :: i, j, k, it
    INTEGER, ALLOCATABLE :: tmpmat(:,:)

    gp = info%gp
    bv_gr = info%bv/REAL(gp)
    it = 0
!   Number of points over which we want our Gaussian to have support in 1D
    P = 9d0
!   Number of standard deviations this should correspond to (m=6, P=10) is good,
!   Based on Lindbo AK Tornberg Paper
    m = 6D0
!   Cutoff distance/Gaussian width
    w = MAXVAL(bv_gr)*P/2D0

!   With all this, we can calculate the new value of eta (smearing parameter)
    info%eta = (2D0*w*info%xi/m)**2
    info%par = (2D0*info%xi**2D0/pi/info%eta)**1.5D0
    info%parexp = -2D0*info%xi*info%xi/info%eta

    ALLOCATE(tmpmat(3,FLOOR(gp*gp*gp)))
    IF(ALLOCATED(info%suppMat)) DEALLOCATE(info%suppMat)

!   Find points that fall within cutoff distance
    DO i = -FLOOR(gp/2),FLOOR(gp/2)
        DO j = -FLOOR(gp/2),FLOOR(gp/2)
            DO k = -FLOOR(gp/2),FLOOR(gp/2)
                rr = bv_gr(:,1)*i + bv_gr(:,2)*j + bv_gr(:,3)*k
                r = NORM2(rr)
                IF(r .lt. w) THEN
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
! This whole thing is extremely dated
SUBROUTINE PeriodicTest(info)
    TYPE(sharedType), INTENT(INOUT), POINTER :: info
    TYPE(YType), POINTER :: Y
    COMPLEX(KIND = 8) :: f(3,2)
    REAL(KIND = 8) ::x0(3,2), Ht(3,info%gp,info%gp,info%gp)
    REAL(KIND = 8), ALLOCATABLE :: x(:,:,:,:), xv(:,:), Jtmp(:,:), x2(:,:), u2(:,:), u3(:,:,:), u4(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fs(:,:,:,:), fv1(:), fv2(:), fv3(:), ft(:,:)
    INTEGER :: i, j, ic, nt, np, tic, toc, rate, pts1, pts

    print *, "Starting tests for periodic functions..."
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
    ! CALL HtCalc(Ht, info, x0, f)
    ! CALL HtCalc(HtT, info, x0, f)
    ! CALL Ewaldint (Ht, info, x2, u2)
    ! print *, u2
    ! print *, ' '
    ! CALL Ewaldint(HtT, info, x2, u3)
    ! print *, u3
    ! print *, ' '
    ! print *, HtT(:,:,16,16,16)
    ! print *, ' '

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

    fs(2,:,:,:) = -fs(1,:,:,:)

    x(2,1,:,:) = x(1,1,:,:)
    x(2,2,:,:) = x(1,2,:,:) + .25D0
    x(2,3,:,:) = x(1,2,:,:)

!   Now the actual process of running the calculations, starting with putting matrices into vectors
    CALL SYSTEM_CLOCK(tic, rate)
    DO ic = 1,2
        ! CALL EwaldPreCalc(fin=fs(ic,1,:,:), f=fv1, Y=Y, x0in=x(ic,:,:,:), x0=xv, Ja=Jtmp)
        ! CALL EwaldPreCalc(fin=fs(ic,2,:,:), f=fv2, Y=Y, Ja=Jtmp)
        ! CALL EwaldPreCalc(fin=fs(ic,3,:,:), f=fv3, Y=Y, Ja=Jtmp)
    ENDDO
    ft(1,:) = fv1
    ft(2,:) = fv2
    ft(3,:) = fv3
    CALL SYSTEM_CLOCK(toc)
    print *, REAL(toc-tic)/REAL(rate)
    
    CALL SYSTEM_CLOCK(tic, rate)
    ! CALL EwaldG(Ht,info, xv, ft)
    CALL SYSTEM_CLOCK(toc)
    print *, REAL(toc-tic)/REAL(rate)

    CALL SYSTEM_CLOCK(tic, rate)
    i = SIZE(xv)
    ALLOCATE(u4(3,i/3))
    ! CALL Ewaldint (Ht, info, xv, u4)
    CALL SYSTEM_CLOCK(toc)
    print *, REAL(toc-tic)/REAL(rate) 
    ! print *, u2

!   Try it with evenly spaced points
    DEALLOCATE(xv, ft)
    pts1 = 40
    pts  = pts1*pts1
    ALLOCATE(xv(3,pts), ft(3,pts))
    ft = 0D0
    ic = 0
    DO i = 1,pts1
        DO j = 1,pts1
            ic = ic+1
            xv(1,ic) = (i-1)/REAL(pts1)
            xv(2,ic) = (j-1)/REAL(pts1)
            ft(1,ic) = CMPLX(MODULO(ic,2)*2-1)
        ENDDO
    ENDDO
    xv(3,:) = 0.5D0

    CALL SYSTEM_CLOCK(tic, rate)
    ! CALL HtCalc(Ht,info, xv, ft)
    CALL SYSTEM_CLOCK(toc)
    ! print *, REAL(toc-tic)/REAL(rate)
    ! CALL Ewaldint (Ht, info, x2, u2)
    ! print *, u2


END SUBROUTINE PeriodicTest

END MODULE PERIODICMOD