MODULE HARMMOD
USE UTILMOD
IMPLICIT NONE

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! The values at the individual nm's and derivatives
TYPE nmType
    INTEGER n, dero
    COMPLEX(KIND = 8), ALLOCATABLE :: v(:,:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: dY(:,:,:,:)

END TYPE nmType

! Full harmonic object evaluated at all values tht and phi
TYPE YType
!   Max order
    INTEGER :: p
!   Parameter space variables (e.g. angles, etc.)
    INTEGER :: np, nt
    REAL(KIND = 8), ALLOCATABLE :: phi(:), tht(:), th(:,:), ph(:,:)

!   Useful, extra items to carry around
    REAL(KIND = 8), ALLOCATABLE :: legs(:,:,:), facs(:), wg(:), ws(:)
    REAL(KIND = 8) :: dphi

!   Individual harmonics
    TYPE(nmType), ALLOCATABLE :: nm(:)

!   Items to keep for the FFT
    REAL(KIND = 8), ALLOCATABLE :: WSAVE(:)
    INTEGER LENSAV 
    

    CONTAINS
    PROCEDURE :: forward  => forwardY
    PROCEDURE :: backward => backwardY
    PROCEDURE :: rotate   => rotateY

END TYPE YType

INTERFACE YType
    PROCEDURE :: newY, newYbare
END INTERFACE YType

INTERFACE nmType
    PROCEDURE :: newnm
END INTERFACE nmType

CONTAINS

!=============================================================================!
!================================= ROUTIUNES =================================!
!=============================================================================!

! Construction routine for whole spheerical harmonic structure
FUNCTION newY(p, dero) RESULT(Y)
    TYPE(YType) Y
    INTEGER, INTENT(IN) :: p, dero

    INTEGER i, j
    REAL(KIND = 8) xs(p + 1)

!   Known constants based off p
    Y%p = p
    Y%nt = p + 1
    Y%np = 2*(p + 1)
    Y%dphi = 2*pi/Y%np

!   Allocate things
    ALLOCATE(Y%facs(2*p + 1), Y%tht(Y%nt), Y%phi(Y%np), &
             Y%ph(Y%nt,Y%np), Y%th(Y%nt, Y%np), Y%wg(Y%nt), Y%ws(Y%nt))

!   Calculate phis and thetas, with meshgrid as well
    CALL lgwt(Y%nt, xs, Y%wg)
    Y%tht = ACOS(xs)

    Y%phi(1) = 0
    DO i = 1,Y%np
        IF(i .gt. 1) Y%phi(i) = Y%phi(i-1) + Y%dphi
        Y%ph(:,i) = Y%phi(i)
        Y%th(:,i) = Y%tht
    ENDDO

!   Factorials, legendres, weights
    Y%facs(1) = 1
    DO i = 1,2*p
        Y%facs(i+1) = Y%facs(i)*i
    ENDDO

!   Since theta is constant across varying phi, only need 1 column of legendres
    ALLOCATE(Y%legs(Y%nt, 1, (p+1)*(p+2)/2))

    Y%legs(:,1,:) = Alegendre(p, xs)

    Y%ws = 0
    DO i = 1,Y%nt
        DO j = 0,p
            Y%ws(i) = Y%ws(i) + legendre(j, xs(i))/cos(Y%tht(i)/2)
        ENDDO
        Y%ws(i) = Y%ws(i)*Y%wg(i)
    ENDDO

!   Get the actual values of the harmonics at the phi and theta locations
    ALLOCATE(Y%nm(p + 1))
    DO i = 0, p
        Y%nm(i+1) = nmType(i, dero, Y)
    ENDDO

!   Initialize for FFT (np values in each transform)
    Y%LENSAV = 2*Y%np + INT(LOG(REAL(Y%np))) + 4
    ALLOCATE(Y%WSAVE(Y%LENSAV))
    CALL zfft1i(Y%np, Y%WSAVE, Y%LENSAV, i)

END FUNCTION newY

!------------------------------------------------------------------!
! Makes a new Y, but takes in thts and phis instead of assuming a 
! plain ole grid
FUNCTION newYbare(th, ph, p) RESULT(Y)
    TYPE(YType) Y

    REAL(KIND = 8) :: th(:,:), ph(:,:)
    INTEGER, OPTIONAL :: p
    INTEGER i, ntp(2)

    IF(present(p)) THEN
        Y%p = p
    ELSE
        Y%p  = INT(sqrt(REAL(size(th)/2)) - 1)
    ENDIF
    ntp  = SHAPE(th)
    Y%nt = ntp(1)
    Y%np = ntp(2)
    ALLOCATE(Y%th(Y%nt,Y%np), Y%ph(Y%nt,Y%np))
    Y%ph = ph
    Y%th = th
    
!   Here we need legendres everywhere, so allcoate based on that !!!!!!!!!! This is prrobably very inefficient
    ALLOCATE(Y%legs(Y%nt,Y%np,(Y%p+1)*(Y%p+2)/2))
    DO i = 1,Y%np
        Y%legs(:,i,:) = Alegendre(Y%p, cos(th(:,i)))
    ENDDO

    ALLOCATE(Y%nm(Y%p + 1))
    DO i = 0, Y%p
        Y%nm(i+1) = nmType(i, 0, Y)
    ENDDO

END FUNCTION newYbare

!------------------------------------------------------------------!
! Construction routine for single order spherical harmonic
FUNCTION newnm(n, dero, Y) RESULT(nm)
    TYPE(nmType) nm
    INTEGER, INTENT(IN) :: n, dero
    TYPE(Ytype), INTENT(IN) :: Y

    INTEGER m, im, il, k, l, sz(3)
    LOGICAL :: unif

    unif = .false.
!   Get some info on n
    nm%n = n
    ALLOCATE(nm%v(2*n + 1, Y%nt, Y%np), nm%dY(dero, 2*n + 1, Y%nt, Y%np))

!   Check if uniform grid    
    sz  = shape(Y%legs)
    IF(sz(2) .eq. 1) unif = .true.

!   Zero-th order
    IF(n .eq. 0) THEN
        nm%v = 0.5D0*ispi
        RETURN
    ENDIF

!   Cycle through individual degrees in a given order
    DO m = n, -n, -1
        im = m + n + 1
!       Negatve order has this special relationship to positive
        IF(m .lt. 0) THEN
            nm%v(im,:,:) = (-1D0)**m*CONJG(nm%v(im - 2*m,:,:))
!       Just the definition of spherical harmonics
        ELSE
            il = n*(n+1)/2 + m + 1
!           Make things cheaper if grid is uniform
            IF(unif) THEN
                DO k = 1,Y%nt
                    nm%v(im,k,:) = Y%legs(k, 1,il)*EXP(ii*m*Y%ph(k,:))
                ENDDO
            ELSE
                DO k = 1,Y%nt
                    DO l = 1,Y%np
                        nm%v(im,k,l) = Y%legs(k, l,il)*EXP(ii*m*Y%ph(k,l))
                    ENDDO
                ENDDO
            ENDIF
        ENDIF
    ENDDO

    IF(dero .gt. 0) THEN
        im = 0
        DO m = n, -n, -1
            im = m + n + 1
!           Calculate derivatives
            DO k = 1,dero
                nm%dY(k,im,:,:) = ThetDer(Y, nm%v, n, m, k)
            ENDDO
        ENDDO
    ENDIF
END FUNCTION newnm

!------------------------------------------------------------------!
! Takes the ord-order derivative at the points with respect to theta
RECURSIVE FUNCTION ThetDer(Y, v, n, m, ord) RESULT(dY)
    TYPE(Ytype), INTENT(IN) :: Y
    INTEGER, INTENT(IN) :: n, m, ord
    COMPLEX(KIND = 8), ALLOCATABLE :: dY(:,:)
    COMPLEX(KIND = 8), INTENT(IN) :: v(:,:,:)

    REAL(KIND = 8) :: nr, mr

    nr = REAL(n)
    mr = REAL(m)
    
    ALLOCATE(dY(Y%nt, Y%np))
    dY = 0D0
!   Recursive: if ord == 1, we're at the bottom. Do last derivative and exit loop
    IF(ord .eq. 1) THEN
!       Two parts added together
        IF(m .gt. -n) THEN
            dY = dY - 0.5D0*sqrt((nr + mr)*(nr - mr + 1D0))&
               * EXP( ii*Y%ph)*v(m + n,:,:)
        ENDIF
        IF(m .lt. n) THEN
            dY = dY + 0.5D0*sqrt((nr - mr)*(nr + mr + 1D0))&
               * EXP(-ii*Y%ph)*v(m + n + 2,:,:)
        ENDIF

!   Higher order derivatives, recursive
    ELSE
!       Two parts added together     
        IF(m .gt. -n) THEN
            dY = dY - 0.5D0*sqrt((nr + mr)*(nr - mr + 1D0))&
               * EXP( ii*Y%ph)*ThetDer(Y, v, n, m-1, ord-1)
        ENDIF
        IF(m .lt. n) THEN
            dY = dY + 0.5D0*sqrt((nr - mr)*(nr + mr + 1D0))&
               * EXP(-ii*Y%ph)*ThetDer(Y, v, n, m+1, ord-1)
        ENDIF
    ENDIF
END FUNCTION ThetDer

!------------------------------------------------------------------!
! Takes a collection of points in physical space and gets the spectral domain
FUNCTION forwardY(Y, f, p) RESULT(fmn)
    REAL(KIND = 8), INTENT(IN) :: f(:,:)
    CLASS(Ytype), TARGET, INTENT(IN) :: Y
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: p

    INTEGER :: pp, it, in, n, m, iserr
    COMPLEX(KIND = 8), ALLOCATABLE :: gm(:,:), gint(:), c(:)
    REAL(KIND = 8), POINTER :: legint(:), wg(:), dphi
    REAL(KIND = 8), ALLOCATABLE :: wrk(:)

    ALLOCATE(gm(Y%nt,Y%np), fmn((Y%p + 1)*(Y%p + 1)), c(Y%np), wrk(2*Y%np))
    wg => Y%wg
    dphi => Y%dphi
!   You can choose the order, but if you don't it'll just do max order    
    IF(PRESENT(p)) THEN
        pp = p
    ELSE
        pp = Y%p
    ENDIF

!   Do FFT along phi for each value of theta
    DO it = 1, Y%nt
        c = f(it,:)
        CALL ZFFT1F(Y%np, 1, c, Y%np, Y%WSAVE, Y%LENSAV, wrk, 2*Y%np, iserr)
        gm(it,:) = c*dphi*Y%np
    ENDDO

!   Now do the integrals across the thetas
    it = 0
    in = 0

    DO n = 0,pp
        in = in + n
        DO m = -n,n
            it = it + 1
!           Get the right Legendres and Fourier constants to sum over,
!           exploiting symmetry of problem
            IF(m .ge. 0)THEN
                legint => Y%legs(:, 1, in + m + 1)
                gint = gm(:, 1 + m)
            ELSE
                legint => Y%legs(:, 1, in - m + 1)
                gint = CONJG(gm(:, 1 - m))*(-1D0)**(m)
            ENDIF

!           Integrate
            fmn(it) = SUM(wg*legint*gint)
        ENDDO
    ENDDO
END FUNCTION forwardY

!------------------------------------------------------------------!
! Takes a collection of points in spectral space and gets the physical domain
FUNCTION backwardY(Y, fmn, p) RESULT(f)
    COMPLEX(KIND = 8), INTENT(IN) :: fmn(:)
    CLASS(Ytype), TARGET, INTENT(IN) :: Y
    REAL(KIND = 8), ALLOCATABLE :: f(:,:)
    INTEGER, OPTIONAL, INTENT(IN) :: p

    INTEGER pp, ih, it, n, m, im
    TYPE(nmType), POINTER :: nm

!   You can choose the order, but if you don't it'll just do max order
    IF(PRESENT(p)) THEN
        pp = p
    ELSE
        pp = Y%p
    ENDIF

!   Performed at all theta and phi in Y
    ALLOCATE(f(Y%nt,Y%np))
    f = 0
    ih = 0
    it = 0

!   Just loop through and perform the sums
    DO n = 0,pp
        ih = ih + 1
        nm => Y%nm(ih)
        im = 0
!       If f is real-valued we really only need to go m = 0,n
        DO m = -n,n
            it = it + 1
            im = im + 1
!           This assumes f is a real-valued function,
!           and that the imaginary values will have canceled anyway
            f = f + REAL(fmn(it)*nm%v(im,:,:))
        ENDDO
    ENDDO
END FUNCTION backwardY


! -------------------------------------------------------------------------!
! Rotate spherical harmonic constants by given Euler angles
FUNCTION rotateY(Y, fmn, a, b, c) RESULT(f)
    REAL(KIND = 8) a, b, c
    COMPLEX(KIND = 8) fmn(:)
    CLASS(YType), TARGET :: Y
    COMPLEX(KIND = 8) f(SIZE(fmn))

    INTEGER n, mp, it, m, im, s, i
    INTEGER, ALLOCATABLE :: frng(:)
    REAL(KIND = 8) Smm, dmms
    COMPLEX(KIND = 8), ALLOCATABLE :: Dmm(:)
    REAL(KIND = 8), POINTER :: facs(:)

    facs => Y%facs

    f = 0
    it = 0

    ALLOCATE(frng(1))
    ALLOCATE(Dmm(1))

!   Loop over harmonic order
    DO n = 0, Y%p
        DEALLOCATE(frng)
        ALLOCATE(frng(2*n + 1))
        frng = (/(i, i=(it+1),(it+2*n+1), 1)/) 
        
!       Loop over harmonic degree we're trying to calculate
        DO mp = -n,n
            DEALLOCATE(Dmm)
            ALLOCATE(Dmm(2*n + 1))
            im = 0
            it = it+1
!           Loop over harmonic degree we're using to calculate
            DO m = -n,n
                Smm = 0D0
                im = im+1
                DO s = MAX(0, m-mp), MIN(n+m, n-mp)
                    Smm = Smm + (-1D0)**s*(COS(b/2D0)**(2D0*(n - s) + m - mp) &
                        * SIN(b/2D0)**(2D0*s - m + mp)) & 
                        / (facs(n+m-s+1)*facs(s+1)*facs(mp-m+s+1)*facs(n-mp-s+1))
                ENDDO
                dmms = (-1D0)**(mp-m)*(facs(n+mp+1)*facs(n-mp+1)*facs(n+m+1)*facs(n-m+1))**0.5D0*Smm
                Dmm(im) = EXP(ii*m*a)*dmms
            ENDDO
            f(it) = f(it) + SUM(fmn(frng)*Dmm)
        ENDDO
    ENDDO
    
!   Rotate back in last angle    
    it = 0
    DO n = 0,Y%p
        DO m = -n,n
            it = it+1
            f(it) = f(it)*EXP(ii*m*c)
        ENDDO
    ENDDO
END FUNCTION rotateY

! -------------------------------------------------------------------------!
! Calculate spherical harmonic coefficients for an RBC
FUNCTION RBCcoeff(Y, ord) RESULT(xmn)
    TYPE(YType),INTENT(IN) :: Y
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:)
    INTEGER, OPTIONAL :: ord

    REAL(KIND = 8), ALLOCATABLE :: x1(:,:), x2(:,:), x3(:,:), xy(:), z(:)
    REAL(KIND = 8) :: alph, a
    INTEGER :: p, i

    ALLOCATE(x1(Y%nt,Y%np), x2(Y%nt,Y%np), x3(Y%nt,Y%np), &
    xy(Y%nt), z(Y%nt), xmn(3,(Y%p + 1)*(Y%p + 1)))

    alph = 1.38581894D0
    a = 1D0

    IF(PRESENT(ord)) THEN
        p = ord
    ELSE
        p = Y%p
    ENDIF

!   x/y coord in a constant phi plane, to be rotated around phi
    xy = a*alph*SIN(Y%tht)
    z = 0.5D0*a*alph*(0.207D0 + 2.003D0*SIN(Y%tht)*SIN(Y%tht) &
      - 1.123D0*SIN(Y%tht)**4D0)*COS(Y%tht)

    DO i = 1,Y%np
        x1(:,i) = xy*COS(Y%phi(i))
        x2(:,i) = xy*SIN(Y%phi(i))
        x3(:,i) = z
    ENDDO

    xmn(1,:) = Y%forward(x1)
    xmn(2,:) = Y%forward(x2)
    xmn(3,:) = Y%forward(x3)

!   Gotta get that max precision
    DO i = 1,INT(size(xmn)/3)
        if(abs(xmn(1,i)).lt.1E-12) xmn(1,i) = 0D0
        if(abs(xmn(2,i)).lt.1E-12) xmn(2,i) = 0D0
        if(abs(xmn(3,i)).lt.1E-12) xmn(3,i) = 0D0
    ENDDO
END FUNCTION RBCcoeff

! -------------------------------------------------------------------------!
! Calculate spherical harmonic coefficients for a Sphere
FUNCTION Spherecoeff(Y, ord) RESULT(xmn)
    TYPE(YType),INTENT(IN) :: Y
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:)
    INTEGER, OPTIONAL :: ord

    REAL(KIND = 8), ALLOCATABLE :: x1(:,:), x2(:,:), x3(:,:)
    INTEGER :: p, i

    ALLOCATE(x1(Y%nt,Y%np), x2(Y%nt,Y%np), x3(Y%nt,Y%np), &
        xmn(3,(Y%p + 1)*(Y%p + 1)))


    IF(PRESENT(ord)) THEN
        p = ord
    ELSE
        p = Y%p
    ENDIF

    x1 = 0.5D0*COS(Y%ph)*SIN(Y%th)
    x2 = 0.5D0*SIN(Y%ph)*SIN(Y%th)
    x3 = COS(Y%th)

    ! x1(1:10,1:10) = x1(1:10,1:10) + .05D0
    ! x2(1:10,1:10) = x2(1:10,1:10) + .08D0
    ! x3(1:3,1:3) = x3(1:3,1:3) + .05D0

    xmn(1,:) = Y%forward(x1)
    xmn(2,:) = Y%forward(x2)
    xmn(3,:) = Y%forward(x3)

    xmn(:,1) = 0D0

!   Gotta get that max precision
    DO i = 1,INT(size(xmn)/3)
        if(abs(xmn(1,i)).lt.1E-12) xmn(1,i) = 0D0
        if(abs(xmn(2,i)).lt.1E-12) xmn(2,i) = 0D0
        if(abs(xmn(3,i)).lt.1E-12) xmn(3,i) = 0D0
    ENDDO
END FUNCTION Spherecoeff

END MODULE HARMMOD