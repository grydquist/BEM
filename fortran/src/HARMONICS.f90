MODULE HARMMOD
USE UTILMOD
USE CMMOD
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

! Object that saves info for rotation
TYPE rotType
    INTEGER :: n
    COMPLEX(KIND = 8), ALLOCATABLE :: Dmm(:,:)
    REAL(KIND = 8), ALLOCATABLE :: dmms(:,:)

END TYPE rotType

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

!   Save rotation info
    LOGICAL :: hasrot = .false.
    TYPE(rotType), ALLOCATABLE :: rot(:,:,:)

    CONTAINS
    PROCEDURE :: forward  => forwardYDirect! forwardY! 
    PROCEDURE :: backward => backwardY
    PROCEDURE :: rotate   => rotateY
    PROCEDURE :: rotcnst   => rotcnstY

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
FUNCTION newY(p, dero, rot) RESULT(Y)
    TYPE(YType) Y
    INTEGER, INTENT(IN) :: p, dero

    INTEGER i, j, n
    REAL(KIND = 8) xs(p + 1)

!   Do we want to store rotation info about each point?
    LOGICAL, INTENT(IN) :: rot

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
            Y%ws(i) = Y%ws(i) + legendre(j, xs(i))*sin(Y%tht(i)/2D0)*2D0
        ENDDO
        Y%ws(i) = Y%ws(i)*Y%wg(i)
    ENDDO

!   Get the actual values of the harmonics at the phi and theta locations
    ALLOCATE(Y%nm(p + 1))
    DO i = 0, p
        Y%nm(i+1) = nmType(i, dero, Y)
    ENDDO

!   Initialize for FFT (np values in each transform)
    Y%LENSAV =  2 * Y%np + int ( log ( real ( Y%np, kind = 8 ) ) ) + 4 !2*Y%np + INT(LOG(REAL(Y%np))) + 4
    ALLOCATE(Y%WSAVE(Y%LENSAV))
    CALL zfft1i(Y%np, Y%WSAVE, Y%LENSAV, i)

!   Rotation info
    Y%hasrot = rot
    IF(rot) THEN
        ALLOCATE(Y%rot(Y%nt,Y%np, Y%p+1))
!       Go to each point and get all the rotation constants at that point
        DO i = 1,Y%nt
            DO j = 1,Y%np
                DO n = 0,Y%p
                    Y%rot(i,j,n + 1) = Y%rotcnst(Y%phi(j), -Y%tht(i), n)
                ENDDO
            ENDDO
        ENDDO
    ENDIF

END FUNCTION newY

!------------------------------------------------------------------!
! Makes a new Y, but takes in thts and phis instead of assuming a 
! plain ole grid
FUNCTION newYbare(th, ph, p, dero) RESULT(Y)
    TYPE(YType) Y

    REAL(KIND = 8) :: th(:,:), ph(:,:)
    INTEGER, OPTIONAL :: p, dero
    INTEGER i, ntp(2)

    IF(PRESENT(p)) THEN
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

!   Get locations and derivatives if present
    ALLOCATE(Y%nm(Y%p + 1))
    IF(PRESENT(dero)) THEN
        DO i = 0, Y%p
            Y%nm(i+1) = nmType(i, dero, Y)
        ENDDO
    ELSE
        DO i = 0, Y%p
            Y%nm(i+1) = nmType(i, 0, Y)
        ENDDO
    ENDIF

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
        nm%dY= 0D0
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
                        nm%v(im,k,l) = Y%legs(k, l,il)*EXP(ii*m*Y%ph(k,l)) !! Most expensive line
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
        IF(iserr.ne.0) THEN
            print *, 'FFT forward error, code: ', iserr
            STOP
        ENDIF

!       Sometimes FFTpack just fails. Check here        
        CALL ZFFT1B(Y%np, 1, c, Y%np, Y%WSAVE, Y%LENSAV, wrk, 2*Y%np, iserr)
        IF(iserr.ne.0) THEN
            print *, 'FFT backward error, code: ', iserr
            STOP
        ENDIF

        IF(MAXVAL(real(c) - f(it,:)) .gt. 1e-10) THEN
            print *, 'FFT failure! Max difference between values:'
            print *, MAXVAL(real(c) - f(it,:))
            STOP
        endif
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
! Takes a collection of points in physical space and gets the spectral domain
! by directly integrating instead of with FFTs
FUNCTION forwardYDirect(Y, f, p) RESULT(fmn)
    REAL(KIND = 8), INTENT(IN) :: f(:,:)
    CLASS(Ytype), TARGET, INTENT(IN) :: Y
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: p
    INTEGER :: i, n, m, im, it, j
    TYPE(nmType), POINTER :: nm

    ALLOCATE(fmn((Y%p + 1)*(Y%p + 1)))
    fmn = 0D0

    IF(PRESENT(p)) THEN
        n = p
    ELSE
        n = Y%p
    ENDIF

    im = 0
    it = 0
!   Just loop through and perform the sums
    DO i = 0,n
        nm => Y%nm(i + 1)
        im = 0
        DO m = -i,i
            it = it + 1
            im = im + 1
            DO j = 1,Y%np
                fmn(it) = fmn(it) + SUM(CONJG(nm%v(im,:,j))*f(:,j)*Y%wg*Y%dphi)
            ENDDO
        ENDDO
    ENDDO

END FUNCTION forwardYDirect
!------------------------------------------------------------------!
! Takes a collection of points in spectral space and gets the physical domain
FUNCTION backwardY(Y, fmn, p) RESULT(f)
    COMPLEX(KIND = 8), INTENT(IN) :: fmn(:)
    CLASS(Ytype), TARGET, INTENT(IN) :: Y
    REAL(KIND = 8), ALLOCATABLE :: f(:,:)
    INTEGER, OPTIONAL, INTENT(IN) :: p

    INTEGER pp, it, n, m, im
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
    it = 0

!   Just loop through and perform the sums
    DO n = 0,pp
        nm => Y%nm(n+1)
        im = n
        it = it + n
!       Exploit symmetry properties
        DO m = 0,n
            it = it + 1
            im = im + 1
!           This assumes f is a real-valued function,
!           and that the imaginary values will have canceled anyway
            IF(m .ne. 0) THEN
                f = f + REAL(2D0*fmn(it)*nm%v(im,:,:))
            ELSE
                f = f + REAL(fmn(it)*nm%v(im,:,:))
            ENDIF
        ENDDO
    ENDDO
END FUNCTION backwardY

! -------------------------------------------------------------------------!
! Get constants for given angles and order
FUNCTION rotcnstY(Y, a, b, n) RESULT(rot)
    REAL(KIND = 8) a, b
    CLASS(YType), TARGET :: Y
    INTEGER m, mp, n, im, s, im2
    TYPE(rotType) :: rot
    REAL(KIND = 8) Smm, dmms
    REAL(KIND = 8), POINTER :: facs(:)

    facs => Y%facs
    rot%n = n
    ALLOCATE(rot%Dmm(2*n + 1, n + 1))
    ALLOCATE(rot%dmms(2*n + 1, 2*n + 1))

!   Loop over harmonic degree we're trying to calculate
    DO mp = 0,n
        im = 0
!       Loop over harmonic degree we're using to calculate
        DO m = -n,n
            Smm = 0D0
            im = im+1
            DO s = MAX(0, m-mp), MIN(n+m, n-mp)
                Smm = Smm + (-1D0)**s*(COS(b/2D0)**(2D0*(n - s) + m - mp) &
                    * SIN(b/2D0)**(2D0*s - m + mp)) & 
                    / (facs(n+m-s+1)*facs(s+1)*facs(mp-m+s+1)*facs(n-mp-s+1))
            ENDDO
            dmms = (-1D0)**(mp-m)*(facs(n+mp+1)*facs(n-mp+1)*facs(n+m+1)*facs(n-m+1))**0.5D0*Smm
            rot%Dmm(im,mp+1) = EXP(ii*m*a)*dmms
        ENDDO
    ENDDO

!!  This is bad, but we need all of the dmms constant
!   Loop over harmonic degree we're trying to calculate
    im2 = 0
    DO mp = -n,n
        im2 = im2+1
        im = 0
!       Loop over harmonic degree we're using to calculate
        DO m = -n,n
            Smm = 0D0
            im = im+1
            DO s = MAX(0, m-mp), MIN(n+m, n-mp)
                Smm = Smm + (-1D0)**s*(COS(b/2D0)**(2D0*(n - s) + m - mp) &
                    * SIN(b/2D0)**(2D0*s - m + mp)) & 
                    / (facs(n+m-s+1)*facs(s+1)*facs(mp-m+s+1)*facs(n-mp-s+1))
            ENDDO
            dmms = (-1D0)**(mp-m)*(facs(n+mp+1)*facs(n-mp+1)*facs(n+m+1)*facs(n-m+1))**0.5D0*Smm
            rot%dmms(im,im2) = dmms
        ENDDO
    ENDDO

END FUNCTION rotcnstY

! -------------------------------------------------------------------------!
! Rotate spherical harmonic constants by given Euler angles
FUNCTION rotateY(Y, fmn, ti, pj, c) RESULT(f)
    REAL(KIND = 8) c
    COMPLEX(KIND = 8) fmn(:)
    CLASS(YType), TARGET :: Y
    COMPLEX(KIND = 8) f(SIZE(fmn))

    INTEGER n, mp, it, m,  i, ti, pj
    INTEGER, ALLOCATABLE :: frng(:)
    COMPLEX(KIND = 8), POINTER :: Dmm(:,:)
    REAL(KIND = 8), POINTER :: facs(:)

!   Need the rotation constants first
    IF(.not. Y%hasrot) THEN
        print *, 'ERROR: Need to calc rotation cnsts before rotation'
    ENDIF

    facs => Y%facs
    f = 0
    it = 0
    ALLOCATE(frng(1))

!   Loop over harmonic order
    DO n = 0, Y%p
        DEALLOCATE(frng)
        ALLOCATE(frng(2*n + 1))
!       Indices of fmn at a given n
        frng = (/(i, i=(it+1),(it+2*n+1), 1)/) 
        it = it+n
!       Loop over harmonic degree we're trying to calculate
        DO mp = 0,n
            it = it+1
            Dmm => Y%rot(ti, pj, n + 1)%Dmm
            f(it) = SUM(fmn(frng)*Dmm(:,mp + 1))
            IF(mp .ne. 0) THEN
                f(it - mp*2) = (-1D0)**mp*CONJG(f(it))
            ENDIF
        ENDDO
    ENDDO
    
!   Rotate back in last angle (if nonzero)
    IF(c .ne. 0) THEN
    it = 0
    DO n = 0,Y%p
        DO m = -n,n
            it = it+1
            f(it) = f(it)*EXP(ii*m*c)
        ENDDO
    ENDDO
    ENDIF
END FUNCTION rotateY

! -------------------------------------------------------------------------!
! Calculate spherical harmonic coefficients for an RBC (Poz 2003)
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
FUNCTION Spherecoeff(Y, obl, ord) RESULT(xmn)
    TYPE(YType),INTENT(IN) :: Y
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:)
    REAL(KIND = 8), OPTIONAL :: obl
    REAL(KIND = 8) :: oblt
    INTEGER, OPTIONAL :: ord

    REAL(KIND = 8), ALLOCATABLE :: x1(:,:), x2(:,:), x3(:,:)
    INTEGER :: p, i

    ALLOCATE(x1(Y%nt,Y%np), x2(Y%nt,Y%np), x3(Y%nt,Y%np), &
        xmn(3,(Y%p + 1)*(Y%p + 1)))

!   How oblate is the sphere (Default perfect sphere)
    IF(PRESENT(obl)) THEN
        oblt = obl
    ELSE
        oblt = 1D0
    ENDIF

    IF(PRESENT(ord)) THEN
        p = ord
    ELSE
        p = Y%p
    ENDIF

    ! x1 = 0.5D0*COS(Y%ph)*SIN(Y%th)
    ! x2 = 0.5D0*SIN(Y%ph)*SIN(Y%th)
    x1 = COS(Y%ph)*SIN(Y%th)
    x2 = SIN(Y%ph)*SIN(Y%th)
    x3 = oblt*COS(Y%th)!0.918D0

    xmn(1,:) = Y%forward(x1)
    xmn(2,:) = Y%forward(x2)
    xmn(3,:) = Y%forward(x3)

    xmn(:,1) = 0D0
    xmn(:,5:(Y%p + 1)*(Y%p + 1)) = 0D0

!   Gotta get that max precision
    DO i = 1,INT(size(xmn)/3)
        if(abs(xmn(1,i)).lt.1E-12) xmn(1,i) = 0D0
        if(abs(xmn(2,i)).lt.1E-12) xmn(2,i) = 0D0
        if(abs(xmn(3,i)).lt.1E-12) xmn(3,i) = 0D0
    ENDDO
END FUNCTION Spherecoeff

! -------------------------------------------------------------------------!
! Calculate spherical harmonic coefficients for a cubish thing (gives flattish surface p = 8)
FUNCTION Cubecoeff(Y, ord) RESULT(xmn)
    TYPE(YType),INTENT(IN) :: Y
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:)
    INTEGER, OPTIONAL :: ord

    REAL(KIND = 8), ALLOCATABLE :: x1(:,:), x2(:,:), x3(:,:)
    INTEGER :: p, i, j

    ALLOCATE(x1(Y%nt,Y%np), x2(Y%nt,Y%np), x3(Y%nt,Y%np), &
        xmn(3,(Y%p + 1)*(Y%p + 1)))


    IF(PRESENT(ord)) THEN
        p = ord
    ELSE
        p = Y%p
    ENDIF

!   If statements for the six faces

!   Cylinder for now
    DO i = 1,Y%nt
        DO j = 1,Y%np
!           First top/bottom faces
            IF(Y%tht(i) .lt. pi/4D0) THEN
                x1(i,j) = TAN(Y%tht(i))*COS(Y%phi(j))
                x2(i,j) = TAN(Y%tht(i))*SIN(Y%phi(j))
                x3(i,j) = 1D0
            ELSEIF (Y%tht(i) .gt.3D0*pi/4D0) THEN
                x1(i,j) = TAN(ABS(Y%tht(i) - pi))*COS(Y%phi(j))!-SQRT(2D0)*ATAN2(Y%tht(i) - pi,1D0)*COS(Y%phi(j))
                x2(i,j) = TAN(ABS(Y%tht(i) - pi))*SIN(Y%phi(j))!-SQRT(2D0)*ATAN2(Y%tht(i) - pi,1D0)*SIN(Y%phi(j))
                x3(i,j) = -1D0
            ELSE
                ! IF(Y%tht(i) .gt. pi/2D0) THEN
                !     thtmp = pi/2D0 - ABS(Y%tht(i) - pi)
                !     x3(i,j) = -TAN(thtmp)
                ! ELSE
                !     thtmp = pi/2D0 - Y%tht(i)
                !     x3(i,j) = TAN(thtmp)
                ! ENDIF
                ! x1(i,j) = COS(Y%phi(j))
                ! x2(i,j) = SIN(Y%phi(j))
                x1(i,j) = sqrt(2D0)*SIN(Y%tht(i))*COS(Y%phi(j))
                x2(i,j) = sqrt(2D0)*SIN(Y%tht(i))*SIN(Y%phi(j))
                x3(i,j) = sqrt(2D0)*COS(Y%tht(i))
            ENDIF
        ENDDO
    ENDDO

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
END FUNCTION Cubecoeff

! -------------------------------------------------------------------------!
! Takes a list of coefficients that go real(x1mn), imag(x1mn), real(x2mn)...
! and outputs them properly
FUNCTION Readcoeff(filen, ord) RESULT(xmn)
    COMPLEX(KIND = 8), ALLOCATABLE :: xmn(:,:)
    REAL(KIND = 8) :: xmnrawind
    REAL(KIND = 8), ALLOCATABLE :: xmnraw(:)
    CHARACTER (len=*) filen
    INTEGER, INTENT(IN) :: ord
    INTEGER stat, p, i, jmp

    ALLOCATE(xmn(3, (ord+1)*(ord+1)))
    p = 0

!   Read once to find file size. Not efficient but probably not a big deal
    OPEN(unit = 13, file = filen, action = 'read')
    DO 
        READ(13, *, iostat = stat) xmnrawind
        IF(stat .ne. 0) EXIT
        p = p+1
    ENDDO
    CLOSE(13)
    
!   Now we know the file size, so we can allocate properly
    ALLOCATE(xmnraw(p))

!   Text file format: all x real, all x imag, all y real, all y imag...
!   And loop through again to get the values into this big array
    OPEN(unit = 13, file = filen, action = 'read')
    DO i = 1,p
        READ(13, *, iostat = stat) xmnraw(i)
    ENDDO
    CLOSE(13)

    xmnrawind = p/6
    p = int(sqrt(xmnrawind)) - 1
    jmp = (p+1)*(p+1)
    
    IF(ord>p) THEN
        print *, 'ERROR: when reading coeffs, desired order higher than supplied by input' !!! I could just pad w/ zeros
        STOP
    ENDIF

!   Fill out complex coefficient matrix!
    xmn(1,:) = xmnraw(1:(ord+1)*(ord+1))
    xmn(1,:) = xmn(1,:) + xmnraw(jmp + 1 : jmp + (ord+1)*(ord+1))*ii 

    xmn(2,:) = xmnraw(2*jmp + 1 : 2*jmp + (ord+1)*(ord+1))
    xmn(2,:) = xmn(2,:) + xmnraw(3*jmp + 1 : 3*jmp + (ord+1)*(ord+1))*ii 

    xmn(3,:) = xmnraw(4*jmp + 1 : 4*jmp + (ord+1)*(ord+1))
    xmn(3,:) = xmn(3,:) + xmnraw(5*jmp + 1 : 5*jmp + (ord+1)*(ord+1))*ii
END FUNCTION Readcoeff

END MODULE HARMMOD