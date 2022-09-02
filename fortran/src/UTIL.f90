MODULE UTILMOD
IMPLICIT NONE
REAL(KIND = 8), PARAMETER :: ispi = 0.56418958354775627928D0, pi = 3.1415926535897932384626433D0
COMPLEX(KIND = 8), PARAMETER :: ii = (0D0,1D0)

!==============================================================================!
!                   This module has a bunch of random, extra,                  !
!                                useful functions                              !
!==============================================================================!

INTERFACE READ_MFS
    PROCEDURE READ_GRINT_CHAR, READ_GRINT_DOUB, READ_GRINT_INT
END INTERFACE

CONTAINS

! -------------------------------------------------------------------------!
! Calculates (all) the asscoiated legendre polynomials at
! location x up to order n using recursive thing from NR.
! No negative m

! IF YOU'RE HAVING PRECISION ISSUES, MAYBE TRY HERE.
! Try comparing a then backward transform of x1mn, and note
! that a couple of the last values are O(10^-8) when they
! should be zero. The Legendres are also off by this much.
! Could be precision of pi??? FFT seems fine
FUNCTION Alegendre(n, x) RESULT(P)
    INTEGER n, i, m, it
    REAL(KIND = 8) x(:), facto, facte
    REAL(KIND = 8) P(size(x), (n+1)*(n+2)/2)

!   0th order, easy
    P(:,1) = 0.5D0*ispi
    IF(n .eq. 0) RETURN

!   Both first orders
    P(:,2) =  ispi*sqrt(3D0*0.25D0)*x
    P(:,3) = -ispi*sqrt(3D0*0.125D0*(1D0-x*x))
    IF(n .eq. 1) RETURN

!   Factorials/indices to keep track of
    facto = 1D0
    facte = 2D0
    it = 3

!   Higher order polynomials, from 2 to n
    DO i = 2,n
        DO m = 0,i
            it = it + 1
!           Analytical form when m == n            
            IF(m .eq. i) THEN
                facto = facto*(2D0*m - 1D0)
                facte = facte*2D0*m*(2D0*m - 1D0)
                P(:,it) = (-1D0)**m * ispi*sqrt((2D0*m + 1D0)/(4D0*facte))*facto &
                        * (1 - x*x)**(0.5D0*m)

!           Recurrence relationship when m == i-1
            ELSEIF(m .eq. i-1) THEN
                P(:,it) = sqrt(2D0*m + 3D0)*x*P(:,it-i)

!           General recurrence relationship
            ELSE
                P(:,it) = sqrt((4D0*i*i - 1D0)/(i*i - m*m))*(x*P(:,it-i) &
                        - sqrt(((i-1D0)*(i-1D0) - m*m)/(4D0*(i-1D0)*(i-1D0) - 1D0))*P(:,it - 2*i + 1))
            ENDIF
        ENDDO
    ENDDO
END FUNCTION Alegendre

! -------------------------------------------------------------------------!
! Calculates Gauss points and weights on interval -1, 1
SUBROUTINE lgwt(N,x,w)
    INTEGER, INTENT(IN) :: N
    REAL(KIND = 8), INTENT(OUT) :: x(N), w(N)
    INTEGER :: i
    REAL(KIND = 8) :: xu(N), L(N,N+1), Lp(N), y0(N)
    
    xu = (/(i, i = 1, N, 1)/)*2D0/N - 1D0
    x = (/(i, i = 1, N, 1)/) - 1D0
    
    x = cos((2D0*x + 1D0)*pi/(2D0*N)) + (0.27D0/N)*sin(pi*xu*(N-1D0)/(N+1D0))
    y0 = 2D0

    L = 0D0

    DO WHILE(MAXVAL(ABS(x-y0)) .gt. 1e-15)
        L(:,1) = 1D0
        L(:,2) = x
        DO i  = 2,N
            L(:,i + 1) = ((2D0*i - 1D0)*x*L(:,i) - (i-1D0)*L(:,i-1))/i
        ENDDO
        Lp = (N+1D0)*(L(:,N) - x*L(:,N+1))/(1D0 - x*x)
        y0 = x
        x = y0 - L(:,N+1)/Lp
    ENDDO
    w = 2D0/((1D0 - x*x)*Lp*Lp)*((real(N)+1D0)/real(N))*((real(N)+1D0)/real(N))
END SUBROUTINE lgwt

! -------------------------------------------------------------------------!
! (Not associated) Legendre polynomials
FUNCTION legendre(n, x) RESULT(P)
    INTEGER n, i
    REAL(KIND = 8) :: x
    REAL(KIND = 8) :: P, ri
    REAL(KIND = 8) :: Pn1, Pn2

!   Start recursive relationship    
    Pn1 = x
    Pn2 = 1

    IF(n .eq. 0) THEN
        P = Pn2
        RETURN
    ENDIF
    
    P = Pn1

    IF(n .eq. 1) THEN
        RETURN
    ENDIF

    DO i = 2,n
        ri = REAL(i)
        P = (2*ri - 1)/ri*x*Pn1 - (ri-1)/ri*Pn2
        Pn2 = Pn1
        Pn1 = P
    ENDDO
END FUNCTION legendre

! -------------------------------------------------------------------------!
! This routine does the cross product for a two given vector of 
! V1 and V2.
PURE FUNCTION CROSS(V1, V2) RESULT(U)
    REAL(KIND=8), INTENT(IN) :: V1(:), V2(:)
    REAL(KIND=8) U(3)
    U(1) = V1(2)*V2(3) - V1(3)*V2(2)
    U(2) = V1(3)*V2(1) - V1(1)*V2(3)
    U(3) = V1(1)*V2(2) - V1(2)*V2(1)
    RETURN
END FUNCTION CROSS

! -------------------------------------------------------------------------!
! This routine does the dot product for a two given vector of 
! V1 and V2.
PURE FUNCTION DOT(V1, V2) RESULT(U)
    REAL(KIND=8), INTENT(IN) :: V1(:), V2(:)
    REAL(KIND=8) U
    U = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)
    RETURN
END FUNCTION DOT

! -------------------------------------------------------------------------!
! This routine does the outer product for a two given vector of 
! V1 and V2.
PURE FUNCTION OUTER(V1, V2) RESULT(U)
    REAL(KIND=8), INTENT(IN) :: V1(:), V2(:)
    REAL(KIND=8) U(3,3)
    U(1,1) = V1(1)*V2(1)
    U(2,1) = V1(2)*V2(1)
    U(3,1) = V1(3)*V2(1)
    U(1,2) = V1(1)*V2(2)
    U(2,2) = V1(2)*V2(2)
    U(3,2) = V1(3)*V2(2)
    U(1,3) = V1(1)*V2(3)
    U(2,3) = V1(2)*V2(3)
    U(3,3) = V1(3)*V2(3)
    RETURN
END FUNCTION OUTER

! -------------------------------------------------------------------------!
! Calculates determinant of 3x3 matrix
PURE FUNCTION DET3(A) RESULT(U)
    REAL(KIND = 8), INTENT(IN) :: A(3,3)
    REAL(KIND = 8) :: U
    U = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
      + A(1,2)*(A(2,3)*A(3,1) - A(2,1)*A(3,3)) &
      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
END FUNCTION DET3

! -------------------------------------------------------------------------!
! Calculates the eigenvals/vecs of a 3x3 symmetric matrix (from wikipedia)
SUBROUTINE EIG3(A, e, ev)
    REAL(KIND = 8), INTENT(IN) :: A(3,3)
    REAL(KIND = 8), INTENT(OUT) :: e(3)
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: ev(3,3)
    REAL(KIND = 8) p1, q, p2, p, B(3,3), r, I(3,3), phi, tmp(3,3), &
                   chkx, chky
    INTEGER :: j

    p1 = A(1,2)*A(1,2) + A(1,3)*A(1,3) + A(2,3)*A(2,3)

!   Identity    
    I = 0D0
    FORALL(j = 1:3) I(j,j) = 1D0
    q = (A(1,1) + A(2,2) + A(3,3))/3D0
    p2 = (A(1,1) - q)*(A(1,1) - q) + (A(2,2) - q)*(A(2,2) - q) &
        + (A(3,3) - q)*(A(3,3) - q) + 2D0*p1
    p = sqrt(p2/6D0)
    B = (1D0/p)*(A - q*I)
    r = DET3(B)*0.5D0

    IF(r .le. -1D0 + 1e-15) THEN
        phi = pi/3D0
    ELSEIF(r .ge. 1D0 - 1e-15) THEN
        phi = 0D0
    ELSE
        phi = ACOS(r)/3D0 !!!!!!!!!!!!!!!!!!!!!! A good bit of precision is lost on this acos
    ENDIF

    e(1) = q + 2D0*p*COS(phi)
    e(3) = q + 2D0*p*COS(phi + 2D0*pi/3D0)
    e(2) = 3D0*q - e(1) - e(3)

!   Eigenvectors if wanted    
    IF(PRESENT(ev)) THEN

!       If it's the first time step, it should be e(1) == e(2) == 1, causing issues.
!       We have an eigen plane instead, which is perpendicular to the normal, and we
!       can just use the surface basis vectors. Alternatively, we can get the normal
!       vector from the 0 eigenvalue, then another eigenvalue, and do the cross product
        IF(ABS(e(1) - e(2)) .lt. 1e-15) THEN
!           0 eigenvector
            tmp = A - e(3)*I
            tmp = RREF3(tmp)
!           Rare case x-value = 0
            IF(ABS(tmp(1,3)).lt.1e-15) THEN
                ev(1,3) = 0D0
                ev(2,3) = 1D0
                ev(3,3) = -1D0/tmp(2,3)
            ELSE
                ev(1,3) = 1D0
                ev(3,3) = -1D0/tmp(1,3)
                ev(2,3) = -ev(3,3)*tmp(2,3)/tmp(2,2)
            ENDIF

!           Arbitrary other Eigenvector
            tmp = A - e(1)*I
            ev(1,1) = 1D0
            ev(2,1) = 1D0
            ev(3,1) = -(tmp(1,1) + tmp(1,2))/tmp(1,3)

            ev(:,2) = CROSS(ev(:,1), ev(:,3))

            ev(:,1) = ev(:,1)/NORM2(ev(:,1))
            ev(:,2) = ev(:,2)/NORM2(ev(:,2))
            ev(:,3) = ev(:,3)/NORM2(ev(:,3))

!       My own rref for the others
        ELSE
!           First let's find the normal eigenvector b/c it usually doesn't have bad behavior
            tmp = A - e(3)*I
            tmp = RREF3(tmp)
!           Rare case x-value = 0
            IF(ABS(tmp(1,3)).lt.1e-15) THEN
                ev(1,3) = 0D0
                ev(2,3) = 1D0
                ev(3,3) = -1D0/tmp(2,3)
            ELSE
                ev(1,3) = 1D0
                ev(3,3) = -tmp(1,1)/tmp(1,3)
                ev(2,3) = -ev(3,3)*tmp(2,3)/tmp(2,2)
            ENDIF
            ev(:,3) = ev(:,3)/NORM2(ev(:,3))

!           Now things can get a little tricky if the x or y coordinate in the following
!           eigenvectors should be zero b/c of precision, so we need to do some careful checks
!           If the Y or X should be zero, the other perpendicular eigenvector WILL have some
!           issues with precision, as the RREF in the respective Y/X column will be O(1e-14)
!           or so, so we want to avoid that one and calc the other one. If either sums are small,
!           Just do the other
            tmp = A-e(1)*I
            tmp = RREF3(tmp)
            chkx = SUM(ABS(tmp(1,:)))
            chky = SUM(ABS(tmp(2,:)))

!           Values too low, don't risk it
            IF(chkx .lt. 1e-10 .or. chky .lt. 1e-10) THEN
                tmp = A-e(2)*I
                tmp = RREF3(tmp)
!               x = 0 case
                IF(ABS(tmp(1,3)).lt.1e-14) THEN
                    ev(1,2) = 0D0
                    ev(2,2) = 1D0
                    ev(3,2) = -tmp(2,2)/tmp(2,3)
                ELSE
                    ev(1,2) = 1D0
                    ev(3,2) = -tmp(1,1)/tmp(1,3)
                    ev(2,2) = -ev(3,2)*tmp(2,3)/tmp(2,2)
                ENDIF
                ev(:,2) = ev(:,2)/NORM2(ev(:,2))
!               Cross for last one
                ev(:,1) = CROSS(ev(:,2), ev(:,3))
                ev(:,1) = ev(:,1)/NORM2(ev(:,1))

            ELSE
!               x = 0 case
                IF(ABS(tmp(1,3)).lt.1e-14) THEN
                    ev(1,1) = 0D0
                    ev(2,1) = 1D0
                    ev(3,1) = -tmp(2,2)/tmp(2,3)
                ELSE
                    ev(1,1) = 1D0
                    ev(3,1) = -tmp(1,1)/tmp(1,3)
                    ev(2,1) = -ev(3,1)*tmp(2,3)/tmp(2,2)
                ENDIF
                ev(:,1) = ev(:,1)/NORM2(ev(:,1))
!               Cross for last one
                ev(:,2) = CROSS(ev(:,1), ev(:,3))
                ev(:,2) = ev(:,2)/NORM2(ev(:,2))
            ENDIF
            
        ENDIF
    ENDIF
END SUBROUTINE EIG3

! -------------------------------------------------------------------------!
! Inner product of 3x1 vector with a 3x3 matrix
PURE FUNCTION INNER3_33(U, V) RESULT(W)
    REAL(KIND = 8), INTENT(IN) :: U(3), V(3,3)
    REAL(KIND = 8) :: W(3)
    W(1) = U(1)*V(1,1) + U(2)*V(2,1) + U(3)*V(3,1)
    W(2) = U(1)*V(1,2) + U(2)*V(2,2) + U(3)*V(3,2)
    W(3) = U(1)*V(1,3) + U(2)*V(2,3) + U(3)*V(3,3)
END FUNCTION
! -------------------------------------------------------------------------!
! Inverts 3x3 Matrix
PURE FUNCTION INVERT33(A) RESULT(B)
    REAL(KIND = 8), INTENT(IN) :: A(3,3)
    REAL(KIND = 8) :: B(3,3), DET

    DET = &
      A(1,1)*A(2,2)*A(3,3)  &
    - A(1,1)*A(2,3)*A(3,2)  &
    - A(1,2)*A(2,1)*A(3,3)  &
    + A(1,2)*A(2,3)*A(3,1)  &
    + A(1,3)*A(2,1)*A(3,2)  &
    - A(1,3)*A(2,2)*A(3,1)
    
    B(1,1) =  A(2,2)*A(3,3) - A(2,3)*A(3,2)
    B(1,2) = -A(2,1)*A(3,3) + A(2,3)*A(3,1)
    B(1,3) =  A(2,1)*A(3,2) - A(2,2)*A(3,1)
    B(2,1) = -A(1,2)*A(3,3) + A(1,3)*A(3,2)
    B(2,2) =  A(1,1)*A(3,3) - A(1,3)*A(3,1)
    B(2,3) = -A(1,1)*A(3,2) + A(1,2)*A(3,1)
    B(3,1) =  A(1,2)*A(2,3) - A(1,3)*A(2,2)
    B(3,2) = -A(1,1)*A(2,3) + A(1,3)*A(2,1)
    B(3,3) =  A(1,1)*A(2,2) - A(1,2)*A(2,1)

    B = TRANSPOSE(B)/DET

END FUNCTION
! -------------------------------------------------------------------------!
! Puts a 3x3 matrix in (unnormalized) reduced row Echelon form to solve eigenvectors
PURE FUNCTION RREF3(U) RESULT(V)
    REAL(KIND = 8), INTENT(IN) :: U(3,3)
    REAL(KIND = 8) :: V(3,3)
    V(1,:) = U(1,:) - U(1,2)*(U(3,:)/U(3,2))
    V(2,:) = U(2,:) - U(2,1)*(U(3,:)/U(3,1))
    V(3,:) = 0D0
END FUNCTION RREF3

! -------------------------------------------------------------------------!
! Takes three points and does quadratic interpolation
PURE FUNCTION QInterp(x, y, xo) RESULT(yo)
    REAL(KIND = 8), INTENT(IN) :: x(3), xo
    REAL(KIND = 8), INTENT(IN), ALLOCATABLE :: y(:,:,:)
    REAL(KIND = 8), ALLOCATABLE :: yo(:,:)
    yo = y(1,:,:)*(xo - x(2))*(xo - x(3))/((x(1) - x(2))*(x(1) - x(3))) &
       + y(2,:,:)*(xo - x(1))*(xo - x(3))/((x(2) - x(1))*(x(2) - x(3))) &
       + y(3,:,:)*(xo - x(1))*(xo - x(2))/((x(3) - x(1))*(x(3) - x(2)))
END FUNCTION QInterp

! -------------------------------------------------------------------------!
! Uses above function to get current velocity gradient
PURE FUNCTION VelInterp(G, t, nts, kfr) RESULT(dU)
    REAL(KIND = 8), INTENT(IN) :: t, kfr
    INTEGER, INTENT(IN) :: nts
    INTEGER :: kts
    REAL(KIND = 8), INTENT(IN) :: G(:,:,:)
    REAL(KIND = 8) :: dU(3,3)
    REAL(KIND = 8) :: xs(3)
    REAL(KIND = 8), ALLOCATABLE :: ys(:,:,:)

    ALLOCATE(ys(3,3,3))

    kts = FLOOR(t/kfr)
    xs(1) = REAL(kts,8)*kfr
    xs(2) = xs(1) + kfr
    ys(1,:,:) = G(kts + 1,:,:)
    ys(2,:,:) = G(kts + 2,:,:)

!       Use two points to right and one left, unless you're at the end
!       Exception if there's just 1 timestep needed
    IF(t + 2D0*kfr .lt. nts*kfr .or. nts .eq. 1) THEN
            xs(3) = xs(2) + kfr
            ys(3,:,:) = G(kts + 3,:,:)
    ELSE
            xs(3) = xs(1) - kfr
            ys(3,:,:) = G(kts - 1,:,:)
    ENDIF
    dU = QInterp(xs,ys,t)
END FUNCTION VelInterp

! -------------------------------------------------------------------------!
! Gets values after specified input

SUBROUTINE READ_GRINT_DOUB(x, filen, srch)
    CHARACTER (len=*), INTENT(IN) :: filen, srch
    CHARACTER (len=1000) :: text
    CHARACTER (len=20) :: word
    REAL(KIND = 8), INTENT(OUT)  :: x
    INTEGER ierr

    OPEN(unit = 88, file = filen, action = 'read')
    DO
        READ (88,"(a)",iostat=ierr) text ! read line into character variable
        IF (ierr /= 0) EXIT
        READ (text,*) word ! read first word of line
        IF (word == srch) THEN ! found search string at beginning of line
           READ (text,*) word,x
           EXIT
        ENDIF
    ENDDO
    CLOSE(88)
END SUBROUTINE READ_GRINT_DOUB

SUBROUTINE READ_GRINT_INT(x, filen, srch)
    CHARACTER (len=*), INTENT(IN) :: filen, srch
    CHARACTER (len=1000) :: text
    CHARACTER (len=20) :: word
    INTEGER, INTENT(OUT)  :: x
    INTEGER ierr

    OPEN(unit = 88, file = filen, action = 'read')
    DO
        READ (88,"(a)",iostat=ierr) text ! read line into character variable
        IF (ierr /= 0) EXIT
        READ (text,*) word ! read first word of line
        IF (word == srch) THEN ! found search string at beginning of line
           READ (text,*) word,x
           EXIT
        ENDIF
    ENDDO
    CLOSE(88)
END SUBROUTINE READ_GRINT_INT

SUBROUTINE READ_GRINT_CHAR(x, filen, srch)
    CHARACTER (len=*), INTENT(IN) :: filen, srch
    CHARACTER (len=1000) :: text
    CHARACTER (len=20) :: word, word2
    CHARACTER (len=:), ALLOCATABLE, INTENT(OUT) :: x
    INTEGER ierr

    OPEN(unit = 88, file = filen, action = 'read')
    DO
        READ (88,"(a)",iostat=ierr) text ! read line into character variable
        IF (ierr /= 0) EXIT
        READ (text,*) word ! read first word of line
        IF (word == srch) THEN ! found search string at beginning of line
           READ (text,*) word,word2
           x = trim(word2)
           EXIT
        ENDIF
    ENDDO
    CLOSE(88)
END SUBROUTINE READ_GRINT_CHAR
! -------------------------------------------------------------------------!
! FFT 1D wrapper
FUNCTION FFT1(c) RESULT(ch)
    COMPLEX(KIND = 8), ALLOCATABLE :: c(:), ch(:)
    INTEGER :: N, ier, LENSAV, LENWRK, LENC, INC
    REAL(KIND = 8), ALLOCATABLE :: WSAVE(:), WORK(:)

    N = SIZE(c)
    INC = 1
    LENC = N
    LENWRK = 2*N
    LENSAV = 2*N + INT(LOG(REAL(N))) + 4
    ALLOCATE(WSAVE(LENSAV), WORK(LENWRK))
    CALL ZFFT1I(N, WSAVE, LENSAV, ier)
    CALL ZFFT1B(N, 1, c, LENC, WSAVE, lensav, work, lenwrk, ier)
    ch = c
    
END FUNCTION FFT1
! -------------------------------------------------------------------------!
! iFFT 1D wrapper
FUNCTION iFFT1(c) RESULT(ch)
    COMPLEX(KIND = 8), ALLOCATABLE :: c(:), ch(:)
    INTEGER :: N, ier, LENSAV, LENWRK, LENC, INC
    REAL(KIND = 8), ALLOCATABLE :: WSAVE(:), WORK(:)

    N = SIZE(c)
    INC = 1
    LENC = N
    LENWRK = 2*N
    LENSAV = 2*N + INT(LOG(REAL(N))) + 4
    ALLOCATE(WSAVE(LENSAV), WORK(LENWRK))
    CALL ZFFT1I(N, WSAVE, LENSAV, ier)
    CALL ZFFT1F(N, 1, c, LENC, WSAVE, lensav, work, lenwrk, ier)
    ch = c
    
END FUNCTION iFFT1

! -------------------------------------------------------------------------!
! FFT 3D cubic wrapper
FUNCTION FFT3(c, WSAVEin) RESULT(ch)
    COMPLEX(KIND = 8), ALLOCATABLE :: c(:,:,:), ch(:,:,:), ctmp(:)
    INTEGER :: N, ier, LENSAV, LENWRK, LENC, INC, i, j, k, N3(3)
    REAL(KIND = 8), ALLOCATABLE :: WORK(:), WSAVE(:)
    REAL(KIND = 8), ALLOCATABLE, OPTIONAL :: WSAVEin(:)

    N3 = SHAPE(c)
    N = N3(1)
    INC = 1
    LENC = N
    LENWRK = 2*N
    LENSAV = 2*N + INT(LOG(REAL(N))) + 4
    ALLOCATE(WSAVE(LENSAV), WORK(LENWRK), ctmp(N))
    IF(.not. present(WSAVEin)) THEN
        CALL ZFFT1I(N, WSAVE, LENSAV, ier)
    ELSE
        WSAVE = WSAVEin
    ENDIF
    
    DO k = 1,N
        DO i = 1,N
            ctmp = c(i,:,k)
            CALL ZFFT1B(N, 1, ctmp, LENC, WSAVE, lensav, work, lenwrk, ier)
            c(i,:,k) = ctmp
        ENDDO
    ENDDO
    DO k = 1,N
        DO j = 1,N
            ctmp = c(:,j,k)
            CALL ZFFT1B(N, 1, ctmp, LENC, WSAVE, lensav, work, lenwrk, ier)
            c(:,j,k) = ctmp
        ENDDO
    ENDDO
    DO i = 1,N
        DO j = 1,N
            ctmp = c(i,j,:)
            CALL ZFFT1B(N, 1, ctmp, LENC, WSAVE, lensav, work, lenwrk, ier)
            c(i,j,:) = ctmp
        ENDDO
    ENDDO
    ch = c
END FUNCTION FFT3
! -------------------------------------------------------------------------!
! iFFT 3D cubic wrapper 
FUNCTION iFFT3(c, WSAVEin) RESULT(ch)
    COMPLEX(KIND = 8), ALLOCATABLE :: c(:,:,:), ch(:,:,:), ctmp(:)
    INTEGER :: N, ier, LENSAV, LENWRK, LENC, INC, i, j, k, N3(3)
    REAL(KIND = 8), ALLOCATABLE :: WORK(:), WSAVE(:)
    REAL(KIND = 8), ALLOCATABLE, OPTIONAL :: WSAVEin(:)

    N3 = SHAPE(c)
    N = N3(1)
    INC = 1
    LENC = N
    LENWRK = 2*N
    LENSAV = 2*N + INT(LOG(REAL(N))) + 4
    ALLOCATE(WSAVE(LENSAV), WORK(LENWRK), ctmp(N))
    IF(.not. present(WSAVEin)) THEN
        CALL ZFFT1I(N, WSAVE, LENSAV, ier)
    ELSE
        WSAVE = WSAVEin
    ENDIF

    DO i = 1,N
        DO j = 1,N
            ctmp = c(i,j,:)
            CALL ZFFT1F(N, 1, ctmp, LENC, WSAVE, lensav, work, lenwrk, ier)
            c(i,j,:) = ctmp
        ENDDO
    ENDDO
    DO k = 1,N
        DO j = 1,N
            ctmp = c(:,j,k)
            CALL ZFFT1F(N, 1, ctmp, LENC, WSAVE, lensav, work, lenwrk, ier)
            c(:,j,k) = ctmp
        ENDDO
    ENDDO
    DO k = 1,N
        DO i = 1,N
            ctmp = c(i,:,k)
            CALL ZFFT1F(N, 1, ctmp, LENC, WSAVE, lensav, work, lenwrk, ier)
            c(i,:,k) = ctmp
        ENDDO
    ENDDO
    ch = c
END FUNCTION iFFT3

! -------------------------------------------------------------------------!
! Shifts 0th wavemode to center for 3D
SUBROUTINE FFTSHIFT(A)
    COMPLEX(KIND = 8), DIMENSION(:,:,:), INTENT(INOUT) :: A
    COMPLEX(KIND = 8), ALLOCATABLE :: B(:,:,:)
    INTEGER :: N3(3), N
    N3 = SHAPE(A)
    N = N3(1)
    ALLOCATE(B(N,N,N))

!   Even case
    IF(MOD(N,2).eq.0) THEN
        B(1:N/2  , 1:N/2  , 1:N/2  ) = A(N/2+1:N, N/2+1:N, N/2+1:N)
        B(N/2+1:N, 1:N/2  , 1:N/2  ) = A(1:N/2  , N/2+1:N, N/2+1:N)
        B(1:N/2  , N/2+1:N, 1:N/2  ) = A(N/2+1:N, 1:N/2  , N/2+1:N)
        B(N/2+1:N, N/2+1:N, 1:N/2  ) = A(1:N/2  , 1:N/2  , N/2+1:N)
        B(1:N/2  , 1:N/2  , N/2+1:N) = A(N/2+1:N, N/2+1:N, 1:N/2  )
        B(N/2+1:N, 1:N/2  , N/2+1:N) = A(1:N/2  , N/2+1:N, 1:N/2  )
        B(1:N/2  , N/2+1:N, N/2+1:N) = A(N/2+1:N, 1:N/2  , 1:N/2  )
        B(N/2+1:N, N/2+1:N, N/2+1:N) = A(1:N/2  , 1:N/2  , 1:N/2  )
    ELSE
        print *, 'Havent done odd fftshift yet' !!!!!!!!!!!!!
        stop
    ENDIF
    A = B
END SUBROUTINE FFTSHIFT

! -------------------------------------------------------------------------!
! Functions to calculate the kernels
FUNCTION Gij(r,eye) RESULT(A)
    REAL(KIND = 8) r(3), A(3,3), eye(3,3), mri
    mri = 1D0/(sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
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
    mri = 1D0/(sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
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
! Cell-cell interatction forces
FUNCTION Morse(r, rn, De, r0, B) RESULT(f)
    REAL(KIND = 8) :: r(3), rn, De, r0, B, f(3)
    f = -2D0*De*B*(exp(B*(r0-rn)) - exp(2D0*B*(r0-rn)))*(r/rn)
END FUNCTION Morse

FUNCTION LJ(r, rn, e, s) RESULT(f)
    REAL(KIND = 8) :: r(3), rn, e, s, f(3)
    f = -24D0*e*(s**6)*(rn**6 - 2D0*s**6)/(rn**13)*(r/rn)
END FUNCTION LJ

! -------------------------------------------------------------------------!
! Given basis vectors and point, gives partial coordinates
FUNCTION bvPartial(bv, x) RESULT(px)
    REAL(KIND = 8) bv(3,3), x(3), px(3), bv1(3), bv2(3), bv3(3), d, pd, n(3)
    
!   Here's how: the cell is composed of 6 planes, or
!   3 pairs of parallel planes. We need to check if the point
!   is between each pair of these planes. For each pair, one of
!   the planes goes through (0,0) which are the ones we'll look at.
!   There are 3 basis vectors, and two basis vectors define a plane.
!   The projection of the third vector onto the normal defines the
!   distance between the pair of parallel planes. We just need to
!   check the points distance to this distance, and make sure it's
!   smaller (paying attention to the sign)

!   Get each basis vector so it's a little easier to work with
    bv1 = bv(:,1)
    bv2 = bv(:,2)
    bv3 = bv(:,3)

!   Get normal vector to first plane
    n = CROSS(bv1,bv2)
    n = n/NORM2(n)
!   Get projected distance of 3rd vector and point
    d  = DOT(n, bv3)
    pd = DOT(n, x)
    px(3) = pd/d

!   And just repeat for the other two directions
    n = CROSS(bv2,bv3)
    n = n/NORM2(n)
    d  = DOT(n, bv1)
    pd = DOT(n, x)
    px(1) = pd/d

    n = CROSS(bv3,bv1)
    n = n/NORM2(n)
    d  = DOT(n, bv2)
    pd = DOT(n, x)
    px(2) = pd/d

END FUNCTION bvPartial

! -------------------------------------------------------------------------!
! Given a point, the basis vectors, and a distance, gives the boxes that
! must be checked for minimum periodic distance
FUNCTION bvBxs(bv, h, x) RESULT(bxs)
    REAL(KIND = 8) bv(3,3), h, x(3), nv(3), d, dpt
    INTEGER, ALLOCATABLE :: bxs(:,:)
    INTEGER fcs(6), bxst(3,27), it

    fcs = 0
    it = 0
!   Do each individual face and find distance, then see if that distance is
!   is less than h
    nv = CROSS(bv(:,1), bv(:,2))
    d = DOT(nv, x)/NORM2(nv)
    IF(d .lt. h) fcs(1) = 1

!   Now do for parallel plane, and the distance must be greater than -h
    dpt = -DOT(nv, bv(:,3))
    d = (DOT(nv, x) + dpt)/NORM2(nv)
    IF(d .gt. -h) fcs(2) = 1

!   And then the other 2 faces
    nv = CROSS(bv(:,2), bv(:,3))
    d = DOT(nv, x)/NORM2(nv)
    IF(d .lt. h) fcs(3) = 1
    
    dpt = -DOT(nv, bv(:,1))
    d = (DOT(nv, x) + dpt)/NORM2(nv)
    IF(d .gt. -h) fcs(4) = 1

    nv = CROSS(bv(:,3), bv(:,1))
    d = DOT(nv, x)/NORM2(nv)
    IF(d .lt. h) fcs(5) = 1
    
    dpt = -DOT(nv, bv(:,2))
    d = (DOT(nv, x) + dpt)/NORM2(nv)
    IF(d .gt. -h) fcs(6) = 1

!   Now we want a list of boxes we need to check.
!   Just basically done via brute force

!   Bottom
    IF(fcs(1) .eq. 1) THEN
        it = it + 1
        bxst(:,it) = (/ 0, 0,-1/)

!       Left
        IF(fcs(3) .eq. 1) THEN
            it = it + 1
            bxst(:,it) = (/-1, 0,-1/)

!           Back
            IF(fcs(5) .eq. 1) THEN
                it = it + 1
                bxst(:,it) = (/-1,-1,-1/)
!           Front
            ELSEIF(fcs(6).eq.1) THEN
                it = it + 1
                bxst(:,it) = (/-1, 1,-1/)
            ENDIF
!       Right
        ELSEIF(fcs(4) .eq. 1) THEN
            it = it + 1
            bxst(:,it) = (/ 1, 0,-1/)

!           Back
            IF(fcs(5) .eq. 1) THEN
                it = it + 1
                bxst(:,it) = (/ 1,-1,-1/)
!           Front
            ELSEIF(fcs(6).eq.1) THEN
                it = it + 1
                bxst(:,it) = (/ 1, 1,-1/)
            ENDIF
        ENDIF

!   Top
    ELSEIF(fcs(2) .eq. 1) THEN
        it = it + 1
        bxst(:,it) = (/ 0, 0, 1/)

!       Left
        IF(fcs(3) .eq. 1) THEN
            it = it + 1
            bxst(:,it) = (/-1, 0, 1/)

!           Back
            IF(fcs(5) .eq. 1) THEN
                it = it + 1
                bxst(:,it) = (/-1,-1, 1/)
!           Front
            ELSEIF(fcs(6).eq.1) THEN
                it = it + 1
                bxst(:,it) = (/-1, 1, 1/)
            ENDIF
!       Right
        ELSEIF(fcs(4) .eq. 1) THEN
            it = it + 1
            bxst(:,it) = (/ 1, 0, 1/)

!           Back
            IF(fcs(5) .eq. 1) THEN
                it = it + 1
                bxst(:,it) = (/ 1,-1, 1/)
!           Front
            ELSEIF(fcs(6).eq.1) THEN
                it = it + 1
                bxst(:,it) = (/ 1, 1, 1/)
            ENDIF
        ENDIF
    ENDIF

!   Left
    If(fcs(3) .eq. 1) THEN
        it = it + 1
        bxst(:,it) = (/-1, 0, 0/)

!       Back
        IF(fcs(5) .eq. 1) THEN
            it = it + 1
            bxst(:,it) = (/-1,-1, 0/)
!       Front
        ELSEIF(fcs(6).eq.1) THEN
            it = it + 1
            bxst(:,it) = (/-1, 1, 0/)
        ENDIF
!   Right
    ELSEIF(fcs(4) .eq. 1) THEN
        it = it + 1
        bxst(:,it) = (/ 1, 0, 0/)

!       Back
        IF(fcs(5) .eq. 1) THEN
            it = it + 1
            bxst(:,it) = (/ 1,-1, 0/)
!       Front
        ELSEIF(fcs(6).eq.1) THEN
            it = it + 1
            bxst(:,it) = (/ 1, 1, 0/)
        ENDIF
    ENDIF

!   Back
    IF(fcs(5) .eq. 1) THEN
        it = it + 1
        bxst(:,it) = (/ 0,-1, 0/)
!   Front
    ELSEIF(fcs(6).eq.1) THEN
        it = it + 1
        bxst(:,it) = (/ 0, 1, 0/)
    ENDIF

!   Last four:
    IF((fcs(1) .eq. 1) .and. (fcs(3) .eq. 0) .and. (fcs(4) .eq. 0) .and. (fcs(5) .eq. 1)) THEN
        it = it + 1
        bxst(:,it) = (/ 0,-1,-1/)
    ENDIF
    IF((fcs(1) .eq. 1) .and. (fcs(3) .eq. 0) .and. (fcs(4) .eq. 0) .and. (fcs(6) .eq. 1)) THEN
        it = it + 1
        bxst(:,it) = (/ 0, 1,-1/)
    ENDIF
    IF((fcs(2) .eq. 1) .and. (fcs(3) .eq. 0) .and. (fcs(4) .eq. 0) .and. (fcs(5) .eq. 1)) THEN
        it = it + 1
        bxst(:,it) = (/ 0,-1, 1/)
    ENDIF
    IF((fcs(2) .eq. 1) .and. (fcs(3) .eq. 0) .and. (fcs(4) .eq. 0) .and. (fcs(6) .eq. 1)) THEN
        it = it + 1
        bxst(:,it) = (/ 0, 1, 1/)
    ENDIF
    
    ALLOCATE(bxs(3,it))
    bxs = bxst(:,1:it)
END FUNCTION bvBxs

! -------------------------------------------------------------------------!
! Periodic functions to calculate the kernels
! Arguments:
!   r:        Distance from eval point to source point in primary cell
!   bxs:      Number of boxes to sum over
!   bv:       Basis vectors representing original cell
!   eye:      Identity matrix
!   n (Tij):  Normal vector
!   wg:       Gauss weight at point
!   wgprim:   Gauss weight at point in primary cell
FUNCTION PGij(r, bxs, bv, eye, xi, fourier) RESULT(A)
    REAL(KIND = 8) :: r(3), bv(3,3), eye(3,3), A(3,3), xi, tau, &
                      rcur(3), kv(3,3), kcur(3)
    INTEGER :: bxs, i, j, k
    LOGICAL, OPTIONAL :: fourier
    LOGICAL :: flag

    flag  = .false.
    IF(PRESENT(fourier)) flag = fourier

    A = 0D0

!   Volume of lattice unit cell (likely could pass as arg, but not huge deal)
    tau = DOT(CROSS(bv(:,1), bv(:,2)), bv(:,3))

    IF(flag) THEN
        kv(:,1) = 2*pi/tau*cross(bv(:,2),bv(:,3))
        kv(:,2) = 2*pi/tau*cross(bv(:,3),bv(:,1))
        kv(:,3) = 2*pi/tau*cross(bv(:,1),bv(:,2))
    ENDIF

!   We do the real and Fourier sums in the same loops
    DO i = -bxs, bxs
        DO j = -bxs, bxs
            DO k = -bxs, bxs

!               Fourier part
                IF(flag .and. .not.((i .eq. 0) .and. (j .eq. 0) .and. (k.eq.0)) ) THEN
                    kcur = i*kv(:,1) + j*kv(:,2) + k*kv(:,3)
                    A = A + FOURIER_G_HAS(kcur, xi, eye)/tau*COS(DOT(kcur,r))
                ENDIF

!               Real part (get current vector first)
                rcur = r + i*bv(:,1) + j*bv(:,2) + k*bv(:,3)

!               Check cutoff in here as well
                IF(.not. flag .and. norm2(rcur)*xi.gt.3.5) CYCLE

!               Cycle if the contribution will be small enough
                A    = A + REAL_G_HAS(rcur, xi, eye)

            ENDDO
        ENDDO
    ENDDO

    CONTAINS
!   Hasimotos
    FUNCTION REAL_G_HAS(r, xi, eye) RESULT(A)
        REAL(KIND = 8) :: r(3), mr, xi, A(3,3), C, D, eye(3,3), er
        mr = SQRT(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
        er = mr*xi
        C = ERFC(er) - 2D0*ispi*er*EXP(-er*er)
        D = ERFC(er) + 2D0*ispi*er*EXP(-er*er)
        A = eye*C/mr + OUTER(r,r)*D/(mr*mr*mr)
    END FUNCTION REAL_G_HAS

    FUNCTION FOURIER_G_HAS(k, xi, eye) RESULT(A)
        REAL(KIND = 8) k(3), xi, eye(3,3), A(3,3), kn, w
        kn = SQRT(k(1)*k(1) + k(2)*k(2) + k(3)*k(3))
        w = kn/xi;
        A = 8D0*PI/(xi*xi*xi*xi)*(1D0/(w*w*w*w) + 0.25D0/(w*w) )*((kn*kn)*eye - OUTER(k,k))*exp(-0.25D0*w*w);
    END FUNCTION FOURIER_G_HAS
END FUNCTION PGij

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
FUNCTION PTij(r, bxs, bv, n, eye, xi, x, fourier) RESULT(A)
    REAL(KIND = 8) :: r(3), bv(3,3), n(3), eye(3,3), A(3,3), xi, tau, &
                      rcur(3), kv(3,3), kcur(3), xr(3)
    INTEGER :: bxs, i, j, k
    LOGICAL, OPTIONAL :: fourier
    REAL(KIND = 8), OPTIONAL :: x(3)
    LOGICAL :: flag

    flag = .false.
    IF(PRESENT(fourier)) flag = fourier
    IF(.not. PRESENT(x)) THEN
        xr = r
    ELSE
        xr = x
    ENDIF 

    A  = 0D0

!   Volume of lattice unit cell (likely could pass as arg, but not huge deal)
    tau = DOT(CROSS(bv(:,1), bv(:,2)), bv(:,3))

    IF(flag) THEN
        kv(:,1) = 2*pi/tau*cross(bv(:,2),bv(:,3))
        kv(:,2) = 2*pi/tau*cross(bv(:,3),bv(:,1))
        kv(:,3) = 2*pi/tau*cross(bv(:,1),bv(:,2))
    ENDIF

!   We do the real and Fourier sums in the same loops
    DO i = -bxs, bxs
        DO j = -bxs, bxs
            DO k = -bxs, bxs

!               Real part (get current vector first)
                rcur = r + i*bv(:,1) + j*bv(:,2) + k*bv(:,3)

!               Check cutoff in here as well
                IF(.not. flag .and. norm2(rcur)*xi.gt.3.5) CYCLE

                ! A   = A + REAL_T_HAS(rcur, xi, n, eye)
                A   = A + REAL_T_MAR(rcur, xi, n, eye)

!               Fourier part
                IF(flag .and. .not.((i .eq. 0) .and. (j .eq. 0) .and. (k.eq.0)) ) THEN
                    kcur = i*kv(:,1) + j*kv(:,2) + k*kv(:,3)
                    ! A = A + FOURIER_T_HAS(kcur, xi, n, eye)/tau*SIN(DOT(kcur,r))
                    A = A - FOURIER_T_MAR(kcur, xi, n, eye)/tau*SIN(DOT(kcur,r))
                ENDIF
            ENDDO
        ENDDO
    ENDDO
!   Non-periodic portion comgin from pressure to balance net force
    A = A - 8D0*PI/tau*OUTER(xr,n) ! Doesn't actually matter r or surface, as constant disappears in integral
    
    CONTAINS
!   Hasimotos
    FUNCTION REAL_T_HAS(r, xi, n, eye) RESULT(A)
        REAL(KIND = 8) :: r(3), mr, xi, n(3), A(3,3), C, D, eye(3,3), er, rh(3), xer, rdn
        mr = SQRT(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
        rh = r/mr
        er = mr*xi
        xer = EXP(-er*er)
        C = -6D0*ERFC(er)/(mr*mr) - xi*ispi/mr*(12D0 + 8D0*er*er - 16D0*er*er*er*er)*xer
        D = 8D0*xi*xi*xi*mr*ispi*(2D0 - er*er)*xer
        rdn = DOT(rh,n)
        A = C*(OUTER(rh,rh)*rdn) + D*(eye*rdn + OUTER(n,rh) + OUTER(rh,n))
    END FUNCTION REAL_T_HAS

    FUNCTION FOURIER_T_HAS(k, xi, n, eye) RESULT(A)
        REAL(KIND = 8) k(3), xi, n(3), eye(3,3), kn, w, kdn, Q2, xer, Q1(3,3), A(3,3)
        kn = SQRT(k(1)*k(1) + k(2)*k(2) + k(3)*k(3))
        w = kn/xi
        kdn = DOT(k,n)
        Q1 = (-2D0/(kn*kn*kn*kn)*OUTER(k,k)*kdn + 1D0/(kn*kn)*(OUTER(k,n) + OUTER(n,k) + eye*kdn))
        Q2 = 8D0 + 2D0*w*w + w*w*w*w
        xer = EXP(-0.25D0*w*w)
        A = -PI*Q1*Q2*xer
    END FUNCTION FOURIER_T_HAS

!   Try with the Marin decomp
    FUNCTION REAL_T_MAR(r, xi, n, eye) RESULT(A)
        REAL(KIND = 8) :: r(3), mr, xi, n(3), A(3,3), C, D, eye(3,3), er, rh(3), xer, rdn
        mr = SQRT(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
        rh = r/mr
        er = mr*xi
        xer = EXP(-er*er)
        C = -6D0*ERFC(er)/(mr*mr) - xi*ispi/mr*(12D0 + 8D0*er*er)*xer
        D = 4D0*xi*xi*ispi*xer*er
        rdn = DOT(rh,n)
        A = C*(OUTER(rh,rh)*rdn) + D*(eye*rdn + OUTER(n, rh) + OUTER(rh, n))
    END FUNCTION REAL_T_MAR

    FUNCTION FOURIER_T_MAR(k, xi, n, eye) RESULT(A)
        REAL(KIND = 8) k(3), xi, n(3), eye(3,3), kn, w, kdn, Q2, xer, Q1(3,3), A(3,3)
        kn = SQRT(k(1)*k(1) + k(2)*k(2) + k(3)*k(3))
        w = kn/xi
        kdn = DOT(k,n)
        xer = EXP(-0.25D0*w*w)
        Q1 = (-2D0/(kn*kn*kn*kn)*OUTER(k,k)*kdn + 1D0/(kn*kn)*(OUTER(k,n) + OUTER(n,k) + eye*kdn))
        Q2 = 8D0*PI*(1 + w*w*0.25D0)
        A = Q1*Q2*xer
    END FUNCTION FOURIER_T_MAR

END FUNCTION PTij

END MODULE UTILMOD