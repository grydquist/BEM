MODULE UTILMOD
IMPLICIT NONE
REAL(KIND = 8), PARAMETER :: ispi = 0.56418958354775627928D0, pi = 3.1415926535897932384626433D0
COMPLEX(KIND = 8), PARAMETER :: ii = (0D0,1D0)

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
    ! IF(ABS(A(1,1) - 0.996933541385099D0) .lt. 1e-10) print *, e

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
! Puts a 3x3 matrix in (unnormalized) reduced row Echelon form to solve eigenvectors
PURE FUNCTION RREF3(U) RESULT(V)
    REAL(KIND = 8), INTENT(IN) :: U(3,3)
    REAL(KIND = 8) :: V(3,3)

    V(1,:) = U(1,:) - U(1,2)*(U(3,:)/U(3,2))
    V(2,:) = U(2,:) - U(2,1)*(U(3,:)/U(3,1))
    V(3,:) = 0D0

END FUNCTION RREF3

END MODULE UTILMOD