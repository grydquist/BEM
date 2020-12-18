PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: x(4),w(4), xv(4,8), a(8), b(8), tmp(3,3),lams(3), ev(3,3),t2(3), &
        tmp2(3,3), pivot(2), aa(3,3), bb(3), cc(102), tic, toc

COMPLEX(KIND = 8) :: A2(3,3), b2(3), ut(3)
COMPLEX :: swrk(9)
REAL(KIND = 8), ALLOCATABLE :: rwrk(:), wrk(:)
COMPLEX(kind = 8), ALLOCATABLE :: fwrk(:), umnt(:,:), xmnt(:,:)
INTEGER :: IPIV(3), iter, info
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
! REAL(KIND = 8), ALLOCATABLE :: tmp(:)
! COMPLEX(KIND = 8) :: x1mn(3,16)

INTEGER ::  iserr, i, rc, argl, stat, n

!================================================================================!
!================================= DECLARATIONS =================================!
!================================================================================!

!-----------------------Material properties/parameters-----------------------!
! REAL(KIND = 8) :: mu, lam, B, C, Eb
! INTEGER, PARAMETER :: nsd = 3, nsb = 4

! !-----------------------Simulation items-----------------------!
INTEGER :: NT

! A2 = 0D0
! A2(1,1) = ii
! A2(2,2) = 1D0
! A2(3,3) = 1D0
! b2 = 1D0
! CALL zcgesv(3, 1, A2, 3, IPIV, b2, 3, ut, 3, wrk, swrk, rwrk, iter, info)
! print *, ut
! print *, MATMUL(A2, ut)
! stop


! n = 5
! ALLOCATE(rwrk(2*n + INT(LOG(REAL(n))) + 4), fwrk(n), wrk(2*n))
! fwrk = (/0D0, 1d0, 2d0, 3d0, 4D0/)

! CALL ZFFT1I(n, rwrk, 2*n + INT(LOG(REAL(n))) + 4, iserr)
! CALL ZFFT1F(n, 1, fwrk, n, rwrk, 2*n + INT(LOG(REAL(n))) + 4, wrk, 2*n, iserr)
! print *, fwrk*n
! stop

! Read in input file
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

! Initialize
CALL cpu_time(tic)
cell = cellType(filein)
print*, 'Initialized'
ALLOCATE(xmnt(3,(cell%p+1)*(cell%p+1)), umnt(3,(cell%p+1)*(cell%p+1)))

! Time step loop
DO i = 1,cell%NT
!       Get surface derivatives, then stress froom deformation, then motion from fluid
        CALL cell%derivs()
        CALL cell%stress()
        CALL cell%fluid()
        cell%umn(:,1) = 0D0

        umnt = cell%umn
        xmnt = cell%xmn

!       Update and output
        cell%cts = cell%cts + 1
        cell%xmn = xmnt + umnt*cell%dt               !!!! make update step
        cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
        
! !       Second part for midpoint
!         CALL cell%derivs()
!         CALL cell%stress()
!         CALL cell%fluid()
!         cell%umn(:,1) = 0D0

!         umnt = 0.5D0*(umnt + cell%umn)
!         cell%xmn = xmnt + umnt*cell%dt

        CALL cell%write()
        print *, i, MAXVAL((ABS(cell%ff))), MAXVAL((ABS(cell%umn)))

ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN