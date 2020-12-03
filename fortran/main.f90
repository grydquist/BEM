PROGRAM MAIN
USE FFTMOD
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: x(4),w(4), xv(4,8), a(8), b(8), tmp(3,3),lams(3), ev(3,3),t2(3), &
        tmp2(3,3), dt, pivot(2), aa(3,3), bb(3), cc(102), tic, toc
TYPE(cellType) :: cell
! REAL(KIND = 8), ALLOCATABLE :: tmp(:)
! COMPLEX(KIND = 8) :: x1mn(3,16)

INTEGER ::  iserr, i, rc

!================================================================================!
!================================= DECLARATIONS =================================!
!================================================================================!

!-----------------------Material properties/parameters-----------------------!
! REAL(KIND = 8) :: mu, lam, B, C, Eb
! INTEGER, PARAMETER :: nsd = 3, nsb = 4

! !-----------------------Simulation items-----------------------!
INTEGER :: NT
! REAL(KIND = 8) :: dU(nsd, nsd)

NT = 1000
CALL cpu_time(tic)
cell = cellType(10, 4)
print*, 'Initialized'
! print *, cell%xmn(3,:)
CALL cell%stress()
call cell%fluid()
DO i = 1,NT
        cell%cts = cell%cts + 1
        print *, i, MAXVAL((ABS(cell%ff)))
        print *, '   '
        cell%xmn = cell%xmn + cell%umn*cell%dt               !!!! make update step
        cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
        call cell%write()

        CALL cell%derivs()
        CALL cell%stress()
        CALL cell%fluid()
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN