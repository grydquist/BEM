PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, gradline(3,3), t
COMPLEX(kind = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:)
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat

! Read in input file
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

! Initialize & allocate
CALL cpu_time(tic)
cell = cellType(filein)
ALLOCATE(xmnt(3,(cell%p+1)*(cell%p+1)), umnt(3,(cell%p+1)*(cell%p+1)))

! Read velocity gradient file
! OPEN(unit = 13, file = cell%gradfile, status = 'old', action = 'read', iostat=stat)

print*, 'Initialized'
t = 0D0
! Time step loop
DO i = 1,cell%NT
!       The time distance between successive rows is 0.1 of the Kolmogorov time scale.
!       To normalize these gradients, you could multiply them by the Kolmogorov time scale, which is 2e-6

!       Get velocity gradient of time step and normalize
        ! READ(13, *) cell%dU
        ! cell%dU = cell%dU*2D-6

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
        t = t + cell%dt
        write(*,'(I,X,F8.4,X,X,F8.4,X,F8.4)'), i, t, MAXVAL((ABS(cell%ff))), MAXVAL((ABS(cell%umn)))
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN