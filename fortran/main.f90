PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, gradline(3,3), t
COMPLEX(kind = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:)
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat

! For a = 1, V = 4.18904795321178, SA = 16.8447913187040, sphere 6.50088174342271

! Read in input file
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

! Initialize & allocate
cell = cellType(filein)
CALL cpu_time(tic)

ALLOCATE(xmnt(3,(cell%p+1)*(cell%p+1)), umnt(3,(cell%p+1)*(cell%p+1)))

! Read velocity gradient file
! OPEN(unit = 7, file = cell%gradfile, status = 'old', action = 'read', iostat=stat)

print*, 'Initialized'
t = 0D0
! Time step loop
DO i = 1,cell%NT
        IF(i.lt.300) THEN
                cell%dU = 0d0
        ELSE
                cell%dU(1,3) = 1d0
        ENDIF
!       The time distance between successive rows is 0.1 of the Kolmogorov time scale.
!       To normalize these gradients, you could multiply them by the Kolmogorov time scale, which is 2e-6

!       Get velocity gradient of time step and normalize
        ! READ(7, *) cell%dU
        ! cell%dU = cell%dU*2D-6

!       Get surface derivatives, then stress froom deformation, then motion from fluid
        CALL cell%derivs()
        CALL cell%stress()
        CALL cell%fluid()
        cell%umn(:,1) = 0D0

        umnt = cell%umn
        xmnt = cell%xmn

        ! IF(cell%vol().gt.4.18904795321178 .and. i.lt.2500) umnt = umnt - 0.1D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))

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
        write(*,'(I,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)'), &
        i, t, MAXVAL((ABS(cell%ff))), MAXVAL((ABS(cell%umn))), cell%vol(), &
        cell%SA(), cell%Eb*cell%intg((2*cell%fab(1,:,:) - cell%C0)*(2*cell%fab(1,:,:) - cell%C0))/2D0
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN