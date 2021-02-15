PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, t, V0, zm, xs(3), kdt, kfr
COMPLEX(kind = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:)
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, kts, nts
REAL(KIND=8), ALLOCATABLE :: G(:,:,:,:), ys(:,:,:)

! For a = 1, V = 4.18904795321178, SA = 16.8447913187040, sphere 6.50088174342271

! Read in input file
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

! Initialize & allocate
cell = cellType(filein)
CALL cpu_time(tic)

ALLOCATE(xmnt(3,(cell%p+1)*(cell%p+1)), umnt(3,(cell%p+1)*(cell%p+1)))

print *, 'Reading in velocity gradient'
!! ============================
OPEN(1,FILE = cell%gradfile, ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
ALLOCATE(G(nts,3,3,i), ys(3,3,3))
READ(1) G
CLOSE(1)
!! ============================

! Info about flow time scales/vel grad file
kdt = READ_GRINT_DOUB(filein, 'Kolm_time')
kfr = READ_GRINT_DOUB(filein, 'Kolm_frac')

! Let shape equilibrate by running a few time steps with no dU
print *, 'Getting initial shape'
!!!!!!!!!!!!!! Not a good way to do this b/c vel coeffs change based
!!!!!!!!!!!!!! on membrane response time, a function of Ca. Should do
!!!!!!!!!!!!!! something like 0.05 * or / Ca
print *,  'Max velocity coefficient (want less than 0.0005): '
cell%dU = 0D0

!! ============================
DO WHILE(MAXVAL(ABS(umnt)) .gt. 0.005)
        CALL cell%derivs()
        CALL cell%stress() 
        CALL cell%fluid()
        cell%umn(:,1) = 0D0
        umnt = cell%umn
        xmnt = cell%xmn
        cell%xmn = xmnt + umnt*cell%dt
        cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:)) 
        cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
        write(*,'(F8.5)'), MAXVAL(ABS(umnt))
ENDDO
!! ============================
!   Write initial configuration
CALL cell%write()

print*, 'Initialized'
t = 0D0
! Time step loop
DO i = 1,cell%NT
!       Old way: now I read entire .bin file
!       Get velocity gradient of time step and normalize
        ! READ(7, *) cell%dU
        ! cell%dU = cell%dU*kdt

!       First, interpolate velocity gradient to current time step
!       Get time steps to interpolate to, first getting kolmogorov timesteps we're between
!! ============================
        kts = FLOOR(t/kfr)
        xs(1) = REAL(kts,8)*kfr
        xs(2) = xs(1) + kfr
        ys(1,:,:) = G(kts + 1,:,:,1)
        ys(2,:,:) = G(kts + 2,:,:,1)

!       Use two points to right and one left, unless you're at the end
        IF(t + 2D0*kfr .lt. nts*kfr) THEN
                xs(3) = xs(2) + kfr
                ys(3,:,:) = G(kts + 3,:,:,1)
        ELSE
                xs(3) = xs(1) - kfr
                ys(3,:,:) = G(kts - 1,:,:,1)
        ENDIF

!       Do interpolation, then normalize by kolm time
        cell%dU = QInterp(xs,ys,t)*kdt
!! ============================
        
        ! ! IF(i.gt.00) THEN
                ! cell%dU = 0d0
        ! ! ELSE
                ! cell%dU(1,3) = 1d0
        ! ! ENDIF

        ! IF(i.gt.00) THEN
                ! cell%dU = 0D0
        ! ENDIF

!       Get surface derivatives, then stress froom deformation, then motion from fluid
        CALL cell%derivs()
        CALL cell%stress()
        CALL cell%fluid()
        cell%umn(:,1) = 0D0
!       Initial volume
        IF(i.eq.1) THEN
                V0 = cell%Vol()
        ENDIF

!       Volume correction: small, inward normal velocity based off current volume/SA/time step
        zm = -(cell%Vol() - V0)/(3D0*cell%SA()*cell%dt)
        cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))

        umnt = cell%umn
        xmnt = cell%xmn

!       Volume reduction (add small inward normal vel every timestep)
        ! IF(cell%vol().gt. 4.22 .and. i.lt.500) umnt = umnt - 0.1D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
        ! IF(cell%vol().lt.4.185 .and. i.lt.500) umnt = umnt + 0.01D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
        ! IF(cell%vol().gt.4.1894 .and. i.lt.500) umnt = umnt - 0.01D0*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))

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
!         zm = -(cell%Vol() - V0)/(3D0*cell%SA()*cell%dt)
!         cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
!         umnt = 0.5D0*(umnt + cell%umn)
!         cell%xmn = xmnt + umnt*cell%dt
!         cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
!         cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
!         cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))

        
        cell%fab(1,:,:) = cell%J

!       Write some output
        IF((cell%cts .eq. 1) .or. (MOD(cell%cts,cell%dtinc)) .eq. 0) THEN
                CALL cell%write()
        ENDIF
        t = t + cell%dt
        write(*,'(I,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)'), &
        i, t, 2D0*MAXVAL((ABS(cell%ff)))/cell%B, MAXVAL((ABS(cell%umn))), cell%vol(), &
        cell%SA()!, cell%intg((2*cell%fab(1,:,:) - cell%C0)*(2*cell%fab(1,:,:) - cell%C0))/2D0
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN