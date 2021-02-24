PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, t = 0D0, xs(3), kdt, kfr, Gfac
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, kts, nts
REAL(KIND=8), ALLOCATABLE :: G(:,:,:), ys(:,:,:), Gtmp(:,:,:,:)

print *, 'Reading in input file...'
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

print *, 'Initializing cell/harmonics...'
cell = cellType(filein)
CALL cpu_time(tic)

! Info about flow time scales/vel grad file
kdt = READ_GRINT_DOUB(filein, 'Kolm_time')
kfr = READ_GRINT_DOUB(filein, 'Kolm_frac')

print *, 'Reading in velocity gradient...'
OPEN(1,FILE = cell%gradfile, ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
! Read into a temporary array so we don't hold onto this big one
ALLOCATE(Gtmp(nts,3,3,i), ys(3,3,3))
READ(1) Gtmp
CLOSE(1)
! How many timesteps from G do we actually need?
Gfac = nts*kfr/(cell%NT*cell%dt)
nts = CEILING(nts/Gfac)
! To prevent only having 2 timesteps
if(nts .eq. 1) THEN
        ALLOCATE(G(3, 3, 3))
        G = Gtmp(1:3, :, :, 1)
ELSE
        ALLOCATE(G(nts+1, 3, 3))
        G = Gtmp(1:nts+1, :, :, 1)
ENDIF
DEALLOCATE(Gtmp)

! Let shape equilibrate by running a few time steps with no dU
print *, 'Getting initial shape -'
print *,  'Max velocity coefficient w.r.t. membrane time (want less than 0.0005*Ca):'
cell%dU = 0D0

!! ============================
cell%umn(1,1) = 1/cell%Ca
DO WHILE(MAXVAL(ABS(cell%umn))*cell%Ca .gt. 0.005)
        CALL cell%derivs()
        CALL cell%stress() 
        CALL cell%fluid()
        cell%umn(:,1) = 0D0
        cell%xmn = cell%xmn + cell%umn*cell%dt
        cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:)) 
        cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
        write(*,'(F8.5)') MAXVAL(ABS(cell%umn))*cell%Ca
ENDDO
!! ============================

!   Write initial configuration
CALL cell%write()

print*, 'Initialized!'

! Time step loop
DO i = 1,cell%NT

!       First, interpolate velocity gradient to current time step
!       Get time steps to interpolate to, first getting kolmogorov timesteps we're between
!! ============================
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

!       Do interpolation, then normalize by kolm time
        cell%dU = QInterp(xs,ys,t)*kdt
!! ============================
        
!       Hardcoded shear flow
        ! cell%dU = 0D0
        ! cell%dU(1,3) = 1d0

!       Get surface derivatives, then stress froom deformation, then motion from fluid
        CALL cell%derivs()
        CALL cell%stress()
        CALL cell%fluid()

!       Initial volume
        IF(i.eq.1) THEN
                cell%V0 = cell%Vol()
        ENDIF

!       Time advancer. Arguments are order of accuracy
!       and if we should do volume reduction routine.
        CALL cell%update(1, .false.)

!       Write and display some output
        IF((cell%cts .eq. 1) .or. (MOD(cell%cts,cell%dtinc)) .eq. 0) THEN
                CALL cell%write()
        ENDIF

!       Check if there's any funny business
        IF(isNaN(MAXVAL((ABS(cell%umn)))) .or. MAXVAL((ABS(cell%umn))) .gt. HUGE(t)) THEN
                print *, 'ERROR: inftys or NaNs'
                stop
        ENDIF

        t = t + cell%dt
        write(*,'(I4,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') &
        i, t, 1D0*2D0*MAXVAL((ABS(cell%ff)))/cell%B, MAXVAL((ABS(cell%umn))), cell%vol(), &
        cell%SA()
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN