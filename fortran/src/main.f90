PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, t = 0D0, xs(3), kdt, kfr, Gfac
TYPE(cellType), ALLOCATABLE :: cell(:)
TYPE(probType) :: prob
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, kts, nts, pthline
REAL(KIND=8), ALLOCATABLE :: G(:,:,:), ys(:,:,:), Gtmp(:,:,:,:)

print *, 'Reading in input file...'
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

print *, 'Initializing cell/harmonics...'
cell = cellType(filein, .false., prob)
CALL cpu_time(tic)

! Info about flow time scales/vel grad file
kdt = READ_GRINT_DOUB(filein, 'Kolm_time')
kfr = READ_GRINT_DOUB(filein, 'Kolm_frac')
pthline = READ_GRINT_INT(filein, 'Path_line')

print *, 'Reading in velocity gradient...'
OPEN(1,FILE = prob%gradfile, ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
! Read into a temporary array so we don't hold onto this big one
ALLOCATE(Gtmp(nts,3,3,i), ys(3,3,3))
READ(1) Gtmp
CLOSE(1)
! How many timesteps from G do we actually need?
Gfac = nts*kfr/(prob%NT*prob%dt)
nts = CEILING(nts/Gfac)
! To prevent only having 2 timesteps
if(nts .eq. 1) THEN
        ALLOCATE(G(3, 3, 3))
        G = Gtmp(1:3, :, :, pthline)
ELSE
        ALLOCATE(G(nts+1, 3, 3))
        G = Gtmp(1:nts+1, :, :, pthline)
ENDIF
DEALLOCATE(Gtmp)

! Let shape equilibrate by running a few time steps with no dU
print *, 'Getting initial shape -'
print *,  'Max velocity coefficient w.r.t. membrane time (want less than 0.0005*Ca):'
prob%dU = 0D0

!! ============================
cell(1)%umn(1,1) = 1/cell(1)%Ca
DO WHILE(MAXVAL(ABS(cell(1)%umn))*cell(1)%Ca .gt. 0.005)
        CALL cell(1)%derivs()
        CALL cell(1)%stress() 
        CALL cell(1)%fluid(prob)
        cell(1)%umn(:,1) = 0D0
        cell(1)%xmn = cell(1)%xmn + cell(1)%umn*prob%dt
        cell(1)%x(1,:,:) = cell(1)%Y%backward(cell(1)%xmn(1,:))
        cell(1)%x(2,:,:) = cell(1)%Y%backward(cell(1)%xmn(2,:)) 
        cell(1)%x(3,:,:) = cell(1)%Y%backward(cell(1)%xmn(3,:))
        write(*,'(F8.6)') MAXVAL(ABS(cell(1)%umn))*cell(1)%Ca
ENDDO
!! ============================

!   Write initial configuration
CALL cell(1)%write(prob)

print*, 'Initialized!'

! Time step loop
DO i = 1,prob%NT

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
        prob%dU = QInterp(xs,ys,t)*kdt
!! ============================
        
!       Hardcoded shear flow
        ! cell%dU = 0D0
        ! cell%dU(1,3) = 1d0

!       Get surface derivatives, then stress froom deformation, then motion from fluid
        CALL cell(1)%derivs()
        CALL cell(1)%stress()
        CALL cell(1)%fluid(prob)

!       Initial volume
        IF(i.eq.1) THEN
                cell(1)%V0 = cell(1)%Vol()
        ENDIF

!       Time advancer. Arguments are order of accuracy
!       and if we should do volume reduction routine.
        CALL cell(1)%update(prob, 1, .false.)

!       Write and display some output
        IF((prob%cts .eq. 1) .or. (MOD(prob%cts,prob%dtinc)) .eq. 0) THEN
                CALL cell(1)%write(prob)
        ENDIF

!       Check if there's any funny business
        IF(isNaN(MAXVAL((ABS(cell(1)%umn)))) .or. MAXVAL((ABS(cell(1)%umn))) .gt. HUGE(t)) THEN
                print *, 'ERROR: inftys or NaNs'
                stop
        ENDIF

        t = t + prob%dt
        write(*,'(I4,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') &
        i, t, 1D0*2D0*MAXVAL((ABS(cell(1)%ff)))/cell%B, MAXVAL((ABS(cell(1)%umn))), cell(1)%vol(), &
        cell(1)%SA()
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN