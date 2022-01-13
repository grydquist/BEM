PROGRAM MAIN
USE PROBMOD
IMPLICIT NONE
REAL(KIND = 8) :: kdt, kfr, Gfac
TYPE(probType) :: prob
TYPE(cmType), TARGET :: cmM
TYPE(sharedType), TARGET :: info
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, nts, pthline, ic
INTEGER(KIND = 8) rate, tic, toc
REAL(KIND = 8), ALLOCATABLE :: G(:,:,:), Gtmp(:,:,:,:)

! MPI communicator startup
CAll MPI_INIT(ic)
cmM = newCm(MPI_COMM_WORLD)

CALL SYSTEM_CLOCK(tic,rate)
IF(cmM%mas()) print *, 'Reading in input file...'
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

IF(cmM%mas()) print *, 'Initializing cell/harmonics...'
! Shared info about the problem
info = sharedType(filein)
prob = probType(filein, .false., cmM, info)

! Info about flow time scales/vel grad file
CALL READ_MFS(kdt, filein, 'Shear_dim_time')
CALL READ_MFS(kfr, filein, 'Gradient_timestep') !!! Check, think about more. Doesn't always end up like this, esp biologically
kfr = kfr/kdt ! Fraction of ts between velgrads
CALL READ_MFS(pthline, filein, 'Path_line')

IF(cm%mas()) print *, 'Reading in velocity gradient...'
OPEN(1,FILE = prob%gradfile, ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
! Read into a temporary array so we don't hold onto this big one
ALLOCATE(Gtmp(nts,3,3,i))
READ(1) Gtmp
CLOSE(1)

! How many timesteps from G do we actually need?
Gfac = nts*kfr/(prob%NT*prob%info%dt) ! Gfac: Ratio of total time in velgrad to total time requested
! If we have fewer velocity gradient time steps than requested, set down to gradient
IF(Gfac .lt. 1) THEN
        prob%NT = FLOOR(prob%NT*Gfac)
        print *, "Warning: Fewer gradient time steps than requested total time steps"
        print *, "Setting total timesteps to "
        print *, prob%NT
ELSE
        nts = CEILING(nts/Gfac)
ENDIF

! To prevent only having 2 timesteps
if(nts .eq. 1) THEN
        ALLOCATE(G(3, 3, 3))
        G = Gtmp(1:3, :, :, pthline)
ELSE
        ALLOCATE(G(nts+1, 3, 3))
        G = Gtmp(1:nts+1, :, :, pthline)
ENDIF
DEALLOCATE(Gtmp)

! Relax cell to get it into an equilibrium state
!! ============================
print *, 'Getting initial shape -'
print *,  'Max velocity coefficient w.r.t. membrane time (want less than 0.005*Ca):'
prob%dU = 0D0
DO ic = 1, prob%NCell
        CALL prob%cell(ic)%relax(0.005D0)
ENDDO
!! ============================

! Get initial volumes
DO ic = 1, prob%NCell
        CALL prob%cell(ic)%derivs()
        CALL prob%cell(ic)%stress()
        prob%cell(ic)%V0 = prob%cell(ic)%Vol()
ENDDO

! Should we continue from where we left off
!!! Not working properly right now (doesn't return same result when restarting)
IF(prob%cont) CALL prob%Continue()

!   Write initial configuration
CALL prob%write()
IF(cm%mas()) print*, 'Initialized!'

! Time step loop
DO i = 1,prob%NT
!! ============================
!       Do interpolation to get current grad tensor, then normalize by kolm time
        prob%info%dU = VelInterp(G,prob%t,nts,kfr)*kdt
!       Hardcoded shear
        ! prob%dU = 0D0
        ! prob%dU(1,3) = 1D0
        ! prob%dU(1,1) =  1D0
        ! prob%dU(3,3) = -1D0
!! ============================

!       Updater
        CALL SYSTEM_CLOCK(tic)
        CALL prob%update(1, .false.)
        CALL SYSTEM_CLOCK(toc)
        ! IF(prob%cm%mas()) print *, REAL(toc - tic)/REAL(rate)

!       Write and display some output
        CALL prob%write()
        CALL prob%output()
ENDDO

END PROGRAM MAIN