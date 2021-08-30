PROGRAM MAIN
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: kdt, kfr, Gfac
TYPE(cellType), ALLOCATABLE, TARGET :: cell(:)
TYPE(probType) :: prob
TYPE(cmType), TARGET :: cm
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, nts, pthline, ic
INTEGER(KIND = 8) rate, tic, toc
REAL(KIND = 8), ALLOCATABLE :: G(:,:,:), Gtmp(:,:,:,:)

! MPI communicator startup
CAll MPI_INIT(ic)
CALL cm%new(MPI_COMM_WORLD)

CALL SYSTEM_CLOCK(tic,rate)
IF(cm%mas()) print *, 'Reading in input file...'
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

IF(cm%mas()) print *, 'Initializing cell/harmonics...'
prob%cm   => cm
cell = cellType(filein, .false., prob)
prob%cell => cell

! Info about flow time scales/vel grad file
kdt = READ_GRINT_DOUB(filein, 'Kolm_time')
kfr = READ_GRINT_DOUB(filein, 'Kolm_frac')
pthline = READ_GRINT_INT(filein, 'Path_line')

IF(cm%mas()) print *, 'Reading in velocity gradient...'
OPEN(1,FILE = prob%gradfile, ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
! Read into a temporary array so we don't hold onto this big one
ALLOCATE(Gtmp(nts,3,3,i))
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

! Relax cell to get it into an equilibrium state
!! ============================
print *, 'Getting initial shape -'
print *,  'Max velocity coefficient w.r.t. membrane time (want less than 0.0005*Ca):'
prob%dU = 0D0
DO ic = 1, prob%NCell
        CALL cell(ic)%relax(prob, 0.005D0)
ENDDO
!! ============================

! Get initial volumes
DO ic = 1, prob%NCell
        CALL cell(ic)%derivs()
        CALL cell(ic)%stress() 
        cell(ic)%V0 = cell(ic)%Vol()
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
        prob%dU = VelInterp(G,prob%t,nts,kfr)*kdt
!       Hardcoded shear
        prob%dU = 0D0
        prob%dU(1,3) = 1D0
        ! prob%dU(1,1) =  1D0
        ! prob%dU(3,3) = -1D0
!! ============================

!       Updater
        CALL SYSTEM_CLOCK(tic)
        CALL prob%update(1, .false.)
        CALL SYSTEM_CLOCK(toc)
        IF(prob%cm%mas()) print *, REAL(toc - tic)/REAL(rate)

!       Write and display some output
        CALL prob%write()
        CALL prob%output()
ENDDO

END PROGRAM MAIN