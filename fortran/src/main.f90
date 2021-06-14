PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, t = 0D0, kdt, kfr, Gfac
TYPE(cellType), ALLOCATABLE, TARGET :: cell(:)
TYPE(probType) :: prob
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, nts, pthline, ic
REAL(KIND = 8), ALLOCATABLE :: G(:,:,:), Gtmp(:,:,:,:)

print *, 'Reading in input file...'
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

print *, 'Initializing cell/harmonics...'
cell = cellType(filein, .false., prob)
prob%cell => cell
CALL cpu_time(tic)

! Info about flow time scales/vel grad file
kdt = READ_GRINT_DOUB(filein, 'Kolm_time')
kfr = READ_GRINT_DOUB(filein, 'Kolm_frac')
pthline = READ_GRINT_INT(filein, 'Path_line')

print *, 'Reading in velocity gradient...'
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

!   Write initial configuration
CALL cell(1)%write(prob)

print*, 'Initialized!'

! Time step loop
DO i = 1,prob%NT
!! ============================
!       Do interpolation to get current grad tensor, then normalize by kolm time
        prob%dU = VelInterp(G,t,nts,kfr)*kdt
!       Hardcoded shear
        prob%dU = 0D0
        prob%dU(1,3) = 1D0
!! ============================

!       Updater
        CALL CPU_TIME(tic)
        CALL prob%update(1, .false.)
        CALL CPU_TIME(toc)
        ! print *, toc - tic

!       Write and display some output
        IF((prob%cts .eq. 1) .or. (MOD(prob%cts,prob%dtinc)) .eq. 0) THEN
                CALL cell(1)%write(prob)
        ENDIF

        t = t + prob%dt
        DO ic = 1, prob%NCell
                write(*,'(I4,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') &
                prob%cts, t, 1D0*2D0*MAXVAL(ABS(cell(ic)%ff))/cell(ic)%B, MAXVAL(ABS(cell(ic)%umn)), cell(ic)%vol(), &
                cell(ic)%SA()
        ENDDO
ENDDO

CALL cpu_time(toc)
print *, toc-tic

END PROGRAM MAIN