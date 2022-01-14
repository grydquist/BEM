PROGRAM MAIN
USE PROBMOD
IMPLICIT NONE
TYPE(probType) :: prob
TYPE(cmType), TARGET :: cmM
TYPE(sharedType), TARGET :: info
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat

! MPI communicator startup
CAll MPI_INIT(i)
cmM = newCm(MPI_COMM_WORLD)

IF(cmM%mas()) print *, 'Reading in input file...'
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

IF(cmM%mas()) print *, 'Initializing cell/harmonics...'
! Shared info about the problem
info = sharedType(filein)
! Problem handler
prob = probType(filein, .false., cmM, info)

IF(cmM%mas()) print*, 'Initialized!'

! Time step loop
DO i = 1,prob%NT

!       Updater
        CALL prob%update(1, .false.)

ENDDO

END PROGRAM MAIN