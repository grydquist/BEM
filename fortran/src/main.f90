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

! The program still runs fine, but there are a few things that should be
! done to make it better/faster:

!!! Needs better parallelization. Real space sums should be able to be split
  ! into both target and integration surfaces instead of integration surfaces
  ! as they are now (i.e., Pcells(2,2)).
  ! Loops inside both Ewald and real can maybe be parallelized
!!! Instead of cycling in real space integral calculations, put the points into
  ! search boxes to only calculate nearest interactions
!!! Can maybe decrease memory in A2f by calculating the matrices one at a time
!!! mult cell not working at all
