PROGRAM MAIN
IMPLICIT NONE
INTEGER ::  i,  nts
REAL(KIND=8), ALLOCATABLE ::  Gtmp(:,:,:,:)
! This takes the bin file I was given and just cuts down the size a bit

print *, 'Reading in velocity gradient'
OPEN(1,FILE = 'GU_1.bin', ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
! Read into a temporary array so we don't hold onto this big one
ALLOCATE(Gtmp(nts,3,3,i))
READ(1) Gtmp
CLOSE(1)

print *, "Writing"
OPEN(1,FILE = 'VG.bin', ACCESS = 'stream', ACTION = 'write')
WRITE(1) 1000, 20
WRITE(1) Gtmp(1:1000,:,:,1:20)
CLOSE(1)

END PROGRAM MAIN