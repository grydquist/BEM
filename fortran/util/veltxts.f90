PROGRAM MAIN
IMPLICIT NONE
INTEGER ::  i,  nts, j, k, z, pth
REAL(KIND=8), ALLOCATABLE ::  Gtmp(:,:,:,:)
REAL(KIND=8) :: G(3,3), t, ts, kfr, kdt
REAL(KIND=8) :: t1(143), J1(143), S1(143)
! This takes the bin file I was given and makes a bunch of individual txt files

print *, 'Reading in velocity gradient'
OPEN(1,FILE = 'VG.bin', ACCESS = 'stream', ACTION = 'read')
READ(1) nts, i
! Read into a temporary array so we don't hold onto this big one
ALLOCATE(Gtmp(nts,3,3,i))
READ(1) Gtmp
CLOSE(1)


OPEN(1,FILE = '../../Turbm20.txt')
DO z=1,143
READ(1,*) t1(z), J1(z), S1(z)
ENDDO
CLOSE(1)

pth=20
ts = 0.35D0
t = 0D0
kfr = 0.1D0
kdt = 2.0D-6
nts = 500
print *, "Writing"
OPEN(1,FILE = 'Turb20.txt', ACTION = 'write')
DO z = 1,142
G = VelInterp(Gtmp(:,:,:,pth),t,nts,kfr)*kdt

! Write the time/stress
WRITE(1,'(F5.3, X)', advance="no") t1(z)
WRITE(1,'(F8.5, X)', advance="no") J1(z)
WRITE(1,'(F8.5, X)', advance="no") S1(z)

DO j=1,3
    DO k =1,3
        IF(j.eq.3 .and. k.eq.3) THEN
            WRITE(1,'(F8.5, X)') G(j,k)!Gtmp(z,j,k,pth)
        ELSE
            WRITE(1,'(F8.5, X)', advance="no") G(j,k)!Gtmp(z,j,k,pth)
        ENDIF
    ENDDO
ENDDO

t = t + ts
ENDDO
CLOSE(1)

CONTAINS
! -------------------------------------------------------------------------!
! Uses above function to get current velocity gradient
FUNCTION VelInterp(G, t, nts, kfr) RESULT(dU)
    REAL(KIND = 8), INTENT(IN) :: t, kfr
    INTEGER, INTENT(IN) :: nts
    INTEGER :: kts
    REAL(KIND = 8), INTENT(IN) :: G(:,:,:)
    REAL(KIND = 8) :: dU(3,3)
    REAL(KIND = 8) :: xs(3)
    REAL(KIND = 8), ALLOCATABLE :: ys(:,:,:)

    ALLOCATE(ys(3,3,3))

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
    dU = ys(1,:,:)*(t - xs(2))*(t - xs(3))/((xs(1) - xs(2))*(xs(1) - xs(3))) &
       + ys(2,:,:)*(t - xs(1))*(t - xs(3))/((xs(2) - xs(1))*(xs(2) - xs(3))) &
       + ys(3,:,:)*(t - xs(1))*(t - xs(2))/((xs(3) - xs(1))*(xs(3) - xs(2)))
END FUNCTION VelInterp

END PROGRAM MAIN
