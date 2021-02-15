PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, t = 0D0, xs(3), kdt, kfr
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat, kts, nts
REAL(KIND=8), ALLOCATABLE :: G(:,:,:,:), ys(:,:,:)

! Read in input file
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

! Initialize & allocate
cell = cellType(filein)
CALL cpu_time(tic)

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

!! ============================
! Let shape equilibrate by running a few time steps with no dU
print *, 'Getting initial shape -'
print *,  'Max velocity coefficient w.r.t. membrane time (want less than 0.0005):'
cell%dU = 0D0

! To get the loop started
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

print*, 'Initialized'

! Time step loop
DO i = 1,cell%NT

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
        
!       Hardcoded shear flow
        cell%dU = 0D0
        cell%dU(1,3) = 1d0

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
        t = t + cell%dt
        write(*,'(I4,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') &
        i, t, 2D0*MAXVAL((ABS(cell%ff)))/cell%B, MAXVAL((ABS(cell%umn))), cell%vol(), &
        cell%SA()
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN