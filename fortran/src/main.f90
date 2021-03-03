PROGRAM MAIN
USE UTILMOD
USE SHAPEMOD
IMPLICIT NONE
REAL(KIND = 8) :: tic, toc, t = 0D0, kdt, kfr
TYPE(cellType) :: cell
CHARACTER(:), ALLOCATABLE :: filein
INTEGER ::  i, argl, stat
! Tweezer stuff: chi = force/(G*a0), dSA = surface area of applied force
REAL(KIND = 8) :: chi, dSA


!! ================================================!!
!! ==                                            ==!!
!! ==           READ INFO BELOW ABOUT            ==!!
!! ==           HOW THIS BRANCH WORKS            ==!!
!! ==                                            ==!!
!! ================================================!!
! Instead of flow field being applied to the cell, a force is applied.
! As such, the capillary number, which relates a flow field to
! cell/fluid properties, is no longer applicable. We need a new
! non-dimensional number that relates the applied force to the cell
! properties. I call this chi, where
!
!       chi = F/(G*a0)
!
! F is the force applied to both ends of the cell. 
! I haven't put this in the input file, so it needs to be hardcoded
! in below! The Sinha/Graham paper has info on the experiments, but
! to match with them the formula is chi = F (in pN)/7.05 !


! Read in input file
CALL get_command_argument(number=1, length=argl)
ALLOCATE(character(argl) :: filein)
CALL get_command_argument(number=1, value=filein, status=stat)

! Initialize & allocate
cell = cellType(filein)
CALL cpu_time(tic)

print *, 'Reading in velocity gradient'
!! ============================
! OPEN(1,FILE = cell%gradfile, ACCESS = 'stream', ACTION = 'read')
! READ(1) nts, i
! ALLOCATE(G(nts,3,3,i), ys(3,3,3))
! READ(1) G
! CLOSE(1)
!! ============================

! Info about flow time scales/vel grad file
kdt = READ_GRINT_DOUB(filein, 'Kolm_time')
kfr = READ_GRINT_DOUB(filein, 'Kolm_frac')

!! ============================
! ! Let shape equilibrate by running a few time steps with no dU
! print *, 'Getting initial shape -'
! print *,  'Max velocity coefficient w.r.t. membrane time (want less than 0.0005):'
! cell%dU = 0D0

! ! To get the loop started
! cell%umn(1,1) = 1/cell%Ca
! DO WHILE(MAXVAL(ABS(cell%umn))*cell%Ca .gt. 0.005)
!         CALL cell%derivs()
!         CALL cell%stress() 
!         CALL cell%fluid()
!         cell%umn(:,1) = 0D0
!         cell%xmn = cell%xmn + cell%umn*cell%dt
!         cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
!         cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:)) 
!         cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
!         write(*,'(F8.5)') MAXVAL(ABS(cell%umn))*cell%Ca
! ENDDO
!! ============================

!   Write initial configuration
CALL cell%write()

print*, 'Initialized'

! New force parameter, chi:
chi = 325D0/7.05D0

! Time step loop
DO i = 1,cell%NT

!       First, interpolate velocity gradient to current time step
!       Get time steps to interpolate to, first getting kolmogorov timesteps we're between
!! ============================
!         kts = FLOOR(t/kfr)
!         xs(1) = REAL(kts,8)*kfr
!         xs(2) = xs(1) + kfr
!         ys(1,:,:) = G(kts + 1,:,:,1)
!         ys(2,:,:) = G(kts + 2,:,:,1)

! !       Use two points to right and one left, unless you're at the end
!         IF(t + 2D0*kfr .lt. nts*kfr) THEN
!                 xs(3) = xs(2) + kfr
!                 ys(3,:,:) = G(kts + 3,:,:,1)
!         ELSE
!                 xs(3) = xs(1) - kfr
!                 ys(3,:,:) = G(kts - 1,:,:,1)
!         ENDIF

! !       Do interpolation, then normalize by kolm time
!         cell%dU = QInterp(xs,ys,t)*kdt
!! ============================
        
!       Hardcoded shear flow
        cell%dU = 0D0
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
        t = t + cell%dt
        write(*,'(I4,X,F8.4,X,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') &
        i, t, 1D0*2D0*MAXVAL((ABS(cell%ff)))/cell%B, MAXVAL((ABS(cell%umn))), cell%vol(), &
        cell%SA()

!       Check if there's any funny business
        IF(isNaN(MAXVAL((ABS(cell%umn)))) .or. MAXVAL((ABS(cell%umn))) .gt. HUGE(t)) THEN
                print *, 'ERROR: inftys or NaNs'
                stop
        ENDIF

!       To keep a constant tweezer force, we need to do some fancy non-dimensionalizing
        
!       Magnitude of tweezer force (note that ff in shape is not non-dimm'd, it's mltiplied by G.
!       However, G is 1 for Ca = 1 in these sims so I guess it is essentially non-dimm'd ¯\_(ツ)_/¯)
!       ff = chi*G/dSA = chi*B/2/dSA
        dSA = cell%intg(REAL(ABS(cell%mark), 8))
        cell%frcmag = chi*cell%B/2D0/dSA

! !       If you want to see where the markers are applied
!         OPEN (UNIT = 88, FILE = 'tmpmark')
!         WRITE(88,*) cell%mark !cell%Yf%backward(cell%fmn(1,:)) !
!         CLOSE(88)
!         OPEN (UNIT = 88, FILE = 'tmpxf')
!         WRITE(88,*) cell%xf(1,:,:) !cell%x(1,:,:) !
!         WRITE(88,*) cell%xf(2,:,:) !cell%x(2,:,:) !
!         WRITE(88,*) cell%xf(3,:,:) !cell%x(3,:,:) !
!         CLOSE(88)

! !       If you want to see the percentage of surface area covered by the markers
!         print *, cell%intg(REAL(ABS(cell%mark), 8)), cell%b/2
!         stop
ENDDO

CALL cpu_time(toc)

print *, toc-tic

END PROGRAM MAIN