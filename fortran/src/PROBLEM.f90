MODULE PROBMOD
USE SHAPEMOD
IMPLICIT NONE

!==============================================================================!
!              The purpose of this module is perform calculations              !
!            that solve the overall math problem, including updating           !
!               the timestep, performing matrix inversions, and                !
!           facilitating a loop to calculate all the cells' responses          !
!==============================================================================!

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! Contains all the miscellaneous info about the problem
TYPE probType
    INTEGER :: cts, NT, dtinc, NCell
    REAL(KIND = 8) :: t
!   Velocity gradient file location
    REAL(KIND = 8) :: dU(3,3)
    CHARACTER(:), ALLOCATABLE :: gradfile

!   The location of the shared container
    TYPE(sharedType), POINTER :: info

!   Pointer to the cells
    TYPE(cellType), POINTER :: cell(:)

!   MPI stuff
    TYPE(cmType), POINTER :: cm
    INTEGER :: PCells(2)

!   Are we continuing from somewhere?
    LOGICAL :: cont = .false.

    CONTAINS
    PROCEDURE :: Update  => UpdateProb
    PROCEDURE :: Write   => WriteProb
    PROCEDURE :: Output  => OutputProb
    PROCEDURE :: Continue=> ContinueProb
END TYPE probType

! -------------------------------------------------------------------------!
INTERFACE probType
    PROCEDURE :: newprob
END INTERFACE probType

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

! Constructs the problem information
FUNCTION newprob(filein, reduce, cm, info) RESULT(prob)
    CHARACTER(len = *), INTENT(IN) :: filein
    LOGICAL, INTENT(IN) :: reduce
    TYPE(cmType), INTENT(IN), TARGET :: cm
    TYPE(sharedType), TARGET :: info
    TYPE(probType), TARGET :: prob
    TYPE(cellType) :: celltmp
    CHARACTER(:), ALLOCATABLE :: fileout
    CHARACTER(len = 3) :: restart, cont
    CHARACTER(len = 30) :: restfile, icC, contfile, cfile2
    REAL(KIND = 8), ALLOCATABLE :: cPt(:,:), ths(:,:), phs(:,:), thts(:), phis(:), xs(:), wg(:)
    REAL(KIND = 8) dphi
    INTEGER :: nt, np, ntf, npf, fali, p, m, ind, n, it, im2, im, ic, stat

    prob%cm => cm
    prob%info => info

!   General problem parameters
    CALL READ_MFS(prob%NT, filein, 'Max_time_steps')
    CALL READ_MFS(prob%NCell, filein, 'Number_cells')
    CALL READ_MFS(prob%dtinc, filein, 'Time_inc')
    CALL READ_MFS(prob%gradfile, filein, 'Gradient_file')
    CALL READ_MFS(fileout, filein, 'Output')
    prob%cts = 0
    prob%t = 0D0

!   Check if there are more processors than cells (can't handle)
    IF(prob%cm%np() .gt. prob%NCell) THEN
        print *, 'ERROR: More processors than cells'
        STOP
    ENDIF
    
!   Should we continue from a spot where we left off? !! Not working right now,
    !! doesn't return the same values for some reason
    IF(cont .eq. "Yes") prob%cont = .true.

!   Make harmonics(order, # of derivs, if we'll rotate or not)
    nt  = info%Y%nt
    np  = info%Y%np
    ntf = info%Yf%nt
    npf = info%Yf%np
    
    ALLOCATE(prob%cell(prob%NCell))

!   Make a master cell of sorts
    celltmp = cellType(filein, reduce, info)

!   Loop and make individual cells
    DO ic = 1, prob%NCell
        prob%cell(ic) = celltmp
        prob%cell(ic)%id = ic
        write(icC, "(I0.1)") ic
        prob%cell(ic)%fileout = TRIM(fileout//icC)
        
        !! Test                                        !!!!!
        ! SELECT CASE(ic)
        ! CASE(1)
        !     cell(ic)%xmn(3,1) =  0D0!  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = -5D0*SQRT(pi)*3D0/5D0
        ! CASE(2)
        !     cell(ic)%xmn(3,1) = 0D0!  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = 5D0*SQRT(pi)*3D0/5D0
        ! CASE(3)
        !     cell(ic)%xmn(3,1) =  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = -5D0*SQRT(pi)*3D0/2D0
        ! CASE(4)
        !     cell(ic)%xmn(3,1) =  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  0D0! -5D0*SQRT(pi)*3D0/3D0
        ! CASE(5)
        !     cell(ic)%xmn(3,1) =  5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  5D0*SQRT(pi)*3D0/2D0
        ! CASE(6)
        !     cell(ic)%xmn(3,1) = -5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) = -5D0*SQRT(pi)*3D0/2D0
        ! CASE(7)
        !     cell(ic)%xmn(3,1) = -5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  0D0! -5D0*SQRT(pi)*3D0/3D0
        ! CASE(8)
        !     cell(ic)%xmn(3,1) = -5D0*SQRT(pi)*3D0/5D0
        !     cell(ic)%xmn(1,1) =  5D0*SQRT(pi)*3D0/2D0
        ! END SELECT
    ENDDO
    
!   Which cells will the current processor handle? Integer div. rounds down
    n = MOD(prob%NCell, prob%cm%np())
    m = prob%NCell/prob%cm%np()

!   Get the top and bottom cell indices for this processor,
!   giving the remainder to the first processors
    IF (prob%cm%id() .lt. n) THEN
        prob%PCells(1) = (prob%cm%id()    )*(m + 1) + 1
        prob%PCells(2) = (prob%cm%id() + 1)*(m + 1)
    ELSE
        prob%PCells(1) = n*(m + 1) + (prob%cm%id()     - n)*m + 1
        prob%PCells(2) = n*(m + 1) + (prob%cm%id() + 1 - n)*m
    ENDIF
END FUNCTION newprob

! -------------------------------------------------------------------------!
! Writes xmn to a text file. Could be better
SUBROUTINE WriteProb(prob)
    TYPE(cellType), POINTER :: cell
    CLASS(probType), INTENT(IN) :: prob
    CHARACTER (LEN = 25) ctsst, datdir, filename
    INTEGER ic
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:)

    IF(prob%cm%slv()) RETURN
!   Don't write if it's not a timestep to write
    IF(.not.((prob%cts .eq. 1) .or. &
         (MOD(prob%cts,prob%dtinc)) .eq. 0)) RETURN

    DO ic = 1,prob%NCell
    cell => prob%cell(ic)
    IF(NOT(ALLOCATED(fmn))) THEN
        ALLOCATE(fmn(3,(cell%info%q+1)*(cell%info%q+1)))
    ENDIF

!   Formatting pain
    write(ctsst, "(I0.5)") prob%cts
    datdir = TRIM('dat/'//cell%fileout//'/')
    filename = TRIM('x_'//ctsst)

!   Write position - first half - real, x,y,z groups, second half imag
    IF(.not. cell%writ) THEN
        CALL MAKEDIRQQ(datdir)
    ENDIF
    OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    WRITE(88,*) REAL(cell%xmn)
    WRITE(88,*) AIMAG(cell%xmn)
    CLOSE(88)

!   Write velocity
    filename = TRIM('u_'//ctsst)
    OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    WRITE(88,*) REAL(cell%umn)!REAL(cell%fmn(:,1:((cell%p+1)*(cell%p+1)))) !
    WRITE(88,*) AIMAG(cell%umn)!AIMAG(cell%fmn(:,1:((cell%p+1)*(cell%p+1))))!
    CLOSE(88)

!   Write force
    ! fmn(1,:) = cell%Yf%forward(cell%ff(1,:,:), cell%q)
    ! fmn(2,:) = cell%Yf%forward(cell%ff(2,:,:), cell%q)
    ! fmn(3,:) = cell%Yf%forward(cell%ff(3,:,:), cell%q)

    ! filename = TRIM('f_'//ctsst)
    ! OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    ! WRITE(88,*) REAL(fmn(:,1:((cell%p+1)*(cell%p+1))))
    ! WRITE(88,*) AIMAG(fmn(:,1:((cell%p+1)*(cell%p+1))))
    ! CLOSE(88)

!   Another file that just has all the material constants for the simulation
    IF(.not. cell%writ) THEN
        filename = 'Params'
        OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
        WRITE(88,*) "p"
        WRITE(88,*) cell%info%p
        WRITE(88,*) "dt"
        WRITE(88,*) prob%info%dt
        WRITE(88,*) "dt_inc"
        WRITE(88,*) prob%dtinc
        WRITE(88,*) "Ed"
        WRITE(88,*) cell%C/cell%B
        WRITE(88,*) "Eb"
        WRITE(88,*) cell%Eb/2D0/cell%B
        WRITE(88,*) "Ca"
        WRITE(88,*) cell%Ca
        WRITE(88,*) "c0"
        WRITE(88,*) cell%c0
        WRITE(88,*) "lambda"
        WRITE(88,*) cell%lam
        CLOSE(88)
    ENDIF
    
!   Not optimal, but just a file for the max timestep
    filename = 'maxdt'
    OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    WRITE(88,*) prob%cts
    CLOSE(88)
    
    IF(.not. cell%writ) cell%writ = .true.
    ENDDO

END SUBROUTINE WriteProb

! -------------------------------------------------------------------------!
! Problem advancement
SUBROUTINE UpdateProb(prob, ord, reduce)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell, celli
    INTEGER, INTENT(IN) :: ord
    LOGICAL, INTENT(IN) :: reduce
    REAL(KIND = 8) :: zm
    COMPLEX(KIND = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:), ut(:), wrk(:), &
                                      A2(:,:), b2(:), A(:,:), b(:)
    INTEGER :: ic, ic2, i, row, col, iter, p_info
    REAL(KIND = 8), ALLOCATABLE :: rwrk(:), utr(:), uti(:)
    COMPLEX, ALLOCATABLE :: swrk(:)
    COMPLEX(KIND = 8), POINTER :: xmn(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER(KIND = 8) tic, toc, rate
    
    ALLOCATE(ut(prob%info%NmatT), &
             uti(prob%info%NmatT), &
             utr(prob%info%NmatT), &
             IPIV(prob%info%NmatT), wrk(prob%info%NmatT), &
             swrk(prob%info%NmatT*(prob%info%NmatT+1)), &
             rwrk(prob%info%NmatT), &
             A(prob%info%NmatT, prob%info%NmatT), &
             b(prob%info%NmatT))
    
    CALL SYSTEM_CLOCK(tic,rate)

    A = 0D0
    b = 0D0
    b2 = 0D0
    A2 = 0D0
!   Construct the matrix
!   First loop: velocity surface
    DO ic = prob%PCells(1), prob%PCells(2)
        row = (ic -1)*prob%info%Nmat + 1
        cell => prob%cell(ic)

!       Second loop: integral surface
        DO ic2 = 1,prob%NCell
            celli => prob%cell(ic2)
            
!           Get geometric info about new state, get stress in new state if first go round
            IF(ic.eq.prob%PCells(1)) THEN
                CALL celli%derivs()
                CALL celli%stress()
!               Characteristic grid spacing
                celli%h = SQRT(celli%SA()/celli%info%Y%nt**2D0)
            ENDIF

            col = (ic2-1)*prob%info%Nmat + 1
!           Get velocity sub-matrix for cell-cell combo (Same or diff cell)
            IF(ic .eq. ic2) THEN
                CALL cell%fluid(A2, b2)
            ELSE
                CALL cell%fluid(A2, b2, celli)
            ENDIF

!           Put sub matrix into big matrix
            A(row:row + prob%info%Nmat - 1, col:col + prob%info%Nmat - 1) = A2
!           Sum over all the integrals
            b(row:row + prob%info%Nmat - 1) = b(row:row + prob%info%Nmat - 1) + b2
        ENDDO
    ENDDO

    DO i = 1,prob%info%NmatT
        A(:,i) = prob%cm%reduce(REAL(A(:,i))) + prob%cm%reduce(AIMAG(A(:,i)))*ii
    ENDDO
    b = prob%cm%reduce(REAL(b)) + prob%cm%reduce(AIMAG(b))*ii

!   Invert big matrix to get a list of all the vel constants of all cells
    ut = 0D0
    IF(prob%cm%mas()) THEN
        CALL zcgesv(prob%info%NmatT, 1, A, prob%info%NmatT, IPIV, b, prob%info%NmatT, &
                    ut, prob%info%NmatT, wrk, swrk, rwrk, iter, p_info)
    ENDIF

!   Broadcast to all processors
    IF(prob%cm%np() .gt. 1) THEN
        utr = REAL(ut)
        uti = AIMAG(ut)
        CALL prob%cm%bcast(utr)
        CALL prob%cm%bcast(uti)
        ut = utr + uti*ii
    ENDIF

!   Advance in time now
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)

!       Reconstruct individual vels
        row = (ic-1)*prob%info%Nmat + 1
        cell%umn = 0D0
        cell%umn(1,1:prob%info%Nmat/3) = ut((/(i, i=row    , row + prob%info%Nmat-2, 3)/))
        cell%umn(2,1:prob%info%Nmat/3) = ut((/(i, i=row + 1, row + prob%info%Nmat-1, 3)/))
        cell%umn(3,1:prob%info%Nmat/3) = ut((/(i, i=row + 2, row + prob%info%Nmat  , 3)/))

!       Volume correction: small, inward normal velocity based off current volume/SA/time step
!       Removed for reduce, because that keeps things at constant volume
        if(.not. reduce) THEN
            zm = -(cell%Vol() - cell%V0)/(cell%SA()*prob%info%dt)
            cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
        ENDIF

        umnt = cell%umn
        xmnt = cell%xmn

!       Volume reduction (add small inward normal vel every timestep)
        IF(reduce) THEN
            IF(cell%Vol().gt. 4.22  ) umnt = umnt - 0.10D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
            IF(cell%Vol().lt. 4.185 ) umnt = umnt + 0.01D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
            IF(cell%Vol().gt. 4.1894) umnt = umnt - 0.01D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
        ENDIF
    
!       Second part for midpoint
        IF(ord .eq. 1) THEN
            ! CALL cell%derivs()
            ! CALL cell%stress()
            ! CALL cell%fluid(prob, A2, b2)
            ! zm = -(cell%Vol() - cell%V0)/(3D0*cell%SA()*info%dt)
            ! cell%umn = cell%umn + zm*cell%nkmn(:,1:((cell%p+1)*(cell%p+1)))
            ! umnt = 0.5D0*(umnt + cell%umn)
            ! cell%xmn = xmnt + umnt*info%dt
            ! cell%x(1,:,:) = cell%Y%backward(cell%xmn(1,:))
            ! cell%x(2,:,:) = cell%Y%backward(cell%xmn(2,:))
            ! cell%x(3,:,:) = cell%Y%backward(cell%xmn(3,:))
        ELSEIF(ord .ne. 1) THEN
            print*, "ERROR: time advancement of order >1 not supported"
            stop
        ENDIF
    
!       Update position and current time step
        cell%xmn = xmnt + umnt*prob%info%dt

! !       Simple periodic
!         IF(REAL(cell%xmn(1,1)).lt.-12D0*sqrt(pi)) THEN
!             cell%xmn(1,1) = cell%xmn(1,1) + 24D0*sqrt(PI)
!         ELSEIF(REAL(cell%xmn(1,1)).gt.12D0*sqrt(pi)) THEN
!             cell%xmn(1,1) = cell%xmn(1,1) - 24D0*sqrt(PI)
!         ENDIF

        cell%x(1,:,:) = cell%info%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = cell%info%Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = cell%info%Y%backward(cell%xmn(3,:))
    ENDDO

    CALL SYSTEM_CLOCK(toc)
    !print *, REAL(toc-tic)/REAL(rate)

    prob%cts = prob%cts + 1
    prob%t = prob%t + prob%info%dt

!   Check if there's any funny business
    IF(isNaN(MAXVAL(ABS(cell%umn))) .or. MAXVAL(ABS(cell%umn)) .gt. HUGE(zm)) THEN
            print *, 'ERROR: inftys or NaNs'
            STOP
    ENDIF
END SUBROUTINE UpdateProb

! -------------------------------------------------------------------------!
! Continue from where we left off at? 
!!! Should make this more robust (read Params, but at the top, etc.)
SUBROUTINE ContinueProb(prob)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell
    CHARACTER (LEN = 25) contfile, cfile2
    INTEGER ic, endt, stat, p, i, jmp
    REAL(KIND = 8), ALLOCATABLE :: xmnraw(:,:)
    LOGICAL :: exist_in

    DO ic = 1,prob%NCell
        cell => prob%cell(ic)
        INQUIRE(FILE = TRIM('dat/'//TRIM(cell%fileout)//'/maxdt'), EXIST = exist_in)
        IF(.not. exist_in) THEN
            print *, "ERROR: Files not found for continue. Have you run previous cases with the same number of cells?"
            stop
        ENDIF

!       Get the timestep of where we left off
        contfile = TRIM('dat/'//TRIM(cell%fileout))
        cfile2 = TRIM('dat/'//TRIM(cell%fileout)//'/maxdt')
        OPEN(unit = 13, file = TRIM(cfile2), action = 'read')
        READ(13, '(I16)', iostat = stat) endt
        CLOSE(13)
        prob%cts = endt

!       Open the file of where we left off and get the shape
        write(cfile2, "(I0.5)") endt
        cfile2 = TRIM('dat/'//TRIM(cell%fileout)//'/x_'//TRIM(cfile2))

!       Read in the files
        p = (cell%info%p + 1)*(cell%info%p + 1)*2
        ALLOCATE(xmnraw(3,p))
        OPEN(unit = 13, file = cfile2, action = 'read')
        DO i = 1,p
            READ(13, *, iostat = stat) xmnraw(:,i)
        ENDDO
        CLOSE(13)

!       Text file format: all real, then imag
        p = cell%info%p
        jmp = (p+1)*(p+1)
    
        cell%xmn = 0D0

        cell%xmn(1,1:(p+1)*(p+1)) = xmnraw(1,1:(p+1)*(p+1))
        cell%xmn(2,1:(p+1)*(p+1)) = xmnraw(2,1:(p+1)*(p+1))
        cell%xmn(3,1:(p+1)*(p+1)) = xmnraw(3,1:(p+1)*(p+1))

        cell%xmn(1,1:(p+1)*(p+1)) = cell%xmn(1,1:(p+1)*(p+1)) + xmnraw(1, jmp+1: 2*jmp)*ii
        cell%xmn(2,1:(p+1)*(p+1)) = cell%xmn(2,1:(p+1)*(p+1)) + xmnraw(2, jmp+1: 2*jmp)*ii
        cell%xmn(3,1:(p+1)*(p+1)) = cell%xmn(3,1:(p+1)*(p+1)) + xmnraw(3, jmp+1: 2*jmp)*ii
        
    ENDDO
    prob%t = prob%cts*prob%info%dt
END SUBROUTINE ContinueProb

! -------------------------------------------------------------------------!
! Some output
SUBROUTINE OutputProb(prob)
    CLASS(probType), INTENT(IN) :: prob
    INTEGER ic
    IF(prob%cm%slv()) RETURN

!   What the output means
    IF(MOD(prob%cts,50).eq.1) THEN
        PRINT *, "  CTS  Time    Cell  Max F     Max umn  Vol     SA"
        PRINT *, "----------------------------------------------------"
    ENDIF

!   The output
    DO ic = 1, prob%NCell
        write(*,'(I5,X,F8.4,X,I5, X, F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') & 
        prob%cts, prob%t, ic, 1D0*2D0*MAXVAL(ABS(prob%cell(ic)%ff))/prob%cell(ic)%B, &
        MAXVAL(ABS(prob%cell(ic)%umn)), prob%cell(ic)%vol(), prob%cell(ic)%SA()
    ENDDO
END SUBROUTINE OutputProb


END MODULE PROBMOD