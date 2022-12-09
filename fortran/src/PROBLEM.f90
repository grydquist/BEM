MODULE PROBMOD
USE SHAPEMOD
USE CMMOD
USE PERIODICMOD
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
    INTEGER :: cts, NT, dtinc, NCell, nts
    REAL(KIND = 8) :: t, kdt, kfr
    REAL(KIND = 8), ALLOCATABLE :: G(:,:,:)

!   The location of the shared container
    TYPE(sharedType), POINTER :: info

!   Pointer to the cells
    TYPE(cellType), POINTER :: cell(:)
    
!   Name of output file
    CHARACTER(:), ALLOCATABLE :: fileout

!   MPI stuff
    TYPE(cmType), POINTER :: cm
    INTEGER :: ictxt, desca(9), descb(9)
    INTEGER, ALLOCATABLE :: map_A(:,:), map_b(:)

    CONTAINS
    PROCEDURE :: Update  => UpdateProb
    PROCEDURE :: Write   => WriteProb
    PROCEDURE :: Output  => OutputProb
    PROCEDURE :: Continue=> ContinueProb
    PROCEDURE :: StrnRt  => StrnRtProb
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
    CHARACTER(:), ALLOCATABLE :: fileout, cont, gfile
    CHARACTER(len = 30) :: icC
    INTEGER :: m, n, ic, pthline, it, sqre, nprow, npcol, nb, p_inf, &
               myrow, mycol, np, nq, nqrhs, ictxt, tot_blocks, &
               xblx, yblx, itb, i, ib, j, jb, strti, strtj, nbx, nby
    REAL(KIND = 8) :: Gfac
    REAL(KIND = 8), ALLOCATABLE :: Gtmp(:,:,:,:)
    LOGICAL :: fl = .false.
    TYPE(cellprops) props

    prob%cm => cm
    prob%info => info

!   General problem parameters
    CALL READ_MFS(prob%NT, filein, 'Max_time_steps')
    CALL READ_MFS(prob%NCell, filein, 'Number_cells')
    CALL READ_MFS(prob%dtinc, filein, 'Time_inc')
    CALL READ_MFS(gfile, filein, 'Gradient_file')
    CALL READ_MFS(fileout, filein, 'Output')
    CALL READ_MFS(cont, filein, 'Continue')
    CALL READ_MFS(prob%kdt, filein, 'Shear_dim_time')
    CALL READ_MFS(prob%kfr, filein, 'Gradient_timestep')
    prob%kfr = prob%kfr/prob%kdt ! Fraction of ts between velgrads
    CALL READ_MFS(pthline, filein, 'Path_line')
    prob%cts = 0
    prob%t = 0D0

!   The rest of the parameters specific to the cell go here
    CALL READ_MFS(props%lam, filein, 'Viscosity_Ratio')
    CALL READ_MFS(props%Ca, filein, 'Capillary')
    CALL READ_MFS(props%C, filein, 'Dilatation_Ratio')
    CALL READ_MFS(props%Eb, filein, 'Bending_Modulus')
    CALL READ_MFS(props%c0, filein, 'Spont_Curvature')
    CALL READ_MFS(props%int_pres, filein, 'Internal_Pressure')

!   Initialize info now
    IF(prob%cm%mas()) print *, 'Done with init file!'
    CALL info%init()

    ALLOCATE(prob%cell(prob%NCell))

!   Make a master cell of sorts
    celltmp = cellType(filein, reduce, info, props)

!   Loop and make individual cells
    DO ic = 1, prob%NCell
        prob%cell(ic) = celltmp
        prob%cell(ic)%id = ic
        write(icC, "(I0.1)") ic
        prob%cell(ic)%fileout = TRIM('cell_'//icC)
    ENDDO
    prob%fileout=TRIM(fileout)
    
!   For parallel, we split the particles up based on theta
!   integration point coordinates

!   If we have fewer processors than cells, each proc must handle multiple cells
!   We only want each cell to be assigned to one proc in this case
!   Otherwise, we want each proc only assigned to one cell
!   PCells(1,:) are starting/ending cell indices,
!   PCells(2,:) are starting/ending theta points (in the cells)

!   Also need to find all the procs corresponding to a cell
    ALLOCATE(info%CProcs(prob%NCell, 2))
    info%CProcs = 0

!   First case: more procs than cells
    IF(prob%cm%np() .ge. prob%NCell) THEN
!       Number of procs in each cell rounded up
        m = prob%cm%np()/prob%NCell + 1
        n = MOD(prob%cm%np(), prob%NCell)
        it = 0
        ic = 0
!       To find the cell for a proc, count up through cell procs
        DO WHILE(it .lt. (prob%cm%id() + 1))
            IF(ic .eq. n) m = m - 1
            
!           Current top processor
            it = it + m
            ic = ic + 1
        ENDDO
        info%PCells(1,1) = ic
        info%PCells(1,2) = ic

!       Also use above to get top and bottom cells of processor
        info%CProcs(ic, 1) = it - m
        info%CProcs(ic, 2) = it - 1

!       Current position in this proc
        it = prob%cm%id() - (it - m)

!       Then allocate thetas in that proc with a similar process as above
        n = MOD(info%Y%nt, m)
        m = info%Y%nt/m
        
        IF (it .lt. n) THEN
            info%PCells(2, 1) = (it    )*(m + 1) + 1
            info%PCells(2, 2) = (it + 1)*(m + 1)
        ELSE
            info%PCells(2, 1) = n*(m + 1) + (it     - n)*m + 1
            info%PCells(2, 2) = n*(m + 1) + (it + 1 - n)*m
        ENDIF
    ELSE
!       Number of cells per proc rounded down
        m = prob%NCell/prob%cm%np()
        n = MOD(prob%NCell, prob%cm%np())
        IF (prob%cm%id() .lt. n) THEN
            info%PCells(1, 1) = (prob%cm%id()    )*(m + 1) + 1
            info%PCells(1, 2) = (prob%cm%id() + 1)*(m + 1)
        ELSE
            info%PCells(1, 1) = n*(m + 1) + (prob%cm%id()     - n)*m + 1
            info%PCells(1, 2) = n*(m + 1) + (prob%cm%id() + 1 - n)*m
        ENDIF

!       Do the whole range of thetas for each cell
        info%PCells(2,1) = 1
        info%PCells(2,2) = info%Y%nt

!       Need to deal with cells per procs
        DO ic = info%PCells(1, 1), info%PCells(1, 2)
            info%CProcs(ic, 1) = prob%cm%id()
            info%CProcs(ic, 2) = prob%cm%id()
        ENDDO
    ENDIF

!   Velocity gradient information
    IF(prob%cm%mas()) print *, 'Reading in velocity gradient...'
    OPEN(1,FILE = gfile, ACCESS = 'stream', ACTION = 'read')
    READ(1) prob%nts, ic
!   Read into a temporary array so we don't hold onto this big one
    ALLOCATE(Gtmp(prob%nts,3,3,ic))
    READ(1) Gtmp
    CLOSE(1)

!   Some pre-processing for the parallel linear solver
    IF(prob%cm%np().gt.1) THEN
        
!       We need to make a grid out of the processors: make as square as possible
        sqre = FLOOR(SQRT(REAL(prob%cm%np()))) + 1
        DO WHILE(MOD(prob%cm%np(), sqre) .ne. 0)
            sqre = sqre - 1
        ENDDO

!       Proc rows and columns for square grid
        nprow = MAX(sqre, prob%cm%np()/sqre)
        npcol = MIN(sqre, prob%cm%np()/sqre)

!       Oranizes data matrix into blocks: make same size as proc grid
        nb = (prob%cm%np() + nprow - 1) / nprow

!       Normal setup things
        CALL blacs_get( -1, 0, prob%ictxt )
        CALL blacs_gridinit( prob%ictxt, 'Row-major', nprow, npcol )
        CALL blacs_gridinfo( prob%ictxt, nprow, npcol, myrow, mycol )

!       number of points this proc is responsible for in x/y
        np    = numroc( info%NmatT, nb, myrow, nprow )
        nq    = numroc( info%NmatT, nb, mycol, npcol )
        nqrhs = numroc( 1, nb, mycol, npcol )

        CALL descinit( prob%desca, info%NmatT, info%NmatT, nb, nb, 0, 0, prob%ictxt, max( 1, np ),p_inf)
        CALL descinit( prob%descb, info%NmatT, 1,          nb, nb, 0, 0, prob%ictxt, max( 1, np ),p_inf)

!       Maps points in A and b to those needed by processor
        ALLOCATE(prob%map_A(2, np*nq))
        IF(nqrhs.gt.0) THEN
            ALLOCATE(prob%map_b(np))
        ENDIF

!       Total blocks in both directions
        tot_blocks = (info%NmatT + nb - 1)/nb

!       Number of blocks this processor gets in each dir
        xblx = tot_blocks/nprow
        yblx = tot_blocks/npcol

!       Distribute extra blocks
        IF(MOD(tot_blocks, nprow) .gt. (myrow)) xblx = xblx + 1
        IF(MOD(tot_blocks, npcol) .gt. (mycol)) yblx = yblx + 1

!       Get the coordinate mapping from the matrix we need to solve for the points we need to solve with this proc
        it = 0
        itb = 0

!       Loops have a very specific order for scalapack
        DO j = 1, yblx
!           Get the indices of the first point in this block
            strtj = (mycol)*nb + (j-1)*nb*npcol + 1
        
!           Check for partial block if this is the last block
            IF(mycol + (j-1)*npcol + 1 .eq. tot_blocks) THEN
                nby = info%NmatT - strtj + 1
            ELSE
                nby = nb
            ENDIF
        
!           Points in this block
            DO jb = 1,nby
        
!               Block loop
                DO i = 1, xblx
                
!                   Get the indices of the first point in this block
                    strti = (myrow)*nb + (i-1)*nb*nprow + 1
                
!                   Check for partial block if this is the last block
                    IF(myrow + (i-1)*nprow + 1 .eq. tot_blocks) THEN
                        nbx = info%NmatT - strti + 1
                    ELSE
                        nbx = nb
                    ENDIF
        
!                   Points in block again
                    DO ib = 1,nbx
                        it = it + 1
                        prob%map_A(1, it) = strti + ib - 1
                        prob%map_A(2, it) = strtj + jb - 1
                        
!                       Also distribute b here
                        IF(nqrhs.gt.0 .and. jb.eq.1 .and. j.eq.1) THEN
                            itb = itb + 1
                            prob%map_b(itb) = strti + ib - 1
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDIF

!   How many timesteps from G do we actually need?
    Gfac = prob%nts*prob%kfr/(prob%NT*prob%info%dt) ! Gfac: Ratio of total time in velgrad to total time requested
!   If we have fewer velocity gradient time steps than requested, set down to gradient
    IF(Gfac .le. 1) THEN
        prob%NT = FLOOR(prob%NT*Gfac)
        print *, "Warning: Fewer gradient time steps than requested total time steps"
        print *, "Setting total timesteps to "
        print *, prob%NT
        fl = .true.
    ELSE
        prob%nts = CEILING(prob%nts/Gfac)
    ENDIF
    
!   To prevent only having 2 timesteps
    IF(prob%nts .eq. 1) THEN
        ALLOCATE(prob%G(3, 3, 3))
        prob%G = Gtmp(1:3, :, :, pthline)
    ELSE
!       When we need to use the whole grad file
        IF(.not. fl) THEN
            ALLOCATE(prob%G(prob%nts+1, 3, 3))
            prob%G = Gtmp(1:prob%nts+1, :, :, pthline)
        ELSE
            ALLOCATE(prob%G(prob%nts, 3, 3))
            prob%G = Gtmp(1:prob%nts, :, :, pthline)
        ENDIF
    ENDIF

!   Inital support point array
    IF(prob%info%periodic) THEN
        CALL SuppPoints(prob%info)
    ENDIF

!   Relax cell to get it into an equilibrium state, get volumes
!!  ============================
    IF(prob%cm%mas()) print *, 'Getting initial shape -'
    IF(prob%cm%mas()) print *, 'Max velocity coefficient w.r.t. membrane time (want less than 0.005*Ca):'
    Gfac = (((4D0*PI/3D0)/((0.95D0**(3D0))*(PI/6D0)))**(1D0/3D0))/2D0
    DO ic = 1, prob%NCell
            CALL prob%cell(ic)%relax(0.005D0)!!!
            CALL prob%cell(ic)%derivs()
            CALL prob%cell(ic)%stress()
            prob%cell(ic)%V0 = prob%cell(ic)%Vol()
            !! Test                                        !!!!!
            SELECT CASE(MOD(ic, 4))!9) )
            CASE(1)
                prob%cell(ic)%xmn(1,1) = (Gfac)/(ispi*0.5D0)!0.50D0/(ispi*0.5D0)!
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!0.50D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (Gfac)/(ispi*0.5D0)!0.50D0/(ispi*0.5D0)!
            CASE(3)
                prob%cell(ic)%xmn(1,1) = (Gfac)/(ispi*0.5D0) !+ info%bvl/(ispi*0.5D0)!3D0/(ispi*0.5D0)! 
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!0.5D0/(ispi*0.5D0)!
            CASE(2)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)!
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!0.5D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!1.5D0/(ispi*0.5D0)!
            CASE(0)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!0.5D0/(ispi*0.5D0)!
                prob%cell(ic)%xmn(2,1) = (-1+2*(ic/4))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)!(1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!
                prob%cell(ic)%xmn(3,1) = (Gfac)/(ispi*0.5D0)!1.5D0/(ispi*0.5D0)!
            CASE(5)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!3.5D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!2.75D0/(ispi*0.5D0)
            CASE(6)
                prob%cell(ic)%xmn(1,1) = (5D0*Gfac)/(ispi*0.5D0)!4.5D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!3D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (3D0*Gfac)/(ispi*0.5D0)!4.25D0/(ispi*0.5D0)
            CASE(7)
                prob%cell(ic)%xmn(1,1) = (Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!0.1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
            CASE(8)
                prob%cell(ic)%xmn(1,1) = (3D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!0.1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
            CASE(4)
                prob%cell(ic)%xmn(1,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1+2*(ic/9))*(Gfac)/(ispi*0.5D0)!0.1D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (5D0*Gfac)/(ispi*0.5D0)!1D0/(ispi*0.5D0)
            END SELECT
            prob%cell(ic)%x(1,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(1,:))
            prob%cell(ic)%x(2,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(2,:))
            prob%cell(ic)%x(3,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(3,:))
    ENDDO
!    ============================
    
!   Should we continue from a spot where we left off? !! Not working right now,
    !! doesn't return the same values for some reason
    IF(cont .eq. "Yes") CALL prob%Continue()
    
!   Write initial configuration
    CALL prob%write()

END FUNCTION newprob

! -------------------------------------------------------------------------!
! Writes xmn to a text file. Could be better
SUBROUTINE WriteProb(prob)
    TYPE(cellType), POINTER :: cell
    CLASS(probType), INTENT(IN) :: prob
    CHARACTER (LEN = 45) ctsst, datdir, filename
    INTEGER ic
    COMPLEX(KIND = 8), ALLOCATABLE :: fmn(:,:)

    IF(prob%cm%slv()) RETURN
!   Don't write if it's not a timestep to write
    IF(.not.((prob%cts .eq. 1) .or. &
         (MOD(prob%cts,prob%dtinc)) .eq. 0)) RETURN

!   Write top-level directory containing all lower level directoriesx
    IF(.not. prob%cell(1)%writ) THEN
        datdir = TRIM('dat/'//prob%fileout//'/')
        CALL MAKEDIRQQ(datdir)
    ENDIF

    DO ic = 1,prob%NCell
    cell => prob%cell(ic)
    IF(NOT(ALLOCATED(fmn))) THEN
        ALLOCATE(fmn(3,(cell%info%q+1)*(cell%info%q+1)))
    ENDIF

!   Formatting pain
    write(ctsst, "(I0.5)") prob%cts
    datdir = TRIM('dat/'//prob%fileout//'/'//cell%fileout//'/')
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
    TYPE(YType), POINTER :: Y
    TYPE(sharedType), POINTER :: info
    REAL(KIND = 8) :: zm
    COMPLEX(KIND = 8), ALLOCATABLE :: umnt(:,:), xmnt(:,:), ut(:), wrk(:), &
                                      A2(:,:), b2(:), A(:,:), b(:), &
                                      fmnR(:,:), xmnR(:,:), &
                                      fv1(:), fv2(:), fv3(:), fv(:,:), &
                                      b2f(:), b2r(:), YVt(:), um(:,:,:,:,:), &
                                      A2f(:,:,:,:,:,:,:,:), A2r(:,:,:,:,:,:), &
                                      A2t(:,:,:,:,:,:), b2t(:), &
                                      v(:,:,:), Ap(:), bp(:)
    INTEGER :: ic, ic2, i, j, row, col, iter, p_info, m, n, im, im2, &
               it, th_st, th_end, tot_gs, curid, inds_proc(2)
    REAL(KIND = 8), ALLOCATABLE :: rwrk(:), utr(:), uti(:), nJt(:,:,:), xv(:,:)
    COMPLEX, ALLOCATABLE :: swrk(:)
    INTEGER, ALLOCATABLE :: IPIV(:), inds(:)
    INTEGER(KIND = 8) tic, toc, rate, tic0, toc0

    Y => prob%info%Y
    info => prob%info

!   If periodic, make variables for the Fourier part
    IF(info%periodic) THEN
        ALLOCATE( &
        fv(3,Y%nt*Y%np*prob%NCell), &
        nJt(3, Y%nt, Y%np), &
        um(3,3,Y%nt,Y%np, prob%NCell), &
        fmnR(3, (info%p+1)*(info%p+1)), &
        xmnR(3, (info%p+1)*(info%p+1)))
    ENDIF


!   Do interpolation to get current grad tensor, then normalize by kolm time
    info%dU = VelInterp(prob%G,prob%t,prob%nts,prob%kfr)*prob%kdt
!! ============================
!   Hardcoded shear
    info%dU = 0D0 !!!
    info%dU(1,3) = 1D0
    ! prob%dU(1,1) =  1D0
    ! prob%dU(3,3) = -1D0
!! ============================

    CALL SYSTEM_CLOCK(tic,rate)
!   First, calculate geometric quantities and stresses in all of the cells
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)

!       For periodic, check if cells are outside primary cell, and put them in if so
        IF(info%periodic) THEN
            CALL cell%inPrim()
        ENDIF

        CALL cell%derivs()
        CALL cell%stress()

!       Construct a large vector containing all the forces/locations in
!       preparation for the fast Ewald summation
        IF(info%periodic) THEN
!           Forces and locations
            fmnR = cell%fmn(:,1:(info%p+1)*(info%p+1))
            xmnR = cell%xmn(:,1:(info%p+1)*(info%p+1))
!           Append vectors
            CALL EwaldPreCalc(fv1, Y, xv, fmn=fmnR(1,:), xmn=xmnR)
            CALL EwaldPreCalc(fv2, Y, fmn=fmnR(2,:))
            CALL EwaldPreCalc(fv3, Y, fmn=fmnR(3,:))
        ENDIF
    ENDDO
    ! CALL prob%cm%barrier()
    CALL SYSTEM_CLOCK(toc)
    IF(prob%cm%mas()) print *, "Stress/derivs: ", REAL(toc-tic)/REAL(rate)

!   Now we need to calculate the SpH velocity coefficients of the cells.
!   For periodic, we calculate LHS matrix and RHS vector in two parts:
!   Long range (via FFTs) and short range (directly, truncated at low r).
!   We calculate these two matrices individually at integration points,
!   add them, and then perform the Galerkin transfo on these summed matrices.

    IF(info%periodic) THEN
        CALL SYSTEM_CLOCK(tic,rate)

!       Configure parallel: Each cell has it's own grid to be sent to the grid for each n,m
!       Split up the total resources among processors, done here
        tot_gs = prob%NCell*(Y%p*Y%p + Y%p)*0.5D0

        ALLOCATE(inds(prob%cm%np()))
        inds = ParallelSplit(tot_gs, prob%cm%np())
        inds_proc(1) = inds(prob%cm%id() + 1)

        IF(prob%cm%id() .ne. prob%cm%np() - 1) THEN
            inds_proc(2) = inds(prob%cm%id() + 2) - 1
        ELSE
            inds_proc(2) = tot_gs
        ENDIF
        ALLOCATE(b2f(3*Y%nt*Y%np*prob%NCell))

!       Calculate RHS vector
        fv(1,:) = fv1 
        fv(2,:) = fv2
        fv(3,:) = fv3
        CALL EwaldG(info=info, x0=xv, f=fv, u1=b2f, full=.true.)

        ALLOCATE(A2f(3,3, 2*(Y%p-1)+1, Y%p, Y%nt, Y%np, prob%NCell, prob%NCell)) !!! Can really clean up memory here
        A2f = 0D0
        it = 0

!       Start with the long range part, as this calculates all points for one integral surface
!       The "forces" here are the spherical harmonics multiplied by the normal vector
        DO ic = 1, prob%NCell
            cell => prob%cell(ic)

!           We need to make a vector for all of the spherical harmonics,
!           multiplied by normal vector, for this integration surface
            fmnR(1,:) = cell%nkt(1,1:(info%p+1)*(info%p+1))
            fmnR(2,:) = cell%nkt(2,1:(info%p+1)*(info%p+1))
            fmnR(3,:) = cell%nkt(3,1:(info%p+1)*(info%p+1))

            nJt(1,:,:) = Y%backward(fmnR(1,:), info%p)
            nJt(2,:,:) = Y%backward(fmnR(2,:), info%p)
            nJt(3,:,:) = Y%backward(fmnR(3,:), info%p)

            DO n = 0, Y%p-1
                ALLOCATE(v(n+1, 3, Y%nt*Y%np))
                im = 0
                DO m = -n,n
                    im = im + 1
                    im2 = m + (Y%p)

                    IF(m .gt. 0) THEN ! Symmetry
                        ! A2f(:,:,im2,n+1,:,:,:, ic) = CONJG(A2f(:,:,im2 - 2*m,n+1,:,:,:, ic))*(-1D0)**m
                        CYCLE
                    ENDIF

!                   Helps with indexing, but could be made slightly more efficient
                    it = it + 1

!                   Cycle if we're not doing this grid in this proc
                    IF( (it .lt. inds_proc(1)) .or. (it .gt. inds_proc(2)) ) CYCLE

                    CALL EwaldPreCalc(YVt, Y, fin=Y%nm(n+1)%v(im,:,:)*nJt(1,:,:))
                    v(im,1,:) = YVt
                    DEALLOCATE(YVt)

                    CALL EwaldPreCalc(YVt, Y, fin=Y%nm(n+1)%v(im,:,:)*nJt(2,:,:))
                    v(im,2,:) = YVt
                    DEALLOCATE(YVt)

                    CALL EwaldPreCalc(YVt, Y, fin=Y%nm(n+1)%v(im,:,:)*nJt(3,:,:))
                    v(im,3,:) = YVt
                    DEALLOCATE(YVt)

!                   Do only the grid spreading this proc is responsible for
                    CALL EwaldT(info=info, x0=xv, f=v(im,:,:), strt=ic, um=um, full=.true.)
                    A2f(:,:, im2, n+1,:,:,:, ic) = um
                ENDDO
                DEALLOCATE(v)
            ENDDO
        ENDDO
        DEALLOCATE(um)

        ! CALL prob%cm%barrier()
        CALL SYSTEM_CLOCK(toc)
        IF(prob%cm%mas()) print *, "Ewald: ", REAL(toc-tic)/REAL(rate)

        CALL SYSTEM_CLOCK(tic,rate)
!       Get this info to all procs
!       Memory issues, so unfortunately just do this piece by piece
!       Reduce only the first half, then use symmetry to get the back half
!       For now, this is the main scaling cost
        DO n = 0, Y%p-1
            DO m = -n,n
                im = m + (Y%p)
                DO ic2 = 1, prob%NCell
                    IF(m .gt. 0) THEN ! Symmetry
                        A2f(:,:,im,n+1,:,:,:, ic2) = CONJG(A2f(:,:,im - 2*m,n+1,:,:,:, ic2))*(-1D0)**m
                        CYCLE
                    ENDIF
                    !!! Memory: don't need all of these, can deallocate some parts
                    A2f(:,:,im,n+1,:,:,:,ic2) = prob%cm%reduce(A2f(:,:,im,n+1,:,:,:,ic2))
                ENDDO
            ENDDO
        ENDDO
        
        ! CALL prob%cm%barrier()
        CALL SYSTEM_CLOCK(toc)
        IF(prob%cm%mas()) print *, "Reduce: ", REAL(toc-tic)/REAL(rate)
    ENDIF

    ALLOCATE( &
        ut(info%NmatT), &
        uti(info%NmatT), &
        utr(info%NmatT), &
        IPIV(info%NmatT), wrk(info%NmatT), &
        swrk(info%NmatT*(info%NmatT+1)), &
        rwrk(info%NmatT), &
        A(info%NmatT, info%NmatT), &
        b(info%NmatT), &
        A2r(3,3, 2*(Y%p-1)+1, Y%p, Y%nt, Y%np), &
        A2t(3,3, 2*(Y%p-1)+1, Y%p, Y%nt, Y%np), &
        b2r(3*Y%nt*Y%np), &
        b2t(3*Y%nt*Y%np), &
        A2(info%Nmat, info%Nmat), &
        b2(info%Nmat))
    A = 0D0
    b = 0D0
    b2 = 0D0
    A2 = 0D0

    CALL SYSTEM_CLOCK(tic,rate)
!   Now the real space part
!   First loop: velocity surface
    DO ic = info%PCells(1,1), info%PCells(1,2)
        row = (ic - 1)*info%Nmat + 1
        cell => prob%cell(ic)

!       Get starting and ending theta current cell for parallel
        th_st  = info%PCells(2,1)
        th_end = info%PCells(2,2)

!       Second loop: integral surface
        DO ic2 = 1,prob%NCell
            col = (ic2-1)*info%Nmat + 1
            celli => prob%cell(ic2)
            
!           Get velocity sub-matrix for cell-cell combo (Same or diff cell)
            IF(ic .eq. ic2) THEN
                CALL cell%fluid(Ao = A2r, bo = b2r, itt1 = th_st, itt2 = th_end)
            ELSE
                CALL cell%fluid(Ao = A2r, bo = b2r, itt1 = th_st, itt2 = th_end, celli = celli)
            ENDIF

!           Here we need to take back in the velocity/integral cell combo info across all processors
!           that had a part in this combo
!           All info is sent to the processor that has the theta = 1 point (skip if one proc has all points)
            IF((prob%cm%id() .eq. info%CProcs(ic,1)) .and. (info%CProcs(ic,1) .ne. info%CProcs(ic,2))) THEN
                DO i = info%CProcs(ic,1) + 1, info%CProcs(ic,2)
                    CALL prob%cm%recv(A2t, i)
                    A2r = A2r + A2t

                    CALL prob%cm%recv(b2t, i)
                    b2r = b2r + b2t
                ENDDO
            ELSEIF(info%CProcs(ic,1) .ne. info%CProcs(ic,2)) THEN
                A2t = A2r
                CALL prob%cm%send(A2t, info%CProcs(ic,1)) 

                b2t = b2r
                CALL prob%cm%send(b2t, info%CProcs(ic,1))
            ENDIF

!           Add long range interactions if periodic
!           We also do after send/receive so it doesn't double add
            IF(info%periodic) THEN
                A2r = A2r + A2f(:,:,:,:,:,:,ic,ic2)*(1D0-cell%lam)/(1D0+cell%lam)

!               b2f is the vector at the velocity points integrated across all surfaces,
!               so we only add in once
                IF(ic .eq. ic2) b2r = b2r + b2f( (ic-1)*(3*Y%nt*Y%np)+1 : (ic)*(3*Y%nt*Y%np))/(1D0 + cell%lam)
            ENDIF

!           Perform the Galerkin
            A2 = 0D0
            b2 = 0D0
            IF(prob%cm%id() .eq. info%CProcs(ic,1)) THEN ! Only if in the processor with the info though
                
                CALL info%Gal(A2r, b2r, A2, b2) !!!!!!!!!!!!! Ideally could be parallelized...
!               Put sub matrix into big matrix
                A(row:row + info%Nmat - 1, col:col + info%Nmat - 1) = A2
!               Sum over all the integrals
                b(row:row + info%Nmat - 1) = b(row:row + info%Nmat - 1) + b2
            ENDIF
        ENDDO
    ENDDO
    ! CALL prob%cm%barrier()
    CALL SYSTEM_CLOCK(toc)
    IF(prob%cm%mas()) print *, 'Fluid: ',REAL(toc-tic)/REAL(rate)

    DO i = 1,info%NmatT
        A(:,i) = prob%cm%reduce(REAL(A(:,i))) + prob%cm%reduce(AIMAG(A(:,i)))*ii
    ENDDO
    b = prob%cm%reduce(REAL(b)) + prob%cm%reduce(AIMAG(b))*ii

!   Invert big matrix to get a list of all the vel constants of all cells
    ut = 0D0
    CALL SYSTEM_CLOCK(tic)
    IF(prob%cm%np() .gt. 1) THEN

!       Submatrix proc is responsible for
        ALLOCATE(Ap(size(prob%map_A)/2))

!       Fill this matrix
        DO i = 1,size(Ap)
            Ap(i) = A(prob%map_A(1,i), prob%map_A(2,i))
        ENDDO

!       Fill b if there's something to fill
        IF(ALLOCATED(prob%map_b)) THEN
            ALLOCATE(bp(size(prob%map_b)))
            bp = b(prob%map_b)
        ELSE
            ALLOCATE(bp(0))
            bp=0
        ENDIF

        CALL pzgesv( info%NmatT, 1, Ap(1), 1, 1, prob%desca, ipiv(1), &
                     bp(1), 1, 1, prob%descb, p_info )

        b = 0D0
        IF(ALLOCATED(prob%map_b)) b(prob%map_b) = bp
        ut = prob%cm%reduce(REAL(b)) + prob%cm%reduce(AIMAG(b))*ii
        
        IF(ALLOCATED(bp)) DEALLOCATE(bp)
        DEALLOCATE(Ap)

    ELSE
        CALL zcgesv(info%NmatT, 1, A, info%NmatT, IPIV, b, info%NmatT, &
                    ut, info%NmatT, wrk, swrk, rwrk, iter, p_info)
    ENDIF
    CALL SYSTEM_CLOCK(toc)
    IF(prob%cm%mas()) print *, 'Invert: ',REAL(toc-tic)/REAL(rate)

!   Advance in time now
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)

!       Reconstruct individual vels
        row = (ic-1)*info%Nmat + 1
        cell%umn = 0D0
        cell%umn(1,1:info%Nmat/3) = ut((/(i, i=row    , row + info%Nmat-2, 3)/))
        cell%umn(2,1:info%Nmat/3) = ut((/(i, i=row + 1, row + info%Nmat-1, 3)/))
        cell%umn(3,1:info%Nmat/3) = ut((/(i, i=row + 2, row + info%Nmat  , 3)/))

!       Volume correction: small, inward normal velocity based off current volume/SA/time step
!       Removed for reduce, because that keeps things at constant volume
        if(.not. reduce) THEN
            zm = -(cell%Vol() - cell%V0)/(cell%SA()*info%dt)
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
            ! cell%x(1,:,:) = Y%backward(cell%xmn(1,:))
            ! cell%x(2,:,:) = Y%backward(cell%xmn(2,:))
            ! cell%x(3,:,:) = Y%backward(cell%xmn(3,:))
        ELSEIF(ord .ne. 1) THEN
            print*, "ERROR: time advancement of order >1 not supported"
            stop
        ENDIF
    
!       Update position and current time step
        cell%xmn = xmnt + umnt*info%dt

        cell%x(1,:,:) = Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = Y%backward(cell%xmn(3,:))
    ENDDO

!   Advance periodic basis vectors
    IF(info%periodic) THEN
        CALL info%bvAdvance()
!       Update array containing support points
        CALL SuppPoints(info)
    ENDIF

    prob%cts = prob%cts + 1
    prob%t = prob%t + info%dt

!   Check if there's any funny business
    IF(ANY(isNaN((ABS(cell%umn)))) .or. ANY((ABS(cell%umn)) .gt. HUGE(zm))) THEN
            print *, 'ERROR: inftys or NaNs'
            STOP
    ENDIF

!   Write and display some output
    CALL prob%write()
    CALL prob%output()
END SUBROUTINE UpdateProb

! -------------------------------------------------------------------------!
! Finds equivalent strain rate
FUNCTION StrnRtProb(prob) RESULT(Sig)
    CLASS(probType), TARGET :: prob
    TYPE(sharedType), POINTER :: info
    TYPE(YType), POINTER :: Y
    TYPE(cellType), POINTER :: cell
    REAL(KIND = 8) :: Sig(3,3), tmp(3,3), lam, tau
    REAL(KIND = 8), ALLOCATABLE :: fMxPuMn(:,:,:,:), u(:,:,:), n(:,:,:)
    INTEGER :: i, j, ic, nt, np, d1, d2

    info => prob%info
    Y => info%Yf
    nt = Y%nt
    np = Y%np
    
    tau = DOT(CROSS(info%bv(:,1), info%bv(:,2)), info%bv(:,3))

    tmp = 0D0

    ALLOCATE(fMxPuMn(3,3,nt, np), u(3,nt, np), n(3,nt, np))

!   One paper says to use x - x_c,
!   but I shouldn't need to as long as net force in cell=0

    DO ic = 1, prob%NCell
        cell => prob%cell(ic)
        lam = cell%lam
        u(1,:,:) = Y%backward(cell%umn(1,:), info%P)
        u(2,:,:) = Y%backward(cell%umn(2,:), info%P)
        u(3,:,:) = Y%backward(cell%umn(3,:), info%P)

        n(1,:,:) = Y%backward(cell%nkmn(1,:))
        n(2,:,:) = Y%backward(cell%nkmn(2,:))
        n(3,:,:) = Y%backward(cell%nkmn(3,:))
!       Get velocity and normalss at integration points
        DO i = 1,nt
            DO j = 1,np
                fMxPuMn(:,:,i,j) = OUTER(cell%ff(:,i,j), cell%xf(:,i,j))  & !! Need Ca here? don't seem to need to.
                                 + (lam - 1D0)*(OUTER(u(:,i,j), n(:,i,j)) &     !! Run with sphere array, compare to Zinchenko
                                 + OUTER(n(:,i,j), u(:,i,j)))
            ENDDO
        ENDDO
        DO d1 = 1,3
        DO d2 = 1,3
            tmp(d1, d2) = tmp(d1, d2) + cell%intg(fMxPuMn(d1,d2,:,:))
        ENDDO
        ENDDO
        ! print *, u(1,(nt+1)/2,1)! - cell%x(3,2,1)
        ! print *, ' '
        ! print *, cell%intg(fMxPuMn(1,3,:,:)), u(1,1,1), SUM(u(1,:,:))/nt/np
    ENDDO
    ! stop

    Sig = tmp/tau + info%dU + TRANSPOSE(info%dU)

END FUNCTION StrnRtProb

! -------------------------------------------------------------------------!
! Continue from where we left off at? 
!!! Should make this more robust (read Params, but at the top, etc.)
SUBROUTINE ContinueProb(prob)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell
    CHARACTER (LEN = 35) contfile, cfile2
    INTEGER ic, endt, stat, p, i, jmp
    REAL(KIND = 8), ALLOCATABLE :: xmnraw(:,:)
    LOGICAL :: exist_in

    DO ic = 1,prob%NCell
        cell => prob%cell(ic)
        INQUIRE(FILE = TRIM('dat/'//prob%fileout//'/'//TRIM(cell%fileout)//'/maxdt'), EXIST = exist_in)
        IF(.not. exist_in) THEN
            print *, "ERROR: Files not found for continue. Have you run previous cases with the same number of cells?"
            stop
        ENDIF

!       Get the timestep of where we left off
        contfile = TRIM('dat/'//prob%fileout//'/'//TRIM(cell%fileout))
        cfile2 = TRIM('dat/'//prob%fileout//'/'//TRIM(cell%fileout)//'/maxdt')
        OPEN(unit = 13, file = TRIM(cfile2), action = 'read')
        READ(13, '(I16)', iostat = stat) endt
        CLOSE(13)
        IF(endt .eq. 0) THEN
!           Starting at 0 anyways, just return
            RETURN
        ENDIF
        prob%cts = endt

!       Open the file of where we left off and get the shape
        write(cfile2, "(I0.5)") endt
        cfile2 = TRIM('dat/'//prob%fileout//'/'//TRIM(cell%fileout)//'/x_'//TRIM(cfile2))

!       Read in the files
        p = (cell%info%p + 1)*(cell%info%p + 1)*2
        IF(.not.ALLOCATED(xmnraw)) ALLOCATE(xmnraw(3,p))
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

        CALL cell%derivs()
        CALL cell%stress()        
    ENDDO

!   Also need to bring the basis vectors back in. Only works for shear flow 1 right now
    prob%info%dU = 0D0
    prob%info%dU(1,3) = 1D0

    prob%info%bv = prob%info%bv + MATMUL(prob%info%dU, prob%info%bv)*(endt - 5)*prob%info%dt
    DO ic = 1,5
        CALL prob%info%bvAdvance()
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
        prob%cts, prob%t, ic, MAXVAL(ABS(prob%cell(ic)%ff))*prob%cell(ic)%Ca, &
        MAXVAL(ABS(prob%cell(ic)%umn)), prob%cell(ic)%vol(), prob%cell(ic)%SA()
    ENDDO
END SUBROUTINE OutputProb

END MODULE PROBMOD