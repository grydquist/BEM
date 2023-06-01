MODULE PROBMOD
USE SHAPEMOD
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
    INTEGER :: ictxt, desca(9), descb(9), thti_f(2), der_n(2), der_m(2)
    INTEGER, ALLOCATABLE :: map_A(:,:), map_b(:)

    CONTAINS
    PROCEDURE :: Update  => UpdateProb
    PROCEDURE :: Write   => WriteProb
    PROCEDURE :: Output  => OutputProb
    PROCEDURE :: Continue=> ContinueProb
    PROCEDURE :: StrnRt  => StrnRtProb
    PROCEDURE :: LHSCmp  => LHSCmpProb
    PROCEDURE :: LHS => LHSGalProb
    PROCEDURE :: RHSGalProb
    PROCEDURE :: RHS => RHSProb
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
               xblx, yblx, itb, i, ib, j, jb, strti, strtj, nbx, nby, ntm, npm
    INTEGER, ALLOCATABLE :: all_thti_f(:)
    REAL(KIND = 8) :: Gfac, zm
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

!   Hardcoded shear
    IF(info%shear) THEN
        info%dU = 0D0
        info%dU(1,3) = 1D0
    ENDIF

!   Hardcoded periodic extension
    IF(info%extens) THEN
        info%dU = ExtensdU()
    ENDIF

!   Hardcoded quiescent
    IF(info%none) THEN
        info%dU = 0D0
    ENDIF

    celltmp = cellType(filein, reduce, info, props)

!   Loop and make individual cells
    DO ic = 1, prob%NCell
        prob%cell(ic) = celltmp
        prob%cell(ic)%id = ic
        write(icC, "(I0.1)") ic
        prob%cell(ic)%fileout = TRIM('cell_'//icC)     
!       Singular integration stuff
        ntm = MAX(info%Y%nt, info%Ys%nt)
        npm = MAX(info%Y%np, info%Ys%np)
        ALLOCATE(prob%cell(ic)%sing_inds(2, info%Y%nt, info%Y%np, prob%NCell), &
                 prob%cell(ic)%cut_box  (3, info%Y%nt, info%Y%np, ntm, npm, prob%NCell))
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
!   Exception for sequential below
    IF(prob%cm%np().eq.1) THEN
        prob%thti_f(1) = 1
        prob%thti_f(2) = info%Yf%nt
        prob%der_n(1)  = 0
        prob%der_n(2)  = info%Y%p
        prob%der_m(1)  = 0
        prob%der_m(2)  = info%Y%p
    ENDIF
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

!       For the parallel part of the stress. This is different b/c we
!       need to divvy this up so that all procs work on a single cell
        ALLOCATE(all_thti_f(prob%cm%np()))
        all_thti_f = ParallelSplit(info%Yf%nt, prob%cm%np())
!       Lower bound of thetas this cell should handle
        prob%thti_f(1) = all_thti_f(prob%cm%id() + 1)

!       Upper bound (with exception for if this is the top proc)
        IF(prob%cm%id() + 1 .ne. prob%cm%np()) THEN
            prob%thti_f(2) = all_thti_f(prob%cm%id() + 2) - 1
        ELSE
            prob%thti_f(2) = info%Yf%nt
        ENDIF
!       Exception for if there's more procs than points
        IF(prob%cm%id() + 1 .gt. info%Yf%nt) prob%thti_f(2) = prob%thti_f(1) - 1

!       Same process for the derivatives, but we need to but the results
!       in terms of a starting/ending p/m
        all_thti_f = ParallelSplit((info%Y%p + 1)*(info%Y%p + 1), prob%cm%np())

!       Lower n/m: cycle til I find the right one
        DO ic = info%Y%p, 0, -1
            IF(all_thti_f(prob%cm%id() + 1) .gt. ic*ic) THEN
                prob%der_n(1) = ic
                prob%der_m(1) = all_thti_f(prob%cm%id() + 1) - ic*ic - ic - 1
                EXIT
            ENDIF
        ENDDO

        IF(prob%cm%id() + 1 .ne. prob%cm%np()) THEN
            DO ic = info%Y%p, 0, -1
                IF(all_thti_f(prob%cm%id() + 2) .gt. ic*ic) THEN
                    prob%der_n(2) = ic
                    prob%der_m(2) = all_thti_f(prob%cm%id() + 2) - ic*ic - ic - 2
                    EXIT
                ENDIF
            ENDDO
        ELSE
            prob%der_n(2) = info%Y%p
            prob%der_m(2) = info%Y%p
        ENDIF
    ENDIF
!   How many timesteps from G do we actually need?
    Gfac = prob%nts*prob%kfr/(prob%NT*prob%info%dt) ! Gfac: Ratio of total time in velgrad to total time requested
!   If we have fewer velocity gradient time steps than requested, set down to gradient
    IF(Gfac .le. 1) THEN
        prob%NT = FLOOR(prob%NT*Gfac)
        IF(prob%cm%mas()) print *, "Warning: Fewer gradient time steps than requested total time steps"
        IF(prob%cm%mas()) print *, "Setting total timesteps to "
        IF(prob%cm%mas()) print *, prob%NT
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
!   Exception here for not continue !!!
!!  ============================
    IF(prob%cm%mas()) print *, 'Getting initial shape -'
    IF(prob%cm%mas()) print *, 'Max velocity coefficient w.r.t. membrane time (want less than 0.005*Ca):'
    Gfac = (((4D0*PI/3D0)/((0.95D0**(3D0))*(PI/6D0)))**(1D0/3D0))/2D0
    IF(.not. info%periodic) info%bvl = 5.1174D0
    DO ic = 1, prob%NCell
            CALL prob%cell(ic)%relax(0.005D0, prob%cm)!!!
            CALL prob%cell(ic)%derivs()
            CALL prob%cell(ic)%stress()
            prob%cell(ic)%V0 = prob%cell(ic)%Vol()
            !! Test                                        !!!!!
            SELECT CASE(ic-1)!MOD(ic, 4))!9) )
            CASE(0)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = .45D0/(ispi*0.5D0)
            CASE(1)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = .45D0/(ispi*0.5D0)
            CASE(2)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + info%bvl/4D0)/(ispi*0.5D0)
            CASE(3)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + info%bvl/4D0)/(ispi*0.5D0)
            CASE(4)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + info%bvl/2D0)/(ispi*0.5D0)
            CASE(5)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + info%bvl/2D0)/(ispi*0.5D0)
            CASE(6)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(7)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(8)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + info%bvl/2D0)/(ispi*0.5D0)
            CASE(9)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(10)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/4D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(11)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(12)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(13)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(14)
                prob%cell(ic)%xmn(1,1) = 1.4D0/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            CASE(15)
                prob%cell(ic)%xmn(1,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(2,1) = (1.4D0 + info%bvl/2D0)/(ispi*0.5D0)
                prob%cell(ic)%xmn(3,1) = (.45D0 + 3D0*info%bvl/4D0)/(ispi*0.5D0)
            END SELECT
            IF(prob%Ncell.eq.1) prob%cell(ic)%xmn(:,1) = 0D0
            prob%cell(ic)%x(1,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(1,:))
            prob%cell(ic)%x(2,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(2,:))
            prob%cell(ic)%x(3,:,:) = prob%cell(ic)%info%Y%backward(prob%cell(ic)%xmn(3,:))
    ENDDO
!    ============================
    
!   Should we continue from a spot where we left off?
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
    REAL(KIND =8), ALLOCATABLE :: ff(:,:), intg1, intg2, intg3

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

!   Write drag coeff
    ! filename = TRIM('K_'//ctsst)
    ! ff = cell%ff(1,:,:)
    ! intg1 = -cell%intg(ff)
    ! ff = cell%ff(2,:,:)
    ! intg2 = -cell%intg(ff)
    ! ff = cell%ff(3,:,:)
    ! intg3 = -cell%intg(ff)
    ! OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
    ! ! WRITE(88,*) -cell%intg(ff)/(6D0*pi*REAL(prob%cell(1)%umn(3,1))*ispi*0.5D0)
    ! WRITE(88,*) SQRT(intg1*intg1+intg2*intg2+intg3D0*intg3)/(6D0*pi*NORM2(REAL(cell%umn(:,1)))*ispi*0.5D0)
    ! CLOSE(88)

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
    WRITE(88,*) prob%t
    CLOSE(88)
    
!   Output strains if requested
    IF(prob%info%OutputStrain) THEN
        filename = TRIM('strains_'//ctsst)
        OPEN (UNIT = 88, FILE = TRIM(datdir)//TRIM(filename))
        WRITE(88,*) cell%avgJ
        WRITE(88,*) cell%maxJ
        WRITE(88,*) cell%avgSh
        WRITE(88,*) cell%maxSh
        CLOSE(88)
    ENDIF
    
    IF(.not. cell%writ) cell%writ = .true.
    ENDDO

END SUBROUTINE WriteProb

! -------------------------------------------------------------------------!
! Problem advancement
SUBROUTINE UpdateProb(prob, ord, reduce)
    CLASS(probType), INTENT(INOUT) :: prob
    TYPE(cellType), POINTER :: cell, celli
    TYPE(YType), POINTER :: Y
    TYPE(sharedType), POINTER :: info
    INTEGER, INTENT(IN) :: ord
    INTEGER :: ic, ic2, i, j, row, it, th_st, th_end, k, &
        im, m, n, mnmatT
    INTEGER(KIND = 8) tic, toc, rate
    LOGICAL, INTENT(IN) :: reduce
    REAL(KIND = 8) :: zm, res, b_n
    REAL(KIND = 8), ALLOCATABLE ::  fv1(:), fv2(:), fv3(:), fv(:,:), &
                                    nv1(:), nv2(:), nv3(:), nv(:,:), &
                                    xv(:,:), e(:)
    COMPLEX(KIND = 8), ALLOCATABLE :: fmnR(:,:), nmnR(:,:), &
                                      Aumn(:), bmn(:), r(:), y_v(:), &
                                      cs(:), sn(:), beta(:), betat(:), &
                                      Ht(:,:), e1(:), Q_k(:),&
                                      Q_kp1(:), Q(:,:), H(:,:)
    COMPLEX(KIND = 8) :: temp

    Y => prob%info%Y
    info => prob%info

    mnmatT = 3*(Y%p + 1)*(Y%p + 1)*info%NCell
    ALLOCATE(&
        r(mnmatT), &
        Q_k(mnmatT), &
        Q_kp1(mnmatT), &
        Q(mnmatT, info%GMRES_it + 1), &
        H(info%GMRES_it + 2, info%GMRES_it+1), &
        cs(info%GMRES_it), &
        sn(info%GMRES_it), &
        e (info%GMRES_it + 1), &
        e1(info%GMRES_it + 1), &
        beta(info%GMRES_it + 1))

!   If periodic, make variables for the Fourier part
    IF(info%periodic) THEN
        ALLOCATE( &
        fv(3,Y%nt*Y%np*prob%NCell), &
        nv(3,Y%nt*Y%np*prob%NCell), &
        fmnR(3, (info%p+1)*(info%p+1)), &
        nmnR(3, (info%p+1)*(info%p+1)))
    ENDIF

    IF(info%shear) THEN
!       Hardcoded shear
        info%dU = 0D0
        info%dU(1,3) = 1D0
    ELSEIF(info%extens) THEN
!       Hardcoded periodic extension
        info%dU = ExtensdU()
    ELSEIF(info%none) THEN
!       Quiescent
        info%dU = 0D0
    ELSE
!       Do interpolation to get current grad tensor, then normalize by kolm time
        info%dU = VelInterp(prob%G,prob%t,prob%nts,prob%kfr)*prob%kdt
    ENDIF

    CALL SYSTEM_CLOCK(tic,rate)
!   First, calculate geometric quantities and stresses in all of the cells
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)

!       For periodic, check if cells are outside primary cell, and put them in if so
        IF(info%periodic) THEN
            CALL cell%inPrim()
        ENDIF

!       These algorithms have been parallelized, but I am not implementing
!       derivs here because the overhead makes it more expensive
        CALL cell%derivs()!prob%der_m, prob%der_n, prob%cm)
        CALL cell%stress()!prob%thti_f, prob%cm)

!       Find which constants should be used in reconstruction
        cell%xrecon(1,:) = Y%reconConsts(cell%xmn(1,:))
        cell%xrecon(2,:) = Y%reconConsts(cell%xmn(2,:))
        cell%xrecon(3,:) = Y%reconConsts(cell%xmn(3,:))
        
        cell%frecon(1,:) = Y%reconConsts(cell%fmn(1,:))
        cell%frecon(2,:) = Y%reconConsts(cell%fmn(2,:))
        cell%frecon(3,:) = Y%reconConsts(cell%fmn(3,:))

        cell%nJrecon(1,:) = Y%reconConsts(cell%nkt(1,:))
        cell%nJrecon(2,:) = Y%reconConsts(cell%nkt(2,:))
        cell%nJrecon(3,:) = Y%reconConsts(cell%nkt(3,:))

!       Construct a large vector containing all the forces/locations in
!       preparation for the fast Ewald summation. I do these outside b/c
!       I don't want to build them everytime
        IF(info%periodic) THEN
!           Forces and locations
            fmnR = cell%fmn(:,1:(info%p+1)*(info%p+1))
            nmnR = cell%nkt(:,1:(info%p+1)*(info%p+1))

!           Append vectors
            CALL EwaldPreCalc(fv1, Y, xv, fmn=fmnR(1,:), x0in=cell%x)
            CALL EwaldPreCalc(fv2, Y, fmn=fmnR(2,:))
            CALL EwaldPreCalc(fv3, Y, fmn=fmnR(3,:))

            CALL EwaldPreCalc(nv1, Y, fmn=nmnR(1,:))
            CALL EwaldPreCalc(nv2, Y, fmn=nmnR(2,:))
            CALL EwaldPreCalc(nv3, Y, fmn=nmnR(3,:))
        ENDIF
    ENDDO

    IF(info%periodic) THEN; 
        fv(1,:) = fv1;  fv(2,:) = fv2;  fv(3,:) = fv3;
        nv(1,:) = nv1;  nv(2,:) = nv2;  nv(3,:) = nv3;
        DEALLOCATE(fv1, fv2, fv3, nv1, nv2, nv3)
    ENDIF

    ! CALL prob%cm%barrier()
    CALL SYSTEM_CLOCK(toc)
    IF(prob%cm%mas()) print *, "Stress/derivs: ", REAL(toc-tic)/REAL(rate)

!   Now we need to calculate the SpH velocity coefficients of the cells.
!   For periodic, we calculate LHS and RHS vectors in two parts:
!   Long range (via FFTs) and short range (directly, truncated at low r).
!   We calculate these two matrices individually at integration points
!   and add them in the below routines

!   Prep for GMRES loop
!   RHS
    bmn = prob%RHSGalProb(xv, fv)
!   LHS with RHS as intial guess
    Aumn = prob%LHS(bmn, xv, nv)

    r = bmn - Aumn
    res = NORM2C(r)

    Q = 0D0
    Q(:,1) = r/res
    H = 0D0
    e1 = 0D0
    e1(1) = 1D0
    k = 0
    cs = 0D0
    sn = 0D0
    b_n = NORM2C(bmn)
    beta = res * e1;

    IF(res.gt.info%GMRES_tol) THEN ! To account for rare case where res=0
!       Now the GMRES loop
        DO k = 1, info%GMRES_it

!           Finding next Q using previous Q in LHS
            CALL SYSTEM_CLOCK(tic,rate)
                Q_k = Q(:,k)
                Q_kp1 = prob%LHS(Q_k, xv, nv)
            CALL SYSTEM_CLOCK(toc)

!           First the Arnoldi iteration
            DO i = 1, k
                H(i, k) = DOT_PRODUCT(Q(:,i), Q_kp1)
                Q_kp1 = Q_kp1 - H(i, k) *Q(:,i)
            ENDDO
            H(k + 1, k) = NORM2C(Q_kp1)
            Q_kp1 = Q_kp1/H(k + 1, k)
            Q(:,k + 1) = Q_kp1

!           Rest of GMRES is just kinda details
            DO i = 1,k - 1
                temp = cs(i)*H(i,k) + sn(i)*H(i+1,k)
                H(i+1,k) = -sn(i) * H(i,k) + cs(i) * H(i + 1,k)
                H(i,k)   = temp
            ENDDO

            temp = SQRT(H(k,k)*H(k,k)+ H(k+1,k)*H(k+1,k))
            cs(k) = H(k,k) / temp
            sn(k) = H(k+1,k) / temp

            H(k, k) = cs(k) * H(k, k) + sn(k) * H(k + 1, k)
            H(k + 1, k) = 0D0

            beta(k + 1) = -sn(k) * beta(k)
            beta(k)     = cs(k) * beta(k)
            res       = ABS(beta(k + 1)) / b_n
            e(k + 1) = res

            IF(prob%cm%mas())print *, k, res, REAL(toc-tic)/REAL(rate)

            IF(res .lt. info%GMRES_tol) exit
        ENDDO

        IF(res .gt. info%GMRES_tol) k = k - 1
        Ht = H(1:k,1:k)
        betat = beta(1:k)

        y_v = GaussElimC(Ht ,betat, k)
    ELSE
        k = 1
        ALLOCATE(y_v(1))
        y_v = 0
        Q = 0
    ENDIF

    row = 1
    cell%umn = 0D0
    DO ic = 1, prob%NCell
        cell => prob%cell(ic)
        it = 1
        im = 0
        DO n  = 0, Y%p
            DO m = -n,n
                im = im + 1
                cell%umn(:,im) = bmn(row + it - 1:row + it + 1) &
                               + MATMUL(Q(row + it - 1:row + it + 1,1:k), y_v)
                IF(n .eq. Y%p) cell%umn(:,im) = 0D0 ! Explicitly set to 0
                IF(n .eq. Y%p) cell%xmn(:,im) = 0D0 ! b/c caused trouble
                Q_k(row + it - 1:row + it + 1) = cell%umn(:,im)
                it = it+3
            ENDDO
        ENDDO

!       Volume correction: small, inward normal velocity based off current volume/SA/time step
!       Removed for reduce, because that keeps things at constant volume
        IF(.not. reduce) THEN
            zm = -(cell%Vol() - cell%V0)/(cell%SA()*info%dt)

!           Changed so that I'm more careful with last term
            cell%umn(:,1:((cell%info%p-1)*(cell%info%p-1))) &
               = cell%umn (:,1:((cell%info%p-1)*(cell%info%p-1))) &
               + zm*cell%nkmn(:,1:((cell%info%p-1)*(cell%info%p-1)))
        ENDIF

!       Volume reduction (add small inward normal vel every timestep)
        IF(reduce) THEN
            IF(cell%Vol().gt. 4.22  ) cell%umn = cell%umn &
                                               - 0.10D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
            IF(cell%Vol().lt. 4.185 ) cell%umn = cell%umn &
                                               + 0.01D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
            IF(cell%Vol().gt. 4.1894) cell%umn = cell%umn &
                                               - 0.01D0*cell%nkmn(:,1:((cell%info%p+1)*(cell%info%p+1)))
        ENDIF

        cell%xmn = cell%xmn + cell%umn*info%dt

        cell%x(1,:,:) = Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = Y%backward(cell%xmn(3,:))
        row = row + mnmatT/info%NCell
    ENDDO

!   Advance periodic basis vectors
    IF(info%periodic) THEN
        CALL info%bvAdvance(prob%t)
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
! Works with the LHSCmpProb to do Galerkin
FUNCTION LHSGalProb(prob, umn, xv, nv, dbg) RESULT(Aumn)
    CLASS(probType), TARGET :: prob
    TYPE(sharedType), POINTER :: info
    REAL(KIND = 8), ALLOCATABLE :: Au(:), xv(:,:), nv(:,:), Au1(:)
    COMPLEX(KIND = 8), ALLOCATABLE :: Aumn(:), umn(:), Aumn1(:)
    LOGICAL, OPTIONAL :: dbg
    LOGICAL :: dbgflg = .false.
    INTEGER :: row, row2, ic, mnmat, th_st, th_end

    info => prob%info

    mnmat = 3*(info%Y%p + 1)*(info%Y%p + 1)
    IF(PRESENT(dbg)) dbgflg = dbg
    Aumn1 = 0d0
    ALLOCATE(Aumn(mnmat*info%NCell), &
             Au1(3*info%Y%nt*info%Y%np))
    Aumn = 0D0

!   Get the vector with the LHS evaluated at all of the intg points
    Au = prob%LHSCmp(umn, xv, nv, dbgflg)

!   Then do the Galerkin over all cells
    row = 1 + (info%PCells(1,1) - 1)*info%Nmat
    row2= 1 + (info%PCells(1,1) - 1)*mnmat

    DO ic = info%PCells(1,1), info%PCells(1,2)

        th_st  = info%PCells(2,1)
        th_end = info%PCells(2,2)

        Au1 = Au(row:row + info%Nmat - 1)
        CALL info%Gal(Au1, Aumn1, th_st, th_end)
        Aumn(row2:row2 + mnmat - 1) = Aumn1

        row = row + info%Nmat
        row2= row2+ mnmat
    ENDDO

    Aumn = prob%cm%reduce(Aumn)

END FUNCTION LHSGalProb

! -------------------------------------------------------------------------!
! For a complex LHS
FUNCTION LHSCmpProb(prob, umn, xv, nv, dbg) RESULT(Au)
    CLASS(probType), TARGET :: prob
    REAL(KIND = 8), ALLOCATABLE :: xv(:,:), nv(:,:), Au(:), vec_t(:), &
                                   ur(:), urm(:,:,:), urmT(:,:,:,:), Au_t(:), &
                                   urbm(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: u3T(:,:,:), u3(:,:)
    COMPLEX(KIND = 8), ALLOCATABLE :: umn(:)
    LOGICAL, OPTIONAL :: dbg
    LOGICAL :: dbgflg=.false.
    TYPE(sharedType), POINTER :: info
    TYPE(cellType), POINTER :: cell, celli
    REAL(KIND = 8) :: lam
    INTEGER :: row, col, ic, ic2, th_st, th_end, i, j, it, n, m, im, dim, it2, mnmat, tic, toc, rate

    IF(PRESENT(dbg)) dbgflg = dbg
    info => prob%info

    mnmat = 3*(info%Y%p + 1)*(info%Y%p + 1)

    lam = prob%cell(1)%lam
    ALLOCATE(u3T(3, (info%Y%p + 1)*(info%Y%p + 1), prob%NCell), &
             u3 (3, (info%Y%p + 1)*(info%Y%p + 1)), Au(info%NmatT), &
             ur(3*info%Y%nt*info%Y%np*info%NCell), &
             urmT(3, info%Y%nt, info%Y%np, info%NCell), &
             urm(3, info%Y%nt, info%Y%np), Au_t(info%NmatT), &
             urbm(3, info%NCell))

    CALL SYSTEM_CLOCK(tic,rate)
!   One argument takes a matrix as an input
    it2 = 1
    u3T = 0D0
    row = 1
    DO ic = 1,prob%NCell
        it = 1
        im = 1
        DO n  = 0, info%Y%p
            DO m = -n,n
                u3T(:,im,ic) = umn(row + it - 1:row + it + 1)

!               The rigid body calculation can be done exactly, so I remove it and
!               add it in later
                IF(n.eq.0) urbm(:, ic) = REAL(u3T(1:3, im, ic))*(0.5D0*ispi)
                IF(n.eq.0) u3T(:,im,ic) = 0D0

                it = it+3
                im = im + 1
            ENDDO
        ENDDO
        row = row + mnmat

!       We need these for periodic, but are also done in RHS_real, 
!       So may as well just do
        urmT(1,:,:,ic) = info%Y%backward(u3T(1,:,ic))
        urmT(2,:,:,ic) = info%Y%backward(u3T(2,:,ic))
        urmT(3,:,:,ic) = info%Y%backward(u3T(3,:,ic))
        DO i = 1, info%Y%nt
            DO j = 1, info%Y%np
                DO dim = 1,3
                    ur(it2) = urmT(dim, i, j, ic)
                    it2 = it2 + 1
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    CALL SYSTEM_CLOCK(toc)
    ! IF(prob%cm%mas()) print *, "Precalcs: ", REAL(toc-tic)/REAL(rate)

    Au = 0D0
    CALL SYSTEM_CLOCK(tic,rate)
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

            u3 = u3T(:,:,ic2)
            urm = urmT(:,:,:,ic2)
            vec_t = 0D0
            IF(ic .eq. ic2) THEN
                CALL cell%LHS_real(v_input_mn = u3, v_input = urm, &
                v = vec_t, itt1 = th_st, itt2 = th_end, dbg=dbgflg)
            ELSE
                CALL cell%LHS_real(v_input_mn = u3, v_input = urm, &
                v = vec_t, itt1 = th_st, itt2 = th_end, celli = celli, dbg=dbgflg)
            ENDIF
            Au(row:row + info%Nmat - 1) = Au(row:row + info%Nmat - 1) + vec_t
        ENDDO
    ENDDO
    CALL SYSTEM_CLOCK(toc)
    ! IF(prob%cm%mas()) print *, "Real: ", REAL(toc-tic)/REAL(rate)

    Au = prob%cm%reduce(Au)

!   And here is where we add the RBM back in
    it=1
    DO ic2 = 1,prob%NCell
        DO i = 1, info%Y%nt
            DO j = 1, info%Y%np
                Au(it:it + 2) = Au(it:it + 2) + urbm(:,ic2)*2D0/(1D0 + lam)
                it = it + 3
            ENDDO
        ENDDO
    ENDDO

!   Calculate Fourier part of double layer at all evaluation points
!   Reduction takes place within this loop
    CALL SYSTEM_CLOCK(tic,rate)

    Au_t = 0D0
    IF(info%periodic) THEN
        CALL EwaldT(info=info, x0=xv, f1=ur, n=nv, strt=1, u1=Au_t, full=.true., cm=prob%cm)
        Au = Au - Au_t*(1D0 - lam)/(4D0*pi*(1D0 + lam))! Most general is to loop over cells here
    ENDIF
    CALL SYSTEM_CLOCK(toc)
    ! IF(prob%cm%mas()) print *, "Ewald: ", REAL(toc-tic)/REAL(rate)

    ! CALL prob%cm%barrier()
    ! stop

END FUNCTION LHSCmpProb

! -------------------------------------------------------------------------!
! Finds equivalent strain rate
FUNCTION RHSGalProb(prob, xv, fv) RESULT(bmn)
    CLASS(probType), TARGET :: prob
    TYPE(sharedType), POINTER :: info
    REAL(KIND = 8), ALLOCATABLE :: xv(:,:), fv(:,:), b(:), b1(:)
    COMPLEX(KIND = 8), ALLOCATABLE :: bmn1(:), bmn(:)
    INTEGER :: ic, row, row2, mnmat, th_st, th_end

    info => prob%info

    mnmat = 3*(info%Y%p + 1)*(info%Y%p + 1)
    bmn1 = 0D0
    ALLOCATE(bmn(mnmat*info%NCell), &
             b1(info%Nmat))
    bmn = 0D0

!   Get the vector with the RHS evaluated at all of the intg points
    b = prob%RHS(xv, fv)

!   Then do the Galerkin over all cells
    row = 1 + (info%PCells(1,1) - 1)*info%Nmat
    row2= 1 + (info%PCells(1,1) - 1)*mnmat
    
    DO ic = info%PCells(1,1), info%PCells(1,2)

        th_st  = info%PCells(2,1)
        th_end = info%PCells(2,2)

        b1 = b(row:row + info%Nmat - 1)
        CALL info%Gal(b1, bmn1, th_st, th_end)
        bmn(row2:row2 + mnmat - 1) = bmn1

        row = row + info%Nmat
        row2= row2+ mnmat
    ENDDO

    bmn = prob%cm%reduce(bmn)

END FUNCTION RHSGalProb

! -------------------------------------------------------------------------!
! EHS calcultor at all int points
FUNCTION RHSProb(prob, xv, fv) RESULT(b)
    CLASS(probType), TARGET :: prob
    REAL(KIND = 8), ALLOCATABLE :: xv(:,:), fv(:,:), b(:), vec_t(:), b_t(:)
    TYPE(sharedType), POINTER :: info
    TYPE(cellType), POINTER :: cell, celli
    REAL(KIND = 8) :: lam
    INTEGER :: row, col, ic, ic2, th_st, th_end, i, j, it

    info => prob%info

    lam = prob%cell(1)%lam
    ALLOCATE(b(info%NmatT), &
           b_t(info%NmatT))

    b = 0D0

!   Now real-space part
    DO ic = info%PCells(1,1), info%PCells(1,2)
        row = (ic - 1)*info%Nmat + 1
        cell => prob%cell(ic)
!       Get starting and ending theta current cell for parallel
        th_st  = info%PCells(2,1)
        th_end = info%PCells(2,2)

!       Second loop: integral surface
        DO ic2 = 1,prob%NCell
            celli => prob%cell(ic2)
            
            vec_t = 0D0
            IF(ic .eq. ic2) THEN
                CALL cell%RHS_real(v_input_mn = cell%fmn, v = vec_t, &
                itt1 = th_st, itt2 = th_end)
            ELSE
                CALL cell%RHS_real(v_input_mn = celli%fmn, v = vec_t, &
                itt1 = th_st, itt2 = th_end, celli = celli)
            ENDIF
            b(row:row + info%Nmat - 1) = b(row:row + info%Nmat - 1) + vec_t
        ENDDO
    ENDDO

    b = prob%cm%reduce(b)

!   Fill out initial RHS vec, starting with periodic portion if present
    b_t = 0D0
    IF(info%periodic) THEN
        CALL EwaldG(info=info, x0=xv, f=fv, strt=1, u1=b_t, full=.true., cm=prob%cm)
        b = b - b_t/(4D0*pi*(1D0 + lam)) ! Most general is to loop over cells here
    ENDIF

END FUNCTION RHSProb

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
    CHARACTER (LEN = 35) contfile, cfile2, trash
    INTEGER ic, endt, stat, p, i, jmp, trash2, repar_num, p_prev
    REAL(KIND = 8) tmpdt, endtime, zm, lamp
    REAL(KIND = 8), ALLOCATABLE :: xmnraw(:,:), dUtmp(:,:), bvt(:,:)
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
        READ(13, *, iostat = stat) endtime
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

!       We may want to try seeing how different p's respond: read in 
!       old p
        contfile = TRIM('dat/'//prob%fileout//'/'//TRIM(cell%fileout)//'/Params')
        OPEN(unit = 13, file = TRIM(contfile), action = 'read')
        READ(13, "(a)", iostat = stat) trash
        READ(13, *, iostat = stat) p_prev
        CLOSE(13)

        p = (p_prev + 1)*(p_prev + 1)*2 !(cell%info%p + 1)*(cell%info%p + 1)*2
        IF(.not.ALLOCATED(xmnraw)) ALLOCATE(xmnraw(3,p))
        xmnraw = 0D0
        OPEN(unit = 13, file = cfile2, action = 'read')
        DO i = 1,p
            READ(13, *, iostat = stat) xmnraw(:,i)
        ENDDO
        CLOSE(13)

!       Text file format: all real, then imag
        p = cell%info%p
        p = MIN(p, p_prev)
        jmp = (p_prev+1)*(p_prev+1)
    
        cell%xmn = 0D0

        cell%xmn(1,1:(p+1)*(p+1)) = xmnraw(1,1:(p+1)*(p+1))
        cell%xmn(2,1:(p+1)*(p+1)) = xmnraw(2,1:(p+1)*(p+1))
        cell%xmn(3,1:(p+1)*(p+1)) = xmnraw(3,1:(p+1)*(p+1))

        cell%xmn(1,1:(p+1)*(p+1)) = cell%xmn(1,1:(p+1)*(p+1)) + xmnraw(1, jmp+1: jmp + (p+1)*(p+1))*ii
        cell%xmn(2,1:(p+1)*(p+1)) = cell%xmn(2,1:(p+1)*(p+1)) + xmnraw(2, jmp+1: jmp + (p+1)*(p+1))*ii
        cell%xmn(3,1:(p+1)*(p+1)) = cell%xmn(3,1:(p+1)*(p+1)) + xmnraw(3, jmp+1: jmp + (p+1)*(p+1))*ii
        
        cell%x(1,:,:) = prob%info%Y%backward(cell%xmn(1,:))
        cell%x(2,:,:) = prob%info%Y%backward(cell%xmn(2,:))
        cell%x(3,:,:) = prob%info%Y%backward(cell%xmn(3,:))

        CALL cell%derivs()
        CALL cell%stress()        
        
        cfile2 = TRIM('dat/'//prob%fileout//'/'//TRIM(cell%fileout)//'/Params')
        OPEN(unit = 13, file = TRIM(cfile2), action = 'read')
        READ(13, "(a)", iostat = stat) trash
        READ(13, '(I16)', iostat = stat) trash2
        READ(13, "(a)", iostat = stat) trash
        READ(13, *, iostat = stat) tmpdt
        CLOSE(13)

    ENDDO

!   Also need to bring the basis vectors back in.
    IF(prob%info%shear) THEN
        prob%info%dU = 0D0
        prob%info%dU(1,3) = 1D0
        prob%info%bv = prob%info%bv + MATMUL(prob%info%dU, prob%info%bv)*(endtime - 5*tmpdt)!)(endt - 5)*tmpdt!prob%info%dt

    ELSEIF(prob%info%extens) THEN
        prob%info%dU = 0D0
        zm = (3D0 + sqrt(9D0 - 4D0))/2D0
        zm = atan2(1D0 - zm, 1D0)
        zm = PI/2D0 + zm
        zm = -zm*2D0
        prob%info%dU(1,1) = cos(zm)
        prob%info%dU(1,3) = sin(zm)
        prob%info%dU(3,1) = sin(zm)
        prob%info%dU(3,3) =-cos(zm)
        prob%info%dU = prob%info%dU*SQRT(2D0)/2D0

        dUtmp = prob%info%dU*(endtime - 5*tmpdt)
        prob%info%bv = MATMUL(prob%info%bv, &
                       EXPM(dUtmp))

!       Now need to get the reparametrized box at this point in time
        lamp = LOG((3D0 + sqrt(9D0 - 4D0))/2D0)
        repar_num = FLOOR((endtime - 5*tmpdt)/lamp + 0.5D0)

        DO ic = 1,repar_num
            bvt = prob%info%bv
            prob%info%bv(1,1) = bvt(1,1) +     bvt(1,3)
            prob%info%bv(1,3) = bvt(1,1) + 2D0*bvt(1,3)
            prob%info%bv(3,1) = bvt(3,1) +     bvt(3,3)
            prob%info%bv(3,3) = bvt(3,1) + 2d0*bvt(3,3)
        ENDDO
    ENDIF

    DO ic = 1,5
        CALL prob%info%bvAdvance(endtime - ((5-ic)*tmpdt))
    ENDDO

    prob%t = endtime
END SUBROUTINE ContinueProb

! -------------------------------------------------------------------------!
! Some output
SUBROUTINE OutputProb(prob)
    CLASS(probType), INTENT(IN) :: prob
    INTEGER ic
    REAL(KIND =8), ALLOCATABLE :: ff(:,:), intg1, intg2, intg3
    IF(prob%cm%slv()) RETURN

!   What the output means
    IF(MOD(prob%cts,50).eq.1) THEN
        PRINT *, "  CTS  Time    Cell  Max F     Max umn  Vol     SA"
        PRINT *, "----------------------------------------------------"
    ENDIF

!   Drag stuff
    ! ff = prob%cell(1)%ff(1,:,:)
    ! intg1 = -prob%cell(1)%intg(ff)
    ! ff = prob%cell(1)%ff(2,:,:)
    ! intg2 = -prob%cell(1)%intg(ff)
    ! ff = prob%cell(1)%ff(3,:,:)
    ! intg3 = -prob%cell(1)%intg(ff)
    ! print *, SQRT(intg1*intg1 + intg2*intg2 + intg3D0*intg3)/(6D0*pi*NORM2(REAL(prob%cell(1)%umn(:,1)))*ispi*0.5D0)

!   The output
    DO ic = 1, prob%NCell
        prob%cell(ic)%umn(:,1) = 0D0
        write(*,'(I5,X,F8.4,X,I5, X, F11.7,X,F11.7,X,F8.4,X,F8.4,X,F8.4,X,F8.4)') & 
        prob%cts, prob%t, ic, MAXVAL(ABS(prob%cell(ic)%ff))*prob%cell(ic)%Ca, &
        MAXVAL(ABS(prob%cell(ic)%umn)), prob%cell(ic)%vol(), prob%cell(ic)%SA()
    ENDDO
END SUBROUTINE OutputProb

END MODULE PROBMOD