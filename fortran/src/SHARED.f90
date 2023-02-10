MODULE SHAREDMOD
USE HARMMOD
IMPLICIT NONE

!==============================================================================!
!             This module is a spot for the container that has all             !
!            info shared between cells, so there's just one pointer            !
!                           to a "master" of sorts                             !
!==============================================================================!

!==============================================================================!
!================================= CONTAINERS =================================!
!==============================================================================!

! Contains all shared info that is the same between cells
! The idea is that there is one of these, and the cells all have a pointer to it
TYPE sharedType

!   Harmonics info
    INTEGER :: p, q, ftot
    TYPE(YType) :: Y, Yf, Ys
    REAL(KIND = 8) :: thtc, k, h

!   Normalized Legendres/exp's at GPs
    COMPLEX(KIND =8), ALLOCATABLE :: es(:,:), esf(:,:)
    REAL(KIND = 8), ALLOCATABLE :: cPmn(:,:,:), cPmnf(:,:,:), xsf(:)

!   Matrix size info
    INTEGER :: Nmat, NmatT

!   Time step
    REAL(KIND = 8) :: dt

!   Velocity gradient
    REAL(KIND = 8) :: dU(3,3)

!   Periodic information (using periodic, basis vecs, volume, etc)
!   Basis vectors are in columns of bv
    LOGICAL :: periodic
    REAL(KIND = 8) :: eye(3,3), bv(3,3), tau, bvl, kv(3,3), xi, eta, par, parexp 
    INTEGER :: gp, suppPoints
    INTEGER, ALLOCATABLE :: suppMat(:,:)

!   3D FFT Parameters
    REAL(KIND = 8), ALLOCATABLE :: WSAVE(:)
    INTEGER :: LENSAVE

!   Cell-cell parameters
    REAL(KIND = 8) :: epsi, r0, D, Beta
    LOGICAL :: CellCell
    INTEGER :: NCell

!   MPI Stuff
    INTEGER :: PCells(2, 2)
    INTEGER, ALLOCATABLE :: CProcs(:,:)

!   Special flow cases (extensional, shear)
    LOGICAL :: shear  = .false.
    LOGICAL :: extens = .false.

!   Number GMRES iterations/tolerance before termination
    INTEGER :: GMRES_it
    REAL(KIND = 8) :: GMRES_tol

    CONTAINS
    PROCEDURE :: bvAdvance => bvAdvanceInfo
    PROCEDURE :: GalInfo
    PROCEDURE :: GalOneInfo
    GENERIC :: Gal => GalInfo, GalOneInfo
    PROCEDURE :: init => initInfo
END TYPE sharedType

TYPE cellprops
    REAL(KIND = 8) lam, Ca, C, Eb, c0, int_pres
END TYPE cellprops
! -------------------------------------------------------------------------!
INTERFACE sharedType
    PROCEDURE :: newinfo
END INTERFACE sharedType

CONTAINS
!=============================================================================!
!================================= ROUTINES ==================================!
!=============================================================================!

FUNCTION newinfo(filein) RESULT (info)
    CHARACTER(len = *), INTENT(IN) :: filein
    CHARACTER(:), ALLOCATABLE :: cellchar
    INTEGER :: NCell
    TYPE(sharedType) :: info
    INTEGER fali, p, ic
    
    CALL READ_MFS(info%dt, filein, 'Time_step')
!   Coarse and fine grids    
    CALL READ_MFS(p, filein, 'Harmonic_order')
    info%p = p
    CALL READ_MFS(fali, filein, 'Refinement_factor')
    info%q = p*fali

!   Number of cells
    CALL READ_MFS(NCell, filein, 'Number_cells')
    info%NCell = NCell

!   Are cell-cell interactions included
    CALL READ_MFS(cellchar, filein, 'CellCell')
    IF(TRIM(cellchar) .eq. "Yes") THEN
        info%CellCell = .true.
    ELSE
        info%CellCell = .false.
    ENDIF

!   Periodic box length
    info%eye = 0D0
    FORALL(ic = 1:3) info%eye(ic,ic) = 1D0
    CALL READ_MFS(info%bvl, filein, 'Periodic_box_size')

    CALL READ_MFS(cellchar, filein, 'Flow_type') 
    IF(TRIM(cellchar) .eq. 'e') info%extens = .true.
    IF(TRIM(cellchar) .eq. 's') info%shear  = .true.

    info%GMRES_it = 25
    info%GMRES_tol = 5D-7
    
END FUNCTION newinfo

! -------------------------------------------------------------------------!
! Initialization
SUBROUTINE initInfo(info)
    CLASS(sharedType), INTENT(INOUT) :: info
    INTEGER p, q, nt, np, ntf, npf, ic, m, ind, it, n, im, im2
    REAL(KIND = 8), ALLOCATABLE :: cPt(:,:), ths(:,:), phs(:,:), thts(:), phis(:), xs(:), wg(:)
    REAL(KIND = 8) :: dphi

    p = info%p
    q = info%q
    IF(info%bvl .eq. 0) THEN
        info%periodic = .false.

!   Initialize periodic components of info
    ELSE
        info%periodic = .true.
!       Start with cube
        info%bv = info%eye*info%bvl
        info%tau = DOT(CROSS(info%bv(:,1), info%bv(:,2)), info%bv(:,3))
        
!       Reciprocal basis vectors
        info%kv(:,1) = 2D0*PI/info%tau*CROSS(info%bv(:,2), info%bv(:,3));
        info%kv(:,2) = 2D0*PI/info%tau*CROSS(info%bv(:,3), info%bv(:,1));
        info%kv(:,3) = 2D0*PI/info%tau*CROSS(info%bv(:,1), info%bv(:,2));

!       Some parameters for the Ewald sum
        info%gp = 2**5
!       Smoothing parameter (Dimensional quantity with unit 1/L, so should be based on box size)
        info%xi = 2.5D0*5D0/info%bvl!SQRT(PI)/(info%tau**(1D0/3D0))!

!       3D FFT Info
        info%LENSAVE = 2*info%gp + INT(LOG(REAL(info%gp))) + 4
        ALLOCATE(info%WSAVE(info%LENSAVE))
        CALL ZFFT1I(info%gp, info%WSAVE, info%LENSAVE, m)
    ENDIF

!   Make harmonics(order, # of derivs, if we'll rotate or not)
    info%Y = YType(p, 1, .true.)
    info%Yf = YType(q, 4, .true., p)
    
!   Stuff needed for calcs
    nt  = info%Y%nt
    np  = info%Y%np
    ntf = info%Yf%nt
    npf = info%Yf%np

!   Harmonics for the singular integration, slightly tricky
!   Construct the singular integration grid, with patch. Essentially 2 grids at once,
!   a fine patch with a sinh transform and a coarse one.
    IF(info%NCell .gt. 1 .or. info%periodic) THEN
        ALLOCATE(thts(nt + ntf), phis(np + npf), &
                 ths (nt + ntf, np + npf), phs(nt + ntf, np + npf))
!       Cutoff theta, can affect accuracy. Depends on spacing, but want consistent. Just hardcode for now
        info%thtc = pi/6D0 !!!

!       Coarser part away from near-singularity
        ALLOCATE(xs(nt), wg(nt))
        CALL lgwt(nt, xs, wg)
        thts(1:nt) = xs*(pi - info%thtc)/2D0 + (pi + info%thtc)/2D0
        DEALLOCATE(xs, wg)

!       Finer part near near-singularity
        ALLOCATE(xs(ntf), wg(ntf))
        CALL lgwt(ntf, xs, wg)
!       The integration rule is based on the separation distance. However, this changes and
!       we don't want to recalculate the grid every time. Instead, just choose a distance
!       where it's accurate (approx spacing/10 here)
        info%h = SQRT(PI)/10D0/nt
!       Sinh normalizing constant (to get -1, 1 range)
        info%k = -0.5D0*LOG(info%thtc/info%h &
                    + sqrt((info%thtc/info%h)**2D0 + 1D0));
        thts(nt + 1:nt + ntf) = info%h*SINH(info%k*(xs - 1D0))
        info%xsf = xs
        DEALLOCATE(xs, wg)

!       Calculate phis and construct mesh grid
        phis(1) = 0D0
        dphi = info%Y%dphi
        DO ic = 1,(np + npf)
            IF(ic .gt. 1) phis(ic) = phis(ic - 1) + dphi

            IF(ic .eq. np + 1) THEN
                dphi = info%Yf%dphi
                phis(ic) = 0D0
            ENDIF
            phs(:,ic) = phis(ic)
            ths(:,ic) = thts
        ENDDO

!       Create a bare harmonics object evaluated at this grid
        info%Ys = YType(ths, phs, p)
    ENDIF

    ALLOCATE(info%es(2*(p-1)+1, np), info%cPmn(p, 2*(p-1)+1, nt), cPt(nt, p*(p+1)/2))
    ALLOCATE(info%esf(2*(p-1)+1, npf + np), info%cPmnf(p, 2*(p-1)+1, ntf + nt))

!   Matrix size
    info%Nmat = 3*info%Y%nt*info%Y%np
    info%NmatT= info%Nmat*info%NCell

!   Exponential calculation part (coarse and fine)
    DO m = -(p-1),(p-1)
        ind = m + p
        info%es(ind,:) = EXP(ii*DBLE(m)*info%Y%phi)*info%Y%dphi
    ENDDO

!   Legendre polynomial calculation part (coarse and fine)
    cPt = Alegendre(p-1,COS(info%Y%tht))
    it = 0
    DO n = 0, p-1
        ind = n+1
        im = 0
        DO m = -(p-1),p-1
            im = im + 1
            IF(ABS(m) .gt. n) THEN
                info%cPmn(ind,im,:) = 0D0
            ELSEIF(m .le. 0) THEN
                it = it + 1
                IF(m.eq.-n) im2 = it
                info%cPmn(ind,im,:) = (-1D0)**m*cPt(:, im2 + abs(m))
            ELSE
                info%cPmn(ind,im,:) = (-1D0)**m*info%cPmn(ind, im - 2*m, :)
            ENDIF
        ENDDO
    ENDDO

    IF(info%NCell .gt. 1 .or. info%periodic) THEN
!       Fine part, essentially done on 2 grids
!       Technically the grid goes up to order q, but we only calculate up to p
        DO m = -(p-1),(p-1)
            ind = m + p
            info%esf(ind,:) = EXP(ii*DBLE(m)*info%Ys%ph(1,:))
        ENDDO

!       Manage the dphi
        DO ic = 1,np + npf
            IF(ic .le. np) THEN
                info%esf(:,ic) = info%esf(:,ic)*info%Y%dphi
            ELSE
                info%esf(:,ic) = info%esf(:,ic)*info%Yf%dphi
            ENDIF
        ENDDO
        DEALLOCATE(cPt)
!       Technically the grid goes up to order q, but we only calculate up to p
        ALLOCATE(cPt(nt + ntf, p*(p+1)/2))
        cPt = Alegendre(p-1,COS(info%Ys%th(:,1)))
        it = 0
        DO n = 0, p-1
            ind = n+1
            im = 0
            DO m = -(p-1),p-1
                im = im + 1
                IF(ABS(m) .gt. n) THEN
                    info%cPmnf(ind,im,:) = 0D0
                ELSEIF(m .le. 0) THEN
                    it = it + 1
                    IF(m.eq.-n) im2 = it
                    info%cPmnf(ind,im,:) = &
                        (-1D0)**m*cPt(:, im2 + abs(m))
                ELSE
                    info%cPmnf(ind,im,:) = &
                        (-1D0)**m*info%cPmnf(ind, im - 2*m, :)
                ENDIF
            ENDDO
        ENDDO
    ENDIF
    
!   Cell-cell interaction paramters (i.e. Morse and Lennard-Jones)
!   The below parameters seem ok, but perhaps could use a bit of tweaking.
    info%D = 0.03D0 ! 0.03 seems to be a good spot
    info%r0 = 0.1D0 ! Tough to know good spot. 0.1-0.15 probably
    info%Beta = 10D0! was at 1.5, but probably need around 10 so it decays fast
    info%epsi = 0.023D0

END SUBROUTINE initInfo

! -------------------------------------------------------------------------!
! Advances the basis vectors in time, reparameterizes if needed
SUBROUTINE bvAdvanceInfo(info, t)
    CLASS(sharedType), INTENT(INOUT) :: info
    REAL(KIND = 8), INTENT(IN), OPTIONAL :: t
    REAL(KIND = 8) :: lamp, bvt(3,3)

    info%bv = info%bv + MATMUL(info%dU,info%bv)*info%dt

!   temporary lees-edwards
    IF(info%shear) THEN
        IF(info%bv(1,3).gt.info%bvl) info%bv(1,3) = info%bv(1,3) - info%bvl

!   This is just for periodic strain
    ELSEIF(info%extens) THEN
        bvt = info%bv
        lamp = LOG((3D0 + sqrt(9D0 - 4D0))/2D0)
        IF(t/lamp + 0.5D0 - floor(t/lamp + 0.5D0) .le. info%dt/lamp) THEN
            info%bv(1,1) = bvt(1,1) +     bvt(1,3) 
            info%bv(1,3) = bvt(1,1) + 2D0*bvt(1,3)
            info%bv(3,1) = bvt(3,1) +     bvt(3,3)
            info%bv(3,3) = bvt(3,1) + 2d0*bvt(3,3)
        ENDIF
    ENDIF

!   Volume (shouldn't really change w/ no dilatation)
    info%tau = DOT(CROSS(info%bv(:,1), info%bv(:,2)), info%bv(:,3))

!   Reciprocal basis vectors
    info%kv(:,1) = 2D0*PI/info%tau*CROSS(info%bv(:,2), info%bv(:,3));
    info%kv(:,2) = 2D0*PI/info%tau*CROSS(info%bv(:,3), info%bv(:,1));
    info%kv(:,3) = 2D0*PI/info%tau*CROSS(info%bv(:,1), info%bv(:,2));

END SUBROUTINE bvAdvanceInfo

! -------------------------------------------------------------------------!
! Takes a matrix (3 X 3 X m X n X nt X np) and vector (3*nt*np) and performs a Galerkin
! projection with spherical hamronic functions to get a (3*m*n X 3*m*n) matrix/vec
SUBROUTINE GalInfo(info, Ai, b, A2, b2)
    CLASS(sharedType), TARGET, INTENT(IN) :: info
    COMPLEX(KIND = 8), ALLOCATABLE, INTENT(IN) :: Ai(:,:,:,:,:,:), b(:)
    COMPLEX(KIND = 8), ALLOCATABLE, INTENT(OUT) :: A2(:,:), b2(:)
    COMPLEX(KIND = 8), POINTER :: vcurn(:,:), es(:,:)
    REAL(KIND = 8), POINTER :: cPmn(:,:,:)
    TYPE(YType), POINTER :: Y
    TYPE(nmType), POINTER :: nm
    INTEGER :: it, n, m, im, col, m2, im2, i2, j2, n2, im3, row, ic, i, j
    COMPLEX(KIND = 8) :: At(3,3), bt(3), tmpsum(3,3)
    COMPLEX(KIND = 8), ALLOCATABLE :: Fi(:,:,:,:)

    es   => info%es
    cPmn => info%cPmn
    Y => info%Y

    ALLOCATE(A2(info%Nmat, info%Nmat), &
             b2(info%Nmat), &
             Fi(3,3, 2*(Y%p-1)+1, Y%nt))

!   Second integral: The outer loops go over the order and degree of the previous integrals
    it = 0
    DO n = 0,Y%p - 1
        nm => Y%nm(n+1)
        ! it = it + n
        DO m = -n,n
            im = m + Y%p
            it = it + 1
            col = 3*it - 2
            ! colm= col - 2*3*m

!           First loop: m2 (Galerkin order), theta, sum phi, 
            DO m2 = -(Y%p-1), Y%p-1
                im2 = m2 + Y%p
                DO i2 = 1,Y%nt
                    tmpsum = 0D0
                    DO j2 = 1,Y%np
                        tmpsum = tmpsum + Ai(1:3,1:3, im, n+1, i2, j2)*CONJG(es(im2,j2))
                    ENDDO
                    Fi(1:3,1:3,im2,i2) = tmpsum
                ENDDO
            ENDDO

!           Second loop: n2, m2, sum theta
            im2 = 0
            DO n2 = 0, Y%p-1
                DO m2 = -n2, n2
                    im3 = m2 + Y%p
                    im2 = im2 + 1
                    row = 3*im2 - 2
                    ! rowm= row - 2*3*m2
                    At = 0D0
                    ! At2 = 0D0
                    DO i2 = 1,Y%nt
                        At  = At  + Fi(1:3,1:3,im3, i2)*cPmn(n2+1,im3,i2)*Y%wg(i2)
                        ! At2 = At2 + CONJG(Fi(1:3,1:3,im3, i2))*cPmn(n2+1,im3,i2)*cell%Y%wg(i2)*(-1D0)**m2
                    ENDDO
                    A2(row:row+2,col:col+2) = At

!                   Can't get symmetry working as before, because I loop over cols, then rows now...
!                   Exploit symmetry (This and the calculation of At2 are a little overkill,
!                   but it's finnicky to get the right if statements so I'll just leave them)
                    ! A2(rowm:rowm+2, col :col +2) = At2
                    ! A2(row :row +2, colm:colm+2) = (-1D0)**(m - m2)*CONJG(At2)
                    ! A2(rowm:rowm+2, colm:colm+2) = (-1D0)**(m + m2)*CONJG(At)
                ENDDO
            ENDDO

            ic = 0
!           Loop over integration points to calc integral
            bt = 0D0
            vcurn => nm%v(m + n + 1,:,:)
            DO i =1,Y%nt
                DO j = 1,Y%np
                ic = ic+1

!               Intg. b (essentially forward transform of RHS!)
                bt = bt + b(3*ic-2:3*ic)*CONJG(vcurn(i,j)) &
                    *Y%wg(i)*Y%dphi
                ENDDO
            ENDDO
            b2(col:col+2) = bt

        ENDDO
    ENDDO
    
    vcurn => NULL()
    es => NULL()
    cPmn => NULL()
    Y => NULL()
    nm => NULL()
END SUBROUTINE GalInfo
    
! -------------------------------------------------------------------------!
! Takes a (real) vector (3*nt*np) and performs a Galerkin
! projection with spherical hamronic functions to get a (3*m*n) matrix/vec
SUBROUTINE GalOneInfo(info, b, b2, itt1, itt2)
    CLASS(sharedType), TARGET, INTENT(IN) :: info
    REAL(KIND = 8), ALLOCATABLE, INTENT(IN) :: b(:)
    COMPLEX(KIND = 8), ALLOCATABLE, INTENT(OUT) ::b2(:)
    COMPLEX(KIND = 8), POINTER :: vcurn(:,:)
    COMPLEX(KIND = 8) :: bt(3)
    TYPE(YType), POINTER :: Y
    TYPE(nmType), POINTER :: nm
    INTEGER, INTENT(IN), OPTIONAL :: itt1, itt2
    INTEGER :: it, n, m, im, col, ic, i, j, it1, it2

    Y => info%Y

    IF(PRESENT(itt1)) THEN
        it1 = itt1
    ELSE
        it1 = 1
    ENDIF

    IF(PRESENT(itt2)) THEN
        it2 = itt2
    ELSE
        it2 = Y%nt
    ENDIF

    ALLOCATE(&
             b2(3*(Y%p + 1)*(Y%p + 1)))
    b2 = 0D0

!   Second integral: The outer loops go over the order and degree of the previous integrals
    it = 0
    DO n = 0,Y%p - 1
        nm => Y%nm(n+1)
        ! it = it + n
        DO m = -n,n !!!!!!!!! Symmetry
            im = m + Y%p
            it = it + 1
            col = 3*it - 2

            ic = (it1 - 1)*Y%np
!           Loop over integration points to calc integral
            bt = 0D0
            vcurn => nm%v(m + n + 1,:,:)
!           Still parallel along
            DO i = it1, it2
                DO j = 1,Y%np
                ic = ic+1

!               Intg. b (essentially forward transform of RHS!)
                bt = bt + b(3*ic-2:3*ic)*CONJG(vcurn(i,j)) &
                    *Y%wg(i)*Y%dphi
                ENDDO
            ENDDO
            b2(col:col+2) = bt

        ENDDO
    ENDDO
    
    vcurn => NULL()
    Y => NULL()
    nm => NULL()
END SUBROUTINE GalOneInfo
    
END MODULE SHAREDMOD