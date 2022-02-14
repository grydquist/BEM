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
    REAL(KIND = 8) :: eye(3,3), bv(3,3), tau, bvl

!   Cell-cell parameters
    REAL(KIND = 8) :: epsi, r0, D, Beta
    LOGICAL :: CellCell

    CONTAINS
    PROCEDURE :: bvAdvance => bvAdvanceInfo
END TYPE sharedType

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
    INTEGER fali, p, nt, np, ntf, npf, ic, m, ind, it, n, im, im2
    REAL(KIND = 8), ALLOCATABLE :: cPt(:,:), ths(:,:), phs(:,:), thts(:), phis(:), xs(:), wg(:)
    REAL(KIND = 8) :: dphi
    
    CALL READ_MFS(info%dt, filein, 'Time_step')
!   Coarse and fine grids    
    CALL READ_MFS(p, filein, 'Harmonic_order')
    info%p = p
    CALL READ_MFS(fali, filein, 'Refinement_factor')
    info%q = p*fali

!   Number of cells
    CALL READ_MFS(NCell, filein, 'Number_cells')

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
    IF(info%bvl .eq. 0) THEN
        info%periodic = .false.
    ELSE
        info%periodic = .true.
!       Start with cube
        info%bv = info%eye*info%bvl
    ENDIF

!   Make harmonics(order, # of derivs, if we'll rotate or not)
    info%Y = YType(p, 1, .true.)
    info%Yf = YType(p*fali, 4, .true., p)
    
!   Stuff needed for calcs
    nt  = info%Y%nt
    np  = info%Y%np
    ntf = info%Yf%nt
    npf = info%Yf%np

!   Harmonics for the singular integration, slightly tricky
!   Construct the singular integration grid, with patch. Essentially 2 grids at once,
!   a fine patch with a sinh transform and a coarse one.
    IF(NCell.gt.1) THEN
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
    info%Nmat = 3*(info%Y%p)*(info%Y%p)
    info%NmatT= info%Nmat*NCell

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

    IF(NCell .gt. 1) THEN
!   Fine part, essentially done on 2 grids
!   Technically the grid goes up to order q, but we only calculate up to p
    DO m = -(p-1),(p-1)
        ind = m + p
        info%esf(ind,:) = EXP(ii*DBLE(m)*info%Ys%ph(1,:))
    ENDDO

!   Manage the dphi
    DO ic = 1,np + npf
        IF(ic .le. np) THEN
            info%esf(:,ic) = info%esf(:,ic)*info%Y%dphi
        ELSE
            info%esf(:,ic) = info%esf(:,ic)*info%Yf%dphi
        ENDIF
    ENDDO
    DEALLOCATE(cPt)
!   Technically the grid goes up to order q, but we only calculate up to p
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
                info%cPmnf(ind,im,:) = (-1D0)**m*cPt(:, im2 + abs(m))
            ELSE
                info%cPmnf(ind,im,:) = (-1D0)**m*info%cPmnf(ind, im - 2*m, :)
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

END FUNCTION newinfo

! -------------------------------------------------------------------------!
! Advances the basis vectors in time, reparameterizes if needed
SUBROUTINE bvAdvanceInfo(info)
    CLASS(sharedType), INTENT(INOUT) :: info

    info%bv = info%bv + MATMUL(info%dU,info%bv)*info%dt

!   temporary lees-edwards
    IF(info%bv(1,3).gt.info%bvl) info%bv(1,3) = info%bv(1,3) - info%bvl

END SUBROUTINE bvAdvanceInfo
    
END MODULE SHAREDMOD