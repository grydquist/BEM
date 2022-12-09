      SUBROUTINE pzlaswp( DIREC, ROWCOL, N, A, IA, JA, DESCA, K1, K2,
     $                    IPIV )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
       CHARACTER          DIREC, ROWCOL
       INTEGER            IA, JA, K1, K2, N
*     ..
*     .. Array Arguments ..
       INTEGER            DESCA( * ), IPIV( * )
       COMPLEX*16         A( * )
*     ..
*
*  Purpose:
*  ========
*
*  PZLASWP performs a series of row or column interchanges on
*  the distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1).  One
*  interchange is initiated for each of rows or columns K1 trough K2 of
*  sub( A ). This routine assumes that the pivoting information has
*  already been broadcast along the process row or column.
*  Also note that this routine will only work for K1-K2 being in the
*  same MB (or NB) block.  If you want to pivot a full matrix, use
*  PZLAPIV.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  DIREC   (global input) CHARACTER
*          Specifies in which order the permutation is applied:
*          = 'F' (Forward)
*          = 'B' (Backward)
*
*  ROWCOL  (global input) CHARACTER
*          Specifies if the rows or columns are permuted:
*          = 'R' (Rows)
*          = 'C' (Columns)
*
*  N       (global input) INTEGER
*          If ROWCOL = 'R', the length of the rows of the distributed
*          matrix A(*,JA:JA+N-1) to be permuted;
*          If ROWCOL = 'C', the length of the columns of the distributed
*          matrix A(IA:IA+N-1,*) to be permuted.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, * ).
*          On entry, this array contains the local pieces of the distri-
*          buted matrix to which the row/columns interchanges will be
*          applied. On exit the permuted distributed matrix.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  K1      (global input) INTEGER
*          The first element of IPIV for which a row or column inter-
*          change will be done.
*
*  K2      (global input) INTEGER
*          The last element of IPIV for which a row or column inter-
*          change will be done.
*
*  IPIV    (local input) INTEGER array, dimension LOCr(M_A)+MB_A for
*          row pivoting and LOCc(N_A)+NB_A for column pivoting.  This
*          array is tied to the matrix A, IPIV(K) = L implies rows
*          (or columns) K and L are to be interchanged.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   lld_, mb_, m_, nb_, n_, rsrc_
       parameter( block_cyclic_2d = 1, dlen_ = 9, dtype_ = 1,
     $                     ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6,
     $                     rsrc_ = 7, csrc_ = 8, lld_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICURCOL, ICURROW, IIA, IP, J, JJA, JP,
     $                   mycol, myrow, npcol, nprow
*     ..
*     .. External Subroutines ..
      EXTERNAL           blacs_gridinfo, infog2l, pzswap
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
      CALL blacs_gridinfo( desca( ctxt_ ), nprow, npcol, myrow, mycol )
*
      IF( lsame( rowcol, 'R' ) ) THEN
         IF( lsame( direc, 'F' ) ) THEN
            CALL infog2l( k1, ja, desca, nprow, npcol, myrow, mycol,
     $                    iia, jja, icurrow, icurcol )
            DO 10 i = k1, k2
               ip = ipiv( iia+i-k1 )
               IF( ip.NE.i )
     $            CALL pzswap( n, a, i, ja, desca, desca( m_ ), a, ip,
     $                         ja, desca, desca( m_ ) )
   10       CONTINUE
         ELSE
            CALL infog2l( k2, ja, desca, nprow, npcol, myrow, mycol,
     $                    iia, jja, icurrow, icurcol )
           DO 20 i = k2, k1, -1
               ip = ipiv( iia+i-k1 )
               IF( ip.NE.i )
     $            CALL pzswap( n, a, i, ja, desca, desca( m_ ), a, ip,
     $                         ja, desca, desca( m_ ) )
   20       CONTINUE
         END IF
      ELSE
         IF( lsame( direc, 'F' ) ) THEN
            CALL infog2l( ia, k1, desca, nprow, npcol, myrow, mycol,
     $                    iia, jja, icurrow, icurcol )
            DO 30 j = k1, k2
               jp = ipiv( jja+j-k1 )
               IF( jp.NE.j )
     $            CALL pzswap( n, a, ia, j, desca, 1, a, ia, jp,
     $                         desca, 1 )
   30       CONTINUE
         ELSE
            CALL infog2l( ia, k2, desca, nprow, npcol, myrow, mycol,
     $                    iia, jja, icurrow, icurcol )
            DO 40 j = k2, k1, -1
               jp = ipiv( jja+j-k1 )
               IF( jp.NE.j )
     $            CALL pzswap( n, a, ia, j, desca, 1, a, ia, jp,
     $                         desca, 1 )
   40       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End PZLASWP
*
      END