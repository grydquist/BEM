      SUBROUTINE infog2l( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW,
     $                    MYCOL, LRINDX, LCINDX, RSRC, CSRC )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            CSRC, GCINDX, GRINDX, LRINDX, LCINDX, MYCOL,
     $                   myrow, npcol, nprow, rsrc
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  INFOG2L computes the starting local indexes LRINDX, LCINDX corres-
*  ponding to the distributed submatrix starting globally at the entry
*  pointed by GRINDX, GCINDX. This routine returns the coordinates in
*  the grid of the process owning the matrix entry of global indexes
*  GRINDX, GCINDX, namely RSRC and CSRC.
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
*  GRINDX    (global input) INTEGER
*            The global row starting index of the submatrix.
*
*  GCINDX    (global input) INTEGER
*            The global column starting index of the submatrix.
*
*  DESC      (input) INTEGER array of dimension DLEN_.
*            The array descriptor for the underlying distributed matrix.
*
*  NPROW     (global input) INTEGER
*            The total number of process rows over which the distributed
*            matrix is distributed.
*
*  NPCOL     (global input) INTEGER
*            The total number of process columns over which the
*            distributed matrix is distributed.
*
*  MYROW     (local input) INTEGER
*            The row coordinate of the process calling this routine.
*
*  MYCOL     (local input) INTEGER
*            The column coordinate of the process calling this routine.
*
*  LRINDX    (local output) INTEGER
*            The local rows starting index of the submatrix.
*
*  LCINDX    (local output) INTEGER
*            The local columns starting index of the submatrix.
*
*  RSRC      (global output) INTEGER
*            The row coordinate of the process that possesses the first
*            row and column of the submatrix.
*
*  CSRC      (global output) INTEGER
*            The column coordinate of the process that possesses the
*            first row and column of the submatrix.
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
      INTEGER            CBLK, GCCPY, GRCPY, RBLK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          mod
*     ..
*     .. Executable Statements ..
*
      grcpy = grindx-1
      gccpy = gcindx-1
*
      rblk = grcpy / desc(mb_)
      cblk = gccpy / desc(nb_)
      rsrc = mod( rblk + desc(rsrc_), nprow )
      csrc = mod( cblk + desc(csrc_), npcol )
*
      lrindx = ( rblk / nprow + 1 ) * desc(mb_) + 1
      lcindx = ( cblk / npcol + 1 ) * desc(nb_) + 1
*
      IF( mod( myrow+nprow-desc(rsrc_), nprow ) .GE.
     $    mod( rblk, nprow ) ) THEN
         IF( myrow.EQ.rsrc )
     $      lrindx = lrindx + mod( grcpy, desc(mb_) )
         lrindx = lrindx - desc(mb_)
      END IF
*
      IF( mod( mycol+npcol-desc(csrc_), npcol ) .GE.
     $    mod( cblk, npcol ) ) THEN
         IF( mycol.EQ.csrc )
     $      lcindx = lcindx + mod( gccpy, desc(nb_) )
         lcindx = lcindx - desc(nb_)
      END IF
*
      RETURN
*
*     End of INFOG2L
*
       END