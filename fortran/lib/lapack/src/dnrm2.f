      DOUBLE PRECISION FUNCTION dnrm2(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER incx,n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION x(*)
*     ..
*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION one,zero
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION absxi,norm,scale,ssq
      INTEGER ix
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
*     ..
      IF (n.LT.1 .OR. incx.LT.1) THEN
          norm = zero
      ELSE IF (n.EQ.1) THEN
          norm = abs(x(1))
      ELSE
          scale = zero
          ssq = one
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 ix = 1,1 + (n-1)*incx,incx
              IF (x(ix).NE.zero) THEN
                  absxi = abs(x(ix))
                  IF (scale.LT.absxi) THEN
                      ssq = one + ssq* (scale/absxi)**2
                      scale = absxi
                  ELSE
                      ssq = ssq + (absxi/scale)**2
                  END IF
              END IF
   10     CONTINUE
          norm = scale*sqrt(ssq)
      END IF
*
      dnrm2 = norm
      RETURN
*
*     End of DNRM2.
*
      END