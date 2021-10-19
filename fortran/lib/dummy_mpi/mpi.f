!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!      
!--------------------------------------------------------------------
!      
!     These are only dummy functions, created for sequential code
!

      SUBROUTINE MPI_ALLREDUCE(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7
      END SUBROUTINE MPI_ALLREDUCE

      SUBROUTINE MPI_BCAST(arg1, arg2, arg3, arg4, arg5, arg6)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6
      END SUBROUTINE MPI_BCAST

      SUBROUTINE MPI_GATHERV(arg1, arg2, arg3, arg4, arg5, arg6, arg7,  &
     &   arg8, arg9, arg10)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9,  &
     &   arg10
      END SUBROUTINE MPI_GATHERV

      SUBROUTINE MPI_SCATTERV(arg1, arg2, arg3, arg4, arg5, arg6, arg7, &
     &   arg8, arg9, arg10)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9,  &
     &   arg10
      END SUBROUTINE MPI_SCATTERV

      SUBROUTINE MPI_SCATTER(arg1, arg2, arg3, arg4, arg5, arg6, arg7, 
     &   arg8, arg9)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
      END SUBROUTINE MPI_SCATTER

      SUBROUTINE MPI_ALLGATHERV(arg1, arg2, arg3, arg4, arg5, arg6,     &
     &   arg7, arg8, arg9)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
      END SUBROUTINE MPI_ALLGATHERV

      SUBROUTINE MPI_GATHER(arg1, arg2, arg3, arg4, arg5, arg6,     
     &   arg7, arg8, arg9)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
      END SUBROUTINE MPI_GATHER

      SUBROUTINE MPI_ALLGATHER(arg1, arg2, arg3, arg4, arg5, arg6, arg7,&
     &   arg8)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      END SUBROUTINE MPI_ALLGATHER
      
      SUBROUTINE MPI_ABORT(arg1, arg2, arg3)
         INTEGER arg1, arg2, arg3
      END SUBROUTINE MPI_ABORT
      
      SUBROUTINE MPI_FINALIZE(arg1)
         INTEGER arg1
      END SUBROUTINE MPI_FINALIZE

      SUBROUTINE MPI_INIT(arg1)
         INTEGER arg1
      END SUBROUTINE MPI_INIT

      SUBROUTINE MPI_COMM_RANK(arg1, arg2, arg3)
         INTEGER arg1, arg2, arg3
      END SUBROUTINE MPI_COMM_RANK

      SUBROUTINE MPI_COMM_SIZE(arg1, arg2, arg3)
         INTEGER arg1, arg2, arg3
      END SUBROUTINE MPI_COMM_SIZE
         
      SUBROUTINE MPI_IRECV(arg1,arg2,arg3, arg4, arg5, arg6, arg7, arg8)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      END SUBROUTINE MPI_IRECV
         
      SUBROUTINE MPI_ISEND(arg1,arg2,arg3, arg4, arg5, arg6, arg7, arg8)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      END SUBROUTINE MPI_ISEND
         
      SUBROUTINE MPI_WAIT(arg1, arg2, arg3)
         INTEGER arg1, arg2, arg3
      END SUBROUTINE MPI_WAIT
      
      SUBROUTINE MPI_RECV(arg1,arg2, arg3, arg4, arg5, arg6, arg7, arg8)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      END SUBROUTINE MPI_RECV
         
      SUBROUTINE MPI_SEND(arg1,arg2, arg3, arg4, arg5, arg6, arg7, arg8)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      END SUBROUTINE MPI_SEND

      SUBROUTINE MPI_TYPE_SIZE(arg1, arg2, arg3)
         INTEGER arg1, arg2, arg3
      END SUBROUTINE MPI_TYPE_SIZE

      SUBROUTINE MPI_COMM_GROUP(arg1, arg2, arg3)
         INTEGER arg1, arg2, arg3
      END SUBROUTINE MPI_COMM_GROUP

      SUBROUTINE MPI_GROUP_INCL(arg1, arg2, arg3, arg4, arg5)
         INTEGER arg1, arg2, arg3, arg4, arg5
      END SUBROUTINE MPI_GROUP_INCL

      SUBROUTINE MPI_GROUP_FREE(arg1, arg2)
         INTEGER arg1, arg2
      END SUBROUTINE MPI_GROUP_FREE

      SUBROUTINE MPI_COMM_CREATE(arg1, arg2, arg3, arg4)
         INTEGER arg1, arg2, arg3, arg4
      END SUBROUTINE MPI_COMM_CREATE

      SUBROUTINE MPI_BARRIER(arg1, arg2)
         INTEGER arg1, arg2
      END SUBROUTINE MPI_BARRIER

      SUBROUTINE MPI_INITIALIZED(arg1, arg2)
         INTEGER arg1, arg2
      END SUBROUTINE MPI_INITIALIZED

      SUBROUTINE SPLIT(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
         INTEGER arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
      END SUBROUTINE SPLIT

