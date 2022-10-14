!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Communicator related module. For usage see TEST_CMMOD at the end
!     of this file.
!      
!--------------------------------------------------------------------

      MODULE CMMOD
      IMPLICIT NONE
      INCLUDE "mpif.h"

!     Size of blocks for openMP communications
      INTEGER, PARAMETER :: mpBs = 1000
!     master is assumed to have zero ID
      INTEGER, PARAMETER :: master = 0
!     Abstracted MPI names
      INTEGER, PARAMETER :: mplog  = MPI_LOGICAL
      INTEGER, PARAMETER :: mpint  = MPI_INTEGER
      INTEGER, PARAMETER :: mpreal = MPI_DOUBLE_PRECISION
      INTEGER, PARAMETER :: mpchar = MPI_CHARACTER
      INTEGER, PARAMETER :: mpcplx = MPI_DOUBLE_COMPLEX

      TYPE cmType
         PRIVATE
!        Communicator handle
         INTEGER :: cHndl = MPI_COMM_NULL
!        Processors ID
         INTEGER :: taskId = 0
!        Number of openMP threads in this cm
         INTEGER nThreads
!        Number of processors
         INTEGER :: nProcs = 1
      CONTAINS 
!        Returns commu handle
         PROCEDURE :: com => comCm
!        Returns processor ID 
         PROCEDURE :: id => idCm
!        Returns number of processors
         PROCEDURE :: np => npCm
!        Returns number of threads
         PROCEDURE :: nT => ntCm
!        Returns processor ID in fortran indexing
         PROCEDURE :: tF => tfCm
!        Returns true if this is master
         PROCEDURE :: mas => masCm
!        Returns true if this is slave
         PROCEDURE :: slv => slvCm
!        Returns true if there is only one processor
         PROCEDURE :: seq => seqCm
!        Forces the MPI communicator to abort
         PROCEDURE :: fStop => fStopCm
!        Finalizes the communicator
         PROCEDURE :: finalize => finalizeCm
!        Broadcasting scaler/vector of logic/integer/real/character
         PROCEDURE :: BCASTLS
         PROCEDURE :: BCASTLV
         PROCEDURE :: BCASTIS
         PROCEDURE :: BCASTIV
         PROCEDURE :: BCASTIM
         PROCEDURE :: BCASTRS
         PROCEDURE :: BCASTRV
         PROCEDURE :: BCASTRM
         PROCEDURE :: BCASTSS
         PROCEDURE :: BCASTSV
         PROCEDURE :: BCASTCS
         PROCEDURE :: BCASTCV
         GENERIC :: bcast => BCASTLS, BCASTLV, BCASTIS, BCASTIV,&
           BCASTIM, BCASTRS, BCASTRV, BCASTRM, BCASTSS, BCASTSV, &
           BCASTCS, BCASTCV
!        Blocking MPI send
         PROCEDURE :: SENDRV
         PROCEDURE :: SENDCV
         PROCEDURE :: SENDCA
         GENERIC :: send => SENDRV, SENDCV, SENDCA
!        Blocking MPI recv
         PROCEDURE :: RECVRV
         PROCEDURE :: RECVCV
         PROCEDURE :: RECVCA
         GENERIC :: recv => RECVRV, RECVCV, RECVCA
!        Non-blocking MPI send
         PROCEDURE :: isend => ISENDRV
!        Non-blocking MPI recv
         PROCEDURE :: irecv => IRECVRV
!        Doing MPI wait on a set of requests
         PROCEDURE, NOPASS :: WAITS
         PROCEDURE :: WAITV
         GENERIC :: wait => WAITS, WAITV
!        Doing MPI reduce on a set of data types
         PROCEDURE :: REDUCEIS
         PROCEDURE :: REDUCEIV
         PROCEDURE :: REDUCERS
         PROCEDURE :: REDUCERV
         PROCEDURE :: REDUCERM
         PROCEDURE :: REDUCECS
         PROCEDURE :: REDUCECV
         GENERIC :: reduce => REDUCEIS, REDUCEIV, REDUCERS, REDUCERV, &
           REDUCERM, REDUCECS, REDUCECV
!        Doing MPI_GATHER on different data types
         PROCEDURE :: GATHERIS
         PROCEDURE :: GATHERIV
         PROCEDURE :: GATHERIM
         PROCEDURE :: GATHERRS
         PROCEDURE :: GATHERRV
         PROCEDURE :: GATHERRM
         PROCEDURE :: GATHERCS
         PROCEDURE :: GATHERCV
         GENERIC :: gather => GATHERIS, GATHERIV, GATHERIM, GATHERRS, &
           GATHERRV, GATHERRM, GATHERCS, GATHERCV
!        Doing MPI_SCATTER on different data types
         PROCEDURE :: SCATRIS
         PROCEDURE :: SCATRIV
         PROCEDURE :: SCATRIM
         PROCEDURE :: SCATRRS
         PROCEDURE :: SCATRRV
         PROCEDURE :: SCATRRM
         GENERIC :: scatter => SCATRIS, SCATRIV, SCATRIM, SCATRRS, &
           SCATRRV, SCATRRM
!        Finding global of a field, provided a specific mapping
         PROCEDURE :: GLOBALIS
         PROCEDURE :: GLOBALRS
         PROCEDURE :: GLOBALRV
         GENERIC :: global => GLOBALIS, GLOBALRS, GLOBALRV
!        Finding local of a field, provided a specific mapping
         PROCEDURE :: LOCALIS
         PROCEDURE :: LOCALRS
         PROCEDURE :: LOCALRV
         GENERIC :: local => LOCALIS, LOCALRS, LOCALRV
!        To read and write in parallel to a binary file
         PROCEDURE :: WRITEIS
         PROCEDURE :: WRITERS
         PROCEDURE :: WRITEIV
         PROCEDURE :: WRITERV
         GENERIC :: write => WRITEIS, WRITERS, WRITEIV, WRITERV
         PROCEDURE :: READIS
         PROCEDURE :: READRS
         PROCEDURE :: READIV
         PROCEDURE :: READRV
         GENERIC :: read => READIS, READRS, READIV, READRV
!        Calls MPI barrier
         PROCEDURE :: barrier => barrierCm
!        Creating a communication handle from a set of processors
         PROCEDURE :: create => createCm
      END TYPE cmType

      INTERFACE cmType
         PROCEDURE :: newCm
      END INTERFACE
!---------------------------------------------------------------------
!     Tag upper bound
      INTEGER, PRIVATE :: tagUB = 32767
!     Sizes of various data types
      INTEGER, PROTECTED :: mplog_size, mpint_size, mpreal_size,&
        mpchar_size, mpcplx_size
!     A generic communicator for everyone use
      TYPE(cmType) :: cm

      CONTAINS 
!####################################################################

!     IMPLEMENTATION

!####################################################################
      FUNCTION newCm(comHandle) RESULT(cm)
      INTEGER, INTENT(IN), OPTIONAL :: comHandle
      TYPE(cmType) :: cm

      LOGICAL flag
      INTEGER ierr

      CALL MPI_INITIALIZED(flag,ierr)
      IF (.NOT.flag) CALL MPI_INIT(ierr)

      cm%cHndl = MPI_COMM_WORLD
      IF (PRESENT(comHandle)) cm%cHndl = comHandle

      CALL MPI_COMM_SIZE(cm%cHndl, cm%nProcs, ierr)
      CALL MPI_COMM_RANK(cm%cHndl, cm%taskId, ierr)

!      CALL MPI_COMM_GET_ATTR(cm%cHndl, MPI_TAG_UB, i, flag, ierr)
!      IF (flag) tagUB = i

!$OMP PARALLEL 
      cm%nThreads = 1!OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

      CALL MPI_TYPE_SIZE(mplog,  mplog_size,  ierr)
      CALL MPI_TYPE_SIZE(mpint,  mpint_size,  ierr)
      CALL MPI_TYPE_SIZE(mpreal, mpreal_size, ierr)
      CALL MPI_TYPE_SIZE(mpchar, mpchar_size, ierr)
      CALL MPI_TYPE_SIZE(mpcplx, mpcplx_size, ierr)

      RETURN
      END FUNCTION newCm
!--------------------------------------------------------------------
      FUNCTION comCm(cm) RESULT(com)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER com

      com = cm%cHndl

      RETURN
      END FUNCTION comCm
!--------------------------------------------------------------------
      FUNCTION idCm(cm) RESULT(taskId)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER taskId

      taskId = cm%taskId

      RETURN
      END FUNCTION idCm
!--------------------------------------------------------------------
      FUNCTION npCm(cm) RESULT(n)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER n

      n = cm%nProcs

      RETURN
      END FUNCTION npCm
!--------------------------------------------------------------------
      FUNCTION ntCm(cm) RESULT(nt)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER nt

      nt = cm%nThreads

      RETURN
      END FUNCTION ntCm
!--------------------------------------------------------------------
      FUNCTION tfCm(cm) RESULT(tf)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER tf

      tf = cm%taskId + 1

      RETURN
      END FUNCTION tfCm
!--------------------------------------------------------------------
      FUNCTION masCm(cm) RESULT(mas)
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL mas

      IF (cm%taskId .NE. master) THEN
         mas = .FALSE.
      ELSE
         mas = .TRUE.
      END IF

      RETURN
      END FUNCTION masCm
!--------------------------------------------------------------------
      FUNCTION slvCm(cm) RESULT(slv)
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL slv

      IF (cm%taskId .NE. master) THEN
         slv = .TRUE.
      ELSE
         slv = .FALSE.
      END IF

      RETURN
      END FUNCTION slvCm
!--------------------------------------------------------------------
      FUNCTION seqCm(cm) RESULT(seq)
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL seq

      seq = .FALSE.
      IF (cm%nProcs .EQ. 1) seq = .TRUE.

      RETURN
      END FUNCTION seqCm
!--------------------------------------------------------------------
      SUBROUTINE fStopCm(cm)
      CLASS(cmType), INTENT(IN) :: cm

      LOGICAL ierr
      
      CALL MPI_ABORT(cm, MPI_ERR_OTHER, ierr)

      RETURN
      END SUBROUTINE fStopCm

!####################################################################
      
      SUBROUTINE BCASTLS(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL, INTENT(INOUT) :: u

      INTEGER ierr

      CALL MPI_BCAST(u, 1, mplog, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTLS
!--------------------------------------------------------------------
      SUBROUTINE BCASTLV(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL, INTENT(INOUT) :: u(:)

      INTEGER m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mplog, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTLV

!--------------------------------------------------------------------
      SUBROUTINE BCASTIS(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(INOUT) :: u

      INTEGER ierr

      CALL MPI_BCAST(u, 1, mpint, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTIS
!--------------------------------------------------------------------
      SUBROUTINE BCASTIV(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(INOUT) :: u(:)

      INTEGER m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpint, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTIV
!--------------------------------------------------------------------
      SUBROUTINE BCASTIM(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(INOUT) :: u(:,:)

      INTEGER m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpint, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTIM
!--------------------------------------------------------------------
      SUBROUTINE BCASTRS(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(INOUT) :: u

      INTEGER ierr

      CALL MPI_BCAST(u, 1, mpreal, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTRS
!--------------------------------------------------------------------
      SUBROUTINE BCASTRV(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(INOUT) :: u(:)

      INTEGER m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpreal, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTRV
!--------------------------------------------------------------------
      SUBROUTINE BCASTRM(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(INOUT) :: u(:,:)

      INTEGER m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpreal, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTRM
!--------------------------------------------------------------------
      SUBROUTINE BCASTSS(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      CHARACTER(LEN=*), INTENT(INOUT) :: u

      INTEGER l, ierr

      l = LEN(u)
      CALL MPI_BCAST(u, l, mpchar, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTSS
!--------------------------------------------------------------------
      SUBROUTINE BCASTSV(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      CHARACTER(LEN=*), INTENT(INOUT) :: u(:)

      INTEGER m, l, ierr
      
      m = SIZE(u)
      l = LEN(u)
      CALL MPI_BCAST(u, l*m, mpchar, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTSV
!--------------------------------------------------------------------
      SUBROUTINE BCASTCS(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(INOUT) :: u

      INTEGER ierr

      CALL MPI_BCAST(u, 1, mpcplx, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTCS
!--------------------------------------------------------------------
      SUBROUTINE BCASTCV(cm, u)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(INOUT) :: u(:)

      INTEGER m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpcplx, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTCV

!####################################################################
      
      SUBROUTINE SENDRV(cm, u, to, tag)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN) :: to
      INTEGER, INTENT(IN), OPTIONAL :: tag

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + cm%id(), tagUB)

      IF (to .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_SEND(u, m, mpreal, to, ftag, cm%com(), ierr)

      RETURN
      END SUBROUTINE SENDRV
!--------------------------------------------------------------------
      SUBROUTINE SENDCV(cm, u, to, tag)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN) :: to
      INTEGER, INTENT(IN), OPTIONAL :: tag

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + cm%id(), tagUB)

      IF (to .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_SEND(u, m, mpcplx, to, ftag, cm%com(), ierr)

      RETURN
      END SUBROUTINE SENDCV
!--------------------------------------------------------------------
      SUBROUTINE SENDCA(cm, u, to, tag)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(IN) :: u(:,:,:,:,:,:)
      INTEGER, INTENT(IN) :: to
      INTEGER, INTENT(IN), OPTIONAL :: tag

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + cm%id(), tagUB)

      IF (to .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_SEND(u, m, mpcplx, to, ftag, cm%com(), ierr)

      RETURN
      END SUBROUTINE SENDCA
!--------------------------------------------------------------------
      SUBROUTINE RECVRV(cm, u, from, tag)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(OUT) :: u(:)
      INTEGER, INTENT(IN) :: from
      INTEGER, INTENT(IN), OPTIONAL :: tag

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + from, tagUB)
      
      IF (from .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_RECV(u, m, mpreal, from, ftag, cm%com(), &
         MPI_STATUS_IGNORE, ierr)

      RETURN
      END SUBROUTINE RECVRV
!--------------------------------------------------------------------
      SUBROUTINE RECVCV(cm, u, from, tag)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(OUT) :: u(:)
      INTEGER, INTENT(IN) :: from
      INTEGER, INTENT(IN), OPTIONAL :: tag

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + from, tagUB)
      
      IF (from .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_RECV(u, m, mpcplx, from, ftag, cm%com(), &
         MPI_STATUS_IGNORE, ierr)

      RETURN
      END SUBROUTINE RECVCV
!--------------------------------------------------------------------
      SUBROUTINE RECVCA(cm, u, from, tag)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(OUT) :: u(:,:,:,:,:,:)
      INTEGER, INTENT(IN) :: from
      INTEGER, INTENT(IN), OPTIONAL :: tag

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + from, tagUB)
      
      IF (from .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_RECV(u, m, mpcplx, from, ftag, cm%com(), &
         MPI_STATUS_IGNORE, ierr)

      RETURN
      END SUBROUTINE RECVCA
!--------------------------------------------------------------------
      FUNCTION ISENDRV(cm, u, to, tag) RESULT(req)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN) :: to
      INTEGER, INTENT(IN), OPTIONAL :: tag
      INTEGER req

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + cm%id(), tagUB)

      req = MPI_REQUEST_NULL
      IF (to .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_ISEND(u, m, mpreal, to, ftag, cm%com(), req, ierr)

      RETURN
      END FUNCTION ISENDRV
!--------------------------------------------------------------------
      FUNCTION IRECVRV(cm, u, from, tag) RESULT(req)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(OUT) :: u(:)
      INTEGER, INTENT(IN) :: from
      INTEGER, INTENT(IN), OPTIONAL :: tag
      INTEGER req

      INTEGER m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag
      ftag = MOD(ftag + from, tagUB)
      
      req = MPI_REQUEST_NULL
      IF (from .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_IRECV(u, m, mpreal, from, ftag, cm%com(), req, ierr)

      RETURN
      END FUNCTION IRECVRV
!--------------------------------------------------------------------
      SUBROUTINE WAITS(iReq)
      INTEGER, INTENT(IN) :: iReq

      INTEGER ierr

      CALL MPI_WAIT(iReq, MPI_STATUS_IGNORE, ierr)

      RETURN
      END SUBROUTINE WAITS
!--------------------------------------------------------------------
      SUBROUTINE WAITV(cm, iReq)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: iReq(:)

      INTEGER i, n

      n = SIZE(iReq)
      DO i=1, n
         CALL cm%wait(iReq(i))
      END DO

      RETURN
      END SUBROUTINE WAITV

!####################################################################
      
      FUNCTION REDUCEIS(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      INTEGER gU

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, 1, mpint, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCEIS
!--------------------------------------------------------------------
      FUNCTION REDUCEIV(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u(:)
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      INTEGER :: gU(SIZE(u))

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, SIZE(u), mpint, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCEIV
!--------------------------------------------------------------------
      FUNCTION REDUCERS(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      REAL(KIND=8) gU

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, 1, mpreal, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCERS
!--------------------------------------------------------------------
      FUNCTION REDUCERV(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      REAL(KIND=8) :: gU(SIZE(u))

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, SIZE(u), mpreal, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCERV
!--------------------------------------------------------------------
      FUNCTION REDUCERM(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:,:)
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      REAL(KIND=8) :: gU(SIZE(u,1),SIZE(u,2))

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, SIZE(u), mpreal, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCERM
!--------------------------------------------------------------------
      FUNCTION REDUCECS(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(IN) :: u
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      COMPLEX(KIND=8) gU

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, 1, mpcplx, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCECS
!--------------------------------------------------------------------
      FUNCTION REDUCECV(cm, u, iOp) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN), OPTIONAL :: iOp
      COMPLEX(KIND=8) :: gU(SIZE(u))

      INTEGER ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, SIZE(u), mpcplx, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCECV

!####################################################################
      
      FUNCTION GATHERIS(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u
      INTEGER, ALLOCATABLE :: gU(:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      IF (cm%seq()) THEN
         ALLOCATE(gU(1))
         gU = u
         RETURN
      END IF

      IF (cm%mas()) THEN
         ALLOCATE(gU(cm%np()))
      ELSE
         ALLOCATE(gU(0))
      END IF
      CALL MPI_GATHER(u, 1, mpint, gU, 1, mpint, master, cm%com(), ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(cm%np()))
         END IF
         CALL cm%bcast(gU)
      END IF

      RETURN
      END FUNCTION GATHERIS
!--------------------------------------------------------------------
      FUNCTION GATHERRS(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u
      REAL(KIND=8), ALLOCATABLE :: gU(:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      IF (cm%seq()) THEN
         ALLOCATE(gU(1))
         gU = u
         RETURN
      END IF

      IF (cm%mas()) THEN
         ALLOCATE(gU(cm%np()))
      ELSE
         ALLOCATE(gU(0))
      END IF
      CALL MPI_GATHER(u, 1, mpreal, gU, 1, mpreal, master, cm%com(), &
         ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(cm%np()))
         END IF
         CALL cm%bcast(gU)
      END IF

      RETURN
      END FUNCTION GATHERRS
!--------------------------------------------------------------------
      FUNCTION GATHERCS(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(IN) :: u
      COMPLEX(KIND=8), ALLOCATABLE :: gU(:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      IF (cm%seq()) THEN
         ALLOCATE(gU(1))
         gU = u
         RETURN
      END IF

      IF (cm%mas()) THEN
         ALLOCATE(gU(cm%np()))
      ELSE
         ALLOCATE(gU(0))
      END IF
      CALL MPI_GATHER(u, 1, mpcplx, gU, 1, mpcplx, master, cm%com(), &
         ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(cm%np()))
         END IF
         CALL cm%bcast(gU)
      END IF

      RETURN
      END FUNCTION GATHERCS
!--------------------------------------------------------------------
      FUNCTION GATHERIV(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u(:)
      INTEGER, ALLOCATABLE :: gU(:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      n = SIZE(u)
      IF (cm%seq()) THEN
         ALLOCATE(gU,SOURCE=u)
         RETURN
      END IF

      sCount = cm%GATHERIS(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(gU(gN), disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(gU(0), disp(0))
      END IF
      CALL MPI_GATHERV(u, n, mpint, gU, sCount, disp, mpint, master, &
         cm%com(), ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         CALL cm%bcast(gN)
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(gN))
         END IF
         CALL cm%bcast(gU)
      END IF
      DEALLOCATE(disp, sCount)

      RETURN
      END FUNCTION GATHERIV
!--------------------------------------------------------------------
      FUNCTION GATHERRV(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:)
      REAL(KIND=8), ALLOCATABLE :: gU(:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      n = SIZE(u)
      IF (cm%seq()) THEN
         ALLOCATE(gU,SOURCE=u)
         RETURN
      END IF

      sCount = cm%GATHERIS(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(gU(gN), disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(gU(0), disp(0))
      END IF
      CALL MPI_GATHERV(u, n, mpreal, gU, sCount, disp, mpreal, master,&
         cm%com(), ierr)
      

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         CALL cm%bcast(gN)
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(gN))
         END IF
         CALL cm%bcast(gU)
      END IF
      DEALLOCATE(disp, sCount)

      RETURN
      END FUNCTION GATHERRV
!--------------------------------------------------------------------
      FUNCTION GATHERCV(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      COMPLEX(KIND=8), INTENT(IN) :: u(:)
      COMPLEX(KIND=8), ALLOCATABLE :: gU(:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      n = SIZE(u)
      IF (cm%seq()) THEN
         ALLOCATE(gU,SOURCE=u)
         RETURN
      END IF

      sCount = cm%GATHERIS(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(gU(gN), disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(gU(0), disp(0))
      END IF
      CALL MPI_GATHERV(u, n, mpcplx, gU, sCount, disp, mpcplx, master,&
         cm%com(), ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         CALL cm%bcast(gN)
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(gN))
         END IF
         CALL cm%bcast(gU)
      END IF
      DEALLOCATE(disp, sCount)

      RETURN
      END FUNCTION GATHERCV
!--------------------------------------------------------------------
      FUNCTION GATHERIM(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u(:,:)
      INTEGER, ALLOCATABLE :: gU(:,:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr, n, m, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      m = SIZE(u,1)
      n = SIZE(u,2)
      IF (cm%seq()) THEN
         ALLOCATE(gU,SOURCE=u)
         RETURN
      END IF

      sCount = cm%GATHERIS(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(gU(m,gN), disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(gU(0,0), disp(0))
      END IF
      sCount = sCount*m
      disp   = disp*m
      CALL MPI_GATHERV(u, n*m, mpint, gU, sCount, disp, mpint, master, &
         cm%com(), ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         CALL cm%bcast(gN)
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(m,gN))
         END IF
         CALL cm%bcast(gU)
      END IF
      DEALLOCATE(disp, sCount)

      RETURN
      END FUNCTION GATHERIM
!--------------------------------------------------------------------
      FUNCTION GATHERRM(cm, u, toAll) RESULT(gU)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:,:)
      REAL(KIND=8), ALLOCATABLE :: gU(:,:)
      LOGICAL, INTENT(IN), OPTIONAL :: toAll

      LOGICAL fToAll
      INTEGER ierr, n, m, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (ALLOCATED(gU)) DEALLOCATE(gU)
      m = SIZE(u,1)
      n = SIZE(u,2)
      IF (cm%seq()) THEN
         ALLOCATE(gU,SOURCE=u)
         RETURN
      END IF

      sCount = cm%GATHERIS(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(gU(m,gN), disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(gU(0,0), disp(0))
      END IF
      sCount = sCount*m
      disp   = disp*m
      CALL MPI_GATHERV(u, n*m, mpreal, gU, sCount, disp, mpreal, master,&
         cm%com(), ierr)

      fToAll = .FALSE.
      IF (PRESENT(toAll)) fToAll = toAll
      IF (fToAll) THEN
         CALL cm%bcast(gN)
         IF (cm%slv()) THEN
            DEALLOCATE(gU)
            ALLOCATE(gU(m,gN))
         END IF
         CALL cm%bcast(gU)
      END IF
      DEALLOCATE(disp, sCount)

      RETURN
      END FUNCTION GATHERRM

!####################################################################
      
      SUBROUTINE SCATRIS(cm, gU, u)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: gU(:)
      INTEGER, INTENT(OUT) :: u

      INTEGER ierr

      IF (cm%seq()) THEN
         u = gU(1)
         RETURN
      END IF

      CALL MPI_SCATTER(gU, 1, mpint, u, 1, mpint, master, cm%com(),ierr)

      RETURN
      END SUBROUTINE SCATRIS
!--------------------------------------------------------------------
      SUBROUTINE SCATRRS(cm, gU, u)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: gU(:)
      REAL(KIND=8), INTENT(OUT) :: u

      INTEGER ierr

      IF (cm%seq()) THEN
         u = gU(1)
         RETURN
      END IF

      CALL MPI_SCATTER(gU, 1, mpreal, u, 1, mpreal, master, cm%com(),&
         ierr)

      RETURN
      END SUBROUTINE SCATRRS
!--------------------------------------------------------------------
      SUBROUTINE SCATRIV(cm, gU, u)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: gU(:)
      INTEGER, INTENT(INOUT) :: u(:)

      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (cm%seq()) THEN
         u = gU
         RETURN
      END IF

      n = SIZE(u)
      sCount = cm%gather(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(disp(0))
      END IF

      CALL MPI_SCATTERV(gU, sCount, disp, mpint, u, n, mpint, master, &
         cm%com(),ierr)

      RETURN
      END SUBROUTINE SCATRIV
!--------------------------------------------------------------------
      SUBROUTINE SCATRRV(cm, gU, u)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: gU(:)
      REAL(KIND=8), INTENT(INOUT) :: u(:)

      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (cm%seq()) THEN
         u = gU
         RETURN
      END IF

      n = SIZE(u)
      sCount = cm%gather(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(disp(0))
      END IF

      CALL MPI_SCATTERV(gU, sCount, disp, mpreal, u, n, mpreal, master, &
         cm%com(),ierr)

      RETURN
      END SUBROUTINE SCATRRV
!--------------------------------------------------------------------
      SUBROUTINE SCATRIM(cm, gU, u)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: gU(:,:)
      INTEGER, INTENT(INOUT) :: u(:,:)

      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (cm%seq()) THEN
         u = gU
         RETURN
      END IF

      n = SIZE(u)
      sCount = cm%gather(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(disp(0))
      END IF

      CALL MPI_SCATTERV(gU, sCount, disp, mpint, u, n, mpint, master, &
         cm%com(),ierr)

      RETURN
      END SUBROUTINE SCATRIM
!--------------------------------------------------------------------
      SUBROUTINE SCATRRM(cm, gU, u)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: gU(:,:)
      REAL(KIND=8), INTENT(INOUT) :: u(:,:)

      INTEGER ierr, n, gN, i
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)

      IF (cm%seq()) THEN
         u = gU
         RETURN
      END IF

      n = SIZE(u)
      sCount = cm%gather(n)
      
      IF (cm%mas()) THEN
         gN = SUM(sCount)
         ALLOCATE(disp(cm%np()))
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(disp(0))
      END IF

      CALL MPI_SCATTERV(gU, sCount, disp, mpreal, u, n, mpreal, master, &
         cm%com(),ierr)

      RETURN
      END SUBROUTINE SCATRRM

!####################################################################
!     Computes gU(map) = u
      SUBROUTINE GLOBALIS(cm, u, gU, map)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u(:), map(:)
      INTEGER, INTENT(OUT) :: gU(:)

      INTEGER, ALLOCATABLE :: tmp(:), gMap(:)

      IF (cm%seq()) THEN
         gU(map) = u
         RETURN
      END IF

      gMap = cm%gather(map)
      tmp  = cm%gather(u)
      IF (cm%slv()) RETURN

      gU(gMap) = tmp
      DEALLOCATE(gMap, tmp)

      RETURN
      END SUBROUTINE GLOBALIS
!---------------------------------------------------------------------
      SUBROUTINE GLOBALRS(cm, u, gU, map)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN) :: map(:)
      REAL(KIND=8), INTENT(OUT) :: gU(:)

      INTEGER, ALLOCATABLE :: gMap(:)
      REAL(KIND=8), ALLOCATABLE :: tmp(:)

      IF (cm%seq()) THEN
         gU(map) = u
         RETURN
      END IF

      gMap = cm%gather(map)
      tmp  = cm%gather(u)
      IF (cm%slv()) RETURN

      gU(gMap) = tmp
      DEALLOCATE(gMap, tmp)

      RETURN
      END SUBROUTINE GLOBALRS
!---------------------------------------------------------------------
      SUBROUTINE GLOBALRV(cm, u, gU, map)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:,:)
      INTEGER, INTENT(IN) :: map(:)
      REAL(KIND=8), INTENT(OUT) :: gU(:,:)

      INTEGER, ALLOCATABLE :: gMap(:)
      REAL(KIND=8), ALLOCATABLE :: tmp(:,:)

      IF (cm%seq()) THEN
         gU(:,map) = u
         RETURN
      END IF

      gMap = cm%gather(map)
      tmp  = cm%gather(u)
      IF (cm%slv()) RETURN

      gU(:,gMap) = tmp
      DEALLOCATE(gMap, tmp)

      RETURN
      END SUBROUTINE GLOBALRV

!####################################################################
!     Computes u = gU(map)
      SUBROUTINE LOCALIS(cm, u, gU, map)
      CLASS(cmType), INTENT(INout) :: cm
      INTEGER, INTENT(IN) :: gU(:), map(:)
      INTEGER, INTENT(OUT) :: u(:)

      INTEGER, ALLOCATABLE :: tmp(:), gMap(:)

      IF (cm%seq()) THEN
         u = gU(map)
         RETURN
      END IF

      gMap = cm%gather(map)
      IF (cm%mas()) THEN
         ALLOCATE(tmp(SIZE(gMap)))
         tmp = gU(gMap)
      ELSE 
         ALLOCATE(tmp(0))
      END IF
      CALL cm%scatter(tmp,u)

      DEALLOCATE(gMap, tmp)

      RETURN
      END SUBROUTINE LOCALIS
!---------------------------------------------------------------------
      SUBROUTINE LOCALRS(cm, u, gU, map)
      CLASS(cmType), INTENT(INout) :: cm
      INTEGER, INTENT(IN) :: map(:)
      REAL(KIND=8), INTENT(IN) :: gU(:)
      REAL(KIND=8), INTENT(OUT) :: u(:)

      INTEGER, ALLOCATABLE :: gMap(:)
      REAL(KIND=8), ALLOCATABLE :: tmp(:)

      IF (cm%seq()) THEN
         u = gU(map)
         RETURN
      END IF

      gMap = cm%gather(map)
      IF (cm%mas()) THEN
         ALLOCATE(tmp(SIZE(gMap)))
         tmp = gU(gMap)
      ELSE 
         ALLOCATE(tmp(0))
      END IF
      CALL cm%scatter(tmp,u)

      DEALLOCATE(gMap, tmp)

      RETURN
      END SUBROUTINE LOCALRS
!---------------------------------------------------------------------
      SUBROUTINE LOCALRV(cm, u, gU, map)
      CLASS(cmType), INTENT(INout) :: cm
      INTEGER, INTENT(IN) :: map(:)
      REAL(KIND=8), INTENT(IN) :: gU(:,:)
      REAL(KIND=8), INTENT(OUT) :: u(:,:)

      INTEGER, ALLOCATABLE :: gMap(:)
      REAL(KIND=8), ALLOCATABLE :: tmp(:,:)

      IF (cm%seq()) THEN
         u = gU(:,map)
         RETURN
      END IF

      gMap = cm%gather(map)
      IF (cm%mas()) THEN
         ALLOCATE(tmp(SIZE(gU,1),SIZE(gMap)))
         tmp = gU(:,gMap)
      ELSE 
         ALLOCATE(tmp(0,0))
      END IF
      CALL cm%scatter(tmp,u)

      DEALLOCATE(gMap, tmp)

      RETURN
      END SUBROUTINE LOCALRV

!####################################################################
!     Receives variable to write, u, file handle, fid, and the latest
!     position in the file, pos, and returns updated position, pos
!     The file must be opened with flags: ACCESS="STREAM"
      SUBROUTINE WRITEIS(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      CALL WRITEIV(cm, (/u/), fid, pos)

      RETURN
      END SUBROUTINE WRITEIS
!---------------------------------------------------------------------
      SUBROUTINE WRITERS(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      CALL WRITERV(cm, (/u/), fid, pos)

      RETURN
      END SUBROUTINE WRITERS
!---------------------------------------------------------------------
      SUBROUTINE WRITEIV(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: u(:)
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      INTEGER n, ipos, fpos
      INTEGER, ALLOCATABLE :: nG(:)

      fpos = 1
      IF (PRESENT(pos)) fpos = pos

      n = SIZE(u)
      IF (cm%seq()) THEN
         WRITE(fid,POS=fpos) mpint, 1, n, u
         fpos = fpos + mpint_size*(3 + n)
         IF (PRESENT(pos)) pos = fpos
         RETURN
      END IF

      nG = cm%gather(n,.TRUE.)
      IF (cm%mas()) THEN
         WRITE(fid,POS=fpos) mpint, cm%np(), nG, u
      ELSE
         ipos = fpos + mpint_size*(2 + cm%np() + SUM(nG(1:cm%id())))
         WRITE(fid,POS=ipos) u
      END IF
      fpos = fpos + mpint_size*(2 + cm%np() + SUM(nG))
      DEALLOCATE(nG)
      IF (PRESENT(pos)) pos = fpos

      RETURN
      END SUBROUTINE WRITEIV
!---------------------------------------------------------------------
      SUBROUTINE WRITERV(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(IN) :: u(:)
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      INTEGER n, ipos, fpos
      INTEGER, ALLOCATABLE :: nG(:)

      fpos = 1
      IF (PRESENT(pos)) fpos = pos
      
      n = SIZE(u)
      IF (cm%seq()) THEN
         WRITE(fid,POS=fpos) mpreal, 1, n, u
         fpos = fpos + mpint_size*3 + mpreal_size*n
         IF (PRESENT(pos)) pos = fpos
         RETURN
      END IF

      nG = cm%gather(n,.TRUE.)
      IF (cm%mas()) THEN
         WRITE(fid,POS=fpos) mpreal, cm%np(), nG, u
      ELSE
         ipos = fpos + mpint_size *(2 + cm%np())&
                     + mpreal_size*SUM(nG(1:cm%id()))
         WRITE(fid,POS=ipos) u
      END IF
      fpos = fpos + mpint_size*(2 + cm%np()) + mpreal_size*SUM(nG)
      DEALLOCATE(nG)
      IF (PRESENT(pos)) pos = fpos

      RETURN
      END SUBROUTINE WRITERV
!---------------------------------------------------------------------
      SUBROUTINE READIS(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(OUT) :: u
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      INTEGER iu(1)

      CALL READIV(cm, iu, fid, pos)
      u = iu(1)

      RETURN
      END SUBROUTINE READIS
!---------------------------------------------------------------------
      SUBROUTINE READRS(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(OUT) :: u
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      REAL(KIND=8) iu(1)

      CALL READRV(cm, iu, fid, pos)
      u = iu(1)

      RETURN
      END SUBROUTINE READRS
!---------------------------------------------------------------------
      SUBROUTINE READIV(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(INOUT) :: u(:)
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      INTEGER n, ipos, np, dType, fpos
      INTEGER, ALLOCATABLE :: nG(:)

      fpos = 1
      IF (PRESENT(pos)) fpos = pos
      
      ALLOCATE(nG(cm%np()))
      IF (cm%mas()) THEN
         READ(fid,POS=fpos) dType, np, nG
         IF (np .NE. cm%np()) ERROR STOP "READIV: incompatible np"
         IF (dType .NE. mpint) ERROR STOP "READIV: incompatible dType"
      END IF
      fpos = fpos + mpint_size*(2 + cm%np())

      CALL cm%bcast(nG)
      n = nG(cm%tf())
      ipos = fpos + mpint_size*SUM(nG(1:cm%id()))
      IF (n .NE. SIZE(u)) ERROR STOP "READIV: incompatible SIZE(u)"
      READ(fid,POS=ipos) u
      fpos = fpos + mpint_size*SUM(nG)
      DEALLOCATE(nG)
      IF (PRESENT(pos)) pos = fpos

      RETURN
      END SUBROUTINE READIV
!---------------------------------------------------------------------
      SUBROUTINE READRV(cm, u, fid, pos)
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=8), INTENT(OUT) :: u(:)
      INTEGER, INTENT(IN) :: fid
      INTEGER, INTENT(INOUT), OPTIONAL :: pos

      INTEGER n, ipos, np, dType, fpos
      INTEGER, ALLOCATABLE :: nG(:)

      fpos = 1
      IF (PRESENT(pos)) fpos = pos
      
      ALLOCATE(nG(cm%np()))
      IF (cm%mas()) THEN
         READ(fid,POS=fpos) dType, np, nG
         IF (np .NE. cm%np()) ERROR STOP "READRV: incompatible np"
         IF (dType .NE. mpreal) ERROR STOP "READRV: incompatible dType"
      END IF
      fpos = fpos + mpint_size*(2 + cm%np())

      CALL cm%bcast(nG)
      n = nG(cm%tf())
      ipos = fpos + mpreal_size*SUM(nG(1:cm%id()))
      IF (n .NE. SIZE(u)) ERROR STOP "READRV: incompatible SIZE(u)"
      READ(fid,POS=ipos) u
      fpos = fpos + mpreal_size*SUM(nG)
      DEALLOCATE(nG)
      IF (PRESENT(pos)) pos = fpos

      RETURN
      END SUBROUTINE READRV

!####################################################################
      
      SUBROUTINE barrierCm(cm)
      CLASS(cmType), INTENT(IN) :: cm

      INTEGER ierr

      CALL MPI_BARRIER(cm%com(),ierr)

      RETURN
      END SUBROUTINE barrierCm

!####################################################################
      
      FUNCTION createCm(cm, pid) RESULT(ch)
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: pid(:)
      INTEGER ch

      INTEGER n, cmGrp, newGrp, ierr

      n = SIZE(pid)
      CALL MPI_COMM_GROUP(cm%com(), cmGrp, ierr)
      CALL MPI_GROUP_INCL(cmGrp, n, pid, newGrp, ierr)
      CALL MPI_GROUP_FREE(cmGrp,ierr)
      CALL MPI_COMM_CREATE(cm%com(), newGrp, ch, ierr)
      CALL MPI_GROUP_FREE(newGrp,ierr)

      RETURN
      END FUNCTION createCm

!####################################################################
      
      SUBROUTINE finalizeCm(cm)
      CLASS(cmType), INTENT(INOUT) :: cm

      INTEGER ierr

      cm%cHndl  = MPI_COMM_NULL
      cm%taskId = 0
      cm%nProcs = 1
      CALL MPI_FINALIZE(ierr)

      RETURN
      END SUBROUTINE finalizeCm

!####################################################################
!     Testing current the implementation. Also this should give you a
!     hint on how to use this module.
      SUBROUTINE TEST_CMMOD()
   
      CHARACTER(LEN=*), PARAMETER :: fName=".dummy.bin"
      INTEGER i, iV(2), iM(2,2), fid, pos, nI, nIV(2)
      REAL(KIND=8) r, rV(2), rM(2,2), nR, nRV(2)
      COMPLEX(KIND=8) c, cV(2)
      CHARACTER(LEN=5) s, sV(2)
      INTEGER, ALLOCATABLE :: gI(:), gIM(:,:)
      REAL(KIND=8), ALLOCATABLE :: gR(:), gRM(:,:)
      COMPLEX(KIND=8), ALLOCATABLE :: gC(:)
      
      IF (cm%com() .NE. MPI_COMM_WORLD) ERROR STOP "Issue with cm%com"
      IF (cm%mas()) THEN
         PRINT *, "Current processor ID: ", cm%id()
         PRINT *, "Total number of processors: ", cm%np()
         PRINT *, "Number of threads: ", cm%nT()
      ELSE
         IF (.NOT.cm%slv()) ERROR STOP "Issue with cm%slv()"
      END IF
      IF (cm%id() .NE. cm%tF()-1) ERROR STOP "Issue with cm%tf()"
      IF (cm%np().NE.1 .AND. cm%seq()) ERROR STOP "Issue with cm%seq()"
      IF (cm%np() .GT. 2) ERROR STOP "Test case is for 1 or procs."
      
      IF (cm%mas()) THEN
         PRINT *, "Testing cm%bcast() ..."
         i  = 1
         iV = 1
         r  = 1D0
         rV = 1D0
         c  = (1D0,1D0)
         cV = c
         s  = "Hello"
         sV = s
      ELSE
         i  = 0
         iV = 0
         r  = 0D0
         rV = 0D0
         c  = 0D0
         cV = 0D0
         s  = ""
         sV = s
      END IF
      CALL cm%bcast(i)
      CALL cm%bcast(iV)
      CALL cm%bcast(r)
      CALL cm%bcast(rV)
      CALL cm%bcast(c)
      CALL cm%bcast(cV)
      CALL cm%bcast(s)
      CALL cm%bcast(sV)
      IF (i.NE.1 .OR. ANY(iV.NE.1) .OR. r.NE.1D0 .OR. ANY(rV.NE.1D0)&
         .OR. c.NE.(1D0,1D0) .OR. ANY(cV.NE.(1D0,1D0)) .OR. s.NE."Hello"&
         .OR. ANY(sV.NE."Hello")) ERROR STOP "Issue with cm%bcast"
      IF (cm%mas()) PRINT *, "(PASSED)"

      IF (.NOT.cm%seq()) THEN
         IF (cm%mas()) THEN
            PRINT *, "Testing cm%send() and cm%recv() ..."
            CALL cm%recv(rV,1)
            IF (ANY(rV .NE. 2D0)) ERROR STOP "Issue with cm%send/recv"
            PRINT *, "(PASSED)"
         ELSE
            rV = 2D0
            CALL cm%send(rV,0)
         END IF
         IF (cm%mas()) THEN
            PRINT *, "Testing cm%isend() and cm%irecv() ..."
            i = cm%irecv(rV,1)
            CALL cm%wait(i)
            IF (ANY(rV .NE. 3D0)) ERROR STOP "Issue with cm%send/recv"
            PRINT *, "(PASSED)"
         ELSE
            rV = 3D0
            i = cm%isend(rV,0)
            CALL cm%wait(i)
         END IF
      END IF

      IF (cm%mas()) PRINT *, "Testing cm%reduce() ..."
      i  = 1
      rV = 1D0
      i  = cm%reduce(i)
      iV = cm%reduce(iV)
      r  = cm%reduce(r)
      rV = cm%reduce(rV)
      c  = cm%reduce(c)
      cV = cm%reduce(cV)
      IF (cm%mas()) THEN
         IF (i.NE.cm%np() .OR. ANY(iV.NE.cm%np()) .OR.&
            r.NE.REAL(cm%np(),8) .OR. ANY(rV.NE.REAL(cm%np(),8)) .OR.&
            c.NE.(1D0,1D0)*cm%np() .OR. ANY(cV.NE.(1D0,1D0)*cm%np()))&
            ERROR STOP "Issue with cm%reduce"
         PRINT *, "(PASSED)"
      END IF

      iM = 1
      rM = 1D0
      IF (cm%mas()) PRINT *, "Testing cm%gather() ..."
      gI = cm%gather(i)
      gR = cm%gather(r)
      gC = cm%gather(c)
      IF (cm%mas()) THEN
         IF (ANY(gI.NE.cm%np()) .OR. ANY(gR.NE.1D0*cm%np()) .OR.&
            ANY(gC.NE.(1D0,1D0)*cm%np())) &
            ERROR STOP "Issue with cm%gather"
      END IF
      gI  = cm%gather(iV)
      gR  = cm%gather(rV)
      gC  = cm%gather(cV)
      gIM = cm%gather(iM)
      gRM = cm%gather(rM)
      IF (cm%mas()) THEN
         IF (ANY(gI.NE.cm%np()) .OR. ANY(gR.NE.1D0*cm%np()) .OR.&
            ANY(gC.NE.(1D0,1D0)*cm%np()) .OR. ANY(gIM.NE.1) .OR. &
            ANY(gRM.NE.1D0)) ERROR STOP "Issue with cm%gather"
         PRINT *, "(PASSED)"
      END IF

      i = 0
      r = 0D0
      IF (cm%mas()) THEN
         iV = (/1,2/)
         rV = (/1D0,2D0/)
      END IF
      IF (cm%mas()) PRINT *, "Testing cm%scatter() ..."
      CALL cm%scatter(iV, i)
      CALL cm%scatter(rV, r)
      IF (cm%mas()) THEN
         IF (i.NE.cm%tf() .OR. r.NE.REAL(i,8)) ERROR STOP &
            "Issue with cm%scatter"
      END IF
      CALL cm%scatter(gI,iV)
      CALL cm%scatter(gR,rV)
      IF (cm%mas()) THEN
         IF (ANY(iV.NE.cm%np()) .OR. ANY(rV.NE.1D0*cm%np()))&
            ERROR STOP "Issue with cm%scatter"
         PRINT *, "(PASSED)"
      END IF

      fid = 74
      pos = 1
      OPEN(fid,FILE=fName,ACCESS='STREAM')
      IF (cm%mas()) PRINT *, "Testing cm%write/read ..."
      CALL cm%write(iV, fid, pos)
      CALL cm%write(rV, fid, pos)
      CALL cm%write(i, fid, pos)
      CALL cm%write(r, fid, pos)
      CLOSE(fid)
      pos = 1
      OPEN(fid,FILE=fName,ACCESS='STREAM')
      CALL cm%read(nIV, fid, pos)
      CALL cm%read(nRV, fid, pos)
      CALL cm%read(nI,  fid, pos)
      CALL cm%read(nR,  fid, pos)
      IF (ANY(iV.NE.nIV) .OR. ANY(rV.NE.nRV) .OR. i.NE.nI .OR. r.NE.nR) &
         ERROR STOP "Issue with cm%write/read"
      IF (cm%mas()) PRINT *, "(PASSED)"
      CLOSE(fid)

      IF (cm%mas()) PRINT *, "Testing cm%barrier ..."
      CALL cm%barrier()
      IF (cm%mas()) THEN
         OPEN(fid,FILE=fName)
         CLOSE(fid,STATUS='DELETE')
         PRINT *, "(PASSED)"
      END IF

      IF (cm%mas()) THEN
         PRINT *, "Testing cm%global ..."
         nIV = (/3,1/)
      ELSE
         nIV = (/2,4/)
      END IF
      iV = nIV + 1
      rV = iV
      rM(1,:) = rV
      rM(2,:) = rV
      DEALLOCATE(gI, gR, gRM)
      ALLOCATE(gI(4), gR(4), gRM(2,4))
      gR = gI
      gRM(1,:) = gR
      gRM(2,:) = gR
      CALL cm%global(iV,gI,nIV)
      CALL cm%global(rV,gR,nIV)
      CALL cm%global(rM,gRM,nIV)
      IF ((cm%np().EQ.2 .AND. ANY(gI.NE.(/2,3,4,5/)) .AND. cm%mas())&
         .OR. ((gI(1).NE.2.OR.gI(3).NE.4).AND.cm%np().EQ.1) .OR. &
         ANY(NINT(gR).NE.gI) .OR. ANY(NINT(gRM(2,:)).NE.gI)) &
         ERROR STOP "Issue with cm%global"
      IF (cm%mas()) PRINT *, "(PASSED)"

      IF (cm%mas()) PRINT *, "Testing cm%local ..."
      iV = 0
      CALL cm%local(iV,gI,nIV)
      IF (ANY(NINT(rV).NE.iV)) ERROR STOP "Issue with cm%local"
      rV = 0
      rM = 0
      CALL cm%local(rV,gR,nIV)
      CALL cm%local(rM,gRM,nIV)
      IF (ANY(NINT(rV).NE.iV) .OR. ANY(NINT(gRM(2,:)).NE.gI)) &
         ERROR STOP "Issue with cm%local"
      IF (cm%mas()) PRINT *, "(PASSED)"

      IF (.NOT.cm%seq()) THEN
         IF (cm%mas()) PRINT *, "Testing cm%create() ..."
         iV = cm%id()
         i  = cm%create(iV)
         IF (cm%mas()) PRINT *, "(PASSED)"
      END IF

      RETURN
      END SUBROUTINE TEST_CMMOD
      END MODULE CMMOD
