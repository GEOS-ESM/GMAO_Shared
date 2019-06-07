      REAL*8 FUNCTION DDOTMPI(N,DX,INCX,DY,INCY,MPICOM)
      INCLUDE 'mpif.h'
C
C     forms the dot product of two conformally distributed vectors
C
      REAL*8 DX(*), DY(*), STEMP, DDOT
      INTEGER IERR,N,INCX,INCY,MPICOM
      EXTERNAL DDOT
C
C     call DDOT on local data to compute partial sum
C
      STEMP = DDOT(N,DX,INCX,DY,INCY)
      IERR = 0
C
C     call MPI communication routine to add the partial sums together
C
      if (mpi_real8.eq.0) then
         CALL MPI_ALLREDUCE(STEMP,DDOTMPI,1,MPI_DOUBLE_PRECISION,
     $        MPI_SUM,MPICOM,IERR)
      else
         CALL MPI_ALLREDUCE(STEMP,DDOTMPI,1,MPI_REAL8,MPI_SUM,
     $        MPICOM,IERR)
      endif
      IF (IERR.NE.MPI_SUCCESS) CALL MPI_ABORT(MPI_COMM_WORLD,IERR)
C
      RETURN
      END
      SUBROUTINE DATX(N,DA,DX,INCX,DY,INCY)
C
C     dy := da*dx
C
      REAL*8 DX(*),DY(*),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0D0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C     
C     unequal increments or equal increments .ne. one
C     
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I = 1,N
         DY(IY) = DA*DX(IX)
         IX = IX+INCX
         IY = IY+INCY
 10   CONTINUE
      RETURN
C     
C     code for both increments equal to 1
C     
 20   M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
         DY(I) = DA*DX(I)
 30   CONTINUE
      IF (N.LT.4) RETURN
 40   MP1 = M+1
      DO 50 I = MP1,N,4
         DY(I) = DA*DX(I)
         DY(I+1) = DA*DX(I+1)
         DY(I+2) = DA*DX(I+2)
         DY(I+3) = DA*DX(I+3)
 50   CONTINUE
      RETURN
      END
C     
