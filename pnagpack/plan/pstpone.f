C
C @(#)pstpone.f	7/31/97 from 3.13 (BNP) 6/3/89
C
      SUBROUTINE pSTPONE(N,astore,R,WRK,ALF,NQ,MSGLVL,MPICOM)
      INTEGER N,NQ(4),MSGLVL,MPICOM
      REAL*8 R(5*N),WRK(N),ALF(1)
      REAL*8 astore(*)
C
C.... This routine performs the first step of the Lanczos algorithm.
C.... It performs a step of extended local re-orthogonalization.
C
C.... N      (local) dimension of the eigenproblem
C.... R      an array containing [r(j),q(j),q(j-1),p(j),p(j-1)/Mr(j)]
C.... ALF    diagonal elements of T
C.... NQ(4)  location pointers for the array R
C.... MPICOM MPI communicator
C
C.... BLAS routines:    DATX,DAXPY,DDOTMPI,DSCAL
C.... subroutines:      STARTV
C.... user-supplied:    OP,OPM,STORE
C
      REAL*8 RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      REAL*8 T,pSTARTV,DDOTMPI
C
      REAL*8 ONE,ZERO
      DATA ONE,ZERO/1.0D0,0.0D0/
C
C.... get initial vector, default is random
C
      RNM = pSTARTV(N,astore,1,R,WRK,NQ,EPS,MSGLVL,MPICOM)
      IF (RNM.EQ.ZERO) RETURN
C
C.... normalize starting vector
C
      T = ONE/RNM
      CALL DATX(N,T,R,1,R(NQ(1)),1)
      CALL DSCAL(N,T,R(NQ(3)),1)
C
C.... take the first step
C
      CALL OP(N,R(NQ(3)),R(NQ(1)),R,MPICOM)
      ALF(1) = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
      CALL DAXPY(N,-ALF(1),R(NQ(1)),1,R,1)
C
C.... restore local orthogonality
C
      T = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
      CALL DAXPY(N,-T,R(NQ(1)),1,R,1)
      ALF(1) = ALF(1)+T
      CALL OPM(N,R,R(NQ(4)),MPICOM)
      RNM = SQRT(DDOTMPI(N,R,1,R(NQ(4)),1,MPICOM))
      ANORM = RNM+ABS(ALF(1))
      TOL = EPSN*ANORM
C
      RETURN
      END
