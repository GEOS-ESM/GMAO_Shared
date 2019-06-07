C
C @(#)ppurge.f	7/31/97 from 3.16 (BNP) 5/11/89
C
      SUBROUTINE pPURGE(N,astore,LL,J,R,Q,RA,QA,WRK,ETA,OLDETA,MSGLVL,MPICOM)
      INTEGER N,LL,J,MSGLVL,MPICOM
      REAL*8 R(N),Q(N),RA(N),QA(N),WRK(N),ETA(J),OLDETA(J)
      REAL*8 Astore(*)
C
C.... This routine examines ETA to decide whether
C.... re-orthogonalization should be performed.
C
C.... N      (local) dimension of the eigenproblem
C.... LL     no. of initial Lanczos vectors in local orthog.
C.... J      current Lanczos step
C.... R      the residual vector to become the next Lanczos vector
C.... Q      the current Lanczos vector
C.... RA     the product of the mass matrix and r
C.... QA     the product of the mass matrix and q
C.... WRK    a temporary vector to hold the previous Lanczos vectors
C.... ETA    state of orthogonality between r and previous Lanczos vectors
C.... OLDETA state of orthogonality between q and previous Lanczos vectors
C.... MPICOM MPI Communicator
C
C.... BLAS routines:    DAXPY,DDOTMPI,IDAMAX
C.... subroutines:      none
C.... user-supplied:    OPM,STORE
C
      INCLUDE 'mpif.h'
      INTEGER RETRQ
      PARAMETER (RETRQ = 2)
C
      REAL*8 RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      INTEGER K,I,LOOP,IDAMAX,MY_ID
      REAL*8 T,TQ,TR,REPS1,DDOTMPI
C
      IF (J.LE.LL+1) RETURN
      MY_ID = -1
      K = IDAMAX(J-(LL+1),ETA(LL+1),1)+LL
C
      IF (ABS(ETA(K)).GT.REPS) THEN
         IF (MSGLVL.GT.10) THEN
            CALL MPI_COMM_RANK(MPICOM, MY_ID, I)
            IF (MY_ID.EQ.0)
     &           PRINT *, 'PPURGE: Omega(', k, j+1, ') = ', abs(eta(k))
         ENDIF
         REPS1 = EPS1/REPS
         DO 55 LOOP = 1,2
            IF (RNM.GT.TOL) THEN
               IF (MY_ID.EQ.0 .AND. MSGLVL.GT.10) THEN
                  PRINT *, 'PPURGE: purging Q_', J, ' and R_', J,
     &                 ' LOOP #', LOOP
               ENDIF
C
C....       Bring in a Lanczos vector t and orthogonalize both r and q
C....       against it
C
               TQ = 0.0D0
               TR = 0.0D0
               DO 50 I = 1,J-1
                  CALL STORE(N,astore,RETRQ,I,WRK)
                  T = -DDOTMPI(N,QA,1,WRK,1,MPICOM)
                  TQ = TQ+ABS(T)
                  CALL DAXPY(N,T,WRK,1,Q,1)
                  T = -DDOTMPI(N,RA,1,WRK,1,MPICOM)
                  TR = TR+ABS(T)
                  CALL DAXPY(N,T,WRK,1,R,1)
50             CONTINUE
               CALL OPM(N,Q,QA,MPICOM)
C
C....          restore local orthogonality
C
               T = -DDOTMPI(N,R,1,QA,1,MPICOM)
               TR = TR+ABS(T)
               CALL DAXPY(N,T,Q,1,R,1)
C
               CALL OPM(N,R,RA,MPICOM)
               RNM = SQRT(DDOTMPI(N,RA,1,R,1,MPICOM))
               IF (TQ.LE.REPS1.AND.TR.LE.REPS1*RNM) GOTO 58
            ENDIF
55       CONTINUE
58       DO 60 I = LL+1,J
            ETA(I) = EPS1
            OLDETA(I) = EPS1
60       CONTINUE
      ENDIF
      RETURN
      END
