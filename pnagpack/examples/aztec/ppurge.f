C
C @(#)shbpurge.f	0.01 (KW) 6/3/97; from purge.f 3.16 5/11/89
C
      SUBROUTINE PPURGE(N,LL,J,R,Q,RA,QA,WRK,ETA,OLDETA,MSGLVL,MPICOM)
      implicit none
      INTEGER N,LL,J,MSGLVL,MPICOM
      REAL*8 R(N),Q(N),RA(N),QA(N),WRK(N),ETA(J),OLDETA(J)
C
C.... This routine examines ETA to decide whether
C.... re-orthogonalization should be performed.
C.... NOTE: this routine uses up to two classic Gram-schmidt process
C.... to restore orthogonality.  If WRK is not large enough, it uses
C.... ETA and OLDETA as workspace.
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
C.... BLAS routines:    IDAMAX
C.... subroutines:      ppurge2 (in C)
C.... user-supplied:    OPM
C
      REAL*8 RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      INTEGER I, K, MY_ID, IDAMAX
      REAL*8 REPS1
C
      IF (J.LE.LL+1) RETURN
      K = IDAMAX(J-(LL+1),ETA(LL+1),1)+LL
C
      IF (ABS(ETA(K)).GT.REPS) THEN
         IF (MSGLVL.GT.10) THEN
            CALL MPI_COMM_RANK(MPICOM, MY_ID, I)
            IF (MY_ID.EQ.0)
     &           PRINT *, 'PPURGE: Omega(', k, j+1, ') = ', abs(eta(k))
         ENDIF
c
         REPS1 = EPS1/REPS
         CALL PPURGE2(N,J,R,Q,RA,QA,WRK,ETA,OLDETA,REPS1,RNM,
     &        MSGLVL,MPICOM)
         DO 60 I = 1,J
            ETA(I) = EPS1
            OLDETA(I) = EPS1
60       CONTINUE
      ENDIF
      RETURN
      END
