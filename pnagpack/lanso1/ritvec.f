C
C @(#)ritvec.f	1.7 (BNP) 5/18/89
C
      SUBROUTINE RITVEC(N,astore,J,EV,KAPPA,Lvec,RITZ,BND,ALF,BET,S,WRK1,WRK2,
     *   IERR,MSGLVL)

      INTEGER N,LANMAX,J,EV,IERR,MSGLVL
      DOUBLE PRECISION KAPPA,Lvec(N,J),RITZ(J),BND(J),ALF(J),BET(J),
     *   S(J,J),WRK1(N),WRK2(N)
C
C.... subroutines:      TQL2,STORE
C.... BLAS routines:    DAXPY,DCOPY,DSCAL
C
      INTEGER RETRQ
      PARAMETER (RETRQ = 2)
C
      INTEGER I,K,NVEC,ICNT

      DOUBLE PRECISION astore(*)
C
      CALL DSCAL(J*J,0.0D0,S,1)
      DO 10 I = 1,J
         S(I,I) = 1.0D0
10    CONTINUE
      CALL DCOPY(J,ALF,1,WRK1,-1)
      IF (J.GT.1) CALL DCOPY(J-1,BET(2),1,WRK2(2),-1)
      CALL TQL2(J,J,WRK1,WRK2,S,IERR)
      IF (IERR.NE.0) RETURN
C
C.... on return WRK1 contains eigenvalues in ascending order
C....       and S contains the corresponding eigenvectors
C
C_JKIM      OPEN(EV,FORM='UNFORMATTED')
C_JKIM      REWIND(EV)
C_JKIM      WRITE(EV)N,J,KAPPA
      nvec= 0
      do 40 k=1,j
         if (bnd(k).le.kappa*abs(ritz(k))) nvec= nvec+1
   40 continue

      icnt=0


      DO 50 K = 1,J
         IF (BND(K).LE.KAPPA*ABS(RITZ(K))) THEN
            CALL DSCAL(N,0.0D0,WRK1,1)
            DO 20 I = 1,J
               CALL STORE(N,astore,RETRQ,I,WRK2)
               CALL DAXPY(N,S(J-I+1,K),WRK2,1,WRK1,1)
20          CONTINUE
C_JKIM            WRITE(EV)RITZ(K),BND(K),(WRK1(I),I=1,N)
          icnt= icnt+1
          do 30 i=1,n
                Lvec(i,icnt)= wrk1(i)
   30 continue

         ENDIF
50    CONTINUE
C_JKIM      CLOSE(EV)
      RETURN
      END
