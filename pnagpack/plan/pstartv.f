C
C @(#)pstartv.f	7/31/97 from 1.10 (BNP) 10/30/90
C
      REAL*8 FUNCTION pSTARTV(N,astore,J,R,WRK,NQ,EPS,MSGLVL,MPICOM)
      include 'mpif.h'
      INTEGER N,J,NQ(4),MSGLVL,MPICOM
      REAL*8 R(5*N),WRK(N),EPS
      REAL*8 astore(*)
C
C.... This routine delivers a starting vector in R and returns |R|;
C.... it returns ZERO if range is spanned or if no starting vector
C.... within range of operator can be found.
C
C.... N      (local) dimension of the eigenproblem
C.... J      starting index for a Lanczos run
C.... R      an array containing [r(j),q(j),q(j-1),p(j),p(j-1)/Mr(j)]
C.... NQ(4)  location pointers for the array R
C.... MPICOM MPI communicator
C
C.... BLAS routines:    DAXPY,DDOTMPI
C.... subroutines:      RANDOM
C.... user-supplied:    OP,OPM,STORE
C
      INTEGER RETRQ
      PARAMETER (RETRQ = 2)
C
      INTEGER IRAND,I,ID,LOOP, doscale
      REAL*8 RNM2,T,RANDOM,DDOTMPI
      SAVE IRAND
C
      REAL*8 ZERO
      DATA ZERO/0.0D0/, IRAND/0/
C
C.... Initialize the sed to be used to generate random numbers
C.... To limit the number of times MPI_COMM_RANK is called, IRAND is
C.... saved across calls to this routine
C
      IF (IRAND.EQ.0) THEN
         CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRAND, I)
         IF (I .NE. MPI_SUCCESS) STOP 'MPI_COMM_RANK FAILED'
         IRAND = 918272+(J+IRAND)*7219
      ENDIF
C
C.... get initial vector, default is random
C
      RNM2 = DDOTMPI(N,R,1,R,1,MPICOM)
      DO 60 ID = 1,3
         IF (ID.GT.1.OR.J.GT.1.OR. (.not. RNM2.GT.ZERO)) THEN
            DO 20 I = 1,N
               R(I) = RANDOM(IRAND)-0.5D0
20          CONTINUE
         ENDIF
C.... put input vector in the range of M, do not multiply with K
C.... better for restarting and allow zero (0) eigenvalue to be computed
         CALL OPM(N,R,R(NQ(3)),MPICOM)
         RNM2 = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
         IF (.NOT. RNM2.GT.ZERO) THEN
            GOTO 60
         ELSE
            if (rnm2.lt.1D-8 .or. rnm2.gt.1D8) then
               doscal = 1
               rnm2 = 1.0D0 / sqrt(rnm2)
               call dscal(n,rnm2,r(NQ(3)),1)
            else
               doscal = 0
            endif
         ENDIF
         CALL OPM(N,R(NQ(3)),R,MPICOM)
         if (doscal .eq. 1) then
            RNM2 = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
            rnm2 = 1.0D0 / sqrt(rnm2)
            call dscal(n,rnm2,r(NQ(3)),1)
         endif
         CALL OPM(N,R,R(NQ(3)),MPICOM)
         RNM2 = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
         IF (J.GT.1) THEN
            DO 50 LOOP = 1,2
               DO 40 I = 1,J-1
                  CALL STORE(N,astore,RETRQ,I,WRK)
                  T = -DDOTMPI(N,R(NQ(3)),1,WRK,1,MPICOM)
                  CALL DAXPY(N,T,WRK,1,R,1)
40             CONTINUE
               CALL OPM(N,R,R(NQ(3)),MPICOM)
               T = DDOTMPI(N,R(NQ(3)),1,R,1,MPICOM)
               IF (T.LE.EPS*RNM2) THEN
                  T = ZERO
                  GOTO 55
               ENDIF
50          CONTINUE
55          RNM2 = T
         ENDIF
         IF (RNM2.GT.ZERO) GOTO 80
60    CONTINUE
80    pSTARTV = SQRT(RNM2)
      RETURN
      END
