C
C @(#)shbstore.f
C     modified from store.f 3.3 (BNP) 3/16/89
C
      SUBROUTINE STORE(N,ISW,J,S)
      implicit none
      INCLUDE 'simphb.h'
      INTEGER N,ISW,J
      REAL*8 S(N)
C
      INTEGER MAXLL
      PARAMETER (MAXLL = 2)
C
      INTEGER STORQ,RETRQ,STORP,RETRP
      PARAMETER (STORQ = 1,RETRQ = 2,STORP = 3,RETRP = 4)
C
      IF (ISW.EQ.STORQ) THEN
         CALL DCOPY(N,S,1,QQ((J+MAXLL-1)*N+1),1)
      ELSE IF (ISW.EQ.RETRQ) THEN
         CALL DCOPY(N,QQ((J+MAXLL-1)*N+1),1,S,1)
      ELSE IF (ISW.EQ.STORP) THEN
         IF (J.GT.MAXLL) STOP 'STORE: (STORP) J.GT.MAXLL'
         CALL DCOPY(N,S,1,QQ((J-1)*N+1),1)
      ELSE IF (ISW.EQ.RETRP) THEN
         IF (J.GT.MAXLL) STOP 'STORE: (RETRP) J.GT.MAXLL'
         CALL DCOPY(N,QQ((J-1)*N+1),1,S,1)
      ENDIF
      RETURN
      END
