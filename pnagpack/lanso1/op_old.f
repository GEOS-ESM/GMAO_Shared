C
C @(#)op.f	3.3 (BNP) 6/3/89; from op.f 2.3 10/13/87
C
      SUBROUTINE OP(N,P,Q,R)
      INTEGER N
      DOUBLE PRECISION P(N),Q(N),R(N)
C
      INTEGER JDIAG
      DOUBLE PRECISION S,D
      COMMON/OPSTIF/JDIAG(1000),S(40000),D(1000)
C
      INTEGER I
C
      DO 100 I = 1,N
         R(I) = P(I)
100   CONTINUE
      CALL ACTCOL(S,R,JDIAG,N,.FALSE.,.TRUE.)
      RETURN
      END
C
      SUBROUTINE OPM(N,A,B)
      INTEGER N
      DOUBLE PRECISION A(N),B(N)
C
      INTEGER JDIAG
      DOUBLE PRECISION S,D
      COMMON/OPSTIF/JDIAG(1000),S(40000),D(1000)
C
      INTEGER I
C
      DO 100 I = 1,N
         B(I) = A(I)*D(I)
100   CONTINUE
      RETURN
      END
C
      SUBROUTINE ACTCOL(A,B,JDIAG,NEQ,AFAC,BACK)
      INTEGER NEQ,JDIAG(NEQ)
      DOUBLE PRECISION A(*),B(*)
      LOGICAL AFAC,BACK
C
      INTEGER NEIGLO,I,JR,J,JD,JH,IS,IE,K,ID,IR,IH
      DOUBLE PRECISION DG,D,DDOT
C
C.... active column profile symmetric equation solver
C.... factor a to ut*d*u, reduce b
C
      NEIGLO = 0
      JR = 0
      DO 600 J = 1,NEQ
      JD = JDIAG(J)
      JH = JD-JR
      IS = J-JH+2
      DG = A(JD)
      IF (JH-2) 550,300,100
100   IF (.NOT.AFAC) GO TO 500
      IE = J-1
      K = JR+2
      ID = JDIAG(IS-1)
C
C.... reduce all equations except diagonal
C
      DO 200 I = IS,IE
      IR = ID
      ID = JDIAG(I)
      IH = MIN0(ID-IR-1,I-IS+1)
      IF (IH.GT.0) A(K) = A(K)-DDOT(IH,A(K-IH),1,A(ID-IH),1)
200   K = K+1
C
C.... reduce diagonal term
C
300   IF (.NOT.AFAC) GO TO 500
      IR = JR+1
      IE = JD-1
      K = J-JD
      DO 400 I = IR,IE
      ID = JDIAG(K+I)
      D = A(I)
      A(I) = A(I)*A(ID)
      DG = DG-D*A(I)
400   CONTINUE
      IF (DG*A(JD).LT.0.0D0) NEIGLO = NEIGLO+1
C
C.... reduce rhs
C
500   IF (BACK) B(J) = B(J)-DDOT(JH-1,A(JR+1),1,B(IS-1),1)
550   IF (DG.NE.0.0D0.AND.AFAC) A(JD) = 1.0D0/DG
600   JR = JD
      IF (.NOT.BACK) RETURN
C
C.... divide by diagonal pivots
C
      DO 700 I = 1,NEQ
      ID = JDIAG(I)
      DG = B(I)
      B(I) = B(I)*A(ID)
C
C.... backsubstitute
C
700   J = NEQ
      JD = JDIAG(J)
800   D = B(J)
      J = J-1
      IF (J.LE.0) GO TO 1100
      JR = JDIAG(J)
      IF (JD-JR.LE.1) GO TO 1000
      IS = J-JD+JR+2
      K = JR-IS+1
      DO 900 I = IS,J
900   B(I) = B(I)-A(I+K)*D
1000  JD = JR
      GO TO 800
1100  RETURN
      END
