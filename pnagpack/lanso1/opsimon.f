      SUBROUTINE OP(N,P,Q,R)
C********************************************************************
C**** MULTIPLICATION BY THE ADJACENCY MATRIX IN BOEING-HARWELL FORMAT
C********************************************************************
C
      INTEGER N
      DOUBLE PRECISION P(N),Q(N),R(N)
      INTEGER XADJ, ADJ
      COMMON /ADJMAT/ XADJ(4000), ADJ(40000)
C
      INTEGER I, J
C
      DO 50 I = 1, N
         R(I) = 0.0D0
   50 CONTINUE
C
      DO 100 I = 1, N
         R(I) =  R(I) + Q(I)
         DO 200 J = XADJ(I) + 1, XADJ(I+1) - 1
            R(I) = R(I) + Q(ADJ(J))
            R(ADJ(J)) = R(ADJ(J)) + Q(I)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE OPM(N,A,B)
      INTEGER N
      DOUBLE PRECISION A(N),B(N)
      INTEGER I
C
      DO 100 I = 1,N
         B(I) = A(I)
100   CONTINUE
      RETURN
      END
