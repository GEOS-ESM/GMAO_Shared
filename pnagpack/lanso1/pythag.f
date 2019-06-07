      DOUBLE PRECISION FUNCTION PYTHAG(A,B) 
      DOUBLE PRECISION A,B 
C 
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW 
C 
      DOUBLE PRECISION P,R,S,T,U 
      P = MAX(ABS(A),ABS(B)) 
      IF (P .EQ. 0.0D0) GO TO 20 
      R = (MIN(ABS(A),ABS(B))/P)**2 
   10 CONTINUE 
         T = 4.0D0 + R 
         IF (T .EQ. 4.0D0) GO TO 20 
         S = R/T 
         U = 1.0D0 + 2.0D0*S 
         P = U*P 
         R = (S/U)**2 * R 
      GO TO 10 
   20 PYTHAG = P 
      RETURN 
      END 
