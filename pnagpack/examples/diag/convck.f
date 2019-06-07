C     @(#)convck.f KJW 4/29/97
      SUBROUTINE CONVCK(J,RITZ,BND,EPS,NEIGL,NEIGR)
      INTEGER I,J,NEIGL,NEIGR
      REAL*8 RITZ(J), BND(J), EPS, TOL
C
C.... This subroutine counts the number of Ritz values with
C.... small error bound from both left and right
C
      NEIGL = 0
      NEIGR = 0
      CALL DSORT2(J,RITZ,BND)
      TOL = 1.6D1*EPS*MAX(ABS(RITZ(1)), ABS(RITZ(J)))
C
C.... count up
C
      I = 1
 10   IF (BND(I).LE.TOL) THEN
         I = I + 1
         IF (I+I.LE.J) GOTO 10
      ENDIF
      NEIGL = I-1
      IF (I.GT.J) RETURN
C
C.... count down
C
      I = J
 20   IF (BND(I).LE.TOL) THEN
         I = I - 1
         IF (I+I.GT.J) GOTO 20
      ENDIF
      NEIGR = J - I
      RETURN
      END
