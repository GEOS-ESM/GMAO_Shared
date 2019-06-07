C
C @(#)blas.f	3.4 (BNP) 12/9/88; from blas.f 2.2 10/13/87
C
C       (standard double precision BLAS)
C
        SUBROUTINE DATX(N,DA,DX,INCX,DY,INCY)
C
C       dy := da*dx
C
        DOUBLE PRECISION DX(*),DY(*),DA
        INTEGER I,INCX,INCY,IX,IY,M,MP1,N
        IF (N.LE.0) RETURN
        IF (DA.EQ.0.0D0) RETURN
        IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C       unequal increments or equal increments .ne. one
C
        IX = 1
        IY = 1
        IF (INCX.LT.0) IX = (-N+1)*INCX+1
        IF (INCY.LT.0) IY = (-N+1)*INCY+1
        DO 10 I = 1,N
          DY(IY) = DA*DX(IX)
          IX = IX+INCX
          IY = IY+INCY
10      CONTINUE
        RETURN
C
C       code for both increments equal to 1
C
20      M = MOD(N,4)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
          DY(I) = DA*DX(I)
30      CONTINUE
        IF (N.LT.4) RETURN
40      MP1 = M+1
        DO 50 I = MP1,N,4
          DY(I) = DA*DX(I)
          DY(I+1) = DA*DX(I+1)
          DY(I+2) = DA*DX(I+2)
          DY(I+3) = DA*DX(I+3)
50      CONTINUE
        RETURN
        END
