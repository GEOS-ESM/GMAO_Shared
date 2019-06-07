      SUBROUTINE LOADEV(EV,N,J,NRITZ,RITZ,BND,RITZV,IERR)
      INTEGER EV,N,J,NRITZ,IERR
      DOUBLE PRECISION RITZ(J),BND(J),RITZV(N,J)
C
C.... Subroutine to load the accepted Ritz values, their associated
C.... error bounds and Ritz vectors back into memory from FORTRAN I/O
C.... channel EV.
C
C.... Input parameters
C
C.... EV    FORTRAN channel containing the accepted Ritz values/vectors
C.... N     dimension of the eigenproblem
C.... J     number of Lanczos steps actually taken
C
C.... Output parameters:
C
C.... NRITZ number of accepted Ritz values/vectors
C.... RITZ  accepted Ritz values
C.... BND   associated (refined) error bounds
C.... RITZV accepted Ritz vectors, NRITZ of them
C.... IERR  error condition as reported by READ
C
      INTEGER I,K
C
      NRITZ = 0
      OPEN(EV,FORM='UNFORMATTED')
      REWIND(EV)
      DO 500 K = 1,J
         READ(EV,IOSTAT = IERR,ERR = 700,END = 600)
     *      RITZ(K),BND(K),(RITZV(I,K),I = 1,N)
500   CONTINUE
600   IERR = 0
      NRITZ = K-1
700   CLOSE(EV)
      RETURN
      END
