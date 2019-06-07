C
C @(#)std.f	1.9 (BNP) 6/3/89
C
      PROGRAM STD
C
C.... The example program uses LANDR to solve the standard eigenvalue
C.... problem [A-lambda*I]x = 0.
C....
C.... For simplicity A here is the famous Wilkinson matrix W(N)+.
C.... (pp.308 Wilkinson AEP)  In this case, /OPSTIF/ SIMPLY stores
C.... M := N/2.
C
      INTEGER M
      COMMON/OPSTIF/M
C....
C.... OP(N,P,Q,R) takes an N-vector Q and should return A*Q in Y.
C.... OPM(N,X,Y)  takes an N-vector X and should return X in Y.
C
C.... The common block "GETPUT" is used solely by the user-supplied
C.... subroutine STORE and holds the Lanczos vectors.
C
      DOUBLE PRECISION A
      COMMON/GETPUT/A(125000)
C
C.... See write-up for details on memory allocation for W and RITZ.
C.... N <= LMTN, MAXPRS <= LMTPRS, LANMAX <= LMTLAN.
C
      INTEGER LMTN,LMTPRS,LMTLAN,LMTNW
      PARAMETER (LMTN = 500,LMTPRS = 250,LMTLAN = LMTPRS)
      PARAMETER (LMTNW = 6*LMTN+4*LMTLAN+1+LMTLAN*LMTLAN)
      INTEGER I,EV,NW,N,LANMAX,MAXPRS,J,NEIG,IERR,MSGLVL
      DOUBLE PRECISION CONDM,ENDL,ENDR,KAPPA,
     *   W(LMTNW),RITZ(LMTN),BND(LMTN)
      INTEGER NCELLS
      NCELLS(N,LANMAX,I) = 6*N+4*LANMAX+1+I*LANMAX*LANMAX
C
C.... NW is the size of the working array W
C
      NW = LMTNW
C
      READ(5,*)N,LANMAX,MAXPRS,ENDL,ENDR,EV,KAPPA,MSGLVL
C
C.... make N an odd number
C
      M = N/2
      N = M+M+1
      WRITE(11,2000)N,LANMAX,MAXPRS,ENDL,ENDR,EV.LE.0,KAPPA
      CONDM = 1.0D0
C
C.... Data validation
C
      IF (ENDL.GE.ENDR) THEN
         WRITE(11,*)'ENDL must be strictly less than ENDR'
         STOP
      ELSE IF (MAXPRS.GT.LANMAX) THEN
         WRITE(11,*)'MAXPRS cannot exceed LANMAX',
     *      '(MAXPRS =',MAXPRS,' LANMAX =',LANMAX,')'
         STOP
      ELSE IF (EV.LE.0.AND.NCELLS(N,LANMAX,0).GT.NW) THEN
         WRITE(11,*)'6*N+4*LANMAX+1 cannot exceed',
     *      NW,'(N =',N,' LANMAX =',LANMAX,')'
         STOP
      ELSE IF (EV.GT.0.AND.NCELLS(N,LANMAX,1).GT.NW) THEN
          WRITE(11,*)'6*N+4*LANMAX+1+LANMAX*LANMAX cannot exceed',
     *      NW,'(N =',N,' LANMAX =',LANMAX,')'
         STOP
      ENDIF
C
C.... To get a random starting vector, the first N cells
C.... must be initialized to zero.
C
      DO 17 I = 1,N
         W(I) = 0.0D0
17    CONTINUE
C
C.... Make a Lanczos run; see LANDR for meanings of parameters.
C
      CALL LANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
C     CALL LANDR(...,J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
     *   J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
      IF (IERR.NE.0) WRITE(11,2999)IERR
C
C.... If you permute eigenvalues, don't forget to permute the
C.... eigenvectors with them.
C
      WRITE(11,9999)J,NEIG,(I,RITZ(I),BND(I),I = 1,J)
      STOP
9999  FORMAT(1X,'...... '
     *   /1X,'...... ',6X,' J =',I3,3X,' NEIG =',I3
     *   /1X,'...... '
     *   /1X,'...... ',3X,3X,'  Computed Ritz values',2X,
     *    '(','Error Bnds',')'
     *   /1X,'...... '
     *   /(1X,'...... ',I3,3X,1PD22.14,2X,'(',1PD10.2,')'))
2000  FORMAT(
     *    1X,'... '
     *   /1X,'... solve the standard eigenproblem'
     *   /1X,'... no. of equations          =',I4
     *   /1X,'... max. no. of Lanczos steps =',I4
     *   /1X,'... max. no. of eigenpairs    =',I4
     *   /1X,'... left  end of the interval =',1PD22.14
     *   /1X,'... right end of the interval =',1PD22.14
     *   /1X,'... eigenvalues only? [T/F]   =',L4
     *   /1X,'... kappa                     =',1PD22.14
     *   /1X,'... ')
2999  FORMAT(1X,'... return flag =',I9,4(' '),'...')
      END
C
      SUBROUTINE OP(N,P,Q,R)
      INTEGER N
      DOUBLE PRECISION P(N),Q(N),R(N)
C
      INTEGER M
      COMMON/OPSTIF/M
C
      INTEGER I
      DOUBLE PRECISION ALPHA
      ALPHA(I) = DBLE(IABS((M+1-I)))
C
      R(1) = ALPHA(1)*Q(1)+Q(1+1)
      DO 100 I = 2,N-1
         R(I) = Q(I-1)+ALPHA(I)*Q(I)+Q(I+1)
100   CONTINUE
      R(N) = Q(N-1)+ALPHA(N)*Q(N)
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
