C
C @(#)diag.f	1.4 (BNP) 6/3/89
C
      PROGRAM DIAG
C
C.... The example program (suggested by nagaab@vax.oxford.ac.uk) uses LANDR
C.... to solve a standard eigenvalue problem [K-lambda*I]x = 0.
C.... ITYPE identifies K.
C....
C.... ITYPE := 1 <=> K(I,I) = 1 for I = 1..N
C....          2 <=> K(I,I) = I for I = 1..N
C....          3 <=> K(I,I) = I*(-1)**(I-1) for I = 1..N
C....          4 <=> K(I,I) = (-1)**(I-1) for I = 1..N
C
      INTEGER ITYPE
      COMMON/OPSTIF/ITYPE
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
      INTEGER I,EV,NW,N,LANMAX,MAXPRS,J,NEIG,IERR,MSGLVL,NMAX,LOOP
      DOUBLE PRECISION CONDM,ENDL,ENDR,KAPPA,T,ERR,EXACT,
     *   W(LMTNW),RITZ(LMTN),BND(LMTN)
      INTEGER NCELLS
      NCELLS(N,LANMAX,I) = 6*N+4*LANMAX+1+I*LANMAX*LANMAX
C
      READ(5,*)ITYPE,N,EV,MSGLVL
C
C.... NW is the size of the working array W
C
      NW = LMTNW
      LANMAX = N
      MAXPRS = N
C
      IF (ITYPE.LT.0.OR.ITYPE.GT.4) THEN
         WRITE(11,*)'ITYPE must be one of 0,1,2,3,4'
         STOP
      ELSE IF (N.LE.0) THEN
         WRITE(11,*)'N must be > 0'
         STOP
      ELSE IF (EV.LE.0.AND.NCELLS(N,LANMAX,0).GT.NW) THEN
         WRITE(11,*)'6*N+4*LANMAX+1 cannot exceed',
     *      NW,'(N =',N,' LANMAX =',LANMAX,')'
         STOP
      ELSE IF (EV.GT.0.AND.NCELLS(N,LANMAX,1).GT.NW) THEN
          WRITE(11,*)'6*N+4*LANMAX+1+LANMAX*LANMAX cannot exceed',
     *      NW,'(N =',N,' LANMAX =',LANMAX,')'
         STOP
      ELSE IF (ITYPE.GT.0) THEN
         ENDL = -1.0D-30
         ENDR = 1.0D-30
         KAPPA = 0.0D0
         CONDM = 1.0D0
         DO 17 I = 1,N
            W(I) = 0.0D0
17       CONTINUE
         WRITE(11,2000)ITYPE,N,LANMAX,MAXPRS,ENDL,ENDR,EV.LE.0,KAPPA
         CALL LANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
     *      J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
         IF (IERR.NE.0) WRITE(11,2999)IERR
         WRITE(11,9999)J,NEIG,(I,RITZ(I),BND(I),I = 1,J)
      ELSE
         NMAX = N
         EV = 0
         MSGLVL = 0
         DO 200 ITYPE = 1,4
            N = 1
            DO 100 LOOP = 1,NMAX
               LANMAX = N
               MAXPRS = N
               ENDL = -1.0D-30
               ENDR = 1.0D-30
               KAPPA = 0.0D0
               CONDM = 1.0D0
               DO 30 I = 1,N
                  W(I) = 0.0D0
30             CONTINUE
               CALL LANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
     *            J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
               ERR = 0.0D0
               DO 50 I = 1,NEIG
                  T = EXACT(ITYPE,N,I)
                  ERR = MAX(ERR,ABS((RITZ(I)-T)/T))
50             CONTINUE
               WRITE(11,3000)ITYPE,N,J,NEIG,ERR
3000           FORMAT(1X,'ITYPE = ',I1,' N = ',I3,' J = ',I3,
     *            ' NEIG = ',I3,' ERR = ',1PE10.2)
               IF (N.GE.NMAX) THEN
                  GOTO 200
               ELSE IF (N.LT.10) THEN
                  N = N+1
               ELSE IF (N.LT.50) THEN
                  N = N+5
               ELSE IF (N.LT.100) THEN
                  N = N+10
               ELSE IF (N.LT.NMAX) THEN
                  N = N+20
               ENDIF
100         CONTINUE
200      CONTINUE
      ENDIF
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
     *   /1X,'... itype                     =',I4
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
      INTEGER ITYPE
      COMMON/OPSTIF/ITYPE
C
      INTEGER I,IS
C
      GOTO(10,20,30,40),ITYPE
10    DO 100 I = 1,N
         R(I) = Q(I)
100   CONTINUE
      RETURN
20    DO 200 I = 1,N
         R(I) = DBLE(I)*Q(I)
200   CONTINUE
      RETURN
30    IS = 1
      DO 300 I = 1,N
         R(I) = DBLE(SIGN(I,IS))*Q(I)
         IS = -IS
300   CONTINUE
      RETURN
40    IS = 1
      DO 400 I = 1,N
         R(I) = DBLE(IS)*Q(I)
         IS = -IS
400   CONTINUE
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
C
      DOUBLE PRECISION FUNCTION EXACT(ITYPE,N,I)
      INTEGER ITYPE,N,I
C
      INTEGER L
C
      GOTO(500,600,700,800),ITYPE
500   EXACT = DBLE(1)
      RETURN
600   EXACT = DBLE(I)
      RETURN
700   L = N/2
      L = L+L
      IF (I.LE.L/2) THEN
         EXACT = DBLE(I+I-L-2)
      ELSE
         EXACT = DBLE(I+I-L-1)
      ENDIF
      RETURN
800   IF (I.LE.N/2) THEN
         EXACT = DBLE(-1)
      ELSE
         EXACT = DBLE(1)
      ENDIF
      RETURN
      END
