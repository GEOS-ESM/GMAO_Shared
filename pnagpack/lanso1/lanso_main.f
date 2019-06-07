C
C @(#)main.f	3.15 (BNP) 6/3/89; from main.f 2.9 6/23/88
C
      PROGRAM MAIN
C
C.... The common block "OPSTIF" is used solely by the user-supplied
C.... subroutines OP, OPM and holds the matrices K & M.
C.... In our implementation, given below, M is diagonal and its
C.... elements are stored in the array D.  K is stored in
C.... active-column form in S; JDIAG goes with S.
C
      INTEGER JDIAG
      DOUBLE PRECISION S,D
      COMMON/OPSTIF/JDIAG(1000),S(40000),D(1000)
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
      INTEGER I,EV,NW,N,LANMAX,MAXPRS,K,J,NEIG,IERR,MSGLVL
      DOUBLE PRECISION CONDM,OLDL,OLDR,SHIFT,ENDL,ENDR,KAPPA,
     *   W(LMTNW),RITZ(LMTN),BND(LMTN)
      INTEGER NCELLS
      NCELLS(N,LANMAX,I) = 6*N+4*LANMAX+1+I*LANMAX*LANMAX
C
C.... NW is the size of W
C
      NW = LMTNW
C
C.... SHIFT is a point in the interval [oldl,oldr], i.e. instead of
C.... [K,M], [M,K-SHIFT*M] problem is solved
C
      READ(1,*)N,LANMAX,MAXPRS,CONDM,OLDL,OLDR,SHIFT,EV,KAPPA,MSGLVL
      WRITE(11,2000)SHIFT,N,LANMAX,MAXPRS,CONDM,OLDL,OLDR,EV.LE.0,KAPPA
C
C.... Data validation
C
      IF (OLDL.GE.OLDR) THEN
         WRITE(11,*)'OLDL must be strictly less than OLDR'
         STOP
      ELSE IF (SHIFT.LE.OLDL.OR.SHIFT.GE.OLDR) THEN
         WRITE(11,*)'SHIFT must be inside [OLDL,OLDR].'
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
C.... Users should replace the following lines with their own code.
C
C.... START of REPLACEMENT CODE.
C
C.... Input matrices M & K in sequence.
C
      READ(5,4000)(JDIAG(I),D(I),I=1,N)
      READ(5,5000)(S(I),I=1,JDIAG(N))
4000  FORMAT((I10,D20.12))
5000  FORMAT((4D20.12))
C
C.... Make the requested shift
C
      DO 10 I = 1,N
         K = JDIAG(I)
         S(K) = S(K)-SHIFT*D(I)
10    CONTINUE
C
C.... ACTCOL factorizes K-SHIFT*M
C
      CALL ACTCOL(S,W,JDIAG,N,.TRUE.,.FALSE.)
C
C.... END of REPLACEMENT CODE
C
      ENDL = 1.0D0/(OLDL-SHIFT)
      ENDR = 1.0D0/(OLDR-SHIFT)
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
C.... User may insert code to recover eigenvalues of [K,M]
C.... using true_eigenvalue(I) := SHIFT+1/RITZ(I).
C.... If you permute eigenvalues, don't forget to permute the
C.... eigenvectors with them.
C
      WRITE(11,9999)J,NEIG,
     *   (I,1.0D0/RITZ(I)+SHIFT,RITZ(I),BND(I),I = 1,J)
      STOP
9999  FORMAT(1X,'...... '
     *   /1X,'...... ',6X,' J =',I3,3X,' NEIG =',I3
     *   /1X,'...... '
     *   /1X,'...... ',3X,3X,' Recovered Ritz values',2X,
     *   '(','  Computed Ritz values',2X,'Error Bnds',')'
     *   /1X,'...... '
     *   /(1X,'...... ',I3,3X,1PD22.14,2X,
     *   '(',1PD22.14,2X,1PD10.2,')'))
2000  FORMAT(
     *    1X,'... '
     *   /1X,'... solve the eigenproblem [ K - shift*M, M ]'
     *   /1X,'... shift                     =',1PD22.14
     *   /1X,'... no. of equations          =',I4
     *   /1X,'... max. no. of Lanczos steps =',I4
     *   /1X,'... max. no. of eigenpairs    =',I4
     *   /1X,'... condition number of M     =',1PD22.14
     *   /1X,'... left  end of the interval =',1PD22.14
     *   /1X,'... right end of the interval =',1PD22.14
     *   /1X,'... eigenvalues only? [T/F]   =',L4
     *   /1X,'... kappa                     =',1PD22.14
     *   /1X,'... ')
2999  FORMAT(1X,'... return flag =',I9,4(' '),'...')
      END
