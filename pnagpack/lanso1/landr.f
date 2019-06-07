C
C @(#)landr.f	3.16 (BNP) 6/3/89; from landr.f 2.16 6/7/88
C
      SUBROUTINE LANDR(N,LANMAX,astore,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
C     SUBROUTINE LANDR(...,J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
     *   J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
      INTEGER N,LANMAX,MAXPRS,EV,J,NEIG,NW,IERR,MSGLVL
      DOUBLE PRECISION CONDM,ENDL,ENDR,KAPPA,
     *   RITZ(LANMAX),BND(LANMAX),W(NW)
C
C.... The program make a Lanczos run using a linear operator that acts
C.... through a user supplied subroutine called OP in order to compute
C.... either all eigenvalues outside an interval [endl,endr] or the
C.... first MAXPRS eigenvalues, whichever occurs first.  The inner
C.... product is defined implicitly by a user supplied subroutine OPM.
C
C            *****************************************
C            *                                       *
C            *         LANCZOS ALORITHM WITH         *
C            *      SELECTIVE ORTHOGONALIZATION      *
C            *            L  A  N  S  O              *
C            *      (Using Simon's Recurrence)       *
C            *                                       *
C            *****************************************
C
C.... inputs
C.... N      dimension of the eigenproblem
C.... LANMAX upper limit to the number of Lanczos steps
C.... MAXPRS upper limit to the number of wanted eigenpairs
C.... CONDM  estimated effective condition number of M
C.... ENDL   left end of the interval containing the UNwanted eigenvalues
C.... ENDR   right end of the interval containing the UNwanted eigenvalues
C.... EV     .LE.0 means eigenvalues only,
C....        .GT.0 means both eigenvalues and eigenvectors are wanted
C....              and EV becomes the output channel for eigenvectors
C.... KAPPA  relative accuracy of Ritz values acceptable as eigenvalues
C....        (KAPPA will not be touched and/or looked at if EV.LE.0)
C.... NW     length of the work array W
C.... W      work array of length NW
C
C.... outputs
C.... J      number of Lanczos steps actually taken
C.... NEIG   number of Ritz values stabilized
C.... RITZ   array to hold the Ritz values
C.... BND    array to hold the error bounds
C.... IERR   error flag
C
C.... subroutines:      MACHAR,LANSO,RITVEC
C
      DOUBLE PRECISION RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      INTEGER MT,I,NQ(4),N1,L2,L3,L4,L5,LSQ6
C
C.... MACHAR specific (EPS declared in COMMON/RDATA/)
C
      INTEGER IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,MAXEXP
      DOUBLE PRECISION EPSNEG,XMIN,XMAX
      DOUBLE PRECISION astore(*)           
C
C.... check input data
C
      MT = 6*N+1+4*LANMAX
      IF (EV.GT.0) MT = MT+LANMAX*LANMAX
      IERR = 0
      IF (N.LE.0) IERR = IERR+1
      IF (LANMAX.LE.0) IERR = IERR+2
      IF (ENDR.LE.ENDL) IERR = IERR+4
      IF (MAXPRS.LE.0) IERR = IERR+8
      IF (MAXPRS.GT.LANMAX) IERR = IERR+16
      IF (LANMAX.GT.N) IERR = IERR+32
      IF (MT.GT.NW) IERR = IERR+64
      IF (IERR.GT.0) RETURN
C
C.... compute machine precision
C
      CALL MACHAR(IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,MAXEXP,
     *   EPS,EPSNEG,XMIN,XMAX)
      EPS1 = EPS*SQRT(DBLE(N))
      EPSN = EPS*DBLE(N)
      RCEPS1 = SQRT(MAX(ABS(CONDM),1.0D0))*EPS1
      REPS = SQRT(EPS)
      EPS34 = REPS*SQRT(REPS)
C
C.... set pointers and initialize
C
      NQ(1) = N+1
      DO 20 I = 2,4
         NQ(I) = NQ(I-1)+N
20    CONTINUE
      N1 = 1+5*N
      L2 = N+N1
      L3 = LANMAX+L2
      L4 = 1+LANMAX+L3
      L5 = LANMAX+L4
      CALL LANSO(N,LANMAX,MAXPRS,ENDL,ENDR,J,NEIG,RITZ,BND,
C     CALL LANSO(...,W,W(N1),W(L2),W(L3),W(L4),W(L5),NQ,IERR,MSGLVL)
     *   W,W(N1),W(L2),W(L3),W(L4),W(L5),NQ,IERR,MSGLVL)
C
C.... compute eigenvectors
C
      IF (EV.GT.0) THEN
         LSQ6 = LANMAX+L5
C
C....    rationalize KAPPA
C
         KAPPA = MAX(ABS(KAPPA),EPS34)
C
C....    W(NQ(4)) becomes working storage from this point on
C
         CALL RITVEC(N,astore,J,EV,KAPPA,RITZ,BND,W(L2),W(L3),
C        CALL RITVEC(...,W(LSQ6),W(NQ(4)),W(N1),IERR,MSGLVL)
     *      W(LSQ6),W(NQ(4)),W(N1),IERR,MSGLVL)
      ENDIF
      RETURN
      END
