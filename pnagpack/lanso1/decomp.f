C*********************************************************************
C*********************************************************************
C***  DECOMPOSITION OF AN ARBTIRARY GRAPH VIA THE EIGENVECTOR ********
C*********************************************************************
C*********************************************************************
C
C     WRITTEN BY HORST D. SIMON AT NASA AMES RESEARCH CENTER
C     PARTS OF THE CODE BASED ON THELANCZOS ALGORITHM IMPLEM-
C     ENTATION OF PARLETT, NOUR-OMID AND LIU
C
C*********************************************************************
C
      PROGRAM  DECOMP
C
C     --------------------------------------------------------------
C     ... THE COMMON BLOCK ADJMAT IS USED FOR STORING THE ADJANCENCY
C         STRUCTURE OF THE MATRIX IN BOEING-HARWELL FORMAT
C     --------------------------------------------------------------
C
      INTEGER XADJ, ADJ
      COMMON  /ADJMAT/ XADJ(4000), ADJ(40000)
C
C     -------------------------------------------------------------
C     ... THE COMMON BLOCK GETPUT IS USED FOR STORING THE LANCZOS
C         VECTORS
C     -------------------------------------------------------------
C
      DOUBLE PRECISION A
      COMMON/GETPUT/A(125000)
C
C     -------------------------------------------------------------
C     ... PARAMETERS
C         LANMAX -- MAXIMUM NUMBER OF LANCZOS STEPS
C         MAXPRS -- MAXIMUM NUMBER OF RITZPAIRS TO BE COMPUTED
C         MSGLVL -- MESSAGE LEVEL
C     -------------------------------------------------------------
C
C.... See write-up for details on memory allocation for W, RITZ, IW and Y.
C.... In our subroutine ENOUGH, 4*LANMAX <= N is silently enforced.
C
      INTEGER LMTN, LMTPRS, LMTLAN, LMTNW, LANMAX, MAXPRS, MSGLVL
C
      PARAMETER (LMTN = 500,LMTPRS = 250,LMTLAN = LMTPRS)
      PARAMETER (LMTNW = 6*LMTN+4*LMTLAN+1+LMTLAN*LMTLAN)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      CHARACTER*80  TITLE
      CHARACTER*3   CODE
      CHARACTER*16  PNTFMT, RINFMT
      CHARACTER*20  NVLFMT
C
      INTEGER NLINE, NPLS, NRILS, NNVLS, NRHSLS,
     1        NROW, NCOL, NIND, NELE
C
      INTEGER I, EV, NW, N, NEIG, IERR, J
C
      DOUBLE PRECISION CONDM,ENDL,ENDR,KAPPA,
     1                 W(6*LMTN+4*LMTLAN+1+LMTLAN*LMTLAN),
     1                 RITZ(LMTN), BND(LMTN)
C
C     -----------------------------------------
C     ... SET NW TO THE LENGTH OF THE WORKSPACE
C     -----------------------------------------
C
      NW = LMTNW
C
C     ------------------------------------------
C     ... READ MATRICES IN BOEING-HARWELL FORMAT
C     ------------------------------------------
C
      READ (5, 9000) TITLE
      READ (5, 9001) NLINE, NPLS, NRILS, NNVLS, NRHSLS
      READ (5, 9002) CODE, NROW, NCOL, NIND, NELE
      READ (5, 9003) PNTFMT, RINFMT, NVLFMT
C
      WRITE (11, 9000) TITLE
      WRITE (11, 9001) NLINE, NPLS, NRILS, NNVLS, NRHSLS
      WRITE (11, 9002) CODE, NROW, NCOL, NIND, NELE
      WRITE (11, 9003) PNTFMT, RINFMT, NVLFMT
C
      READ (5, PNTFMT) (XADJ(I) , I = 1, NROW+1)
      READ (5, RINFMT) (ADJ(I)  , I = 1, NIND)
C
C     -----------------------------------------------
C     ... WRITE STRUCTURAL INFORMATION TO OUTPUT FILE
C     -----------------------------------------------
C
      WRITE (9, 9000) TITLE
      WRITE (9, 9001) NLINE, NPLS, NRILS, NNVLS, NRHSLS
      WRITE (9, 9002) CODE, NROW, NCOL, NIND, NELE
      WRITE (9, 9003) PNTFMT, RINFMT, NVLFMT
C
      WRITE (9, PNTFMT) (XADJ(I) , I = 1, NROW+1)
      WRITE (9, RINFMT) (ADJ(I)  , I = 1, NIND)
C
C     ----------------------------------------
C     ... INITIALIZE VARIABLES FOR LANCZOS RUN
C     ----------------------------------------
C
      N   = NROW
C     ENDL = 0.0D0 
C     ENDR = 0.0D0
C
C     DO 200 I = 1, N
C        ENDL = MAX (ENDL, DFLOAT(2*(XADJ(I+1)-XADJ(I))-1))
C 200 CONTINUE
C     ENDL = - ENDL
C
C.... EV.LE.0 means "eigenvalues only"
C
      CONDM = 1.0D0
      ENDR = 1.0D-30
      ENDL = -ENDR
      EV = 0
      KAPPA = 0.0D0
      LANMAX = 112
      MAXPRS = 112
      MSGLVL = 3
      WRITE (11,2000) N, LANMAX, MAXPRS, ENDL, ENDR, EV.LE.0
C
C     --------------------------------
C     ... INITIALIZE WORKSPACE TO ZERO
C     --------------------------------
C
      DO 300 I = 1,N
         W(I) = 0.0D0
300   CONTINUE
C
C     ----------------------
C     ... MAKE A LANCZOS RUN
C     ----------------------
C
      CALL LANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
     1   J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
C
      IF (IERR.NE.0) WRITE(11,2999) IERR
C
C     ---------------------
C     ... PRINT EIGENVALUES
C     ---------------------
C
      WRITE(11,9999)J,NEIG,(I,RITZ(I),BND(I),I = 1,J)
      STOP
C    
C     -----------
C     ... FORMATS
C     -----------
C
9999  FORMAT(
     *   /1X,'...... '
     *   /1X,'...... ',6X,' J =',I3,3X,' NEIG =',I3
     *   /1X,'...... '
     *   /1X,'...... ',3X,3X,
     *   '  Computed Ritz values',2X,'Error Bnds'
     *   /1X,'...... '
     *   /(1X,'...... ',I3,3X,1PD22.14,2X,1PD10.2))
C
2000  FORMAT(
     *    1X,'... '
     *   /1X,'... no. of equations          =',I4
     *   /1X,'... max. no. of Lanczos steps =',I4
     *   /1X,'... max. no. of eigenpairs    =',I4
     *   /1X,'... left  end of the interval =',1PD12.4
     *   /1X,'... right end of the interval =',1PD12.4
     *   /1X,'... eigenvalues only? [T/F]   =',L4
     *   /1X,'... ')
2999  FORMAT(1X,'... return flag =',I9,4(' '),'...')
C
9000  FORMAT ( A80 )
C
9001  FORMAT ( 5I14 )
C
9002  FORMAT ( A3, 11X, 4I14 )
C
9003  FORMAT ( 2A16, 2A20 )
C
      END
