C
C @(#)diag.f	1.4 (BNP) 6/3/89
C
      PROGRAM DIAG
      include 'mpif.h'
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
      INTEGER ITYPE, OFFSET
      COMMON/OPSTIF/ITYPE, OFFSET
C
C.... The common block "GETPUT" is used solely by the user-supplied
C.... subroutine STORE and holds the Lanczos vectors.
C
      REAL*8 A, Astore(126000)
      COMMON/GETPUT/A(126000)
C
C.... See write-up for details on memory allocation for W and RITZ.
C.... N <= LMTN, MAXPRS <= LMTPRS, LANMAX <= LMTLAN.
C
      INTEGER LMTN,LMTPRS,LMTLAN,LMTNW
      PARAMETER (LMTN = 500,LMTPRS = 250,LMTLAN = LMTPRS)
      PARAMETER (LMTNW = 6*LMTN+5*LMTLAN+1+LMTLAN*LMTLAN)
      INTEGER I,EV,NW,N,LANMAX,MAXPRS,J,NEIG,IERR,MSGLVL,NMAX
      INTEGER MPICOM, MPIRANK, NTOTAL, NPROC, NEIGL, NEIGR
      REAL*8 CONDM,ENDL,ENDR,KAPPA,T,ERR,EXACT,
     *   W(LMTNW),RITZ(LMTLAN),BND(LMTLAN)
      INTEGER NCELLS
      NCELLS(N,LANMAX,I) = 5*N+4*LANMAX+1+I*LANMAX*LANMAX+
     *     MAX(N,LANMAX+1)
C
C.... NW is the size of the working array W
C
      NW = LMTNW
      NMAX = LMTN
      N = NMAX
      EV = 0
      MSGLVL = 15
C
      IF (EV.LE.0.AND.NCELLS(NMAX,NMAX,0).GT.NW) THEN
         WRITE(*,*)'6*N+4*LANMAX+1 cannot exceed',
     *      NW,'(N =',NMAX,' LANMAX =',NMAX,')'
         STOP
      ELSE IF (EV.GT.0.AND.NCELLS(NMAX,NMAX,1).GT.NW) THEN
          WRITE(*,*)'6*N+4*LANMAX+1+LANMAX*LANMAX cannot exceed',
     *      NW,'(N =',NMAX,' LANMAX =',NMAX,')'
         STOP
      ENDIF
C
C     START MPI
C
      CALL MPI_INIT(IERR)
      IF (IERR .NE. MPI_SUCCESS) STOP 'FAILED TO INIT MPI'
      MPICOM = MPI_COMM_WORLD
      CALL MPI_COMM_SIZE(MPICOM, NPROC, IERR)
      CALL MPI_COMM_RANK(MPICOM, MPIRANK, IERR)
C
      DO 200 ITYPE = 1,4
         N = 1
 100     IF (N.LE.NMAX) THEN
            NTOTAL = N * NPROC
            OFFSET = N * MPIRANK
            LANMAX = MIN(NTOTAL, LMTLAN)
            MAXPRS = MAX(1, LANMAX / 2)
            ENDL = -1.0D-30
            ENDR = 1.0D-30
            KAPPA = 0.0D0
            CONDM = 1.0D0
C           initial guess to LANCZOS routine (zero ==> a random vector)
            DO 30 I = 1,N
               W(I) = 1.0D0
 30         CONTINUE
            IF (MSGLVL.GT.2) THEN
               WRITE(10+MPIRANK,2000)ITYPE,NTOTAL,LANMAX,MAXPRS,ENDL,
     &              ENDR,EV.LE.0,KAPPA
            ENDIF
!            CALL pLANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
!     &           J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL,MPICOM)
            CALL pLANDR(N,LANMAX,astore,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
     &           J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL,MPICOM)
            IF (IERR.NE.0) WRITE(10+MPIRANK,2999)IERR
            IF (MSGLVL.GE.0) THEN
               IF (MSGLVL.GT.3)
     &              WRITE(10+MPIRANK,9999)J,NEIG,
     &              (I,RITZ(I),BND(I),I = 1,J)
               ERR = 1.1D-16
               CALL CONVCK(J,RITZ,BND,ERR,NEIGL,NEIGR)
               ERR = 0.0D0
               DO 50 I = 1,NEIGL
                  T = EXACT(ITYPE,NTOTAL,I)
                  ERR = MAX(ERR,ABS((RITZ(I)-T)/T))
 50            CONTINUE
               DO 60 I = J,J-NEIGR+1,-1
                  T = EXACT(ITYPE,NTOTAL,NTOTAL-J+I)
                  ERR = MAX(ERR,ABS((RITZ(I)-T)/T))
 60            CONTINUE
               IF (NEIGL+NEIGR .EQ. 0) THEN
                  T = EXACT(ITYPE,NTOTAL,1)
                  ERR = ABS((RITZ(1) - T)/T)
                  T = EXACT(ITYPE,NTOTAL,NTOTAL)
                  ERR = MIN(ERR, ABS((RITZ(J)-T)/T))
               ENDIF
               WRITE(10+MPIRANK,3000)ITYPE,NTOTAL,J,NEIGL,NEIGR,ERR
               IF (MPIRANK.EQ.0)
     &              WRITE(*,3000)ITYPE,NTOTAL,J,NEIGL,NEIGR,ERR
            ENDIF
            IF (N.LT.10) THEN
               N = N+1
            ELSE IF (N.LT.50) THEN
               N = N+10
            ELSE IF (N.LT.100) THEN
               N = N+25
            ELSE IF (N.LE.NMAX) THEN
               N = N+100
            ENDIF
            GOTO 100
         ENDIF
 200  CONTINUE
C
      CALL MPI_FINALIZE(IERR)
      STOP
 9999 FORMAT(1X,'...... '
     *   /1X,'...... ',6X,' J =',I3,3X,' NEIG =',I3
     *   /1X,'...... '
     *   /1X,'...... ',3X,3X,'  Computed Ritz values',2X,
     *    '(','Error Bnds',')'
     *   /1X,'...... '
     *   /(1X,'...... ',I3,3X,1PD22.14,2X,'(',1PD10.2,')'))
 2000 FORMAT(
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
 2999 FORMAT(1X,'... return flag =',I9,4(' '),'...')
 3000 FORMAT(1X,'ITYPE= ',I1,'  N= ',I6,'  J= ',I4,
     *            '  NEIG= ',I3,I3,'  ACT ERR= ',1PE10.2)
      END
C
      SUBROUTINE OP(N,P,Q,R,MPICOM)
      INTEGER N,MPICOM
      REAL*8 P(N),Q(N),R(N)
C
      INTEGER ITYPE, OFFSET
      COMMON/OPSTIF/ITYPE, OFFSET
C
      INTEGER I,IS
C
      GOTO(10,20,30,40),ITYPE
 10   DO 100 I = 1,N
         R(I) = Q(I)
 100  CONTINUE
      RETURN
 20   DO 200 I = 1,N
         R(I) = DBLE(OFFSET+I)*Q(I)
 200  CONTINUE
      RETURN
 30   IS = (-1)**OFFSET
      DO 300 I = 1,N
         R(I) = DBLE(SIGN(OFFSET+I,IS))*Q(I)
         IS = -IS
 300  CONTINUE
      RETURN
 40   IF (MOD(OFFSET,2).EQ.0) THEN
         IS = 1
      ELSE
         IS = -1
      ENDIF
      DO 400 I = 1,N
         R(I) = DBLE(IS)*Q(I)
         IS = -IS
 400  CONTINUE
      RETURN
      END
C
      SUBROUTINE OPM(N,A,B,MPICOM)
      INTEGER N,MPICOM
      REAL*8 A(N),B(N)
      INTEGER I
C
      DO 100 I = 1,N
         B(I) = A(I)
 100  CONTINUE
      RETURN
      END
C
      REAL*8 FUNCTION EXACT(ITYPE,N,I)
      INTEGER ITYPE,N,I
C
      INTEGER L
C
      GOTO(500,600,700,800),ITYPE
 500  EXACT = DBLE(1)
      RETURN
 600  EXACT = DBLE(I)
      RETURN
 700  L = N/2
      L = L+L
      IF (I.LE.L/2) THEN
         EXACT = DBLE(I+I-L-2)
      ELSE
         EXACT = DBLE(I+I-L-1)
      ENDIF
      RETURN
 800  IF (I.LE.N/2) THEN
         EXACT = DBLE(-1)
      ELSE
         EXACT = DBLE(1)
      ENDIF
      RETURN
      END
