C @(#)plandr2.f	7/31/97; from landr.f 3.16 (BNP) 6/3/89
C
C     SUBROUTINE LANDR(...,J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
      SUBROUTINE pLANDR2(N,LANMAX,astore,MAXPRS,LOHI,CONDM,KAPPA,
     *   J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL,MPICOM)
      INTEGER N,LANMAX,MAXPRS,LOHI,J,NEIG,NW,IERR,MSGLVL,MPICOM
      REAL*8 CONDM,KAPPA,
     *   RITZ(LANMAX),BND(LANMAX),W(NW)
      INCLUDE 'mpif.h'
C
c     this modified version is intended to compute a number of
c     extreme eigenvalues of given OPerator.
c     The alpha and beta arrays are also copied to the beginning
c     of the workspace (W) before returning.
c
C.... The program make a Lanczos run using a linear operator that acts
C.... through a user supplied subroutine called OP in order to compute
C.... first MAXPRS smallest or largest eigenvalues or till maximum
C.... Lanczos steps are taken.  The inner product is defined implicitly
C.... by a user supplied subroutine OPM.
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
C.... N      (local) dimension of the eigenproblem
C.... LANMAX upper limit to the number of Lanczos steps
C.... MAXPRS upper limit to the number of wanted eigenpairs
C.... LOHI   which end of spectrum to compute
C....        .LE.0 lower end of the spectrum, otherwise the high end
C.... CONDM  estimated effective condition number of M
C.... KAPPA  relative accuracy of Ritz values acceptable as eigenvalues
C.... NW     length of the work array W, no less than
C....        5*N+4*LANMAX+max(N,LANMAX+1)+1
C
C.... outputs
C.... J      number of Lanczos steps actually taken
C.... NEIG   number of Ritz values stabilized
C.... RITZ   array to hold the Ritz values
C.... BND    array to hold the error bounds
C.... IERR   error flag
C
C.... Workspace
C.... W      work array of length NW
C....        On entry, first N elements of W is the initial guess
C....        On return, first 2*j elements are alpha, beta
C....         *****     w(1:j) == alpha(1:j)
C....         *****     w(j+1) == 0, w(j+2:j+j)=beta(2:j)
C
C.... NOTE: the j-1 elements produced from the Lanczos algorithm are
C.... stored as beta(2:j) NOT beta(1:j-1) ******
C
C.... subroutines:      MACHAR,LANSO,RITVEC
C
      REAL*8 RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      INTEGER MT,I,NQ(4),N1,L2,L3,L4,L5
C
C.... MACHAR specific (EPS declared in COMMON/RDATA/)
C
      INTEGER IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,MAXEXP,NTOTAL
      REAL*8 EPSNEG,XMIN,XMAX
      DOUBLE PRECISION astore(*)
C
C.... check input data
C
      CALL MPI_ALLREDUCE(N,NTOTAL,1,MPI_INTEGER,MPI_SUM,MPICOM,IERR)
      if (ierr .ne. MPI_SUCCESS) call MPI_ABORT(mpicom, ierr)
      IERR = 0
      MT = 5*N+1+4*LANMAX+MAX(N,LANMAX+1)
      IF (NTOTAL.LE.0) IERR = IERR+1
      IF (LANMAX.LE.0) IERR = IERR+2
      IF (MAXPRS.LE.0) IERR = IERR+8
      IF (MAXPRS.GT.LANMAX) IERR = IERR+16
      IF (LANMAX.GT.NTOTAL) IERR = IERR+32
      IF (MT.GT.NW) IERR = IERR+64
      IF (IERR.GT.0) RETURN
C
C.... compute machine precision
C
      CALL MACHAR(IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,MAXEXP,
     &            EPS,EPSNEG,XMIN,XMAX)
      EPS1 = EPS*SQRT(DBLE(NTOTAL))
      EPSN = EPS*DBLE(NTOTAL)
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
      L2 = 1+5*N
      L3 = LANMAX+L2
      L4 = LANMAX+1+L3
      L5 = LANMAX+L4
      N1 = LANMAX+L5
C     CALL LANSO(...,W,W(N1),W(L2),W(L3),W(L4),W(L5),NQ,IERR,MSGLVL)
      CALL pLANSO2(N,LANMAX,astore,MAXPRS,LOHI,KAPPA,J,NEIG,RITZ,BND,W,W(N1),
     *           W(L2),W(L3),W(L4),W(L5),NQ,IERR,MSGLVL,MPICOM)
      CALL DCOPY(j, w(l2), -1, w, -1)
      CALL DCOPY(j, w(l3), -1, w(j+1), -1)
C
      RETURN
      END
c-----------------------------------------------------------------------
C     SUBROUTINE LANSO(...,R,WRK,ALF,BET,ETA,OLDETA,NQ,IERR,MSGLVL)
      SUBROUTINE pLANSO2(N,LANMAX,astore,MAXPRS,LOHI,KAPPA,J,NEIG,RITZ,BND,
     &   R,WRK,ALF,BET,ETA,OLDETA,NQ,IERR,MSGLVL,MPICOM)
      IMPLICIT NONE
      INTEGER N,LANMAX,MAXPRS,LOHI,J,NEIG,NQ(4),IERR,MSGLVL,MPICOM
      REAL*8 astore(*)
      REAL*8 KAPPA,R(5*N),WRK(*),
     &  ALF(LANMAX),BET(LANMAX+1),RITZ(LANMAX),BND(LANMAX),
     &  ETA(LANMAX),OLDETA(LANMAX)
C
c     This is modified from planso to compute a few extreme eigenvalues
c     of a given operator.  An eigenpair is declared converged if
c     error bound < kappa * max(|Ritz value|, eps^2/3)
c
C.... inputs
C.... N      (local) dimension of the eigenproblem
C.... LANMAX upper limit to the number of Lanczos steps
C.... MAXPRS upper limit to the number of wanted eigenpairs
C.... LOHI   <=0 compute smallest eigenvalue, >0 compute
C....        largest eigenvalues
C.... KAPPA  Relative tolerance on error of Ritz values
C.... MPICOM the MPI communicator for the group of processor
C....        participating in the Lanczos process
C
C.... work space
C.... R      holds 5 vectors of length n. see the text for details.
C.... NQ(4)  contains the pointers to the begining of each vector in R.
C.... ALF    array to hold diagonal of the tridiagonal T
C.... BET    array to hold off-diagonal of T
C.... ETA    orthogonality estimate of lanczos vectors at step j
C.... OLDETA orthogonality estimate of lanczos vectors at step j-1
C.... WRK    miscellaneous usage, size max(N, LANMAX+1)
C
C.... outputs
C.... J      actual number of Lanczos steps taken
C.... NEIG   number of converged Ritz-pairs
C.... RITZ   array to hold the computed Ritz values
C.... BND    array to hold the error bounds
C.... IERR   error flag
C
C.... BLAS routines:    DATX,DAXPY,DCOPY,DDOTMPI,DSCAL
C.... subroutines:      TQLB,ORTBND,pPURGE,pSTARTV,pSTPONE,
C....                   RITZSORT,LANCONV
C.... user-supplied:    OP,OPM,STORE
C
      INTEGER STORQ,RETRQ,STORP,RETRP
      PARAMETER (STORQ = 1,RETRQ = 2,STORP = 3,RETRP = 4)
C
      INTEGER MAXLL
      REAL*8 FOUR
      PARAMETER (MAXLL = 2,FOUR = 4.0D0)
C
      REAL*8 RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      LOGICAL ENOUGH
      INTEGER LL,I,L,FIRST,LAST,MID,ID2,ID3
      REAL*8 T,DDOTMPI,pSTARTV
C
      REAL*8 ONE,ZERO
      DATA ONE,ZERO/1.0D0,0.0D0/
C
      CALL pSTPONE(N,astore,R,WRK,ALF,NQ,MSGLVL,MPICOM)
      J = 1
cMPI  special case no longer applicable
c$$$      IF (N.EQ.1) THEN
c$$$         NEIG = 1
c$$$         RITZ(1) = ALF(1)
c$$$         BND(1) = ZERO
c$$$      ENDIF
C
      ETA(1) = EPS1
      OLDETA(1) = EPS1
C
      LL = 0
      FIRST = 2
      LAST = MIN(MAXPRS+MAX(8,MAXPRS),LANMAX)
      ENOUGH = .FALSE.
 100  IF (FIRST.LE.LANMAX .AND. .NOT. ENOUGH) THEN
         IF (RNM.LE.TOL) RNM = ZERO
C
C....    lanczos loop
C
         DO 10 J = FIRST,LAST
            MID = NQ(2)
            NQ(2) = NQ(1)
            NQ(1) = MID
            MID = NQ(3)
            NQ(3) = NQ(4)
            NQ(4) = MID
            CALL STORE(N,astore,STORQ,J-1,R(NQ(2)))
            IF (J-1.LE.MAXLL) CALL STORE(N,astore,STORP,J-1,R(NQ(4)))
            BET(J) = RNM
C
C...        Restart if invariant subspace is found
C
            IF (BET(J).EQ.ZERO) THEN
               RNM = pSTARTV(N,astore,J,R,WRK,NQ,EPS,MSGLVL,MPICOM)
               ENOUGH = RNM.EQ.ZERO
               IF (ENOUGH) GOTO 15
            ENDIF
C
C....       take a Lanczos step
C
            T = ONE/RNM
            CALL DATX(N,T,R,1,R(NQ(1)),1)
            CALL DSCAL(N,T,R(NQ(3)),1)
C
            CALL OP(N,R(NQ(3)),R(NQ(1)),R,MPICOM)
            IF (BET(J).GT.ZERO) CALL DAXPY(N,-RNM,R(NQ(2)),1,R,1)
C
            ALF(J) = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
            CALL DAXPY(N,-ALF(J),R(NQ(1)),1,R,1)
C
C....       orthogonalize against initial Lanczos vectors
C
            IF (J.LE.MAXLL+1.AND.ABS(ALF(J-1)).GT.FOUR*ABS(ALF(J))) THEN
               LL = J-1
            ENDIF
            DO 5 I = 1,MIN(LL,J-2)
               CALL STORE(N,astore,RETRP,I,WRK)
               T = DDOTMPI(N,WRK,1,R,1,MPICOM)
               CALL STORE(N,astore,RETRQ,I,WRK)
               CALL DAXPY(N,-T,WRK,1,R,1)
               ETA(I) = EPS1
               OLDETA(I) = EPS1
5           CONTINUE
C
C....       extended local reorthogonalization
C
            T = DDOTMPI(N,R,1,R(NQ(4)),1,MPICOM)
            CALL DAXPY(N,-T,R(NQ(2)),1,R,1)
            IF (BET(J).GT.ZERO) BET(J) = BET(J)+T
            T = DDOTMPI(N,R,1,R(NQ(3)),1,MPICOM)
            CALL DAXPY(N,-T,R(NQ(1)),1,R,1)
            ALF(J) = ALF(J)+T
            CALL OPM(N,R,R(NQ(4)),MPICOM)
            RNM = SQRT(DDOTMPI(N,R,1,R(NQ(4)),1,MPICOM))
            ANORM = BET(J)+ABS(ALF(J))+RNM
            TOL = EPSN*ANORM
C
C....       update the orthogonality bounds
C
            CALL ORTBND(J,ALF,BET,ETA,OLDETA)
C
C....       restore the orthogonality state when needed
C
C           CALL PURGE(...,ETA,OLDETA,MSGLVL)
            CALL pPURGE(N,astore,LL,J,R,R(NQ(1)),R(NQ(4)),R(NQ(3)),WRK,
     &         ETA,OLDETA,MSGLVL,MPICOM)
            IF (RNM.LE.TOL) RNM = ZERO
            IF (MSGLVL.GT.10) THEN
               call mpi_comm_rank(mpicom, i, ierr)
               if (i.eq.0) PRINT *, J, ALF(J), BET(J)
            ENDIF
10       CONTINUE
         J = LAST
15       IF (ENOUGH) J = J-1
         FIRST = J+1
         BET(J+1) = RNM
C
C....    Now analyze T
C
         L = 1
         DO 40 ID2 = 1,J
            IF (L.GT.J) GOTO 50
            DO 20 I = L,J
               IF (BET(I+1).EQ.ZERO) GOTO 30
20          CONTINUE
            I = J
C....       Now i is at the end of an unreduced submatrix
C
30          CALL DCOPY(I-L+1,ALF(L),1,RITZ(L),-1)
            IF (I.GT.L) CALL DCOPY(I-L,BET(L+1),1,WRK(L+1),-1)
            CALL TQLB(I-L+1,RITZ(L),WRK(L),BND(L),IERR)
            IF (IERR.NE.0) THEN
               PRINT *,' TQLB failed to converge (ierr =',IERR,')'
               PRINT *,' L =',L,' I =',I
               PRINT *,(ID3,RITZ(ID3),WRK(ID3),BND(ID3),ID3 = L,I)
            ENDIF
            DO 35 ID3 = L,I
               BND(ID3) = RNM*ABS(BND(ID3))
35          CONTINUE
            L = I+1
40       CONTINUE
C
C....    Sort eigenvalues so the wanted ones are in the front
C
50       CALL SORTRITZ(LOHI,J,RITZ,BND)
         call lanconv(j, ritz, bnd, maxprs, eps, kappa, neig)
         IF (MSGLVL.GT.5) THEN
            call mpi_comm_rank(mpicom, i, ierr)
            if (i.eq.0) then
               PRINT *, 'LANSO: STEP ', j, ';  NEIG ', neig,
     &              '   RES ', BND(1)
               if (msglvl.gt.10) then
                  do 60 i = 1, min(msglvl, j)
                     print *, i, ritz(i), bnd(i)
 60               continue
               endif
            endif
         ENDIF
         enough = (neig.ge.maxprs)
C
C....    Should we stop?
C
         IF (.NOT.ENOUGH) THEN
            IF (NEIG.EQ.0) THEN
               LAST = FIRST + MAX(8, MAXPRS)
            ELSE
               LAST = FIRST + MIN(128, MAX(MAXPRS-NEIG,
     &              NINT(0.25*J/NEIG)))
            ENDIF
            IF (3+LAST .ge. LANMAX) LAST = LANMAX
         ELSE
            ENOUGH = .TRUE.
         ENDIF
         ENOUGH = ENOUGH.OR.FIRST.GT.LANMAX
         GOTO 100
      ENDIF
      CALL STORE(N,astore,STORQ,J,R(NQ(1)))
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine sortritz(lohi,j,ritz,bnd)
      implicit none
      integer i, j, lohi
      real*8 ritz(j), bnd(j)
c
      if (lohi.le.0) then
         call dsort2(j,ritz,bnd)
      else
         do 10 i = 1, j
            ritz(i) = -ritz(i)
 10      continue
         call dsort2(j,ritz,bnd)
         do 20 i = 1, j
            ritz(i) = -ritz(i)
 20      continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine lanconv(j, ritz, bnd, maxeig, eps, kappa, neig)
      implicit none
      integer j, neig, i, maxeig
      real*8 ritz(j), bnd(j), eps, kappa, gap, gapl, eps23
c
c     an eigenpair is declared converged if
c     bnd(i) < kappa*(eps^2/3, |ritz(i)|)
c
      neig = 0
      eps23 = exp(log(eps)*0.66666667D0)
      gapl = abs(ritz(j) - ritz(1))
      do 10 i = 1, min(j-1, maxeig)
         gap = gapl
         gapl = abs(ritz(i+1) - ritz(i))
         gap = min(gap, gapl)
         if (gap.gt.bnd(i)) bnd(i) = bnd(i)*(bnd(i)/gap)
         if (bnd(i) .le. kappa*max(eps23, abs(ritz(i)))) neig = neig+1
 10   continue
      return
      end
