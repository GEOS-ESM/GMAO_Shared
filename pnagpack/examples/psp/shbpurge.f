C
C @(#)shbpurge.f	0.01 (KW) 6/3/97; from purge.f 3.16 5/11/89
C
      SUBROUTINE pPURGE(N,LL,J,R,Q,RA,QA,WRK,ETA,OLDETA,MSGLVL,MPICOM)
      implicit none
      include 'mpif.h'
      include 'simphb.h'
      INTEGER N,LL,J,MSGLVL,MPICOM
      REAL*8 R(N),Q(N),RA(N),QA(N),WRK(N),ETA(J),OLDETA(J)
C
C.... This routine examines ETA to decide whether
C.... re-orthogonalization should be performed.
C.... NOTE: this routine uses up to two classic Gram-schmidt process
C.... to restore orthogonality.  If WRK is not large enough, it uses
C.... ETA and OLDETA as workspace.
C
C.... N      (local) dimension of the eigenproblem
C.... LL     no. of initial Lanczos vectors in local orthog.
C.... J      current Lanczos step
C.... R      the residual vector to become the next Lanczos vector
C.... Q      the current Lanczos vector
C.... RA     the product of the mass matrix and r
C.... QA     the product of the mass matrix and q
C.... WRK    a temporary vector to hold the previous Lanczos vectors
C.... ETA    state of orthogonality between r and previous Lanczos vectors
C.... OLDETA state of orthogonality between q and previous Lanczos vectors
C.... MPICOM MPI Communicator
C
C.... BLAS routines:    DAXPY,DDOTMPI,IDAMAX
C.... subroutines:      NONE
C.... user-supplied:    OPM
C
      REAL*8 RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      INTEGER K,I,LOOP,IDAMAX,MY_ID
      REAL*8 T,TQ,TR,REPS1,tmptime,zero,one,DDOTMPI
      parameter (zero=0.0D0, one=1.0D0)
C
      IF (J.LE.LL+1) RETURN
      MY_ID = -1
      K = IDAMAX(J-(LL+1),ETA(LL+1),1)+LL
      IF (MSGLVL.GT.10) THEN
         CALL MPI_COMM_RANK(MPICOM, MY_ID, I)
         IF (MY_ID.EQ.0)
     &        PRINT *, 'PPURGE: Omega(', k, j+1, ') = ', abs(eta(k))
      ENDIF
C
      IF (ABS(ETA(K)).GT.REPS) THEN
c     subpurge - tracking
         tmptime = mpi_wtime()
c
         REPS1 = EPS1/REPS
         DO 55 LOOP = 1,2
            IF (RNM.GT.TOL) THEN
               nreorth = nreorth + 1
               IF (MY_ID.EQ.0 .AND. MSGLVL.GT.10) THEN
                  PRINT *, 'PPURGE: purging Q_', J, ' and R_', J,
     &                 ' LOOP #', LOOP
               ENDIF
C
C....       Bring in a Lanczos vector t and orthogonalize both r and q
C....       against it
C
               TQ = 0.0D0
               TR = 0.0D0
c
c     try to put the two part sums into WRK if possible, otherwise
c     use eta and old eta as workspace
c
               if (n.ge.4*(j-1)) then
                  call dgemv('T',n,j-1,one,qq(1+n+n),n,qa,1,zero,
     &                 wrk(j+j-1),1)
                  call dgemv('T',n,j-1,one,qq(1+n+n),n,ra,1,zero,
     &                 wrk(3*j-2),1)
                  if (mpi_real8.eq.0) then
                     call mpi_allreduce(wrk(j+j-1),wrk,j+j-2,
     &                    mpi_double_precision,mpi_sum,mpicom,i)
                  else
                     call mpi_allreduce(wrk(j+j-1),wrk,j+j-2,mpi_real8,
     &                    mpi_sum,mpicom,i)
                  endif
                  if (I .NE. MPI_SUCCESS)
     &                 CALL MPI_ABORT(MPI_COMM_WORLD,I)
                  call dgemv('N',n,j-1,-one,qq(1+N+N),n,wrk,1,
     &                 one,q,1)
                  call dgemv('N',n,j-1,-one,qq(1+N+N),n,wrk(j),1,
     &                 one,r,1)
                  do 50 i = 1, j-1
                     tq = tq + abs(wrk(i))
                     tr = tr + abs(wrk(j-1+i))
 50               continue
               else
                  call dgemv('T',n,j-1,one,qq(1+n+n),n,qa,1,zero,eta,1)
                  if (mpi_real8.eq.0) then
                     call mpi_allreduce(eta,oldeta,j-1,
     &                    mpi_double_precision,mpi_sum,mpicom,i)
                  else
                     call mpi_allreduce(eta,oldeta,j-1,mpi_real8,
     &                    mpi_sum,mpicom,i)
                  endif
                  if (I .NE. MPI_SUCCESS)
     &                 CALL MPI_ABORT(MPI_COMM_WORLD,I)
                  call dgemv('N',n,j-1,-one,qq(1+N+N),n,oldeta,1,
     &                 one,q,1)
                  do 51 i = 1, j-1
                     tq = tq + abs(oldeta(i))
 51               continue
                  call dgemv('T',n,j-1,one,qq(1+n+n),n,ra,1,zero,eta,1)
                  if (mpi_real8.eq.0) then
                     call mpi_allreduce(eta,oldeta,j-1,
     &                    mpi_double_precision,mpi_sum,mpicom,i)
                  else
                     call mpi_allreduce(eta,oldeta,j-1,mpi_real8,
     &                    mpi_sum,mpicom,i)
                  endif
                  if (I .NE. MPI_SUCCESS)
     &                 CALL MPI_ABORT(MPI_COMM_WORLD,I)
                  call dgemv('N',n,j-1,-one,qq(1+N+N),n,oldeta,1,
     &                 one,r,1)
                  do 52 i = 1, j-1
                     tr = tr + abs(oldeta(i))
 52               continue
               endif
               CALL OPM(N,Q,QA,MPICOM)
C
C....          restore local orthogonality
C
               T = -DDOTMPI(N,R,1,QA,1,MPICOM)
               TR = TR+ABS(T)
               CALL DAXPY(N,T,Q,1,R,1)
C
               CALL OPM(N,R,RA,MPICOM)
               RNM = SQRT(DDOTMPI(N,RA,1,R,1,MPICOM))
               IF (TQ.LE.REPS1.AND.TR.LE.REPS1*RNM) GOTO 80
            ENDIF
55       CONTINUE
80       DO 60 I = 1,J
            ETA(I) = EPS1
            OLDETA(I) = EPS1
60       CONTINUE
c     subpurge tracking
         reorthtime = reorthtime + mpi_wtime() - tmptime
c
      ENDIF
      RETURN
      END
