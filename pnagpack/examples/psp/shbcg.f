      program shbcg
      implicit none
c
c     This program performs the following tasks:
c     (1) It constructs a simple symmetric eigenvalue problem by
c     reading a matrix stored in Boeing-Harwell format.  The matrix
c     is read by one processor and distributed to all processors in
c     the group.  If the matrix is RSA type, it is expanded before
c     distribution, if it is RUA type, the symmetric part of the matrix
c     is used.
c     (2) It calls LANDR once using MEMSIZE memory as workspace to
c     compute as much eigenvalues as possible using shift-and-invert
c     scheme.  Eigenvectors are not computed, Lanczos vectors are
c     stored only to perform reorthogonalization.
c     (3) The invert operator is solved using preconditioned CG from
c     SPARSKIT. The matrix-vector multiplication routine from PSPARSLIB
c     is used.
c     *** NOTE that the matrix (A - delta I) must be symmetric ***
c     *** positive definite!                                   ***
c     Dynamic memory allocation is used through-out the preprocessing
c     stages.
c
      include 'mpif.h'
c
c     the arrays for the distrubited MATVEC (amxdis) is stored in
c     the following header file
c
      include 'shbcg.h'
c     .. Parameters ..
      integer            memsize
      parameter          (memsize = 2000000)
      double precision   one
      parameter          (one = 1.0D0)
c     ..
c     .. Local Scalars ..
      character*128      delstr, hbfile
      integer            i, ierr, iou, j, maxlan, maxnloc, maxnzloc,
     &                   mpicom, neig, nloc, nw
      double precision   condm, endl, endr, kappa, lanso_time,
     &                   setup_time
c     ..
c     .. Local Arrays ..
      integer            ord(*)
      double precision   bnd(*), ritz(*), wrk(*)
      pointer            (pwrk, wrk), (pord,ord), (pritz, ritz),
     &                   (pbnd, bnd)
c     ..
c     .. External Functions ..
      integer            iargc
      external           iargc
c     ..
c     .. External Subroutines ..
      external           dmalloc, dmfree, genprec, hbread1,
     &                   plandr2, psp_setup
c     ..
c     .. Statement Functions ..
      integer            ncells
c     ..
c     .. Statement Function definitions ..
      ncells(nloc, i) = 5*nloc + 4*i + 1 + max(nloc, i+1)
c     ..
c     .. Executable Statements ..
c
c.... A few constants
c
      verbose = 3
      locfact = 2
      kjac = 5
      iou = 17
      amxtype = 29
      ptrn = 1
      kappa = 1.0D-8
      delta = 0.0D0
      lusize = 0
c
      paloc = 0
      pjaloc = 0
      pialoc = 0
      pix = 0
      pord = 0
      pqq = 0
      pwrk = 0
      pritz = 0
      pbnd = 0
      palu = 0
      pjlu = 0
      pju = 0
      kmvcnt = 0
      mmvcnt = 0
      nreorth = 0
      kmvtime = 0
      mmvtime = 0
      reorthtime = 0
c
c
c     start MPI
c
      call mpi_init(ierr)
      if (ierr.ne.mpi_success)stop 'FAILED TO INIT MPI'
      mpicom = mpi_comm_world
      call mpi_comm_size(mpicom, totprocs, ierr)
      call mpi_comm_rank(mpicom, myproc, ierr)
      if (totprocs.gt.max_procs) then
         print *, 'Increase MAX_PROCS in file shbcg.h and PSPARSLIB.'
         call mpi_abort(mpicom, ierr)
      endif
      if (myproc.eq.0)print *, 'SHBCG'
c
c     the file name can be given as a command line argument
c
      if (iargc().gt.0) then
         call getarg(1, hbfile)
         if (iargc().gt.1) then
            call getarg(2, delstr)
            read (delstr, FMT = *, iostat = ierr)delta
            if (ierr.ne.0) then
               if (myproc.eq.0)print *,
     &             'The second argument must be a number.'
               delta = 0.0D0
            endif
         else
            delta = 0.0D0
         endif
         if (iargc().gt.2) then
            call getarg(3, delstr)
            read (delstr, FMT = *, iostat = ierr)kjac
            if (ierr.ne.0) then
               if (myproc.eq.0)print *,
     &             'The third argument must be a number.'
               kjac = 5
            endif
         else
            kjac = 5
         endif
      else
         print *, 'USAGE: shbcg hbfilename delta kjac'
         stop
      endif
c
c     read the given file and distribute the matrix
c
      setup_time = mpi_wtime()
      call hbread1(hbfile, 0, myproc, totprocs, nloc, split, paloc,
     &             pjaloc, pialoc, maxnzloc, maxnloc, iou, mpicom, ierr,
     &             verbose)
      if (ierr.ne.0) then
         print *, 'HBREAD failed with error code ', ierr
         call mpi_abort(mpicom, ierr)
      endif
c
c     prepare the matrix for distributed matrix-vector multiplication
c
      call psp_setup(maxnloc, maxnzloc, nloc, myproc, totprocs, split,
     &               aloc, jaloc, ialoc, pord, nbnd, nproc, proc, pix,
     &               ipr, mpicom, ierr)
      if (ierr.ne.0) then
         print *, 'PSP_SETUP failed with error code ', ierr
         call mpi_abort(mpicom, ierr)
      endif
      if (verbose.gt.15) then
         write (10+myproc, *) 'Local components on PE', myproc
         call dump(1, nloc, .true., aloc, jaloc, ialoc, 10+myproc)
         write (10+myproc, *) 'External edges:'
         call dump(nloc+1, nloc+nloc-nbnd+1, .true., aloc, jaloc, ialoc,
     &        10+myproc)
      endif
c
c     generate preconditioner for CG
c
      call genprec(nloc, ierr)
      if (ierr.ne.0) then
         print *, 'GENPREC returns error ', ierr
         print *, '        NO preconditioning will be used.'
         if (palu.ne.0) call dmfree(palu)
         if (pjlu.ne.0) call dmfree(pjlu)
         if (pju.ne.0) call dmfree(pju)
         palu = 0
         pjlu = 0
         pju = 0
         locfact = -1
      endif
c
c     the workspace for LANSO is allocated to fill up to
c     the maximum memory size available
c
      i = 4*maxnzloc + 10*maxnloc + 8*nloc + lusize
      if (totprocs.gt.1 .and. pix.ne.0) then
         i = i + 3*(nloc+max(ipr(nproc+1),ipr(nproc+2)))
      endif
      endl = dble(memsize-i) / dble(nloc+5)
      endl = max(2.0D1, endl)
      if (mpi_real8.eq.0) then
         call mpi_allreduce(endl, endr, 1, mpi_double_precision,
     &                      mpi_min, mpicom, ierr)
      else
         call mpi_allreduce(endl, endr, 1, mpi_real8, mpi_min, mpicom,
     &                      ierr)
      endif
      if (ierr.ne.mpi_success) then
         print *, 'SHBCG: failed to compute MAXLAN.'
         call mpi_abort(mpicom, ierr)
      endif
      maxlan = min(int(endr), split(totprocs+1)-1)
      nw = ncells(nloc, maxlan)
      call dmalloc(pqq, nloc*(maxlan+2), 8)
      call dmalloc(pwrk, nw, 8)
      call dmalloc(pritz, maxlan, 8)
      call dmalloc(pbnd, maxlan, 8)
      endl = -1.0D-30
      endr = 1.0D-30
      condm = one
      if (myproc.eq.0 .and. verbose.gt.0) then
         write(*, FMT = '(/A$)')'SHBCG: '
         print *, 'NPROCS = ', totprocs, '      MAXLAN = ', maxlan,
     &      '    Shift = ', delta, '    Kjac = ', kjac
         if (verbose.gt.10) then
            print *, 'Matrix distribution:', split(1),
     &         ('--', split(i), i = 2, totprocs+1)
         endif
      endif
      if (pqq.eq.0 .or. pwrk.eq.0 .or. pritz.eq.0 .or. pbnd.eq.0) then
         print *, 'SHBCG failed to allocate workspace on PE', myproc
         call MPI_abort(mpicom, ierr)
      endif
c
c     the initial guess for the eigenvector is always [1,1,...,1]^T
c
      do 10 i = 1, nloc
         wrk(i) = one
 10   continue
      lanso_time = mpi_wtime()
      setup_time = lanso_time - setup_time
c
c     main loop to compute the smallest eigenvalue
c
      j = 0
      call plandr2(nloc, maxlan, 5, 1, condm, kappa, j, neig, ritz, bnd,
     &             wrk, nw, ierr, verbose, mpicom)
      if (ierr.ne.0) then
         print *, 'LANDR2 returned with error code = ', ierr
         if (j.le.0)call mpi_abort(mpi_comm_world, ierr)
      endif
c
      lanso_time = mpi_wtime() - lanso_time
      call dmfree(pwrk)
      call dmfree(pqq)
      if (mpi_real8.eq.0) then
         call mpi_reduce(setup_time, endl, 1, mpi_double_precision,
     &                   mpi_sum, 0, mpicom, ierr)
         call mpi_reduce(lanso_time, endr, 1, mpi_double_precision,
     &                   mpi_sum, 0, mpicom, ierr)
      else
         call mpi_reduce(setup_time, endl, 1, mpi_real8, mpi_sum, 0,
     &                   mpicom, ierr)
         call mpi_reduce(lanso_time, endr, 1, mpi_real8, mpi_sum, 0,
     &                   mpicom, ierr)
      endif
      condm = max(abs(ritz(1)), abs(ritz(j)))
      call dsort2a(j, ritz, bnd)
c
      if (myproc.eq.0) then
         setup_time = endl / totprocs
         lanso_time = endr / totprocs
         print *, 'LANCZOS STEPS = ', j, '    NEIG = ',
     &      neig, '    || OP || = ', condm
         print *, 'Setup time = ', setup_time, 'sec'
         print *, 'LANSO time = ', lanso_time, 'sec'
         print *, 'MATVEC = ', kmvcnt, mmvcnt, '  time ', kmvtime,
     &      mmvtime
         print *, 'N reorth = ', nreorth, '  time ', reorthtime
         if (verbose.gt.3) then
            condm = condm*kappa
            write(*, FMT = 9999)
            write(*, FMT = 9989)1, delta+one/ritz(1), bnd(1)
            do 20 i = 2, j - 1
               if (bnd(i).lt.condm)write(*, FMT = 9989)i,
     &              delta+one/ritz(i), bnd(i)
 20         continue
            write(*, FMT = 9989)j, delta+one/ritz(j), bnd(j)
         endif
      endif
 9999 format(1X, '   INDEX       Lambda(LANSO)       ErrorBound')
 9989 format(1X, I8, G25.17, 2X, 1P, E10.2)
c
      call mpi_finalize(ierr)
c     end of program shbcg
      end
c-----------------------------------------------------------------------
c     modified from store.f 3.3 (BNP) 3/16/89
c
      subroutine store(n, isw, j, s)
      implicit none
      INCLUDE 'shbcg.h'
c     .. Parameters ..
      integer            maxll
      parameter          (maxll = 2)
      integer            storq, retrq, storp, retrp
      parameter          (storq = 1, retrq = 2, storp = 3, retrp = 4)
c     ..
c     .. Scalar Arguments ..
      integer            isw, j, n
c     ..
c     .. Array Arguments ..
      double precision   s(n)
c     ..
c     .. External Subroutines ..
      external           dcopy
c     ..
c     .. Executable Statements ..
c
      if (isw.eq.storq) then
         call dcopy(n, s, 1, qq((j+maxll-1)*n+1), 1)
      else if (isw.eq.retrq) then
         call dcopy(n, qq((j+maxll-1)*n+1), 1, s, 1)
      else if (isw.eq.storp) then
         if (j.gt.maxll)stop 'STORE: (STORP) J.GT.MAXLL'
         call dcopy(n, s, 1, qq((j-1)*n+1), 1)
      else if (isw.eq.retrp) then
         if (j.gt.maxll)stop 'STORE: (RETRP) J.GT.MAXLL'
         call dcopy(n, qq((j-1)*n+1), 1, s, 1)
      endif
      return
      end
c-----------------------------------------------------------------------
c @(#)shbpurge.f        0.01 (KW) 6/3/97; from purge.f 3.16 5/11/89
c
      subroutine ppurge(n, ll, j, r, q, ra, qa, wrk, eta, oldeta,
     &                  msglvl, mpicom)
      implicit none
      include 'mpif.h'
      include 'shbcg.h'
c
c.... This routine examines ETA to decide whether
c.... re-orthogonalization should be performed.
c.... NOTE: this routine uses up to two classic Gram-schmidt process
c.... to restore orthogonality.  If WRK is not large enough, it uses
c.... ETA and OLDETA as workspace.
c
c.... N      (local) dimension of the eigenproblem
c.... LL     no. of initial Lanczos vectors in local orthog.
c.... J      current Lanczos step
c.... R      the residual vector to become the next Lanczos vector
c.... Q      the current Lanczos vector
c.... RA     the product of the mass matrix and r
c.... QA     the product of the mass matrix and q
c.... WRK    a temporary vector to hold the previous Lanczos vectors
c.... ETA    orthogonality between r and previous Lanczos vectors
c.... OLDETA orthogonality between q and previous Lanczos vectors
c.... MPICOM MPI Communicator
c
c.... BLAS routines:    DAXPY,DDOTMPI,IDAMAX
c.... subroutines:      DAXPY, DGEMV
c.... user-supplied:    OPM
c
c     .. Parameters ..
      double precision   zero, one
      parameter          (zero = 0.0D0, one = 1.0D0)
c     ..
c     .. Scalar Arguments ..
      integer            j, ll, mpicom, msglvl, n
c     ..
c     .. Array Arguments ..
      double precision   eta(j), oldeta(j), q(n), qa(n), r(n), ra(n),
     &                   wrk(n)
c     ..
c     .. Scalars in Common ..
      double precision   anorm, eps, eps1, eps34, epsn, rceps1, reps,
     &                   rnm, tol
c     ..
c     .. Local Scalars ..
      integer            i, k, loop, my_id
      double precision   reps1, t, tmptime, tq, tr
c     ..
c     .. External Functions ..
      integer            idamax
      double precision   ddotmpi
      external           idamax, ddotmpi
c     ..
c     .. External Subroutines ..
      external           daxpy, dgemv, opm
c     ..
c     .. Common blocks ..
      common             / rdata / rnm, anorm, tol, eps, eps1, rceps1,
     &                   epsn, reps, eps34
c     ..
c     .. Executable Statements ..
c
      if (j.le.ll+1)return
      my_id = -1
      k = idamax(j-(ll+1), eta(ll+1), 1) + ll
c
      if (abs(eta(k)).gt.reps) then
         if (msglvl.gt.10) then
            call mpi_comm_rank(mpicom, my_id, i)
            if (my_id.eq.0)print *, 'PPURGE: Omega(', k, j + 1, ') = ',
     &          abs(eta(k))
         endif
c     subpurge - tracking
         tmptime = mpi_wtime()
c
         reps1 = eps1 / reps
         do 40 loop = 1, 2
            if (rnm.gt.tol) then
               nreorth = nreorth + 1
               if (my_id.eq.0 .and. msglvl.gt.10) then
                  print *, 'PPURGE: purging Q_', j, ' and R_', j,
     &               ' LOOP #', loop
               endif
c
c....       Bring in a Lanczos vector t and orthogonalize both r and q
c....       against it
c
               tq = 0.0D0
               tr = 0.0D0
c
c     try to put the two part sums into WRK if possible, otherwise
c     use eta and old eta as workspace
c
               if (n.ge.4*(j-1)) then
                  call dgemv('T', n, j-1, one, qq(1+n+n), n, qa, 1,
     &                       zero, wrk(j+j-1), 1)
                  call dgemv('T', n, j-1, one, qq(1+n+n), n, ra, 1,
     &                       zero, wrk(3*j-2), 1)
                  if (mpi_real8.eq.0) then
                     call mpi_allreduce(wrk(j+j-1), wrk, j+j-2,
     &                                  mpi_double_precision, mpi_sum,
     &                                  mpicom, i)
                  else
                     call mpi_allreduce(wrk(j+j-1), wrk, j+j-2,
     &                                  mpi_real8, mpi_sum, mpicom, i)
                  endif
                  if (i.ne.mpi_success)call mpi_abort(mpi_comm_world, i)
                  call dgemv('N', n, j-1, -one, qq(1+n+n), n, wrk, 1,
     &                       one, q, 1)
                  call dgemv('N', n, j-1, -one, qq(1+n+n), n, wrk(j), 1,
     &                       one, r, 1)
                  do 10 i = 1, j - 1
                     tq = tq + abs(wrk(i))
                     tr = tr + abs(wrk(j-1+i))
 10               continue
               else
                  call dgemv('T', n, j-1, one, qq(1+n+n), n, qa, 1,
     &                       zero, eta, 1)
                  if (mpi_real8.eq.0) then
                     call mpi_allreduce(eta, oldeta, j-1,
     &                                  mpi_double_precision, mpi_sum,
     &                                  mpicom, i)
                  else
                     call mpi_allreduce(eta, oldeta, j-1, mpi_real8,
     &                                  mpi_sum, mpicom, i)
                  endif
                  if (i.ne.mpi_success)call mpi_abort(mpi_comm_world, i)
                  call dgemv('N', n, j-1, -one, qq(1+n+n), n, oldeta, 1,
     &                       one, q, 1)
                  do 20 i = 1, j - 1
                     tq = tq + abs(oldeta(i))
 20               continue
                  call dgemv('T', n, j-1, one, qq(1+n+n), n, ra, 1,
     &                       zero, eta, 1)
                  if (mpi_real8.eq.0) then
                     call mpi_allreduce(eta, oldeta, j-1,
     &                                  mpi_double_precision, mpi_sum,
     &                                  mpicom, i)
                  else
                     call mpi_allreduce(eta, oldeta, j-1, mpi_real8,
     &                                  mpi_sum, mpicom, i)
                  endif
                  if (i.ne.mpi_success)call mpi_abort(mpi_comm_world, i)
                  call dgemv('N', n, j-1, -one, qq(1+n+n), n, oldeta, 1,
     &                       one, r, 1)
                  do 30 i = 1, j - 1
                     tr = tr + abs(oldeta(i))
 30               continue
               endif
               call opm(n, q, qa, mpicom)
c
c....          restore local orthogonality
c
               t = -ddotmpi(n, r, 1, qa, 1, mpicom)
               tr = tr + abs(t)
               call daxpy(n, t, q, 1, r, 1)
c
               call opm(n, r, ra, mpicom)
               rnm = sqrt(ddotmpi(n,ra,1,r,1,mpicom))
               if (tq.le.reps1 .and. tr.le.reps1*rnm) go to 50
            endif
 40      continue
 50      do 60 i = 1, j
            eta(i) = eps1
            oldeta(i) = eps1
 60      continue
c     subpurge tracking
         reorthtime = reorthtime + mpi_wtime() - tmptime
c
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine genprec(nloc, ierr)
      implicit none
      include 'shbcg.h'
c     .. Scalar Arguments ..
      integer            ierr, nloc
c     ..
c     .. Local Scalars ..
      integer            i, ii, j, knz
c     ..
c     .. Local Arrays ..
      integer            iwk(*)
      pointer            (piwk, iwk)
c     ..
c     .. External Subroutines ..
      external           csort, dmalloc, dmfree, ilu0, slu_fact
c     ..
c     .. Executable Statements ..
      ierr = 0
c
c     call SuperLU
c
      if (locfact.gt.1) then
         locfact = 2
         call slu_fact(nloc, aloc, jaloc, ialoc, ierr)
         if (ierr.gt.0) then
            lusize = ierr
            ierr = 0
            go to 40
         else
            print *, 'GENPREC: SuperLU failed on PE', myproc,
     &         ' Use ILU0 as local factorization.'
            locfact = 1
         endif
      endif
c
c     generate ilu0
c
      if (locfact.eq.0) then
         knz = nloc+1
      else
         knz = ialoc(nloc+1) - ialoc(1) + nloc
      endif
      call dmalloc(palu, knz, 8)
      call dmalloc(pjlu, knz, 1)
      call dmalloc(pju, nloc, 1)
      if (palu.eq.0 .or. pjlu.eq.0 .or. pju.eq.0) then
         print *, 'GENPREC failed to allocate space for LU on PE',
     &        myproc
         stop
      endif
      lusize= knz + knz + nloc
c
      if (locfact.eq.0) then
         locfact = 0
         j = nloc+2
         do 20 i = 1, nloc
            alu(i) = 0.0D0
            jlu(i) = j
            ju(i) = j
            do 10 ii = ialoc(i), ialoc(i+1)-1
               if (jaloc(ii).eq.i .and. aloc(ii).ne.0.0D0)
     &              alu(i) = 1.0D0/aloc(ii)
 10         continue
 20      continue
         i = nloc + 1
         alu(i) = 0.0D0
         jlu(i) = j
      else
         locfact = 1
         call dmalloc(piwk, nloc, 1)
         if (piwk.eq.0) then
            print *, 'GENPREC failed to allocate workspace for ILU0',
     &           ' on PE', myproc
            stop
         endif
         call ilu0(nloc, aloc, jaloc, ialoc, alu, jlu, ju, iwk, ierr)
         call dmfree(piwk)
         if (ierr.ne.0) then
            print *, 'ILU0 encountered zero pivot at step ', ierr
            call dmfree(palu)
            call dmfree(pjlu)
            call dmfree(pju)
            lusize = 0
         endif
      endif
c
 9999 format(1X, 'Local factorization: ', A, ' on PE', I3)
 40   if (verbose.gt.2) then
         if (locfact.eq.0) then
            print 9999, 'diagonal', myproc
         else if (locfact.eq.1) then
            print 9999, 'ILU0', myproc
         else if (locfact.eq.2) then
            print 9999, '(Super)LU', myproc
         else
            print 9999, 'Unknown', myproc
         endif
      endif
c
      return
c     end subroutine genprec
      end
c-----------------------------------------------------------------------
C @(#) modified from dsort2.f	3.2 (BNP) 12/9/88
C
      SUBROUTINE DSORT2A(N,ARRAY1,ARRAY2)
      INTEGER N
      REAL*8 ARRAY1(0:N-1),ARRAY2(0:N-1)
C
C.... Sort array1 and array2 in decreasing order for array1's magnitude
C
      INTEGER IGAP,I,J
      REAL*8 TEMP
C
      IGAP = N/2
 10   IF (IGAP.GT.0) THEN
        DO 200 I = IGAP,N-1
          J = I-IGAP
 50       IF (abs(ARRAY1(J)).LT.abs(ARRAY1(J+IGAP))) THEN
              TEMP = ARRAY1(J)
              ARRAY1(J) = ARRAY1(J+IGAP)
              ARRAY1(J+IGAP) = TEMP
              TEMP = ARRAY2(J)
              ARRAY2(J) = ARRAY2(J+IGAP)
              ARRAY2(J+IGAP) = TEMP
          ELSE
            GO TO 200
          ENDIF
          IF (J.GE.IGAP) THEN
             J = J-IGAP
             GO TO 50
          ENDIF
 200    CONTINUE
      ELSE
        RETURN
      ENDIF
      IGAP = IGAP/2
      GO TO 10
      END
c-----------------------------------------------------------------------
      double precision    function distdot(n, x, ix, y, iy)
c     distdot is required by CG from SPARSKIT
      include 'mpif.h'
c     .. Scalar Arguments ..
      integer            ix, iy, n
c     ..
c     .. Array Arguments ..
      double precision   x(*), y(*)
c     ..
c     .. Local Scalars ..
      integer            ierr
c     ..
c     .. Local Arrays ..
      double precision   tmp(2)
c     ..
c     .. External Functions ..
      double precision   ddot
      external           ddot
c     ..
c     .. Executable Statements ..
      tmp(2) = ddot(n, x, ix, y, iy)
      if (mpi_real8.gt.0) then
         call mpi_allreduce(tmp(2), tmp(1), 1, mpi_real8, mpi_sum,
     &                      mpi_comm_world, ierr)
      else
         call mpi_allreduce(tmp(2), tmp(1), 1, mpi_double_precision,
     &                      mpi_sum, mpi_comm_world, ierr)
      endif
      if (ierr.ne.mpi_success) then
         call mpi_abort(mpi_comm_world, ierr)
      endif
      distdot = tmp(1)
      return
c     end distdot
      end
c-----------------------------------------------------------------------
      subroutine opm(n, x, y, mpicom)
      implicit none
      include 'mpif.h'
      include 'shbcg.h'
c     .. Scalar Arguments ..
      integer            mpicom, n
c     ..
c     .. Array Arguments ..
      double precision   x(n), y(n)
c     ..
c     .. Local Scalars ..
      integer            i
      double precision   tmptime
c     ..
c     .. Executable Statements ..
      mmvcnt = mmvcnt + 1
      tmptime = mpi_wtime()
      do 10 i = 1, n
         y(i) = x(i)
 10   continue
      mmvtime = mmvtime + mpi_wtime() - tmptime
      return
c     end subroutine opm
      end
c-----------------------------------------------------------------------
      subroutine op(n, p, q, r, mpicom)
      implicit none
      include 'mpif.h'
      include 'shbcg.h'
c     .. Scalar Arguments ..
      integer            mpicom, n
c     ..
c     .. Array Arguments ..
      double precision   p(n), q(n), r(n)
c     ..
c     .. Local Scalars ..
      integer            i, ii, j, length
      double precision   flops, tmptime, res
c     ..
c     .. Local Arrays ..
      integer            ipar(16)
      double precision   fpar(16), itwrk(*), xtmp(*), ytmp(*)
      pointer            (pxtmp, xtmp), (pytmp, ytmp), (pitwrk, itwrk)
c     ..
c     .. External Subroutines ..
      double precision   ddot
      external           amux, amxdis, cg, daxpy, dcopy, dmalloc,
     &                   lusol, ddot
c     ..
c     .. Save statement ..
      save               pitwrk, pxtmp, pytmp, length, flops
c     ..
c     .. Data statements ..
      data               pitwrk, pxtmp, pytmp, length, flops
     &                   / 4*0, 0.0D0 /
c     ..
c     .. Executable Statements ..
c
c     this routine performs the following solve
c     (A - delta I) r = q
c     always use zero as initial guess
c
      tmptime = mpi_wtime()
c
c     allocte enough temporary workspace for amxdis and CG
c
      if (pitwrk.eq.0) then
         call dmalloc(pitwrk, 5*n, 8)
         if (pitwrk.eq.0) then
            print *, 'OP failed to allocate workspace for CG on PE',
     &           myproc
            call MPI_abort(mpicom, i)
         endif
      endif         
      if (totprocs.gt.1 .and. pxtmp.eq.0) then
         length = n + max(ipr(nproc+1), ipr(nproc+2))
         call dmalloc(pxtmp, length, 8)
         call dmalloc(pytmp, length, 8)
         if (pxtmp.eq.0 .or. pytmp.eq.0) then
            print *, 'OP failed to allocate workspace for AMXDIS on PE',
     &           myproc
            call MPI_abort(mpicom, i)
         endif
      endif
c
      ipar(1) = 0
      if ((locfact.eq.2 .or. (palu.ne.0 .and. pjlu.ne.0 .and. pju.ne.0))
     &     .and. kjac.gt.0) then
         ipar(2) = 2
         call MPI_Allreduce(ipar(2), ipar(3), 1, mpi_integer, mpi_sum,
     &        mpicom, i)
         if (ipar(3).ne.2*totprocs) then
            ipar(2) = 0
            locfact = -1
            if (palu.ne.0) then
               call dmfree(palu)
               call dmfree(pjlu)
               call dmfree(pju)
               lusize=0
            endif
         endif
      else
         ipar(2) = 0
      endif
      ipar(3) = 0
      ipar(4) = 5*n
      ipar(5) = 10
      ipar(6) = 100
      ipar(16) = 0
      fpar(1) = 1.0D-8
      fpar(2) = 0.0D0
      fpar(11) = 0.0D0
      do 30 i = 1, n
         r(i) = 0.0D0
 30   continue
 10   call cg(n, q, r, ipar, fpar, itwrk)
      if (ipar(1).eq.1) then
         if (totprocs.gt.1) then
            call dcopy(n, itwrk(ipar(8)), 1, xtmp, 1)
            call amxdis(n, nbnd, xtmp, ytmp, aloc, jaloc, ialoc, nproc,
     &                  proc, ix, ipr, amxtype, length, ptrn)
            if (delta.ne.0.0D0) then
               ipar(8) = ipar(8) - 1
               ipar(9) = ipar(9) - 1
c$DIR unroll=16
               do 20 i = 1, n
                  itwrk(ipar(9)+i) = ytmp(i) - delta*itwrk(ipar(8)+i)
 20            continue
               ipar(8) = ipar(8) + 1
               ipar(9) = ipar(9) + 1
            else
               call dcopy(n, ytmp, 1, itwrk(ipar(9)), 1)
            endif
            if (xtmp(length+1).ne.0 .or. ytmp(length+1).ne.0) then
               print *, 'Guard elements of xtmp and ytmp are:',
     &            xtmp(length+1), ytmp(length+1)
            endif
         else
            call amux(n, itwrk(ipar(8)), itwrk(ipar(9)), aloc, jaloc,
     &                ialoc)
            if (delta.ne.0.0D0) call daxpy(n, -delta,
     &           itwrk(ipar(8)), 1, itwrk(ipar(9)), 1)
         endif
         fpar(11) = fpar(11) + ialoc(n+1) + ialoc(n+1)
         kmvcnt = kmvcnt + 1
         if (ipar(7).eq.ipar(6)-1) then
            if (fpar(5).gt.fpar(4) .and. fpar(5).lt.fpar(3)) then
               ipar(6) = min(ipar(6)+ipar(6), split(totprocs+1))
            endif
            if (verbose.gt.3 .and. myproc.eq.0) then
               print 9999, ipar(7), fpar(5), fpar(5)/fpar(3)
            endif
         endif
         go to 10
      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
         if (locfact.gt.1) then
            call dcopy(n, itwrk(ipar(8)), 1, itwrk(ipar(9)), 1)
            call SLU_SOLVE(n, 1, n, itwrk(ipar(9)), fpar(11))
         else
            call lusol(n, itwrk(ipar(8)), itwrk(ipar(9)), alu, jlu, ju)
            fpar(11) = fpar(11) + jlu(n+1) + jlu(n+1)
         endif
         if (kjac.le.1 .or. nbnd.gt.n) go to 10
         j = 1
 100     call dcopy(n, itwrk(ipar(9)), 1, xtmp, 1)
         call msg_bdx_sendi(n, xtmp, ytmp, nproc, proc, ix, ipr, ptrn,
     &        i)
         if (i.ne.mpi_success)call mpi_abort(mpi_comm_world, i)
c
c     check the convergence so far to see whether it is necessary to
c     increase kjac
c
         if (j+1.eq.kjac) then
            call amux(n, xtmp, ytmp, aloc, jaloc, ialoc)
            call amux1(n-nbnd+1, xtmp, ytmp(nbnd), aloc, jaloc,
     &           ialoc(n+1))
            ipar(8) = ipar(8) - 1
            do 40 i = 1, n
               ytmp(i) = ytmp(i) - itwrk(ipar(8)+i) - delta*xtmp(i)
 40         continue
            ipar(8) = ipar(8) + 1
            fpar(15) = ddot(n, itwrk(ipar(8)), 1, itwrk(ipar(8)), 1)
            fpar(16) = ddot(n, ytmp, 1, ytmp, 1)
            call mpi_allreduce(fpar(15), ytmp, 2, mpi_real8, mpi_sum,
     &           mpicom, i)
            if (ytmp(2).lt.ytmp(1)) then
               i = kjac
               if (ytmp(1).gt.0.0D0 .and. ytmp(2).gt.1D-2*ytmp(1)
     &              .and. ipar(16).eq.0) then
                  kjac = 1+int(j*log(1D-2)/log(ytmp(2)/ytmp(1)))
                  kjac = min(100, kjac)
               endif
               if (myproc.eq.0 .and. i.ne.kjac .and. verbose.gt.1) then
                  print *, 'OP: modified KJAC from ', i, ' to ', kjac
               endif
            else
               ipar(16) = ipar(16) - 1
               if (verbose.gt.2 .and. myproc.eq.0) then
                  print 9989, sqrt(ytmp(1)), j, sqrt(ytmp(2)), ipar(7)
               endif
               if (ipar(16).le.nint(-0.8*ipar(7)) .and.
     &              ipar(16).lt.-5) then
                  if (verbose.gt.0 .and. myproc.eq.0) then
                     print *, 'OP: restart CG without preconditioning.'
                  endif
                  ipar(1) = 0
                  ipar(2) = 0
                  goto 10
               endif
            endif
         endif
 9989    format(1X, 'Block-Jacobi failing: res_0', 1PE12.3,
     &        ', res_', I3.3, 1PE12.3, ', CG step', I8)
         do 70 i = 1, n
            res = itwrk(ipar(8)+i-1)
            if (i.ge.nbnd) then
               do 60 ii = ialoc(n+i-nbnd+1), ialoc(n+i-nbnd+2)-1
                  res = res - aloc(ii)*xtmp(jaloc(ii))
 60            continue
            endif
            itwrk(ipar(9)+i-1) = res
 70      continue
         fpar(11) = fpar(11) + 2*(ialoc(n+n-nbnd+2)-ialoc(n+1))
         if (locfact.gt.1) then
            call SLU_SOLVE(n, 1, n, itwrk(ipar(9)), fpar(11))
         else
            call lusol(n, itwrk(ipar(9)), itwrk(ipar(9)), alu, jlu, ju)
            fpar(11) = fpar(11) + jlu(n+1) + jlu(n+1)
         endif
         j = j + 1
         if (j.lt.kjac) then
            goto 100
         else
            goto 10
         endif
      else if (ipar(1).ne.0 .and. myproc.eq.0) then
         print *, 'OP: iterative method returned with error code',
     &      ipar(1)
      endif
c
 9999 format(1X, 'CG: MATVEC', I6, '   res=', 1PE10.2,
     &     '   reduction=', 1PE10.2)
      if (verbose.gt.1 .and. myproc.eq.0) then
         print 9999, ipar(7), fpar(5), fpar(5)/fpar(3)
      endif
      if (ipar(1).lt.-1 .or. fpar(5).gt.fpar(3)) then
         if (myproc.eq.0) then
            print *, 'NOTE that OP is using CG as the solver.'
            print *, 'The matrix (A - delta I) must be SPD!'
         endif
         call mpi_abort(mpicom, ipar(1))
      endif
c
      kmvtime = kmvtime + mpi_wtime() - tmptime
      flops = flops + fpar(11)
      if (ipar(16).lt.-ipar(7)/2) then
         if (verbose.gt.1 .and. myproc.eq.0) then
            print *, 'Block-Jacobi has faild ', -ipar(16), ' out of ',
     &           ipar(7), ' times.'
            print *, 'It will be disabled.'
         endif
         kjac = -1
      else if (ipar(16).le.0) then
         kjac = min(kjac+kjac, 100)
      endif
      return
c     end subroutine op
      end
