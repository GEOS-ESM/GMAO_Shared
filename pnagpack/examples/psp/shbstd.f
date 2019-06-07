      program shbstd
      implicit none
c
c     This program performs the following tasks:
c     (1) It constructs a generalized symmetric eigenvalue problem by
c     reading two matrices stored in Boeing-Harwell format.  The matrices
c     are read by one processor and distributed to all processors in
c     the group.  The second matrix is assumed to be the mass matrix
c     and only its diagonal is used.  The generalized eigenvalue problem
c     is solved in standard form.
c     OP = M^{-1} K, OPM = M.
c     The matrix K is stored in the data structure named A
c     (Aloc, ...), the matrix M is stored as a vector call BLOC,
c     its inverse is stored in an array called BINV.
c     (2) It calls PLANDR once using MEMSIZE memory as workspace to
c     compute as many eigenvalues as possible.  Eigenvectors are
c     not computed, Lanczos vectors are stored only to perform
c     reorthogonalization.  The Lanczos vectors are stored in common
c     block GETPUT defined in simphb.h
c     The matrix-vector multiplication routine from PSPARSLIB is used.
c     Dynamic memory allocation is used through-out the preprocessing
c     stages.
c
      include 'mpif.h'
c
c     the arrays for the distrubited MATVEC (amxdis) is stored in
c     the following header file
c
      include 'simphb.h'
c
c     the size of machine's main memory in real*8 words
c     20Megawords will be allocated to the workspace of this program
c
      integer memsize
      parameter (memsize=2000000)
c
c     locally declared variables
c
      character*128 mfile, kfile
      integer iargc,mpicom,ierr,iou,msglvl
      integer i,ev,nw,nloc,j,neig,ncells,maxlan
      integer maxnzloc,maxnloc,ord(*), split(MAX_PROCS+1)
      real*8 condm,endl,endr,kappa,setup_time,lanso_time,one
      real*8 wrk(*),ritz(*),bnd(*)
      pointer (pwrk, wrk), (pord,ord), (pritz, ritz), (pbnd, bnd)
      parameter(one=1.0D0)
      NCELLS(Nloc,i) = 5*Nloc+4*i+1+MAX(Nloc,i+1)
C
C.... A few constants
C
      ev = 0
      iou = 17
      msglvl = 11
      amxtype = 30
      ptrn = 1
      kappa = 1.0D-6
c
      pbloc = 0
      pbinv = 0
      paloc = 0
      pjaloc = 0
      pialoc = 0
      pix = 0
      pord = 0
      pqq = 0
      pwrk = 0
      pritz = 0
      pbnd = 0
      kmvcnt = 0
      mmvcnt = 0
      nreorth = 0
      kmvtime = 0
      mmvtime = 0
      reorthtime = 0
c
c     start MPI
c
      call mpi_init(ierr)
      IF (IERR .NE. MPI_SUCCESS) STOP 'FAILED TO INIT MPI'
      MPICOM = MPI_COMM_WORLD
      call mpi_comm_size(mpicom, totprocs, ierr)
      call mpi_comm_rank(mpicom, myproc, ierr)
      if (totprocs.gt.MAX_PROCS) then
         print *, 'Increase MAX_PROCS in file shbstd.h and PSPARSLIB.'
         call mpi_abort(mpicom, ierr)
      endif
      if (myproc.eq.0) print *, 'SHBSTD'
c
c     the file name can be given as a command line argument
c
      if (iargc().gt.1) then
         call getarg(1, kfile)
         call getarg(2, mfile)
cCRAY         call pxfgetarg(1, kfile, i, ierr)
cCRAY         call pxfgetarg(2, mfile, i, ierr)
      else if (myproc.eq.0) then
         print *, 'Please input the name of the stiffness matrix:'
         read '(A)', kfile
         print *, 'Please input the name of the stiffness matrix:'
         read '(A)', mfile
      endif
c
c     read the given file and distribute the matrix
c
      setup_time = mpi_wtime()
      call hbread2(kfile, mfile, 0, myproc, totprocs, nloc, split,
     &     paloc, pjaloc, pialoc, pbloc, maxnzloc,
     &     maxnloc, iou, mpicom, ierr, msglvl)
      if (ierr.ne.0) then
         print *, 'HBREAD2 failed with error code ', ierr
         call mpi_abort(mpicom, ierr)
      endif
c
c     prepare the matrix for distributed matrix-vector multiplication
c
      call psp_setup(maxnloc, maxnzloc, nloc, myproc, totprocs, split,
     &     aloc, jaloc, ialoc, pord, nbnd, nproc, proc, pix, ipr,
     &     mpicom, ierr)
      if (ierr.ne.0) then
         print *, 'PSP_SETUP failed with error code ', ierr
         call mpi_abort(mpicom, ierr)
      endif
c
c     reorder the M-matrix after psp_setup has reordered K-matrix
c
      call dmalloc(pbinv, nloc, 8)
      if (pbinv.eq.0) call MPI_abort(mpicom,i)
      if (pord .ne. 0) then
         do 20 i = 1, nloc
            binv(i) = bloc(ord(nloc+i))
 20      continue
      else
         do 25 i = 1, nloc
            binv(i) = bloc(i)
 25      continue
      endif
      do 30 i = 1, nloc
         if (binv(i).gt.0.0D0 .or. binv(i).lt.0.0D0) then
            bloc(i) = binv(i)
            binv(i) = 1.0D0/binv(i)
         else
            bloc(i) = 0.0D0
            binv(i) = 0.0D0
         endif
 30   continue
c
c     the workspace for LANSO is allocated to fill up to
c     the maximum memory size available
c
      endl = max(2.0D1, dble(memsize - 8*nloc)/dble(nloc+5))
      if (mpi_real8.eq.0) then
         call MPI_ALLREDUCE(endl, endr, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MIN, mpicom, ierr)
      else
         call MPI_ALLREDUCE(endl, endr, 1, MPI_REAL8, MPI_MIN,
     &        mpicom, ierr)
      endif
      if (ierr.ne.MPI_SUCCESS) then
         print *, 'SHBSTD: failed to compute MAXLAN.'
         call mpi_abort(mpicom, ierr)
      endif
      maxlan = min(int(endr), split(totprocs+1)-1)
      nw = ncells(nloc, maxlan)
      call dmalloc(pqq, nloc*(maxlan+2), 8)
      call dmalloc(pwrk, nw, 8)
      call dmalloc(pritz, maxlan, 8)
      call dmalloc(pbnd, maxlan, 8)
      if (pqq.eq.0 .or. pwrk.eq.0 .or. pritz.eq.0 .or. pbnd.eq.0) then
         print *, 'SHBSTD: failed to allocate workspace.'
         call mpi_abort(mpicom, ierr)
      endif
      endl = -1.0D-30
      endr = 1.0D-30
      condm = one
c
c     the initial guess for the eigenvector is always [1,1,...,1]^T
c
      do 10 i = 1, nloc
         wrk(i) = one
 10   continue
      if (myproc.eq.0 .and. msglvl.gt.0) then
         write (*, '(/A$)') 'SHBSTD: '
         print *, 'NPROCS = ', totprocs, '	MAXLAN = ', maxlan
         if (msglvl.gt.5) then
            print *, 'Matrix distribution:', split(1),
     &           ('--', split(i), i=2,totprocs+1)
         endif
      endif
      lanso_time = mpi_wtime()
      setup_time = lanso_time - setup_time
c
c     main loop to compute the smallest eigenvalue
c
      j = 0
      call pLANDR(Nloc,MAXLAN,MAXLAN,CONDM,ENDL,ENDR,EV,KAPPA,
     &           J,NEIG,RITZ,BND,WRK,NW,IERR,MSGLVL,MPICOM)
      if (ierr.ne.0) then
         print *, 'LANDR returned with error code = ', ierr
         call mpi_abort(mpi_comm_world,ierr)
      endif
c
      lanso_time = mpi_wtime() - lanso_time
      call dmfree(pwrk)
      call dmfree(pqq)
      if (mpi_real8.eq.0) then
         call mpi_reduce(setup_time, endl, 1, MPI_double_precision,
     &         MPI_SUM, 0, mpicom, ierr)
         call mpi_reduce(lanso_time, endr, 1, MPI_double_precision,
     &         MPI_SUM, 0, mpicom, ierr)
      else
         call mpi_reduce(setup_time, endl, 1, MPI_REAL8, MPI_SUM, 0,
     &        mpicom, ierr)
         call mpi_reduce(lanso_time, endr, 1, MPI_REAL8, MPI_SUM, 0,
     &        mpicom, ierr)
      endif
c
      condm = kappa*max(abs(ritz(1)), abs(ritz(j)))
      if (myproc .eq. 0) then
         setup_time = endl / totprocs
         lanso_time = endr / totprocs
         print *, 'LANCZOS STEPS = ', j
         print *, 'Eigenpairs converged = ', neig, condm
         print *, 'Setup time = ', setup_time, 'sec'
         print *, 'LANSO time = ', lanso_time, 'sec'
         print *, 'MATVEC = ', kmvcnt, mmvcnt,
     &        '  time ', kmvtime, mmvtime
         print *, 'N reorth = ', nreorth, '  time ', reorthtime
         write (*, 9000)
         write (*, 9010) 1, ritz(1), bnd(1)
         do i = 2, j-1
            if (msglvl.gt.3 .or. bnd(i).le.condm)
     &           write (*, 9010) i, ritz(i), bnd(i)
         enddo
         write (*, 9010) j, ritz(j), bnd(j)
      endif
 9000 format(1X, '   INDEX       Lambda(LANSO)       ErrorBound')
 9010 format(1X, I8, G25.17, 2X,1PE10.2)
c
      call mpi_finalize(ierr)
c$$$      stop 'DONE'
      end
c-----------------------------------------------------------------------
      subroutine opm(n,x,y,mpicom)
      implicit none
      include 'mpif.h'
      include 'simphb.h'
      integer i, n,mpicom
      real*8 x(n), y(n), tmptime
      mmvcnt = mmvcnt + 1
      tmptime = MPI_Wtime()
      do 10 i = 1, n
         y(i) = bloc(i) * x(i)
 10   continue
      mmvtime = mmvtime + MPI_Wtime() - tmptime
      return
c     end subroutine opm
      end
c-----------------------------------------------------------------------
      subroutine op(n,p,q,r,mpicom)
      implicit none
      integer i, n, length, lentmp, mpicom
      real*8 p(n), q(n), r(n), xtmp(*), ytmp(*), tmptime
      pointer (pxtmp, xtmp), (pytmp, ytmp)
      include 'mpif.h'
      include 'simphb.h'
      save pxtmp, pytmp, length
      data pxtmp, pytmp, length/3*0/
c
      kmvcnt = kmvcnt + 1
      tmptime = MPI_Wtime()
      if (totprocs.gt.1) then
c
c     allocte enough temporary workspace before calling amxdis
c
         lentmp = n+max(ipr(nproc+1), ipr(nproc+2))
         if (length.lt.lentmp) then
            if (length.gt.0) then
               call dmfree(pxtmp)
               call dmfree(pytmp)
            endif
            call dmalloc(pxtmp,lentmp,8)
            call dmalloc(pytmp,lentmp,8)
            if (pxtmp.eq.0 .or. pytmp.eq.0) call MPI_abort(mpicom,i)
            length = lentmp
         endif
         call dcopy(n,q,1,xtmp,1)
         call amxdis(n,nbnd,xtmp,ytmp,aloc,jaloc, ialoc, nproc, proc,
     &     ix, ipr, amxtype, length, ptrn)
CDIR$UNROLL 16
         do 10 i = 1, n
            r(i) = binv(i) * ytmp(i)
 10      continue
         if (xtmp(length+1).ne.0 .or. ytmp(length+1).ne.0) then
            print *, 'Guard elements of xtmp and ytmp are:',
     &           xtmp(length+1), ytmp(length+1)
         endif
      else
         call amux(n,q,r,aloc,jaloc,ialoc)
CDIR$UNROLL 16
         do 20 i = 1, n
            r(i) = binv(i) * r(i)
 20      continue
      endif
c
      kmvtime = kmvtime + MPI_Wtime() - tmptime
      return
c     end subroutine op
      end
