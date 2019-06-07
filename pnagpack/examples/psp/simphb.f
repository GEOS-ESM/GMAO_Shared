      program simphb
      implicit none
c
c     This program performs the following tasks:
c     (1) It constructs a simple symmetric eigenvalue problem by
c     reading a matrix stored in Boeing-Harwell format.  The matrix
c     is read by one processor and distributed to all processors in
c     the group.  If the matrix is RSA type, it is expanded to full
c     matrix and it is RUA type, its symmetric part will be used.
c     (2) It calls LANDR once using MEMSIZE memory as workspace to
c     compute as many eigenvalues as possible.  Eigenvectors are
c     not computed, Lanczos vectors are stored only to perform
c     reorthogonalization.  The Lanczos vectors are stored in the common
c     block GETPUT defined in file simphb.h
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
      character*128 hbfile
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
      msglvl = 1
      amxtype = 29
      ptrn = 1
      kappa = 1.0D-6
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
         print *, 'Increase MAX_PROCS in file simphb.h and PSPARSLIB.'
         call mpi_abort(mpicom, ierr)
      endif
      if (myproc.eq.0) print *, 'SIMPHB'
c
c     the file name can be given as a command line argument
c
      if (iargc().gt.0) then
         call getarg(1, hbfile)
cNON-CRAY         call pxfgetarg(1, hbfile, i, ierr)
      else if (myproc.eq.0) then
         print *, 'Please input the Boeing-Harwell file name:'
         read *, hbfile
      endif
c
c     read the given file and distribute the matrix
c
      setup_time = mpi_wtime()
      call hbread1(hbfile, 0, myproc, totprocs, nloc, split,
     &     paloc, pjaloc, pialoc, maxnzloc,
     &     maxnloc, iou, mpicom, ierr, msglvl)
      if (ierr.ne.0) then
         print *, 'HBREAD1 failed with error code ', ierr
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
         print *, 'SIMPHB: failed to compute MAXLAN.'
         call mpi_abort(mpicom, ierr)
      endif
      maxlan = min(int(endr), split(totprocs+1)-1)
      nw = ncells(nloc, maxlan)
      call dmalloc(pqq, nloc*(maxlan+2), 8)
      call dmalloc(pwrk, nw, 8)
      call dmalloc(pritz, maxlan, 8)
      call dmalloc(pbnd, maxlan, 8)
      if (pqq.eq.0 .or. pwrk.eq.0 .or. pritz.eq.0 .or. pbnd.eq.0) then
         print *, 'SIMPHB: failed to allocate workspace.'
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
         write (*, '(/A$)') 'SIMPHB: '
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
      integer i, n,mpicom
      real*8 x(n), y(n), tmptime
      include 'mpif.h'
      include 'simphb.h'
      mmvcnt = mmvcnt + 1
      tmptime = MPI_Wtime()
      do 10 i = 1, n
         y(i) = x(i)
 10   continue
      mmvtime = mmvtime + MPI_Wtime() - tmptime
      return
c     end subroutine opm
      end
c-----------------------------------------------------------------------
      subroutine op(n,p,q,r,mpicom)
      implicit none
      integer n, length, lentmp, mpicom
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
         call dcopy(n,ytmp,1,r,1)
         if (xtmp(length+1).ne.0 .or. ytmp(length+1).ne.0) then
            print *, 'Guard elements of xtmp and ytmp are:',
     &           xtmp(length+1), ytmp(length+1)
         endif
      else
         call amux(n,q,r,aloc,jaloc,ialoc)
      endif
c
      kmvtime = kmvtime + MPI_Wtime() - tmptime
      return
c     end subroutine op
      end
