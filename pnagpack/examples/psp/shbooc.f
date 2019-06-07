      program simphb
      implicit none
c
c     This program performs the following tasks:
c     (1) It constructs a simple symmetric eigenvalue problem by
c     reading a matrix stored in Boeing-Harwell format.  The matrix
c     is read by one processor and distributed to all processors in
c     the group.  If the matirx is RSA type, it is expanded to contain
c     both upper and lower triangular parts before distributed.  If the
c     matrix is RUA, its symmetric part is used.
c     (2) It computes a few smallest eigenvalues and corresponding
c     eigenvectors by calling PLANDR2.
c     The matrix-vector multiplication routine from PSPARSLIB is used.
c     Dynamic memory allocation is used through-out the preprocessing
c     stages.
c     The disk I/O is performed by dskstore.f.
c
      include 'mpif.h'
c
c     the arrays for the distrubited MATVEC (amxdis) is stored in
c     the following header file
c
      include 'simphb.h'
c
c     locally declared variables
c
      character*128 hbfile
      integer iargc,mpicom,ierr,iou,msglvl,ned,lohi
      integer i,nw,nloc,j,neig,ncells,maxlan,maxprs
      integer maxnzloc,maxnloc,ord(*), split(MAX_PROCS+1)
      real*8 condm,endl,endr,kappa,setup_time,lanso_time
      real*8 ddotmpi,one
      real*8 wrk(*),ritz(*),bnd(*)
      pointer (pwrk, wrk), (pord,ord), (pritz, ritz), (pbnd, bnd)
      parameter(one=1.0D0)
c
c     a common block to store IO time
c
      real*8 iotime, iosize
      common /dskstore/iotime, iosize
c
      NCELLS(Nloc,i) = 5*Nloc+4*i+1+MAX(Nloc,i+1)
C
C.... A few constants
C
      iou = 17
      ned = 2
      lohi = -1
      msglvl = 1
      amxtype = 29
      ptrn = 1
      kappa = 1.0D-8
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
      iotime = 0
      iosize = 0
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
      if (msglvl.ge.0) then
         condm = ddotmpi(ialoc(nloc+1)-ialoc(1), aloc, 1, aloc, 1,
     &        mpicom)
         if (myproc.eq.0) print *, 'Matrix Frobenius norm = ',
     &        sqrt(condm)
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
c     the workspace for LANSO is allocated here
c
      maxlan = split(totprocs+1) - split(1)
      maxprs = min(ned, lmtprs)
      nw = ncells(nloc, maxlan)
      call dmalloc(pwrk, nw, 8)
      call dmalloc(pritz, maxlan, 8)
      call dmalloc(pbnd, maxlan, 8)
      condm = one
c
c     the initial guess for the eigenvector is always [1,1,...,1]^T
c
      do 10 i = 1, nloc
         wrk(i) = one
 10   continue
      if (myproc.eq.0 .and. msglvl.gt.0) then
         write (*, '(/A$)') 'SHBOOC: '
         print *, 'NPROCS = ', totprocs, '	MAXLAN = ', maxlan,
     &        '	LOHI =', lohi, '	Kappa =', kappa
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
      call PLANDR2(Nloc,MAXLAN,MAXPRS,LOHI,CONDM,KAPPA,
     &           J,NEIG,RITZ,BND,WRK,NW,IERR,MSGLVL,MPICOM)
      if (ierr.ne.0) then
         print *, 'LANDR returned with error code = ', ierr
         call mpi_abort(mpi_comm_world,ierr)
      endif
c
      lanso_time = mpi_wtime() - lanso_time
      call dmfree(pwrk)
      call dmfree(pqq)
      if (MPI_REAL8.eq.0) then
         call mpi_allreduce(setup_time, endl, 1, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, mpicom, ierr)
         call mpi_allreduce(lanso_time, endr, 1, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, mpicom, ierr)
         setup_time = endl / totprocs
         lanso_time = endr / totprocs
         call mpi_allreduce(iotime, endl, 1, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, mpicom, ierr)
         iotime = endl / totprocs
      else
         call mpi_allreduce(setup_time, endl, 1, MPI_REAL8, MPI_SUM,
     &        mpicom, ierr)
         call mpi_allreduce(lanso_time, endr, 1, MPI_REAL8, MPI_SUM,
     &        mpicom, ierr)
         setup_time = endl / totprocs
         lanso_time = endr / totprocs
         call mpi_allreduce(iotime, endl, 1, MPI_REAL8, MPI_SUM,
     &        mpicom, ierr)
         iotime = endl / totprocs
      endif
c
      if (myproc .eq. 0) then
         print *, 'LANSO steps = ', j
         print *, 'Eigenpairs: wanted', ned, '  computed', neig
         print *, 'Setup time = ', setup_time, 'sec'
         print *, 'LANSO time = ', lanso_time, 'sec'
         print *, 'Eigenairs wanted =', maxprs, ' computed = ', neig
         print *, 'MATVEC = ', kmvcnt, mmvcnt,
     &        '  time ', kmvtime, mmvtime
         print *, 'N reorth = ', nreorth, '  time ', reorthtime
         print *, 'Time used in IO = ', iotime, iosize, iosize/iotime
         write (*, 9000)
         do i = 1, min(maxprs+4, j)
            write (*, 9010) ritz(i), bnd(i)
         enddo
      endif
 9000 format(1X, '       Lambda(LANSO)       ',
     &     'ErrorBound')
 9010 format(1X, G25.17, 2X, 1PE10.2)
c
      call dmfree(pritz)
      call dmfree(pbnd)
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
