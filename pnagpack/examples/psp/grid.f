c filename: grid3d.f
c
c     An example of using PSPARSLIB to find extreme eigenvalues of
c     3-D Lapalacian operator.
c
      program grid3d
      implicit none
      include 'simphb.h'
      include 'mpif.h'
c     .. Parameters ..
      integer            memsize
      parameter          (memsize=2000000)
      double precision   one
      parameter          (one = 1.0D0)
c     ..
c     .. Local Scalars ..
      character*128      dimstr
      integer            i, ierr, j, maxlan, maxnloc, maxnzloc,
     &                   mpicom, mpx, msglvl, neig, ngridx, nloc, nw
      double precision   condm, endl, endr, kappa, lanso_time,
     &                   setup_time
c     ..
c     .. Local Arrays ..
      integer            riord(*), split(MAX_PROCS+1)
      double precision   bnd(*), ritz(*), wrk(*)
      pointer (pwrk, wrk), (pritz, ritz), (pbnd, bnd), (priord, riord)
c     ..
c     .. External Functions ..
      integer            iargc, ncells
      external           iargc
c     ..
c     .. External Subroutines ..
      external           dmalloc, dmfree, g3dmat, plandr2,
     &                   psp_setup
c     ..
      NCELLS(Nloc,i) = 5*Nloc+4*i+1+MAX(Nloc,i+1)
c     .. Executable Statements ..
      msglvl = 5
      amxtype = 30
      ptrn = 1
      kappa = 1.0D-8
c
      paloc = 0
      pjaloc = 0
      pialoc = 0
      pix = 0
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
      call mpi_init(ierr)
      if (ierr.ne.mpi_success)stop 'FAILED TO INIT MPI'
      mpicom = mpi_comm_world
      call mpi_comm_size(mpi_comm_world, totprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, myproc, ierr)
c
      if (iargc().gt.0) then
         call getarg(1, dimstr)
cCRAY         call pxfgetarg(1, dimstr, i, ierr)
         read(dimstr, *, iostat=ierr) ngridx
         if (ierr.ne.0)call mpi_abort(mpicom, ierr)
      else
         ngridx = 32
      endif
c
      mpx = nint(exp(log(dble(totprocs))/3.0D0))
      if (mpx*mpx*mpx.ne.totprocs) then
         if (myproc.eq.0)print *,
     &       'GRID3D: number of processors must be cubic.'
         call mpi_abort(mpicom, ierr)
      endif
      setup_time = mpi_wtime()
c
c     generate the local portion of the matrix and setup the data
c     structure required for matrix-vector multiplication
c
      call g3dmat(ngridx, ngridx, ngridx, mpx, mpx, mpx, myproc,
     &            maxnloc, maxnzloc, nloc, paloc, pjaloc, pialoc, split)
      call psp_setup(maxnloc, maxnzloc, nloc, myproc, totprocs, split,
     &               aloc, jaloc, ialoc, priord, nbnd, nproc, proc, pix,
     &               ipr, mpicom, ierr)
      if (ierr.ne.0) then
         print *, 'PSP_SETUP failed with error code ', ierr
         call mpi_abort(mpicom, ierr)
      endif
c
c     the workspace for LANSO is allocated to fill up to
c     the maximum memory size available
c
      i = 2*maxnzloc + 4*maxnloc + 8*nloc
      if (totprocs.gt.1 .and. pix.ne.0) then
         i = i + 3*(nloc+max(ipr(nproc+1),ipr(nproc+2)))
      endif
      endl = dble(memsize-i) / dble(nloc+5)
      if (mpi_real8.eq.0) then
         call mpi_allreduce(endl, endr, 1, mpi_double_precision,
     &                      mpi_min, mpicom, ierr)
      else
         call mpi_allreduce(endl, endr, 1, mpi_real8, mpi_min, mpicom,
     &                      ierr)
      endif
      if (ierr.ne.mpi_success) then
         print *, 'GRID: failed to compute MAXLAN.'
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
c
c     the initial guess for the eigenvector is always [1,1,...,1]^T
c
      do 10 i = 1, nloc
         wrk(i) = one
 10   continue
      if (myproc.eq.0 .and. msglvl.gt.0) then
         write(*, FMT = '(/A$)')'GRID: '
         print *, 'NPROCS = ', totprocs, '      MAXLAN = ', maxlan
         if (msglvl.gt.5) then
            print *, 'Matrix distribution:', split(1),
     &         ('--', split(i), i = 2, totprocs+1)
         endif
      endif
      lanso_time = mpi_wtime()
      setup_time = lanso_time - setup_time
c
c     main loop to compute the smallest eigenvalue
c
      j = 0
c$$$      call plandr2(nloc, maxlan, 5, 1, condm, kappa, j, neig, ritz, bnd,
c$$$     &             wrk, nw, ierr, msglvl, mpicom)
      call pLANDR(Nloc,MAXLAN,MAXLAN,CONDM,ENDL,ENDR,0,KAPPA,
     &           J,NEIG,RITZ,BND,WRK,NW,IERR,MSGLVL,MPICOM)
      if (ierr.ne.0) then
         print *, 'LANDR returned with error code = ', ierr
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
c
      condm = max(abs(ritz(1)), abs(ritz(j)))
      if (myproc.eq.0) then
         setup_time = endl / totprocs
         lanso_time = endr / totprocs
         print *, 'STEPS= ', j, '    NEIG = ',
     &      neig, '    ||A||=', condm
         print *, 'Setup time = ', setup_time, 'sec'
         print *, 'LANSO time = ', lanso_time, 'sec'
         print *, 'MATVEC = ', kmvcnt, mmvcnt, '  time ', kmvtime,
     &      mmvtime
         print *, 'N reorth = ', nreorth, '  time ', reorthtime
         if (msglvl.gt.3) then
            condm = kappa*condm
            write(*, FMT = 9999)
            write(*, FMT = 9989)1, ritz(1), bnd(1)
            do 20 i = 2, j-1
               if (bnd(i).le.condm)
     &              write(*, FMT = 9989)i, ritz(i), bnd(i)
 20         continue
            write(*, FMT = 9989)j, ritz(j), bnd(j)
         endif
      endif
 9999 format(1X, '   INDEX       Lambda(LANSO)       ErrorBound')
 9989 format(1X, I8, G25.17, 2X, 1P, E10.2)
c
      call mpi_finalize(ierr)
c
      end
c-----------------------------------------------------------------------
      subroutine opm(n, x, y, mpicom)
      implicit none
      include 'mpif.h'
      include 'simphb.h'
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
      include 'simphb.h'
c     .. Scalar Arguments ..
      integer            mpicom, n
c     ..
c     .. Array Arguments ..
      double precision   p(n), q(n), r(n)
c     ..
c     .. Local Scalars ..
      integer            length, lentmp
      double precision   tmptime
c     ..
c     .. Local Arrays ..
      double precision   xtmp(*), ytmp(*)
      pointer (pxtmp, xtmp), (pytmp, ytmp)
c     ..
c     .. External Subroutines ..
      external           amux, amxdis, dcopy, dmalloc, dmfree
c     ..
c     .. Save statement ..
      save               pxtmp, pytmp, length
c     ..
c     .. Data statements ..
      data               pxtmp, pytmp, length / 3*0 /
c     ..
c     .. Executable Statements ..
c
      kmvcnt = kmvcnt + 1
      tmptime = mpi_wtime()
      if (totprocs.gt.1) then
c
c     allocte enough temporary workspace before calling amxdis
c
         lentmp = n + max(ipr(nproc+1), ipr(nproc+2))
         if (length.lt.lentmp) then
            if (length.gt.0) then
               call dmfree(pxtmp)
               call dmfree(pytmp)
            endif
            call dmalloc(pxtmp, lentmp, 8)
            call dmalloc(pytmp, lentmp, 8)
            length = lentmp
         endif
         call dcopy(n, q, 1, xtmp, 1)
         call amxdis(n, nbnd, xtmp, ytmp, aloc, jaloc, ialoc, nproc,
     &               proc, ix, ipr, amxtype, length, ptrn)
         call dcopy(n, ytmp, 1, r, 1)
         if (xtmp(length+1).ne.0 .or. ytmp(length+1).ne.0) then
            print *, 'Guard elements of xtmp and ytmp are:',
     &         xtmp(length+1), ytmp(length+1)
         endif
      else
         call amux(n, q, r, aloc, jaloc, ialoc)
      endif
c
      kmvtime = kmvtime + mpi_wtime() - tmptime
      return
c     end subroutine op
      end
c-----------------------------------------------------------------------
      subroutine g3dmat(nx, ny, nz, mpx, mpy, mpz, myproc, nmax, nzmax,
     &                  nloc, paval, pja, pia, split)
c     .. Scalar Arguments ..
      integer            mpx, mpy, mpz, myproc, nloc, nmax, nx, ny, nz,
     &                   nzmax
c     ..
c     .. Array Arguments ..
      integer            split(mpx*mpy*mpz+1)
      integer            ia(*), ja(*)
      double precision   aval(*)
      pointer (paval, aval), (pja, ja), (pia, ia)
c     ..
c
c     G3DMAT: creats the part of the matrix on a subcube to be used by
c     PSPARSLIB
c
c     .. Local Scalars ..
      integer            i, ioffx, ipx, ipy, ipz, irow, ix, jcol,
     &                   joffy, jy, knz, koffz, kz, nxy
      double precision   dim
c     ..
c     .. External Subroutines ..
      external           dmalloc
c     ..
c     .. Intrinsic Functions ..
      intrinsic          int
c     ..
c     .. Executable Statements ..
      if (nx.gt.1 .and. ny.gt.1 .and. nz.gt.1) then
         dim = 6
      else if ((nx.gt.1 .and. ny.gt.1) .or. (nz.gt.1 .and. ny.gt.1) .or.
     &         (nx.gt.1 .and. nz.gt.1)) then
         dim = 4
      else
         dim = 2
      endif
      nxy = nx*ny
      nloc = nxy*nz
      ioffx = nloc - nx + 1
      joffy = nloc*mpx - (ny-1)*nx
      koffz = nloc*mpx*mpy - (nz-1)*nxy
      ipz = myproc / (mpx*mpy)
      ipy = (myproc-ipz*mpx*mpy) / mpx
      ipx = myproc - (ipz*mpy+ipy)*mpx
      split(1) = 1
      do 10 i = 1, mpx*mpy*mpz
         split(i+1) = split(i) + nloc
 10   continue
      nmax = nloc+nloc
      nzmax = int(dim+1)*nloc
      call dmalloc(paval, nzmax, 8)
      call dmalloc(pja, nzmax, 1)
      call dmalloc(pia, nmax+1, 1)
c
      irow = 0
      jcol = 0
      knz = 1
      ia(1) = 1
      do 40 kz = 1, nz
         do 30 jy = 1, ny
            do 20 ix = 1, nx
               jcol = split(myproc+1) + irow
               irow = irow + 1
               aval(knz) = dim
               ja(knz) = jcol
               knz = knz + 1
c     +/- x
               if (ipx.gt.0 .or. ix.gt.1) then
                  aval(knz) = -1
                  if (ix.gt.1) then
                     ja(knz) = jcol - 1
                  else
                     ja(knz) = jcol - ioffx
                  endif
                  knz = knz + 1
               endif
               if (ipx.lt.mpx-1 .or. ix.lt.nx) then
                  aval(knz) = -1
                  if (ix.lt.nx) then
                     ja(knz) = jcol + 1
                  else
                     ja(knz) = jcol + ioffx
                  endif
                  knz = knz + 1
               endif
c     +/- y
               if (ipy.gt.0 .or. jy.gt.1) then
                  aval(knz) = -1
                  if (jy.gt.1) then
                     ja(knz) = jcol - nx
                  else
                     ja(knz) = jcol - joffy
                  endif
                  knz = knz + 1
               endif
               if (ipy.lt.mpy-1 .or. jy.lt.ny) then
                  aval(knz) = -1
                  if (jy.lt.ny) then
                     ja(knz) = jcol + nx
                  else
                     ja(knz) = jcol + joffy
                  endif
                  knz = knz + 1
               endif
c     +/- z
               if (ipz.gt.0 .or. kz.gt.1) then
                  aval(knz) = -1
                  if (kz.gt.1) then
                     ja(knz) = jcol - nxy
                  else
                     ja(knz) = jcol - koffz
                  endif
                  knz = knz + 1
               endif
               if (ipz.lt.mpz-1 .or. kz.lt.nz) then
                  aval(knz) = -1
                  if (kz.lt.nz) then
                     ja(knz) = jcol + nxy
                  else
                     ja(knz) = jcol + koffz
                  endif
                  knz = knz + 1
               endif
c
               ia(irow+1) = knz
 20         continue
 30      continue
 40   continue
c
      return
c     end g3dmat
      end
