c
c @(#)Converted from store.f 3.3 (BNP) 3/16/89
c
c     This store function stores P vectors in memory and stores Q
c     vectors into a seriours of direct access unformated temporary
c     files.  Each process will manage its own temporary file which
c     is located in the current working directory.  The user should
c     run this program in a directory with large available space.
c
c     This program uses Fortran IO unit 99, user should
c     avoid using this unit number!
c
      subroutine store(n, isw, j, s)
      implicit none
      include 'mpif.h'
      include 'simphb.h'
c     .. Scalar Arguments ..
      integer            isw, j, n
c     ..
c     .. Array Arguments ..
      real*8             s(n)
c     ..
c     .. Parameters ..
      integer            maxll, iou
      parameter          (maxll = 2, iou = 99)
      integer            storq, retrq, storp, retrp
      parameter          (storq = 1, retrq = 2, storp = 3, retrp = 4)
      character*8        fprefix
      parameter          (fprefix = '_LANVEC_')
c     ..
c     .. Local Scalars ..
      logical            first
      character*12       filename
      integer            i, ierr, my_id
      real*8             tmp
c     ..
c
c     a common block to store IO time
c
      real*8 iotime, iosize
      common /dskstore/iotime, iosize
c     .. External Subroutines ..
      external           dcopy, dmalloc
c     ..
c     .. Save statement ..
      save               first, my_id
c     ..
c     .. Data statements ..
      data               first / .true. /
c     ..
c     .. Executable Statements ..
      tmp = MPI_Wtime()
c
c     First call to this routine
c
      if (first) then
         first = .false.
         call dmalloc(pqq, maxll*n, 8)
         call mpi_comm_rank(mpi_comm_world, my_id, ierr)
         write(filename, FMT = '(A8, I4.4)')fprefix, my_id
         open(iou, file = filename, access = 'DIRECT',
     &       form = 'UNFORMATTED', recl = 262144, iostat = ierr)
         if (ierr.ne.0) then
            print *, 'DSKSTORE: PE', my_id, 'can not open file ',
     &         filename, ' on I/O unit', iou
            call mpi_abort(mpi_comm_world, ierr)
         endif
      endif
c
c     regular operations
c
      if (isw.eq.storq) then
         write(iou, rec = j, iostat = ierr)(s(i), i = 1, n)
         if (ierr.ne.0) then
            print *, 'DSKSTORE: PE', my_id,
     &         'failed to store Lanczos vector', j
            call mpi_abort(mpi_comm_world, ierr)
         endif
      else if (isw.eq.retrq) then
         read (iou, rec = j, iostat = ierr)(s(i), i = 1, n)
         if (ierr.ne.0) then
            print *, 'DSKSTORE: PE', my_id,
     &         'failed to read Lanczos vector', j
            call mpi_abort(mpi_comm_world, ierr)
         endif
      else if (isw.eq.storp) then
         if (j.gt.maxll)stop 'STORE: (STORP) J.GT.MAXLL'
         call dcopy(n, s, 1, qq((j-1)*n+1), 1)
      else if (isw.eq.retrp) then
         if (j.gt.maxll)stop 'STORE: (RETRP) J.GT.MAXLL'
         call dcopy(n, qq((j-1)*n+1), 1, s, 1)
      endif
      iosize = iosize + 8*n
      iotime = iotime + MPI_wtime() - tmp
      return
c     end subroutine dskstore
      end
