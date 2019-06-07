      subroutine psp_setup(nmax, nzmax, nloc, my_id, npe, split,
     &     aval, ja, ia, priord, nbnd, nproc, proc, pix, ipr,
     &     mpicom, ierr)
      implicit none
      integer nloc, nbnd, nproc, npe, mpicom, ierr, nmax, nzmax
      integer split(npe+1), ja(nzmax), ia(nmax+1),
     &     ipr(npe+1), proc(npe+npe), riord(*), ix(*)
      real*8 aval(nzmax)
      pointer (priord, riord), (pix, ix)
      include 'mpif.h'
c
c     This is a driver for BDRY and SETUP1 of the PSPARSLIB.
c     It allocates workspace as needed. The output is just enough
c     run the matrix-vector multiplication routine from PSARPSLIB.
c
c NOTE ON MATRIX DISTRIBUTION:
c     the matrix is assumed to have been reordered already before
c     entering into this routine.  each processor in the group
c     is to have a consecutive block of rows of the original
c     matrix where the global row indices are stored in array
c     'split'.  the column indices in ja are still the global
c     indices -- this setup routine will change them to local
c     ordering and produce ordering in 'riord'.
c
c     INPUT:
c     nmax   -- the size of the array ia
c     nloc   -- number of rows in the current processor
c     split  -- the (global) starting row number on each processor
c     aval,ja,ia -- the submatrix (block row) in CSR format
c     mpicom -- MPI communicator
c
c     OUTPUT:
c     aval,ja,ia -- reordered with interior nodes first
c     riord  -- (through priord) integer array of size nloc+nloc
c               where riord(i) is global row number of ith local
c               row before permutation, riord(nloc+i) is the ith
c               row's index when entering this routine. After
c               returning from this routine, the ith row's global
c               row index is riord(riord(nloc+i)).
c               If there is only one processor, priord=0.
c     nbnd   -- beginning of the interface nodes in ia
c     nproc  -- number of neighboring processors
c     proc   -- list of neighboring processors (NPE*2)
c     ix     -- (through pix) local interface points found (NLOC)
c     ipr    -- ipr(j) points to the beginning of the list of interface
c               points that are coupled with processor proc(j)
c     ierr   -- error flag
c
      real*8 cval(*), bval(*)
      integer my_id, i, j, iou, nzloc, stype, iprow, ntot, nix
      integer jamap(*), iwk(*), link(*),
     &     ic(*), jc(*), ib(*), jb(*)
      pointer (pjamap, jamap), (piwk, iwk), (plink, link),
     &     (pbval, bval), (pib, ib), (pjb, jb),
     &     (pcval, cval), (pic, ic), (pjc, jc)
      data pjamap, piwk, plink, pcval, pic, pjc, pib, pjb,
     &     pbval /9*0/
c
      if (npe.eq.1) then
         nbnd = nloc+1
         nproc = 0
         proc(1) = 0
         ipr(1) = 1
         ipr(2) = 1
         priord = 0
         ierr = 0
         return
      endif
c
      iou = 6
      stype = 71
      nzloc = ia(nloc+1)-ia(1)
      if (nzloc.gt.nzmax) then
         ierr = -1
         return
      else
         ierr = 0
      endif
c     *** NOTE ***
c     mpicom can only be MPI_COMM_WORLD because PSPARSLIB only uses
c     MPI_COMM_WORLD communicator!
      call mpi_comm_compare(MPI_COMM_WORLD, mpicom, i, ierr)
      if (ierr .ne. mpi_success) then
         print *, 'PSP_SETUP failed MPI_COMM_COMPARE.'
         call MPI_ABORT(MPI_COMM_WORLD, ierr)
      else if (i .eq. MPI_UNEQUAL) then
         print *, 'PSPARSLIB was designed to use MPI_COMM_WORLD.'
         call MPI_ABORT(MPI_COMM_WORLD, ierr)
      endif
c
      call MPI_ALLREDUCE(nloc, ntot, 1, MPI_INTEGER, MPI_SUM, mpicom,
     &     ierr)
      if (ierr.ne.MPI_SUCCESS) then
         print *, 'MPI_ALLREDUCE failed to compute NTOT on processor ',
     &        my_id
         call MPI_ABORT(MPI_COMM_WORLD, ierr)
      endif
c
      nix = max(max(nloc, npe), min(nzloc, (npe-1)*nloc))
      call dmalloc(priord, 2*nloc, 1)
      call dmalloc(pix, nix, 1)
      call dmalloc(piwk, ntot, 1)
      call dmalloc(pjamap, nzloc, 1)
      call dmalloc(plink, ntot, 1)
      if (priord.eq.0 .or. pix.eq.0 .or. piwk.eq.0 .or. pjamap.eq.0
     &     .or. plink.eq.0) then
         call MPI_abort(mpicom, ierr)
      endif
c
      j =  split(my_id+1) - 1
      do 10 i = 1, nloc
         riord(i) = j + i
 10   continue
      jamap(1) = iprow(npe, split, ja(1), 1)
      do 30 i = 2, nzloc
         jamap(i) = iprow(npe, split, ja(i), jamap(i-1))
 30   continue
c$$$      print *, 'PE', my_id, ' entering into BDRY..'
c$$$      write (10+my_id, *) 'RIORD:', (riord(i), i=1, nloc)
c$$$      write (10+my_id, '(1X,2I10)') (ja(i), jamap(i), i=1,nzloc)
      call bdry(nloc, npe, my_id, ia, jamap, riord, iwk,
     &     link, nbnd, nproc, ix, ipr, proc, iou)
c
c     check to make sure enough workspace was allocated to IX
c
      i = ipr(nproc+1)-ipr(1)
      if (nix .lt. i) then
         print *, 'PE', my_id, ': array IX(',nix,
     &        ') smaller than needed', i
         call MPI_ABORT(MPI_COMM_WORLD, ierr)
      endif
c
c     sanity checks on the dynamic arrays
c
      if (riord(nloc+nloc+1).ne.0 .or. iwk(ntot+1).ne.0 .or.
     &     jamap(nzloc+1).ne.0 .or. link(ntot+1).ne.0) then
         print *, 'PE', my_id, ': Guard elemements of riord, iwk,',
     &        ' jamap and link: ', riord(nloc+nloc+1), iwk(ntot+1),
     &        jamap(nzloc+1), link(ntot+1)
      endif
      call dmfree(pjamap)
      call dmfree(plink)
      my_id = my_id - 1
c$$$      write (10+my_id, *) 'BDRY: Nloc =', nloc, '   Nbnd =', nbnd,
c$$$     &     '   Nproc =', nproc, '   NNZloc =', ia(nloc+1)-ia(1)
c$$$      write (10+my_id, *) 'RIORD:', (riord(i), i=1, nloc)
c$$$      write (10+my_id, *) 'IPR ', (ipr(i), i=1, nproc+1)
c$$$      write (10+my_id, *) 'PROC', (proc(i), i=1,nproc)
c$$$      write (10+my_id, *) 'IX  ', (ix(i), i=1,ipr(nproc+1)-1)
c$$$      call dump(1,nloc,.true.,aval,ja,ia,10+my_id)
c
      call dmalloc(pbval, max(nzloc, ntot), 8)
      call dmalloc(pcval, nzloc, 8)
      call dmalloc(pjb, max(nzloc, ntot), 1)
      call dmalloc(pjc, nzloc, 1)
      call dmalloc(pib, ntot+1, 1)
      call dmalloc(pic, nloc+1, 1)
      if (pbval.eq.0 .or. pcval.eq.0 .or. pjb.eq.0 .or. pjc.eq.0
     &     .or. pib.eq.0 .or. pic.eq.0) then
         call MPI_ABORT(MPI_COMM_WORLD, ierr)
      endif
      do 40 i = 1, nloc+1
         ic(i) = ia(i)
 40   continue
      do 50 i = 1, nzloc
         jc(i) = ja(i)
         cval(i) = aval(i)
 50   continue
c$$$      print *, 'PE', my_id, ' entering into SETUP1..'
      call setup1(nloc, nbnd, cval, jc, ic, nproc, proc, ix, ipr,
     &     aval, ja, ia, bval, jb, ib, iwk, riord, stype, nzloc,
     &     j, iou, nzloc, nmax+1)
      i = max(nzloc, ntot) + 1
      if (bval(i).ne.0 .or. jb(i).ne.0 .or. cval(nzloc+1).ne.0 .or.
     &     jc(nzloc+1).ne.0) then
         print *, 'PE', my_id, ': Guard elemements of jb, bval,',
     &        ' jc and cval: ', jb(i), bval(i),
     &        jc(nzloc+1), cval(nzloc+1)
      endif
      if (ib(ntot+2).ne.0 .or. ic(nloc+2).ne.0 .or.
     &     riord(nloc+nloc+1).ne.0 .or. iwk(ntot+1).ne.0) then
         print *, 'PE', my_id, ': Guard elemements of riord, iwk,',
     &        ' ib and ic: ', riord(nloc+nloc+1), iwk(ntot+1),
     &        ib(ntot+2), ic(nloc+2)
      endif
c$$$      write (10+my_id, *) 'SETUP1: interface nodes =', nloc-nbnd+1,
c$$$     &     ' Number of nonzero columns =', j
c$$$      write (10+my_id, *) (('(',i,riord(nloc+i),')'), i=1,nloc)
c$$$      call dump(1,nloc+nloc-nbnd+1,.true.,aval,ja,ia,10+my_id)
c$$$      call flush(10+my_id)
      call dmfree(pic)
      call dmfree(pib)
      call dmfree(pjc)
      call dmfree(pjb)
      call dmfree(pcval)
      call dmfree(pbval)
      call dmfree(piwk)
c
      j = nloc - nbnd + 1
      if (j.gt.nmax) then
         print *, 'PSP_SETUP (PE',my_id,') No. of boundry nodes=', j,
     &        ' nloc=', nloc, ' nmax=', nmax
         ierr = -2
      endif
c
      j = 0
      ipr(nproc+2) = 0
      do 60 i = 1, nproc
         if (ipr(i).gt.ipr(i+1) .or. proc(i).lt.0 .or.
     &        proc(nproc+i).lt.0) j = j + 1
         ipr(nproc+2) = ipr(nproc+2) + proc(nproc+i)
 60   continue
      if (j.gt.0) then
         print *, 'PSP_SETUP (PE',my_id,') array IPR is corrupt.'
         print *, (ipr(i), i=1,nproc+2)
      endif
c
c     realloc array IX to reduce its size
c
      nix = ipr(nproc+1)-1
      call dmalloc(piwk, nix, 1)
      if (piwk.eq.0) call MPI_Abort(MPI_COMM_WORLD, ierr)
      do 70 i = 1, nix
         iwk(i) = ix(i)
 70   continue
      call dmfree(pix)
      call dmalloc(pix, nix, 1)
      if (pix.eq.0) call MPI_ABORT(MPI_COMM_WORLD, ierr)
      do 80 i = 1, nix
         ix(i) = iwk(i)
 80   continue
      call dmfree(piwk)
c
c     sort the elements in each row in column order
c
      call dmalloc(piwk, nzloc+nzloc, 1)
      if (piwk.eq.0) call MPI_Abort(MPI_COMM_WORLD, ierr)
      call csort(nloc, aval, ja, ia, iwk, .true.)
      call dmfree(piwk)
      return
c     end subroutine psp_setup
      end
c-----------------------------------------------------------------------
      integer function iprow(nprocs, split, index, ip0)
      implicit none
      integer nprocs, index, ip0, split(nprocs+1)
c
c     given a row number, this function computes the processor it
c     belongs to by consulting array SPLIT.
c
      integer iup, imed
      iprow = ip0
      if (index.ge.split(ip0) .and. index.lt.split(ip0+1)) return
      if (index.lt.1 .or. index.ge.split(nprocs+1)) then
         print *, 'IPROW: ', index, ip0, split(ip0), split(nprocs+1)
         stop
      endif
c
      if (nprocs.gt.10) then
c
c     binary search
c
         iprow = 1
         iup = nprocs+1
 10      imed = (iprow+iup) / 2
         if (iup.gt.iprow+1) then
            if (index .ge. split(imed)) then
               iprow = imed
            else if (index .lt. split(imed)) then
               iup = imed
            endif
            goto 10
         endif
      else if (nprocs.gt.2) then
c
c     linear search starting from ip0
c
         if (index.ge.split(ip0+1)) then
 20         iprow = iprow + 1
            if (index.ge.split(iprow+1)) goto 20
         else
 30         iprow = iprow - 1
            if (index .lt. split(iprow)) goto 30
         endif
      else
         iprow = 3 - ip0
      endif
c
      return
c     end function iprow
      end
