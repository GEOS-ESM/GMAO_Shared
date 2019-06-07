c
c     This subroutine read a matrix in Harwell-Boeing format
c     and distribute it to all members of MPICOM group.
c     It reads the whole matrix on one processor, performs
c     reordering and distributes the reordered matrix as
c     row-blocks.
c     If the matrix is RSA type, it expands into full matrix.
c     If the matrix is RUA type, it takes the symmetric part of the
c     matrix.
c
      subroutine hbread2(kfile, mfile, mapopt, my_id, nprocs, nrow,
     &     split, paval, pja, pia, pbdiag, nzmax, nmax, iou,
     &     mpicom, ierr, msglvl)
      implicit none
      include 'mpif.h'
      integer nrow, nzmax, nmax, ierr, mapopt, msglvl, mpicom, iou
      integer nprocs, my_id
      integer ia(*), ja(*), split(nprocs+1)
      real*8  aval(*), bdiag(*)
      character kfile*(*), mfile*(*)
      pointer (paval, aval), (pbdiag, bdiag), (pja, ja), (pia, ia)
c
c     This routine reads two Harwell-Boeing matrices from files and
c     distribute the matrix as NPROCS consecutive row blocks,
c     where NPROCS is the number of processors in the MPICOM
c     group.
c
c     Only the diagonal part of the second matrix is used (bdiag).
c     The first matrix is scaled by bdiag and the zero rows are removed
c     from storage (matrix unchanged).  This setup is usful for solve
c     a generalized eigenvalue problem in standard form.
c     OP = M^{-1} K, OPM = M.
c
c     The variable 'map' should be passed back to the caller.
c
c     INPUT:
c     filename  -- name of the Boeing-Harwell file
c     my_id     -- MPI_RANK
c     nprocs    -- number of processors in current group
c     iou       -- the fortran I/O unit number to be used when read file
c     mpicom    -- MPI communicator (group)
c     msglvl    -- more message is printed if msglvl is higher
c
c     OUTPUT:
c     nrow      -- number rows distributed to this processor
c     split     -- the starting point of each row block distributed
c     aval,ja,ia-- the CSR structure to store the submatrix (row block)
c     bdiag     -- the diagonal elements of mass matrix
c     nzmax     -- size of aval and ja
c     nmax      -- size of ia
c     ierr      -- error flag
c
      integer itmp(*), jtmp(*)
      real*8  vtmp(*), diag(*)
      pointer (pitmp, itmp), (pjtmp, jtmp), (pvtmp, vtmp), (pdiag, diag)
c     .. Local Scalars ..
      character*8        key
      character*72       title
      integer            totrow, totnz, i
c
c     arrays for partitioning the matrix
c
      integer map(*)
      pointer (pmap, map)
c     ..
      data pitmp, pjtmp, pvtmp, pmap, pdiag /5*0/
c
      totrow = 1
      totnz = 1
      if (my_id.eq.0) then
c
c     read in the mass matrix first since it require less temporary
c     storage
c
         totrow = 0
         totnz = 0
         call hbresy(msglvl, iou, mfile, key, title, totrow,
     &        totnz, pvtmp, pjtmp, pitmp, ierr)
         if (ierr.ne.0) call MPI_ABORT(MPI_COMM_WORLD, ierr)
         nrow = totrow
         call dmalloc(pdiag, totrow, 8)
         call dmalloc(pmap, totrow, 1)
         call getdia(totrow, totrow, 0, vtmp, jtmp, itmp, i, diag,
     &        map, 0)
         call hbresy(msglvl, iou, kfile, key, title, totrow,
     &        totnz, pvtmp, pjtmp, pitmp, ierr)
         if (ierr.ne.0) call MPI_ABORT(MPI_COMM_WORLD, ierr)
         totnz = itmp(totrow+1) - itmp(1)
         if (nrow .ne. totrow) then
            print *,
     &           'The size of the two input matrices must be the same.'
            call MPI_ABORT(MPI_COMM_WORLD, ierr)
         endif
c
c     .. decide what rows to send to each processor
c     this is a partitioning routine which probably also reorder
c     the matrix.  the input matrix is to be reordered so that
c     each process will receive one consecutive block rows of
c     the reordered matrix
c     
         if (mapopt.eq.0) then
            call dummap(totrow, jtmp, itmp, nprocs, map, split)
         else if (mapopt.eq.-1) then
            call glimap(totrow, jtmp, itmp, nprocs, map, split)
         else
            call dummap(totrow, jtmp, itmp, nprocs, map, split)
         endif
         call permpair(totrow, vtmp, jtmp, itmp, diag, map)
         call dmfree(pmap)
      endif
c
c     .. call distribution routine
c
      call scatter_pair(totrow, totnz, pvtmp, pjtmp, pitmp, pdiag,
     &     my_id, nprocs, split, nrow, paval, pja, pia, pbdiag,
     &     nzmax, nmax, mpicom, ierr)
c$$$      call dump(1,nrow,.true.,aval,ja,ia,10+my_id)
c
      if (my_id .eq. 0 .and. totrow.gt.nrow) then
         call dmfree(pitmp)
         call dmfree(pjtmp)
         call dmfree(pvtmp)
         call dmfree(pdiag)
      endif
      return
c     end of subroutine hbread2
      end
c-----------------------------------------------------------------------
      subroutine permpair(nrow, aval, ja, ia, bdiag, map)
      implicit none
      integer nrow, ja(*), ia(nrow+1), map(nrow)
      real*8 aval(*), bdiag(nrow)
c
c     permute a matrix stored in CSR format (symmetric permutation)
c     a,ja,ia- (INPUT/OUTPUT) permuted matrix is tored in the same array
c              internally, another copy is produced.
c     bdiag  - (INPUT/OUTPUT) permuted diagonal of associated matrix
c     map    - (INPUT) the ith row in the new matrix is the map(i)th
c              row of the old one.
c
      logical isnulmap
      real*8  tmp(*)
      integer i, j, k, nnz, pntr(*), indx(*), invmap(*)
      pointer (ppntr, pntr), (pindx, indx), (pinvmap, invmap),
     $     (ptmp, tmp)
      data ppntr, pindx, pinvmap, ptmp/4*0/
c
c     don't do anything if the map is a trivial map
c
      if (isnulmap(nrow, map)) return
c
c     really have to change the matrix
c
      nnz = ia(nrow+1)-ia(1)
      call dmalloc(ppntr, nrow+1, 1)
      call dmalloc(pinvmap, nrow, 1)
      call dmalloc(pindx, nnz, 1)
c
      do 10 i = 1, nrow+1
         pntr(i) = ia(i)
 10   continue
      do 20 i = 1, nnz
         indx(i) = ja(i)
 20   continue
      do 30 i = 1, nrow
         invmap(map(i)) = i
 30   continue
c
c     compute the new ia, (ja is temporarily used as workspace)
c
      do 40 i = 1, nrow
         ja(i) = ia(i+1) - ia(i)
 40   continue
      ia(1) = 1
      do 50 i = 1, nrow
         ia(i+1) = ia(i) + ja(map(i))
 50   continue
c
c     shuffling the array ja
c
      do 70 i = 1, nrow
         k = pntr(map(i))
         do 60 j = ia(i), ia(i+1)-1
            ja(j) = invmap(indx(k))
            k = k + 1
 60      continue
 70   continue
c
      call dmfree(pinvmap)
      call dmfree(pindx)
c
c     permute the values
c
      call dmalloc(ptmp, max(nnz, nrow), 8)
      do 80 i = 1, nnz
         tmp(i) = aval(i)
 80   continue
      do 100 i = 1, nrow
         k = pntr(map(i))
         do 90 j = ia(i), ia(i+1)-1
            aval(j) = tmp(k)
            k = k + 1
 90      continue
 100  continue
c
      call dcopy(nrow, bdiag, 1, tmp, 1)
      do 110 i = 1, nrow
         bdiag(map(i)) = tmp(i)
 110  continue
c
      call dmfree(ptmp)
      call dmfree(ppntr)
      return
c     end subroutine permpair
      end
c-----------------------------------------------------------------------
      subroutine scatter_pair(ntot, nztot, pag, pjg, pig, pdg, my_id,
     &     nprocs, split, nloc, pal, pjl, pil, pdl, nzmax, nmax,
     &     mpicom, ierr)
      implicit none
      integer ntot, nztot, nprocs,nloc, nzmax, nmax, mpicom, ierr,
     &     il(*), jl(*), split(nprocs+1), jg(*), ig(*)
      real*8 ag(*), al(*), dg(*), dl(*)
      pointer (pag, ag), (pjg, jg), (pig, ig), (pdg, dg)
      pointer (pal, al), (pjl, jl), (pil, il), (pdl, dl)
      include 'mpif.h'
c
c     distribute a CSR matrix on PE 0 to everyone in the group MPICOM
c     according to the array SPLIT
c
c     INPUT:
c     ntot   -- the global matrix size
c     nztot  -- total number of nonzero elements
c     ag,jg,ig -- the whole matrix in CSR format
c     dg     -- the diagonal part of the second matrix (global)
c     my_id  -- MPI_RANK
c     nprocs -- number of processes 
c     split  -- the global index of the blocks
c     mpicom -- MPI communicator for the group
c
c     OUTPUT:
c     nloc   -- numver of rows in this processor
c     al,jl,il -- the submatrix (row block) on this processor
c               in CSR format
c     dl     -- the local part of the matrix
c     nmax   -- size of array IL
c     nzmax  -- size of array JL and AL
c     ierr   -- error flag
c
c     local variables
c
      integer i, my_id, nzloc, mpitype, nelm(*), disp(*)
      pointer (pnelm, nelm), (pdisp, disp)
      data pnelm, pdisp/2*0/
c
c     if there is only one processor, simply copy CSR
c
      if (nprocs.eq.1) then
         nloc = ntot
         nmax = ntot
         nzmax = nztot
         pal = pag
         pjl = pjg
         pil = pig
         pdl = pdg
         return
      endif
      if (MPI_REAL8.eq.0) then
         mpitype = MPI_Double_precision
      else
         mpitype = MPI_REAL8
      endif
c
c     broadcast how the rows are split
c
      call MPI_BCAST(split, nprocs+1, MPI_INTEGER, 0, mpicom, ierr)
      if (ierr .ne. MPI_SUCCESS) then
         print *, 'MPI_BCAST(split) error flag = ', ierr
         ierr = -2
         return
      endif
c
c     number of rows in this block is nloc
c
      nloc = split(my_id+2) - split(my_id+1)
      nmax = nloc+nloc
      call dmalloc(pil, nmax+1, 1)
      call dmalloc(pdl, nloc, 8)
c
c     send the ia array, need to convert it to row element count
c     processor 0 needs to compute the element counts and displacements
c     the displacements are C-style indices
c
      if (my_id .eq. 0) then
         call dmalloc(pnelm, nprocs, 1)
         call dmalloc(pdisp, nprocs, 1)
         do 10 i = 1, nprocs
            nelm(i) = split(i+1) - split(i)
            disp(i) = ig(split(i)) - 1
            split(i) = split(i) - 1
 10      continue
         do 20 i = ntot+1, 2, -1
            ig(i) = ig(i) - ig(i-1)
 20      continue
c$$$         write (10+my_id, *) ntot, (ig(i), i=2, ntot+1)
      endif
      call MPI_SCATTERV(dg, nelm, split, mpitype,
     &     dl, nloc, mpitype, 0, mpicom, ierr)
      if (ierr .ne. mpi_success) then
         print *, 'MPI_SCATTERV failed (bdiag). Error flag = ', ierr
         ierr = -3
         return
      endif
      call MPI_SCATTERV(ig(2), nelm, split, MPI_INTEGER,
     &     il(2), nloc, MPI_INTEGER, 0, mpicom, ierr)
      if (ierr .ne. mpi_success) then
         print *, 'MPI_SCATTERV failed (ia). Error flag = ', ierr
         ierr = -4
         return
      endif
c$$$      write (10+my_id, *) nloc, (il(i), i=2, nloc+1)
      il(1) = 1
      do 30 i = 2, nloc+1
         il(i) = il(i) + il(i-1)
 30   continue
c
c     send out the ja array
c
      nzloc = il(nloc+1) - 1
      nzmax = nzloc
      call dmalloc(pjl, nzmax, 1)
      call dmalloc(pal, nzmax, 8)
      if (my_id .eq. 0) then
         do 40 i = 1, nprocs-1
            nelm(i) = disp(i+1)-disp(i)
            split(i) = split(i) + 1
 40      continue
         split(nprocs) = split(nprocs) + 1
         nelm(nprocs) = nztot - disp(nprocs)
c$$$         write (10+my_id, *) 'nelm = ', (nelm(i), i=1,nprocs)
      endif
c$$$      write (10+my_id, *) 'nzloc = ', nzloc
c$$$      call flush(10+my_id)
      call MPI_SCATTERV(jg, nelm, disp, MPI_INTEGER,
     &     jl, nzloc, MPI_INTEGER, 0, mpicom, ierr)
      if (ierr .ne. mpi_success) then
         print *, 'MPI_SCATTERV failed (ja). Error flag = ', ierr
         ierr = -5
         return
      endif
c
c     send out the nonzero value array
c
      call MPI_SCATTERV(ag, nelm, disp, mpitype,
     &        al, nzloc, MPI_DOUBLE_PRECISION, 0, mpicom, ierr)
      if (ierr .ne. mpi_success) then
         print *, 'MPI_SCATTERV failed (a). Error flag = ', ierr
         ierr = -6
         return
      endif
c
      if (my_id .eq. 0) then
         call dmfree(pdisp)
         call dmfree(pnelm)
      endif
c
      return
c     end subroutine scatter_pair
      end
