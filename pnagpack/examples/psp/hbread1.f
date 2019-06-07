c
c     This subroutine read a matrix in Harwell-Boeing format
c     and distribute it to all members of MPICOM group.
c     It reads the whole matrix on one processor, performs
c     reordering and distribution the reordered matrix as
c     row-blocks.
c     If the matrix is RSA type, it expands into full matrix.
c     If the matrix is RUA type, it takes the symmetric part of the
c     matrix.
c
      subroutine hbread1(filename, mapopt, my_id, nprocs, nrow,
     &     split, paval, pja, pia, nzmax, nmax, iou, mpicom, ierr,
     &     msglvl)
      implicit none
      include 'mpif.h'
      integer nrow, nzmax, nmax, ierr, mapopt, msglvl, mpicom, iou
      integer nprocs, my_id
      integer ia(*), ja(*), split(nprocs+1)
      real*8  aval(*)
      character filename*(*)
      pointer (paval, aval), (pja, ja), (pia, ia)
c
c     This routine reads a Harwell-Boeing matrix from file and
c     distribute the matrix as NPROCS consecutive row blocks,
c     where NPROCS is the number of processors in the MPICOM
c     group.
c
c***  The matrix could be potentially reordered in this routine.
c***  however, there is no reordering information recorded.
c***  need to rework this if reordering is found to be actually
c***  desired.
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
c     nzmax     -- size of aval and ja
c     nmax      -- size of ia
c     ierr      -- error flag
c
      integer itmp(*), jtmp(*)
      real*8  vtmp(*)
      pointer (pitmp, itmp), (pjtmp, jtmp), (pvtmp, vtmp)
c     .. Local Scalars ..
      character*8        key
      character*72       title
      integer            totrow, totnz
c
c     arrays for partitioning the matrix
c
      integer map(*)
      pointer (pmap, map)
c     ..
      data pitmp, pjtmp, pvtmp, pmap/4*0/
c
      totrow = 1
      totnz = 1
      if (my_id.eq.0) then
c
c     read in the whole matrix and symmetrize it as needed
c
         totrow = 0
         totnz = 0
         call hbresy(msglvl, iou, filename, key, title, totrow,
     &        totnz, pvtmp, pjtmp, pitmp, ierr)
         if (ierr.ne.0) call MPI_ABORT(MPI_COMM_WORLD, ierr)
         totnz = itmp(totrow+1) - itmp(1)
c$$$         call dump(1,totrow,.true.,vtmp,jtmp,itmp,10+my_id)
c
c     .. decide what rows to send to each processor
c     this is a partitioning routine which probably also reorder
c     the matrix.  the input matrix is to be reordered so that
c     each process will receive one consecutive block rows of
c     the reordered matrix
c     
         call dmalloc(pmap, totrow, 1)
         if (pmap.eq.0) then
            print *, 'HBRAD1 failed to allocate map array.'
            call mpi_abort(MPI_COMM_WORLD, ierr)
         endif
         if (mapopt.eq.0) then
            call dummap(totrow, jtmp, itmp, nprocs, map, split)
         else if (mapopt.eq.-1) then
            call glimap(totrow, jtmp, itmp, nprocs, map, split)
         else
            call dummap(totrow, jtmp, itmp, nprocs, map, split)
         endif
         call permcsr(totrow, vtmp, jtmp, itmp, map)
         call dmfree(pmap)
      endif
c
c     .. call distribution routine
c
      call scatter_csr(totrow, totnz, pvtmp, pjtmp, pitmp, my_id,
     &     nprocs, split, nrow, paval, pja, pia, nzmax, nmax, mpicom,
     &     ierr)
c$$$      call dump(1,nrow,.true.,aval,ja,ia,10+my_id)
c
      if (my_id .eq. 0 .and. totrow.gt.nrow) then
         call dmfree(pitmp)
         call dmfree(pjtmp)
         call dmfree(pvtmp)
      endif
      return
c     end of subroutine hbread_mpi
      end
c-----------------------------------------------------------------------
      subroutine dummap(nrow, ja, ia, ndom, map, split)
      implicit none
      integer nrow, ndom, ia(nrow+1), split(ndom+1), map(nrow), ja(*)
c
c     a dummy partition routine which simply uniformly divide the the
c     matrix into NDOM row-blocks
c
      integer domsize, i, remainder
      do 10 i = 1, nrow
         map(i) = i
 10   continue
      split(1) = 1
      split(max(2,ndom+1)) = nrow+1
c
c     trivial case of 1 domain(processor) only
c
      if (ndom.le.1) return
c
c     more than one domain (processor)
c
      domsize = nrow/ndom
      remainder = nrow - domsize*ndom
      do 20 i = 1, ndom-1
         if (i .le. remainder) then
            split(i+1) = split(i) + domsize + 1
         else
            split(i+1) = split(i) + domsize
         endif
 20   continue
c
      return
c     end subroutine dummap
      end
c-----------------------------------------------------------------------
      subroutine glimap(nrow, ja, ia, ndom, map, split)
      implicit none
      integer nrow, ndom, ia(nrow+1), split(ndom+1), map(nrow), ja(*)
c
c     a map scheme that only leave the last process with different
c     size (required by PSGECSRMV by G. Li at Cray)
c     An example of bad mapping: Nrow=13, ndom=6, GLIMAP=(3 3 3 3 1 0).
c
      integer domsize, i
      do 10 i = 1, nrow
         map(i) = i
 10   continue
      split(1) = 1
      split(max(2,ndom+1)) = nrow+1
c
c     trivial case of 1 domain(processor) only
c
      if (ndom.le.1) return
c
c     more than one domain (processor)
c
      domsize = nrow/ndom
      if (domsize*ndom.ne.nrow) domsize = domsize+1
      do 20 i = 1, ndom-1
         split(i+1) = split(i) + domsize
 20   continue
c
      if (domsize*(ndom-1).gt.nrow) then
         print *, 'GLIMAP: too few rows per processor. ',
     &        'DO NOT use GLIMAP.'
         do 30 i = nrow/domsize+1, nrow
            split(i) = nrow+1
 30      continue
      endif
c
      return
c     end subroutine glimap
      end
c-----------------------------------------------------------------------
      logical function isnulmap(nn, map)
      integer nn, map(nn), i
c
c     verify if the given map is a trivial map
c
      isnulmap = .true.
      do 10 i = 1, nn
         if (map(i) .ne. i .and. map(i).ne.0) then
            isnulmap = .false.
            return
         endif
 10   continue
      return
c     end of function isnulmap
      end
c-----------------------------------------------------------------------
      subroutine permcsr(nrow, aval, ja, ia, map)
      implicit none
      integer nrow, ja(*), ia(nrow+1), map(nrow)
      real*8 aval(*)
c
c     permute a matrix stored in CSR format (symmetric permutation)
c     a,ja,ia- (INPUT/OUTPUT) permuted matrix is tored in the same array
c              internally, another copy is produced.
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
      if (ppntr.eq.0 .or. pinvmap.eq.0 .or. pindx.eq.0) then
         print *, 'PERMCSR failed to allcoate workspace.'
         stop
      endif
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
      call dmalloc(ptmp, nnz, 8)
      if (ptmp.eq.0) then
         print *, 'PERMCSR failed to allcoate workspace.'
         stop
      endif
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
      call dmfree(ptmp)
      call dmfree(ppntr)
      return
c     end subroutine permcsr
      end
c-----------------------------------------------------------------------
      subroutine permindx(nrow, ja, ia, map)
      implicit none
      integer nrow, ja(*), ia(nrow+1), map(nrow)
c
c     permute the indeces according to the give map
c     ja, ia - (INPUT/OUTPUT) permuted matrix is tored in the same array
c              internally, another copy is produced.
c     map    - (INPUT) the ith row in the new matrix is the map(i)th
c              row of the old one.
c
      logical isnulmap
      integer i, j, k, nnz, pntr(*), indx(*), invmap(*)
      pointer (ppntr, pntr), (pindx, indx), (pinvmap, invmap)
      data ppntr, pindx, pinvmap/3*0/
c
c     don't do anything if the map is a trivial map
c
      if (isnulmap(nrow, map)) return
c
c     really have to change the matrix
c
      nnz = ia(nrow+1)-ia(1)
      call dmalloc(pindx, nnz, 1)
      call dmalloc(ppntr, nrow+1, 1)
      call dmalloc(pinvmap, nrow, 1)
      if (pindx.eq.0 .or. ppntr.eq.0 .or. pinvmap.eq.0) then
         print *, 'PERMINDX failed to allcoate workspace.'
         stop
      endif
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
      call dmfree(ppntr)
      call dmfree(pindx)
      return
c     end subroutine permindx
      end
c-----------------------------------------------------------------------
      subroutine permval(nrow, aval, ia, map)
      implicit none
      integer nrow, map(nrow), ia(nrow+1)
      real*8 aval(*)
c
c     This is the companion routine for permindx. It permutes an array
c     that is supposed to be the nonzero value of a sparse matrix in
c     CSR format.
c
      logical isnulmap
      integer i,j,k,nnz, pntr(*)
      real*8 tmp(*)
      pointer (ptmp, tmp), (ppntr, pntr)
      data ptmp, ppntr/2*0/
c
c     nothing to do if the map is a trivial map
c
      if (isnulmap(nrow, map)) return
c
c     allocate temporary memory for workspace
c
      nnz = ia(nrow+1) - ia(1)
      call dmalloc(ptmp, nnz, 8)
      call dmalloc(ppntr, nrow+1, 1)
      if (ptmp.eq.0 .or. ppntr.eq.0) then
         print *, 'PERMVAL failed to allcoate workspace.'
         stop
      endif
c
c     the array ia is in the new ordering, needs to reconstruct the
c     old pointer array
c
      do 10 i = 1, nnz
         tmp(i) = aval(i)
 10   continue
      do 20 i = 1, nrow
         pntr(map(i)+1) = ia(i+1) - ia(i)
 20   continue
      pntr(1) = 1
      do 30 i = 2, nrow+1
         pntr(i) = pntr(i-1) + pntr(i)
 30   continue
      do 50 i = 1, nrow
         k = pntr(map(i))
         do 40 j = ia(i), ia(i+1)-1
            aval(j) = tmp(k)
            k = k + 1
 40      continue
 50   continue
c
      call dmfree(ppntr)
      call dmfree(ptmp)
      return
c     end subroutine permval
      end
c-----------------------------------------------------------------------
      subroutine scatter_csr(ntot, nztot, pag, pjg, pig, my_id, nprocs,
     &     split, nloc, pal, pjl, pil, nzmax, nmax, mpicom, ierr)
      implicit none
      integer ntot, nztot, nprocs,nloc, nzmax, nmax, mpicom, ierr,
     &     il(*), jl(*), split(nprocs+1), jg(*), ig(*)
      real*8 ag(*), al(*)
      pointer (pag, ag), (pjg, jg), (pig, ig)
      pointer (pal, al), (pjl, jl), (pil, il)
      include 'mpif.h'
c
c     distribute a CSR matrix on PE 0 to everyone in the group MPICOM
c     according to the array SPLIT
c
c     INPUT:
c     ntot   -- the global matrix size
c     nztot  -- total number of nonzero elements
c     ag,jg,ig -- the whole matrix in CSR format
c     my_id  -- MPI_RANK
c     nprocs -- number of processes 
c     split  -- the global index of the blocks
c     mpicom -- MPI communicator for the group
c
c     OUTPUT:
c     nloc   -- numver of rows in this processor
c     al,jl,il -- the submatrix (row block) on this processor
c               in CSR format
c     nmax   -- size of array IL
c     nzmax  -- size of array JL and AL
c     ierr   -- error flag
c
c     local variables
c
      integer i, my_id, nzloc, nelm(*), disp(*)
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
         return
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
      if (pil.eq.0) then
         print *, 'SCATTER_CSR failed to allcoate IL array.'
         stop
      endif
c
c     send the ia array, need to convert it to row element count
c     processor 0 needs to compute the element counts and displacements
c     the displacements are C-style indices
c
      if (my_id .eq. 0) then
         call dmalloc(pnelm, nprocs, 1)
         call dmalloc(pdisp, nprocs, 1)
         if (pnelm.eq.0 .or. pdisp.eq.0) then
            print *, 'SCATTER_CSR failed to allcoate workspace.'
            stop
         endif
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
      call MPI_SCATTERV(ig(2), nelm, split, MPI_INTEGER,
     &     il(2), nloc, MPI_INTEGER, 0, mpicom, ierr)
      if (ierr .ne. mpi_success) then
         print *, 'MPI_SCATTERV failed (ia). Error flag = ', ierr
         ierr = -3
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
      if (pjl.eq.0 .or. pal.eq.0) then
         print *, 'SCATTER_CSR failed to allcoate workspace.'
         stop
      endif
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
         ierr = -4
         return
      endif
c
c     send out the nonzero value array
c
      if (MPI_REAL8.eq.0) then
         call MPI_SCATTERV(ag, nelm, disp, MPI_DOUBLE_PRECISION,
     &        al, nzloc, MPI_DOUBLE_PRECISION, 0, mpicom, ierr)
      else
         call MPI_SCATTERV(ag, nelm, disp, MPI_REAL8,
     &        al, nzloc, MPI_REAL8, 0, mpicom, ierr)
      endif
      if (ierr .ne. mpi_success) then
         print *, 'MPI_SCATTERV failed (a). Error flag = ', ierr
         ierr = -5
         return
      endif
c
      if (my_id .eq. 0) then
         call dmfree(pdisp)
         call dmfree(pnelm)
      endif
c
      return
c     end subroutine scatter_csr
      end
c.......................................................................
c\Name: hbresy (HB read and symmetrize)
c
c\Description:
c     read in a matrix from Harwell-Boeing file given, return the
c     symmetric part of the matrix in CSR format. i.e., both upper and
c     lower triangular part of the matrix are stored.
c
c\Arguments:
c     On input the arguments nnz indicates the size of the
c     arrays pointed by PARY,PJA
c     New memory are allocated if the old memory is not sufficient for
c     the new matrix
c     PIA is always deallocated and reallocated.
c
c     iprt -- if > 0 the name and the size of the matrix will be printed
c             else, only the error message is printed
c     iou  -- I/O unit to be used
c     fname-- the HB file name
c     n    -- number of columns in the matrix
c     nnz  -- the array size of ary and ja
c     pary -- pointer to ary (the nonzero values)
c     pja  -- pointer to ja (the row indices)
c     pia  -- pointer to ia (the pointer to beginning of each row)
c     ierr -- the error flag, an error message will be printed
c
c.......................................................................
      subroutine hbresy(iprt, iou, fname, key, title, n, nnz,
     &     pary, pja, pia, ierr)
      implicit none
c     .. Scalar Arguments ..
      character*(*)      fname, title, key
      integer            ierr, iou, iprt, n, nnz
      pointer            (pary,ary), (pja,ja), (pia, ia)
c     ..
c     .. Local Scalars ..
      character*2        guesol
      character*3        type
      integer            job, m, ncol, nrhs, nrow, nz
c     ..
c     .. Local Arrays ..
      integer            ia(*), iau(*), iwk(*), ja(*)
      double precision   ary(*)
      pointer            (piau,iau), (piwk,iwk)
c     ..
c     .. External Subroutines ..
      external           dmalloc, dmfree, symmat, readmt
c     ..
c     .. Intrinsic Functions ..
      intrinsic          abs, mod
c     ..
      data piau, piwk /2*0/
c     .. Executable Statements ..
      if (iou.le.0) then
         iou = mod(abs(iou)+n, 97) + 3
      endif
c
      open(iou, file = fname, status = 'old', ERR = 10)
c
      m = n
      nz = nnz
      nrhs = 0
c     .. read the header only
      job = 0
      call readmt(n, nnz, job, iou, ary, ja, ia, ary, nrhs, guesol,
     &            nrow, ncol, nnz, title, key, type, ierr)
      if (ierr.ne.0) then
         print *, 'HBRESY: File ', fname,
     &      ' does not have  correct header.'
         return
      endif
      if (iprt.gt.0) then
         write(*, FMT = *)
         write(*, FMT = 9989)key, nrow, ncol, nnz
         write(*, FMT = *)title
      else
         write(*, FMT = 9999)key
      endif
 9999 format(1X, 'HBRESY: ', A8)
 9989 format(1X, 'HBRESY: ', A8, '  Size = ', I7, ' X ', I7,
     &      '  NNZ  = ', I10)
c
c     .. allocate desired memory
c
      nnz = nnz + nnz
      if (nz.lt.nnz) then
         if (nz.ne.0) then
            call dmfree(pary)
            call dmfree(pja)
         endif
         call dmalloc(pary, nnz, 8)
         call dmalloc(pja, nnz, 1)
      endif
      if (pia.ne.0) call dmfree(pia)
      call dmalloc(pia, nrow+1, 1)
      if (pary.eq.0 .or. pja.eq.0 .or. pia.eq.0) then
         print *, 'HBRESY failed allocate space for the matrix.'
         stop
      endif
c
c     .. rewind the I/O unit and read the matrix
c
      job = 2
      nrhs = 0
      rewind (iou)
      call readmt(nrow, nnz, job, iou, ary, ja, ia, ary, nrhs, guesol,
     &            n, m, nz, title, key, type, ierr)
      close (iou)
      if (job.lt.2)ierr = 6
      if (ierr.ne.0)print *, 'READMT: returned with error code ', ierr
      if (ierr.eq.1) then
         print *, 'The matrix was not read due to nmax < ncol :', nrow,
     &      m
      else if (ierr.eq.2) then
         print *, 'The matrix was not read due to nzmax < nnz :', nnz,
     &      nz
      else if (ierr.eq.3) then
         print *, 'The matrix was not read due to both (nmax < ncol)',
     &      ' and (nzmax < nnz)', nrow, m, nnz, nz
      else if (ierr.eq.4) then
         print *, 'Can''t read sparse vectors, e.g., RHS,',
     &      ' guess or solution(s).'
      else if (ierr.eq.5) then
         print *, 'RHS/guess/solution was not read due to insurficient',
     &      ' space in array RHS.'
      else if (ierr.eq.6) then
         print *, 'Matrix file does not contain values.'
      else if (ierr.ne.0) then
         print *, '*** Unspecified error!'
      endif
      if (ierr.ne.0) go to 20
c
c     if the matrix is not symmetric, get the symmtric part of it
c
      nz = ia(n+1) + ia(n+1)
      call dmalloc(piau, n, 1)
      call dmalloc(piwk, n+1, 1)
      if (piau.eq.0 .or. piwk.eq.0) then
         print *, 'HBRESY failed to allcoate workspace.'
         stop
      endif
      if (type.eq.'RSA' .or. type.eq.'rsa') then
         call ssrcsr(3, 1, n, ary, ja, ia, nz, ary, ja, ia, iau,
     &        iwk, ierr)
      else if (type.eq.'RUA' .or. type.eq.'rua') then
         call symmat(n, nz, ary, ja, ia, ary, ja, ia, iau, iwk, ierr)
      else
         print *, 'HBRESY can only handle RSA and RUA inputs.'
         STOP
      endif
      call dmfree(piau)
      call dmfree(piwk)
      if (ierr.ne.0) then
         print *, 'HBSYMM error = ', ierr
      endif
c
      go to 20
c
 10   print *, 'HBRESY: Unable to open file ', fname, ' on I/O unit ',
     &   iou, '.'
      print *, '        Please make sure the named file exist,'
      print *, '        and the desired I/O unit is free.'
      ierr = -9999
 20   return
c     end subroutine hbresy
      end
c-----------------------------------------------------------------------
      subroutine symmat(nrow, nzmax, a, ja, ia, ao, jao, iao, indu, iwk,
     &                  ierr)
      implicit none
c     .. Scalar Arguments ..
      integer            ierr, nrow, nzmax
c     ..
c     .. Array Arguments ..
      integer            ia(nrow+1), iao(nrow+1), indu(nrow),
     &                   iwk(nrow+1), ja(*), jao(nzmax)
      real*8             a(*), ao(nzmax)
c     ..
c-----------------------------------------------------------------------
c     Symmetric part of a CSR matrix (In-place)
c-----------------------------------------------------------------------
c     This subroutine computes the symmetric part of a matrix in
c     CSR format by computing Ao = (A + A')/2, where A' is A
c     transpose.
c
c     Typically this routine is used to expand the SSR matrix of
c     Harwell Boeing matrices, or to obtain a symmetrized graph of
c     unsymmetric matrices.
c
c     This routine is inplace, i.e., (Ao,jao,iao) may be same as
c     (a,ja,ia).
c
c     It is possible to input an arbitrary CSR matrix to this routine,
c     since there is no syntactical difference between CSR and SSR
c     format. It also removes duplicate entries and perform a partial
c     ordering. The output matrix has an order of lower half, main
c     diagonal and upper half after the partial ordering.
c-----------------------------------------------------------------------
c on entry:
c---------
c
c nrow  = column dimension of inout matrix
c a,
c ia,
c ja    = matrix in compressed sparse row format.
c
c nzmax = size of arrays ao and jao. SYMMAT will abort if the storage
c         provided in ao, jao is not sufficient to store A. See ierr.
c         It should be no less than NNZ(A) + number of off-diagonal
c         nonzero elements in A, where A is the input matrix.
c
c on return:
c----------
c ao, jao, iao
c       = output matrix in compressed sparse row format. The resulting
c         matrix is symmetric and is equal to A+A'. ao, jao, iao,
c         can be the same as a, ja, ia in the calling sequence.
c
c indu  = integer array of length nrow. INDU will contain pointers
c         to the beginning of upper traigular part.
c
c iwk   = integer work space (size nrow+1).
c
c ierr  = integer. Serving as error message. If the length of the arrays
c         ao, jao exceeds nzmax, ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c
c-----------------------------------------------------------------------
c     .. Local Scalars ..
      integer            i, ipos, j, k, kfirst, klast, ko, kosav, nnz
      real*8             tmp
c     ..
c     .. Executable Statements ..
      ierr = 0
      do 10 i = 1, nrow
         indu(i) = 0
         iwk(i) = 0
 10   continue
      iwk(nrow+1) = 0
c
c     .. compute number of elements in each row of (A'-D)
c     put result in iwk(i+1)  for row i.
c
      do 30 i = 1, nrow
         do 20 k = ia(i), ia(i+1) - 1
            j = ja(k)
            if (j.ne.i)
     &         iwk(j+1) = iwk(j+1) + 1
 20      continue
 30   continue
c
c     .. find addresses of first elements of ouput matrix. result in iwk
c
      iwk(1) = 1
      do 40 i = 1, nrow
         indu(i) = iwk(i) + ia(i+1) - ia(i)
         iwk(i+1) = iwk(i+1) + indu(i)
         indu(i) = indu(i) - 1
 40   continue
c.....Have we been given enough storage in ao, jao ?
      nnz = iwk(nrow+1) - 1
      if (nnz.gt.nzmax) then
         ierr = nnz
         return
      endif
c
c     .. copy the existing matrix (backwards).
c
      kosav = iwk(nrow+1)
      do 60 i = nrow, 1, -1
         klast = ia(i+1) - 1
         kfirst = ia(i)
         iao(i+1) = kosav
         kosav = iwk(i)
         ko = iwk(i) - kfirst
         iwk(i) = ko + klast + 1
         do 50 k = klast, kfirst, -1
            ao(k+ko) = a(k)
            jao(k+ko) = ja(k)
 50      continue
 60   continue
      iao(1) = 1
c
c     now add A' to A. Go through the structure of ao, jao, iao
c     that has already been copied. iwk(i) is the address
c     of the next free location in row i for ao, jao.
c
      do 80 i = 1, nrow
         do 70 k = iao(i), indu(i)
            j = jao(k)
            if (j.ne.i) then
               ipos = iwk(j)
               ao(ipos) = ao(k)
               jao(ipos) = i
               iwk(j) = ipos + 1
            else
               ao(k) = ao(k) + ao(k)
            endif
 70      continue
 80   continue
c
c     scale every nonzero value by half
c
      tmp = 0.5D0
      do 85 i = 1, iao(nrow+1)-1
         ao(i) = tmp*ao(i)
 85   continue
c
c     .. eliminate duplicate entries --
c     array INDU is used as marker for existing indices, it is also the
c     location of the entry.
c     IWK is used to stored the old IAO array.
c     matrix is copied to squeeze out the space taken by the duplicated
c     entries.
c
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = iao(i)
 90   continue
      iwk(nrow+1) = iao(nrow+1)
      k = 1
      do 120 i = 1, nrow
         iao(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = jao(ipos)
            if (indu(j).eq.0) then
c     .. new entry ..
               indu(j) = k
               jao(k) = jao(ipos)
               ao(k) = ao(ipos)
               k = k + 1
            else
c     .. duplicate entry ..
               ao(indu(j)) = ao(indu(j)) + ao(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
c     .. remove marks before working on the next row ..
         do 110 ipos = iao(i), k - 1
            indu(jao(ipos)) = 0
 110     continue
 120  continue
      iao(nrow+1) = k
c
c     .. partial ordering ..
c     split the matrix into strict upper/lower triangular
c     parts, INDU points to the beginning of the upper part.
c
      do 140 i = 1, nrow
         klast = iao(i+1) - 1
         kfirst = iao(i)
 130     if (klast.gt.kfirst) then
            if (jao(klast).le.i .and. jao(kfirst).ge.i) then
c     .. swap klast with kfirst ..
               j = jao(klast)
               jao(klast) = jao(kfirst)
               jao(kfirst) = j
               tmp = ao(klast)
               ao(klast) = ao(kfirst)
               ao(kfirst) = tmp
            endif
            if (jao(klast).gt.i)
     &         klast = klast - 1
            if (jao(kfirst).lt.i)
     &         kfirst = kfirst + 1
            go to 130
         endif
c
         if (jao(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
c
c     .. order the entries according to column indices
c     bubble-sort is used
c
      do 190 i = 1, nrow
         do 160 ipos = iao(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  tmp = ao(k)
                  ao(k) = ao(j)
                  ao(j) = tmp
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i)+1, iao(i+1)-1
            do 170 j = iao(i+1)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  tmp = ao(k)
                  ao(k) = ao(j)
                  ao(j) = tmp
               endif
 170        continue
 180     continue
 190  continue
c
      return
c---- end of symmat ----------------------------------------------------
c-----------------------------------------------------------------------
      end
