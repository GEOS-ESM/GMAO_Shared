      subroutine hbread3(filename, mapopt, my_id, nprocs, nrow, split,
     &                   paval, pja, pia, nzmax, nmax, iou, mpicom,
     &                   ierr, msglvl)
      implicit none
      include 'mpif.h'
c     .. Scalar Arguments ..
      character*(*)      filename
      integer            ierr, iou, mapopt, mpicom, msglvl, my_id, nmax,
     &                   nprocs, nrow, nzmax
c     ..
c     .. Array Arguments ..
      integer            split(nprocs+1), ia(*), ja(*)
      double precision   aval(*)
      pointer            (paval,aval), (pia, ia), (pja, ja)
c     ..
c
c     This routine reads a Harwell-Boeing matrix from file and
c     distribute the matrix as NPROCS consecutive row blocks,
c     where NPROCS is the number of processors in the MPICOM
c     group.  The matrix must contain nonzero values.
c
c     The matrix is read on the last processor, as each piece of
c     data is read, it is sent to its destination.
c
C     SUGGESTION FOR READING RSA MATRICES --> expand the RSA matrices
c     into RUA format, then use this routine to read the symmetric
c     matrices in RUA format.
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
c     .. Local Scalars ..
      character          ch
      character*3        rhstyp, type
      character*8        key
      character*16       indfmt, ptrfmt
      character*20       rhsfmt, valfmt
      character*72       title
      integer            i, indcrd, ipe, iroot, j, ncol,
     &                   nelm, neltvl, nnz, nrhs, nrwindx, ptrcrd,
     &                   rhscrd, totcrd, totnz, totrow, valcrd
c     ..
c     .. Local Arrays ..
      integer            itmp(*), status(mpi_status_size)
      pointer            (pitmp, itmp)
c     ..
c     .. External Subroutines ..
      external           dmalloc, dmfree
c     ..
c     .. Intrinsic Functions ..
      intrinsic          ichar, index, lge, lle, max, mod
c     ..
c     .. Statement Functions ..
      logical            isnumber
c     ..
c     .. Data statements ..
      data               pitmp / 0 /
c     ..
c     .. Statement Function definitions ..
      isnumber(ch) = (lge(ch,'0') .and. lle(ch,'9'))
c     ..
c     .. Executable Statements ..
      ierr = 0
      totnz = 1
      totrow = 1
      iroot = nprocs - 1
      if (my_id.eq.iroot) then
c
c     open the given file
c
         open(iou, file = filename, status = 'OLD',
     &       iostat = ierr)
         if (ierr.ne.0) then
            print *, 'HBREAD3: can not open file ', filename, ierr
            call mpi_abort(mpicom, ierr)
         endif
c
c     read the header
c
         read (iou, FMT = 9999)title, key, totcrd, ptrcrd, indcrd,
     &      valcrd, rhscrd, type, totrow, ncol, totnz, neltvl, ptrfmt,
     &      indfmt, valfmt, rhsfmt
 9999    format(A72, A8, /5I14, /A3, 11X, 4I14, /2A16, 2A20)
c
         if (rhscrd.gt.0)read (iou, FMT = 9989)rhstyp, nrhs, nrwindx
 9989    format(A3, 11X, I14, I14)
         if (msglvl.gt.0) then
            write (*, '(A $)') 'HBREAD3:'
            i = index(filename, ' ') - 1
            if (i.lt.1) i = len(filename)
            if (msglvl.gt.2)print *, 'Opened Harwell-Boeing file ',
     &          filename(1:i)
            print *, key, '(', type, '),        Size: ', totrow, 'X',
     &         ncol, ',       NNZ: ', totnz
            if (msglvl.gt.1)print *, title
         endif
c
c     make sure that the matrix is assembled
c
         if (type(3: 3).ne.'A' .and. type(3: 3).ne.'a') then
            print *, 'HBREAD3 expects the matrix to be assembled.'
            call mpi_abort(mpicom, ierr)
         else if (valcrd.le.0) then
            print *, 'HBREAD3 expects the matrix to contain values.'
            call mpi_abort(mpicom, ierr)
         endif
c
c     figure out how many numbers are in each row.  If the format
c     does not start with an integer, one (1) is assumed.
c
         ptrcrd = 0
         i = index(ptrfmt, '(') + 1
 30      if (isnumber(ptrfmt(i: i))) then
            ptrcrd = 10*ptrcrd + ichar(ptrfmt(i: i)) - ichar('0')
            i = i + 1
            go to 30
         else if (ptrcrd.eq.0 .and. ptrfmt(i:i).eq.' ') then
            i = i + 1
            go to 30
         endif
         if (ptrcrd.eq.0)ptrcrd = 1
         indcrd = 0
         i = index(indfmt, '(') + 1
 40      if (isnumber(indfmt(i: i))) then
            indcrd = 10*indcrd + ichar(indfmt(i: i)) - ichar('0')
            i = i + 1
            go to 40
         else if (indcrd.eq.0 .and. indfmt(i:i).eq.' ') then
            i = i + 1
            go to 40
         endif
         if (indcrd.eq.0)indcrd = 1
         valcrd = 0
         i = index(valfmt, '(') + 1
 50      if (isnumber(valfmt(i: i))) then
            valcrd = 10*valcrd + ichar(valfmt(i: i)) - ichar('0')
            i = i + 1
            go to 50
         else if (valfmt(i:i).eq.'p' .or. valfmt(i:i).eq.'P') then
            valcrd = 0
            i = i + 1
            go to 50
         else if (valcrd.eq.0 .and. valfmt(i:i).eq.' ') then
            i = i + 1
            go to 50
         endif
         if (valcrd.eq.0)valcrd = 1
c
c     two ways to split the matrix
c
         neltvl = totrow / nprocs
         split(1) = 1
         if (mapopt.ne.-1) then
            neltvl = neltvl + 1
            do 60 i = 1, mod(totrow, nprocs)
               split(i+1) = split(i) + neltvl
 60         continue
            neltvl = neltvl - 1
            do 70 i = mod(totrow, nprocs) + 1, nprocs
               split(i+1) = split(i) + neltvl
 70         continue
            neltvl = neltvl + 1
         else
            neltvl = neltvl + 1
            do 80 i = 1, totrow / neltvl
               split(i+1) = split(i) + neltvl
 80         continue
            do 90 i = totrow / neltvl, nprocs
               split(i+1) = totrow + 1
 90         continue
         endif
      endif
c
c     broadcast split array to all processors
c
      call mpi_bcast(split, nprocs+1, mpi_integer, iroot, mpicom, ierr)
      if (ierr.ne.mpi_success) then
         print *, 'HBREAD3: PE', my_id, ' failed MPI_BCAST. IERR=', ierr
         call mpi_abort(mpicom, ierr)
      endif
c
c     the number of rows belong to this processor
c
      nrow = split(my_id+2) - split(my_id+1)
c
c     allocate enough space to store pointer array
c
      nmax = nrow + nrow
      if (my_id.eq.iroot) then
         nmax = max(neltvl+neltvl, neltvl+ptrcrd)
      endif
      call dmalloc(pia, nmax+1, 1)
      if (pia.eq.0) then
         print *, 'HBREAD3 failed to allocate IA array on PE', my_id
         call MPI_Abort(mpicom, ierr)
      endif
c
c     read the data file and send the pointer array segments
c
      j = 0
      do 110 ipe = 0, nprocs - 2
         if (my_id.eq.iroot) then
            nelm = split(ipe+2) - split(ipe+1) + 1
            if (nelm.gt.j) then
               if (mod(nelm-j,ptrcrd).gt.0) then
                  neltvl = ((nelm-j)/ptrcrd+1)*ptrcrd + j
               else
                  neltvl = nelm
               endif
               if (neltvl+split(ipe+1).gt.totrow)neltvl = totrow -
     &             split(ipe+1) + 2
               read (iou, FMT = ptrfmt)(ia(i), i = j+1, neltvl)
            else
               neltvl = j
            endif
            if (msglvl.gt.10)print *, 'PE', my_id, ' sending ', nelm,
     &          'IA elements to', ipe
            call mpi_send(ia, nelm, mpi_integer, ipe, 1, mpicom, ierr)
            nelm = nelm - 1
            j = neltvl - nelm
            do 100 i = 1, j
               ia(i) = ia(i+nelm)
 100        continue
         else if (my_id.eq.ipe) then
            if (msglvl.gt.10)print *, 'PE', my_id, ' waiting for ',
     &          nrow + 1, 'IA elements from', iroot
            call mpi_recv(ia, nrow+1, mpi_integer, iroot, 1, mpicom,
     &                    status, ierr)
            call mpi_get_count(status, mpi_integer, nelm, i)
            if (nelm.ne.nrow+1) then
               print *, 'HBREAD3: PE ', ipe, ' received ', nelm,
     &            'IA elements while expecting ', nrow + 1
               if (ierr.eq.mpi_success)ierr = -1
            endif
         endif
         if (ierr.ne.mpi_success) then
            print *, 'HBREAD3 failed while sending IA to PE ', ipe,
     &         '. Ierr=', ierr
            call mpi_abort(mpicom, ierr)
         endif
 110  continue
      if (my_id.eq.iroot .and. j.le.nrow) then
         read (iou, FMT = ptrfmt)(ia(i), i = j+1, nrow+1)
      endif
c
c     convert the pointer array to start with 1
c
      neltvl = ia(1) - 1
      do 120 i = 1, nrow + 1
         ia(i) = ia(i) - neltvl
 120  continue
c
c     number of nonzero elements on this processor
c
      nnz = ia(nrow+1) - ia(1)
      if (my_id.eq.iroot) then
         call dmalloc(pitmp, nprocs+nprocs, 1)
         if (pitmp.eq.0) then
            print *, 'HBREAD3: unable to allocate workspace on PE',
     &           my_id
            call mpi_abort(mpicom, ierr)
         endif
      endif
      call mpi_gather(nnz, 1, mpi_integer, itmp, 1, mpi_integer, iroot,
     &                mpicom, ierr)
      if (my_id.eq.iroot) then
         nzmax = itmp(1)
         itmp(nprocs+1) = 1
         do 130 i = 2, nprocs
            nzmax = max(nzmax, itmp(i))
            itmp(nprocs+i) = itmp(nprocs+i-1) + itmp(i-1)
 130     continue
         nzmax = nzmax + indcrd
      else
         nzmax = nnz
      endif
      call dmalloc(pja, nzmax, 1)
      call dmalloc(paval, nzmax, 8)
      if (pja.eq.0 .or. paval.eq.0) then
         print *, 'HBREAD3 failed to allocate JA/A array on PE', my_id
         call mpi_abort(mpicom, ierr)
      endif
c
c     read and distribute the index array
c
      j = 0
      do 150 ipe = 0, nprocs - 2
         if (my_id.eq.iroot) then
            nelm = itmp(ipe+1)
            if (nelm.gt.j) then
               if (mod(nelm-j,indcrd).gt.0) then
                  neltvl = j + ((nelm-j)/indcrd+1)*indcrd
               else
                  neltvl = nelm
               endif
               if (itmp(nprocs+ipe+1)+neltvl.gt.totnz)neltvl = totnz -
     &             itmp(nprocs+ipe+1) + 1
               read (iou, FMT = indfmt)(ja(i), i = j+1, neltvl)
            else
               neltvl = j
            endif
            if (msglvl.gt.10)print *, 'PE', my_id, ' sending ', nelm,
     &          'JA elements to', ipe
            call mpi_send(ja, nelm, mpi_integer, ipe, 2, mpicom, ierr)
            j = neltvl - nelm
            do 140 i = 1, j
               ja(i) = ja(i+nelm)
 140        continue
         else if (my_id.eq.ipe) then
            if (msglvl.gt.10)print *, 'PE', my_id, ' waiting for ', nnz,
     &          'JA elements from', iroot
            call mpi_recv(ja, nzmax, mpi_integer, iroot, 2, mpicom,
     &                    status, ierr)
            call mpi_get_count(status, mpi_integer, nelm, i)
            if (nelm.ne.nnz) then
               print *, 'HBREAD3: PE ', my_id, ' received ', nelm,
     &            'JA elements while expecting ', nnz
               if (ierr.eq.mpi_success)ierr = -1
            endif
         endif
         if (ierr.ne.mpi_success) then
            print *, 'HBREAD3 failed while sending JA to PE ', ipe,
     &         '. IERR=', ierr
            call mpi_abort(mpicom, ierr)
         endif
 150  continue
      if (my_id.eq.iroot .and. j.lt.nnz) then
         read (iou, FMT = indfmt)(ja(i), i = 1+j, nnz)
      endif
c
c     read and distribute the value array
c
      j = 0
      do 170 ipe = 0, nprocs - 2
         if (my_id.eq.iroot) then
            nelm = itmp(ipe+1)
            if (nelm.gt.j) then
               if (mod(nelm-j,valcrd).gt.0) then
                  neltvl = ((nelm-j)/valcrd+1)*valcrd + j
               else
                  neltvl = nelm
               endif
               if (itmp(nprocs+ipe+1)+neltvl.gt.totnz)neltvl = totnz -
     &             itmp(nprocs+ipe+1) + 1
               read (iou, FMT = valfmt)(aval(i), i = 1+j, neltvl)
            else
               neltvl = j
            endif
            if (msglvl.gt.10)print *, 'PE', my_id, ' sending ', nelm,
     &          'AVAL elements to', ipe
            call mpi_send(aval, nelm, mpi_real8, ipe, 3, mpicom, ierr)
            j = neltvl - nelm
            do 160 i = 1, j
               aval(i) = aval(i+nelm)
 160        continue
         else if (my_id.eq.ipe) then
            if (msglvl.gt.10)print *, 'PE', my_id, ' waiting for ', nnz,
     &          'AVAL elements from', iroot
            call mpi_recv(aval, nzmax, mpi_real8, iroot, 3, mpicom,
     &                    status, ierr)
            call mpi_get_count(status, mpi_real8, nelm, i)
            if (nelm.ne.nnz) then
               print *, 'HBREAD3: PE ', my_id, ' received ', nelm,
     &            'AVAL elements while expecting ', nnz
               if (ierr.eq.mpi_success)ierr = -1
            endif
         endif
         if (ierr.ne.mpi_success) then
            print *, 'HBREAD3 failed while sending AVAL to PE ', ipe,
     &         '. IERR=', ierr
            call mpi_abort(mpicom, ierr)
         endif
 170  continue
      if (my_id.eq.iroot) then
         if (j.lt.nnz) then
            read (iou, FMT = valfmt)(aval(i), i = 1+j, nnz)
         endif
         call dmfree(pitmp)
         close (iou)
         if (msglvl.gt.2) write (*, '(A)') 'HBREAD3: completed.'
      endif
c
      return
c     end of subroutine hbread3
      end
