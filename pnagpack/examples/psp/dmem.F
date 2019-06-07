c     file name: dmem.F
c
c     This is a fortran version of dynamic memory (de)allocation routine.
c
c     One must compile this file with -D option to define either
c     SUN or CRAY. Only SUN and CRAY dynamic memory allocation
c     scheme are implemented here!
c
c     to generate dmem.o on SUN/SGI
c     f77 -DSUN [your other options] -c dmem.F
c
c     to generate dmem.o on CRAY
c     cf77 -DCRAY [your other options] -c dmem.F
c
c     ** DO NOT change to .f extension **
c     !! .F extension tells the compilers to put the code through
c     !! a cpp type of preprocessor first.
c
c     The routine is intended to dynamically (de)allocate arrays (one
c     at a time) for fortran programs.
c
c     two routines are defined:
c     DMALLOC -- dynamic memory allocation
c     DMFREE -- deallocation routine
c
c     interface of DMALLOC:
c
c     subroutine dmalloc(ptr, nelm, type)
c
c     PTR: POINTER (OUTPUT)
c     if the routine was successful in allocating the desired memory,
c     it is the start position of the memory allocated, else it returns
c     zero (0).
c
c     NELM: INTEGER (INPUT)
c     number of elements desired.
c
c     TYPE: INTEGER (INPUT)
c     element type
c     data type convention (as used in TOOLPACK by Cowell and Garbow)
c     1 integer
c     2 real
c     3 logical
c     4 complex
c     5 double precision
c     6 character
c     7 double complex
c     8 8-byte elements (whatever type it is, typically real*8)
c
c
c     interface of DMFREE:
c
c     subroutine dmfree(ptr)
c
c     PTR: POINTER (INPUT)
c     the pointer of the array allocated through DMALLOC
c
#if defined(SUN)
      subroutine dmalloc(ptr,nelm,type)
      integer nelm,type
      character dummy(*)
      pointer (ptr,dummy)
      integer size(8), nbytes, i, malloc
      external malloc
c     .. definition of data size on 4-byte word workstations ..
      data size/4,4,4,8,8,1,16,8/
c
c     -- allocate one extra element for potential debug use
c     -- always try to allocate multiple of 8 bytes, an attempt to
c        optimize data alignment
c
      if (nelm.le.0 .or. type.le.0 .or. type.gt.8) then
         print *, 'DMALLOC: illegal arguments. NELM=',nelm,
     &        '  TYPE=',type
         ptr = 0
      endif
c
      nbytes = (nelm+1)*size(type)
      if (mod(nbytes, 8).ne.0) then
         nbytes = 8*(nbytes/8 + 1)
      endif
      ptr = malloc(nbytes)
      if (ptr.eq.0) then
         print *, 'DMALLOC: malloc failed to allocate ', nbytes,
     $        ' bytes.'
         return
      endif
c
c     clear the extra bytes to zero
c
      do 10 i = nelm*size(type)+1, nbytes
         dummy(i) = char(0)
 10   continue
c
      return
      end
c
      subroutine dmfree(ptr)
      integer dummy(1)
      pointer (ptr, dummy)
      if (ptr.ne.0) then
         call free(ptr)
         ptr = 0
      endif
      return
      end
#endif

#if defined(CRAY)
      subroutine dmalloc(ptr,nelm,type)
      integer nelm,type,dummy(*)
      pointer (ptr,dummy)
      integer i, nwords, size(8)
c     .. definition of data size on Cray machines ..
      data size/1,1,1,2,2,1,4,1/
c
c     -- allocate one extra element for possible debug use
c
      if (nelm.le.0 .or. type.le.0 .or. type.gt.8) then
         print *, 'DMALLOC: illegal arguments. NELM=',nelm,
     &        '  TYPE=',type
         ptr = 0
      endif
c
      nwords = (nelm+1)*size(type)
      call hpalloc(ptr, nwords, i, 0)
      if (i.ne.0) then
         print *, 'DMALLOC: hpalloc failed to allocate ',
     $        nwords, ' words.'
         print *, '         hpalloc error code = ', i
         return
      endif
c
c     clear the extra words to zero
c
      do 10 i = 1+nelm*size(type), nwords
         dummy(i) = 0
 10   continue
c
      return
      end
c
      subroutine dmfree(ptr)
      integer dummy(1),ierr
      pointer (ptr, dummy)
      if (ptr.ne.0) then
         call hpdeallc(ptr,ierr,0)
         ptr = 0
         if (ierr.ne.0) then
            print *, 'DMFREE: hpdeallc returns with error code = ', ierr
         endif
      endif
      return
      end
#endif
