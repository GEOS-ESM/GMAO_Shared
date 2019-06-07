c
c     This head file is used to define the common block for storing
c     the matrix information required to perform matrix-vector
c     multiplications with PSPARSLIB amxdis.
c
      INTEGER LMTLAN,LMTprs
      PARAMETER (LMTLAN = 25, lmtprs=10)
      integer MAX_PROCS
      parameter (MAX_PROCS=512)
c
c     This is intended to run preconditioned iterative method to
c     solve with the following operator
c     (A - delta I)
c     where the local part of A is stored in aloc, jaloc, ialoc
c     and an ILU type preconditioner is stored in alu, jlu, ju.
c     The edges connected to external nodes are listed in the
c     data structure aloc, jaloc, ialoc(nloc+1) [nloc-nbnd+1 rows].
c
c     locfact controls how the local factorization is done in
c     the preconditioner for iterative method
c     0 -- diagonal
c     1 -- ILU0 (from PSPARSLIB)
c     2 -- Full factorization (SuperLU, alu, jlu, ju structure not used)
c     kjac is the maxium allowed number of global (block) Jacobian
c     iterations
      integer locfact, kjac, lusize
c
c     The following is the header block of the matrix-vector
c     multiplication routine
c$$$  subroutine amxdis(nloc,nbnd,x,y,aloc,jaloc,ialoc,nproc,proc,ix,
c$$$ *     ipr,type,xlen,ptrn)
c$$$  implicit none
c$$$  integer nloc,nbnd,nproc,type,xlen,ptrn,myproc
c$$$  integer jaloc(*),ix(*),ialoc(*),ipr(nproc+1),proc(2*nproc)
c$$$  real*8 aloc(*),x(xlen),y(*)
      integer nbnd, nproc, amxtype, ptrn, myproc, totprocs, verbose
      integer PROC(4*MAX_PROCS+7), IPR(MAX_PROCS+1), SPLIT(MAX_PROCS+1)
      integer jaloc(*), ialoc(*), ix(*), jlu(*), ju(*)
      real*8 aloc(*), alu(*), delta
      pointer (paloc, aloc), (pjaloc, jaloc), (pialoc, ialoc),
     &     (pix, ix), (palu, alu), (pjlu, jlu), (pju, ju)
c
      common /pspmat/delta, paloc, pjaloc, pialoc, pix,
     &     palu, pjlu, pju, proc, ipr, split,
     &     nbnd, nproc, amxtype, ptrn, myproc, totprocs, verbose,
     &     locfact, lusize, kjac
C
C.... The common block "GETPUT" is used solely by the user-supplied
C.... subroutine STORE and holds the Lanczos vectors.
C
      REAL*8 QQ(*)
      pointer (pqq, qq)
      COMMON/GETPUT/pqq
C
C.... statistics of time and matvec count
C
      real*8 mmvtime, kmvtime, reorthtime
      integer mmvcnt, kmvcnt, nreorth
      common/shb_statistics/mmvtime, kmvtime, reorthtime,
     &     kmvcnt, mmvcnt, nreorth
