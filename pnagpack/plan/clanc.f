c-------------------------------------------------------------------------
c !ROUTINE: clanc --- call for the lanczos driver plandr and post process
c Interface to main planso driver and post-process the computed eigen vectors
c
c  23jul2007 Kim    - Initial implementation
c---------------------------------------------------------------------------
      subroutine clanc(n,lanmax,ev,kappa,
     *          j,neig,Lvec,ritz,bnd,w,nw,nvec,indev,indvec,ierr,mpicom)

      implicit none

c !INPUT PARAMETERS:
      integer, intent(in)  :: n, lanmax, ev, nw, mpicom
c n      dimension of the eigenproblem
c lanmax upper limit to the number of lanczos steps 
c ev     .le.0 means eigenvalues only,
c        .gt.0 means both eigenvalues and eigenvectors are wanted
c              and ev becomes the output channel for eigenvectors
c nw     length of the work array w
c mpicom mpi communicator

      real*8, intent(in) :: kappa 
c relative accuracy of ritz values acceptable as eigenvalues
c (kappa will not be touched and/or looked at if ev.le.0)

      integer,  intent(inout) :: j, neig, nvec, ierr
c j      number of lanczos steps actually taken
c neig   number of ritz values stabilized
c nvec   number of accepted ritz values
c ierr   error flag

      real*8, intent(inout) :: Lvec(n,lanmax), ritz(lanmax), bnd(lanmax), w(nw)
c Lvec   array to hold the eigen vectors
c ritz   array to hold the ritz values
c bnd    array to hold the error bounds
c w      work array of length nw

      integer, intent(inout)  :: indev(lanmax) 
c convergence indicatgor for computed eigen values
      integer, intent(inout)  :: indvec(lanmax) 
c convergence indicatgor for computed eigen vectors

c !LOCAL VARIALBLES
      integer i, maxprs, msglvl

      real*8 eps, condm, endl, endr

      real*8,  allocatable :: astore(:)

      maxprs = lanmax           
      condm  = 1.000            
      endl   = 0.000            
      endr   = 0.1e-04
          
	              
      msglvl = 15         
	              
		      

c
c Storage for the computed Lanczos vectors
c ----------------------------------------
      allocate ( astore(n*(lanmax+2)) )

c Lanczos package main driver
c---------------------------------------------------------------
      call plandr(n,lanmax,astore,maxprs,condm,endl,endr,ev,kappa, 
     *         j,neig,lvec,ritz,bnd,w,nw,ierr,msglvl,mpicom)

c on successful exit from landr do three things: 
c---------------------------------------------------------------
c  (1) change order of eigenvalues and bounds from largest
c      to smallest. 
c  (2) determine the number of eigenvectors written out.
c      this could be done more directly 
c      in ritvec but is done here for 
c      minimal changes in the landr-package. this information 
c      is needed for routine vretr.
c  (3) take care of the index vectors indev and indvec

      if (ierr.eq.0) then 

         do i = 1, lanmax
            w(i)        = ritz(i)
            w(i+lanmax) = bnd (i)
         enddo
         do i = 1 , lanmax
            ritz(i) = w(lanmax-i+1)
            bnd(i)  = w(2*lanmax-i+1)
         enddo
         do i = 1, lanmax 
            indev(i)  = 0
            indvec(i) = 0 
         enddo

         nvec = 0 
         if (ev.gt.0) then 
            do i = 1, lanmax
               if (bnd(i).le.kappa*abs(ritz(i))) then 
                  nvec = nvec + 1
                  indvec(i) = 1
               endif
            enddo
         endif

         call deteps (eps)
         do i = 1, lanmax
            if (bnd(i).le.16.000*eps*abs(ritz(i))) indev(i) = 1
         enddo

      endif

      deallocate ( astore )
      return
      end subroutine clanc
