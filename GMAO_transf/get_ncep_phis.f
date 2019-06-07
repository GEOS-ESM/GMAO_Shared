      subroutine get_ncep_phis ( filename,phis,im,jm,ierr )

c this routine reads a file containg the NCEP topography
c the topography was extracted from an NCEP first guess file     
c using ss2fv.x and printed out immediately after gg2ll
!
! 01Aug2003 Todling  Replaced lu=30 by luvail call

      use m_ioutil, only : luavail
      implicit none
      integer     im,jm,i,j
      integer     lu
      integer     ierr,ios
      real phis(im,jm)
      real dum(im,jm)

      character(len=*) filename

c surface geopotential

!     lu = 30
      ierr = 0
      lu = luavail()
      open (lu,file=trim(filename),form='unformatted',
     &      access='sequential',iostat=ios)
      if (ios .ne. 0) then
        ierr = 41
        print *,'error opening ',filename
        return
      endif
      read (lu,iostat=ios) ((dum(i,j),i=1,im),j=1,jm)
      if (ios .ne. 0) then
        ierr = 42
        return
      endif
      close(lu)

      phis = dum

      return
      end
