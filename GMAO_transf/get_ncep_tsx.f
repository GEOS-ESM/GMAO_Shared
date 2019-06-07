      subroutine get_ncep_tsx( filename,ts,im,jm,nymd,nhms )
c program to read ts from ncep pressure level analysis file
c ---------------------------------------------------------
      use m_hinterp, only : hinterp
      implicit none
      integer  im,jm,nymd,nhms
      real  ts(im,jm)
      character*80 filename

      real,    allocatable :: dum(:,:)
      real*4,  allocatable :: bum(:,:)

      integer, parameter :: imax = 360
      integer, parameter :: jmax = 181
      integer, parameter :: nrec = 17+5*26  ! 17 surface fields + 5 upper air fields

      real*4, parameter :: undef4 = 9.999E+20
      real    undef
      logical defined
      integer L,i,j,m,n

      open (30,file=trim(filename),form='unformatted',access='direct',recl=imax*jmax*4)

      allocate ( bum(imax,jmax) )
      allocate ( dum(imax,jmax) )

      undef = undef4

c surface temperature
c -------------------
      m = nhms/60000
      read(30,rec=11+m*nrec) ((bum(i,j),i=1,imax),j=1,jmax)
      dum(:,:) = bum(:,:)

      if( imax.ne.im .or. jmax.ne.jm ) call hinterp ( dum,imax,jmax,ts,im,jm,1,undef,1,.false. )

      deallocate ( dum,bum )

      close (30)
      return
      end
