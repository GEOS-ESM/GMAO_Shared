      module m_hinterp
      private
      public hinterp

      contains

      subroutine hinterp ( qin,iin,jin,qout,iout,jout,mlev,undef,msgn,check )
      implicit   none
      integer    iin,jin,       iout,jout, mlev,msgn
      real   qin(iin,jin,mlev), qout(iout,jout,mlev)
      real undef,pi,dlin,dpin,dlout,dpout
      real dlam(iin), lons(iout*jout), lon
      real dphi(jin), lats(iout*jout), lat
      integer i,j,loc
      logical check

      pi = 4.0*atan(1.0)
      dlin = 2*pi/iin
      dpin = pi/(jin-1)
      dlam(:) = dlin
      dphi(:) = dpin

      dlout = 2*pi/iout
      dpout = pi/(jout-1)
      
      loc = 0
      do j=1,jout
      do i=1,iout
      loc = loc + 1
      lon = -pi + (i-1)*dlout
      lons(loc) = lon
      enddo
      enddo

      loc = 0
      do j=1,jout
      lat = -pi/2.0 + (j-1)*dpout
      do i=1,iout
      loc = loc + 1
      lats(loc) = lat
      enddo
      enddo

      call interp_h ( qin,iin,jin,mlev,
     .                dlam,dphi,0.0,90.0,0.0,
     .                qout,iout*jout,lons,lats,
     .                msgn,3,check,undef )

      return
      end subroutine hinterp
      end module m_hinterp
