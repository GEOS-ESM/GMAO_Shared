   module m_fv2fv
   use m_const, only: kappa
   implicit none
   PRIVATE
   PUBLIC fv2fv
   interface fv2fv
   module procedure fv2fv_
   end interface
   CONTAINS
   subroutine fv2fv_ ( ps1,dp1,u1,v1,thv1,q1,o1,phis1,lm1, &
                       ps2,dp2,u2,v2,thv2,q2,o2,phis2,lm2,im,jm, &
                       w1, w2 )  ! optional arguments

!***********************************************************************
!  
!  Purpose
!     Driver for remapping of fv model levels at different suface ps
!
!  Argument Description
!
!     output:
!
!     ps1 ...... model surface  pressure
!     dp1 ...... model pressure thickness
!     u1 ....... model zonal      wind
!     v1 ....... model meridional wind
!     thv1 ..... model virtual potential  temperature
!     q1 ....... model specific   humidity
!     o1 ....... model ozone
!     w1 ....... model cloud water fraction (optional)
!     phis1 .... model surface geopotential
!     lm1 ...... model vertical   dimension
!
!     input:
!
!     ps2 ...... model surface  pressure
!     dp2 ...... model pressure thickness
!     u2 ....... model zonal      wind
!     v2 ....... model meridional wind
!     t2 . ..... model dry-bulb temperature
!     q2 ....... model specific   humidity
!     o2 ....... model ozone
!     w2 ....... model cloud water fraction (optional)
!     phis2 .... model surface geopotential
!     lm2 ...... model vertical   dimension
!
!     im ....... zonal      dimension
!     jm ....... meridional dimension
!
! 05Nov2003 Todling OpenMP in the last loop here is deactivated on the
!                   Compaq: version 551 of compiler does not like it.
! 14May2004 Todling Moved interface form m_xform here.
! 12Jan2005 Todling Replaced getcon by m_const
! 23Jan2005 B Zhang Updated for total cloud water fraction
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************
   use m_gmap, only : gmap
   use util, only : myout
   implicit none
   character(len=*), parameter :: myname = 'fv2fv'
   integer  im,jm,lm1,lm2
   logical same

! output 
! ------
   real      dp1(im,jm,lm1)
   real       u1(im,jm,lm1)
   real       v1(im,jm,lm1)
   real     thv1(im,jm,lm1)
   real       q1(im,jm,lm1)
   real       o1(im,jm,lm1)
   real,optional:: w1(im,jm,lm1)
   real      ps1(im,jm)

   real   phis1(im,jm)

! input 
! -----
   real     dp2(im,jm,lm2)
   real      u2(im,jm,lm2)
   real      v2(im,jm,lm2)
   real      t2(im,jm,lm2)
   real    thv2(im,jm,lm2)
   real      q2(im,jm,lm2)
   real      o2(im,jm,lm2)
   real,optional:: w2(im,jm,lm1)
   real     ps2(im,jm)
   real   phis2(im,jm)

! Local variables
! ---------------
   real   pe1(im,jm,lm1+1)
   real   pe2(im,jm,lm2+1)
   real   pk (im,jm,lm2  )
   real  pke1(im,jm,lm1+1)
   real  pke2(im,jm,lm2+1)

   real, allocatable, dimension(:) :: ak,bk
   integer i,j,L

! Size check
! ----------
      if (present(w1) .and. present(w2) ) then
          if (size(w1,1)/=im .and. size(w1,2)/=jm .and. size(w1,3)/=lm1) then
              print *, 'fv2fv: bad dimension for w1 '
              print *, 'fv2fv: program aborting ... '
              stop
          endif
          if (size(w2,1)/=im .and. size(w2,2)/=jm .and. size(w2,3)/=lm2) then
              print *, 'fv2fv: bad dimension for w2 '
              print *, 'fv2fv: program aborting ... '
              stop
          endif
      endif

! start

   same = .true.

! Compute edge-level pressures
! ----------------------------
   pe1(:,:,lm1+1) = ps1
   do L=lm1,1,-1
      pe1(:,:,L) = pe1(:,:,L+1)-dp1(:,:,L)    
   enddo                                  

#if   (openmp)
!$omp  parallel do &
!$omp& default (shared) &
!$omp& private (i,j,L)
#endif
   do L=1,lm1+1
      do j=1,jm
        do i=1,im
          pke1(i,j,L) = pe1(i,j,L)**kappa
        enddo
      enddo
   enddo
      
! Compute edge-level pressures
! ----------------------------
   pe2(:,:,lm1+1) = ps2
   do L=lm2,1,-1
      pe2(:,:,L) = pe2(:,:,L+1)-dp2(:,:,L)    
   enddo                                  
      
#if   (openmp)
!$omp  parallel do &
!$omp& default (shared) &
!$omp& private (i,j,L)
#endif
   do L=1,lm2+1
      do j=1,jm
        do i=1,im
          pke2(i,j,L) = pe2(i,j,L)**kappa
        enddo
      enddo
   enddo
      
! Map target analysis onto fv grid
! --------------------------------
   call myout ( 6, myname, 'map target grid onto new grid' )
   if (present(w1) .and. present(w2) ) then
   call gmap ( im,jm,1, kappa, &
               lm2,  pke2,  pe2, u2,  v2,  thv2,  q2, &
               lm1,  pke1,  pe1, u1,  v1,  thv1,  q1, &
               o_m=o2, o_n=o1, &
               w_m=w2, w_n=w1 )                       ! <= Optional I/O cloud water frac
   else
   call gmap ( im,jm,1, kappa, &
               lm2,  pke2,  pe2, u2,  v2,  thv2,  q2, &
               lm1,  pke1,  pe1, u1,  v1,  thv1,  q1, &
               o_m=o2, o_n=o1 )
   endif

   if ( lm1+1 > lm2+1 ) then
        print *, 'fv2fv: loop following this will have problems '
        print *, 'fv2fv: program aborting ... '
        stop
   endif

#ifdef sgi
#if   (openmp)
!$omp  parallel do &
!$omp& default (shared) &
!$omp& private (i,j,L)
#endif
#endif
col: do j=1,jm
row:   do i=1,im
lev:     do L=1,lm1+1
           if(pe1(i,j,L).ne.pe2(i,j,L)) then
             same = .false.
             exit col
            end if
       end do lev
     end do row
   end do col
   if(same) then
     u1 = u2
     v1 = v2
     thv1 = thv2
     q1 = q2
     o1 = o2
    if (present(w1) .and. present(w2) ) then
     w1 = w2
    endif
   end if

   end subroutine fv2fv_
   end module m_fv2fv
