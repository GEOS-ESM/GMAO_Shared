   module util

   implicit none

   interface myalloc
     module procedure alloc1d, alloc2d, alloc3d, alloc4d
   end interface

   contains

!-------------------------------------------------------------------------
   subroutine minmax_2d(name,data,l,m,n)
!------------------------------------------------------------------------------
!
!  This subroutine finds the min-max and the
!  mean of a 2D array. 
!
!  Input     :
!
!  data      : 2D array to be searched
!  l ,m      : First and second dimension of data
!
!  Output    :
!
!  amin,amax : minimun and maximum of the field
!------------------------------------------------------------------------------

   implicit none
   character(len=*), intent(in) :: name
   integer, intent(in) :: l,m
   integer, intent(in), optional :: n
!
   real :: data(l,m)
   integer :: k,i,j,imax,jmax,imin,jmin
   real :: amax, amin,amean
   real, parameter :: big = -9.99e33

   amax      = data(1,1)
   amin      = data(1,1)
   imin=1;jmin=1;imax=1;jmax=1
   amean     = 0.0
   do j = 1, m
     do i = 1, l
       if( abs(data(i,j)) < big ) then
         if (data(i,j) > amax) then
           imax=i; jmax=j
           amax = data(i,j)
         end if
         if (data(i,j) < amin) then
           imin=i; jmin=j
           amin = data(i,j)
         end if
       end if
       amean  = amean + data(i,j)/float(l*m)
     end do
   end do
   if(present(n)) then
     write(6,*)name,' min = ',amin,' at i = ',imin,' j = ',jmin,' level = ', n
   else
     write(6,*)name,' min = ',amin,' at i = ',imin,' j = ',jmin
   end if
   write(6,*)name,' max = ',amax,' at i = ',imax,' j = ',jmax,' mean  = ',amean

   end subroutine

!-------------------------------------------------------------------------
   SUBROUTINE myexit(routine, msg, status)
!-------------------------------------------------------------------------
   implicit none
   integer , intent(in) :: status
   character *(*) , intent(in) :: routine
   character *(*) , intent(in) :: msg

   if(status==0) then
      write (*,'(1x,a,a,a)') ' # ',routine,' executed successfully'
   else
      write (*,'(1x,a,a)') ' *** Error in routine ',routine
      write(* ,*)' *** rc = ',status,' -> ',msg
      call exit (status)
   end if

   end SUBROUTINE myexit

!-------------------------------------------------------------------------
   SUBROUTINE mywarn(routine, msg )
!-------------------------------------------------------------------------
   implicit none
   character *(*) , intent(in) :: routine
   character *(*) , intent(in) :: msg

   write (*,'(1x,a,a,a)') ' *** ',routine,' : non-zero return code'
   write(* ,*)' *** Reason : ',msg

   end SUBROUTINE mywarn

!-------------------------------------------------------------------------
   SUBROUTINE myout(lu, routine, msg )
!-------------------------------------------------------------------------
   implicit none
   integer , intent(in) :: lu
   character *(*) , intent(in) :: routine
   character *(*) , intent(in) :: msg

   write (lu,'(1x,a,a,a)') routine,' : ',msg

   end SUBROUTINE myout


!-------------------------------------------------------------------------
   SUBROUTINE alloc1d( arr, d1 )
!-------------------------------------------------------------------------
   implicit none
   integer , intent(in)        :: d1
   real, dimension(:), pointer :: arr
   integer                     :: err
   character(len=*), parameter :: myname = 'alloc1d'
 
   if(.not.associated(arr)) then
     allocate (arr(d1),stat=err)  
     if (err/=0) call myexit(myname,'no more memory',2)
     arr = 0.0
   else
     call mywarn(myname,'pointer is associated')
   end if
     
   end SUBROUTINE alloc1d

!-------------------------------------------------------------------------
   SUBROUTINE alloc2d( arr, d1, d2 )
!-------------------------------------------------------------------------
   implicit none
   integer , intent(in)        :: d1, d2
   real, dimension(:,:), pointer :: arr
   integer                     :: err
   character(len=*), parameter :: myname = 'alloc2d'
 
   if(.not.associated(arr)) then
     allocate (arr(d1,d2),stat=err)  
     if (err/=0) call myexit(myname,'no more memory',2)
     arr = 0.0
   else
     call mywarn(myname,'pointer is associated')
   end if
     
   end SUBROUTINE alloc2d

!-------------------------------------------------------------------------
   SUBROUTINE alloc3d( arr, d1, d2, d3 )
!-------------------------------------------------------------------------
   implicit none
   integer , intent(in)        :: d1, d2, d3
   real, dimension(:,:,:), pointer :: arr
   integer                     :: err
   character(len=*), parameter :: myname = 'alloc3d'
 
   if(.not.associated(arr)) then
     allocate (arr(d1,d2,d3),stat=err)  
     if (err/=0) call myexit(myname,'no more memory',2)
     arr = 0.0
   else
     call mywarn(myname,'pointer is associated')
   end if
     
   end SUBROUTINE alloc3d

!-------------------------------------------------------------------------
   SUBROUTINE alloc4d( arr, d1, d2, d3, d4 )
!-------------------------------------------------------------------------
   implicit none
   integer , intent(in)        :: d1, d2, d3, d4
   real, dimension(:,:,:,:), pointer :: arr
   integer                     :: err
   character(len=*), parameter :: myname = 'alloc4d'
 
   if(.not.associated(arr)) then
     allocate (arr(d1,d2,d3,d4),stat=err)  
     if (err/=0) call myexit(myname,'no more memory',2)
     arr = 0.0
   else
     call mywarn(myname,'pointer is associated')
   end if
     
   end SUBROUTINE alloc4d

!-------------------------------------------------------------------------
   SUBROUTINE open_file (unit, filename, ff)
!-------------------------------------------------------------------------
! open formatted or unformatted (sequential) files

   implicit none
   integer :: unit, error_number, ff
   character (len = 132) :: error_message, ffc, stat
   character (*) :: filename
   logical :: fatal

   if(ff==1) then
      ffc='formatted'
   else
      ffc='unformatted'
   end if

   open (unit   = unit        , &
         file   = filename    , &
         form   = ffc         , &
         iostat = error_number)

   error_open : if (error_number .ne. 0) then
      error_message = ' *** open_file: Error OPENing '//filename
      fatal = .true.
      call myexit(  filename, error_message, error_number )
   endif error_open
   write(*,'(1x,a,a)')' # Opened ',filename

   end SUBROUTINE open_file

!-------------------------------------------------------------------------
   subroutine writit (q,im,jm,lm,ku)
!-------------------------------------------------------------------------

   implicit none
   integer im,jm,lm,ku,L
   real   q (im,jm,lm)
   real*4 q2(im,jm)
   do L=lm,1,-1
   q2(:,:) = q(:,:,L)
   write(ku) q2
   enddo

   end subroutine

!-------------------------------------------------------------------------
   subroutine write_grads(field,lon,lat,lu)
!-------------------------------------------------------------------------
   implicit none
   integer lon,lat,lu
   real field(lon,lat)
   real*4 my_f(lon,lat)
   my_f = field
   write(lu) my_f
   end subroutine

   end module


