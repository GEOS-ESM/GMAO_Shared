   program getgfiodim

   implicit none

   integer, parameter :: READ_ONLY = 1
   character(len=255) fname, argv
   integer  iarg, argc
   integer fid, err
   integer im, jm, km, lm, nvars, ngatts

   argc =  command_argument_count()
   if ( argc .lt. 1 ) then 
       print *
       print *, "Usage: getgfiodim.x gfiofilename"
       print *
       error stop 1
   endif
   iarg = 1
   call get_command_argument ( iArg, argv )
   fname = trim(argv)

!  Open the file
!  -------------
   call GFIO_Open ( trim(fname), READ_ONLY, fid, err )
   if ( err .ne. 0 ) then
      print *, 'cannot open file '
      error stop 2
   endif
 
!  Get dimensions
!  --------------
   call GFIO_DimInquire ( fid, im, jm, km, lm, nvars, ngatts, err)
   if ( err .ne. 0 ) then
      print *, 'cannot extra dimensions '
      error stop 3
   end if

!   Close GFIO file
!   ---------------
    call GFIO_close ( fid, err )

    write(6,'(i5,1x,i5,1x,i5)') im, jm, km
    end
