!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsnxtime:  read next date/time from ODS (or other diag) file
!
! !INTERFACE:

      subroutine odsnxtime ( ODSFile, nymd, nhms )

! !USES:
      use netcdf
      implicit none

      include 'ods_stdio.h'

! !INPUT PARAMETERS:

      character(*), intent(in)     ::  ODSFile

! !INPUT/OUTPUT PARAMETERS:

      integer      , intent(in out) ::  nymd
      integer      , intent(in out) ::  nhms
      
!
! !REVISION HISTORY:
!       ????????? - Redder   - Initial code.
!       ????????? - Todling  - Imported here from old iolib from C.Redder
!       15Apr2004 - Todling  - Added prologue; added support for diag_conv
!       20Dec2004 - Dee      - Support for diag_sat
!       03Mar2005 - Dee      - Fixed a bug introduced with 20Dec2004 revision
!       15Apr2019 - Sienkiewicz  - modify check for binary diag, add netcdf4 diag 
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      character(len=*), parameter :: myname_ = 'ods_next_time'

      integer     id, ierr, ios
      integer     idate, mobs 
      character*16  satype
      logical     conv
      integer     hour, ifrq, nsyn
      integer     first_jday, first_jhour
      integer     first_nymd, first_nhms
      integer     latest_jday, latest_hour, latest_jhour
      integer     latest_nymd, latest_nhms
      integer     jday, jhour, jhour2
      integer     status, ncid, nymdhh
      save        jhour
      data        jhour / -1 /
      logical     diff_file
      character * ( 255 )
     .            PrevODSFile
      save        PrevODSFile
      data        PrevODSFile / '!@#$%^&**()_+|' /   ! just garbage

      integer     ODS_CalDat
      external    ODS_CalDat

      integer     ODS_Handle
      external    ODS_Handle

!     Nothing to do, return
!     ---------------------
      if ( nymd .ne. -1 ) return

!     try opening as NetCDF (for nc4-diag or ODS)
!     -------------------------------------------
      status =  nf90_open(path=trim(ODSFile), mode = nf90_nowrite, ncid = ncid )
      
      if (status /= nf90_noerr) then ! not a netCDF file - treat as diag_bin
         call ods_dcscan ( .true., ODSFile, idate, mobs, conv, satype, ierr )
         if(ierr==0)then
            nymd = idate / 100
            nhms = 10000 * (idate - 100 * nymd)
            return
         else
            write(stderr,'(4a,i5)') myname_, 
     .                             ': Cannot open file ', trim(ODSFile), 
     .                             ' error code: ',ierr
            return
         end if
      else
!
!     determine if nc4/hdf file is nc4-diag by checking for date_time metadata
!     ------------------------------------------------------------------------
         status = nf90_inquire_attribute(ncid,NF90_GLOBAL,'date_time')
         if ( status == nf90_noerr ) then
!     
!     treat as nc4-diag file, otherwise continue to ODS processing
            status = nf90_get_att(ncid,NF90_GLOBAL,'date_time',nymdhh)
            nhms = mod(nymdhh,100)*10000
            nymd = int(nymdhh/100)
            status =  nf90_close(ncid)
            return
         end if
         status =  nf90_close(ncid)      !  ODS file, continue processing
      end if

!     open the ODS file:
!     -----------------
      call ODS_Open ( id, ODSFile, 'r', ierr )
      if ( ierr .ne. 0 ) then
      
         call ods_dcscan ( .true., ODSFile, idate, mobs, conv, satype, ierr )
         if(ierr==0)then
               nymd = idate / 100
               nhms = 10000 * (idate - 100 * nymd)
              return
         else
               write(stderr,'(4a,i5)') myname_, 
     .                      ': Cannot open file ', trim(ODSFile), 
     .                      ' error code: ',ierr
               return
         endif
	 
      end if

!     Get first Julian day on file
!     ----------------------------
      call ODS_IGet ( id, 'syn_beg:first_julian_day',
     .                     first_jday, ierr )
      if ( ierr .ne. 0 ) then
         write(stderr,'(4a,i5)') myname_,
     .                ': Error in first Julian day ',
     .                trim(ODSFile),
     .                ' error code: ',ierr
         call ODS_Close ( id, ' ', ierr )
         return
      end if
      first_jhour = 24 * ( first_jday - 1 )
      first_nymd  = ODS_CalDat ( first_jday )
      first_nhms  = 0

!     Get latest Julian day on file
!     -----------------------------
      call ODS_IGet ( id, 'syn_beg:latest_julian_day',
     .                     latest_jday, ierr )
      if ( ierr .ne. 0 ) then
         write(stderr,'(4a,i5)') myname_, 
     .                ': Error in latest Julian day ',
     .                 trim(ODSFile),
     .                ' error code: ',ierr
         call ODS_Close ( id, ' ', ierr )
         return
      end if
      latest_nymd = ODS_CalDat ( latest_jday )

!     Get synoptic hour of latest Julian day on file
!     ----------------------------------------------
      call ODS_IGet ( id, 'syn_beg:latest_synoptic_hour',
     .                     latest_hour, ierr )
      if ( ierr .ne. 0 ) then
         write(stderr,'(4a,i5)') myname_, 
     .                ': Error in latest synoptic hour ',
     .                 trim(ODSFile),
     .                ' error code: ', ierr
         call ODS_Close ( id, ' ', ierr )
         return
      end if
      latest_jhour = 24 * ( latest_jday - 1 ) + latest_hour
      latest_nhms  = latest_hour * 10000


!     Compare previous and current file names to
!     determine if the file is different
!     ------------------------------------------
      diff_file = .false.
      if ( trim(PrevODSFile) .ne. trim(ODSFile) ) then
         diff_file   = .true.
         PrevODSFile = trim(ODSFile)
      end if

!     Inquire frequency of data in file
!     ---------------------------------
      call ODS_IGet ( id, 'nsyn', nsyn, ierr )
      if ( ierr .ne. 0 ) then 
         nsyn = 4  ! set default
         write(stderr,'(3a,i5)') myname_, 
     .                ': Error cannot determined nsyn ',
     .                ' using default: ', nsyn
      end if
      ifrq = 24 / nsyn

!     Set Julian hour time tag for current set of data
!     ------------------------------------------------
      if      ( jhour .eq. -1 ) then ! Just opened file
         jhour = first_jhour         ! ----------------

      else if ( diff_file )     then ! Different file just opened
         jhour = first_jhour         ! --------------------------

      else
         jhour = jhour + ifrq        ! File already opened
                                     ! -------------------
      end if

!     Check to determine if more data exist.  More data exists if the
!     latest time tag on file is more recent than the current tag
!     ---------------------------------------------------------------
      if ( jhour .gt. latest_jhour ) then
         jhour = -1
         call ODS_Close ( id, ' ', ierr )
         return
      end if

!     Determine important tags for date and time
!     ------------------------------------------
      jday = jhour / 24 + 1
      hour = mod ( jhour, 24 )
      nymd = ODS_CalDat ( jday )
      nhms = 10000 * hour

      call ODS_Close ( id, ' ', ierr )

      return
      end subroutine odsnxtime
