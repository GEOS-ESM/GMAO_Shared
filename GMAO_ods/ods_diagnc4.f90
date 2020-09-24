module diagnc4_read
  use netcdf
  implicit none
  
! Integer types
  integer, parameter, public  :: i_byte  = selected_int_kind(1)      ! byte  integer
  integer, parameter, public  :: i_short = selected_int_kind(4)      ! short integer
  integer, parameter, public  :: i_long  = selected_int_kind(8)      ! long  integer
  integer, parameter, private :: llong_t = selected_int_kind(16)     ! llong integer
  integer, parameter, public  :: i_llong = max( llong_t, i_long )

! Expected 8-bit byte sizes of the integer kinds
  integer, parameter, public :: num_bytes_for_i_byte  = 1
  integer, parameter, public :: num_bytes_for_i_short = 2
  integer, parameter, public :: num_bytes_for_i_long  = 4
  integer, parameter, public :: num_bytes_for_i_llong = 8

! Define arrays for default definition
  integer, parameter, private :: num_i_kinds = 4
  integer, parameter, dimension( num_i_kinds ), private :: integer_types = (/ &
       i_byte, i_short, i_long,  i_llong  /)
  integer, parameter, dimension( num_i_kinds ), private :: integer_byte_sizes = (/ &
       num_bytes_for_i_byte, num_bytes_for_i_short, &
       num_bytes_for_i_long, num_bytes_for_i_llong  /)
  integer, parameter, private :: default_integer = 3  ! 1=byte,
                                                      ! 2=short,
                                                      ! 3=long,
                                                      ! 4=llong
  integer, parameter, public  :: i_kind = integer_types( default_integer )
  integer, parameter, public  :: r_single = selected_real_kind(6)  ! single precision
  integer, parameter, public  :: r_double = selected_real_kind(15) ! double precision

  interface get_1d_var
     module procedure get_1d_double_var_
     module procedure get_1d_real_var_
     module procedure get_1d_int_var_
     module procedure get_1d_char_var_
  end interface get_1d_var

  
contains

  integer function get_dim(ncid,name,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    integer,           intent(out):: rc
    integer                       :: dimid
    integer                       :: status
    integer                       :: length

    rc = 0
    status = nf90_inq_dimid(ncid,name,dimid)
    if (status /= nf90_noerr) then
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_inquire_dimension(ncid,dimid,len=length)
    if (status /= nf90_noerr) then
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    get_dim = length
    return
  end function get_dim
    
  subroutine get_1d_double_var_(ncid,name,values,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    real(r_double)    ,intent(out):: values(:)
    integer           ,intent(out):: rc
    integer                       :: status
    integer                       :: varid
    rc = 0
    status = nf90_inq_varid(ncid,name,varid)
    if (status /= nf90_noerr) then
       print *,'inquire error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_get_var(ncid, varid, values)
    if (status /= nf90_noerr) then
       print *,'error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    return
  end subroutine get_1d_double_var_
  
  subroutine get_1d_real_var_(ncid,name,values,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    real(r_single)    ,intent(out):: values(:)
    integer           ,intent(out):: rc
    integer                       :: status
    integer                       :: varid
    rc = 0
    status = nf90_inq_varid(ncid,name,varid)
    if (status /= nf90_noerr) then
       print *,'inquire error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_get_var(ncid, varid, values)
    if (status /= nf90_noerr) then
       print *,'error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    return
  end subroutine get_1d_real_var_
  
  subroutine get_1d_int_var_(ncid,name,ivalues,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    integer(i_kind)   ,intent(out):: ivalues(:)
    integer           ,intent(out):: rc
    integer                       :: status
    integer                       :: varid
    rc = 0
    status = nf90_inq_varid(ncid,name,varid)
    if (status /= nf90_noerr) then
       print *,'inquire error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_get_var(ncid, varid, ivalues)
    if (status /= nf90_noerr) then
       print *,'error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    return
  end subroutine get_1d_int_var_
  
  subroutine get_1d_char_var_(ncid,name,cvalues,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    character(len = *),intent(out):: cvalues(:)
    integer           ,intent(out):: rc
    integer                       :: status
    integer                       :: varid
    rc = 0
    status = nf90_inq_varid(ncid,name,varid)
    if (status /= nf90_noerr) then
       print *,'inquire error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_get_var(ncid, varid, cvalues)
    if (status /= nf90_noerr) then
       print *,'error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    return
  end subroutine get_1d_char_var_

  subroutine get_2d_int_var(ncid,name,values,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    integer(i_kind)   ,intent(out):: values(:,:)
    integer           ,intent(out):: rc
    integer                       :: status
    integer                       :: varid
    rc = 0
    status = nf90_inq_varid(ncid,name,varid)
    if (status /= nf90_noerr) then
       print *,'inquire error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_get_var(ncid, varid, values)
    if (status /= nf90_noerr) then
       print *,'error reading ', name,', status= ',status
       rc = 2
       status = nf90_close(ncid)
       return
    end if
  end subroutine get_2d_int_var

  subroutine get_2d_double_var(ncid,name,values,rc)
    integer,           intent(in) :: ncid
    character(len = *),intent(in) :: name
    real(r_double)    ,intent(out):: values(:,:)
    integer           ,intent(out):: rc
    integer                       :: status
    integer                       :: varid
    rc = 0
    status = nf90_inq_varid(ncid,name,varid)
    if (status /= nf90_noerr) then
       print *,'inquire error reading ', name
       rc = 2
       status = nf90_close(ncid)
       return
    end if
    status = nf90_get_var(ncid, varid, values)
    if (status /= nf90_noerr) then
       print *,'error reading ', name,', status=',status 
       rc = 2
       status = nf90_close(ncid)
       return
    end if
  end subroutine get_2d_double_var

end module diagnc4_read

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ods_diagnc4:  get data from diag\_nc4 file and convert to ODS
!
! !INTERFACE:

subroutine ods_diagnc4(fname, nymd, nhms, ods, rc)

! !USES
  
  use m_odsmeta
  use m_ods
  use m_odsxsup, only : getodsmeta
!  use netcdf
  use diagnc4_read
  use m_Sndx,    only : setSndx
  use m_ods_obsdiags, only : ods_obsdiags_getparam, ods_obsdiags


  implicit none

! !INPUT PARAMETERS:
 
  character(len=*), intent(in)   :: fname   ! GSI diag_ file name
  integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
  integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

! !OUTPUT PARAMETERS:
 
  type(ods_vect), intent(inout)  :: ods     ! ODS vector
  
  integer, intent(out)           :: rc      ! Error return code:
  
! !DESCRIPTION: get data from GSI diag\_nc4 files and convert to ODS
!
! !REVISION HISTORY:
!   2019-04-16 - Sienkiewicz - Initial code from standalone nc4diag2ods
!   2019-04-22 - Sienkiewicz - added 'spd' and 'sst' processing
!   2019-05-21 - Sienkiewicz - replace nc_diag_read routines
!   2019-08-28 - Sienkiewicz - polymorphic get_1d_var, initial implementation
!                               of sensitivity processing 
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  integer input_id
  real(r_single)  zero_single,tiny_single
  real small_num
  logical verbose
  logical satknown
  real, parameter :: undef = 1.e15
  integer nchan_dim, nspot, i, j, jj, n, knt, isat, kidsat, ii
  integer clen, stlen
  integer statks
  integer status
  integer idimid, varid
  character(len=5)                         :: statid

  logical :: lobsdiagsave, lobssens, ladjsigo
  integer miter

  character(len=20) isis
  character(len=10) dplat,satype,diag_type

  integer iiuse, inldpt, itldpt, iobssen
  integer iunldpt, ivnldpt, iutldpt, ivtldpt, iuobssen, ivobssen
  integer ndiag, ioff
  integer, allocatable :: obs_iuse(:,:)
  real(r_single)     :: nlomx(2), tlomx(2), obimp(2)
  real(r_single), allocatable :: pdata(:,:)
  real(r_double), allocatable :: u_obs_nldepart(:,:), u_obs_tldepart(:,:),  &
       u_obs_obssen(:,:), v_obs_nldepart(:,:), v_obs_tldepart(:,:), v_obs_obssen(:,:), &
       obs_nldepart(:,:), obs_tldepart(:,:), obs_obssen(:,:)
  logical passed

  integer, parameter :: nnln = 8
  character(len=*), parameter :: nlnqct(nnln)=(/ ' dw', ' ps', 'tcp', '  q', &
       'spd', '  t', ' uv', 'sst'/)
  integer, parameter :: noz = 15
  character(len=*), parameter :: oztype(noz)=(/'sbuv2', 'omi', 'mls',    & 
       'mls20', 'mls22', 'mls30', 'mls55', 'tomseff', 'omieff', 'o3lev', &
       'gome','ompslpuv','ompslpvis','ompsnm','ompsnp'/)

!     ODS variables
!     -------------
  integer nobs_ods, nobs, nobs1, nobs2
  integer ierr
  integer(i_kind), allocatable :: ivals(:),iuse(:)
  real(r_single), allocatable :: rvals(:),rvals2(:),sigo(:),o3pres(:)
  real(r_double), allocatable :: varchn(:)
  character(len=3)  obsclass
  character(len=:), dimension(:), allocatable :: ob_class
  character(len=:), dimension(:), allocatable :: station
  integer nymdh,nymdf,nhmsf
  integer iks

  tiny_single = tiny(zero_single)
  small_num = 10. * tiny_single
  ladjsigo = .true.
  verbose = .true.

  nymdh = 0

  rc = 0


! Inquire whether this is diag file w/ history of omx, osens, ...
! --------------------------------------------------------------- 
  call ods_obsdiags_getparam ( 'lobsdiagsave', lobsdiagsave )
  call ods_obsdiags_getparam ( 'lobssens', lobssens )
  call ods_obsdiags_getparam ( 'miter', miter )
  call ods_obsdiags_getparam ( 'ladjsigo', ladjsigo )

! Initialize the nc_diag reader
! -----------------------------

  status = nf90_open(path= trim(fname), mode= nf90_nowrite,  &
       ncid = input_id)
  if (status /= nf90_noerr) then
     rc=1
     print *,'Unable to read file ',trim(fname)
     return
  end if


! Get the date_time value and check against the 'requested'
! ---------------------------------------------------------
  status = nf90_inquire_attribute(input_id,NF90_GLOBAL,'date_time')
  if ( status /= nf90_noerr ) then
     print *,'Unable to find date/time in file ',trim(fname)
     rc = 1
     status = nf90_close(input_id)
     return
  end if
  status = nf90_get_att(input_id,NF90_GLOBAL,'date_time',nymdh)
  if ( status /= nf90_noerr ) then
     print *,'Unable to read date/time from file ',trim(fname)
     rc = 1
     status = nf90_close(input_id)
     return
  end if
  
  nymdf = nymdh / 100
  nhmsf = mod(nymdh,100)*10000
  if (nymdf/=nymd .or. nhmsf/=nhms) then
     rc = 1
     print*, 'Date/Time requested not found in nc4 file'
     print*, ' nymd , nhms =',nymd ,nhms
     print*, ' nymdf, nhmsf=',nymdf,nhmsf
     status = nf90_close(input_id)
     return
  end if

! Start reading the data from the file
! ------------------------------------  
  nobs = get_dim(input_id,'nobs',rc)
  if (rc /= 0) return
  
! Check for the observation type.
! -------------------------------
  status = nf90_inquire_attribute(input_id,NF90_GLOBAL,'Satellite')
  if (status == nf90_noerr) then
     print *,'found Satellite attribute'

! Processing radiance or ozone type
! ---------------------------------
     status = nf90_get_att(input_id,NF90_GLOBAL,'Observation_type', satype)
     if (status /= nf90_noerr) then
        print *,'Unable to read Observation_type metadata'
        rc = 1
        status = nf90_close(input_id)
        return
     end if
     status = nf90_get_att(input_id,NF90_GLOBAL,'Satellite_Sensor', isis)
     if (status /= nf90_noerr) then
        print *,'Unable to read Satellite_Sensor metadata'
        rc = 1
        status = nf90_close(input_id)
        return
     end if
     status = nf90_get_att(input_id,NF90_GLOBAL,'Satellite', dplat)
     if (status /= nf90_noerr) then
        print *,'Unable to read Satellite metadata'
        rc = 1
        status = nf90_close(input_id)
        return
     end if

     if( any(trim(satype) == oztype )) then
        diag_type = 'ozone'
        call ozone_getisat_(satype,dplat,isis,isat)
        print *,'ozone: ',satype,dplat,isis,isat
     else
        diag_type = 'radiance'
        call getsatid_(isat)
     endif
     nobs_ods = nobs
     print *,' diag type = ',diag_type, ' sensor = ',isis,' nobs = ',nobs_ods

     satknown = .false.
     kidsat = isat
     do n = 1,nsats
        if (trim(sats(n))==trim(satype)) then
           satknown = .true.
           kidsat = kidsat + idsats(n)
           exit
        else
           cycle
        endif
     enddo
     if (.not. satknown .or. kidsat == 0) then
        print *,'Cannot identify satellite type:'
        print *,'isis = ',isis,' dplat = ',dplat,' satype = ',satype,' isat =',isat
        status = nf90_close(input_id)
        rc=3
        return
     end if
        
  else if (  nf90_inq_dimid(input_id,'Observation_Class_maxstrlen',idimid) == nf90_noerr ) then

!
! We -ought- to be able to get the type from a global (assuming the diag file
! is homogeneous) but for now will take it from the first element of the 
! observation class array.  
!!TO_DO - add global variable with observation class to non-radiance diag files?


     status = nf90_inquire_dimension(input_id,idimid,len=clen)
     if (status /= nf90_noerr) then
        print *, 'problem getting obs class string length'
        rc=1
        status = nf90_close(input_id)
        return
     end if
        
     allocate(character(len=clen) :: ob_class(nobs))
     status = nf90_inq_varid(input_id,'Observation_Class',varid)
     if (status /= nf90_noerr) then
        print *,'error getting Observation Class'
        rc = 1
        status = nf90_close(input_id)
        return
     end if
     status = nf90_get_var(input_id, varid, ob_class)
     if (status /= nf90_noerr) then
        print *,'error getting Observation Class'
        rc = 1
        status = nf90_close(input_id)
        return
     end if
     obsclass = ob_class(1)(clen-2:clen)
     if ( obsclass == ' uv') then 
        diag_type = 'uvconv'
        nobs1 = nobs+1
        nobs2 = nobs*2
        nobs_ods = nobs2
     else 
        diag_type = 'conv'
        nobs_ods = nobs
     end if
     print *,' diag type = ',diag_type,'ob_class = ',ob_class(1),' nobs = ',nobs_ods
  else
     print *,'missing info, can''t convert obstype, ',trim(fname)
     rc = 4
     status = nf90_close(input_id)
     return
  end if

  call ODS_Init (ods, nobs_ods, ierr)
  if (ierr /= 0) then
     print *, 'error initializing ODS vector, ierr = ',ierr
     rc = 7
     status = nf90_close(input_id)
     return
  end if

  allocate(ivals(nobs),rvals(nobs),rvals2(nobs),sigo(nobs),stat=ierr)
  if (ierr /= 0) then
     print *, 'error alllocating arrays for reading data, ierr= ',ierr
     rc = 7
     status = nf90_close(input_id)
     return
  end if

!  Fields common to all data types - latitude, longitude, obs, omf
!  ---------------------------------------------------------------  
  call get_1d_var(input_id,'Latitude',rvals,rc)
  if (rc /= 0) return
  ods%data%lat(1:nobs) = rvals
  call get_1d_var(input_id,'Longitude',rvals,rc)
  if (rc /= 0) return
  where( rvals > 180)  rvals = rvals - 360.
  ods%data%lon(1:nobs) = rvals

  if (diag_type == 'uvconv') then
     call get_1d_var(input_id,'u_Observation',rvals,rc)
     if (rc /= 0) return
     ods%data%obs(1:nobs) = rvals
     call get_1d_var(input_id,'u_Obs_Minus_Forecast_adjusted',rvals,rc)
     if (rc /= 0) return
     ods%data%omf(1:nobs) = rvals
     call get_1d_var(input_id,'v_Observation',rvals,rc)
     if (rc /= 0) return     
     ods%data%obs(nobs1:nobs_ods) = rvals
     call get_1d_var(input_id,'v_Obs_Minus_Forecast_adjusted',rvals,rc)
     if (rc /= 0) return
     ods%data%omf(nobs1:nobs_ods) = rvals

     ods%data%lat(nobs1:nobs_ods) = ods%data%lat(1:nobs)
     ods%data%lon(nobs1:nobs_ods) = ods%data%lon(1:nobs)
  else
     call get_1d_var(input_id,'Observation',rvals,rc)
     if (rc /= 0) return
     ods%data%obs = rvals
     call get_1d_var(input_id,'Obs_Minus_Forecast_adjusted',rvals,rc)
     if (rc /= 0) return
     ods%data%omf = rvals
  end if
  ods%data%oma = undef
  ods%data%qchist = 0

! Processing for radiance types
! -----------------------------
  if ( diag_type == 'radiance' ) then

     nchan_dim = get_dim(input_id,'nchans',rc)
     if (rc /= 0) return

     call get_1d_var(input_id,'Obs_Time',rvals,rc)
     if (rc /= 0) return
     ods%data%time = int(rvals * 60.)            ! use 'int' to match diag_bin time
     call get_1d_var(input_id,'Channel_Index',ivals,rc)
     if (rc /= 0) return
     ods%data%lev = float(ivals)
     call get_1d_var(input_id,'Obs_Minus_Forecast_unadjusted',rvals,rc)
     if (rc /= 0) return
     ods%data%xm=rvals-ods%data%omf

     call get_1d_var(input_id,'Inverse_Observation_Error',rvals,rc)
     if (rc /= 0) return

     sigo = undef
     
! Check if original or adjusted (final) sigO value is requested
! -------------------------------------------------------------
     if (ladjsigo) then
        where (rvals > small_num)  sigo = 1.0/rvals
     else
        allocate(varchn(nchan_dim))
        call get_1d_var(input_id,'error_variance',varchn,rc)
        if (rc /= 0) return

        varchn = sqrt(varchn)
        where (rvals > small_num) sigo = varchn(ods%data%lev)
        deallocate(varchn)
     end if
     ods%data%Xvec = sigo       ! radiance Xvec
     
     allocate(iuse(nchan_dim))
     call get_1d_var(input_id,'use_flag',iuse,rc)
     if (rc /= 0) return
     ivals = 0
     where(rvals <= small_num) ivals = 2
     call get_1d_var(input_id,'QC_Flag',rvals,rc)
     if (rc /= 0) return
     where( rvals < 0 ) ivals  = 1
     do i = 1,nchan_dim
        if (iuse(i) == -1) where(ods%data%lev == i) ivals=1
     end do
     ods%data%qcexcl = ivals
     deallocate(iuse)
     
     nspot = nobs/nchan_dim
     knt = 0
     do i = 1,nspot
        do j = 1,nchan_dim
           knt = knt + 1
           ods%data%ks(knt) = i
        end do
     end do

     ods%data%kt = 40
     ods%data%kx = kidsat

! Processing for ozone types
! --------------------------
  else if ( diag_type == 'ozone') then

     ods%data%kx = kidsat

     call get_1d_var(input_id,'Time',rvals,rc)
     if (rc /= 0) return
     ods%data%time = int(rvals * 60.)

     call get_1d_var(input_id,'Reference_Pressure',rvals,rc)
     if (rc /= 0) return
     ods%data%lev = rvals
    

     select case(satype)
        
! Level ozone observation processing
! ----------------------------------
     case('o3lev','mls','mls20','mls22','mls30','mls55')
        ods%data%kt = kto3mx
        do i = 1,nobs
           ods%data%ks(i) = i
        end do

! Layer (and total) ozone observation processing
! ----------------------------------------------
     case default

        ods%data%kt = kto3
        iks = 1
        ods%data%xm(1) = 0.0
        do i = 1,nobs
           ods%data%ks(i) = iks
           if (ods%data%lev(i) == 0.0) then
              ods%data%kt(i) = kttco3             ! total ozone
              ods%data%lev(i) = undef
              ods%data%xm(i) = 0.0
              iks = iks + 1
           else
              ods%data%kt(i) = kto3               ! layer ozone
              if (i < nobs) ods%data%xm(i+1) = ods%data%lev(i)
           end if
        end do

     end select

     sigo = undef
     ivals = 0
     call get_1d_var(input_id,'Inverse_Observation_Error',rvals,rc)
     if (rc /= 0) return
     where(rvals > small_num)
        sigo = 1./rvals
     elsewhere
        ivals = 2
     endwhere
     ods%data%Xvec = sigo       ! ozone Xvec
     ods%data%qcexcl = ivals



  else if ( diag_type == 'uvconv' .or. diag_type == 'conv' ) then

! Process conventional data types (including GPS)
! ------------------------------------------------
     stlen = get_dim(input_id, 'Station_ID_maxstrlen',rc)
     if (rc /= 0) return
     allocate(character(len=stlen) :: station(nobs))
     call get_1d_var(input_id,'Station_ID',station,rc)
     if (rc /= 0) return

     call get_1d_var(input_id,'Time',rvals,rc)
     if (rc /= 0) return
     ods%data%time(1:nobs) = int(rvals * 60,)

     call get_1d_var(input_id,'Observation_Type',ivals,rc)
     if (rc /= 0) return
     ods%data%kx(1:nobs) = ivals

     call get_1d_var(input_id,'Pressure',rvals,rc)
     if (rc /= 0) return
     ods%data%lev(1:nobs) = rvals
 
     ivals = 0
     call get_1d_var(input_id, 'Analysis_Use_Flag',rvals,rc)
     if (rc /= 0) return
     where( rvals < 0 ) ivals = X_PASSIVE

     sigo = undef
     call get_1d_var(input_id,'Errinv_Final',rvals,rc)
     if (rc /= 0) return
     call get_1d_var(input_id,'Nonlinear_QC_Rel_Wgt',rvals2,rc)
     if (rc /= 0) return
!
! use nonlin QC mark where rel weight < 1  (overrides X_PASSIVE)
! --------------------------------------------------------------
     if( any(obsclass == nlnqct) ) where(rvals2 < 1 )  ivals = X_NCEP_NLNQC

! Check if original or adjusted (final) sigO value is requested
! -------------------------------------------------------------
     if (ladjsigo) then                ! use adjusted sigo
        rvals2 = rvals
     else                              ! use input sigo
        call get_1d_var(input_id,'Errinv_Input',rvals2,rc)
        if (rc /= 0) return
     end if

     where (rvals > small_num)
        sigo = 1./rvals2
     elsewhere
        sigo = undef
        ivals = 2                        ! rejected by QC
     endwhere
     ods%data%qcexcl(1:nobs) = ivals
     ods%data%Xvec(1:nobs) = sigo           !  conventional Xvec

     if (diag_type == 'uvconv') then
        ods%data%time(nobs1:nobs_ods)   = ods%data%time(1:nobs)
        ods%data%kx(nobs1:nobs_ods)     = ods%data%kx(1:nobs)
        ods%data%lev(nobs1:nobs_ods)    = ods%data%lev(1:nobs)
        ods%data%qcexcl(nobs1:nobs_ods) = ods%data%qcexcl(1:nobs)
        ods%data%Xvec(nobs1:nobs_ods)   = ods%data%Xvec(1:nobs)    ! v-component Xvec
     end if

     select case( obsclass )

     case ('gps') 
!!TO_DO - need GPS rdiag(19) 'hob' vertical grid location, to use 
!!        for adjusting error with pressure values == 0.0
        call get_1d_var(input_id,'GPS_Type',rvals,rc)
        if (rc /= 0) return
        
        if(rvals(1) == 0) then               ! assuming homogeneous file
           ods%data%kt= ktGPSr
           call get_1d_var(input_id,'Model_Elevation',rvals,rc)
           if (rc /= 0) return
           ods%data%xm = rvals
        else
           ods%data%kt= ktGPSb
           call get_1d_var(input_id,'Height',rvals,rc)
           if (rc /= 0) return
           ods%data%xm = rvals
        endif

!  fix for pressure = 0.
        if (any(ods%data%lev == 0.0)) then
           status = nf90_inq_varid(input_id,'Vertical_Grid_Location',varid)
           if (status == nf90_noerr) then
              call get_1d_var(input_id,'Vertical_Grid_Location',rvals,rc)
              if (rc /= 0) return
              do i = 1,nobs
                 if (ods%data%lev(i) == 0.0) then
                    if (rvals(i) < 1.) then
                       ods%data%lev(i) = 1050.       ! obs below model sfc
                    else
                       ods%data%lev(i) = 0.01        ! obs above model top
                    endif
                 endif
              end do
           end if
        end if
        
     case (' ps', 'tcp')  ! surface pressure
        ods%data%kt= ktps2m
        ods%data%lev = undef
        call get_1d_var(input_id,'Height',rvals,rc)
        if (rc /= 0) return
        ods%data%xm = rvals

     case('  q')
        ods%data%kt= ktqq        ! specific humidity
        ods%data%obs = 1.e3*ods%data%obs
        ods%data%omf = 1.e3*ods%data%omf
        call get_1d_var(input_id,'Forecast_Saturation_Spec_Hum',rvals,rc)
        if (rc /= 0) return
        ods%data%xm = 1.e3*rvals
        where(sigo /= undef) ods%data%Xvec = 1.e3*sigo      ! adjust Q Xvec
        
     case ('spd')
        ods%data%kt= ktus10       ! define it as 10m speeds
        ods%data%lev = undef
        call get_1d_var(input_id,'Height',rvals,rc)
        if (rc /= 0) return
        ods%data%xm = rvals

     case('sst')
        ods%data%kt= ktSST
        ods%data%lev = undef
! (use Station_Elevation to match diag_bin processing)
        call get_1d_var(input_id,'Station_Elevation',rvals,rc)
        if (rc /= 0) return
        ods%data%xm = rvals
        
     case ('  t')    ! virtual temperature
        ods%data%kt= ktTv
        where(ods%data%kx==311) ods%data%kx=304
        call get_1d_var(input_id,'Obs_Minus_Forecast_unadjusted',rvals,rc)
        if (rc /= 0) return        
        ods%data%xm=rvals-ods%data%omf   ! fill in bias correction as xm

     case(' uv')  ! vector wind
        ods%data%kt(1:nobs)          = ktuu
        ods%data%kt(nobs1:nobs_ods)  = ktvv
        call get_1d_var(input_id,'Height',rvals,rc)
        if (rc /= 0) return        
        ods%data%xm(1:nobs)         = rvals
        ods%data%xm(nobs1:nobs_ods) = rvals

     case(' pw')    ! total column water
        ods%data%kt      = ktTPW
        call get_1d_var(input_id,'Prep_QC_Mark',rvals,rc)
        if (rc /= 0) return       
        ods%data%xm = rvals
        
     case default
        print *,'can''t handle var = ',trim(ob_class(1)), nobs, ' observations'
        rc = 4
        status = nf90_close(input_id)
        return
     end select

!    Set sounding index as in diag_bin processing
!    --------------------------------------------
     
     call setsndx (ods%data%ks(1:nobs),ods%data%kx(1:nobs),station(1:nobs))

!    Set ks to actual station id for radiosondes
!    -------------------------------------------
     do i = 1, nobs
        if (ods%data%kx(i)==120 .OR. ods%data%kx(i)==220  .or.      &  ! radiosondes
             ods%data%kx(i)==901 .OR. ods%data%kx(i)==902) then  ! lagragian ballon data

!          Normally the station id is a 5-digit number ...
!          -------------------------------------------
           statid = station(i)(1:5)
           read(statid,'(i5)',err=999) statks
           if (statks > 0 ) then
              ods%data%ks(i) = statks
           else         !  if invalid numeric id, treat like non-numeric id
                        !  RT: code will never be here!
              do ii =1,len(station(i))
                 ods%data%ks(i) = ods%data%ks(i)+ii*iachar(station(i)(ii:ii)) !  arbitrary scheme to use a number better than 99999
              enddo
              print*, 'ods_diagnc4: YES CODE PASSES HERE SOMETIMES!'
           endif
           cycle

!          ... but if not, make sure this ks is distinct from any other station id
!          -----------------------------------------------------------------------
999        continue

           do ii =1,len(station(i))
              ods%data%ks(i) = ods%data%ks(i)+ii*iachar(station(i)(ii:ii)) !  arbitrary scheme to use a number better than 99999
           enddo
           if (verbose) print *,'Non-numeric station id ', statid, ' for kx = ', ods%data%kx(i), 'assigned ks =', ods%data%ks(i)

        end if
     end do
     
     deallocate(station)

     deallocate(ob_class)

  end if

  if ( diag_type == 'uvconv' ) then
     ods%data%ks(nobs1:nobs_ods)      = ods%data%ks(1:nobs)
  end if

  call getodsmeta( ods )

! The initial implementation of sensitivity calculation attempts to put the
! sensitivity data into an array arranged in the same way as the data in the
! binary files (except of course with only the sensitivity information) and
! calls the same m_ods_obsdiags routines in the same way to produce the same
! result (hopefully) as for the binary diag files.  Later on more efficient
! routines to process the data can be worked out.

  if (lobsdiagsave) then

     iiuse = get_dim(input_id,'ObsDiagSave_iuse_arr_dim',rc)
     if (rc /= 0) then
        print *,'problem getting iuse dimension, exiting'
        return
     end if

     if (iiuse /= miter) then
        print *,'iuse array inconsistent with miter, exiting'
        rc = 10
        return
     end if
     
     if (diag_type == 'uvconv') then
        iunldpt = get_dim(input_id,'u_ObsDiagSave_nldepart_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting nldepart dimension, exiting'
           return
        end if
        iutldpt = get_dim(input_id,'u_ObsDiagSave_tldepart_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting tldepart dimension, exiting'
           return
        end if
        iuobssen = get_dim(input_id,'u_ObsDiagSave_obssen_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting obssen dimension, exiting'
           return
        end if
        
        ivnldpt = get_dim(input_id,'v_ObsDiagSave_nldepart_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting nldepart dimension, exiting'
           return
        end if
        ivtldpt = get_dim(input_id,'v_ObsDiagSave_tldepart_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting tldepart dimension, exiting'
           return
        end if
        ivobssen = get_dim(input_id,'v_ObsDiagSave_obssen_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting obssen dimension, exiting'
           return
        end if
        
! allocate arrays for reading in sensitivity information
        allocate(obs_iuse(iiuse,nobs),u_obs_nldepart(iunldpt,nobs),   &
             u_obs_tldepart(iutldpt,nobs),u_obs_obssen(iuobssen,nobs), &
             v_obs_nldepart(ivnldpt,nobs),v_obs_tldepart(ivtldpt,nobs),&
             v_obs_obssen(ivobssen,nobs),pdata(7*miter+2,nobs),stat=ierr)
        if (ierr /= 0) then
           print *, 'error alllocating arrays for reading data, ierr= ',ierr
           rc = 7
           status = nf90_close(input_id)
           return
        end if
! read arrays with sensitivity information
        call get_2d_int_var(input_id,'ObsDiagSave_iuse',obs_iuse,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'u_ObsDiagSave_nldepart',u_obs_nldepart,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'u_ObsDiagSave_tldepart',u_obs_tldepart,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'u_ObsDiagSave_obssen',u_obs_obssen,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'v_ObsDiagSave_nldepart',v_obs_nldepart,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'v_ObsDiagSave_tldepart',v_obs_tldepart,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'v_ObsDiagSave_obssen',v_obs_obssen,rc)
        if (rc /= 0) return

        ioff = 0
        do jj = 1,miter
           ioff = ioff + 1
           pdata(ioff,1:nobs) = obs_iuse(jj,1:nobs)
        enddo
        do jj = 1,miter+1
           ioff = ioff + 1
           pdata(ioff,1:nobs) = u_obs_nldepart(jj,1:nobs)
           ioff = ioff + 1
           pdata(ioff,1:nobs) = v_obs_nldepart(jj,1:nobs)
        end do
        do jj = 1,miter
           ioff = ioff + 1
           pdata(ioff,1:nobs) = u_obs_tldepart(jj,1:nobs)
           ioff = ioff + 1
           pdata(ioff,1:nobs) = v_obs_tldepart(jj,1:nobs)
        end do
        do jj = 1,miter
           ioff = ioff + 1
           pdata(ioff,1:nobs) = u_obs_obssen(jj,1:nobs)
           ioff = ioff + 1
           pdata(ioff,1:nobs) = v_obs_obssen(jj,1:nobs)
        end do

        do i = 1,nobs
           ods%data%Xvec(i) = 0.0
           ods%data%Xvec(i+nobs) = 0.0
           call ods_obsdiags(nlomx, tlomx, obimp, pdata, 0, i,  &
                7*miter+2, nobs, undef, passed)
           if (passed) then
              ods%data%Xvec(i) = obimp(1)
              ods%data%Xvec(i+nobs) = obimp(2)
              ods%data%qcexcl(i) = 0
              ods%data%qcexcl(i+nobs) = 0
           end if
        end do
        
        deallocate(obs_iuse, u_obs_nldepart, u_obs_tldepart,  &
             u_obs_obssen, v_obs_nldepart, v_obs_tldepart,    &
             v_obs_obssen, pdata)

     else
! get dimensions for sensitivity variables
        inldpt = get_dim(input_id,'ObsDiagSave_nldepart_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting nldepart dimension, exiting'
           return
        end if
        itldpt = get_dim(input_id,'ObsDiagSave_tldepart_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting tldepart dimension, exiting'
           return
        end if
        iobssen = get_dim(input_id,'ObsDiagSave_obssen_arr_dim',rc)
        if (rc /= 0) then
           print *,'problem getting obssen dimension, exiting'
           return
        end if

! allocate arrays for reading in sensitivity information
        allocate(obs_iuse(iiuse,nobs),obs_nldepart(inldpt,nobs),   &
             obs_tldepart(itldpt,nobs),obs_obssen(iobssen,nobs),   &
             pdata(4*miter+1,nobs),stat=ierr)
        if (ierr /= 0) then
           print *, 'error alllocating arrays for reading data, ierr= ',ierr
           rc = 7
           status = nf90_close(input_id)
           return
        end if
! read arrays with sensitivity information
        call get_2d_int_var(input_id,'ObsDiagSave_iuse',obs_iuse,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'ObsDiagSave_nldepart',obs_nldepart,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'ObsDiagSave_tldepart',obs_tldepart,rc)
        if (rc /= 0) return
        call get_2d_double_var(input_id,'ObsDiagSave_obssen',obs_obssen,rc)
        if (rc /= 0) return

        ioff = 0
        do jj = 1,miter
           ioff = ioff + 1
           pdata(ioff,1:nobs) = obs_iuse(jj,1:nobs)
        enddo
        do jj = 1,miter+1
           ioff = ioff + 1
           pdata(ioff,1:nobs) = obs_nldepart(jj,1:nobs)
        end do
        do jj = 1,miter
           ioff = ioff + 1
           pdata(ioff,1:nobs) = obs_tldepart(jj,1:nobs)
        end do
        do jj = 1,miter
           ioff = ioff + 1
           pdata(ioff,1:nobs) = obs_obssen(jj,1:nobs)
        end do

        do i = 1,nobs
           ods%data%Xvec(i) = 0.0
           call ods_obsdiags(nlomx(1), tlomx(1), obimp(1), pdata, 0, i,  &
                4*miter+1, nobs, undef, passed)
           if (passed) then
              ods%data%Xvec(i) = obimp(1)
              ods%data%qcexcl(i) = 0
           end if
        end do        
        
        deallocate(obs_iuse, obs_nldepart, obs_tldepart, obs_obssen, pdata)

     end if
  end if
  
  status = nf90_close(input_id)

  deallocate(ivals, rvals, rvals2, sigo)

  return

contains
  
  subroutine getsatid_(myisat)
! RT: some heck that needs more work
! RT: this needs serious attention as it is becoming a huge heck now (3/30/09)
    implicit none
    integer, intent(out) :: myisat
    integer ios
    
    myisat     = 0  ! take fixed sat index as in idsats
    
!           select case( trim(ladjust(dplat)) )
!           case ('aura')
!              myisat = 999
!             return
!           end select

! first handle precip types
    i = index('pcp',isis(1:3))
    if (i>0) then
       i = index('trmm',dplat(1:4))
       if (i>0) then
          if (dplat(6:8) == 'lnd') then
             myisat = 1
          else if (dplat(6:8) == 'ocn') then
             myisat = 2
          else
             myisat = 0
          end if
          return
       else
          read(dplat(5:6),'(i2)',iostat=ios)myisat
          return
       end if
    end if

! Need to distinguish between AMSUA from AQUA and METOP for example
! (NOAA sats are already distinguished)
    i = index('metop',dplat(1:5))
    if(i>0)then
       myisat = 25 + iachar(dplat(7:7)) - iachar('a')
       return
    endif
    i = index('tiros',dplat(1:5))
    if(i>0)then
       myisat = 5
       return
    endif
!           i = index('aqua',dplat(1:4))
!           if(i>0)then
!              return
!           endif
! the "word" fgnm stands for the platforms dmsp/goes/noaa/meteosat this needs generalization
! e.g. dmsp -> f15   goes -> g12   noaa -> n18  meteosat -> m09
    i = index('fgnm',dplat(1:1))
    if(i>0)then
       read(dplat(2:3),'(i2)',iostat=ios)myisat
    endif
       
  end subroutine getsatid_

  subroutine ozone_getisat_(dtype,dplat,dsis,myisat)
    implicit none
    character(len=*),intent(in):: dtype,dplat,dsis
    integer,intent(out):: myisat
    
    myisat=0        ! for a default isat value.
    
    if(len_trim(dplat)==3 .and.                            &
         verify(dplat(1:1),'fgnm')==0 .and.                &
         verify(dplat(2:3),'0123456789')==0) then
! the "word" fgnm stands for the platforms dmsp/goes/noaa/meteosat this needs generalization
! e.g. dmsp -> f15   goes -> g12   noaa -> n18  meteosat -> m09
       read(dplat(2:3),'(i2)') myisat
       
    else
       select case(dplat)
       case('aura')
          myisat = 0
       case('nim07')
          myisat = 1
       case('ep')
          myisat = 2
       case('metop-a')
          myisat = 25
       case('metop-b')
          myisat = 26
       case default
          if(dplat(1:6)=='metop-') myisat = 25 + iachar(dplat(7:7)) - iachar('a')
       end select
    endif
  end subroutine ozone_getisat_

end subroutine ods_diagnc4
