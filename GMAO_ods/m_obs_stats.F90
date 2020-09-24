module m_obs_stats

use m_die, only: die

implicit none

private

public :: obs_stats_init
public :: obs_stats_collect
public :: obs_stats_final
public :: obs_stats_write

public :: obs_stats

! All internal from here down

character(len=*), parameter :: myname = "obs_stats"

integer, parameter :: mychar = 30

type obs_stats

  character(len=mychar) :: expid
  character(len=mychar) :: forecast
  character(len=mychar) :: source
  character(len=mychar) :: verification

  integer :: ns
  integer, pointer ::               counter(:) => null()
  integer, pointer ::                nymdhh(:) => null()
  character(len=mychar), pointer :: domain(:)  => null()   ! global,tropics,n.hem,etc
  real, pointer ::               boundary(:,:) => null()   ! end points of domain east,north,south,west
  real, pointer ::                    level(:) => null()   ! level/channel
  character(len=mychar), pointer :: levtype(:) => null()   ! level or channel: pl or ch
  character(len=mychar), pointer :: varname(:) => null()   ! variable name: tv, u, tb, bend, etc
  integer, pointer ::                 step(:)  => null()   ! hours: 0 from analysis, but could be 24 if omf from 24 hr fcsts
  character(len=mychar), pointer ::  sname(:)  => null()   ! rms,efficiency,etc
  real, pointer ::                   svalue(:) => null()   ! actual value of statistic

end type obs_stats

interface obs_stats_init
   module procedure init_
end interface
interface obs_stats_collect
   module procedure collect_
end interface
interface obs_stats_final
   module procedure final_
end interface
interface obs_stats_write
   module procedure write_
end interface

contains

  subroutine init_ ( ns, stats, expid, fcst, source, verif )

  implicit none

  integer, intent(in)            :: ns ! number of stats entries
  type(obs_stats), intent(inout) :: stats
  character(len=*),intent(in)    :: expid
  character(len=*),intent(in)    :: fcst
  character(len=*),intent(in)    :: source
  character(len=*),intent(in)    :: verif

  stats%expid = expid 
  stats%forecast = fcst 
  stats%source = source 
  stats%verification = verif 

  stats%ns=ns
  allocate(stats%counter(ns))
  allocate(stats%nymdhh(ns))
  allocate(stats%domain(ns))
  allocate(stats%boundary(ns,4))
  allocate(stats%level(ns))
  allocate(stats%levtype(ns))
  allocate(stats%varname(ns))
  allocate(stats%step(ns))
  allocate(stats%sname(ns))
  allocate(stats%svalue(ns))

  end subroutine init_

  subroutine final_ ( stats )

  implicit none

  type(obs_stats), intent(inout) :: stats

  stats%expid = "null" 
  stats%forecast = "null"
  stats%source = "null"
  stats%verification = "null"

  deallocate(stats%svalue)
  deallocate(stats%sname)
  deallocate(stats%step)
  deallocate(stats%varname)
  deallocate(stats%levtype)
  deallocate(stats%level)
  deallocate(stats%boundary)
  deallocate(stats%domain)
  deallocate(stats%nymdhh)
  deallocate(stats%counter)

  end subroutine final_

  subroutine collect_(nymd,nhms,hour,varname,sfctype,sname,levtype,level, &
                      regnames,regbounds,&
                      svalue,stats, &
                      fakevars)
  implicit none

  integer,intent(in) :: nymd, nhms
  integer,intent(in) :: hour
  character(len=*),intent(in) :: varname(:)
  character(len=*),intent(in) :: sname(:)
  character(len=*),intent(in) :: levtype(:)
  character(len=*),intent(in) :: regnames(:)
  logical,intent(in) :: sfctype(:)
  real,intent(in) :: level(:)
  real,intent(in) :: regbounds(:,:)
  real,intent(in) :: svalue(:,:,:,:)

  character(len=*),intent(in), optional :: fakevars(:)

  character(len=*), parameter :: myname_ = myname//"*collect_"
  type(obs_stats), intent(inout) :: stats

  integer nlevs,kk,nv,ir,ns,ic,nhmdhh,nfake

!count|date|domain_name|east|expver|forecast|level|levtype|north|variable|source|south|statistic|step|type|value|verify|west
!0.0|2020031900|global|180.0000|sapoes00rtXEC.21z|gmao|1000.000|pl|90.00000|p|gmao|-90.00000|cor|0|fc|0.9695317|gmao|-180.0000
! fixed for now:
  stats%counter = 0
  nfake = 0
  if(present(fakevars)) then
     nfake = size(fakevars)
  endif

  if(size(svalue,1)/=size(regnames)) then
      call die (myname_,': inconsistent region dims, aborting ', 99)
  endif
  if(size(svalue,3)/=size(varname)) then
      call die (myname_,': inconsistent varname dims, aborting ', 99)
  endif
  if(size(svalue,4)/=size(sname)) then
      call die (myname_,': inconsistent statistic dims, aborting ', 99)
  endif
  if(size(sfctype)/=size(varname)) then
      call die (myname_,': inconsistent sfc/varname dims, aborting ', 99)
  endif

! from arg list:
  stats%nymdhh    = nymd*100 + nhms/10000

  nlevs = size(svalue,2)
  ic=0
  do ns = 1,size(svalue,4)           ! number of statistic types
     do nv = 1,size(svalue,3)        ! variables
        do ir = 1,size(svalue,1)     ! regions
           do kk = 1,nlevs           ! levels
              if(sfctype(nv).and.kk>1) cycle
              ic = ic + 1
              if ( ic > stats%ns ) then
                  call die (myname_,': inconsistent dims, aborting ', 99)
              endif

              stats%step(ic)    = hour
              stats%sname(ic)   = sname(ns)
              stats%varname(ic) = varname(nv)
              stats%svalue(ic)  = svalue(ir,kk,nv,ns)
              stats%level(ic)   = level(kk)
              stats%levtype(ic) = levtype(1)       ! to be revisited
              stats%domain(ic)  = regnames(ir)     ! to be revisited
              stats%boundary(ic,1) = 180.0         ! east
              stats%boundary(ic,2) =-180.0         ! west
              stats%boundary(ic,3) = regbounds(2,ir) ! north
              stats%boundary(ic,4) = regbounds(1,ir) ! souhh

           enddo
        enddo
     enddo
  enddo
  if (nfake>0) then ! place h as fake variable to satisfy wired scorecard plotting program
     nlevs = size(svalue,2)
     do ns = 1,size(svalue,4)           ! number of statistic types
        do nv = 1,nfake                 ! number of fake variables
           do ir = 1,size(svalue,1)     ! regions
              do kk = 1,nlevs           ! levels
                 ic = ic + 1
                 if ( ic > stats%ns ) then
                     call die (myname_,': inconsistent dims, aborting ', 99)
                 endif

                 stats%step(ic)    = hour
                 stats%sname(ic)   = sname(ns)
                 stats%varname(ic) = fakevars(nv)
                 stats%svalue(ic)  = 1e15
                 stats%level(ic)   = level(kk)
                 stats%levtype(ic) = levtype(1)
                 stats%domain(ic)  = regnames(ir)
                 stats%boundary(ic,1) = 180.0
                 stats%boundary(ic,2) =-180.0
                 stats%boundary(ic,3) = regbounds(2,ir)
                 stats%boundary(ic,4) = regbounds(1,ir)
              enddo
           enddo
        enddo
     enddo
  endif
  
  end subroutine collect_
 
  subroutine write_(outfile,stats)

  use m_ioutil, only: luavail
  implicit none
  
  character(len=*) outfile
  type(obs_stats), intent(in) :: stats

  character(len=256)  record
  character(len=mychar)  cdate,cstat,czlev,cstep,ceast,cwest,cnorth,csouth
  integer iii,luout

  luout=luavail()
  open(luout,file=trim(outfile),form='formatted')
                       ! the following is in alphabetic order
  write(luout,'(5a)') 'count|date|domain_name|east|',   &
                      'expver|forecast|level|levtype|', &
                      'north|variable|source|south|',   &
                      'statistic|step|type|value|',     &
                      'verify|west'
   do iii=1,stats%ns
!     if (abs(stats%svalue(iii))>1.e14) cycle
      write( cstat,*) stats%svalue(iii)
!     write( cdate,'(i10.10)') stats%nymdhh(iii)
      write( cdate,*) stats%nymdhh(iii)
      write( czlev,*) stats%level(iii)
      write( cstep,*) stats%step(iii)
      write( ceast,*) stats%boundary(iii,1)
      write( cwest,*) stats%boundary(iii,2)
      write(cnorth,*) stats%boundary(iii,3)
      write(csouth,*) stats%boundary(iii,4)
#ifdef _DEBUG_
print*,                 trim(adjustl(cdate))
print*,     trim(adjustl(stats%domain(iii)))
print*,                 trim(adjustl(ceast))
print*,           trim(adjustl(stats%expid))
print*,          trim(adjustl(stats%source))
print*,                 trim(adjustl(czlev))
print*,    trim(adjustl(stats%levtype(iii)))
print*,                trim(adjustl(cnorth))
print*,    trim(adjustl(stats%varname(iii)))
print*,          trim(adjustl(stats%source))
print*,                trim(adjustl(csouth))
print*,      trim(adjustl(stats%sname(iii)))
print*,                 trim(adjustl(cstep))
print*,                                 'fc'
print*,                 trim(adjustl(cstat))
print*,    trim(adjustl(stats%verification))
print*,                 trim(adjustl(cwest))          ! west
#endif
      write(record,'(f3.1,34a)') 0.0, '|',  &  ! count (dummy for this code)
                 trim(adjustl(cdate)),'|',  &  ! date
     trim(adjustl(stats%domain(iii))),'|',  &  ! domain_name
                 trim(adjustl(ceast)),'|',  &  ! east
           trim(adjustl(stats%expid)),'|',  &  ! expver
          trim(adjustl(stats%source)),'|',  &  ! forecast
                 trim(adjustl(czlev)),'|',  &  ! level
    trim(adjustl(stats%levtype(iii))),'|',  &  ! levtype
                trim(adjustl(cnorth)),'|',  &  ! north
    trim(adjustl(stats%varname(iii))),'|',  &  ! variable
          trim(adjustl(stats%source)),'|',  &  ! source=forecast
                trim(adjustl(csouth)),'|',  &  ! south
      trim(adjustl(stats%sname(iii))),'|',  &  ! statistic
                 trim(adjustl(cstep)),'|',  &  ! step
                                 'fc','|',  &  ! type
                 trim(adjustl(cstat)),'|',  &  ! value
    trim(adjustl(stats%verification)),'|',  &  ! verify
                 trim(adjustl(cwest))          ! west
      write(luout,'(a)') trim(adjustl(record))

  enddo
  close(luout)

  end subroutine write_


end module m_obs_stats
