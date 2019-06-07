!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! NASA/GSFC, Global Modeling and Assimilation Office, 900.3, GEOS/DAS  !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_checksums - Show local check sum values
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_checksums
      implicit none
      private	! except

      public :: checksums_show		! data structure

    interface checksums_show; module procedure	&
      showlv_,	&
      showgd_; end interface

! !REVISION HISTORY:
! 	04Feb05	- Jing Guo <guo@gmao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_checksums'
  logical,parameter :: CHECKSPV=.true.
#include "assert.H"
contains
subroutine showgd_(a,atype,header,lu,ref,inc,undef)
  use m_undef,only : undef_ssi
  use m_stdio,only : stdout
  use m_die,only : assert_
  implicit none
  real,dimension(:,:,:),intent(in) :: a
  character(len=*),intent(in) :: atype
  character(len=*),intent(in) :: header
  integer,optional,intent(in) :: lu
  real   ,optional,dimension(:),intent(in) :: ref
  integer,optional,intent(in) :: inc
  real   ,optional,intent(in) :: undef

  integer :: lu_,inc_
  real,dimension(size(a,3)) :: ref_
  real :: undef_
  integer :: k,mlev

  lu_=stdout
  if(present(lu)) lu_=lu

  mlev=size(a,3)
  do k=1,mlev
    ref_(k)=k
  end do

  inc_=1
  if(present(inc)) inc_=inc

  undef_=undef_ssi()
  if(present(undef)) undef_=undef

  if(present(ref)) then
    ASSERT(size(ref_)==size(ref))
    ref_=ref
  endif
  call show3d(lu_,size(a,1),size(a,2),size(a,3),	&
	a,ref_,atype,'INDX', undef_,header,inc_)
end subroutine showgd_

subroutine showlv_(a,atype,flag,lu,ref,undef)
  use m_stdio,only : stdout
  use m_undef,only : undef_ssi
  implicit none
  real,dimension(:,:),intent(in) :: a
  character(len=*),intent(in) :: atype
  character(len=*),intent(in) :: flag
  integer,optional,intent(in) :: lu
  real   ,optional,intent(in) :: ref
  real   ,optional,intent(in) :: undef

  integer :: lu_
  real :: ref_
  real :: undef_
  lu_=stdout
  if(present(lu)) lu_=lu
  ref_=0.
  if(present(ref)) ref_=ref
  undef_=undef_ssi()
  if(present(undef)) undef_=undef
  call showlv(lu_,size(a,1),size(a,2),a,ref_,atype,'INDX', &
	undef_,flag)
end subroutine showlv_


	subroutine showlv (lu,mx,my,a,h,atype,htype,amiss,flag)

!	..Print statistics of one 2-d variable
	implicit none

	integer,intent(in) :: lu		! Output unit
	integer,intent(in) :: mx,my		! Array sizes
	real,intent(in) :: a(mx,my)		! The array
	real,intent(in) :: h			! The argument
	character(len=*),intent(in) :: atype	! Type of the variable(array)
	character(len=*),intent(in) :: htype	! Type of the level(argument)
	real,intent(in) :: amiss		! missing value flag of a
	character(len=*),intent(in) :: flag

	integer i,j
	integer imx,imn,jmx,jmn
	integer knt
	real amx,amn
	real avg,dev,d
	logical first

!	..A practical value for the magnitude of the fraction of a real
!	number.

	real rfrcval
	parameter(rfrcval=1.e-5)

!	..function

	logical spv
	real aspv
	spv(aspv)=abs((aspv-amiss)/amiss).le.rfrcval .and. CHECKSPV

	knt=0
	avg=0.
	do j=1,my
	  do i=1,mx
	    if(.not.spv(a(i,j))) then
	      knt=knt+1
	      avg=avg+(a(i,j)-avg)/knt
	    endif
	  end do
	end do

	knt=0
	dev=0.
	do j=1,my
	  do i=1,mx
	    if(.not.spv(a(i,j))) then
	      d=a(i,j)-avg
	      knt=knt+1
	      dev=dev+(d*d-dev)/knt
	    endif
	  end do
	end do
	dev=sqrt(dev*knt/max(1,knt-1))

	amx=a(1,1)
	amn=a(1,1)
	first=.true.
	do j=1,my
	  do i=1,mx
	    if(.not.spv(a(i,j))) then
	      if(first) then
		imx=i
		imn=i
		jmx=j
		jmn=j
		amx=a(imx,jmx)
		amn=a(imn,jmn)
		first=.false.
	      else
		if(a(i,j).gt.amx) then
		  amx=a(i,j)
		  imx=i
		  jmx=j
		endif
		if(a(i,j).lt.amn) then
		  amn=a(i,j)
		  imn=i
		  jmn=j
		endif
	      endif
	    endif
	  end do
	end do

	if(atype.eq.'RELH') then
	  avg=avg*100.
	  dev=dev*100.
	  amx=amx*100.
	  amn=amn*100.
	endif

	if(htype.eq.'NULL'.or.htype.eq.'SRFC') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')	&
     &	      flag,':',knt,nint(avg),'+',nint(dev),	&
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',	&
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.	&
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.	&
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,i6,f7.1,a,f6.1,2(a,f7.1,a,i3,a,i3,a))')	&
     &	      flag,':',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,i6,f7.3,a,f6.3,2(a,f7.3,a,i3,a,i3,a))')	&
     &	      flag,':',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,i6,e10.3,a,e9.3,'//	&
     &	      '2(a,e9.3,a,i3,a,i3,a))')	&
     &	      flag,':',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	elseif(htype.eq.'PRES'.or.htype.eq.'INDX') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,i4,a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,nint(avg),'+',nint(dev),	&
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',	&
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.	&
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.	&
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,i4,a,i6,f7.1,a,f6.1,'//	&
     &	      '2(a,f7.1,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,i4,a,i6,f7.3,a,f6.3,'//	&
     &	      '2(a,f7.3,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,i4,a,i6,e10.3,a,e9.3,'//	&
     &	      '2(a,e10.3,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	elseif(htype.eq.'HGHT') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,i5,a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,nint(avg),'+',nint(dev),	&
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',	&
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.	&
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.	&
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,i5,a,i6,f7.1,a,f6.1,'//	&
     &	      '2(a,f7.1,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,i5,a,i6,f7.3,a,f6.3,'//	&
     &	      '2(a,f7.3,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,i5,a,i6,e10.3,a,e9.3,'//	&
     &	      '2(a,e10.3,a,i3,a,i3,a))')	&
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	elseif(htype.eq.'TEMP') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,f7.2,a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')	&
     &	      flag,'(',h,')',knt,nint(avg),'+',nint(dev),	&
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',	&
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.	&
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.	&
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,f7.2,a,i6,f7.1,a,f6.1,'//	&
     &	      '2(a,f7.1,a,i3,a,i3,a))')	&
     &	      flag,'(',h,')',knt,avg,'+',dev,	&
     &	      ' mx=',amx,'(',imx,',',jmx,')',	&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,f7.2,a,i6,f7.3,a,f6.3,'//		&
     &	      '2(a,f7.3,a,i3,a,i3,a))')				&
     &	      flag,'(',h,')',knt,avg,'+',dev,			&
     &	      ' mx=',amx,'(',imx,',',jmx,')',			&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,f7.2,a,i6,e10.3,a,e9.3,'//	&
     &	      '2(a,e10.3,a,i3,a,i3,a))')			&
     &	      flag,'(',h,')',knt,avg,'+',dev,			&
     &	      ' mx=',amx,'(',imx,',',jmx,')',			&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	else
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,e10.3,a,i6,i7,a,i6,'//	&
     &	      '2(a,i7,a,i3,a,i3,a))')				&
     &	      flag,'(',h,')',knt,nint(avg),'+',nint(dev),	&
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',		&
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.		&
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.		&
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,e10.3,a,i6,f7.1,a,f6.1,'//	&
     &	      '2(a,f7.1,a,i3,a,i3,a))')				&
     &	      flag,'(',h,')',knt,avg,'+',dev,			&
     &	      ' mx=',amx,'(',imx,',',jmx,')',			&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,e10.3,a,i6,f7.3,a,f6.3,'//	&
     &	      '2(a,f7.3,a,i3,a,i3,a))')				&
     &	      flag,'(',h,')',knt,avg,'+',dev,			&
     &	      ' mx=',amx,'(',imx,',',jmx,')',			&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,e10.3,a,i6,e10.3,a,e9.3,'//&
     &	      '2(a,e10.3,a,i3,a,i3,a))')			&
     &	      flag,'(',h,')',knt,avg,'+',dev,			&
     &	      ' mx=',amx,'(',imx,',',jmx,')',			&
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	endif

	end subroutine showlv

	subroutine show3d(lu,mx,my,mz,a,h,atype,htype,amiss,header,inc)

!	..Print statistics of one 3-d variable
	implicit none

	integer,intent(in) :: lu		! Output unit
	integer,intent(in) :: mx,my,mz	! Array sizes
	real,intent(in) :: a(mx,my,mz)	! The array
	real,intent(in) :: h(mz)		! The argument(levels)
	character(len=*),intent(in) :: atype	! Type of the variable
	character(len=*),intent(in) :: htype	! Typf of the levels
	real,intent(in) :: amiss		! missing value flag of a
	character(len=*),intent(in) :: header	! A header message
	integer,intent(in) :: inc		! order of the listing

	integer i,j,k
	integer kfr,kto,kinc
	integer imx,imn,jmx,jmn
	integer knt
	real amx,amn
	real avg,dev,d
	logical first

!	..A practical value for the magnitude of the fraction of a real
!	number.

	real rfrcval
	parameter(rfrcval=1.e-5)

!	..function

	logical spv
	real aspv
	spv(aspv)=abs((aspv-amiss)/amiss).le.rfrcval .and. CHECKSPV

	write(lu,'(/a,3(2x,i4,a))') header,mx,' X ',my,' X ',mz
	if(htype.eq.'PRES') then
	  write(lu,'(a,3x,a,2x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','mbar', &
     &	    'count','mean','stdv','maxi','mini'
	elseif(htype.eq.'HGHT') then
	  write(lu,'(a,2x,a,2x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','meter',&
     &	    'count','mean','stdv','maxi','mini'
	elseif(htype.eq.'TEMP') then
	  write(lu,'(a,4x,a,4x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','K',    &
     &	    'count','mean','stdv','maxi','mini'
	else
	  write(lu,'(a,4x,a4,1x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl',htype, &
     &	    'count','mean','stdv','maxi','mini'
	endif

!	..Check the order of the listing, increase or decrease
	if(inc.ge.0) then
	  kfr=1
	  kto=mz
	  kinc=1
	else
	  kfr=mz
	  kto=1
	  kinc=-1
	endif

	do k=kfr,kto,kinc
	  knt=0
	  avg=0.
	  do j=1,my
	    do i=1,mx
	      if(.not.spv(a(i,j,k))) then
		knt=knt+1
	      	avg=avg+(a(i,j,k)-avg)/knt
	      endif
	    end do
	  end do

	  dev=0.
	  knt=0
	  do j=1,my
	    do i=1,mx
	      if(.not.spv(a(i,j,k))) then
		d=a(i,j,k)-avg
		knt=knt+1
		dev=dev+(d*d-dev)/knt
	      endif
	    end do
	  end do
	  dev=sqrt(dev*knt/max(1,knt-1))

	  amx=a(1,1,k)
	  amn=a(1,1,k)
	  first=.true.
	  do j=1,my
	    do i=1,mx
	      if(.not.spv(a(i,j,k))) then
		if(first) then
		  imx=i
		  imn=i
		  jmx=j
		  jmn=j
		  amx=a(imx,jmx,k)
		  amn=a(imn,jmn,k)
		  first=.false.
		else
		  if(a(i,j,k).gt.amx) then
		    amx=a(i,j,k)
		    imx=i
		    jmx=j
		  endif
	      	  if(a(i,j,k).lt.amn) then
		    amn=a(i,j,k)
		    imn=i
		    jmn=j
		  endif
		endif
	      endif
	    end do
	  end do

	  if(atype.eq.'RELH') then
	    avg=avg*100.
	    dev=dev*100.
	    amx=amx*100.
	    amn=amn*100.
	  endif

	  if(htype.eq.'PRES'.or.htype.eq.'HGHT'.or.htype.eq.'INDX') then
	    if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	      write(lu,'(i3,2i7,2i10,'//			&
     &		'2(i10,a,i3,a,i3,a))')				&
     &		k,nint(h(k)),knt,nint(avg),nint(dev),		&
     &		nint(amx),'(',imx,',',jmx,')',			&
     &		nint(amn),'(',imn,',',jmn,')'
	    elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or.	&
     &	      atype.eq.'WIND'.or.atype.eq.'%REH'.or.		&
     &	      atype.eq.'RELH'.or.atype.eq.'MIXR') then
	      write(lu,'(i3,2i7,2f10.2,'//			&
     &		'2(f10.2,a,i3,a,i3,a))')			&
     &		k,nint(h(k)),knt,avg,dev,			&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    elseif(atype.eq.'NORM') then
	      write(lu,'(i3,2i7,2f10.4,'//			&
     &		'2(f10.4,a,i3,a,i3,a))')			&
     &		k,nint(h(k)),knt,avg,dev,			&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    else
	      write(lu,'(i3,2i7,2(1x,e9.3),'//		&
     &		'2(1x,e9.3,a,i3,a,i3,a))')		&
     &		k,nint(h(k)),knt,avg,dev,			&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    endif

	  elseif(htype.eq.'TEMP') then
	    if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	      write(lu,'(i3,f7.2,i7,2i10,'//			&
     &		'2(i10,a,i3,a,i3,a))')				&
     &		k,h(k),knt,nint(avg),nint(dev),			&
     &		nint(amx),'(',imx,',',jmx,')',			&
     &		nint(amn),'(',imn,',',jmn,')'
	    elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or.	&
     &	      atype.eq.'WIND'.or.atype.eq.'%REH'.or.		&
     &	      atype.eq.'RELH'.or.atype.eq.'MIXR') then
	      write(lu,'(i3,f7.2,i7,2f10.2,'//			&
     &		'2(f10.2,a,i3,a,i3,a))')			&
     &		k,h(k),knt,avg,dev,				&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    elseif(atype.eq.'NORM') then
	      write(lu,'(i3,f7.2,i7,2f10.4,'//			&
     &		'2(f10.4,a,i3,a,i3,a))')			&
     &		k,h(k),knt,avg,dev,				&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    else
	      write(lu,'(i3,f7.2,i7,2(1x,e9.3),'//	&
     &		'2(1x,e9.3,a,i3,a,i3,a))')		&
     &		k,h(k),knt,avg,dev,				&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    endif

	  else
	    if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	      write(lu,'(i3,e10.3,i7,2i10,'//		&
     &		'2(i10,a,i3,a,i3,a))')				&
     &		k,h(k),knt,nint(avg),nint(dev),			&
     &		nint(amx),'(',imx,',',jmx,')',			&
     &		nint(amn),'(',imn,',',jmn,')'
	    elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or.	&
     &	      atype.eq.'WIND'.or.atype.eq.'%REH'.or.		&
     &	      atype.eq.'RELH'.or.atype.eq.'MIXR') then
	      write(lu,'(i3,e10.3,i7,2f10.2,'//		&
     &		'2(f10.2,a,i3,a,i3,a))')			&
     &		k,h(k),knt,avg,dev,				&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    elseif(atype.eq.'NORM') then
	      write(lu,'(i3,e10.3,i7,2f10.4,'//		&
     &		'2(f10.4,a,i3,a,i3,a))')			&
     &		k,h(k),knt,avg,dev,				&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    else
	      write(lu,'(i3,1x,e9.3,i7,2(1x,e9.3),'// &
     &		'2(1x,e9.3,a,i3,a,i3,a))')		&
     &		k,h(k),knt,avg,dev,				&
     &		amx,'(',imx,',',jmx,')',			&
     &		amn,'(',imn,',',jmn,')'
	    endif

	  endif
	end do		! k=kfr,kto,kinc
	
	end subroutine show3d	! gdstat()


end module m_checksums

