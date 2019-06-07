!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: ssstats: spectral statistics

   program ssstats

!
! !USAGE: ssstats  ssfile1 [ssfile2]
!
! !USES:
!
   use m_ss
   use util
   implicit NONE

! !REVISION HISTORY:
!
!! 2003.06.10  C. Cruz:   Initial Code
!
!-------------------------------------------------------------------------
!EOP

   character(len=*), parameter :: myname = 'ssstats'

   type(ss_vect)  :: ss,ss2     ! analysis vector in spectral space

!  Grid

!  Locals

   integer :: n,iargc,rc,k
   integer :: nymd, nhms, ndt
   integer :: im,jm,km,lm,tt,jcap12,jcap,nlath
   real :: wna1,wna2,rmsad,rmspd

   logical :: do_diff_met

!  File names

   character(len=255) :: ssifile, ssifile2

! start

   do_diff_met=.false.
   call init ( ssifile, ssifile2, do_diff_met )

   call myout ( 6, myname, 'read meta data')
   call ss_init ( ssifile, ss, rc )
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot read from file', rc )
   end if

   if(do_diff_met) then
     call myout ( 6, myname, 'read 2nd ss vector')
     call ss_init ( ssifile2, ss2, rc )
     if ( rc .ne. 0 ) then
       call myexit ( myname,'cannot read from file', rc )
     end if
     if(ss%meta%jcap /= ss2%meta%jcap) then
        print *,' Both spectral files must have the same spectral resolution'
        stop
     end if
   end if

   jcap=ss%meta%jcap
   jcap12=(jcap+1)*(jcap+2)
   nlath=ss%meta%latf/2
 
   call myout ( 6, myname, 'Get ss vector')
   call ss_get ( ssifile, ss, rc )
   if ( rc .ne. 0 ) then
     call myexit ( myname,'cannot get ss vector', rc )
   end if

   if(do_diff_met) then
     call myout ( 6, myname, 'Get 2nd ss vector')
     call ss_get ( ssifile2, ss2, rc )
     if ( rc .ne. 0 ) then
       call myexit ( myname,'cannot get ss vector', rc )
     end if
   end if

   write(*,*) ' ### spectral descriptive metrics ###'
   write(*,*)
   write(*,*) ' * pressure *'
   write(*,*)
   write(*,996) 'lvl','MWN1 ','MWN2'
   call desc_spectral_stats(ss%ps,jcap,nlath,wna1)
   if(do_diff_met) call desc_spectral_stats(ss2%ps,jcap,nlath,wna2)
   k=1 ; write(*,997) k,wna1,wna2
   write(*,*)
   write(*,*) ' * virtual temp *'
   write(*,*)
   write(*,996) 'lvl','MWN1 ','MWN2'
   do k=1,ss%meta%nsig
   call desc_spectral_stats(ss%t(:,k),jcap,nlath,wna1)
   if(do_diff_met) call desc_spectral_stats(ss2%t(:,k),jcap,nlath,wna2)
   write(*,997) k,wna1,wna2
   end do
   write(*,*)
   write(*,*) ' * divergence *'
   write(*,*)
   write(*,996) 'lvl','MWN1 ','MWN2'
   do k=1,ss%meta%nsig
   call desc_spectral_stats(ss%d(:,k),jcap,nlath,wna1)
   if(do_diff_met) call desc_spectral_stats(ss2%d(:,k),jcap,nlath,wna2)
   write(*,997) k,wna1,wna2
   end do
   write(*,*)
   write(*,*) ' * vorticity *'
   write(*,*)
   write(*,996) 'lvl','MWN1 ','MWN2'
   do k=1,ss%meta%nsig
   call desc_spectral_stats(ss%z(:,k),jcap,nlath,wna1)
   if(do_diff_met) call desc_spectral_stats(ss2%z(:,k),jcap,nlath,wna2)
   write(*,997) k,wna1,wna2
   end do
   write(*,*)
   write(*,*) ' * mixing ratio *'
   write(*,*)
   write(*,996) 'lvl','MWN1 ','MWN2'
   do k=1,ss%meta%nsig
   call desc_spectral_stats(ss%q(:,k,1),jcap,nlath,wna1)
   if(do_diff_met) call desc_spectral_stats(ss2%q(:,k,1),jcap,nlath,wna2)
   write(*,997) k,wna1,wna2
   end do
   
   if(do_diff_met) then
     write(*,*)
     write(*,*) ' ### spectral difference metrics ###'
     write(*,*)
     write(*,*) ' * pressure *'
     write(*,*)
     write(*,996) 'lvl','RMSAD','RMSPD'
     call diff_spectral_stats(ss%ps,ss2%ps,jcap,nlath,rmsad,rmspd)
     k=1 ; write(*,997) k,rmsad,rmspd
     write(*,*)
     write(*,*)
     write(*,*) ' * virtual temp *'
     write(*,*)
     write(*,996) 'lvl','RMSAD','RMSPD'
     do k=1,ss%meta%nsig
     call desc_spectral_stats(ss%t(:,k),jcap,nlath,wna1)
     if(do_diff_met) call desc_spectral_stats(ss2%t(:,k),jcap,nlath,wna2)
     write(*,997) k,wna1,wna2
     end do
     write(*,*)
     write(*,*)
     write(*,*) ' * divergence *'
     write(*,*)
     write(*,996) 'lvl','RMSAD','RMSPD'
     do k=1,ss%meta%nsig
     call desc_spectral_stats(ss%d(:,k),jcap,nlath,wna1)
     if(do_diff_met) call desc_spectral_stats(ss2%d(:,k),jcap,nlath,wna2)
     write(*,997) k,wna1,wna2
     end do
     write(*,*)
     write(*,*)
     write(*,*) ' * vorticity *'
     write(*,*)
     write(*,996) 'lvl','RMSAD','RMSPD'
     do k=1,ss%meta%nsig
     call desc_spectral_stats(ss%z(:,k),jcap,nlath,wna1)
     if(do_diff_met) call desc_spectral_stats(ss2%z(:,k),jcap,nlath,wna2)
     write(*,997) k,wna1,wna2
     end do
     write(*,*)
     write(*,*)
     write(*,*) ' * mixing ratio *'
     write(*,*)
     write(*,996) 'lvl','RMSAD','RMSPD'
     do k=1,ss%meta%nsig
     call desc_spectral_stats(ss%q(:,k,1),jcap,nlath,wna1)
     if(do_diff_met) call desc_spectral_stats(ss2%q(:,k,1),jcap,nlath,wna2)
     write(*,997) k,wna1,wna2
     end do
   end if
   
! clean up

   call myout ( 6, myname, 'Clean up' )
   call ss_clean ( ss )
   call ss_clean ( ss2 )
   
   call myout ( 6, myname, '-- SSSTATS has successfully ended --' )
   call exit(0)

996 format(a,4x,a,7x,a,8x)
997 format(i3,2x,f10.6,2x,f10.6)

   CONTAINS

!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: init --- 
!
! !DESCRIPTION: parses command line.
!
! !INTERFACE:
!
   subroutine init ( ssifile, ssifile2, do_diff_met )
   implicit NONE

   character*255, intent(out) :: ssifile, ssifile2
   logical, intent(inout) :: do_diff_met

!EOP

!BOC

   character*4, parameter :: myname = 'init'
   character*255 :: res
   integer :: n,iargc,nargs
   character*255 :: argv
   character*255, allocatable :: arg(:)

   nargs =  iargc()
   if( nargs.lt.1 ) then 
     print *,' Usage: ssstats.x ssfile [-file2 ssfile2]'
     stop
   end if
   allocate ( arg(nargs) )

! process options

   do n=1,nargs
     call getarg(n,arg(n))
   enddo
   do n=1,nargs
     if( trim(arg(n)).eq.'-file2'   ) then
       ssifile2 = trim(arg(n+1))
       do_diff_met = .true.
     end if
   end do
   ssifile=trim(arg(1))

   print *, ' File #1 :',trim(ssifile)
   if(do_diff_met)print *, ' File #2 :',trim(ssifile2)

   deallocate(arg)
   rc = 0

   end subroutine init

!******************************************************************************
   SUBROUTINE desc_spectral_stats(spec_coeff1,jcap,nlath,wna1)
!******************************************************************************
! input

   integer :: jcap,nlath
   real :: wna1
   real :: spec_coeff1(0:jcap,0:jcap)

! local

   real :: top,b,den,xlog
   integer :: l,m,mp,mpi,lii,mii
   complex :: spcoeff1_compx_pair

! start

   top=0.
   b=0.
   den=0.
   xlog=0.

! flag=1 -> do descriptive and difference stats (both spec_coeff1 and
! spec_coeff2 must contain data) else do only descriptive stats


   ! do l=0 first

   l=0
   do m=0,jcap-l,2
      ! first even/odd terms (m=0,2,...)
!       write(78,66)   ' m(re) =',m, ' val =', spec_coeff1(m,l)
      ! now odd/even (m=1,3,...)
      if(m+1.le.jcap-l) then
         mp=m+1
!         write(78,66)' m(im) =',mp,' val =', spec_coeff1(mp,l)
      end if
   end do
   ! now do l > 0
   do l=1,jcap
      lii=jcap+1-l
!      write(78,*) ' l(re) = ',l,' l(im) = ',lii
      do m=0,jcap-l,2
      ! first even/odd terms (m=0,2,...)
         mii=jcap-m
!         write(78,66)' m(re)=',m,   ' val(re) =', spec_coeff1(m,l)
!         write(78,66)' m(im)= ',mii,' val(im) =', spec_coeff1(mii,lii)
         spcoeff1_compx_pair=cmplx(spec_coeff1(m,l),spec_coeff1(mii,lii))
         b=spcoeff1_compx_pair*conjg(spcoeff1_compx_pair)
         xlog=log((l*l+m*m)**(0.5))
         top=top+b*xlog
         den=den+b
!         write(78,*)' l*l = ',l*l,' m*m = ',m*m, &
!                 ' log=',log((l*l+m*m)**(0.5)),' top = ',top,' b = ',b
         ! now odd/even (m=1,3,...)
         if(m+1.le.jcap-l)then
            mp=m+1
            mpi=jcap-mp
!            write(78,66)' m(re)=',mp,  ' val(re) =', spec_coeff1(mp,l)
!            write(78,66)' m(im)=',mpi, ' val(im) =', spec_coeff1(mpi,lii)
             spcoeff1_compx_pair=cmplx(spec_coeff1(mp,l),spec_coeff1(mpi,lii))
             b=spcoeff1_compx_pair*conjg(spcoeff1_compx_pair)
             xlog=log((l*l+m*m)**(0.5))
             top=top+b*xlog
             den=den+b
!             write(78,*)' l*l = ',l*l,' m*m = ',m*m, &
!                        ' log=',log((l*l+m*m)**(0.5)),' top = ',top,' b = ',b
         END IF
      end do
   end do

   wna1=exp(top/den)

66 format(a,i4,a,f16.10)
67 format(a,i4,a,i4,a,f20.10,a,f20.10,a,f20.10)

   end SUBROUTINE desc_spectral_stats

!******************************************************************************
   SUBROUTINE diff_spectral_stats(spec_coeff1,spec_coeff2, &
                                  jcap,nlath,rmsad,rmspd)
!******************************************************************************
! input

   integer :: jcap,nlath
   real :: rmsad,rmspd
   real :: spec_coeff1(0:jcap,0:jcap)
   real :: spec_coeff2(0:jcap,0:jcap)

! local

   real :: absp,argp,b1,b2,den,wlm,abs1,abs2,arg1,arg2
   integer :: l,m,mp,mpi,lii,mii
   complex :: spcoeff1_compx_pair,spcoeff2_compx_pair

! start

   absp=0.
   argp=0.
   b1=0.;b2=0.
   abs1=0.
   abs2=0.
   arg1=0.
   arg2=0.
   den=0.
   wlm=0.

   ! do l=0 first

   l=0
   do m=0,jcap-l,2
      ! first even/odd terms (m=0,2,...)
!       write(78,66)   ' m(re) =',m, ' val =', spec_coeff1(m,l)
      ! now odd/even (m=1,3,...)
      if(m+1.le.jcap-l) then
         mp=m+1
!         write(78,66)' m(im) =',mp,' val =', spec_coeff1(mp,l)
      end if
   end do
   ! now do l > 0
   do l=1,jcap
      lii=jcap+1-l
      do m=0,jcap-l,2
      ! first even/odd terms (m=0,2,...)
         mii=jcap-m
!         write(78,*)' l,m,lii,mii :',l,m,lii,mii
!         write(78,'(2(f16.12,1x))')spec_coeff1(m,l),spec_coeff1(mii,lii)
!         write(78,'(2(f16.12,1x))')spec_coeff2(m,l),spec_coeff2(mii,lii)

         spcoeff1_compx_pair=cmplx(spec_coeff1(m,l),spec_coeff1(mii,lii))
         spcoeff2_compx_pair=cmplx(spec_coeff2(m,l),spec_coeff2(mii,lii))
         b1=spcoeff1_compx_pair*conjg(spcoeff1_compx_pair)
         b2=spcoeff2_compx_pair*conjg(spcoeff2_compx_pair)
         wlm=(b1**0.5)*(b2**0.5)
!         write(78,'(a,3(f16.12,1x))')' b1,b2,wlm :',b1,b2,wlm

         abs1=(spec_coeff1(m,l)**2+spec_coeff1(mii,lii)**2)**0.5
         abs2=(spec_coeff2(m,l)**2+spec_coeff2(mii,lii)**2)**0.5
         absp=absp+wlm*(abs1-abs2)**2
!         write(78,'(a,3(f16.12,1x))')'abs1,abs2,absp :',abs1,abs2,absp

         if(spec_coeff1(m,l).ne.0.0.and.spec_coeff2(m,l).ne.0) then
            arg1=atan(spec_coeff1(mii,lii)/spec_coeff1(m,l))
            arg2=atan(spec_coeff2(mii,lii)/spec_coeff2(m,l))
            argp=argp+wlm*(arg1-arg2)**2
!            write(78,'(a,3(f16.12,1x))')' arg1,arg2,argp :',arg1,arg2,argp
         end if
         den=den+wlm

         ! now odd/even (m=1,3,...)
         if(m+1.le.jcap-l)then
            mp=m+1
            mpi=jcap-mp
!         write(78,*)' mpi :',mpi
!         write(78,'(2(f16.12,1x))')spec_coeff1(mp,l),spec_coeff1(mpi,lii)
!         write(78,'(2(f16.12,1x))')spec_coeff2(mp,l),spec_coeff2(mpi,lii)

            spcoeff1_compx_pair=cmplx(spec_coeff1(mp,l),spec_coeff1(mpi,lii))
            spcoeff2_compx_pair=cmplx(spec_coeff2(mp,l),spec_coeff2(mpi,lii))
            b1=spcoeff1_compx_pair*conjg(spcoeff1_compx_pair)
            b2=spcoeff2_compx_pair*conjg(spcoeff2_compx_pair)
            wlm=(b1**0.5)*(b2**0.5)
!         write(78,'(a,3(f16.12,1x))')' b1,b2,wlm :',b1,b2,wlm

            abs1=(spec_coeff1(mp,l)**2+spec_coeff1(mpi,lii)**2)**0.5
            abs2=(spec_coeff2(mp,l)**2+spec_coeff2(mpi,lii)**2)**0.5
            absp=absp+wlm*(abs1-abs2)**2
!         write(78,'(a,3(f16.12,1x))')'abs1,abs2,absp :',abs1,abs2,absp

            if(spec_coeff1(mp,l).ne.0.0.and.spec_coeff2(mp,l).ne.0) then
               arg1=atan(spec_coeff1(mpi,lii)/spec_coeff1(mp,l))
               arg2=atan(spec_coeff2(mpi,lii)/spec_coeff2(mp,l))
               argp=argp+wlm*(arg1-arg2)**2
            end if
!         write(78,'(a,3(f16.12,1x))')' arg1,arg2,argp :',arg1,arg2,argp

            den=den+wlm

         END IF
      end do
   end do
   rmsad = (absp/den)**0.5
   rmspd = (argp/den)**0.5

66 format(a,i4,a,f16.10)
67 format(a,i4,a,i4,a,f20.10,a,f20.10,a,f20.10)

   end SUBROUTINE diff_spectral_stats

   end program ssstats
