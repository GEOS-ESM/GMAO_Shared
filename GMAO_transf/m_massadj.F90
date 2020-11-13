!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_massadj --- Implements dry mass conservation for analysis
! 
! !INTERFACE:
!
      MODULE  m_massadj
            
! !USES:
 
      use m_realkinds, only : r_kind => kind_r8

      Implicit NONE

      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  madj_init
      PUBLIC  madj_pre
      PUBLIC  madj_post
      PUBLIC  madj_clean

!
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  29Jan2007  Todling   Created module, based on da Silva's GEOS-4 code.
!
!EOP
!-------------------------------------------------------------------------

    interface madj_init ; module procedure init_ ; end interface
    interface madj_pre  ; module procedure pre_  ; end interface
    interface madj_post ; module procedure post_ ; end interface
    interface madj_clean; module procedure clean_; end interface

    logical, save :: initilized = .false.
    logical, save :: DEBUG      = .false.

    real(r_kind)                              :: gps_f
    real(r_kind), allocatable, dimension(:)   :: gqcol_f
    real(r_kind), allocatable, dimension(:)   :: dbk
    real(r_kind), allocatable, dimension(:,:) :: qcol_f

!   Hardwire mass fixer method
!   --------------------------
    integer, parameter :: NO_FIX = 0, BK_BASED = 1, SIMPLE_RATIO = 2
    integer            :: method = SIMPLE_RATIO

    CONTAINS

    subroutine init_ ( im, jm, km, nq )

    integer, intent(in) :: im, jm, km, nq

    if( .not. initilized ) then
!      if(im==0.or.jm==0.or.km==0) then
!        allocate(gqcol_f(0), dbk(0), qcol_f(0,0) )
!      else
         allocate(gqcol_f(nq), dbk(km), qcol_f(im,jm) )
!      endif
       initilized = .true.
       if(DEBUG) print*, 'madj_init: im,jm,km,lm=',im,jm,km,nq
    endif

    end subroutine init_

    subroutine clean_ ( )

    if( initilized ) then
       deallocate ( gqcol_f, dbk, qcol_f )
       initilized = .false.
       if(DEBUG) print*, 'madj_clean: done'
    endif
                                                                                                                                             
    end subroutine clean_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  pre_ --- Calculate total dry mass of background
! 
! !INTERFACE:
!

     subroutine pre_ ( im, jm, km, nq, ak, bk,  &
                       ps, delp, q, rc )

! !USES:

! !INPUT PARAMETERS: 
                                               ! First guess dimensions
      integer, intent(in)  :: im               !  zonal
      integer, intent(in)  :: jm               !  meridional
      integer, intent(in)  :: km               !  vertical
      integer, intent(in)  :: nq               !  total # of tracers

      real(r_kind),    intent(in)  ::  ak(km+1)        ! vertical grid a-coefficient
      real(r_kind),    intent(in)  ::  bk(km+1)        ! vertical grid b-coefficient

! !INPUT/OUTPUT PARAMETERS:
                                               ! First Guess/Analysis fields:
      real(r_kind),    intent(in)  ::  ps(im,jm)      
                                               ! Surface pressure (Pa)
      real(r_kind),    intent(in)  ::  delp(im,jm,km) 
                                               ! Delta pressure (Pa)
      real(r_kind),    intent(in)  ::  q(im,jm,km,nq)
                                               ! Specific humdity (kg/kg)
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out)  :: rc              ! Error return code:
                                                 !  0   all is well
                                                 !  1   ...

! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  11Oct2006  da Silva  Initial code (GEOS-4 fvchem)
!  29Jan2007  Todling   Slight adaptation to work in transf/gsi mode.
!
!EOP
!-------------------------------------------------------------------------

!     Local Variables
!     ---------------
      integer :: i, j, lm               ! index

!     Constituent mass conservation
!     -----------------------------
      integer :: n, k
      real(r_kind)    :: alpha, glambda         ! global values
      real(r_kind)    :: qmax, qmin, gps_a
      real(r_kind)    :: lambda(im,jm)
      real(r_kind)    :: ps_corr, delp_min

      rc = 0
      lm = min ( 2, nq )  ! Note: Needs 2 to properly replay GEOS-4.0.3
                          !       ana.eta files. However, in replay mode
                          !       tracers are not updated.

!     if ( im==0 .or. jm==0 .or. km==0 ) return
       if(DEBUG) print*, 'madj_pre: im,jm,km,lm=',im,jm,km,nq

!     Start by computing constituent column mass before analysis
!     ----------------------------------------------------------
      gqcol_f = 0.0
      do n = 1, nq
         qcol_f = 0.0
         do k = 1, km              ! yes, whole column for all tracers
            do j = 1, jm
               do i = 1, im
                  qcol_f(i,j) = qcol_f(i,j) + q(i,j,k,n) * delp(i,j,k)
               end do
            end do
         end do
         gqcol_f(n) = gmean_ ( im, jm, qcol_f )
              gps_f = gmean_ ( im, jm, ps     )
      end do

!    Layer mean bk for mass fixer
!    ----------------------------
      do k = 1, km
         dbk(k) = bk(k+1) - bk(k)  ! normalized layer thickness
      end do

      call pmaxmin_('delpbkg_bot',delp(:,:,km), &
                     qmin, qmax, im*jm,1, 1._r_kind )

      call pmaxmin_('     ps_bkg', ps, &
                     qmin, qmax, im*jm,1, 1._r_kind )

       if(DEBUG) print*, 'madj_pre: done'
    end subroutine pre_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  post_ --- Make dry mass is not changed by analysis
!
! !INTERFACE:
!
    subroutine post_ ( im, jm, km, nq, gid, ptop, &
                       ps, delp, q, rc, delpmin )

! !USES:

! !INPUT PARAMETERS: 

      integer, intent(in)  :: im               !  zonal
      integer, intent(in)  :: jm               !  meridional
      integer, intent(in)  :: km               !  vertical
      integer, intent(in)  :: nq               !  total # of tracers
      integer, intent(in)  :: gid              !  mpi-processor number

      real(r_kind), intent(in)  :: ptop        ! top of atmos pressure (Pa)

! !OUTPUT PARAMETERS: 
                                               ! First Guess/Analysis fields:
      real(r_kind),    intent(inout)  ::  ps(im,jm)      
                                               ! Surface pressure (Pa)
      real(r_kind),    intent(inout)  ::  delp(im,jm,km) 
                                               ! Delta pressure (Pa)
      real(r_kind),    intent(inout)  ::  q(im,jm,km,nq)
                                               ! Specific humdity (kg/kg)

      integer, intent(out) :: rc               ! return error code

      real(r_kind),    optional, intent(out)    ::  delpmin

! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!
!  11Oct2006  da Silva  Initial code (GEOS-4 fvchem)
!  29Jan2007  Todling   Slight adaptation to work in transf/gsi mode.
!
!EOP
!-------------------------------------------------------------------------

!   Local Variables
!   ---------------
    integer :: i, j, lm               ! index
                                                                                                                                           
!   Constituent mass conservation
!   -----------------------------
    integer :: n, k
    real(r_kind)    :: qcol_a(im,jm)
    real(r_kind)    :: gqcol_a                ! global means
    real(r_kind)    :: alpha, glambda         ! global values
    real(r_kind)    :: qmax, qmin, gps_a
    real(r_kind)    :: lambda(im,jm)
    real(r_kind)    :: ps_corr, delp_min
    real(r_kind)    :: delpmin_

    rc = 0
!   if ( im==0 .or. jm==0 .or. km==0 ) return
    if(DEBUG) print*, 'madj_post: im,jm,km,lm=',im,jm,km,nq

!   Dry mass adjustment: make sure analysis preserves dry mass
!   ----------------------------------------------------------
    qcol_a = 0.0
    do k = 1, km
       do j = 1, jm
          do i = 1, im
             qcol_a(i,j) = qcol_a(i,j) + q(i,j,k,1) * delp(i,j,k)
          end do
       end do
    end do
    gqcol_a = gmean_ ( im, jm, qcol_a )
      gps_a = gmean_ ( im, jm, ps     )

   ps_corr = ( gps_a - gps_f ) - ( gqcol_a - gqcol_f(1) )

!  Distribute the mass correction proportionally to delp 
!  Notes: 
!     1) Since correction is proportional to delp, it should work when
!        lowest layer is shaved by the analysis.
!        We apply this correction only for those layers dependent on "ps".
!     2) If one could assume that the analysis is in eta-coordinates, then
!        a more sensible distribution of ps_corr would be:
!           delp_corr(i,j,k) = dbk(k) * ps_corr
!      
!  -------------------------------------------------------------------------

   alpha = 1.0 - ps_corr / ( gps_a - ptop )

   delp(:,:,1:km) = alpha * delp(:,:,1:km) ! no longer eta coords 
               ps = ps - ps_corr

!  Apply dry mass correction to SPHU only
!  --------------------------------------
   q(:,:,:,1) = q(:,:,:,1) / alpha
              
!  This will signal remapping if so needed
!   DO NOT COMMENT!!!!
!  ---------------------------------------
   call pmaxmin_('delpana_bot',delp(:,:,km), &
                  delpmin_, qmax, im*jm,1, 1._r_kind )

   call pmaxmin_('   NO_PRINT', ps, &
                  qmin, qmax, im*jm,1, 1._r_kind )

   if(present(delpmin)) delpmin = delpmin_

   if ( gid .eq. 0 ) print 10, 'ps_dry  ', &
                     qmin, gps_f, qmax, ps_corr, 100.*(alpha-1.)

!    At this point the dry mass correction has been applied
!     and as long we do not replay moisture, the total air mass
!     should be conserved as well.

!    For each constituent (except water), apply column mass fixer
!    Note: most of of the mass non-conservation could be coming 
!          from the shaving of the lowest layer by the analysis
!    ------------------------------------------------------------
     if ( method /= NO_FIX ) then

      do n = 2, nq

!       Compute after analysis mass
!       ---------------------------
        qcol_a = 0.0
        lambda = 0.0
        do k = 1, km
           do j = 1, jm
              do i = 1, im
                 qcol_a(i,j) = qcol_a(i,j) +          q(i,j,k,n) * delp(i,j,k)
                 lambda(i,j) = lambda(i,j) + dbk(k) * q(i,j,k,n) * delp(i,j,k)
              end do
           end do
        end do
        gqcol_a = gmean_ ( im, jm, qcol_a )
        glambda = gmean_ ( im, jm, lambda )

!       Compute final lambda
!       --------------------
        if ( glambda /= 0 ) then
             glambda = ( gqcol_f(n) - gqcol_a ) / glambda
        end if

!       Modify constituent q to ensure mass conservation
!       ------------------------------------------------
        do k = 1, km

           if ( method == BK_BASED ) then
                                                    ! <0 should never happen
              alpha = max ( 0._r_kind, 1._r_kind + glambda * dbk(k) ) 

           else if ( method == SIMPLE_RATIO ) then
              alpha = gqcol_f(n) / gqcol_a

           else
              alpha = 1._r_kind  ! no mass fixer, otherwise

           end if
              
           q(:,:,k,n) = alpha * q(:,:,k,n) 
              
        end do

        call pmaxmin_('NO_PRINT', alpha*qcol_a, &
                      qmin, qmax, im*jm,1,1._r_kind )


        if ( gid .eq. 0 ) print 10, 'Species nuumber: ', n, &
                                    qmin, alpha*gqcol_a, qmax, &
                                    (alpha-1.)*gqcol_a, 100*(alpha-1.)

      end do  ! constituent

   end if     ! method /= NO_FIX

10 format(a8,' min/ave/max:',3(1pE10.2),' - Mass Fix:',2(1pE10.2),'%')


   if(DEBUG) print*, 'madj_post: done'
   end subroutine post_


!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean_ --- Calculate the mean of a 2D field
!
! !INTERFACE:

      real(r_kind) function gmean_(im, jm, q)

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer  im, jm                  ! Horizontal dimensions
      real(r_kind) q(im,jm)            ! 2D field 

! !DESCRIPTION: Calculate global mean of a 2D field
!
!EOP
!---------------------------------------------------------------------

      integer i, j
      real(r_kind) sine(jm),cosp(jm),sinp(jm),cose(jm)
      real(r_kind), allocatable, save :: gw(:)
      real(r_kind) dl, dp, xsum(jm)
      logical first
      data first /.true./

      if(first) then
         call setrig_(im,jm,dp,dl,cosp,cose,sinp,sine)

         allocate( gw(jm) )
            gw(1) = 1 + sine(2)
         do j=2,jm-1
            gw(j) = sine(j+1) - sine(j)
         enddo
            gw(jm) = 1 - sine(jm)
 
        first = .false.
      endif

          xsum = 0.
      do j=1,jm
        do i=1,im
          xsum(j) = xsum(j) + q(i,j)
        enddo
          xsum(j) = xsum(j)*gw(j)
      enddo

      gmean_ = sum(xsum) / (2*im)

      end function gmean_


      subroutine setrig_(im, jm, dp, dl, cosp, cose, sinp, sine)

      implicit none

      integer im, jm
      integer j, jm1
      real(r_kind) sine(jm),cosp(jm),sinp(jm),cose(jm)
      real(r_kind) dp, dl
      double precision pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)
      dp  = pi/dble(jm1)

      do 10 j=2,jm
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
10    sine(j) = dsin(ph5)

      cosp( 1) =  0.
      cosp(jm) =  0.

      do 80 j=2,jm1
80    cosp(j) = (sine(j+1)-sine(j)) / dp

! Define cosine at edges..

      do 90 j=2,jm
90    cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do 100 j=2,jm1
100   sinp(j) = 0.5 * (sine(j) + sine(j+1))

      end subroutine setrig_


! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin_( qname, a, pmin, pmax, im, jt, fac )

#if defined( SPMD )
      use mod_comm, only : gid, mp_reduce_max
#define CPP_PRT_PREFIX  if(gid.eq.0)
#else
#define CPP_PRT_PREFIX
#endif

      implicit none

      character*(*)  qname
      integer im, jt
      integer i, j
      real(r_kind) a(im,jt)
      real(r_kind) qmin(jt), qmax(jt)
      real(r_kind) pmax, pmin
      real(r_kind) fac                     ! multiplication factor
      real(r_kind) pm1(2)

!$omp parallel do private(i, j, pmax, pmin)

      do j=1,jt
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
      pmax = qmax(1)
      pmin = qmin(1)
      do j=2,jt
         pmax = max(pmax, qmax(j))
         pmin = min(pmin, qmin(j))
      enddo


#if defined( SPMD )
      pm1(1) = pmax
      call mp_reduce_max(1, pm1)
      pmax=pm1(1)
      pm1(1) = -pmin
      call mp_reduce_max(1, pm1)
      pmin=-pm1(1)
#endif

      CPP_PRT_PREFIX write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

      end subroutine pmaxmin_

     end module m_massadj
