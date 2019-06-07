 module get_pattern 

   use stoch_data

   implicit none 
   public get_pattern_sppt,get_pattern_skeb,           & 
          scalarspect_to_grid,scalargrid_to_spect,     &    
          vrtdivspect_to_uvgrid,uvgrid_to_vrtdivspect
   Contains 
 
   !========================================================================

   subroutine get_pattern_sppt(sppt2d,PREF,ls_node,ls_nodes,max_ls_nodes,  &
                               lats_nodes_r,global_lats_r,lonsperlar,      &
                               plnev_r,plnod_r)
   ! generate random pattern for sppt.

     use mod_param
     use patterngenerator

     implicit none

     real(kind=kind_evod),intent(out) :: sppt2d(lonr,lats_node_r)
     real(kind=kind_evod),intent( in) :: plnev_r(len_trie_ls,latr2),  & 
                                         plnod_r(len_trio_ls,latr2)
     real*4,              intent( in) :: PREF(levs+1)
     integer, intent(in)              :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
                                         max_ls_nodes(nodes),lats_nodes_r(nodes), & 
                                         global_lats_r(latr),lonsperlar(latr)
     integer i,j,k,l,lat,ierr,n

     real(kind_evod), dimension(lonr,lats_node_r):: wrk2d

     if (.not. allocated (vfact_sppt)) call get_vfact(PREF)

     sppt2d = 0.
     wrk2d = 0.
     do n=1,nsppt
       call patterngenerator_advance(spec_sppt_e(:,:,n),spec_sppt_o(:,:,n), & 
                                     rpattern_sppt(n))
       call scalarspect_to_grid(spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),wrk2d,&
                                ls_node,ls_nodes,max_ls_nodes,              &
                                lats_nodes_r,global_lats_r,lonsperlar,      &
                                plnev_r,plnod_r,1)
      sppt2d = sppt2d + wrk2d
     enddo
   ! print *,'min/max sppt2d',minval(sppt2d),maxval(sppt2d)
     if (sppt_logit) sppt2d = (2./(1.+exp(sppt2d)))-1.
     
   end subroutine get_pattern_sppt
 !========================================================================
   subroutine get_pattern_skeb(u3d,v3d,PREF,                          & 
                               ls_node,ls_nodes,max_ls_nodes,         &
                               lats_nodes_r,global_lats_r,lonsperlar, & 
                               epsedn,epsodn,snnp1ev,snnp1od,         & 
                               plnew_r,plnow_r,plnev_r,plnod_r)

   ! generate random patterns for skeb (stochastic kinetic energy backscatter).
   ! output arrays skeb3d_u,skeb3d_v contains u and v patterns for latitudes on this task.

     use mod_param
     use mod_param, only : rerth=>con_rerth
     use mod_param, only : kind_io4,r_kind => kind_io4, kind_evod
     use patterngenerator

     implicit none

     real(kind=kind_evod), intent(inout) :: u3d(lonr,latr,levs),v3d(lonr,latr,levs)
!     real(kind=kind_evod), intent(  out) :: ptrn3d(lonr,latr,levs)
     real*4,               intent(in   ) :: PREF(levs+1)
     real(kind=kind_evod), intent(in   ) :: plnev_r(len_trie_ls,latr2)
     real(kind=kind_evod), intent(in   ) :: plnod_r(len_trio_ls,latr2)
     real(kind=kind_evod), intent(in   ) :: plnew_r(len_trie_ls,latr2)
     real(kind=kind_evod), intent(in   ) :: plnow_r(len_trio_ls,latr2)
     real(kind=kind_evod), intent(in   ) :: epsedn(len_trie_ls),epsodn(len_trio_ls)
     real(kind=kind_evod), intent(in   ) :: snnp1ev(len_trie_ls),snnp1od(len_trio_ls)
     integer             , intent(in   ) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes), & 
                                            max_ls_nodes(nodes), lats_nodes_r(nodes), & 
                                            global_lats_r(latr),lonsperlar(latr)

     ! locals
     real(kind=kind_evod), dimension(lonr,latr,levs)        :: udiffg,vdiffg,dissrate
     real(kind=kind_evod), dimension(len_trie_ls,2,levs)    :: vrtspec_e,divspec_e
     real(kind=kind_evod), dimension(len_trio_ls,2,levs)    :: vrtspec_o,divspec_o
     real(kind=kind_evod), dimension(len_trie_ls,2,levs)    :: vrtdiffspec_e,divdiffspec_e
     real(kind=kind_evod), dimension(len_trio_ls,2,levs)    :: vrtdiffspec_o,divdiffspec_o
     real(kind=kind_evod), dimension(len_trie_ls)           :: smoothfact_e,kenorm_e,wavenumsq_e
     real(kind=kind_evod), dimension(len_trio_ls)           :: smoothfact_o,kenorm_o,wavenumsq_o
     complex(r_kind),      dimension((jcap+1)*(jcap+2)/2)   :: workspec
     integer i,j,k,l,n,nn,locl,indev,indod,ierr,jbasev,jbasod,indlsod,indlsev,lat
     include 'function_indlsod.h'
     include 'function_indlsev.h'

     real(r_kind) rnn1,rnn0,rnnmax,disscoef
     integer,parameter          :: lu=20

     !--------------------------------------------------------------

     if (.not. allocated (vfact_skeb)) call get_vfact(PREF)

     rnn1 = skeb_diss_smooth          ! use namelist parameter to define smoothing scale.
     rnn0 = rnn1*(rnn1+1.)
     smoothfact_e=1.; smoothfact_o=1. ! used to smooth dissipation estimate.
     kenorm_e=0.; kenorm_o=0.         ! used to convert forcing pattern to wind field.
     do locl=1,ls_max_node
       l = ls_node(locl,1)
       jbasev = ls_node(locl,2)
       indev = indlsev(l,l)
       jbasod = ls_node(locl,3)
       indod = indlsod(l+1,l)
       do n=l,jcap,2
         rnn1 = n*(n+1.)
         smoothfact_e(indev) = exp(-(rnn1/rnn0))
         kenorm_e(indev) = sqrt(rnn1)/rerth
         indev = indev + 1
       enddo
       do n=l+1,jcap,2
         rnn1 = n*(n+1.)
         smoothfact_o(indod) = exp(-(rnn1/rnn0))
         kenorm_o(indod) = sqrt(rnn1)/rerth
         indod = indod + 1
       enddo
     enddo
     ! set the even and odd (n-l) terms of the top row to zero
     do locl=1,ls_max_node
       l = ls_node(locl,1)
       jbasev = ls_node(locl,2)
       jbasod = ls_node(locl,3)
       if (mod(l,2) .eq. mod(jcap+1,2)) then
         smoothfact_e(indlsev(jcap+1,l)) = 0.
         kenorm_e(indlsev(jcap+1,l)) = 0.
       endif
       if (mod(l,2) .ne. mod(jcap+1,2)) then
         smoothfact_o(indlsod(jcap+1,l)) = 0.
         kenorm_o(indlsod(jcap+1,l)) = 0.
       endif
     enddo
     wavenumsq_e = ((kenorm_e*rerth)**2) ! n*(n+1)
     wavenumsq_o = ((kenorm_o*rerth)**2) 

     ! generate random streamfunction forcing patterns.
     u3d=0; v3d=0.
     do n=1,nskeb
       do k=1,levs
         call patterngenerator_advance_skeb(spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),&
                                            rpattern_skeb(n))
       enddo
           
       ! don't modify spectral arrays used to evolve pattern
       vrtdiffspec_e = spec_skeb_e(:,:,:,n)
       vrtdiffspec_o = spec_skeb_o(:,:,:,n)

       ! apply successive applications of 1-2-1 filter in vertical 
       ! to introduce vertical correlations.

       if (skeb_vfilt(n) > 0) then
         do nn=1,skeb_vfilt(n)
           do k=2,levs-1
             divdiffspec_e(:,:,k) = vrtdiffspec_e(:,:,k+1)+ 2.*vrtdiffspec_e(:,:,k)+ & 
                                    vrtdiffspec_e(:,:,k-1)
             divdiffspec_o(:,:,k) = vrtdiffspec_o(:,:,k+1)+ 2.*vrtdiffspec_o(:,:,k)+ & 
                                    vrtdiffspec_o(:,:,k-1)
           enddo
           divdiffspec_e(:,:,   1)= (1.+1./3.)*vrtdiffspec_e(:,:,2)+      & 
                                     2.*(1.+1./3.)*vrtdiffspec_e(:,:,1)
           divdiffspec_e(:,:,levs)= (1.+1./3.)*vrtdiffspec_e(:,:,levs-1)+ &
                                     2.*(1.+1./3.)*vrtdiffspec_e(:,:,levs)
           divdiffspec_o(:,:,   1)= (1.+1./3.)*vrtdiffspec_o(:,:,2)+      &
                                     2.*(1.+1./3.)*vrtdiffspec_o(:,:,1)
           divdiffspec_o(:,:,levs)= (1.+1./3.)*vrtdiffspec_o(:,:,levs-1)+ & 
                                     2.*(1.+1./3.)*vrtdiffspec_o(:,:,levs)
           vrtdiffspec_e = 0.25*divdiffspec_e
           vrtdiffspec_o = 0.25*divdiffspec_o
         enddo
       end if

       !ptrn3d =0.
       !call scalarspect_to_grid(vrtdiffspec_e,vrtdiffspec_o,ptrn3d,   &
       !                         ls_node,ls_nodes,max_ls_nodes,        &
       !                         lats_nodes_r,global_lats_r,lonsperlar,&
       !                         plnev_r,plnod_r,levs)
       !open (lu,file="pattern.grd",form='unformatted',access='sequential',convert='little_endian')
       !call writout8(lu,ptrn3d,lonr,lats_node_r,levs)
       !close(lu)

       ! ke norm (convert streamfunction forcing to vorticity forcing)
       divdiffspec_e = 0; divdiffspec_o = 0.
       do k=1,levs
         do nn=1,2
           vrtdiffspec_e(:,nn,k) = -kenorm_e**2*vrtdiffspec_e(:,nn,k)*vfact_skeb(k)
           vrtdiffspec_o(:,nn,k) = -kenorm_o**2*vrtdiffspec_o(:,nn,k)*vfact_skeb(k)
         enddo
       enddo

       call vrtdivspect_to_uvgrid(vrtdiffspec_e,vrtdiffspec_o, divdiffspec_e,divdiffspec_o, &
                                  vdiffg,udiffg,ls_node,ls_nodes,max_ls_nodes,lats_nodes_r, & 
                                  global_lats_r,lonsperlar,epsedn,epsodn,snnp1ev,snnp1od,plnev_r,plnod_r,levs)

       ! disscoef = 7.75   ! sqrt(r * D * delt/delE) with rD = 6.e-3 m2s-3 and delE/delt = 1.e-4 m2s-3  
                           ! This should be a tunable parameter
       ! u3d = u3d + disscoef*udiffg; v3d = v3d + disscoef*vdiffg
       u3d = u3d + udiffg; v3d = v3d + vdiffg
     enddo

   end subroutine get_pattern_skeb

!============================================================================
! the routines below are spectral transform routines used internally.
   subroutine scalarspect_to_grid(trie_ls,trio_ls,datag,                &
                                  ls_node,ls_nodes,max_ls_nodes,        &
                                  lats_nodes_r,global_lats_r,lonsperlar,&
                                  plnev_r,plnod_r,nlevs)

      use mod_four
      use mod_param

      implicit none
      real(kind=kind_evod), intent(in ) :: trie_ls(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: trio_ls(len_trio_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: plnev_r(len_trie_ls,latr2),plnod_r(len_trio_ls,latr2)
      real(kind=kind_evod), intent(out) :: datag(lonr,lats_node_r,nlevs)
      integer             , intent(in ) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),      &
                                           nlevs,max_ls_nodes(nodes),lats_nodes_r(nodes), &
                                           global_lats_r(latr),lonsperlar(latr)
      ! local vars
      real(kind=kind_evod) for_gr_r_1(lonrx,nlevs,lats_dim_r)
      real(kind=kind_evod) for_gr_r_2(lonr,nlevs,lats_dim_r)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat
      real (kind=kind_evod) tx1

      call sumfln_slg_gg(trie_ls,trio_ls,lat1s_r,              &
                         plnev_r,plnod_r,                      &
                         nlevs,ls_node,latr2,                  &
                         lats_dim_r,nlevs,for_gr_r_1,          &
                         ls_nodes,max_ls_nodes,                &
                         lats_nodes_r,global_lats_r,           &
                         lats_node_r,ipt_lats_node_r,lon_dim_r,&
                         lonsperlar,lonrx,latr,0)

      do lan=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         call four_to_grid(for_gr_r_1(1,1,lan),for_gr_r_2(1,1,lan),&
                           lonrx,lonr,lons_lat,nlevs)
      enddo  

      datag = 0.
      do lan=1,lats_node_r
        lat      = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)
        do k=1,nlevs
          do i=1,lons_lat
            datag(i,lan,k) = for_gr_r_2(i,k,lan)
          enddo
        enddo
      enddo

      return

   end subroutine scalarspect_to_grid  
   !==============================================================================
   subroutine scalargrid_to_spect(trie_ls,trio_ls,datag,                & 
                                  ls_node,ls_nodes,max_ls_nodes,        &
                                  lats_nodes_r,global_lats_r,lonsperlar,&
                                  plnew_r,plnow_r,nlevs)

      use mod_four
      use mod_param

      implicit none
      real(kind=kind_evod), intent(out) :: trie_ls(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(out) :: trio_ls(len_trio_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: datag(lonr,lats_node_r,nlevs)
      real(kind=kind_evod), intent(in ) :: plnew_r(len_trie_ls,latr2),plnow_r(len_trio_ls,latr2)
      integer             , intent(in ) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes), &
                                           max_ls_nodes(nodes),lats_nodes_r(nodes),  & 
                                           global_lats_r(latr),lonsperlar(latr),nlevs
      ! local vars
      real(kind=kind_evod) for_gr_r_1(lonrx,nlevs,lats_dim_r)
      real(kind=kind_evod) for_gr_r_2(lonr,nlevs,lats_dim_r)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat
      real (kind=kind_evod) tx1

      trie_ls = 0.; trio_ls = 0.

      do lan=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)
        do k=1,nlevs
          do i=1,lons_lat
            for_gr_r_2(i,k,lan) = datag(i,lan,k)
          enddo
        enddo
      enddo

      do lan=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         call grid_to_four(for_gr_r_2(1,1,lan),for_gr_r_1(1,1,lan),&
                           lonr,lonrx,lons_lat,nlevs)
      enddo

      call four2fln_gg(lats_dim_r,nlevs,nlevs,for_gr_r_1,&
                    ls_nodes,max_ls_nodes,&
                    lats_nodes_r,global_lats_r,lon_dim_r,&
                    lats_node_r,ipt_lats_node_r,&
                    lat1s_r,lonrx,latr,latr2,&
                    trie_ls(1,1,1), trio_ls(1,1,1),&
                    plnew_r, plnow_r,&
                    ls_node,0,&
                    nlevs,nlevs)

      return
    end subroutine scalargrid_to_spect 

    !==============================================================================

    subroutine vrtdivspect_to_uvgrid(trie_di,trio_di,trie_ze,trio_ze,       &
                                     uug,vvg,ls_node,ls_nodes,max_ls_nodes, &
                                     lats_nodes_r,global_lats_r,lonsperlar, &
                                     epsedn,epsodn,snnp1ev,snnp1od,plnev_r,plnod_r,nlevs)

      use mod_four
      use mod_param
      use mod_uvdz
      use initialize, only: coslat_r

      implicit none
      real(kind=kind_evod), intent(out) :: uug(lonr,lats_node_r,nlevs)
      real(kind=kind_evod), intent(out) :: vvg(lonr,lats_node_r,nlevs)
      real(kind=kind_evod), intent(in ) :: trie_di(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: trio_di(len_trio_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: trie_ze(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: trio_ze(len_trio_ls,2,nlevs)
      real(kind=kind_evod), intent(in ) :: epsedn(len_trie_ls),epsodn(len_trio_ls),       &
                                           snnp1ev(len_trie_ls),snnp1od(len_trio_ls),     &
                                           plnev_r(len_trie_ls,latr2),plnod_r(len_trio_ls,latr2)
      integer             , intent(in ) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),      &
                                           nlevs,max_ls_nodes(nodes),lats_nodes_r(nodes), & 
                                           global_lats_r(latr),lonsperlar(latr)
      ! local vars
      real(kind=kind_evod) trie_ls(len_trie_ls,2,2*nlevs)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,2*nlevs)
      real(kind=kind_evod) for_gr_r_1(lonrx,2*nlevs,lats_dim_r)
      real(kind=kind_evod) for_gr_r_2(lonr,2*nlevs,lats_dim_r)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat
      real (kind=kind_evod) tx1

      do k=1,nlevs
        call dezouv(trie_di(1,1,k),       trio_ze(1,1,k),&
                    trie_ls(1,1,k), trio_ls(1,1,nlevs+k),&
                    epsedn,epsodn,snnp1ev,snnp1od,ls_node)
        call dozeuv(trio_di(1,1,k),       trie_ze(1,1,k),&
                    trio_ls(1,1,k), trie_ls(1,1,nlevs+k),&
                    epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo

      call sumfln_slg_gg(trie_ls,trio_ls,lat1s_r,       &
                  plnev_r,plnod_r,                      &
                  2*nlevs,ls_node,latr2,                &
                  lats_dim_r,2*nlevs,for_gr_r_1,        &
                  ls_nodes,max_ls_nodes,                &
                  lats_nodes_r,global_lats_r,           &
                  lats_node_r,ipt_lats_node_r,lon_dim_r,&
                  lonsperlar,lonrx,latr,0)

      do lan=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         call four_to_grid(for_gr_r_1(1,1,lan),for_gr_r_2(1,1,lan),&
                           lonrx,lonr,lons_lat,2*nlevs)
      enddo  

      uug = 0.; vvg = 0.
      do lan=1,lats_node_r
        lat      = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)
        ! this is not correct  tx1  = 1./(coslat_r(lat)*rcs2_r(min(lat,latg-lat+1)))
        ! the following give exactly the same result 
        ! tx1      = coslat_r(lat)*rcs2_r(min(lat,latg-lat+1))
        tx1      = 1./ coslat_r(lat)

        do k=1,nlevs
          do i=1,lons_lat
            uug(i,lan,k) = for_gr_r_2(i,k,lan) * tx1
            vvg(i,lan,k) = for_gr_r_2(i,nlevs+k,lan) * tx1
          enddo
        enddo
      enddo

      return
   end subroutine  vrtdivspect_to_uvgrid 

   !==============================================================================

   subroutine uvgrid_to_vrtdivspect(trie_di,trio_di,trie_ze,trio_ze,      &
                                    uug,vvg,ls_node,ls_nodes,max_ls_nodes,&
                                    lats_nodes_r,global_lats_r,lonsperlar,&
                                    epse,epso,snnp1ev,snnp1od,plnew_r,plnow_r,nlevs)

      use mod_four
      use mod_param
      use initialize, only: coslat_r,rcs2_r
      use mod_uvdz

      implicit none
      real(kind=kind_evod), intent(out) :: trie_di(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(out) :: trio_di(len_trio_ls,2,nlevs)
      real(kind=kind_evod), intent(out) :: trie_ze(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(out) :: trio_ze(len_trio_ls,2,nlevs)
      real(kind=kind_evod),  intent(in) :: uug(lonr,lats_node_r,nlevs)
      real(kind=kind_evod),  intent(in) :: vvg(lonr,lats_node_r,nlevs)
      real(kind=kind_evod),  intent(in) :: epse(len_trie_ls),epso(len_trio_ls),      & 
                                           snnp1ev(len_trie_ls),snnp1od(len_trio_ls),&
                                           plnew_r(len_trie_ls,latr2),plnow_r(len_trio_ls,latr2)
      integer             ,  intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),      &
                                           nlevs,max_ls_nodes(nodes),lats_nodes_r(nodes), & 
                                           global_lats_r(latr),lonsperlar(latr)
! local vars
      real(kind=kind_evod) trie_ls(len_trie_ls,2,2*nlevs)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,2*nlevs)
      real(kind=kind_evod) for_gr_r_1(lonrx,2*nlevs,lats_dim_r)
      real(kind=kind_evod) for_gr_r_2(lonr,2*nlevs,lats_dim_r)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat
      real (kind=kind_evod) tx1

      do lan=1,lats_node_r 
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
!         tx1  = sqrt(rcs2_r(min(lat,latg-lat+1)))
         tx1  = coslat_r(lat)*rcs2_r(min(lat,latg-lat+1))
         do k=1,nlevs
           do i=1,lons_lat
             for_gr_r_2(i,k,lan)=uug(i,lan,k)*tx1
             for_gr_r_2(i,nlevs+k,lan)=vvg(i,lan,k)*tx1
           enddo
         enddo
         call grid_to_four(for_gr_r_2(1,1,lan),for_gr_r_1(1,1,lan),&
                           lonr,lonrx,lons_lat,2*nlevs)
      enddo  

      call four2fln_gg(lats_dim_r,2*nlevs,2*nlevs,for_gr_r_1,&
                    ls_nodes,max_ls_nodes,&
                    lats_nodes_r,global_lats_r,lon_dim_r,&
                    lats_node_r,ipt_lats_node_r,&
                    lat1s_r,lonrx,latr,latr2,&
                    trie_ls(1,1,1), trio_ls(1,1,1),&
                    plnew_r, plnow_r,&
                    ls_node,0,&
                    2*nlevs,2*nlevs)
!this is incorrect      1,2*levs)


      do k=1,nlevs
         call uveodz(trie_ls(1,1,k), trio_ls(1,1,nlevs+k),&
                     trie_di(1,1,k), trio_ze(1,1,k),&
                     epse,epso,ls_node)
         call uvoedz(trio_ls(1,1,k), trie_ls(1,1,nlevs+k),&
                     trio_di(1,1,k), trie_ze(1,1,k),&
                     epse,epso,ls_node)
      enddo

      return
   end subroutine uvgrid_to_vrtdivspect
 end module get_pattern 
