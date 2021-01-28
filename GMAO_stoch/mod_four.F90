 module mod_four 

   ! Module containing Fourier transform routines
   ! 

   use mod_param, only: kind_evod, kind_mpi
   use mod_param, only: me,jcap,ls_dim,nodes,len_trie_ls,len_trio_ls,ls_max_node

   implicit none 
   public four2fln_gg,sumfln_slg_gg,four_to_grid,grid_to_four
   Contains
   subroutine four2fln_gg(workdim,nvarsdim,nvars,four_gr,      &
                          ls_nodes,max_ls_nodes,               &
                          lats_nodes,global_lats,lon_dims,     &
                          lats_node,ipt_lats_node,             &
                          lat1s,londi,latl,latl2,              &
                          flnev,flnod,                         &
                          plnev,plnod,ls_node,nvars_0,         &
                          nvar_zero_top_1,nvar_zero_top_2)
!
!
      use num_parthd
      implicit none
!
      integer              nvarsdim,latl2,latl
      integer              nvars,nvars_0
      integer              nvar_zero_top_1
      integer              nvar_zero_top_2
      integer              workdim
      integer              londi
      integer              lat1s(0:jcap)
      integer              lats_node,ipt_lats_node
      integer              lon_dims   ! NOT USED ! lon_dims(latgd)
!
      real(kind=kind_evod) four_gr(londi,nvarsdim,workdim)
!
      integer              ls_nodes(ls_dim,nodes)
      integer                 max_ls_nodes(nodes)
      integer                   lats_nodes(nodes)
      integer              global_lats(latl)
!
!$$$      real(kind=kind_mpi) works(2,nvars,ls_dim*workdim,nodes)
!$$$      real(kind=kind_mpi) workr(2,nvars,ls_dim*workdim,nodes)
      real(kind=kind_mpi) ,allocatable ::works(:,:,:,:),workr(:,:,:,:)
      
!
      integer                    kpts(1+jcap)
      integer                    kptr(1+jcap)
      integer              sendcounts(1+jcap)
      integer              recvcounts(1+jcap)
      integer                 sdispls(1+jcap)
!
      integer              ierr,ilat,ipt_ls
      integer              lval,node,nvar,ifin
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      real(kind=kind_evod) fp(latl2,2,nvars)
      real(kind=kind_evod) fm(latl2,2,nvars)
!
      real(kind=kind_evod) flnev(len_trie_ls,2,nvars)
      real(kind=kind_evod) flnod(len_trio_ls,2,nvars)
!
      real(kind=kind_evod) plnev(len_trie_ls,latl2)
      real(kind=kind_evod) plnod(len_trio_ls,latl2)
!
      integer              ls_node(ls_dim,3)
!
!    local scalars
!    -------------
!
      integer              j,k,l,n
      integer              lat,lat1
      integer              indev1,indev2
      integer              indod1,indod2
      integer              num_threads
      integer              nvar_thread_max
      integer              nvar_1,nvar_2
      integer              thread
      integer              ipt_wr(4,latl2,ls_max_node)
!
!     statement functions
!     -------------------
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
!
      include 'function_indlsev.h'
      include 'function_indlsod.h'
!
      real(kind=kind_evod) cons0     !constant
      real(kind=kind_evod) cons1     !constant

!AA      integer num_parthds,ip1,ip2,ip3,ip4
      integer ip1,ip2,ip3,ip4

      allocate (works(2,nvars,ls_dim*workdim,nodes))
      allocate (workr(2,nvars,ls_dim*workdim,nodes))

!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      cons0 = 0.d0     !constant
      cons1 = 1.d0     !constant
!
!!
      ifin=lats_node
!!
      kpts   = 0
!$omp parallel do private(node,l,lval,j,lat,nvar)
      do node=1,nodes
        do l=1,max_ls_nodes(node)
          lval = ls_nodes(l,node) + 1
          do j=1,ifin
            lat = global_lats(ipt_lats_node-1+j)
            if ( min(lat,latl-lat+1) >= lat1s(lval-1) ) then
              kpts(node) = kpts(node) + 1
              do nvar=1,nvars
!
                 works(1,nvar,kpts(node),node) = four_gr(2*lval-1,nvars_0+nvar,j)
!
                 works(2,nvar,kpts(node),node) = four_gr(2*lval,nvars_0+nvar,j)
!
              enddo
            endif
          enddo
        enddo
      enddo
!
!
      kptr   = 0
      do l=1,ls_max_node
        ilat   = 1
        do node=1,nodes
          ifin=lats_nodes(node)
!jfe      do j=1,lats_nodes_ext(node)
          do j=1,ifin
            lat    = global_lats(ilat)
            ipt_ls = min(lat,latl-lat+1)
            if( ipt_ls >= lat1s(ls_nodes(l,me+1)) ) then
              kptr(node) = kptr(node) + 1
              if ( lat <= latl2 ) then
                ipt_wr(1,ipt_ls,l) = kptr(node)
                ipt_wr(2,ipt_ls,l) =      node
              else
                ipt_wr(3,ipt_ls,l) = kptr(node)
                ipt_wr(4,ipt_ls,l) =      node
              endif
            endif
             ilat = ilat + 1
          enddo
        enddo
      enddo
!
!
      do node=1,nodes
         sendcounts(node) = kpts(node) * 2 * nvars
         recvcounts(node) = kptr(node) * 2 * nvars
            sdispls(node) = (node-1) * 2*ls_dim*workdim*nvars
      end do

!
!      call mpi_barrier (mc_comp,ierr)
!
!AA      call mpi_alltoallv(works,sendcounts,sdispls,mpi_r_mpi,
!AA     x                   workr,recvcounts,sdispls,mpi_r_mpi,
!AA     x                   mc_comp,ierr)
      workr = works 
!
      num_threads     = min(num_parthds(),nvars)
      nvar_thread_max = (nvars+num_threads-1)/num_threads
!
!    -------------------------------------------------
!    compute the coefficients of the expansion
!    in spherical harmonics of the field at each level
!    -------------------------------------------------
!
!
      do j = 1, ls_max_node   ! start of j loop ########################
!
              l=ls_node(j,1)
         jbasev=ls_node(j,2)
         jbasod=ls_node(j,3)
!
         lat1 = lat1s(l)
!
         indev1 = indlsev(l,l)
!
         indod1 = indlsod(l+1,l)
         if (mod(l,2) == mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l)
            indod2 = indlsod(jcap  ,l)
         else
            indev2 = indlsev(jcap  ,l)
            indod2 = indlsod(jcap+1,l)
         endif
!$omp parallel do shared(fm,fp,workr,ipt_wr) &
!$omp shared(plnev,plnod,flnev,flnod) &
!$omp shared(indev1,indev2,indod1,indod2,jbasev,jbasod) &
!$omp shared(j,l,lat1,nvar_thread_max) &
!$omp private(thread,k,lat,nvar_1,nvar_2,ip1,ip2,ip3,ip4)
!
         do thread=1,num_threads   ! start of thread loop ..............
           nvar_1 = (thread-1)*nvar_thread_max+1
           nvar_2 = min(nvar_1+nvar_thread_max-1,nvars)
!
           do k = nvar_1,nvar_2
             do lat = lat1, latl2
               ip1 = ipt_wr(1,lat,j)
               ip2 = ipt_wr(2,lat,j)
               ip3 = ipt_wr(3,lat,j)
               ip4 = ipt_wr(4,lat,j)
!AA               fp(lat,1,k) = workr(1,k,ip1,ip2) + workr(1,k,ip3,ip4)
!AA               fp(lat,2,k) = workr(2,k,ip1,ip2) + workr(2,k,ip3,ip4)
!AA               fm(lat,1,k) = workr(1,k,ip1,ip2) - workr(1,k,ip3,ip4)
!AA               fm(lat,2,k) = workr(2,k,ip1,ip2) - workr(2,k,ip3,ip4)
               fp(lat,1,k) = works(1,k,ip1,ip2) + works(1,k,ip3,ip4)
               fp(lat,2,k) = works(2,k,ip1,ip2) + works(2,k,ip3,ip4)
               fm(lat,1,k) = works(1,k,ip1,ip2) - works(1,k,ip3,ip4)
               fm(lat,2,k) = works(2,k,ip1,ip2) - works(2,k,ip3,ip4)
               enddo
            enddo
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
            if ( kind_evod == 8 ) then !------------------------------
!
!          compute even real      expansion coefficients
!          compute even imaginary expansion coefficients
!
            call dgemm ('n', 'n', indev2-indev1+1, 2*(nvar_2-nvar_1+1), & 
                        latl2-lat1+1, cons1,                            & !constant
                        plnev(indev1,lat1), len_trie_ls,                &
                        fp(lat1,1,nvar_1), latl2, cons0,                & !constant
                        flnev(indev1,1,nvar_1), len_trie_ls)
!
!          compute odd real      expansion coefficients
!          compute odd imaginary expansion coefficients
!
            call dgemm ('n', 'n', indod2-indod1+1, 2*(nvar_2-nvar_1+1), &
                        latl2-lat1+1, cons1,                            & !constant
                        plnod(indod1,lat1), len_trio_ls,                &
                        fm(lat1,1,nvar_1), latl2, cons0,                & !constant
                        flnod(indod1,1,nvar_1), len_trio_ls)
            else !------------------------------------------------------
!
!          compute even real      expansion coefficients
!          compute even imaginary expansion coefficients
!
            call sgemm ('n', 'n', indev2-indev1+1, 2*(nvar_2-nvar_1+1), & 
                        latl2-lat1+1, cons1,                            & !constant
                        plnev(indev1,lat1), len_trie_ls,                &
                        fp(lat1,1,nvar_1), latl2, cons0,                & !constant
                        flnev(indev1,1,nvar_1), len_trie_ls)
!
!          compute odd real      expansion coefficients
!          compute odd imaginary expansion coefficients
!
            call sgemm ('n', 'n', indod2-indod1+1, 2*(nvar_2-nvar_1+1), &
                        latl2-lat1+1, cons1,                            & !constant
                        plnod(indod1,lat1), len_trio_ls,                & 
                        fm(lat1,1,nvar_1), latl2, cons0,                & !constant
                        flnod(indod1,1,nvar_1), len_trio_ls)
            endif !-----------------------------------------------------
!
            if (mod(l,2) == mod(jcap+1,2)) then
!             set the even (n-l) terms of the top row to zero
               do k = max(nvar_1,nvar_zero_top_1),min(nvar_2,nvar_zero_top_2) 
                  flnev(indev2,1,k) = cons0     !constant
                  flnev(indev2,2,k) = cons0     !constant
               end do
            else
!             set the  odd (n-l) terms of the top row to zero
               do k = max(nvar_1,nvar_zero_top_1), min(nvar_2,nvar_zero_top_2) 
                  flnod(indod2,1,k) = cons0     !constant
                  flnod(indod2,2,k) = cons0     !constant
               end do
            endif
!
         end do   ! end of thread loop .................................
!
      end do   ! end of do j loop ######################################
!
      deallocate(workr,works)
      return
      end

      !-----------------------------------------------------------------------------------!
      !-----------------------------------------------------------------------------------!
      subroutine four_to_grid(syn_gr_a_1,syn_gr_a_2,lon_dim_coef, &
                              lon_dim_grid,lons_lat,lot)
      use num_parthd
      implicit none
!!
      real(kind=kind_evod)     syn_gr_a_1(lon_dim_coef,lot)
      real(kind=kind_evod)     syn_gr_a_2(lon_dim_grid,lot)
      integer                  lon_dim_coef
      integer                  lon_dim_grid
      integer                  lons_lat
      integer                  lot
!________________________________________________________
      real(kind=kind_evod) aux1crs(42002)
      real(kind=kind_evod)     scale_ibm
      integer                  ibmsign
      integer                  init
      integer                  lot_thread
      integer                  num_threads
      integer                  nvar_thread_max
      integer                  nvar_1
      integer                  nvar_2
      integer                  thread
!ML      integer                  num_parthds
!________________________________________________________
      num_threads=min(num_parthds(),lot)
 
      nvar_thread_max=(lot+num_threads-1)/num_threads

 
      if ( kind_evod .eq. 8 ) then !------------------------------------
!$omp parallel do shared(syn_gr_a_1,syn_gr_a_2,lons_lat) &
!$omp shared(lon_dim_coef,lon_dim_grid) &
!$omp shared(lot,num_threads,nvar_thread_max) &
!$omp shared(ibmsign,scale_ibm) &
!$omp private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
 
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1=(thread-1)*nvar_thread_max+1
            nvar_2=min(nvar_1+nvar_thread_max-1,lot)
            lot_thread=nvar_2 - nvar_1 +1

            init=1
            ibmsign=-1
            scale_ibm=1.0d0
            call dcrft(init,                                     &     
                       syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,   &
                       syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,     &
                       lons_lat,lot_thread,ibmsign,scale_ibm,    &
                       aux1crs,22000,                            &
                       aux1crs(22001),20000)
            init=0
            call dcrft(init,                                     &
                       syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,   &
                       syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,     &
                       lons_lat,lot_thread,ibmsign,scale_ibm,    &
                       aux1crs,22000,                            &
                       aux1crs(22001),20000)
 
         enddo  ! fin thread loop ......................................
      else !------------------------------------------------------------
!$omp parallel do shared(syn_gr_a_1,syn_gr_a_2,lons_lat) &
!$omp shared(lon_dim_coef,lon_dim_grid) &
!$omp shared(lot,num_threads,nvar_thread_max) &
!$omp shared(ibmsign,scale_ibm) &
!$omp private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
 
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1=(thread-1)*nvar_thread_max+1
            nvar_2=min(nvar_1+nvar_thread_max-1,lot)
            lot_thread=nvar_2 - nvar_1 +1
 
            init=1
            ibmsign=-1
            scale_ibm=1.0d0
            call scrft(init,                                &
                   syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2,  &
                   syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,    &
                   lons_lat,lot_thread,ibmsign,scale_ibm,   &
                   aux1crs,22000,                           &
                   aux1crs(22001),20000,                    &
                   aux1crs(22001),0)
            init=0
            call scrft(init,                                &
                   syn_gr_a_1(1,nvar_1)   ,lon_dim_coef/2, &
                   syn_gr_a_2(1,nvar_1)   ,lon_dim_grid,   &
                   lons_lat,lot_thread,ibmsign,scale_ibm,  &
                   aux1crs,22000,                          &
                   aux1crs(22001),20000,                   &
                   aux1crs(22001),0)
 
         enddo  ! fin thread loop ......................................
      endif !-----------------------------------------------------------
!!
      return
      end

      !-----------------------------------------------------------------------------------!
      !-----------------------------------------------------------------------------------!
      subroutine grid_to_four(anl_gr_a_2,anl_gr_a_1,lon_dim_grid,lon_dim_coef,lons_lat,lot)

      use num_parthd
      implicit none
!!
      real(kind=kind_evod)     anl_gr_a_2(lon_dim_grid,lot)
      real(kind=kind_evod)     anl_gr_a_1(lon_dim_coef,lot)
      integer                  lon_dim_grid
      integer                  lon_dim_coef
      integer                  lons_lat
      integer                  lot
!________________________________________________________
!      real(kind=kind_evod) aux1crs(42002)
      real(kind=kind_evod) aux1crs(44002)
      real(kind=kind_evod)     scale_ibm,rone
      integer                  ibmsign
      integer                  init
      integer                  lot_thread
      integer                  num_threads
      integer                  nvar_thread_max
      integer                  nvar_1,nvar_2
      integer                  thread
!AA      integer                  num_parthds
!________________________________________________________
      num_threads=min(num_parthds(),lot)
 
      nvar_thread_max=(lot+num_threads-1)/num_threads
 
      if ( kind_evod .eq. 8 ) then !------------------------------------
!$omp parallel do shared(anl_gr_a_1,anl_gr_a_2,lons_lat) &
!$omp shared(lon_dim_coef,lon_dim_grid) &
!$omp shared(lot,num_threads,nvar_thread_max) &
!$omp shared(ibmsign,scale_ibm,rone) &
!$omp private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
 
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1=(thread-1)*nvar_thread_max+1
            nvar_2=min(nvar_1+nvar_thread_max-1,lot)
            lot_thread=nvar_2 - nvar_1 +1
            init=1
            ibmsign=1
            rone=1.0d0
            scale_ibm=rone/lons_lat
            call drcft(init,                               &  
                   anl_gr_a_2(1,nvar_1),   lon_dim_grid,   &
                   anl_gr_a_1(1,nvar_1),   lon_dim_coef/2, & 
                   lons_lat,lot_thread,ibmsign,scale_ibm,  &
                   aux1crs,22000,                          &
                   aux1crs(22001),20000)
            init=0
            call drcft(init,                               &
                   anl_gr_a_2(1,nvar_1),   lon_dim_grid,   &
                   anl_gr_a_1(1,nvar_1),   lon_dim_coef/2, &
                   lons_lat,lot_thread,ibmsign,scale_ibm,  &
                   aux1crs,22000,                          &
                   aux1crs(22001),20000)
 
         enddo  ! fin thread loop ......................................
      else !------------------------------------------------------------
!$omp parallel do shared(anl_gr_a_1,anl_gr_a_2,lons_lat) &
!$omp shared(lon_dim_coef,lon_dim_grid) &
!$omp shared(lot,num_threads,nvar_thread_max) &
!$omp shared(ibmsign,scale_ibm,rone) &
!$omp private(thread,nvar_1,nvar_2,lot_thread,init,aux1crs)
 
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1=(thread-1)*nvar_thread_max+1
            nvar_2=min(nvar_1+nvar_thread_max-1,lot)
            lot_thread=nvar_2 - nvar_1 +1
 
            init=1
            ibmsign=1
            rone=1.0d0
            scale_ibm=rone/lons_lat
            call srcft(init,                               &
                   anl_gr_a_2(1,nvar_1),   lon_dim_grid,   &
                   anl_gr_a_1(1,nvar_1),   lon_dim_coef/2, &
                   lons_lat,lot_thread,ibmsign,scale_ibm,  &
                   aux1crs,22000,                          &
                   aux1crs(22001),20000,                   &
                   aux1crs(22001),0)
            init=0
            call srcft(init,                               &
                   anl_gr_a_2(1,nvar_1),   lon_dim_grid,   &
                   anl_gr_a_1(1,nvar_1),   lon_dim_coef/2, &
                   lons_lat,lot_thread,ibmsign,scale_ibm,  &
                   aux1crs,22000,                          &
                   aux1crs(22001),20000,                   &
                   aux1crs(22001),0) 
 
         enddo  ! fin thread loop ......................................
      endif !-----------------------------------------------------------
!!
      return
      end

      subroutine sumfln_slg_gg(flnev,flnod,lat1s,plnev,plnod,   &
                              nvars,ls_node,latl2,              &
                              workdim,nvarsdim,four_gr,         &
                              ls_nodes,max_ls_nodes,            &
                              lats_nodes,global_lats,           &
                              lats_node,ipt_lats_node,lon_dims, &
                              lons_lat,londi,latl,nvars_0)
!
      use num_parthd
      implicit none
!
      integer lat1s(0:jcap),latl2
!
      integer              nvars,nvars_0
      real(kind=kind_evod) flnev(len_trie_ls,2*nvars)
      real(kind=kind_evod) flnod(len_trio_ls,2*nvars)
!
      real(kind=kind_evod) plnev(len_trie_ls,latl2)
      real(kind=kind_evod) plnod(len_trio_ls,latl2)
!
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of l
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
!    local scalars
!    -------------
!
      integer              j, k, l, lat, lat1, n, kn, n2,indev,indod
!
!    local arrays
!    ------------
!
      real(kind=kind_evod) apev(nvars*2,latl2)
      real(kind=kind_evod) apod(nvars*2,latl2)
      integer              num_threads
      integer              nvar_thread_max
      integer              nvar_1,nvar_2
      integer              thread
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      integer              nvarsdim,latl
      integer              workdim
      integer              londi
      integer lats_node,ipt_lats_node
      integer lon_dims        ! NOT USED legacy lon_dims(latgd)
!
!
      real(kind=kind_evod) four_gr(londi,nvarsdim,workdim)
!
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer                lats_nodes(nodes)
!jfe  integer        global_lats(latg+2*jintmx+2*nypt*(nodes-1))
      integer        global_lats(latl)
      integer               lons_lat(latl)
!
      real(kind=kind_mpi) works(2,nvars,ls_dim*workdim,nodes)
      real(kind=kind_mpi) workr(2,nvars,ls_dim*workdim,nodes)
!
      integer                    kpts(1+jcap)
      integer                    kptr(1+jcap)
      integer              sendcounts(1+jcap)
      integer              recvcounts(1+jcap)
      integer                 sdispls(1+jcap)
!
      integer              ierr,ilat,ipt_ls
      integer              lmax,lval,i,jj
      integer              node,nvar
 
!    for omp buffer copy
      integer ilat_list(nodes)
!
!    statement functions
!    -------------------
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
!
      include 'function_indlsev.h'
      include 'function_indlsod.h'
!
      real(kind=kind_evod), parameter ::  cons0=0.0d0, cons1=1.0d0
!AA      integer num_parthds
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      num_threads     = min(num_parthds(),nvars)
      nvar_thread_max = (nvars+num_threads-1)/num_threads
      kpts   = 0

!
      do j = 1, ls_max_node   ! start of do j loop #####################
!
              l = ls_node(j,1)
         jbasev = ls_node(j,2)
         jbasod = ls_node(j,3)

         indev  = indlsev(l,l)
         indod  = indlsod(l+1,l)
!
         lat1 = lat1s(l)
      if ( kind_evod == 8 ) then !------------------------------------

!$omp parallel do private(thread,nvar_1,nvar_2,n2)
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1 = (thread-1)*nvar_thread_max+1
            nvar_2 = min(nvar_1+nvar_thread_max-1,nvars)
            n2     = 2*(nvar_2-nvar_1+1)

!           compute the even and odd components of the fourier coefficients
!
!           compute the sum of the even real      terms for each level
!           compute the sum of the even imaginary terms for each level
!
!           call dgemm('t','n',latl2-lat1+1, 2*(nvar_2-nvar_1+1),
!    &                 (jcap+3-l)/2,cons1,     !constant
!    &                 plnev(indev,lat1),len_trie_ls,
!    &                 flnev(indev,1,nvar_1),len_trie_ls,cons0,
!    &                 apev(lat1,1,nvar_1), latl2)
            call dgemm('t','n',n2,latl2-lat1+1,      &  
                       (jcap+3-l)/2,                 &
                       cons1,                        &
                       flnev(indev,2*nvar_1-1),      &
                       len_trie_ls,                  &
                       plnev(indev,lat1),            &
                       len_trie_ls,                  &
                       cons0,                        & 
                       apev(2*nvar_1-1,lat1),        &
                       2*nvars                       &
                       )
!
!           compute the sum of the odd real      terms for each level
!           compute the sum of the odd imaginary terms for each level
!
!           call dgemm('t','n',latl2-lat1+1, 2*(nvar_2-nvar_1+1),
!    &                 (jcap+2-l)/2,cons1,     !constant
!    &                 plnod(indod,lat1), len_trio_ls,
!    &                 flnod(indod,1,nvar_1),len_trio_ls,cons0,
!    &                 apod(lat1,1,nvar_1), latl2)
            call dgemm(                            &
                       't',                        &   
                       'n',                        &
                       n2,                         &
                       latl2-lat1+1,               &
                       (jcap+2-l)/2,               &
                       cons1,                      &
                       flnod(indod,2*nvar_1-1),    &
                       len_trio_ls,                &
                       plnod(indod,lat1),          &
                       len_trio_ls,                &
                       cons0,                      &
                       apod(2*nvar_1-1,lat1),      &
                       2*nvars                     &
                       )
!
         enddo   ! end of thread loop ..................................
      else !------------------------------------------------------------
!$omp parallel do private(thread,nvar_1,nvar_2)
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1=(thread-1)*nvar_thread_max+1
            nvar_2=min(nvar_1+nvar_thread_max-1,nvars)
         enddo   ! end of thread loop ..................................
      endif !-----------------------------------------------------------
!
!ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       compute the fourier coefficients for each level
!       -----------------------------------------------
!
         ilat_list(1) = 0
         do node = 1, nodes - 1
           ilat_list(node+1) = ilat_list(node) + lats_nodes(node)
         end do
 
!$omp parallel do private(node,jj,ilat,lat,ipt_ls,nvar,kn,n2)
         do node=1,nodes
           do jj=1,lats_nodes(node)
             ilat   = ilat_list(node) + jj
             lat    = global_lats(ilat)
             ipt_ls = min(lat,latl-lat+1)
             if ( ipt_ls >= lat1s(ls_nodes(j,me+1)) ) then
               kpts(node) = kpts(node) + 1
               kn = kpts(node)
!
               if ( lat <= latl2 ) then
!                                                northern hemisphere
                 do nvar=1,nvars
                   n2 = nvar + nvar
                   works(1,nvar,kn,node) = apev(n2-1,ipt_ls) + apod(n2-1,ipt_ls)
                   works(2,nvar,kn,node) = apev(n2,  ipt_ls) + apod(n2,  ipt_ls)
                 enddo
               else
!                                                southern hemisphere
                 do nvar=1,nvars
                   n2 = nvar + nvar
                   works(1,nvar,kn,node) = apev(n2-1,ipt_ls) - apod(n2-1,ipt_ls)
                   works(2,nvar,kn,node) = apev(n2,  ipt_ls) - apod(n2,  ipt_ls)
                 enddo
                endif
             endif
           enddo
         enddo
!
      enddo   ! end of do j loop #######################################
!
      kptr = 0
      do node=1,nodes
         do l=1,max_ls_nodes(node)
            lval = ls_nodes(l,node)+1
            do j=1,lats_node
               lat = global_lats(ipt_lats_node-1+j)
               if ( min(lat,latl-lat+1) >= lat1s(lval-1) ) then
                  kptr(node) = kptr(node) + 1
               endif
            enddo
         enddo
      enddo
!
!
      n2 = nvars + nvars
!$omp parallel do private(node)
      do node=1,nodes
         sendcounts(node) = kpts(node) * n2
         recvcounts(node) = kptr(node) * n2
            sdispls(node) = (node-1)   * n2 * ls_dim * workdim
      end do
!
!_AA   call mpi_alltoallv(works,sendcounts,sdispls,mpi_r_mpi, 
!_AA .                    workr,recvcounts,sdispls,mpi_r_mpi,
!_AA .                    mc_comp,ierr)
       workr=works
       
!
!$omp parallel do private(j,lat,lmax,nvar,lval,n2)
      do j=1,lats_node
         lat  = global_lats(ipt_lats_node-1+j)
         lmax = min(jcap,lons_lat(lat)/2)
         n2   = lmax + lmax + 3
         if ( n2 <= lons_lat(lat)+2 ) then
           do nvar=1,nvars
             do lval = n2, lons_lat(lat)+2
               four_gr(lval,nvars_0+nvar,j) = cons0
             enddo
           enddo
         endif
      enddo
!
      kptr = 0
!!
!$omp parallel do private(node,l,lval,j,lat,nvar,kn,n2)
      do node=1,nodes
        do l=1,max_ls_nodes(node)
          lval = ls_nodes(l,node)+1
          n2   = lval + lval
          do j=1,lats_node
            lat = global_lats(ipt_lats_node-1+j)
            if ( min(lat,latl-lat+1) >= lat1s(lval-1) ) then
              kptr(node) = kptr(node) + 1
              kn = kptr(node)

              do nvar=1,nvars
                four_gr(n2-1,nvars_0+nvar,j) = workr(1,nvar,kn,node)
                four_gr(n2,  nvars_0+nvar,j) = workr(2,nvar,kn,node)
              enddo
            endif
          enddo
        enddo
      enddo
!
      return
      end

 end module mod_four
