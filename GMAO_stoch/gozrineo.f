      subroutine gozrineo_r(plnev_r,plnod_r,
     x                      pddev_r,pddod_r,
     x                      plnew_r,plnow_r,
     x                      epse,epso,rcs2_r,wgt_r,ls_node,num_lat)
cc
      use mod_param, only : kind_evod,jcap,latr2
      use mod_param, only : ls_dim,ls_max_node,len_trie_ls,len_trio_ls
      use mod_param, only : rerth => con_rerth
      implicit none
cc
      real(kind=kind_evod) plnev_r(len_trie_ls,latr2)
      real(kind=kind_evod) plnod_r(len_trio_ls,latr2)
      real(kind=kind_evod) pddev_r(len_trie_ls,latr2)
      real(kind=kind_evod) pddod_r(len_trio_ls,latr2)
      real(kind=kind_evod) plnew_r(len_trie_ls,latr2)
      real(kind=kind_evod) plnow_r(len_trio_ls,latr2)
cc
      real(kind=kind_evod)    epse(len_trie_ls)
      real(kind=kind_evod)    epso(len_trio_ls)
cc
      real(kind=kind_evod) rcs2_r(latr2)
      real(kind=kind_evod)  wgt_r(latr2)
cc
      integer                  ls_node(ls_dim,3)
cc
      integer                  num_lat
cc
cc
      integer                  l,lat,locl,n
cc
      integer                  indev,indev1,indev2
      integer                  indod,indod1,indod2
      integer                  inddif
cc
      real(kind=kind_evod) rn,rnp1,wcsa
cc
      real(kind=kind_evod) cons0     !constant
      real(kind=kind_evod) cons2     !constant
cc
      integer                  indlsev,jbasev
      integer                  indlsod,jbasod
cc
      include 'function_indlsev.h'
      include 'function_indlsod.h'
cc
cc
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
cc
cc
      do lat=1,num_lat
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            jbasod=ls_node(locl,3)
cc
            rn=l
cc
cc
            pddev_r(indlsev(l,l),lat) = -epso(indlsod(l+1,l))
     x                              * plnod_r(indlsod(l+1,l),lat) * rn
            indev1 = indlsev(l,l) + 1
            if (mod(l,2).eq.mod(jcap+1,2)) then
               indev2 = indlsev(jcap-1,l)
            else
               indev2 = indlsev(jcap  ,l)
            endif
            indod1 = indlsod(l+1,l)
            inddif = indev1 - indod1
cc
            rn  =l+2
            rnp1=l+2+1
            do indev = indev1 , indev2
cc
               pddev_r(indev,lat) = epse(indev)
     x                         * plnod_r(indev-inddif  ,lat) * rnp1
     x                            - epso(indev-inddif+1)
     x                         * plnod_r(indev-inddif+1,lat) * rn
cc
               rn   = rn   + cons2     !constant
               rnp1 = rnp1 + cons2     !constant
            enddo
cc
cc......................................................................
            indev1 = indlsev(l,l)
            if (mod(l,2).eq.mod(jcap+1,2)) then
               indev2 = indlsev(jcap+1,l) - 1
            else
               indev2 = indlsev(jcap  ,l) - 1
            endif
            indod1 = indlsod(l+1,l)
            inddif = indev1 - indod1
cc
            rn  =l+1
            rnp1=l+1+1
            do indev = indev1 , indev2
cc
               pddod_r(indev-inddif,lat) = epso(indev-inddif)
     x                                * plnev_r(indev  ,lat) * rnp1
     x                                   - epse(indev+1)
     x                                * plnev_r(indev+1,lat) * rn
cc
               rn   = rn   + cons2     !constant
               rnp1 = rnp1 + cons2     !constant
            enddo
cc
         enddo
cc
cc......................................................................
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            jbasod=ls_node(locl,3)
cc
            if (mod(l,2).eq.mod(jcap+1,2)) then
cc
cc             set the even (n-l) terms of the top row to zero
               pddev_r(indlsev(jcap+1,l),lat) = cons0     !constant
cc
            else
cc
cc             set the  odd (n-l) terms of the top row to zero
               pddod_r(indlsod(jcap+1,l),lat) = cons0     !constant
cc
            endif
cc
         enddo
cc
cc......................................................................
cc
         wcsa=rcs2_r(lat)/rerth
cc
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            jbasod=ls_node(locl,3)
            indev1 = indlsev(l,l)
            indod1 = indlsod(l+1,l)
            if (mod(l,2).eq.mod(jcap+1,2)) then
               indev2 = indlsev(jcap+1,l)
               indod2 = indlsod(jcap  ,l)
            else
               indev2 = indlsev(jcap  ,l)
               indod2 = indlsod(jcap+1,l)
            endif
            do indev = indev1 , indev2
cc
               pddev_r(indev,lat) = pddev_r(indev,lat) * wcsa
               plnew_r(indev,lat) = plnev_r(indev,lat) * wgt_r(lat)
cc
            enddo
cc
            do indod = indod1 , indod2
cc
               pddod_r(indod,lat) = pddod_r(indod,lat) * wcsa
               plnow_r(indod,lat) = plnod_r(indod,lat) * wgt_r(lat)
cc
            enddo
cc
         enddo
cc
      enddo
cc
      return
      end
