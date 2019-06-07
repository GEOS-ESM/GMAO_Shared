 module mod_uvdz

   use mod_param, only : len_trie_ls,len_trio_ls,ls_dim,ls_max_node,jcap
   use mod_param, rerth => con_rerth
   implicit none
   public uveodz,uvoedz,dezouv,dozeuv 
   Contains

   subroutine uveodz(ulnev,vlnod,flnev,flnod,epse,epso,ls_node)
!
      implicit none

      real(kind=kind_evod) ulnev(len_trie_ls,2)
      real(kind=kind_evod) vlnod(len_trio_ls,2)
      real(kind=kind_evod) flnev(len_trie_ls,2)
      real(kind=kind_evod) flnod(len_trio_ls,2)
!
      real(kind=kind_evod)  epse(len_trie_ls)
      real(kind=kind_evod)  epso(len_trio_ls)
!
      integer              ls_node(ls_dim,3)
!
      integer              l,locl,n
!
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif
!
      real(kind=kind_evod) rl,rn,rnp1
!
      real(kind=kind_evod) cons0     !constant
      real(kind=kind_evod) cons2     !constant
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
!
      include 'function_indlsev.h'
      include 'function_indlsod.h'
!
!......................................................................
!
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
!
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
           rl = l
           rn = l
         do indev = indev1 , indev2
            flnev(indev,1) = -rl*ulnev(indev,2) + rn*epso(indev-inddif)*vlnod(indev-inddif,1)
!
            flnev(indev,2) = rl*ulnev(indev,1) + rn*epso(indev-inddif)*vlnod(indev-inddif,2)
            rn = rn + cons2     !constant
         end do
!
      end do
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         rl   = l
         rnp1 = l+2+1
         do indev = indev1 , indev2
            flnev(indev,1) = flnev(indev,1) - rnp1*epse(indev)*vlnod(indev-inddif,1)
!
            flnev(indev,2) = flnev(indev,2) - rnp1*epse(indev)*vlnod(indev-inddif,2)
            rnp1 = rnp1 + cons2     !constant
         end do
!
      end do
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l) - 1
         else
            indev2 = indlsev(jcap  ,l) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         rl   = l
         rn   = l+1
         rnp1 = l+1+1
         do indev = indev1 , indev2
            flnod(indev-inddif,1) = -rl  * vlnod(indev-inddif,2)                 &
                                    -rn  * epse(indev+1     ) * ulnev(indev+1,1) &
                                    +rnp1* epso(indev-inddif) * ulnev(indev  ,1)
!
            flnod(indev-inddif,2) =  rl  * vlnod(indev-inddif,1)                 &
                                    -rn  * epse(indev+1     ) * ulnev(indev+1,2) &
                                    +rnp1* epso(indev-inddif) * ulnev(indev  ,2)
            rn   = rn   + cons2     !constant
            rnp1 = rnp1 + cons2     !constant
         end do
!
      end do
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
!
         if (mod(l,2).eq.mod(jcap+1,2)) then
!          set the even (n-l) terms of the top row to zero
            flnev(indlsev(jcap+1,l),1) = cons0     !constant
            flnev(indlsev(jcap+1,l),2) = cons0     !constant
         else
!          set the  odd (n-l) terms of the top row to zero
            flnod(indlsod(jcap+1,l),1) = cons0     !constant
            flnod(indlsod(jcap+1,l),2) = cons0     !constant
         endif
!
      enddo
!
!......................................................................
!
!
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
            flnev(indev,1)=flnev(indev,1)/rerth
            flnev(indev,2)=flnev(indev,2)/rerth
         enddo
!
         do indod = indod1 , indod2
            flnod(indod,1)=flnod(indod,1)/rerth
            flnod(indod,2)=flnod(indod,2)/rerth
         enddo
!
!
      enddo
!
      return
      end

      subroutine uvoedz(ulnod,vlnev,flnod,flnev,epse,epso,ls_node)
!
      implicit none
!
      real(kind=kind_evod) ulnod(len_trio_ls,2)
      real(kind=kind_evod) vlnev(len_trie_ls,2)
      real(kind=kind_evod) flnod(len_trio_ls,2)
      real(kind=kind_evod) flnev(len_trie_ls,2)
!
      real(kind=kind_evod)  epse(len_trie_ls)
      real(kind=kind_evod)  epso(len_trio_ls)
!
      integer              ls_node(ls_dim,3)
!
!
      integer              l,locl,n
!
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif
!
      real(kind=kind_evod) rl,rn,rnp1
!
      real(kind=kind_evod) cons0     !constant
      real(kind=kind_evod) cons2     !constant
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
!
      include 'function_indlsev.h'
      include 'function_indlsod.h'
!
!......................................................................
!
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
!
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l) - 1
         else
            indev2 = indlsev(jcap  ,l) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         rl   = l
         rn   = l+1
         rnp1 = l+1+1
         do indev = indev1 , indev2
!
            flnod(indev-inddif,1) = -rl  * ulnod(indev-inddif,2)                  &
                                    +rn  * epse(indev+1     ) * vlnev(indev+1 ,1) & 
                                    -rnp1* epso(indev-inddif) * vlnev(indev   ,1)
!
            flnod(indev-inddif,2) =  rl  * ulnod(indev-inddif,1)                  &
                                    +rn  * epse(indev+1     ) * vlnev(indev+1 ,2) & 
                                    -rnp1* epso(indev-inddif) * vlnev(indev   ,2)
!
            rn   = rn   + cons2     !constant
            rnp1 = rnp1 + cons2     !constant
         end do
!
      end do
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
           rl = l
           rn = l
         do indev = indev1 , indev2
!
            flnev(indev,1) = -rl * vlnev(indev,2) - rn * epso(indev-inddif) * ulnod(indev-inddif,1)
!
            flnev(indev,2) =  rl * vlnev(indev,1) - rn * epso(indev-inddif) * ulnod(indev-inddif,2)
!
            rn = rn + cons2     !constant
         end do
!
      end do
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         rnp1 = l+2+1
         do indev = indev1 , indev2
!
            flnev(indev ,1) = flnev(indev ,1) + rnp1 * epse(indev) * ulnod(indev-inddif,1)
!
            flnev(indev ,2) = flnev(indev ,2) + rnp1 * epse(indev) * ulnod(indev-inddif,2)
!
            rnp1 = rnp1 + cons2     !constant
         end do
!
      end do
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
!
         if (mod(l,2).eq.mod(jcap+1,2)) then
!          set the even (n-l) terms of the top row to zero
            flnev(indlsev(jcap+1,l),1) = cons0     !constant
            flnev(indlsev(jcap+1,l),2) = cons0     !constant
         else
!          set the  odd (n-l) terms of the top row to zero
            flnod(indlsod(jcap+1,l),1) = cons0     !constant
            flnod(indlsod(jcap+1,l),2) = cons0     !constant
         endif
!
!
      enddo
!
!......................................................................
!
!
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
            flnev(indev,1)=flnev(indev,1)/rerth
            flnev(indev,2)=flnev(indev,2)/rerth
         enddo
!
         do indod = indod1 , indod2
            flnod(indod,1)=flnod(indod,1)/rerth
            flnod(indod,2)=flnod(indod,2)/rerth
         enddo
!
!
      enddo
!
      return
      end

      subroutine dezouv(dev,zod,uev,vod,epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!
      implicit none
!
      real(kind=kind_evod)     dev(len_trie_ls,2)
      real(kind=kind_evod)     zod(len_trio_ls,2)
      real(kind=kind_evod)     uev(len_trie_ls,2)
      real(kind=kind_evod)     vod(len_trio_ls,2)
!
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
!
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
!
      integer              ls_node(ls_dim,3)
!
      integer              l,locl,n
!
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif
!
      real(kind=kind_evod) rl
!
      real(kind=kind_evod) cons0     !constant
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
!
      include 'function_indlsev.h'
      include 'function_indlsod.h'
!
!
!......................................................................
!
!
      cons0 = 0.d0     !constant
!
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
!
         uev(indlsev(l,l),1) = cons0     !constant
         uev(indlsev(l,l),2) = cons0     !constant
!
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            uev(indev,1) = -epsedn(indev) * zod(indev-inddif,1)

            uev(indev,2) = -epsedn(indev) * zod(indev-inddif,2)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            vod(indev-inddif,1) = epsodn(indev-inddif) * dev(indev,1)
!
            vod(indev-inddif,2) = epsodn(indev-inddif) * dev(indev,2)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indev = indev1 , indev2
!             u(l,n)=-i*l*d(l,n)/(n*(n+1))
!
               uev(indev,1) = uev(indev,1) + rl * dev(indev,2) / snnp1ev(indev)
!
               uev(indev,2) = uev(indev,2) - rl * dev(indev,1) / snnp1ev(indev)
!
            enddo
         endif
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod1 = indlsod(l+1,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indod2 = indlsod(jcap  ,l)
         else
            indod2 = indlsod(jcap+1,l) - 1
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indod = indod1 , indod2
!             u(l,n)=-i*l*d(l,n)/(n*(n+1))
!
               vod(indod,1) = vod(indod,1) + rl * zod(indod,2) / snnp1od(indod)
!
               vod(indod,2) = vod(indod,2) - rl * zod(indod,1) / snnp1od(indod)
!
            enddo
         endif
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l) - 1
         else
            indev2 = indlsev(jcap  ,l) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            uev(indev,1) = uev(indev ,1) + epsodn(indev-inddif) * zod(indev-inddif,1)
!
            uev(indev,2) = uev(indev ,2)+ epsodn(indev-inddif) * zod(indev-inddif,2)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            vod(indev-inddif,1) = vod(indev-inddif,1) - epsedn(indev) * dev(indev ,1)
!
            vod(indev-inddif,2) = vod(indev-inddif,2) - epsedn(indev) * dev(indev  ,2)
!
         enddo
!
      enddo
!
!......................................................................
!
!
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
!
            uev(indev,1) = uev(indev,1) * rerth
            uev(indev,2) = uev(indev,2) * rerth
!
         enddo
!
         do indod = indod1 , indod2
!
            vod(indod,1) = vod(indod,1) * rerth
            vod(indod,2) = vod(indod,2) * rerth
!
         enddo
!
      enddo
!
      return
      end

      subroutine dozeuv(dod,zev,uod,vev,epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!
      implicit none
!
      real(kind=kind_evod)     dod(len_trio_ls,2)
      real(kind=kind_evod)     zev(len_trie_ls,2)
      real(kind=kind_evod)     uod(len_trio_ls,2)
      real(kind=kind_evod)     vev(len_trie_ls,2)
!
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
!
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
!
      integer              ls_node(ls_dim,3)
!
      integer              l,locl,n
!
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif
!
      real(kind=kind_evod) rl
!
      real(kind=kind_evod) cons0     !constant
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
!
      include 'function_indlsev.h'
      include 'function_indlsod.h'
!
!
!......................................................................
!
!
      cons0 = 0.d0     !constant
!
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
!
         vev(indlsev(l,l),1) = cons0     !constant
         vev(indlsev(l,l),2) = cons0     !constant
!
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            uod(indev-inddif,1) = -epsodn(indev-inddif) * zev(indev,1)
!
            uod(indev-inddif,2) = -epsodn(indev-inddif) * zev(indev,2)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            vev(indev,1) = epsedn(indev) * dod(indev-inddif,1)
!
            vev(indev,2) = epsedn(indev) * dod(indev-inddif,2)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod1 = indlsod(l+1,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indod2 = indlsod(jcap  ,l)
         else
            indod2 = indlsod(jcap+1,l) - 1
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indod = indod1 , indod2
!             u(l,n)=-i*l*d(l,n)/(n*(n+1))
!
               uod(indod,1) = uod(indod,1) + rl * dod(indod,2) / snnp1od(indod)
!
               uod(indod,2) = uod(indod,2) - rl * dod(indod,1) / snnp1od(indod)
!
            enddo
         endif
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         if ( l .ge. 1 ) then
              rl = l
            do indev = indev1 , indev2
!             u(l,n)=-i*l*d(l,n)/(n*(n+1))
!
               vev(indev,1) = vev(indev,1) + rl * zev(indev,2) / snnp1ev(indev)
!
               vev(indev,2) = vev(indev,2) - rl * zev(indev,1) / snnp1ev(indev)
!
            enddo
         endif
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,l)
         else
            indev2 = indlsev(jcap  ,l)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            uod(indev-inddif,1) = uod(indev-inddif,1) + epsedn(indev) * zev(indev ,1)
!
            uod(indev-inddif,2) = uod(indev-inddif,2) + epsedn(indev) * zev(indev ,2)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,l) - 1
         else
            indev2 = indlsev(jcap  ,l) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
!
         do indev = indev1 , indev2
!
            vev(indev,1) = vev(indev ,1) - epsodn(indev-inddif) * dod(indev-inddif,1)
!
            vev(indev,2) = vev(indev ,2) - epsodn(indev-inddif) * dod(indev-inddif,2)
!
         enddo
!
      enddo
!
!......................................................................
!
!
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
         do indod = indod1 , indod2
!
            uod(indod,1) = uod(indod,1) * rerth
            uod(indod,2) = uod(indod,2) * rerth
!
         enddo
!
         do indev = indev1 , indev2
!
            vev(indev,1) = vev(indev,1) * rerth
            vev(indev,2) = vev(indev,2) * rerth
!
         enddo
!
      enddo
!
      return
      end
  end module mod_uvdz
