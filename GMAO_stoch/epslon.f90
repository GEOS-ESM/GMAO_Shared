  subroutine epslon(epse,epso,epsedn,epsodn,ls_node)

      use mod_param, only : jcap,ls_max_node,kind_evod
      use mod_param, only : len_trie_ls,len_trio_ls,ls_dim
      implicit none

      real(kind=kind_evod)   epse(len_trie_ls)
      real(kind=kind_evod)   epso(len_trio_ls)
      real(kind=kind_evod) epsedn(len_trie_ls)
      real(kind=kind_evod) epsodn(len_trio_ls)
      integer                  ls_node(ls_dim,3)
      integer                  l,locl,n
      integer                  indev
      integer                  indod
      real(kind=kind_evod) f1,f2,rn,val
      real(kind=kind_evod) cons0     !constant
       integer                  indlsev,jbasev
       integer                  indlsod,jbasod

      include 'function_indlsev.h'
      include 'function_indlsod.h'

      cons0=0.0d0     !constant
!!......................................................................

      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev=indlsev(l,l)
         epse  (indev)=cons0     !constant
         epsedn(indev)=cons0     !constant
          indev=indev+1

         do n=l+2,jcap+1,2
            rn=n
            f1=n*n-l*l
            f2=4*n*n-1
            val=sqrt(f1/f2)
            epse  (indev)=val
            epsedn(indev)=val/rn
            indev=indev+1
         enddo

      enddo

!!......................................................................

      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod=indlsod(l+1,l)

         do n=l+1,jcap+1,2
            rn=n
            f1=n*n-l*l
            f2=4*n*n-1
            val=sqrt(f1/f2)
            epso  (indod)=val
            epsodn(indod)=val/rn
            indod=indod+1
         enddo

      enddo

      return
  end
