      subroutine setlats_r(lats_nodes_r,global_lats_r,iprint,lonsperlar)
!
      use mod_param, only: nodes,latr,lonr,icolor,liope
      implicit none
!
      integer              lats_nodes_r(nodes)
      integer              global_lats_r(latr)
      integer              iprint,opt,ifin,nodesio
      integer              lonsperlar(latr)
      integer              jcount,jpt,lat,lats_sum,node,i
      integer              ilatpe,ngrptg,ngrptl,ipe,irest,idp
!
      integer,allocatable :: lats_hold(:,:)
      allocate ( lats_hold(latr,nodes) )
!
!AA      opt = 1
      opt = 2
      if (opt == 2) lonsperlar = lonr  ! full grid
      lats_nodes_r = 0
      if (liope .and. icolor == 2) then
        nodesio = 1
      else
        nodesio = nodes
      endif
!
      ngrptg = 0
      do lat=1,latr
         do i=1,lonsperlar(lat)
           ngrptg = ngrptg + 1
         enddo
      enddo
!
!     distribution of the grid
!
      ilatpe = ngrptg/nodesio
      ngrptl = 0
      ipe    = 0
      irest  = 0
      idp    = 1
      do lat=1,latr
        ifin   = lonsperlar(lat)
        ngrptl = ngrptl + ifin
        if (ngrptl*nodesio <= ngrptg+irest) then
           lats_nodes_r(ipe+1)  = lats_nodes_r(ipe+1) + 1
           lats_hold(idp,ipe+1) = lat
           idp = idp + 1
        else
          ipe = ipe + 1
          if (ipe <= nodesio) lats_hold(1,ipe+1) = lat
          idp    = 2
          irest  = irest + ngrptg - (ngrptl-ifin)*nodesio
          ngrptl = ifin
          lats_nodes_r(ipe+1) = lats_nodes_r(ipe+1) + 1
        endif
      enddo
!!
      jpt = 0
      do node=1,nodesio
        if ( lats_nodes_r(node) > 0 ) then
          do jcount=1,lats_nodes_r(node)
            global_lats_r(jpt+jcount) = lats_hold(jcount,node)
          enddo
        endif
        jpt = jpt+lats_nodes_r(node)
      enddo
!
      deallocate (lats_hold)
      if ( iprint .ne. 1 ) return
!
      jpt=0
      do node=1,nodesio
         if ( lats_nodes_r(node) .gt. 0 ) then
            print 600
            lats_sum=0
            do jcount=1,lats_nodes_r(node)
               lats_sum=lats_sum + lonsperlar(global_lats_r(jpt+jcount))
               print 700, node-1,
     x                    node,    lats_nodes_r(node),
     x                    jpt+jcount, global_lats_r(jpt+jcount),
     x                     lonsperlar(global_lats_r(jpt+jcount)),
     x                    lats_sum
            enddo
         endif
         jpt=jpt+lats_nodes_r(node)
      enddo
!
      print 600
!
  600 format ( ' ' )
!
  700 format (  'setlats_r  me=', i4,
     x          '  lats_nodes_r(',  i4, ' )=', i4,
     x          '  global_lats_r(', i4, ' )=', i4,
     x          '  lonsperlar=', i5,
     x          '  lats_sum=',   i6 )
!
      return
      end
