 
      subroutine dcrft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)
 
      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,nz,i,j
      real x(2*ldx,*),y(ldy,*),scale,table(44002),wrk,z
 
      if (init.ne.0) then
        call rffti(n,table)
      else
!ocl novrec
        do j=1,m
          y(1,j)=x(1,j)
          do i=2,n
            y(i,j)=x(i+1,j)
          enddo
          call rfftb(n,y(1,j),table)
          do i=1,n
            y(i,j)=scale*y(i,j)
          enddo
        enddo
      endif
 
      return
      end
      subroutine scrft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)

      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,nz,i,j
      real x(2*ldx,*),y(ldy,*),scale,table(44002),wrk,z

      if (init.ne.0) then
        call rffti(n,table)
      else
!ocl novrec
        do j=1,m
          y(1,j)=x(1,j)
          do i=2,n
            y(i,j)=x(i+1,j)
          enddo
          call rfftb(n,y(1,j),table)
          do i=1,n
            y(i,j)=scale*y(i,j)
          enddo
        enddo
      endif

      return
      end
c
c***********************************************************************
c
      subroutine csfft(isign,n,scale,x,y,table,work,isys)
 
      implicit none
      integer isign,n,isys,i
      real scale,x(*),y(*),table(*),work(*)
 
      if (isign.eq.0) then
        call rffti(n,table)
      endif
      if (isign.eq.1) then
        y(1)=x(1)
        do i=2,n
          y(i)=x(i+1)
        enddo
        call rfftb(n,y,table)
        do i=1,n
          y(i)=scale*y(i)
        enddo
      endif
 
      return
      end
c
c***********************************************************************
c
      subroutine drcft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)
 
      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,nz,i,j
      real x(ldx,*),y(2*ldy,*),scale,table(44002),wrk,z
 
      if (init.ne.0) then
        call rffti(n,table)
      else
        do j=1,m
          do i=1,n
            y(i,j)=x(i,j)
          enddo
          call rfftf(n,y(1,j),table)
          do i=1,n
            y(i,j)=scale*y(i,j)
          enddo
          do i=n,2,-1
            y(i+1,j)=y(i,j)
          enddo
          y(2,j)   = 0.
          y(n+2,j) = 0.
        enddo
      endif
 
      return
      end
c
c***********************************************************************
c
      subroutine srcft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)

      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,nz,i,j
      real x(ldx,*),y(2*ldy,*),scale,table(44002),wrk,z

      if (init.ne.0) then
        call rffti(n,table)
      else
        do j=1,m
          do i=1,n
            y(i,j)=x(i,j)
          enddo
          call rfftf(n,y(1,j),table)
          do i=1,n
            y(i,j)=scale*y(i,j)
          enddo
          do i=n,2,-1
            y(i+1,j)=y(i,j)
          enddo
          y(2,j)   = 0.
          y(n+2,j) = 0.
        enddo
      endif

      return
      end
c
c***********************************************************************
c
      subroutine scfft(isign,n,scale,x,y,table,work,isys)
 
      implicit none
      integer isign,n,isys,i
      real scale,x(*),y(*),table(*),work(*)
 
      if (isign.eq.0) then
        call rffti(n,table)
      endif
      if (isign.eq.-1) then
        do i=1,n
          y(i)=x(i)
        enddo
        call rfftf(n,y,table)
        do i=1,n
          y(i)=scale*y(i)
        enddo
        do i=n,2,-1
          y(i+1)=y(i)
        enddo
        y(2)=0.
      endif
 
      return
      end
