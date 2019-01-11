      subroutine ludcmp(ain,n,indx,d,a)
c--- Adapted from Numerical Recipes
c--- extended so that original matrix a is not destroyed
      implicit none
      integer n,indx(n)
      double precision d,a(n,n),ain(n,n),TINY
      parameter (tiny=1d-20)
      integer i,imax,j,k
      double precision aamax,dum,sum,vv(n)
      a=ain
      d=1d0
      do i=1,n
         aamax=0d0
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo                  ! 11
         if (aamax.eq.0d0) stop 'singular matrix.'
         vv(i)=1d0/aamax
      enddo                     ! 12
      do j=1,n
         if (j.gt.1) then
            do i=1,j-1
               sum=a(i,j)
               if (i.gt.1)then
                  do k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
                  enddo         ! 13
                  a(i,j)=sum
               endif
            enddo               ! 14
         endif
         aamax=0d0
         do i=j,n
            sum=a(i,j)
            if (j.gt.1)then
               do k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
               enddo            ! 15
               a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo                  ! 16
         if (j.ne.imax)then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo               ! 17
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0d0) a(j,j)=tiny
         if(j.ne.n)then
            dum=1d0/a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            enddo               ! 18
         endif
      enddo                     ! 19
      return
      end
