      subroutine zlubksb(a,n,indx,bin,b)
c--- Adapted from Numerical Recipes
C--- back substitution for real a and complex b
c--- extended so that original vector b is not destroyed
      implicit none
      integer n,indx(n)
      double precision a(n,n)
      double complex b(n),bin(n),sum
      integer i,ii,j,ll
      b=bin
      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-dcmplx(a(i,j))*b(j)
            enddo               ! 11
         else if (sum.ne.dcmplx(0d0,0d0)) then
            ii=i
         endif
         b(i)=sum
      enddo                     ! 12
      do i=n,1,-1
         sum=b(i)
         if(i.lt.n)then
            do j=i+1,n
               sum=sum-dcmplx(a(i,j))*b(j)
            enddo               ! 13
         endif
         b(i)=sum/dcmplx(a(i,i))
      enddo                     ! 14
      return
      end
