      subroutine zsvbksb(u,w,v,m,n,b,x) 
c--- Adapted from Numerical Recipes
C--- back substitution for real u,w,v and complex b
      implicit none
      integer i,j,jj,m,n,nmax
      parameter (nmax=100) 
      double precision u(m,n),w(n),v(n,n) 
      double complex b(m),x(n),s,tmp(nmax)

      do j=1,n 
        s=dcmplx(0d0,0d0)
        if(w(j).ne.0d0)then 
          do i=1,m 
            s=s+u(i,j)*b(i) 
          enddo
          s=s/w(j) 
        endif 
        tmp(j)=s 
      enddo 

      do j=1,n 
        s=dcmplx(0d0,0d0)
        do jj=1,n 
          s=s+v(j,jj)*tmp(jj) 
        enddo
        x(j)=s 
      enddo
      
      return 
      end 
