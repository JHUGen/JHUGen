      subroutine mprove(a,alud,n,indx,b,x)  
c--- Adapted from Numerical Recipes
c--- (extension to complex b)
      implicit none
      integer n,indx(n)  
      double precision a(n,n),alud(n,n)
      double complex b(n),x(n),sdp,r(n),rout(n)
CU    USES zlubksb  
      integer i,j 
       
      do i=1,n  
        sdp=-b(i)  
        do j=1,n  
          sdp=sdp+a(i,j)*x(j) 
        enddo
        r(i)=sdp  
      enddo
      call zlubksb(alud,n,indx,r,rout)  
      do i=1,n  
        x(i)=x(i)-rout(i)  
      enddo

      return  
      end  
 

