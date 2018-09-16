      subroutine integratehisto(N)
c--- Compute cumulative integral of histogram N and use it to
c--- replace current histogram N 
      implicit none
      include 'histo.f'
      integer N,j
      double precision xint
      
      xint=0d0
      
      do j=1,NBIN(N)
      xint=xint+HIST(N,j)*HDEL(N)
      HIST(N,j)=xint
      enddo
      
      j=index(title(N),'+INTEGRAL+')
      title(N)(j:j+9)='Cumulative'
      
      return
      end
      
      
