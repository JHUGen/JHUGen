      subroutine integratehisto(N)
      implicit none
      include 'types.f'
c--- Compute cumulative integral of histogram N and use it to
c--- replace current histogram N 
      
      include 'histo.f'
      integer:: N,j
      real(dp):: xint
      
      xint=0._dp
      
      do j=1,NBIN(N)
      xint=xint+HIST(N,j)*HDEL(N)
      HIST(N,j)=xint
      enddo
      
      j=index(title(N),'+INTEGRAL+')
      title(N)(j:j+9)='Cumulative'
      
      return
      end
      
      
