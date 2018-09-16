      subroutine FBhisto(N)
c--- Compute FB asymmetry in histogram N,
c--- assuming F is currently in N-2 and B is in N-1
      implicit none
      include 'types.f'
      include 'histo.f'
      integer N,j
      real(dp):: xrel
      
      do j=1,NBIN(N)
      xrel=(HIST(N-2,j)-HIST(N-1,j))/(HIST(N-2,j)+HIST(N-1,j))
      HIST(N,j)=xrel
      enddo
      
      j=index(title(N),'+FB+')
      title(N)(j:j+3)=' FB '
      
      return
      end
