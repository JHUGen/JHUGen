      subroutine getptilde(nd,p)
      include 'constants.f'
      include 'npart.f'
      include 'ptilde.f'
      integer nd,i,j
      double precision p(mxpart,4)
      
      do j=1,4
        do i=1,npart+2
        p(i,j)=ptilde(nd,i,j)
        enddo
      enddo
      
      return
      end
      
