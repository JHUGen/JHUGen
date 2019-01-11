      subroutine storeptilde(nd,p)
      include 'constants.f'
      include 'ptilde.f'
      integer nd,i,j
      double precision p(mxpart,4)
      
      do j=1,4
        do i=1,mxpart
        ptilde(nd,i,j)=p(i,j)
        enddo
      enddo
      
      return
      end
      
