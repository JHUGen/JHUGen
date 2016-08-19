      subroutine clearcoeffs(coeff)
c--- routine to zero out all integral coefficients in Wbbm calculation   
      implicit none
      include 'constants.f'
      include 'Wbbmlabels.f'
      integer j,k
      
      do j=0,4
      do k=1,20
      coeff(j,k)=czip
      enddo
      enddo
      
      return
      end
      
