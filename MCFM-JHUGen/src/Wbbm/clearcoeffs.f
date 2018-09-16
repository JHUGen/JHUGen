      subroutine clearcoeffs(coeff)
      implicit none
      include 'types.f'
c--- routine to zero out all integral coefficients in Wbbm calculation   
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'Wbbmlabels.f'
      integer:: j,k
      
      do j=0,4
      do k=1,20
      coeff(j,k)=czip
      enddo
      enddo
      
      return
      end
      
