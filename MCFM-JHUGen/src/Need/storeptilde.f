      subroutine storeptilde(nd,p)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      integer:: nd,i,j
      real(dp):: p(mxpart,4)
      
      do j=1,4
        do i=1,mxpart
        ptilde(nd,i,j)=p(i,j)
        enddo
      enddo
      
      return
      end
      
