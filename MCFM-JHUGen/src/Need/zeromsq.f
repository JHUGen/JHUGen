      subroutine zeromsq(msq,msqv)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)

      do j=-nf,nf
      do k=-nf,nf
        msq(j,k)=0._dp
        msqv(j,k)=0._dp
      enddo      
      enddo      

      return
      end
            
