      subroutine qqb_dirgam_g_swap(pin,msq)
      implicit none
      include 'types.f'
c--- this is just a wrapper routine to qqb_dirgam_g,
c--- that interchanges p4 and p5 before the call
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pin(mxpart,4),msq(-nf:nf,-nf:nf),
     & pswap(mxpart,4)
     
      pswap(:,:)=pin(:,:)
      pswap(5,:)=pin(4,:)
      pswap(4,:)=pin(5,:)
      
      call qqb_dirgam_g(pswap,msq)
      
      return
      end
      
      
