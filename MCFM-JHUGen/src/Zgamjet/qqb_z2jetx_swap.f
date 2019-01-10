      subroutine qqb_z2jetx_swap(pin,msq,msqd1,msqx,msqd2)
c--- this is just a wrapper routine to qqb_dirgam_g,
c--- that interchanges p4 and p5 before the call
      implicit none
      include 'constants.f'
      integer nu
      double precision pin(mxpart,4),msq(-nf:nf,-nf:nf),
     & pswap(mxpart,4)
      double precision msqd1(0:2,-nf:nf,-nf:nf),msqd2(0:2,-nf:nf,-nf:nf)
     & ,msqx(0:2,-nf:nf,-nf:nf)
     
      pswap(:,:)=pin(:,:)
      pswap(5,:)=pin(6,:)
      pswap(6,:)=pin(5,:)
      
      call qqb_z2jetx(pswap,msq,msqd1,msqx,msqd2)
      
      return
      end
      
