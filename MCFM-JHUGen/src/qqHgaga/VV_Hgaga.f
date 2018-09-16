      subroutine VV_Hgaga(p,msq)
      implicit none
      include 'types.f'
      
c--- Weak Bosion Fusion : sums up WW and ZZ contributions
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: msq_ww(-nf:nf,-nf:nf),msq_zz(-nf:nf,-nf:nf)
      integer:: j,k
 
      call WW_Hgaga(p,msq_ww)
      call ZZ_Hgaga(p,msq_zz)
      
      do j=-nf,nf
      do k=-nf,nf
        msq(j,k)=msq_ww(j,k)+msq_zz(j,k)
      enddo
      enddo      
      
      return
      end
      

