      subroutine dkqqb_ww_v(p,msq)
      implicit none
      include 'types.f'

c----Virtual corrections in decay
C----averaged over initial colours and spins
c     q(-p1)+q~(-p2)-->  W^+ (nu(p3)+e^+(p4))+W^- (s(p5)+cbar(p6))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'masses.f'
      include 'qcdcouple.f'
c      include 'ewcouple.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s,c0,ct
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      call qqb_ww(p,msq)
      call coefswdk(s(5,6),ct,c0)
      msq(:,:)=msq(:,:)*ason2pi*CF*(c0+ct)
      return
      end
