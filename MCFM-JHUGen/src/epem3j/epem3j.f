      subroutine epem3j(p,msq)
      implicit none
      include 'types.f'

c--- simple modification of qqb_w_g.f: permuted 1 and 4, 2 and 3
c--- to switch leptons with quarks and added a factor of Nc

c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     u(-p1)+dbar(-p2)--> W^+(n(p3)+e^+(p4))   + g(p5)
C For nwz=-1
c     d(-p1)+ubar(-p2)--> W^-(e^-(p3)+nbar(p4))+ g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq,w1jet

c--- initialize
      msq(:,:)=zip

      call dotem(5,p,s)
c---calculate the propagator
      fac=gwsq**2*gsq*V*xn

      qqbWg= +aveqq*fac*w1jet(4,3,2,1,5)

c--- put result in msq(0,0) element for LO and real, fill the other
c--- elements too to make sure that the virtual CT's work as well
      msq(0,0)=qqbWg
      msq(1,0)=msq(0,0)
      msq(0,1)=msq(0,0)

      return
      end

