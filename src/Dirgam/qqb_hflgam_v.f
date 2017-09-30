      subroutine qqb_hflgam_v(p,msq)
      implicit none
c----Matrix element for gamma + heavy quark production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->gamma(p3)+c/b(p4)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'heavyflav.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qaggam,qg,gq
      integer j,k

      fac=4d0*V*gsq*esq*ason2pi

      scheme='tH-V'

      call dotem(4,p,s)
c      qa=+qaggam(1,2,4)*fac*aveqq
c      aq=+qaggam(2,1,4)*fac*aveqq
c      aq=qa
      qg=-qaggam(1,4,2)*fac*aveqg
c      ag=-qaggam(4,1,2)*fac*aveqg
c      ag=qg

      gq=-qaggam(2,4,1)*fac*aveqg
c      ga=-qaggam(4,2,1)*fac*aveqg
c      ga=gq

c--set msq=0 to initalize
      msq(:,:)=0d0

      msq( flav,0)=Q(flav)**2*qg
c      msq(-flav,0)=Q(j)**2*ag      
      msq(0, flav)=Q(flav)**2*gq
c      msq(0,-flav)=Q(k)**2*ga

      return
      end


