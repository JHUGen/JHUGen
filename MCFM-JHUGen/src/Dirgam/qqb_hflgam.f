      subroutine qqb_hflgam(p,msq)
      implicit none
      include 'types.f'
C----- Matrix element for f(-p1)+f(-p2)->gamma(p3)+Q(p4)
c-----  where Q is a heavy quark: b (flav=5) or c (flav=4)
c-----  treated in the massless approximation
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'heavyflav.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qg,gq,ag,ga

c--- initalize to zero
      msq(:,:)=0._dp
      call dotem(3,p,s)
      fac=4._dp*V*gsq*esq

      qg=-fac*aveqg*(s(1,3)/s(1,2)+s(1,2)/s(1,3))
      ag=qg
      gq=-fac*aveqg*(s(1,2)/s(2,3)+s(2,3)/s(1,2))
      ga=gq

      msq( flav,0)=Q(flav)**2*qg
c      msq(-flav,0)=Q(flav)**2*ag
      msq(0, flav)=Q(flav)**2*gq
c      msq(0,-flav)=Q(flav)**2*ga
      
      return
      end


      
