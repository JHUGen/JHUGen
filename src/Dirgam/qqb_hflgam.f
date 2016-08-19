      subroutine qqb_hflgam(p,msq)
C----- Matrix element for f(-p1)+f(-p2)->gamma(p3)+Q(p4)
c-----  where Q is a heavy quark: b (flav=5) or c (flav=4)
c-----  treated in the massless approximation
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'heavyflav.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qg,gq,ag,ga

c--- initalize to zero
      msq(:,:)=0d0
      call dotem(3,p,s)
      fac=4d0*V*gsq*esq

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


      
