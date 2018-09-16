      subroutine gg_ZZ(p,msqgg)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, September 2013
c--- Matrix element squared for the process gg->ZZ
c--- Effects of massive quarks in the third generation may be included
c--- (default: included)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'docheck.f'
      include 'interference.f'
      logical:: includegens1and2,includebottom,includetop
      integer:: h1,h2,h34,h56
      real(dp):: p(mxpart,4),pswap(mxpart,4),msqgg,fac,oprat
      complex(dp)::
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Sloop_uptype(2,2,2,2),Sloop_dntype(2,2,2,2),
     & Sloop_bquark(2,2,2,2),Sloop_tquark(2,2,2,2),
     & Mamp,Samp

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.
c--- set this to true to include massive bottom quark
      includebottom=.true.
c--- set this to true to include massive top quark
      includetop=.true.

c--- if set, performs check against numerical results at specific PS point
      docheck=.false.

c--- compute all gg->ZZ amplitudes
      call getggZZamps(p,includegens1and2,includebottom,includetop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)

      if (interference) then
c--- for interference, compute amplitudes after 4<->6 swap
       pswap(1,:)=p(1,:)
       pswap(2,:)=p(2,:)
       pswap(3,:)=p(3,:)
       pswap(4,:)=p(6,:)
       pswap(5,:)=p(5,:)
       pswap(6,:)=p(4,:)
       call getggZZamps(pswap,includegens1and2,includebottom,includetop,
     &  Sloop_uptype,Sloop_dntype,Sloop_bquark,Sloop_tquark)
      endif

      msqgg=0._dp
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2

      Mamp=
     &  +two*Mloop_uptype(h1,h2,h34,h56)
     &  +two*Mloop_dntype(h1,h2,h34,h56)
     &      +Mloop_bquark(h1,h2,h34,h56)
     &      +Mloop_tquark(h1,h2,h34,h56)

      if (interference .eqv. .false.) then
c--- normal case
        msqgg=msqgg+abs(Mamp)**2
      else
c--- with interference
        Samp=
     &  +two*Sloop_uptype(h1,h2,h34,h56)
     &  +two*Sloop_dntype(h1,h2,h34,h56)
     &      +Sloop_bquark(h1,h2,h34,h56)
     &      +Sloop_tquark(h1,h2,h34,h56)
        if (h34 == h56) then
          oprat=1._dp-two*real(conjg(Mamp)*Samp)
     &                 /(abs(Mamp)**2+abs(Samp)**2)
        else
          oprat=1._dp
        endif
        if (bw34_56) then
          msqgg=msqgg+two*abs(Mamp)**2*oprat
        else
          msqgg=msqgg+two*abs(Samp)**2*oprat
        endif
      endif

      enddo
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggZZamps.f)
      fac=avegg*V*(4._dp*esq*gsq/(16._dp*pisq)*esq)**2

      msqgg=msqgg*fac*vsymfact

      return
      end




