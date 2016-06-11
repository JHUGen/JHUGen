      subroutine gg_hZZ_tb(p,msq)
      implicit none
c--- Author: J. M. Campbell, September 2013
c--- Matrix element squared for gg -> H -> ZZ signal process
c--- The exact result for massive bottom and top quark loops is included
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'qlfirst.f'
      include 'interference.f'
      integer h1,h2,h34,h56
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac,
     & pswap(mxpart,4),oprat
      double complex ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),Ahiggs,
     & ggH_bquark_swap(2,2,2,2),ggH_tquark_swap(2,2,2,2),Ahiggs_swap
      double complex ggH2_bquark(2,2,2,2),ggH2_tquark(2,2,2,2),
     & ggH2_bquark_swap(2,2,2,2),ggH2_tquark_swap(2,2,2,2)

      if (qlfirst) then
        qlfirst=.false. 
        call qlinit
      endif
      
      msq(:,:)=0d0
      
      call getggHZZamps(p,ggH_bquark,ggH_tquark)
      call getggH2ZZamps(p,ggH2_bquark,ggH2_tquark)
      
      if (interference) then
c--- for interference, compute amplitudes after 4<->6 swap
       pswap(1,:)=p(1,:)
       pswap(2,:)=p(2,:)
       pswap(3,:)=p(3,:)
       pswap(4,:)=p(6,:)
       pswap(5,:)=p(5,:)
       pswap(6,:)=p(4,:)
       call getggHZZamps(pswap,ggH_bquark_swap,ggH_tquark_swap)
       call getggH2ZZamps(pswap,ggH2_bquark_swap,ggH2_tquark_swap)
      endif
      
      msqgg=0d0
      do h1=1,2
      h2=h1
      do h34=1,2
      do h56=1,2
      
c--- compute total Higgs amplitude
      AHiggs=
     &  +ggH_bquark(h1,h2,h34,h56)
     &  +ggH_tquark(h1,h2,h34,h56)
     &  +ggH2_bquark(h1,h2,h34,h56)
     &  +ggH2_tquark(h1,h2,h34,h56)
     
      if (interference .eqv. .false.) then
c--- normal case
        msqgg=msqgg+cdabs(AHiggs)**2
      else
c--- with interference
        AHiggs_swap=
     &  +ggH_bquark_swap(h1,h2,h34,h56)
     &  +ggH_tquark_swap(h1,h2,h34,h56)
     &  +ggH2_bquark_swap(h1,h2,h34,h56)
     &  +ggH2_tquark_swap(h1,h2,h34,h56)
        if (h34 .eq. h56) then
          oprat=1d0-2d0*dble(dconjg(AHiggs)*AHiggs_swap)
     &                 /(cdabs(AHiggs)**2+cdabs(AHiggs_swap)**2)
        else
          oprat=1d0
        endif

        msqgg=msqgg+2d0*cdabs(AHiggs)**2*oprat
     &  +2d0*cdabs(AHiggs_swap)**2*oprat

      endif
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggHZZamps.f)
      fac=avegg*V*(4d0*esq*gsq/(16d0*pisq)*esq)**2
      
      msq(0,0)=msqgg*fac*vsymfact

      return
      end
      
      
