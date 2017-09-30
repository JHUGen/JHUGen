      subroutine gg_ZZ_all(p,msq)
      implicit none
c--- Author: J. M. Campbell, September 2013
c--- Total of gg -> H -> ZZ signal process
c--- and gg -> ZZ NNLO contribution to continuum background
c--- added at the amplitude level, i.e. correctly including interference effects
c--- The effect of massive bottom and top quark loops is included
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'noglue.f'
      include 'qlfirst.f'
      include 'interference.f'
      integer h1,h2,h34,h56
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac,
     & pswap(mxpart,4),oprat
      double complex 
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Sloop_uptype(2,2,2,2),Sloop_dntype(2,2,2,2),
     & Sloop_bquark(2,2,2,2),Sloop_tquark(2,2,2,2),
     & ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),Acont,Ahiggs,
     & ggH_bquark_swap(2,2,2,2),ggH_tquark_swap(2,2,2,2),Ahiggs_swap,
     & Acont_swap,Mamp,Samp
      logical includegens1and2,includebottom,includetop

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.      
c--- set this to true to include massive bottom quark
      includebottom=.true.
c--- set this to true to include massive top quark
      includetop=.true.

      if (qlfirst) then
        qlfirst=.false. 
        call qlinit
      endif

c--- if neither contribution is included print warning message and stop
      if ((includegens1and2 .eqv. .false.) .and.
     &    (includebottom    .eqv. .false.) .and.
     &    (includetop       .eqv. .false.)) then
         write(6,*) 'Box loop is set to zero, please edit gg_ZZ_int.f'
         stop
      endif
c--- if noglue print warning message and stop
      if (noglue) then
         write(6,*) 'Please set noglue .false. in input file'
         stop
      endif
            
      msq(:,:)=0d0
      
c      if (pttwo(3,4,p) .lt. 7d0) return ! Kauer gg2VV cut on |H+C|^2

      call getggZZamps(p,includegens1and2,includebottom,includetop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)

      call getggHZZamps(p,ggH_bquark,ggH_tquark)
      
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
       call getggHZZamps(pswap,ggH_bquark_swap,ggH_tquark_swap)
      endif
      
      msqgg=0d0
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
      
c--- compute total continuum amplitude 
      Acont=     
     &  +2d0*Mloop_uptype(h1,h2,h34,h56)
     &  +2d0*Mloop_dntype(h1,h2,h34,h56)
     &      +Mloop_bquark(h1,h2,h34,h56)
     &      +Mloop_tquark(h1,h2,h34,h56)
c--- compute total Higgs amplitude   
      AHiggs=
     &  +ggH_bquark(h1,h2,h34,h56)   
     &  +ggH_tquark(h1,h2,h34,h56)   

c---- This accumulates all contributions
      Mamp=Acont+AHiggs
      
      if (interference .eqv. .false.) then
c--- normal case
        msqgg=msqgg+cdabs(Mamp)**2
      else
c--- with interference
        Acont_swap=
     &  +2d0*Sloop_uptype(h1,h2,h34,h56)
     &  +2d0*Sloop_dntype(h1,h2,h34,h56)
     &      +Sloop_bquark(h1,h2,h34,h56)
     &      +Sloop_tquark(h1,h2,h34,h56)
        AHiggs_swap=
     &  +ggH_bquark_swap(h1,h2,h34,h56)
     &  +ggH_tquark_swap(h1,h2,h34,h56)
        Samp=Acont_swap+AHiggs_swap
        if (h34 .eq. h56) then
          oprat=1d0-2d0*dble(dconjg(Mamp)*Samp)
     &                 /(cdabs(Mamp)**2+cdabs(Samp)**2)
        else
          oprat=1d0
        endif
        if (bw34_56) then
          msqgg=msqgg+2d0*cdabs(Mamp)**2*oprat
        else
          msqgg=msqgg+2d0*cdabs(Samp)**2*oprat
        endif
      endif

      enddo
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggZZamps.f and getggHZZamps.f )
      fac=avegg*V*(4d0*esq*gsq/(16d0*pisq)*esq)**2
      
      msq(0,0)=msqgg*fac*vsymfact

      return
      end
      
      
