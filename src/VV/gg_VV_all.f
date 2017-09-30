      subroutine gg_VV_all(p,msq)
      implicit none
c--- Author: J. M. Campbell, April 2014
c--- Matrix element squared for gg -> e- e+ ve ve~  Higgs signal process
c--- and gg -> ZZ NNLO contribution to continuum background
c--- added at the amplitude level, i.e. correctly including interference effects
c--- Effects of massive quarks in the third generation may be included
c--- with contributions from both ZZ and WW intermediate states
c--- (default: included)
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'docheck.f'
      include 'runstring.f'
      include 'nuflav.f'
      logical includegens1and2,includebottom,includetop,includegen3
      integer h1,h2,h34,h56
      double precision p(mxpart,4),pswap(mxpart,4),msq(fn:nf,fn:nf),
     & msqgg,fac,pttwo
      double complex 
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Mlight(2,2,2,2),Mgen3(2,2,2,2),
     & Mamp,faczz,facww
      double complex ggHZZ_bquark(2,2,2,2),ggHZZ_tquark(2,2,2,2),
     & ggHWW_bquark(2,2,2,2),ggHWW_tquark(2,2,2,2)

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.      
c--- set this to true to include massive bottom quark
      includebottom=.true.
c--- set this to true to include massive top quark
      includetop=.true.
c--- set this to true to include 3rd generation of quarks (t,b)
      includegen3=.true.

c--- if set, performs check against numerical results at specific PS point
      docheck=.false.
      
      ggHWW_bquark(:,:,:,:)=czip
      ggHWW_tquark(:,:,:,:)=czip

c--- this is the swap to get to the right configuration for the WW call
      pswap(1:2,:)=p(1:2,:)
      pswap(3,:)=p(5,:)
      pswap(4,:)=p(4,:)
      pswap(5,:)=p(3,:)
      pswap(6,:)=p(6,:)
      
c--- overall factor from getggZZamps.f that is not common
      faczz=4d0*esq**2
c--- overall factor from getggWWamps.f that is not common
      facww=gwsq**2

      msq(:,:)=0d0
      
      Mlight(:,:,:,:)=czip
      Mgen3(:,:,:,:)=czip
      
c--- compute all gg->WW and gg->ZZ amplitudes
      call getggWWamps(pswap,includegens1and2,includegen3,
     & Mlight(:,:,1,1),Mgen3(:,:,1,1))
      call getggZZamps(p,includegens1and2,includebottom,includetop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)
            
c--- compute all gg->H->WW and gg->H->ZZ amplitudes
      call getggHWWamps(pswap,
     & ggHWW_bquark(:,:,1,1),ggHWW_tquark(:,:,1,1))
      call getggHZZamps(p,ggHZZ_bquark,ggHZZ_tquark)

c--- pt cut to reproduce Kauer paper
      if (pttwo(3,4,p) .lt. 1d0) then
        Mloop_uptype=czip
        Mloop_dntype=czip
        Mloop_bquark=czip
        Mloop_tquark=czip
        Mloop_uptype=czip
      endif
      if (pttwo(4,5,p) .lt. 1d0) then
        Mlight=czip
        Mgen3=czip
      endif
      
c--- Hack to isolate different pieces depending on runstring
      if(index(runstring,'_ww') .gt. 0) faczz=czip
      if(index(runstring,'_zz') .gt. 0) facww=czip
      
      msqgg=0d0
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
      
c--- compute total amplitude: note relative sign in interference
      Mamp=
     &  (2d0*Mloop_uptype(h1,h2,h34,h56)
     &  +2d0*Mloop_dntype(h1,h2,h34,h56)
     &      +Mloop_bquark(h1,h2,h34,h56)
     &      +Mloop_tquark(h1,h2,h34,h56)
     &      +ggHZZ_bquark(h1,h2,h34,h56)
     &      +ggHZZ_tquark(h1,h2,h34,h56))*faczz

c--- add contribution from vmu,vtau if necessary
      if (nuflav .gt. 1) then
        msqgg=msqgg+dfloat(nuflav-1)*cdabs(Mamp)**2        
      endif
     
      Mamp=Mamp
     & -(2d0*Mlight(h1,h2,h34,h56)
     &      +Mgen3(h1,h2,h34,h56)
     &      +ggHWW_bquark(h1,h2,h34,h56)
     &      +ggHWW_tquark(h1,h2,h34,h56))*facww
          
      msqgg=msqgg+cdabs(Mamp)**2

      enddo
      enddo
      enddo
      enddo
      
c--- overall common factor extracted (c.f. getggZZamps.f and getggWWamps.f)
      fac=avegg*V*(gsq/(16d0*pisq))**2
      
      msqgg=msqgg*fac
      
      msq(0,0)=msqgg
      
      return
      end
      
      


