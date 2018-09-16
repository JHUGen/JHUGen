      subroutine gg_hVV_tb(p,msq)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, April 2014
c--- Matrix element squared for gg -> H -> e- e+ ve ve~ signal process
c--- The exact result for massive bottom and top quark loops is included,
c--- with contributions from both ZZ and WW intermediate states
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'qlfirst.f'
      include 'runstring.f'
      include 'nuflav.f'
c      include 'interference.f'
      integer:: h1,h2,h34,h56
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac,
     & pswap(mxpart,4),oprat,pttwo
      complex(dp):: ggHZZ_bquark(2,2,2,2),ggHZZ_tquark(2,2,2,2),
     & ggHWW_bquark(2,2,2,2),ggHWW_tquark(2,2,2,2),Ahiggs,
     & faczz,facww

      if (qlfirst) then
        qlfirst=.false.
        call qlinit
      endif

      ggHWW_bquark(:,:,:,:)=czip
      ggHWW_tquark(:,:,:,:)=czip

c--- debug to check phase: checked by JC 4/1/14 (really)
c      call getggHZZamps(p,ggHZZ_bquark,ggHZZ_tquark)
c      call getggHWWamps(p,ggHWW_bquark(:,:,1,1),ggHWW_tquark(:,:,1,1))
c      do h1=1,2
c      do h2=1,2
c      write(6,*) ggHZZ_bquark(h1,h2,1,1)/ggHWW_bquark(h1,h2,1,1)
c      write(6,*) ggHZZ_tquark(h1,h2,1,1)/ggHWW_tquark(h1,h2,1,1)
c      enddo
c      enddo
c      pause

c--- this is the swap to get to the right configuration for the WW call
      pswap(1:2,:)=p(1:2,:)
      pswap(3,:)=p(5,:)
      pswap(4,:)=p(4,:)
      pswap(5,:)=p(3,:)
      pswap(6,:)=p(6,:)

      msq(:,:)=0._dp

      call getggHWWamps(pswap,
     & ggHWW_bquark(:,:,1,1),ggHWW_tquark(:,:,1,1))
      call getggHZZamps(p,ggHZZ_bquark,ggHZZ_tquark)

c--- overall factor from getggHZZamps.f that is not common
      faczz=4._dp*esq**2
c--- overall factor from getggHWWamps.f that is not common
      facww=gwsq**2

c--- pt cut to reproduce Kauer paper
c      if (pttwo(3,4,p) < 1._dp) then
c        faczz=czip
c      endif
c      if (pttwo(4,5,p) < 1._dp) then
c        facww=czip
c      endif

c--- Hack to isolate different pieces depending on runstring
      if(index(runstring,'_ww') > 0) faczz=czip
      if(index(runstring,'_zz') > 0) facww=czip

      msqgg=0._dp
      do h1=1,2
      h2=h1
      do h34=1,2
      do h56=1,2

c--- compute total Higgs amplitude: note relative sign in interference
      AHiggs=
     &  +ggHZZ_bquark(h1,h2,h34,h56)*faczz
     &  +ggHZZ_tquark(h1,h2,h34,h56)*faczz

c--- add contribution from vmu,vtau if necessary
      if (nuflav > 1) then
        msqgg=msqgg+real(nuflav-1,dp)*abs(AHiggs)**2
      endif

      AHiggs=AHiggs
     &  -ggHWW_bquark(h1,h2,h34,h56)*facww
     &  -ggHWW_tquark(h1,h2,h34,h56)*facww

      msqgg=msqgg+abs(AHiggs)**2
      enddo
      enddo
      enddo

c--- overall common factor extracted (c.f. getggHZZamps.f and getggHWWamps.f)
      fac=avegg*V*(gsq/(16._dp*pisq))**2

      msq(0,0)=msqgg*fac

      return
      end


