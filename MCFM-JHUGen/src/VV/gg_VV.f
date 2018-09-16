      subroutine gg_VV(p,msqgg)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, April 2014
c--- Matrix element squared for gg -> e- e+ ve ve~ box loops
c--- Effects of massive quarks in the third generation may be included
c--- with contributions from both ZZ and WW intermediate states
c--- (default: included)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'docheck.f'
      include 'runstring.f'
      include 'nuflav.f'
      logical:: includegens1and2,includebottom,includetop,includegen3
      integer:: h1,h2,h34,h56
      real(dp):: p(mxpart,4),pswap(mxpart,4),msqgg,fac,pttwo
      complex(dp)::
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Mlight(2,2,2,2),Mgen3(2,2,2,2),
     & Mamp,faczz,facww

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

c--- this is the swap to get to the right configuration for the WW call
      pswap(1:2,:)=p(1:2,:)
      pswap(3,:)=p(5,:)
      pswap(4,:)=p(4,:)
      pswap(5,:)=p(3,:)
      pswap(6,:)=p(6,:)

c--- overall factor from getggZZamps.f that is not common
      faczz=4._dp*esq**2
c--- overall factor from getggWWamps.f that is not common
      facww=gwsq**2

      Mlight(:,:,:,:)=czip
      Mgen3(:,:,:,:)=czip

c--- compute all gg->WW and gg->ZZ amplitudes
      call getggWWamps(pswap,includegens1and2,includegen3,
     & Mlight(:,:,1,1),Mgen3(:,:,1,1))
      call getggZZamps(p,includegens1and2,includebottom,includetop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)

c--- pt cut to reproduce Kauer paper
      if (pttwo(3,4,p) < 1._dp) then
        faczz=czip
      endif
      if (pttwo(4,5,p) < 1._dp) then
        facww=czip
      endif

c--- Hack to isolate different pieces depending on runstring
      if(index(runstring,'_ww') > 0) faczz=czip
      if(index(runstring,'_zz') > 0) facww=czip

      msqgg=0._dp
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2

      Mamp=
     &  (2._dp*Mloop_uptype(h1,h2,h34,h56)
     &  +2._dp*Mloop_dntype(h1,h2,h34,h56)
     &      +Mloop_bquark(h1,h2,h34,h56)
     &      +Mloop_tquark(h1,h2,h34,h56))*faczz

c--- add contribution from vmu,vtau if necessary
      if (nuflav > 1) then
        msqgg=msqgg+real(nuflav-1,dp)*abs(Mamp)**2
      endif

      Mamp=Mamp
     & -(2._dp*Mlight(h1,h2,h34,h56)
     &      +Mgen3(h1,h2,h34,h56))*facww

      msqgg=msqgg+abs(Mamp)**2

      enddo
      enddo
      enddo
      enddo

c--- overall common factor extracted (c.f. getggZZamps.f and getggWWamps.f)
      fac=avegg*V*(gsq/(16._dp*pisq))**2

      msqgg=msqgg*fac

      return
      end




