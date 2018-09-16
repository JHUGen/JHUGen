      subroutine ATLAS_hww2013(p,failed_cuts)
      implicit none
      include 'types.f'
c--- Implementation of H->WW cuts in ATLAS-CONF-2013-030 (Njet=0)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'runstring.f'
      real(dp):: p(mxpart,4)
      logical:: failed_cuts
      real(dp):: pt,etarap
      integer:: i,j,idlept
      real(dp):: ran2,tmp,dot,ptminl1,ptminl2,etamaxmu,etamaxe,
     & etaegapmin,etaegapmax,metrelmin,mll,phimin,ptl1,ptl2,etal1,etal2,
     & mllmin,mllmax,phimax,mwindow,phi,etvec(2),missinget,ptll,ptllmin,
     & Etll,MET,mtrans
      logical:: first,l1muon,l2muon,highMT
      data first/.true./
      save first,highMT,mllmax,phimax

C---- optional cut on mtrans
c      etvec(1)=p(3,1)+p(6,1)
c      etvec(2)=p(3,2)+p(6,2)
c--- transverse mass
c      Etll=sqrt(max(zip,(p(4,4)+p(5,4))**2-(p(4,3)+p(5,3))**2))
c      MET=sqrt(max(zip,etvec(1)**2+etvec(2)**2))
c      mtrans=MET+Etll

c      if(mtrans<300._dp) then
c        failed_cuts=.true.
c         return
c      endif

c--- flag for (e,mu + mu,e) or (e,e) or (mu,mu)
c      idlept=1    ! (e,mu)
c      idlept=2    ! (e,e)
c      idlept=3    ! (mu,mu)

      idlept=1

      ptminl1=25._dp
      ptminl2=15._dp

      etamaxmu=2.5_dp
      etamaxe=2.47_dp
      etaegapmin=1.37_dp
      etaegapmax=1.52_dp

      ptllmin=30._dp

      if     (idlept == 1) then
        mllmin=10._dp
        mwindow=zip
        metrelmin=25._dp
        if (ran2() < 0.5_dp) then
          l1muon=.true.
          l2muon=.false.
        else
          l1muon=.false.
          l2muon=.true.
        endif
      elseif (idlept == 2) then
        mllmin=12._dp
        mwindow=15._dp
        metrelmin=45._dp
        l1muon=.false.
        l2muon=.false.
      elseif (idlept == 3) then
        mllmin=12._dp
        mwindow=15._dp
        metrelmin=45._dp
        l1muon=.true.
        l2muon=.true.
      else
        write(6,*) 'idlept set wrong in ATLAS_hww2013: ',idlept
        stop
      endif


      if(first) then
      mllmax=50._dp
      phimax=1.8_dp
      if (index(runstring,'nomll') > 0) then
        mllmax=1d4
      endif
      if (index(runstring,'basic') > 0) then
        mllmax=1d4
        phimax=1d4
      endif
      if (index(runstring,'highMT') > 0) then
        highMT=.true.
      else
        highMT=.false.
      endif
      first=.false.
      write(6,*)  '******** ATLAS H->WW Search cuts ******'
      write(6,*)  '*                                     *'
      if     (idlept == 1) then
      write(6,*)  '*     Treating W decays as (e,mu)     *'
      elseif (idlept == 2) then
      write(6,*)  '*     Treating W decays as (e,e)      *'
      elseif (idlept == 3) then
      write(6,*)  '*     Treating W decays as (mu,mu)    *'
      endif
      write(6,*)  '*                                     *'
      write(6,99) '*     pt(l1) >         ', ptminl1,     ' GeV     *'
      write(6,99) '*     pt(l2) >         ', ptminl2,     ' GeV     *'
      write(6,99) '*     |eta_mu| <       ', etamaxmu,    '         *'
      write(6,99) '*     |eta_el| <       ', etamaxe,     '         *'
      write(6,99) '*     |eta_el| <       ', etaegapmin,  '         *'
      write(6,99) '*     |eta_el| >       ', etaegapmax,  '         *'
      write(6,99) '*     mll >            ',mllmin,       ' GeV     *'
      write(6,99) '*     mll <            ',mllmax,       ' GeV     *'
      write(6,99) '*     |mll-mZ| >       ',mwindow,      ' GeV     *'
      write(6,99) '*     METrel >         ',metrelmin,    ' GeV     *'
      write(6,99) '*     dphill <         ',phimax,       '         *'
      write(6,99) '*     ptll >           ',ptllmin,      ' GeV     *'
      if (highMT) then
      write(6,*)  '*                                     *'
      write(6,*)  '*     transverse mass > 250 GeV       *'
      write(6,*)  '*                                     *'
      endif
      write(6,*)  '***************************************'
      endif

      failed_cuts = .false.

c--- reject event if transverse mass too small, if requested
      if (highMT) then
        etvec(1:2)=p(3,1:2)+p(6,1:2)
        Etll=sqrt(max(zip,(p(4,4)+p(5,4))**2-(p(4,3)+p(5,3))**2))
        MET=sqrt(max(zip,etvec(1)**2+etvec(2)**2))
        mtrans=MET+Etll
        if (mtrans < 250._dp) then
          failed_cuts=.true.
          return
        endif
      endif

      mll=sqrt(max(zip,two*dot(p,4,5)))

c--- minimum dilepton invariant mass
      if (mll < mllmin) then
        failed_cuts=.true.
        return
      endif

c--- maximum dilepton invariant mass
      if (mll > mllmax) then
        failed_cuts=.true.
        return
      endif

c--- veto on dilepton invariant mass around Z mass
      if (abs(mll-zmass) < mwindow) then
        failed_cuts=.true.
        return
      endif

      ptl1=pt(4,p)
      ptl2=pt(5,p)

      if (ptl2 > ptl1) then
        tmp=ptl1
        ptl1=ptl2
        ptl2=tmp
        etal1=etarap(5,p)
        etal2=etarap(4,p)
      else
        etal1=etarap(4,p)
        etal2=etarap(5,p)
      endif

c--- pt cut on hardest lepton
      if (ptl1 < ptminl1) then
        failed_cuts=.true.
        return
      endif
c--- pt cut on softest lepton
      if (ptl2 < ptminl2) then
        failed_cuts=.true.
        return
      endif

c--- rapidity cut on hardest lepton
      if (l1muon) then
        if (abs(etal1) > etamaxmu) then
          failed_cuts=.true.
          return
        endif
      else
        if ((abs(etal1) > etamaxe) .or.
     &   (abs(etal1)>etaegapmin) .and.(abs(etal1)<etaegapmax))then
          failed_cuts=.true.
          return
        endif
      endif

c--- rapidity cut on softest lepton
      if (l2muon) then
        if (abs(etal2) > etamaxmu) then
          failed_cuts=.true.
          return
        endif
      else
        if ((abs(etal2) > etamaxe) .or.
     &   (abs(etal2)>etaegapmin) .and.(abs(etal2)<etaegapmax))then
          failed_cuts=.true.
          return
        endif
      endif

c--- missing et "relative"
      do j=1,2
        etvec(j)=p(3,j)+p(6,j)
      enddo
      missinget=sqrt(etvec(1)**2+etvec(2)**2)

      phimin=99._dp
      do i=4,5
        phi=(p(i,1)*etvec(1)+p(i,2)*etvec(2))
     &      /sqrt((p(i,1)**2+p(i,2)**2)*missinget**2)
        if (phi > +0.9999999_dp) phi=+1._dp
        if (phi < -0.9999999_dp) phi=-1._dp
        phi=acos(phi)
        if (phi < phimin) phimin=phi
      enddo

      if (phi < pi/two) then
        missinget=missinget*sin(phi)
      endif

      if (missinget < metrelmin) then
        failed_cuts=.true.
        return
      endif

c--- dilepton transverse momentum
      ptll=sqrt((p(4,1)+p(5,1))**2+(p(4,2)+p(5,2))**2)
      if (ptll < ptllmin) then
        failed_cuts=.true.
        return
      endif

c--- dilepton phi
      phi=(p(4,1)*p(5,1)+p(4,2)*p(5,2))
     &    /sqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (phi > +0.9999999_dp) phi=+1._dp
      if (phi < -0.9999999_dp) phi=-1._dp
      phi=acos(phi)
      if (phi > phimax) then
        failed_cuts=.true.
        return
      endif

 99   format(1x,a23,f6.2,a10)

      return
      end




