      subroutine ATLAS_hww2013(p,failed_cuts) 
c--- Implementation of H->WW cuts in ATLAS-CONF-2013-030 (Njet=0)
      implicit none 
      include 'constants.f' 
      include 'masses.f'
      include 'runstring.f'
      double precision p(mxpart,4)
      logical failed_cuts
      double precision pt,etarap
      integer i,j,idlept
      double precision ran2,tmp,dot,ptminl1,ptminl2,etamaxmu,etamaxe,
     & etaegapmin,etaegapmax,metrelmin,mll,phimin,ptl1,ptl2,etal1,etal2,
     & mllmin,mllmax,phimax,mwindow,phi,etvec(2),missinget,ptll,ptllmin,
     & Etll,MET,mtrans
      logical first,l1muon,l2muon,highMT
      data first/.true./
      save first,highMT,mllmax,phimax
 
C---- optional cut on mtrans 
c      etvec(1)=p(3,1)+p(6,1)
c      etvec(2)=p(3,2)+p(6,2)
c--- transverse mass
c      Etll=dsqrt(max(0d0,(p(4,4)+p(5,4))**2-(p(4,3)+p(5,3))**2))
c      MET=dsqrt(max(0d0,etvec(1)**2+etvec(2)**2))
c      mtrans=MET+Etll

c      if(mtrans.lt.300d0) then 
c        failed_cuts=.true. 
c         return 
c      endif

c--- flag for (e,mu + mu,e) or (e,e) or (mu,mu)
c      idlept=1    ! (e,mu)
c      idlept=2    ! (e,e)
c      idlept=3    ! (mu,mu)
      
      idlept=1

      ptminl1=25d0
      ptminl2=15d0
      
      etamaxmu=2.5d0
      etamaxe=2.47d0
      etaegapmin=1.37d0
      etaegapmax=1.52d0
      
      ptllmin=30d0

      if     (idlept .eq. 1) then
        mllmin=10d0
        mwindow=0d0
        metrelmin=25d0
        if (ran2() .lt. 0.5d0) then
          l1muon=.true.
          l2muon=.false.
        else
          l1muon=.false.
          l2muon=.true.
        endif
      elseif (idlept .eq. 2) then
        mllmin=12d0
        mwindow=15d0
        metrelmin=45d0
        l1muon=.false.
        l2muon=.false.
      elseif (idlept .eq. 3) then
        mllmin=12d0
        mwindow=15d0
        metrelmin=45d0
        l1muon=.true.
        l2muon=.true.
      else
        write(6,*) 'idlept set wrong in ATLAS_hww2013: ',idlept
        stop
      endif

         
      if(first) then 
      mllmax=50d0
      phimax=1.8d0      
      if (index(runstring,'nomll') .gt. 0) then
        mllmax=1d4
      endif
      if (index(runstring,'basic') .gt. 0) then
        mllmax=1d4
        phimax=1d4
      endif      
      if (index(runstring,'highMT') .gt. 0) then
        highMT=.true.
      else
        highMT=.false.
      endif      
      first=.false.
      write(6,*)  '******** ATLAS H->WW Search cuts ******'
      write(6,*)  '*                                     *'
      if     (idlept .eq. 1) then
      write(6,*)  '*     Treating W decays as (e,mu)     *'
      elseif (idlept .eq. 2) then
      write(6,*)  '*     Treating W decays as (e,e)      *'
      elseif (idlept .eq. 3) then
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
        Etll=dsqrt(max(0d0,(p(4,4)+p(5,4))**2-(p(4,3)+p(5,3))**2))
        MET=dsqrt(max(0d0,etvec(1)**2+etvec(2)**2))
        mtrans=MET+Etll
        if (mtrans .lt. 250d0) then
          failed_cuts=.true.
          return
        endif
      endif

      mll=sqrt(max(0d0,2d0*dot(p,4,5)))

c--- minimum dilepton invariant mass
      if (mll .lt. mllmin) then
        failed_cuts=.true.
        return
      endif

c--- maximum dilepton invariant mass
      if (mll .gt. mllmax) then
        failed_cuts=.true.
        return
      endif

c--- veto on dilepton invariant mass around Z mass
      if (abs(mll-zmass) .lt. mwindow) then
        failed_cuts=.true.
        return
      endif

      ptl1=pt(4,p)
      ptl2=pt(5,p)
      
      if (ptl2 .gt. ptl1) then
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
      if (ptl1 .lt. ptminl1) then
        failed_cuts=.true.
        return
      endif 
c--- pt cut on softest lepton
      if (ptl2 .lt. ptminl2) then
        failed_cuts=.true.
        return
      endif 

c--- rapidity cut on hardest lepton
      if (l1muon) then
        if (abs(etal1) .gt. etamaxmu) then
          failed_cuts=.true.
          return
        endif
      else
        if ((abs(etal1) .gt. etamaxe) .or.
     &   (abs(etal1).gt.etaegapmin) .and.(abs(etal1).lt.etaegapmax))then
          failed_cuts=.true.
          return
        endif
      endif

c--- rapidity cut on softest lepton
      if (l2muon) then
        if (abs(etal2) .gt. etamaxmu) then
          failed_cuts=.true.
          return
        endif
      else
        if ((abs(etal2) .gt. etamaxe) .or.
     &   (abs(etal2).gt.etaegapmin) .and.(abs(etal2).lt.etaegapmax))then
          failed_cuts=.true.
          return
        endif
      endif

c--- missing et "relative"
      do j=1,2
        etvec(j)=p(3,j)+p(6,j)
      enddo
      missinget=dsqrt(etvec(1)**2+etvec(2)**2)

      phimin=99d0
      do i=4,5
        phi=(p(i,1)*etvec(1)+p(i,2)*etvec(2))
     &      /dsqrt((p(i,1)**2+p(i,2)**2)*missinget**2)
        if (phi .gt. +0.9999999D0) phi=+1D0
        if (phi .lt. -0.9999999D0) phi=-1D0
        phi=dacos(phi)
        if (phi .lt. phimin) phimin=phi
      enddo

      if (phi .lt. pi/2d0) then
        missinget=missinget*dsin(phi)
      endif

      if (missinget .lt. metrelmin) then
        failed_cuts=.true.
        return
      endif

c--- dilepton transverse momentum
      ptll=dsqrt((p(4,1)+p(5,1))**2+(p(4,2)+p(5,2))**2)
      if (ptll .lt. ptllmin) then
        failed_cuts=.true.
        return
      endif

c--- dilepton phi
      phi=(p(4,1)*p(5,1)+p(4,2)*p(5,2))
     &    /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (phi .gt. +0.9999999D0) phi=+1D0
      if (phi .lt. -0.9999999D0) phi=-1D0
      phi=dacos(phi)
      if (phi .gt. phimax) then
        failed_cuts=.true.
        return
      endif

 99   format(1x,a23,f6.2,a10)
      
      return 
      end


         

