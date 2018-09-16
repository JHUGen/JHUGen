      subroutine VBS(p,failed_cuts) 
c--- a generic set of vector boson scattering (VBS) cuts, based
c--- on the cuts used for like-sign WWjj production in
c--- 1405.6241 (ATLAS) and 1410.6315 (CMS)
      implicit none 
      include 'constants.f' 
      include 'masses.f'
      include 'runstring.f'
      include 'jetcuts.f'
      double precision p(mxpart,4)
      integer ijet,ilep,id(4),idj(4),j,i1(6),i2(6),iperms,j1,j2,l1,l2,
     & ilomomenta
      logical failed_cuts,is_lepton,is_hadronic,doszleper
      double precision pt,etarap,ptlep,etalep,mjj,etaj1,etaj2,mlj
      double precision ptparton(3:9),etaparton(3:9),etmiss,mljmin,
     & mll,mllmin,metvec(4),met,metmin,mjets,ygap,etajmin,etajmax
      logical first
      common/ilomomenta/ilomomenta
      data first/.true./
      data i1/1,1,2,1,2,3/
      data i2/2,3,3,4,4,4/
      save first,i1,i2,doszleper

      ptlep=20d0 
      etalep=2.5d0
      metmin=40d0
      mllmin=10d0
      mljmin=200d0

c--- additional cuts for VBF topology
      ygap=2.5d0
      mjets=500d0
         
      if(first) then 
      first=.false.

      write(6,*)  '************** VBS cuts ***************'
      write(6,*)  '*                                     *'
      write(6,99) '*     pt(lep) >          ', ptlep,         ' GeV   *'
      write(6,99) '*     |eta_lep| <        ', etalep,        '       *'
      write(6,99) '*     m(lep1,lep2) >     ',mllmin,         ' GeV   *'
      write(6,99) '*     MET >              ',metmin,         ' GeV   *'
      write(6,99) '*     |y(j1)-y(j2))| >   ', ygap,          '       *'
      write(6,99) '*     m(j1,j2) >         ', mjets,         ' GeV   *'
      write(6,99) '*     m(l1,j2),m(l2,j1)> ', mljmin,        ' GeV   *'
      write(6,*)  '*     + jets in opposite hemispheres  *'
      write(6,*)  '*     + lepton rap. between jets      *'
      write(6,*)  '***************************************'

      if(index(runstring,'szleper') .gt. 0) then
        doszleper=.true.
        write(6,*) '*         + Szleper Rpt > 2.0         *'
      else
        doszleper=.false.
      endif

      endif
 
      failed_cuts = .false. 

      do j=3,8
      ptparton(j)=pt(j,p)
      etaparton(j)=etarap(j,p)
      enddo
      
      if (ilomomenta .gt. 9) then
        write(6,*) 'VBS cuts routine does not accept > 9 momenta'
        stop
      endif
      
      ilep=0
      ijet=0
      do j=3,ilomomenta
        if (is_lepton(j)) then      ! lepton cuts
          ilep=ilep+1
          id(ilep)=j
          if (ptparton(j) .lt. ptlep)  then 
             failed_cuts=.true. 
             return 
          endif
          if (abs(etaparton(j)) .gt. etalep)  then 
             failed_cuts=.true. 
             return 
          endif
        elseif (is_hadronic(j)) then      ! jet cuts
          ijet=ijet+1
          idj(ijet)=j
        endif
      enddo
      
      if ((ilep .lt. 2) .or. (ilep .gt. 4)) then
        write(6,*) 'Error in cutting routine:'
        write(6,*) 'did not find 2 or 4 leptons: ',ilep
        stop
      endif

      if (ijet .lt. 2) then
        write(6,*) 'Error in cutting routine:'
        write(6,*) 'did not find 2 jets: ',ijet
        stop
      write(6,*) 'id',id
      write(6,*) 'idj',idj
      endif

c--- if there are three jets, determine two with largest gap
      if (ijet .eq. 3) then
        call etagapsort(idj,etaparton)
c--- veto remaining jet if in acceptance
c        if ((ptparton(idj(3)) .gt. ptjetmin) .and. 
c     &      (abs(etaparton(idj(3))) .lt. etajetmax)) then
c          failed_cuts=.true. 
c          return 
c        endif
      endif

c--- ensure a rapidity gap of at least ygap between the tagging jets
      etaj1=etaparton(idj(1))
      etaj2=etaparton(idj(2))      
      if (abs(etaj1-etaj2) .lt. ygap) then
         failed_cuts=.true. 
         return 
      endif
      
c--- ensure the tagging jets lie in opposite hemispheres
      if (etaj1*etaj2 .ge. 0d0) then
         failed_cuts=.true. 
         return 
      endif

c--- lepton rapidities between jets
      etajmin=min(etaj1,etaj2)
      etajmax=max(etaj1,etaj2)
      
      do j=1,ilep
      if ( (etaparton(id(j)) .lt. etajmin)
     & .or.(etaparton(id(j)) .gt. etajmax)) then
         failed_cuts=.true. 
         return 
      endif
      enddo
      
c--- ensure the tagging jets have an invariant mass larger than mjets
      mjj=dsqrt(max(0d0,(p(idj(1),4)+p(idj(2),4))**2
     & -(p(idj(1),1)+p(idj(2),1))**2
     & -(p(idj(1),2)+p(idj(2),2))**2
     & -(p(idj(1),3)+p(idj(2),3))**2))
      if (mjj .lt. mjets) then
         failed_cuts=.true. 
         return 
      endif
 
c--- minimum value of m_ll
      if     (ilep .eq. 2) then
        iperms=1
      elseif (ilep .eq. 3) then
        iperms=3
      elseif (ilep .eq. 4) then
        iperms=6
      endif
      do j=1,iperms
        mll=(p(id(i1(j)),4)+p(id(i2(j)),4))**2
     &     -(p(id(i1(j)),1)+p(id(i2(j)),1))**2
     &     -(p(id(i1(j)),2)+p(id(i2(j)),2))**2
     &     -(p(id(i1(j)),3)+p(id(i2(j)),3))**2
        if (mll .lt. mllmin**2) then
          failed_cuts=.true. 
          return 
        endif
      enddo

c--- only apply missing ET cut if <4 leptons
      if (ilep .lt. 4) then
        met=etmiss(p,metvec)
        if (met .lt. metmin) then
          failed_cuts=.true.
          return
        endif
      endif

c--- cut on invariant mass of m(l1,j2) and (ml2,j1)
c--- where (j1,j2) and (l1,l2) are pt-ordered
      if (ilep .eq. 2) then
        if (ptparton(id(1)) .gt. ptparton(id(2))) then
          l1=id(1)
          l2=id(2)
        else
          l1=id(2)
          l2=id(1)
        endif
        if (ptparton(idj(1)) .gt. ptparton(idj(2))) then
          j1=idj(1)
          j2=idj(2)
        else
          j1=idj(2)
          j2=idj(1)
        endif
        mlj=(p(l1,4)+p(j2,4))**2
     &     -(p(l1,1)+p(j2,1))**2
     &     -(p(l1,2)+p(j2,2))**2
     &     -(p(l1,3)+p(j2,3))**2
        if (mlj .lt. mljmin**2) then
          failed_cuts=.true. 
          return 
        endif
        mlj=(p(l2,4)+p(j1,4))**2
     &     -(p(l2,1)+p(j1,1))**2
     &     -(p(l2,2)+p(j1,2))**2
     &     -(p(l2,3)+p(j1,3))**2
        if (mlj .lt. mljmin**2) then
          failed_cuts=.true. 
          return 
        endif
      endif


c--- Szleper Rpt cut
      if (doszleper) then
        if (ptparton(id(1))*ptparton(id(2))
     &     /(ptparton(idj(1))*ptparton(idj(2))) .lt. 2.0d0) then
          failed_cuts=.true.
          return
         endif
      endif
      
 99   format(1x,a25,f6.2,a8)
      
      return 
      end
      
      
      
      subroutine etagapsort(idj,etaparton)
c--- given array ptparton, sort indices idj such that
c--- |eta(idj(1))-eta(idj(2))|>|eta(idj(1))-eta(idj(3))|>|eta(idj(2))-eta(idj(3))| 
      implicit none
      integer idj(4),idjout(4)
      double precision etaparton(3:9),etaj1,etaj2,etaj3
      
      idjout(:)=idj(:)

      etaj1=etaparton(idj(1))
      etaj2=etaparton(idj(2))
      etaj3=etaparton(idj(3))
      
      if(abs(etaj1-etaj2).gt.max(abs(etaj1-etaj3),abs(etaj1-etaj3)))then
        idjout(1)=idj(1)
        idjout(2)=idj(2)
        idjout(3)=idj(3)
      else
     &if(abs(etaj1-etaj3).gt.max(abs(etaj1-etaj2),abs(etaj2-etaj3)))then
        idjout(1)=idj(1)
        idjout(2)=idj(3)
        idjout(3)=idj(2)
      else
        idjout(1)=idj(2)
        idjout(2)=idj(3)
        idjout(3)=idj(1)
      endif

c--- return in-place
      idj(:)=idjout(:)
      
      return
      end
      
