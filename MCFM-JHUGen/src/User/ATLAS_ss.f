      subroutine ATLAS_ss(p,failed_cuts)
      implicit none
      include 'types.f'
c--- a set of vector boson scattering (VBS) cuts, to mimic
c--- the cuts used for like-sign WWjj production in 1405.6241 (ATLAS)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'runstring.f'
      include 'jetcuts.f'
      real(dp):: p(mxpart,4)
      integer:: ijet,ilep,id(4),idj(4),j,i1(6),i2(6),iperms,j1,j2,l1,l2,
     & ilomomenta,i
      logical:: failed_cuts,is_lepton,is_hadronic,doszleper
      real(dp):: pt,etarap,ptlep,etalep,mjj,etaj1,etaj2,mlj
      real(dp):: ptparton(3:9),etaparton(3:9),etmiss,mljmin,
     & mll,mllmin,metvec(4),met,metmin,mjets,ygap,etajmin,etajmax,
     & Rllmin,Rljmin,R
      logical:: first
      common/ilomomenta/ilomomenta
      data first/.true./
      data i1/1,1,2,1,2,3/
      data i2/2,3,3,4,4,4/
      save first,i1,i2,doszleper

      ptlep=25._dp
      etalep=2.5_dp
      metmin=40._dp
      mllmin=20._dp
      Rllmin=0.3_dp
      Rljmin=0.3_dp

c--- additional cuts for VBF topology
      ygap=2.4_dp
      mjets=500._dp

      if(first) then
      first=.false.

      write(6,*)  '************** VBS cuts ***************'
      write(6,*)  '*                                     *'
      write(6,99) '*     pt(lep) >          ', ptlep,         ' GeV   *'
      write(6,99) '*     |eta_lep| <        ', etalep,        '       *'
      write(6,99) '*     m(lep1,lep2) >     ',mllmin,         ' GeV   *'
      write(6,99) '*     R(lep1,lep2) >     ', Rllmin,        '       *'
      write(6,99) '*     R(lep,jet) >       ', Rljmin,        '       *'
      write(6,99) '*     MET >              ',metmin,         ' GeV   *'
      write(6,99) '*     |y(j1)-y(j2))| >   ', ygap,          '       *'
      write(6,99) '*     m(j1,j2) >         ', mjets,         ' GeV   *'
      write(6,*)  '***************************************'

      if(index(runstring,'szleper') > 0) then
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

      if (ilomomenta > 9) then
        write(6,*) 'VBS cuts routine does not accept > 9 momenta'
        stop
      endif

      ilep=0
      ijet=0
      do j=3,ilomomenta
        if (is_lepton(j)) then      ! lepton cuts
          ilep=ilep+1
          id(ilep)=j
          if (ptparton(j) < ptlep)  then
             failed_cuts=.true.
             return
          endif
          if (abs(etaparton(j)) > etalep)  then
             failed_cuts=.true.
             return
          endif
        elseif (is_hadronic(j)) then      ! jet cuts
          ijet=ijet+1
          idj(ijet)=j
        endif
      enddo

      if ((ilep < 2) .or. (ilep > 4)) then
        write(6,*) 'Error in cutting routine:'
        write(6,*) 'did not find 2 or 4 leptons: ',ilep
        stop
      endif

      if (ijet < 2) then
        write(6,*) 'Error in cutting routine:'
        write(6,*) 'did not find 2 jets: ',ijet
        stop
      write(6,*) 'id',id
      write(6,*) 'idj',idj
      endif

c--- if there are three jets, determine two with largest gap
      if (ijet == 3) then
        call etagapsort(idj,etaparton)
c--- veto remaining jet if in acceptance
c        if ((ptparton(idj(3)) > ptjetmin) .and.
c     &      (abs(etaparton(idj(3))) < etajetmax)) then
c          failed_cuts=.true.
c          return
c        endif
      endif

c--- ensure a rapidity gap of at least ygap between the tagging jets
      etaj1=etaparton(idj(1))
      etaj2=etaparton(idj(2))
      if (abs(etaj1-etaj2) < ygap) then
         failed_cuts=.true.
         return
      endif

c--- ensure the tagging jets have an invariant mass larger than mjets
      mjj=sqrt(max(zip,(p(idj(1),4)+p(idj(2),4))**2
     & -(p(idj(1),1)+p(idj(2),1))**2
     & -(p(idj(1),2)+p(idj(2),2))**2
     & -(p(idj(1),3)+p(idj(2),3))**2))
      if (mjj < mjets) then
         failed_cuts=.true.
         return
      endif

c--- minimum value of m_ll
      if     (ilep == 2) then
        iperms=1
      elseif (ilep == 3) then
        iperms=3
      elseif (ilep == 4) then
        iperms=6
      endif
      do j=1,iperms
        mll=(p(id(i1(j)),4)+p(id(i2(j)),4))**2
     &     -(p(id(i1(j)),1)+p(id(i2(j)),1))**2
     &     -(p(id(i1(j)),2)+p(id(i2(j)),2))**2
     &     -(p(id(i1(j)),3)+p(id(i2(j)),3))**2
        if (mll < mllmin**2) then
          failed_cuts=.true.
          return
        endif
      enddo

c--- only apply missing ET cut if <4 leptons
      if (ilep < 4) then
        met=etmiss(p,metvec)
        if (met < metmin) then
          failed_cuts=.true.
          return
        endif
      endif

c--- cut on Delta-R(l,l)
      do i=1,ilep
      do j=i+1,ilep
      if (R(p,id(i),id(j)) < Rllmin) then
        failed_cuts=.true.
        return
      endif
      enddo
      enddo

c--- cut on Delta-R(l,j)
      do i=1,ilep
      do j=1,ijet
      if (R(p,id(i),idj(j)) < Rljmin) then
        failed_cuts=.true.
        return
      endif
      enddo
      enddo

 99   format(1x,a25,f6.2,a8)

      return
      end
