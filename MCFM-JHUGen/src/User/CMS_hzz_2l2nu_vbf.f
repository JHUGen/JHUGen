      subroutine CMS_hzz_2l2nu_vbf(p,failed_cuts)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'runstring.f'
      real(dp):: p(mxpart,4)
      integer:: ijet,ilep,id(4),idj(4),j
      logical:: failed_cuts,is_electron,is_muon,is_hadronic
      real(dp):: pt,etarap,ptmu,etamu,ptel,etael,mjj,etaj1,etaj2
      real(dp):: ptparton(3:8),etaparton(3:8),mllmax,etmiss,
     & mll,mllwidth,metvec(4),met,metmin,mjets,
     & ptleading,ptnexttoleading,etajmax,etajmin,
     & ygap

      logical:: first,nomll
      data first/.true./
      save first,nomll

      ptmu=20._dp
      ptel=20._dp
      etamu=2.4_dp
      etael=2.4_dp
      metmin=80._dp ! this is the CMS value
c      metmin=25._dp ! this is me relaxing the cut to see the H peak
      mllwidth=30._dp

      ptleading=20._dp
      ptnexttoleading=10._dp

c--- additional cuts for VBF topology
      ygap=2.4_dp
      mjets=500._dp


      if(first) then
         first=.false.

        if (index(runstring,'nomll') > 0) then
          nomll=.true.
        else
          nomll=.false.
        endif

      write(6,*)  '********** CMS (2l,2nu) cuts **********'
      write(6,*)  '*              + VBF topology         *'
      write(6,*)  '*                                     *'
      write(6,99) '*     pt(mu) >         ', ptmu,          ' GeV     *'
      write(6,99) '*     pt(el) >         ', ptel,          ' GeV     *'
      write(6,99) '*     |eta_mu| <       ', etamu,         '         *'
      write(6,99) '*     |eta_el| <       ', etael,         '         *'
      write(6,99) '*     MET >            ',metmin,         ' GeV     *'
      if (nomll .eqv. .false.) then
      write(6,99) '*     |mll-mz| <       ',mllwidth/2._dp,   ' GeV     *'
      endif
      write(6,99) '*     |y(j1)-y(j2))| > ', ygap, '     *'
      write(6,99) '*     m(j1,j2) >       ', mjets,' GeV *'
      write(6,*)  '*     + jets in opposite hemispheres  *'
      write(6,*)  '*     + lept rapidities between jets  *'
      write(6,*)  '***************************************'

      endif

      failed_cuts = .false.

      do j=3,8
      ptparton(j)=pt(j,p)
      etaparton(j)=etarap(j,p)
      enddo

      ilep=0
      ijet=0
      do j=3,8
        if (is_electron(j)) then      ! electron cuts
          ilep=ilep+1
          id(ilep)=j
          if (ptparton(j) < ptel)  then
             failed_cuts=.true.
             return
          endif
          if (abs(etaparton(j)) > etael)  then
             failed_cuts=.true.
             return
          endif

        elseif (is_muon(j)) then      ! muon cuts
          ilep=ilep+1
          id(ilep)=j
          if (ptparton(j) < ptmu)  then
             failed_cuts=.true.
             return
          endif
          if (abs(etaparton(j)) > etamu)  then
             failed_cuts=.true.
             return
          endif

        elseif (is_hadronic(j)) then      ! jet cuts
          ijet=ijet+1
          idj(ijet)=j
        endif
      enddo

      if (ilep .ne. 2) then
        write(6,*) 'Error in cutting routine:'
        write(6,*) 'did not find 2 leptons: ',ilep
        stop
      endif

      if (ijet < 2) then
        write(6,*) 'Error in cutting routine:'
        write(6,*) 'did not find 2 jets: ',ijet
        stop
      write(6,*) 'id',id
      write(6,*) 'idj',idj
      endif
c      goto 101
c--- ensure a rapidity gap of at least ygap between the tagging jets
      etaj1=etaparton(idj(1))
      etaj2=etaparton(idj(2))
      if (abs(etaj1-etaj2) < ygap) then
         failed_cuts=.true.
         return
      endif

c--- ensure the tagging jets lie in opposite hemispheres
      if (etaj1*etaj2 >= zip) then
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


c--- lepton rapidities between jets
      etajmin=min(etaj1,etaj2)
      etajmax=max(etaj1,etaj2)

      if  ((etaparton(id(1)) < etajmin)
     & .or.(etaparton(id(1)) > etajmax)) then
      failed_cuts=.true.
         return
      endif
      if  ((etaparton(id(2)) < etajmin)
     & .or.(etaparton(id(2)) > etajmax)) then
      failed_cuts=.true.
         return
      endif

 101  continue

      if (nomll .eqv. .false.) then
c--- m_ll cut
      mll=(p(id(1),4)+p(id(2),4))**2
     &   -(p(id(1),1)+p(id(2),1))**2
     &   -(p(id(1),2)+p(id(2),2))**2
     &   -(p(id(1),3)+p(id(2),3))**2
      mll=sqrt(max(mll,zip))

      if (abs(mll-zmass) > mllwidth/2._dp) then
        failed_cuts=.true.
        return
      endif

      endif

c--- missing ET cut
      met=etmiss(p,metvec)
      if (met < metmin) then
        failed_cuts=.true.
        return
      endif

 99   format(1x,a23,f6.2,a10)

      return
      end

