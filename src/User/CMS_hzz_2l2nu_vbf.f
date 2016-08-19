      subroutine CMS_hzz_2l2nu_vbf(p,failed_cuts) 
      implicit none 
      include 'constants.f' 
      include 'masses.f'
      include 'runstring.f'
      double precision p(mxpart,4)
      integer ijet,ilep,id(4),idj(4),j
      logical failed_cuts,is_electron,is_muon,is_hadronic 
      double precision pt,etarap,ptmu,etamu,ptel,etael,mjj,etaj1,etaj2
      double precision ptparton(3:8),etaparton(3:8),mllmax,etmiss,
     & mll,mllwidth,metvec(4),met,metmin,mjets,
     & ptleading,ptnexttoleading,etajmax,etajmin,
     & ygap

      logical first,nomll
      data first/.true./
      save first,nomll

      ptmu=20d0 
      ptel=20d0 
      etamu=2.4d0
      etael=2.4d0
      metmin=80d0 ! this is the CMS value
c      metmin=25d0 ! this is me relaxing the cut to see the H peak
      mllwidth=30d0

      ptleading=20d0
      ptnexttoleading=10d0

c--- additional cuts for VBF topology
      ygap=2.4d0
      mjets=500d0

         
      if(first) then 
         first=.false.

        if (index(runstring,'nomll') .gt. 0) then
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
      write(6,99) '*     |mll-mz| <       ',mllwidth/2d0,   ' GeV     *'
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
          if (ptparton(j) .lt. ptel)  then 
             failed_cuts=.true. 
             return 
          endif
          if (abs(etaparton(j)) .gt. etael)  then 
             failed_cuts=.true. 
             return 
          endif

        elseif (is_muon(j)) then      ! muon cuts
          ilep=ilep+1
          id(ilep)=j
          if (ptparton(j) .lt. ptmu)  then 
             failed_cuts=.true. 
             return 
          endif
          if (abs(etaparton(j)) .gt. etamu)  then 
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

      if (ijet .lt. 2) then
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
      if (abs(etaj1-etaj2) .lt. ygap) then
         failed_cuts=.true. 
         return 
      endif
      
c--- ensure the tagging jets lie in opposite hemispheres
      if (etaj1*etaj2 .ge. 0d0) then
         failed_cuts=.true. 
         return 
      endif

c--- ensure the tagging jets have an invariant mass larger than mjets
      mjj=dsqrt(max(0d0,(p(idj(1),4)+p(idj(2),4))**2
     & -(p(idj(1),1)+p(idj(2),1))**2
     & -(p(idj(1),2)+p(idj(2),2))**2
     & -(p(idj(1),3)+p(idj(2),3))**2))
      if (mjj .lt. mjets) then
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
      mll=sqrt(max(mll,0d0))
      
      if (abs(mll-zmass) .gt. mllwidth/2d0) then
        failed_cuts=.true. 
        return 
      endif
      
      endif

c--- missing ET cut      
      met=etmiss(p,metvec)
      if (met .lt. metmin) then
        failed_cuts=.true.
        return
      endif
      
 99   format(1x,a23,f6.2,a10)
      
      return 
      end

