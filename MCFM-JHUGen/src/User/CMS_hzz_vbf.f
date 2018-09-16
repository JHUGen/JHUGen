      subroutine CMS_hzz_vbf(p,failed_cuts) 
c--- Same cuts as CMS_hzz, with the addition of a few to select VBF topology
      implicit none 
      include 'constants.f' 
      include 'masses.f'
      double precision p(mxpart,4)
      logical failed_cuts,is_electron,is_muon 
      double precision pt,etarap,ptmu,etamu,ptel,etael,m4lmin,m2lmin,
     & etaj1,etaj2,m78,etajmin,etajmax
      double precision ptlep(3:6),etalep(3:6),mllmax,  
     & mllnearmin,mllfarmin,
     & orderedpt(4),ptleading,ptnexttoleading,m34,m56,m36,m45,m3456,
     & ygap,mjets
      integer i,j
      logical first 
      data first/.true./

      ptmu=5d0 
      ptel=7d0 
      etamu=2.4d0
      etael=2.5d0
      m4lmin=100d0
      m2lmin=4d0
      ptleading=20d0
      ptnexttoleading=10d0
      mllnearmin=40d0
      mllfarmin=12d0
      mllmax=120d0

c--- additional cuts for VBF topology
      ygap=2.5d0
      mjets=500d0
         
      if(first) then 
         first=.false.
      write(6,*)  '*********** CMS Higgs Search cuts******'
      write(6,*)  '*              + VBF topology         *'
      write(6,*)  '*                                     *'
      write(6,99) '*     pt(mu) >         ', ptmu,          ' GeV     *'
      write(6,99) '*     pt(el) >         ', ptel,          ' GeV     *'
      write(6,99) '*     |eta_mu| <       ', etamu,         '         *'
      write(6,99) '*     |eta_el| <       ', etael,         '         *'
      write(6,99) '*     m4lmin >         ',m4lmin,         ' GeV     *'
      write(6,99) '*     m2lmin >         ',m2lmin,         ' GeV     *'
      write(6,99) '*     ptleading >      ',ptleading,      ' GeV     *'
      write(6,99) '*     ptnexttoleading >',ptnexttoleading,' GeV     *'
      write(6,99) '*     mllnear >        ',mllnearmin,     ' GeV     *'
      write(6,99) '*     mllfar >         ',mllfarmin,      ' GeV     *'
      write(6,99) '*     mllmax >         ',mllmax,         ' GeV     *'
      write(6,*)  '*                                     *'
      write(6,99) '*     |y(j1)-y(j2))| > ', ygap, '     *'
      write(6,99) '*     m(j1,j2) >       ', mjets,' GeV *'
      write(6,*)  '*     + jets in opposite hemispheres  *'
      write(6,*)  '*     + lept rapidities between jets  *'
      write(6,*)  '***************************************'
      endif
      
      failed_cuts = .false. 

c--- ensure a rapidity gap of at least ygap between the tagging jets
      etaj1=etarap(7,p)
      etaj2=etarap(8,p)      
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
      m78=dsqrt(max(0d0,(p(7,4)+p(8,4))**2
     & -(p(7,1)+p(8,1))**2-(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2))
      if (m78 .lt. mjets) then
         failed_cuts=.true. 
         return 
      endif

      do j=3,6
      ptlep(j)=pt(j,p)
      etalep(j)=etarap(j,p)
      orderedpt(j-2)=ptlep(j)
      enddo
      
      do j=3,6
      if (is_electron(j)) then      ! electron cuts
      if (ptlep(j) .lt. ptel)  then 
         failed_cuts=.true. 
         return 
      endif
      if (abs(etalep(j)) .gt. etael)  then 
         failed_cuts=.true. 
         return 
      endif
      elseif (is_muon(j)) then      ! muon cuts
      if (ptlep(j) .lt. ptmu)  then 
         failed_cuts=.true. 
         return 
      endif
      if (abs(etalep(j)) .gt. etamu)  then 
         failed_cuts=.true. 
         return 
      endif
      endif
      enddo

c--- lepton rapidities between jets
      etajmin=min(etaj1,etaj2)
      etajmax=max(etaj1,etaj2)
      
      if ((etalep(3) .lt. etajmin) .or. (etalep(3) .gt. etajmax)) then
         failed_cuts=.true. 
         return 
      endif
      if ((etalep(4) .lt. etajmin) .or. (etalep(4) .gt. etajmax)) then
         failed_cuts=.true. 
         return 
      endif
      if ((etalep(5) .lt. etajmin) .or. (etalep(5) .gt. etajmax)) then
         failed_cuts=.true. 
         return 
      endif
      if ((etalep(6) .lt. etajmin) .or. (etalep(6) .gt. etajmax)) then
         failed_cuts=.true. 
         return 
      endif
      
      call piksrt(4,orderedpt)

      if (orderedpt(4) .lt. ptleading)  then 
         failed_cuts=.true. 
         return 
      endif
      
      if (orderedpt(3) .lt. ptnexttoleading)  then 
         failed_cuts=.true. 
         return 
      endif
      


!--- m_ll cut 
      i=4
      m34=(p(3,i)+p(4,i))**2
      m56=(p(5,i)+p(6,i))**2
      m36=(p(3,i)+p(6,i))**2
      m45=(p(4,i)+p(5,i))**2
      m3456=(p(3,i)+p(4,i)+p(5,i)+p(6,i))**2
      do i=1,3
            m34=m34-(p(3,i)+p(4,i))**2 
            m56=m56-(p(5,i)+p(6,i))**2 
            m36=m36-(p(3,i)+p(6,i))**2 
            m45=m45-(p(4,i)+p(5,i))**2 
            m3456=m3456-(p(3,i)+p(4,i)+p(5,i)+p(6,i))**2 
      enddo
      m34=sqrt(max(m34,0d0))
      m56=sqrt(max(m56,0d0))
      m36=sqrt(max(m36,0d0))
      m45=sqrt(max(m45,0d0))
      m3456=sqrt(max(m3456,0d0))

      if(m34 .lt. m2lmin) then 
         failed_cuts=.true.
         return 
      endif

      if(m56  .lt. m2lmin) then 
         failed_cuts=.true.
         return 
      endif

      if(m36  .lt. m2lmin) then 
         failed_cuts=.true.
         return 
      endif

      if(m45  .lt. m2lmin) then 
         failed_cuts=.true.
         return 
      endif

      if(m3456 .lt. m4lmin) then 
         failed_cuts=.true.
         return 
      endif

      if (abs(m34-zmass) .lt. abs(m56-zmass)) then
      if ((m34 .lt. mllnearmin) .or. (m34 .gt. mllmax)) then
         failed_cuts=.true.
         return 
      endif
      if ((m56 .lt. mllfarmin) .or. (m56 .gt. mllmax)) then
         failed_cuts=.true.
         return 
      endif
      elseif (abs(m34-zmass) .ge. abs(m56-zmass)) then
      if ((m56 .lt. mllnearmin) .or. (m56 .gt. mllmax)) then
         failed_cuts=.true.
         return 
      endif
      if ((m34 .lt. mllfarmin) .or. (m34 .gt. mllmax)) then
         failed_cuts=.true.
         return 
      endif
      endif

 99   format(1x,a23,f6.2,a10)
      
      return 
      end


