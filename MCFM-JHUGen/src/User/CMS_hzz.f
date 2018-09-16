      subroutine CMS_hzz(p,failed_cuts)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: p(mxpart,4)
      logical:: failed_cuts,is_electron,is_muon
      real(dp):: pt,etarap,ptmu,etamu,ptel,etael,m4lmin,m2lmin
      real(dp):: ptlep(3:6),etalep(3:6),mllmax,
     & mllnearmin,mllfarmin,
     & orderedpt(4),ptleading,ptnexttoleading,m34,m56,m36,m45,m3456
      integer:: i,j
      logical:: first
      data first/.true./

      ptmu=5._dp
      ptel=7._dp
      etamu=2.4_dp
      etael=2.5_dp
      m4lmin=100._dp
      m2lmin=4._dp
      ptleading=20._dp
      ptnexttoleading=10._dp
      mllnearmin=40._dp
      mllfarmin=12._dp
      mllmax=120._dp


      if(first) then
         first=.false.
      write(6,*)  '*********** CMS Higgs Search cuts******'
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
      write(6,*)  '***************************************'
      endif

      failed_cuts = .false.

      do j=3,6
      ptlep(j)=pt(j,p)
      etalep(j)=etarap(j,p)
      orderedpt(j-2)=ptlep(j)
      enddo

      do j=3,6
      if (is_electron(j)) then      ! electron cuts
      if (ptlep(j) < ptel)  then
         failed_cuts=.true.
         return
      endif
      if (abs(etalep(j)) > etael)  then
         failed_cuts=.true.
         return
      endif
      elseif (is_muon(j)) then      ! muon cuts
      if (ptlep(j) < ptmu)  then
         failed_cuts=.true.
         return
      endif
      if (abs(etalep(j)) > etamu)  then
         failed_cuts=.true.
         return
      endif
      endif
      enddo

      call piksrt(4,orderedpt)

      if (orderedpt(4) < ptleading)  then
         failed_cuts=.true.
         return
      endif

      if (orderedpt(3) < ptnexttoleading)  then
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
      m34=sqrt(max(m34,zip))
      m56=sqrt(max(m56,zip))
      m36=sqrt(max(m36,zip))
      m45=sqrt(max(m45,zip))
      m3456=sqrt(max(m3456,zip))

      if(m34 < m2lmin) then
         failed_cuts=.true.
         return
      endif

      if(m56  < m2lmin) then
         failed_cuts=.true.
         return
      endif

      if(m36  < m2lmin) then
         failed_cuts=.true.
         return
      endif

      if(m45  < m2lmin) then
         failed_cuts=.true.
         return
      endif

      if(m3456 < m4lmin) then
         failed_cuts=.true.
         return
      endif

      if (abs(m34-zmass) < abs(m56-zmass)) then
      if ((m34 < mllnearmin) .or. (m34 > mllmax)) then
         failed_cuts=.true.
         return
      endif
      if ((m56 < mllfarmin) .or. (m56 > mllmax)) then
         failed_cuts=.true.
         return
      endif
      elseif (abs(m34-zmass) >= abs(m56-zmass)) then
      if ((m56 < mllnearmin) .or. (m56 > mllmax)) then
         failed_cuts=.true.
         return
      endif
      if ((m34 < mllfarmin) .or. (m34 > mllmax)) then
         failed_cuts=.true.
         return
      endif
      endif

 99   format(1x,a23,f6.2,a10)

      return
      end




      SUBROUTINE piksrt(n,arr)
      include 'types.f'
      integer:: n
      real(dp):: arr(n)
      integer:: i,j
      real(dp):: a
      do 12 j=2,n
        a=arr(j)
        do 11 i=j-1,1,-1
          if(arr(i)<=a)goto 10
          arr(i+1)=arr(i)
 11             continue
        i=0
 10           arr(i+1)=a
 12               continue
      return
      END
