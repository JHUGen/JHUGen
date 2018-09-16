      subroutine cms_higgsWW(p,failed_cuts)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: p(mxpart,4)
      logical:: failed_cuts
      real(dp):: pt,etarap,m45
      real(dp):: pt_l_hard,pt_l_soft
      real(dp):: pts,pth,phimax,mllmax,etmiss,etmiss_min
      real(dp):: et_vec(4),r2,delphi,eta_max_h
      real(dp):: eta_max_s,eta_hard,eta_soft
      integer:: i
      logical:: first
      data first/.true./
      common/HWW_Cuts/pts,pth,mllmax,phimax,etmiss_min,
     &     eta_max_h,eta_max_s


!------- ETA DOESNT CHANGE -------------
c      etmiss_min=20._dp
c      eta_max_h=2.5_dp
c      eta_max_s=2.5_dp
!----------------------------------------

      !-- cut params change to change cuts
!      pts=25._dp
!      pth=30._dp
!      mllmax=50._dp
!      phimax=60._dp

!------- ETA DOESNT CHANGE -------------
      etmiss_min=20._dp
      eta_max_h=2.5_dp
      eta_max_s=2.5_dp
!----------------------------------------

!      if (hmass < 160._dp) then
!        pth=25._dp
!        pts=20._dp
!        mllmax=45._dp
!        phimax=60._dp
!      endif

!      if (hmass > 180._dp) then
!        pth=40._dp
!        pts=25._dp
!        mllmax=90._dp
!        phimax=100._dp
!      endif

!      if (hmass > 380._dp) then
!        pth=90._dp
!        pts=25._dp!
!        mllmax=300._dp
!        phimax=175._dp
!      endif

      call HWW_cuts_params(hmass)

      if(first) then
         first=.false.
      write(6,*)  '**************** Higgs Search cuts  ****************'
      write(6,*)  '*                                                  *'
      write(6,99)  '*     pt(lep)_max>  ', pth, '                    *'
      write(6,99)  '*     pt(lep)_min >  ', pts, '                    *'
      write(6,99)  '*     mll  <  ', mllmax, '                     *'
      write(6,99)  '*     phi(lep,lep) <  ', phimax, '               *'
      write(6,99)  '*     |eta_l | <  ', eta_max_h, '                 *'
      write(6,99) '*       ET_miss > ',etmiss_min, '                *'
      write(6,*)  '****************************************************'
      endif


!---- Subroutine for calculting CMS Higgs (m_H=160) cuts namely
!---- pt_l(max) > 30 pt_l(min) > 25 m_ll < 50 delta_phi_ll < 60 degrees
!---- eta_l < 2.5, ET_miss = 20
!-----Inclusive in the jet
! '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))'
c      if(kcase.ne.kWWqqbr) then
c         write(6,*)'Attempted to apply Higgs search cuts to process'
c         write(6,*)'Check Runstring '
c         stop
c         return
c      endif


      failed_cuts = .false.
      do i=1,4
         et_vec(i)=zip
      enddo

      if(pt(4,p)>pt(5,p)) then
         pt_l_hard=pt(4,p)
         eta_hard=etarap(4,p)
         pt_l_soft=pt(5,p)
         eta_soft=etarap(5,p)
      else
         pt_l_hard=pt(5,p)
         eta_hard=etarap(5,p)
         pt_l_soft=pt(4,p)
         eta_soft=etarap(4,p)
      endif



!---- pt cuts
      if((pt_l_hard<pth).or.(pt_l_soft<pts)) then
         failed_cuts=.true.
         return
      endif

!Binoth et al cuts
!      if((pt_l_hard>50._dp).or.(pt_l_hard<35._dp)) then
!         failed_cuts=.true.
!         return
!      elseif((pt_l_soft>25._dp).or.(pt_l_soft<20._dp)) then
!         failed_cuts=.true.
!         return
!      endif

      if(abs(etmiss(p,et_vec))<etmiss_min) then
         failed_cuts=.true.
         return
      endif

!---- eta_cuts
      if((abs(eta_hard)>eta_max_h).or.
     &     ((abs(eta_soft)>eta_max_s))) then
         failed_cuts=.true.
!         write(6,*) eta_hard,eta_soft
         return
      endif


!--- m_ll cut
      m45=zip
      do i=1,4
         if(i.ne.4) then
            m45=m45-(p(4,i)+p(5,i))**2
         else
            m45=m45+(p(4,i)+p(5,i))**2
         endif
      enddo

      if(sqrt(max(m45,zip))>mllmax) then
         failed_cuts=.true.
         return
      endif

!---- phi_ll
!--- convert phi max to radians
      phimax=phimax/(360._dp)*2._dp*pi

       r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     &     /sqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (r2 > +0.9999999_dp) r2=+one
      if (r2 < -0.9999999_dp) r2=-one
      delphi=acos(r2)
      if(delphi>phimax) then
         failed_cuts=.true.
         return
      endif

 99   format(1x,a29,f6.2,a17)

      return
      end


      subroutine HWW_cuts_params(in_mh)
      implicit none
      include 'types.f'

      real(dp):: mh
      real(dp):: pts,pth,mllmax,phimax
      real(dp):: etmiss_min,eta_max_h,eta_max_s
      common/HWW_Cuts/pts,pth,mllmax,phimax,etmiss_min,
     &     eta_max_h,eta_max_s
      real(dp):: tiny,in_mh
      tiny=1.e-4_dp
!----- Ensure that an input of 130 etc is treated properly (i.e. always boosted into correct cuts
      mh=in_mh+tiny

!---- HWW_cuts from CMS 7 TeV
      if(mh<130._dp) then
         pth=20._dp
         pts=20._dp
         mllmax=40._dp
         phimax=60._dp
      elseif((mh>130._dp).and.(mh<150._dp)) then
         pth=25._dp
         pts=20._dp
         mllmax=45._dp
         phimax=60._dp

      elseif((mh>150._dp).and.(mh<160._dp)) then
         pth=27._dp
         pts=25._dp
         mllmax=50._dp
         phimax=60._dp


      elseif((mh>160._dp).and.(mh<170._dp)) then
         pth=30._dp
         pts=25._dp
         mllmax=50._dp
         phimax=60._dp


      elseif((mh>170._dp).and.(mh<180._dp)) then
         pth=34._dp
         pts=25._dp
         mllmax=50._dp
         phimax=60._dp

      elseif((mh>180._dp).and.(mh<190._dp)) then
         pth=36._dp
         pts=25._dp
         mllmax=60._dp
         phimax=70._dp

      elseif((mh>190._dp).and.(mh<200._dp)) then
         pth=38._dp
         pts=25._dp
         mllmax=80._dp
         phimax=90._dp

      elseif((mh>200._dp).and.(mh<210._dp)) then
         pth=40._dp
         pts=25._dp
         mllmax=90._dp
         phimax=100._dp

      elseif((mh>210._dp).and.(mh<220._dp)) then
         pth=44._dp
         pts=25._dp
         mllmax=110._dp
         phimax=110._dp

      elseif((mh>220._dp).and.(mh<230._dp)) then
         pth=48._dp
         pts=25._dp
         mllmax=120._dp
         phimax=120._dp

      elseif((mh>230._dp).and.(mh<250._dp)) then
         pth=52._dp
         pts=25._dp
         mllmax=130._dp
         phimax=130._dp

      elseif((mh>250._dp).and.(mh<300._dp)) then
         Pth=55._dp
         pts=25._dp
         mllmax=150._dp
         phimax=140._dp

      elseif((mh>300._dp)) then
         pth=70._dp
         pts=25._dp
         mllmax=200._dp
         phimax=175._dp

      endif
      return

      end
