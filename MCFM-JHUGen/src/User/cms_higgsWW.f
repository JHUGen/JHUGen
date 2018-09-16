      subroutine cms_higgsWW(p,failed_cuts) 
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision p(mxpart,4)
      logical failed_cuts 
      double precision pt,etarap,m45
      double precision pt_l_hard,pt_l_soft
      double precision pts,pth,phimax,mllmax,etmiss,etmiss_min
      double precision et_vec(4),r2,delphi,eta_max_h
      double precision eta_max_s,eta_hard,eta_soft
      integer i 
      logical first 
      data first/.true./
      common/HWW_Cuts/pts,pth,mllmax,phimax,etmiss_min,
     &     eta_max_h,eta_max_s


!------- ETA DOESNT CHANGE ------------- 
c      etmiss_min=20d0
c      eta_max_h=2.5d0
c      eta_max_s=2.5d0
!----------------------------------------   

      !-- cut params change to change cuts
!      pts=25d0
!      pth=30d0
!      mllmax=50d0
!      phimax=60d0

!------- ETA DOESNT CHANGE ------------- 
      etmiss_min=20d0
      eta_max_h=2.5d0
      eta_max_s=2.5d0
!----------------------------------------      

!      if (hmass .lt. 160d0) then
!        pth=25d0
!        pts=20d0
!        mllmax=45d0
!        phimax=60d0
!      endif
      
!      if (hmass .gt. 180d0) then
!        pth=40d0
!        pts=25d0
!        mllmax=90d0
!        phimax=100d0
!      endif
      
!      if (hmass .gt. 380d0) then
!        pth=90d0
!        pts=25d0!
!        mllmax=300d0
!        phimax=175d0
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
c      if(case.ne.'WWqqbr') then 
c         write(6,*)'Attempted to apply Higgs search cuts to process' 
c         write(6,*)'Check Runstring ' 
c         stop
c         return 
c      endif

      
      failed_cuts = .false. 
      do i=1,4
         et_vec(i)=0d0
      enddo

      if(pt(4,p).gt.pt(5,p)) then 
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
      if((pt_l_hard.lt.pth).or.(pt_l_soft.lt.pts)) then 
         failed_cuts=.true. 
         return 
      endif

!Binoth et al cuts 
!      if((pt_l_hard.gt.50d0).or.(pt_l_hard.lt.35d0)) then 
!         failed_cuts=.true. 
!         return
!      elseif((pt_l_soft.gt.25d0).or.(pt_l_soft.lt.20d0)) then 
!         failed_cuts=.true.
!         return
!      endif
     
      if(dabs(etmiss(p,et_vec)).lt.etmiss_min) then 
         failed_cuts=.true.
         return 
      endif

!---- eta_cuts
      if((dabs(eta_hard).gt.eta_max_h).or.
     &     ((dabs(eta_soft).gt.eta_max_s))) then 
         failed_cuts=.true. 
!         write(6,*) eta_hard,eta_soft
         return
      endif

     
!--- m_ll cut 
      m45=0d0
      do i=1,4
         if(i.ne.4) then 
            m45=m45-(p(4,i)+p(5,i))**2
         else
            m45=m45+(p(4,i)+p(5,i))**2 
         endif
      enddo

      if(dsqrt(max(m45,0d0)).gt.mllmax) then 
         failed_cuts=.true.
         return 
      endif

!---- phi_ll 
!--- convert phi max to radians
      phimax=phimax/(360d0)*2d0*pi

       r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     .     /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1D0
      if (r2 .lt. -0.9999999D0) r2=-1D0
      delphi=dacos(r2)
      if(delphi.gt.phimax) then 
         failed_cuts=.true. 
         return 
      endif

 99   format(1x,a29,f6.2,a17)

      return 
      end
      

      subroutine HWW_cuts_params(in_mh) 
      implicit none 
      double precision mh
      double precision pts,pth,mllmax,phimax
      double precision etmiss_min,eta_max_h,eta_max_s
      common/HWW_Cuts/pts,pth,mllmax,phimax,etmiss_min,
     &     eta_max_h,eta_max_s
      double precision tiny,in_mh
      tiny=1d-4 
!----- Ensure that an input of 130 etc is treated properly (i.e. always boosted into correct cuts
      mh=in_mh+tiny
      
!---- HWW_cuts from CMS 7 TeV
      if(mh.lt.130d0) then 
         pth=20d0
         pts=20d0
         mllmax=40d0
         phimax=60d0
      elseif((mh.gt.130d0).and.(mh.lt.150d0)) then 
         pth=25d0
         pts=20d0
         mllmax=45d0
         phimax=60d0
         
      elseif((mh.gt.150d0).and.(mh.lt.160d0)) then 
         pth=27d0
         pts=25d0
         mllmax=50d0
         phimax=60d0
         
         
      elseif((mh.gt.160d0).and.(mh.lt.170d0)) then 
         pth=30d0
         pts=25d0
         mllmax=50d0
         phimax=60d0

         
      elseif((mh.gt.170d0).and.(mh.lt.180d0)) then 
         pth=34d0
         pts=25d0
         mllmax=50d0
         phimax=60d0
         
      elseif((mh.gt.180d0).and.(mh.lt.190d0)) then 
         pth=36d0
         pts=25d0
         mllmax=60d0
         phimax=70d0
         
      elseif((mh.gt.190d0).and.(mh.lt.200d0)) then 
         pth=38d0
         pts=25d0
         mllmax=80d0
         phimax=90d0

      elseif((mh.gt.200d0).and.(mh.lt.210d0)) then 
         pth=40d0
         pts=25d0
         mllmax=90d0
         phimax=100d0
    
      elseif((mh.gt.210d0).and.(mh.lt.220d0)) then 
         pth=44d0
         pts=25d0
         mllmax=110d0
         phimax=110d0
         
      elseif((mh.gt.220d0).and.(mh.lt.230d0)) then 
         pth=48d0
         pts=25d0
         mllmax=120d0
         phimax=120d0

      elseif((mh.gt.230d0).and.(mh.lt.250d0)) then 
         pth=52d0
         pts=25d0
         mllmax=130d0
         phimax=130d0

      elseif((mh.gt.250d0).and.(mh.lt.300d0)) then 
         Pth=55d0
         pts=25d0
         mllmax=150d0
         phimax=140d0
         
      elseif((mh.gt.300d0)) then 
         pth=70d0
         pts=25d0
         mllmax=200d0
         phimax=175d0

      endif
      return 
         
      end
