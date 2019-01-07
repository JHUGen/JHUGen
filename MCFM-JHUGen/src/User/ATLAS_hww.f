!---- Subroutine which calculates the ATLAS Higgs search cuts of 1106.2748
!---- C. Williams July 11
      subroutine ATLAS_hww(p,failed_cuts)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'first.f'
      double precision p(mxpart,4)
      logical failed_cuts
      double precision pt,etarap,m45
      double precision pt_l_hard,pt_l_soft,pttwo
      double precision pts,pth,phimax,mllmax,etmiss,etmiss_min
      double precision et_vec(4),r2,delphi,eta_max_h
      double precision eta_max_s,eta_hard,eta_soft
      double precision mtmin,mtmax,mt45,ptsq,p3456(2)
      integer i
      common/HWW_Cuts/pts,pth,mllmax,phimax,etmiss_min,
     &     eta_max_h,eta_max_s
      double precision tiny
      tiny=1d-4

!------- ETA DOESNT CHANGE -------------
      eta_max_h=2.5d0
      eta_max_s=2.5d0
!----------------------------------------

!------- PT of leptons doesn't change either

      pts=15d0
      pth=20d0

!------- Neither does ETMiss
      etmiss_min=30d0

!----- m_ll_max delphi depend on m_h
      hmass=hmass-tiny
      if(hmass.lt.170d0) then
         mllmax=50d0
         phimax=1.3d0 !-- in radians
      else
         mllmax=60d0
         phimax=1.8d0
      endif
!----- Reset hmass
      hmass=hmass+tiny

      mtmin=0.75d0*hmass
      mtmax=hmass





      if(first) then
         first=.false.
      write(6,*)  '********** ATLAS Higgs Search cuts  ****************'
      write(6,*)  '*                                                  *'
      write(6,99)  '*     pt(lep)_max>  ', pth, '                    *'
      write(6,99)  '*     pt(lep)_min >  ', pts, '                    *'
      write(6,99)  '*     mll  <  ', mllmax, '                     *'
      write(6,99)  '*     phi(lep,lep) <  ', phimax, '               *'
      write(6,99)  '*     |eta_l | <  ', eta_max_h, '                 *'
      write(6,99) '*       ET_miss > ',etmiss_min, '                *'
      write(6,99) '*       mT  > ',mtmin,'                           * '
      write(6,99) '*       mT  < ',mtmax,'                           *'
      write(6,*)  '****************************************************'
      endif

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
!----- Phi cut
      r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     .     /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1D0
      if (r2 .lt. -0.9999999D0) r2=-1D0
      delphi=dacos(r2)
      if(delphi.gt.phimax) then
         failed_cuts=.true.
         return
      endif


!------mt cut (use m45=mt)
      ptsq=0d0
      do i=1,2
         p3456(i)=0d0
         p3456(i)=p(3,i)+p(4,i)+p(5,i)+p(6,i)
         ptsq=ptsq+p3456(i)**2
      enddo


      mt45=0d0
      mt45=(dsqrt(dsqrt(pttwo(4,5,p)**2+m45)+etmiss(p,et_vec))**2)
!      write(6,*) mt45
 !     mt45=dsqrt(max(0d0,dsqrt(mt45-ptsq)))
 !     write(6,*) mt45
 !     pause
!------ Mt cuts
      if((mt45.lt.mtmin).or.(mt45.gt.mtmax)) then
         failed_cuts=.true.
         return
      endif


 99   format(1x,a29,f6.2,a17)

      return
      end subroutine


