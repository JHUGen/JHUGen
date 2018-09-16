!---- Subroutine which calculates the ATLAS Higgs search cuts of 1106.2748
!---- C. Williams July 11
      subroutine ATLAS_hww(p,failed_cuts)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'first.f'
      real(dp):: p(mxpart,4)
      logical:: failed_cuts
      real(dp):: pt,etarap,m45
      real(dp):: pt_l_hard,pt_l_soft,pttwo
      real(dp):: pts,pth,phimax,mllmax,etmiss,etmiss_min
      real(dp):: et_vec(4),r2,delphi,eta_max_h
      real(dp):: eta_max_s,eta_hard,eta_soft
      real(dp):: mtmin,mtmax,mt45,ptsq,p3456(2)
      integer:: i
      common/HWW_Cuts/pts,pth,mllmax,phimax,etmiss_min,
     &     eta_max_h,eta_max_s
      real(dp):: tiny
      tiny=1.e-4_dp

!------- ETA DOESNT CHANGE -------------
      eta_max_h=2.5_dp
      eta_max_s=2.5_dp
!----------------------------------------

!------- PT of leptons doesn't change either

      pts=15._dp
      pth=20._dp

!------- Neither does ETMiss
      etmiss_min=30._dp

!----- m_ll_max delphi depend on m_h
      hmass=hmass-tiny
      if(hmass<170._dp) then
         mllmax=50._dp
         phimax=1.3_dp !-- in radians
      else
         mllmax=60._dp
         phimax=1.8_dp
      endif
!----- Reset hmass
      hmass=hmass+tiny

      mtmin=0.75_dp*hmass
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
!----- Phi cut
      r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     &     /sqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (r2 > +0.9999999_dp) r2=+one
      if (r2 < -0.9999999_dp) r2=-one
      delphi=acos(r2)
      if(delphi>phimax) then
         failed_cuts=.true.
         return
      endif


!------mt cut (use m45=mt)
      ptsq=zip
      do i=1,2
         p3456(i)=zip
         p3456(i)=p(3,i)+p(4,i)+p(5,i)+p(6,i)
         ptsq=ptsq+p3456(i)**2
      enddo


      mt45=zip
      mt45=(sqrt(sqrt(pttwo(4,5,p)**2+m45)+etmiss(p,et_vec))**2)
!      write(6,*) mt45
 !     mt45=sqrt(max(zip,sqrt(mt45-ptsq)))
 !     write(6,*) mt45
 !     pause
!------ Mt cuts
      if((mt45<mtmin).or.(mt45>mtmax)) then
         failed_cuts=.true.
         return
      endif


 99   format(1x,a29,f6.2,a17)

      return
      end subroutine


