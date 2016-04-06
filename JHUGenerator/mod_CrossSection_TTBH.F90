MODULE ModCrossSection_TTBH
IMPLICIT NONE

integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)


 CONTAINS

 
 

FUNCTION EvalWeighted_TTBH(yRnd,VgsWgt)
use ModKinematics
use ModParameters
use ModTTBHiggs
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
real(8) :: yRnd(1:16),VgsWgt, EvalWeighted_TTBH
real(8) :: pdf(-6:6,1:2)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi, PDFFac1,PDFFac2,xRnd
real(8) :: MomExt(1:4,1:13),MomOffShell(1:4,1:13),PSWgt,PSWgt2,PSWgt3,PSWgt4
real(8) :: LO_Res_Unpol,LO_Res_GG_Unpol,LO_Res_QQB_Unpol,PreFac  !, MG_MOM(0:3,1:5),MadGraph_tree
real(8) :: WdecayKfactor
integer :: NBin(1:NumHistograms),NHisto,iPartChannel,PartChannelAvg
integer :: MY_IDUP(1:13),ICOLUP(1:2,1:13),DK_IDUP(1:6),DK_ICOLUP(1:2,3:6)
logical :: applyPSCut
integer, parameter :: NumPartonicChannels=6
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13
EvalWeighted_TTBH = 0d0
WdecayKfactor = 1d0

   iPartChannel = int(yRnd(1) * NumPartonicChannels)
   if( PChannel.eq.0 ) then
      iPartChannel = 0
      PartChannelAvg = 1
   elseif( PChannel.eq.1 ) then
      iPartChannel = int(yRnd(1) * (NumPartonicChannels-1d0))+1
      PartChannelAvg = NumPartonicChannels - 1
   else
      PartChannelAvg = NumPartonicChannels
   endif
   if( iPartChannel.eq.0 ) then
      iPart_sel = 0
      jPart_sel = 0
      if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return
   elseif( iPartChannel.eq.1 ) then
      iPart_sel = Up_
      jPart_sel = AUp_
      if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return
   elseif( iPartChannel.eq.2 ) then
      iPart_sel = Dn_
      jPart_sel = ADn_
      if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return
   elseif( iPartChannel.eq.3 ) then
      iPart_sel = Chm_
      jPart_sel = AChm_
      if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return
   elseif( iPartChannel.eq.4 ) then
      iPart_sel = Str_
      jPart_sel = AStr_
      if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return
   elseif( iPartChannel.eq.5 ) then
      iPart_sel = Bot_
      jPart_sel = ABot_
      if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return
   endif

   
   
   call PDFMapping(2,yRnd(2:3),eta1,eta2,Ehat,sHatJacobi)  
   if (EHat.lt.2*M_Top+M_Reso) return
   call EvalPhasespace_2to3M(EHat,(/M_Reso,M_Top,M_Top/),yRnd(4:8),MomExt(1:4,1:5),PSWgt)! a(1)b(2)-->H(3)+tbar(4)+t(5)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
  
   call VVBranchings(DK_IDUP(1:6),DK_ICOLUP(1:2,3:6),700) ! Do not assign MY_IDUP yet
   if( TOPDECAYS.NE.0 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,tbar),yRnd(09:12),MomExt(1:4,06:08),PSWgt2)    ! ATop 
      MomExt(1:4,bbar) = MomExt(1:4,06)
      MomExt(1:4,nubar)= MomExt(1:4,08)
      MomExt(1:4,lepM) = MomExt(1:4,07)
      MomExt(1:4,Wm)   = MomExt(1:4,lepM) + MomExt(1:4,nubar)
      
      call EvalPhasespace_TopDecay(MomExt(1:4,t),yRnd(13:16),MomExt(1:4,10:12),PSWgt3)    !  Top
      MomExt(1:4,b)   = MomExt(1:4,10)
      MomExt(1:4,nu)  = MomExt(1:4,12)
      MomExt(1:4,lepP)= MomExt(1:4,11)
      MomExt(1:4,Wp)  = MomExt(1:4,lepP) + MomExt(1:4,nu)
      PSWgt = PSWgt * PSWgt2*PSWgt3

      call TTbar_OffShellProjection(MomExt,MomOffShell,PSWgt4)
      MomOffShell(1:4,1:3) = MomExt(1:4,1:3)      
!       PSWgt = PSWgt * PSWgt4       ! not using the Jacobian because the mat.el. don't have BW-propagators

      WdecayKfactor = (ScaleFactor( convertLHE(DK_IDUP(3)),convertLHE(DK_IDUP(4)) ))*(ScaleFactor( convertLHE(DK_IDUP(5)),convertLHE(DK_IDUP(6)) ))
   endif
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(17:18),MomExt(1:4,12:13),PSWgt5)
!    PSWgt = PSWgt * PSWgt5 


   call Kinematics_TTBH(MomOffShell,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * PartChannelAvg * WdecayKfactor
   call SetRunningScales( (/ MomExt(1:4,Hbos),MomExt(1:4,tbar),MomExt(1:4,t) /) , (/ Not_a_particle_,ATop_,Top_,Not_a_particle_ /) )
   call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)
   
   
   
   if( iPartChannel.eq.0 ) then
      call EvalAmp_GG_TTBH(MomExt,LO_Res_Unpol)    
      PDFFac1 = pdf(iPart_sel,1)*pdf(jPart_sel,2)
   else
      call EvalAmp_QQB_TTBH(MomExt,LO_Res_Unpol) 
      PDFFac1 = pdf(iPart_sel,1)*pdf(jPart_sel,2) + pdf(iPart_sel,2)*pdf(jPart_sel,1)
   endif
   EvalWeighted_TTBH = LO_Res_Unpol * PDFFac1 * PreFac
   
   
   
   
!    if( PChannel.eq.0 .or. PChannel.eq.2 ) then
!       call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
!       PDFFac1 = pdf(0,1)*pdf(0,2)
!       EvalWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1
!    endif
!    if( PChannel.eq.1 .or. PChannel.eq.2 ) then
!       call EvalAmp_QQB_TTBH(MomExt,LO_Res_QQB_Unpol)
!       PDFFac1 = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
!               + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
!               + pdf(Bot_,1)*pdf(ABot_,2)                            
!       PDFFac2 = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
!               + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
!               + pdf(Bot_,2)*pdf(ABot_,1)
!       EvalWeighted_TTBH = EvalWeighted_TTBH + LO_Res_QQB_Unpol * ( PDFFac1 + PDFFac2 )  
!    endif
!    EvalWeighted_TTBH = EvalWeighted_TTBH * PreFac
   
! print *, "checker",eta1,eta2,ehat,Mu_Fact,FluxFac,pdf(0,1)*pdf(0,2),( PDFFac1 + PDFFac2 )   
! print *, "gg",LO_Res_GG_Unpol         *FluxFac*pdf(0,1)*pdf(0,2)
! print *, "qq",LO_Res_QQB_Unpol        *FluxFac*( PDFFac1 + PDFFac2 )   
! pause
!       LO_Res_QQB_Unpol = LO_Res_QQB_Unpol/alphas**2*0.13d0**2
!       LO_Res_GG_Unpol = LO_Res_GG_Unpol/alphas**2*0.13d0**2
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBH(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alphas,m_top,m_z
!       print *, "My gg tree:         ", LO_Res_GG_Unpol/(100d0)**2
!       print *, "MadGraph gg hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_GG_Unpol)/(100d0)**2)
!       
!       call SUUB_TTBH(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alphas,m_top,m_z
!       print *, "My qqb tree:         ", LO_Res_QQB_Unpol/(100d0)**2
!       print *, "MadGraph qqb hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_QQB_Unpol)/(100d0)**2)
!       pause



if( unweighted ) then 

  if( warmup ) then
  
      CrossSec(iPart_sel,jPart_sel) = CrossSec(iPart_sel,jPart_sel) + EvalWeighted_TTBH*VgsWgt
      CrossSecMax(iPart_sel,jPart_sel) = max(CrossSecMax(iPart_sel,jPart_sel),EvalWeighted_TTBH*VgsWgt)
   
  else  
  
     ICOLUP(:,:) = 000
     ICOLUP(1:2,1) = (/501,510/)
     ICOLUP(1:2,2) = (/510,502/)   
     ICOLUP(1:2,tbar) = (/000,502/)
     ICOLUP(1:2,t)    = (/501,000/)
     MY_IDUP(b)    = Bot_;        ICOLUP(1:2,b) = (/501,00/)
     MY_IDUP(Wp)   = DK_IDUP(1);  ICOLUP(1:2,Wp)   = (/000,000/)
     MY_IDUP(lepP) = DK_IDUP(3);  ICOLUP(1:2,lepP) = DK_ICOLUP(1:2,3)
     MY_IDUP(nu)   = DK_IDUP(4);  ICOLUP(1:2,nu)   = DK_ICOLUP(1:2,4)  
     MY_IDUP(bbar) = ABot_;       ICOLUP(1:2,bbar) = (/000,502/)
     MY_IDUP(Wm)   = DK_IDUP(2);  ICOLUP(1:2,Wm)   = (/000,000/)             
     MY_IDUP(lepM) = DK_IDUP(6);  ICOLUP(1:2,lepM) = DK_ICOLUP(1:2,6)
     MY_IDUP(nubar)= DK_IDUP(5);  ICOLUP(1:2,nubar)= DK_ICOLUP(1:2,5)  
     
     call random_number(xRnd) 
     if( EvalWeighted_TTBH*VgsWgt.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",EvalWeighted_TTBH*VgsWgt, CrossSecMax(iPart_sel,jPart_sel)
!          write(io_stdout, "(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",EvalWeighted_TTBH*VgsWgt, CrossSecMax(iPart_sel,jPart_sel)
         AlertCounter = AlertCounter + 1
     elseif( EvalWeighted_TTBH*VgsWgt .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPart_sel,jPart_sel) = AccepCounter_part(iPart_sel,jPart_sel) + 1
         call WriteOutEvent_TTBH(MomOffShell,MY_IDUP(1:13),ICOLUP(1:2,1:13))
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
     endif
     
  endif! warmup
  
  
else


   if( writeWeightedLHE ) then 
        call Error("WriteLHE not yet supported for ttb+H")
   endif
   do NHisto=1,NumHistograms
       call intoHisto(NHisto,NBin(NHisto),EvalWeighted_TTBH*VgsWgt)
   enddo

   
endif! unweighted
   
   
RETURN
END FUNCTION 





FUNCTION EvalUnWeighted_TTBH(yRnd,genEvt,iPartons,RES)
use ModKinematics
use ModParameters
use ModTTBHiggs
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
real(8) :: yRnd(1:16),VgsWgt, EvalUnWeighted_TTBH
real(8) :: pdf(-6:6,1:2),RES(-5:5,-5:5)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi,CS_Max,DKRnd
real(8) :: MomExt(1:4,1:13),MomOffShell(1:4,1:13),PSWgt,PSWgt2,PSWgt3,PSWgt4
real(8) :: LO_Res_GG_Unpol,LO_Res_QQB_Unpol,PreFac,PDFFac1,PDFFac2
real(8) :: WdecayKfactor
integer :: NBin(1:NumHistograms),NHisto,iPartons(1:2),DKFlavor
integer :: MY_IDUP(1:13),ICOLUP(1:2,1:13),nparton,DK_IDUP(1:6),DK_ICOLUP(1:2,3:6)
logical :: applyPSCut,genEvt
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13
include 'csmaxvalue.f'  
EvalUnWeighted_TTBH = 0d0
WdecayKfactor = 1d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   
   if (EHat.lt.2*M_Top+M_Reso) return
   call EvalPhasespace_2to3M(EHat,(/M_Reso,M_Top,M_Top/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)! a(1)b(2)-->H(3)+tbar(4)+t(5)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   ICOLUP(1:2,Hbos) = (/000,000/)
   ICOLUP(1:2,tbar) = (/000,502/)
   ICOLUP(1:2,t)    = (/501,000/)
   if( TOPDECAYS.NE.0 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,tbar),yRnd(08:11),MomExt(1:4,06:08),PSWgt2)    ! ATop 
      MomExt(1:4,bbar) = MomExt(1:4,06)
      MomExt(1:4,nubar)= MomExt(1:4,08)
      MomExt(1:4,lepM) = MomExt(1:4,07)
      MomExt(1:4,Wm)   = MomExt(1:4,lepM) + MomExt(1:4,nubar)
      
      call EvalPhasespace_TopDecay(MomExt(1:4,t),yRnd(12:15),MomExt(1:4,10:12),PSWgt3)    !  Top
      MomExt(1:4,b)   = MomExt(1:4,10)
      MomExt(1:4,nu)  = MomExt(1:4,12)
      MomExt(1:4,lepP)= MomExt(1:4,11)
      MomExt(1:4,Wp)  = MomExt(1:4,lepP) + MomExt(1:4,nu)
      PSWgt = PSWgt * PSWgt2*PSWgt3
      
      call TTbar_OffShellProjection(MomExt,MomOffShell,PSWgt4)
      MomOffShell(1:4,1:3) = MomExt(1:4,1:3)  
!       PSWgt = PSWgt * PSWgt4       ! not using the Jacobian because the mat.el. don't have BW-propagators

      if( RandomizeVVdecays ) then 
         call random_number(DKRnd)
         if( DKRnd.lt.0.5d0 ) call swapi(DecayMode1,DecayMode2)
      endif
      call VVBranchings(DK_IDUP(1:6),DK_ICOLUP(1:2,3:6),700)
      MY_IDUP(b)    = Bot_;        ICOLUP(1:2,b) = (/501,00/)
      MY_IDUP(Wp)   = DK_IDUP(1);  ICOLUP(1:2,Wp)   = (/000,000/)
      MY_IDUP(lepP) = DK_IDUP(3);  ICOLUP(1:2,lepP) = DK_ICOLUP(1:2,3)
      MY_IDUP(nu)   = DK_IDUP(4);  ICOLUP(1:2,nu)   = DK_ICOLUP(1:2,4)  
      MY_IDUP(bbar) = ABot_;       ICOLUP(1:2,bbar) = (/000,502/)
      MY_IDUP(Wm)   = DK_IDUP(2);  ICOLUP(1:2,Wm)   = (/000,000/)             
      MY_IDUP(lepM) = DK_IDUP(6);  ICOLUP(1:2,lepM) = DK_ICOLUP(1:2,6)
      MY_IDUP(nubar)= DK_IDUP(5);  ICOLUP(1:2,nubar)= DK_ICOLUP(1:2,5)

      WdecayKfactor = (ScaleFactor( convertLHE(DK_IDUP(3)),convertLHE(DK_IDUP(4)) ))*(ScaleFactor( convertLHE(DK_IDUP(5)),convertLHE(DK_IDUP(6)) ))
   else
      MY_IDUP(6:11)=-9999
   endif
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
!    PSWgt = PSWgt * PSWgt4 
   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * WdecayKfactor

   call Kinematics_TTBH(MomOffShell,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return

   call SetRunningScales( (/ MomExt(1:4,Hbos),MomExt(1:4,tbar),MomExt(1:4,t) /) , (/ Not_a_particle_,ATop_,Top_,Not_a_particle_ /) )
   call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)

   
IF( GENEVT ) THEN   
      
      if( iPartons(1).eq.0 .and. iPartons(2).eq.0 ) then
          call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
          PDFFac1 = pdf(0,1)*pdf(0,2)
          EvalUnWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1 * PreFac
          MY_IDUP(1:5) = (/Glu_,Glu_,Hig_,ATop_,Top_/)
          ICOLUP(1:2,1) = (/501,510/)
          ICOLUP(1:2,2) = (/510,502/)    
      else
          call EvalAmp_QQB_TTBH(MomExt,LO_Res_QQB_Unpol)
          PDFFac1 = pdf( LHA2M_pdf(iPartons(1)),1) * pdf( LHA2M_pdf(iPartons(2)),2)
          EvalUnWeighted_TTBH = LO_Res_QQB_Unpol * PDFFac1 * PreFac 
          MY_IDUP(1:5) = (/ LHA2M_pdf(iPartons(1)),LHA2M_pdf(iPartons(2)),Hig_,ATop_,Top_/)
          if( iPartons(1).gt.0 ) then
             ICOLUP(1:2,1) = (/501,000/)
             ICOLUP(1:2,2) = (/000,502/)
          else
             ICOLUP(1:2,1) = (/000,502/)
             ICOLUP(1:2,2) = (/501,000/)
          endif
      endif
      
      CS_max = CSmax(iPartons(1),iPartons(2))
      if( EvalUnWeighted_TTBH .gt. CS_max) then
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_TTBH, CS_max
         AlertCounter = AlertCounter + 1
      elseif( EvalUnWeighted_TTBH .gt. yRnd(16)*CS_max ) then
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPartons(1),iPartons(2)) = AccepCounter_part(iPartons(1),iPartons(2))+1
         call WriteOutEvent_TTBH(MomOffShell,MY_IDUP(1:13),ICOLUP(1:2,1:13))
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
      endif
      EvalCounter = EvalCounter + 1 
      


ELSE! NOT GENEVT

      
      if( PChannel.eq.0 .or. PChannel.eq.2 ) then
          call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
          PDFFac1 = pdf(0,1)*pdf(0,2)
          EvalUnWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1 * PreFac
          
          RES(0,0) = EvalUnWeighted_TTBH
          if (EvalUnWeighted_TTBH.gt.csmax(0,0)) then
              CSmax(0,0) = EvalUnWeighted_TTBH
          endif
      endif
      
      if( PChannel.eq.1 .or. PChannel.eq.2 ) then
          call EvalAmp_QQB_TTBH(MomExt,LO_Res_QQB_Unpol)       
          do nparton = -5,5
            if (nparton.eq.-5) then
              PDFFac1 = pdf(Bot_,2)*pdf(ABot_,1)
            elseif(nparton.eq.-4) then
              PDFFac1 = pdf(Chm_,2)*pdf(AChm_,1)
            elseif(nparton.eq.-3) then
              PDFFac1 = pdf(Str_,2)*pdf(AStr_,1)
            elseif(nparton.eq.-2) then
              PDFFac1 = pdf(Up_,2) *pdf(AUp_,1)
            elseif(nparton.eq.-1) then
              PDFFac1 = pdf(Dn_,2) *pdf(ADn_,1)
            elseif (nparton.eq.0) then
              cycle
            elseif (nparton.eq.1) then
              PDFFac1 = pdf(Dn_,1) *pdf(ADn_,2)
            elseif (nparton.eq.2) then
              PDFFac1 = pdf(Up_,1) *pdf(AUp_,2)
            elseif(nparton.eq.3) then
              PDFFac1 = pdf(Str_,1)*pdf(AStr_,2)
            elseif(nparton.eq.4) then
              PDFFac1 = pdf(Chm_,1)*pdf(AChm_,2)
            elseif(nparton.eq.5) then
              PDFFac1 = pdf(Bot_,1)*pdf(ABot_,2)
            endif
            EvalUnWeighted_TTBH = LO_Res_QQB_Unpol * PreFac *PDFFac1
            RES(nparton,-nparton) = EvalUnWeighted_TTBH
            if (EvalUnWeighted_TTBH.gt.csmax(nparton,-nparton)) CSmax(nparton,-nparton) = EvalUnWeighted_TTBH
          enddo
      endif


ENDIF! GENEVT 


RETURN
END FUNCTION






END MODULE ModCrossSection_TTBH







