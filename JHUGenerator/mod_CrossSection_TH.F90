MODULE ModCrossSection_TH
! Author: R. Rontsch, June 2015

 integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)




 contains



FUNCTION EvalWeighted_TH(yRnd,VgsWgt)
! Routine for production of H(p3)+t(p4)+jet(p5)
! Top decays taken from implementation in MCFM, see hep-ph:/1204.1513
use ModKinematics
use ModParameters
use ModMisc
use ModTHiggs
implicit none
real(8) :: EvalWeighted_TH,yRnd(1:11),VgsWgt
real(8) :: Ehat,MH_Inv,eta1,eta2,ISFac,sHatJacobi,PreFac,FluxFac,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:9),MomOffShell(1:4,1:9),LO_Res_Unpol(-6:6,-6:6),MuFac,pdf(-6:6,1:2)
integer :: NBin(1:NumHistograms),NHisto
logical :: applyPSCut
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9


    
   EvalWeighted_TH = 0d0   
   EvalCounter = EvalCounter + 1

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   MH_Inv = M_Reso
   if( EHat.le.m_Top+MH_Inv ) then
      EvalWeighted_TH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3ArbMass(EHat,(/MH_Inv,M_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   MuFac=(M_Top + M_Reso)/4d0
   call setPDFs(eta1,eta2,MuFac,pdf)
   
   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),MomExt(1:4,6:8),PSWgt2)
!       MomExt(1:4,b)  = MomExt(1:4,6)
      MomExt(1:4,nu) = MomExt(1:4,8)
      MomExt(1:4,lep)= MomExt(1:4,7)
      MomExt(1:4,W)  = MomExt(1:4,lep) + MomExt(1:4,nu)
      PSWgt = PSWgt * PSWgt2

      call Top_OffShellProjection(MomExt,MomOffShell,PSWgt3)  
      MomOffShell(1:4,1:3) = MomExt(1:4,1:3)            
!       PSWgt = PSWgt * PSWgt3       ! not using the Jacobian because the mat.el. don't have BW-propagators
      call Kinematics_TH(MomOffShell,applyPSCut,NBin)
   ELSE
      call Kinematics_TH(MomExt,applyPSCut,NBin)
   ENDIF
   
   
   if( applyPSCut ) then
      EvalWeighted_TH = 0d0
      return
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt
   
   IF( PROCESS.EQ.110 ) THEN        ! t-channel production of t H
      call EvalAmp_QB_TH(MomExt,LO_Res_Unpol)
      EvalWeighted_TH = &
                      + LO_Res_Unpol(Up_,Bot_)   * ( pdf(Up_,1) *pdf(Bot_,2)  +  pdf(Chm_,1) *pdf(Bot_,2) )   &
                      + LO_Res_Unpol(Bot_,Up_)   * ( pdf(Bot_,1)*pdf(Up_,2)   +  pdf(Bot_,1) *pdf(Chm_,2) )   &
                      + LO_Res_Unpol(ADn_,Bot_)  * ( pdf(ADn_,1)*pdf(Bot_,2)  +  pdf(AStr_,1)*pdf(Bot_,2) )   &
                      + LO_Res_Unpol(Bot_,ADn_)  * ( pdf(Bot_,1)*pdf(ADn_,2)  +  pdf(Bot_,1) *pdf(AStr_,2))
   ELSEIF( PROCESS.EQ.111 ) THEN      ! t-channel production of tbar H
      call EvalAmp_QbarBbar_TH(MomExt,LO_Res_Unpol)
      EvalWeighted_TH = &
                      + LO_Res_Unpol(Dn_,ABot_)  * ( pdf(Dn_,1)*pdf(ABot_,2)  + pdf(Str_,1)*pdf(ABot_,2) )    &
                      + LO_Res_Unpol(ABot_,Dn_)  * ( pdf(ABot_,1)*pdf(Dn_,2)  + pdf(ABot_,1)*pdf(Str_,2) )    &
                      + LO_Res_Unpol(AUp_,ABot_) * ( pdf(AUp_,1)*pdf(ABot_,2) + pdf(AChm_,1)*pdf(ABot_,2))    &
                      + LO_Res_Unpol(ABot_,AUp_) * ( pdf(ABot_,1)*pdf(AUp_,2) + pdf(ABot_,1)*pdf(AChm_,2))
   ELSEIF (PROCESS .EQ. 112) THEN      ! s-channel production of t H
      call EvalAmp_QQB_THBBAR(MomExt,LO_Res_Unpol)
      EvalWeighted_TH = &
           + LO_Res_Unpol(Up_,ADn_)  * ( pdf(Up_,1)*pdf(ADn_,2)  + pdf(Chm_,1)*pdf(AStr_,2) )    &
           + LO_Res_Unpol(ADn_,Up_)  * ( pdf(ADn_,1)*pdf(Up_,2)  + pdf(AStr_,1)*pdf(Chm_,2) )    
   ELSEIF (PROCESS .EQ. 113) THEN      ! s-channel production of tbar H
      call EvalAmp_QQB_TBARHB(MomExt,LO_Res_Unpol)
      EvalWeighted_TH = &
           + LO_Res_Unpol(Dn_,AUp_)  * ( pdf(Dn_,1)*pdf(AUp_,2)  + pdf(Str_,1)*pdf(AChm_,2) )    &
           + LO_Res_Unpol(AUp_,Dn_)  * ( pdf(AUp_,1)*pdf(Dn_,2)  + pdf(AChm_,1)*pdf(Str_,2) )                 
   ENDIF
   EvalWeighted_TH = EvalWeighted_TH * PreFac


   if( writeWeightedLHE ) then 
        call Error("WriteLHE not yet supported for t+H")
   endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalWeighted_TH*VgsWgt)                
   enddo
   AccepCounter=AccepCounter+1


 RETURN
 END FUNCTION
 












FUNCTION EvalUnWeighted_TH(yRnd,genEvt,iPartons,RES)
use ModKinematics
use ModParameters
use ModTHiggs
use ModMisc
implicit none
real(8) :: yRnd(1:16),VgsWgt, EvalUnWeighted_TH
real(8) :: pdf(-6:6,1:2),RES(-5:5,-5:5)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi,MH_Inv,MuFac
real(8) :: MomExt(1:4,1:9),MomOffShell(1:4,1:9),PSWgt,PSWgt2,PSWgt3
real(8) :: LO_Res_Unpol(-6:6,-6:6),PreFac,PDFFac1,CS_Max,DKRnd
integer :: NBin(1:NumHistograms),NHisto,iPartons(1:2),DKFlavor
integer :: MY_IDUP(1:9),ICOLUP(1:2,1:9),nparton,DK_IDUP(1:6),DK_ICOLUP(1:2,3:6)
logical :: applyPSCut,genEvt
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'csmaxvalue.f'  
EvalUnWeighted_TH = 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   MH_Inv = M_Reso
   if( EHat.le.m_Top+MH_Inv ) return

   call EvalPhaseSpace_2to3ArbMass(EHat,(/MH_Inv,M_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)!  inLeft, inRight, Higgs, Top, Quark
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),MomExt(1:4,6:8),PSWgt2)
      MomExt(1:4,b)  = MomExt(1:4,6)
      MomExt(1:4,nu) = MomExt(1:4,8)
      MomExt(1:4,lep)= MomExt(1:4,7)
      MomExt(1:4,W)  = MomExt(1:4,lep) + MomExt(1:4,nu)
      PSWgt = PSWgt * PSWgt2
      
      call Top_OffShellProjection(MomExt,MomOffShell,PSWgt3)  
      MomOffShell(1:4,1:3) = MomExt(1:4,1:3)            
!       PSWgt = PSWgt * PSWgt3        ! not using the Jacobian because the mat.el. don't have BW-propagators
      
      call VVBranchings(DK_IDUP(1:6),DK_ICOLUP(1:2,3:6))
      if( PROCESS.EQ.110 ) then
          ICOLUP(1:2,Hbos) = (/000,000/)
          ICOLUP(1:2,t)    = (/501,000/)
          if( iPartons(inLeft)*iPartons(inRight).gt.0 ) then
              ICOLUP(1:2,qout) = (/502,000/)
          else
              ICOLUP(1:2,qout) = (/000,502/)
          endif
          if( iPartons(inLeft).gt.0 ) then
              ICOLUP(1:2,inLeft) = (/502,000/)
          else
              ICOLUP(1:2,inLeft) = (/000,502/)
          endif
          if( iPartons(inRight).gt.0 ) then
              ICOLUP(1:2,inRight)= (/501,000/)
          else
              ICOLUP(1:2,inRight)= (/000,501/)
          endif
          MY_IDUP(b)   = Bot_;       ICOLUP(1:2,b)   = (/501,000/)
          MY_IDUP(W)   = DK_IDUP(1); ICOLUP(1:2,W)   = (/000,000/)
          MY_IDUP(lep) = DK_IDUP(3); ICOLUP(1:2,lep) = DK_ICOLUP(1:2,3)
          MY_IDUP(nu)  = DK_IDUP(4); ICOLUP(1:2,nu)  = DK_ICOLUP(1:2,4)  
      elseif( PROCESS.EQ.111 ) then
          ICOLUP(1:2,Hbos) = (/000,000/)
          ICOLUP(1:2,t)    = (/000,501/)
          if( iPartons(inLeft)*iPartons(inRight).gt.0 ) then
              ICOLUP(1:2,qout) = (/000,502/)
          else
              ICOLUP(1:2,qout) = (/502,000/)
          endif
          if( iPartons(inLeft).gt.0 ) then
              ICOLUP(1:2,inLeft) = (/502,000/)
          else
              ICOLUP(1:2,inLeft) = (/000,502/)
          endif
          if( iPartons(inRight).gt.0 ) then
              ICOLUP(1:2,inRight)= (/501,000/)
          else
              ICOLUP(1:2,inRight)= (/000,501/)
          endif
          MY_IDUP(b)   = ABot_;      ICOLUP(1:2,b)  = (/000,501/)
          MY_IDUP(W)   = DK_IDUP(2); ICOLUP(1:2,W)  = (/000,000/)             
          MY_IDUP(lep) = DK_IDUP(6); ICOLUP(1:2,lep)= DK_ICOLUP(1:2,6)
          MY_IDUP(nu)  = DK_IDUP(5); ICOLUP(1:2,nu) = DK_ICOLUP(1:2,5)  
      elseif( PROCESS.EQ.112 ) then
          ICOLUP(1:2,Hbos) = (/000,000/)
          ICOLUP(1:2,t)    = (/502,000/)
          ICOLUP(1:2,qout) = (/000,502/)
          if( iPartons(inLeft).gt.0 ) then
              ICOLUP(1:2,inLeft) = (/501,000/)
          else
              ICOLUP(1:2,inLeft) = (/000,501/)
          endif
          if( iPartons(inRight).gt.0 ) then
              ICOLUP(1:2,inRight)= (/501,000/)
          else
              ICOLUP(1:2,inRight)= (/000,501/)
          endif
          MY_IDUP(b)   = Bot_;       ICOLUP(1:2,b)   = (/502,00/)
          MY_IDUP(W)   = DK_IDUP(1); ICOLUP(1:2,W)   = (/000,000/)
          MY_IDUP(lep) = DK_IDUP(3); ICOLUP(1:2,lep) = DK_ICOLUP(1:2,3)
          MY_IDUP(nu)  = DK_IDUP(4); ICOLUP(1:2,nu)  = DK_ICOLUP(1:2,4)  
      elseif( PROCESS.EQ.113 ) then   
          ICOLUP(1:2,Hbos) = (/000,000/)
          ICOLUP(1:2,t)    = (/000,502/)
          ICOLUP(1:2,qout) = (/502,000/)
          if( iPartons(inLeft).gt.0 ) then
              ICOLUP(1:2,inLeft) = (/501,000/)
          else
              ICOLUP(1:2,inLeft) = (/000,501/)
          endif
          if( iPartons(inRight).gt.0 ) then
              ICOLUP(1:2,inRight)= (/501,000/)
          else
              ICOLUP(1:2,inRight)= (/000,501/)
          endif
          MY_IDUP(b)   = ABot_;      ICOLUP(1:2,b)  = (/000,502/)
          MY_IDUP(W)   = DK_IDUP(2); ICOLUP(1:2,W)  = (/000,000/)             
          MY_IDUP(lep) = DK_IDUP(6); ICOLUP(1:2,lep)= DK_ICOLUP(1:2,6)
          MY_IDUP(nu)  = DK_IDUP(5); ICOLUP(1:2,nu) = DK_ICOLUP(1:2,5)  
      endif

   ELSE

      MomOffShell(1:4,1:5) = MomExt(1:4,1:5)
!       if( PROCESS.EQ.110 ) then
!           ICOLUP(1:2,Hbos) = (/000,000/)
!           ICOLUP(1:2,t)    = (/501,000/)
!           ICOLUP(1:2,qout) = (/502,000/)
!           ICOLUP(1:2,inLeft) = ICOLUP(1:2,qout)
!           ICOLUP(1:2,inRight)= ICOLUP(1:2,t) 
!       elseif( PROCESS.EQ.111 ) then
!           ICOLUP(1:2,Hbos) = (/000,000/)
!           ICOLUP(1:2,t)    = (/000,501/)
!           ICOLUP(1:2,qout) = (/000,502/)
!           ICOLUP(1:2,inLeft) = ICOLUP(1:2,qout)
!           ICOLUP(1:2,inRight)= ICOLUP(1:2,t) 
! ! Markus, please check this: 
!       elseif( PROCESS.EQ.112 ) then
!           ICOLUP(1:2,Hbos) = (/000,000/)
!           ICOLUP(1:2,t)    = (/501,000/)
!           ICOLUP(1:2,qout) = (/000,502/)
!           ICOLUP(1:2,inLeft) = ICOLUP(1:2,qout)
!           ICOLUP(1:2,inRight)= ICOLUP(1:2,t) 
!       elseif( PROCESS.EQ.113 ) then
!           ICOLUP(1:2,Hbos) = (/000,000/)
!           ICOLUP(1:2,t)    = (/000,501/)
!           ICOLUP(1:2,qout) = (/502,000/)
!           ICOLUP(1:2,inLeft) = ICOLUP(1:2,qout)
!           ICOLUP(1:2,inRight)= ICOLUP(1:2,t) 
!       endif
!       MY_IDUP(6:9)=-9999
   ENDIF

   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt

   call Kinematics_TH(MomOffShell,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return

   MuFac=(M_Top + M_Reso)/4d0
   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   LO_Res_Unpol = 0d0


IF( GENEVT ) THEN   

          IF( PROCESS.EQ.110 ) THEN
              call EvalAmp_QB_TH(MomExt,LO_Res_Unpol)
              MY_IDUP(1:5) = (/ LHA2M_pdf(iPartons(1)),LHA2M_pdf(iPartons(2)),Hig_,SU2flip(LHA2M_pdf(iPartons(2))),SU2flip(LHA2M_pdf(iPartons(1))) /)
              if( abs(iPartons(1)).eq.5 ) then
                  call swapi( MY_IDUP(4),MY_IDUP(5) )
                  if( iPartons(1)*iPartons(2).gt.0 ) then
                      call swapi( ICOLUP(1,inLeft),ICOLUP(1,inRight) )
                      call swapi( ICOLUP(2,inLeft),ICOLUP(2,inRight) )
                  else
                      call swapi( ICOLUP(1,inLeft),ICOLUP(2,inRight) )
                      call swapi( ICOLUP(2,inLeft),ICOLUP(1,inRight) )
                  endif
              endif              
          ELSEIF( PROCESS.EQ.111 ) THEN
              call EvalAmp_QbarBbar_TH(MomExt,LO_Res_Unpol)
              MY_IDUP(1:5) = (/ LHA2M_pdf(iPartons(1)),LHA2M_pdf(iPartons(2)),Hig_,SU2flip(LHA2M_pdf(iPartons(2))),SU2flip(LHA2M_pdf(iPartons(1))) /)
              if( abs(iPartons(1)).eq.5 ) then
                  call swapi( MY_IDUP(4),MY_IDUP(5) )
                  if( iPartons(1)*iPartons(2).gt.0 ) then
                      call swapi( ICOLUP(1,inLeft),ICOLUP(1,inRight) )
                      call swapi( ICOLUP(2,inLeft),ICOLUP(2,inRight) )
                  else
                      call swapi( ICOLUP(1,inLeft),ICOLUP(2,inRight) )
                      call swapi( ICOLUP(2,inLeft),ICOLUP(1,inRight) )
                  endif
              endif              
          ELSEIF (PROCESS .EQ. 112) THEN      
              MY_IDUP(1:5) = (/ LHA2M_pdf(iPartons(1)),LHA2M_pdf(iPartons(2)),Hig_,Top_,ABot_ /)
              call EvalAmp_QQB_THBBAR(MomExt,LO_Res_Unpol)
          ELSEIF (PROCESS .EQ. 113) THEN      
              MY_IDUP(1:5) = (/ LHA2M_pdf(iPartons(1)),LHA2M_pdf(iPartons(2)),Hig_,ATop_,Bot_ /)
              call EvalAmp_QQB_TBARHB(MomExt,LO_Res_Unpol)
          ENDIF


          
          PDFFac1 = pdf( LHA2M_pdf(iPartons(1)),1) * pdf( LHA2M_pdf(iPartons(2)),2)
          EvalUnWeighted_TH = LO_Res_Unpol(LHA2M_pdf(iPartons(1)),LHA2M_pdf(iPartons(2))) * PDFFac1 * PreFac 
      
          CS_max = CSmax(iPartons(1),iPartons(2))
          if( EvalUnWeighted_TH .gt. CS_max) then
            write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_TH, CS_max
            AlertCounter = AlertCounter + 1
          elseif( EvalUnWeighted_TH .gt. yRnd(16)*CS_max ) then
            AccepCounter = AccepCounter + 1
            AccepCounter_part(iPartons(1),iPartons(2)) = AccepCounter_part(iPartons(1),iPartons(2))+1
            call WriteOutEvent_TH(MomOffShell,MY_IDUP(1:9),ICOLUP(1:2,1:9))
            do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),1d0)
            enddo
          endif
          EvalCounter = EvalCounter + 1 
      

ELSE! NOT GENEVT


          IF( PROCESS.EQ.110 ) THEN
              call EvalAmp_QB_TH(MomExt,LO_Res_Unpol)
              do nparton = -5,5   ! LHE conventions
                PDFFac1 = pdf( LHA2M_pdf(nparton),1) * pdf(Bot_,2)
                EvalUnWeighted_TH = LO_Res_Unpol(LHA2M_pdf(nparton),Bot_) * PreFac *PDFFac1
                RES(nparton,+5) = EvalUnWeighted_TH
                if (EvalUnWeighted_TH.gt.csmax(nparton,+5)) CSmax(nparton,+5) = EvalUnWeighted_TH

                PDFFac1 = pdf( LHA2M_pdf(nparton),2) * pdf(Bot_,1)
                EvalUnWeighted_TH = LO_Res_Unpol(Bot_,LHA2M_pdf(nparton)) * PreFac *PDFFac1
                RES(+5,nparton) = EvalUnWeighted_TH
                if (EvalUnWeighted_TH.gt.csmax(+5,nparton)) CSmax(+5,nparton) = EvalUnWeighted_TH
              enddo

          ELSEIF( PROCESS.EQ.111 ) THEN
              call EvalAmp_QbarBbar_TH(MomExt,LO_Res_Unpol)
              do nparton = -5,5
                PDFFac1 = pdf( LHA2M_pdf(nparton),1) * pdf(ABot_,2)
                EvalUnWeighted_TH = LO_Res_Unpol(LHA2M_pdf(nparton),ABot_) * PreFac *PDFFac1
                RES(nparton,-5) = EvalUnWeighted_TH
                if (EvalUnWeighted_TH.gt.csmax(nparton,-5)) CSmax(nparton,-5) = EvalUnWeighted_TH

                PDFFac1 = pdf( LHA2M_pdf(nparton),2) * pdf(ABot_,1)
                EvalUnWeighted_TH = LO_Res_Unpol(ABot_,LHA2M_pdf(nparton)) * PreFac *PDFFac1
                RES(-5,nparton) = EvalUnWeighted_TH
                if (EvalUnWeighted_TH.gt.csmax(-5,nparton)) CSmax(-5,nparton) = EvalUnWeighted_TH
             enddo
          ELSEIF (PROCESS .EQ. 112) THEN      ! s-channel production of t H
               call EvalAmp_QQB_THBBAR(MomExt,LO_Res_Unpol)
               RES(2,-1) = LO_Res_Unpol(Up_,ADn_)   * pdf(Up_,1)*pdf(ADn_,2)   * PreFac
               RES(4,-3) = LO_Res_Unpol(Chm_,AStr_) * pdf(Chm_,1)*pdf(AStr_,2) * PreFac
               RES(-1,2) = LO_Res_Unpol(ADn_,Up_)   * pdf(Up_,2)*pdf(ADn_,1)   * PreFac
               RES(-3,4) = LO_Res_Unpol(AStr_,Chm_) * pdf(Chm_,2)*pdf(AStr_,1) * PreFac
                
               if (RES(2,-1).gt.csmax(2,-1)) CSmax(2,-1) = RES(2,-1)
               if (RES(4,-3).gt.csmax(4,-3)) CSmax(4,-3) = RES(4,-3)
               if (RES(-1,2).gt.csmax(-1,2)) CSmax(-1,2) = RES(-1,2)
               if (RES(-3,4).gt.csmax(-3,4)) CSmax(-3,4) = RES(-3,4)

          ELSEIF( PROCESS.EQ.113 ) THEN      ! t-channel production of tbar H
              call EvalAmp_QQB_TBARHB(MomExt,LO_Res_Unpol)
               RES(-2,1) = LO_Res_Unpol(AUp_,Dn_)   * pdf(AUp_,1)*pdf(Dn_,2)   * PreFac
               RES(-4,3) = LO_Res_Unpol(AChm_,Str_) * pdf(AChm_,1)*pdf(Str_,2) * PreFac
               RES(1,-2) = LO_Res_Unpol(Dn_,AUp_)   * pdf(AUp_,2)*pdf(Dn_,1)   * PreFac
               RES(3,-4) = LO_Res_Unpol(Str_,AChm_) * pdf(AChm_,2)*pdf(Str_,1) * PreFac
                
               if (RES(-2,1).gt.csmax(-2,1)) CSmax(-2,1) = RES(-2,1)
               if (RES(-4,3).gt.csmax(-4,3)) CSmax(-4,3) = RES(-4,3)
               if (RES(1,-2).gt.csmax(1,-2)) CSmax(1,-2) = RES(1,-2)
               if (RES(3,-4).gt.csmax(3,-4)) CSmax(3,-4) = RES(3,-4)
          ENDIF


ENDIF! GENEVT 


RETURN
END FUNCTION













   
END MODULE

      


