MODULE ModCrossSection_TH
! Author: R. Rontsch, June 2015


 implicit none



 contains



FUNCTION EvalWeighted_TH(yRnd,VgsWgt)
! Routine for production of H(p3)+t(p4)+jet(p5)
! Top decays taken from implementation in MCFM, see hep-ph:/1204.1513
use ModKinematics
use ModParameters
use ModMisc
use ModTHiggs
#if compiler==1
    use ifport
#endif
implicit none
real(8) :: EvalWeighted_TH,yRnd(1:11),VgsWgt
real(8) :: Ehat,MH_Inv,eta1,eta2,ISFac,sHatJacobi,PreFac,FluxFac,PSWgt,PSWgt2,pdf(-6:6,1:2)
real(8) :: MomExt(1:4,1:8),LO_Res_Unpol(-6:6,-6:6),MuFac
integer :: NBin(1:NumHistograms),NHisto
logical :: applyPSCut

    
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
      PSWgt = PSWgt * PSWgt2
   ENDIF

   call Kinematics_TH(MomExt,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalWeighted_TH = 0d0
      return
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt

   IF( PROCESS.EQ.110 ) THEN
      call EvalAmp_QB_TH(MomExt,LO_Res_Unpol)
      EvalWeighted_TH = &
                      + LO_Res_Unpol(Up_,Bot_)   * ( pdf(Up_,1) *pdf(Bot_,2)  +  pdf(Chm_,1) *pdf(Bot_,2) )   &
                      + LO_Res_Unpol(Bot_,Up_)   * ( pdf(Bot_,1)*pdf(Up_,2)   +  pdf(Bot_,1) *pdf(Chm_,2) )   &
                      + LO_Res_Unpol(ADn_,Bot_)  * ( pdf(ADn_,1)*pdf(Bot_,2)  +  pdf(AStr_,1)*pdf(Bot_,2) )   &
                      + LO_Res_Unpol(Bot_,ADn_)  * ( pdf(Bot_,1)*pdf(ADn_,2)  +  pdf(Bot_,1) *pdf(AStr_,2))
   ELSEIF( PROCESS.EQ.111 ) THEN
      call EvalAmp_QbarBbar_TH(MomExt,LO_Res_Unpol)
      EvalWeighted_TH = &
                      + LO_Res_Unpol(Dn_,ABot_)  * ( pdf(Dn_,1)*pdf(ABot_,2)  + pdf(Str_,1)*pdf(ABot_,2) )    &
                      + LO_Res_Unpol(ABot_,Dn_)  * ( pdf(ABot_,1)*pdf(Dn_,2)  + pdf(ABot_,1)*pdf(Str_,2) )    &
                      + LO_Res_Unpol(AUp_,ABot_) * ( pdf(AUp_,1)*pdf(ABot_,2) + pdf(AChm_,1)*pdf(ABot_,2))    &
                      + LO_Res_Unpol(ABot_,AUp_) * ( pdf(ABot_,1)*pdf(AUp_,2) + pdf(ABot_,1)*pdf(AChm_,2))
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
 

   
END MODULE

      


