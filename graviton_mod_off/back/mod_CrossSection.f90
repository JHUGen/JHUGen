MODULE ModCrossSection
implicit none

contains


Function EvalCS(yRnd,VgsWgt)
use ModKinematics
use ModParameters
use ModGraviton
use ModHiggs
use ModZprime
use ModMisc
use ifport
implicit none
real(8) :: EvalCS,LO_Res_Unpol_old,LO_Res_Unpol,yRnd(1:22),VgsWgt
real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:11),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4)
logical :: applyPSCut
include 'csmaxvalue.f'

   EvalCS = 0d0
   call PDFMapping(11,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   EvalCounter = EvalCounter+1


   call EvalPhaseSpace_2to2Z(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   call EvalPhasespace_ZDecay(MomExt(1:4,3),yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
   call EvalPhasespace_ZDecay(MomExt(1:4,4),yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3

   call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)

   if( applyPSCut ) then
      EvalCS = 0d0
      return
   endif

   call setPDFs(eta1,eta2,m_Grav,pdf)
   FluxFac = 1d0/(2d0*EHat**2)




   if (PChannel.eq.0.or.PChannel.eq.2) then 
      PDFFac = pdf(0,1) * pdf(0,2)
      if (Process.eq.0) then 
      call EvalAmp_gg_H_ZZ( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      elseif(Process.eq.2) then
      call EvalAmp_gg_G_ZZ( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
!     call EvalAmp_gg_G_ZZ_old( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      endif

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2

   PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
   EvalCS = LO_Res_Unpol * PreFac

   

    if (EvalCS.gt.csmax_gg) csmax_gg = EvalCS

    endif 

   if (PChannel.eq.1.or.PChannel.eq.2) then 
      PDFFac = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
             + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
             + pdf(Bot_,1)*pdf(ABot_,2)                             &
             + pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
             + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
             + pdf(Bot_,2)*pdf(ABot_,1)

      if (Process.eq.1) then 
      call EvalAmp_qqb_Zprime_ZZ((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      elseif(Process.eq.2) then 
      call EvalAmp_qqb_G_ZZ((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      endif 

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2

   PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
   EvalCS = LO_Res_Unpol * PreFac

    if (EvalCS.gt.csmax_qq) csmax_qq = EvalCS

   endif




      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalCS*VgsWgt)
      enddo





RETURN
END FUNCTION


FUNCTION EvalCS_LO_ppllll(yRnd)
use ModKinematics
use ModParameters
use ModHiggs
use ModZprime
use ModGraviton
use ModMisc
use ifport
implicit none
real(8) :: EvalCS_LO_ppllll,LO_Res_Unpol_old,LO_Res_Unpol,yRnd(1:22),VgsWgt
real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:11),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4)
logical :: applyPSCut
real(8) :: CS_max, channel_ratio
real(8) :: oneovervolume
integer :: gluon,quark
include 'vegas_common.f'
include 'csmaxvalue.f'

   gluon = 0 
   quark = 0
   oneovervolume = one

   EvalCS_LO_ppllll = 0d0
   call PDFMapping(11,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   EvalCounter = EvalCounter+1


   call EvalPhaseSpace_2to2Z(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   call EvalPhasespace_ZDecay(MomExt(1:4,3),yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
   call EvalPhasespace_ZDecay(MomExt(1:4,4),yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3

   call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
   if( applyPSCut ) then
      EvalCS_LO_ppllll = 0d0
      return
   endif

   call setPDFs(eta1,eta2,m_Grav,pdf)
   FluxFac = 1d0/(2d0*EHat**2)


    if (fix_channels_ratio .and. Process.eq.2) then
       channel_ratio = adj_par*csmax_qq/(csmax_qq*adj_par+csmax_gg)  ! fix qq/ total (gg + qq) 
    else
      channel_ratio = adj_par*csmax_qq/(csmax_qq*adj_par+csmax_gg)
    endif  
   

   if( yRnd(10).gt.channel_ratio ) then

      gluon = 1

      CS_max = csmax_gg


      PDFFac = pdf(0,1) * pdf(0,2)

      if (Process.eq.0) then 
   call EvalAmp_gg_H_ZZ( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      elseif(Process.eq.2) then 
   call EvalAmp_gg_G_ZZ( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      endif 
!       call EvalAmp_gg_G_ZZ_old( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
 
   else

      quark = 1

      CS_max = csmax_qq


      PDFFac = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
             + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
             + pdf(Bot_,1)*pdf(ABot_,2)                             &
             + pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
             + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
             + pdf(Bot_,2)*pdf(ABot_,1)
      if (Process.eq.1) then 
  call EvalAmp_qqb_Zprime_ZZ((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)

      elseif(Process.eq.2) then 
      call EvalAmp_qqb_G_ZZ((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
      endif 

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2
   endif

   PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
   EvalCS_LO_ppllll = LO_Res_Unpol * PreFac




      if( EvalCS_LO_ppllll.gt. CS_max) then
          print *, "CS_max is too small. Adjust CS_max!",EvalCS_LO_ppllll, CS_max
          stop
      endif



      if( EvalCS_LO_ppllll .gt. yRnd(9)*CS_max ) then
           do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
           enddo
          AccepCounter = AccepCounter + 1
          if (gluon.eq.1)           AccepCounter_g = AccepCounter_g + 1
          if (quark.eq.1)           AccepCounter_q = AccepCounter_q + 1
!          call WriteOutEvent(eta1,eta2,(/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/))
      else
          RejeCounter = RejeCounter + 1
      endif

RETURN
END FUNCTION



END MODULE ModCrossSection


