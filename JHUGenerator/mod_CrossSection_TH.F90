MODULE ModCrossSection_TH
! Author: R. Rontsch, June 2015
  use ModTopDecay
  implicit none


integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)  


contains

  FUNCTION EvalWeighted_TH(yRnd,VgsWgt)
! Routine for production of H(p3)+t(p4)+jet(p5)
! Top decays taken from implementation in MCFM, see hep-ph:/1204.1513
    use ModKinematics
    use ModParameters
!    use ModTTBHiggs
    use ModMisc
#if compiler==1
    use ifport
#endif
    implicit none
    real(8) :: EvalWeighted_TH,yRnd(1:11),VgsWgt
    real(8) :: Ehat,MH_Inv,eta1,eta2,ISFac,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles,PSWgt,PSWgt2,PSWgt3,pdf(-6:6,1:2)
    real(8) :: MomExt(1:4,1:11),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol,s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4),MuFac
    complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),CoupFac,decay_amp(1:2)
    real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
    integer :: NBin(1:NumHistograms),NHisto,j
    logical :: applyPSCut
    
    EvalWeighted_TH = 0d0

    call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    MH_Inv = M_Reso
   if( EHat.le.m_Top+MH_Inv ) then
      EvalWeighted_TH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3ArbMass(EHat,(/MH_Inv,M_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

!   ISFac = MomCrossing(MomExt)
!   IsFac=1/4d0/9d0
   IsFac=SpinAvg*QuarkColAvg**2

   MuFac=(M_Top + M_Reso)/4d0
   call setPDFs(eta1,eta2,MuFac,pdf)


! couplings
  CoupFac=2d0*gwsq/vev*ci
!   CoupFac=1d0

   ColFac=9d0
   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),MomExt(1:4,6:8),PSWgt2)
      PSWgt = PSWgt * PSWgt2
! usual decay top(p4) --> b(p6) + e+(p7) + nu(p8)
!      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
   ENDIF

   call Kinematics_TH(MomExt(1:4,1:11),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalWeighted_TH = 0d0
      return
   endif


!  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt



! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,5)*MomExt(1,4)-MomExt(2,5)*MomExt(2,4)-MomExt(3,5)*MomExt(3,4)-MomExt(4,5)*MomExt(4,4)
      p4Dp7=MomExt(1,7)*MomExt(1,4)-MomExt(2,7)*MomExt(2,4)-MomExt(3,7)*MomExt(3,4)-MomExt(4,7)*MomExt(4,4)
      p2Dp3=MomExt(1,2)*MomExt(1,3)-MomExt(2,2)*MomExt(2,3)-MomExt(3,2)*MomExt(3,3)-MomExt(4,2)*MomExt(4,3)
      MomExtFlat(1,1:4)=MomExt(1:4,1)
      MomExtFlat(2,1:4)=MomExt(1:4,2)
      MomExtFlat(3,1:4)=m_Reso**2/2d0/p2Dp3*MomExt(1:4,2)
      MomExtFlat(4,1:4)=MomExt(1:4,3)-MomExtFlat(3,1:4)
      MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,5)/2d0/p4Dp5
      MomExtFlat(6,1:4)=MomExt(1:4,4)-MomExtFlat(5,1:4)
      MomExtFlat(7,1:4)=MomExt(1:4,5)

    ! use different flattened momenta for top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,7)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,4)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,6)
         MomExtFlatDK(9,1:4)=MomExt(1:4,7)
         MomExtFlatDK(10,1:4)=MomExt(1:4,8)
      ENDIF
      

! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call tdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
         call ubhtdamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Up_,Bot_,1:2))
         call ubhtdamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Bot_,Up_,1:2))        
         call ubhtdamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,LOAmp(ADn_,Bot_,1:2))
         call ubhtdamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,LOAmp(Bot_,ADn_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp=LOAmp*CoupFac


      LO_Res_UnPol= &
           + (abs(LOAmp(Up_,Bot_,1))**2+abs(LOAmp(Up_,Bot_,2))**2)   * (pdf(Up_,1)*pdf(Bot_,2)  + pdf(Chm_,1)*pdf(Bot_,2)) &
           + (abs(LOAmp(Bot_,Up_,1))**2+abs(LOAmp(Bot_,Up_,2))**2)   * (pdf(Bot_,1)*pdf(Up_,2)  + pdf(Bot_,1)*pdf(Chm_,2)) &
           + (abs(LOAmp(ADn_,Bot_,1))**2+abs(LOAmp(ADn_,Bot_,2))**2) * (pdf(ADn_,1)*pdf(Bot_,2) + pdf(AStr_,1)*pdf(Bot_,2)) &
           + (abs(LOAmp(Bot_,ADn_,1))**2+abs(LOAmp(Bot_,ADn_,2))**2) * (pdf(Bot_,1)*pdf(ADn_,2) + pdf(Bot_,1)*pdf(AStr_,2))

      LO_Res_UnPol = LO_Res_UnPol * ColFac * ISFac
      EvalWeighted_TH = LO_Res_Unpol * PreFac


   AccepCounter=AccepCounter+1

   if( IsNan(EvalWeighted_TH) ) then
        print *, "NAN:",EvalWeighted_TH
        print *, yRnd(:)
        print *, LO_Res_UnPol

        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalWeighted_TH = 0d0
!         pause                                                                                                                                 
        return
   endif

   if( writeWeightedLHE ) then 
        call Error("WriteLHE not yet supported for t+H")
     endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalWeighted_TH*VgsWgt)                
   enddo
   EvalCounter = EvalCounter + 1



 end FUNCTION EvalWeighted_TH



 FUNCTION EvalWeighted_TBH(yRnd,VgsWgt)
! Routine for production of H(p3)+tbar(p4)+jet(p5)
! Top decays taken from implementation in MCFM, see hep-ph:/1204.1513
    use ModParameters
    use ModKinematics
    use ModMisc
    implicit none
    real(8) :: EvalWeighted_TBH,yRnd(1:11),VgsWgt
    real(8) :: Ehat,MH_Inv,eta1,eta2,ISFac,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles,PSWgt,PSWgt2,PSWgt3,pdf(-6:6,1:2)
    real(8) :: MomExt(1:4,1:11),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol,s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4),MuFac
    complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),CoupFac,decay_amp(1:2)
    real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
    integer :: NBin(1:NumHistograms),NHisto,j,k
    logical :: applyPSCut



    EvalWeighted_TBH = 0d0

    call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    MH_Inv = M_Reso

   if( EHat.le.m_Top+MH_Inv ) then
      EvalWeighted_TBH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3ArbMass(EHat,(/MH_Inv,M_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

!   ISFac = MomCrossing(MomExt)
   IsFac=SpinAvg*QuarkColAvg**2
   MuFac=(M_Top + M_Reso)/4d0
   call setPDFs(eta1,eta2,MuFac,pdf)

!couplings
   CoupFac=2d0*gwsq/vev*ci
!   CoupFac=1d0

   ColFac=9d0

   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),MomExt(1:4,6:8),PSWgt2)
      PSWgt = PSWgt * PSWgt2
! usual decay atop(p4) --> ab(p6) + e-(p7) + nubar(p8)
!      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
   ENDIF



   call Kinematics_TH(MomExt(1:4,1:11),applyPSCut,NBin)
   if( applyPSCut ) then
      EvalWeighted_TBH = 0d0
      return
   endif
   

  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt 

! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,5)*MomExt(1,4)-MomExt(2,5)*MomExt(2,4)-MomExt(3,5)*MomExt(3,4)-MomExt(4,5)*MomExt(4,4)
   p4Dp7=MomExt(1,7)*MomExt(1,4)-MomExt(2,7)*MomExt(2,4)-MomExt(3,7)*MomExt(3,4)-MomExt(4,7)*MomExt(4,4)
   p2Dp3=MomExt(1,2)*MomExt(1,3)-MomExt(2,2)*MomExt(2,3)-MomExt(3,2)*MomExt(3,3)-MomExt(4,2)*MomExt(4,3)
   MomExtFlat(1,1:4)=MomExt(1:4,1)
   MomExtFlat(2,1:4)=MomExt(1:4,2)
   MomExtFlat(3,1:4)=m_Reso**2/2d0/p2Dp3*MomExt(1:4,2)
   MomExtFlat(4,1:4)=MomExt(1:4,3)-MomExtFlat(3,1:4)
   MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,5)/2d0/p4Dp5
   MomExtFlat(6,1:4)=MomExt(1:4,4)-MomExtFlat(5,1:4)
   MomExtFlat(7,1:4)=MomExt(1:4,5)

    ! use different flattened momenta for anti-top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,7)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,4)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,6)
         MomExtFlatDK(9,1:4)=MomExt(1:4,7)
         MomExtFlatDK(10,1:4)=MomExt(1:4,8)

      ENDIF



! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call atdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
      call dbbarhtbaruamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Dn_,ABot_,1:2))
      call dbbarhtbaruamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(ABot_,Dn_,1:2))
      call dbbarhtbaruamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,LOAmp(AUp_,ABot_,1:2))
      call dbbarhtbaruamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,LOAmp(ABot_,AUp_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp=LOAmp*CoupFac


      LO_Res_UnPol= &
           + (abs(LOAmp(Dn_,ABot_,1))**2+abs(LOAmp(Dn_,ABot_,2))**2)   * (pdf(Dn_,1)*pdf(ABot_,2)  + pdf(Str_,1)*pdf(ABot_,2)) &
           + (abs(LOAmp(ABot_,Dn_,1))**2+abs(LOAmp(ABot_,Dn_,2))**2)   * (pdf(ABot_,1)*pdf(Dn_,2)  + pdf(ABot_,1)*pdf(Str_,2)) &
           + (abs(LOAmp(AUp_,ABot_,1))**2+abs(LOAmp(AUp_,ABot_,2))**2) * (pdf(AUp_,1)*pdf(ABot_,2) + pdf(AChm_,1)*pdf(ABot_,2)) &
           + (abs(LOAmp(ABot_,AUp_,1))**2+abs(LOAmp(ABot_,AUp_,2))**2) * (pdf(ABot_,1)*pdf(AUp_,2) + pdf(ABot_,1)*pdf(AChm_,2))

      LO_Res_UnPol = LO_Res_UnPol * ColFac * ISFac 
      EvalWeighted_TBH = LO_Res_Unpol * PreFac

   
   if( IsNan(EvalWeighted_TBH) ) then
        print *, "NAN:",EvalWeighted_TBH
        print *, yRnd(:)
        print *, LO_Res_UnPol

        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalWeighted_TBH = 0d0
!         pause                                                                                                                                                                 
        return
   endif


   AccepCounter=AccepCounter+1
   if( writeWeightedLHE ) then 
        call Error("WriteLHE not yet supported for tb+H")
   endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalWeighted_TBH*VgsWgt)
   enddo
   EvalCounter = EvalCounter + 1



 end FUNCTION EvalWeighted_TBH



   
      subroutine ubhtdamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production u(p1)+b(p2)->H(p3)+t(p4)+d(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        use ModParameters
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s24,s34,s15,mt,mw

        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top

! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)        
        
        ampw(1) = 1/( - mw**2 + s24)/( - mw**2 + s15)*za(p5,e4)*zb(p1,&
     & p2)*mw**2 + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)/(zb(k4,e4)&
     & )*za(p5,e3)*zb(p2,k4)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + s24)&
     & /( - mw**2 + s15)/(zb(k4,e4))*za(p5,k3)*zb(p2,k4)*zb(k3,p1)*&
     & mt**2
        
        ampw(2) = 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)*za(p5,e3)*&
     & zb(p2,e4)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + &
     & s15)*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*mt + 1/( - mw**2 + s24)/( - &
     & mw**2 + s15)/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*mw**2*mt
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p5,e4&
     & )*zb(p1,p2)*mt*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15&
     & )/(zb(k4,e4))*za(p5,p1)*zb(p1,p2)*zb(p1,k4)*mt*vev*KL - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s15)/(zb(k4,e4))*za(p5,p2)*zb(p1,&
     & p2)*zb(p2,k4)*mt*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p5,p1&
     & )*zb(p1,p2)*zb(p1,e4)*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s15)*za(p5,p2)*zb(p1,p2)*zb(p2,e4)*vev*KL - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s15)/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*&
     & mt**2*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      end subroutine ubhtdamp
        
    
      subroutine dbbarhtbaruamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production d(p1)+bbar(p2)->H(p3)+tbar(p4)+u(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        use ModParameters
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s24,s34,s15,mt,mw
    
        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)
        
        ampw(1) = 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)*za(p2,e4)*&
     & za(p5,e3)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + &
     & s15)*za(p2,e4)*za(p5,k3)*zb(k3,p1)*mt - 1/( - mw**2 + s24)/( - &
     & mw**2 + s15)/(zb(k4,e4))*za(p2,p5)*zb(p1,k4)*mw**2*mt
        
        ampw(2) =  - 1/( - mw**2 + s24)/( - mw**2 + s15)*za(p2,p5)*zb(&
     & p1,e4)*mw**2 + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)/(za(k4,&
     & e4))*za(p2,k4)*za(p5,e3)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + &
     & s24)/( - mw**2 + s15)/(za(k4,e4))*za(p2,k4)*za(p5,k3)*zb(k3,p1)*&
     & mt**2
        
        ampt(1) = 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p2,p5)*&
     & za(e4,p2)*zb(p2,p1)*vev*KR + 1d0/2d0/( - mt**2 + s34)/( - mw**2&
     &  + s15)*za(p2,p5)*za(e4,p5)*zb(p5,p1)*vev*KR + 1d0/2d0/( - mt**2&
     &  + s34)/( - mw**2 + s15)/(zb(k4,e4))*za(p2,p5)*zb(p1,k4)*mt**2*&
     & vev*KL
        
        ampt(2) = 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p2,p5)*&
     & zb(p1,e4)*mt*vev*KL + 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)&
     & /(za(k4,e4))*za(p2,p5)*za(k4,p2)*zb(p2,p1)*mt*vev*KR + 1d0/2d0/(&
     &  - mt**2 + s34)/( - mw**2 + s15)/(za(k4,e4))*za(p2,p5)*za(k4,p5)&
     & *zb(p5,p1)*mt*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      end subroutine dbbarhtbaruamp
    
    
    
    SUBROUTINE TDECAY(k4,e4,b,ep,nu,za,zb,dkamp)
! top decay routine, taken from MCFM, see hep-ph:/1204.1513
       use ModParameters
       implicit none
       integer :: k4,e4,b,ep,nu
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp
     
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = za(b,nu)*zb(ep,e4)
       dkamp(2) = m_top * za(b,nu)*zb(ep,k4)/za(e4,k4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * g2_weak
   
   
     end SUBROUTINE TDECAY
   
    
    SUBROUTINE ATDECAY(k4,e4,bbar,em,nubar,za,zb,dkamp)
! anti-top decay routine, taken from MCFM, see hep-ph:/1204.1513
       use ModParameters
       implicit none
       integer :: k4,e4,bbar,em,nubar
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
       complex(8) :: WProp
     
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = -m_top * zb(bbar,nubar)*za(em,k4)/zb(e4,k4)
       dkamp(2) = -zb(bbar,nubar)*za(em,e4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * g2_weak
   
   
     end SUBROUTINE ATDECAY
   


end MODULE ModCrossSection_TH

      


