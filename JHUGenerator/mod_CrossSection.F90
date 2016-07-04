MODULE ModCrossSection
implicit none

contains


 Function EvalWeighted(yRnd,VgsWgt)    ! this is a function which is only for computations
 use ModKinematics                     ! with weighted events
 use ModParameters
 use ModGraviton
 use ModHiggs
 use ModZprime
 use ModMisc
#if compiler==1
 use ifport
#endif
 implicit none
 real(8) :: EvalWeighted,LO_Res_Unpol_old,LO_Res_Unpol,yRnd(1:22),VgsWgt
 real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
 real(8) :: pdf(-6:6,1:2)
 integer :: NBin(1:11),NHisto,i,MY_IDUP(1:9), ICOLUP(1:2,1:9)
 real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
 real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4)
 real(8) :: EZ,  EZ1, EZ2, pz, xmax, xx1, xx2
 real(8) :: pz12, MomExt_f(1:4,1:4), MomDK_f(1:4,1:4)
 real(8) :: MZ1,MZ2,EZ_max,dr, MG, yz1,yz2
 real(8) :: offzchannel
 logical :: applyPSCut
 include 'csmaxvalue.f'

    EvalWeighted = 0d0
    if( OffShellReson ) then
      call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    else
!       call PDFMapping(11,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
      call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    endif
    EvalCounter = EvalCounter+1




!    particle associations:
!    
!    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor

    call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9))


  yz1 = yRnd(10)
  yz2 = yRnd(11)
  offzchannel = yRnd(12) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell 
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

            EZ_max = EHat - MZ1*0.99
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

        elseif(M_Reso.lt.2d0*M_V) then
            if (offzchannel.le.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
                MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
                MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
                sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/(  &
                1d0/((MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )     &
                + 1d0/((MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
            elseif(offzchannel.gt.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
                MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
                MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
                sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/( &
                1d0/((MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )    &
                + 1d0/((MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
            endif
       endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
        MZ1 = M_V
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat - MZ1*0.99
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
        else
            MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
        endif

  elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
        MZ2 = M_V
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat - MZ2*0.99
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
         else
            MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
        endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
        MZ1 = M_V
        MZ2 = M_V
  endif


    if( EHat.lt.MZ1+MZ2 ) then
      EvalWeighted = 0d0
      return
    endif


    call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
    call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
    if( DecayMode1.ne.7 ) then ! don't decay the photon
        call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
        call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
        PSWgt = PSWgt * PSWgt2*PSWgt3
    else
        MomDK(1:4,1) = MomExt(1:4,3)
        MomDK(1:4,2) = 0d0
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
    endif


    if( (OffShellV1).or.(OffShellV2).or.(DecayMode1.eq.7) ) then
        call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    else
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalWeighted = 0d0
      return
    endif

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)

   if (PChannel.eq.0.or.PChannel.eq.2) then
      PDFFac = pdf(0,1) * pdf(0,2)
      if (Process.eq.0) then
            call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
            call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      endif

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalWeighted = LO_Res_Unpol * PreFac
    endif



   if (PChannel.eq.1.or.PChannel.eq.2) then
      PDFFac = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
             + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
             + pdf(Bot_,1)*pdf(ABot_,2)                             &
             + pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
             + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
             + pdf(Bot_,2)*pdf(ABot_,1)

      if (Process.eq.1) then
         call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
         call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      endif
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2

      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalWeighted = LO_Res_Unpol * PreFac
   endif


      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalWeighted*VgsWgt)
      enddo

RETURN
END FUNCTION









FUNCTION EvalUnWeighted(yRnd,genEvt,RES)
use ModKinematics
use ModParameters
use ModHiggs
use ModZprime
use ModGraviton
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
real(8) :: RES(-5:5,-5:5)
real(8) :: EvalUnWeighted,LO_Res_Unpol_old,LO_Res_Unpol,yRnd(1:22),VgsWgt
real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:11),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomExt_f(1:4,1:4),MomDK_f(1:4,1:4)
logical :: applyPSCut,genEvt
real(8) :: CS_max, channel_ratio
real(8) :: oneovervolume, bound(1:11), sumtot,yz1,yz2,EZ_max,dr,MZ1,MZ2
integer :: parton(-5:5,-5:5), i1, ifound, i2, MY_IDUP(1:9), ICOLUP(1:2,1:9)
real(8)::ntRnd,ZMass(1:2)
real(8) :: offzchannel
include 'vegas_common.f'
include 'csmaxvalue.f'


   parton = 0
   oneovervolume = one
   ICOLUP(1:2,1:9) = 0
   EvalUnWeighted = 0d0
   RES = 0d0

   if(OffShellReson) then
      call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   else
      call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   endif
   EvalCounter = EvalCounter+1

   MY_IDUP(3)= 0

!    particle associations:
!    
!    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor
   call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9))


  yz1 = yRnd(10)
  yz2 = yRnd(11)
  offzchannel = yRnd(15) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

            EZ_max = EHat - MZ1*0.99
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

        elseif(M_Reso.lt.2d0*M_V) then
            if (offzchannel.le.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
                MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
                MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
                sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/(  &
                1d0/((MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )     &
                + 1d0/((MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
            elseif(offzchannel.gt.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
                MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
                MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
                sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/( &
                1d0/((MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )    &
                + 1d0/((MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
            endif
       endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
        MZ1 = M_V
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat - MZ1*0.99
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
        else
            MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
	endif

  elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
        MZ2 = M_V
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat - MZ2*0.99
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
         else
            MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
	endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
        MZ1 = M_V
        MZ2 = M_V
  endif


    if( MZ1+MZ2.gt.EHat ) then
      EvalUnWeighted = 0d0
      RejeCounter = RejeCounter + 1
      return
    endif


   call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   if( DecayMode1.ne.7 ) then ! don't decay the photon
      call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
      call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
      PSWgt = PSWgt * PSWgt2*PSWgt3
    else
        MomDK(1:4,1) = MomExt(1:4,3)
        MomDK(1:4,2) = 0d0
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
   endif


    if( (OffShellV1).or.(OffShellV2).or.(DecayMode1.eq.7) ) then
        call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    else
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalUnWeighted = 0d0
      return
    endif

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)

IF( GENEVT ) THEN
    sumtot = csmax(0,0)+csmax(-5,5)+csmax(-4,4)+csmax(-3,3)+csmax(-2,2)+csmax(-1,1)  &
    +csmax(1,-1)+csmax(2,-2)+csmax(3,-3)+csmax(4,-4)+csmax(5,-5)

    bound(1)  = csmax(0,0)/sumtot
    bound(2)  = bound(1) + csmax(-5,5)/sumtot
    bound(3)  = bound(2) + csmax(-4,4)/sumtot
    bound(4)  = bound(3) + csmax(-3,3)/sumtot
    bound(5)  = bound(4) + csmax(-2,2)/sumtot
    bound(6)  = bound(5) + csmax(-1,1)/sumtot
    bound(7)  = bound(6) + csmax(1,-1)/sumtot
    bound(8)  = bound(7) + csmax(2,-2)/sumtot
    bound(9)  = bound(8) + csmax(3,-3)/sumtot
    bound(10) = bound(9) + csmax(4,-4)/sumtot
    bound(11) = one

    if (yRnd(13).lt.bound(1)) ifound = 1
    do i1=1,10
       if (yRnd(13).gt.bound(i1).and.yRnd(13).lt.bound(i1+1)) ifound = i1+1
    enddo
!    if (fix_channels_ratio .and. Process.eq.2) then
!       channel_ratio = adj_par*csmax_qq/(csmax_qq*adj_par+csmax_gg)  ! fix qq/ total (gg + qq)
!    else
!      channel_ratio = adj_par*csmax_qq/(csmax_qq*adj_par+csmax_gg)
!    endif


   if(ifound.eq.1 ) then
      parton(0,0) = 1
      CS_max = csmax(0,0)*adj_par   ! this is necessary here for the follow up
      MY_IDUP(1:2)=(/Glu_,Glu_/)
      ICOLUP(1:2,1) = (/501,502/)
      ICOLUP(1:2,2) = (/502,501/)
      PDFFac = pdf(0,1) * pdf(0,2)

      if (Process.eq.0) then
            call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
            call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      endif
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2

   else

      if (ifound.eq.2) then
          PDFFac = pdf(Bot_,2)*pdf(ABot_,1)
          i2 = -5
          MY_IDUP(1:2)=(/ABot_,Bot_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
      elseif(ifound.eq.3) then
          PDFFac = pdf(Chm_,2)*pdf(AChm_,1)
          i2 = -4
          MY_IDUP(1:2)=(/AChm_,Chm_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
      elseif(ifound.eq.4) then
          PDFFac = pdf(Str_,2)*pdf(AStr_,1)
          i2 = -3
          MY_IDUP(1:2)=(/AStr_,Str_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
      elseif(ifound.eq.5) then
          PDFFac = pdf(Up_,2) *pdf(AUp_,1)
          i2 = -2
          MY_IDUP(1:2)=(/AUp_,Up_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
      elseif(ifound.eq.6) then
          PDFFac = pdf(Dn_,2) *pdf(ADn_,1)
          i2 = -1
          MY_IDUP(1:2)=(/ADn_,Dn_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
      elseif (ifound.eq.7) then
          PDFFac = pdf(Dn_,1) *pdf(ADn_,2)
          i2 = 1
          MY_IDUP(1:2)=(/Dn_,ADn_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
      elseif(ifound.eq.8) then
          PDFFac = pdf(Up_,1) *pdf(AUp_,2)
          i2 = 2
          MY_IDUP(1:2)=(/Up_,AUp_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
      elseif(ifound.eq.9) then
          PDFFac = pdf(Str_,1)*pdf(AStr_,2)
          i2 = 3
          MY_IDUP(1:2)=(/Str_,AStr_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
      elseif(ifound.eq.10) then
          PDFFac = pdf(Chm_,1)*pdf(AChm_,2)
          i2 = 4
          MY_IDUP(1:2)=(/Chm_,AChm_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
       elseif(ifound.eq.11) then
          PDFFac = pdf(Bot_,1)*pdf(ABot_,2)
          i2 = 5
          MY_IDUP(1:2)=(/Bot_,ABot_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
       endif
       parton(i2,-i2) = 1
       CS_max = csmax(i2,-i2)

      if (Process.eq.1) then
         call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
         call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      endif

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2
   endif

   PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
   if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
   if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
   EvalUnWeighted = LO_Res_Unpol * PreFac

      if( EvalUnWeighted.gt. CS_max) then
          print *, "CS_max is too small. Adjust CS_max!",EvalUnWeighted, CS_max
!          stop
      endif

      if( EvalUnWeighted .gt. yRnd(14)*CS_max ) then
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1
         AccepCounter_part = AccepCounter_part  + parton
         if( (OffShellV1).or.(OffShellV2).or.(DecayMode1.eq.7) ) then
              call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
          else
              call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
         endif
      else
          RejeCounter = RejeCounter + 1
      endif


ELSE! GENEVT


   if (PChannel.eq.0.or.PChannel.eq.2) then
      PDFFac = pdf(0,1) * pdf(0,2)
      if (Process.eq.0) then
         call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
         call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      endif
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2

      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc

      EvalUnWeighted = LO_Res_Unpol * PreFac
      RES(0,0) = EvalUnWeighted

      if (EvalUnWeighted.gt.csmax(0,0)) then
          csmax(0,0) = EvalUnWeighted
      endif
   endif

   if (PChannel.eq.1.or.PChannel.eq.2) then
      if (Process.eq.1) then
         call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
         call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
      endif

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2
      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt *   SymmFac

      do i1 = -5,5
         if (i1.eq.-5) then
          PDFFac = pdf(Bot_,2)*pdf(ABot_,1)
         elseif(i1.eq.-4) then
          PDFFac = pdf(Chm_,2)*pdf(AChm_,1)
         elseif(i1.eq.-3) then
          PDFFac = pdf(Str_,2)*pdf(AStr_,1)
         elseif(i1.eq.-2) then
          PDFFac = pdf(Up_,2) *pdf(AUp_,1)
         elseif(i1.eq.-1) then
          PDFFac = pdf(Dn_,2) *pdf(ADn_,1)
         elseif (i1.eq.0) then
          PDFFac = 0d0
         elseif (i1.eq.1) then
          PDFFac = pdf(Dn_,1) *pdf(ADn_,2)
         elseif (i1.eq.2) then
          PDFFac = pdf(Up_,1) *pdf(AUp_,2)
         elseif(i1.eq.3) then
          PDFFac = pdf(Str_,1)*pdf(AStr_,2)
         elseif(i1.eq.4) then
          PDFFac = pdf(Chm_,1)*pdf(AChm_,2)
          elseif(i1.eq.5) then
          PDFFac = pdf(Bot_,1)*pdf(ABot_,2)
          endif

          EvalUnWeighted = LO_Res_Unpol * PreFac *PDFFac
          RES(i1,-i1) = EvalUnWeighted
          if (EvalUnWeighted.gt.csmax(i1,-i1)) csmax(i1,-i1) = EvalUnWeighted
      enddo
   endif

ENDIF! genEvt


RETURN
END FUNCTION

















! subroutine EvalCS_un(yRnd,RES)
! use ModKinematics
! use ModParameters
! use ModGraviton
! use ModHiggs
! use ModZprime
! use ModMisc
! #if compiler==1
! use ifport
! #endif
! implicit none
! real(8) :: RES(-5:5,-5:5)
! real(8) :: EvalCS,LO_Res_Unpol_old,LO_Res_Unpol,yRnd(1:22),VgsWgt
! real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
! real(8) :: pdf(-6:6,1:2)
! integer :: NBin(1:11),NHisto,i, i1,Z1DKFlavor,Z2DKFlavor
! real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
! real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4)
! real(8) :: MomExt_f(1:4,1:4), MomDK_f(1:4,1:4),yz1,yz2,EZ_max,dr,MZ1,MZ2
! real(8) :: offzchannel
! logical :: applyPSCut
! include 'csmaxvalue.f'
!
!    RES = 0d0
!
!
!    EvalCS = 0d0
!    if(OffShellReson) then
!       call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    else
!       call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    endif
!    EvalCounter = EvalCounter+1
!
!    if( DecayMode.eq.0 ) then
!       Z1DKFlavor = ElM_
!       Z2DKFlavor = MuM_
!       yz1 = yRnd(10)
!       yz2 = yRnd(11)
!    offzchannel = yRnd(12) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
!    elseif( DecayMode.eq.1 ) then
!       Z1DKFlavor = ElM_
!       Z2DKFlavor = ZQBranching( yRnd(12) )
!       yz1 = yRnd(10)
!       yz2 = yRnd(11)
!     offzchannel = yRnd(13) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
!    endif
!
!
! !---- new stuff
!
!
!
!     if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
!
!         if(M_Reso.gt.2d0*M_V) then
!
!
!         EZ_max = EHat
!         dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
!         MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
!         sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
!
!         EZ_max = EHat - MZ1*0.99
!         dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
!         MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
!         sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
!
!         elseif(M_Reso.lt.2d0*M_V) then
!
!
!         if (offzchannel.le.0.5d0) then
!
!         EZ_max = EHat
!         dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
!         MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
!       MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
!
!         sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/(  &
!         1d0/((MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )     &
!         + 1d0/((MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
!
!          sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
!
!         elseif(offzchannel.gt.0.5d0) then
!
!         EZ_max = EHat
!         dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
!         MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
!
!       MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz1)))
!
!
!         sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/( &
!         1d0/((MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )    &
!         + 1d0/((MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
!
!          sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
!
!
!       endif
!       endif
!
!         elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
!
!         MZ1 = M_V
!
!         if(M_Reso.gt.2d0*M_V) then
!         EZ_max = EHat - MZ1*0.99
!         dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
!         MZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
!         sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
!          else
!        MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
!         sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
!
!       endif
!
!         elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
!
!         MZ2 = M_V
!
!         if(M_Reso.gt.2d0*M_V) then
!         EZ_max = EHat - MZ2*0.99
!         dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
!         MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
!         sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)*( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )
!
!          else
!        MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
!         sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
!
!       endif
!
!        elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
!
!                 MZ1 = M_V
!                 MZ2 = M_V
!
!
!          endif
!
!
!
!
!     if( MZ1+MZ2.gt.EHat ) then!NEW
!       EvalCS = 0d0
!       return
!     endif
!
!
!
!    call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
!    call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
!    call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
!    call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
!    PSWgt = PSWgt * PSWgt2*PSWgt3
!
!     if( OffShellV1.or.OffShellV2 ) then!NEW
!         call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
!     else
!         call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
!         call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
!     endif
!
!
!    if( applyPSCut ) then
!       EvalCS = 0d0
!       return
!    endif
!
!    call setPDFs(eta1,eta2,Mu_Fact,pdf)
!    FluxFac = 1d0/(2d0*EHat**2)
!
!
!    if (PChannel.eq.0.or.PChannel.eq.2) then
!       PDFFac = pdf(0,1) * pdf(0,2)
!       if (Process.eq.0) then
!       call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),Z1DKFlavor,Z2DKFlavor,LO_Res_Unpol)
!       elseif(Process.eq.2) then
!       call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),Z1DKFlavor,Z2DKFlavor,LO_Res_Unpol)
! !     call EvalAmp_gg_G_VV_old( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),LO_Res_Unpol)
!       endif
!
!       LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
!
!    PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
!
!     EvalCS = LO_Res_Unpol * PreFac
!     RES(0,0) = EvalCS
!
!     if (EvalCS.gt.csmax(0,0)) then
!         csmax(0,0) = EvalCS
! !     print *, offzchannel,MZ1, MZ2, LO_Res_Unpol,sHatJacobi,  csmax(0,0)
! !     pause
!     endif
!
!     endif
!
!    if (PChannel.eq.1.or.PChannel.eq.2) then
!
!
!       if (Process.eq.1) then
!       call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),Z1DKFlavor,Z2DKFlavor,LO_Res_Unpol)
!       elseif(Process.eq.2) then
!       call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),Z1DKFlavor,Z2DKFlavor,LO_Res_Unpol)
!       endif
!
!       LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2
!       PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt *   SymmFac
!
!       do i1 = -5,5
!
!          if (i1.eq.-5) then
!           PDFFac = pdf(Bot_,2)*pdf(ABot_,1)
!          elseif(i1.eq.-4) then
!           PDFFac = pdf(Chm_,2)*pdf(AChm_,1)
!          elseif(i1.eq.-3) then
!           PDFFac = pdf(Str_,2)*pdf(AStr_,1)
!          elseif(i1.eq.-2) then
!           PDFFac = pdf(Up_,2) *pdf(AUp_,1)
!          elseif(i1.eq.-1) then
!           PDFFac = pdf(Dn_,2) *pdf(ADn_,1)
!          elseif (i1.eq.0) then
!           PDFFac = 0d0
!          elseif (i1.eq.1) then
!           PDFFac = pdf(Dn_,1) *pdf(ADn_,2)
!          elseif (i1.eq.2) then
!           PDFFac = pdf(Up_,1) *pdf(AUp_,2)
!          elseif(i1.eq.3) then
!           PDFFac = pdf(Str_,1)*pdf(AStr_,2)
!          elseif(i1.eq.4) then
!           PDFFac = pdf(Chm_,1)*pdf(AChm_,2)
!           elseif(i1.eq.5) then
!           PDFFac = pdf(Bot_,1)*pdf(ABot_,2)
!           endif
!
!    EvalCS = LO_Res_Unpol * PreFac *PDFFac
!    RES(i1,-i1) = EvalCS
!
!     if (EvalCS.gt.csmax(i1,-i1)) csmax(i1,-i1) = EvalCS
!
!     enddo
!
!    endif
!
! ! if(EvalCS.lt.minCS) minCS=EvalCS
! ! if(EvalCS.gt.maxCS) maxCS=EvalCS
! ! avgCS = avgCS + EvalCS
! RETURN
! END subroutine

END MODULE ModCrossSection



