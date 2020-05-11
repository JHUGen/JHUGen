MODULE ModCrossSection
use ModHashCollection
implicit none
integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

 CONTAINS






 ! This function handles weighted and unweighted events for PP-->spin0,1,2-->VV
 FUNCTION EvalWeighted(yRnd,VgsWgt)
 use ModKinematics
 use ModParameters
 use ModGraviton
 use ModHiggs
 use ModZprime
 use ModMisc
#if compiler==1
 use ifport
#endif
 implicit none
 real(8) :: EvalWeighted,LO_Res_Unpol,yRnd(1:22),VgsWgt
 real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac,m1ffwgt,m2ffwgt,FinalStateWeight
 real(8) :: pdf(-6:6,1:2)
 integer :: NBin(1:NumHistograms),NHisto,i,MY_IDUP(1:9),ICOLUP(1:2,1:9),xBin(1:4)
 integer, pointer :: ijSel(:,:)
 integer :: NumPartonicChannels,iPartChannel,PartChannelAvg,flav1,flav2,ID_DK(6:9)
 real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,ML1,ML2,ML3,ML4,MZ1,MZ2
 real(8) :: MomExt(1:4,1:8),MomDK(1:4,1:4),xRnd
 logical :: applyPSCut
 include 'vegas_common.f'

    MomExt(:,:)=0d0
    MomDK(:,:)=0d0
    LO_Res_Unpol = 0d0
    EvalWeighted = 0d0
    m1ffwgt = 1d0 ! Multiplicative factor
    m2ffwgt = 1d0 ! Multiplicative factor

    if( OffShellReson ) then
      call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    else
      call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    endif
    EvalCounter = EvalCounter+1

   if( Process.eq.0 .or. PChannel.eq.0) then
      NumPartonicChannels = 1
      iPart_sel = 0
      jPart_sel = 0
   elseif( Process.eq.1 .or. PChannel.eq.1 ) then
      NumPartonicChannels = 10
      iPartChannel = int(yRnd(12) * (NumPartonicChannels)) +2 ! this runs from 2..11
      call getRef_PPXchannelHash(ijSel)
      iPart_sel = ijSel(iPartChannel,1)
      jPart_sel = ijSel(iPartChannel,2)
   elseif( Process.eq.2 .or. PChannel.eq.2 ) then
      NumPartonicChannels = 11
      iPartChannel = int(yRnd(12) * (NumPartonicChannels)) +1 ! this runs from 1..11
      call getRef_PPXchannelHash(ijSel)
      iPart_sel = ijSel(iPartChannel,1)
      jPart_sel = ijSel(iPartChannel,2)
   endif
   PartChannelAvg = NumPartonicChannels

   if( unweighted .and. .not.warmup .and.  sum(AccepCounter_part(:,:)) .eq. sum(RequEvents(:,:)) ) then
      stopvegas=.true.
   endif
   if( unweighted .and. .not. warmup .and. AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel)  ) return

!    particle associations:
!
!    IDUP(6)  -->  MomExt(:,6) == MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomExt(:,5) == MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomExt(:,8) == MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomExt(:,7) == MomDK(:,3)  -->  ubar-spinor
    ICOLUP(1:2,1:9) = 0
    call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9),FinalStateWeight)
    MY_IDUP(1:3) = 0
    if( RandomizeVVdecays .and. OffShellV1.eqv.OffShellV2) then
       call random_number(yrnd(17:18))
       if( (MY_IDUP(6).ne.MY_IDUP(8)) .and. (IsAZDecay(DecayMode1)) .and. (IsAZDecay(DecayMode2)) ) then
         if( (yrnd(18).le.0.5d0) ) then
           call swapi(MY_IDUP(4),MY_IDUP(5))
           call swapi(MY_IDUP(6),MY_IDUP(8))
           call swapi(MY_IDUP(7),MY_IDUP(9))
           call swapi(ICOLUP(1,6),ICOLUP(1,8))
           call swapi(ICOLUP(1,7),ICOLUP(1,9))
           call swapi(ICOLUP(2,6),ICOLUP(2,8))
           call swapi(ICOLUP(2,7),ICOLUP(2,9))
         endif
       elseif( (IsAWDecay(DecayMode1)) .and. (IsAWDecay(DecayMode2)) ) then
         if( (yrnd(18).le.0.5d0) ) then
           MY_IDUP(4) = ChargeFlip(MY_IDUP(4))
           MY_IDUP(5) = ChargeFlip(MY_IDUP(5))
           MY_IDUP(6) = ChargeFlip(MY_IDUP(6))
           MY_IDUP(7) = ChargeFlip(MY_IDUP(7))
           MY_IDUP(8) = ChargeFlip(MY_IDUP(8))
           MY_IDUP(9) = ChargeFlip(MY_IDUP(9))
           ! if there's a charge flip then the order of particle and anti-particles needs to be flipped, too
           call swapi(MY_IDUP(6),MY_IDUP(7))
           call swapi(MY_IDUP(8),MY_IDUP(9))
         endif
         if( (yrnd(17).le.0.5d0) ) then
           call swapi(MY_IDUP(4),MY_IDUP(5))
           call swapi(MY_IDUP(6),MY_IDUP(8))
           call swapi(MY_IDUP(7),MY_IDUP(9))
           call swapi(ICOLUP(1,6),ICOLUP(1,8))
           call swapi(ICOLUP(1,7),ICOLUP(1,9))
           call swapi(ICOLUP(2,6),ICOLUP(2,8))
           call swapi(ICOLUP(2,7),ICOLUP(2,9))
         endif
       endif ! Do not swap any state with no two identical VV decays! This would break mass determination in EvalPhasespace_H4f.
    endif
    ID_DK(6:9)=MY_IDUP(6:9) ! Decay ids to pass into the ME
    if(MY_IDUP(4).eq.Pho_ .and. ID_DK(6).eq.Not_a_particle_) then
      ID_DK(6)=MY_IDUP(4)
    endif
    if(MY_IDUP(5).eq.Pho_ .and. ID_DK(8).eq.Not_a_particle_) then
      ID_DK(8)=MY_IDUP(5)
    endif

    if(abs(EHat-M_Reso).ge.20.0d0*Ga_Reso) return ! for some reason this removes some of the cross section, but significantly improves speed for OffXVV=111 !!
    if( any(yRnd(4:5).gt.0.99d0) .or. EHat.lt.5d0*GeV ) return ! the cut at 0.99 is required for EvalPhasespace_H4f when interference is turned on. Otherwise, it becomes unstable.

    call EvalPhasespace_H4f(yRnd(3),yRnd(4:11),EHat,MomExt(1:4,1:8),ID_DK(6:9),PSWgt)

    ! Line above can handle Zgam and gamgam as well
    !if( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
    !   call EvalPhasespace_HVga(yRnd(3:4),EHat,MomExt(1:4,1:8),PSWgt)  !   ga-ga
    !elseif( .not. IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
    !   call EvalPhasespace_HVga(yRnd(3:7),EHat,MomExt(1:4,1:8),PSWgt)  !    Z-ga
    !else
    !   call EvalPhasespace_H4f(yRnd(3),yRnd(4:11),EHat,MomExt(1:4,1:8),PSWgt)  ! VV-->4l
    !endif

    if(.not.IsAPhoton(DecayMode1)) then
       ML1 = getMass(MY_IDUP(6))
       ML2 = getMass(MY_IDUP(7))
       MZ1 = get_MInv(MomExt(1:4,5)+MomExt(1:4,6))
       if( MZ1.lt.(ML1+ML2) ) then
          EvalWeighted = 0d0
          return
       endif
    endif
    if(.not.IsAPhoton(DecayMode2)) then
       ML3 = getMass(MY_IDUP(8))
       ML4 = getMass(MY_IDUP(9))
       MZ2 = get_MInv(MomExt(1:4,7)+MomExt(1:4,8))
       if( MZ2.lt.(ML3+ML4) ) then
          EvalWeighted = 0d0
          return
       endif
    endif

    call boost2Lab(eta1,eta2,8,MomExt(1:4,1:8))

    call Kinematics(4,MomExt,MomExt(1:4,5:8),applyPSCut,NBin)
    if( applyPSCut  .or. PSWgt.eq.0d0  ) then
       EvalWeighted = 0d0
       return
    endif

   call SetRunningScales( (/ (MomExt(1:4,3)+MomExt(1:4,4)),Mom_Not_a_particle(1:4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) )
   call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)

   FluxFac = 1d0/(2d0*EHat**2)
   PDFFac = pdf(LHA2M_pdf(iPart_sel),1)  *  pdf(LHA2M_pdf(jPart_sel),2)
   PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * PDFFac * VgsWgt * PartChannelAvg * FinalStateWeight

   !print *,"Mom1: ",MomExt(1:4,1)
   !print *,"Mom2: ",MomExt(1:4,2)
   !print *,"Mom3: ",MomExt(1:4,3)
   !print *,"Mom4: ",MomExt(1:4,4)
   !print *,"Mom5: ",MomExt(1:4,5)
   !print *,"Mom6: ",MomExt(1:4,6)
   !print *,"Mom7: ",MomExt(1:4,7)
   !print *,"Mom8: ",MomExt(1:4,8)
   !print *,"iPart, jPart: ",iPart_sel,jPart_sel

   if ( PChannel.eq.0 .or.(PChannel.eq.2 .and. iPart_sel.eq.0 .and. jPart_sel.eq.0)) then
      if (Process.eq.0) then
          call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8)/),ID_DK(6:9),LO_Res_Unpol)
      elseif(Process.eq.2) then
          call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8)/),ID_DK(6:9),LO_Res_Unpol)
      endif
      MY_IDUP(1:2)=(/Glu_,Glu_/)
      ICOLUP(1:2,1) = (/501,502/)
      ICOLUP(1:2,2) = (/502,501/)

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
      !print *,"LO_Res_Unpol gg:",LO_Res_Unpol
   elseif (PChannel.eq.1.or.PChannel.eq.2) then
      if( Process.eq.1 .and. iPart_sel.gt.0d0 ) then
          call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8)/),ID_DK(6:9),LO_Res_Unpol)
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/000,501/)
      elseif( Process.eq.1 .and. iPart_sel.lt.0d0 ) then
          call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8)/),ID_DK(6:9),LO_Res_Unpol)
          ICOLUP(1:2,1) = (/000,502/)
          ICOLUP(1:2,2) = (/502,000/)
      elseif( Process.eq.2 .and. iPart_sel.gt.0d0 ) then
          call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8)/),ID_DK(6:9),LO_Res_Unpol)
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/000,501/)
      elseif( Process.eq.2 .and. iPart_sel.lt.0d0 ) then
          call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,8)/),ID_DK(6:9),LO_Res_Unpol)
          ICOLUP(1:2,1) = (/000,502/)
          ICOLUP(1:2,2) = (/502,000/)
      endif
      MY_IDUP(1:2)=(/LHA2M_pdf(iPart_sel),LHA2M_pdf(jPart_sel)/)

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2
      !print *,"LO_Res_Unpol qq:",LO_Res_Unpol
   endif

   if(.not.IsAPhoton(DecayMode1)) then
      call ShiftMass(MomExt(1:4,5),MomExt(1:4,6), GetMass(MY_IDUP(6)),GetMass(MY_IDUP(7)),  MomDK(1:4,1),MomDK(1:4,2), MassWeight=m1ffwgt)
      PreFac = PreFac*m1ffwgt
   else
      MomDK(:,1)=MomExt(:,5)
      MomDK(:,2)=MomExt(:,6)
   endif
   if(.not.IsAPhoton(DecayMode2)) then
      call ShiftMass(MomExt(1:4,7),MomExt(1:4,8), GetMass(MY_IDUP(8)),GetMass(MY_IDUP(9)),  MomDK(1:4,3),MomDK(1:4,4), MassWeight=m2ffwgt)
      PreFac = PreFac*m2ffwgt
   else
      MomDK(:,3)=MomExt(:,7)
      MomDK(:,4)=MomExt(:,8)
   endif


   EvalWeighted = LO_Res_Unpol * PreFac
   if( WidthScheme.ne.2 .and. Process.eq.0 ) EvalWeighted = EvalWeighted * ReweightBWPropagator( Get_MInv2( MomExt(1:4,3)+MomExt(1:4,4) ) )


   if( unweighted ) then
     if( warmup ) then

       CrossSec(iPart_sel,jPart_sel) = CrossSec(iPart_sel,jPart_sel) + EvalWeighted
       CrossSecMax(iPart_sel,jPart_sel) = max(CrossSecMax(iPart_sel,jPart_sel),EvalWeighted)

     else! not warmup

       call random_number(xRnd)
       if( EvalWeighted.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",EvalWeighted, CrossSecMax(iPart_sel,jPart_sel)
         write(io_stdout, "(2X,A,1PE13.6,1PE13.6,1PE13.6,I3,I3)") "CrossSecMax is too small.",EvalWeighted, CrossSecMax(iPart_sel,jPart_sel),EvalWeighted/CrossSecMax(iPart_sel,jPart_sel),iPart_sel,jPart_sel
         AlertCounter = AlertCounter + 1

         CrossSecMax(iPart_sel,jPart_sel)=CrossSecMax(iPart_sel,jPart_sel) * 1.1d0

       elseif( EvalWeighted .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then

         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPart_sel,jPart_sel) = AccepCounter_part(iPart_sel,jPart_sel) + 1
         do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo

         call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))

         if( abs(MY_IDUP(7)).eq.abs(ElM_) ) then
             flav1 = 1
         elseif( abs(MY_IDUP(7)).eq.abs(MuM_) ) then
             flav1 = 2
         elseif( abs(MY_IDUP(7)).eq.abs(TaM_) ) then
             flav1 = 3
         elseif( IsANeutrino(MY_IDUP(7)) ) then
             flav1 = 4
         elseif( IsAQuark(MY_IDUP(7)) ) then
             flav1 = 5
         else! for photons
             flav1 = 1
         endif
         if( abs(MY_IDUP(9)).eq.abs(ElM_) ) then
             flav2 = 1
         elseif( abs(MY_IDUP(9)).eq.abs(MuM_) ) then
             flav2 = 2
         elseif( abs(MY_IDUP(9)).eq.abs(TaM_) ) then
             flav2 = 3
         elseif( IsANeutrino(MY_IDUP(9)) ) then
             flav2 = 4
         elseif( IsAQuark(MY_IDUP(9)) ) then
             flav2 = 5
         else! for photons
             flav2 = 1
         endif
         Br_counter(flav1,flav2) = Br_counter(flav1,flav2) + 1

       endif

     endif! warmup

   else! weighted


      if( EvalWeighted.ne.0d0 ) then
        AccepCounter=AccepCounter+1
        if( writeWeightedLHE .and. (.not. warmup) ) then
            call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=EvalWeighted)
        endif
        do NHisto=1,NumHistograms-3
            call intoHisto(NHisto,NBin(NHisto),EvalWeighted)
        enddo
      endif
   endif! unweighted



   EvalWeighted = EvalWeighted/VgsWgt
   RETURN


END FUNCTION

FUNCTION EvalUnWeighted_DecayToVV(yRnd,genEvt,EHat,Res,AcceptedEvent,MY_IDUP,ICOLUP)
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
real(8) :: Res!  .ne.0: accepted event,  .eq.0: reject event,   .eq.-1: reject event and exit the loop over 'tries'
real(8) :: EvalUnWeighted_DecayToVV,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
real(8) :: tau,x1,x2,sHatJacobi,PreFac,FinalStateWeight
integer :: NBin(1:NumHistograms),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomExt_f(1:4,1:4),MomDK_f(1:4,1:4),MomDK_massless(1:4,1:4)
logical :: applyPSCut,genEvt
real(8) :: CS_max,eta1,eta2
real(8) :: oneovervolume, bound(1:11), sumtot,yz1,yz2,EZ_max,dr,MZ1,MZ2,ML1,ML2,ML3,ML4
integer :: i1, i2, MY_IDUP(1:9), ICOLUP(1:2,1:9),OSPair,OSSFPair,LeptInEvent_tmp(0:8),JetsInEvent_tmp(0:8),ordered_Lept(1:8), ID_DK(6:9)
real(8)::  ntRnd,ZMass(1:2),AcceptedEvent(1:4,1:4)
real(8) :: offzchannel
include 'vegas_common.f'
include 'csmaxvalue.f'


   oneovervolume = one
   ICOLUP(1:2,1:9) = 0
   EvalUnWeighted_DecayToVV = 0d0
   Res = 0d0
   EvalCounter = EvalCounter+1



   MY_IDUP(3)= Hig_
!    particle associations:
!
!    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor
   call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9),FinalStateWeight)
!    print *, MY_IDUP(4:9);pause

   if( (RandomizeVVdecays.eqv..true.) ) then
      if( (MY_IDUP(6).ne.MY_IDUP(8)) .and. (IsAZDecay(DecayMode1)) .and. (IsAZDecay(DecayMode2)) ) then
        if( (yrnd(16).le.0.5d0) ) then
         call swapi(MY_IDUP(4),MY_IDUP(5))
         call swapi(MY_IDUP(6),MY_IDUP(8))
         call swapi(MY_IDUP(7),MY_IDUP(9))
         call swapi(ICOLUP(1,6),ICOLUP(1,8))
         call swapi(ICOLUP(1,7),ICOLUP(1,9))
         call swapi(ICOLUP(2,6),ICOLUP(2,8))
         call swapi(ICOLUP(2,7),ICOLUP(2,9))
        endif
     elseif( (IsAWDecay(DecayMode1)) .and. (IsAWDecay(DecayMode2)) ) then
        if( (yrnd(16).le.0.5d0) ) then
         MY_IDUP(4) = ChargeFlip(MY_IDUP(4))
         MY_IDUP(5) = ChargeFlip(MY_IDUP(5))
         MY_IDUP(6) = ChargeFlip(MY_IDUP(6))
         MY_IDUP(7) = ChargeFlip(MY_IDUP(7))
         MY_IDUP(8) = ChargeFlip(MY_IDUP(8))
         MY_IDUP(9) = ChargeFlip(MY_IDUP(9))
         ! if there's a charge flip then the order of particle and anti-particles needs to be flipped, too
         call swapi(MY_IDUP(6),MY_IDUP(7))
         call swapi(MY_IDUP(8),MY_IDUP(9))
        endif
        if( (yrnd(17).le.0.5d0) ) then
         call swapi(MY_IDUP(4),MY_IDUP(5))
         call swapi(MY_IDUP(6),MY_IDUP(8))
         call swapi(MY_IDUP(7),MY_IDUP(9))
         call swapi(ICOLUP(1,6),ICOLUP(1,8))
         call swapi(ICOLUP(1,7),ICOLUP(1,9))
         call swapi(ICOLUP(2,6),ICOLUP(2,8))
         call swapi(ICOLUP(2,7),ICOLUP(2,9))
        endif
     endif
  endif
  ID_DK(6:9)=MY_IDUP(6:9) ! Decay ids to pass into the ME
  if(MY_IDUP(4).eq.Pho_ .and. ID_DK(6).eq.Not_a_particle_) then
    ID_DK(6)=MY_IDUP(4)
  endif
  if(MY_IDUP(5).eq.Pho_ .and. ID_DK(8).eq.Not_a_particle_) then
    ID_DK(8)=MY_IDUP(5)
  endif

  eta1=1d0; eta2=1d0
  sHatJacobi = 1d0

  yz1 = yRnd(10)
  yz2 = yRnd(11)
  offzchannel = yRnd(15) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * ( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )

            EZ_max = EHat - MZ1*0.99d0
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )

        else   !M_Reso.le.2d0*M_V_ps
            if (offzchannel.le.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
                MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
                MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
                sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/(  &
                1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )     &
                + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
            elseif(offzchannel.gt.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
                MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
                MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
                sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/( &
                1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )    &
                + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
            endif
       endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
        MZ1 = getMass(MY_IDUP(4))
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat - MZ1*0.99
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
        else
            MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
        endif

  elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
        MZ2 = getMass(MY_IDUP(5))
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat - MZ2*0.99
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
         else
            MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
        endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
        MZ1 = getMass(MY_IDUP(4))
        MZ2 = getMass(MY_IDUP(5))
  endif


    if( MZ1+MZ2.gt.EHat ) then
      EvalUnWeighted_DecayToVV = 0d0
      RejeCounter = RejeCounter + 1
      return
    endif




   call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   if( .not.IsAPhoton(DecayMode1) .and. .not.IsAPhoton(DecayMode2) ) then ! don't decay the photon
      ML1 = getMass(MY_IDUP(6))
      ML2 = getMass(MY_IDUP(7))
      ML3 = getMass(MY_IDUP(8))
      ML4 = getMass(MY_IDUP(9))
      if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
          EvalUnWeighted_DecayToVV = 0d0
          RejeCounter = RejeCounter + 1
          return
      endif
      call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
      call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,ML3,ML4,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
      PSWgt = PSWgt * PSWgt2*PSWgt3

      if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then! introduce this momentum flip to allow proper mapping of integrand with Z-poles at MZ2=(p2+p3)^2 and MZ2=(p1+p4)^2
          if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK(1:4,1),MomDK(1:4,3) )
!           PSWgt = PSWgt * 2d0
      endif
    elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        ML1=0d0; ML2=0d0; ML3=0d0; ML4=0d0
        MomDK(1:4,1) = MomExt(1:4,3)
        MomDK(1:4,2) = 0d0
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
    elseif( .not.IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        ML1 = getMass(MY_IDUP(6))
        ML2 = getMass(MY_IDUP(7))
        ML3=0d0; ML4=0d0
        if( (MZ1.lt.ML1+ML2) ) then
            EvalUnWeighted_DecayToVV = 0d0
            return
        endif
        call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
        PSWgt = PSWgt * PSWgt2
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
   endif


    if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode1)) ) then
        call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    else
        call Error("This should not happen anymore")
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalUnWeighted_DecayToVV = 0d0
      return
    endif


    call SetRunningScales( (/ (MomExt(1:4,3)+MomExt(1:4,4)),Mom_Not_a_particle(1:4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) ) ! Call anyway

IF( GENEVT ) THEN

      MY_IDUP(1:2)=(/Glu_,Glu_/)
      ICOLUP(1:2,1) = (/501,502/)
      ICOLUP(1:2,2) = (/502,501/)

      if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode1)) ) then
           AcceptedEvent(1:4,1:4) = MomDK(1:4,1:4)
      else
           AcceptedEvent(1:4,1:4) = MomDK_f(1:4,1:4)
      endif

      if (Process.eq.0) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            else
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            endif
      endif


      PreFac = 2d0 * hbarc2XsecUnit * sHatJacobi * PSWgt * SymmFac * FinalStateWeight
      EvalUnWeighted_DecayToVV = LO_Res_Unpol * PreFac

      CS_max = csmax(0,0)

      if( EvalUnWeighted_DecayToVV .gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_DecayToVV, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_DecayToVV, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted_DecayToVV .gt. yRnd(14)*CS_max ) then

         if( any(RequestNLeptons.ge.0) .or. any(RequestNJets.ge.0) ) then! lepton or jet filter
                LeptInEvent_tmp(0:8) = LeptInEvent(0:8)
                JetsInEvent_tmp(0:8) = JetsInEvent(0:8)
    !             print *, ""
                do i1=6,9
                    if( IsALepton(MY_IDUP(i1)) ) then
                      LeptInEvent_tmp(0) = LeptInEvent_tmp(0)+1
                      LeptInEvent_tmp( LeptInEvent_tmp(0) ) = ConvertLHE(MY_IDUP(i1))
                    endif
                    if( IsAJet(MY_IDUP(i1)) ) then
                      JetsInEvent_tmp(0) = JetsInEvent_tmp(0)+1
                      JetsInEvent_tmp( JetsInEvent_tmp(0) ) = ConvertLHE(MY_IDUP(i1))
                    endif
                enddo
! print *, "leptons in event: ",LeptInEvent_tmp(1: LeptInEvent_tmp(0))
                ordered_Lept(1:8) = (/1,2,3,4,5,6,7,8/)! order leptons for tau decay associations
                call BubleSort(LeptInEvent_tmp(0),dabs(dble(LeptInEvent_tmp(1:LeptInEvent_tmp(0)))), ordered_Lept(1:LeptInEvent_tmp(0)))
! do i1=1,LeptInEvent_tmp(0)
!   print *, "new order:",LeptInEvent_tmp( ordered_Lept(i1) )
! enddo
! pause
                if( (RequestNLeptons(1).ge.0 .and. LeptInEvent_tmp(0) .lt. RequestNLeptons(1)) .or. (RequestNLeptons(2).ge.0 .and. LeptInEvent_tmp(0) .gt. RequestNLeptons(2)) ) then
    !                 print *,"not enough leptons, reject!" !,LeptInEvent_tmp(1: LeptInEvent_tmp(0))
                    Res = -1d0
                    return
                endif
                if( (RequestNJets(1).ge.0 .and. JetsInEvent_tmp(0) .lt. RequestNJets(1)) .or. (RequestNJets(2).ge.0 .and. JetsInEvent_tmp(0) .gt. RequestNJets(2)) ) then
    !                 print *,"not enough jets, reject!" !,JetsInEvent_tmp(1: JetsInEvent_tmp(0))
                    Res = -1d0
                    return
                endif
                if( any(RequestOS.ge.0) .or. any(RequestOSSF.ge.0) ) then
                    OSPair = 0
                    OSSFPair = 0
!                     do i1=1,LeptInEvent_tmp(0)-1
!                         do i2=i1+1,LeptInEvent_tmp(0)
                    do i1=LeptInEvent_tmp(0),2,-1
                        do i2=i1-1,1,-1
                            if(      ( LeptInEvent_tmp(i1)+LeptInEvent_tmp(i2).eq.0                                                     )    &     ! found a l+ l- pair
                                .OR. ( CountTauAsAny .AND. LeptInEvent_tmp(i1).ne.-999 .AND. LeptInEvent_tmp(i2).ne.-999 .AND. ( &
                                            ( abs(LeptInEvent_tmp(i1)).eq.ConvertLHE(TaM_) .and. LeptInEvent_tmp(i1)*LeptInEvent_tmp(i2).lt.0  )    &     ! found l tau pair
                                       .OR. ( abs(LeptInEvent_tmp(i2)).eq.ConvertLHE(TaM_) .and. LeptInEvent_tmp(i1)*LeptInEvent_tmp(i2).lt.0  )    &     ! found l tau pair
                                     )                                                                                         ) &
                            ) then
                              LeptInEvent_tmp(i1) = -999! remove from list
                              LeptInEvent_tmp(i2) = -999! remove from list
                              OSPair = OSPair + 1
                              OSSFPair = OSSFPair + 1
                              exit
                            endif
                        enddo
                    enddo
                    do i1=LeptInEvent_tmp(0),2,-1
                        do i2=i1-1,1,-1
                            if(      ( LeptInEvent_tmp(i1)*LeptInEvent_tmp(i2).lt.0 )    &     ! found a l+ l'- pair
                               .and. ( (LeptInEvent_tmp(i1).ne.-999) .and. (LeptInEvent_tmp(i2).ne.-999) ) &
                            ) then
                              LeptInEvent_tmp(i1) = -999! remove from list
                              LeptInEvent_tmp(i2) = -999! remove from list
                              OSPair = OSPair + 1
                              exit
                            endif
                        enddo
                    enddo
!                     print *, "found ",OSPair," OS pairs"
!                     print *, "found ",OSSFPair," OSSF pairs"
                    if( (RequestOS(1).ge.0 .and. OSPair.lt.RequestOS(1)) .or. (RequestOS(2).ge.0 .and. OSPair.gt.RequestOS(2)) ) then
!                         print *,"no OS pair, reject!" !,LeptInEvent_tmp(1: LeptInEvent_tmp(0))
                        Res = -1d0
                        return
                    endif
                    if( (RequestOSSF(1).ge.0 .and. OSSFPair.lt.RequestOSSF(1)) .or. (RequestOSSF(2).ge.0 .and. OSSFPair.gt.RequestOSSF(2)) ) then
!                         print *,"no OSSF pair, reject!" !,LeptInEvent_tmp(1: LeptInEvent_tmp(0))
                        Res = -1d0
                        return
                    endif
                endif
!                 print *, "accept event"
!                 pause
         endif! lepton filter

         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1
         Res = 1d0

      else
          RejeCounter = RejeCounter + 1
          Res = 0d0
      endif


ELSE! NOT GENEVT


      if (Process.eq.0) then
         MomExt(1:4,1) = (/EHat,0d0,0d0,0d0/)
         MomExt(1:4,2) = (/0d0,0d0,0d0,0d0/)

         if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
         else
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
         endif
      endif

     PreFac = 2d0 * hbarc2XsecUnit * sHatJacobi * PSWgt * SymmFac * FinalStateWeight
     EvalUnWeighted_DecayToVV = LO_Res_Unpol * PreFac
     Res = EvalUnWeighted_DecayToVV


      if (EvalUnWeighted_DecayToVV.gt.csmax(0,0)) then
          csmax(0,0) = EvalUnWeighted_DecayToVV
      endif


ENDIF! genEvt


RETURN
END FUNCTION


FUNCTION EvalWeighted_gg4f_fullproddec(yRnd,VgsWgt)
#if linkMELA==1
use ModKinematics
use ModParameters
use ModMisc
use ModMCFMWrapper
!use ModHiggs
#if compiler==1
use ifport
#endif
implicit none
integer,parameter :: inTop=1, inBot=2, V1=3, V2=4, Lep1P=5, Lep1M=6, Lep2P=7, Lep2M=8, NUP=8
real(8) :: yRnd(1:10),VgsWgt, EvalWeighted_gg4f_fullproddec
real(8) :: pdf(-6:6,1:2),me2(-5:5,-5:5)
real(8) :: eta1, eta2, FluxFac, Ehat
real(8) :: MomExt(1:4,1:NUP),MomShifted(1:4,1:NUP),PSWgt,FinalStateWeight,m1ffwgt,m2ffwgt
real(8) :: p_MCFM(mxpart,1:4),msq_MCFM(-5:5,-5:5),msq_VgsWgt(-5:5,-5:5),Wgt_Ratio_Interf,originalprobability
integer :: id_MCFM(mxpart),MY_IDUP(1:NUP),ICOLUP(1:2,1:NUP),NBin(1:NumHistograms),NHisto,ipart,jpart
integer :: i,j,k
real(8) :: PreFac,VegasWeighted_fullproddec,xRnd,LeptonAndVegasWeighted_fullproddec!,LO_Res_Unpol
logical :: applyPSCut,swap34_56
include 'vegas_common.f'
include 'maxwt.f'

   EvalWeighted_gg4f_fullproddec = 0d0
   Wgt_Ratio_Interf = 1d0
   m1ffwgt = 1d0
   m2ffwgt = 1d0
   iPart_sel = 0
   jPart_sel = 0

   !   throwing random number for accept-reject
   if(.not. warmup) then
       call random_number(xRnd)
   endif

   if( unweighted .and. .not. warmup .and. sum(AccepCounter_part(:,:)) .eq. sum(RequEvents(:,:)) ) then
      stopvegas=.true.
   endif
   if( unweighted .and. .not. warmup .and. AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel)  ) return


   call VVBranchings(MY_IDUP(V1:NUP),ICOLUP(1:2,V2+1:NUP),FinalStateWeight,700)
   call swap(MY_IDUP(Lep1P),MY_IDUP(Lep1M))
   call swap(MY_IDUP(Lep2P),MY_IDUP(Lep2M))
   MY_IDUP(1:2) = Glu_
   id_MCFM(1:2) = MY_IDUP(1:2)
   id_MCFM(3:6) = MY_IDUP(5:8)

   if (IsNaN(VgsWgt)) then
      write(6,*) "VegasWgt is NaN!"
      !pause
   endif

   call EvalPhasespace_gg4f(yRnd(1:10),eta1,eta2,EHat,MomExt(1:4,1:NUP),PSWgt,id_MCFM(1:6),swap34_56)
   PSWgt = PSWgt * FinalStateWeight
   if( PSWgt.lt.1d-33 ) then
      return
   endif
   call boost2Lab(eta1,eta2,NUP,MomExt(1:4,1:NUP))

   call Kinematics_gg4f_fullproddec(MomExt,id_MCFM,applyPSCut,NBin)
   if( applyPSCut ) then
      return
   endif

   call SetRunningScales( &
      (/ (MomExt(1:4,V1)+MomExt(1:4,V2)), Mom_Not_a_particle(1:4), Mom_Not_a_particle(1:4) /), &
      (/ Not_a_particle_, Not_a_particle_, Not_a_particle_, Not_a_particle_ /) &
      )
   call EvalAlphaS()
   !write(6,*) "setPDFs args:",eta1,eta2,alphas,alphas_mz
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)
   !pause

   ! GeV conversion is now done inside EvalAmp_qqVVqq
   call convert_to_MCFM(-MomExt(1:4,inTop), p_MCFM(1,1:4))
   call convert_to_MCFM(-MomExt(1:4,inBot), p_MCFM(2,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep1P), p_MCFM(3,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep1M), p_MCFM(4,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2P), p_MCFM(5,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2M), p_MCFM(6,1:4))
   msq_MCFM(:,:) = 0d0

   MomShifted = MomExt
   if (swap34_56) then
      if(.not.IsAPhoton(DecayMode1)) then
         call ShiftMass(MomExt(1:4,Lep1P),MomExt(1:4,Lep2M), GetMass(MY_IDUP(Lep1P)),GetMass(MY_IDUP(Lep2M)),MomShifted(1:4,Lep1P),MomShifted(1:4,Lep2M),MassWeight=m1ffwgt)
      endif
      if(.not.IsAPhoton(DecayMode2)) then
         call ShiftMass(MomExt(1:4,Lep2P),MomExt(1:4,Lep1M), GetMass(MY_IDUP(Lep2P)),GetMass(MY_IDUP(Lep1M)),MomShifted(1:4,Lep2P),MomShifted(1:4,Lep1M),MassWeight=m2ffwgt)
      endif
   else
      if(.not.IsAPhoton(DecayMode1)) then
         call ShiftMass(MomExt(1:4,Lep1P),MomExt(1:4,Lep1M), GetMass(MY_IDUP(Lep1P)),GetMass(MY_IDUP(Lep1M)),MomShifted(1:4,Lep1P),MomShifted(1:4,Lep1M),MassWeight=m1ffwgt)
      endif
      if(.not.IsAPhoton(DecayMode2)) then
         call ShiftMass(MomExt(1:4,Lep2P),MomExt(1:4,Lep2M), GetMass(MY_IDUP(Lep2P)),GetMass(MY_IDUP(Lep2M)),MomShifted(1:4,Lep2P),MomShifted(1:4,Lep2M),MassWeight=m2ffwgt)
      endif
   endif
   PSWgt = PSWgt * m1ffwgt * m2ffwgt


   call EvalAmp_gg4f(id_MCFM, p_MCFM, msq_MCFM)
   !msq_MCFM(iPart_sel,jPart_sel) = 1d0
   !call EvalAmp_gg_H_VV( &
   !   (/-MomExt(1:4,inTop),-MomExt(1:4,inBot),MomExt(1:4,Lep1P),MomExt(1:4,Lep1M),MomExt(1:4,Lep2P),MomExt(1:4,Lep2M)/), &
   !   id_MCFM(3:6), &
   !   LO_Res_Unpol &
   !   )
   !do i=1,6
   !   print *,"id(",i,"):",convertLHE(id_MCFM(i))
   !enddo
   !do i=1,2
   !   print *,"Mom",i,"=",MomExt(1:4,i)/GeV
   !enddo
   !do i=Lep1P,Lep2M
   !   print *,"Mom",i,"=",MomExt(1:4,i)/GeV
   !enddo
   !print *,"MCFM ME = ",msq_MCFM(0,0)
   !!print *,"JHU ME = ",LO_Res_Unpol
   !write(6,*) "Mu_Fact =",Mu_Fact
   !write(6,*) "Mu_Ren =",Mu_Ren
   !write(6,*) "alphas =",alphas
   !write(6,*) "alphas_mz =",alphas_mz
   !write(6,*) "pdf1 =",pdf(LHA2M_pdf(iPart_sel),1)
   !write(6,*) "pdf2 =",pdf(LHA2M_pdf(jPart_sel),2)
   !pause


   originalprobability = msq_MCFM(iPart_sel,jPart_sel)

   PreFac = hbarc2XsecUnit / (GeV**4)  ! adjust msq_MCFM for GeV units of MCFM mat.el.
   msq_MCFM = msq_MCFM * PreFac

   if ( &
      msq_MCFM(iPart_sel,jPart_sel) .le. 0d0 .or. &
      pdf(LHA2M_pdf(iPart_sel),1) .le. 0d0 .or. &
      pdf(LHA2M_pdf(jPart_sel),2) .le. 0d0 .or. &
      IsNaN(msq_MCFM(iPart_sel,jPart_sel)) .or. &
      IsNaN(pdf(LHA2M_pdf(iPart_sel),1)) .or. &
      IsNaN(pdf(LHA2M_pdf(jPart_sel),2)) &
      ) then
      write(6,*) "Mu_Fact =",Mu_Fact
      write(6,*) "Mu_Ren =",Mu_Ren
      write(6,*) "alphas =",alphas
      write(6,*) "alphas_mz =",alphas_mz
      write(6,*) "msq_MCFM(",iPart_sel,",",jPart_sel,") =",msq_MCFM(iPart_sel,jPart_sel)
      write(6,*) "pdf1 =",pdf(LHA2M_pdf(iPart_sel),1)
      write(6,*) "pdf2 =",pdf(LHA2M_pdf(jPart_sel),2)
      do jpart=1,6
         write(6,*) "P_MCFM(",convertLHE(id_MCFM(jpart)),")=",p_MCFM(jpart,:)
      enddo
      if ( &
            IsNaN(msq_MCFM(iPart_sel,jPart_sel)) .or. &
            IsNaN(pdf(LHA2M_pdf(iPart_sel),1)) .or. &
            IsNaN(pdf(LHA2M_pdf(jPart_sel),2)) &
         ) then
         pause
      endif
      return
    endif

   !EvalWeighted_gg4f_fullproddec = msq_MCFM(iPart_sel,jPart_sel) * pdf(LHA2M_pdf(iPart_sel),1)*pdf(LHA2M_pdf(jPart_sel),2)
   EvalWeighted_gg4f_fullproddec = PSWgt * FluxFac * msq_MCFM(iPart_sel,jPart_sel) * pdf(LHA2M_pdf(iPart_sel),1)*pdf(LHA2M_pdf(jPart_sel),2)
   VegasWeighted_fullproddec = EvalWeighted_gg4f_fullproddec * VgsWgt
   !if (EvalWeighted_gg4f_fullproddec.eq.0d0) then
   !   write(6,*) "EvalWeighted_gg4f_fullproddec==0. Ids:",id_MCFM
   !endif
   !write(6,*) "originalprobability,EvalWeighted_gg4f_fullproddec,VgsWgt=",originalprobability,EvalWeighted_gg4f_fullproddec,VgsWgt
   !pause


   if( unweighted ) then

     if( warmup ) then

!        if( VegasWeighted_fullproddec.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
!            print *, "New max",iPart_sel,jPart_sel,VegasWeighted_fullproddec
!        endif

       CrossSec(iPart_sel,jPart_sel) = CrossSec(iPart_sel,jPart_sel) + VegasWeighted_fullproddec
       CrossSecMax(iPart_sel,jPart_sel) = max(CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_fullproddec)

       do NHisto=1,NumHistograms
         call intoHisto(NHisto,NBin(NHisto),VegasWeighted_fullproddec)
       enddo

       !if (FindCrossSectionWithWeights) then
       !  LeptonAndVegasWeighted_fullproddec = VegasWeighted_fullproddec * ReweightLeptonInterference(id_MCFM, p_MCFM, originalprobability)
       !  CrossSectionWithWeights = CrossSectionWithWeights + LeptonAndVegasWeighted_fullproddec
       !  CrossSectionWithWeightsErrorSquared = CrossSectionWithWeightsErrorSquared + LeptonAndVegasWeighted_fullproddec**2
       !endif

     else! not warmup

       EvalCounter = EvalCounter+1

       if( VegasWeighted_fullproddec.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",VegasWeighted_fullproddec, CrossSecMax(iPart_sel,jPart_sel)
          write(io_stdout, "(2X,A,1PE13.6,1PE13.6,1PE13.6,I4,I4)") "CrossSecMax is too small.",VegasWeighted_fullproddec, CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_fullproddec/CrossSecMax(iPart_sel,jPart_sel),iPart_sel,jPart_sel
          AlertCounter = AlertCounter + 1

!          This dynamically increases the maximum in case it is exceeded
          CrossSecMax(iPart_sel,jPart_sel) = VegasWeighted_fullproddec
          write(io_LogFile,"(2X,A,1PE13.6)") "Increasing CrossSecMax to ",VegasWeighted_fullproddec
          write(io_stdout, "(2X,A,1PE13.6)") "Increasing CrossSecMax to ",VegasWeighted_fullproddec

       elseif( VegasWeighted_fullproddec .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then
          AccepCounter = AccepCounter + 1
          AccepCounter_part(iPart_sel,jPart_sel) = AccepCounter_part(iPart_sel,jPart_sel) + 1

          !Wgt_Ratio_Interf = ReweightLeptonInterference(id_MCFM, p_MCFM, originalprobability)

          call WriteOutEvent_gg4f_fullproddec(MomShifted,MY_IDUP,ICOLUP,EventWeight=Wgt_Ratio_Interf)

          do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),1d0)
          enddo
       else
          RejeCounter=RejeCounter+1
       endif

     endif! warmup

   else! weighted

      if( VegasWeighted_fullproddec.ne.0d0 ) then
        AccepCounter=AccepCounter+1
        if( writeWeightedLHE .and. (.not. warmup) ) then
            call WriteOutEvent_gg4f_fullproddec(MomShifted,MY_IDUP,ICOLUP,EventWeight=VegasWeighted_fullproddec)
        endif
        do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),VegasWeighted_fullproddec)
        enddo

      endif

   endif! unweighted

#else
implicit none
real(8) :: yRnd(1:10),VgsWgt, EvalWeighted_gg4f_fullproddec
   EvalWeighted_gg4f_fullproddec = 0d0
   print *, "To use this process, please set linkMELA=Yes in the makefile and recompile."
   print *, "You will also need to have a compiled JHUGenMELA in the directory specified by JHUGenMELADir in the makefile."
   stop 1
#endif
RETURN
END FUNCTION


!  FUNCTION EvalWeighted(yRnd,VgsWgt)    ! this is a function which is only for computations
!  use ModKinematics                     ! with weighted events
!  use ModParameters
!  use ModGraviton
!  use ModHiggs
!  use ModZprime
!  use ModMisc
! #if compiler==1
!  use ifport
! #endif
!  implicit none
!  real(8) :: EvalWeighted,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
!  real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac1,PDFFac2
!  real(8) :: pdf(-6:6,1:2)
!  integer :: NBin(1:NumHistograms),NHisto,i,MY_IDUP(1:9), ICOLUP(1:2,1:9),xBin(1:4)
!  real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,FinalStateWeight
!  real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomDK_massless(1:4,1:4)
!  real(8) :: EZ,  EZ1, EZ2, pz, xmax, xx1, xx2
!  real(8) :: pz12, MomExt_f(1:4,1:4), MomDK_f(1:4,1:4)
!  real(8) :: MZ1,MZ2,ML1,ML2,ML3,ML4,EZ_max,dr, MG, yz1,yz2
!  real(8) :: offzchannel
!  logical :: applyPSCut
!  include 'csmaxvalue.f'
!
!     EvalWeighted = 0d0
!     if( OffShellReson ) then
!       call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!     else
! !       call PDFMapping(11,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!       call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!     endif
!     EvalCounter = EvalCounter+1
!
!
! !    particle associations:
! !
! !    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
! !    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
! !    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
! !    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor
!     ICOLUP(1:2,1:9) = 0
!     call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9),FinalStateWeight)
!     MY_IDUP(1:3) = 0
!     yz1 = yRnd(10)
!     yz2 = yRnd(11)
!     offzchannel = yRnd(12) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
!   if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
!         if(M_Reso.gt.2d0*M_V_ps) then
!             EZ_max = EHat
!             dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
!             MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
!             sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * ( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
!
!             EZ_max = EHat - MZ1*0.99d0
!             dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
!             MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
!             sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
!
!         elseif(M_Reso.lt.2d0*M_V_ps) then
!             if (offzchannel.le.0.5d0) then
!                 EZ_max = EHat
!                 dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
!                 MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
!                 MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
!                 sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/(  &
!                 1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )     &
!                 + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
!                 sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
!             elseif(offzchannel.gt.0.5d0) then
!                 EZ_max = EHat
!                 dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
!                 MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
!                 MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
!                 sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/( &
!                 1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )    &
!                 + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
!                 sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
!             endif
!        endif
!
!   elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
!         MZ1 = getMass(MY_IDUP(4))
!         if(M_Reso.gt.2d0*M_V_ps) then
!             EZ_max = EHat - MZ1*0.99
!             dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
!             MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
!             sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
!         else
!             MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
!             sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
!         endif
!
!   elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
!         MZ2 = getMass(MY_IDUP(5))
!         if(M_Reso.gt.2d0*M_V_ps) then
!             EZ_max = EHat - MZ2*0.99
!             dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
!             MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
!             sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
!          else
!             MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
!             sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
!         endif
!
!   elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
!         MZ1 = getMass(MY_IDUP(4))
!         MZ2 = getMass(MY_IDUP(5))
!   endif
!
!     if( EHat.lt.MZ1+MZ2 ) then
!       EvalWeighted = 0d0
!       return
!     endif
!
!     call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
!     call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
!     if( .not.IsAPhoton(DecayMode1) .and. .not.IsAPhoton(DecayMode2) ) then ! decay ZZ's and WW's
!         ML1 = getMass(MY_IDUP(7))
!         ML2 = getMass(MY_IDUP(6))
!         ML3 = getMass(MY_IDUP(9))
!         ML4 = getMass(MY_IDUP(8))
!         if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
!             EvalWeighted = 0d0
!             return
!         endif
!         call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
!         call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,ML3,ML4,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
!         PSWgt = PSWgt * PSWgt2*PSWgt3
!         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then! introduce this momentum flip to allow proper mapping of integrand with Z-poles at MZ2=(p2+p3)^2 and MZ2=(p1+p4)^2
!             if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK(1:4,1),MomDK(1:4,3) )
!             PSWgt = PSWgt * 2d0
!         endif
!         if( (includeInterference.eqv..false.) .and. (OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!             PSWgt = PSWgt * 1d0/2d0
!         endif
!     elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
!         ML1=0d0; ML2=0d0; ML3=0d0; ML4=0d0
!         MomDK(1:4,1) = MomExt(1:4,3)
!         MomDK(1:4,2) = 0d0
!         MomDK(1:4,3) = MomExt(1:4,4)
!         MomDK(1:4,4) = 0d0
!     elseif( .not.IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
!         ML1 = getMass(MY_IDUP(7))
!         ML2 = getMass(MY_IDUP(6))
!         ML3=0d0; ML4=0d0
!         if( (MZ1.lt.ML1+ML2) ) then
!             EvalWeighted = 0d0
!             return
!         endif
!         call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
!         PSWgt = PSWgt * PSWgt2
!         MomDK(1:4,3) = MomExt(1:4,4)
!         MomDK(1:4,4) = 0d0
!     endif
!
!
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! !     call EvalPhasespace_H4f(yRnd(12),yRnd(3:11),EHat,MomExt(1:4,1:8),PSWgt)
! !     MomDK(1:4,1) = MomExt(1:4,5)
! !     MomDK(1:4,2) = MomExt(1:4,6)
! !     MomDK(1:4,3) = MomExt(1:4,7)
! !     MomDK(1:4,4) = MomExt(1:4,8)
! !
! !     MZ1 = get_MInv(MomExt(1:4,3))
! !     MZ2 = get_MInv(MomExt(1:4,4))
! !     if( EHat.lt.MZ1+MZ2  ) then
! !       EvalWeighted = 0d0
! !       return
! !     endif
! !         ML1 = getMass(MY_IDUP(7)) *0d0
! !         ML2 = getMass(MY_IDUP(6)) *0d0
! !         ML3 = getMass(MY_IDUP(9)) *0d0
! !         ML4 = getMass(MY_IDUP(8)) *0d0
! !         if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
! !             EvalWeighted = 0d0
! !             return
! !         endif
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!
!
!     if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
!         call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
!     else
!         call Error("This should not happen anymore")
!         call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
!         call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
!     endif
!     if( applyPSCut ) then
!       EvalWeighted = 0d0
!       return
!     endif
!
!    call SetRunningScales( (/ (MomExt(1:4,3)+MomExt(1:4,4)),Mom_Not_a_particle(1:4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) )
!    call EvalAlphaS()
!    call setPDFs(eta1,eta2,pdf)
!    FluxFac = 1d0/(2d0*EHat**2)
!
!    if (PChannel.eq.0.or.PChannel.eq.2) then
!       PDFFac = pdf(0,1) * pdf(0,2)
!
!       if (Process.eq.0) then
!             if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
!                call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
!                call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
!                if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!                   if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
!                endif
!                call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
!             else
!                call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
!             endif
!
!       elseif(Process.eq.2) then
!             if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
!                call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
!                call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
!                if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!                   if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
!                endif
!                call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
!             else
!                call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
!             endif
!       endif
!
!       LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
!       PreFac = 2d0 * hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac * FinalStateWeight
!       EvalWeighted = LO_Res_Unpol * PreFac
!
! ! EvalWeighted = PreFac  ! for PS output   (only run 1 iteration without vegas adaptation)
!
!     endif
!
!    if (PChannel.eq.1.or.PChannel.eq.2) then
!       PDFFac1 = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
!               + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
!               + pdf(Bot_,1)*pdf(ABot_,2)
!       PDFFac2 = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
!               + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
!               + pdf(Bot_,2)*pdf(ABot_,1)
!
!
!       if (Process.eq.1) then
!             if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
!                call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
!                call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
!                if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!                   if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
!                endif
!                call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
!                call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
!             else
!                call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
!                call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
!             endif
!
!       elseif(Process.eq.2) then
!             if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
!                call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
!                call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
!                if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!                   if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
!                endif
!                call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
!                call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
!             else
!                call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
!                call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
!             endif
!       endif
!
!       LO_Res_Unpol1 = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac1
!       LO_Res_Unpol2 = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac2
!       LO_Res_Unpol = LO_Res_Unpol1 + LO_Res_Unpol2
!       PreFac = 2d0 * hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * SymmFac * FinalStateWeight
!
!       EvalWeighted = LO_Res_Unpol * PreFac
!    endif
!
!    if( WidthScheme.ne.2 ) EvalWeighted = EvalWeighted * ReweightBWPropagator( Get_MInv2( MomExt(1:4,3)+MomExt(1:4,4) ) )
! !    print *, ReweightToCPS( Get_MInv( MomExt(1:4,3)+MomExt(1:4,4) ) );pause
!
!    if( writeWeightedLHE .and. (.not. warmup) ) then
!       if (PChannel.eq.0) then
!          My_IDUP(1) = Glu_
!          My_IDUP(2) = Glu_
!          ICOLUP(1:2,1) = (/501,502/)
!          ICOLUP(1:2,2) = (/502,501/)
!       endif
!       if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2))   ) then
!             call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=EvalWeighted*VgsWgt)
!       else
!             call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=EvalWeighted*VgsWgt)
!       endif
!    endif
!
!
!       do NHisto=1,NumHistograms-3
!           call intoHisto(NHisto,NBin(NHisto),EvalWeighted*VgsWgt)
!       enddo
!
!
!       call intoHisto(18,NBin(18),EvalWeighted*VgsWgt)
!
!
! RETURN
! END FUNCTION

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
real(8) :: EvalUnWeighted,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumHistograms),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,FinalStateWeight
real(8) :: MomExt(1:4,1:8),MomDK(1:4,1:4),MomDK_massless(1:4,1:4),MomExt_f(1:4,1:4),MomDK_f(1:4,1:4)
logical :: applyPSCut,genEvt
real(8) :: CS_max, channel_ratio
real(8) :: oneovervolume, bound(1:11), sumtot,yz1,yz2,EZ_max,dr,MZ1,MZ2,ML1,ML2,ML3,ML4
integer :: parton(-5:5,-5:5), i1, ifound, i2, MY_IDUP(1:9), ICOLUP(1:2,1:9),flav1,flav2, ID_DK(6:9)
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
   call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9),FinalStateWeight)
! if(  my_idup(6).eq.my_idup(8)) return! for 2e2mu   ! noi+ 0.92591989,     i+ 0.92775679,   i- 0.14072385,     ix+ 3.2770660  ix-  0.88181581E-02
! if(  my_idup(6).ne.my_idup(8)) return! for 4mu/4e  ! noi+ 0.92608512,     i+ 1.01761060,   i- 0.12915384,     ix+  3.5925481 ix-  0.80721147E-02
!                                                            50/50%              48/52%         52/48%               48%                52%


   if( RandomizeVVdecays .and. OffShellV1.eqv.OffShellV2 ) then
      if( (MY_IDUP(6).ne.MY_IDUP(8)) .and. (IsAZDecay(DecayMode1)) .and. (IsAZDecay(DecayMode2)) ) then
        if( (yrnd(18).le.0.5d0) ) then
         call swapi(MY_IDUP(4),MY_IDUP(5))
         call swapi(MY_IDUP(6),MY_IDUP(8))
         call swapi(MY_IDUP(7),MY_IDUP(9))
         call swapi(ICOLUP(1,6),ICOLUP(1,8))
         call swapi(ICOLUP(1,7),ICOLUP(1,9))
         call swapi(ICOLUP(2,6),ICOLUP(2,8))
         call swapi(ICOLUP(2,7),ICOLUP(2,9))
        endif
     elseif( (IsAWDecay(DecayMode1)) .and. (IsAWDecay(DecayMode2)) ) then
        if( (yrnd(18).le.0.5d0) ) then
         MY_IDUP(4) = ChargeFlip(MY_IDUP(4))
         MY_IDUP(5) = ChargeFlip(MY_IDUP(5))
         MY_IDUP(6) = ChargeFlip(MY_IDUP(6))
         MY_IDUP(7) = ChargeFlip(MY_IDUP(7))
         MY_IDUP(8) = ChargeFlip(MY_IDUP(8))
         MY_IDUP(9) = ChargeFlip(MY_IDUP(9))
         ! if there's a charge flip then the order of particle and anti-particles needs to be flipped, too
         call swapi(MY_IDUP(6),MY_IDUP(7))
         call swapi(MY_IDUP(8),MY_IDUP(9))
        endif
        if( (yrnd(17).le.0.5d0) ) then
         call swapi(MY_IDUP(4),MY_IDUP(5))
         call swapi(MY_IDUP(6),MY_IDUP(8))
         call swapi(MY_IDUP(7),MY_IDUP(9))
         call swapi(ICOLUP(1,6),ICOLUP(1,8))
         call swapi(ICOLUP(1,7),ICOLUP(1,9))
         call swapi(ICOLUP(2,6),ICOLUP(2,8))
         call swapi(ICOLUP(2,7),ICOLUP(2,9))
        endif
     endif ! Do not swap for VV'
  endif
  ID_DK(6:9)=MY_IDUP(6:9) ! Decay ids to pass into the ME
  if(MY_IDUP(4).eq.Pho_ .and. ID_DK(6).eq.Not_a_particle_) then
    ID_DK(6)=MY_IDUP(4)
  endif
  if(MY_IDUP(5).eq.Pho_ .and. ID_DK(8).eq.Not_a_particle_) then
    ID_DK(8)=MY_IDUP(5)
  endif

  yz1 = yRnd(10)
  yz2 = yRnd(11)
  offzchannel = yRnd(15) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * ( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )

            EZ_max = EHat - MZ1*0.99d0
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )

        elseif(M_Reso.lt.2d0*M_V_ps) then
            if (offzchannel.le.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
                MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
                MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
                sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/(  &
                1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )     &
                + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ1*0.999d0)**2
            elseif(offzchannel.gt.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
                MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
                MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
                sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/( &
                1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )    &
                + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ2*0.999d0)**2
            endif
       endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
        MZ1 = getMass(MY_IDUP(4))
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat - MZ1*0.99d0
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
        else
            MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
  endif

  elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
        MZ2 = getMass(MY_IDUP(5))
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat - MZ2*0.99d0
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
         else
            MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ2*0.999d0)**2
  endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
        MZ1 = getMass(MY_IDUP(4))
        MZ2 = getMass(MY_IDUP(5))
  endif


    if( MZ1+MZ2.gt.EHat ) then
      EvalUnWeighted = 0d0
      RejeCounter = RejeCounter + 1
      return
    endif



   call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
    if( .not.IsAPhoton(DecayMode1) .and. .not.IsAPhoton(DecayMode2) ) then ! don't decay the photon
      ML1 = getMass(MY_IDUP(6))  !*0d0
      ML2 = getMass(MY_IDUP(7))  !*0d0
      ML3 = getMass(MY_IDUP(8))  !*0d0
      ML4 = getMass(MY_IDUP(9))  !*0d0
      if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
          EvalUnWeighted = 0d0
          RejeCounter = RejeCounter + 1
          return
      endif
      call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
      call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,ML3,ML4,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
      PSWgt = PSWgt * PSWgt2*PSWgt3

      if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9))  ) then! introduce this momentum flip to allow proper mapping of integrand with Z-poles at MZ2=(p2+p3)^2 and MZ2=(p1+p4)^2
          if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK(1:4,1),MomDK(1:4,3) )
!           PSWgt = PSWgt * 2d0
      endif
    elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        ML1=0d0; ML2=0d0; ML3=0d0; ML4=0d0
        MomDK(1:4,1) = MomExt(1:4,3)
        MomDK(1:4,2) = 0d0
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
    elseif( .not.IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        ML1 = getMass(MY_IDUP(6))
        ML2 = getMass(MY_IDUP(7))
        ML3=0d0; ML4=0d0
        if( (MZ1.lt.ML1+ML2) ) then
            EvalUnWeighted = 0d0
            return
        endif
        call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
        PSWgt = PSWgt * PSWgt2
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
   endif



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!
!     call EvalPhasespace_H4f(yRnd(12),yRnd(3:11),EHat,MomExt(1:4,1:8),PSWgt)
!     MomDK(1:4,1) = MomExt(1:4,5)
!     MomDK(1:4,2) = MomExt(1:4,6)
!     MomDK(1:4,3) = MomExt(1:4,7)
!     MomDK(1:4,4) = MomExt(1:4,8)
!
!     MZ1 = get_MInv(MomExt(1:4,3))
!     MZ2 = get_MInv(MomExt(1:4,4))
!     if( EHat.lt.MZ1+MZ2  ) then
!       EvalUnWeighted = 0d0
!       return
!     endif
!         ML1 = getMass(MY_IDUP(7)) *0d0
!         ML2 = getMass(MY_IDUP(6)) *0d0
!         ML3 = getMass(MY_IDUP(9)) *0d0
!         ML4 = getMass(MY_IDUP(8)) *0d0
!         if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
!             EvalUnWeighted = 0d0
!             return
!         endif
!
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !



    if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
        call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    else
        call Error("This should not happen anymore")
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalUnWeighted = 0d0
      return
    endif

   call SetRunningScales( (/ (MomExt(1:4,3)+MomExt(1:4,4)),Mom_Not_a_particle(1:4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) )
   call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)



IF( GENEVT ) THEN
!   csmax(0,0) is rescaled by par_adj
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



!   yrnd(13) selects the partonic channel according to its relative contribution to the total cross section
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
      CS_max = csmax(0,0)*adj_par   ! this is necessary here for the follow up = undo rescaling after partonic channel has been chosen
      MY_IDUP(1:2)=(/Glu_,Glu_/)
      ICOLUP(1:2,1) = (/501,502/)
      ICOLUP(1:2,2) = (/502,501/)
      PDFFac = pdf(0,1) * pdf(0,2)

      if (Process.eq.0) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            endif
      else
            LO_Res_Unpol = 0d0
      endif
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2 * PDFFac

   else

      if (Process.eq.1) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            endif
      else
         LO_Res_Unpol1 = 0d0
         LO_Res_Unpol2 = 0d0
      endif


      if (ifound.eq.2) then
          PDFFac = pdf(Bot_,2)*pdf(ABot_,1)
          i2 = -5
          MY_IDUP(1:2)=(/ABot_,Bot_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
          LO_Res_Unpol = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.3) then
          PDFFac = pdf(Chm_,2)*pdf(AChm_,1)
          i2 = -4
          MY_IDUP(1:2)=(/AChm_,Chm_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
          LO_Res_Unpol = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.4) then
          PDFFac = pdf(Str_,2)*pdf(AStr_,1)
          i2 = -3
          MY_IDUP(1:2)=(/AStr_,Str_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
          LO_Res_Unpol = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.5) then
          PDFFac = pdf(Up_,2) *pdf(AUp_,1)
          i2 = -2
          MY_IDUP(1:2)=(/AUp_,Up_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
          LO_Res_Unpol = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.6) then
          PDFFac = pdf(Dn_,2) *pdf(ADn_,1)
          i2 = -1
          MY_IDUP(1:2)=(/ADn_,Dn_/)
          ICOLUP(1:2,1) = (/0,502/)
          ICOLUP(1:2,2) = (/502,0/)
          LO_Res_Unpol = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif (ifound.eq.7) then
          PDFFac = pdf(Dn_,1) *pdf(ADn_,2)
          i2 = 1
          MY_IDUP(1:2)=(/Dn_,ADn_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
          LO_Res_Unpol = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.8) then
          PDFFac = pdf(Up_,1) *pdf(AUp_,2)
          i2 = 2
          MY_IDUP(1:2)=(/Up_,AUp_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
          LO_Res_Unpol = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.9) then
          PDFFac = pdf(Str_,1)*pdf(AStr_,2)
          i2 = 3
          MY_IDUP(1:2)=(/Str_,AStr_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
          LO_Res_Unpol = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac
      elseif(ifound.eq.10) then
          PDFFac = pdf(Chm_,1)*pdf(AChm_,2)
          i2 = 4
          MY_IDUP(1:2)=(/Chm_,AChm_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
          LO_Res_Unpol = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac
       elseif(ifound.eq.11) then
          PDFFac = pdf(Bot_,1)*pdf(ABot_,2)
          i2 = 5
          MY_IDUP(1:2)=(/Bot_,ABot_/)
          ICOLUP(1:2,1) = (/501,0/)
          ICOLUP(1:2,2) = (/0,501/)
          LO_Res_Unpol = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac
       endif
       parton(i2,-i2) = 1
       CS_max = csmax(i2,-i2)
   endif

   PreFac = 2d0 * hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * SymmFac * FinalStateWeight
   EvalUnWeighted = LO_Res_Unpol * PreFac
   if( WidthScheme.ne.2 ) EvalUnWeighted = EvalUnWeighted * ReweightBWPropagator( Get_MInv2( MomExt(1:4,3)+MomExt(1:4,4) ) )

      if( EvalUnWeighted.gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted .gt. yRnd(14)*CS_max ) then
         do NHisto=1,NumHistograms-3
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo

         AccepCounter = AccepCounter + 1
         AccepCounter_part(-5:+5,-5:5) = AccepCounter_part(-5:+5,-5:5)  + parton(-5:+5,-5:5)
         call BranchingCounter(MY_IDUP(6:9))
         if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
              call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
          else
              call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
         endif


          if( abs(MY_IDUP(7)).eq.abs(ElM_) ) then
              flav1 = 1
          elseif( abs(MY_IDUP(7)).eq.abs(MuM_) ) then
              flav1 = 2
          elseif( abs(MY_IDUP(7)).eq.abs(TaM_) ) then
              flav1 = 3
          elseif( IsANeutrino(MY_IDUP(7)) ) then
              flav1 = 4
          elseif( IsAQuark(MY_IDUP(7)) ) then
              flav1 = 5
          else! for photons
              flav1 = 1
          endif
          if( abs(MY_IDUP(9)).eq.abs(ElM_) ) then
              flav2 = 1
          elseif( abs(MY_IDUP(9)).eq.abs(MuM_) ) then
              flav2 = 2
          elseif( abs(MY_IDUP(9)).eq.abs(TaM_) ) then
              flav2 = 3
          elseif( IsANeutrino(MY_IDUP(9)) ) then
              flav2 = 4
          elseif( IsAQuark(MY_IDUP(9)) ) then
              flav2 = 5
          else! for photons
              flav2 = 1
          endif
          Br_counter(flav1,flav2) = Br_counter(flav1,flav2) + 1

      else
          RejeCounter = RejeCounter + 1
      endif

ELSE! NOT GENEVT


   if (PChannel.eq.0.or.PChannel.eq.2) then
      PDFFac = pdf(0,1) * pdf(0,2)
      if (Process.eq.0) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol)
            endif
      else
         LO_Res_Unpol = 0d0
      endif
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2


      PreFac = 2d0 * hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac * FinalStateWeight
      EvalUnWeighted = LO_Res_Unpol * PreFac
      if( WidthScheme.ne.2 ) EvalUnWeighted = EvalUnWeighted * ReweightBWPropagator( Get_MInv2( MomExt(1:4,3)+MomExt(1:4,4) ) )
      RES(0,0) = EvalUnWeighted

      if (EvalUnWeighted.gt.csmax(0,0)) then
          csmax(0,0) = EvalUnWeighted
      endif
   endif

   if (PChannel.eq.1.or.PChannel.eq.2) then
      if (Process.eq.1) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),ID_DK(6:9),LO_Res_Unpol2)
            endif
      else
         LO_Res_Unpol1 = 0d0
         LO_Res_Unpol2 = 0d0
      endif

      LO_Res_Unpol1 = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2
      LO_Res_Unpol2 = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2
      PreFac = 2d0 * hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt *   SymmFac

      do i1 = -5,5
         if (i1.eq.-5) then
          PDFFac = pdf(Bot_,2)*pdf(ABot_,1)
          EvalUnWeighted = LO_Res_Unpol2 * PreFac *PDFFac
         elseif(i1.eq.-4) then
          PDFFac = pdf(Chm_,2)*pdf(AChm_,1)
          EvalUnWeighted = LO_Res_Unpol2 * PreFac *PDFFac
         elseif(i1.eq.-3) then
          PDFFac = pdf(Str_,2)*pdf(AStr_,1)
          EvalUnWeighted = LO_Res_Unpol2 * PreFac *PDFFac
         elseif(i1.eq.-2) then
          PDFFac = pdf(Up_,2) *pdf(AUp_,1)
          EvalUnWeighted = LO_Res_Unpol2 * PreFac *PDFFac
         elseif(i1.eq.-1) then
          PDFFac = pdf(Dn_,2) *pdf(ADn_,1)
          EvalUnWeighted = LO_Res_Unpol2 * PreFac *PDFFac
         elseif (i1.eq.0) then
          PDFFac = 0d0
         elseif (i1.eq.1) then
          PDFFac = pdf(Dn_,1) *pdf(ADn_,2)
          EvalUnWeighted = LO_Res_Unpol1 * PreFac *PDFFac
         elseif (i1.eq.2) then
          PDFFac = pdf(Up_,1) *pdf(AUp_,2)
          EvalUnWeighted = LO_Res_Unpol1 * PreFac *PDFFac
         elseif(i1.eq.3) then
          PDFFac = pdf(Str_,1)*pdf(AStr_,2)
          EvalUnWeighted = LO_Res_Unpol1 * PreFac *PDFFac
         elseif(i1.eq.4) then
          PDFFac = pdf(Chm_,1)*pdf(AChm_,2)
          EvalUnWeighted = LO_Res_Unpol1 * PreFac *PDFFac
         elseif(i1.eq.5) then
          PDFFac = pdf(Bot_,1)*pdf(ABot_,2)
          EvalUnWeighted = LO_Res_Unpol1 * PreFac *PDFFac
         endif

          RES(i1,-i1) = EvalUnWeighted
          if (EvalUnWeighted.gt.csmax(i1,-i1)) csmax(i1,-i1) = EvalUnWeighted
      enddo
   endif

ENDIF! genEvt


RETURN
END FUNCTION


Function EvalWeighted_HJ(yRnd,VgsWgt)
 use ModKinematics
 use ModParameters
 use ModHiggsJ
 use ModMisc
 !use ModMisc
#if compiler==1
 use ifport
#endif
    implicit none
    real(8) :: yRnd(1:7)
    real(8) :: VgsWgt, EvalWeighted_HJ
    real(8) :: pdf(-6:6,1:2)
    real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
    real(8) :: MomExt(1:4,1:4), Masses(1:2), PSWgt, PSWgt2
    real(8) :: me2(-5:5,-5:5), lheweight(-5:5,-5:5)
    integer :: i,j,MY_IDUP(1:4),ICOLUP(1:2,1:4),NBin(1:NumHistograms),NHisto
    real(8) :: LO_Res_Unpol, PreFac
    logical :: applyPSCut

    EvalWeighted_HJ=0d0

    call PDFMapping(16,yrnd(1:2),eta1,eta2,Ehat,sHatJacobi) !!!!efficiency improvement as ~1/s?
    if (EHat.lt.M_Reso) return

    EvalCounter = EvalCounter+1

    Masses(1) = M_Reso
    Masses(2) = zero
    call EvalPhasespace_2to2(EHat,Masses,yRnd(3:4),MomExt,PSWgt)
    call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
    call Kinematics_HJ(4,MomExt,applyPSCut,NBin)
    if( applyPSCut .or. PSWgt.eq.zero ) return

    MY_IDUP(1:4)  = (/Glu_,Glu_,Hig_,Glu_/) ! weighted events, dominated by gg>Hg
    ICOLUP(1:2,1) = (/502,501/)!weighted events, dominated by gg>Hg
    ICOLUP(1:2,2) = (/503,502/)!weighted events, dominated by gg>Hg
    ICOLUP(1:2,3) = (/000,000/)!weighted events, dominated by gg>Hg
    ICOLUP(1:2,4) = (/503,501/)!weighted events, dominated by gg>Hg

    call SetRunningScales( (/ MomExt(1:4,3),MomExt(1:4,4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,MY_IDUP(4),Not_a_particle_,Not_a_particle_ /) )
    call EvalAlphaS()
    call setPDFs(eta1,eta2,pdf)
    call EvalAmp_HJ(MomExt,me2)

    LO_Res_Unpol = 0d0
    do i = -5,5
      do j = -5,5
         LO_Res_Unpol = LO_Res_Unpol + me2(i,j)*pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2)
      enddo
    enddo
!print *, me2(0,0)*pdf(LHA2M_pdf(0),1)*pdf(LHA2M_pdf(0),2), me2(0,1)*pdf(LHA2M_pdf(0),1)*pdf(LHA2M_pdf(1),2)
    FluxFac = 1d0/(2d0*(EHat)**2)
    PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt
    EvalWeighted_HJ = LO_Res_Unpol * PreFac

    AccepCounter=AccepCounter+1

    if( writeWeightedLHE ) then

      call WriteOutEvent_HJ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),MY_IDUP(1:4),ICOLUP(1:2,1:4),EventWeight=EvalWeighted_HJ*VgsWgt)
    endif

    do NHisto = 1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalWeighted_HJ*VgsWgt)
    enddo
!print *, EvalWeighted_HJ, EvalWeighted_HJ*VgsWgt
!if(EvalWeighted_HJ.lt.0.0d0)then
!print *, EvalWeighted_HJ
!print *, EvalWeighted_HJ, eta1, eta2, PSWgt
!print *, PSWgt
!print *, Ehat, LHC_Energy
!print *, MomExt
!print *, "-------------"
!print *, pdf
!print *, "-------------"
!print *, me2
!print *, "-------------"
!endif

   RETURN

 end Function EvalWeighted_HJ

FUNCTION EvalUnWeighted_HJ(yRnd,genEvt,RES)
 use ModKinematics
 use ModParameters
 use ModHiggsJ
 use ModMisc
#if compiler==1
 use ifport
#endif
implicit none
real(8) :: yRnd(:),VgsWgt, EvalUnWeighted_HJ,RES(-5:5,-5:5)
real(8) :: pdf(-6:6,1:2)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:5), Masses(1:2), PSWgt!,MomDK(1:4,1:4)
real(8) :: me2(-5:5,-5:5),bound(0:121)
integer :: i,j,k,ifound,jfound
integer :: MY_IDUP(1:4),ICOLUP(1:2,1:4),NBin(1:NumHistograms),NHisto
real(8) :: LO_Res_Unpol, PreFac, CS_max, sumtot
logical :: applyPSCut,genEVT
include 'csmaxvalue.f'

   EvalUnWeighted_HJ = 0d0

   call PDFMapping(16,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)

   if (EHat.lt.M_Reso) return
   Masses(1) = M_Reso
   Masses(2) = zero
   call EvalPhasespace_2to2(EHat,Masses,yRnd(3:4),MomExt,PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   call Kinematics_HJ(4,MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return

   FluxFac = 1d0/(2d0*EHat**2)
   EvalCounter = EvalCounter+1


IF( GENEVT ) THEN

      sumtot = 0d0
      do i = -5,5
         do j = -5,5
            sumtot = sumtot + csmax(i,j)
         enddo
      enddo
      !print *,csmax(0,0), csmax(0,1)
      k=0; bound(0)=0d0
      do i = -5,5
         do j = -5,5
            k=k+1
            bound(k) = bound(k-1) + csmax(i,j)/sumtot
            if( yRnd(8).gt.bound(k-1) .and. yRnd(8).lt.bound(k)  ) then
               ifound=i; jfound=j;
               goto 1313
            endif
         enddo
      enddo
 1313 continue
      !print *,ifound, jfound


      MY_IDUP(1:3) = (/LHA2M_ID(ifound),LHA2M_ID(jfound),Hig_/)
      ICOLUP(1:2,3) = (/000,000/)! H
! color assignments may need correction/randomization
      if( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).eq.Glu_ ) then! gg->Hg
          if(yRnd(10).gt.0.5d0)then
            ICOLUP(1:2,1) = (/502,501/)
            ICOLUP(1:2,2) = (/503,502/)
            ICOLUP(1:2,4) = (/503,501/)
          else
            ICOLUP(1:2,1) = (/501,502/)
            ICOLUP(1:2,2) = (/502,503/)
            ICOLUP(1:2,4) = (/501,503/)
          endif
          MY_IDUP(4) = Glu_
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).gt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qg->Hq
          ICOLUP(1:2,1) = (/502,000/)
          ICOLUP(1:2,2) = (/501,502/)
          ICOLUP(1:2,4) = (/501,000/)
          MY_IDUP(4) = MY_IDUP(1)
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).lt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qbg->Hqb
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,2) = (/501,502/)
          ICOLUP(1:2,4) = (/000,502/)
          MY_IDUP(4) = MY_IDUP(1)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).gt.0 ) then! gq->Hq
          ICOLUP(1:2,1) = (/501,502/)
          ICOLUP(1:2,2) = (/502,000/)
          ICOLUP(1:2,4) = (/501,000/)
          MY_IDUP(4) = MY_IDUP(2)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).lt.0 ) then! gqb->Hqb
          ICOLUP(1:2,1) = (/501,502/)
          ICOLUP(1:2,2) = (/000,501/)
          ICOLUP(1:2,4) = (/000,502/)
          MY_IDUP(4) = MY_IDUP(2)
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(1).ne.Glu_ .and. MY_IDUP(2).lt.0 ) then! qqb->Hg
          ICOLUP(1:2,1) = (/502,000/)
          ICOLUP(1:2,2) = (/000,501/)
          ICOLUP(1:2,4) = (/502,501/)
          MY_IDUP(4) = Glu_
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).gt.0 .and. MY_IDUP(2).ne.Glu_ ) then! qbq->Hg
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,2) = (/502,000/)
          ICOLUP(1:2,4) = (/502,501/)
          MY_IDUP(4) = Glu_
      endif

      call SetRunningScales( (/ MomExt(1:4,3),MomExt(1:4,4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Glu_,Not_a_particle_,Not_a_particle_ /) ) ! Here, pass Glu_ since all partons are treated as massless
      call EvalAlphaS()
      call setPDFs(eta1,eta2,pdf)
      call EvalAmp_HJ(MomExt,me2)

      PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt

      LO_Res_Unpol =  me2(ifound,jfound) * pdf(LHA2M_pdf(ifound),1)*pdf(LHA2M_pdf(jfound),2)
      EvalUnWeighted_HJ = LO_Res_Unpol * PreFac

      if( ifound.eq.0 .and. jfound.eq.0 ) then
          CS_max = csmax(ifound,jfound) * adj_par
      else
          CS_max = csmax(ifound,jfound)
      endif

      if( EvalUnWeighted_HJ.gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_HJ, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_HJ, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted_HJ .gt. yRnd(9)*CS_max ) then
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1
         call WriteOutEvent_HJ((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4)/),MY_IDUP(1:4),ICOLUP(1:2,1:4))
      else
          RejeCounter = RejeCounter + 1
      endif
!print *, ifound, jfound

ELSE! NOT GENEVT


   call SetRunningScales( (/ MomExt(1:4,3),MomExt(1:4,4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Glu_,Not_a_particle_,Not_a_particle_ /) )
   call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)
   call EvalAmp_HJ(MomExt,me2)
   PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt


   LO_Res_Unpol = 0d0
   do i = -5,5
      do j = -5,5

         LO_Res_Unpol = me2(i,j) * PreFac *pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2)

         EvalUnWeighted_HJ = EvalUnWeighted_HJ  + LO_Res_Unpol

          RES(i,j) = LO_Res_Unpol
          if (LO_Res_Unpol.gt.csmax(i,j)) then
              csmax(i,j) = LO_Res_Unpol
          endif
      enddo
   enddo
!print *,csmax(0,0), csmax(0,1)


ENDIF! GENEVT


 RETURN
 END FUNCTION EvalUnWeighted_HJ






Function EvalWeighted_VHiggs(yRnd,VgsWgt)
 use ModKinematics
 use ModParameters
 use ModVHiggs
 use ModMisc
 !use ModMisc
#if compiler==1
 use ifport
#endif
    implicit none
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5): flavor in Z/W decay mode
! yrnd(6:7): cos(theta_4) and phi_4 in the CM frame of Z*(3)
! yrnd(8:9): cos(theta_6) and phi_6 in the CM frame of decay product of Z(4)
! yrnd(10:11): cos(theta_8) and phi_8 in the CM frame of decay product of H(5)
! yRnd(12): inv_mass(4)
! yRnd(13): inv_mass(5)
! yrnd(14:15): PDF mapping
! yrnd(16): partonic channel in unweighted events
! yRnd(17): accept or reject in unweighted events
    real(8) :: yRnd(1:20),VgsWgt, EvalWeighted_VHiggs
    real(8) :: pdf(-6:6,1:2)
    real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
    real(8) :: MomExt(1:4,1:9), PSWgt, PSWgt2,NewME2(1:6)
    real(8) :: me2!, lheweight(-6:6,-6:6)
    integer :: i,j,k,h,NBin(1:NumHistograms),NHisto
    real(8) :: DKWgt, LO_Res_Unpol, PreFac
    logical :: applyPSCut !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!phase space cuts?
    real(8) :: cyRnd(4)

    real(8) :: inv_mass(9),mass(9,2)
    !double precision beam_momentum(2,4), four_momentum(7,4),inv_mass(7),mass(7,2)
    real(8) :: helicity(9)!, beam_h(2) !helicities
    integer :: id(9), id2(9)!, beam_id(2)
    integer :: tmp_idup(1:3), tmp_icolup(1:2,1:2)
    integer :: nVhels


    EvalWeighted_VHiggs=0d0
    EvalCounter = EvalCounter+1

    DKWgt=1
    id(:)=0
    helicity(:)=0
    mass(1:2,1:2)=0d0
    mass(3,1)=M_V
    mass(3,2)=Ga_V
    mass(4,1)=M_V
    mass(4,2)=Ga_V
    nVhels=2
    if(IsAPhoton(DecayMode1))then
       if (includeVprime .and. Ga_Vprime.gt.0d0) then
          mass(3,1)=(M_V+M_Vprime)/2d0
          mass(3,2)=sqrt(abs((M_V-M_Vprime)/2d0)**2+Ga_V**2+Ga_Vprime**2)
       endif
       mass(4,1)=getMass(Pho_)
       mass(4,2)=getDecayWidth(Pho_)
    else
       if (includeVprime .and. Ga_Vprime.gt.0d0) then
          mass(3,1)=(M_V+M_Vprime)/2d0
          mass(3,2)=sqrt(abs((M_V-M_Vprime)/2d0)**2+Ga_V**2+Ga_Vprime**2)
          mass(4,:)=mass(3,:)
       endif
    endif
    mass(5,1)=M_Reso
    mass(5,2)=Ga_Reso
    mass(6:9,1:2)=0d0


    id(5)=convertLHE(Hig_)
    id(8)=convertLHE(Bot_)
    id(9)=-id(8)
!toss coin and decide beam A helicity
        if (yRnd(1).lt.(0.5d0+POL_A/200d0))then
          helicity(1)=1d0
        else
          helicity(1)=-1d0
        endif
!toss coin and decide beam B helicity
        if (yRnd(2).lt.(0.5d0+POL_B/200d0))then
          helicity(2)=1d0
        else
          helicity(2)=-1d0
        endif
!toss coin and decide particle 8,9 helicities
        if (yRnd(4).gt.0.5d0)then
          helicity(8)=1d0
        else
          helicity(8)=-1d0
        endif
        helicity(9)=helicity(8)
!toss coin and decide particle 6,7 helicities
        if (yRnd(3).gt.0.5d0)then
          helicity(6)=1d0
        else
          helicity(6)=-1d0
        endif
        helicity(7)=-helicity(6)

    call VBranching(DecayMode1, tmp_idup, tmp_icolup, DKWgt)
    id(3)=convertLHE(tmp_idup(1))
    id(4)=convertLHE(tmp_idup(1))
    id(6:7)=Not_a_particle_
    if (tmp_idup(2).ne.Not_a_particle_)then
      id(6)=convertLHE(tmp_idup(2))
    elseif(DecayMode1.eq.7) then
      id(3)=convertLHE(Z0_)
      id(6)=convertLHE(Pho_)
    endif
    if (tmp_idup(3).ne.Not_a_particle_)then
      id(7)=convertLHE(tmp_idup(3))
    endif
    if (id(6).eq.-id(7) .and. id(6).lt.0) then
      call swap(id(6),id(7))
    endif

    if((IsAWDecay(DecayMode1) .or. IsALHENeutrino(id(6)) .or. IsALHENeutrino(id(7))) .and. .not. includeVprime) then
      helicity(6)=sign(1d0,-dble(id(6)))
      helicity(7)=-helicity(6)
      nVhels=1
    endif
    if(H_DK) then
      DkWgt = DKWgt*6d0 ! H->bb decay
    endif

    DkWgt = DkWgt * dble(nVhels) ! 2 possible helicities

!print*, id(6:7),helicity(6:7)

if( IsAZDecay(DecayMode1) .or. IsAPhoton(DecayMode1) ) then
!if pp collider
    if(Collider.eq.1)then
      call PDFMapping(14,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      MomExt(1,3)=EHat
      MomExt(2,3)=0d0
      MomExt(3,3)=0d0
      MomExt(4,3)=0d0

      MomExt(1,1)=EHat/2d0
      MomExt(2,1)=0d0
      MomExt(3,1)=0d0
      MomExt(4,1)=MomExt(1,1)
      MomExt(1,2)=EHat/2d0
      MomExt(2,2)=0d0
      MomExt(3,2)=0d0
      MomExt(4,2)=-MomExt(1,2)
!print*,includeGammaStar
      call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,IsAPhoton(DecayMode1))
      call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,useAonshell=IsAPhoton(DecayMode1))
      if( applyPSCut .or. PSWgt.eq.zero ) return

      if(IsAZDecay(DecayMode1))then
        call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      else !IsAPhoton(DecayMode1)
        call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
      endif
      call setPDFs(eta1,eta2,pdf)
      FluxFac = 1d0/(2d0*EHat**2)
      PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * DKWgt
      LO_Res_Unpol=0d0
      EvalWeighted_VHiggs=0d0
      do i = -6,6
        j = -i
        id(1:2) = (/LHA2M_PDF(i),LHA2M_PDF(j)/)
        if (abs(LHA2M_PDF(i)).ne.6   .and.   abs(LHA2M_PDF(j)).ne.6.  .and.  i.ne.0)then
          call EvalAmp_VHiggs(id,helicity,MomExt,me2)
          if(IsNaN(me2).or.(me2.eq.0d0))return
          !if(IsNaN(me2))return
          if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90


!     if( me2.lt.1d-8 ) cycle
!     NewME2(1:6) = cdabs(TreeAmp_qqb_ZH((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,6),MomExt(1:4,7)/),int((/helicity(1),helicity(2),helicity(6),helicity(7)/))))**2
!     print *, "heli",helicity(1),helicity(2),helicity(6),helicity(7)
!     print *, "Mes", me2,NewME2(up_)
!     print *, "ratio", me2/NewME2(up_)
!     pause




        else
          me2=0d0
        endif
          LO_Res_Unpol = me2/3d0*pdf(i,1)*pdf(j,2)* PreFac
          EvalWeighted_VHiggs = EvalWeighted_VHiggs + LO_Res_Unpol
          !lheweight(i,j)=LO_Res_Unpol

      enddo

!if e+ e- collider
    else if(Collider.eq.0)then
      MomExt(1,3)=ILC_Energy
      MomExt(2,3)=0d0
      MomExt(3,3)=0d0
      MomExt(4,3)=0d0

      MomExt(1,1)=ILC_Energy/2d0
      MomExt(2,1)=0d0
      MomExt(3,1)=0d0
      MomExt(4,1)=MomExt(1,1)
      MomExt(1,2)=ILC_Energy/2d0
      MomExt(2,2)=0d0
      MomExt(3,2)=0d0
      MomExt(4,2)=-MomExt(1,2)

      call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,IsAPhoton(DecayMode1))
      call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,useAonshell=IsAPhoton(DecayMode1))
      if( applyPSCut .or. PSWgt.eq.zero ) return

      FluxFac = 1d0/(2d0*ILC_Energy**2)
      PreFac = hbarc2XsecUnit * FluxFac * PSWgt * DKWgt
      LO_Res_Unpol=0d0
      EvalWeighted_VHiggs=0d0
      id(2)=convertLHE(ElM_)
      id(1)=-id(2)
      call EvalAmp_VHiggs(id,helicity,MomExt,me2)
      if(IsNaN(me2).or.(me2.eq.0d0))return
      if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
      LO_Res_Unpol =me2 * PreFac
      EvalWeighted_VHiggs = LO_Res_Unpol

    endif

elseif( IsAWDecay(DecayMode1) ) then
      call PDFMapping(15,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      MomExt(1,3)=EHat
      MomExt(2,3)=0d0
      MomExt(3,3)=0d0
      MomExt(4,3)=0d0

      MomExt(1,1)=EHat/2d0
      MomExt(2,1)=0d0
      MomExt(3,1)=0d0
      MomExt(4,1)=MomExt(1,1)
      MomExt(1,2)=EHat/2d0
      MomExt(2,2)=0d0
      MomExt(3,2)=0d0
      MomExt(4,2)=-MomExt(1,2)

      call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,IsAPhoton(DecayMode1))

      call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut)
      if( applyPSCut .or. PSWgt.eq.zero ) return

      call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      call setPDFs(eta1,eta2,pdf)
      FluxFac = 1d0/(2d0*EHat**2)
      PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * DKWgt
      LO_Res_Unpol=0d0
      EvalWeighted_VHiggs=0d0
      do i = -5,5
      do j = -5,5
         if (i.eq.0 .or. j.eq.0) cycle

         id2=id
         id2(1:2) = (/i,j/)

         if (CoupledVertex(id2(1:2),-1) .ne. CoupledVertex(id2(6:7),-1)) cycle

         call EvalAmp_VHiggs(id2,helicity,MomExt,me2)
         if(IsNaN(me2).or.(me2.eq.0d0))return
         if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90

         LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac
         EvalWeighted_VHiggs = EvalWeighted_VHiggs+LO_Res_Unpol
      enddo
      enddo

endif

   cyRnd(1)=yRnd(9)
   cyRnd(2)=yRnd(8)
   if(.not.IsAPhoton(DecayMode1)) then
      call EvalPhasespace_VDecay(MomExt(1:4,4),inv_mass(4),getMass(convertLHEreverse(id(6))),getMass(convertLHEreverse(id(7))),cyRnd(1:2),MomExt(1:4,6:7),PSWgt2)
   endif
   if(Collider.eq.1)then
      call boost2Lab(eta1,eta2,9,MomExt(1:4,1:9))
   endif
   do i=6,7
      inv_mass(i)=dsqrt(dabs(MomExt(1:4,i).dot.MomExt(1:4,i)))
   enddo

   AccepCounter=AccepCounter+1

   if( writeWeightedLHE ) then
      if( IsAZDecay(DecayMode1) .or. IsAPhoton(DecayMode1) ) then
         if(Collider.eq.1)then
            id(1:2) = (/convertLHE(Up_),convertLHE(AUp_)/)
!if e+ e- collider
         else if(Collider.eq.0)then
            id(1:2) = (/convertLHE(ElP_),convertLHE(ElM_)/)
         endif
      elseif( IsAWDecay(DecayMode1) ) then
         id(1:2) = (/convertLHE(Up_),convertLHE(ADn_)/)
      endif
      if( IsNaN(EvalWeighted_VHiggs) ) then
         EvalWeighted_VHiggs = 0d0
         return
      endif
      if(EvalWeighted_VHiggs.ne.0d0)then
         call WriteOutEvent_VHiggs(id,helicity,MomExt,inv_mass,EventWeight=EvalWeighted_VHiggs*VgsWgt)
      endif
   endif

   do NHisto = 1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalWeighted_VHiggs*VgsWgt)
   enddo

   RETURN

 end Function EvalWeighted_VHiggs

Function EvalUnWeighted_VHiggs(yRnd,genEvt,RES)
use ModKinematics
 use ModParameters
 use ModVHiggs
 use ModMisc
#if compiler==1
 use ifport
#endif
implicit none
real(8) :: yRnd(1:21), EvalUnWeighted_VHiggs, RES(-5:5,-5:5)
real(8) :: pdf(-6:6,1:2)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:9), PSWgt, PSWgt2
real(8) :: me2,bound(0:121)
integer :: i,j,k,ifound,jfound
integer :: NBin(1:NumHistograms),NHisto
real(8) :: DKWgt, LO_Res_Unpol, PreFac, CS_max, sumtot
logical :: applyPSCut,genEVT
real(8) :: cyRnd(4)
real(8) :: inv_mass(9),mass(9,2)
!real(8) :: beam_momentum(2,4), four_momentum(7,4),inv_mass(7),mass(7,2)
real(8) :: helicity(9) !helicities
integer :: id(9), id2(9)
integer :: tmp_idup(1:3), tmp_icolup(1:2,1:2)
include 'csmaxvalue.f'

EvalUnWeighted_VHiggs = 0d0

DKWgt=1d0
id(:)=0
helicity(:)=0

mass(1:2,1:2)=0d0
mass(3,1)=M_V
mass(3,2)=Ga_V
mass(4,1)=M_V
mass(4,2)=Ga_V
if(IsAPhoton(DecayMode1))then
  if (includeVprime .and. Ga_Vprime.gt.0d0) then
    mass(3,1)=(M_V+M_Vprime)/2d0
    mass(3,2)=sqrt(abs((M_V-M_Vprime)/2d0)**2+Ga_V**2+Ga_Vprime**2)
  endif
  mass(4,1)=getMass(Pho_)
  mass(4,2)=getDecayWidth(Pho_)
else
  if (includeVprime .and. Ga_Vprime.gt.0d0) then
    mass(3,1)=(M_V+M_Vprime)/2d0
    mass(3,2)=sqrt(abs((M_V-M_Vprime)/2d0)**2+Ga_V**2+Ga_Vprime**2)
    mass(4,:)=mass(3,:)
  endif
endif
mass(5,1)=M_Reso
mass(5,2)=Ga_Reso
mass(6:9,1:2)=0d0

id(5)=convertLHE(Hig_)
id(8)=convertLHE(Bot_)
id(9)=-id(8)
!toss coin and decide beam A helicity
if (yRnd(1).lt.(0.5d0+POL_A/200d0))then
  helicity(1)=1d0
else
  helicity(1)=-1d0
endif
!toss coin and decide beam B helicity
if (yRnd(2).lt.(0.5d0+POL_B/200d0))then
  helicity(2)=1d0
else
  helicity(2)=-1d0
endif
!toss coin and decide particle 8,9 helicities
if (yRnd(4).gt.0.5d0)then
  helicity(8)=1d0
else
  helicity(8)=-1d0
endif
helicity(9)=helicity(8)
!toss coin and decide particle 6,7 helicities
if (yRnd(3).gt.0.5d0)then
  helicity(6)=1d0
else
  helicity(6)=-1d0
endif
helicity(7)=-helicity(6)

call VBranching(DecayMode1, tmp_idup, tmp_icolup, DKWgt)
id(3)=convertLHE(tmp_idup(1))
id(4)=convertLHE(tmp_idup(1))
id(6:7)=Not_a_particle_
if (tmp_idup(2).ne.Not_a_particle_)then
   id(6)=convertLHE(tmp_idup(2))
elseif(DecayMode1.eq.7) then
   id(3)=convertLHE(Z0_)
   id(6)=convertLHE(Pho_)
endif
if (tmp_idup(3).ne.Not_a_particle_)then
   id(7)=convertLHE(tmp_idup(3))
endif
if (id(6).eq.-id(7) .and. id(6).lt.0) then
   call swap(id(6),id(7))
endif
if(IsAWDecay(DecayMode1)) then
   helicity(6)=sign(1d0,-dble(id(6)))
   helicity(7)=-helicity(6)
endif
if(H_DK) then
   DkWgt = DKWgt*6d0 ! H->bb decay
endif



if( IsAZDecay(DecayMode1) .or. IsAPhoton(DecayMode1)) then
!if pp collider
    if(Collider.eq.1)then
      call PDFMapping(14,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      MomExt(1,3)=EHat
      MomExt(2,3)=0d0
      MomExt(3,3)=0d0
      MomExt(4,3)=0d0

      MomExt(1,1)=EHat/2d0
      MomExt(2,1)=0d0
      MomExt(3,1)=0d0
      MomExt(4,1)=MomExt(1,1)
      MomExt(1,2)=EHat/2d0
      MomExt(2,2)=0d0
      MomExt(3,2)=0d0
      MomExt(4,2)=-MomExt(1,2)

      call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,IsAPhoton(DecayMode1))
      call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,useAonshell=IsAPhoton(DecayMode1))
      if( applyPSCut .or. PSWgt.eq.zero ) return

      if(IsAZDecay(DecayMode1))then
        call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      else !IsAPhoton(DecayMode1)
        call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
      endif
      call setPDFs(eta1,eta2,pdf)
      FluxFac = 1d0/(2d0*EHat**2)
      PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * DKWgt

!if e+ e- collider
    else if(Collider.eq.0)then
      MomExt(1,3)=ILC_Energy
      MomExt(2,3)=0d0
      MomExt(3,3)=0d0
      MomExt(4,3)=0d0

      MomExt(1,1)=ILC_Energy/2d0
      MomExt(2,1)=0d0
      MomExt(3,1)=0d0
      MomExt(4,1)=MomExt(1,1)
      MomExt(1,2)=ILC_Energy/2d0
      MomExt(2,2)=0d0
      MomExt(3,2)=0d0
      MomExt(4,2)=-MomExt(1,2)

      call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,IsAPhoton(DecayMode1))
      call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,useAonshell=IsAPhoton(DecayMode1))
      if( applyPSCut .or. PSWgt.eq.zero ) return

      FluxFac = 1d0/(2d0*ILC_Energy**2)
      PreFac = hbarc2XsecUnit * FluxFac * PSWgt * DKWgt
    endif

elseif( IsAWDecay(DecayMode1) ) then
      call PDFMapping(15,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      MomExt(1,3)=EHat
      MomExt(2,3)=0d0
      MomExt(3,3)=0d0
      MomExt(4,3)=0d0

      MomExt(1,1)=EHat/2d0
      MomExt(2,1)=0d0
      MomExt(3,1)=0d0
      MomExt(4,1)=MomExt(1,1)
      MomExt(1,2)=EHat/2d0
      MomExt(2,2)=0d0
      MomExt(3,2)=0d0
      MomExt(4,2)=-MomExt(1,2)

      call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,IsAPhoton(DecayMode1))
      call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut)
      if( applyPSCut .or. PSWgt.eq.zero ) return

      call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      call setPDFs(eta1,eta2,pdf)
      FluxFac = 1d0/(2d0*EHat**2)
      PreFac = hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * DKWgt
endif

EvalCounter = EvalCounter+1

IF( GENEVT ) THEN

   sumtot = 0d0
   do i = -5,5
      do j = -5,5
         sumtot = sumtot + csmax(i,j)
      enddo
   enddo

   k=0; bound(0)=0d0
   do i = -5,5
      do j = -5,5
         k=k+1
         bound(k) = bound(k-1) + csmax(i,j)/sumtot
         if( yRnd(16).gt.bound(k-1) .and. yRnd(16).lt.bound(k)  ) then
            ifound=i; jfound=j;
            goto 1313
         endif
      enddo
   enddo
1313 continue


if( IsAZDecay(DecayMode1) .or. IsAPhoton(DecayMode1) ) then
!if pp collider
  if(Collider.eq.1)then
    id(1:2) = (/ifound,jfound/)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    if(IsNaN(me2).or.(me2.eq.0d0))return
    if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2) * PreFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol
!if e+ e- collider
  else if(Collider.eq.0)then
    ifound=0
    jfound=0
    id(2)=convertLHE(ElM_)
    id(1)=-id(2)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    if(IsNaN(me2).or.(me2.eq.0d0))return
    if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
    LO_Res_Unpol = me2 * PreFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol
  endif

elseif( IsAWDecay(DecayMode1) ) then
!pp>WH
    id(1:2) = (/ifound,jfound/)
    if( ((id(1).eq.convertLHE(AUp_).or.id(1).eq.convertLHE(AChm_)) .and. &
     (id(2).eq.convertLHE(Dn_) .or. id(2).eq.convertLHE(Str_) .or. id(2).eq.convertLHE(Bot_))) .or. &
    ((id(2).eq.convertLHE(AUp_).or. id(2).eq.convertLHE(AChm_)) .and. &
     (id(1).eq.convertLHE(Dn_) .or. id(1).eq.convertLHE(Str_) .or. id(1).eq.convertLHE(Bot_)))   )then
      id(3)=-id(3)
      id(4)=-id(4)
      id(6)=-id(6)
      id(7)=-id(7)
      helicity(6)=sign(1d0,-dble(id(6)))
      helicity(7)=-helicity(6)
    endif
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    if(IsNaN(me2).or.(me2.eq.0d0))return
    if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2) * PreFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol
endif



  CS_max = csmax(ifound,jfound)
  if( EvalUnWeighted_VHiggs.gt. CS_max) then
    write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_VHiggs, CS_max
    write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_VHiggs, CS_max
    AlertCounter = AlertCounter + 1
    Res = 0d0
  elseif( EvalUnWeighted_VHiggs .gt. yRnd(17)*CS_max ) then
    call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,useAonshell=(IsAPhoton(DecayMode1)))
    do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
    enddo
    AccepCounter = AccepCounter + 1
    cyRnd(1)=yRnd(9)
    cyRnd(2)=yRnd(8)
    if(.not.IsAPhoton(DecayMode1)) then
      call EvalPhasespace_VDecay(MomExt(1:4,4),inv_mass(4),getMass(convertLHEreverse(id(6))),getMass(convertLHEreverse(id(7))),cyRnd(1:2),MomExt(1:4,6:7),PSWgt2)
    endif
    if(Collider.eq.1)then
      call boost2Lab(eta1,eta2,9,MomExt(1:4,1:9))
    endif
    do i=6,7
      inv_mass(i)=dsqrt(dabs(MomExt(1:4,i).dot.MomExt(1:4,i)))
    enddo
    call WriteOutEvent_VHiggs(id,helicity,MomExt,inv_mass,EventWeight=1d0)
  else
    RejeCounter = RejeCounter + 1
  endif


ELSE! NOT GENEVT

if( IsAZDecay(DecayMode1) .or. IsAPhoton(DecayMode1) ) then
!if pp collider
  if(Collider.eq.1)then
  do i = -5,5
    j = -i
    id(1:2) = (/i,j/)
    if (abs(i).ne.0)then
      call EvalAmp_VHiggs(id,helicity,MomExt,me2)
      if(IsNaN(me2).or.(me2.eq.0d0))return
      if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
    else
      me2=0d0
    endif
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac
    EvalUnWeighted_VHiggs = EvalUnWeighted_VHiggs+LO_Res_Unpol
    RES(i,j) = LO_Res_Unpol
    if (LO_Res_Unpol.gt.csmax(i,j)) then
      csmax(i,j) = LO_Res_Unpol
    endif
  enddo
!if e+ e- collider
  else if(Collider.eq.0)then
    id(2)=convertLHE(ElM_)
    id(1)=-id(2)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    if(IsNaN(me2).or.(me2.eq.0d0))return
    if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
    LO_Res_Unpol = me2 * PreFac
    EvalUnWeighted_VHiggs = EvalUnWeighted_VHiggs + LO_Res_Unpol
    RES(0,0) = LO_Res_Unpol
    if (LO_Res_Unpol.gt.csmax(0,0)) then
      csmax(0,0) = LO_Res_Unpol
    endif
  endif
elseif( IsAWDecay(DecayMode1) ) then
!pp>WH
  do i = -5,5
  do j = -5,5
     id2=id
     id2(1:2) = (/i,j/)
     if    ( ((id2(1).eq.convertLHE(Up_).or.id2(1).eq.convertLHE(Chm_)) .and. &
      (id2(2).eq.convertLHE(ADn_) .or. id2(2).eq.convertLHE(AStr_) .or. id2(2).eq.convertLHE(ABot_))) .or. &
     ((id2(2).eq.convertLHE(Up_).or.id2(2).eq.convertLHE(Chm_)) .and. &
      (id2(1).eq.convertLHE(ADn_) .or. id2(1).eq.convertLHE(AStr_) .or. id2(1).eq.convertLHE(ABot_)))   )then
           helicity(6)=sign(1d0,-dble(id2(6)))
           helicity(7)=-helicity(6)
           call EvalAmp_VHiggs(id2,helicity,MomExt,me2)
           if(IsNaN(me2).or.(me2.eq.0d0))return
           if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
     elseif( ((id2(1).eq.convertLHE(AUp_).or.id2(1).eq.convertLHE(AChm_)) .and. &
      (id2(2).eq.convertLHE(Dn_) .or. id2(2).eq.convertLHE(Str_) .or. id2(2).eq.convertLHE(Bot_))) .or. &
     ((id2(2).eq.convertLHE(AUp_).or. id2(2).eq.convertLHE(AChm_)) .and. &
      (id2(1).eq.convertLHE(Dn_) .or. id2(1).eq.convertLHE(Str_) .or. id2(1).eq.convertLHE(Bot_)))   )then
           id2(3)=-id2(3)
           id2(4)=-id2(4)
           id2(6)=-id2(6)
           id2(7)=-id2(7)
           helicity(6)=sign(1d0,-dble(id2(6)))
           helicity(7)=-helicity(6)
           call EvalAmp_VHiggs(id2,helicity,MomExt,me2)
           if(IsNaN(me2).or.(me2.eq.0d0))return
           if(H_DK.eqv..false.)me2=me2*(M_Reso*Ga_Reso)**2!remove erroneous H propagator with stable H in mod_VHiggs.F90
     else
           me2=0d0
     endif
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac
    EvalUnWeighted_VHiggs = EvalUnWeighted_VHiggs+LO_Res_Unpol
    RES(i,j) = LO_Res_Unpol
    if (LO_Res_Unpol.gt.csmax(i,j)) then
      csmax(i,j) = LO_Res_Unpol
    endif
  enddo
  enddo
endif


ENDIF! GENEVT


 RETURN
 end Function EvalUnWeighted_VHiggs



! with stable H, MATRIXELEMENT0() in mod_VHiggs.H90 mistakenly multiples an H propagator.
! This propagator is calculated here to help remove it, without changing MATRIXELEMENT0().
!function Remove_HPropagator(s)
!use ModParameters
!implicit none
!   Remove_HPropagator = (M_Reso*Ga_Reso)**2

!   return Remove_HPropagator

!end function







 FUNCTION EvalWeighted_tautau(yRnd,VgsWgt)
 use ModKinematics
 use ModParameters
 use ModGraviton
 use ModHiggs
 use ModZprime
 use ModMisc
#if compiler==1
 use ifport
#endif
 implicit none
 real(8) :: EvalWeighted_tautau,LO_Res_Unpol,yRnd(1:22),VgsWgt
 real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac1,PDFFac2
 real(8) :: pdf(-6:6,1:2)
 integer :: NBin(1:NumHistograms),NHisto,i,MY_IDUP(1:9), ICOLUP(1:2,1:9)
 real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
 logical :: applyPSCut
 real(8) :: pHiggs(1:4),Mom(1:4,1:12)



    EvalWeighted_tautau = 0d0
    if( OffShellReson ) then
      call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    else
!       call PDFMapping(11,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
      call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
    endif
    EvalCounter = EvalCounter+1

    MY_IDUP(1:2) = (/ElP_,MuM_/)
    pHiggs(2:4) = (/110d0, -60d0,95d0  /) ! What is this?
    pHiggs(1) = dsqrt( M_Reso**2 + pHiggs(2)**2+ pHiggs(3)**2+ pHiggs(4)**2 )
    call EvalPhasespace_tautau(yRnd(1:12),pHiggs,MY_IDUP,Mom,PSWgt)


    call Kinematics_Htautau(Mom,applyPSCut,NBin)
    if( applyPSCut ) then
      EvalWeighted_tautau = 0d0
      return
    endif


   call SetRunningScales( (/ pHiggs(1:4),Mom_Not_a_particle(1:4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) )
   call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)

   PDFFac = pdf(0,1) * pdf(0,2)

!       call matrix element here
      LO_Res_Unpol = 1d0


      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
      PreFac = 2d0 * hbarc2XsecUnit * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalWeighted_tautau = LO_Res_Unpol * PreFac

      do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),EvalWeighted_tautau*VgsWgt)
      enddo



RETURN
END FUNCTION

FUNCTION EvalUnWeighted_DecayToTauTau(yRnd,genEvt,Ehat,Res,AcceptedEvent,MY_IDUP,ICOLUP)
use ModKinematics
use ModParameters
use ModMisc
use ModHiggs
#if compiler==1
use ifport
#endif
implicit none
real(8) :: EvalUnWeighted_DecayToTauTau,Ehat,Res,AcceptedEvent(:,:),DKRnd
real(8) :: yRnd(:),Mom(1:4,1:13),pHiggs(1:4),PSWgt,PreFac,LO_Res_Unpol,CS_max,m1ffwgt,m2ffwgt,FinalStateWeight
integer :: MY_IDUP(:),ICOLUP(:,:),NBin(1:NumHistograms),NHisto=1,DK_IDUP(1:6),DK_ICOLUP(1:2,3:6)
logical :: applyPSCut,genEvt
integer, parameter :: inLeft=1, inRight=2, Hig=3, tauP=4, tauM=5, Wp=6, Wm=7,   nu=8, nubar_tau=9, lepP=10,   lepM=11, nu_tau=12, nubar=13
include 'vegas_common.f'
include 'csmaxvalue.f'
EvalUnWeighted_DecayToTauTau = 0d0
m1ffwgt=1d0
m2ffwgt=1d0


  ICOLUP(1:2,tauP) = (/000,000/)
  ICOLUP(1:2,tauM) = (/000,000/)
  MY_IDUP(tauP) = TaP_;
  MY_IDUP(tauM) = TaM_;
!   if( RandomizeVVdecays ) then
!      call random_number(DKRnd)
!      if( DKRnd.lt.0.5d0 ) call swapi(DecayMode1,DecayMode2)
!   endif
  call VVBranchings(DK_IDUP(1:6),DK_ICOLUP(1:2,3:6),FinalStateWeight,800)
  MY_IDUP(nu_tau)   = NuT_;        ICOLUP(1:2,nu_tau)   = (/000,000/)
  MY_IDUP(Wp)       = DK_IDUP(1);  ICOLUP(1:2,Wp)       = (/000,000/)
  MY_IDUP(lepP)     = DK_IDUP(3);  ICOLUP(1:2,lepP)     = DK_ICOLUP(1:2,3)
  MY_IDUP(nu)       = DK_IDUP(4);  ICOLUP(1:2,nu)       = DK_ICOLUP(1:2,4)
  MY_IDUP(nubar_tau)= ANuT_;       ICOLUP(1:2,nubar_tau)= (/000,000/)
  MY_IDUP(Wm)       = DK_IDUP(2);  ICOLUP(1:2,Wm)       = (/000,000/)
  MY_IDUP(lepM)     = DK_IDUP(6);  ICOLUP(1:2,lepM)     = DK_ICOLUP(1:2,6)
  MY_IDUP(nubar)    = DK_IDUP(5);  ICOLUP(1:2,nubar)    = DK_ICOLUP(1:2,5)

  pHiggs(1:4) = (/Ehat,0d0,0d0,0d0/)
  if( TauDecays.eq.0 ) then
     call genps(2,Ehat,yRnd(1:2),(/m_tau,m_tau/),Mom(1:4,tauP:tauM),PSWgt)
     PSWgt = PSWgt/(2d0*Pi)**2/(4d0*Pi)
  else
     call EvalPhasespace_tautau(yRnd(1:12),pHiggs,(/DK_IDUP(3),DK_IDUP(6)/),Mom,PSWgt)
  endif
  call Kinematics_Htautau(Mom,applyPSCut,NBin)
  if( applyPSCut ) then
      EvalUnWeighted_DecayToTauTau = 0d0
      return
  endif
  PreFac = hbarc2XsecUnit * PSWgt * FinalStateWeight

  call SetRunningScales( (/ pHiggs(1:4),Mom_Not_a_particle(1:4),Mom_Not_a_particle(1:4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) ) ! Call anyway

  if( TauDecays.eq.0 ) then
!     if (genevt) then
!        call printMom(Mom(1:4,tauP:tauM))
!     endif
     call EvalAmp_H_FF(Mom(1:4,tauP:tauM),m_tau,LO_Res_Unpol)
  else
     call EvalAmp_H_TT_decay((/Mom(1:4,lepM),Mom(1:4,nubar),Mom(1:4,nu_tau),Mom(1:4,nu),Mom(1:4,lepP),Mom(1:4,nubar_tau)/),m_tau,ga_tau,LO_Res_Unpol)
  endif

  if( TauDecays.ne.0 ) then
     call ShiftMass(Mom(1:4,LepP),Mom(1:4,Wp)-Mom(1:4,LepP), GetMass(MY_IDUP(LepP)),0d0,  AcceptedEvent(1:4,LepP),AcceptedEvent(1:4,Nu), MassWeight=m1ffwgt )
     call ShiftMass(Mom(1:4,LepM),Mom(1:4,Wm)-Mom(1:4,LepM), GetMass(MY_IDUP(LepM)),0d0,  AcceptedEvent(1:4,LepM),AcceptedEvent(1:4,Nubar), MassWeight=m2ffwgt )
     PreFac = PreFac*m1ffwgt*m2ffwgt
  endif

  EvalUnWeighted_DecayToTauTau = LO_Res_Unpol * PreFac

IF( GENEVT ) THEN

      CS_max = csmax(0,0)

      AcceptedEvent(:,:) = Mom(:,:)

      if( EvalUnWeighted_DecayToTauTau .gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_DecayToTauTau, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_DecayToTauTau, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0
      elseif( EvalUnWeighted_DecayToTauTau .gt. yRnd(14)*CS_max ) then
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1

         Res = 1d0
      else
         RejeCounter = RejeCounter + 1
         Res = 0d0
      endif
      EvalCounter = EvalCounter + 1


ELSE! NOT GENEVT

      Res = EvalUnWeighted_DecayToTauTau
      if( EvalUnWeighted_DecayToTauTau.gt.csmax(0,0) ) then
          csmax(0,0) = EvalUnWeighted_DecayToTauTau
      endif

ENDIF! GENEVT



RETURN
END FUNCTION


 FUNCTION EvalOnlyPS(yRnd,VgsWgt)    ! this is a function which generates only the PS events
 use ModKinematics
 use ModParameters
 use ModMisc
#if compiler==1
 use ifport
#endif
 implicit none
 real(8) :: EvalOnlyPS,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
 real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac1,PDFFac2
 real(8) :: pdf(-6:6,1:2)
 integer :: NBin(1:NumHistograms),NHisto,i,MY_IDUP(1:9), ICOLUP(1:2,1:9),xBin(1:4)
 real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,FinalStateWeight
 real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomDK_massless(1:4,1:4)
 real(8) :: EZ,  EZ1, EZ2, pz, xmax, xx1, xx2
 real(8) :: pz12, MomExt_f(1:4,1:4), MomDK_f(1:4,1:4)
 real(8) :: MZ1,MZ2,ML1,ML2,ML3,ML4,EZ_max,dr, MG, yz1,yz2
 real(8) :: offzchannel
 logical :: applyPSCut
 include 'csmaxvalue.f'

    EvalOnlyPS = 0d0
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
    call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9),FinalStateWeight)
    MY_IDUP(1:3) = 0
    yz1 = yRnd(10)
    yz2 = yRnd(11)
    offzchannel = yRnd(12) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * ( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )

            EZ_max = EHat - MZ1*0.99d0
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )

        else   !M_Reso.le.2d0*M_V_ps
            if (offzchannel.le.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
                MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
                MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
                sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/(  &
                1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )     &
                + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
            elseif(offzchannel.gt.0.5d0) then
                EZ_max = EHat
                dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
                MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
                MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
                sHatJacobi = sHatJacobi * dr/(Ga_V_ps*M_V_ps) * 1d0/( &
                1d0/((MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )    &
                + 1d0/((MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 ) )
                sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
            endif
       endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..true.)) then
        MZ1 = getMass(MY_IDUP(4))
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat - MZ1*0.99
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ2 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz2-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ2**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
        else
            MZ2 = abs(EHat - MZ1*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ1*0.999)**2
        endif

  elseif((OffShellV1.eqv..true.).and.(OffShellV2.eqv..false.)) then
        MZ2 = getMass(MY_IDUP(5))
        if(M_Reso.gt.2d0*M_V_ps) then
            EZ_max = EHat - MZ2*0.99
            dr = datan((EZ_max**2-M_V_ps**2)/(Ga_V_ps*M_V_ps)) + datan(M_V_ps/Ga_V_ps)
            MZ1 = dsqrt( M_V_ps*Ga_V_ps * dtan(dr*yz1-datan(M_V_ps/Ga_V_ps)) + M_V_ps**2 )
            sHatJacobi = sHatJacobi*dr/(Ga_V_ps*M_V_ps)*( (MZ1**2 - M_V_ps**2)**2 + M_V_ps**2*Ga_V_ps**2 )
         else
            MZ1 = abs(EHat - MZ2*0.999999999999999d0)*dsqrt(abs(dble(yz2)))
            sHatJacobi = sHatJacobi *(EHat - MZ2*0.999)**2
        endif

  elseif((OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.)) then
        MZ1 = getMass(MY_IDUP(4))
        MZ2 = getMass(MY_IDUP(5))
  endif

    if( EHat.lt.MZ1+MZ2 ) then
      EvalOnlyPS = 0d0
      return
    endif

    call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
    call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
    if( .not.IsAPhoton(DecayMode1) .and. .not.IsAPhoton(DecayMode2) ) then ! decay ZZ's and WW's
        ML1 = getMass(MY_IDUP(7))
        ML2 = getMass(MY_IDUP(6))
        ML3 = getMass(MY_IDUP(9))
        ML4 = getMass(MY_IDUP(8))
        if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
            EvalOnlyPS = 0d0
            return
        endif
        call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
        call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,ML3,ML4,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
        PSWgt = PSWgt * PSWgt2*PSWgt3
        if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then! introduce this momentum flip to allow proper mapping of integrand with Z-poles at MZ2=(p2+p3)^2 and MZ2=(p1+p4)^2
            if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK(1:4,1),MomDK(1:4,3) )
            PSWgt = PSWgt * 2d0
        endif
        if( (includeInterference.eqv..false.) .and. (OffShellV1.eqv..false.).and.(OffShellV2.eqv..false.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
            PSWgt = PSWgt * 1d0/2d0
        endif
    elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        ML1=0d0; ML2=0d0; ML3=0d0; ML4=0d0
        MomDK(1:4,1) = MomExt(1:4,3)
        MomDK(1:4,2) = 0d0
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
    elseif( .not.IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        ML1 = getMass(MY_IDUP(7))
        ML2 = getMass(MY_IDUP(6))
        ML3=0d0; ML4=0d0
        if( (MZ1.lt.ML1+ML2) ) then
            EvalOnlyPS = 0d0
            return
        endif
        call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
        PSWgt = PSWgt * PSWgt2
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
    endif



    if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
        call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    else
        call Error("This should not happen anymore")
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalOnlyPS = 0d0
      return
    endif


    EvalOnlyPS = PSWgt * FinalStateWeight

    if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2))   ) then
          call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=1d0)
    else
          call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=1d0)
    endif


    do NHisto=1,NumHistograms-3
       call intoHisto(NHisto,NBin(NHisto),1d0)
    enddo


RETURN
END FUNCTION



END MODULE ModCrossSection







