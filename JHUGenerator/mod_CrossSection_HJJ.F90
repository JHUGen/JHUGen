MODULE ModCrossSection_HJJ
implicit none
integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

 CONTAINS




! since the me2(:,:) array is defined from -5..+5, the iPart_sel,jPart_sel have to follow the LHE numbering convention
FUNCTION EvalWeighted_HJJ_fulldecay(yRnd,VgsWgt)
use ModKinematics
use ModParameters
use ModHiggsjj
use ModHiggs
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
integer,parameter :: mxpart=14 ! this has to match the MCFM parameter
real(8) :: yRnd(1:18),VgsWgt, EvalWeighted_HJJ_fulldecay
real(8) :: pdf(-6:6,1:2)           ,me2(-5:5,-5:5),me2_tmpzz(-5:5,-5:5),me2_tmpww(-5:5,-5:5),me2_hdk,me2_prop
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:10),PSWgt
real(8) :: p_MCFM(mxpart,1:4),msq_MCFM(-5:5,-5:5)
complex(8) :: HZZcoupl(1:32),HWWcoupl(1:32)
integer :: MY_IDUP(1:10),ICOLUP(1:2,1:10),NBin(1:NumHistograms),NHisto
integer :: iPartChannel,PartChannelAvg,NumPartonicChannels,ijSel(1:121,1:3),i,j,iflip
real(8) :: LO_Res_Unpol, PreFac,VegasWeighted_HJJ_fulldecay,xRnd
logical :: applyPSCut
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10
real(8) :: s13,s14,s15,s16,s23,s24,s25,s26,s34,s35,s36,s45,s46,s56,s78,s910,s710,s89
include 'vegas_common.f'
EvalWeighted_HJJ_fulldecay = 0d0



   NumPartonicChannels = 121

! selecting only u-d channels for checks
NumPartonicChannels = 2
! after removing also repair get_GENchannelHash !!

   iPartChannel = int(yRnd(18) * (NumPartonicChannels)) +1 ! this runs from 1..121
   call get_GENchannelHash(ijSel)
   iPart_sel = ijSel(iPartChannel,1)
   jPart_sel = ijSel(iPartChannel,2)
   PartChannelAvg = NumPartonicChannels
   if( unweighted .and. .not.warmup .and.  sum(AccepCounter_part(:,:)) .eq. sum(RequEvents(:,:)) ) then
      stopvegas=.true.
   endif
   if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return


   call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   call EvalPhasespace_VBF_H4f(yRnd(3),yRnd(4:17),EHat,MomExt(1:4,1:10),PSWgt)

!       call genps(6,EHat,yRnd(3:16),(/0d0,0d0,0d0,0d0,0d0,0d0/),MomExt(1:4,3:8),PSWgt)
!       MomExt(1:4,1)=(/Ehat,0d0,0d0,+Ehat/)/2d0
!       MomExt(1:4,2)=(/Ehat,0d0,0d0,-Ehat/)/2d0
!       MomExt(1:4,10)=MomExt(1:4,8)
!       MomExt(1:4,9) =MomExt(1:4,7)
!       MomExt(1:4,8) =MomExt(1:4,6)
!       MomExt(1:4,7) =MomExt(1:4,5)
!       MomExt(1:4,5) = MomExt(1:4,7)+MomExt(1:4,8)
!       MomExt(1:4,6) = MomExt(1:4,9)+MomExt(1:4,10)
!       PSWgt = PSWgt * (2d0*Pi)**(4-(6)*3) * (4d0*Pi)**((6)-1)

!
!       EvalWeighted_HJJ_fulldecay=PSWgt*sHatJacobi  * ( MomExt(1:4,3).dot.MomExt(1:4,7) ) * ( MomExt(1:4,4).dot.MomExt(1:4,10) ) * ( MomExt(1:4,8).dot.MomExt(1:4,9) ) * ( MomExt(1:4,7).dot.MomExt(1:4,10) ) / EHat**8
!       return

!    call EvalPhasespace_VBF_NEW2(yRnd(17),yRnd(3:7),EHat,MomExt(1:4,1:10),PSWgt)! stable Higgs
!       call genps(3,EHat,yRnd(3:7),(/0d0,0d0,M_Reso/),MomExt(1:4,3:5),PSWgt)
!       PSWgt = PSWgt * (2d0*Pi)**(4-(3)*3) * (4d0*Pi)**((3)-1)
   call boost2Lab(eta1,eta2,10,MomExt(1:4,1:10))
   PSWgt = PSWgt * (100d0)**8 * PartChannelAvg ! adjust PSWgt for GeV units of MCFM mat.el.


   call Kinematics_HVBF_fulldecay(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.lt.1d-12 ) return
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)
   EvalCounter = EvalCounter+1


   MY_IDUP(1:10) = (/Up_,Up_,Up_,Up_, Z0_,Z0_, ElM_,ElP_,MuM_,MuP_/)
   ICOLUP(1:2,1) = (/501,0/)
   ICOLUP(1:2,2) = (/0,502/)
   ICOLUP(1:2,3) = (/501,0/)
   ICOLUP(1:2,4) = (/0,502/)
   ICOLUP(1:2,5:10) = 0

   call convert_to_MCFM(-MomExt(1:4,inTop)*100d0, p_MCFM(1,1:4))
   call convert_to_MCFM(-MomExt(1:4,inBot)*100d0, p_MCFM(2,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep1P)*100d0, p_MCFM(3,1:4))! check f fbar assignment
   call convert_to_MCFM(+MomExt(1:4,Lep1M)*100d0, p_MCFM(4,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2P)*100d0, p_MCFM(5,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2M)*100d0, p_MCFM(6,1:4))
   call convert_to_MCFM(+MomExt(1:4,outTop)*100d0,p_MCFM(7,1:4))
   call convert_to_MCFM(+MomExt(1:4,outBot)*100d0,p_MCFM(8,1:4))

   HZZcoupl(1) = 1d0  !ghz1
   HZZcoupl(2) = ghz2
   HZZcoupl(3) = ghz3
   HZZcoupl(4) = ghz4
   HZZcoupl(5:) = (0d0,0d0)

   HWWcoupl(1) = 1d0! HZZcoupl(1)!  ghw1  ! this is actually wrong
   HWWcoupl(2) = ghw2
   HWWcoupl(3) = ghw3
   HWWcoupl(4) = ghw4
   HWWcoupl(5:) = (0d0,0d0)

   do i=1,10
   print *,MomExt(1:4,i),(MomExt(1:4,i)).dot.(MomExt(1:4,i))
   enddo
    print *, "---"
   do i=1,8
   print *,p_MCFM(i,1:4)
   enddo
!    pause



 msq_MCFM(:,:) = 0d0
#if linkMELA==1
 call qq_ZZqq(p_MCFM,msq_MCFM,HZZcoupl,HWWcoupl,Lambda*100d0,Lambda_Q*100d0,(/Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4/)*100d0)!  q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8)
#else
 print *, "To use this process, please set linkMELA=Yes in the makefile and recompile."
 print *, "You will also need to have a compiled JHUGenMELA in the directory specified by JHUGenMELADir in the makefile."
 stop 1
#endif
 ! large overhead here because MCFM computes all partonic channels and we only use msq_MCFM(iPart_sel,jPart_sel)  !


!   CHECKS:
   print *,"msq_MCFM:"
   do j=-5,5
   print *, msq_MCFM(-5,j), msq_MCFM(-4,j), msq_MCFM(-3,j), msq_MCFM(-2,j), msq_MCFM(-1,j), msq_MCFM(0,j), msq_MCFM(1,j), msq_MCFM(2,j), msq_MCFM(3,j), msq_MCFM(4,j), msq_MCFM(5,j)
   enddo
   print *,""

   me2(:,:)=0d0
   do i=-5,5
   do j=-5,5
   me2_tmpzz(:,:)=0d0
   me2_tmpww(:,:)=0d0
   call EvalAmp_WBFH_UnSymm_SA_Select( (/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)+MomExt(1:4,6)/),i,j,.true.,iflip,me2_tmpzz) ! calling on-shell VBF with stable Higgs
   if (iflip.eq.2) then
   call swap(me2_tmpzz(i,j),me2_tmpzz(j,i))
   endif
   call EvalAmp_WBFH_UnSymm_SA_Select( (/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)+MomExt(1:4,6)/),i,j,.true.,iflip,me2_tmpww) ! calling on-shell VBF with stable Higgs
   if (iflip.eq.1) then
   call swap(me2_tmpww(i,j),me2_tmpww(j,i))
   endif
   me2 = me2 + me2_tmpzz + me2_tmpww
   enddo
   enddo
!   call EvalAmp_WBFH_UnSymm_SA(MomExt(1:4,1:5),me2)
! !   msq_MCFM(:,:) = me2(:,:)

   call EvalAmp_H_VV( (/MomExt(1:4,5)+MomExt(1:4,6),(/0d0,0d0,0d0,0d0/),MomExt(1:4,7),MomExt(1:4,8),MomExt(1:4,9),MomExt(1:4,10)/),(/ElM_,ElP_,MuM_,MuP_/),me2_hdk)                     ! adding higgs decay
   print *,"me2_hdk:",me2_hdk

   me2_prop = cdabs(1d0/( ((MomExt(1:4,5)+MomExt(1:4,6)).dot.(MomExt(1:4,5)+MomExt(1:4,6))) - m_Reso**2 + (0d0,1d0)*m_Reso*Ga_Reso ))**2
   print *,"me2_prop:",me2_prop

   me2(:,:) = me2(:,:) * me2_hdk * me2_prop !  adding higgs propagator
   print *,"me2:"
   do j=-5,5
   print *, me2(-5,j), me2(-4,j), me2(-3,j), me2(-2,j), me2(-1,j), me2(0,j), me2(1,j), me2(2,j), me2(3,j), me2(4,j), me2(5,j)
   enddo
   print *,""
   !print *, "rat", msq_MCFM(j,i)/me2(i,j)
   pause


   LO_Res_Unpol = msq_MCFM(iPart_sel,jPart_sel)  *  pdf(LHA2M_pdf(iPart_sel),1) * pdf(LHA2M_pdf(jPart_sel),2)

   PreFac = fbGeV2 * FluxFac * PSWgt * sHatJacobi
   EvalWeighted_HJJ_fulldecay = LO_Res_Unpol * PreFac
   VegasWeighted_HJJ_fulldecay = EvalWeighted_HJJ_fulldecay*VgsWgt



   if( unweighted ) then

     if( warmup ) then

       CrossSec(iPart_sel,jPart_sel) = CrossSec(iPart_sel,jPart_sel) + VegasWeighted_HJJ_fulldecay
       CrossSecMax(iPart_sel,jPart_sel) = max(CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_HJJ_fulldecay)

     else! not warmup

       call random_number(xRnd)
       if( VegasWeighted_HJJ_fulldecay.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",VegasWeighted_HJJ_fulldecay, CrossSecMax(iPart_sel,jPart_sel)
         write(io_stdout, "(2X,A,1PE13.6,1PE13.6,1PE13.6,I3,I3)") "CrossSecMax is too small.",VegasWeighted_HJJ_fulldecay, CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_HJJ_fulldecay/CrossSecMax(iPart_sel,jPart_sel),iPart_sel,jPart_sel
         AlertCounter = AlertCounter + 1
       elseif( VegasWeighted_HJJ_fulldecay .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPart_sel,jPart_sel) = AccepCounter_part(iPart_sel,jPart_sel) + 1
!         call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=1d0)
         call WriteOutEvent_HJJ_fulldecay(MomExt,MY_IDUP,ICOLUP)
         do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
       endif

     endif! warmup

   else! weighted

      if( VegasWeighted_HJJ_fulldecay.ne.0d0 ) then
        AccepCounter=AccepCounter+1
        if( writeWeightedLHE .and. (.not. warmup) ) then
            call WriteOutEvent_HJJ_fulldecay(MomExt,MY_IDUP,ICOLUP,EventWeight=VegasWeighted_HJJ_fulldecay)
        endif
        do NHisto=1,9 !NumHistograms
          call intoHisto(NHisto,NBin(NHisto),VegasWeighted_HJJ_fulldecay)
        enddo



! compute all invariants to check PS performance

 s13=get_MInv(MomExt(1:4,1)-MomExt(1:4,3))
 s14=get_MInv(MomExt(1:4,1)-MomExt(1:4,4))
 s15=get_MInv(MomExt(1:4,1)-MomExt(1:4,5))
 s16=get_MInv(MomExt(1:4,1)-MomExt(1:4,6))

 s23=get_MInv(MomExt(1:4,2)-MomExt(1:4,3))
 s24=get_MInv(MomExt(1:4,2)-MomExt(1:4,4))
 s25=get_MInv(MomExt(1:4,2)-MomExt(1:4,5))
 s26=get_MInv(MomExt(1:4,2)-MomExt(1:4,6))

 s34=get_MInv(MomExt(1:4,3)+MomExt(1:4,4))
 s35=get_MInv(MomExt(1:4,3)+MomExt(1:4,5))
 s36=get_MInv(MomExt(1:4,3)+MomExt(1:4,6))

 s45=get_MInv(MomExt(1:4,4)+MomExt(1:4,5))
 s46=get_MInv(MomExt(1:4,4)+MomExt(1:4,6))

 s56=get_MInv(MomExt(1:4,5)+MomExt(1:4,6))

 s78=get_MInv(MomExt(1:4,7)+MomExt(1:4,8))
 s910=get_MInv(MomExt(1:4,9)+MomExt(1:4,10))
 s710=get_MInv(MomExt(1:4,7)+MomExt(1:4,10))
 s89=get_MInv(MomExt(1:4,8)+MomExt(1:4,9))

 LO_Res_Unpol = LO_Res_Unpol * VgsWgt

 call intoHisto(10,WhichBin(NHisto,dsqrt(s13)),LO_Res_Unpol)
 call intoHisto(11,WhichBin(NHisto,dsqrt(s13)),1d0/PSWgt)

 call intoHisto(12,WhichBin(NHisto,dsqrt(s14)),LO_Res_Unpol)
 call intoHisto(13,WhichBin(NHisto,dsqrt(s14)),1d0/PSWgt)

 call intoHisto(14,WhichBin(NHisto,dsqrt(s15)),LO_Res_Unpol)
 call intoHisto(15,WhichBin(NHisto,dsqrt(s15)),1d0/PSWgt)

 call intoHisto(16,WhichBin(NHisto,dsqrt(s16)),LO_Res_Unpol)
 call intoHisto(17,WhichBin(NHisto,dsqrt(s16)),1d0/PSWgt)

 call intoHisto(18,WhichBin(NHisto,dsqrt(s23)),LO_Res_Unpol)
 call intoHisto(19,WhichBin(NHisto,dsqrt(s23)),1d0/PSWgt)

 call intoHisto(20,WhichBin(NHisto,dsqrt(s24)),LO_Res_Unpol)
 call intoHisto(21,WhichBin(NHisto,dsqrt(s24)),1d0/PSWgt)

 call intoHisto(22,WhichBin(NHisto,dsqrt(s25)),LO_Res_Unpol)
 call intoHisto(23,WhichBin(NHisto,dsqrt(s25)),1d0/PSWgt)

 call intoHisto(24,WhichBin(NHisto,dsqrt(s26)),LO_Res_Unpol)
 call intoHisto(25,WhichBin(NHisto,dsqrt(s26)),1d0/PSWgt)

 call intoHisto(26,WhichBin(NHisto,dsqrt(s34)),LO_Res_Unpol)
 call intoHisto(27,WhichBin(NHisto,dsqrt(s34)),1d0/PSWgt)

 call intoHisto(28,WhichBin(NHisto,dsqrt(s35)),LO_Res_Unpol)
 call intoHisto(29,WhichBin(NHisto,dsqrt(s35)),1d0/PSWgt)

 call intoHisto(30,WhichBin(NHisto,dsqrt(s36)),LO_Res_Unpol)
 call intoHisto(31,WhichBin(NHisto,dsqrt(s36)),1d0/PSWgt)

 call intoHisto(32,WhichBin(NHisto,dsqrt(s45)),LO_Res_Unpol)
 call intoHisto(33,WhichBin(NHisto,dsqrt(s45)),1d0/PSWgt)

 call intoHisto(34,WhichBin(NHisto,dsqrt(s46)),LO_Res_Unpol)
 call intoHisto(35,WhichBin(NHisto,dsqrt(s46)),1d0/PSWgt)

 call intoHisto(36,WhichBin(NHisto,dsqrt(s56)),LO_Res_Unpol)
 call intoHisto(37,WhichBin(NHisto,dsqrt(s56)),1d0/PSWgt)

 call intoHisto(38,WhichBin(NHisto,dsqrt(s78)),LO_Res_Unpol)
 call intoHisto(39,WhichBin(NHisto,dsqrt(s78)),1d0/PSWgt)

 call intoHisto(40,WhichBin(NHisto,dsqrt(s910)),LO_Res_Unpol)
 call intoHisto(41,WhichBin(NHisto,dsqrt(s910)),1d0/PSWgt)

 call intoHisto(42,WhichBin(NHisto,dsqrt(s710)),LO_Res_Unpol)
 call intoHisto(43,WhichBin(NHisto,dsqrt(s710)),1d0/PSWgt)

 call intoHisto(44,WhichBin(NHisto,dsqrt(s89)),LO_Res_Unpol)
 call intoHisto(45,WhichBin(NHisto,dsqrt(s89)),1d0/PSWgt)



      endif

   endif! unweighted


RETURN
END FUNCTION





FUNCTION EvalUnWeighted_HJJ_fulldecay(yRnd,genEvt,iPartons,RES)
use ModKinematics
use ModParameters
use ModHiggsjj
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
integer,parameter :: mxpart=14 ! this has to match the MCFM parameter
real(8) :: yRnd(1:17),VgsWgt, EvalUnWeighted_HJJ_fulldecay
real(8) :: pdf(-6:6,1:2)   ,me2(-5:5,-5:5),RES(-5:5,-5:5)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:10),PSWgt,CS_Max
real(8) :: p_MCFM(mxpart,1:4),msq_MCFM(-5:5,-5:5)
complex(8) :: HZZcoupl(1:32),HWWcoupl(1:32)
integer :: i,j,MY_IDUP(1:10),ICOLUP(1:2,1:10),NBin(1:NumHistograms),NHisto,iPartons(1:2)
real(8) :: LO_Res_Unpol,LO_Res_Pol, PreFac,BWJacobi
logical :: applyPSCut,genEvt
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10
include 'csmaxvalue.f'
EvalUnWeighted_HJJ_fulldecay = 0d0


   call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    if( ehat.lt. max(m_reso+2*pTjetcut,VBF_4ml_minmax(1)) ) return ! 30GeV = pTcut on fwd jets
   if( ehat.lt. m_reso+2*pTjetcut ) return ! 30GeV = pTcut on fwd jets
!   if( ehat.lt.300d0*GeV+30d0*GeV ) return ! 30GeV = pTcut on fwd jets

   call EvalPhasespace_VBF_H4f(yRnd(17),yRnd(3:16),EHat,MomExt(1:4,1:10),PSWgt)
   call boost2Lab(eta1,eta2,11,MomExt(1:4,1:10))
   PSWgt = PSWgt * (100d0)**8

   call Kinematics_HVBF_fulldecay(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.lt.1d-12 ) return

   call SetRunningScales( (/ (MomExt(1:4,Lep1P)+MomExt(1:4,Lep1M)+MomExt(1:4,Lep2P)+MomExt(1:4,Lep2M)),MomExt(1:4,outTop),MomExt(1:4,outBot) /) , (/ Not_a_particle_,Up_,AUp_,Not_a_particle_ /) ) ! Implement like this for now, has to change with deterministic flavors!
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)


   MY_IDUP(1:10) = (/Up_,Up_,Up_,Up_, Z0_,Z0_, ElM_,ElP_,MuM_,MuP_/)
   ICOLUP(1:2,1) = (/501,0/)
   ICOLUP(1:2,2) = (/0,502/)
   ICOLUP(1:2,3) = (/501,0/)
   ICOLUP(1:2,4) = (/0,502/)
   ICOLUP(1:2,5:10) = 0


   call convert_to_MCFM(-MomExt(1:4,inTop)*100d0, p_MCFM(1,1:4))
   call convert_to_MCFM(-MomExt(1:4,inBot)*100d0, p_MCFM(2,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep1P)*100d0, p_MCFM(3,1:4))! check f fbar assignment
   call convert_to_MCFM(+MomExt(1:4,Lep1M)*100d0, p_MCFM(4,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2P)*100d0, p_MCFM(5,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2M)*100d0, p_MCFM(6,1:4))
   call convert_to_MCFM(+MomExt(1:4,outTop)*100d0,p_MCFM(7,1:4))
   call convert_to_MCFM(+MomExt(1:4,outBot)*100d0,p_MCFM(8,1:4))

   HZZcoupl(1) = (1d0,0d0)
   HWWcoupl(:) = HZZcoupl(:)
   msq_MCFM(:,:) = 0d0




IF( GENEVT ) THEN

#if linkMELA==1
          call qq_ZZqq(p_MCFM,msq_MCFM,HZZcoupl,HWWcoupl,Lambda*100d0,Lambda_Q*100d0,(/Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4/)*100d0)!  q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8)
#else
          print *, "To use this process, please set linkMELA=Yes in the makefile and recompile."
          print *, "You will also need to have a compiled JHUGenMELA in the directory specified by JHUGenMELADir in the makefile."
          stop 1
#endif
          LO_Res_Unpol = 0d0
!           do i = -5,5
!               do j = -5,5
do i = 1,1
do j = 2,2
                LO_Res_Unpol = LO_Res_Unpol + msq_MCFM(i,j) * pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2)
              enddo
          enddo
          PreFac = fbGeV2 * FluxFac * PSWgt * sHatJacobi
          EvalUnWeighted_HJJ_fulldecay = LO_Res_Unpol * PreFac

      CS_max = CSmax(iPartons(1),iPartons(2))

! print *, "check",iPartons(1:2)
! print *, "check",EvalUnWeighted_HJJ_fulldecay ,yRnd(16)*CS_max,CS_max;pause

      if( EvalUnWeighted_HJJ_fulldecay .gt. CS_max) then
         write(*,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_HJJ_fulldecay, CS_max
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_HJJ_fulldecay, CS_max
         AlertCounter = AlertCounter + 1
      elseif( EvalUnWeighted_HJJ_fulldecay .gt. yRnd(16)*CS_max ) then
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPartons(1),iPartons(2)) = AccepCounter_part(iPartons(1),iPartons(2))+1
         call WriteOutEvent_HJJ_fulldecay(MomExt,MY_IDUP,ICOLUP)
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
      endif
      EvalCounter = EvalCounter + 1


ELSE! NOT GENEVT


#if linkMELA==1
   call qq_ZZqq(p_MCFM,msq_MCFM,HZZcoupl,HWWcoupl,Lambda*100d0,Lambda_Q*100d0,(/Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4/)*100d0)!  q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8)
#else
   print *, "To use this process, please set linkMELA=Yes in the makefile and recompile."
   print *, "You will also need to have a compiled JHUGenMELA in the directory specified by JHUGenMELADir in the makefile."
   stop 1
#endif
   PreFac = fbGeV2 * FluxFac * PSWgt * sHatJacobi

   LO_Res_Unpol = 0d0
!    do i = -5,5
!       do j = -5,5
do i = 1,1
do j = 2,2
         LO_Res_Pol = msq_MCFM(i,j) * pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2) * PreFac

         RES(i,j) = LO_Res_Pol
         if( RES(i,j).gt.CSmax(i,j) ) then
           CSmax(i,j) = RES(i,j)
         endif

         LO_Res_Unpol = LO_Res_Unpol + LO_Res_Pol

      enddo
   enddo
   EvalUnWeighted_HJJ_fulldecay = LO_Res_Unpol



ENDIF! GENEVT


RETURN
END FUNCTION





 ! test VBF/HJJ function for process==-60/-61
FUNCTION EvalWeighted_HJJ_test(yRnd,VgsWgt)
use ModKinematics
use ModParameters
use ModHiggsjj
use ModMisc
#if compiler==1
use ifport
#endif
   implicit none
   real(8) :: yRnd(1:9),VgsWgt, EvalWeighted_HJJ_test
   real(8) :: VegasWeighted_HJJ
   real(8) :: pdf(-6:6,1:2)
   real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
   real(8) :: MomExt(1:4,1:5), PSWgt
   real(8) :: me2(-5:5,-5:5)
   real(8) :: me2_testhvv
   integer :: i,j,MY_IDUP(1:5),ICOLUP(1:2,1:5),NBin(1:NumHistograms),NHisto
   integer :: iPartChannel,PartChannelAvg,NumPartonicChannels,ijSel(1:121,1:3),flavor_tag
   real(8) :: LO_Res_Unpol, PreFac,xRnd,partonic_flip,outgoing_flip
   integer :: partbound(0:121), part_multiplicity, part_tracker
   real(8) :: part_detect
   logical :: applyPSCut
   integer, parameter :: ij_neg_offset=6, ij_max=+5
   integer, parameter :: ij_num=ij_max+ij_neg_offset
   integer :: rPart_sel, sPart_sel, iStore, jStore
   real(8) :: me2_swap(-5:5,-5:5)
   include 'vegas_common.f'

   me2(:,:) = 0d0
   me2_swap(:,:) = 0d0
   EvalWeighted_HJJ_test = 0d0
   VegasWeighted_HJJ = 0d0

!   print *, "Begin EvalWeighted_HJJ_test"

   ! Determine which partonic channel to generate
   if( Process.eq.60 ) then!  assuming everywhere that i>j  (apart from the LHE writeout)
      call get_VBFchannelHash_nosplit(ijSel,NumPartonicChannels)

      part_tracker=0; partbound(part_tracker)=0 ! These are here to make sure WW:ZZ ratio in the id phase space is (3*3):1 due to subdivision into each CKM partner
      do part_tracker = 1,NumPartonicChannels
         flavor_tag = ijSel(part_tracker,3)
         if(flavor_tag.eq.2) then ! ZZ-WW fusion
            partbound(part_tracker) = partbound(part_tracker-1) + 9
         elseif(flavor_tag.eq.0) then ! WW-only fusion
            partbound(part_tracker) = partbound(part_tracker-1) + 9
         else ! ZZ fusion
         partbound(part_tracker) = partbound(part_tracker-1) + 1
         endif
      enddo
      part_detect = int(yRnd(8) * partbound(NumPartonicChannels)) ! [0...268)
      do part_tracker = 1,NumPartonicChannels
         iPartChannel = part_tracker
         if(part_detect.ge.partbound(part_tracker-1) .and. part_detect.lt.partbound(part_tracker)) exit ! part_detect==268 also goes into the last bin (==68) this way.
      enddo

!     print *, "partbound: ",partbound
      PartChannelAvg = partbound(NumPartonicChannels)
   elseif( Process.eq.61 ) then
      call get_HJJchannelHash(ijSel)
      NumPartonicChannels = 77
      iPartChannel = int(yRnd(8) * (NumPartonicChannels)) +1 ! this runs from 1..77
      PartChannelAvg = NumPartonicChannels
   endif
   iStore = ijSel(iPartChannel,1)
   jStore = ijSel(iPartChannel,2)
   flavor_tag = ijSel(iPartChannel,3)

   if( unweighted .and. .not.warmup .and.  sum(AccepCounter_part(:,:)) .eq. sum(RequEvents(:,:)) ) stopvegas=.true.
   if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iStore,jStore) .ge. RequEvents(iStore,jStore))  ) return

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if (EHat.lt.M_Reso) return
   ! Outgoing ids~incoming ids up to swapping and flavor changes, and main.F90 protexts from errors due to swapping (ie. massive mJ_mJ-fixed scenario)
   call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,3),MomExt(1:4,4) /) , (/ Not_a_particle_,LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel),Not_a_particle_ /) )
   if ( Process.eq.61 ) call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)

   ! Assign the stored variables to the working ones now
   iPart_sel = iStore
   jPart_sel = jStore
   call random_number(partonic_flip)
   if( partonic_flip.gt.0.5d0 ) call swapi(iPart_sel,jPart_sel)

   if( Process.eq.60 ) then
      MY_IDUP(1:2)= (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel)/)

      if( MY_IDUP(1).gt.0 ) then ! quark
         ICOLUP(1:2,1) = (/501,000/)
      else! anti-quark
         ICOLUP(1:2,1) = (/000,501/)
      endif
      if( MY_IDUP(2).gt.0 ) then! quark
         ICOLUP(1:2,2) = (/502,000/)
      else! anti-quark
         ICOLUP(1:2,2) = (/000,502/)
      endif

      if( flavor_tag.eq.1 ) then ! ZZ H
         MY_IDUP(3:4)= (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel)/)
         ICOLUP(1:2,3) = ICOLUP(1:2,1)
         ICOLUP(1:2,4) = ICOLUP(1:2,2)
      elseif( flavor_tag.eq.2 ) then ! ZZ/WW H
         MY_IDUP(3) = -GetCKMPartner_flat( LHA2M_ID(iPart_sel) )
         MY_IDUP(4) = -GetCKMPartner_flat( LHA2M_ID(jPart_sel) )
         if( abs(MY_IDUP(3)).eq.Top_ ) return
         if( abs(MY_IDUP(4)).eq.Top_ ) return
         ICOLUP(1:2,3) = ICOLUP(1:2,1)
         ICOLUP(1:2,4) = ICOLUP(1:2,2)
      else ! WW H
         MY_IDUP(3) = -GetCKMPartner_flat( LHA2M_ID(iPart_sel) )
         MY_IDUP(4) = -GetCKMPartner_flat( LHA2M_ID(jPart_sel) )
         if( abs(MY_IDUP(3)).eq.Top_ ) return
         if( abs(MY_IDUP(4)).eq.Top_ ) return
         ICOLUP(1:2,3) = ICOLUP(1:2,1)
         ICOLUP(1:2,4) = ICOLUP(1:2,2)
      endif

      rPart_sel=convertToPartIndex(MY_IDUP(3))
      sPart_sel=convertToPartIndex(MY_IDUP(4))
      MY_IDUP(5)  = Hig_
      ICOLUP(1:2,5) = (/000,000/)

      call EvalPhasespace_VBF_deterministic(yRnd(9),yRnd(3:7),EHat,iPart_sel,jPart_sel,rPart_Sel,sPart_sel,MomExt,PSWgt)
      call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
      call Kinematics_HVBF(5,MomExt,applyPSCut,NBin)
      if( applyPSCut .or. PSWgt.eq.zero ) return
      EvalCounter = EvalCounter+1
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt  * PartChannelAvg
      if( iStore.ne.jStore ) PreFac = PreFac*2d0

      call EvalAmp_WBFH_UnSymm_SA_Select_exact( MomExt,iPart_sel,jPart_sel,rPart_Sel,sPart_Sel,me2(iPart_sel,jPart_sel))
      !call wrapHVV(MomExt,iPart_sel,jPart_sel,rPart_Sel,sPart_Sel,me2_testhvv)
      !write(6,*) me2(iPart_sel,jPart_sel)
      !write(6,*) me2_testhvv
      !pause

   elseif( Process.eq.61 ) then
      MY_IDUP(1:5) = (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel),LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel),Hig_/)! flavor default is out3=in1 out4=in2

      if( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).eq.Glu_ ) then! gg->gg/qqb
         ICOLUP(1:2,1) = (/501,502/)
         ICOLUP(1:2,2) = (/503,501/)
         if( flavor_tag.eq.2 ) then! gg->qqb
            call random_number(xRnd)
            ICOLUP(1:2,3) = (/503,000/)
            ICOLUP(1:2,4) = (/000,502/)
            if( xRnd.lt.1d0/5d0 ) then
               MY_IDUP(3:4) = (/Up_,AUp_/)
            elseif( xRnd.lt.2d0/5d0 ) then
               MY_IDUP(3:4) = (/Dn_,ADn_/)
            elseif( xRnd.lt.3d0/5d0 ) then
               MY_IDUP(3:4) = (/Chm_,AChm_/)
            elseif( xRnd.lt.4d0/5d0 ) then
               MY_IDUP(3:4) = (/Str_,AStr_/)
            elseif( xRnd.lt.5d0/5d0 ) then
               MY_IDUP(3:4) = (/Bot_,ABot_/)
            endif
         else! gg->gg
            ICOLUP(1:2,3) = (/504,502/)
            ICOLUP(1:2,4) = (/503,504/)
         endif
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).gt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qg->qg
         ICOLUP(1:2,1) = (/501,000/)
         ICOLUP(1:2,2) = (/502,501/)
         ICOLUP(1:2,3) = (/503,000/)
         ICOLUP(1:2,4) = (/502,503/)
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).lt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qbg->qbg
         ICOLUP(1:2,1) = (/000,501/)
         ICOLUP(1:2,2) = (/501,502/)
         ICOLUP(1:2,3) = (/000,503/)
         ICOLUP(1:2,4) = (/503,502/)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).gt.0 ) then! gq->gq
         ICOLUP(1:2,2) = (/501,000/)
         ICOLUP(1:2,1) = (/502,501/)
         ICOLUP(1:2,4) = (/503,000/)
         ICOLUP(1:2,3) = (/502,503/)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).lt.0 ) then! gqb->gqb
         ICOLUP(1:2,2) = (/000,501/)
         ICOLUP(1:2,1) = (/501,502/)
         ICOLUP(1:2,4) = (/000,503/)
         ICOLUP(1:2,3) = (/503,502/)
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0 ) then! qqb->qqb, qqb'->qqb'
         ICOLUP(1:2,1) = (/501,000/)
         ICOLUP(1:2,2) = (/000,502/)
         ICOLUP(1:2,3) = (/501,000/)
         ICOLUP(1:2,4) = (/000,502/)
         if( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.1 ) then! qqb->gg
            MY_IDUP(3:4) = (/Glu_,Glu_/)
            ICOLUP(1:2,1) = (/501,000/)
            ICOLUP(1:2,2) = (/000,501/)
            ICOLUP(1:2,3) = (/502,503/)
            ICOLUP(1:2,4) = (/503,502/)
         elseif( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.3 ) then! qqb->q' qbar'
            ICOLUP(1:2,1) = (/501,000/)
            ICOLUP(1:2,2) = (/000,501/)
            ICOLUP(1:2,3) = (/502,000/)
            ICOLUP(1:2,4) = (/000,502/)
            do while (.true.) ! infinite loop, sorry bad programming...
               call random_number(xRnd)
               if( xRnd.lt.1d0/5d0 ) then
                  MY_IDUP(3:4) = (/Up_,AUp_/)
               elseif( xRnd.lt.2d0/5d0 ) then
                  MY_IDUP(3:4) = (/Dn_,ADn_/)
               elseif( xRnd.lt.3d0/5d0 ) then
                  MY_IDUP(3:4) = (/Chm_,AChm_/)
               elseif( xRnd.lt.4d0/5d0 ) then
                  MY_IDUP(3:4) = (/Str_,AStr_/)
               elseif( xRnd.lt.5d0/5d0 ) then
                  MY_IDUP(3:4) = (/Bot_,ABot_/)
               endif
               if( abs(MY_IDUP(3)).ne.abs(MY_IDUP(1)) ) exit
            enddo
         endif
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).gt.0 ) then! qbq->qbq
         ICOLUP(1:2,2) = (/501,000/)
         ICOLUP(1:2,1) = (/000,501/)
         ICOLUP(1:2,4) = (/502,000/)
         ICOLUP(1:2,3) = (/000,502/)
         if( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.1 ) then! qbq->gg
            MY_IDUP(3:4) = (/Glu_,Glu_/)
            ICOLUP(1:2,2) = (/501,000/)
            ICOLUP(1:2,1) = (/000,501/)
            ICOLUP(1:2,4) = (/502,503/)
            ICOLUP(1:2,3) = (/503,502/)
         elseif( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.3 ) then! qbq->qbar'q'
            ICOLUP(1:2,2) = (/501,000/)
            ICOLUP(1:2,1) = (/000,501/)
            ICOLUP(1:2,4) = (/502,000/)
            ICOLUP(1:2,3) = (/000,502/)
            do while (.true.) ! infinite loop, sorry bad programming...
               call random_number(xRnd)
               if( xRnd.lt.1d0/5d0 ) then
                   MY_IDUP(3:4) = (/AUp_,Up_/)
               elseif( xRnd.lt.2d0/5d0 ) then
                   MY_IDUP(3:4) = (/ADn_,Dn_/)
               elseif( xRnd.lt.3d0/5d0 ) then
                   MY_IDUP(3:4) = (/AChm_,Chm_/)
               elseif( xRnd.lt.4d0/5d0 ) then
                   MY_IDUP(3:4) = (/AStr_,Str_/)
               elseif( xRnd.lt.5d0/5d0 ) then
                   MY_IDUP(3:4) = (/ABot_,Bot_/)
               endif
               if( abs(MY_IDUP(3)).ne.abs(MY_IDUP(1)) ) exit
            enddo
         endif
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).gt.0 ) then! qq->qq
         ICOLUP(1:2,1) = (/501,000/)
         ICOLUP(1:2,2) = (/502,000/)
         ICOLUP(1:2,3) = (/501,000/)
         ICOLUP(1:2,4) = (/502,000/)
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).lt.0 ) then! qbqb->qbqb
         ICOLUP(1:2,1) = (/000,501/)
         ICOLUP(1:2,2) = (/000,502/)
         ICOLUP(1:2,3) = (/000,501/)
         ICOLUP(1:2,4) = (/000,502/)
      endif
      ICOLUP(1:2,5) = (/000,000/)

      rPart_sel=convertToPartIndex(MY_IDUP(3))
      sPart_sel=convertToPartIndex(MY_IDUP(4))

      call EvalPhasespace_VBF_or_HJJ(yRnd(9),yRnd(3:7),EHat,MomExt,PSWgt)
      call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
      call Kinematics_HJJ(5,MomExt,applyPSCut,NBin)
      if( applyPSCut .or. PSWgt.eq.zero ) return
      EvalCounter = EvalCounter+1
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt  * PartChannelAvg
      if( iStore.ne.jStore ) PreFac = PreFac*2d0

      call EvalAmp_SBFH_UnSymm_SA_Select_exact(MomExt,iPart_sel,jPart_sel,rPart_sel,sPart_sel,me2(iPart_sel,jPart_sel))
      !write(6,*) me2(iPart_sel,jPart_sel)
      !pause
   endif

   LO_Res_Unpol = me2(iPart_sel,jPart_sel) * pdf(LHA2M_pdf(iPart_sel),1)*pdf(LHA2M_pdf(jPart_sel),2)
   EvalWeighted_HJJ_test = LO_Res_Unpol * PreFac
   VegasWeighted_HJJ = EvalWeighted_HJJ_test*VgsWgt

   ! iPart_sel,jPart_sel are no longer used in parton id determination or ME calculations, use iStore and jStore
   if( unweighted ) then

      if( warmup ) then

         CrossSec(iStore,jStore) = CrossSec(iStore,jStore) + VegasWeighted_HJJ
         CrossSecMax(iStore,jStore) = max(CrossSecMax(iStore,jStore),VegasWeighted_HJJ)

      else! not warmup

         call random_number(xRnd)
         if( VegasWeighted_HJJ.gt.CrossSecMax(iStore,jStore) ) then
            write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",VegasWeighted_HJJ, CrossSecMax(iStore,jStore)
            write(io_stdout, "(2X,A,1PE13.6,1PE13.6,1PE13.6,I3,I3)") "CrossSecMax is too small.",VegasWeighted_HJJ, CrossSecMax(iStore,jStore),VegasWeighted_HJJ/CrossSecMax(iStore,jStore),iStore,jStore
            AlertCounter = AlertCounter + 1
         elseif( VegasWeighted_HJJ .gt. xRnd*CrossSecMax(iStore,jStore) ) then
            AccepCounter = AccepCounter + 1
            AccepCounter_part(iStore,jStore) = AccepCounter_part(iStore,jStore) + 1
            call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=1d0)
            do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)
            enddo
         endif

      endif! warmup

   else! weighted

      if( VegasWeighted_HJJ.ne.0d0 ) then
         AccepCounter=AccepCounter+1
         if( writeWeightedLHE ) then
            call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=VegasWeighted_HJJ)
         endif
         do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),VegasWeighted_HJJ)
         enddo
      endif

   endif! unweighted

   !print *, "End EvalWeighted_HJJ_test"

   RETURN
END FUNCTION





 ! since the me2(:,:) array is defined from -5..+5, the iPart_sel,jPart_sel have to follow the LHE numbering convention
 FUNCTION EvalWeighted_HJJ(yRnd,VgsWgt)
 use ModKinematics
 use ModParameters
 use ModHiggsjj
 use ModMisc
#if compiler==1
 use ifport
#endif
   implicit none
   real(8) :: yRnd(1:9),VgsWgt, EvalWeighted_HJJ
   real(8) :: VegasWeighted_HJJ
   real(8) :: pdf(-6:6,1:2)
   real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
   real(8) :: MomExt(1:4,1:5), PSWgt
   real(8) :: me2(-5:5,-5:5)
   integer :: i,j,MY_IDUP(1:5),ICOLUP(1:2,1:5),NBin(1:NumHistograms),NHisto,iflip
   integer :: iPartChannel,PartChannelAvg,NumPartonicChannels,ijSel(1:121,1:3),flavor_tag
   real(8) :: LO_Res_Unpol, PreFac,xRnd,partonic_flip
   logical :: applyPSCut,ZZ_Fusion
   include 'vegas_common.f'

   EvalWeighted_HJJ = 0d0
   VegasWeighted_HJJ = 0d0

   if( Process.eq.60 ) then!  assuming everywhere that i>j  (apart from the LHE writeout)
      NumPartonicChannels = 71
      iPartChannel = int(yRnd(8) * (NumPartonicChannels)) +1 ! this runs from 1..71
      call get_VBFchannelHash(ijSel)
      iPart_sel = ijSel(iPartChannel,1)
      jPart_sel = ijSel(iPartChannel,2)
      ZZ_Fusion = .false.
      if( ijSel(iPartChannel,3).eq.1 ) ZZ_Fusion = .true.
   elseif( Process.eq.61 ) then
      NumPartonicChannels = 77
      iPartChannel = int(yRnd(8) * (NumPartonicChannels)) +1 ! this runs from 1..77
      call get_HJJchannelHash(ijSel)
      iPart_sel = ijSel(iPartChannel,1)
      jPart_sel = ijSel(iPartChannel,2)
      flavor_tag= ijSel(iPartChannel,3)
   endif
   PartChannelAvg = NumPartonicChannels


   if( unweighted .and. .not.warmup .and.  sum(AccepCounter_part(:,:)) .eq. sum(RequEvents(:,:)) ) then
      stopvegas=.true.
   endif
   if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) return

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   if (EHat.lt.M_Reso) return
   if( Process.eq.60 ) call EvalPhasespace_VBF_NEW2(yRnd(9),yRnd(3:7),EHat,MomExt,PSWgt)
   if( Process.eq.61 ) call EvalPhasespace_VBF_NEW2(yRnd(9),yRnd(3:7),EHat,MomExt,PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   if( Process.eq.60 ) call Kinematics_HVBF(5,MomExt,applyPSCut,NBin)
   if( Process.eq.61 ) call Kinematics_HJJ(5,MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   EvalCounter = EvalCounter+1

   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt  * PartChannelAvg
   call random_number(partonic_flip)
   if( partonic_flip.gt.0.5d0 ) call swapi(iPart_sel,jPart_sel)
   if( iPart_sel.ne.jPart_sel ) PreFac = PreFac*2d0

   ! Outgoing ids~incoming ids up to swapping and flavor changes, and main.F90 protexts from errors due to swapping (ie. massive mJ_mJ-fixed scenario)
   call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,3),MomExt(1:4,4) /) , (/ Not_a_particle_,LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel),Not_a_particle_ /) )
   if ( Process.eq.61 ) call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)

   if( Process.eq.60 ) then
      MY_IDUP(1:2)= (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel)/)

      if( MY_IDUP(1).gt.0 ) then ! quark
          ICOLUP(1:2,1) = (/501,000/)
      else! anti-quark
          ICOLUP(1:2,1) = (/000,501/)
      endif
      if( MY_IDUP(2).gt.0 ) then! quark
          ICOLUP(1:2,2) = (/502,000/)
      else! anti-quark
          ICOLUP(1:2,2) = (/000,502/)
      endif

      call EvalAmp_WBFH_UnSymm_SA_Select( MomExt,iPart_sel,jPart_sel,zz_fusion,iflip,me2)

      if( ZZ_Fusion ) then

          if( iflip.eq.2 ) then ! wrong configuration --> swap 3 and 4
             MY_IDUP(3:4)= (/LHA2M_ID(jPart_sel),LHA2M_ID(iPart_sel)/)
             ICOLUP(1:2,4) = ICOLUP(1:2,1)
             ICOLUP(1:2,3) = ICOLUP(1:2,2)
          else!
             MY_IDUP(3:4)= (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel)/)
             ICOLUP(1:2,3) = ICOLUP(1:2,1)
             ICOLUP(1:2,4) = ICOLUP(1:2,2)
          endif

      else! WW fusion

          if( iflip.eq.1 ) then ! wrong configuration --> swap 3 and 4  (opposite to zz case)
             MY_IDUP(3) = -GetCKMPartner( LHA2M_ID(jPart_sel) )
             MY_IDUP(4) = -GetCKMPartner( LHA2M_ID(iPart_sel) )

             if( abs(MY_IDUP(3)).eq.Top_ ) return !MY_IDUP(3) = sign(1,MY_IDUP(3))*Chm_
             if( abs(MY_IDUP(4)).eq.Top_ ) return !MY_IDUP(4) = sign(1,MY_IDUP(4))*Chm_
             ICOLUP(1:2,4) = ICOLUP(1:2,1)
             ICOLUP(1:2,3) = ICOLUP(1:2,2)
          else
             MY_IDUP(3) = -GetCKMPartner( LHA2M_ID(iPart_sel) )
             MY_IDUP(4) = -GetCKMPartner( LHA2M_ID(jPart_sel) )
             if( abs(MY_IDUP(3)).eq.Top_ ) return !MY_IDUP(3) = sign(1,MY_IDUP(3))*Chm_
             if( abs(MY_IDUP(4)).eq.Top_ ) return !MY_IDUP(4) = sign(1,MY_IDUP(4))*Chm_
             ICOLUP(1:2,3) = ICOLUP(1:2,1)
             ICOLUP(1:2,4) = ICOLUP(1:2,2)
          endif

      endif
      MY_IDUP(5)  = Hig_
      ICOLUP(1:2,5) = (/000,000/)


   elseif( Process.eq.61 ) then

      call EvalAmp_SBFH_UnSymm_SA_Select(MomExt,iPart_sel,jPart_sel,flavor_tag,iflip,me2)

      MY_IDUP(1:5) = (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel),LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel),Hig_/)! flavor default is out3=in1 out4=in2

      if( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).eq.Glu_ ) then! gg->?
          ICOLUP(1:2,1) = (/501,502/)
          ICOLUP(1:2,2) = (/503,501/)
          if( flavor_tag.eq.2 ) then! gg->qqb
             call random_number(xRnd)
             ICOLUP(1:2,3) = (/503,000/)
             ICOLUP(1:2,4) = (/000,502/)
             if( xRnd.lt.1d0/5d0 ) then
                MY_IDUP(3:4) = (/Up_,AUp_/)
             elseif( xRnd.lt.2d0/5d0 ) then
                MY_IDUP(3:4) = (/Dn_,ADn_/)
             elseif( xRnd.lt.3d0/5d0 ) then
                MY_IDUP(3:4) = (/Chm_,AChm_/)
             elseif( xRnd.lt.4d0/5d0 ) then
                MY_IDUP(3:4) = (/Str_,AStr_/)
             elseif( xRnd.lt.5d0/5d0 ) then
                MY_IDUP(3:4) = (/Bot_,ABot_/)
             endif
          else! gg->gg
             ICOLUP(1:2,3) = (/504,502/)
             ICOLUP(1:2,4) = (/503,504/)
          endif
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).gt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qg->qg
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/502,501/)
          ICOLUP(1:2,3) = (/503,000/)
          ICOLUP(1:2,4) = (/502,503/)
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).lt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qbg->qbg
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,2) = (/501,502/)
          ICOLUP(1:2,3) = (/000,503/)
          ICOLUP(1:2,4) = (/503,502/)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).gt.0 ) then! gq->gq
          ICOLUP(1:2,2) = (/501,000/)
          ICOLUP(1:2,1) = (/502,501/)
          ICOLUP(1:2,4) = (/503,000/)
          ICOLUP(1:2,3) = (/502,503/)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).lt.0 ) then! gqb->gqb
          ICOLUP(1:2,2) = (/000,501/)
          ICOLUP(1:2,1) = (/501,502/)
          ICOLUP(1:2,4) = (/000,503/)
          ICOLUP(1:2,3) = (/503,502/)
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0 ) then! qqb->qqb
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/000,502/)
          ICOLUP(1:2,3) = (/501,000/)
          ICOLUP(1:2,4) = (/000,502/)
          if( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.1 ) then! qqb->gg
             MY_IDUP(3:4) = (/Glu_,Glu_/)
             ICOLUP(1:2,1) = (/501,000/)
             ICOLUP(1:2,2) = (/000,501/)
             ICOLUP(1:2,3) = (/502,503/)
             ICOLUP(1:2,4) = (/503,502/)
          elseif( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.3 ) then! qqb->q' qbar'
             ICOLUP(1:2,1) = (/501,000/)
             ICOLUP(1:2,2) = (/000,501/)
             ICOLUP(1:2,3) = (/502,000/)
             ICOLUP(1:2,4) = (/000,502/)
             do while (.true.) ! infinite loop, sorry bad programming...
                call random_number(xRnd)
                if( xRnd.lt.1d0/5d0 ) then
                    MY_IDUP(3:4) = (/Up_,AUp_/)
                elseif( xRnd.lt.2d0/5d0 ) then
                    MY_IDUP(3:4) = (/Dn_,ADn_/)
                elseif( xRnd.lt.3d0/5d0 ) then
                    MY_IDUP(3:4) = (/Chm_,AChm_/)
                elseif( xRnd.lt.4d0/5d0 ) then
                    MY_IDUP(3:4) = (/Str_,AStr_/)
                elseif( xRnd.lt.5d0/5d0 ) then
                    MY_IDUP(3:4) = (/Bot_,ABot_/)
                endif
                if( abs(MY_IDUP(3)).ne.abs(MY_IDUP(1)) ) exit
             enddo
          endif
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).gt.0 ) then! qbq->qbq
          ICOLUP(1:2,2) = (/501,000/)
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,4) = (/502,000/)
          ICOLUP(1:2,3) = (/000,502/)
          if( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.1 ) then! qbq->gg
             MY_IDUP(3:4) = (/Glu_,Glu_/)
             ICOLUP(1:2,2) = (/501,000/)
             ICOLUP(1:2,1) = (/000,501/)
             ICOLUP(1:2,4) = (/502,503/)
             ICOLUP(1:2,3) = (/503,502/)
          elseif( MY_IDUP(1).eq.-MY_IDUP(2) .and. flavor_tag.eq.3 ) then! qbq->qbar'q'
             ICOLUP(1:2,2) = (/501,000/)
             ICOLUP(1:2,1) = (/000,501/)
             ICOLUP(1:2,4) = (/502,000/)
             ICOLUP(1:2,3) = (/000,502/)
             do while (.true.) ! infinite loop, sorry bad programming...
                call random_number(xRnd)
                if( xRnd.lt.1d0/5d0 ) then
                    MY_IDUP(3:4) = (/AUp_,Up_/)
                elseif( xRnd.lt.2d0/5d0 ) then
                    MY_IDUP(3:4) = (/ADn_,Dn_/)
                elseif( xRnd.lt.3d0/5d0 ) then
                    MY_IDUP(3:4) = (/AChm_,Chm_/)
                elseif( xRnd.lt.4d0/5d0 ) then
                    MY_IDUP(3:4) = (/AStr_,Str_/)
                elseif( xRnd.lt.5d0/5d0 ) then
                    MY_IDUP(3:4) = (/ABot_,Bot_/)
                endif
                if( abs(MY_IDUP(3)).ne.abs(MY_IDUP(1)) ) exit
             enddo
          endif
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).gt.0 ) then! qq->qq
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/502,000/)
          ICOLUP(1:2,3) = (/501,000/)
          ICOLUP(1:2,4) = (/502,000/)
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).lt.0 ) then! qbqb->qbqb
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,2) = (/000,502/)
          ICOLUP(1:2,3) = (/000,501/)
          ICOLUP(1:2,4) = (/000,502/)
      endif
      if( iflip.eq.2 ) then
        call swapi(MY_IDUP(3),MY_IDUP(4))
        call swapi(ICOLUP(1,3),ICOLUP(1,4))
        call swapi(ICOLUP(2,3),ICOLUP(2,4))
      endif
      ICOLUP(1:2,5) = (/000,000/)

   endif

   LO_Res_Unpol = me2(iPart_sel,jPart_sel) * pdf(LHA2M_pdf(iPart_sel),1)*pdf(LHA2M_pdf(jPart_sel),2)
   EvalWeighted_HJJ = LO_Res_Unpol * PreFac
   VegasWeighted_HJJ = EvalWeighted_HJJ*VgsWgt

   if( jPart_Sel.gt.iPart_sel ) call swapi(iPart_sel,jPart_sel) ! iPar,jPart are no longer used in parton id determination or ME calculations
   if( unweighted ) then

     if( warmup ) then

       CrossSec(iPart_sel,jPart_sel) = CrossSec(iPart_sel,jPart_sel) + VegasWeighted_HJJ
       CrossSecMax(iPart_sel,jPart_sel) = max(CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_HJJ)

     else! not warmup

       call random_number(xRnd)
       if( VegasWeighted_HJJ.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",VegasWeighted_HJJ, CrossSecMax(iPart_sel,jPart_sel)
         write(io_stdout, "(2X,A,1PE13.6,1PE13.6,1PE13.6,I3,I3)") "CrossSecMax is too small.",VegasWeighted_HJJ, CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_HJJ/CrossSecMax(iPart_sel,jPart_sel),iPart_sel,jPart_sel
         AlertCounter = AlertCounter + 1
       elseif( VegasWeighted_HJJ .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPart_sel,jPart_sel) = AccepCounter_part(iPart_sel,jPart_sel) + 1
         call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=1d0)
         do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
       endif

     endif! warmup

   else! weighted

      if( VegasWeighted_HJJ.ne.0d0 ) then
        AccepCounter=AccepCounter+1
        if( writeWeightedLHE ) then
          call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=VegasWeighted_HJJ)
        endif
        do NHisto=1,NumHistograms
          call intoHisto(NHisto,NBin(NHisto),VegasWeighted_HJJ)
        enddo
      endif

   endif! unweighted

   RETURN

 END FUNCTION EvalWeighted_HJJ








 FUNCTION EvalUnWeighted_HJJ(yRnd,genEvt,iPartons,RES)
 use ModKinematics
 use ModParameters
 use ModHiggsjj
 use ModMisc
#if compiler==1
 use ifport
#endif
implicit none
real(8) :: yRnd(:),VgsWgt, EvalUnWeighted_HJJ,RES(-5:5,-5:5)
real(8) :: pdf(-6:6,1:2)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:5), PSWgt
real(8) :: me2(-5:5,-5:5)
integer :: i,j,k,iPartons(1:2)
integer :: MY_IDUP(1:5),ICOLUP(1:2,1:5),NBin(1:NumHistograms),NHisto
real(8) :: LO_Res_Unpol, PreFac, CS_max, sumtot
logical :: applyPSCut,genEVT,zz_fusion
include 'csmaxvalue.f'

   EvalUnWeighted_HJJ = 0d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)


   if (EHat.lt.M_Reso) return
   if( Process.eq.60 ) call EvalPhaseSpace_VBF(EHat,M_Reso,yRnd(3:7),MomExt,PSWgt)
   if( Process.eq.61 ) call EvalPhaseSpace_VBF(EHat,M_Reso,yRnd(3:7),MomExt,PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))


   if( Process.eq.60 ) call Kinematics_HVBF(5,MomExt,applyPSCut,NBin)
   if( Process.eq.61 ) call Kinematics_HJJ(5,MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return


   ! Outgoing ids~incoming ids up to swapping and flavor changes, and main.F90 protexts from errors due to swapping (ie. massive mJ_mJ-fixed scenario)
   call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,3),MomExt(1:4,4) /) , (/ Not_a_particle_,LHA2M_ID(iPartons(1)),LHA2M_ID(iPartons(2)),Not_a_particle_ /) )
   if ( Process.eq.61 ) call EvalAlphaS()
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)
   EvalCounter = EvalCounter+1

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt



IF( GENEVT ) THEN

!    sumtot = 0d0
!    do i = -5,5
!       do j = -5,5
!          sumtot = sumtot + csmax(i,j)
!       enddo
!    enddo
!    k=0; bound(0)=0d0
!    do i = -5,5
!       do j = -5,5
!          k=k+1
!          bound(k) = bound(k-1) + csmax(i,j)/sumtot
!          if( yRnd(8).gt.bound(k-1) .and. yRnd(8).lt.bound(k)  ) then
!             ifound=i; jfound=j;
!             goto 1313
!          endif
!       enddo
!    enddo
! 1313 continue



   if( Process.eq.60 ) then
      call EvalAmp_WBFH_UnSymm_SA(MomExt,me2)

      MY_IDUP(1:2)= (/LHA2M_ID(iPartons(1)),LHA2M_ID(iPartons(2))/)
      if( MY_IDUP(1).gt.0 ) then ! quark
          ICOLUP(1:2,1) = (/501,000/)
      else! anti-quark
          ICOLUP(1:2,1) = (/000,501/)
      endif
      if( MY_IDUP(2).gt.0 ) then! quark
          ICOLUP(1:2,2) = (/502,000/)
      else! anti-quark
          ICOLUP(1:2,2) = (/000,502/)
      endif

      ZZ_fusion=.false.
      if( MY_IDUP(1).eq.MY_IDUP(2) ) ZZ_fusion=.true.
      if( any(MY_IDUP(1).eq.(/ Up_, Chm_,ADn_,AStr_,ABot_/)) .and. any(MY_IDUP(2).eq.(/ Up_, Chm_,ADn_,AStr_,ABot_/)) ) ZZ_fusion=.true.
      if( any(MY_IDUP(1).eq.(/AUp_,AChm_, Dn_, Str_, Bot_/)) .and. any(MY_IDUP(2).eq.(/AUp_,AChm_, Dn_, Str_, Bot_/)) ) ZZ_fusion=.true.

      if( ZZ_Fusion ) then
          if( (MomExt(4,1)*MomExt(4,3).lt.0d0) .and. (MomExt(4,2)*MomExt(4,4).lt.0d0) ) then ! wrong configuration --> swap 3 and 4
             MY_IDUP(3:4)= (/LHA2M_ID(iPartons(2)),LHA2M_ID(iPartons(1))/)
             ICOLUP(1:2,4) = ICOLUP(1:2,1)
             ICOLUP(1:2,3) = ICOLUP(1:2,2)
          else!
             MY_IDUP(3:4)= (/LHA2M_ID(iPartons(1)),LHA2M_ID(iPartons(2))/)
             ICOLUP(1:2,3) = ICOLUP(1:2,1)
             ICOLUP(1:2,4) = ICOLUP(1:2,2)
          endif
      else! WW fusion
          if( (MomExt(4,1)*MomExt(4,3).lt.0d0) .and. (MomExt(4,2)*MomExt(4,4).lt.0d0) ) then ! wrong configuration --> swap 3 and 4
             MY_IDUP(3:4)= (/SU2flip(LHA2M_ID(iPartons(2))),SU2flip(LHA2M_ID(iPartons(1)))/)
             if( abs(MY_IDUP(3)).eq.Top_ ) MY_IDUP(3) = sign(1,MY_IDUP(3))*Chm_
             if( abs(MY_IDUP(4)).eq.Top_ ) MY_IDUP(4) = sign(1,MY_IDUP(4))*Chm_
             ICOLUP(1:2,4) = ICOLUP(1:2,1)
             ICOLUP(1:2,3) = ICOLUP(1:2,2)
          else
             MY_IDUP(3:4)= (/SU2flip(LHA2M_ID(iPartons(1))),SU2flip(LHA2M_ID(iPartons(2)))/)
             if( abs(MY_IDUP(3)).eq.Top_ ) MY_IDUP(3) = sign(1,MY_IDUP(3))*Chm_
             if( abs(MY_IDUP(4)).eq.Top_ ) MY_IDUP(4) = sign(1,MY_IDUP(4))*Chm_
             ICOLUP(1:2,3) = ICOLUP(1:2,1)
             ICOLUP(1:2,4) = ICOLUP(1:2,2)
          endif
      endif
      MY_IDUP(5)  = Hig_
      ICOLUP(1:2,5) = (/000,000/)


   elseif( Process.eq.61 ) then
      call EvalAmp_SBFH_UnSymm_SA(MomExt,me2)
      MY_IDUP(1:5)  = (/LHA2M_ID(iPartons(1)),LHA2M_ID(iPartons(2)),LHA2M_ID(iPartons(1)),LHA2M_ID(iPartons(2)),Hig_/)

      if( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).eq.Glu_ ) then! gg->gg
          ICOLUP(1:2,1) = (/501,502/)
          ICOLUP(1:2,2) = (/503,501/)
          ICOLUP(1:2,3) = (/504,502/)
          ICOLUP(1:2,4) = (/503,504/)
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).gt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qg->qg
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/502,501/)
          ICOLUP(1:2,3) = (/503,000/)
          ICOLUP(1:2,4) = (/502,503/)
      elseif( MY_IDUP(1).ne.Glu_ .and. MY_IDUP(1).lt.0 .and. MY_IDUP(2).eq.Glu_ ) then! qbg->qbg
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,2) = (/501,502/)
          ICOLUP(1:2,3) = (/000,503/)
          ICOLUP(1:2,4) = (/503,502/)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).gt.0 ) then! gq->gq
          ICOLUP(1:2,2) = (/501,000/)
          ICOLUP(1:2,1) = (/502,501/)
          ICOLUP(1:2,4) = (/503,000/)
          ICOLUP(1:2,3) = (/502,503/)
      elseif( MY_IDUP(1).eq.Glu_ .and. MY_IDUP(2).ne.Glu_ .and. MY_IDUP(2).lt.0 ) then! gqb->gqb
          ICOLUP(1:2,2) = (/000,501/)
          ICOLUP(1:2,1) = (/501,502/)
          ICOLUP(1:2,4) = (/000,503/)
          ICOLUP(1:2,3) = (/503,502/)
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0 ) then! qqb->qqb
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/000,501/)
          ICOLUP(1:2,3) = (/502,000/)
          ICOLUP(1:2,4) = (/000,502/)
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).gt.0 ) then! qq->qq
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/502,000/)
          ICOLUP(1:2,3) = (/501,000/)
          ICOLUP(1:2,4) = (/502,000/)
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).gt.0 ) then! qbq->qbq
          ICOLUP(1:2,2) = (/501,000/)
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,4) = (/502,000/)
          ICOLUP(1:2,3) = (/000,502/)
      elseif( MY_IDUP(1).lt.0 .and. MY_IDUP(2).lt.0 ) then! qbqb->qbqb
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,2) = (/000,502/)
          ICOLUP(1:2,3) = (/000,501/)
          ICOLUP(1:2,4) = (/000,502/)
      endif
      if( (MomExt(4,1)*MomExt(4,3).lt.0d0) .and. (MomExt(4,2)*MomExt(4,4).lt.0d0) ) then
        call swapi(MY_IDUP(3),MY_IDUP(4))
        call swapi(ICOLUP(1,3),ICOLUP(1,4))
        call swapi(ICOLUP(2,3),ICOLUP(2,4))
      endif

      ICOLUP(1:2,5) = (/000,000/)
   endif

   LO_Res_Unpol =  me2(iPartons(1),iPartons(2)) * pdf(LHA2M_pdf(iPartons(1)),1)*pdf(LHA2M_pdf(iPartons(2)),2)
   EvalUnWeighted_HJJ = LO_Res_Unpol * PreFac

!    if( iPartons(1).eq.0 .and. iPartons(2).eq.0 ) then
!        CS_max = csmax(iPartons(1),iPartons(2)) * adj_par
!    else
!        CS_max = csmax(iPartons(1),iPartons(2))
!    endif

      CS_max = CSmax(iPartons(1),iPartons(2))
      if( EvalUnWeighted_HJJ.gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_HJJ, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_HJJ, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted_HJJ .gt. yRnd(14)*CS_max ) then
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPartons(1),iPartons(2)) = AccepCounter_part(iPartons(1),iPartons(2))+1
         call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5))
      else
          RejeCounter = RejeCounter + 1
      endif
      EvalCounter = EvalCounter + 1




ELSE! NOT GENEVT




   if( Process.eq.60 ) then
      call EvalAmp_WBFH_UnSymm_SA(MomExt,me2)
   elseif( Process.eq.61 ) then
      call EvalAmp_SBFH_UnSymm_SA(MomExt,me2)
   endif


   LO_Res_Unpol = 0d0
   do i = -5,5
      do j = -5,5

          LO_Res_Unpol = me2(i,j)*pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2) * PreFac
          EvalUnWeighted_HJJ = EvalUnWeighted_HJJ  + LO_Res_Unpol

          RES(i,j) = LO_Res_Unpol
          if (LO_Res_Unpol.gt.CSmax(i,j)) then
              CSmax(i,j) = LO_Res_Unpol
          endif
      enddo
   enddo



ENDIF! GENEVT


RETURN
END FUNCTION EvalUnWeighted_HJJ











END MODULE ModCrossSection_HJJ







