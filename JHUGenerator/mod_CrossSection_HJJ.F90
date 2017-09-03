MODULE ModCrossSection_HJJ
use ModHashCollection
implicit none
integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

 CONTAINS




! since the me2(:,:) array is defined from -5..+5, the iPart_sel,jPart_sel have to follow the LHE numbering convention
FUNCTION EvalWeighted_HJJ_fulldecay(yRnd,VgsWgt)
use ModKinematics
use ModParameters
#if linkMELA==1
use ModMCFMWrapper
#endif
use ModHiggsjj
use ModHiggs
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
integer,parameter :: mxpart=14 ! this has to match the MCFM parameter
real(8) :: yRnd(1:18),VgsWgt, EvalWeighted_HJJ_fulldecay
real(8) :: pdf(-6:6,1:2)           ,me2(-5:5,-5:5)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:10),PSWgt
real(8) :: p_MCFM(mxpart,1:4),msq_MCFM(-5:5,-5:5)
integer :: id_MCFM(mxpart),MY_IDUP(1:10),ICOLUP(1:2,1:10),NBin(1:NumHistograms),NHisto,jpart,ichan
integer, pointer :: ijSel(:,:)
integer :: iPartChannel,PartChannelAvg,NumPartonicChannels,iflip,i,j,k,flavor_tag  !,ijSel(1:121,1:3)
real(8) :: LO_Res_Unpol, PreFac,VegasWeighted_HJJ_fulldecay,xRnd,me2_hdk,me2_prop,me2_tmpzz(-5:5,-5:5),me2_tmpww(-5:5,-5:5),me2_zzcontr(-5:5,-5:5),me2_wwcontr(-5:5,-5:5)
logical :: applyPSCut
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10
real(8) :: s13,s14,s15,s16,s23,s24,s25,s26,s34,s35,s36,s45,s46,s56,s78,s910,s710,s89
include 'vegas_common.f'
include 'maxwt.f'
EvalWeighted_HJJ_fulldecay = 0d0

   call getRef_MCFM_qqVVqq_GenHash(ijSel)
   NumPartonicChannels=Hash_MCFM_qqVVqq_Gen_Size
   do ichan=1,Hash_MCFM_qqVVqq_Gen_Size
      if(ijSel(ichan,1).eq.Not_a_particle_ .or. ijSel(ichan,2).eq.Not_a_particle_ .or. ijSel(ichan,3).eq.Not_a_particle_ .or. ijSel(ichan,4).eq.Not_a_particle_) then
         NumPartonicChannels=NumPartonicChannels-1
      endif
   enddo
   if (NumPartonicChannels.eq.0) return
   PartChannelAvg = NumPartonicChannels

   iPartChannel = int(yRnd(18) * (NumPartonicChannels)) +1 ! this runs from 1..100
   iPart_sel = ijSel(iPartChannel,1)
   jPart_sel = ijSel(iPartChannel,2)

   if(.not. warmup) then
       call random_number(xRnd)!   throwing random number for accept-reject
       if( (ThisDmax.gt.0d0) .and. (ThisDmax .lt. xRnd*CrossSecMax(iPart_sel,jPart_sel)) ) then
         RejeCounter_part(iPart_sel,jPart_sel) = RejeCounter_part(iPart_sel,jPart_sel) + 1
         RejeCounter=RejeCounter+1
         return
       endif
   endif

   if( unweighted .and. .not.warmup .and.  sum(AccepCounter_part(:,:)) .eq. sum(RequEvents(:,:)) ) then
      stopvegas=.true.
      print *, "NumPartonicChannels=",NumPartonicChannels
   endif
   if( (unweighted) .and. (.not. warmup) .and. (AccepCounter_part(iPart_sel,jPart_sel) .ge. RequEvents(iPart_sel,jPart_Sel))  ) then
      call removeOffshellChannelFromHashRef(ijSel,iPartChannel,NumPartonicChannels,4)
      print *, "REMOVED CHANNEL ",iPart_sel,jPart_sel," FROM HASH: ",AccepCounter_part(iPart_sel,jPart_sel),RequEvents(iPart_sel,jPart_Sel)
      NumPartonicChannels = NumPartonicChannels-1
      print *, NumPartonicChannels," CHANNELS REMAINING"
      return
   endif

   call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi,EhatMin=dmax1(m4l_minmax(1),0d0)+mJJcut)
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
!       EvalWeighted_HJJ_fulldecay=PSWgt !*sHatJacobi  * ( MomExt(1:4,3).dot.MomExt(1:4,7) ) * ( MomExt(1:4,4).dot.MomExt(1:4,10) ) * ( MomExt(1:4,8).dot.MomExt(1:4,9) ) * ( MomExt(1:4,7).dot.MomExt(1:4,10) ) / EHat**8
!       return

   call boost2Lab(eta1,eta2,10,MomExt(1:4,1:10))
   PSWgt = PSWgt * PartChannelAvg


   call Kinematics_HVBF_fulldecay(MomExt,applyPSCut,NBin)
   DebugCounter(9) = DebugCounter(9) + 1
   if( applyPSCut .or. PSWgt.lt.1d-33 ) then
      return
   endif
   DebugCounter(10) = DebugCounter(10) + 1
   call SetRunningScales( (/MomExt(1:4,5)+MomExt(1:4,6),MomExt(1:4,3),MomExt(1:4,4) /) , (/ Not_a_particle_,Not_a_particle_,Not_a_particle_,Not_a_particle_ /) )
   call setPDFs(eta1,eta2,pdf)
   FluxFac = 1d0/(2d0*EHat**2)


   ! GeV conversion is now done inside EvalAmp_qqVVqq
   call convert_to_MCFM(-MomExt(1:4,inTop), p_MCFM(1,1:4))
   call convert_to_MCFM(-MomExt(1:4,inBot), p_MCFM(2,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep1P), p_MCFM(3,1:4))! check f fbar assignment
   call convert_to_MCFM(+MomExt(1:4,Lep1M), p_MCFM(4,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2P), p_MCFM(5,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2M), p_MCFM(6,1:4))
   call convert_to_MCFM(+MomExt(1:4,outTop),p_MCFM(7,1:4))
   call convert_to_MCFM(+MomExt(1:4,outBot),p_MCFM(8,1:4))


 msq_MCFM(:,:) = 0d0
#if linkMELA==1

   ! FIXME: TEMPORARY ASSIGNMENT OF I J R S
   do jpart=1,2
      id_MCFM(jpart) = ijSel(iPartChannel,jpart)
      if(id_MCFM(jpart) .ne. 0) id_MCFM(jpart)=convertLHE(id_MCFM(jpart))
   enddo
   do jpart=3,4
      id_MCFM(jpart+4) = ijSel(iPartChannel,jpart)
      if(id_MCFM(jpart+4) .ne. 0) id_MCFM(jpart+4)=convertLHE(id_MCFM(jpart+4))
   enddo
   id_MCFM(3:6) = (/ ElM_,ElP_,MuM_,MuP_ /)

   call EvalAmp_qqVVqq(id_MCFM, p_MCFM, 1, msq_MCFM) ! 1 for ZZ decay, 2 for WW decay, 3 for ZZ+WW mixture
!    call qq_ZZqq(p_MCFM,msq_MCFM)


   msq_MCFM = msq_MCFM * (100d0)**8  ! adjust msq_MCFM for GeV units of MCFM mat.el.

#else
 print *, "To use this process, please set linkMELA=Yes in the makefile and recompile."
 print *, "You will also need to have a compiled JHUGenMELA in the directory specified by JHUGenMELADir in the makefile."
 stop 1
#endif


   LO_Res_Unpol = msq_MCFM(iPart_sel,jPart_sel)  *  pdf(LHA2M_pdf(iPart_sel),1) * pdf(LHA2M_pdf(jPart_sel),2)
   PreFac = fbGeV2 * FluxFac * PSWgt * sHatJacobi
   EvalWeighted_HJJ_fulldecay = LO_Res_Unpol * PreFac
   VegasWeighted_HJJ_fulldecay = EvalWeighted_HJJ_fulldecay*VgsWgt

   if( unweighted ) then

     if( warmup ) then

!        if( VegasWeighted_HJJ_fulldecay.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
!            print *, "New max",iPart_sel,jPart_sel,VegasWeighted_HJJ_fulldecay
!        endif

       CrossSec(iPart_sel,jPart_sel) = CrossSec(iPart_sel,jPart_sel) + VegasWeighted_HJJ_fulldecay
       CrossSecMax(iPart_sel,jPart_sel) = max(CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_HJJ_fulldecay)

     else! not warmup

       EvalCounter = EvalCounter+1

       if( VegasWeighted_HJJ_fulldecay .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then
         RejeCounter_part(iPart_sel,jPart_sel) = RejeCounter_part(iPart_sel,jPart_sel) + 1
         RejeCounter=RejeCounter+1
       endif

       if( VegasWeighted_HJJ_fulldecay.gt.CrossSecMax(iPart_sel,jPart_sel) ) then
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CrossSecMax is too small.",VegasWeighted_HJJ_fulldecay, CrossSecMax(iPart_sel,jPart_sel)
         write(io_stdout, "(2X,A,1PE13.6,1PE13.6,1PE13.6,I3,I3)") "CrossSecMax is too small.",VegasWeighted_HJJ_fulldecay, CrossSecMax(iPart_sel,jPart_sel),VegasWeighted_HJJ_fulldecay/CrossSecMax(iPart_sel,jPart_sel),iPart_sel,jPart_sel
         AlertCounter = AlertCounter + 1

!          This dynamically increases the maximum in case it is exceeded
!          CrossSecMax(iPart_sel,jPart_sel) = VegasWeighted_HJJ_fulldecay
!          write(io_LogFile,"(2X,A,1PE13.6)") "Increasing CrossSecMax to ",VegasWeighted_HJJ_fulldecay
!          write(io_stdout, "(2X,A,1PE13.6)") "Increasing CrossSecMax to ",VegasWeighted_HJJ_fulldecay

       elseif( VegasWeighted_HJJ_fulldecay .gt. xRnd*CrossSecMax(iPart_sel,jPart_sel) ) then
         AccepCounter = AccepCounter + 1
         AccepCounter_part(iPart_sel,jPart_sel) = AccepCounter_part(iPart_sel,jPart_sel) + 1

         MY_IDUP(1:2)= (/LHA2M_ID(iPart_sel),LHA2M_ID(jPart_sel)/)
         MY_IDUP(3:4)= (/LHA2M_ID(id_MCFM(7)),LHA2M_ID(id_MCFM(8))/)
         call WriteOutEvent_HJJ_fulldecay(MomExt,MY_IDUP,ICOLUP)

         do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
       else
         RejeCounter_part(iPart_sel,jPart_sel) = RejeCounter_part(iPart_sel,jPart_sel) + 1
         RejeCounter=RejeCounter+1
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

      endif

   endif! unweighted


RETURN
END FUNCTION






FUNCTION EvalUnWeighted_HJJ_fulldecay(yRnd,genEvt,iPartons,RES)
implicit none
real(8) :: yRnd(1:17),EvalUnWeighted_HJJ_fulldecay,RES(-5:5,-5:5)
integer :: iPartons(1:2)
logical :: genEvt
EvalUnWeighted_HJJ_fulldecay = 0d0


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
   integer, pointer :: ijSel(:,:)
   integer :: iPartChannel,PartChannelAvg,NumPartonicChannels,flavor_tag
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
      call getRef_VBFchannelHash_nosplit(ijSel,NumPartonicChannels)

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
      call getRef_HJJchannelHash_nosplit(ijSel,NumPartonicChannels)
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
   integer, pointer :: ijSel(:,:)
   integer :: iPartChannel,PartChannelAvg,NumPartonicChannels,flavor_tag
   real(8) :: LO_Res_Unpol, PreFac,xRnd,partonic_flip
   logical :: applyPSCut,ZZ_Fusion
   include 'vegas_common.f'

   EvalWeighted_HJJ = 0d0
   VegasWeighted_HJJ = 0d0

   if( Process.eq.60 ) then!  assuming everywhere that i>j  (apart from the LHE writeout)
      NumPartonicChannels = 71
      iPartChannel = int(yRnd(8) * (NumPartonicChannels)) +1 ! this runs from 1..71
      call getRef_VBFchannelHash(ijSel)
      iPart_sel = ijSel(iPartChannel,1)
      jPart_sel = ijSel(iPartChannel,2)
      ZZ_Fusion = .false.
      if( ijSel(iPartChannel,3).eq.1 ) ZZ_Fusion = .true.
   elseif( Process.eq.61 ) then
      NumPartonicChannels = 77
      iPartChannel = int(yRnd(8) * (NumPartonicChannels)) +1 ! this runs from 1..77
      call getRef_HJJchannelHash(ijSel)
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







