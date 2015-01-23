MODULE ModCrossSection_TTBH
IMPLICIT NONE

integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)


 CONTAINS


!  ./JHUGen  Process=80 TopDK=1 PChannel=0 Unweighted=0  1.4725086

! ./JHUGen  Process=80 TopDK=1 PChannel=1 Unweighted=0   0.70816989

!   ./JHUGen  Process=80 TopDK=1 PChannel=2 Unweighted=0  2.1803083

! ratio gg/qqb = 2.08

! unweighting:
!   Acceptance  Counter_part:            0                   620
!   Acceptance  Counter_part:            1                   380

!   Acceptance  Counter_part:            0                  6417
!   Acceptance  Counter_part:            1                  3583




FUNCTION EvalWeighted_TTBH(yRnd,VgsWgt)
use ModKinematics
use ModParameters
use ModTTBHiggs
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
real(8) :: yRnd(1:15),VgsWgt, EvalWeighted_TTBH
real(8) :: pdf(-6:6,1:2)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi, PDFFac1,PDFFac2
real(8) :: MomExt(1:4,1:13), PSWgt,PSWgt2,PSWgt3
real(8) :: LO_Res_GG_Unpol,LO_Res_QQB_Unpol, PreFac, MG_MOM(0:3,1:5),MadGraph_tree
integer :: NBin(1:NumHistograms),NHisto
logical :: applyPSCut
EvalWeighted_TTBH = 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   
   if (EHat.lt.2*M_Top+M_Reso) return
   call EvalPhasespace_2to3M(EHat,M_Reso,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! a(1)b(2)-->H(3)+tbar(4)+t(5)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   
   if( TOPDECAYS.NE.0 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(08:11),MomExt(1:4,06:08),PSWgt2)    ! ATop 
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),MomExt(1:4,09:11),PSWgt3)    !  Top
      PSWgt = PSWgt * PSWgt2*PSWgt3
   endif 
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
!    PSWgt = PSWgt * PSWgt4 
   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt

!  MomExt(1:4,1)=(/   3.2772203957555925d0,        0.0000000000000000d0,        0.0000000000000000d0,        3.2772203957555925d0     /)
!  MomExt(1:4,2)=(/     2.8653204768734621d0,        0.0000000000000000d0,        0.0000000000000000d0,       -2.8653204768734621d0         /)
!  MomExt(1:4,3)=(/     2.0695031298922770d0,       -1.3995194639860060d0,       0.4309336522215461d0,      -0.74896506056107459d0     /)
!  MomExt(1:4,4)=(/     1.9107949309969310d0,       0.80690401116821620d0,      -3.12533731972295947d-002,  7.85265034749671742d-002    /)
!  MomExt(1:4,5)=(/    2.1622428117398473d0,       0.59261545281778993d0,      -0.39968027902431658d0,        1.0823384759682382d0         /)
!  MomExt(1:4,6)=(/    1.0151002823876631d0,       0.96306893538507010d0,       0.32070427775393112d0,      -8.69340152707267361d-003    /)
!  MomExt(1:4,7)=(/    0.26916769467215163d0,      -7.41472049268541988d-002,  0.24741337553749249d0,      -7.57631933183877115d-002    /)
!  MomExt(1:4,8)=(/    0.62652695393711610d0,      -8.20177192899996244d-002, -0.59937102648865326d0,       0.16298309832042757d0     /)
!  MomExt(1:4,9)=(/    0.63056238611088244d0,       0.18199675333260884d0,       0.53765938388191470d0,       0.27460606598900683d0         /)  
!  MomExt(1:4,10)=(/    0.34622959545354920d0,      -2.44789896017240799d-002, -0.32712245562490633d0,      -0.11075473290987667d0         /)
!  MomExt(1:4,11)=(/     1.1854508301754154d0,       0.43509768908690516d0,      -0.61021720728132500d0,       0.91848714288910793d0        /)  
  
   call Kinematics_TTBH(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   
   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   if( PChannel.eq.0 .or. PChannel.eq.2 ) then
      call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
      PDFFac1 = pdf(0,1)*pdf(0,2)
      EvalWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1
   endif
   if( PChannel.eq.1 .or. PChannel.eq.2 ) then
      call EvalAmp_QQB_TTBH(MomExt,LO_Res_QQB_Unpol)
      PDFFac1 = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
              + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
              + pdf(Bot_,1)*pdf(ABot_,2)                            
      PDFFac2 = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
              + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
              + pdf(Bot_,2)*pdf(ABot_,1)
      EvalWeighted_TTBH = EvalWeighted_TTBH + LO_Res_QQB_Unpol * ( PDFFac1 + PDFFac2 )  
   endif
   EvalWeighted_TTBH = EvalWeighted_TTBH * PreFac
   
! print *, "gg",LO_Res_GG_Unpol
! print *, "qq",LO_Res_QQB_Unpol
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


   AccepCounter=AccepCounter+1
   if( writeWeightedLHE ) then 
        call Error("WriteLHE not yet supported for ttb+H")
!        call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=EvalWeighted_TTBH*VgsWgt)
   endif
   do NHisto=1,NumHistograms
       call intoHisto(NHisto,NBin(NHisto),EvalWeighted_TTBH*VgsWgt)
   enddo
   EvalCounter = EvalCounter+1



RETURN
END FUNCTION 







FUNCTION EvalUnWeighted_TTBH(yRnd,genEvt,RES)
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
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi,bound(0:121)
real(8) :: MomExt(1:4,1:13), PSWgt,PSWgt2,PSWgt3,CS_Max,sumtot
real(8) :: LO_Res_GG_Unpol,LO_Res_QQB_Unpol,PreFac,PDFFac1,PDFFac2
integer :: NBin(1:NumHistograms),NHisto
integer :: MY_IDUP(1:11),ICOLUP(1:2,1:11)
integer :: nparton,k,ifound,jfound
logical :: applyPSCut,genEvt
include 'csmaxvalue.f'  
EvalUnWeighted_TTBH = 0d0
! RES(:,:) = 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   
   if (EHat.lt.2*M_Top+M_Reso) return
   call EvalPhasespace_2to3M(EHat,M_Reso,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! a(1)b(2)-->H(3)+tbar(4)+t(5)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   
   if( TOPDECAYS.NE.0 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(08:11),MomExt(1:4,06:08),PSWgt2)    ! ATop 
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),MomExt(1:4,09:11),PSWgt3)    !  Top
      PSWgt = PSWgt * PSWgt2*PSWgt3
      if( TOPDECAYS.EQ.1 ) MY_IDUP(6:11)=(/ABot_,ElM_,ANuE_,Bot_,MuP_,NuM_/)
      if( TOPDECAYS.EQ.2 ) MY_IDUP(6:11)=(/ABot_,Dn_,AUp_,Bot_,ADn_,Up_/)
      if( TOPDECAYS.EQ.3 ) MY_IDUP(6:11)=(/ABot_,ElM_,ANuE_,Bot_,ADn_,Up_/)
      if( TOPDECAYS.EQ.4 ) MY_IDUP(6:11)=(/ABot_,Dn_,AUp_,Bot_,MuP_,NuM_/)
   endif 
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
!    PSWgt = PSWgt * PSWgt4 
   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt

   call Kinematics_TTBH(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   

IF( GENEVT ) THEN   

      sumtot = 0d0
      do nparton = -5,5
!           do j = -5,5
!             sumtot = sumtot + csmax(nparton,-nparton)!   WHY IS THIS CSMAX AND NOT CS ? ! AND SHOULDNT IT RUN FROM -6..6
            sumtot = sumtot + VG(nparton,-nparton)!   WHY IS THIS CSMAX AND NOT CS ? ! AND SHOULDNT IT RUN FROM -6..6
!           enddo
      enddo
      
      
      k=0; bound(0)=0d0
      do nparton = -5,5
!           do j = -5,5
            k=k+1
!             bound(k) = bound(k-1) + csmax(nparton,-nparton)/sumtot
            bound(k) = bound(k-1) + VG(nparton,-nparton)/sumtot
            if( yRnd(8).gt.bound(k-1) .and. yRnd(8).lt.bound(k)  ) then
                ifound=nparton; jfound=-nparton;
                exit
            endif
! !           enddo
      enddo

!       put in counter here to see if bound really generates the correct fraction
!       print *, "Xx ",bound(:)
!       pause
      
      
      call setPDFs(eta1,eta2,Mu_Fact,pdf)

      if( PChannel.eq.2 ) then
          PDFFac1 = pdf(ifound,1)*pdf(jfound,2)
          CS_max = CSmax(ifound,jfound)
          MY_IDUP(1:5) = (/LHA2M_ID(ifound),LHA2M_ID(jfound),Hig_,ATop_,Top_/)
          ICOLUP(1:2,1:11) = 0
          if( ifound.eq.0 ) then
              call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
              EvalUnWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1 * PreFac
              CS_max = CSmax(0,0)
          else
              call EvalAmp_QQB_TTBH(MomExt,LO_Res_QQB_Unpol)
              EvalUnWeighted_TTBH = LO_Res_QQB_Unpol * PDFFac1 * PreFac
              CS_max = CSmax(ifound,jfound)
          endif
          
      elseif( PChannel.eq.0 ) then
          call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
          PDFFac1 = pdf(0,1)*pdf(0,2)
          EvalUnWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1 * PreFac
          CS_max = CSmax(0,0)
          MY_IDUP(1:5) = (/Glu_,Glu_,Hig_,ATop_,Top_/)
          ICOLUP(1:2,1:11) = 0
          
      elseif( PChannel.eq.1 ) then
          call EvalAmp_QQB_TTBH(MomExt,LO_Res_QQB_Unpol)
          PDFFac1 = pdf(ifound,1)*pdf(jfound,2)
          EvalUnWeighted_TTBH = LO_Res_QQB_Unpol * PDFFac1 * PreFac 
          CS_max = CSmax(ifound,jfound)
          MY_IDUP(1:5) = (/Up_,AUp_,Hig_,ATop_,Top_/)
          ICOLUP(1:2,1:11) = 0
      endif

      if( EvalUnWeighted_TTBH .gt. CS_max) then
         write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_TTBH, CS_max
         write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_TTBH, CS_max
         AlertCounter = AlertCounter + 1
!          RES = 0d0
      elseif( EvalUnWeighted_TTBH .gt. yRnd(16)*CS_max ) then
         AccepCounter = AccepCounter + 1
         if( ifound.eq.0 ) then
           AccepCounter_part(0,0)=AccepCounter_part(0,0)+1
         else
           AccepCounter_part(1,-1)=AccepCounter_part(1,-1)+1
         endif
         call WriteOutEvent_TTBH(MomExt,MY_IDUP(1:11),ICOLUP(1:2,1:11))
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)
         enddo
      else
         RejeCounter = RejeCounter + 1
      endif
      

ELSE! NOT GENEVT

      call setPDFs(eta1,eta2,Mu_Fact,pdf)

      
      if( PChannel.eq.0 .or. PChannel.eq.2 ) then
          call EvalAmp_GG_TTBH(MomExt,LO_Res_GG_Unpol)
          PDFFac1 = pdf(0,1)*pdf(0,2)
          EvalUnWeighted_TTBH = LO_Res_GG_Unpol * PDFFac1 * PreFac
          
          RES(0,0) = EvalUnWeighted_TTBH
          if (EvalUnWeighted_TTBH.gt.csmax(0,0)) then
              csmax(0,0) = EvalUnWeighted_TTBH
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
            if (EvalUnWeighted_TTBH.gt.csmax(nparton,-nparton)) csmax(nparton,-nparton) = EvalUnWeighted_TTBH
          enddo
      endif


ENDIF! GENEVT 


 RETURN
 END FUNCTION






END MODULE ModCrossSection_TTBH







