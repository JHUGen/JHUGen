MODULE ModCrossSection_TTBH
IMPLICIT NONE

integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)


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
   call EvalPhasespace_2to3M(EHat,(/M_Reso,M_Top,M_Top/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)! a(1)b(2)-->H(3)+tbar(4)+t(5)
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

   call Kinematics_TTBH(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return

!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,1)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,2)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,3)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,4)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,5)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,6)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,7)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,8)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,9)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,10)
!    write(*,"(PE21.14,PE21.14,PE21.14,PE21.14)") MomExt(1:4,11)

   Mu_Fact = 0.5d0*( 2d0*M_top + M_Reso )
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
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:13), PSWgt,PSWgt2,PSWgt3,CS_Max
real(8) :: LO_Res_GG_Unpol,LO_Res_QQB_Unpol,PreFac,PDFFac1,PDFFac2
integer :: NBin(1:NumHistograms),NHisto,iPartons(1:2)
integer :: MY_IDUP(1:11),ICOLUP(1:2,1:11),nparton
logical :: applyPSCut,genEvt
include 'csmaxvalue.f'  
EvalUnWeighted_TTBH = 0d0


   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   
   if (EHat.lt.2*M_Top+M_Reso) return
   call EvalPhasespace_2to3M(EHat,(/M_Reso,M_Top,M_Top/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)! a(1)b(2)-->H(3)+tbar(4)+t(5)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   ICOLUP(1:2,3) = (/000,000/)
   ICOLUP(1:2,4) = (/000,502/)
   ICOLUP(1:2,5) = (/501,000/)
   if( TOPDECAYS.NE.0 ) then
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(08:11),MomExt(1:4,06:08),PSWgt2)    ! ATop 
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),MomExt(1:4,09:11),PSWgt3)    !  Top
      PSWgt = PSWgt * PSWgt2*PSWgt3
      if( TOPDECAYS.EQ.1 ) then
          MY_IDUP(6:11)=(/ABot_,ElM_,ANuE_,Bot_,MuP_,NuM_/)
          ICOLUP(1:2,6) = (/000,502/)
          ICOLUP(1:2,7) = (/000,000/)
          ICOLUP(1:2,8) = (/000,000/)
          ICOLUP(1:2,9) = (/501,000/)
          ICOLUP(1:2,10)= (/000,000/)
          ICOLUP(1:2,11)= (/000,000/)             
      elseif( TOPDECAYS.EQ.2 ) then
          MY_IDUP(6:11)=(/ABot_,Dn_,AUp_,Bot_,ADn_,Up_/)
          ICOLUP(1:2,6) = (/000,502/)
          ICOLUP(1:2,7) = (/503,000/)
          ICOLUP(1:2,8) = (/000,503/)
          ICOLUP(1:2,9) = (/501,000/)
          ICOLUP(1:2,10)= (/000,504/)
          ICOLUP(1:2,11)= (/504,000/)   
      elseif( TOPDECAYS.EQ.3 ) then
          MY_IDUP(6:11)=(/ABot_,ElM_,ANuE_,Bot_,ADn_,Up_/)
          ICOLUP(1:2,6) = (/000,502/)
          ICOLUP(1:2,7) = (/000,000/)
          ICOLUP(1:2,8) = (/000,000/)
          ICOLUP(1:2,9) = (/501,000/)
          ICOLUP(1:2,10)= (/000,504/)
          ICOLUP(1:2,11)= (/504,000/)   
      elseif( TOPDECAYS.EQ.4 ) then
          MY_IDUP(6:11)=(/ABot_,Dn_,AUp_,Bot_,MuP_,NuM_/)
          ICOLUP(1:2,6) = (/000,502/)
          ICOLUP(1:2,7) = (/503,000/)
          ICOLUP(1:2,8) = (/000,503/)
          ICOLUP(1:2,9) = (/501,000/)
          ICOLUP(1:2,10)= (/000,000/)
          ICOLUP(1:2,11)= (/000,000/)
      endif
      
   else
      MY_IDUP(6:11)=-9999
   endif 
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
!    PSWgt = PSWgt * PSWgt4 
   FluxFac = 1d0/(2d0*EHat**2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt

   call Kinematics_TTBH(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return

   Mu_Fact = 0.5d0*( 2d0*M_top + M_Reso )   
   call setPDFs(eta1,eta2,Mu_Fact,pdf)



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
         call WriteOutEvent_TTBH(MomExt,MY_IDUP(1:11),ICOLUP(1:2,1:11))
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







