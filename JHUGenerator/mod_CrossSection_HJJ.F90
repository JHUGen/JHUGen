MODULE ModCrossSection_HJJ
implicit none
integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

 CONTAINS



 
 
FUNCTION EvalWeighted_HJJ_fulldecay(yRnd,VgsWgt)
use ModKinematics
use ModParameters
use ModHiggsjj
use ModMisc
#if compiler==1
use ifport
#endif
implicit none
integer,parameter :: mxpart=14 ! this has to match the MCFM parameter
real(8) :: yRnd(1:17),VgsWgt, EvalWeighted_HJJ_fulldecay
real(8) :: pdf(-6:6,1:2),me2(-5:5,-5:5)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:11),PSWgt,PSWgt1,PSWgt2,PSWgt3,MInvH,MInvZ1,MInvZ2
real(8) :: yz1,yz2,offzchannel,EZ_max,dr,ML1,ML2,ML3,ML4
real(8) :: p_MCFM(mxpart,1:4),msq_MCFM(-5:5,-5:5),HZZcoupl(1:32),HWWcoupl(1:32)
integer :: i,j,MY_IDUP(1:5),ICOLUP(1:2,1:5),NBin(1:NumHistograms),NHisto
real(8) :: LO_Res_Unpol, PreFac,BWJacobi
logical :: applyPSCut
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, Higgs=5, V1=6, V2=7, Lep1P=8, Lep1M=9, Lep2P=10, Lep2M=11
   
   EvalWeighted_HJJ_fulldecay = 0d0
   
   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)   
   MInvH = M_Reso

       BWJacobi = 1d0
!      call SmearExternal(yRnd(17),M_Reso,Ga_Reso,M_Reso/4d0,M_Reso*4d0,MInvH,BWJacobi)
       

   
   if (EHat.lt.MInvH) return
   call EvalPhaseSpace_VBF(EHat,MInvH,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))
   

    yz1 = yRnd(14)
    yz2 = yRnd(15)
    offzchannel = yRnd(16) 
! yz1=0.2d0
! yz2=0.5d0
! offzchannel=0.8d0

         if(MInvH.gt.2d0*M_V) then
              EZ_max = EHat
              dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
              MInvZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
              sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MInvZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

              EZ_max = EHat - MInvZ1*0.99999d0
              dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
              MInvZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
              sHatJacobi = sHatJacobi*dr/(Ga_V*M_V)  *  ( (MInvZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

          elseif(MInvH.lt.2d0*M_V) then
              if (offzchannel.le.0.5d0) then
                  EZ_max = EHat
                  dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
                  MInvZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
                  MInvZ2 = abs(EHat - MInvZ1*0.999999999999999d0)*dsqrt(dabs(dble(yz2)))
                  sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/(  &
                  1d0/((MInvZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )   )
                  sHatJacobi = sHatJacobi *(EHat - MInvZ1*0.99999d0)**2
              elseif(offzchannel.gt.0.5d0) then
                  EZ_max = EHat
                  dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
                  MInvZ2 = dsqrt( M_V*Ga_V * dtan(dr*yz2-datan(M_V/Ga_V)) + M_V**2 )
                  MInvZ1 = abs(EHat - MInvZ2*0.999999999999999d0)*dsqrt(dabs(dble(yz1)))
                  sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * 1d0/(  &
                  + 1d0/((MInvZ2**2 - M_V**2)**2 + M_V**2*Ga_V**2 ) )
                  sHatJacobi = sHatJacobi *(EHat - MInvZ2*0.99999d0)**2
              endif
        endif

    if( MInvZ1+MInvZ2.gt.MInvH ) then
      EvalWeighted_HJJ_fulldecay = 0d0
      RejeCounter = RejeCounter + 1
      return
    endif
    
    call EvalPhaseSpace_2(MInvH,(/MInvZ1,MInvZ2/),yRnd(8:9),MomExt(1:4,6:7),PSWgt1)
    ML1 = 0d0
    ML2 = 0d0
    ML3 = 0d0
    ML4 = 0d0
    call boost(MomExt(1:4,6),MomExt(1:4,5),MInvH)
    call boost(MomExt(1:4,7),MomExt(1:4,5),MInvH)
    
    
    call EvalPhasespace_VDecay(MomExt(1:4,6),MInvZ1,ML1,ML2,yRnd(10:11),MomExt(1:4,8:9),PSWgt2)
    call EvalPhasespace_VDecay(MomExt(1:4,7),MInvZ2,ML3,ML4,yRnd(12:13),MomExt(1:4,10:11),PSWgt3)
    PSWgt = PSWgt * PSWgt1*PSWgt2*PSWgt3 * BWJacobi

!      if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then! introduce this momentum flip to allow proper mapping of integrand with Z-poles at MInvZ2=(p2+p3)^2 and MInvZ2=(p1+p4)^2
!          if( yrnd(16).gt.0.5d0 ) call swapmom( MomExt(1:4,8),MomExt(1:4,10) )
!      endif

   call Kinematics_HVBF_fulldecay(MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   
   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)
   


   call convert_to_MCFM(-MomExt(1:4,inTop), p_MCFM(1,1:4))
   call convert_to_MCFM(-MomExt(1:4,inBot), p_MCFM(2,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep1P), p_MCFM(3,1:4))! check f fbar assignment
   call convert_to_MCFM(+MomExt(1:4,Lep1M), p_MCFM(4,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2P), p_MCFM(5,1:4))
   call convert_to_MCFM(+MomExt(1:4,Lep2M), p_MCFM(6,1:4))
   call convert_to_MCFM(+MomExt(1:4,outTop),p_MCFM(7,1:4))
   call convert_to_MCFM(+MomExt(1:4,outBot),p_MCFM(8,1:4))
   HZZcoupl(:) = 0d0
   HZZcoupl(1) = 0d0
   HWWcoupl(:) = HZZcoupl(:)

!    print *, (MomExt(1:4,1)).dot.(MomExt(1:4,1))
!    print *, (MomExt(1:4,2)).dot.(MomExt(1:4,2))
!    print *, (MomExt(1:4,8)).dot.(MomExt(1:4,8))
!    print *, (MomExt(1:4,9)).dot.(MomExt(1:4,9))
!    print *, (MomExt(1:4,10)).dot.(MomExt(1:4,10))
!    print *, (MomExt(1:4,11)).dot.(MomExt(1:4,11))
!    print *, (MomExt(1:4,3)).dot.(MomExt(1:4,3))
!    print *, (MomExt(1:4,4)).dot.(MomExt(1:4,4))
!    print *, "---"
!    print *, p_MCFM(1,1:4)+p_MCFM(2,1:4)  +p_MCFM(3,1:4)+p_MCFM(4,1:4)  +p_MCFM(5,1:4)+p_MCFM(6,1:4)  +p_MCFM(7,1:4)+p_MCFM(8,1:4)
!    print *, p_MCFM(1,4)**2-p_MCFM(1,1)**2-p_MCFM(1,2)**2-p_MCFM(1,3)**2
!    print *, p_MCFM(2,4)**2-p_MCFM(2,1)**2-p_MCFM(2,2)**2-p_MCFM(2,3)**2
!    print *, p_MCFM(3,4)**2-p_MCFM(3,1)**2-p_MCFM(3,2)**2-p_MCFM(3,3)**2
!    print *, p_MCFM(4,4)**2-p_MCFM(4,1)**2-p_MCFM(4,2)**2-p_MCFM(4,3)**2
!    print *, p_MCFM(5,4)**2-p_MCFM(5,1)**2-p_MCFM(5,2)**2-p_MCFM(5,3)**2
!    print *, p_MCFM(6,4)**2-p_MCFM(6,1)**2-p_MCFM(6,2)**2-p_MCFM(6,3)**2
!    print *, p_MCFM(7,4)**2-p_MCFM(7,1)**2-p_MCFM(7,2)**2-p_MCFM(7,3)**2
!    print *, p_MCFM(8,4)**2-p_MCFM(8,1)**2-p_MCFM(8,2)**2-p_MCFM(8,3)**2
!    pause

!   call qq_ZZqq(p_MCFM,msq_MCFM,HZZcoupl,HWWcoupl,Lambda,Lambda_Q,Lambda_z1)!  q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8)
   
   LO_Res_Unpol = 0d0
   do i = -5,5
      do j = -5,5
         LO_Res_Unpol = LO_Res_Unpol + msq_MCFM(i,j)*pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2)
      enddo
   enddo


   PreFac = fbGeV2 * FluxFac * PSWgt * sHatJacobi
   EvalWeighted_HJJ_fulldecay = LO_Res_Unpol * PreFac

   AccepCounter=AccepCounter+1
   if( writeWeightedLHE ) then 
       call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=EvalWeighted_HJJ_fulldecay*VgsWgt)
   endif
   do NHisto=1,NumHistograms
       call intoHisto(NHisto,NBin(NHisto),EvalWeighted_HJJ_fulldecay*VgsWgt)
   enddo
   EvalCounter = EvalCounter+1




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
real(8) :: yRnd(:),VgsWgt, EvalUnWeighted_HJJ_fulldecay,RES(-5:5,-5:5)
real(8) :: pdf(-6:6,1:2)
real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
real(8) :: MomExt(1:4,1:5), PSWgt
real(8) :: me2(-5:5,-5:5)
integer :: i,j,k,iPartons(1:2)
integer :: MY_IDUP(1:5),ICOLUP(1:2,1:5),NBin(1:NumHistograms),NHisto
real(8) :: LO_Res_Unpol, PreFac, CS_max, sumtot
logical :: applyPSCut,genEvt,zz_fusion
include 'csmaxvalue.f'

 EvalUnWeighted_HJJ_fulldecay = 0d0
 
RETURN
END FUNCTION




 


 FUNCTION EvalWeighted_HJJ(yRnd,VgsWgt)
 use ModKinematics
 use ModParameters
 use ModHiggsjj
#if compiler==1
 use ifport
#endif
   implicit none
   real(8) :: yRnd(1:7),VgsWgt, EvalWeighted_HJJ
   real(8) :: pdf(-6:6,1:2)
   real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
   real(8) :: MomExt(1:4,1:5), PSWgt
   real(8) :: me2(-5:5,-5:5)
   integer :: i,j,MY_IDUP(1:5),ICOLUP(1:2,1:5),NBin(1:NumHistograms),NHisto
   real(8) :: LO_Res_Unpol, PreFac
   logical :: applyPSCut
   
   EvalWeighted_HJJ = 0d0

   call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   EvalCounter = EvalCounter+1

   if (EHat.lt.M_Reso) return
   if( Process.eq.60 ) call EvalPhaseSpace_VBF(EHat,M_Reso,yRnd(3:7),MomExt,PSWgt)
   if( Process.eq.61 ) call EvalPhaseSpace_VBF(EHat,M_Reso,yRnd(3:7),MomExt,PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))


   if( Process.eq.60 ) call Kinematics_HVBF(5,MomExt,applyPSCut,NBin)
   if( Process.eq.61 ) call Kinematics_HJJ(5,MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)

   if (process.eq.60) then
      call EvalAmp_WBFH_UnSymm_SA(MomExt,(/ghz1,ghz2,ghz3,ghz4/),(/ghw1,ghw2,ghw3,ghw4/),me2)
      MY_IDUP(1:5)  = (/Up_,Up_,Up_,Up_,Hig_/)
      ICOLUP(1:2,1) = (/501,000/)
      ICOLUP(1:2,2) = (/502,000/)
      ICOLUP(1:2,3) = (/501,000/)
      ICOLUP(1:2,4) = (/502,000/)
      ICOLUP(1:2,5) = (/000,000/)
   elseif (process.eq.61) then
      call EvalAmp_SBFH_UnSymm_SA(MomExt,(/ghg2,ghg3,ghg4/),me2)
      me2 = me2 * (2d0/3d0*alphas**2)**2 
      MY_IDUP(1:5)  = (/Up_,Up_,Up_,Up_,Hig_/)
      ICOLUP(1:2,1) = (/501,000/)
      ICOLUP(1:2,2) = (/502,000/)
      ICOLUP(1:2,3) = (/501,000/)
      ICOLUP(1:2,4) = (/502,000/)
      ICOLUP(1:2,5) = (/000,000/)
   endif
   
   LO_Res_Unpol = 0d0
   do i = -5,5
      do j = -5,5
         LO_Res_Unpol = LO_Res_Unpol + me2(i,j)*pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2)
      enddo
   enddo

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt 
   EvalWeighted_HJJ = LO_Res_Unpol * PreFac


   AccepCounter=AccepCounter+1
   if( writeWeightedLHE ) then 
       call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5),EventWeight=EvalWeighted_HJJ*VgsWgt)
   endif
   do NHisto=1,NumHistograms
       call intoHisto(NHisto,NBin(NHisto),EvalWeighted_HJJ*VgsWgt)
   enddo




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
   

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
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
      call EvalAmp_WBFH_UnSymm_SA(MomExt,(/ghz1,ghz2,ghz3,ghz4/),(/ghw1,ghw2,ghw3,ghw4/),me2)

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
      call EvalAmp_SBFH_UnSymm_SA(MomExt,(/ghg2,ghg3,ghg4/),me2)
      me2 = me2 * (2d0/3d0*alphas**2)**2 !-- (alphas/sixpi gs^2)^2
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
      call EvalAmp_WBFH_UnSymm_SA(MomExt,(/ghz1,ghz2,ghz3,ghz4/),(/ghw1,ghw2,ghw3,ghw4/),me2)
   elseif( Process.eq.61 ) then
      call EvalAmp_SBFH_UnSymm_SA(MomExt,(/ghg2,ghg3,ghg4/),me2)
      me2 = me2 * (2d0/3d0*alphas**2)**2
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







