MODULE ModCrossSection
implicit none
integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

contains


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
   real(8) :: MomExt(1:4,1:5), PSWgt,MomDK(1:4,1:4)
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


   if( Process.eq.60 ) call Kinematics_HVBF(5,MomExt,MomDK,applyPSCut,NBin)
   if( Process.eq.61 ) call Kinematics_HJJ(5,MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)

   if (process.eq.60) then
      call EvalAmp_WBFH_UnSymm_SA(MomExt,(/ghz1,ghz2,ghz3,ghz4/),me2)
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



 




 FUNCTION EvalUnWeighted_HJJ(yRnd,genEvt,RES)
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
real(8) :: MomExt(1:4,1:5), PSWgt,MomDK(1:4,1:4)
real(8) :: me2(-5:5,-5:5),bound(0:121)
integer :: i,j,k,ifound,jfound
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


   if( Process.eq.60 ) call Kinematics_HVBF(5,MomExt,MomDK,applyPSCut,NBin)
   if( Process.eq.61 ) call Kinematics_HJJ(5,MomExt,applyPSCut,NBin)
   if( applyPSCut .or. PSWgt.eq.zero ) return
   

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)
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
         if( yRnd(8).gt.bound(k-1) .and. yRnd(8).lt.bound(k)  ) then
            ifound=i; jfound=j;
            goto 1313
         endif
      enddo
   enddo
1313 continue

   if( Process.eq.60 ) then
      call EvalAmp_WBFH_UnSymm_SA(MomExt,(/ghz1,ghz2,ghz3,ghz4/),me2)

      MY_IDUP(1:2)= (/LHA2M_ID(ifound),LHA2M_ID(jfound)/)
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
      if( any(MY_IDUP(1).eq.(/Up_,Chm_/))   .and. any(MY_IDUP(2).eq.(/ADn_,AStr_,ABot_/)) ) ZZ_fusion=.true.
      if( any(MY_IDUP(1).eq.(/AUp_,AChm_/)) .and. any(MY_IDUP(2).eq.(/Dn_,Str_,Bot_/))    ) ZZ_fusion=.true.
      if( any(MY_IDUP(2).eq.(/Up_,Chm_/))   .and. any(MY_IDUP(1).eq.(/ADn_,AStr_,ABot_/)) ) ZZ_fusion=.true.
      if( any(MY_IDUP(2).eq.(/AUp_,AChm_/)) .and. any(MY_IDUP(1).eq.(/Dn_,Str_,Bot_/))    ) ZZ_fusion=.true.

      if( ZZ_Fusion ) then
          if( (MomExt(4,1)*MomExt(4,3).lt.0d0) .and. (MomExt(4,2)*MomExt(4,4).lt.0d0) ) then ! wrong configuration --> swap 3 and 4
             MY_IDUP(3:4)= (/LHA2M_ID(jfound),LHA2M_ID(ifound)/)
             ICOLUP(1:2,4) = ICOLUP(1:2,1)
             ICOLUP(1:2,3) = ICOLUP(1:2,2)
          else! 
             MY_IDUP(3:4)= (/LHA2M_ID(ifound),LHA2M_ID(jfound)/)
             ICOLUP(1:2,3) = ICOLUP(1:2,1)
             ICOLUP(1:2,4) = ICOLUP(1:2,2)
          endif
      else! WW fusion
          if( (MomExt(4,1)*MomExt(4,3).lt.0d0) .and. (MomExt(4,2)*MomExt(4,4).lt.0d0) ) then ! wrong configuration --> swap 3 and 4
             MY_IDUP(3:4)= (/SU2flip(LHA2M_ID(jfound)),SU2flip(LHA2M_ID(ifound))/)
             if( abs(MY_IDUP(3)).eq.Top_ ) MY_IDUP(3) = sign(1,MY_IDUP(3))*Chm_
             if( abs(MY_IDUP(4)).eq.Top_ ) MY_IDUP(4) = sign(1,MY_IDUP(4))*Chm_ 
             ICOLUP(1:2,4) = ICOLUP(1:2,1)
             ICOLUP(1:2,3) = ICOLUP(1:2,2)
          else
             MY_IDUP(3:4)= (/SU2flip(LHA2M_ID(ifound)),SU2flip(LHA2M_ID(jfound))/)
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
      MY_IDUP(1:5)  = (/LHA2M_ID(ifound),LHA2M_ID(jfound),LHA2M_ID(ifound),LHA2M_ID(jfound),Hig_/)

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
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0 ) then! qq->qq
          ICOLUP(1:2,1) = (/501,000/)
          ICOLUP(1:2,2) = (/502,000/)
          ICOLUP(1:2,3) = (/501,000/)
          ICOLUP(1:2,4) = (/502,000/) 
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0 ) then! qbq->qbq
          ICOLUP(1:2,2) = (/501,000/)
          ICOLUP(1:2,1) = (/000,501/)
          ICOLUP(1:2,4) = (/502,000/)
          ICOLUP(1:2,3) = (/000,502/) 
      elseif( MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0 ) then! qbqb->qbqb
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
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt 

   LO_Res_Unpol =  me2(ifound,jfound) * pdf(LHA2M_pdf(ifound),1)*pdf(LHA2M_pdf(jfound),2)
   EvalUnWeighted_HJJ = LO_Res_Unpol * PreFac

   if( ifound.eq.0 .and. jfound.eq.0 ) then
       CS_max = csmax(ifound,jfound) * adj_par 
   else
       CS_max = csmax(ifound,jfound)   
   endif

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
         call WriteOutEvent_HVBF((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,3),MomExt(1:4,4),MomExt(1:4,5)/),MY_IDUP(1:5),ICOLUP(1:2,1:5))
      else
          RejeCounter = RejeCounter + 1
      endif


ELSE! NOT GENEVT


   if( Process.eq.60 ) then
      call EvalAmp_WBFH_UnSymm_SA(MomExt,(/ghz1,ghz2,ghz3,ghz4/),me2)

   elseif( Process.eq.61 ) then
      call EvalAmp_SBFH_UnSymm_SA(MomExt,(/ghg2,ghg3,ghg4/),me2)
      me2 = me2 * (2d0/3d0*alphas**2)**2

   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt 


   LO_Res_Unpol = 0d0
   do i = -5,5
      do j = -5,5

         LO_Res_Unpol = me2(i,j)*pdf(LHA2M_pdf(i),1)*pdf(LHA2M_pdf(j),2) * PreFac

         EvalUnWeighted_HJJ = EvalUnWeighted_HJJ  + LO_Res_Unpol 

          RES(i,j) = LO_Res_Unpol
          if (LO_Res_Unpol.gt.csmax(i,j)) then
              csmax(i,j) = LO_Res_Unpol
          endif
      enddo
   enddo



ENDIF! GENEVT


 RETURN
 END FUNCTION EvalUnWeighted_HJJ














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
    real(8) :: yRnd(1:20),VgsWgt, EvalWeighted_VHiggs
    real(8) :: pdf(-6:6,1:2)
    real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
    real(8) :: MomExt(1:4,1:9), PSWgt, PSWgt2
    real(8) :: me2, lheweight(-6:6,-6:6)
    integer :: i,j,k,NBin(1:NumHistograms),NHisto
    real(8) :: LO_Res_Unpol, PreFac
    logical :: applyPSCut !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!phase space cuts?
    real(8) :: cyRnd(4)
  
    double precision beam_momentum(2,4), four_momentum(7,4),inv_mass(7),mass(7,2)
    double precision helicity(7), beam_h(2) !helicities
    integer id(7), beam_id(2), id2(7)

    EvalWeighted_VHiggs=0d0
    EvalCounter = EvalCounter+1

    mass(1,1)=M_V*1d2
    mass(2,1)=M_V*1d2
    mass(1,2)=Ga_V*1d2
    mass(2,2)=Ga_V*1d2
    mass(3,1)=M_Reso*1d2
    mass(3,2)=Ga_Reso*1d2
    mass(4,1)=0d0
    mass(5,1)=0d0
    mass(6,1)=0d0
    mass(7,1)=0d0
    mass(4,2)=0d0
    mass(5,2)=0d0
    mass(6,2)=0d0
    mass(7,2)=0d0

    id(3)=convertLHE(Hig_)
    id(6)=convertLHE(Bot_)
    id(7)=-id(6)
!toss coin and decide beam A helicity
        if (yRnd(1).lt.(0.5d0+POL_A/200d0))then
          beam_h(1)=1d0
        else
          beam_h(1)=-1d0
        endif
!toss coin and decide beam B helicity
        if (yRnd(2).lt.(0.5d0+POL_B/200d0))then
          beam_h(2)=1d0
        else
          beam_h(2)=-1d0
        endif
!toss coin and decide particle 6,7 helicities
        if (yRnd(4).gt.0.5d0)then
          helicity(6)=1d0
        else
          helicity(6)=-1d0
        endif
        helicity(7)=helicity(6)
!toss coin and decide particle 4,5 helicities
        if (yRnd(3).gt.0.5d0)then
          helicity(4)=1d0
        else
          helicity(4)=-1d0
        endif
        helicity(5)=-helicity(4)

    if(DecayMode1.eq.0)then
      id(1)=convertLHE(Z0_)
      id(2)=convertLHE(Z0_)
      if(yRnd(5).lt.0.5d0)then
        id(4)=convertLHE(MuM_)
        id(5)=-id(4)
      else
        id(4)=convertLHE(ElM_)
        id(5)=-id(4)
      endif

    elseif(DecayMode1.eq.1)then
      id(1)=convertLHE(Z0_)
      id(2)=convertLHE(Z0_)
      if(yRnd(5).lt.Brhadr_Z_uu)then
        id(4)=convertLHE(Up_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Brhadr_Z_uu+Brhadr_Z_cc))then
        id(4)=convertLHE(Chm_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Brhadr_Z_uu+Brhadr_Z_cc+Brhadr_Z_dd))then
        id(4)=convertLHE(Dn_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Brhadr_Z_uu+Brhadr_Z_cc+Brhadr_Z_dd+Brhadr_Z_ss))then
        id(4)=convertLHE(Str_)
        id(5)=-id(4)
      else
        id(4)=convertLHE(Bot_)
        id(5)=-id(4)  
      endif

    elseif(DecayMode1.eq.2)then
      id(1)=convertLHE(Z0_)
      id(2)=convertLHE(Z0_)
      id(4)=convertLHE(TaM_)
      id(5)=-id(4)  

    elseif(DecayMode1.eq.3)then
      id(1)=convertLHE(Z0_)
      id(2)=convertLHE(Z0_)
      if(yRnd(5).lt.0.33333333333333333333d0)then
        id(4)=convertLHE(NuE_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.0.66666666666666666667d0)then
        id(4)=convertLHE(NuM_)
        id(5)=-id(4)
      else
        id(4)=convertLHE(NuT_)
        id(5)=-id(4)
      endif
    helicity(4)=sign(1d0,-dble(id(4)))
    helicity(5)=-helicity(4)

    elseif (DecayMode1.eq.4) then
      id(1)=convertLHE(Wp_)
      id(2)=convertLHE(Wp_)
      if(yRnd(5).lt.0.5d0)then
        id(5)=convertLHE(ElP_)
        id(4)=convertLHE(NuE_)
      else
        id(5)=convertLHE(MuP_)
        id(4)=convertLHE(NuM_)
      endif
      helicity(4)=sign(1d0,-dble(id(4)))
      helicity(5)=-helicity(4)

    elseif(DecayMode1.eq.5)then
      id(1)=convertLHE(Wp_)
      id(2)=convertLHE(Wp_)
      if(yRnd(5).lt.0.5d0)then
      id(4)=convertLHE(Up_)
      id(5)=convertLHE(Adn_)
    else
      id(4)=convertLHE(Chm_)
      id(5)=convertLHE(AStr_)
    endif
      helicity(4)=sign(1d0,-dble(id(4)))
      helicity(5)=-helicity(4)

    elseif(DecayMode1.eq.6)then
      id(1)=convertLHE(Wp_)
      id(2)=convertLHE(Wp_)
      id(5)=convertLHE(TaP_)
      id(4)=convertLHE(NuT_)
      helicity(4)=sign(1d0,-dble(id(4)))
      helicity(5)=-helicity(4)

    elseif(DecayMode1.eq.7)then
      print *, "invalid final states for V > VH"
      stop

    elseif(DecayMode1.eq.8)then
      id(1)=convertLHE(Z0_)
      id(2)=convertLHE(Z0_)
      if(yRnd(5).lt.0.33333333333333333333d0)then
        id(4)=convertLHE(ElM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.0.66666666666666666667d0)then
        id(4)=convertLHE(MuM_)
        id(5)=-id(4)
      else
        id(4)=convertLHE(TaM_)
        id(5)=-id(4)
      endif

    elseif(DecayMode1.eq.9)then
      id(1)=convertLHE(Z0_)
      id(2)=convertLHE(Z0_)
      if(yRnd(5).lt.Br_Z_uu)then
        id(4)=convertLHE(Up_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc))then
        id(4)=convertLHE(Chm_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd))then
        id(4)=convertLHE(Dn_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss))then
        id(4)=convertLHE(Str_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb))then
        id(4)=convertLHE(Bot_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee))then
        id(4)=convertLHE(ElM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm))then
        id(4)=convertLHE(MuM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt))then
        id(4)=convertLHE(TaM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt+Br_Z_nn))then
        id(4)=convertLHE(NuE_)
        id(5)=-id(4)
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)
      elseif(yRnd(5).lt.(Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt+Br_Z_nn+Br_Z_nn))then
        id(4)=convertLHE(NuM_)
        id(5)=-id(4)
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)
      else
        id(4)=convertLHE(NuT_)
        id(5)=-id(4)
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)
      endif

      elseif(DecayMode1.eq.10)then
        id(1)=convertLHE(Wp_)
        id(2)=convertLHE(Wp_)
        if(yRnd(5).lt.0.33333333333333333333d0)then
          id(5)=convertLHE(ElP_)
          id(4)=convertLHE(NuE_)
        elseif(yRnd(5).lt.0.66666666666666666667d0)then
          id(5)=convertLHE(MuP_)
          id(4)=convertLHE(NuM_)
        else
          id(5)=convertLHE(TaP_)
          id(4)=convertLHE(NuT_)
        endif
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)

      elseif(DecayMode1.eq.11)then
        id(1)=convertLHE(Wp_)
        id(2)=convertLHE(Wp_)
        if(yRnd(5).lt.Br_W_en)then
          id(5)=convertLHE(ElP_)
          id(4)=convertLHE(NuE_)
        elseif(yRnd(5).lt.(Br_W_en+Br_W_mn))then
          id(5)=convertLHE(MuP_)
          id(4)=convertLHE(NuM_)
        elseif(yRnd(5).lt.(Br_W_en+Br_W_mn+Br_W_tn))then
          id(5)=convertLHE(TaP_)
          id(4)=convertLHE(NuT_)
        elseif(yRnd(5).lt.(Br_W_en+Br_W_mn+Br_W_tn+Br_W_ud))then
          id(4)=convertLHE(Up_)
          id(5)=convertLHE(Adn_)
        else
          id(4)=convertLHE(Chm_)
          id(5)=convertLHE(AStr_)
        endif
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)

      else
        print *, "invalid final states"
        stop

      endif


if( IsAZDecay(DecayMode1) ) then
!if pp collider
    if(Collider.eq.1)then
      call PDFMapping(14,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      four_momentum(1,1)=EHat*1d2
      four_momentum(1,2)=0d0
      four_momentum(1,3)=0d0
      four_momentum(1,4)=0d0
 
      beam_momentum(1,1)=EHat*1d2/2d0
      beam_momentum(1,2)=0d0
      beam_momentum(1,3)=0d0
      beam_momentum(1,4)=beam_momentum(1,1)
      beam_momentum(2,1)=EHat*1d2/2d0
      beam_momentum(2,2)=0d0
      beam_momentum(2,3)=0d0
      beam_momentum(2,4)=-beam_momentum(2,1)

      call PHASESPACEGEN(yRnd,four_momentum,inv_mass,mass,PSWgt)
      call setPDFs(eta1,eta2,Mu_Fact,pdf)
      FluxFac = 1d0/(2d0*(EHat*1d2)**2)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt/1d4 *6d0 !2 for e and mu, 3 for colors of b
      LO_Res_Unpol=0d0
      EvalWeighted_VHiggs=0d0
      do i = -6,6
        j = -i
        beam_id = (/LHA2M_PDF(i),LHA2M_PDF(j)/)
        if (abs(LHA2M_PDF(i)).ne.6   .and.   abs(LHA2M_PDF(j)).ne.6.  .and.  i.ne.0)then
          call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
        else
          me2=0d0
        endif
          LO_Res_Unpol = me2/3d0*pdf(i,1)*pdf(j,2)* PreFac
          EvalWeighted_VHiggs = EvalWeighted_VHiggs + LO_Res_Unpol 
          lheweight(i,j)=LO_Res_Unpol
      enddo

!if e+ e- collider
    else if(Collider.eq.0)then
      four_momentum(1,1)=ILC_Energy*1d2
      four_momentum(1,2)=0d0
      four_momentum(1,3)=0d0
      four_momentum(1,4)=0d0 

      beam_momentum(1,1)=ILC_Energy*1d2/2d0
      beam_momentum(1,2)=0d0
      beam_momentum(1,3)=0d0
      beam_momentum(1,4)=beam_momentum(1,1)
      beam_momentum(2,1)=ILC_Energy*1d2/2d0
      beam_momentum(2,2)=0d0
      beam_momentum(2,3)=0d0
      beam_momentum(2,4)=-beam_momentum(2,1)

      call PHASESPACEGEN(yRnd,four_momentum,inv_mass,mass,PSWgt)
      FluxFac = 1d0/(2d0*(ILC_Energy*1d2)**2)
      PreFac = fbGeV2 * FluxFac * PSWgt/1d4 *6d0 !2 for e and mu, 3 for colors of b
      LO_Res_Unpol=0d0
      EvalWeighted_VHiggs=0d0
      beam_id(2)=convertLHE(ElM_)
      beam_id(1)=-beam_id(2)
      call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
      LO_Res_Unpol =me2 * PreFac     
      EvalWeighted_VHiggs = EvalWeighted_VHiggs + LO_Res_Unpol
    endif

elseif( IsAWDecay(DecayMode1) ) then
      call PDFMapping(15,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      four_momentum(1,1)=EHat*1d2
      four_momentum(1,2)=0d0
      four_momentum(1,3)=0d0
      four_momentum(1,4)=0d0
 
      beam_momentum(1,1)=EHat*1d2/2d0
      beam_momentum(1,2)=0d0
      beam_momentum(1,3)=0d0
      beam_momentum(1,4)=beam_momentum(1,1)
      beam_momentum(2,1)=EHat*1d2/2d0
      beam_momentum(2,2)=0d0
      beam_momentum(2,3)=0d0
      beam_momentum(2,4)=-beam_momentum(2,1)

      call PHASESPACEGEN(yRnd,four_momentum,inv_mass,mass,PSWgt)
      call setPDFs(eta1,eta2,Mu_Fact,pdf)
      FluxFac = 1d0/(2d0*(EHat*1d2)**2)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt/1d4 *6d0 !2 for e and mu, 3 for colors of b
      LO_Res_Unpol=0d0
      EvalWeighted_VHiggs=0d0
      do i = -6,6
      do j = -6,6
        beam_id = (/LHA2M_PDF(i),LHA2M_PDF(j)/)
        if    ( ((beam_id(1).eq.convertLHE(Up_).or.beam_id(1).eq.convertLHE(Chm_)) .and. &
         (beam_id(2).eq.convertLHE(ADn_) .or. beam_id(2).eq.convertLHE(AStr_) .or. beam_id(2).eq.convertLHE(ABot_))) .or. & 
        ((beam_id(2).eq.convertLHE(Up_).or.beam_id(2).eq.convertLHE(Chm_)) .and. &
         (beam_id(1).eq.convertLHE(ADn_) .or. beam_id(1).eq.convertLHE(AStr_) .or. beam_id(1).eq.convertLHE(ABot_)))   )then
              helicity(4)=sign(1d0,-dble(id(4)))
              helicity(5)=-helicity(4)
              call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
        elseif( ((beam_id(1).eq.convertLHE(AUp_).or.beam_id(1).eq.convertLHE(AChm_)) .and. &
         (beam_id(2).eq.convertLHE(Dn_) .or. beam_id(2).eq.convertLHE(Str_) .or. beam_id(2).eq.convertLHE(Bot_))) .or. & 
        ((beam_id(2).eq.convertLHE(AUp_).or.beam_id(2).eq.convertLHE(AChm_)) .and. &
         (beam_id(1).eq.convertLHE(Dn_) .or. beam_id(1).eq.convertLHE(Str_) .or. beam_id(1).eq.convertLHE(Bot_)))   )then
              id2=id
              id2(2)=-id(2)
              id2(4)=-id(4)
              id2(5)=-id(5)
              helicity(4)=sign(1d0,-dble(id2(4)))
              helicity(5)=-helicity(4)
              call EvalAmp_VHiggs(yRnd,beam_id,id2,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
        else
              me2=0d0
        endif
          LO_Res_Unpol = me2/3d0*pdf(i,1)*pdf(j,2) * PreFac
          EvalWeighted_VHiggs = EvalWeighted_VHiggs+LO_Res_Unpol
          lheweight(i,j)=LO_Res_Unpol
      enddo
      enddo

elseif( IsAPhoton(DecayMode1) ) then
  print *, "invalid process"
  stop
endif


   call Kinematics_VHiggs(beam_momentum,four_momentum,inv_mass,NBin,applyPSCut)
! boost to the lab frame before writing .lhe
   
   do i=1,4
   do j=1,7
     MomExt(i,j)=four_momentum(j,i)
   enddo
   enddo
   do i=1,4
   MomExt(i,8)=beam_momentum(1,i)
   MomExt(i,9)=beam_momentum(2,i)
   enddo


!print *, MomExt(:,6)

   cyRnd(1)=yRnd(9)
   cyRnd(2)=yRnd(8)
!   cyRnd(3)=yRnd(11)
!   cyRnd(4)=yRnd(10)
   call EvalPhasespace_VDecay(MomExt(1:4,2),inv_mass(2),getMass(convertLHEreverse(id(4)))*1d2,getMass(convertLHEreverse(id(5)))*1d2,cyRnd(1:2),MomExt(1:4,4:5),PSWgt2)
!   call EvalPhasespace_VDecay(MomExt(1:4,3),inv_mass(3),getMass(convertLHEreverse(id(6)))*1d2,getMass(convertLHEreverse(id(7)))*1d2,cyRnd(3:4),MomExt(1:4,6:7),PSWgt2)

!print *, MomExt(:,6)
!print *, "end"
!  pause



   if(Collider.eq.1)then
     call boost2Lab(eta1,eta2,9,MomExt(1:4,1:9))
   endif

   do i=1,4
   do j=1,7
     four_momentum(j,i)=MomExt(i,j)
   enddo
   enddo
   do i=1,4
   beam_momentum(1,i)=MomExt(i,8)
   beam_momentum(2,i)=MomExt(i,9)
   enddo

   do i=4,5
     !inv_mass(i)=getMass(convertLHEreverse(id(i)))*1d2
     inv_mass(i)=dsqrt(dabs(four_momentum(i,:).dot.four_momentum(i,:)))
   enddo
!print *, inv_mass(4),inv_mass(5)!,inv_mass(6),inv_mass(7)
   AccepCounter=AccepCounter+1

  if( writeWeightedLHE ) then
! temporary solution enabling parton shower
    if( IsAZDecay(DecayMode1) ) then
      if(Collider.eq.1)then
        do i = -6,6
          j = -i
          beam_id = (/LHA2M_PDF(i),LHA2M_PDF(j)/)
          if (abs(LHA2M_PDF(i)).ne.6   .and.   abs(LHA2M_PDF(j)).ne.6.  .and.  i.ne.0)then
            if(lheweight(i,j).ne.0d0)then
              call WriteOutEvent_VHiggs(beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,EventWeight=lheweight(i,j)*VgsWgt)
            endif
          endif
        enddo
!if e+ e- collider
      else if(Collider.eq.0)then
        if(EvalWeighted_VHiggs.ne.0d0)then
          call WriteOutEvent_VHiggs(beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,EventWeight=EvalWeighted_VHiggs*VgsWgt)
        endif
      endif
    elseif( IsAWDecay(DecayMode1) ) then
      do i = -6,6
      do j = -6,6
       if(lheweight(i,j).ne.0d0)then
        beam_id = (/LHA2M_PDF(i),LHA2M_PDF(j)/)
        if    ( ((beam_id(1).eq.convertLHE(Up_).or.beam_id(1).eq.convertLHE(Chm_)) .and. &
         (beam_id(2).eq.convertLHE(ADn_) .or. beam_id(2).eq.convertLHE(AStr_) .or. beam_id(2).eq.convertLHE(ABot_))) .or. & 
        ((beam_id(2).eq.convertLHE(Up_).or.beam_id(2).eq.convertLHE(Chm_)) .and. &
         (beam_id(1).eq.convertLHE(ADn_) .or. beam_id(1).eq.convertLHE(AStr_) .or. beam_id(1).eq.convertLHE(ABot_)))   )then
              helicity(4)=sign(1d0,-dble(id(4)))
              helicity(5)=-helicity(4)
              call WriteOutEvent_VHiggs(beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,EventWeight=lheweight(i,j)*VgsWgt)
        elseif( ((beam_id(1).eq.convertLHE(AUp_).or.beam_id(1).eq.convertLHE(AChm_)) .and. &
         (beam_id(2).eq.convertLHE(Dn_) .or. beam_id(2).eq.convertLHE(Str_) .or. beam_id(2).eq.convertLHE(Bot_))) .or. & 
        ((beam_id(2).eq.convertLHE(AUp_).or.beam_id(2).eq.convertLHE(AChm_)) .and. &
         (beam_id(1).eq.convertLHE(Dn_) .or. beam_id(1).eq.convertLHE(Str_) .or. beam_id(1).eq.convertLHE(Bot_)))   )then
              id2=id
              id2(2)=-id(2)
              id2(4)=-id(4)
              id2(5)=-id(5)
              helicity(4)=sign(1d0,-dble(id2(4)))
              helicity(5)=-helicity(4)
              call WriteOutEvent_VHiggs(beam_id,id2,beam_h,helicity,beam_momentum,four_momentum,inv_mass,EventWeight=lheweight(i,j)*VgsWgt)
        endif
       endif
      enddo
      enddo

    elseif( IsAPhoton(DecayMode1) ) then
      print *, "invalid final states"
      stop
    endif
! temporary solution enabling parton shower END
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
real(8) :: LO_Res_Unpol, PreFac, CS_max, sumtot
logical :: applyPSCut,genEVT
real(8) :: cyRnd(4)
real(8) :: beam_momentum(2,4), four_momentum(7,4),inv_mass(7),mass(7,2)
real(8) :: helicity(7), beam_h(2) !helicities
integer :: id(7), beam_id(2), id2(7)
include 'csmaxvalue.f'

EvalUnWeighted_VHiggs = 0d0

mass(1,1)=M_V*1d2
mass(2,1)=M_V*1d2
mass(1,2)=Ga_V*1d2
mass(2,2)=Ga_V*1d2
mass(3,1)=M_Reso*1d2
mass(3,2)=Ga_Reso*1d2
mass(4,1)=0d0
mass(5,1)=0d0
mass(6,1)=0d0
mass(7,1)=0d0
mass(4,2)=0d0
mass(5,2)=0d0
mass(6,2)=0d0
mass(7,2)=0d0

id(3)=convertLHE(Hig_)
id(6)=convertLHE(Bot_)
id(7)=-id(6)
!toss coin and decide beam A helicity
if (yRnd(1).lt.(0.5d0+POL_A/200d0))then
  beam_h(1)=1d0
else
  beam_h(1)=-1d0
endif
!toss coin and decide beam B helicity
if (yRnd(2).lt.(0.5d0+POL_B/200d0))then
  beam_h(2)=1d0
else
  beam_h(2)=-1d0
endif
!toss coin and decide particle 6,7 helicities
if (yRnd(4).gt.0.5d0)then
  helicity(6)=1d0
else
  helicity(6)=-1d0
endif
helicity(7)=helicity(6)
!toss coin and decide particle 4,5 helicities
if (yRnd(3).gt.0.5d0)then
  helicity(4)=1d0
else
  helicity(4)=-1d0
endif
helicity(5)=-helicity(4)

if(DecayMode1.eq.0)then
  id(1)=convertLHE(Z0_)
  id(2)=convertLHE(Z0_)
  if(yRnd(5).lt.0.5d0)then
id(4)=convertLHE(MuM_)
id(5)=-id(4)
  else
id(4)=convertLHE(ElM_)
id(5)=-id(4)
endif

elseif(DecayMode1.eq.1)then
  id(1)=convertLHE(Z0_)
  id(2)=convertLHE(Z0_)
  id(4)=convertLHE(ZQuaBranching(yRnd(5)))
  id(5)=-id(4)

elseif(DecayMode1.eq.2)then
  id(1)=convertLHE(Z0_)
  id(2)=convertLHE(Z0_)
  id(4)=convertLHE(TaM_)
  id(5)=-id(4)  

elseif(DecayMode1.eq.3)then
  id(1)=convertLHE(Z0_)
  id(2)=convertLHE(Z0_)
  if(yRnd(5).lt.0.33333333333333333333d0)then
    id(4)=convertLHE(NuE_)
    id(5)=-id(4)
  elseif(yRnd(5).lt.0.66666666666666666667d0)then
    id(4)=convertLHE(NuM_)
    id(5)=-id(4)
  else
    id(4)=convertLHE(NuT_)
    id(5)=-id(4)
  endif
  helicity(4)=sign(1d0,-dble(id(4)))
  helicity(5)=-helicity(4)

elseif (DecayMode1.eq.4) then
  id(1)=convertLHE(Wp_)
  id(2)=convertLHE(Wp_)
  if(yRnd(5).lt.0.5d0)then
    id(5)=convertLHE(ElP_)
    id(4)=convertLHE(NuE_)
  else
    id(5)=convertLHE(MuP_)
    id(4)=convertLHE(NuM_)
  endif
  helicity(4)=sign(1d0,-dble(id(4)))
  helicity(5)=-helicity(4)

elseif(DecayMode1.eq.5)then
  id(1)=convertLHE(Wp_)
  id(2)=convertLHE(Wp_)
  if(yRnd(5).lt.0.5d0)then
    id(4)=convertLHE(Up_)
    id(5)=convertLHE(Adn_)
  else
    id(4)=convertLHE(Chm_)
    id(5)=convertLHE(AStr_)
  endif
  helicity(4)=sign(1d0,-dble(id(4)))
  helicity(5)=-helicity(4)

elseif(DecayMode1.eq.6)then
  id(1)=convertLHE(Wp_)
  id(2)=convertLHE(Wp_)
  id(5)=convertLHE(TaP_)
  id(4)=convertLHE(NuT_)
  helicity(4)=sign(1d0,-dble(id(4)))
  helicity(5)=-helicity(4)

elseif(DecayMode1.eq.7)then
  print *, "invalid final states for V > VH"
  stop

elseif(DecayMode1.eq.8)then
  id(1)=convertLHE(Z0_)
  id(2)=convertLHE(Z0_)
  if(yRnd(5).lt.0.33333333333333333333d0)then
    id(4)=convertLHE(ElM_)
    id(5)=-id(4)
  elseif(yRnd(5).lt.0.66666666666666666667d0)then
    id(4)=convertLHE(MuM_)
    id(5)=-id(4)
  else
    id(4)=convertLHE(TaM_)
    id(5)=-id(4)
  endif

elseif(DecayMode1.eq.9)then
  id(1)=convertLHE(Z0_)
  id(2)=convertLHE(Z0_)
  id(4)=convertLHE(ZAnyBranching(yRnd(5)))
  id(5)=-id(4)
  if(id(4).eq.convertLHE(NuE_) .or. id(4).eq.convertLHE(NuM_) .or. id(4).eq.convertLHE(NuT_))then
    helicity(4)=sign(1d0,-dble(id(4)))
    helicity(5)=-helicity(4)
  endif

elseif(DecayMode1.eq.10)then
  id(1)=convertLHE(Wp_)
  id(2)=convertLHE(Wp_)
  if(yRnd(5).lt.6d0/39d0)then
        id(4)=convertLHE(Up_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(12d0/39d0))then
        id(4)=convertLHE(Chm_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(18d0/39d0))then
        id(4)=convertLHE(Dn_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(24d0/39d0))then
        id(4)=convertLHE(Str_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(30d0/39d0))then
        id(4)=convertLHE(Bot_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(32d0/39d0))then
        id(4)=convertLHE(ElM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(34d0/39d0))then
        id(4)=convertLHE(MuM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(36d0/39d0))then
        id(4)=convertLHE(TaM_)
        id(5)=-id(4)
      elseif(yRnd(5).lt.(37d0/39d0))then
        id(4)=convertLHE(NuE_)
        id(5)=-id(4)
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)
      elseif(yRnd(5).lt.(38d0/39d0))then
        id(4)=convertLHE(NuM_)
        id(5)=-id(4)
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)
      else
        id(4)=convertLHE(NuT_)
        id(5)=-id(4)
        helicity(4)=sign(1d0,-dble(id(4)))
        helicity(5)=-helicity(4)
      endif

elseif(DecayMode1.eq.11)then
  id(1)=convertLHE(Wp_)
  id(2)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/9d0)then
    id(5)=convertLHE(ElP_)
    id(4)=convertLHE(NuE_)
  elseif(yRnd(5).lt.(2d0/9d0))then
    id(5)=convertLHE(MuP_)
    id(4)=convertLHE(NuM_)
  elseif(yRnd(5).lt.(3d0/9d0))then
    id(5)=convertLHE(TaP_)
    id(4)=convertLHE(NuT_)
  elseif(yRnd(5).lt.(6d0/9d0))then
    id(4)=convertLHE(Up_)
    id(5)=convertLHE(Adn_)
  else
    id(4)=convertLHE(Chm_)
    id(5)=convertLHE(AStr_)
  endif
  helicity(4)=sign(1d0,-dble(id(4)))
  helicity(5)=-helicity(4)

else
  print *, "invalid final states"
  stop

endif


if( IsAZDecay(DecayMode1) ) then
!if pp collider
    if(Collider.eq.1)then
      call PDFMapping(14,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      four_momentum(1,1)=EHat*1d2
      four_momentum(1,2)=0d0
      four_momentum(1,3)=0d0
      four_momentum(1,4)=0d0
 
      beam_momentum(1,1)=EHat*1d2/2d0
      beam_momentum(1,2)=0d0
      beam_momentum(1,3)=0d0
      beam_momentum(1,4)=beam_momentum(1,1)
      beam_momentum(2,1)=EHat*1d2/2d0
      beam_momentum(2,2)=0d0
      beam_momentum(2,3)=0d0
      beam_momentum(2,4)=-beam_momentum(2,1)

      call PHASESPACEGEN(yRnd,four_momentum,inv_mass,mass,PSWgt)
      call setPDFs(eta1,eta2,Mu_Fact,pdf)
      FluxFac = 1d0/(2d0*(EHat*1d2)**2)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt/1d4/3d0! *6d0 !2 for e and mu, 3 for colors of b

!if e+ e- collider
    else if(Collider.eq.0)then
      four_momentum(1,1)=ILC_Energy*1d2
      four_momentum(1,2)=0d0
      four_momentum(1,3)=0d0
      four_momentum(1,4)=0d0 

      beam_momentum(1,1)=ILC_Energy*1d2/2d0
      beam_momentum(1,2)=0d0
      beam_momentum(1,3)=0d0
      beam_momentum(1,4)=beam_momentum(1,1)
      beam_momentum(2,1)=ILC_Energy*1d2/2d0
      beam_momentum(2,2)=0d0
      beam_momentum(2,3)=0d0
      beam_momentum(2,4)=-beam_momentum(2,1)

      call PHASESPACEGEN(yRnd,four_momentum,inv_mass,mass,PSWgt)
      FluxFac = 1d0/(2d0*(ILC_Energy*1d2)**2)
      PreFac = fbGeV2 * FluxFac * PSWgt/1d4! *6d0 !2 for e and mu, 3 for colors of b
    endif

elseif( IsAWDecay(DecayMode1) ) then
      call PDFMapping(15,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)

      four_momentum(1,1)=EHat*1d2
      four_momentum(1,2)=0d0
      four_momentum(1,3)=0d0
      four_momentum(1,4)=0d0
 
      beam_momentum(1,1)=EHat*1d2/2d0
      beam_momentum(1,2)=0d0
      beam_momentum(1,3)=0d0
      beam_momentum(1,4)=beam_momentum(1,1)
      beam_momentum(2,1)=EHat*1d2/2d0
      beam_momentum(2,2)=0d0
      beam_momentum(2,3)=0d0
      beam_momentum(2,4)=-beam_momentum(2,1)

      call PHASESPACEGEN(yRnd,four_momentum,inv_mass,mass,PSWgt)
      call setPDFs(eta1,eta2,Mu_Fact,pdf)
      FluxFac = 1d0/(2d0*(EHat*1d2)**2)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt/1d4/3d0! *6d0 !2 for e and mu, 3 for colors of b

elseif( IsAPhoton(DecayMode1) ) then
  print *, "invalid process"
  stop
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


if( IsAZDecay(DecayMode1) ) then
!if pp collider
  if(Collider.eq.1)then
    beam_id = (/ifound,jfound/)

    call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)

    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2) * PreFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol

!if e+ e- collider
  else if(Collider.eq.0)then
    ifound=0
    jfound=0
    beam_id(2)=convertLHE(ElM_)
    beam_id(1)=-beam_id(2)
    call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
    LO_Res_Unpol = me2 * PreFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol
  endif

elseif( IsAWDecay(DecayMode1) ) then
!pp>WH
    beam_id = (/ifound,jfound/)

    if( ((beam_id(1).eq.convertLHE(AUp_).or.beam_id(1).eq.convertLHE(AChm_)) .and. &
     (beam_id(2).eq.convertLHE(Dn_) .or. beam_id(2).eq.convertLHE(Str_) .or. beam_id(2).eq.convertLHE(Bot_))) .or. & 
    ((beam_id(2).eq.convertLHE(AUp_).or.beam_id(2).eq.convertLHE(AChm_)) .and. &
     (beam_id(1).eq.convertLHE(Dn_) .or. beam_id(1).eq.convertLHE(Str_) .or. beam_id(1).eq.convertLHE(Bot_)))   )then
      id(1)=-id(1)
      id(2)=-id(2)
      id(4)=-id(4)
      id(5)=-id(5)
      helicity(4)=sign(1d0,-dble(id(4)))
      helicity(5)=-helicity(4)
    endif
    call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)

    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2) * PreFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol

elseif( IsAPhoton(DecayMode1) ) then
  print *, "invalid process"
  stop
endif

   
  CS_max = csmax(ifound,jfound)
  if( EvalUnWeighted_VHiggs.gt. CS_max) then
    write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_VHiggs, CS_max
    write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_VHiggs, CS_max
    AlertCounter = AlertCounter + 1
    Res = 0d0
  elseif( EvalUnWeighted_VHiggs .gt. yRnd(17)*CS_max ) then

    call Kinematics_VHiggs(beam_momentum,four_momentum,inv_mass,NBin,applyPSCut)

    do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
    enddo
    AccepCounter = AccepCounter + 1
    do i=1,4
    do j=1,7
      MomExt(i,j)=four_momentum(j,i)
    enddo
    enddo
    do i=1,4
      MomExt(i,8)=beam_momentum(1,i)
      MomExt(i,9)=beam_momentum(2,i)
    enddo
    cyRnd(1)=yRnd(9)
    cyRnd(2)=yRnd(8)
    call EvalPhasespace_VDecay(MomExt(1:4,2),inv_mass(2),getMass(convertLHEreverse(id(4)))*1d2,getMass(convertLHEreverse(id(5)))*1d2,cyRnd(1:2),MomExt(1:4,4:5),PSWgt2)
    if(Collider.eq.1)then
      call boost2Lab(eta1,eta2,9,MomExt(1:4,1:9))
    endif

    do i=1,4
    do j=1,7
      four_momentum(j,i)=MomExt(i,j)
    enddo
    enddo
    do i=1,4
    beam_momentum(1,i)=MomExt(i,8)
    beam_momentum(2,i)=MomExt(i,9)
    enddo

    do i=4,5
      inv_mass(i)=dsqrt(dabs(four_momentum(i,:).dot.four_momentum(i,:)))
    enddo

    call WriteOutEvent_VHiggs(beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,EventWeight=1d0)
  
  else
    RejeCounter = RejeCounter + 1
  endif


ELSE! NOT GENEVT

if( IsAZDecay(DecayMode1) ) then
!if pp collider
  if(Collider.eq.1)then
  do i = -5,5
    j = -i
    beam_id = (/i,j/)
    if (abs(i).ne.0)then
      call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
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
    beam_id(2)=convertLHE(ElM_)
    beam_id(1)=-beam_id(2)
    call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
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
    beam_id = (/i,j/)
    if    ( ((beam_id(1).eq.convertLHE(Up_).or.beam_id(1).eq.convertLHE(Chm_)) .and. &
     (beam_id(2).eq.convertLHE(ADn_) .or. beam_id(2).eq.convertLHE(AStr_) .or. beam_id(2).eq.convertLHE(ABot_))) .or. & 
    ((beam_id(2).eq.convertLHE(Up_).or.beam_id(2).eq.convertLHE(Chm_)) .and. &
     (beam_id(1).eq.convertLHE(ADn_) .or. beam_id(1).eq.convertLHE(AStr_) .or. beam_id(1).eq.convertLHE(ABot_)))   )then
      helicity(4)=sign(1d0,-dble(id(4)))
      helicity(5)=-helicity(4)
      call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
    elseif( ((beam_id(1).eq.convertLHE(AUp_).or.beam_id(1).eq.convertLHE(AChm_)) .and. &
     (beam_id(2).eq.convertLHE(Dn_) .or. beam_id(2).eq.convertLHE(Str_) .or. beam_id(2).eq.convertLHE(Bot_))) .or. & 
    ((beam_id(2).eq.convertLHE(AUp_).or.beam_id(2).eq.convertLHE(AChm_)) .and. &
     (beam_id(1).eq.convertLHE(Dn_) .or. beam_id(1).eq.convertLHE(Str_) .or. beam_id(1).eq.convertLHE(Bot_)))   )then
      id(1)=-id(1)
      id(2)=-id(2)
      id(4)=-id(4)
      id(5)=-id(5)
      helicity(4)=sign(1d0,-dble(id(4)))
      helicity(5)=-helicity(4)
      call EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
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

elseif( IsAPhoton(DecayMode1) ) then
  print *, "invalid process"
  stop
endif


ENDIF! GENEVT


 RETURN
 end Function EvalUnWeighted_VHiggs















 FUNCTION EvalWeighted(yRnd,VgsWgt)    ! this is a function which is only for computations
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
 real(8) :: EvalWeighted,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
 real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac,PDFFac1,PDFFac2
 real(8) :: pdf(-6:6,1:2)
 integer :: NBin(1:NumHistograms),NHisto,i,MY_IDUP(1:9), ICOLUP(1:2,1:9),xBin(1:4)
 real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
 real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomDK_massless(1:4,1:4)
 real(8) :: EZ,  EZ1, EZ2, pz, xmax, xx1, xx2
 real(8) :: pz12, MomExt_f(1:4,1:4), MomDK_f(1:4,1:4)
 real(8) :: MZ1,MZ2,ML1,ML2,ML3,ML4,EZ_max,dr, MG, yz1,yz2
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
    MY_IDUP(1:3) = 0
    yz1 = yRnd(10)
    yz2 = yRnd(11)
    offzchannel = yRnd(12) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell 
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

            EZ_max = EHat - MZ1*0.99d0
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
        MZ1 = getMass(MY_IDUP(4))
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
        MZ2 = getMass(MY_IDUP(5))
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
        MZ1 = getMass(MY_IDUP(4))
        MZ2 = getMass(MY_IDUP(5))
  endif

    if( EHat.lt.MZ1+MZ2 ) then
      EvalWeighted = 0d0
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
            EvalWeighted = 0d0
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
            EvalWeighted = 0d0
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
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif
      endif

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2
      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalWeighted = LO_Res_Unpol * PreFac

! EvalWeighted = PreFac  ! for PS output   (only run 1 iteration without vegas adaptation)

    endif

   if (PChannel.eq.1.or.PChannel.eq.2) then
      PDFFac1 = pdf(Up_,1) *pdf(AUp_,2)  + pdf(Dn_,1) *pdf(ADn_,2)   &
              + pdf(Chm_,1)*pdf(AChm_,2) + pdf(Str_,1)*pdf(AStr_,2)  &
              + pdf(Bot_,1)*pdf(ABot_,2)                            
      PDFFac2 = pdf(Up_,2) *pdf(AUp_,1)  + pdf(Dn_,2) *pdf(ADn_,1)   &
              + pdf(Chm_,2)*pdf(AChm_,1) + pdf(Str_,2)*pdf(AStr_,1)  &
              + pdf(Bot_,2)*pdf(ABot_,1)


      if (Process.eq.1) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(13).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif
      endif

      LO_Res_Unpol1 = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2 * PDFFac1
      LO_Res_Unpol2 = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2 * PDFFac2
      LO_Res_Unpol = LO_Res_Unpol1 + LO_Res_Unpol2
      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * SymmFac

      if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalWeighted = LO_Res_Unpol * PreFac
   endif


   if( writeWeightedLHE .and. (.not. warmup) ) then
      if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2))   ) then
            call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=EvalWeighted*VgsWgt)
      else
            call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventWeight=EvalWeighted*VgsWgt)
      endif
   endif


      do NHisto=1,NumHistograms-7
          call intoHisto(NHisto,NBin(NHisto),EvalWeighted*VgsWgt)
      enddo

!       this is for z decays
!       if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(7)).eq.ElP_ .and. abs(MY_IDUP(8)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(12,NBin(12),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.MuP_ .and. abs(MY_IDUP(7)).eq.MuP_ .and. abs(MY_IDUP(8)).eq.MuP_ .and. abs(MY_IDUP(9)).eq.MuP_ ) call intoHisto(13,NBin(13),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(7)).eq.taP_ .and. abs(MY_IDUP(8)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(14,NBin(14),EvalWeighted*VgsWgt)
! 
!       if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(7)).eq.ElP_ .and. abs(MY_IDUP(8)).eq.muP_ .and. abs(MY_IDUP(9)).eq.muP_ ) call intoHisto(15,NBin(15),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.muP_ .and. abs(MY_IDUP(7)).eq.muP_ .and. abs(MY_IDUP(8)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(15,NBin(15),EvalWeighted*VgsWgt)
! 
!       if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(7)).eq.ElP_ .and. abs(MY_IDUP(8)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(16,NBin(16),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(7)).eq.taP_ .and. abs(MY_IDUP(8)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(16,NBin(16),EvalWeighted*VgsWgt)
! 
!       if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(7)).eq.taP_ .and. abs(MY_IDUP(8)).eq.MuP_ .and. abs(MY_IDUP(9)).eq.MuP_ ) call intoHisto(17,NBin(17),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.MuP_ .and. abs(MY_IDUP(7)).eq.MuP_ .and. abs(MY_IDUP(8)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(17,NBin(17),EvalWeighted*VgsWgt)

!       this is for w decays
!       if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(12,NBin(12),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.MuP_ .and. abs(MY_IDUP(9)).eq.MuP_ ) call intoHisto(13,NBin(13),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(14,NBin(14),EvalWeighted*VgsWgt)
! 
!       if( abs(MY_IDUP(6)).eq.ElP_.and. abs(MY_IDUP(9)).eq.muP_) call intoHisto(15,NBin(15),EvalWeighted*VgsWgt)
!       if( abs(MY_IDUP(6)).eq.muP_.and. abs(MY_IDUP(9)).eq.ElP_) call intoHisto(15,NBin(15),EvalWeighted*VgsWgt)


      call intoHisto(18,NBin(18),EvalWeighted*VgsWgt)


      xBin(1) = WhichXBin(1,yrnd(1))
      xBin(2) = WhichXBin(2,yrnd(2))
!      xBin(3) = WhichXBin(3,yrnd(10))
!      xBin(4) = WhichXBin(4,yrnd(11))
      if( EvalWeighted .gt. globalMax ) globalMax = EvalWeighted
      if( EvalWeighted .lt. globalMin ) globalMin = EvalWeighted
      if( EvalWeighted .gt. PartitionMax(xBin(1),xBin(2)) ) PartitionMax(xBin(1),xBin(2)) = EvalWeighted
!       if( EvalWeighted .gt. PartitionMax(xBin(1),xBin(2),xBin(3),xBin(4)) ) PartitionMax(xBin(1),xBin(2),xBin(3),xBin(4)) = EvalWeighted


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
real(8) :: EvalUnWeighted,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumHistograms),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomDK_massless(1:4,1:4),MomExt_f(1:4,1:4),MomDK_f(1:4,1:4)
logical :: applyPSCut,genEvt
real(8) :: CS_max, channel_ratio
real(8) :: oneovervolume, bound(1:11), sumtot,yz1,yz2,EZ_max,dr,MZ1,MZ2,ML1,ML2,ML3,ML4
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
! if(  my_idup(6).eq.my_idup(8)) return! for 2e2mu   ! noi+ 0.92591989,     i+ 0.92775679,   i- 0.14072385,     ix+ 3.2770660  ix-  0.88181581E-02
! if(  my_idup(6).ne.my_idup(8)) return! for 4mu/4e  ! noi+ 0.92608512,     i+ 1.01761060,   i- 0.12915384,     ix+  3.5925481 ix-  0.80721147E-02
!                                                            50/50%              48/52%         52/48%               48%                52%

  yz1 = yRnd(10)
  yz2 = yRnd(11)
  offzchannel = yRnd(15) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

            EZ_max = EHat - MZ1*0.99d0
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
        MZ1 = getMass(MY_IDUP(4))
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
        MZ2 = getMass(MY_IDUP(5))
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
      ML1 = getMass(MY_IDUP(7))
      ML2 = getMass(MY_IDUP(6))
      ML3 = getMass(MY_IDUP(9))
      ML4 = getMass(MY_IDUP(8))
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
        ML1 = getMass(MY_IDUP(7))
        ML2 = getMass(MY_IDUP(6))
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



    if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
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
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
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
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
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

   PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * SymmFac
   EvalUnWeighted = LO_Res_Unpol * PreFac

      if( EvalUnWeighted.gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted .gt. yRnd(14)*CS_max ) then
         do NHisto=1,NumHistograms-7
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo

!       this is for z decays
! 	if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(7)).eq.ElP_ .and. abs(MY_IDUP(8)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(12,1,1d0)
! 	if( abs(MY_IDUP(6)).eq.MuP_ .and. abs(MY_IDUP(7)).eq.MuP_ .and. abs(MY_IDUP(8)).eq.MuP_ .and. abs(MY_IDUP(9)).eq.MuP_ ) call intoHisto(13,1,1d0)
! 	if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(7)).eq.taP_ .and. abs(MY_IDUP(8)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(14,1,1d0)
! 
! 	if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(7)).eq.ElP_ .and. abs(MY_IDUP(8)).eq.muP_ .and. abs(MY_IDUP(9)).eq.muP_ ) call intoHisto(15,1,1d0)
! 	if( abs(MY_IDUP(6)).eq.muP_ .and. abs(MY_IDUP(7)).eq.muP_ .and. abs(MY_IDUP(8)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(15,1,1d0)
! 
! 	if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(7)).eq.ElP_ .and. abs(MY_IDUP(8)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(16,1,1d0)
! 	if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(7)).eq.taP_ .and. abs(MY_IDUP(8)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(16,1,1d0)
! 
! 	if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(7)).eq.taP_ .and. abs(MY_IDUP(8)).eq.MuP_ .and. abs(MY_IDUP(9)).eq.MuP_ ) call intoHisto(17,1,1d0)
! 	if( abs(MY_IDUP(6)).eq.MuP_ .and. abs(MY_IDUP(7)).eq.MuP_ .and. abs(MY_IDUP(8)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(17,1,1d0)
!       this is for w decays
!       if( abs(MY_IDUP(6)).eq.ElP_ .and. abs(MY_IDUP(9)).eq.ElP_ ) call intoHisto(12,1,1d0)
!       if( abs(MY_IDUP(6)).eq.MuP_ .and. abs(MY_IDUP(9)).eq.MuP_ ) call intoHisto(13,1,1d0)
!       if( abs(MY_IDUP(6)).eq.taP_ .and. abs(MY_IDUP(9)).eq.taP_ ) call intoHisto(14,1,1d0)
! 
!       if( abs(MY_IDUP(6)).eq.ElP_.and. abs(MY_IDUP(9)).eq.muP_) call intoHisto(15,1,1d0)
!       if( abs(MY_IDUP(6)).eq.muP_.and. abs(MY_IDUP(9)).eq.ElP_) call intoHisto(15,1,1d0)

      call intoHisto(18,NBin(18),1d0)


! debugcounter(0)=debugcounter(0)+1
! if(  abs(MY_IDUP(6)).ge.7 .and.  abs(MY_IDUP(6)).le.16  ) debugcounter(1)=debugcounter(1)+1
! if(  abs(MY_IDUP(8)).ge.7 .and.  abs(MY_IDUP(8)).le.16  ) debugcounter(2)=debugcounter(2)+1
! if(  (abs(MY_IDUP(8)).ge.7 .and.  abs(MY_IDUP(8)).le.16) .and.  (abs(MY_IDUP(6)).ge.7 .and.  abs(MY_IDUP(6)).le.16)  ) debugcounter(3)=debugcounter(3)+1

! if( (MY_IDUP(6)).eq.(MY_IDUP(8)) ) debugcounter(1)=debugcounter(1)+1

         AccepCounter = AccepCounter + 1
         AccepCounter_part = AccepCounter_part  + parton
         if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
              call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
          else
              call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
         endif
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
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif
      else
         LO_Res_Unpol = 0d0           
      endif
      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * GluonColAvg**2


      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
      EvalUnWeighted = LO_Res_Unpol * PreFac
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
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif
      else
         LO_Res_Unpol1 = 0d0  
         LO_Res_Unpol2 = 0d0  
      endif

      LO_Res_Unpol1 = LO_Res_Unpol1 * SpinAvg * QuarkColAvg**2
      LO_Res_Unpol2 = LO_Res_Unpol2 * SpinAvg * QuarkColAvg**2
!       if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
!       if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt *   SymmFac

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











FUNCTION EvalUnWeighted_BETA(yRnd)
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
real(8) :: EvalUnWeighted_BETA,LO_Res_Unpol,LO_Res_Unpol1,LO_Res_Unpol2,yRnd(1:22),VgsWgt
real(8) :: eta1,eta2,tau,x1,x2,sHatJacobi,PreFac,FluxFac,PDFFac
real(8) :: pdf(-6:6,1:2)
integer :: NBin(1:NumHistograms),NHisto,i,XBin(1:4)
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomExt_f(1:4,1:4),MomDK_f(1:4,1:4),MomDK_massless(1:4,1:4)
logical :: applyPSCut,genEvt
real(8) :: CS_max, channel_ratio
real(8) :: oneovervolume, bound(1:11), sumtot,yz1,yz2,EZ_max,dr,MZ1,MZ2,ML1,ML2,ML3,ML4
integer :: parton(-5:5,-5:5), i1, ifound, i2, MY_IDUP(1:9), ICOLUP(1:2,1:9)
real(8)::ntRnd,ZMass(1:2)
real(8) :: offzchannel
include 'vegas_common.f'
include 'csmaxvalue.f'


   parton = 0
   oneovervolume = one
   ICOLUP(1:2,1:9) = 0
   EvalUnWeighted_BETA = 0d0
   RES = 0d0

   if(OffShellReson) then
      call PDFMapping(10,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   else
      call PDFMapping(12,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
   endif
   EvalCounter = EvalCounter+1

   MY_IDUP(3)= 0

   xBin(1) = WhichXBin(1,yrnd(1))
   xBin(2) = WhichXBin(2,yrnd(2))
!    xBin(3) = WhichXBin(3,yrnd(10))
!    xBin(4) = WhichXBin(4,yrnd(11))
   if( yRnd(14)*globalMax .gt. PartitionMax(xBin(1),xBin(2)) ) then! reject
!    if( yRnd(14)*globalMax .gt. PartitionMax(xBin(1),xBin(2),xBin(3),xBin(4)) ) then! reject
          RejeCounter = RejeCounter + 1
!print *, yRnd(14)*globalMax,globalmax,PartitionMax(xBin(1),xBin(2)),xBin(1),xBin(2);pause
          return
   endif




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
      EvalUnWeighted_BETA = 0d0
      RejeCounter = RejeCounter + 1
      return
    endif




   call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   call boost2Lab(eta1,eta2,4,MomExt(1:4,1:4))
   if( .not. IsAPhoton(DecayMode1) ) then ! don't decay the photon
      ML1 = getMass(MY_IDUP(7))
      ML2 = getMass(MY_IDUP(6))
      ML3 = getMass(MY_IDUP(9))
      ML4 = getMass(MY_IDUP(8))
      if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
          EvalUnWeighted_BETA = 0d0
          RejeCounter = RejeCounter + 1
          return
      endif
      call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,ML1,ML2,yRnd(5:6),MomDK(1:4,1:2),PSWgt2)
      call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,ML3,ML4,yRnd(7:8),MomDK(1:4,3:4),PSWgt3)
      PSWgt = PSWgt * PSWgt2*PSWgt3
    else
        ML1=0d0; ML2=0d0; ML3=0d0; ML4=0d0
        MomDK(1:4,1) = MomExt(1:4,3)
        MomDK(1:4,2) = 0d0
        MomDK(1:4,3) = MomExt(1:4,4)
        MomDK(1:4,4) = 0d0
   endif

    if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then! introduce this momentum flip to allow proper mapping of integrand with Z-poles at MZ2=(p2+p3)^2 and MZ2=(p1+p4)^2
       if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK(1:4,1),MomDK(1:4,3) )
       !             PSWgt = PSWgt * 2d0
    endif

    if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode1)) ) then
        call Kinematics(4,MomExt,MomDK,applyPSCut,NBin)
    else
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalUnWeighted_BETA = 0d0
      return
    endif

   call setPDFs(eta1,eta2,Mu_Fact,pdf)
   FluxFac = 1d0/(2d0*EHat**2)


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


  ifound=1 !!!!!   FIX ifound to one --> this subroutine works only for gg initial states!!

   if(ifound.eq.1 ) then
      parton(0,0) = 1
      CS_max = csmax(0,0)*adj_par   ! this is necessary here for the follow up
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
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_H_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_gg_G_VV( (/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif
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
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_Zprime_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif

      elseif(Process.eq.2) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            else
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,1),-MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol1)
               call EvalAmp_qqb_G_VV((/-MomExt(1:4,2),-MomExt(1:4,1),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol2)
            endif
      endif

      LO_Res_Unpol = LO_Res_Unpol * SpinAvg * QuarkColAvg**2
   endif

   PreFac = 2d0 * fbGeV2 * FluxFac * sHatJacobi * PSWgt * PDFFac * SymmFac
   EvalUnWeighted_BETA = LO_Res_Unpol * PreFac

      if( EvalUnWeighted_BETA .gt. globalMax) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_BETA, globalmax
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_BETA, globalmax
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted_BETA .gt. yRnd(14)*globalMax ) then
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1
         AccepCounter_part = AccepCounter_part  + parton
         if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode1)) ) then
              call WriteOutEvent((/MomExt(1:4,1),MomExt(1:4,2),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
          else
              call WriteOutEvent((/MomExt_f(1:4,1),MomExt_f(1:4,2),MomDK_f(1:4,1),MomDK_f(1:4,2),MomDK_f(1:4,3),MomDK_f(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9))
         endif
      else
          RejeCounter = RejeCounter + 1
      endif




RETURN
END FUNCTION







FUNCTION EvalUnWeighted_withoutProduction(yRnd,genEvt,Ehat,Res,AcceptedEvent,MY_IDUP,ICOLUP)
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
real(8) :: Res
real(8) :: EvalUnWeighted_withoutProduction,LO_Res_Unpol,yRnd(1:22),VgsWgt,LO_Res_Unpol1,LO_Res_Unpol2
real(8) :: tau,x1,x2,sHatJacobi,PreFac
integer :: NBin(1:NumHistograms),NHisto,i
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3
real(8) :: MomExt(1:4,1:4),MomDK(1:4,1:4),MomExt_f(1:4,1:4),MomDK_f(1:4,1:4),MomDK_massless(1:4,1:4)
logical :: applyPSCut,genEvt
real(8) :: CS_max,eta1,eta2
real(8) :: oneovervolume, bound(1:11), sumtot,yz1,yz2,EZ_max,dr,MZ1,MZ2,ML1,ML2,ML3,ML4
integer :: i1, ifound, i2, MY_IDUP(1:9), ICOLUP(1:2,1:9)
real(8)::  ntRnd,ZMass(1:2),AcceptedEvent(1:4,1:4)
real(8) :: offzchannel
include 'vegas_common.f'
include 'csmaxvalue.f'


   oneovervolume = one
   ICOLUP(1:2,1:9) = 0
   EvalUnWeighted_withoutProduction = 0d0
   Res = 0d0
   EvalCounter = EvalCounter+1



   MY_IDUP(3)= Hig_
!    particle associations:
!    
!    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor
   call VVBranchings(MY_IDUP(4:9),ICOLUP(1:2,6:9))


  eta1=1d0; eta2=1d0
  sHatJacobi = 1d0

  yz1 = yRnd(10)
  yz2 = yRnd(11)
  offzchannel = yRnd(15) ! variable to decide which Z is ``on''- and which Z is off- the mass-shell
  if ((OffShellV1.eqv..true.).and.(OffShellV2.eqv..true.)) then
        if(M_Reso.gt.2d0*M_V) then
            EZ_max = EHat
            dr = datan((EZ_max**2-M_V**2)/(Ga_V*M_V)) + datan(M_V/Ga_V)
            MZ1 = dsqrt( M_V*Ga_V * dtan(dr*yz1-datan(M_V/Ga_V)) + M_V**2 )
            sHatJacobi = sHatJacobi * dr/(Ga_V*M_V) * ( (MZ1**2 - M_V**2)**2 + M_V**2*Ga_V**2 )

            EZ_max = EHat - MZ1*0.99d0
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
        MZ1 = getMass(MY_IDUP(4))
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
        MZ2 = getMass(MY_IDUP(5))
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
        MZ1 = getMass(MY_IDUP(4))
        MZ2 = getMass(MY_IDUP(5))
  endif


    if( MZ1+MZ2.gt.EHat ) then
      EvalUnWeighted_withoutProduction = 0d0
      RejeCounter = RejeCounter + 1
      return
    endif




   call EvalPhaseSpace_2to2(EHat,(/MZ1,MZ2/),yRnd(3:4),MomExt(1:4,1:4),PSWgt)
   if( .not.IsAPhoton(DecayMode1) .and. .not.IsAPhoton(DecayMode2) ) then ! don't decay the photon
      ML1 = getMass(MY_IDUP(7))
      ML2 = getMass(MY_IDUP(6))
      ML3 = getMass(MY_IDUP(9))
      ML4 = getMass(MY_IDUP(8))
      if( (MZ1.lt.ML1+ML2) .or. (MZ2.lt.ML3+ML4) ) then
          EvalUnWeighted_withoutProduction = 0d0
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
        ML1 = getMass(MY_IDUP(7))
        ML2 = getMass(MY_IDUP(6))
        ML3=0d0; ML4=0d0
        if( (MZ1.lt.ML1+ML2) ) then
            EvalUnWeighted_withoutProduction = 0d0
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
        call AdjustKinematics(eta1,eta2,MomExt,MomDK,yRnd(9),yRnd(10),yRnd(11),MomExt_f,MomDK_f)
        call Kinematics(4,MomExt_f,MomDK_f,applyPSCut,NBin)
    endif
    if( applyPSCut ) then
      EvalUnWeighted_withoutProduction = 0d0
      return
    endif



IF( GENEVT ) THEN


      MY_IDUP(1:2)=(/Glu_,Glu_/)
      ICOLUP(1:2,1) = (/501,502/)
      ICOLUP(1:2,2) = (/502,501/)

      if (Process.eq.0) then
            if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            else
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
            endif
      endif


      PreFac = 2d0 * fbGeV2 * sHatJacobi * PSWgt * SymmFac
!       if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
!       if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalUnWeighted_withoutProduction = LO_Res_Unpol * PreFac

      CS_max = csmax(0,0)

      if( EvalUnWeighted_withoutProduction .gt. CS_max) then
          write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_withoutProduction, CS_max
          write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_withoutProduction, CS_max
          AlertCounter = AlertCounter + 1
          Res = 0d0

      elseif( EvalUnWeighted_withoutProduction .gt. yRnd(14)*CS_max ) then
         do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
         enddo
         AccepCounter = AccepCounter + 1
         if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode1)) ) then
              AcceptedEvent(1:4,1:4) = MomDK(1:4,1:4)
          else
              AcceptedEvent(1:4,1:4) = MomDK_f(1:4,1:4)
         endif
         Res = 1d0

      else
          RejeCounter = RejeCounter + 1
          Res = 0d0
      endif


ELSE! NOT GENEVT


      if (Process.eq.0) then
         MomExt(1:4,1) = (/M_Reso,0d0,0d0,0d0/)
         MomExt(1:4,2) = (/0d0,0d0,0d0,0d0/)

         if( ML1.gt.1d-6 .or. ML2.gt.1d-6 .or. ML3.gt.1d-6 .or. ML4.gt.1d-6 ) then
               call EvalPhasespace_VDecay(MomExt(1:4,3),MZ1,0d0,0d0,yRnd(5:6),MomDK_massless(1:4,1:2),PSWgt2)
               call EvalPhasespace_VDecay(MomExt(1:4,4),MZ2,0d0,0d0,yRnd(7:8),MomDK_massless(1:4,3:4),PSWgt3)
               if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                  if( yrnd(16).gt.0.5d0 ) call swapmom( MomDK_massless(1:4,1),MomDK_massless(1:4,3) )
               endif
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK_massless(1:4,1),MomDK_massless(1:4,2),MomDK_massless(1:4,3),MomDK_massless(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
         else
               call EvalAmp_H_VV( (/-MomExt(1:4,1)-MomExt(1:4,2),(/0d0,0d0,0d0,0d0/),MomDK(1:4,1),MomDK(1:4,2),MomDK(1:4,3),MomDK(1:4,4)/),MY_IDUP(6:9),LO_Res_Unpol)
         endif
      endif

     PreFac = 2d0 * fbGeV2 * sHatJacobi * PSWgt * SymmFac
!       if( abs(MY_IDUP(6)).ge.1 .and. abs(MY_IDUP(6)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
!       if( abs(MY_IDUP(8)).ge.1 .and. abs(MY_IDUP(8)).le.6 ) PreFac = PreFac * 3d0 ! =Nc
      EvalUnWeighted_withoutProduction = LO_Res_Unpol * PreFac
      Res = EvalUnWeighted_withoutProduction


      if (EvalUnWeighted_withoutProduction.gt.csmax(0,0)) then
          csmax(0,0) = EvalUnWeighted_withoutProduction
      endif


ENDIF! genEvt


RETURN
END FUNCTION






END MODULE ModCrossSection







