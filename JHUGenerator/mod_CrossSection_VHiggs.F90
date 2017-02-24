MODULE ModCrossSection_VHiggs
implicit none
integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

 CONTAINS





Function EvalWeighted_VHiggs(yRnd,VgsWgt)
 use ModKinematics
 use ModKinematics_VHiggs
 use ModParameters
 use ModVHiggs
 use ModMisc
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
    real(8) :: MomExt(1:4,1:9), PSWgt, PSWgt2
    real(8) :: me2!, lheweight(-6:6,-6:6)
    integer :: i,j,k,NBin(1:NumHistograms),NHisto
    real(8) :: LO_Res_Unpol, PreFac, PostFac
    logical :: applyPSCut !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!phase space cuts?
    real(8) :: cyRnd(4)

    real(8) :: inv_mass(9),mass(9,2)
    !double precision beam_momentum(2,4), four_momentum(7,4),inv_mass(7),mass(7,2)
    real(8) :: helicity(9)!, beam_h(2) !helicities
    integer id(9), id2(9)!, beam_id(2)
    logical :: PhoOnshell

    EvalWeighted_VHiggs=0d0
    EvalCounter = EvalCounter+1

    id(:)=0
    helicity(:)=0
    mass(1:2,1:2)=0d0
    if(IsAPhoton(DecayMode1))then
       mass(3,1)=M_Z
       mass(3,2)=Ga_Z
       mass(4,1)=getMass(Pho_)
       mass(4,2)=getDecayWidth(Pho_)
    else
       mass(3,1)=M_V
       mass(3,2)=Ga_V
       mass(4,1)=M_V
       mass(4,2)=Ga_V
    endif
    mass(5,1)=M_Reso
    mass(5,2)=Ga_Reso
    mass(6:9,1:2)=0d0


    id(5)=convertLHE(Hig_)
    if(HbbDecays.eqv..true.)then
      id(8)=convertLHE(Bot_)
      id(9)=-id(8)
    else
      id(8)=Not_a_particle_
      id(9)=Not_a_particle_
    endif
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

if(DecayMode1.eq.0)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.0.5d0)then
    id(6)=convertLHE(MuM_)
  else
    id(6)=convertLHE(ElM_)
  endif
  id(7)=-id(6)
  PostFac=4d0 ! of which 2 is for summing Z > ff~ helicities
!elseif(DecayMode1.eq.1)then
!  id(3)=convertLHE(Z0_)
!  id(4)=convertLHE(Z0_)
!  id(6)=convertLHE(ZQuaBranching_flat(yRnd(5)))
!  id(7)=-id(6)
  
elseif(DecayMode1.eq.1)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.0.2d0)then
    id(6)=convertLHE(Up_)
  elseif(yRnd(5).lt.0.4d0)then
    id(6)=convertLHE(Chm_)
  elseif(yRnd(5).lt.0.6d0)then
    id(6)=convertLHE(Dn_)
  elseif(yRnd(5).lt.0.8d0)then
    id(6)=convertLHE(Str_)
  else
    id(6)=convertLHE(Bot_)
  endif
  id(7)=-id(6)
  PostFac=30d0 ! of which 2 is for summing Z > ff~ helicities and 3 for summing colors

elseif(DecayMode1.eq.2)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  id(6)=convertLHE(TaM_)
  id(7)=-id(6)
  PostFac=2d0  ! of which 2 is for summing Z > ff~ helicities
  
elseif(DecayMode1.eq.3)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.1d0/3d0)then
    id(6)=convertLHE(NuE_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(6)=convertLHE(NuM_)
  else
    id(6)=convertLHE(NuT_)
  endif
  id(7)=-id(6)
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=3d0
  
elseif (DecayMode1.eq.4) then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.0.5d0)then
    id(7)=convertLHE(ElP_)
    id(6)=convertLHE(NuE_)
  else
    id(7)=convertLHE(MuP_)
    id(6)=convertLHE(NuM_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=2d0 
  
elseif(DecayMode1.eq.5)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/6d0)then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.1d0/3d0)then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.0.5d0)then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.5d0/6d0)then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Abot_)
  else
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Abot_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=18d0 ! of which 3 is for summing Z > ff~ colors
  
elseif(DecayMode1.eq.6)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  id(7)=convertLHE(TaP_)
  id(6)=convertLHE(NuT_)
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=1d0

elseif(DecayMode1.eq.7)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Pho_)
  id(6)=convertLHE(Pho_)
  id(7)=Not_a_particle_
  PostFac=2d0  ! of which 2 is for summing photon helicities
  
elseif(DecayMode1.eq.8)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.1d0/3d0)then
    id(6)=convertLHE(ElM_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(6)=convertLHE(MuM_)
  else
    id(6)=convertLHE(TaM_)
  endif
  id(7)=-id(6)
  PostFac=6d0 ! of which 2 is for summing Z > ff~ helicities
  
elseif(DecayMode1.eq.9)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.6d0/39d0)then
    id(6)=convertLHE(Up_)
  elseif(yRnd(5).lt.(12d0/39d0))then
    id(6)=convertLHE(Chm_)
  elseif(yRnd(5).lt.(18d0/39d0))then
    id(6)=convertLHE(Dn_)
  elseif(yRnd(5).lt.(24d0/39d0))then
    id(6)=convertLHE(Str_)
  elseif(yRnd(5).lt.(30d0/39d0))then
    id(6)=convertLHE(Bot_)
  elseif(yRnd(5).lt.(32d0/39d0))then
    id(6)=convertLHE(ElM_)
  elseif(yRnd(5).lt.(34d0/39d0))then
    id(6)=convertLHE(MuM_)
  elseif(yRnd(5).lt.(36d0/39d0))then
    id(6)=convertLHE(TaM_)
  elseif(yRnd(5).lt.(37d0/39d0))then
    id(6)=convertLHE(NuE_)
    helicity(6)=sign(1d0,-dble(id(6)))
    helicity(7)=-helicity(6)
  elseif(yRnd(5).lt.(38d0/39d0))then
    id(6)=convertLHE(NuM_)
    helicity(6)=sign(1d0,-dble(id(6)))
    helicity(7)=-helicity(6)
  else
    id(6)=convertLHE(NuT_)
    helicity(6)=sign(1d0,-dble(id(6)))
    helicity(7)=-helicity(6)
  endif
  id(7)=-id(6)
  PostFac=39d0
  
elseif(DecayMode1.eq.10)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/3d0)then
    id(7)=convertLHE(ElP_)
    id(6)=convertLHE(NuE_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(7)=convertLHE(MuP_)
    id(6)=convertLHE(NuM_)
  else
    id(7)=convertLHE(TaP_)
    id(6)=convertLHE(NuT_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=3d0
  
elseif(DecayMode1.eq.11)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/21d0)then
    id(7)=convertLHE(ElP_)
    id(6)=convertLHE(NuE_)
  elseif(yRnd(5).lt.(2d0/21d0))then
    id(7)=convertLHE(MuP_)
    id(6)=convertLHE(NuM_)
  elseif(yRnd(5).lt.(3d0/21d0))then
    id(7)=convertLHE(TaP_)
    id(6)=convertLHE(NuT_)
  elseif(yRnd(5).lt.(6d0/21d0))then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.(9d0/21d0))then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.(12d0/21d0))then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.(15d0/21d0))then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.(18d0/21d0))then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Abot_)
  else
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Abot_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=21d0

else
  print *, "invalid final states"
  stop

endif

if(HbbDecays) PostFac=PostFac*6 !of which 2 is for summing H > bb~ helicities and 3 for summing bb~ colors


if( IsAZDecay(DecayMode1).or.IsAPhoton(DecayMode1) ) then
!if pp collider
  if(Collider.eq.1.or.Collider.eq.2)then

    if( IsAZDecay(DecayMode1) )then
      call PDFMapping(14,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
    elseif( IsAPhoton(DecayMode1) )then
      call PDFMapping(16,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
    endif      
    if( Ehat .ge. 2d0*m_top ) return !without complex top mass in ZH amplitude
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

    call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1),ZAinterference=includeGammaStar)
    call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))

    if( applyPSCut .or. PSWgt.eq.zero ) return
    if( IsAZDecay(DecayMode1) )then
      call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
    elseif( IsAPhoton(DecayMode1) )then
      call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
    endif

    call setPDFs(eta1,eta2,pdf)
    FluxFac = 1d0/(2d0*EHat**2)
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt
    LO_Res_Unpol=0d0
    EvalWeighted_VHiggs=0d0
    if((VHiggs_PC.eq."qq".or.VHiggs_PC.eq."lo").and.PChannel.eq.1)then
      do i = -6,6
        j = -i
        id(1:2) = (/LHA2M_PDF(i),LHA2M_PDF(j)/)
        if (abs(LHA2M_PDF(i)).ne.6   .and.   abs(LHA2M_PDF(j)).ne.6.  .and.  i.ne.0)then
          call EvalAmp_VHiggs(id,helicity,MomExt,me2)
        else
          me2=0d0
        endif
          LO_Res_Unpol = me2*pdf(i,1)*pdf(j,2)* PreFac *PostFac
          EvalWeighted_VHiggs = EvalWeighted_VHiggs + LO_Res_Unpol
      enddo
    elseif((VHiggs_PC.eq."gg".or.VHiggs_PC.eq."bo".and.VHiggs_PC.eq."tr").or.PChannel.eq.0)then
      if( Ehat .ge. 2d0*m_top ) return !without complex top mass in ZH amplitude
      id(1)=0
      id(2)=0
      i=0
      j=0
      call EvalAmp_VHiggs(id,helicity,MomExt,me2)
      LO_Res_Unpol = me2*pdf(i,1)*pdf(j,2)* PreFac *PostFac
    !  print *, me2,pdf(i,1),pdf(j,2)
      EvalWeighted_VHiggs = EvalWeighted_VHiggs + LO_Res_Unpol
    endif

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
  
    call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1),ZAinterference=includeGammaStar)
    call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
  
    if( applyPSCut .or. PSWgt.eq.zero ) return
  
    FluxFac = 1d0/(2d0*ILC_Energy**2)
    PreFac = fbGeV2 * FluxFac * PSWgt
    LO_Res_Unpol=0d0
    EvalWeighted_VHiggs=0d0
    id(2)=convertLHE(ElM_)
    id(1)=-id(2)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    LO_Res_Unpol =me2 * PreFac *PostFac
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

  call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HbbDecays)
  call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays)
  if( applyPSCut .or. PSWgt.eq.zero ) return
  call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
  call setPDFs(eta1,eta2,pdf)

  FluxFac = 1d0/(2d0*EHat**2)
  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt ! *6d0 !2 for e and mu, 3 for colors of b
  LO_Res_Unpol=0d0
  EvalWeighted_VHiggs=0d0
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
     else
         me2=0d0
     endif
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac *PostFac
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
 use ModKinematics_VHiggs
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
real(8) :: LO_Res_Unpol, PreFac, CS_max, sumtot, PostFac
logical :: applyPSCut,genEVT
real(8) :: cyRnd(4)
real(8) :: inv_mass(9),mass(9,2)
!real(8) :: beam_momentum(2,4), four_momentum(7,4),inv_mass(7),mass(7,2)
real(8) :: helicity(9) !helicities
integer :: id(9), id2(9)
logical :: PhoOnshell
include 'csmaxvalue.f'

EvalUnWeighted_VHiggs = 0d0
id(:)=0
helicity(:)=0

mass(1:2,1:2)=0d0
if(IsAPhoton(DecayMode1))then
  mass(3,1)=M_Z
  mass(3,2)=Ga_Z
  mass(4,1)=getMass(Pho_)
  mass(4,2)=getDecayWidth(Pho_)
else
  mass(3,1)=M_V
  mass(3,2)=Ga_V
  mass(4,1)=M_V
  mass(4,2)=Ga_V
endif
mass(5,1)=M_Reso
mass(5,2)=Ga_Reso
mass(6:9,1:2)=0d0

id(5)=convertLHE(Hig_)
if(HbbDecays.eqv..true.)then
  id(8)=convertLHE(Bot_)
  id(9)=-id(8)
else
  id(8)=Not_a_particle_
  id(9)=Not_a_particle_
endif
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

if(DecayMode1.eq.0)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.0.5d0)then
    id(6)=convertLHE(MuM_)
  else
    id(6)=convertLHE(ElM_)
  endif
  id(7)=-id(6)
  PostFac=4d0 ! of which 2 is for summing Z > ff~ helicities
!elseif(DecayMode1.eq.1)then
!  id(3)=convertLHE(Z0_)
!  id(4)=convertLHE(Z0_)
!  id(6)=convertLHE(ZQuaBranching_flat(yRnd(5)))
!  id(7)=-id(6)
  
elseif(DecayMode1.eq.1)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.0.2d0)then
    id(6)=convertLHE(Up_)
  elseif(yRnd(5).lt.0.4d0)then
    id(6)=convertLHE(Chm_)
  elseif(yRnd(5).lt.0.6d0)then
    id(6)=convertLHE(Dn_)
  elseif(yRnd(5).lt.0.8d0)then
    id(6)=convertLHE(Str_)
  else
    id(6)=convertLHE(Bot_)
  endif
  id(7)=-id(6)
  PostFac=30d0 ! of which 2 is for summing Z > ff~ helicities and 3 for summing colors

elseif(DecayMode1.eq.2)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  id(6)=convertLHE(TaM_)
  id(7)=-id(6)
  PostFac=2d0  ! of which 2 is for summing Z > ff~ helicities
  
elseif(DecayMode1.eq.3)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.1d0/3d0)then
    id(6)=convertLHE(NuE_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(6)=convertLHE(NuM_)
  else
    id(6)=convertLHE(NuT_)
  endif
  id(7)=-id(6)
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=3d0
  
elseif (DecayMode1.eq.4) then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.0.5d0)then
    id(7)=convertLHE(ElP_)
    id(6)=convertLHE(NuE_)
  else
    id(7)=convertLHE(MuP_)
    id(6)=convertLHE(NuM_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=2d0 
  
elseif(DecayMode1.eq.5)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/6d0)then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.1d0/3d0)then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.0.5d0)then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.5d0/6d0)then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Abot_)
  else
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Abot_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=18d0 ! of which 3 is for summing Z > ff~ colors
  
elseif(DecayMode1.eq.6)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  id(7)=convertLHE(TaP_)
  id(6)=convertLHE(NuT_)
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=1d0

elseif(DecayMode1.eq.7)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Pho_)
  id(6)=convertLHE(Pho_)
  id(7)=Not_a_particle_
  PostFac=2d0  ! of which 2 is for summing photon helicities
  
elseif(DecayMode1.eq.8)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.1d0/3d0)then
    id(6)=convertLHE(ElM_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(6)=convertLHE(MuM_)
  else
    id(6)=convertLHE(TaM_)
  endif
  id(7)=-id(6)
  PostFac=6d0 ! of which 2 is for summing Z > ff~ helicities
  
elseif(DecayMode1.eq.9)then
  id(3)=convertLHE(Z0_)
  id(4)=convertLHE(Z0_)
  if(yRnd(5).lt.6d0/39d0)then
    id(6)=convertLHE(Up_)
  elseif(yRnd(5).lt.(12d0/39d0))then
    id(6)=convertLHE(Chm_)
  elseif(yRnd(5).lt.(18d0/39d0))then
    id(6)=convertLHE(Dn_)
  elseif(yRnd(5).lt.(24d0/39d0))then
    id(6)=convertLHE(Str_)
  elseif(yRnd(5).lt.(30d0/39d0))then
    id(6)=convertLHE(Bot_)
  elseif(yRnd(5).lt.(32d0/39d0))then
    id(6)=convertLHE(ElM_)
  elseif(yRnd(5).lt.(34d0/39d0))then
    id(6)=convertLHE(MuM_)
  elseif(yRnd(5).lt.(36d0/39d0))then
    id(6)=convertLHE(TaM_)
  elseif(yRnd(5).lt.(37d0/39d0))then
    id(6)=convertLHE(NuE_)
    helicity(6)=sign(1d0,-dble(id(6)))
    helicity(7)=-helicity(6)
  elseif(yRnd(5).lt.(38d0/39d0))then
    id(6)=convertLHE(NuM_)
    helicity(6)=sign(1d0,-dble(id(6)))
    helicity(7)=-helicity(6)
  else
    id(6)=convertLHE(NuT_)
    helicity(6)=sign(1d0,-dble(id(6)))
    helicity(7)=-helicity(6)
  endif
  id(7)=-id(6)
  PostFac=39d0
  
elseif(DecayMode1.eq.10)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/3d0)then
    id(7)=convertLHE(ElP_)
    id(6)=convertLHE(NuE_)
  elseif(yRnd(5).lt.2d0/3d0)then
    id(7)=convertLHE(MuP_)
    id(6)=convertLHE(NuM_)
  else
    id(7)=convertLHE(TaP_)
    id(6)=convertLHE(NuT_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=3d0
  
elseif(DecayMode1.eq.11)then
  id(3)=convertLHE(Wp_)
  id(4)=convertLHE(Wp_)
  if(yRnd(5).lt.1d0/21d0)then
    id(7)=convertLHE(ElP_)
    id(6)=convertLHE(NuE_)
  elseif(yRnd(5).lt.(2d0/21d0))then
    id(7)=convertLHE(MuP_)
    id(6)=convertLHE(NuM_)
  elseif(yRnd(5).lt.(3d0/21d0))then
    id(7)=convertLHE(TaP_)
    id(6)=convertLHE(NuT_)
  elseif(yRnd(5).lt.(6d0/21d0))then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.(9d0/21d0))then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.(12d0/21d0))then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(AStr_)
  elseif(yRnd(5).lt.(15d0/21d0))then
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Adn_)
  elseif(yRnd(5).lt.(18d0/21d0))then
    id(6)=convertLHE(Up_)
    id(7)=convertLHE(Abot_)
  else
    id(6)=convertLHE(Chm_)
    id(7)=convertLHE(Abot_)
  endif
  helicity(6)=sign(1d0,-dble(id(6)))
  helicity(7)=-helicity(6)
  PostFac=21d0

else
  print *, "invalid final states"
  stop

endif

if(HbbDecays) PostFac=PostFac*6 !of which 2 is for summing H > bb~ helicities and 3 for summing colors


if( IsAZDecay(DecayMode1) .or. IsAPhoton(DecayMode1) ) then
!if pp collider
  if(Collider.eq.1.or.Collider.eq.2)then
    
    if( IsAZDecay(DecayMode1) )then
      call PDFMapping(14,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
    elseif( IsAPhoton(DecayMode1) )then
      call PDFMapping(16,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
    endif      
    if( Ehat .ge. 2d0*m_top ) return !without complex top mass in ZH amplitude
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
    
    call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1),ZAinterference=includeGammaStar)
    call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))

    if( applyPSCut .or. PSWgt.eq.zero ) return
    
    if( IsAZDecay(DecayMode1) )then
      call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
    elseif( IsAPhoton(DecayMode1) )then
      call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
    endif

    call setPDFs(eta1,eta2,pdf)

    FluxFac = 1d0/(2d0*EHat**2)
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt !/3d0! *6d0 !2 for e and mu, 3 for colors of b

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
  
    call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1),ZAinterference=includeGammaStar)
    call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))

    if( applyPSCut .or. PSWgt.eq.zero ) return

    FluxFac = 1d0/(2d0*ILC_Energy**2)
    PreFac = fbGeV2 * FluxFac * PSWgt
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

  call EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HbbDecays)
  call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays)
  if( applyPSCut .or. PSWgt.eq.zero ) return
  call SetRunningScales( (/ MomExt(1:4,5),MomExt(1:4,6),MomExt(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
  call setPDFs(eta1,eta2,pdf)

  FluxFac = 1d0/(2d0*EHat**2)
  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt !/3d0! *6d0 !2 for e and mu, 3 for colors of qqb

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
  if(Collider.eq.1.or.Collider.eq.2)then
    id(1:2) = (/ifound,jfound/)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2) * PreFac *PostFac
    EvalUnWeighted_VHiggs = LO_Res_Unpol
!if e+ e- collider
  else if(Collider.eq.0)then
    ifound=0
    jfound=0
    id(2)=convertLHE(ElM_)
    id(1)=-id(2)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    LO_Res_Unpol = me2 * PreFac *PostFac
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
  LO_Res_Unpol = me2 *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2) * PreFac *PostFac
  EvalUnWeighted_VHiggs = LO_Res_Unpol
endif


  CS_max = csmax(ifound,jfound)
  if( EvalUnWeighted_VHiggs.gt. CS_max) then
    write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_VHiggs, CS_max
    write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_VHiggs, CS_max
    AlertCounter = AlertCounter + 1
    Res = 0d0
  elseif( EvalUnWeighted_VHiggs .gt. yRnd(17)*CS_max ) then
    call Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
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
  if(Collider.eq.1.or.Collider.eq.2)then
    if((VHiggs_PC.eq."qq".or.VHiggs_PC.eq."lo").and.PChannel.eq.1)then
      do i = -5,5
        j = -i
        id(1:2) = (/i,j/)
        if (abs(i).ne.0)then
          call EvalAmp_VHiggs(id,helicity,MomExt,me2)
        else
          me2=0d0
        endif
        LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac *PostFac
        EvalUnWeighted_VHiggs = EvalUnWeighted_VHiggs+LO_Res_Unpol
        RES(i,j) = LO_Res_Unpol
        if (LO_Res_Unpol.gt.csmax(i,j)) then
          csmax(i,j) = LO_Res_Unpol
        endif
      enddo
    elseif((VHiggs_PC.eq."gg".or.VHiggs_PC.eq."bo".or.VHiggs_PC.eq."tr").and.PChannel.eq.0)then
      me2=0d0
      if( Ehat .ge. 2d0*m_top ) return !without complex top mass in ZH amplitude
      id(1)=0
      id(2)=0
      i=0
      j=0
      call EvalAmp_VHiggs(id,helicity,MomExt,me2)
      LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac *PostFac
      EvalUnWeighted_VHiggs = EvalUnWeighted_VHiggs+LO_Res_Unpol
      RES(i,j) = LO_Res_Unpol
      if (LO_Res_Unpol.gt.csmax(i,j)) then
        csmax(i,j) = LO_Res_Unpol
      endif
    endif
!if e+ e- collider
  else if(Collider.eq.0)then
    id(2)=convertLHE(ElM_)
    id(1)=-id(2)
    call EvalAmp_VHiggs(id,helicity,MomExt,me2)
    LO_Res_Unpol = me2 * PreFac *PostFac
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
     else
           me2=0d0
     endif
    LO_Res_Unpol = me2 *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac *PostFac
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




END MODULE ModCrossSection_VHiggs







