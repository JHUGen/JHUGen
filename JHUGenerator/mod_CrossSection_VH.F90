MODULE ModCrossSection_VH
  use ModVHaux
  implicit none
  integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
  integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)

 CONTAINS





Function EvalWeighted_VH(yRnd,VgsWgt)
  use ModKinematics
  use ModKinematics_VH
  use ModParameters
  use ModMisc
  use ModVHLO
  use ModVHreal
  use ModVHvirtual
  use ModVHdipole
  use ModVHgg
  use ModVHqg
  use ModPhasespace
#if compiler==1
 use ifport
#endif
  implicit none
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5): flavor in Z/W decay mode
! yrnd(6:7): phi_4 and cos(theta_4) in the CM frame of Z*(3)
! yrnd(8:9): phi_6 and cos(theta_6) in the CM frame of decay product of Z(4)
! yrnd(10:11): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
! yRnd(12): inv_mass(4)
! yRnd(13): inv_mass(5)
! yRnd(14:16): real emission phase space
! yRnd(17): NLO integration for + distribution for PDF renormalization
! yrnd(18:19): PDF mapping
! yrnd(20): flavor of j in gq > WH+j
! yrnd(21): partonic channel in unweighted events
! yRnd(22): accept or reject in unweighted events
  real(8) :: yRnd(1:19),VgsWgt, EvalWeighted_VH
  real(8) :: pdf(-6:6,1:2), pdf_ren(-6:6,1:2,1:2), pdf_dip(-6:6,1:2,1:2), pdf_real(-6:6,1:2)
  real(8) :: eta1, eta2, eta1_ren(1:2), eta2_ren(1:2), FluxFac, FluxFac_ren, Ehat, Ehat_ren, shat, shat_ren, sHatJacobi, sHatJacobi_ren
  real(8) :: Mom(1:4,1:10), Mom_real(1:4,1:10), PSWgt, PSWgt_ren, PSWgt_tilde(1:2), PSWgt2
  complex(8) :: amp_dummy,amp_dummy_dummy,amp_virture_finite
  real(8) :: me2real,me2lo,dip(1:2),me2gg,me2gq,me2sub,me2sup,me2_dummy,me2_dummy_dummy,me2_virture_finite,me2_PDF_ren(1:2)!, lheweight(-6:6,-6:6)
  integer :: i_dipole, i_PDF_ren
  real(8) :: ptilde(1:4,1:9,1:2),p_PDF_ren(1:4,1:9,1:2)
  integer :: i,j,k,l,p,q,NBin(1:NumHistograms),NHisto
  real(8) :: PreFac, PreFac_ren, PostFac, Mu_Fact_ren(1:2)
  real(8) :: alphas_ren(1:2), alphas_dip(1:2)
  logical :: applyPSCut, applyPSCut_ren(1:2)
  real(8) :: mass(1:10,1:2)
  real(8) :: helicity(10)
  integer id(10), id2(10)
  real(8) :: finite_factor
  real(8) :: Mu_Fact_real, Mu_ren_real, alphas_real, PSWgt_real, PreFac_real
  real(8) :: Mom_save(1:4,1:9)
! for tests!!!!!!!!!!!!!!
real(8) :: MomExt1(1:4,1:10),MomExt2(1:4,1:10),MomExt3(1:4,1:10),MomExt4(1:4,1:10),MomExt1t(1:4,1:9),MomExt2t(1:4,1:9),MomExt3t(1:4,1:9)
! for tests!!!!!!!!!!!!!!

  EvalWeighted_VH=0d0
  EvalCounter = EvalCounter+1

! initialization and decaying mode related
  id(:)=0
  helicity(:)=0
  mass(1:2,1:2)=0d0
  Mom=0d0
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
  mass(6:10,1:2)=0d0

  id(5)=convertLHE(Hig_)
  if(HbbDecays.eqv..true.)then
    id(8)=convertLHE(Bot_)
    id(9)=-id(8)
  else
    id(8)=Not_a_particle_
    id(9)=Not_a_particle_
  endif

  id(10)=convertLHE(Glu_)!This may need to be overwritten later for gq/qg!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    PostFac=2d0  ! of which 2 is for summing photon polarizations
      
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
    stop 1
    
  endif

  if(HbbDecays) PostFac=PostFac*6d0 !of which 2 is for summing H > bb~ helicities and 3 for summing bb~ colors

! end initialization and decaying mode related

    if(VH_PC.ne."ee".and. &! "ee" ( = e+ e- @LO)
       VH_PC.ne."qq".and. &! "qq" ( = q qbar @LO)
       VH_PC.ne."lo".and. &! "lo" ( = q qbar @LO)
       VH_PC.ne."tr".and. &! "tr" ( = triangles of gg)
       VH_PC.ne."bo".and. &! "bo" ( = boxes of gg)
       VH_PC.ne."gg".and. &! "gg" ( = triangles + boxes of gg)
       VH_PC.ne."sp".and. &! "sp" ( = virtual + dipoles, for development only)
       VH_PC.ne."sb".and. &! "sb" ( = real - dipoles, for development only)
       VH_PC.ne."gq".and. &! "qg" or "gq" ( = qg + gq)
       VH_PC.ne."qg".and. &! "qg" or "gq" ( = qg + gq)
       VH_PC.ne."nl")then  !"nl" ( = full oneloop = q qbar @LO + NLO + gg + gq)
      print*,"invalid VH_PC ", VH_PC
      stop 1
    endif

! begin event

  if(((Collider.eq.1.or.Collider.eq.2).and.VH_PC.eq."ee") .or. &
      (Collider.eq.0                  .and.VH_PC.ne."ee"))then
    print*,"e+ e- collisions with Collider=0 only."
    print*,"VH_PC =",VH_PC," Collider =", Collider
    stop 1
  endif

!if e+ e- collider
  if(Collider.eq.0.and.VH_PC.eq."ee")then  
    call EvalPhasespace_VH(yRnd(6:13),ILC_Energy,Mom(:,1:9),id(6:9),PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
    Mom_save(1:4,1:9)=Mom(1:4,1:9)
    call Kinematics_VH(id,Mom,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
    if( applyPSCut .or. PSWgt.eq.zero ) return    
    FluxFac = 1d0/(2d0*ILC_Energy**2)
    PreFac = fbGeV2 * FluxFac * PSWgt
    EvalWeighted_VH=0d0
    id(1:2)=(/convertLHE(ElP_),convertLHE(ElM_)/)
    call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_dummy)
    EvalWeighted_VH = dble(amp_dummy*dconjg(amp_dummy)) *PreFac *PostFac

!if pp/ppbar collider
  elseif(Collider.eq.1.or.Collider.eq.2)then
    !PDFMapping
    call PDFMapping(1,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!Z/W
!print*,Ehat
    !needed for PDFMapping=1
    if(Ehat.le.(M_Z+M_Reso))return!!!!may be improved
!marginally more efficient
!    if( IsAZDecay(DecayMode1) )then
!      call PDFMapping(14,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!Z      
!    elseif( IsAWDecay(DecayMode1) )then
!      call PDFMapping(15,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!W
!    elseif( IsAPhoton(DecayMode1) )then
!      call PDFMapping(16,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!Z/gamma with gamma in the final state
!    endif
!marginally more efficient

    !Phase space, scales, and PDF's
    if(VH_PC.eq."qq".or.VH_PC.eq."lo".or.VH_PC.eq."tr".or.VH_PC.eq."bo".or.VH_PC.eq."gg".or.VH_PC.eq."sp".or.VH_PC.eq."gq".or.VH_PC.eq."qg".or.VH_PC.eq."nl")then
      call EvalPhasespace_VH(yRnd(6:13),Ehat,Mom(:,1:9),id(6:9),PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
      Mom(:,10)=0d0 ! no QCD particle emitted
      Mom_save(:,1:9)=Mom(:,1:9)
      !boost from center of mass frame to lab frame
      call boost2Lab(eta1,eta2,10,Mom)
!print*,"=========="
!print*,eta1,eta2
!print*,eta1*yRnd(17),eta2*yRnd(17)
      !Kinematics and cuts
      call Kinematics_VH(id,Mom,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
      if( applyPSCut .or. PSWgt.eq.zero )return
      !Set running scales
      if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) )then
        call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      elseif( IsAPhoton(DecayMode1) )then
        call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
      endif
      ! do not forget to set scales in command lines
      call EvalAlphaS()
      call setPDFs(eta1,eta2,pdf)
      if(VH_PC.eq."sp".or.VH_PC.eq."gq".or.VH_PC.eq."qg".or.VH_PC.eq."nl")then
        Ehat_ren = Ehat * dsqrt(yRnd(17))
        if(Ehat_ren.le.(M_Z+M_Reso))return!!!!may be improved
        shat = Ehat**2
        shat_ren = shat * yRnd(17)
        call EvalPhasespace_VH(yRnd(6:13),Ehat_ren,p_PDF_ren(:,:,1),id(6:9),PSWgt_ren,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
        p_PDF_ren(:,:,2)=p_PDF_ren(:,:,1)
        eta1_ren(1) = yRnd(17)*eta1
        eta2_ren(1) = eta2
        eta1_ren(2) = eta1
        eta2_ren(2) = yRnd(17)*eta2
        do i_PDF_ren=1,2
          call boost2Lab(eta1_ren(i_PDF_ren),eta2_ren(i_PDF_ren),9,p_PDF_ren(:,:,i_PDF_ren))
          call Kinematics_VH(id,p_PDF_ren(:,:,i_PDF_ren),NBin,applyPSCut_ren(i_PDF_ren),HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
          !Set running scales
          if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) )then
            call SetRunningScales( (/ p_PDF_ren(1:4,5,i_PDF_ren),p_PDF_ren(1:4,6,i_PDF_ren),p_PDF_ren(1:4,7,i_PDF_ren) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
          elseif( IsAPhoton(DecayMode1) )then
            call SetRunningScales( (/ p_PDF_ren(1:4,5,i_PDF_ren),p_PDF_ren(1:4,6,i_PDF_ren),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
          endif
          Mu_Fact_ren(i_PDF_ren)=Mu_Fact
          ! do not forget to set scales in command lines
          call EvalAlphaS()
          alphas_ren(i_PDF_ren) = alphas
          call setPDFs(eta1_ren(i_PDF_ren),eta2_ren(i_PDF_ren),pdf_ren(:,:,i_PDF_ren))
        enddo!i_PDF_ren
        FluxFac_ren = 1d0/(2d0*Ehat_ren**2)
        sHatJacobi_ren = sHatJacobi * 1d0!!!!may change if PDFMapping gets improved.
        PreFac_ren = fbGeV2 * FluxFac_ren * sHatJacobi_ren * PSWgt_ren            
      endif
    endif

    if(VH_PC.eq."sb".or.VH_PC.eq."nl")then
!print*,yRnd(6:16)
!print*,"=================="
      call EvalPhasespace_VHglu(yRnd(6:16),Ehat,Mom_real,id(6:9),PSWgt_real,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
!if(isNan(ehat))print*,eta1,eta2,Ehat,yrnd(18:19)
!print*,Mom_real
!call EvalPhasespace_VHglu_singular(yRnd(6:16),Ehat,Mom_real,id(6:9),PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
!print*,"IR",(Mom(:,1).dot.Mom(:,10))/Ehat**2,(Mom(:,2).dot.Mom(:,10))/Ehat**2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( ((Mom_real(:,1).dot.Mom_real(:,10))/Ehat**2 .lt. IRmin) .or. ((Mom_real(:,2).dot.Mom_real(:,10))/Ehat**2 .lt. IRmin) )return !for numerical stibility, should not be physical

!if(PSWgt.eq.0d0.or.Mom_real(1,5).lt.0d0)then
!print*,Mom_real(:,1)
!print*,Mom_real(:,2)
!print*,Mom_real(:,3)
!print*,Mom_real(:,4)
!print*,Mom_real(:,5)
!print*,Mom_real(:,6)
!print*,Mom_real(:,7)
!print*,Mom_real(:,8)
!print*,Mom_real(:,9)
!print*,Mom_real(:,10)
!call boost2Lab(eta1,eta2,10,Mom_real)
!print*,"----------------"
!print*,Mom_real(:,1).dot.Mom_real(:,1),Mom_real(1,1),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,2).dot.Mom_real(:,2),Mom_real(1,2),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,3).dot.Mom_real(:,3),Mom_real(1,3),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,4).dot.Mom_real(:,4),Mom_real(1,4),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,5).dot.Mom_real(:,5),Mom_real(1,5),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,6).dot.Mom_real(:,6),Mom_real(1,6),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,7).dot.Mom_real(:,7),Mom_real(1,7),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,8).dot.Mom_real(:,8),Mom_real(1,8),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,9).dot.Mom_real(:,9),Mom_real(1,9),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,10).dot.Mom_real(:,10),Mom_real(1,10),(2d0*getMass(convertLHEreverse(id(8))))**2
  !print*,"==================="
  !pause!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!endif
      !boost from center of mass frame to lab frame
      call boost2Lab(eta1,eta2,10,Mom_real)
      !Kinematics and cuts
      call Kinematics_VH(id,Mom_real,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
      if( applyPSCut .or. PSWgt_real.eq.zero )return
!if(Get_PT(Mom_real(:,10)).le.20d0*GeV)return !jet pt cut such that real contribution .is finite
      !Set running scales
      if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) )then
        call SetRunningScales( (/ Mom_real(1:4,5),Mom_real(1:4,6),Mom_real(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      elseif( IsAPhoton(DecayMode1) )then
        call SetRunningScales( (/ Mom_real(1:4,5),Mom_real(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
      endif
      Mu_Fact_real = Mu_Fact
      Mu_ren_real = mu_Ren
      ! do not forget to set scales in command lines
      call EvalAlphaS()
      alphas_real = alphas
      call setPDFs(eta1,eta2,pdf_real)
      call p_tilde(Mom_real,ptilde)
      ptilde(:,3,:) = ptilde(:,1,:)+ptilde(:,2,:)
      ptilde(:,4,:) = ptilde(:,6,:)+ptilde(:,7,:)
      do i_dipole=1,2
        if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) )then
          call SetRunningScales( (/ ptilde(1:4,5,i_dipole),ptilde(1:4,6,i_dipole),ptilde(1:4,7,i_dipole) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
        elseif( IsAPhoton(DecayMode1) )then
          call SetRunningScales( (/ ptilde(1:4,5,i_dipole),ptilde(1:4,6,i_dipole),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
        endif
        call EvalAlphaS()
        alphas_dip(i_dipole) = alphas
        call setPDFs(eta1,eta2,pdf_dip(-6:6,1:2,i_dipole))
      enddo
      FluxFac = 1d0/(2d0*Ehat**2)
      PreFac_real = fbGeV2 * FluxFac * sHatJacobi * PSWgt_real
    endif

    !Set running scales BACK. Everything needed has been stored... (I hope so!)
    if(VH_PC.ne."sb")then
      if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) )then
        call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
      elseif( IsAPhoton(DecayMode1) )then
        call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
      endif
      call EvalAlphaS()
      ! do not forget to set scales in command lines
      FluxFac = 1d0/(2d0*EHat**2)
      PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt
    endif

    EvalWeighted_VH=0d0

!Momext1(:,1)=(/6.4524836219780077d0,        0.0000000000000000d0,    0.0000000000000000d0,        6.4524836219780077d0 /)
!Momext1(:,2)=(/4.8184722322140914d0,        0.0000000000000000d0,    0.0000000000000000d0,       -4.8184722322140914d0 /)
!Momext1(:,6)=(/0.64327343583692165d0,      -0.39848130152182004d0,    -2.0663649089824071d-002, -0.50456553508657986d0 /)
!Momext1(:,7)=(/4.5071972333345549d0,       -3.8828400769136620d0,    -1.4140898576572971d0,       -1.7996471075637406d0 /)
!Momext1(:,5)=(/6.1204851850206214d0,        4.2813213784354822d0,    1.4347535067471213d0,        3.9382240324142379d0 /)
!MomExt1(:,3) = MomExt1(:,1)+MomExt1(:,2)
!MomExt1(:,4) = MomExt1(:,6)+MomExt1(:,7)
!MomExt1(:,8) = 0d0
!MomExt1(:,9) = 0d0
!MomExt1(:,10) = 0d0
!
!Momext2(:,1)=(/0.97480053689893575d0,        0.0000000000000000d0, 0.0000000000000000d0,       0.97480053689893575d0 /)
!Momext2(:,2)=(/1.4350050018754614d0,        0.0000000000000000d0, 0.0000000000000000d0,       -1.4350050018754614d0 /)
!Momext2(:,6)=(/0.68336592278878705d0,       -5.6516824645623021d-002, -0.65522446557556058d0,      -0.18567641926794989d0 /)
!Momext2(:,7)=(/0.37430166987959523d0,       0.27974435028580369d0, 0.23955804829684374d0,       -6.6759119631958977d-002 /)
!Momext2(:,5)=(/1.3521379461060148d0,      -0.22322752564018064d0, 0.41566641727871684d0,      -0.20776892607661673d0 /)
!MomExt2(:,3) = MomExt2(:,1)+MomExt2(:,2)
!MomExt2(:,4) = MomExt2(:,6)+MomExt2(:,7)
!MomExt2(:,8) = 0d0
!MomExt2(:,9) = 0d0
!MomExt2(:,10) = 0d0
!
!
!
!Momext3(:,1) = (/ 6.4524836219780077d0,        0.0000000000000000d0, 0.0000000000000000d0,        6.4524836219780077d0 /)
!Momext3(:,2) = (/ 4.8184722322140914d0,        0.0000000000000000d0, 0.0000000000000000d0,       -4.8184722322140914d0 /)
!Momext3(:,6) = (/ 0.12578560564158095d0,        1.3118396634171864d-002, -0.10162933196370739d0,        7.2947961869941402d-002 /)
!Momext3(:,7) = (/ 0.88779616338491874d0,       -6.6710343680087852d-002, -0.42302981787951177d0,       0.77767443763562971d0 /)
!Momext3(:,5) = (/ 5.3808667464113640d0,       0.67712353606180309d0, -2.5375058271648707d0,        4.5270072715482916d0 /)
!Momext3(:,10) = (/ 4.8765073387542346d0,      -0.62353158901588701d0, 3.0621649770080897d0,       -3.7436182812899457d0 /)
!MomExt3(:,3) = MomExt3(:,1)+MomExt3(:,2)-MomExt3(:,10)
!MomExt3(:,4) = MomExt3(:,6)+MomExt3(:,7)
!MomExt3(:,8) = 0d0
!MomExt3(:,9) = 0d0
!
!
!Momext4(:,1) = (/  0.97480053689893575d0,        0.0000000000000000d0, 0.0000000000000000d0,       0.97480053689893575d0 /)
!Momext4(:,2) = (/   1.4350050018754614d0,        0.0000000000000000d0, 0.0000000000000000d0,       -1.4350050018754614d0 /)
!Momext4(:,6) = (/  0.20774249642309217d0,      -0.14249354751410651d0, 8.7199473183037501d-005, -0.15117052005333384d0 /)
!Momext4(:,7) = (/  0.13192679274531913d0,       -9.7621490917663784d-002, -6.2189942016411301d-002,  -6.3301929410360641d-002 /)
!Momext4(:,5) = (/   1.3761528216957035d0,      -0.40227650161239736d0, -0.16529074650324624d0,      -0.37702675506114647d0 /)
!Momext4(:,10) = (/  0.69398342791028245d0,       0.64239154004416754d0, 0.22739348904647450d0,       0.13129473954831528d0 /)
!MomExt4(:,3) = MomExt4(:,1)+MomExt4(:,2)-MomExt4(:,10)
!MomExt4(:,4) = MomExt4(:,6)+MomExt4(:,7)
!MomExt4(:,8) = 0d0
!MomExt4(:,9) = 0d0


!MomExt1(:,1)  = (/    3.6993467144712677d0,        0.0000000000000000d0,    0.0000000000000000d0,        3.6993467144712677d0 /)
!MomExt1(:,2)  = (/    3.6993467144712677d0,        0.0000000000000000d0,    0.0000000000000000d0,       -3.6993467144712677d0 /)
!MomExt1(:,5)  = (/    2.3219815450459187d0,      -0.44511328965007180d0,    0.46855893187333120d0,        1.8470043265440523d0 /)
!MomExt1(:,6)  = (/    1.1469531245593170d0,      -0.17902363116310704d0,    0.17891445071840281d0,        1.1186785189435293d0 /)
!MomExt1(:,7)  = (/    0.58938520841107445d0,       0.22417316451991204d0,    -0.42356728395178866d0,       0.34309192961741136d0 /)
!MomExt1(:,10) = (/    3.3403735509262251d0,       0.39996375629326680d0,    -0.22390609863994532d0,       -3.3087747751049927d0 /)
!MomExt1(:,3) = MomExt1(:,1)+MomExt1(:,2)-MomExt1(:,10)
!MomExt1(:,4) = MomExt1(:,6)+MomExt1(:,7)
!
!
!
!MomExt2(:,1)  = (/    5.6071343647042067d0,        0.0000000000000000d0, 0.0000000000000000d0,        5.6071343647042067d0 /)
!MomExt2(:,2)  = (/    5.6071343647042067d0,        0.0000000000000000d0, 0.0000000000000000d0,       -5.6071343647042067d0 /)
!MomExt2(:,5)  = (/    4.3435236991819917d0,        1.3651075645974082d0, -2.8073784259512808d0,        2.7493282517622109d0 /)
!MomExt2(:,6)  = (/    1.2972599389283468d0,       0.51286299876669883d0, -0.76013317363117228d0,       0.91763416021299538d0 /)
!MomExt2(:,7)  = (/    9.7568396212771191d-002,  -5.6646395025800364d-002, -7.8799491502101748d-002,  -1.0070650878840348d-002 /)
!MomExt2(:,10) = (/    5.4759166950853038d0,       -1.8213241683383068d0, 3.6463110910845549d0,       -3.6568917610963658d0 /)
!MomExt2(:,3) = MomExt2(:,1)+MomExt2(:,2)-MomExt2(:,10)
!MomExt2(:,4) = MomExt2(:,6)+MomExt2(:,7)
!
!
!
!
!vev=2.4621845810181631d0
!
!MomExt1t(:,1) = (/  1.4612636158434757d0  ,        0.0000000000000000d0,        0.0000000000000000d0     ,   1.4612636158434757d0     /)
!MomExt1t(:,2) = (/  1.4612636158434757d0  ,        0.0000000000000000d0,        0.0000000000000000d0     ,  -1.4612636158434757d0     /)
!MomExt1t(:,3) = (/  2.9225272316869515d0  ,        0.0000000000000000d0,        0.0000000000000000d0     ,   0.0000000000000000d0     /)
!MomExt1t(:,4) = (/  1.3351579764951791d0  ,       0.23912248810916401d0,      -0.49887803363686017d0     , -0.80695867029432466d0     /)
!MomExt1t(:,5) = (/  1.5873692551917724d0  ,      -0.23912248810916401d0,       0.49887803363686017d0     ,  0.80695867029432466d0     /)
!MomExt1t(:,6) = (/ 0.32428361334240258d0  ,      -0.30847990233180261d0,       -4.0824295707895482d-002  ,  -9.1287395733039817d-002  /)
!MomExt1t(:,7) = (/  1.0108743631527766d0  ,       0.54760239044096659d0,      -0.45805373792896470d0     , -0.71567127456128488d0     /)
!MomExt1t(:,8) = (/  0.0000000000000000d0  ,        0.0000000000000000d0,        0.0000000000000000d0     ,   0.0000000000000000d0     /)
!MomExt1t(:,9) = (/  0.0000000000000000d0  ,        0.0000000000000000d0,        0.0000000000000000d0     ,   0.0000000000000000d0     /)
!! ======================
!MomExt2t(:,1) = (/   1.4993889013229147d0,        0.0000000000000000d0     ,   0.0000000000000000d0,        1.4993889013229147d0     /)
!MomExt2t(:,2) = (/   1.4993889013229147d0,        0.0000000000000000d0     ,   0.0000000000000000d0,       -1.4993889013229147d0     /)
!MomExt2t(:,3) = (/   2.9987778026458294d0,        0.0000000000000000d0     ,   0.0000000000000000d0,        0.0000000000000000d0     /)
!MomExt2t(:,4) = (/   1.3795905093925491d0,       0.37240043559138541d0     , -0.90237656049098836d0,      -0.32603979708109893d0     /)
!MomExt2t(:,5) = (/   1.6191872932532803d0,      -0.37240043559138541d0     ,  0.90237656049098836d0,       0.32603979708109893d0     /)
!MomExt2t(:,6) = (/  0.69987628038380456d0,       0.41499286401781232d0     , -0.51432374828037830d0,       0.23038839513522705d0     /)
!MomExt2t(:,7) = (/  0.67971422900874456d0,       -4.2592428426426909d-002  , -0.38805281221061005d0,      -0.55642819221632600d0     /)
!MomExt2t(:,8) = (/   0.0000000000000000d0,        0.0000000000000000d0     ,   0.0000000000000000d0,        0.0000000000000000d0     /)
!MomExt2t(:,9) = (/   0.0000000000000000d0,        0.0000000000000000d0     ,   0.0000000000000000d0,        0.0000000000000000d0     /)
!! ======================
!MomExt3t(:,1) = (/   4.9249999999999998d0,        0.0000000000000000d0,        0.0000000000000000d0,        4.9249999999999998d0 /)
!MomExt3t(:,2) = (/   4.9249999999999998d0,        0.0000000000000000d0,        0.0000000000000000d0,       -4.9249999999999998d0 /)
!MomExt3t(:,3) = (/   9.8499999999999996d0,        0.0000000000000000d0,        0.0000000000000000d0,        0.0000000000000000d0 /)
!MomExt3t(:,4) = (/   4.8874276696732233d0,       -1.6191114510912707d0,       -4.2507436371871847d0,       -1.5408701352102065d0 /)
!MomExt3t(:,5) = (/   4.9625723303267764d0,        1.6191114510912707d0,        4.2507436371871847d0,        1.5408701352102065d0 /)
!MomExt3t(:,6) = (/   3.7143691570810575d0,       -1.4939929156902811d0,       -3.0908500572228248d0,       -1.4181570176493015d0 /)
!MomExt3t(:,7) = (/   1.1730585125921660d0,      -0.12511853540098977d0,       -1.1598935799643604d0,      -0.12271311756090503d0 /)
!MomExt3t(:,8) = (/   0.0000000000000000d0,        0.0000000000000000d0,        0.0000000000000000d0,        0.0000000000000000d0 /)
!MomExt3t(:,9) = (/   0.0000000000000000d0,        0.0000000000000000d0,        0.0000000000000000d0,        0.0000000000000000d0 /)


!!id(1:2)=(/convertLHE(up_),convertLHE(Aup_)/)
!id(1:2)=(/convertLHE(dn_),convertLHE(Adn_)/)
!
!print*,"========================="
!Mom=MomExt1
!me2lo=0d0
!call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!do l=0,1
!do p=0,1
!!print*,"helicities",(l*2-1),(p*2-1)
!call amp_VH_LO(Mom(:,1:9),mass(3:5,:),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9)/),id(1:9),amp_dummy)
!me2lo = me2lo + dble(amp_dummy*dconjg(amp_dummy))
!!print*,"MomExt1 amptd = ",amp_dummy
!enddo
!enddo
!me2lo = me2lo * aveqq * 3d0
!print*,"MomExt1 me2qq = ",me2lo
!
!Mom=MomExt2
!me2lo=0d0
!call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!do l=0,1
!do p=0,1
!!print*,"helicities",(l*2-1),(p*2-1)
!call amp_VH_LO(Mom(:,1:9),mass(3:5,:),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9)/),id(1:9),amp_dummy)
!me2lo = me2lo + dble(amp_dummy*dconjg(amp_dummy))
!!print*,"MomExt1 amptd = ",amp_dummy
!enddo
!enddo
!me2lo = me2lo * aveqq * 3d0
!print*,"MomExt2 me2qq = ",me2lo
!
!Mom=MomExt3
!me2real=0d0
!call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!alphas_real=alphas
!do l=0,1
!do p=0,1
!do q=0,1
!  call amp_VH_real(Mom,mass(3:5,1:2),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9),dble(q*2-1)/),id,amp_dummy)
!  me2real = me2real + dble(amp_dummy*dconjg(amp_dummy))
!enddo
!enddo
!enddo
!me2real = me2real *aveqq *Cf *3d0 !I think I understand now.
!me2real = me2real * (4d0 * pi * alphas_real)! gs^2
!print*,"MomExt3 me2real = ",me2real
!
!Mom=MomExt4
!me2real=0d0
!call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!alphas_real=alphas
!do l=0,1
!do p=0,1
!do q=0,1
!  call amp_VH_real(Mom,mass(3:5,1:2),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9),dble(q*2-1)/),id,amp_dummy)
!  me2real = me2real + dble(amp_dummy*dconjg(amp_dummy))
!enddo
!enddo
!enddo
!me2real = me2real *aveqq *Cf *3d0 !I think I understand now.
!me2real = me2real * (4d0 * pi * alphas_real)! gs^2
!print*,"MomExt4 me2real = ",me2real
!
!pause
!print*,"========================="
!



!print*,"========================="
!do l=0,1
!do p=0,1
!do q=0,1
!print*,"helicities",(l*2-1),(p*2-1),(q*2-1)
!call SetRunningScales( (/ MomExt1t(1:4,5),MomExt1t(1:4,6),MomExt1t(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!print *, "alphas = ",alphas
!call amp_VH_gg(MomExt1t(:,1:9),mass(3:5,:),(/dble(l*2-1),dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id(1:9),amp_dummy)
!me2gg = dble(amp_dummy*dconjg(amp_dummy)) *GluonColAvg**2 *8d0 !8 = summing delta(a,b)*delta(a,b)
!print*,"MomExt1t amptd = ",amp_dummy
!print*,"MomExt1t me2gg = ",me2gg
!enddo
!enddo
!enddo
!do l=0,1
!do p=0,1
!do q=0,1
!print*,"helicities",(l*2-1),(p*2-1),(q*2-1)
!call SetRunningScales( (/ MomExt2t(1:4,5),MomExt2t(1:4,6),MomExt2t(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!print *, "alphas = ",alphas
!call amp_VH_gg(MomExt2t(:,1:9),mass(3:5,:),(/dble(l*2-1),dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id(1:9),amp_dummy)
!me2gg = dble(amp_dummy*dconjg(amp_dummy)) *GluonColAvg**2 *8d0 !8 = summing delta(a,b)*delta(a,b)
!print*,"MomExt2t amptd = ",amp_dummy
!print*,"MomExt2t me2gg = ",me2gg
!enddo
!enddo
!enddo
!do l=0,1
!do p=0,1
!do q=0,1
!print*,"helicities",(l*2-1),(p*2-1),(q*2-1)
!call SetRunningScales( (/ MomExt3t(1:4,5),MomExt3t(1:4,6),MomExt3t(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()
!print *, "alphas = ",alphas
!call amp_VH_gg(MomExt3t(:,1:9),mass(3:5,:),(/dble(l*2-1),dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id(1:9),amp_dummy)
!me2gg = dble(amp_dummy*dconjg(amp_dummy)) *GluonColAvg**2 *8d0 !8 = summing delta(a,b)*delta(a,b)
!print*,"MomExt3t amptd = ",amp_dummy
!print*,"MomExt3t me2gg = ",me2gg
!enddo
!enddo
!enddo
!pause
!print*,"========================="

!Mom(:,1:9)=MomExt1t(:,1:9)
!call SetRunningScales( (/ MomExt1t(1:4,5),MomExt1t(1:4,6),MomExt1t(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
!call EvalAlphaS()

    !gg
    me2gg=0d0
    if(VH_PC.eq."gg".or.VH_PC.eq."bo".or.VH_PC.eq."tr")then
      if(DecayMode1.eq.4.or.DecayMode1.eq.5.or.DecayMode1.eq.6.or.DecayMode1.eq.10.or.DecayMode1.eq.11)then
        print*,"DecayMode1 = ",DecayMode1," which is a W decay, not compatible with gg"
        stop 1
      endif
      id(1:2)=(/convertLHE(Glu_),convertLHE(Glu_)/)
      call amp_VH_gg(Mom(:,1:9),mass(3:5,1:2),helicity,id(1:9),amp_dummy)
!===================This section average helicities, which results the same==========================
!me2_dummy=0d0
!do l=0,1
!do p=0,1
!do q=0,1
!  call amp_VH_gg(Mom(:,1:9),mass(3:5,1:2),(/dble(l*2-1),dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id(1:9),amp_dummy_dummy)
!  me2_dummy=me2_dummy+dble(amp_dummy_dummy*dconjg(amp_dummy_dummy))
!enddo
!enddo
!enddo
!me2_dummy=me2_dummy/8d0
!!print*,"me2_dummy = ",me2_dummy*pdf(0,1)*pdf(0,2)*PreFac
!me2gg = me2_dummy *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0
!===================This section average helicities, which results the same==========================
      me2gg = dble(amp_dummy*dconjg(amp_dummy)) *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0
      !8 = summing delta(a,b)*delta(a,b)
    endif

    !lo/qq
    me2lo=0d0
    if(VH_PC.eq."qq".or.VH_PC.eq."lo".or.VH_PC.eq."nl")then
      !Z/gamma
      if( IsAZDecay(DecayMode1).or.IsAPhoton(DecayMode1) )then
        do i = -5,5
          j = -i
          id(1:2) = (/i,j/)
          if (i.ne.0)then
            call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_dummy)
!print*,dble(amp_dummy*dconjg(amp_dummy))
          else
            amp_dummy=0d0
          endif
          me2lo = me2lo + dble(amp_dummy*dconjg(amp_dummy)) *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2)
        enddo! parton
        me2lo = me2lo *PreFac *PostFac * QuarkColAvg**2 * 3d0
        !summing 3 colors in intial qq, no factor from spins because they are casted randomly, not summed.
!print*,me2lo
!print*,Mom(:,1)
!print*,Mom(:,2)
!print*,Mom(:,3)
!print*,Mom(:,4)
!print*,Mom(:,5)
!print*,Mom(:,6)
!print*,Mom(:,7)
!print*,Mom(:,8)
!print*,Mom(:,9)
!print*,"-------------"
!pause

      !W
      elseif( IsAWDecay(DecayMode1) )then
        do i = -5,5 !gluon = 0, otherwise PDG codes
        do j = -5,5

          id2=id
          id2(1:2) = (/i,j/)
          !W+
          if( CouplToLHEWp(id2(1:2)) )then
            helicity(6)=sign(1d0,-dble(id2(6)))
            helicity(7)=-helicity(6)
          !W-
          elseif( CouplToLHEWm(id2(1:2)) )then
            id2(3)=-id2(3)
            id2(4)=-id2(4)
            id2(6)=-id2(6)
            id2(7)=-id2(7)
            helicity(6)=sign(1d0,-dble(id2(6)))
            helicity(7)=-helicity(6)
          else
            cycle !skip this i,j combination
          endif
          call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id2(1:9),amp_dummy)
          me2lo = me2lo + dble(amp_dummy*dconjg(amp_dummy))*pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2)
         
        enddo! parton
        enddo! parton
        me2lo = me2lo *PreFac *PostFac *QuarkColAvg**2 *3d0
        !summing 3 colors in intial qq, no factor from spins because they are casted randomly, not summed.

      endif
    endif


!for checking cancellation approaching singularities, set FacScheme=-1 RenScheme=-1 for this test
!call EvalPhasespace_VHglu_singular(yRnd(6:16),Ehat,Mom_real,id(6:9),PSWgt_real,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
!alphas_real=alphas
!call p_tilde(Mom_real,ptilde)
!print*,"IR",(Mom_real(:,1).dot.Mom_real(:,10))/Ehat**2,(Mom_real(:,2).dot.Mom_real(:,10))/Ehat**2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !NLO: real - dipoles
    me2sub=0d0
    if(VH_PC.eq."nl".or.VH_PC.eq."sb")then
      !Z/gamma
      if( IsAZDecay(DecayMode1).or.IsAPhoton(DecayMode1) )then
        do i = -5,5
        !i=1
          j = -i

          id(1:2) = (/i,j/)
          me2real=0d0
          if (i.eq.0)cycle !skip gg
          do l=0,1
          do p=0,1
          do q=0,1
            call amp_VH_real(Mom_real,mass(3:5,1:2),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9),dble(q*2-1)/),id,amp_dummy)
            me2real = me2real + dble(amp_dummy*dconjg(amp_dummy))
          enddo
          enddo
          enddo
          me2real = me2real *pdf_real(LHA2M_PDF(i),1)*pdf_real(LHA2M_PDF(j),2) *PreFac_real *PostFac *aveqq *3d0 !I think I understand now.
          me2real = me2real * Cf * (4d0 * pi * alphas_real)! gs^2

          call Kinematics_VH(id,Mom_real,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
          do NHisto = 1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),me2real*VgsWgt)
          enddo

          dip=0d0
          do i_dipole=1,2
            me2_dummy=0d0
            do p=0,1
            do q=0,1
              call amp_VH_LO(ptilde(:,:,i_dipole),mass(3:5,1:2),(/dble(p*2-1),-dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id,amp_dummy)
              me2_dummy = me2_dummy + dble(amp_dummy*dconjg(amp_dummy))
            enddo
            enddo
            me2_dummy = me2_dummy *pdf_dip(LHA2M_PDF(i),1,i_dipole)*pdf_dip(LHA2M_PDF(j),2,i_dipole) *PreFac_real *PostFac * aveqq * 3d0
            dip(i_dipole) = me2_dummy * (-1d0) / (2d0*Mom_real(:,i_dipole).dot.Mom_real(:,10)) * (-1d0) & !-1 for color (arXiv:0709.2881 EQ. 5.136)
                                      * split_qiqi(alphas_dip(i_dipole),Mom_real(:,1),Mom_real(:,2),Mom_real(:,10))

            call Kinematics_VH(id,ptilde(:,:,i_dipole),NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
            do NHisto = 1,NumHistograms
              call intoHisto(NHisto,NBin(NHisto),-dip(i_dipole)*VgsWgt)
            enddo

          enddo!i_dipole
          me2sub = me2sub + me2real - dip(1) - dip(2)

!if(isnan(me2real).or.isnan(dip(1)).or.isnan(dip(2))) then
! print*,me2real,dip(1),dip(2), PSWgt_real
!print*,Mom(:,1)
!print*,Mom(:,2)
!print*,Mom(:,3)
!print*,Mom(:,4)
!print*,Mom(:,5)
!print*,Mom(:,6)
!print*,Mom(:,7)
!print*,Mom(:,8)
!print*,Mom(:,9)
!print*,Mom(:,10)
!print*,dsqrt(Mom(:,4).dot.Mom(:,4)),dsqrt(Mom(:,5).dot.Mom(:,5))
!print *,"===================="
!endif

!print*,"IR",(Mom_real(:,1).dot.Mom_real(:,10))/Ehat**2,(Mom_real(:,2).dot.Mom_real(:,10))/Ehat**2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!pause
        enddo! parton

!if(isNan(me2real))then
!  print*,me2real,dip(1),dip(2),PSWgt
!  !print*,Mom_real(:,5).dot.Mom_real(:,5),getMass(convertLHEreverse(id(8)))**2,getMass(convertLHEreverse(id(9)))**2
!print*,Mom_real(:,1).dot.Mom_real(:,1),Mom_real(1,1),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,2).dot.Mom_real(:,2),Mom_real(1,2),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,3).dot.Mom_real(:,3),Mom_real(1,3),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,4).dot.Mom_real(:,4),Mom_real(1,4),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,5).dot.Mom_real(:,5),Mom_real(1,5),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,6).dot.Mom_real(:,6),Mom_real(1,6),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,7).dot.Mom_real(:,7),Mom_real(1,7),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,8).dot.Mom_real(:,8),Mom_real(1,8),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,9).dot.Mom_real(:,9),Mom_real(1,9),(2d0*getMass(convertLHEreverse(id(8))))**2
!print*,Mom_real(:,10).dot.Mom_real(:,10),Mom_real(1,10),(2d0*getMass(convertLHEreverse(id(8))))**2
!  print*,"==================="
!  pause!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!endif

      !W
      elseif( IsAWDecay(DecayMode1) )then
        do i = -5,5 !gluon = 0, otherwise PDG codes
        do j = -5,5

          id2=id
          id2(1:2) = (/i,j/)
          !W+
          if( CouplToLHEWp(id2(1:2)) )then
            helicity(6)=sign(1d0,-dble(id2(6)))
            helicity(7)=-helicity(6)
          !W-
          elseif( CouplToLHEWm(id2(1:2)) )then
            id2(3)=-id2(3)
            id2(4)=-id2(4)
            id2(6)=-id2(6)
            id2(7)=-id2(7)
            helicity(6)=sign(1d0,-dble(id2(6)))
            helicity(7)=-helicity(6)
          else
            cycle !skip this i,j combination
          endif

          me2real=0d0
          do l=0,1
          do p=0,1
          do q=0,1
            call amp_VH_real(Mom,mass(3:5,1:2),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9),dble(q*2-1)/),id,amp_dummy)
            me2real = me2real + dble(amp_dummy*dconjg(amp_dummy))
          enddo
          enddo
          enddo
          me2real = me2real *pdf_real(LHA2M_PDF(i),1)*pdf_real(LHA2M_PDF(j),2) *PreFac_real *PostFac *aveqq *Cf *3d0 !I think I understand now.
          me2real = me2real * (4d0 * pi * alphas_real)! gs^2

          !call Kinematics_VH(id,Mom_real,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
          do NHisto = 1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),me2real*VgsWgt)
          enddo

          dip=0d0
          do i_dipole=1,2
            me2_dummy=0d0
            do p=0,1
            do q=0,1
              call amp_VH_LO(ptilde(:,:,i_dipole),mass(3:5,1:2),(/dble(p*2-1),-dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id,amp_dummy)
              me2_dummy = me2_dummy + dble(amp_dummy*dconjg(amp_dummy))
            enddo
            enddo
            me2_dummy = me2_dummy *pdf_dip(LHA2M_PDF(i),1,i_dipole)*pdf_dip(LHA2M_PDF(j),2,i_dipole) *PreFac_real *PostFac * aveqq * 3d0
            dip(i_dipole) = me2_dummy * (-1d0) / (2d0*Mom_real(:,i_dipole).dot.Mom_real(:,10)) * (-1d0) & !-1 for color (arXiv:0709.2881 EQ. 5.136)
                                      * split_qiqi(alphas_dip(i_dipole),Mom_real(:,1),Mom_real(:,2),Mom_real(:,10))

            call Kinematics_VH(id,ptilde(:,:,i_dipole),NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
            do NHisto = 1,NumHistograms
              call intoHisto(NHisto,NBin(NHisto),-dip(i_dipole)*VgsWgt)
            enddo

          enddo!i_dipole
          me2sub = me2sub + me2real - dip(1) - dip(2)
         
        enddo! parton
        enddo! parton
      endif

    endif

    !NLO: virtual + dipoles + PDF renormalization
    !interference of virtual and LO, plus "I bit" "K bit" and "P bit" from arXiv:0709.2881
    me2sup=0d0
    if(VH_PC.eq."nl".or.VH_PC.eq."sp")then
      if( IsAZDecay(DecayMode1).or.IsAPhoton(DecayMode1) )then
        do i = -5,5
          j = -i
  
          id(1:2) = (/i,j/)
  
          if (i.eq.0)cycle !skip gg
          !the finite part of the virtual diagram + other finite parts from I, K, and P bits of arXiv:0709.2881
          me2_dummy=0d0
          do p=0,1
          do q=0,1
            call amp_VH_LO(Mom(:,1:9),mass(3:5,1:2),(/dble(p*2-1),-dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id,amp_dummy)
            me2_dummy = me2_dummy + dble(amp_dummy*dconjg(amp_dummy))
          enddo
          enddo

          me2_virture_finite = me2_dummy * fac_ZH_virtual(Ehat)
          me2_virture_finite = me2_virture_finite*pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2) * PreFac * PostFac * aveqq * 3d0

          call Kinematics_VH(id,Mom,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
          do NHisto = 1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),me2_virture_finite*VgsWgt)
          enddo

          if(1d0-yRnd(17).le.yRndmin) return!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !integration over x for PDF renormalization, except for the parts included in the above (due to Dirac delta function)
          me2_PDF_ren=0d0
          do i_PDF_ren=1,2

            me2_dummy_dummy=0d0
            amp_dummy_dummy=0d0
            !print*,Ehat - M_Reso - M_V, applyPSCut
            if(applyPSCut_ren(i_PDF_ren) .or. PSWgt.eq.0d0 .or. PSWgt_ren.eq.0d0)then!may need a guard against Ehat < m_H+m_V
              cycle
            else
              do p=0,1
              do q=0,1
                call amp_VH_LO(p_PDF_ren(1:4,1:9,i_PDF_ren),mass(3:5,1:2),(/dble(p*2-1),-dble(p*2-1),helicity(3:5),dble(q*2-1),-dble(q*2-1),helicity(8:9)/),id,amp_dummy_dummy)
                me2_dummy_dummy = me2_dummy_dummy + dble(amp_dummy_dummy*dconjg(amp_dummy_dummy))
              enddo
              enddo

              me2_PDF_ren(i_PDF_ren) = IPK_ZH(Mu_Fact_ren(i_PDF_ren),alphas_ren(i_PDF_ren),me2_dummy,me2_dummy_dummy,shat,shat_ren,pdf(LHA2M_PDF(i),1:2),pdf_ren(LHA2M_PDF(i),1:2,i_PDF_ren),PreFac,PreFac_ren,yRnd(17))

              me2_PDF_ren(i_PDF_ren) = me2_PDF_ren(i_PDF_ren) * PostFac * aveqq * 3d0
!print*,"xs"
!if(i_PDF_ren.eq.2)print*,me2_PDF_ren(i_PDF_ren)
!print*,me2_PDF_ren(i_PDF_ren),pdf_ren(LHA2M_PDF(i),1:2,i_PDF_ren),i_PDF_ren
!print*,dsqrt(shat),dsqrt(shat_ren),dsqrt(shat*yRnd(17))

!              me2_PDF_ren(i_PDF_ren) = me2_PDF_ren(i_PDF_ren) &
!                                     + (me2_dummy_dummy *pdf_ren(LHA2M_PDF(i),1)*pdf_ren(LHA2M_PDF(j),2) * PreFac_ren &
!                                      - me2_dummy       *pdf(LHA2M_PDF(i),1)    *pdf(LHA2M_PDF(j),2)     * PreFac)    &
!                                     * 2d0 / (1d0-yRnd(17)) * dlog((1d0-yRnd(17))**2/yRnd(17))
!
!              me2_PDF_ren(i_PDF_ren) = me2_PDF_ren(i_PDF_ren) &
!                                     - (me2_dummy_dummy *pdf_ren(LHA2M_PDF(i),1)*pdf_ren(LHA2M_PDF(j),2) * PreFac_ren *dlog(4d0*pi*Mu_Fact**2/shat_ren) &
!                                      - me2_dummy       *pdf(LHA2M_PDF(i),1)    *pdf(LHA2M_PDF(j),2)     * PreFac     *dlog(4d0*pi*Mu_Fact**2/shat)) &
!                                     * (1d0+yRnd(17)**2) / (1d0-yRnd(17))
!
!              me2_PDF_ren(i_PDF_ren) = me2_PDF_ren(i_PDF_ren) &
!                                     +  me2_dummy_dummy *pdf_ren(LHA2M_PDF(i),1)*pdf_ren(LHA2M_PDF(j),2) * PreFac_ren &
!                                     * (-(1d0+yRnd(17))*dlog((1d0-yRnd(17))**2/yRnd(17)) + 1d0 - yRnd(17))
!
!              me2_PDF_ren(i_PDF_ren) = me2_PDF_ren(i_PDF_ren) * (alphas*Cf/2d0/pi) * PostFac * aveqq * 3d0

!print *, me2_dummy, me2_dummy_dummy

              do NHisto = 1,NumHistograms
                call intoHisto(NHisto,NBin(NHisto),me2_PDF_ren(i_PDF_ren)*VgsWgt)
              enddo

            endif
            
          enddo! x*p1 or x*p2 i_PDF_ren
          me2sup = me2sup + me2_virture_finite + me2_PDF_ren(1) + me2_PDF_ren(2)

!if(me2_PDF_ren(1).lt.-300d0)
!if(me2sup.gt.1000000d0)print*, me2sup, me2_PDF_ren(1) , me2_PDF_ren(2), me2_virture_finite, yRnd(17)



        enddo! parton
!print*,"xs"
!print*, me2sup
      elseif( IsAWDecay(DecayMode1) )then

      endif

    endif


    !gq/qg
    me2gq=0d0
    if(VH_PC.eq."qg".or.VH_PC.eq."gq".or.VH_PC.eq."nl")then
      do l=0,1
      do p=0,1
      do q=0,1
        call amp_VH_qg(Mom_real,mass(3:5,1:2),(/dble(l*2-1),-dble(l*2-1),helicity(3:5),dble(p*2-1),-dble(p*2-1),helicity(8:9),dble(q*2-1)/),id,amp_dummy)
        print*,amp_dummy
      enddo
      enddo
      enddo
    endif
pause

      !summing event weights
    EvalWeighted_VH = me2lo + me2sub + me2sup + me2gg + me2gq
    !print*,EvalWeighted_VH, me2lo , me2sub , me2sup , me2gg , me2gq
    !print*,"==================="
!    if( VH_PC.eq."nl" )then
!      EvalWeighted_VH = me2lo + me2sub + me2sup + me2gg + me2gq
!    elseif( VH_PC.eq."gg" .or. VH_PC.eq."tr" .or. VH_PC.eq."bo")then
!      EvalWeighted_VH = me2gg
!    elseif( VH_PC.eq."lo" .or. VH_PC.eq."qq" .or. VH_PC.eq."ee")then
!      EvalWeighted_VH = me2lo
!    elseif( VH_PC.eq."gq" .or. VH_PC.eq."qg")then
!      EvalWeighted_VH = me2gq
!    else
!      print *, "invalid VH_PC, VH_PC =", VH_PC
!      return
!    endif
  endif

  if(EvalWeighted_VH.lt.1d-12)return

  !JHUGen requires final state fermions being massive
  !Z/W > f f~
  if(.not.IsAPhoton(DecayMode1) .and. writeWeightedLHE) then
    PSWgt2 = s_channel_decay(Mom_save(:,4),getMass(convertLHEreverse(id(6)))**2,getMass(convertLHEreverse(id(7)))**2,yRnd(8:9),Mom_save(:,6),Mom_save(:,7))
    if(Collider.eq.1.or.Collider.eq.2)call boost2Lab(eta1,eta2,2,Mom_save(:,6:7))
    Mom(:,6:7)=Mom_save(:,6:7)
  endif

  !H > b b~
  if(HbbDecays .and. writeWeightedLHE) then
    PSWgt2 = s_channel_decay(Mom_save(:,5),getMass(convertLHEreverse(id(8)))**2,getMass(convertLHEreverse(id(9)))**2,yRnd(10:11),Mom_save(:,8),Mom_save(:,9))
    if(Collider.eq.1.or.Collider.eq.2)call boost2Lab(eta1,eta2,2,Mom_save(:,8:9))
    Mom(:,8:9)=Mom_save(:,8:9)
  endif

!  !boost from center of mass frame to lab frame
!  if(Collider.eq.1.or.Collider.eq.2)then
!    call boost2Lab(eta1,eta2,10,Mom)
!  endif
  do NHisto = 1,NumHistograms
    call intoHisto(NHisto,NBin(NHisto),(me2lo+me2gg+me2sup+me2gq)*VgsWgt)
    !me2sub was filled on the run.
  enddo

  !write events
  if(writeWeightedLHE .and. EvalWeighted_VH.gt.1d-12)then
    !consult Markus and Ulascan how id(1:2) and jet information are written!!!!!!!!!!!!!!!!!!!!!!!
    call WriteOutEvent_VH(id,helicity,Mom(:,1:9),EventWeight=EvalWeighted_VH*VgsWgt)
    !events with a jet raidated by partons may carry negative or large weight...
  endif

  AccepCounter=AccepCounter+1

  RETURN

end Function EvalWeighted_VH

Function EvalUnWeighted_VH(yRnd,genEvt,RES)
  use ModKinematics
  use ModKinematics_VH
  use ModParameters
  use ModMisc
  use ModVHLO
  use ModVHreal
  use ModVHdipole
  use ModVHgg
  use ModVHqg
  use ModPhasespace
#if compiler==1
 use ifport
#endif
  implicit none
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5): flavor in Z/W decay mode
! yrnd(6:7): phi_4 and cos(theta_4) in the CM frame of Z*(3)
! yrnd(8:9): phi_6 and cos(theta_6) in the CM frame of decay product of Z(4)
! yrnd(10:11): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
! yRnd(12): inv_mass(4)
! yRnd(13): inv_mass(5)
! yRnd(14:16): real emission phase space
! yRnd(17): NLO integration for + distribution for PDF renormalization
! yrnd(18:19): PDF mapping
! yrnd(20): flavor of j in gq > WH+j
! yrnd(21): partonic channel in unweighted events
! yRnd(22): accept or reject in unweighted events
  real(8) :: yRnd(1:22), EvalUnWeighted_VH, RES(-5:5,-5:5)
  real(8) :: pdf(-6:6,1:2)
  real(8) :: eta1, eta2, FluxFac, Ehat, sHatJacobi
  real(8) :: Mom(1:4,1:10), Mom_save(1:4,1:10), PSWgt, PSWgt2
  complex(8) :: amp_dummy
  real(8) :: me2real,me2lo,dip(1:2),me2gg,me2gq,me2sub,me2sup,me2_dummy!, lheweight(-6:6,-6:6)
  real(8) :: bound(0:121)
  integer :: i_dipole
  real(8) :: ptilde(1:4,1:9,1:2)
  integer :: i,j,ifound,jfound,k,l,p,q
  integer :: NBin(1:NumHistograms),NHisto
  real(8) :: PreFac, PostFac
  real(8) :: CS_max, sumtot
  logical :: applyPSCut
  logical :: genEVT
  real(8) :: mass(10,2)
  real(8) :: helicity(10)
  integer id(10), id2(10), idj

  include 'csmaxvalue.f'
  EvalUnWeighted_VH = 0d0
  id(:)=0
  helicity(:)=0
  mass(1:2,1:2)=0d0
  Mom=0d0
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
  mass(6:10,1:2)=0d0

  id(5)=convertLHE(Hig_)
  if(HbbDecays.eqv..true.)then
    id(8)=convertLHE(Bot_)
    id(9)=-id(8)
  else
    id(8)=Not_a_particle_
    id(9)=Not_a_particle_
  endif

  !id(10)=convertLHE(Glu_)!This may need to be overwritten later for gq!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    stop 1
    
  endif

  if(HbbDecays) PostFac=PostFac*6 !of which 2 is for summing H > bb~ helicities and 3 for summing bb~ colors

! end initialization and decaying mode related
! begin event

  if(VH_PC.ne."ee".and.VH_PC.ne."qq".and.VH_PC.ne."lo".and.VH_PC.ne."tr".and.VH_PC.ne."bo".and.VH_PC.ne."gg")then
    print*,"VH @NLO in development"
    stop 1
  endif


  if(((Collider.eq.1.or.Collider.eq.2).and.VH_PC.eq."ee") .or. &
      (Collider.eq.0                  .and.VH_PC.ne."ee"))then
    print*,"e+ e- collisions with Collider=0 only."
    print*,"VH_PC = ",VH_PC," Collider = ", Collider
    stop 1
  endif

!if e+ e- collider
  if(Collider.eq.0.and.VH_PC.eq."ee")then
    call EvalPhasespace_VH(yRnd(6:13),ILC_Energy,Mom(:,1:9),id(6:9),PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
    Mom_save=Mom
    call Kinematics_VH(id,Mom,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
    if( applyPSCut .or. PSWgt.eq.zero ) return    

    FluxFac = 1d0/(2d0*ILC_Energy**2)
    PreFac = fbGeV2 * FluxFac * PSWgt
!    EvalWeighted_VH=0d0
!    id(2)=convertLHE(ElM_)
!    id(1)=-id(2)
!    call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_dummy)
!    EvalWeighted_VH = dble(amp_dummy*dconjg(amp_dummy)) *PreFac *PostFac

!if pp/ppbar collider
  elseif(Collider.eq.1.or.Collider.eq.2)then
    !PDFMapping
    if( IsAZDecay(DecayMode1) )then
      call PDFMapping(14,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!Z
    elseif( IsAWDecay(DecayMode1) )then
      call PDFMapping(15,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!W
    elseif( IsAPhoton(DecayMode1) )then
      call PDFMapping(16,yrnd(18:19),eta1,eta2,Ehat,sHatJacobi)!Z/gamma with gamma in the final state
    endif

    !Phase space
    call EvalPhasespace_VH(yRnd(6:13),Ehat,Mom(:,1:9),id(6:9),PSWgt,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
    Mom_save=Mom
    !boost from center of mass frame to lab frame
    if(Collider.eq.1.or.Collider.eq.2)then
      call boost2Lab(eta1,eta2,9,Mom)
    endif

    !Kinematics and cuts
    call Kinematics_VH(id,Mom,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
    !p~ will have their own runs of the fuction

    if( applyPSCut .or. PSWgt.eq.zero ) return

    !Set running scales
    if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) )then
      call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
    elseif( IsAPhoton(DecayMode1) )then
      call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom_Not_a_particle(1:4) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),Not_a_particle_,convertLHEreverse(id(4)) /) )
    endif
    ! do not forget to set scales in command lines
    call EvalAlphaS()

    !PDF
    call setPDFs(eta1,eta2,pdf)
    FluxFac = 1d0/(2d0*EHat**2)
    PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt

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
      if( yRnd(21).gt.bound(k-1) .and. yRnd(21).lt.bound(k)  ) then
        ifound=i; jfound=j;
        goto 1313
      endif
    enddo
    enddo
1313 continue


    id(1:2)=(/ifound,jfound/)!to be replaced in (/0,0/)
    !ee or gg
    if(ifound.eq.0.and.jfound.eq.0)then
      if(Collider.eq.0.and.VH_PC.eq."ee")then!ee
        id(1:2)=(/convertLHE(ElP_),convertLHE(ElM_)/)
        call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_dummy)
        EvalUnweighted_VH = dble(amp_dummy*dconjg(amp_dummy)) *PreFac *PostFac

      elseif(VH_PC.eq."tr".or.VH_PC.eq."bo".or.VH_PC.eq."gg")then!gg
        id(1:2)=(/convertLHE(Glu_),convertLHE(Glu_)/)
        call amp_VH_gg(Mom(:,1:9),mass(3:5,:),helicity,id(1:9),amp_dummy)
        EvalUnweighted_VH = dble(amp_dummy*dconjg(amp_dummy)) *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0

      else
        print*,"invalid parton combination ",ifound,jfound,"for VH_PC=", VH_PC
        stop 1
      endif

    elseif(VH_PC.eq."lo".or.VH_PC.eq."qq")then
      !Z/A or W
      if( (ifound.eq.-jfound) .or. CouplToLHEWp((/ifound,jfound/)) .or. CouplToLHEWm((/ifound,jfound/)) )then
        id2=id
        !if W-, reverse sign of 3,4,6,7
        if( CouplToLHEWm((/ifound,jfound/)) )then
          id2(3)=-id2(3)
          id2(4)=-id2(4)
          id2(6)=-id2(6)
          id2(7)=-id2(7)
        endif
        helicity(6)=sign(1d0,-dble(id2(6)))
        helicity(7)=-helicity(6)
        call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id2(1:9),amp_dummy)
        me2lo = dble(amp_dummy*dconjg(amp_dummy)) *pdf(LHA2M_PDF(ifound),1)*pdf(LHA2M_PDF(jfound),2)
        EvalUnweighted_VH = me2lo *PreFac *PostFac * QuarkColAvg**2 * 3d0
        !summing 3 colors in intial qq, no factor from spins because they are casted randomly, not summed.
      else
        print*,"invalid parton combination ",ifound,jfound,"for VH_PC=", VH_PC
        stop 1
      endif

    else
      print*,"invalid parton combination or in development",ifound,jfound,"for VH_PC=", VH_PC
      stop 1
    endif




    CS_max = csmax(ifound,jfound)
    if( EvalUnWeighted_VH.gt. CS_max) then
      write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_VH, CS_max
      write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_VH, CS_max
      AlertCounter = AlertCounter + 1
      Res = 0d0
    elseif( EvalUnWeighted_VH .gt. yRnd(22)*CS_max ) then
      call Kinematics_VH(id,Mom,NBin,applyPSCut,HbbDecays,PhoOnshell=IsAPhoton(DecayMode1))
      do NHisto=1,NumHistograms
        call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
      enddo
      
      !JHUGen requires final state fermions being massive
      !Z/W > f f~
      if(.not.IsAPhoton(DecayMode1)) then
        PSWgt2 = s_channel_decay(Mom_save(:,4),getMass(convertLHEreverse(id(6)))**2,getMass(convertLHEreverse(id(7)))**2,yRnd(8:9),Mom_save(:,6),Mom_save(:,7))
      endif
      !H > b b~
      if(HbbDecays) then
        PSWgt2 = s_channel_decay(Mom_save(:,5),getMass(convertLHEreverse(id(8)))**2,getMass(convertLHEreverse(id(9)))**2,yRnd(10:11),Mom_save(:,8),Mom_save(:,9))
      endif

      if(Collider.eq.1.or.Collider.eq.2)call boost2Lab(eta1,eta2,4,Mom_save(:,6:9))

      Mom(:,6:9)=Mom_save(:,6:9)

      !write events
      if(EvalUnweighted_VH.ne.0d0)then
        !consult Markus and Ulascan how id(1:2) and jet information are written!!!!!!!!!!!!!!!!!!!!!!!
        call WriteOutEvent_VH(id,helicity,Mom(:,1:9),1d0)
        !events with a jet raidated by partons may carry negative or large weight...
      endif

      AccepCounter=AccepCounter+1

    else
      RejeCounter = RejeCounter + 1
    endif







  ELSE! NOT GENEVT =================Looking for max. weights in each partonic channel===================

    if(Collider.eq.0 .and. VH_PC.eq."ee")then! ee > ZH

      id(1:2)=(/convertLHE(ElP_),convertLHE(ElM_)/)
      call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_dummy)
      EvalUnweighted_VH = dble(amp_dummy*dconjg(amp_dummy)) *PreFac *PostFac
      !update max weight
      if (EvalUnweighted_VH.gt.csmax(0,0)) then
        csmax(0,0) = EvalUnweighted_VH
      endif

    else!hadron collisions
    
      do i=-5,5 !partonic channels
      do j=-5,5 !0 = glu, otherwise PDG code. NOT JHUGen code.
  
        !gg
        me2gg=0d0
        if((VH_PC.eq."gg".or.VH_PC.eq."bo".or.VH_PC.eq."tr").and.i.eq.0.and.j.eq.0)then
          call amp_VH_gg(Mom(:,1:9),mass(3:5,:),helicity,id(1:9),amp_dummy)
          me2gg = dble(amp_dummy*dconjg(amp_dummy)) *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0
          !update max weight
          if (me2gg.gt.csmax(i,j)) then
            csmax(i,j) = me2gg
          endif
        endif

        !lo/qq
        me2lo=0d0
        if((VH_PC.eq."qq".or.VH_PC.eq."lo"))then
          !Z/gamma
          if( (IsAZDecay(DecayMode1).or.IsAPhoton(DecayMode1)) .and. (i.eq.-j) .and. (i.ne.0) )then
            id(1:2) = (/i,j/)
            call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_dummy)
            me2lo = dble(amp_dummy*dconjg(amp_dummy)) *pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2)
            me2lo = me2lo *PreFac *PostFac * QuarkColAvg**2 * 3d0
            !summing 3 colors in intial qq, no factor from spins because they are casted randomly, not summed.
    
          !W
          elseif( IsAWDecay(DecayMode1) .and. ( CouplToLHEWp((/i,j/)).or.CouplToLHEWm((/i,j/)) ) )then
  
    
            id2=id
            id2(1:2) = (/i,j/)
            !W+
            if( CouplToLHEWp(id2(1:2)) )then
              helicity(6)=sign(1d0,-dble(id2(6)))
              helicity(7)=-helicity(6)
            !W-
            elseif( CouplToLHEWm(id2(1:2)) )then
              id2(3)=-id2(3)
              id2(4)=-id2(4)
              id2(6)=-id2(6)
              id2(7)=-id2(7)
              helicity(6)=sign(1d0,-dble(id2(6)))
              helicity(7)=-helicity(6)
            else
              print *, "invalid initial states for WH @LO. id(1:2) =", id2(1:2)
              stop 1
            endif
    
            call amp_VH_LO(Mom(:,1:9),mass(3:5,:),helicity(1:9),id2(1:9),amp_dummy)
            me2lo = dble(amp_dummy*dconjg(amp_dummy))*pdf(LHA2M_PDF(i),1)*pdf(LHA2M_PDF(j),2)
            me2lo = me2lo *PreFac *PostFac *QuarkColAvg**2 *3d0
            !summing 3 colors in intial qq, no factor from spins because they are casted randomly, not summed.
    
          endif!W
  
          !update max weight
          if (me2lo.gt.csmax(i,j)) then
            csmax(i,j) = me2lo
          endif
  
        endif!lo/qq

        !summing event weights
        EvalUnWeighted_VH = me2lo + me2gg
    
      enddo
      enddo!partonic channels

    endif!hadron collisions

  ENDIF! GENEVT

  RETURN

end Function EvalUnWeighted_VH




END MODULE ModCrossSection_VH







