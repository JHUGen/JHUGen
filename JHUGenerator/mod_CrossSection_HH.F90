MODULE ModCrossSection_HH
  use ModVHaux
  use ModHH
  implicit none
  integer, parameter,private :: LHA2M_pdf(-6:6) = (/-5,-6,-3,-4,-1,-2,0 ,2,1,4,3,6,5/)
  integer, parameter,private :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)  

 CONTAINS





Function EvalWeighted_HH(yRnd,VgsWgt)
  use ModKinematics
  use ModKinematics_HH
  use ModParameters
  use ModMisc
  use ModPhasespace
#if compiler==1
 use ifport
#endif
  implicit none
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5:6): phi_4 and cos(theta_4) in the CM frame of 1+2(3)
! yrnd(7:8): phi_6 and cos(theta_6) in the CM frame of decay product of H(4)
! yrnd(9:10): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
! yRnd(11): inv_mass(4)
! yRnd(12): inv_mass(5)
! yRnd(13): swap momenta in PS for stability
! yrnd(14:15): PDF mapping
! yrnd(16): partonic channel in unweighted events
! yRnd(17): accept or reject in unweighted events
  real(8) :: yRnd(1:15),VgsWgt, EvalWeighted_HH
  real(8) :: pdf(-6:6,1:2)
  real(8) :: eta1, eta2, FluxFac, Ehat, shat, sHatJacobi
  real(8) :: Mom(1:4,1:9), Mom_swap(1:4,1:9), Mom_save(1:4,1:9),PSWgt, PSWgt2
  complex(8) :: amp, amp_swap
  real(8) :: me2
  integer :: i,j,k,l
  integer :: NBin(1:NumHistograms),NHisto
  real(8) :: PreFac, PostFac
  logical :: applyPSCut
  real(8) :: mass(1:9,1:2)
  real(8) :: helicity(9)
  integer id(9)
  
  EvalWeighted_HH=0d0
  EvalCounter = EvalCounter+1

! initialization and decaying mode related
  id(3)=convertLHE(Hig_)
  id(4)=convertLHE(Hig_)
  id(5)=convertLHE(Hig_)
  id(6)=convertLHE(Bot_)
  id(7)=convertLHE(Bot_)
  id(8)=convertLHE(Bot_)
  id(9)=convertLHE(Bot_)
  helicity(:)=0
  Mom=0d0
  mass(1:2,1:2)=0d0
  mass(3,1)=M_Reso
  mass(3,2)=Ga_Reso
  mass(4,1)=M_Reso
  mass(4,2)=Ga_Reso
  mass(5,1)=M_Reso
  mass(5,2)=Ga_Reso
  mass(6:9,1:2)=0d0


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
  helicity(9)=-helicity(8)
!toss coin and decide particle 6,7 helicities
  if (yRnd(3).gt.0.5d0)then
    helicity(6)=1d0
  else
    helicity(6)=-1d0
  endif
  helicity(7)=-helicity(6)

  PostFac=18d0 !of which 3^2 for summing bb~ colors, and 2 is for summing H > bb~ helicities (not 2^2 beause h(67)=h(89) or not are separate cases.)

  call PDFMapping(17,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
!call PDFMapping(1,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
  if(Ehat.le.(2d0*M_Reso-10d0*Ga_Reso))return!!!!may be improved
  call EvalPhasespace_HH(yRnd(5:13),Ehat,Mom(:,1:9),id(6:9),PSWgt)
  !call EvalPhasespace_HH_old(yRnd(5:13),Ehat,Mom(:,1:9),mass(3:5,1:2),PSWgt)
  Mom_save=Mom
  call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
  call EvalAlphaS()
  call setPDFs(eta1,eta2,pdf)
!print*,pdf(0,1),pdf(0,2)
  call boost2Lab(eta1,eta2,9,Mom)
  call Kinematics_HH(id,Mom,NBin,applyPSCut)
  if( applyPSCut .or. PSWgt.eq.zero ) return    
  FluxFac = 1d0/(2d0*Ehat**2)
  PreFac = hbarc2XsecUnit * FluxFac * PSWgt

  id(1:2)=(/convertLHE(Glu_),convertLHE(Glu_)/)
  amp=0d0
  amp_swap=0d0
!print*,PSWgt
  call amp_HH(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp)
!print*,amp
  if(includeInterference.eqv..true.)then
    if((helicity(6).eq.helicity(8)).and.(helicity(7).eq.helicity(9)))then
      Mom_swap(:,1:9)=Mom
      Mom_swap(:,6)=Mom(:,6)
      Mom_swap(:,7)=Mom(:,9)
      Mom_swap(:,8)=Mom(:,8)
      Mom_swap(:,9)=Mom(:,7)
      call amp_HH(Mom_swap(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_swap)
      amp=amp+amp_swap
    endif
  endif

!amp=1d5
  EvalWeighted_HH = dble(amp*dconjg(amp)) *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0
!EvalWeighted_HH = dble(amp*dconjg(amp)) *FluxFac*PSWgt

!if(EvalWeighted_HH.gt.10000000000d0)then
!if(EvalWeighted_HH*VgsWgt.gt.100d0)then
!print*,FluxFac,PSWgt,VgsWgt
!print*,Get_MInv(Mom(:,4)),Get_MInv(Mom(:,5))
!print*,EvalWeighted_HH*VgsWgt
!print*,"=============="
!pause
!endif
  if(EvalWeighted_HH.lt.1d-12)return

  do NHisto = 1,NumHistograms
    call intoHisto(NHisto,NBin(NHisto),EvalWeighted_HH*VgsWgt)
  enddo

  !write events
  !JHUGen requires final state fermions being massive
  if(writeWeightedLHE)then
    PSWgt2 = s_channel_decay(Mom_save(:,4),getMass(convertLHEreverse(id(6)))**2,getMass(convertLHEreverse(id(7)))**2,yRnd(7:8),Mom_save(:,6),Mom_save(:,7))
    PSWgt2 = s_channel_decay(Mom_save(:,5),getMass(convertLHEreverse(id(8)))**2,getMass(convertLHEreverse(id(9)))**2,yRnd(9:10),Mom_save(:,8),Mom_save(:,9))
    call boost2Lab(eta1,eta2,4,Mom_save(:,6:9))
    Mom(:,6:9)=Mom_save(:,6:9)
    call WriteOutEvent_HH(id,helicity,Mom(:,1:9),EventWeight=EvalWeighted_HH*VgsWgt)
  endif  

  AccepCounter=AccepCounter+1

  RETURN

end Function EvalWeighted_HH

Function EvalUnWeighted_HH(yRnd,genEvt,RES)
  use ModKinematics
  use ModKinematics_HH
  use ModParameters
  use ModMisc
  use ModPhasespace
#if compiler==1
 use ifport
#endif
  implicit none
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5:6): phi_4 and cos(theta_4) in the CM frame of 1+2(3)
! yrnd(7:8): phi_6 and cos(theta_6) in the CM frame of decay product of H(4)
! yrnd(9:10): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
! yRnd(11): inv_mass(4)
! yRnd(12): inv_mass(5)
! yRnd(13): swap momenta in PS for stability
! yrnd(14:15): PDF mapping
! yrnd(16): partonic channel in unweighted events
! yRnd(17): accept or reject in unweighted events
  real(8) :: yRnd(1:17),VgsWgt, EvalUnWeighted_HH, RES(-5:5,-5:5)
  real(8) :: pdf(-6:6,1:2)
  real(8) :: eta1, eta2, FluxFac, Ehat, shat, sHatJacobi
  real(8) :: Mom(1:4,1:9), Mom_swap(1:4,1:9), Mom_save(1:4,1:9),PSWgt, PSWgt2
  complex(8) :: amp, amp_swap
  real(8) :: me2
  integer :: i,j,k,l,ifound,jfound
  real(8) :: bound(0:121)
  integer :: NBin(1:NumHistograms),NHisto
  real(8) :: PreFac, PostFac
  logical :: applyPSCut
  real(8) :: mass(1:9,1:2)
  real(8) :: helicity(9)
  integer id(9)
  real(8) :: CS_max, sumtot
  logical :: genEVT

  include 'csmaxvalue.f'
  EvalUnWeighted_HH = 0d0
  EvalCounter = EvalCounter+1

! initialization and decaying mode related
  id(3)=convertLHE(Hig_)
  id(4)=convertLHE(Hig_)
  id(5)=convertLHE(Hig_)
  id(6)=convertLHE(Bot_)
  id(7)=convertLHE(Bot_)
  id(8)=convertLHE(Bot_)
  id(9)=convertLHE(Bot_)
  helicity(:)=0
  Mom=0d0
  mass(1:2,1:2)=0d0
  mass(3,1)=M_Reso
  mass(3,2)=Ga_Reso
  mass(4,1)=M_Reso
  mass(4,2)=Ga_Reso
  mass(5,1)=M_Reso
  mass(5,2)=Ga_Reso
  mass(6:9,1:2)=0d0


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

  PostFac=18d0 !of which 3^2 for summing bb~ colors, and 2 is for summing H > bb~ helicities (not 2^2 beause h(67)=h(89) or not are separate cases.)

  call PDFMapping(1,yrnd(14:15),eta1,eta2,Ehat,sHatJacobi)
  !if(Ehat.le.(2d0*M_Reso))return!!!!may be improved
  call EvalPhasespace_HH(yRnd(5:13),Ehat,Mom(:,1:9),id(6:9),PSWgt)
  Mom_save=Mom
  call SetRunningScales( (/ Mom(1:4,5),Mom(1:4,6),Mom(1:4,7) /) , (/ convertLHEreverse(id(3)),convertLHEreverse(id(6)),convertLHEreverse(id(7)),convertLHEreverse(id(4)) /) )
  call EvalAlphaS()
  call setPDFs(eta1,eta2,pdf)
  call boost2Lab(eta1,eta2,9,Mom)
  call Kinematics_HH(id,Mom,NBin,applyPSCut)
  if( applyPSCut .or. PSWgt.eq.zero ) return    
  FluxFac = 1d0/(2d0*Ehat**2)
  PreFac = hbarc2XsecUnit * FluxFac * PSWgt

  
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


  id(1:2)=(/convertLHE(Glu_),convertLHE(Glu_)/)
  amp=0d0
  amp_swap=0d0
  call amp_HH(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp)
  if((helicity(6).eq.helicity(8)).and.(helicity(7).eq.helicity(9)))then
    Mom_swap(:,1:9)=Mom
    Mom_swap(:,6)=Mom(:,8)
    Mom_swap(:,7)=Mom(:,9)
    Mom_swap(:,8)=Mom(:,6)
    Mom_swap(:,9)=Mom(:,7)
    call amp_HH(Mom_swap(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_swap)
    amp=amp+amp_swap
  endif
  EvalUnWeighted_HH = dble(amp*dconjg(amp)) *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0

  if(EvalUnWeighted_HH.lt.1d-12)return

  CS_max = csmax(ifound,jfound)
  if( EvalUnWeighted_HH.gt. CS_max) then
    write(io_stdout,"(2X,A,1PE13.6,1PE13.6)")  "CS_max is too small.",EvalUnWeighted_HH, CS_max
    write(io_LogFile,"(2X,A,1PE13.6,1PE13.6)") "CS_max is too small.",EvalUnWeighted_HH, CS_max
    AlertCounter = AlertCounter + 1
    Res = 0d0
  
  elseif( EvalUnWeighted_HH .gt. yRnd(17)*CS_max ) then
    call Kinematics_HH(id,Mom,NBin,applyPSCut)
    do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),1d0)  ! CS_Max is the integration volume
    enddo
    !JHUGen requires final state fermions being massive
    PSWgt2 = s_channel_decay(Mom_save(:,4),getMass(convertLHEreverse(id(6)))**2,getMass(convertLHEreverse(id(7)))**2,yRnd(7:8),Mom_save(:,6),Mom_save(:,7))
    PSWgt2 = s_channel_decay(Mom_save(:,5),getMass(convertLHEreverse(id(8)))**2,getMass(convertLHEreverse(id(9)))**2,yRnd(9:10),Mom_save(:,8),Mom_save(:,9))
    call boost2Lab(eta1,eta2,4,Mom_save(:,6:9))
    Mom(:,6:9)=Mom_save(:,6:9)
    call WriteOutEvent_HH(id,helicity,Mom(:,1:9),1d0)
    AccepCounter=AccepCounter+1
  
  else
    RejeCounter = RejeCounter + 1
  endif







  ELSE! NOT GENEVT =================Looking for max. weights in each partonic channel===================


  id(1:2)=(/convertLHE(Glu_),convertLHE(Glu_)/)
  call amp_HH(Mom(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp)
  if((helicity(6).eq.helicity(8)).and.(helicity(7).eq.helicity(9)))then
    Mom_swap(:,1:9)=Mom
    Mom_swap(:,6)=Mom(:,8)
    Mom_swap(:,7)=Mom(:,9)
    Mom_swap(:,8)=Mom(:,6)
    Mom_swap(:,9)=Mom(:,7)
    call amp_HH(Mom_swap(:,1:9),mass(3:5,:),helicity(1:9),id(1:9),amp_swap)
    amp=amp+amp_swap
  endif
  EvalUnWeighted_HH = dble(amp*dconjg(amp)) *pdf(0,1)*pdf(0,2) *PreFac *PostFac *GluonColAvg**2 *8d0
  !update max weight
  if (EvalUnweighted_HH.gt.csmax(0,0)) then
    csmax(0,0) = EvalUnweighted_HH
  endif

   

  ENDIF! GENEVT

  RETURN

end Function EvalUnWeighted_HH




END MODULE ModCrossSection_HH







