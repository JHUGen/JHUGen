!--YaofuZhou-----------------------------------------
module ModVHLO
  use ModParameters
  use ModMisc
  use ModVHaux
!  use ModKinematics
  implicit none
  
  public :: amp_VH_LO
  private :: qqbffbHa1
  private :: qqbffbHa2
  private :: qqbffbHg4
  private :: qqbAHa1
  private :: qqbAHa2
  private :: qqbAHg4

contains

subroutine amp_VH_LO(Mom,mass,helicity,id,amp)

  implicit none
  real(8), intent(in) :: Mom(1:4,1:9)
  real(8), intent(in) :: mass(3:5,1:2)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: amp
  integer :: id(9)

  complex(8) :: PROP3, PROP4, PROP5
  complex(8) :: gFFZ, gFFW, gFFA
  complex(8) :: gVVS1(5), gVVS2(5), gVVP(5) ! 1 = Z>ZH, 2 = A>ZH, 3 = Z>AH, 4 = A>AH, 5 = W>WH
  complex(8) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
  real(8) :: q3_q4,q3_q3,q4_q4,q5_q5
  complex(8) :: qqffbHa1, qqffbHa2, qqffbHg4, qqAHa1, qqAHa2, qqAHg4, qqZZ, qqZA, qqAZ, qqAA, qqWW
  real(8) :: sprod(4,4)
  complex(8) :: Spaa(4,4), Spbb(4,4)
  complex(8) :: za_bb(2,2), zb_bb(2,2)
  real(8) :: s_bb(2,2)!fot H>bb~

! SM couplings for fermion currents
  gFFZ = ci*dsqrt(couplZffsq) ! = i * sqrt[ gwsq/4d0/(1.0_dp-xw) ]
  gFFW = ci*dsqrt(couplWffsq) ! = i sqrt[ gwsq/2.0_dp ]
  gFFA = -ci*dsqrt(couplAffsq) ! = -i sqrt[ gwsq*xw ]

  q3_q4 = -scr(Mom(:,3),Mom(:,4))
  q3_q3 = scr(Mom(:,3),Mom(:,3))
  q4_q4 = scr(Mom(:,4),Mom(:,4))
  q5_q5 = scr(Mom(:,5),Mom(:,5))

!HVV vertex
  gVVS1=0d0
  gVVS2=0d0

  if(abs(id(3)).eq.convertLHE(Wp_))then !WH
    ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q5_q5)
    ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q5_q5)
    ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q5_q5)
    ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q5_q5)
    gVVS1(5) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVS2(5) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVP(5) = -2d0*ghz4_dyn
  else !Z/A
    ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q5_q5)
    ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q5_q5)
    ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q5_q5)
    ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q5_q5)
    gVVS1(1) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVS2(1) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVP(1) = -2d0*ghz4_dyn

    ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q3_q3,q5_q5)
    ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q3_q3,q5_q5)
    ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q3_q3,q5_q5)
    ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q3_q3,q5_q5)
    gVVS1(2) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVS2(2) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVP(2) = -2d0*ghz4_dyn

    ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q5_q5)
    ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q5_q5)
    ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q5_q5)
    ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q5_q5)
    gVVS1(3) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVS2(3) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVP(3) = -2d0*ghz4_dyn

    ghz1_dyn = czero
    ghz2_dyn = HVVSpinZeroDynamicCoupling(9,q3_q3,q4_q4,q5_q5)
    ghz3_dyn = HVVSpinZeroDynamicCoupling(10,q3_q3,q4_q4,q5_q5)
    ghz4_dyn = HVVSpinZeroDynamicCoupling(11,q3_q3,q4_q4,q5_q5)
    gVVS1(4) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVS2(4) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
    gVVP(4) = -2d0*ghz4_dyn
  endif

  gVVS1=ci*gVVS1/vev
  gVVS2=ci*gVVS2/vev
  gVVP=ci*gVVP/vev
  !gVVP=ci*gVVP * dsqrt(4d0*pi*alpha_QED) /sitW/2d0/M_W !i * gVVSP / vev
  !gVVS1=ci*gVVS1 * dsqrt(4d0*pi*alpha_QED) /sitW/2d0/M_W !i * gVVS1 / vev
  !gVVS2=ci*gVVS2 * dsqrt(4d0*pi*alpha_QED) /sitW/2d0/M_W !i * gVVS2 / vev

  qqffbHa1=0d0
  qqffbHa2=0d0
  qqffbHg4=0d0
  qqAHa1=0d0
  qqAHa2=0d0
  qqAHg4=0d0
  qqZZ=0d0
  qqAZ=0d0
  qqZA=0d0
  qqAA=0d0
  qqWW=0d0
  amp=0d0

! amplitudes w/o HVV or VFF couplings to calculate
  if(id(6).ne.convertLHE(Pho_))then !Z/W/A decays
    if( id(1).gt.0 .and. id(2).lt.0 )then!q q~
      call spinoru2(4,(/Mom(1:4,1),Mom(1:4,2),Mom(1:4,6),Mom(1:4,7)/),Spaa,Spbb,sprod)
    elseif( id(1).lt.0 .and. id(2).gt.0 )then!q~ q
      call spinoru2(4,(/Mom(1:4,2),Mom(1:4,1),Mom(1:4,6),Mom(1:4,7)/),Spaa,Spbb,sprod)
    else
      print*,"no qq / q~q~ initiated processes for VH"
      return
    endif
    if(gVVS1(1).ne.0d0.or.gVVS1(2).ne.0d0.or.gVVS1(3).ne.0d0.or.gVVS1(4).ne.0d0.or.gVVS1(5).ne.0d0)then
      call qqbffbHa1(Spaa,Spbb,sprod,id,helicity,qqffbHa1)
    endif
    if(gVVS2(1).ne.0d0.or.gVVS2(2).ne.0d0.or.gVVS2(3).ne.0d0.or.gVVS2(4).ne.0d0.or.gVVS2(5).ne.0d0)then
      call qqbffbHa2(Spaa,Spbb,sprod,id,helicity,qqffbHa2)
    endif
    if(gVVP(1).ne.0d0.or.gVVP(2).ne.0d0.or.gVVP(3).ne.0d0.or.gVVP(4).ne.0d0.or.gVVP(5).ne.0d0)then
      call qqbffbHg4(Spaa,Spbb,sprod,id,helicity,qqffbHg4)
    endif
  else ! A is final state
    if( id(1).gt.0 .and. id(2).lt.0 )then!q q~
      call spinoru2(4,(/Mom(1:4,1),Mom(1:4,2),Mom(1:4,2),Mom(1:4,4)/),Spaa,Spbb,sprod)!reference momentum to be set!!!!!!!!!!!!!!!
    elseif( id(1).lt.0 .and. id(2).gt.0 )then!q~ q
      call spinoru2(4,(/Mom(1:4,2),Mom(1:4,1),Mom(1:4,2),Mom(1:4,4)/),Spaa,Spbb,sprod)!reference momentum to be set!!!!!!!!!!!!!!!
    else
      print*,"no qq / q~q~ initiated processes for AH"
      return
    endif
    if(gVVS1(3).ne.0d0.or.gVVS1(4).ne.0d0)then
      call qqbAHa1(Spaa,Spbb,sprod,id,helicity,qqAHa1)
    endif
    if(gVVS2(3).ne.0d0.or.gVVS2(4).ne.0d0)then
      call qqbAHa2(Spaa,Spbb,sprod,id,helicity,qqAHa2)
    endif
    if(gVVP(3).ne.0d0.or.gVVP(4).ne.0d0)then
      call qqbAHg4(Spaa,Spbb,sprod,id,helicity,qqAHg4)
    endif
  endif

  if(abs(id(3)).eq.convertLHE(Wp_))then !WH
    !incoming current w/ W coupling
    qqWW = qqffbHa1*gVVS1(5) + qqffbHa2*gVVS2(5) + qqffbHg4*gVVP(5)
    if((id(1)*helicity(1)).le.0d0.and.(id(6)*helicity(6)).le.0d0)then
      qqWW=qqWW*gFFW**2*CKMbare(id(1),id(2))*CKM(id(6),id(7))
    else
      qqWW=0d0
    endif
    PROP3 = PROPAGATOR(q3_q3,mass(3,1),mass(3,2))
    PROP4 = PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
    amp = qqWW * PROP3 * PROP4

  else ! Z or A
  ! HVV vertex
    if(id(6).ne.convertLHE(Pho_))then !Z/A decays
      qqZZ = qqffbHa1*gVVS1(1) + qqffbHa2*gVVS2(1) + qqffbHg4*gVVP(1)
      qqAZ = qqffbHa1*gVVS1(2) + qqffbHa2*gVVS2(2) + qqffbHg4*gVVP(2)
      qqZA = qqffbHa1*gVVS1(3) + qqffbHa2*gVVS2(3) + qqffbHg4*gVVP(3)
      qqAA = qqffbHa1*gVVS1(4) + qqffbHa2*gVVS2(4) + qqffbHg4*gVVP(4)
      if(qqZZ.ne.0d0)then
        PROP3 = PROPAGATOR(q3_q3,mass(3,1),mass(3,2))
        PROP4 = PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
        qqZZ=qqZZ*ZFFbare(id(1),id(2),helicity(1),helicity(2))*ZFF(id(6),id(7),helicity(6),helicity(7))
        qqZZ=qqZZ*PROP3*PROP4*gFFZ**2
      endif
      if(qqZA.ne.0d0)then
        PROP3 = PROPAGATOR(q3_q3,mass(3,1),mass(3,2))
        PROP4 = PROPAGATOR(q4_q4,0d0,0d0)
        qqZA=qqZA*ZFFbare(id(1),id(2),helicity(1),helicity(2))*AFF(id(6),id(7),helicity(6),helicity(7))
        qqZA=qqZA*PROP3*PROP4*gFFZ*gFFA
      endif
      if(qqAZ.ne.0d0)then
        PROP3 = PROPAGATOR(q3_q3,0d0,0d0)
        PROP4 = PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
        qqAZ=qqAZ*AFF(id(1),id(2),helicity(1),helicity(2))*ZFF(id(6),id(7),helicity(6),helicity(7))
        qqAZ=qqAZ*PROP3*PROP4*gFFA*gFFZ
      endif
      if(qqAA.ne.0d0)then
        PROP3 = PROPAGATOR(q3_q3,0d0,0d0)
        PROP4 = PROPAGATOR(q4_q4,0d0,0d0)
        qqAA=qqAA*AFF(id(1),id(2),helicity(1),helicity(2))*AFF(id(6),id(7),helicity(6),helicity(7))
        qqAA=qqAA*PROP3*PROP4*gFFA*gFFA
      endif
    else !A is final state
      qqZA = qqAHa1*gVVS1(3) + qqAHa2*gVVS2(3) + qqAHg4*gVVP(3)
      qqAA = qqAHa1*gVVS1(4) + qqAHa2*gVVS2(4) + qqAHg4*gVVP(3)
      if(qqAZ.ne.0d0)then
        PROP3 = PROPAGATOR(q3_q3,0d0,0d0)
        PROP4 = PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
        qqAZ=qqAZ*AFF(id(1),id(2),helicity(1),helicity(2))
        qqAZ=qqAZ*PROP3*PROP4*gFFA
      endif
      if(qqAA.ne.0d0)then
        PROP3 = PROPAGATOR(q3_q3,0d0,0d0)
        PROP4 = PROPAGATOR(q4_q4,0d0,0d0)
        qqAA=qqAA*AFF(id(1),id(2),helicity(1),helicity(2))
        qqAA=qqAA*PROP3*PROP4*gFFA
      endif
    endif

    amp = qqZZ+qqAZ+qqZA+qqAA

  endif


! assemble everything and get iM
  if(id(8).ne.Not_a_particle_) then
    if(get_minv(Mom(:,5)).gt.0d0)then
      PROP5 = -PROPAGATOR(q5_q5,mass(5,1),mass(5,2))
      amp=amp *PROP5 &
      *(kappa*FFS(id(8), Mom(:,8), helicity(8), id(9), Mom(:,9), helicity(9)) &
        +ci*kappa_tilde*FFP(id(8), Mom(:,8), helicity(8), id(9), Mom(:,9), helicity(9)))&
      *(-ci/vev*massfrun(getMass(convertLHEreverse(id(8))),get_minv(Mom(:,5))))
    else
      amp=0d0
    endif
  endif ! else H does not decay

  return
END subroutine amp_VH_LO



subroutine qqbffbHa1(Spaa,Spbb,sprod,id,helicity,qqffbHa1)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: qqffbHa1

  !if(id(1)*helicity(1).lt.0d0)then!J1 is left
!  if(dble(id(1))*dble(id(2))*helicity(1)*helicity(2).lt.0d0)then
!    qqffbHa1 = 0d0
!  elseif((id(1).gt.0 .and. helicity(1).lt.0d0 .and. id(2).lt.0 .and. helicity(2).gt.0d0) .or. &
!         (id(2).gt.0 .and. helicity(2).gt.0d0 .and. id(1).lt.0 .and. helicity(1).lt.0d0) )then!J1 is left
!    if(id(6)*helicity(6).lt.0d0)then!J2 is left
!      qqffbHa1 = - 2d0*Spaa(2,3)*Spbb(1,4)!LL
!    else!J2 is right
!      qqffbHa1 = - 2d0*Spaa(2,4)*Spbb(1,3)!LR
!    endif
!  elseif((id(2).gt.0 .and. helicity(2).lt.0d0 .and. id(1).lt.0 .and. helicity(1).gt.0d0) .or. &
!         (id(1).gt.0 .and. helicity(1).gt.0d0 .and. id(2).lt.0 .and. helicity(2).lt.0d0) )then!J1 is right
!    if(id(6)*helicity(6).lt.0d0)then!J2 is left
!      qqffbHa1 = - 2d0*Spaa(1,3)*Spbb(2,4)!RL
!    else!J2 is right
!      qqffbHa1 = - 2d0*Spaa(1,4)*Spbb(2,3)!RR
!    endif
!  endif

  if(    id(1)*helicity(1).lt.0d0 .and. id(2)*helicity(2).lt.0d0)then!J1 is left
    if(id(6)*helicity(6).lt.0d0)then!J2 is left
      qqffbHa1 = - 2d0*Spaa(2,3)*Spbb(1,4)!LL
    else!J2 is right
      qqffbHa1 = - 2d0*Spaa(2,4)*Spbb(1,3)!LR
    endif
  elseif(id(1)*helicity(1).gt.0d0 .and. id(2)*helicity(2).gt.0d0)then!J1 is right
    if(id(6)*helicity(6).lt.0d0)then!J2 is left
      qqffbHa1 = - 2d0*Spaa(1,3)*Spbb(2,4)!RL
    else!J2 is right
      qqffbHa1 = - 2d0*Spaa(1,4)*Spbb(2,3)!RR
    endif
  else
    qqffbHa1 = 0d0
  endif

!  if(    (id(1).gt.0 .and. helicity(1).lt.0d0 .and. id(2).lt.0 .and. helicity(2).gt.0d0) .or. &
!         (id(2).gt.0 .and. helicity(2).gt.0d0 .and. id(1).lt.0 .and. helicity(1).lt.0d0) )then!J1 is left
!    if(id(6)*helicity(6).lt.0d0)then!J2 is left
!      qqffbHa1 = - 2d0*Spaa(2,3)*Spbb(1,4)!LL
!    else!J2 is right
!      qqffbHa1 = - 2d0*Spaa(2,4)*Spbb(1,3)!LR
!    endif
!  elseif((id(2).gt.0 .and. helicity(2).lt.0d0 .and. id(1).lt.0 .and. helicity(1).gt.0d0) .or. &
!         (id(1).gt.0 .and. helicity(1).gt.0d0 .and. id(2).lt.0 .and. helicity(2).lt.0d0) )then!J1 is right
!    if(id(6)*helicity(6).lt.0d0)then!J2 is left
!      qqffbHa1 = - 2d0*Spaa(1,3)*Spbb(2,4)!RL
!    else!J2 is right
!      qqffbHa1 = - 2d0*Spaa(1,4)*Spbb(2,3)!RR
!    endif
!  else
!    qqffbHa1 = 0d0
!  endif

  !qqffbHa1=1d0

!  if(id(1).gt.0 .and. id(2).lt.0)then
!    if( helicity(1).lt.0d0 .and. helicity(2).gt.0d0 )then!J1 is left
!      if(id(6)*helicity(6).lt.0d0)then!J2 is left
!        qqffbHa1 = - 2d0*Spaa(2,3)*Spbb(1,4)!LL
!      else!J2 is right
!        qqffbHa1 = - 2d0*Spaa(2,4)*Spbb(1,3)!LR
!      endif
!    elseif( helicity(1).gt.0d0 .and. helicity(2).lt.0d0 )then!J1 is right
!      if(id(6)*helicity(6).lt.0d0)then!J2 is left
!        qqffbHa1 = - 2d0*Spaa(1,3)*Spbb(2,4)!RL
!      else!J2 is right
!        qqffbHa1 = - 2d0*Spaa(1,4)*Spbb(2,3)!RR
!      endif
!    else
!      qqffbHa1 = 0d0
!    endif
!  elseif(id(1).lt.0 .and. id(2).gt.0)then
!    if( helicity(1).lt.0d0 .and. helicity(2).gt.0d0 )then!J1 is left
!      if(id(6)*helicity(6).lt.0d0)then!J2 is left
!        qqffbHa1 = - 2d0*Spaa(1,3)*Spbb(2,4)!LL
!      else!J2 is right
!        qqffbHa1 = - 2d0*Spaa(2,3)*Spbb(1,4)!LR
!      endif
!    elseif( helicity(1).gt.0d0 .and. helicity(2).lt.0d0 )then!J1 is right
!      if(id(6)*helicity(6).lt.0d0)then!J2 is left
!        qqffbHa1 = - 2d0*Spaa(1,4)*Spbb(2,3)!RL
!      else!J2 is right
!        qqffbHa1 = - 2d0*Spaa(2,4)*Spbb(1,3)!RR
!      endif
!    else
!      qqffbHa1 = 0d0
!    endif
!  endif

  return
end subroutine qqbffbHa1



subroutine qqbffbHa2(Spaa,Spbb,sprod,id,helicity,qqffbHa2)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: qqffbHa2

  if(    id(1)*helicity(1).lt.0d0 .and. id(2)*helicity(2).lt.0d0)then!J1 is left
    if(id(6)*helicity(6).lt.0d0)then!J2 is left
      qqffbHa2 = - Spaa(1,3)*Spaa(2,3)*Spbb(1,3)*Spbb(1,4) &
                 - Spaa(1,3)*Spaa(2,4)*Spbb(1,4)**2 &
                 - Spaa(2,3)**2*Spbb(1,3)*Spbb(2,4) &
                 - Spaa(2,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,4)!LL
    else!J2 is right
      qqffbHa2 = - Spaa(1,4)*Spaa(2,3)*Spbb(1,3)**2 &
                 - Spaa(1,4)*Spaa(2,4)*Spbb(1,3)*Spbb(1,4) &
                 - Spaa(2,3)*Spaa(2,4)*Spbb(1,3)*Spbb(2,3) &
                 - Spaa(2,4)**2*Spbb(1,4)*Spbb(2,3) !LR
    endif
  elseif(id(1)*helicity(1).gt.0d0 .and. id(2)*helicity(2).gt.0d0)then!J1 is right
    if(id(6)*helicity(6).lt.0d0)then!J2 is left
      qqffbHa2 = - Spaa(1,3)**2*Spbb(1,4)*Spbb(2,3) &
                 - Spaa(1,3)*Spaa(1,4)*Spbb(1,4)*Spbb(2,4) &
                 - Spaa(1,3)*Spaa(2,3)*Spbb(2,3)*Spbb(2,4) &
                 - Spaa(1,4)*Spaa(2,3)*Spbb(2,4)**2!RL
    else!J2 is right
      qqffbHa2 = - Spaa(1,3)*Spaa(1,4)*Spbb(1,3)*Spbb(2,3) &
                 - Spaa(1,3)*Spaa(2,4)*Spbb(2,3)**2 &
                 - Spaa(1,4)**2*Spbb(1,3)*Spbb(2,4) &
                 - Spaa(1,4)*Spaa(2,4)*Spbb(2,3)*Spbb(2,4)!RR
    endif
  else
    qqffbHa2 = 0d0
  endif

  return
end subroutine qqbffbHa2


subroutine qqbffbHg4(Spaa,Spbb,sprod,id,helicity,qqffbHg4)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: qqffbHg4

  if(    id(1)*helicity(1).lt.0d0 .and. id(2)*helicity(2).lt.0d0)then!J1 is left
    if(id(6)*helicity(6).lt.0d0)then!J2 is left
      qqffbHg4 = - Spaa(1,3)*Spaa(2,4)*Spbb(1,4)**2 &
                 + Spaa(1,4)*Spaa(2,3)*Spbb(1,4)**2 &
                 + Spaa(2,3)**2*Spbb(1,3)*Spbb(2,4) &
                 - Spaa(2,3)**2*Spbb(1,4)*Spbb(2,3)!LL
    else!J2 is right
      qqffbHg4 =   Spaa(1,3)*Spaa(2,4)*Spbb(1,3)**2 &
                 - Spaa(1,4)*Spaa(2,3)*Spbb(1,3)**2 &
                 - Spaa(2,4)**2*Spbb(1,3)*Spbb(2,4) &
                 + Spaa(2,4)**2*Spbb(1,4)*Spbb(2,3)!LR
    endif
  elseif(id(1)*helicity(1).gt.0d0 .and. id(2)*helicity(2).gt.0d0)then!J1 is right
    if(id(6)*helicity(6).lt.0d0)then!J2 is left
      qqffbHg4 =   Spaa(1,2)*Spaa(3,4)*Spbb(2,4)**2 &
                 - Spaa(1,3)**2*Spbb(1,2)*Spbb(3,4)!RL
    else!J2 is right
      qqffbHg4 = - Spaa(1,3)*Spaa(2,4)*Spbb(2,3)**2 &
                 + Spaa(1,4)**2*Spbb(1,3)*Spbb(2,4) &
                 - Spaa(1,4)**2*Spbb(1,4)*Spbb(2,3) &
                 + Spaa(1,4)*Spaa(2,3)*Spbb(2,3)**2!RR
    endif
  else
    qqffbHg4 =0d0
  endif

  qqffbHg4 = qqffbHg4 * ci!so that it is in phase with Process=50

  return
end subroutine qqbffbHg4



subroutine qqbAHa1(Spaa,Spbb,sprod,id,helicity,qqAHa1)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: qqAHa1

  qqAHa1=0d0

  return
end subroutine qqbAHa1



subroutine qqbAHa2(Spaa,Spbb,sprod,id,helicity,qqAHa2)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: qqAHa2

  qqAHa2=0d0

  return
end subroutine qqbAHa2


subroutine qqbAHg4(Spaa,Spbb,sprod,id,helicity,qqAHg4)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: qqAHg4

  qqAHg4=0d0

  return
end subroutine qqbAHg4





end module ModVHLO
!!--YaofuZhou-----------------------------------------
