!--YaofuZhou-----------------------------------------
module ModVHreal
  use ModParameters
  use ModMisc
  use ModVHaux
  implicit none
  
  public :: amp_VH_real
  private :: qqbffbHa1
  private :: qqbffbHa2
  private :: qqbffbHg4
  private :: qqbAHa1
  private :: qqbAHa2
  private :: qqbAHg4

contains

subroutine amp_VH_real(MomExt,mass,helicity,id,amp)

  implicit none
  real(8), intent(in) :: MomExt(1:4,1:10)
  real(8), intent(in) :: mass(3:5,1:2)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: amp
  integer :: id(10)

  complex(8) :: PROP3, PROP4, PROP5
  complex(8) :: gFFZ, gFFW, gFFA
  complex(8) :: gVVS1(5), gVVS2(5), gVVP(5) ! 1 = Z>ZH, 2 = A>ZH, 3 = Z>AH, 4 = A>AH, 5 = W>WH
  complex(8) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
  real(8) :: q3_q4,q3_q3,q4_q4,q5_q5
  complex(8) :: qqffbHa1, qqffbHa2, qqffbHg4, qqAHa1, qqAHa2, qqAHg4, qqZZ, qqZA, qqAZ, qqAA, qqWW
  real(8) :: sprod(5,5)
  complex(8) :: Spaa(5,5), Spbb(5,5)

! SM couplings for fermion currents
  gFFZ = ci*dsqrt(couplZffsq) ! = i * sqrt[ gwsq/4d0/(1.0_dp-xw) ]
  gFFW = ci*dsqrt(couplWffsq) ! = i sqrt[ gwsq/2.0_dp ]
  gFFA = -ci*dsqrt(couplAffsq) ! = -i sqrt[ gwsq*xw ]

  q3_q4 = -scr(MomExt(:,3),MomExt(:,4))
  q3_q3 = scr(MomExt(:,3),MomExt(:,3))
  q4_q4 = scr(MomExt(:,4),MomExt(:,4))
  q5_q5 = scr(MomExt(:,5),MomExt(:,5))

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
    call spinoru2(5,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,6),MomExt(1:4,7),MomExt(1:4,10)/),Spaa,Spbb,sprod)
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
    call spinoru2(5,(/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,1)+MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,10)/),Spaa,Spbb,sprod)!!!!!!!!!!to be determined depending on reference momentum
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
    PROP4 = PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
    amp = qqWW * PROP4

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

  amp = amp * gs * ci ! and color factor/spin factors have not been dealt with

! assemble everything and get iM
  if(id(8).ne.Not_a_particle_) then
    PROP5 = -PROPAGATOR(q5_q5,mass(5,1),mass(5,2))
    amp = amp *PROP5 &
    *(kappa*FFS(id(8), MomExt(:,8), helicity(8), id(9), MomExt(:,9), helicity(9)) &
      +kappa_tilde*FFP(id(8), MomExt(:,8), helicity(8), id(9), MomExt(:,9), helicity(9)))&
    *(-ci/vev*getMass(convertLHEreverse(id(8))))
  endif ! else H does not decay

  return
END subroutine amp_VH_real






subroutine qqbffbHa1(Spaa,Spbb,sprod,id,helicity,qqffbHa1)
  implicit none
  complex(8), intent(in) :: Spaa(1:5,1:5),Spbb(1:5,1:5)
  real(8), intent(in) :: sprod(1:5,1:5)
  integer, intent(in) :: id(10)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: qqffbHa1

  if(id(1).gt.0d0)then!q q~
    if(id(1)*helicity(1).lt.0d0)then!J1 is left
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        if(helicity(10).lt.0d0)then
          qqffbHa1 = 4d0/dsqrt(2d0)/(Spbb(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(2,3)*Spaa(2,5)*Spbb( &
                     1,2)*Spbb(1,4) + 4d0/dsqrt(2d0)/(Spbb(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(2,5 &
                     )*Spaa(3,5)*Spbb(1,4)*Spbb(1,5)!LL-
        else
          qqffbHa1 = - 4d0/(dsqrt(2d0))/(Spaa(1,5))/(Spaa(1,5)*Spbb(1,5))*Spaa(1,5)*Spaa(2,3)* &
                     Spbb(1,5)*Spbb(4,5) + 4d0/(dsqrt(2d0))/(Spaa(1,5))/(Spaa(2,5)*Spbb(2,5))* &
                     Spaa(1,2)*Spaa(2,3)*Spbb(1,4)*Spbb(2,5)!LL+
        endif
      else!J2 is right
        if(helicity(10).lt.0d0)then
          qqffbHa1 = 4d0/dsqrt(2d0)/(Spbb(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(2,4)*Spaa(2,5)*Spbb( &
                     1,2)*Spbb(1,3) + 4d0/dsqrt(2d0)/(Spbb(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(2,5 &
                     )*Spaa(4,5)*Spbb(1,3)*Spbb(1,5)!LR-
        else
          qqffbHa1 = - 4d0/(dsqrt(2d0))/(Spaa(1,5))/(Spaa(1,5)*Spbb(1,5))*Spaa(1,5)*Spaa(2,4)* &
                     Spbb(1,5)*Spbb(3,5) + 4d0/(dsqrt(2d0))/(Spaa(1,5))/(Spaa(2,5)*Spbb(2,5))* &
                     Spaa(1,2)*Spaa(2,4)*Spbb(1,3)*Spbb(2,5)!RL+
        endif
      endif
    else!J1 is right
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        if(helicity(10).lt.0d0)then
          qqffbHa1 = - 4d0/(dsqrt(2d0))/(Spbb(1,5))/(Spaa(1,5)*Spbb(1,5))*Spaa(1,5)*Spaa(3,5)* &
                     Spbb(1,5)*Spbb(2,4) + 4d0/(dsqrt(2d0))/(Spbb(1,5))/(Spaa(2,5)*Spbb(2,5))* &
                     Spaa(1,3)*Spaa(2,5)*Spbb(1,2)*Spbb(2,4)!RL-
        else
          qqffbHa1 = 4d0/dsqrt(2d0)/(Spaa(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(1,2)*Spaa(1,3)*Spbb( &
                     2,4)*Spbb(2,5) + 4d0/dsqrt(2d0)/(Spaa(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(1,3 &
                     )*Spaa(1,5)*Spbb(2,5)*Spbb(4,5)!RL+
        endif
      else!J2 is right
        if(helicity(10).lt.0d0)then
          qqffbHa1 = - 4d0/(dsqrt(2d0))/(Spbb(1,5))/(Spaa(1,5)*Spbb(1,5))*Spaa(1,5)*Spaa(4,5)* &
                     Spbb(1,5)*Spbb(2,3) + 4d0/(dsqrt(2d0))/(Spbb(1,5))/(Spaa(2,5)*Spbb(2,5))* &
                      Spaa(1,4)*Spaa(2,5)*Spbb(1,2)*Spbb(2,3)!RR-

        else
          qqffbHa1 = 4d0/dsqrt(2d0)/(Spaa(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(1,2)*Spaa(1,4)*Spbb( &
                     2,3)*Spbb(2,5) + 4d0/dsqrt(2d0)/(Spaa(1,5))/(Spaa(2,5)*Spbb(2,5))*Spaa(1,4 &
                     )*Spaa(1,5)*Spbb(2,5)*Spbb(3,5)!RR+
        endif
      endif
    endif
  else!q~ q
    if(id(1)*helicity(1).lt.0d0)then!J1 is left
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        if(helicity(10).lt.0d0)then
          qqffbHa1 = 0d0!LL-
        else
          qqffbHa1 = 0d0!LL+
        endif
      else!J2 is right
        if(helicity(10).lt.0d0)then
          qqffbHa1 = 0d0!LR-
        else
          qqffbHa1 = 0d0!LR+
        endif
      endif
    else !J1 is right
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        if(helicity(10).lt.0d0)then
          qqffbHa1 = 0d0!RL-
        else
          qqffbHa1 = 0d0!RL+
        endif
      else!J2 is right
        if(helicity(10).lt.0d0)then
          qqffbHa1 = 0d0!RR-
        else
          qqffbHa1 = 0d0!RR+
        endif
      endif
    endif
  endif

  return
end subroutine qqbffbHa1



subroutine qqbffbHa2(Spaa,Spbb,sprod,id,helicity,qqffbHa2)
  implicit none
  complex(8), intent(in) :: Spaa(1:5,1:5),Spbb(1:5,1:5)
  real(8), intent(in) :: sprod(1:5,1:5)
  integer, intent(in) :: id(10)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: qqffbHa2

  if(id(1).gt.0d0)then!q q~
    if(id(1)*helicity(1).lt.0d0)then!J1 is left
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHa2 = 0d0!LL
      else!J2 is right
        qqffbHa2 = 0d0!LR
      endif
    else!J1 is right
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHa2 = 0d0!RL
      else!J2 is right
        qqffbHa2 = 0d0!RR
      endif
    endif
  else!q~ q
    if(id(1)*helicity(1).lt.0d0)then!J1 is left
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHa2 = 0d0!LL
      else!J2 is right
        qqffbHa2 = 0d0!LR
      endif
    else!J1 is right
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHa2 = 0d0!RL
      else!J2 is right
        qqffbHa2 = 0d0!RR
      endif
    endif
  endif

  return
end subroutine qqbffbHa2


subroutine qqbffbHg4(Spaa,Spbb,sprod,id,helicity,qqffbHg4)
  implicit none
  complex(8), intent(in) :: Spaa(1:5,1:5),Spbb(1:5,1:5)
  real(8), intent(in) :: sprod(1:5,1:5)
  integer, intent(in) :: id(10)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: qqffbHg4

  if(id(1).gt.0d0)then!q q~
    if(id(1)*helicity(1).lt.0d0)then!J1 is left
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHg4 = 0d0!LL
      else!J2 is right
        qqffbHg4 = 0d0!LR
      endif
    else!J1 is right
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHg4 = 0d0!RL
      else!J2 is right
        qqffbHg4 = 0d0!RR
      endif
    endif
  else!q~ q
    if(id(1)*helicity(1).lt.0d0)then!J1 is left
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHg4 = 0d0!LL
      else!J2 is right
        qqffbHg4 = 0d0!LR
      endif
    else!J1 is right
      if(id(6)*helicity(6).lt.0d0)then!J2 is left
        qqffbHg4 = 0d0!RL
      else!J2 is right
        qqffbHg4 = 0d0!RR
      endif
    endif
  endif

  return
end subroutine qqbffbHg4



subroutine qqbAHa1(Spaa,Spbb,sprod,id,helicity,qqAHa1)
  implicit none
  complex(8), intent(in) :: Spaa(1:5,1:5),Spbb(1:5,1:5)
  real(8), intent(in) :: sprod(1:5,1:5)
  integer, intent(in) :: id(10)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: qqAHa1
  qqAHa1=0d0
  return
end subroutine qqbAHa1



subroutine qqbAHa2(Spaa,Spbb,sprod,id,helicity,qqAHa2)
  implicit none
  complex(8), intent(in) :: Spaa(1:5,1:5),Spbb(1:5,1:5)
  real(8), intent(in) :: sprod(1:5,1:5)
  integer, intent(in) :: id(10)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: qqAHa2
  qqAHa2=0d0
  return
end subroutine qqbAHa2


subroutine qqbAHg4(Spaa,Spbb,sprod,id,helicity,qqAHg4)
  implicit none
  complex(8), intent(in) :: Spaa(1:5,1:5),Spbb(1:5,1:5)
  real(8), intent(in) :: sprod(1:5,1:5)
  integer, intent(in) :: id(10)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: qqAHg4
  qqAHg4=0d0
  return
end subroutine qqbAHg4







end module ModVHreal
!!--YaofuZhou-----------------------------------------
