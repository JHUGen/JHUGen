!--YaofuZhou-----------------------------------------
module ModVHgg
  use ModParameters
!  use ModKinematics
  use ModMisc
  use ModVHaux
  use ModffbHBox7_8mmm
  use ModffbHBox7_8ppm
  use ModffbHBox7_8mmp
  use ModffbHBox7_8ppp
  use ModffbHBox9_10mmm
  use ModffbHBox9_10ppm
  use ModffbHBox9_10mmp
  use ModffbHBox9_10ppp
  use ModffbHBox11_12mmm
  use ModffbHBox11_12ppm
  use ModffbHBox11_12mmp
  use ModffbHBox11_12ppp
  use ModffbHBox7_8mpm
  use ModffbHBox7_8mpp
  use ModffbHBox9_10mpm
  use ModffbHBox9_10mpp
  use ModffbHBox11_12mpm
  use ModffbHBox11_12mpp
  use ModffbHBox7_8pmm
  use ModffbHBox7_8pmp
  use ModffbHBox9_10pmm
  use ModffbHBox9_10pmp
  use ModffbHBox11_12pmm
  use ModffbHBox11_12pmp

  use ModffbH5Box7_8mmm
  use ModffbH5Box7_8ppm
  use ModffbH5Box7_8mmp
  use ModffbH5Box7_8ppp
  use ModffbH5Box9_10mmm
  use ModffbH5Box9_10ppm
  use ModffbH5Box9_10mmp
  use ModffbH5Box9_10ppp
  use ModffbH5Box11_12mmm
  use ModffbH5Box11_12ppm
  use ModffbH5Box11_12mmp
  use ModffbH5Box11_12ppp
  use ModffbH5Box7_8mpm
  use ModffbH5Box7_8mpp
  use ModffbH5Box9_10mpm
  use ModffbH5Box9_10mpp
  use ModffbH5Box11_12mpm
  use ModffbH5Box11_12mpp
  use ModffbH5Box7_8pmm
  use ModffbH5Box7_8pmp
  use ModffbH5Box9_10pmm
  use ModffbH5Box9_10pmp
  use ModffbH5Box11_12pmm
  use ModffbH5Box11_12pmp

  use Collier
  implicit none
  
  public :: amp_VH_gg

contains

!#if linkMELA==1
subroutine amp_VH_gg(Mom,mass,helicity,id,VHmode,amp)
  implicit none
  character(len=2), intent(in) :: VHmode
  real(8), intent(in) :: Mom(1:4,1:9)
  real(8), intent(in) :: mass(3:5,1:2)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: amp
!      integer, intent(in) :: id(9)
  integer  id(9)

  complex(8) PROP3, PROP4, PROP5
  complex(8) gFFZ, gZAFF
  complex(8) gVVS1(3), gVVS2(3) !1 = Z > ZH, 3 = Z > AH, 2 is not used
  complex(8) ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
  real(8) q3_q4,q3_q3,q4_q4,q5_q5

  complex(8) TriffbHa1, TriffbHa2, TriAHa1, TriAHa2, TriZZ, TriZA
  complex(8) BoxffbH, BoxffbH5, BoxZ
  complex(8) Zffb, Affb
  real(8) :: sprod(4,4)
  complex(8) :: Spaa(4,4), Spbb(4,4)

! SM couplings for fermion currents
  gFFZ = ci*dsqrt(couplZffsq) ! = i * sqrt[ gwsq/4d0/(1.0_dp-xw) ]
  gZAFF = ci*couplAZff ! = (loop > Z) * (A > ff) (minus from A > ff included)

  q3_q4 = -scr(Mom(:,3),Mom(:,4))
  q3_q3 = scr(Mom(:,3),Mom(:,3))
  q4_q4 = scr(Mom(:,4),Mom(:,4))
  q5_q5 = scr(Mom(:,5),Mom(:,5))

!HVV vertex
  gVVS1=0d0
  gVVS2=0d0

  ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q5_q5)
  ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q5_q5)
  ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q5_q5)
  ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q5_q5)
  gVVS1(1) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
  gVVS2(1) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )

  ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q5_q5)
  ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q5_q5)
  ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q5_q5)
  ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q5_q5)
  gVVS1(3) = ghz1_dyn*(mass(3,1)**2) + q3_q4 * ( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )
  gVVS2(3) = -( 2d0*ghz2_dyn + ghz3_dyn*q3_q4/Lambda**2 )

  gVVS1=ci*gVVS1/vev
  gVVS2=ci*gVVS2/vev
!gVVS1=0d0
!gVVS2=0d0
!gVVS2=0d0
!gVVS1(1)=ci*1d0*m_Z**2/vev !a1, for debugging use
!gVVS2(1)=ci*1d0/vev !a2, for debugging use
!gVVS1=ci*gVVS1 * dsqrt(4d0*pi*alpha_QED) /sitW/2d0/M_W !i * gVVS1 / vev
!gVVS2=ci*gVVS2 * dsqrt(4d0*pi*alpha_QED) /sitW/2d0/M_W !i * gVVS1 / vev

  PROP3 = PROPAGATOR(q3_q3,mass(3,1),mass(3,2))

  TriffbHa1=0d0
  TriffbHa2=0d0
  TriAHa1=0d0
  TriAHa2=0d0
  BoxffbH=0d0
  BoxffbH5=0d0
! amplitudes w/o HVV or VFF couplings to calculate
  if(id(6).ne.convertLHE(Pho_))then !Z/A decays
    call spinoru2(4,(/Mom(1:4,1),Mom(1:4,2),Mom(1:4,6),Mom(1:4,7)/),Spaa,Spbb,sprod) !!!correct one
    !call spinoru2(4,(/-Mom(1:4,1),-Mom(1:4,2),Mom(1:4,7),Mom(1:4,6)/),Spaa,Spbb,sprod) !!!test only
    if(VHmode.eq."tr".or.VHmode.eq."gg".or.VHmode.eq."nl")then
      if(gVVS1(1).ne.0d0.or.gVVS1(3).ne.0d0)then
        call ggTriffbHa1(Spaa,Spbb,sprod,helicity,TriffbHa1)
        TriffbHa1 = TriffbHa1*PROP3
      endif
      if(gVVS2(1).ne.0d0.or.gVVS2(3).ne.0d0)then
        call ggTriffbHa2(Spaa,Spbb,sprod,helicity,TriffbHa2)
        TriffbHa2 = TriffbHa2*PROP3

      endif
    endif
    if(VHmode.eq."bo".or.VHmode.eq."gg".or.VHmode.eq."nl")then
      if(kappa.ne.(0d0,0d0))then
        call ggBoxffbH(Mom,Spaa,Spbb,sprod,helicity,BoxffbH)!SM
      endif
      if(kappa_tilde.ne.(0d0,0d0))then
        call ggBoxffbH5(Mom,Spaa,Spbb,sprod,helicity,BoxffbH5)!gamma5 in Hff coupling
      endif
    endif

  else ! A is final state
    call spinoru2(4,(/Mom(1:4,1),Mom(1:4,2),-Mom(1:4,1)-Mom(1:4,2),Mom(1:4,4)/),Spaa,Spbb,sprod)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!???????
    if(VHmode.eq."tr".or.VHmode.eq."gg".or.VHmode.eq."nl")then
      if(gVVS1(3).ne.0d0)then
        call ggTriAHa1(Spaa,Spbb,sprod,helicity,TriAHa1)
        TriAHa1 = TriAHa1*PROP3
      endif
      if(gVVS2(3).ne.0d0)then
        call ggTriAHa2(Spaa,Spbb,sprod,helicity,TriAHa2)
        TriAHa2 = TriAHa2*PROP3
      endif
    endif
    if(VHmode.eq."bo")then!gg > box couples to axial vector only, which photon is not.
      print *, "gg > box > photon + H is 0!"
    endif
  endif
!print*, "TriffbHa2 = ",TriffbHa2
  TriZZ=0d0
  TriZA=0d0
  BoxZ=0d0
! HVV vertex
  if(id(6).ne.convertLHE(Pho_))then !Z/A decays
    TriZZ = TriffbHa1*gVVS1(1) + TriffbHa2*gVVS2(1)
    TriZA = TriffbHa1*gVVS1(3) + TriffbHa2*gVVS2(3)
    if((BoxffbH.ne.0d0).or.(BoxffbH5.ne.0d0))then
      BoxZ = (BoxffbH*kappa + BoxffbH5*ci*kappa_tilde) * (-ci) * m_top / vev
    endif
  else !A is final state
    amp = TriAHa1*gVVS1(3) + TriAHa2*gVVS2(3)
  endif
!print*, "TriffbHa1*gVVS1(1) = ",TriffbHa1*gVVS1(1)
!print*, "TriffbHa2*gVVS2(1) = ",TriffbHa2*gVVS2(1)
!print*, "TriZZ = ",TriZZ
  Zffb=0d0
  Affb=0d0
! triangle/box > Z couplings, and Z/A > f f~ couplings. i.e. all the EW couplings w/o H decay.
  if(id(6).ne.convertLHE(Pho_))then !Z/A decays
    if((gVVS1(1).ne.0d0).or.(gVVS2(1).ne.0d0).or.(BoxZ.ne.0d0))then !Z>ff~ couplings
      PROP4 = PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
      Zffb = (TriZZ+BoxZ) * ZFF(id(6),id(7),helicity(6),helicity(7))
      Zffb = Zffb * gFFZ**2 / 4d0 ! 1/4 is for axial couplings of t and b to Z*, which is (1/2)*(1/2)
      Zffb = Zffb * PROP4
    endif
    if(gVVS1(3).ne.0d0.or.gVVS2(3).ne.0d0)then !A>ff~ couplings
      PROP4 = PROPAGATOR(q4_q4,0d0,0d0)
      Affb = TriZA * AFF(id(6),id(7),helicity(6),helicity(7))
      Affb = Affb * gZAFF / 4d0 ! 1/4 is for axial couplings of t and b to Z, which is (1/2)*(1/2)
      Affb = Affb * PROP4
    endif
    amp = Affb + Zffb
  endif

! gg > triangle/box couplings
  amp = amp * (-1d0) * gs**2 !-1 = i^2 from g_s each

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
END subroutine amp_VH_gg

! gg > ZH
subroutine ggTriffbHa1(Spaa,Spbb,sprod,helicity,Tri)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  complex(8) C0, qli3
  complex(8) Tri

  Tri=0d0

  if(helicity(1)*helicity(2).lt.0d0)then
    Tri=0d0
  elseif(helicity(1).gt.0d0)then
    if(helicity(6).gt.0d0)then
      Tri=-Spbb(1,2)*(Spaa(1,4)*Spbb(1,3) + Spaa(2,4)*Spbb(2,3))/Spaa(1,2) !+++
    else
      Tri=-Spbb(1,2)*(Spaa(1,3)*Spbb(1,4) + Spaa(2,3)*Spbb(2,4))/Spaa(1,2) !++-
    endif
  else
    if(helicity(6).gt.0d0)then
      Tri= Spaa(1,2)*(Spaa(1,4)*Spbb(1,3) + Spaa(2,4)*Spbb(2,3))/Spbb(1,2) !--+
    else
      Tri= Spaa(1,2)*(Spaa(1,3)*Spbb(1,4) + Spaa(2,3)*Spbb(2,4))/Spbb(1,2) !---
    endif
  endif
      
  !C0 = qlI3(0d0,0d0,sprod(1,2),M_Top**2,M_Top**2,M_Top**2,Mu_Ren**2,0)
  call C0_cll(C0,dcmplx(0d0),dcmplx(sprod(1,2)),dcmplx(0d0),dcmplx(M_Top**2),dcmplx(M_Top**2),dcmplx(M_Top**2))

  C0 = C0 * ci * pisq / (2d0*pi)**4

  Tri = Tri * (sprod(1,2)/M_Z**2 - 1d0)

  Tri = ci * 8d0 * Tri * C0 * m_top**2 ! i = -1 (fermion loop) i^3 (3 fermion propagators)

  return
END subroutine ggTriffbHa1


subroutine ggTriffbHa2(Spaa,Spbb,sprod,helicity,Tri)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  complex(8) C0!, qli3, B0_00, B0_0m, B0_s0, B0_sm
  complex(8) Tri

  Tri=0d0

  !C0 = qlI3(0d0,0d0,sprod(1,2),M_Top**2,M_Top**2,M_Top**2,Mu_Ren**2,0)
  call C0_cll(C0,dcmplx(0d0),dcmplx(sprod(1,2)),dcmplx(0d0),dcmplx(M_Top**2),dcmplx(M_Top**2),dcmplx(M_Top**2))
  C0 = C0 * ci * pisq / (2d0*pi)**4

  if(helicity(1)*helicity(2).lt.0d0)then
    Tri=0d0
  elseif(helicity(1).gt.0d0)then
    if(helicity(6).gt.0d0)then
      Tri=(Spbb(1,2)*(m_Z**2 + Spaa(1,2)*Spbb(1,2))* &
    (Spaa(1,4)*Spbb(1,3) + Spaa(2,4)*Spbb(2,3))* &
    (Spaa(1,3)*Spbb(1,3) + Spaa(1,4)*Spbb(1,4) + Spaa(2,3)*Spbb(2,3) +  &
      Spaa(2,4)*Spbb(2,4)))/(m_Z**2*Spaa(1,2)) !+++
    else
      Tri=(Spbb(1,2)*(m_Z**2 + Spaa(1,2)*Spbb(1,2))* &
    (Spaa(1,3)*Spbb(1,4) + Spaa(2,3)*Spbb(2,4))* &
    (Spaa(1,3)*Spbb(1,3) + Spaa(1,4)*Spbb(1,4) + Spaa(2,3)*Spbb(2,3) +  &
      Spaa(2,4)*Spbb(2,4)))/(m_Z**2*Spaa(1,2)) !++-
    endif
  else
    if(helicity(6).gt.0d0)then
      Tri=(-Spaa(1,2)*(m_Z**2 + Spaa(1,2)*Spbb(1,2))* &
    (Spaa(1,4)*Spbb(1,3) + Spaa(2,4)*Spbb(2,3))* &
    (Spaa(1,3)*Spbb(1,3) + Spaa(1,4)*Spbb(1,4) + Spaa(2,3)*Spbb(2,3) +  &
      Spaa(2,4)*Spbb(2,4)))/(m_Z**2*Spbb(1,2)) !--+
    else
      Tri=(-Spaa(1,2)*(m_Z**2 + Spaa(1,2)*Spbb(1,2))* &
    (Spaa(1,3)*Spbb(1,4) + Spaa(2,3)*Spbb(2,4))* &
    (Spaa(1,3)*Spbb(1,3) + Spaa(1,4)*Spbb(1,4) + Spaa(2,3)*Spbb(2,3) +  &
      Spaa(2,4)*Spbb(2,4)))/(m_Z**2*Spbb(1,2)) !---
    endif
  endif

  Tri = ci * 4d0 * Tri * C0 * m_top**2! i = -1 (fermion loop) i^3 (3 fermion propagators)

  return
END subroutine ggTriffbHa2


subroutine ggTriAHa1(Spaa,Spbb,sprod,helicity,Tri)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  complex(8) C0, qli3
  complex(8) Tri

  Tri=0d0

  if(helicity(1)*helicity(2).lt.0d0)then
    Tri=0d0
  elseif(helicity(1).gt.0d0)then
    if(helicity(6).gt.0d0)then
      Tri=0d0 !+++
    else
      Tri=0d0 !++-
    endif
  else
    if(helicity(6).gt.0d0)then
      Tri=0d0 !--+
    else
      Tri=0d0 !---
    endif
  endif

  call C0_cll(C0,dcmplx(0d0),dcmplx(-sprod(1,2)),dcmplx(0d0),dcmplx(M_Top**2),dcmplx(M_Top**2),dcmplx(M_Top**2))

  C0 = C0 / ci / pisq

  Tri = 0d0

  return
END subroutine ggTriAHa1

subroutine ggTriAHa2(Spaa,Spbb,sprod,helicity,Tri)
  implicit none
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  complex(8) C0, qli3
  complex(8) Tri

  Tri=0d0

  if(helicity(1)*helicity(2).lt.0d0)then
    Tri=0d0
  elseif(helicity(1).gt.0d0)then
    if(helicity(6).gt.0d0)then
      Tri=0d0 !+++
    else
      Tri=0d0 !++-
    endif
  else
    if(helicity(6).gt.0d0)then
      Tri=0d0 !--+
    else
      Tri=0d0 !---
    endif
  endif
      
  call C0_cll(C0,dcmplx(0d0),dcmplx(-sprod(1,2)),dcmplx(0d0),dcmplx(M_Top**2),dcmplx(M_Top**2),dcmplx(M_Top**2))

  C0 = C0 / ci / pisq

  Tri = 0d0

  return
END subroutine ggTriAHa2

subroutine ggBoxffbH(Mom,Spaa,Spbb,sprod,helicity,Box)
  implicit none
  real(8), intent(in) :: Mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  complex(8) Box, Box7_8, Box9_10, Box11_12

  Box = 0d0

  if(helicity(1).gt.0d0.and.helicity(2).lt.0d0)then
    if(helicity(6).gt.0d0)then !+-+
      call  ffbHBox7_8pmp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10pmp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12pmp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !+--
      call  ffbHBox7_8pmm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10pmm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12pmm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  elseif(helicity(1).lt.0d0.and.helicity(2).gt.0d0)then
    if(helicity(6).gt.0d0)then !-++
      call  ffbHBox7_8mpp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10mpp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12mpp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !-+-
      call  ffbHBox7_8mpm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10mpm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12mpm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  elseif(helicity(1).gt.0d0.and.helicity(2).gt.0d0)then
    if(helicity(6).gt.0d0)then !+++
      call  ffbHBox7_8ppp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10ppp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12ppp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !++-
      call  ffbHBox7_8ppm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10ppm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12ppm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  else
    if(helicity(6).gt.0d0)then !--+
      call  ffbHBox7_8mmp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10mmp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12mmp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !---
      call  ffbHBox7_8mmm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbHBox9_10mmm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbHBox11_12mmm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  endif

  Box = Box7_8 + Box9_10 + Box11_12
  Box = -Box * ci * pisq / (2d0*pi)**4 ! -1 = -1 (fermion loop) i^4 (4 fermion propagators)
  Box = -Box!reason not identified but it makes the box-triangle interference correct...

  return
  END subroutine ggBoxffbH



subroutine ggBoxffbH5(Mom,Spaa,Spbb,sprod,helicity,Box)
  implicit none
  real(8), intent(in) :: Mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  complex(8) Box, Box7_8, Box9_10, Box11_12

  Box = 0d0

  if(helicity(1).gt.0d0.and.helicity(2).lt.0d0)then
    if(helicity(6).gt.0d0)then !+-+
      call  ffbH5Box7_8pmp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10pmp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12pmp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !+--
      call  ffbH5Box7_8pmm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10pmm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12pmm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  elseif(helicity(1).lt.0d0.and.helicity(2).gt.0d0)then
    if(helicity(6).gt.0d0)then !-++
      call  ffbH5Box7_8mpp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10mpp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12mpp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !-+-
      call  ffbH5Box7_8mpm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10mpm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12mpm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  elseif(helicity(1).gt.0d0.and.helicity(2).gt.0d0)then
    if(helicity(6).gt.0d0)then !+++
      call  ffbH5Box7_8ppp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10ppp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12ppp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !++-
      call  ffbH5Box7_8ppm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10ppm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12ppm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  else
    if(helicity(6).gt.0d0)then !--+
      call  ffbH5Box7_8mmp(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10mmp(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12mmp(Mom,Spaa,Spbb,sprod,Box11_12)
    else !---
      call  ffbH5Box7_8mmm(Mom,Spaa,Spbb,sprod,Box7_8)
      call  ffbH5Box9_10mmm(Mom,Spaa,Spbb,sprod,Box9_10)
      call ffbH5Box11_12mmm(Mom,Spaa,Spbb,sprod,Box11_12)
    endif
  endif

  Box = Box7_8 + Box9_10 + Box11_12
  Box = -Box * ci * pisq / (2d0*pi)**4 ! -1 = -1 (fermion loop) i^4 (4 fermion propagators)
  Box = -Box!reason not identified but it makes the box-triangle interference correct...

  return
  END subroutine ggBoxffbH5





!#endif







end module ModVHgg
!!--YaofuZhou-----------------------------------------
