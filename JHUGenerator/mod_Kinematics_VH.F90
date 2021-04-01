MODULE ModKinematics_VH
implicit none
save


contains








SUBROUTINE WriteOutEvent_VH(id,helicity,MomExt,EventWeight)
use ModKinematics
use ModParameters
use ModMisc
implicit none
integer, intent(in) :: id(9)
double precision, intent(in) :: helicity(9)
double precision, intent(in) :: MomExt(1:4,1:9)
real(8), intent(in) :: EventWeight
real(8) :: MomDummy(1:4,1:9), MassDummy(9)
integer :: ICOLUP(4,2)
real(8) :: Spin
integer :: i
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1F3.0)"



MomDummy = MomExt/GeV
do i=1,9
  MassDummy(i) = Get_MInv(MomDummy(:,i))
enddo

IDPRUP=Process
SCALUP=Mu_Fact/GeV
AQEDUP=alpha_QED
AQCDUP=alphas

call getHiggsDecayLength(HiggsDKLength)

NUP = 4
if(HbbDecays) NUP=NUP+2
if(.not.IsAPhoton(DecayMode1)) NUP=NUP+2

    XWGTUP=EventWeight

write(io_LHEOutFile,"(A)") "<event>"
if( .not. ReadLHEFile ) write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

Spin = 0.1d0

if(COLLIDER.eq.0)then
  ICOLUP=0
else
  if(id(1).ne.convertLHE(Glu_).and.id(2).ne.convertLHE(Glu_))then
    if(id(1).gt.0)then
      ICOLUP(1,1)=503
      ICOLUP(1,2)=0
      ICOLUP(2,1)=0
      ICOLUP(2,2)=503
    else
      ICOLUP(1,1)=0
      ICOLUP(1,2)=503
      ICOLUP(2,1)=503
      ICOLUP(2,2)=0
    endif
  elseif(id(1).eq.convertLHE(Glu_).and.id(2).eq.convertLHE(Glu_))then
    ICOLUP(1,1)=503
    ICOLUP(1,2)=504
    ICOLUP(2,1)=504
    ICOLUP(2,2)=503
  else!gq cases
    ICOLUP(1:2,:)=0!to be configured when gq is developed
  endif
endif

if((id(6).eq.convertLHE(Up_)).or.(id(6).eq.convertLHE(Dn_)).or.(id(6).eq.convertLHE(Str_)).or.(id(6).eq.convertLHE(Chm_)).or.(id(6).eq.convertLHE(Bot_)))then
    ICOLUP(3,1)=502
    ICOLUP(3,2)=0
    ICOLUP(4,1)=0
    ICOLUP(4,2)=502
elseif((id(6).eq.convertLHE(AUp_)).or.(id(6).eq.convertLHE(ADn_)).or.(id(6).eq.convertLHE(AStr_)).or.(id(6).eq.convertLHE(AChm_)).or.(id(6).eq.convertLHE(ABot_)))then
    ICOLUP(3,1)=0
    ICOLUP(3,2)=502
    ICOLUP(4,1)=502
    ICOLUP(4,2)=0
else
    ICOLUP(3:4,1:2)=0
endif

do i=1,2
    write(io_LHEOutFile,fmt1) id(i), -1,0,0,ICOLUP(i,1),ICOLUP(i,2),MomDummy(2:4,i), MomDummy(1,i), 0.0d0, 0.0d0, Spin
enddo

if(IsAPhoton(DecayMode1)) then
  write(io_LHEOutFile,fmt1) id(4), 1,1,2,0,0,MomDummy(2:4,4), MomDummy(1,4), MassDummy(4), 0d0, Spin
else
  write(io_LHEOutFile,fmt1) id(4), 2,1,2,0,0,MomDummy(2:4,4), MomDummy(1,4), MassDummy(4), 0d0, Spin
endif

if(HbbDecays.eqv..true.)then
  write(io_LHEOutFile,fmt1) id(5), 2,1,2,0,0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength, Spin
else
  write(io_LHEOutFile,fmt1) id(5), 1,1,2,0,0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength, Spin
endif

if(.not.IsAPhoton(DecayMode1)) then
  write(io_LHEOutFile,fmt1) id(6), 1,3,3,ICOLUP(3,1),ICOLUP(3,2),MomDummy(2:4,6), MomDummy(1,6), MassDummy(6), 0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(7), 1,3,3,ICOLUP(4,1),ICOLUP(4,2),MomDummy(2:4,7), MomDummy(1,7), MassDummy(7), 0.0d0, Spin
endif

if(HbbDecays.eqv..true.)then
  write(io_LHEOutFile,fmt1) id(8), 1,4,4,501,0,MomDummy(2:4,8), MomDummy(1,8), MassDummy(8), 0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(9), 1,4,4,0,501,MomDummy(2:4,9), MomDummy(1,9), MassDummy(9), 0.0d0, Spin
endif

write(io_LHEOutFile,"(A)") "</event>"


END SUBROUTINE








!
!
!SUBROUTINE EvalPhaseSpace_VH(yRnd,MomExt,inv_mass,mass,PSWgt,HDecays,PhoOnshell,ZAinterference)
!use ModParameters
!implicit none
!
!      logical, intent(in) :: HDecays
!      logical, intent(in), optional :: PhoOnshell
!      logical, optional :: ZAinterference
!      real(8), intent(in) :: yRnd(1:20),mass(9,2)
!      real(8) :: phi, beta, gamma
!      real(8) :: temp_vector(4), temp_boost(4)
!      integer :: i
!      real(8), parameter :: Twopi = 2d0 * Pi
!      real(8) :: MomExt(1:4,1:9)
!      real(8) :: MomDummy(1:4)
!      !double precision four_momentum(7,4)
!      real(8), intent(out) :: PSWgt,inv_mass(9)
!! 1=E, 2,3,4=p_x,y,z
!!psg_mass(mass, width)
!      real(8) :: cm_abs3p(9)
!      real(8) :: cm_sin_theta(9), cm_cos_theta(9)
!      real(8) :: cm_sin_phi(9), cm_cos_phi(9)
!!use Cauchy distribution for Breit-Wigner distribution for the invariant mass of 2 and 3?
!!      logical, parameter :: breit_wigner = .true.
!      real(8) :: jacobian4, jacobian5
!      logical :: hasAonshell, hasInterference
!
!      hasAonshell = .false.
!      if(present(PhoOnshell)) then
!         hasAonshell=PhoOnshell
!      endif
!      hasInterference = .false.
!      if(present(ZAinterference)) then
!         hasInterference=ZAinterference
!      endif
!
!      MomExt(:,4:9)=0d0
!      inv_mass(4:9)=0d0
!
!!333333333
!      inv_mass(3) = dsqrt((MomExt(1,3)+MomExt(4,3))*(MomExt(1,3)-MomExt(4,3)))
!
!!555555555
!      if(HDecays)then
!        inv_mass(5) = dsqrt(dabs(bw_sq(yRnd(13),mass(5,1), mass(5,2), inv_mass(3)**2, jacobian5)))
!        jacobian5 = jacobian5 /16d0 /Pi**2 != (ds5/2pi)*(1/8pi)
!      else
!        inv_mass(5) = mass(5,1)
!        jacobian5=1d0
!      endif
!
!!4444444444
!      if(hasAonshell) then
!        inv_mass(4) = mass(4,1)
!        jacobian4=1d0
!      elseif(.not.hasInterference)then
!        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
!        jacobian4 = jacobian4 /16d0 /Pi**2 != (ds4/2pi)*(1/8pi)
!      else
!        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
!        !inv_mass(4) = dsqrt(one_over_s_sq(yRnd(12), MPhotonCutoff**2, (inv_mass(3)-inv_mass(5))**2, jacobian4))
!        jacobian4 = jacobian4 /16d0 /Pi**2
!      endif
!
!!444444444444
!!energy of 4 in the CM frame of 3
!      MomExt(1,4)=(inv_mass(3)**2+(inv_mass(4)+inv_mass(5))*(inv_mass(4)-inv_mass(5)))/2d0/inv_mass(3)
!!|3-momentum| of 4 in the CM frame of 3
!      cm_abs3p(4) = dsqrt((MomExt(1,4)+inv_mass(4)) * (MomExt(1,4)-inv_mass(4)))
!!generating cos(theta_4) and phi_4 in the CM frame of 3
!      cm_cos_theta(4) = yRnd(7)
!      cm_cos_theta(4) = cm_cos_theta(4)*2d0-1d0
!      cm_sin_theta(4) = dsqrt((1d0+cm_cos_theta(4))  *(1d0-cm_cos_theta(4)))
!      phi = yRnd(6)
!      phi=Twopi*phi
!      cm_cos_phi(4) = dcos(phi)
!      cm_sin_phi(4) = dsin(phi)
!!3-momentum of 4 in the CM frame of 3
!      MomExt(2,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_cos_phi(4)
!      MomExt(3,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_sin_phi(4)
!      MomExt(4,4)=cm_abs3p(4)*cm_cos_theta(4)
!!boost the 4-momentum of 4 to the lab frame
!!x and y components stay the same
!!z component
!      beta = MomExt(4,3)/MomExt(1,3)
!      gamma = 1d0/dsqrt((1d0+beta)*(1d0-beta))
!      MomDummy(1:4)=MomExt(1:4,4)
!      MomExt(4,4)=(MomDummy(4)+MomDummy(1)*beta) *gamma
!!energy
!      MomExt(1,4)=(MomDummy(1)+MomDummy(4)*beta) *gamma
!
!!555555555555555555
!      MomExt(1:4,5) = MomExt(1:4,3) - MomExt(1:4,4)
!
!!666666666666666666
!      if(.not.hasAonshell) then
!!invariant mass of 6
!         inv_mass(6)=0d0
!!energy of 6 in the CM frame of 4
!         MomExt(1,6)=inv_mass(4)/2d0
!!|3-momentum| of 6 in the CM frame of 4
!         cm_abs3p(6)=MomExt(1,6)
!!generating cos(theta_6) and phi_6 in the CM frame of 4
!!z-axis is along the boost of 2
!         cm_cos_theta(6) = yRnd(9)
!         cm_cos_theta(6) = cm_cos_theta(6)*2d0-1d0
!         cm_sin_theta(6) = dsqrt((1d0+cm_cos_theta(6)) *(1d0-cm_cos_theta(6)))
!         cm_cos_phi(6) = yRnd(8)
!         phi=Twopi*cm_cos_phi(6)
!         cm_cos_phi(6) = dcos(phi)
!         cm_sin_phi(6) = dsin(phi)
!!3-momentum of 6 in the CM frame of 4
!         MomExt(2,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_cos_phi(6)
!         MomExt(3,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_sin_phi(6)
!         MomExt(4,6)=cm_abs3p(6)*cm_cos_theta(6)
!!boost the 4-momentum of 6 to the lab frame
!         temp_vector = MomExt(1:4,6)
!         temp_boost = MomExt(1:4,4)
!         call LORENTZ(temp_vector, temp_boost)
!         MomExt(1:4,6) = temp_vector
!      else
!         MomExt(1:4,6)=MomExt(1:4,4)
!      endif
!
!!7777777777777777777777
!!invariant mass of 7
!      inv_mass(7)=0d0
!!4-momentum of 7 (lab frame) by energy-momentum conservation
!      MomExt(1:4,7)=MomExt(1:4,4)-MomExt(1:4,6)
!
!!8888888888888888888888
!      if(HDecays)then
!!invariant mass of 8
!        inv_mass(8)=0d0
!!energy of 8 in the CM frame of 5
!        MomExt(1,8)=inv_mass(5)/2d0
!!|3-momentum| of 8 in the CM frame of 5
!        cm_abs3p(8)=MomExt(1,8)
!!generating cos(theta_8) and phi_8 in the CM frame of 5
!!z-axis is along the boost of 5
!        cm_cos_theta(8) = yRnd(11)
!        cm_cos_theta(8) = cm_cos_theta(8)*2d0-1d0
!        cm_sin_theta(8) = dsqrt((1d0+cm_cos_theta(8)) *(1d0-cm_cos_theta(8)))
!        cm_cos_phi(8) = yRnd(10)
!        phi=Twopi*cm_cos_phi(8)
!        cm_cos_phi(8) = dcos(phi)
!        cm_sin_phi(8) = dsin(phi)
!!3-momentum of 8 in the CM frame of 5
!!x and y components are not necessary, yet
!        MomExt(2,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_cos_phi(8)
!        MomExt(3,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_sin_phi(8)
!        MomExt(4,8)=cm_abs3p(8)*cm_cos_theta(8)
!!boost the 4-momentum of 6 to the lab frame
!        temp_vector = MomExt(1:4,8)
!        temp_boost = MomExt(1:4,5)
!        call LORENTZ(temp_vector, temp_boost)
!        MomExt(1:4,8) = temp_vector
!!9999999999999999999999
!!4-momentum of 9 (lab frame) by energy-momentum conservation
!        MomExt(1:4,9)=MomExt(1:4,5)-MomExt(1:4,8)
!      endif
!
!      PSWgt = jacobian4*jacobian5 * cm_abs3p(4)/(4d0*pi)/inv_mass(3)
!      !print *,  "()",inv_mass, jacobian4, jacobian5, cm_abs3p(4), PSWgt, "()"
!
!!do i=4,7
!!print *, dsqrt(dabs(four_momentum(i,:).dot.four_momentum(i,:)))
!!enddo
!!pause
!
!RETURN
!END SUBROUTINE
!
!
!
!
!
!
!
!! Breit-Wigner mass^2
!function bw_sq(x, m, ga, smax, jacobian)
!implicit none
!real(8) :: bw_sq
!real(8), intent(in) :: m, ga, smax,x
!real(8) :: xmin, xmax,xprime
!real(8), intent(out) :: jacobian
!
!xmin=-datan(m/ga)/m/ga
!xmax=-datan((-smax+m**2)/ga/m)/ga/m
!xprime=x*(xmax-xmin)+xmin
!bw_sq=m**2+dtan(xprime*ga*m)*ga*m
!jacobian=(ga*m)**2 * (1d0+dtan(ga*m*xprime)**2) * (xmax-xmin)
!
!return
!end function bw_sq
!
!!generate s according to 1/s^2 distribution
!function one_over_s_sq(x, smin, smax, jacobian)
!implicit none
!real(8) :: one_over_s_sq
!real(8), intent(in) :: x, smin, smax
!real(8), intent(out) :: jacobian
!
!one_over_s_sq=smin/(1d0-x*(smax-smin)/smax)
!jacobian=smin*smax*(smax-smin)/((x-1d0)*smax-x*smin)**2
!
!return
!end function one_over_s_sq
!
!
!!LORENTZ.F
!!VERSION 20130123
!!
!!A subroutine that performs a general boost to a four vector
!!(vector) based on another four vector (boost). The primed and
!!unprimed frames have their axes in parallel to one another.
!!Rotation is not performed by this subroutine.
!
!      subroutine LORENTZ(vector, boost)
!
!      implicit none
!
!      double precision vector(4), boost(4)
!      double precision lambdaMtrx(4,4), vector_copy(4)
!      double precision beta(2:4), beta_sq, gamma
!      integer i,j
!      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
!
!!      double precision KRONECKER_DELTA
!!      external KRONECKER_DELTA
!
!      do i=2,4
!        beta(i) = boost(i)/boost(1)
!      enddo
!
!      beta_sq = beta(2)**2+beta(3)**2+beta(4)**2
!
!  if(beta_sq.ge.epsilon)then
!
!      gamma = 1d0/dsqrt(1d0-beta_sq)
!
!      lambdaMtrx(1,1) = gamma
!
!      do i=2,4
!        lambdaMtrx(1,i) = gamma*beta(i)
!        lambdaMtrx(i,1) = lambdaMtrx(1,i)
!      enddo
!
!      do i=2,4
!      do j=2,4
!        lambdaMtrx(i,j) = (gamma-1d0)*beta(i)*beta(j)/beta_sq + KRONECKER_DELTA(i,j)
!      enddo
!      enddo
!
!
!!apply boost to vector1
!      vector_copy = vector
!      vector = 0d0
!      do i=1,4
!      do j=1,4
!        vector(i) = vector(i) + lambdaMtrx(i,j)*vector_copy(j)
!      enddo
!      enddo
!  endif
!
!      return
!      END subroutine LORENTZ
!
!
!
!! KRONECKER_DELTA.F
!!
!! KRONECKER_DELTA(i,j)
!! A function that returns 1 if i=j, and 0 otherwise.
!      double precision function KRONECKER_DELTA(i,j)
!      integer i,j
!      if(i.eq.j)then
!        KRONECKER_DELTA = 1d0
!      else
!        KRONECKER_DELTA = 0d0
!      endif
!
!      return
!      end function KRONECKER_DELTA



!moved to mod_Parameters.F90
!SUBROUTINE getHiggsDecayLength(ctau)
!use ModParameters
!implicit none
!real(8) :: x,xp,xpp,Len0,propa,ctau,ctau0
!integer :: loop
!
!
!     ctau  = 0d0
!     ctau0 = HiggsDecayLengthMM
!     if( ctau0.lt.1d-16 ) RETURN
!
!     do loop=1,4000000!  4Mio. tries otherwise return zero
!          call random_number(x)
!          xp = 10*x*ctau0             ! scan between 0..10*ctau0
!
!          propa = dexp( -xp/(ctau0) ) ! the max of propa is 1.0
!          call random_number(xpp)
!          if( xpp.lt.propa ) then!   accept
!                ctau = xp
!                RETURN
!          endif
!     enddo
!
!
!RETURN
!END SUBROUTINE



SUBROUTINE EvalPhasespace_VH(xRnd,Energy,Mom,id,Jac,HDecays,PhoOnshell)
  use ModParameters
  use ModPhasespace
  use ModMisc
  implicit none
  real(8) :: Mom(1:4,1:9)
!  real(8) :: xRnd(6:13),Energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8), intent(in) :: xRnd(6:13),Energy
  integer, intent(in) :: id(6:9)
  real(8), intent(out) :: Jac
  real(8) :: Jac1,Jac2,Jac3,Jac4,Jac5
  real(8) :: s67,s89
  logical, intent(in), optional :: PhoOnshell
  logical, intent(in) :: HDecays
  logical :: hasAonshell

  hasAonshell = .false.
  if(present(PhoOnshell)) then
    hasAonshell=PhoOnshell
  endif

  Mom(:,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/) ! beam 1
  Mom(:,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/) ! beam 2
  Mom(:,3) = Energy * (/1d0,0d0,0d0,0d0/) ! V*

! masses
! 555555555 Higgs
  if(HDecays)then
    Jac3 = s_channel_propagator(M_Reso**2,Ga_Reso,(getMass(convertLHEreverse(id(8)))+getMass(convertLHEreverse(id(9))))**2,Energy**2,xRnd(13),s89) ! H --> 89
  else
    s89 = M_Reso**2 ! if H onshell
    Jac3 = 1d0
  endif

! 444444444 V on-shell
  if(hasAonshell) then
    s67 = 0d0 ! if photon in the final state
    Jac2 = 1d0
  else ! if V decays
    Jac2 = s_channel_propagator(M_V**2,Ga_V,(getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))**2,(Energy-dsqrt(s89))**2,xRnd(12),s67) ! V --> 67
  endif

! splittings
  Jac1 = s_channel_decay(Mom(:,3),s67,s89,xRnd(6:7),Mom(:,4),Mom(:,5)) ! V* --> VH

! 444444444 V on-shell
  if(hasAonshell) then
    Mom(:,6) = Mom(:,4)!for JHUGen, if associated photo does not decay, its "decay product 6" carrys the same info as the photon.
    Mom(:,7) = 0d0
    Jac4 = 1d0
  else ! if V decays
    Jac4 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(8:9),Mom(:,6),Mom(:,7)) !  V --> 67
  endif

! 555555555 Higgs
  if(HDecays)then
    Jac5 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(10:11),Mom(:,8),Mom(:,9)) !  H --> 89
  else
    Jac5 = 1d0
    Mom(:,8:9) = 0d0
  endif  

  Jac = Jac1*Jac2*Jac3*Jac4*Jac5

  if((.not.HbbDecays).and.hasAonshell)then   ! if 2 particle final state
    Jac = Jac * PSNorm2
  elseif(HbbDecays.and.(.not.hasAonshell))then ! if 4 particle final state
    Jac = Jac * PSNorm4
  else                                       ! if 3 particle final state
    Jac = Jac * PSNorm3
  endif


  if( isNan(jac).or.(Jac.le.0d0).or.(Jac5.le.0d0) .or.isNan(mom(1,8)) ) then
     Jac = 0d0
     Mom=0d0
  endif

  if(HDecays.and.((Jac5.le.0d0).or.isnan(Jac5).or.isNan(mom(1,8)).or.(mom(1,8).le.0d0)).or.isNan(mom(1,5)).or.(mom(1,5).le.0d0) )then
    Jac = 0d0
    Mom=0d0
  endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_VHglu(xRnd,Energy,Mom,id,Jac,HDecays,PhoOnshell)! ordering 1,2:in  3:Z*  4:Z   5:Higgs  67:l-l+   89:bb~  10:glu
  use ModParameters
  use ModPhasespace
  use ModMisc
  implicit none
  real(8) :: Mom(1:4,1:10)
  real(8), intent(in) :: xRnd(6:16),Energy
  integer, intent(in) :: id(6:9)
  real(8), intent(out) :: Jac
  real(8) :: Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7
  real(8) :: s6789,s67,s89
  logical, intent(in), optional :: PhoOnshell
  logical, intent(in) :: HDecays
  logical :: hasAonshell
  real(8) :: Momab(1:4)  

  hasAonshell = .false.
  if(present(PhoOnshell)) then
    hasAonshell=PhoOnshell
  endif

  Mom(:,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/) ! beam 1
  Mom(:,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/) ! beam 2
  Momab = Energy * (/1d0,0d0,0d0,0d0/) ! 678910

! masses
  if((.not.HbbDecays).and.hasAonshell)then   ! if H, V
    Jac1 = s_channel_propagator(m_V**2,0d0,(M_Reso+M_V)**2,Energy**2,xRnd(14),s6789)
  elseif(HbbDecays.and.(.not.hasAonshell))then ! if H > bb, V > ff
    Jac1 = s_channel_propagator(m_V**2,0d0,(M_Reso+M_V-5d0*Ga_V-5d0*Ga_Reso)**2,Energy**2,xRnd(14),s6789)
  elseif((.not.HbbDecays).and.(.not.hasAonshell))then  ! if H, V > ff
    Jac1 = s_channel_propagator(m_V*2,0d0,(M_Reso+M_V-5d0*Ga_V+getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))**2,Energy**2,xRnd(14),s6789)
  else                                                 ! if H > bb, V
    Jac1 = s_channel_propagator(m_V**2,0d0,(M_Reso+M_V-5d0*Ga_V-5d0*Ga_Reso+getMass(convertLHEreverse(id(8)))+getMass(convertLHEreverse(id(9))))**2,Energy**2,xRnd(14),s6789)
  endif

  if((.not.HbbDecays).and.hasAonshell)then   ! if H, V
    s67 = m_V**2
    s89 = M_Reso**2
    Jac2 = 1d0
    Jac3 = 1d0
  elseif(HbbDecays.and.(.not.hasAonshell))then ! if H > bb, V > ff
    Jac3 = s_channel_propagator(M_Reso**2,Ga_Reso,(getMass(convertLHEreverse(id(8)))+getMass(convertLHEreverse(id(9))))**2,s6789,xRnd(13),s89)
    Jac2 = s_channel_propagator(M_V**2,Ga_V,(M_V-5d0*Ga_V)**2,(dsqrt(s6789)-dsqrt(s89))**2,xRnd(12),s67)
  elseif((.not.HbbDecays).and.(.not.hasAonshell))then  ! if H, V > ff
    Jac2 = s_channel_propagator(M_V**2,Ga_V,(getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))**2,(dsqrt(s6789)-m_Reso)**2,xRnd(12),s67)
    s89 = M_Reso**2
    Jac3 = 1d0
  else                                                 ! if H > bb, V
    s67 = m_V**2
    Jac2 = 1d0
    Jac3 = s_channel_propagator(M_Reso**2,Ga_Reso,(M_Reso-5d0*Ga_Reso)**2,(dsqrt(s6789)-dsqrt(s67))**2,xRnd(13),s89)
  endif


!  splittings
  Jac4 = s_channel_decay(Momab,s6789,0d0,xRnd(15:16),Mom(:,3),Mom(:,10)) ! Z* + glu
  Jac5 = s_channel_decay(Mom(:,3),s67,s89,xRnd(6:7),Mom(:,4),Mom(:,5)) ! Z + H

  if(hasAonshell) then
    Jac6 = 1d0
    Mom(:,6) = Mom(:,4)!for JHUGen, if associated photo does not decay, its "decay product 6" carrys the same info as the photon.
  else ! if V decays
    Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(8:9),Mom(:,6),Mom(:,7)) !  Z --> 67
  endif

  if(HDecays)then
    Jac7 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(10:11),Mom(:,8),Mom(:,9)) !  H --> 89
  else
    Jac7 = 1d0
    Mom(:,8:9)=0d0
  endif

  Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6*Jac7

  if((.not.HbbDecays).and.hasAonshell)then   ! if 3 particle final state
    Jac = Jac * PSNorm3
  elseif(HbbDecays.and.(.not.hasAonshell))then ! if 5 particle final state
    Jac = Jac * PSNorm5
  else                                       ! if 4 particle final state
    Jac = Jac * PSNorm4
  endif


  if( isNan(jac).or.(jac.le.0d0)) then
    Jac = 0d0
    Mom=0d0
  endif

  if(HDecays.and.((jac7.le.0d0).or.isnan(jac7).or.isNan(mom(1,8)).or.(mom(1,8).le.0d0)).or.isNan(mom(1,5)).or.(mom(1,5).le.0d0) )then
    Jac = 0d0
    Mom=0d0
  endif


RETURN
END SUBROUTINE







!SUBROUTINE EvalPhasespace_VHglu(xRnd,Energy,Mom,id,Jac,HDecays,PhoOnshell)! ordering 1,2:in  3,4:Higgs  56:Z  7:glu
!use ModParameters
!use ModPhasespace
!use ModMisc
!implicit none
!real(8) :: xchannel
!real(8) :: s45,s345,Mom_DummyX(1:4),Mom_DummyZ(1:4),Mom_DummyZ2(1:4)
!real(8) :: SingDepth
!integer :: Pcol1,Pcol2,Steps
!real(8) :: Mom(1:4,1:10)
!real(8), intent(in) :: xRnd(6:16),Energy
!integer, intent(in) :: id(6:9)
!real(8), intent(out) :: Jac
!real(8) :: Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7
!real(8) :: s6789,s67,s89
!logical, intent(in), optional :: PhoOnshell
!logical, intent(in) :: HDecays
!logical :: hasAonshell
!real(8) :: Momab(1:4)  
!
!
!   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
!   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)
!
!!  stable Higgs,  decay not yet implemented but momenta 3+4 are allocated 
!   ! masses
!   Jac1 = s_channel_propagator(M_V**2,0d0,M_Reso**2,Energy**2,xRnd(14),s345)     ! Z*+H
!!print*,M_V**2,M_Reso**2,Energy**2,xRnd(14),s345
!   Jac2 = s_channel_propagator(M_V**2,Ga_V,(getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))**2,(dsqrt(s345)-m_Reso)**2,xRnd(12),s45)  ! Z-->56
!
!
!!  splittings
!   Mom_DummyX(1:4) = (/Energy,0d0,0d0,0d0/)   
!   Jac3 = s_channel_decay(Mom_DummyX(1:4),s345,0d0,xRnd(15:16),Mom(:,3),Mom(:,10)) ! Z* + glu
!   Jac4 = s_channel_decay(Mom(:,3),s45,M_Reso**2,xRnd(6:7),Mom(:,4),Mom(:,5)) ! Z + H
!   !Mom(1:4,4) = 0d0
!   Jac5 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(8:9),Mom(:,6),Mom(:,7)) !  Z --> 56
!   Mom(:,8)=0d0
!   Mom(:,9)=0d0
!   Jac = Jac1*Jac2*Jac3*Jac4*Jac5* PSNorm4
!
!
!   !this is for the soft/collinear limit checks only!!
!!         Pcol1= 2 -1
!!         Pcol2= 6 -1
!!         SingDepth = 1e-13
!!         Steps = 10
!!         call gensing(4,Energy,(/M_Reso,0d0,0d0,0d0/),Mom(1:4,4:7),Pcol1,Pcol2,SingDepth,Steps)
!!         Mom(1:4,3) = Mom(1:4,4)
!!         Mom(1:4,4) = 0d0
!!         print *, "generating singular phase space"
!!         Jac=1d0        
!   !end check
!
!!print*,dsqrt(s345),dsqrt(s45),m_reso
!!Print*,Jac1
!!print*,Mom(:,1)
!!print*,Mom(:,2)
!!print*,Mom(:,3)
!!print*,Mom(:,4)
!!print*,Mom(:,5)
!!print*,Mom(:,6)
!!print*,Mom(:,7)
!!print*,Mom(:,8)
!!print*,Mom(:,9)
!!print*,Mom(:,10)
!!pause
!
!!   if( isNan(jac) ) then
!!      print *, "EvalPhasespace_VHglu NaN"
!!      print *, Energy
!!      print *, s345,s45
!!      print *, Jac1,Jac2,Jac3,Jac4,Jac5
!!      print *, xRnd  
!!      Jac = 0d0
!!   endif   
!
!RETURN
!END SUBROUTINE









SUBROUTINE EvalPhasespace_VHglu_singular(xRnd,Energy,Mom,id,Jac,HDecays,PhoOnshell)
!  this is for the soft/collinear limit checks only!!
!  modified from Markus' email...
  use ModParameters
  use ModPhasespace
  use ModMisc
  implicit none
  real(8) :: Mom(1:4,1:10)
  real(8), intent(in) :: xRnd(6:16),Energy
  integer, intent(in) :: id(6:9)
  real(8), intent(out) :: Jac
  real(8) :: Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7
  real(8) :: s6789,s67,s89
  logical, intent(in), optional :: PhoOnshell
  logical, intent(in) :: HDecays
  logical :: hasAonshell
  real(8) :: Momab(1:4)

!for singular PS check
  real(8) :: SingDepth
  integer :: Pcol1,Pcol2,Steps
  real(8) :: MomC(1:4,1:4)
!for singular PS check

  print *, "generating singular phase space"
  print *, "Energy = ", Energy
  Pcol1= 6 -1  !6 6 is soft glu, 1 6 = glu // p1, 2 6 = glu // p2
  Pcol2= 6 -1
  SingDepth = 1e-13
  Steps = 30

  Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
  Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

  Mom(1:4,8) = 0d0
  Mom(1:4,9) = 0d0

  call gensing(4,Energy,(/M_Reso,0d0,0d0,0d0/),MomC,Pcol1,Pcol2,SingDepth,Steps)

  Mom(1:4,5)=MomC(1:4,1)
  Mom(1:4,6)=MomC(1:4,2)
  Mom(1:4,7)=MomC(1:4,3)
  Mom(1:4,10)=MomC(1:4,4)

  Mom(1:4,4) = Mom(1:4,6) + Mom(1:4,7)
  Mom(1:4,3) = Mom(1:4,4) + Mom(1:4,5)
  Jac=1d0

!print *,"===================="
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
!print *,"===================="

  RETURN
END SUBROUTINE



END MODULE
