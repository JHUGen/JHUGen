MODULE ModKinematics_HH
implicit none
save


contains




SUBROUTINE WriteOutEvent_HH(id,helicity,Mom,EventWeight)
  use ModKinematics
  use ModParameters
  use ModMisc
  implicit none
  integer, intent(in) :: id(9)
  double precision, intent(in) :: helicity(9)
  double precision, intent(in) :: Mom(1:4,1:9)
  real(8), intent(in) :: EventWeight
  real(8) :: MomDummy(1:4,1:9), MassDummy(9)
  integer :: ICOLUP(4,2)
  real(8) :: Spin
  integer :: i
  integer :: NUP,IDPRUP
  real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
  character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1F3.0)"
  
  
  MomDummy = Mom/GeV
  MassDummy = Get_MInv(MomDummy)
  
  IDPRUP=Process
  SCALUP=Mu_Fact/GeV
  AQEDUP=alpha_QED
  AQCDUP=alphas

  call getHiggsDecayLength(HiggsDKLength)

  NUP = 8

  XWGTUP=EventWeight

  write(io_LHEOutFile,"(A)") "<event>"
  if( .not. ReadLHEFile ) write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

  Spin = 0.1d0

  write(io_LHEOutFile,fmt1) id(1), -1,0,0,503,504,MomDummy(2:4,1), MomDummy(1,1),        0.0d0,         0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(2), -1,0,0,504,503,MomDummy(2:4,2), MomDummy(1,2),        0.0d0,         0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(4),  2,1,2,  0,  0,MomDummy(2:4,4), MomDummy(1,4), MassDummy(4), HiggsDKLength, Spin
  write(io_LHEOutFile,fmt1) id(5),  2,1,2,  0,  0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength, Spin
  write(io_LHEOutFile,fmt1) id(6),  1,3,3,502,  0,MomDummy(2:4,6), MomDummy(1,6), MassDummy(6),         0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(7),  1,3,3,  0,502,MomDummy(2:4,7), MomDummy(1,7), MassDummy(7),         0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(8),  1,4,4,501,  0,MomDummy(2:4,8), MomDummy(1,8), MassDummy(8),         0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(9),  1,4,4,  0,501,MomDummy(2:4,9), MomDummy(1,9), MassDummy(9),         0.0d0, Spin

  write(io_LHEOutFile,"(A)") "</event>"


END SUBROUTINE






SUBROUTINE EvalPhasespace_HH(xRnd,Energy,Mom,id,Jac)
  use ModParameters
  use ModPhasespace
  use ModMisc
  implicit none
  real(8) :: Mom(1:4,1:9)
  real(8), intent(in) :: xRnd(5:13),Energy
! yrnd(5:6): phi_4 and cos(theta_4) in the CM frame of 1+2(3) > 45
! yrnd(7:8): phi_6 and cos(theta_6) in the CM frame of decay product of H(4) > 67
! yrnd(9:10): phi_8 and cos(theta_8) in the CM frame of decay product of H(5) > 89
! yRnd(11): inv_mass(4)
! yRnd(12): inv_mass(5)
! yRnd(13): swap momenta in PS for stability
  integer, intent(in) :: id(6:9)
  real(8), intent(out) :: Jac
  real(8) :: Jac1,Jac2,Jac3,Jac4,Jac5
  real(8) :: s67,s89

  Mom(:,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/) ! beam 1
  Mom(:,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/) ! beam 2
  Mom(:,3) = Energy * (/1d0,0d0,0d0,0d0/) ! 1+2=3

  !Jac2 = s_channel_propagator(M_Reso**2,Ga_Reso,(getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))**2,(Energy-M_Reso+5d0*Ga_Reso)**2,xRnd(11),s67) ! H --> 67
  !Jac3 = s_channel_propagator(M_Reso**2,Ga_Reso,(getMass(convertLHEreverse(id(8)))+getMass(convertLHEreverse(id(9))))**2,(Energy-dsqrt(s67))**2,xRnd(12),s89) ! H --> 89
  Jac2 = s_channel_propagator(M_Reso**2,Ga_Reso,(M_Reso-5d0*Ga_Reso)**2,(Energy-M_Reso+5d0*Ga_Reso)**2,xRnd(11),s67) ! H --> 67
  Jac3 = s_channel_propagator(M_Reso**2,Ga_Reso,(M_Reso-5d0*Ga_Reso)**2,(Energy-dsqrt(s67))**2,xRnd(12),s89) ! H --> 89
  !Jac2 = s_channel_propagator(0d0,0d0,0d0,(Energy)**2,xRnd(11),s67) ! H --> 67
  !Jac3 = s_channel_propagator(0d0,0d0,0d0,(Energy-dsqrt(s67))**2,xRnd(12),s89) ! H --> 89

  Jac1 = s_channel_decay(Mom(:,3),s67,s89,xRnd(5:6),Mom(:,4),Mom(:,5)) ! 3 --> HH

  if(includeInterference.eqv..false.)then
    Jac4 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,6),Mom(:,7)) !  H4 --> 67
    Jac5 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(9:10),Mom(:,8),Mom(:,9)) !  H5 --> 89
  else
    if( xRnd(13).lt.0.5d0 ) then
      Jac4 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,6),Mom(:,7)) !  H4 --> 67
      Jac5 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(9:10),Mom(:,8),Mom(:,9)) !  H5 --> 89
    else
      Jac4 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,6),Mom(:,9)) !  H4 --> 69
      Jac5 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(9:10),Mom(:,8),Mom(:,7)) !  H5 --> 87
    endif
  endif

  Jac = Jac1*Jac2*Jac3*Jac4*Jac5*PSNorm4

!  print*,"Jac2 ",Jac2
!  print*,"Jac3 ",Jac3
!  print*,"Jac4 ",Jac4
!  print*,"Jac5 ",Jac5
!  print*,"PSNorm4 ",PSNorm4
!print*,dsqrt(s67),dsqrt(s89),jac
  if( isNan(jac) ) then
     print *, "EvalPhasespace_HH NaN"
     print *, Energy
     print *, s67,s89
     print *, Jac1,Jac2,Jac3,Jac4,Jac5
     print *, xRnd
     Jac = 0d0
     print*,Mom(:,1)
     print*,Mom(:,2)
     print*,Mom(:,3)
     print*,Mom(:,4)
     print*,Mom(:,5)
     print*,Mom(:,6)
     print*,Mom(:,7)
     print*,Mom(:,8)
     print*,Mom(:,9)
     pause
  endif
!if(jac.gt.10000d0)then
!print*,dsqrt(s67),dsqrt(s89),jac
!print *,Jac1,Jac2,Jac3,Jac4,Jac5
!print*,Mom(:,1)
!print*,Mom(:,2)
!print*,Mom(:,3)
!print*,Mom(:,4)
!print*,Mom(:,5)
!print*,Mom(:,6)
!print*,Mom(:,7)
!print*,Mom(:,8)
!print*,Mom(:,9)
!print *,"===================="
!endif
RETURN
END SUBROUTINE




!SUBROUTINE EvalPhasespace_HH_old(yRnd,Energy,Mom,mass,PSWgt)
!use ModParameters
!implicit none
!! yrnd(5:6): phi_4 and cos(theta_4) in the CM frame of 1+2(3)
!! yrnd(7:8): phi_6 and cos(theta_6) in the CM frame of decay product of H(4)
!! yrnd(9:10): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
!! yRnd(11): inv_mass(4)
!! yRnd(12): inv_mass(5)
!      real(8), intent(in) :: yRnd(5:12),mass(3:5,1:2),Energy
!      real(8) :: phi, beta, gamma
!      real(8) :: temp_vector(4), temp_boost(4)
!      integer :: i
!      real(8), parameter :: Twopi = 2d0 * Pi
!      real(8) :: Mom(1:4,1:9)
!      real(8) :: MomDummy(1:4)
!      !double precision four_momentum(7,4)
!      real(8), intent(out) :: PSWgt
!      real(8) :: inv_mass(9)
!! 1=E, 2,3,4=p_x,y,z
!!psg_mass(mass, width)
!      real(8) :: cm_abs3p(9)
!      real(8) :: cm_sin_theta(9), cm_cos_theta(9)
!      real(8) :: cm_sin_phi(9), cm_cos_phi(9)
!!use Cauchy distribution for Breit-Wigner distribution for the invariant mass of 2 and 3?
!!      logical, parameter :: breit_wigner = .true.
!      real(8) :: jacobian4, jacobian5
!
!  Mom(:,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/) ! beam 1
!  Mom(:,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/) ! beam 2
!  Mom(:,3) = Energy * (/1d0,0d0,0d0,0d0/) ! 1+2=3
!
!
!      Mom(:,4:9)=0d0
!      inv_mass(4:9)=0d0
!
!!333333333
!      inv_mass(3) = dsqrt((Mom(1,3)+Mom(4,3))*(Mom(1,3)-Mom(4,3)))
!
!!555555555
!        inv_mass(5) = dsqrt(dabs(bw_sq(yRnd(12),mass(5,1), mass(5,2), inv_mass(3)**2, jacobian5)))
!        jacobian5 = jacobian5 /16d0 /Pi**2 != (ds5/2pi)*(1/8pi)
!
!
!!4444444444
!
!        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(11),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
!        jacobian4 = jacobian4 /16d0 /Pi**2 != (ds4/2pi)*(1/8pi)
!
!!444444444444
!!energy of 4 in the CM frame of 3
!      Mom(1,4)=(inv_mass(3)**2+(inv_mass(4)+inv_mass(5))*(inv_mass(4)-inv_mass(5)))/2d0/inv_mass(3)
!!|3-momentum| of 4 in the CM frame of 3
!      cm_abs3p(4) = dsqrt((Mom(1,4)+inv_mass(4)) * (Mom(1,4)-inv_mass(4)))
!!generating cos(theta_4) and phi_4 in the CM frame of 3
!      cm_cos_theta(4) = yRnd(6)
!      cm_cos_theta(4) = cm_cos_theta(4)*2d0-1d0
!      cm_sin_theta(4) = dsqrt((1d0+cm_cos_theta(4))  *(1d0-cm_cos_theta(4)))
!      phi = yRnd(5)
!      phi=Twopi*phi
!      cm_cos_phi(4) = dcos(phi)
!      cm_sin_phi(4) = dsin(phi)
!!3-momentum of 4 in the CM frame of 3
!      Mom(2,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_cos_phi(4)
!      Mom(3,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_sin_phi(4)
!      Mom(4,4)=cm_abs3p(4)*cm_cos_theta(4)
!!boost the 4-momentum of 4 to the lab frame
!!x and y components stay the same
!!z component
!      beta = Mom(4,3)/Mom(1,3)
!      gamma = 1d0/dsqrt((1d0+beta)*(1d0-beta))
!      MomDummy(1:4)=Mom(1:4,4)
!      Mom(4,4)=(MomDummy(4)+MomDummy(1)*beta) *gamma
!!energy
!      Mom(1,4)=(MomDummy(1)+MomDummy(4)*beta) *gamma
!
!!555555555555555555
!      Mom(1:4,5) = Mom(1:4,3) - Mom(1:4,4)
!
!!666666666666666666
!!invariant mass of 6
!         inv_mass(6)=0d0
!!energy of 6 in the CM frame of 4
!         Mom(1,6)=inv_mass(4)/2d0
!!|3-momentum| of 6 in the CM frame of 4
!         cm_abs3p(6)=Mom(1,6)
!!generating cos(theta_6) and phi_6 in the CM frame of 4
!!z-axis is along the boost of 2
!         cm_cos_theta(6) = yRnd(8)
!         cm_cos_theta(6) = cm_cos_theta(6)*2d0-1d0
!         cm_sin_theta(6) = dsqrt((1d0+cm_cos_theta(6)) *(1d0-cm_cos_theta(6)))
!         cm_cos_phi(6) = yRnd(7)
!         phi=Twopi*cm_cos_phi(6)
!         cm_cos_phi(6) = dcos(phi)
!         cm_sin_phi(6) = dsin(phi)
!!3-momentum of 6 in the CM frame of 4
!         Mom(2,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_cos_phi(6)
!         Mom(3,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_sin_phi(6)
!         Mom(4,6)=cm_abs3p(6)*cm_cos_theta(6)
!!boost the 4-momentum of 6 to the lab frame
!         temp_vector = Mom(1:4,6)
!         temp_boost = Mom(1:4,4)
!         call LORENTZ(temp_vector, temp_boost)
!         Mom(1:4,6) = temp_vector
!
!
!!7777777777777777777777
!!invariant mass of 7
!      inv_mass(7)=0d0
!!4-momentum of 7 (lab frame) by energy-momentum conservation
!      Mom(1:4,7)=Mom(1:4,4)-Mom(1:4,6)
!
!!8888888888888888888888
!!invariant mass of 8
!        inv_mass(8)=0d0
!!energy of 8 in the CM frame of 5
!        Mom(1,8)=inv_mass(5)/2d0
!!|3-momentum| of 8 in the CM frame of 5
!        cm_abs3p(8)=Mom(1,8)
!!generating cos(theta_8) and phi_8 in the CM frame of 5
!!z-axis is along the boost of 5
!        cm_cos_theta(8) = yRnd(10)
!        cm_cos_theta(8) = cm_cos_theta(8)*2d0-1d0
!        cm_sin_theta(8) = dsqrt((1d0+cm_cos_theta(8)) *(1d0-cm_cos_theta(8)))
!        cm_cos_phi(8) = yRnd(9)
!        phi=Twopi*cm_cos_phi(8)
!        cm_cos_phi(8) = dcos(phi)
!        cm_sin_phi(8) = dsin(phi)
!!3-momentum of 8 in the CM frame of 5
!!x and y components are not necessary, yet
!        Mom(2,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_cos_phi(8)
!        Mom(3,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_sin_phi(8)
!        Mom(4,8)=cm_abs3p(8)*cm_cos_theta(8)
!!boost the 4-momentum of 6 to the lab frame
!        temp_vector = Mom(1:4,8)
!        temp_boost = Mom(1:4,5)
!        call LORENTZ(temp_vector, temp_boost)
!        Mom(1:4,8) = temp_vector
!!9999999999999999999999
!!4-momentum of 9 (lab frame) by energy-momentum conservation
!        Mom(1:4,9)=Mom(1:4,5)-Mom(1:4,8)
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
!!!generate s according to 1/s^2 distribution
!!function one_over_s_sq(x, smin, smax, jacobian)
!!implicit none
!!real(8) :: one_over_s_sq
!!real(8), intent(in) :: x, smin, smax
!!real(8), intent(out) :: jacobian
!!
!!one_over_s_sq=smin/(1d0-x*(smax-smin)/smax)
!!jacobian=smin*smax*(smax-smin)/((x-1d0)*smax-x*smin)**2
!!
!!return
!!end function one_over_s_sq
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

END MODULE
