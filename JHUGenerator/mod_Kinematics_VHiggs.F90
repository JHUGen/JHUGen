MODULE ModKinematics_VHiggs
implicit none
save


contains








SUBROUTINE WriteOutEvent_VHiggs(id,helicity,MomExt,inv_mass,EventWeight)
use ModParameters
implicit none
double precision, intent(in) :: MomExt(1:4,1:9), inv_mass(9)
!double precision, intent(in) :: beam_h(2)
double precision helicity(9)
real(8) :: MomDummy(1:4,1:9), MassDummy(9)
integer :: id(9)
integer :: ICOLUP(4,2)
real(8) :: EventWeight
real(8) :: Spin
integer :: i
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"



MomDummy = MomExt/GeV
MassDummy = inv_mass/GeV

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

helicity(3:5)=0d0
Spin = 0.1d0

if(COLLIDER.eq.0)then
  ICOLUP=0
else
  if(id(1).ne.0.and.id(2).ne.0)then
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
  elseif(id(1).eq.0.and.id(2).eq.0)then
    id(1:2) = 21
    ICOLUP(1,1)=501
    ICOLUP(1,2)=502
    ICOLUP(2,1)=502
    ICOLUP(2,2)=501
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










SUBROUTINE EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,HDecays,PhoOnshell,ZAinterference)
use ModParameters
implicit none

      logical, intent(in) :: HDecays
      logical, intent(in), optional :: PhoOnshell
      logical, optional :: ZAinterference
      real(8), intent(in) :: yRnd(1:20),mass(9,2)
      real(8) :: phi, beta, gamma
      real(8) :: temp_vector(4), temp_boost(4)
      integer :: i
      real(8), parameter :: Twopi = 2d0 * Pi
      real(8) :: MomExt(1:4,1:9)
      real(8) :: MomDummy(1:4)
      !double precision four_momentum(7,4)
      real(8), intent(out) :: PSWgt,inv_mass(9)
! 1=E, 2,3,4=p_x,y,z
!psg_mass(mass, width)
      real(8) :: cm_abs3p(9)
      real(8) :: cm_sin_theta(9), cm_cos_theta(9)
      real(8) :: cm_sin_phi(9), cm_cos_phi(9)
!use Cauchy distribution for Breit-Wigner distribution for the invariant mass of 2 and 3?
!      logical, parameter :: breit_wigner = .true.
      real(8) :: jacobian4, jacobian5
      logical :: hasAonshell, hasInterference

      hasAonshell = .false.
      if(present(PhoOnshell)) then
         hasAonshell=PhoOnshell
      endif
      hasInterference = .false.
      if(present(ZAinterference)) then
         hasInterference=ZAinterference
      endif

      MomExt(:,4:9)=0d0
      inv_mass(4:9)=0d0

!333333333
      inv_mass(3) = dsqrt((MomExt(1,3)+MomExt(4,3))*(MomExt(1,3)-MomExt(4,3)))

!555555555
      if(HDecays)then
        inv_mass(5) = dsqrt(dabs(bw_sq(yRnd(13),mass(5,1), mass(5,2), inv_mass(3)**2, jacobian5)))
        jacobian5 = jacobian5 /16d0 /Pi**2 != (ds5/2pi)*(1/8pi)
      else
        inv_mass(5) = mass(5,1)
        jacobian5=1d0
      endif

!4444444444
      if(hasAonshell) then
        inv_mass(4) = mass(4,1)
        jacobian4=1d0
      elseif(.not.hasInterference)then
        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
        jacobian4 = jacobian4 /16d0 /Pi**2 != (ds4/2pi)*(1/8pi)
      else
        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
        !inv_mass(4) = dsqrt(one_over_s_sq(yRnd(12), MPhotonCutoff**2, (inv_mass(3)-inv_mass(5))**2, jacobian4))
        jacobian4 = jacobian4 /16d0 /Pi**2
      endif

!444444444444
!energy of 4 in the CM frame of 3
      MomExt(1,4)=(inv_mass(3)**2+(inv_mass(4)+inv_mass(5))*(inv_mass(4)-inv_mass(5)))/2d0/inv_mass(3)
!|3-momentum| of 4 in the CM frame of 3
      cm_abs3p(4) = dsqrt((MomExt(1,4)+inv_mass(4)) * (MomExt(1,4)-inv_mass(4)))
!generating cos(theta_4) and phi_4 in the CM frame of 3
      cm_cos_theta(4) = yRnd(6)
      cm_cos_theta(4) = cm_cos_theta(4)*2d0-1d0
      cm_sin_theta(4) = dsqrt((1d0+cm_cos_theta(4))  *(1d0-cm_cos_theta(4)))
      phi = yRnd(7)
      phi=Twopi*phi
      cm_cos_phi(4) = dcos(phi)
      cm_sin_phi(4) = dsin(phi)
!3-momentum of 4 in the CM frame of 3
      MomExt(2,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_cos_phi(4)
      MomExt(3,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_sin_phi(4)
      MomExt(4,4)=cm_abs3p(4)*cm_cos_theta(4)
!boost the 4-momentum of 4 to the lab frame
!x and y components stay the same
!z component
      beta = MomExt(4,3)/MomExt(1,3)
      gamma = 1d0/dsqrt((1d0+beta)*(1d0-beta))
      MomDummy(1:4)=MomExt(1:4,4)
      MomExt(4,4)=(MomDummy(4)+MomDummy(1)*beta) *gamma
!energy
      MomExt(1,4)=(MomDummy(1)+MomDummy(4)*beta) *gamma

!555555555555555555
      MomExt(1:4,5) = MomExt(1:4,3) - MomExt(1:4,4)

!666666666666666666
      if(.not.hasAonshell) then
!invariant mass of 6
         inv_mass(6)=0d0
!energy of 6 in the CM frame of 4
         MomExt(1,6)=inv_mass(4)/2d0
!|3-momentum| of 6 in the CM frame of 4
         cm_abs3p(6)=MomExt(1,6)
!generating cos(theta_6) and phi_6 in the CM frame of 4
!z-axis is along the boost of 2
         cm_cos_theta(6) = yRnd(8)
         cm_cos_theta(6) = cm_cos_theta(6)*2d0-1d0
         cm_sin_theta(6) = dsqrt((1d0+cm_cos_theta(6)) *(1d0-cm_cos_theta(6)))
         cm_cos_phi(6) = yRnd(9)
         phi=Twopi*cm_cos_phi(6)
         cm_cos_phi(6) = dcos(phi)
         cm_sin_phi(6) = dsin(phi)
!3-momentum of 6 in the CM frame of 4
         MomExt(2,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_cos_phi(6)
         MomExt(3,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_sin_phi(6)
         MomExt(4,6)=cm_abs3p(6)*cm_cos_theta(6)
!boost the 4-momentum of 6 to the lab frame
         temp_vector = MomExt(1:4,6)
         temp_boost = MomExt(1:4,4)
         call LORENTZ(temp_vector, temp_boost)
         MomExt(1:4,6) = temp_vector
      else
         MomExt(1:4,6)=MomExt(1:4,4)
      endif

!7777777777777777777777
!invariant mass of 7
      inv_mass(7)=0d0
!4-momentum of 7 (lab frame) by energy-momentum conservation
      MomExt(1:4,7)=MomExt(1:4,4)-MomExt(1:4,6)

!8888888888888888888888
      if(.not.HDecays)then
!invariant mass of 8
        inv_mass(8)=0d0
!energy of 8 in the CM frame of 5
        MomExt(1,8)=inv_mass(5)/2d0
!|3-momentum| of 8 in the CM frame of 5
        cm_abs3p(8)=MomExt(1,8)
!generating cos(theta_8) and phi_8 in the CM frame of 5
!z-axis is along the boost of 5
        cm_cos_theta(8) = yRnd(10)
        cm_cos_theta(8) = cm_cos_theta(8)*2d0-1d0
        cm_sin_theta(8) = dsqrt((1d0+cm_cos_theta(8)) *(1d0-cm_cos_theta(8)))
        cm_cos_phi(8) = yRnd(11)
        phi=Twopi*cm_cos_phi(8)
        cm_cos_phi(8) = dcos(phi)
        cm_sin_phi(8) = dsin(phi)
!3-momentum of 8 in the CM frame of 5
!x and y components are not necessary, yet
        MomExt(2,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_cos_phi(8)
        MomExt(3,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_sin_phi(8)
        MomExt(4,8)=cm_abs3p(8)*cm_cos_theta(8)
!boost the 4-momentum of 6 to the lab frame
        temp_vector = MomExt(1:4,8)
        temp_boost = MomExt(1:4,5)
        call LORENTZ(temp_vector, temp_boost)
        MomExt(1:4,8) = temp_vector
!9999999999999999999999
!4-momentum of 9 (lab frame) by energy-momentum conservation
        MomExt(1:4,9)=MomExt(1:4,5)-MomExt(1:4,8)
      endif

      PSWgt = jacobian4*jacobian5 * cm_abs3p(4)/(4d0*pi)/inv_mass(3)
      !print *,  "()",inv_mass, jacobian4, jacobian5, cm_abs3p(4), PSWgt, "()"

!do i=4,7
!print *, dsqrt(dabs(four_momentum(i,:).dot.four_momentum(i,:)))
!enddo
!pause

RETURN
END SUBROUTINE







! Breit-Wigner mass^2
function bw_sq(x, m, ga, smax, jacobian)
implicit none
real(8) :: bw_sq
real(8), intent(in) :: m, ga, smax,x
real(8) :: xmin, xmax,xprime
real(8), intent(out) :: jacobian

xmin=-datan(m/ga)/m/ga
xmax=-datan((-smax+m**2)/ga/m)/ga/m
xprime=x*(xmax-xmin)+xmin
bw_sq=m**2+dtan(xprime*ga*m)*ga*m
jacobian=(ga*m)**2 * (1d0+dtan(ga*m*xprime)**2) * (xmax-xmin)

return
end function bw_sq

!generate s according to 1/s^2 distribution
function one_over_s_sq(x, smin, smax, jacobian)
implicit none
real(8) :: one_over_s_sq
real(8), intent(in) :: x, smin, smax
real(8), intent(out) :: jacobian

one_over_s_sq=smin/(1d0-x*(smax-smin)/smax)
jacobian=smin*smax*(smax-smin)/((x-1d0)*smax-x*smin)**2

return
end function one_over_s_sq


!LORENTZ.F
!VERSION 20130123
!
!A subroutine that performs a general boost to a four vector
!(vector) based on another four vector (boost). The primed and
!unprimed frames have their axes in parallel to one another.
!Rotation is not performed by this subroutine.

      subroutine LORENTZ(vector, boost)

      implicit none

      double precision vector(4), boost(4)
      double precision lambdaMtrx(4,4), vector_copy(4)
      double precision beta(2:4), beta_sq, gamma
      integer i,j
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

!      double precision KRONECKER_DELTA
!      external KRONECKER_DELTA

      do i=2,4
        beta(i) = boost(i)/boost(1)
      enddo

      beta_sq = beta(2)**2+beta(3)**2+beta(4)**2

  if(beta_sq.ge.epsilon)then

      gamma = 1d0/dsqrt(1d0-beta_sq)

      lambdaMtrx(1,1) = gamma

      do i=2,4
        lambdaMtrx(1,i) = gamma*beta(i)
        lambdaMtrx(i,1) = lambdaMtrx(1,i)
      enddo

      do i=2,4
      do j=2,4
        lambdaMtrx(i,j) = (gamma-1d0)*beta(i)*beta(j)/beta_sq + KRONECKER_DELTA(i,j)
      enddo
      enddo


!apply boost to vector1
      vector_copy = vector
      vector = 0d0
      do i=1,4
      do j=1,4
        vector(i) = vector(i) + lambdaMtrx(i,j)*vector_copy(j)
      enddo
      enddo
  endif

      return
      END subroutine LORENTZ



! KRONECKER_DELTA.F
!
! KRONECKER_DELTA(i,j)
! A function that returns 1 if i=j, and 0 otherwise.
      double precision function KRONECKER_DELTA(i,j)
      integer i,j
      if(i.eq.j)then
        KRONECKER_DELTA = 1d0
      else
        KRONECKER_DELTA = 0d0
      endif

      return
      end function KRONECKER_DELTA




SUBROUTINE getHiggsDecayLength(ctau)
use ModParameters
implicit none
real(8) :: x,xp,xpp,Len0,propa,ctau,ctau0
integer :: loop


     ctau  = 0d0
     ctau0 = HiggsDecayLengthMM
     if( ctau0.lt.1d-16 ) RETURN

     do loop=1,4000000!  4Mio. tries otherwise return zero
          call random_number(x)
          xp = 10*x*ctau0             ! scan between 0..10*ctau0

          propa = dexp( -xp/(ctau0) ) ! the max of propa is 1.0
          call random_number(xpp)
          if( xpp.lt.propa ) then!   accept
                ctau = xp
                RETURN
          endif
     enddo


RETURN
END SUBROUTINE







END MODULE
