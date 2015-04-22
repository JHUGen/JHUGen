! For the amplitude squared for the following process,
! q/e+(1) q'/e-(2) > Z/W(3) > Z/W(4) H(5), Z/W(4) > l-/nu(6) l+/nu~(7), H(5) > b(8) b~(9)
! call EvalAmp_VHiggs(id,helicity,MomExt,mass,me2)
! id = PDG code in integer type.
! helicity = +/-1 in double precision type.
! Example:
! id(1:9) = (/-11, 11, 23, 23, 25, 11, -11, 5, -5/)
! helicity(1:9) = (/1d0, -1d0, 0d0, 0d0, 0d0, 1d0, -1d0, 1d0, 1d0/)
! for a e+ e- > Z > Z H, Z > e- e+, H > b b~

module ModVHiggs
  implicit none
  private
    !-- general definitions, to be merged with Markus final structure
  integer, parameter  :: dp = selected_real_kind(15)

  include './variables.F90'

  real(dp), parameter :: pi =3.141592653589793238462643383279502884197d0
  real(dp), parameter :: T3lL= -0.5d0
  real(dp), parameter :: T3lR=  0d0
  real(dp), parameter :: T3nL=  0.5d0
  real(dp), parameter :: T3uL= 0.5d0
  real(dp), parameter :: T3uR= 0d0
  real(dp), parameter :: T3dL= -0.5d0
  real(dp), parameter :: T3dR= 0d0
  real(dp), parameter :: QlL = -1d0
  real(dp), parameter :: QlR = -1d0
  real(dp), parameter :: QnL =  0d0
  real(dp), parameter :: QuL = 0.66666666666666666666666666667d0
  real(dp), parameter :: QuR = 0.66666666666666666666666666667d0
  real(dp), parameter :: QdL = -0.33333333333333333333333333333d0
  real(dp), parameter :: QdR = -0.33333333333333333333333333333d0

  !spin-0 couplings
  real(dp), parameter :: gFFS=1d0
  real(dp), parameter :: gFFP=0d0

  !----- notation for subroutines
  public :: EvalAmp_VHiggs

contains



subroutine EvalAmp_VHiggs(id,helicity,MomExt,vvcoupl,mass,me2)
      real(dp), intent(in) :: MomExt(1:4,1:9) !beam_momentum(2,4),four_momentum(7,4)
      complex(dp) :: vvcoupl(32)
      !real(dp) :: MomExt(1:4,1:9)
      real(dp) :: inv_mass(9)
      !real(dp) :: inv_mass(9)
      real(dp), intent(out) :: me2
      real(dp), intent(in) :: helicity(9)!, beam_h(2)
      real(dp), intent(in) :: mass(9,2) !(mass, width)
      integer, intent(in) :: id(9)
      integer :: i
      complex(dp) amplitude
      inv_mass=0d0
      do i=3,5
        inv_mass(i) = dsqrt(MomExt(1,i)**2 - MomExt(2,i)**2 - MomExt(3,i)**2 - MomExt(4,i)**2)
      enddo

      amplitude=MATRIXELEMENT0(MomExt,inv_mass,mass,helicity,id,vvcoupl)
      me2=dble(amplitude*dconjg(amplitude))
!print *, me2, helicity(1:2), helicity(6:7)
    return

end subroutine EvalAmp_VHiggs










!MATRIXELEMENT0.F
!VERSION 20130710

!

      complex(dp) function MATRIXELEMENT0(MomExt,inv_mass,mass,helicity,id,vvcoupl)

      implicit none

      complex(dp) dMATRIXELEMENT
      real(dp), intent(in) :: MomExt(1:4,1:9) !,four_momentum(7,4)
      complex(8), intent(in) :: vvcoupl(32)
      real(dp), intent(in) :: inv_mass(9)
      real(dp) ::  mass(9,2)
      !real(dp), intent(in) ::  beam_momentum(2,4)
      real(dp), intent(in) ::  helicity(9)!, beam_h(2) !helicities
      integer, intent(in) ::  id(9)!, beam_id(2)

      integer mu1,mu2,mu3,mu4,lambda1,lambda2
      complex(dp) PVVX0P      
      complex(dp) Vcurrent1(4), Acurrent1(4), current1(4), Vcurrent2(4)
      complex(dp) Acurrent2(4), current2(4),POL1(3,4), POL2(3,4)
      complex(dp) g_mu_nu(4,4), pp(4,4), epp(4,4)
      complex(dp) VVX0(4,4)
      complex(dp) PROP1, PROP2, PROP3, qq, gVVP, gVVS1, gVVS2, gFFZ, gFFW
      complex(dp) ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      real(dp) q3_q3,q4_q4
      complex(dp) ghz1, ghz2, ghz3, ghz4
      
      if(id(3).eq.23)then
        mass(3,1)=M_Z
        mass(4,1)=M_Z
        mass(3,2)=Ga_Z
        mass(4,2)=Ga_Z
      else
        mass(3,1)=M_W
        mass(4,1)=M_W
        mass(3,2)=Ga_W
        mass(4,2)=Ga_W
      endif

      ghz1 = vvcoupl(1)
      ghz2 = vvcoupl(2)
      ghz3 = vvcoupl(3)
      ghz4 = vvcoupl(4)

      ghz1_prime = vvcoupl(5) 
      ghz1_prime2= vvcoupl(6) 
      ghz1_prime3= vvcoupl(7) 
      ghz1_prime4= vvcoupl(8)
      ghz1_prime5= vvcoupl(9)

      ghz2_prime = vvcoupl(10) 
      ghz2_prime2= vvcoupl(11)
      ghz2_prime3= vvcoupl(12)
      ghz2_prime4= vvcoupl(13)
      ghz2_prime5= vvcoupl(14)

      ghz3_prime = vvcoupl(15)
      ghz3_prime2= vvcoupl(16)
      ghz3_prime3= vvcoupl(17)
      ghz3_prime4= vvcoupl(18)
      ghz3_prime5= vvcoupl(19)

      ghz4_prime = vvcoupl(20)
      ghz4_prime2= vvcoupl(21)
      ghz4_prime3= vvcoupl(22)
      ghz4_prime4= vvcoupl(23)
      ghz4_prime5= vvcoupl(24)

      ghz1_prime6= vvcoupl(25)
      ghz1_prime7= vvcoupl(26)

      ghz2_prime6= vvcoupl(27)
      ghz2_prime7= vvcoupl(28)

      ghz3_prime6= vvcoupl(29)
      ghz3_prime7= vvcoupl(30)

      ghz4_prime6= vvcoupl(31)
      ghz4_prime7= vvcoupl(32)


      gFFZ = (0d0,1d0)*dsqrt(4d0*pi*alpha_QED/(1d0-sitW**2))/sitW
      gFFW = (0d0,1d0)*dsqrt(2d0*pi*alpha_QED)/sitW
!qq = s in the paper      
      qq=-MomExt(1,3)*MomExt(1,4) +MomExt(2,3)*MomExt(2,4) +MomExt(3,3)*MomExt(3,4) +MomExt(4,3)*MomExt(4,4)
!narrow-width approximation
!      qq=(H_mass**2-Z_mass**2-(inv_mass(1)**2))/2d0
!     print *,qq

      PROP1 = PROPAGATOR(inv_mass(3),mass(3,1),mass(3,2))
      PROP2 = PROPAGATOR(inv_mass(4),mass(4,1),mass(4,2))
      PROP3 = PROPAGATOR(inv_mass(5),mass(5,1),mass(5,2))

      if(id(1).gt.0)then
        call FFV(id(2), MomExt(:,2), helicity(2), id(1), MomExt(:,1), helicity(1), Vcurrent1)
        if((id(1)+id(2)).eq.0)then
          call FFA(id(2), MomExt(:,2), helicity(2), id(1), MomExt(:,1), helicity(1), Acurrent1)
        endif
      else
        call FFV(id(1), MomExt(:,1), helicity(1), id(2), MomExt(:,2), helicity(2), Vcurrent1)
        if((id(1)+id(2)).eq.0)then
          call FFA(id(1), MomExt(:,1), helicity(1), id(2), MomExt(:,2), helicity(2), Acurrent1)
        endif
      endif

      if(id(6).gt.0)then
        call FFV(id(6), MomExt(:,6), helicity(6), id(7), MomExt(:,7), helicity(7), Vcurrent2)
        if((id(6)+id(7)).eq.0)then
          call FFA(id(6), MomExt(:,6), helicity(6), id(7), MomExt(:,7), helicity(7), Acurrent2)
        endif
      else
        call FFV(id(7), MomExt(:,7), helicity(7), id(6), MomExt(:,6), helicity(6), Vcurrent2)
        if((id(6)+id(7)).eq.0)then
          call FFA(id(7), MomExt(:,7), helicity(7), id(6), MomExt(:,6), helicity(6), Acurrent2)
        endif
      endif

!WH
      if((id(1)+id(2)).ne.0)then
        if((id(1)*helicity(1)).le.0d0)then
          current1=Vcurrent1*gFFW*CKM(abs(id(1)),abs(id(2)))
        else
          current1=0d0
        endif
        current2=Vcurrent2*gFFW*CKM(abs(id(6)),abs(id(7)))

!ZH
      else if((abs(id(1)).eq.11).or.(abs(id(1)).eq.13))then
!        print *, beam_id, beam_h
!e+ e- Z vertex for incoming states
        if((id(1)*helicity(1)).gt.0d0)then
          current1=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent1 -(0.5d0*T3lR)*Acurrent1
        else
          current1=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent1 -(0.5d0*T3lL)*Acurrent1
        endif
        current1=current1*gFFZ

!u u~ Z vertex for incoming states
      else if((abs(id(1)).eq.2).or.(abs(id(1)).eq.4))then
        if((id(1)*helicity(1)).gt.0d0)then
          current1=(0.5d0*T3uR - QuR*sitW**2) *Vcurrent1 -(0.5d0*T3uR)*Acurrent1
        else
          current1=(0.5d0*T3uL - QuL*sitW**2) *Vcurrent1 -(0.5d0*T3uL)*Acurrent1
        endif
        current1=current1*gFFZ
!d d~ Z vertex for incoming states
      else if((abs(id(1)).eq.1).or.(abs(id(1)).eq.3).or.(abs(id(1)).eq.5))then
        if((id(1)*helicity(1)).gt.0d0)then
          current1=(0.5d0*T3dR - QdR*sitW**2) *Vcurrent1 -(0.5d0*T3dR)*Acurrent1
        else
          current1=(0.5d0*T3dL - QdL*sitW**2) *Vcurrent1 -(0.5d0*T3dL)*Acurrent1
        endif
        current1=current1*gFFZ

      else
      current1=0d0
      print *, "invalid incoming state"

      endif

!f+ f- Z vertex for final states
      if((id(6)+id(7)).eq.0)then
!l+ l- Z vertex for final state
        if((abs(id(6)).eq.11).or.(abs(id(6)).eq.13).or.(abs(id(6)).eq.15))then
          if((id(6)*helicity(6)).gt.0d0)then
            current2=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent2 -(0.5d0*T3lR)*Acurrent2
          else
            current2=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent2 -(0.5d0*T3lL)*Acurrent2
          endif
          current2=current2*gFFZ*dsqrt(scale_alpha_Z_ll)

!u u~ Z vertex for final state
        else if((abs(id(6)).eq.2).or.(abs(id(6)).eq.4))then
          if((id(6)*helicity(6)).gt.0d0)then
            current2=(0.5d0*T3uR - QuR*sitW**2) *Vcurrent2 -(0.5d0*T3uR)*Acurrent2
          else
            current2=(0.5d0*T3uL - QuL*sitW**2) *Vcurrent2 -(0.5d0*T3uL)*Acurrent2
          endif
          current2=current2*gFFZ*dsqrt(scale_alpha_Z_uu)

!d d~ Z vertex for final state
        else if((abs(id(6)).eq.1).or.(abs(id(6)).eq.3).or.(abs(id(6)).eq.5))then
          if((id(6)*helicity(6)).gt.0d0)then
            current2=(0.5d0*T3dR - QdR*sitW**2) *Vcurrent2 -(0.5d0*T3dR)*Acurrent2
          else
            current2=(0.5d0*T3dL - QdL*sitW**2) *Vcurrent2 -(0.5d0*T3dL)*Acurrent2
          endif
          current2=current2*gFFZ*dsqrt(scale_alpha_Z_dd)

!nu nu~ Z vertex for final state        
        else if((abs(id(6)).eq.12).or.(abs(id(6)).eq.14).or.(abs(id(6)).eq.16))then
          current2=(0.5d0*T3nL - QnL*sitW**2) *Vcurrent2 -(0.5d0*T3nL)*Acurrent2
          current2=current2*gFFZ*dsqrt(scale_alpha_Z_nn)

        else
        current2=0d0
        print *, "invalid final state"!, id(4:5)
        stop
        endif

      endif
!f+ f- Z vertex for final states END

      call POLARIZATION(MomExt(:,3), POL1)
      call POLARIZATION(MomExt(:,4), POL2)
   

!ZZX vertex
      q3_q3 = inv_mass(3)**2
      q4_q4 = inv_mass(4)**2
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))  &
                        +   ghz1_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z1**2                                   &
                        +   ghz1_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z1**2                                   &
                        +   ghz1_prime4* inv_mass(5)**2 / Lambda_z1**2                                            &
                        +   ghz1_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z1**4                             &
                        +   ghz1_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z1**4                             &
                        +   ghz1_prime7* ( abs(q3_q3)*abs(q4_q4) )/Lambda_z1**4
      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))  &
                        +   ghz2_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z2**2                                   &
                        +   ghz2_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z2**2                                   &
                        +   ghz2_prime4* inv_mass(5)**2 / Lambda_z2**2                                            &
                        +   ghz2_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z2**4                             &
                        +   ghz2_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z2**4                             &
                        +   ghz2_prime7* ( abs(q3_q3)*abs(q4_q4) )/Lambda_z4**4
      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))  &
                        +   ghz3_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z3**2                                   &
                        +   ghz3_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z3**2                                   &
                        +   ghz3_prime4* inv_mass(5)**2 / Lambda_z3**2                                            &
                        +   ghz3_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z3**4                             &
                        +   ghz3_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z3**4                             &
                        +   ghz3_prime7* ( abs(q3_q3)*abs(q4_q4) )/Lambda_z3**4
      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))  &
                        +   ghz4_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z4**2                                   &
                        +   ghz4_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z4**2                                   &
                        +   ghz4_prime4* inv_mass(5)**2 / Lambda_z4**2                                            &
                        +   ghz4_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z4**4                             &
                        +   ghz4_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z4**4                             &
                        +   ghz4_prime7* ( abs(q3_q3)*abs(q4_q4) )/Lambda_z4**4

      gVVS1 = ghz1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghz2_dyn + ghz3_dyn*qq/Lambda )
      gVVS2 = -( 2d0*ghz2_dyn + ghz3_dyn*qq/Lambda )
      gVVP = -2d0*ghz4_dyn
      
      VVX0 = 0d0
      if(gVVS1.ne.0d0)then
        call VVS1(g_mu_nu)
        VVX0 = VVX0 + gVVS1*g_mu_nu
      endif

      if(gVVS2.ne.0d0)then
        call VVS2(MomExt(:,5),MomExt(:,5),pp)
        VVX0 = VVX0 + gVVS2*pp
      endif

      if(gVVP.ne.0d0)then
        call VVP(-MomExt(:,3),MomExt(:,4),epp)
        VVX0 = VVX0 + gVVP*epp
      endif

      VVX0 = (0d0,1d0)/vev*VVX0


! assemble everything and get iM
      MATRIXELEMENT0=(0d0,0d0)
      dMATRIXELEMENT=(0d0,0d0)
      do lambda1=1,3
      do lambda2=1,3  
      dMATRIXELEMENT=current1(1)*dconjg(POL1(lambda1, 1)) -current1(2)*dconjg(POL1(lambda1, 2)) -current1(3)*dconjg(POL1(lambda1, 3))-current1(4)*dconjg(POL1(lambda1, 4))
      dMATRIXELEMENT=dMATRIXELEMENT *(current2(1)*POL2(lambda2, 1) -current2(2)*POL2(lambda2, 2) -current2(3)*POL2(lambda2, 3) -current2(4)*POL2(lambda2, 4))
      PVVX0P=(0d0,0d0)
      do mu3=1,4
      do mu4=1,4
      PVVX0P=PVVX0P +POL1(lambda1,mu3)*VVX0(mu3,mu4)*dconjg(POL2(lambda2,mu4))
      enddo !mu4
      enddo !mu3
      dMATRIXELEMENT=dMATRIXELEMENT*PVVX0P
      MATRIXELEMENT0=MATRIXELEMENT0+dMATRIXELEMENT
      enddo !lambda2
      enddo !lambda1

      if(H_DK.eqv..false.)then
        MATRIXELEMENT0=MATRIXELEMENT0 *PROP1*PROP2*PROP3
      else
        MATRIXELEMENT0=MATRIXELEMENT0 *PROP1*PROP2*PROP3 &
        *(gFFS*FFS(id(8), MomExt(:,8), helicity(8), id(9), MomExt(:,9), helicity(9)) &
         +gFFP*FFP(id(8), MomExt(:,8), helicity(8), id(9), MomExt(:,9), helicity(9)))&
        *(0d0,-1d0)*m_bot/vev
      endif

      return
      END function MATRIXELEMENT0







!ANGLES.F
!VERSION 20130531
!
!A subroutine that calculates the polar and azimuthal angles of a 
!given vector(4) in terms of their sin and cos, which will be 
!returned by the array sincos(4).

      subroutine ANGLES(sincos, vector)
      implicit none
!     real(dp) Pi
      real(dp) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      real(dp) sincos(4), vector(4), abs3p, phi
!sincos(1)=cos(theta)
!sincos(2)=sin(theta)
!sincos(3)=cos(phi)
!sincos(4)=sin(phi)

!|3-momentum|
      abs3p = dsqrt(vector(2)**2+vector(3)**2+vector(4)**2)

!if |3-momentum|=0
      if(abs3p.lt.epsilon)then
        sincos(1)=1d0
        sincos(2)=0d0
      else
        sincos(1)=vector(4)/abs3p
        sincos(2)=dsqrt((1d0+sincos(1))*(1d0-sincos(1)))
      endif

!if colinear
      if(dabs(vector(3)).lt.epsilon)then
        phi=0d0
      else
        if(dabs(vector(2)).lt.epsilon)then
           phi=(TwoPi/2d0)/2d0 * dsign(1d0,vector(3))
        else
           phi=datan(vector(3)/vector(2))
        endif
      endif
!shift phi so that 0 < phi < 2Pi
      if(vector(2).lt.0d0)then
            phi=phi+Pi
      endif
      if(phi.lt.0d0)then
            phi=phi+Twopi
      endif
!     print *,phi
      sincos(3)=dcos(phi)
      sincos(4)=dsin(phi)

      return
      END subroutine ANGLES








!ANTISYMMETRIC2.F
!VERSION 20130702

!in epp(_mu, _nu) returns the 

      subroutine ANTISYMMETRIC2(p1,p2,epp)

      implicit none
!     real(dp) Pi
      real(dp) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision

      complex(dp) p1(4), p2(4)
      complex(dp) epp(4,4)
!      real(dp) ANTISYMMETRIC
      integer i,j,k,l

!      external ANTISYMMETRIC

!     do i=1,4
!     do j=1,4
      epp(i,j)=0d0
!     enddo
!     enddo

      do i=1,4
      do j=1,4
      do k=1,4
      do l=1,4
          epp(i,j)=epp(i,j)+ANTISYMMETRIC(i,j,k,l)*p1(k)*p2(l)
      enddo
      enddo
      enddo
      enddo

      return
      END subroutine ANTISYMMETRIC2









!ANTISYMMETRIC.F
!VERSION 20130618

!returns the element of the rank-4 COVARIANT total antysymmetric 
!tensor.
!ANTISYMMETRIC(0,1,2,3)=1.

      real(dp) function ANTISYMMETRIC(i,j,k,l)

      implicit none
!     include '../COMMON.INI'

      integer i,j,k,l

      ANTISYMMETRIC=-dble((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l))/12d0

      return
      END function ANTISYMMETRIC







!CONTRA_FIELD_TENSOR.F
!VERSION 20130630

!in T_mu_nu(lambda, ^mu, ^nu) returns the contrvariant field tensors
!for given polarization vectors.

      subroutine CONTRA_FIELD_TENSOR(POL, T_mu_nu)

      implicit none
!     real(dp) Pi
      real(dp) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      complex(dp) epep(4,4),emem(4,4),epe0(4,4),eme0(4,4),e0e0(4,4)
      complex(dp) epem(4,4),e0ep(4,4),e0em(4,4),emep(4,4)
      complex(dp) POL(3,4), T_mu_nu(5,4,4)
      
      call CONTRA_OUTER(POL(1,:), POL(1,:), epep)
      call CONTRA_OUTER(POL(2,:), POL(2,:), emem)
      call CONTRA_OUTER(POL(1,:), POL(3,:), epe0)
      call CONTRA_OUTER(POL(3,:), POL(1,:), e0ep)
      call CONTRA_OUTER(POL(2,:), POL(3,:), eme0)
      call CONTRA_OUTER(POL(3,:), POL(2,:), e0em)
      call CONTRA_OUTER(POL(3,:), POL(3,:), e0e0)
      call CONTRA_OUTER(POL(1,:), POL(2,:), epem)
      call CONTRA_OUTER(POL(2,:), POL(1,:), emep)

!lambda = +2
      T_mu_nu(1,:,:)=epep
!lambda = -2
      T_mu_nu(2,:,:)=emem
!lambda = +3
      T_mu_nu(3,:,:)=(epe0+e0ep)/dsqrt(2d0)
!lambda = -3
      T_mu_nu(4,:,:)=(eme0+e0em)/dsqrt(2d0)
!lambda = 0
      T_mu_nu(5,:,:)=(epem+emep)/dsqrt(6d0) + e0e0/dsqrt(1.5d0)


      return
      END subroutine CONTRA_FIELD_TENSOR









!CONTRA_OUTER.F
!VERSION 20130620

!in pp(_mu, _nu) returns the CONTRAVARIANT tensor of
!p1 p2 outer product.

      subroutine CONTRA_OUTER(p1,p2,pp)

      implicit none
!     include '../COMMON.INI'
      complex(dp) p1(4), p2(4)
      complex(dp) pp(4,4)
      integer mu, nu

      do mu=1,4
        do nu=1,4
          pp(mu,nu)=p1(mu)*p2(nu)
        enddo
      enddo

      return
      END subroutine CONTRA_OUTER









!COVARIANT_FIELD_TENSOR.F
!VERSION 20130630

!in T_mu_nu(lambda, _mu, _nu) returns the covariant field tensors
!for given polarization vectors.

      subroutine COVARIANT_FIELD_TENSOR(POL, T_mu_nu)

      implicit none
!     real(dp) Pi
      real(dp) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      complex(dp) epep(4,4),emem(4,4),epe0(4,4),eme0(4,4),e0e0(4,4)
      complex(dp) epem(4,4),e0ep(4,4),e0em(4,4),emep(4,4)
      complex(dp) POL(3,4), T_mu_nu(5,4,4)
      
      call COVARIANT_OUTER(POL(1,:), POL(1,:), epep)
      call COVARIANT_OUTER(POL(2,:), POL(2,:), emem)
      call COVARIANT_OUTER(POL(1,:), POL(3,:), epe0)
      call COVARIANT_OUTER(POL(3,:), POL(1,:), e0ep)
      call COVARIANT_OUTER(POL(2,:), POL(3,:), eme0)
      call COVARIANT_OUTER(POL(3,:), POL(2,:), e0em)
      call COVARIANT_OUTER(POL(3,:), POL(3,:), e0e0)
      call COVARIANT_OUTER(POL(1,:), POL(2,:), epem)
      call COVARIANT_OUTER(POL(2,:), POL(1,:), emep)

!lambda = +2
      T_mu_nu(1,:,:)=epep
!lambda = -2
      T_mu_nu(2,:,:)=emem
!lambda = +3
      T_mu_nu(3,:,:)=(epe0+e0ep)/dsqrt(2d0)
!lambda = -3
      T_mu_nu(4,:,:)=(eme0+e0em)/dsqrt(2d0)
!lambda = 0
      T_mu_nu(5,:,:)=(epem+emep)/dsqrt(6d0) + e0e0/dsqrt(1.5d0)

      return
      END subroutine COVARIANT_FIELD_TENSOR








!COVARIANT_OUTER.F
!VERSION 20130620

!in pp(_mu, _nu) returns the COVARIANT tensor of p1 p2
!outer product.

      subroutine COVARIANT_OUTER(p1,p2,pp)

      implicit none
!     include '../COMMON.INI'
      complex(dp) p1(4), p2(4)
      complex(dp) pp(4,4)
      integer mu, nu

      do mu=1,4
        do nu=1,4
          pp(mu,nu)=p1(mu)*p2(nu)
          if( ( (mu.ne.1).and.(nu.eq.1) ).or. &
              ( (mu.eq.1).and.(nu.ne.1) ) )then
            pp(mu,nu)=-pp(mu,nu)
          endif
        enddo
      enddo

      return
      END subroutine COVARIANT_OUTER








!COVARIANT_VECTOR.F
!VERSION 20130703

!returns the component of the COVARIANT vector for given 4-vector
!and Lorentz index.

      complex(dp) function COVARIANT_VECTOR(p,mu)

      implicit none

      complex(dp) p(4)
      integer mu

      if(mu.ne.1)then
        COVARIANT_VECTOR = -p(mu)
      else
        COVARIANT_VECTOR = p(mu)
      endif

      return
      END function COVARIANT_VECTOR










!FFP.A
!VERSION 20130522

!returns i.Psi~(p1,s1).gamma5.Psi(p2,s2) for massless states

      complex(dp) function FFP(pdg_code1, p1, h1, pdg_code2, p2, h2)

      implicit none
      real(dp), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      real(dp) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(dp) sqrt_pp1Dpp2

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).gt.0d0)then
        FFP=0d0

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. (   dabs( p2(1)+p2(4) ).lt.epsilon ) )then
        FFP=0d0
      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then
        FFP=-dsqrt(2d0*p1(1)*(p2(1)+p2(4)))
      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then
        FFP= dsqrt(2d0*p2(1)*(p1(1)+p1(4)))
      else
        sqrt_pp1Dpp2 = dsqrt((p1(1)+p1(4))/(p2(1)+p2(4)))
        FFP=(p2(2)-(0d0,1d0)*p2(3))*sqrt_pp1Dpp2- (p1(2)-(0d0,1d0)*p1(3))/sqrt_pp1Dpp2
      endif

      FFP=FFP*(0d0,-1d0)
      
      if( (dble(pdg_code1)*h1) .lt. 0d0)then
        FFP=-dconjg(FFP)
        
      endif    

      return
      END function FFP









!FFA.F
!VERSION 20130523

!in Acurrent(4) returns Psi~(p1,s1).gamma^mu.gamma5.Psi(p2,s2) for massless
!states.

      subroutine FFA(pdg_code1, p1, h1, pdg_code2, p2, h2, Acurrent)

      implicit none
      real(dp), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      real(dp) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(dp) sqrt_pp1Dpp2, sqrt_pp1Xpp2
      complex(dp) Acurrent(4)
      integer mu

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).lt.0d0)then
        do mu=1,4
          Acurrent(mu)=0d0
        enddo

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. &
               ( dabs( p2(1)+p2(4) ).lt.epsilon ) )then

        Acurrent(1)= 2d0*dsqrt(p1(1)*p2(1))
        Acurrent(2)= 0d0
        Acurrent(3)= 0d0
        Acurrent(4)=-Acurrent(1)

      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then

        Acurrent(1)= dsqrt( 2d0*p1(1) / ( p2(1)+p2(4) ) ) &
                    *( p2(2) + (0d0,1d0)*p2(3) )
        Acurrent(2)= dsqrt( 2d0*p1(1) * ( p2(1)+p2(4) ) )
        Acurrent(3)= (0d0,1d0)*Acurrent(2)
        Acurrent(4)=-Acurrent(1)

      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then

        Acurrent(1)= dsqrt( 2d0*p2(1) / ( p1(1)+p1(4) ) ) &
                    *( p1(2) - (0d0,1d0)*p1(3) )
        Acurrent(2)= dsqrt( 2d0*p2(1) * ( p1(1)+p1(4) ) )
        Acurrent(3)=-(0d0,1d0)*Acurrent(2)
        Acurrent(4)=-Acurrent(1)

      else

        sqrt_pp1Dpp2= dsqrt( (p1(1)+p1(4)) / (p2(1)+p2(4)) )
        sqrt_pp1Xpp2= dsqrt( (p1(1)+p1(4)) * (p2(1)+p2(4)) )
        Acurrent(1)= sqrt_pp1Xpp2         &
                    +( p1(2) - (0d0,1d0)*p1(3) )     &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
        Acurrent(2)= ( p1(2) - (0d0,1d0)*p1(3) )/sqrt_pp1Dpp2   &
                    +( p2(2) + (0d0,1d0)*p2(3) )*sqrt_pp1Dpp2
        Acurrent(3)= ( (0d0,1d0)*p1(2) + p1(3) )/sqrt_pp1Dpp2   &
                    -( (0d0,1d0)*p2(2) - p2(3) )*sqrt_pp1Dpp2
        Acurrent(4)=sqrt_pp1Xpp2          &
                    -( p1(2) - (0d0,1d0)*p1(3) )                &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
      endif

!     print *, pdg_code1,h1,dble(pdg_code1)*h1,'!'
      if( (dble(pdg_code1)*h1) .lt. 0d0)then
!       print *, Acurrent
        do mu=1,4
          Acurrent(mu)=-dconjg(Acurrent(mu))
        enddo       
      endif    

      return
      END subroutine FFA










!FFS.F
!VERSION 20130522

!returns Psi~(p1,s1).Psi(p2,s2) for massless states

      complex(dp) function FFS(pdg_code1, p1, h1, pdg_code2, p2, h2)

      implicit none
      real(dp), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      real(dp) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(dp) sqrt_pp1Dpp2

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).gt.0d0)then
        FFS=0d0

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. (   dabs( p2(1)+p2(4) ).lt.epsilon ) )then
        FFS=0d0
      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then
        FFS=-dsqrt(2d0*p1(1)*(p2(1)+p2(4)))
      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then
        FFS= dsqrt(2d0*p2(1)*(p1(1)+p1(4)))
      else
        sqrt_pp1Dpp2 = dsqrt((p1(1)+p1(4))/(p2(1)+p2(4)))
        FFS=(p2(2)-(0d0,1d0)*p2(3))*sqrt_pp1Dpp2- (p1(2)-(0d0,1d0)*p1(3))/sqrt_pp1Dpp2
      endif

      if( (dble(pdg_code1)*h1) .lt. 0d0)then
        FFS=-dconjg(FFS)
        
      endif    

      return
      END function FFS










!FFV.F
!VERSION 20130523

!in Vcurrent(4) returns Psi~(p1,s1).gamma^mu.Psi(p2,s2) for massless
!states.

      subroutine FFV(pdg_code1, p1, h1, pdg_code2, p2, h2, Vcurrent)

      implicit none
      real(dp), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      real(dp) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(dp) sqrt_pp1Dpp2, sqrt_pp1Xpp2
      complex(dp) Vcurrent(4)
      integer mu

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).lt.0d0)then
        do mu=1,4
          Vcurrent(mu)=0d0
        enddo

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. ( dabs( p2(1)+p2(4) ).lt.epsilon ) )then

        Vcurrent(1)= 2d0*dsqrt(p1(1)*p2(1))
        Vcurrent(2)= 0d0
        Vcurrent(3)= 0d0
        Vcurrent(4)=-Vcurrent(1)

      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then

        Vcurrent(1)= dsqrt( 2d0*p1(1) / ( p2(1)+p2(4) ) ) *( p2(2) + (0d0,1d0)*p2(3) )
        Vcurrent(2)= dsqrt( 2d0*p1(1) * ( p2(1)+p2(4) ) )
        Vcurrent(3)= (0d0,1d0)*Vcurrent(2)
        Vcurrent(4)=-Vcurrent(1)

      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then

        Vcurrent(1)= dsqrt( 2d0*p2(1) / ( p1(1)+p1(4) ) ) *( p1(2) - (0d0,1d0)*p1(3) )
        Vcurrent(2)= dsqrt( 2d0*p2(1) * ( p1(1)+p1(4) ) )
        Vcurrent(3)=-(0d0,1d0)*Vcurrent(2)
        Vcurrent(4)=-Vcurrent(1)

      else

        sqrt_pp1Dpp2= dsqrt( (p1(1)+p1(4)) / (p2(1)+p2(4)) )
        sqrt_pp1Xpp2= dsqrt( (p1(1)+p1(4)) * (p2(1)+p2(4)) )
        Vcurrent(1)= sqrt_pp1Xpp2  &
                    +( p1(2) - (0d0,1d0)*p1(3) )   &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
        Vcurrent(2)= ( p1(2) - (0d0,1d0)*p1(3) )/sqrt_pp1Dpp2  &
                    +( p2(2) + (0d0,1d0)*p2(3) )*sqrt_pp1Dpp2
        Vcurrent(3)= ( (0d0,1d0)*p1(2) + p1(3) )/sqrt_pp1Dpp2  &
                    -( (0d0,1d0)*p2(2) - p2(3) )*sqrt_pp1Dpp2
        Vcurrent(4)=sqrt_pp1Xpp2  &
                    -( p1(2) - (0d0,1d0)*p1(3) )   &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
      endif

      if( (dble(pdg_code1)*h1) .lt. 0d0)then
        do mu=1,4
          Vcurrent(mu)=dconjg(Vcurrent(mu))
        enddo       
      endif    

      return
      END subroutine FFV










!INNER.F
!VERSION 20130620

!returns the (complex) inner product of p1 and p2, p1_mu p2^nu.

      complex(dp) function INNER(p1,p2)

      implicit none
      complex(dp) p1(4), p2(4)

      INNER = p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3) - p1(4)*p2(4)

      return
      END function INNER











!INV_LORENTZ.F
!VERSION 20130602
!
!A subroutine that performs a general inverse boost to a four vector
!(vector) based on another four vector (boost). The primed and
!unprimed frames have their axes in parallel to one another.
!Rotation is not performed by this subroutine.

      subroutine INV_LORENTZ(vector, boost)

      implicit none

      real(dp) vector(4), boost(4) 
      real(dp) lambda(4,4), vector_copy(4)
      real(dp) beta(2:4), beta_sq, gamma
      integer i,j

      do i=2,4
        beta(i) = -boost(i)/boost(1)
      enddo

      beta_sq = beta(2)**2+beta(3)**2+beta(4)**2

      gamma = 1d0/dsqrt(1d0-beta_sq)

      lambda(1,1) = gamma

      do i=2,4
        lambda(1,i) = gamma*beta(i)
        lambda(i,1) = lambda(1,i)
      enddo

      do i=2,4
      do j=2,4
        lambda(i,j) = (gamma-1d0)*beta(i)*beta(j)/beta_sq + KRONECKER_DELTA(i,j)
      enddo
      enddo

!apply boost to vector1
      vector_copy = vector
      vector = 0d0
      do i=1,4
      do j=1,4
        vector(i) = vector(i) + lambda(i,j)*vector_copy(j)
      enddo
      enddo

      return
      END subroutine INV_LORENTZ











!KRONECKER_DELTA.F

!KRONECKER_DELTA(i,j)
!A function that returns 1 if i=j, and 0 otherwise.
      real(dp) function KRONECKER_DELTA(i,j)
      integer i,j
      if(i.eq.j)then
        KRONECKER_DELTA = 1d0
      else
        KRONECKER_DELTA = 0d0
      endif

      return
      end function KRONECKER_DELTA











!METRIC.F
!VERSION 20130524

!in METRIC returns the element of the Minkovski metric, with
!signature (1,-1,-1,-1), for given (_mu, _nu) or given (^mu, ^nu).

      real(dp) function METRIC(mu,nu)

      implicit none
      integer mu, nu

      if(mu.ne.nu)then
        METRIC=0d0
      else if(mu.eq.1)then
        METRIC=1d0
      else
        METRIC=-1d0
      endif

      return
      END function METRIC










!`.F
!VERSION 20130524

!in POL(lambda,^mu) returns the polarization vectors for given
!4-momentum p

      subroutine POLARIZATION(p, POL)

      implicit none
      real(dp) p(4), sincos(4), inv_mass, abs3p
      complex(dp) POL(3,4)
!     integer lambda, mu
      
      call ANGLES(sincos, p)

!lambda = +1
      POL(1,1)= 0d0
      POL(1,2)= (-sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
      POL(1,3)= (-sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
      POL(1,4)= sincos(2)/dsqrt(2d0)
!lambda = -1
      POL(2,1)= 0d0
      POL(2,2)= ( sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
      POL(2,3)= ( sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
      POL(2,4)= -sincos(2)/dsqrt(2d0)
!|3-momentum|
      abs3p = dsqrt(p(2)**2+p(3)**2+p(4)**2)
!invariant mass     
      inv_mass= dsqrt(p(1)**2-abs3p**2)
!lambda = L
      POL(3,1)= abs3p/inv_mass
      POL(3,2)= sincos(3)*sincos(2)*p(1)/inv_mass
      POL(3,3)= sincos(4)*sincos(2)*p(1)/inv_mass
      POL(3,4)= sincos(1)*p(1)/inv_mass

      return
      END subroutine POLARIZATION








!POLARIZATIONA.F
!VERSION 20130529

!in POL(lambda,^mu) returns the photon polarization vectors for given
!4-momentum p

      subroutine POLARIZATIONA(p, POL)

      implicit none
      real(dp) p(4), sincos(4), inv_mass, abs3p
      complex(dp) POL(2,4)
      
      call ANGLES(sincos, p)

!lambda = +1
      POL(1,1)= 0d0
      POL(1,2)= (-sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
      POL(1,3)= (-sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
      POL(1,4)= sincos(2)/dsqrt(2d0)
!lambda = -1
      POL(2,1)= 0d0
      POL(2,2)= ( sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
      POL(2,3)= ( sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
      POL(2,4)= -sincos(2)/dsqrt(2d0)


      return
      END subroutine POLARIZATIONA









!POLARIZATIONX.F
!VERSION 20130529

!in POL(lambda,^mu) returns the polarization vectors for given
!4-momentum p, with POL(L, ^mu) = (0, 0, 0, 1) when p is along the
!Z direction.

      subroutine POLARIZATIONX(p, POL)

      implicit none
      real(dp) p(4), sincos(4), inv_mass, abs3p
      complex(dp) POL(3,4)
!     integer lambda, mu
      
      call ANGLES(sincos, p)

!lambda = +1
      POL(1,1)= 0d0
      POL(1,2)= (-sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
      POL(1,3)= (-sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
      POL(1,4)= sincos(2)/dsqrt(2d0)
!lambda = -1
      POL(2,1)= 0d0
      POL(2,2)= ( sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
      POL(2,3)= ( sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
      POL(2,4)= -sincos(2)/dsqrt(2d0)
!|3-momentum|
      abs3p = dsqrt(p(2)**2+p(3)**2+p(4)**2)
!invariant mass     
      inv_mass= dsqrt(p(1)**2-abs3p**2)
!lambda = L
      POL(3,1)= 0d0
      POL(3,2)= sincos(3)*sincos(2)
      POL(3,3)= sincos(4)*sincos(2)
      POL(3,4)= sincos(1)

      return
      END subroutine POLARIZATIONX









!PROPAGATOR.F
!VERSION 20130522
!
!PROPAGATOR() returns the generi!complex-valued propagator
!without tensor structure (numerator), given mass, invariant mass
!and width.

      complex(dp) function PROPAGATOR(inv_mass, mass, width)
      implicit none

      real(dp) inv_mass, mass, width

!not assuming auto-conversion
!     PROPAGATOR = (0d0,1d0)/(dcmplx(inv_mass**2,0d0)
!    &            -dcmplx(mass**2,0d0)+
!    &             (0d0,1d0)*dcmplx(mass,0d0)*dcmplx(width,0d0))

!assuming auto-conversion. works with gfortran
      PROPAGATOR = (0d0,1d0) / ( inv_mass**2 - mass**2 + (0d0,1d0)*mass*width )
!     print *, PROPAGATOR

      return
      END function PROPAGATOR







!VVP.F
!VERSION 20130524

!in epp(_mu, _nu) returns the 

      subroutine VVP(p1,p2,epp)

      implicit none
      real(dp) p1(4), p2(4)
      complex(dp) epp(4,4)
!      real(dp) ANTISYMMETRIC
      integer i,j,k,l

!      external ANTISYMMETRIC

      do i=1,4
      do j=1,4
        epp(i,j)=0d0
      enddo
      enddo

      do i=1,4
      do j=1,4
      do k=1,4
      do l=1,4
          epp(i,j)=epp(i,j)+ANTISYMMETRIC(i,j,k,l)*p1(k)*p2(l)
      enddo
      enddo
      enddo
      enddo

      return
      END subroutine VVP







!VVS1.F
!VERSION 20130524

!in g_mu_nu(_mu, _nu) returns the VVS1 4*4 tensor.

      subroutine VVS1(g_mu_nu)

      implicit none
      complex(dp) g_mu_nu(4,4)

      g_mu_nu = 0d0

      g_mu_nu(1,1) =  1d0
      g_mu_nu(2,2) = -1d0
      g_mu_nu(3,3) = -1d0
      g_mu_nu(4,4) = -1d0

      return
      END subroutine VVS1








!VVS2.F
!VERSION 20130524

!in pp(_mu, _nu) returns the tensor of p1_mu * p2_nu.

      subroutine VVS2(p1,p2,pp)

      implicit none
      real(dp) p1(4), p2(4)
      complex(dp) pp(4,4)
      integer mu, nu

      do mu=1,4
        do nu=1,4
          pp(mu,nu)=p1(mu)*p2(nu)
          if( ( (mu.ne.1).and.(nu.eq.1) ).or. &
              ( (mu.eq.1).and.(nu.ne.1) ) )then
            pp(mu,nu)=-pp(mu,nu)
          endif
        enddo
      enddo

      return
      END subroutine VVS2






FUNCTION CKM(id1in,id2in)
implicit none
real(8) :: CKM
integer :: id1, id2, id1in, id2in
id1 = abs(id1in)
id2 = abs(id2in)
if((id1.eq.2  .and.  id2.eq.1)  .or.  (id1.eq.1  .and.  id2.eq.2))then
  CKM= 0.97425d0 * dsqrt(scale_alpha_W_ud)
elseif((id1.eq.2  .and.  id2.eq.3)  .or.  (id1.eq.3  .and.  id2.eq.2))then
  CKM= 0.2253d0
elseif((id1.eq.2  .and.  id2.eq.1)  .or.  (id1.eq.1  .and.  id2.eq.2))then
  CKM= 0.00413d0
elseif((id1.eq.4  .and.  id2.eq.1)  .or.  (id1.eq.1  .and.  id2.eq.4))then
  CKM= 0.225d0
elseif((id1.eq.4  .and.  id2.eq.3)  .or.  (id1.eq.3  .and.  id2.eq.4))then
  CKM= 0.986d0 * dsqrt(scale_alpha_W_cs)
elseif((id1.eq.4  .and.  id2.eq.1)  .or.  (id1.eq.1  .and.  id2.eq.4))then
  CKM= 0.0411d0
!elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Top_)))then
!  CKM= 0.0084d0
!elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.3)  .or.  (id1.eq.3  .and.  id2.eq.convertLHE(Top_)))then
!  CKM= 0.0400d0
!elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Top_)))then
!  CKM= 1.021
else
  CKM= 1d0! * dsqrt(scale_alpha_W_ln)
endif

END FUNCTION







end module ModVHiggs
!!--YaofuZhou-----------------------------------------