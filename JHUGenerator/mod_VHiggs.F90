!--YaofuZhou-----------------------------------------
module ModVHiggs
  use ModParameters
  implicit none
  private
  double precision, parameter :: T3lL= -0.5d0
  double precision, parameter :: T3lR=  0d0
  double precision, parameter :: T3nL=  0.5d0
  double precision, parameter :: T3uL= 0.5d0
  double precision, parameter :: T3uR= 0d0
  double precision, parameter :: T3dL= -0.5d0
  double precision, parameter :: T3dR= 0d0
  double precision, parameter :: QlL = -1d0
  double precision, parameter :: QlR = -1d0
  double precision, parameter :: QnL =  0d0
  double precision, parameter :: QuL = 0.66666666666666666666666666667d0
  double precision, parameter :: QuR = 0.66666666666666666666666666667d0
  double precision, parameter :: QdL = -0.33333333333333333333333333333d0
  double precision, parameter :: QdR = -0.33333333333333333333333333333d0

  !spin-0 couplings
  double precision, parameter :: gFFS=1d0
  double precision, parameter :: gFFP=0d0
  double precision, parameter :: b_Yukawa=4.18d0

  !----- notation for subroutines
  public :: EvalAmp_VHiggs

contains



  subroutine EvalAmp_VHiggs(yRnd,beam_id,id,beam_h,helicity,beam_momentum,four_momentum,inv_mass,mass,me2)
      real(8), intent(in) :: yRnd(1:20)
      real(8), intent(in) :: beam_momentum(2,4),four_momentum(7,4)
      real(8), intent(in) :: inv_mass(7)
      real(8), intent(out) :: me2
      double precision, intent(in) :: helicity(7), beam_h(2)
      double precision, intent(in) ::  mass(7,2) !(mass, width)
      integer, intent(in) ::  beam_id(2), id(7)
      double complex amplitude

      amplitude=MATRIXELEMENT0(four_momentum,inv_mass,mass,beam_momentum,helicity,beam_h,id,beam_id)
      me2=dble(amplitude*dconjg(amplitude))

    return

  end subroutine EvalAmp_VHiggs










!MATRIXELEMENT0.F
!VERSION 20130710

!

      double complex function MATRIXELEMENT0(four_momentum,inv_mass,mass,beam_momentum,helicity,beam_h,id,beam_id)

      implicit none

      double complex dMATRIXELEMENT
      double precision, intent(in) :: four_momentum(7,4)
      double precision, intent(in) :: inv_mass(7)
      double precision, intent(in) ::  mass(7,2)
      double precision, intent(in) ::  beam_momentum(2,4)
      double precision, intent(in) ::  helicity(7), beam_h(2) !helicities
      integer, intent(in) ::  id(7), beam_id(2)

      integer mu1,mu2,mu3,mu4,lambda1,lambda2
      double complex PVVX0P      
      double complex Vcurrent1(4), Acurrent1(4), current1(4), Vcurrent2(4)
      double complex Acurrent2(4), current2(4),POL1(3,4), POL2(3,4)
      double complex g_mu_nu(4,4), pp(4,4), epp(4,4)
      double complex VVX0(4,4)
      double complex PROP1, PROP2, PROP3, qq, gVVP, gVVS1, gVVS2, gFFZ, gFFW



      gFFZ = (0d0,1d0)*dsqrt(4d0*pi*alpha_QED/(1d0-sitW*2))/sitW
      gFFW = (0d0,1d0)*dsqrt(2d0*pi*alpha_QED)/sitW
!qq = s in the paper      
      qq=-four_momentum(1,1)*four_momentum(2,1) +four_momentum(1,2)*four_momentum(2,2) +four_momentum(1,3)*four_momentum(2,3) +four_momentum(1,4)*four_momentum(2,4)
!narrow-width approximation
!      qq=(H_mass**2-Z_mass**2-(inv_mass(1)**2))/2d0
!     print *,qq

      PROP1 = PROPAGATOR(inv_mass(1),mass(1,1),mass(1,2))
      PROP2 = PROPAGATOR(inv_mass(2),mass(2,1),mass(2,2))
      PROP3 = PROPAGATOR(inv_mass(3),mass(3,1),mass(3,2))

      if(beam_id(1).gt.0)then
        call FFV(beam_id(2), beam_momentum(2,:), beam_h(2), beam_id(1), beam_momentum(1,:), beam_h(1), Vcurrent1)
        if((beam_id(1)+beam_id(2)).eq.0)then
          call FFA(beam_id(2), beam_momentum(2,:), beam_h(2), beam_id(1), beam_momentum(1,:), beam_h(1), Acurrent1)
        endif
      else
        call FFV(beam_id(1), beam_momentum(1,:), beam_h(1), beam_id(2), beam_momentum(2,:), beam_h(2), Vcurrent1)
        if((beam_id(1)+beam_id(2)).eq.0)then
          call FFA(beam_id(1), beam_momentum(1,:), beam_h(1), beam_id(2), beam_momentum(2,:), beam_h(2), Acurrent1)
        endif
      endif

      if(id(4).gt.0)then
        call FFV(id(4), four_momentum(4,:), helicity(4), id(5), four_momentum(5,:), helicity(5), Vcurrent2)
        if((id(4)+id(5)).eq.0)then
          call FFA(id(4), four_momentum(4,:), helicity(4), id(5), four_momentum(5,:), helicity(5), Acurrent2)
        endif
      else
        call FFV(id(5), four_momentum(5,:), helicity(5), id(4), four_momentum(4,:), helicity(4), Vcurrent2)
        if((id(4)+id(5)).eq.0)then
          call FFA(id(5), four_momentum(5,:), helicity(5), id(4), four_momentum(4,:), helicity(4), Acurrent2)
        endif
      endif



!WH
      if((beam_id(1)+beam_id(2)).ne.0)then
        if((beam_id(1)*beam_h(1)).le.0d0)then
          current1=Vcurrent1*gFFW*CKM(abs(beam_id(1)),abs(beam_id(2)))
        else
          current1=0d0
        endif
        current2=Vcurrent2*gFFW*CKM(abs(id(4)),abs(id(5)))

!ZH
      else if((abs(beam_id(1)).eq.11).or.(abs(beam_id(1)).eq.13))then
!        print *, beam_id, beam_h
!e+ e- Z vertex for incoming states
        if((beam_id(1)*beam_h(1)).gt.0d0)then
          current1=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent1 -(0.5d0*T3lR)*Acurrent1
        else
          current1=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent1 -(0.5d0*T3lL)*Acurrent1
        endif
        current1=current1*gFFZ
!u u~ Z vertex for incoming states
      else if((abs(beam_id(1)).eq.2).or.(abs(beam_id(1)).eq.4))then
        if((beam_id(1)*beam_h(1)).gt.0d0)then
          current1=(0.5d0*T3uR - QuR*sitW**2) *Vcurrent1 -(0.5d0*T3uR)*Acurrent1
        else
          current1=(0.5d0*T3uL - QuL*sitW**2) *Vcurrent1 -(0.5d0*T3uL)*Acurrent1
        endif
        current1=current1*gFFZ
!d d~ Z vertex for incoming states
      else if((abs(beam_id(1)).eq.1).or.(abs(beam_id(1)).eq.3).or.(abs(beam_id(1)).eq.5))then
        if((beam_id(1)*beam_h(1)).gt.0d0)then
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
      if((id(4)+id(5)).eq.0)then
!l+ l- Z vertex for final state
        if((abs(id(4)).eq.11).or.(abs(id(4)).eq.13).or.(abs(id(4)).eq.15))then
          if((id(4)*helicity(4)).gt.0d0)then
            current2=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent2 -(0.5d0*T3lR)*Acurrent2
          else
            current2=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent2 -(0.5d0*T3lL)*Acurrent2
          endif
          current2=current2*gFFZ
!u u~ Z vertex for final state
        else if((abs(id(4)).eq.2).or.(abs(id(4)).eq.4))then
          if((id(4)*helicity(4)).gt.0d0)then
            current2=(0.5d0*T3uR - QuR*sitW**2) *Vcurrent2 -(0.5d0*T3uR)*Acurrent2
          else
            current2=(0.5d0*T3uL - QuL*sitW**2) *Vcurrent2 -(0.5d0*T3uL)*Acurrent2
          endif
          current2=current2*gFFZ

!d d~ Z vertex for final state
        else if((abs(id(4)).eq.1).or.(abs(id(4)).eq.3).or.(abs(id(4)).eq.5))then
          if((id(4)*helicity(4)).gt.0d0)then
            current2=(0.5d0*T3dR - QdR*sitW**2) *Vcurrent2 -(0.5d0*T3dR)*Acurrent2
          else
            current2=(0.5d0*T3dL - QdL*sitW**2) *Vcurrent2 -(0.5d0*T3dL)*Acurrent2
          endif
          current2=current2*gFFZ

!nu nu~ Z vertex for final state        
        else if((abs(id(4)).eq.12).or.(abs(id(4)).eq.14).or.(abs(id(4)).eq.16))then
          current2=(0.5d0*T3nL - QnL*sitW**2) *Vcurrent2 -(0.5d0*T3nL)*Acurrent2
          current2=current2*gFFZ

        else
        current2=0d0
        print *, "invalid final state"!, id(4:5)
        stop
        endif

      endif
!f+ f- Z vertex for final states END

      call POLARIZATION(four_momentum(1,:), POL1)
      call POLARIZATION(four_momentum(2,:), POL2)
   

!ZZX vertex
      gVVS1 = ghz1*(mass(1,1)**2) + qq * ( 2d0*ghz2 + ghz3*qq/(Lambda*1d2) )

      gVVS2 = -( 2d0*ghz2 + ghz3*qq/(Lambda*1d2) )

      gVVP = -2d0*ghz4

      VVX0 = 0d0
      if(gVVS1.ne.0d0)then
        call VVS1(g_mu_nu)
        VVX0 = VVX0 + gVVS1*g_mu_nu
      endif

      if(gVVS2.ne.0d0)then
        call VVS2(four_momentum(3,:),four_momentum(3,:),pp)
        VVX0 = VVX0 + gVVS2*pp
      endif

      if(gVVP.ne.0d0)then
        call VVP(-four_momentum(1,:),four_momentum(2,:),epp)
        VVX0 = VVX0 + gVVP*epp
      endif

      VVX0 = (0d0,1d0)/(vev*1d2)*VVX0


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

      MATRIXELEMENT0=MATRIXELEMENT0 *PROP1*PROP2*PROP3 &
      *(gFFS*FFS(id(6), four_momentum(6,:), helicity(6), id(7), four_momentum(7,:), helicity(7)) &
       +gFFP*FFP(id(6), four_momentum(6,:), helicity(6), id(7), four_momentum(7,:), helicity(7)))&
        *(0d0,-1d0)*b_Yukawa/(vev*1d2)


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
!     double precision Pi
      double precision Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      double precision sincos(4), vector(4), abs3p, phi
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
        phi=datan(vector(3)/vector(2))
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
!     double precision Pi
      double precision Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision

      double complex p1(4), p2(4)
      double complex epp(4,4)
!      double precision ANTISYMMETRIC
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

      double precision function ANTISYMMETRIC(i,j,k,l)

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
!     double precision Pi
      double precision Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      double complex epep(4,4),emem(4,4),epe0(4,4),eme0(4,4),e0e0(4,4)
      double complex epem(4,4),e0ep(4,4),e0em(4,4),emep(4,4)
      double complex POL(3,4), T_mu_nu(5,4,4)
      
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
      double complex p1(4), p2(4)
      double complex pp(4,4)
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
!     double precision Pi
      double precision Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      double complex epep(4,4),emem(4,4),epe0(4,4),eme0(4,4),e0e0(4,4)
      double complex epem(4,4),e0ep(4,4),e0em(4,4),emep(4,4)
      double complex POL(3,4), T_mu_nu(5,4,4)
      
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
      double complex p1(4), p2(4)
      double complex pp(4,4)
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

      double complex function COVARIANT_VECTOR(p,mu)

      implicit none

      double complex p(4)
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

      double complex function FFP(pdg_code1, p1, h1, pdg_code2, p2, h2)

      implicit none
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      double precision p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      double precision sqrt_pp1Dpp2

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
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      double precision p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      double precision sqrt_pp1Dpp2, sqrt_pp1Xpp2
      double complex Acurrent(4)
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

      double complex function FFS(pdg_code1, p1, h1, pdg_code2, p2, h2)

      implicit none
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      double precision p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      double precision sqrt_pp1Dpp2

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
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
      
      double precision p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      double precision sqrt_pp1Dpp2, sqrt_pp1Xpp2
      double complex Vcurrent(4)
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

      double complex function INNER(p1,p2)

      implicit none
      double complex p1(4), p2(4)

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

      double precision vector(4), boost(4) 
      double precision lambda(4,4), vector_copy(4)
      double precision beta(2:4), beta_sq, gamma
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
      double precision function KRONECKER_DELTA(i,j)
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

      double precision function METRIC(mu,nu)

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
      double precision p(4), sincos(4), inv_mass, abs3p
      double complex POL(3,4)
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
      double precision p(4), sincos(4), inv_mass, abs3p
      double complex POL(2,4)
      
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
      double precision p(4), sincos(4), inv_mass, abs3p
      double complex POL(3,4)
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

      double complex function PROPAGATOR(inv_mass, mass, width)
      implicit none

      double precision inv_mass, mass, width

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
      double precision p1(4), p2(4)
      double complex epp(4,4)
!      double precision ANTISYMMETRIC
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
      double complex g_mu_nu(4,4)

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
      double precision p1(4), p2(4)
      double complex pp(4,4)
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

























end module ModVHiggs
!!--YaofuZhou-----------------------------------------