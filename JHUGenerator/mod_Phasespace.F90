MODULE ModPhasespace 
!    This implementation closely follows the notation in Stefan Kallweit's Master Thesis
!
!    NOTE: (1) global factors of pi are excluded, i.e. multiply by (2*pi)^(4-3*N) for N-particle phase space (see PSNormi below)
!          (2) no external dependencies required



  public :: s_channel_prop_decay
  public :: t_channel_prop_decay
  public :: s_channel_propagator
  public :: s_channel_decay

  public  :: k_s
  public  :: k_t
  public  :: k_l
  public  :: k_BreitWigner 

  public :: get_minmax_s
  public :: get_minmax_t

    
  private :: rotate3D_phi_theta
  private :: boost_from_CMS_to_RefMom
  private :: sqrt_lambda
  private :: g_s
  private :: g_s_BreitWigner
  private :: g_d
  private :: g_p
  private :: h
  private :: h_l
  private :: h_BreitWigner 



!   private :: remap_s_channel_propagator
!   private :: generate_s_channel_phasespace
!   private :: generate_t_channel_phasespace

!   public :: Propagator_S_channel
!   public :: DecayPS_S_channel
!   public :: DecayPS_T_channel
  
!   type :: Propagator_S_channel
!       real(8) :: RandomVar
!       real(8) :: PropMass_sq
!       real(8) :: Width
!       real(8) :: Power = 2d0
!       real(8) :: InvMinMax_sq(1:2) = (/0d0,14000d0**2/)      
!       contains
!         procedure :: Remap => remap_s_channel_propagator
!   end type


!   type :: DecayPS_S_channel
!       real(8) :: RandomVars(1:2)
!       real(8) :: p0(1:4)
!       real(8) :: Mass1_sq
!       real(8) :: Mass2_sq
!       contains
!         procedure :: generate => generate_s_channel_phasespace
!   end type


!   type :: DecayPS_T_channel
!       real(8) :: RandomVars(1:2)
!       real(8) :: pa(1:4)
!       real(8) :: pb(1:4)
!       real(8) :: Mass1_sq
!       real(8) :: Mass2_sq
!       real(8) :: PropMass_sq
!       real(8) :: Power = 2d0
!       real(8) :: InvMinMax_sq(1:2) = (/0d0,14000d0**2/)      
!       contains
!         procedure :: generate => generate_t_channel_phasespace
!   end type


  real(8), private, parameter :: pi = 3.1415926535897932384626433d0

  real(8), public, parameter :: PSNorm2 =  (2*pi)**(4-3*2)
  real(8), public, parameter :: PSNorm3 =  (2*pi)**(4-3*3)
  real(8), public, parameter :: PSNorm4 =  (2*pi)**(4-3*4)
  real(8), public, parameter :: PSNorm5 =  (2*pi)**(4-3*5)
  real(8), public, parameter :: PSNorm6 =  (2*pi)**(4-3*6)
  real(8), public, parameter :: PSNorm7 =  (2*pi)**(4-3*7)
  real(8), public, parameter :: PSNorm8 =  (2*pi)**(4-3*8)
  real(8), public, parameter :: PSNorm9 =  (2*pi)**(4-3*9)


  
  
 CONTAINS
 
 
 
 
  ! note: if partx(1)=partx(3)=partx(4) (i.e. sqrt(smin)=sqrt(smax)=mass) then we assume a stable internal particle with narrow-width jacobian
  ! note: if 0=partx(3)=partx(4) (i.e. sqrt(smin)=sqrt(smax)=0) then we assume a stable external particle
  ! note: the number of xRnd's has to be chosen correctly
  FUNCTION s_channel_prop_decay(p0,part1,part2,xRnd,Mom1,Mom2,PropPower) 
  implicit none
  real(8) :: s_channel_prop_decay
  real(8) :: p0(:),xRnd(:),Mom1(1:4),Mom2(1:4),Jac1,Jac2,Jac3,Minvsq_1,Minvsq_2,Power
  real(8) :: part1(1:4),part2(1:4) ! 1=mass, 2=width, 3=sqrt(smin), 4=sqrt(smax)
  real(8),optional :: PropPower
  integer :: iRnd
  
      if( present(PropPower) ) then
         Power = PropPower
      else
         Power = 2d0
      endif
      iRnd = 1

      if( part1(3).eq.0d0 .and. part1(4).eq.0d0 ) then! on-shell decay
          Minvsq_1 = part1(1)**2
          Jac1 = 1d0
      elseif( part1(1).eq.part1(3) .and. part1(1).eq.part1(4) ) then! narrow-width decay
          Minvsq_1 = part1(1)**2
          Jac1 = pi/(part1(1)*part1(2))
      else! propagator
          Jac1 = s_channel_propagator( part1(1)**2,part1(2), part1(3)**2,part1(4)**2, xRnd(iRnd),Minvsq_1,Power )
          iRnd = iRnd+1
      endif
      

      if( part2(3).eq.0d0 .and. part2(4).eq.0d0 ) then! on-shell decay
          Minvsq_2 = part2(1)**2
          Jac2 = 1d0
      elseif( part2(1).eq.part2(3) .and. part2(1).eq.part2(4) ) then! narrow-width decay
          Minvsq_2 = part2(1)**2
          Jac2 = pi/(part2(1)*part2(2))
      else! propagator
          Jac2 = s_channel_propagator( part2(1)**2,part2(2), part2(3)**2,part2(4)**2, xRnd(iRnd),Minvsq_2,Power )
          iRnd = iRnd+1          
      endif
      
      
      Jac3 = s_channel_decay( p0,Minvsq_1,Minvsq_2,xRnd(iRnd:iRnd+1),Mom1,Mom2 )       
      s_channel_prop_decay = Jac1*Jac2*Jac3
 
  RETURN
  END FUNCTION

 
 
 
 
! t-channel phase space
  FUNCTION t_channel_prop_decay(pa,pb,PropMass_sq,Mass1_sq,Mass2_sq,xRnd,Mom1,Mom2,PropPower) 
  implicit none
  real(8) :: t_channel_prop_decay
  real(8) :: pa(1:4),pb(1:4),p0(1:4),xRnd(:),Mom1(1:4),Mom2(1:4)
  real(8) :: phi,Theta_hat,CosTheta,p0_sq,Mandelstam_S,Mandelstam_T,pa_sq,pb_sq,Eab,E1,p1z,T_MinMax_sq(1:2)
  real(8) :: Mass1_sq,Mass2_sq,PropMass_sq,phi_hat,CosTheta_hat,pa_hat(1:4),p_boost(1:4),Power
  real(8),optional :: PropPower
  integer :: iRnd
  
       if( present(PropPower) ) then
         Power = PropPower
       else
         Power = 2d0
       endif
       iRnd = 1
       
       p0(1:4) = pa(1:4) + pb(1:4)
       pa_sq = pa(1)**2 - pa(2)**2 - pa(3)**2 - pa(4)**2
       pb_sq = pb(1)**2 - pb(2)**2 - pb(3)**2 - pb(4)**2
       Mandelstam_S = p0(1)*p0(1) -p0(2)*p0(2) -p0(3)*p0(3) -p0(4)*p0(4)
       call get_minmax_t(Mandelstam_S,pa_sq,pb_sq,Mass1_sq,Mass2_sq,T_MinMax_sq)       
!        Mandelstam_T = -h(xRnd(2), -PropMass_sq, Power, -T_MinMax_sq(1), -T_MinMax_sq(2))
!        t_channel_prop_decay = 1d0/g_p( Mandelstam_S, pa_sq, pb_sq, Mandelstam_T, PropMass_sq, power, T_MinMax_sq(1), T_MinMax_sq(2))
       t_channel_prop_decay = k_t(xRnd(2),Mandelstam_S, pa_sq, pb_sq, PropMass_sq, Power, T_MinMax_sq(1), T_MinMax_sq(2),Mandelstam_T)
       
       E1 = ( Mandelstam_S + Mass1_sq - Mass2_sq )/2d0/dsqrt(Mandelstam_S)
       p1z= sqrt_lambda(Mandelstam_S,Mass1_sq,Mass2_sq)/2d0/dsqrt(Mandelstam_S)
       Mom1(1:4) = (/ E1,0d0,0d0,p1z /)

       phi   = 2d0*pi * xRnd(1)
       CosTheta = 1d0/sqrt_lambda(Mandelstam_S,Mass1_sq,Mass2_sq)/sqrt_lambda(Mandelstam_S,pa_sq,pb_sq)    &
                * ( (Mandelstam_S+Mass1_sq-Mass2_sq)*(Mandelstam_S+pa_sq-pb_sq) + 2d0*Mandelstam_S*(Mandelstam_T-Mass1_sq-pa_sq) )       
       call rotate3D_phi_theta(phi,dacos(CosTheta),Mom1(:))
       
       Eab = Mandelstam_S/dsqrt(Mandelstam_S)
       Mom2(1:4) = (/Eab,0d0,0d0,0d0/) - Mom1(1:4)


       
       pa_hat(1:4) = pa(1:4)
       p_boost(1:4) = (/ pa(1)+pb(1), -(pa(2:4)+pb(2:4)) /)       
       call boost_from_CMS_to_RefMom(p_boost,pa_hat)
       CosTheta_hat = pa_hat(4)/dsqrt( pa_hat(2)**2+pa_hat(3)**2+pa_hat(4)**2 )
       Theta_hat = dacos(CosTheta_hat)       
       phi_hat = datan( pa_hat(3)/(pa_hat(2)+1d-15) )
       if( pa_hat(2).lt.0d0 ) phi_hat = phi_hat + pi

       
! !checks
! print *,""
! E1 = ( Mandelstam_S + pa_sq - pb_sq )/2d0/dsqrt(Mandelstam_S)
! p1z= sqrt_lambda(Mandelstam_S,pa_sq,pb_sq)/2d0/dsqrt(Mandelstam_S)
! pa_hat(1:4) = (/ E1,0d0,0d0,p1z /)! = pa_star
! call rotate3D_phi_theta(-phi_hat,dacos(CosTheta_hat),pa_hat)
! print *, "R * pa_star=",pa_hat
! 
! p_boost(1:4) = (/ pa(1)+pb(1), -(pa(2:4)+pb(2:4)) /)
! call boost_from_CMS_to_RefMom(p_boost,pa(:))
! print *, "B * pa=",pa
! print *,""
! !checks
       
       call rotate3D_phi_theta(-phi_hat,Theta_hat,Mom1(:),Mom2(:))
!        call rotate3D_phi_theta(-phi_hat,Theta_hat,Mom2(:))
       p_boost(1:4) = pa(1:4)+pb(1:4)
       call boost_from_CMS_to_RefMom(p_boost,Mom1(:),Mom2(:))
!        call boost_from_CMS_to_RefMom(p_boost,Mom2(:))
       

  RETURN
  END FUNCTION


 
 
 

 
  FUNCTION s_channel_propagator(PropMass_sq,Width,sMin,sMax,xRnd,InvMass_sq,PropPower) 
  implicit none
  real(8) :: s_channel_propagator
  real(8) :: PropMass_sq,Width,sMin,sMax,xRnd,InvMass_sq,Power
  real(8), optional :: PropPower

      if( present(PropPower) ) then
         Power = PropPower
      else
         Power = 2d0
      endif

      if( Width.lt.1d-15 ) then
           InvMass_sq = h(xRnd, PropMass_sq, Power, sMin, sMax)
           s_channel_propagator = 1d0/g_s(InvMass_sq, PropMass_sq, Power, sMin, sMax)
print *, "check old",InvMass_sq,s_channel_propagator
           s_channel_propagator = k_s(xRnd, PropMass_sq, Power, sMin, sMax,InvMass_sq)
print *, "check new",InvMass_sq,s_channel_propagator
pause
                      
      else
           if( dabs(sMin + sMax -2d0*PropMass_sq).lt.1d-10 .and.PropMass_sq.gt.1d-14 ) then! this is the narrow-width mapping (i.e. no integration)
              InvMass_sq = PropMass_sq
              s_channel_propagator = pi/(dsqrt(PropMass_sq)*Width)
           else ! this is the normal s-channel mapping  
              if( sMax.gt.PropMass_sq-Width**2 ) then
!                   InvMass_sq = h_BreitWigner(xRnd, PropMass_sq, Width, Power, sMin, sMax)
!                   s_channel_propagator = 1d0/g_s_BreitWigner(InvMass_sq, PropMass_sq, Width, Power, sMin, sMax)
                  s_channel_propagator = k_BreitWigner(xRnd, PropMass_sq, Width, Power, sMin, sMax,InvMass_sq)
              else
!                   InvMass_sq = sMin + (sMax-sMin) * xRnd
!                   s_channel_propagator = sMax-sMin
                  s_channel_propagator = k_l(xRnd,sMin,sMax,InvMass_sq)
              endif
           endif
      endif

  RETURN
  END FUNCTION



  
  
  

! s-channel phase space
  FUNCTION s_channel_decay( p0,Mass1_sq,Mass2_sq,xRnd,Mom1,Mom2  ) 
  implicit none
  real(8) :: s_channel_decay
  real(8) :: p0(1:4),Mass1_sq,Mass2_sq,xRnd(1:2),Mom1(1:4),Mom2(1:4)
  real(8) :: E1,p1z,phi,theta,Jac,Mandelstam_S

        phi   = 2d0*pi * xRnd(1)
        theta = 1d0*pi * xRnd(2)
        Mandelstam_S = p0(1)**2 - p0(2)**2 - p0(3)**2 - p0(4)**2 

        E1 = ( Mandelstam_S + Mass1_sq - Mass2_sq )/2d0/dsqrt(Mandelstam_S)
        p1z= sqrt_lambda(Mandelstam_S,Mass1_sq,Mass2_sq)/2d0/dsqrt(Mandelstam_S)
        Mom1(1:4) = (/ E1,0d0,0d0,p1z /)
!         call rotate3D_phi_theta(phi,theta,Mom1)
        
        Mom2(1)   = dsqrt(Mandelstam_S) - Mom1(1)
        Mom2(2:4) = - Mom1(2:4)
        call boost_from_CMS_to_RefMom(p0,Mom1,Mom2)
!         call boost_from_CMS_to_RefMom(p0,Mom2)
 
        s_channel_decay = 1d0/g_d(Mandelstam_S, Mass1_sq,Mass2_sq)  *  dSin(pi*xRnd(2))

  RETURN
  END FUNCTION
 

 
 


! ----- transformations ---------

  SUBROUTINE rotate3D_phi_theta(phi,theta,mom,mom2)
  implicit none
  real(8) :: phi,theta,mom(1:4)
  real(8), optional :: Mom2(1:4)
  real(8) :: cos_phi,sin_phi,cos_theta,sin_theta,mom_tmp(2:4),mom_tmp2(2:4)
  
      
      cos_phi = dCos(phi)
      sin_phi = dSin(phi)
      cos_theta = dCos(theta)
      sin_theta = dSin(theta)

      mom_tmp(2:4) = mom(2:4)      
      mom(2) = mom_tmp(2) * cos_phi * cos_theta + mom_tmp(3) * sin_phi + mom_tmp(4) * cos_phi * sin_theta
      mom(3) = mom_tmp(3) * cos_phi - mom_tmp(2) * cos_theta * sin_phi - mom_tmp(4) * sin_phi * sin_theta
      mom(4) = mom_tmp(4) * cos_theta - mom_tmp(2) * sin_theta
      
      if( present(mom2) ) then
         mom_tmp2(2:4) = mom2(2:4)
         mom2(2) = mom_tmp2(2) * cos_phi * cos_theta + mom_tmp2(3) * sin_phi + mom_tmp2(4) * cos_phi * sin_theta
         mom2(3) = mom_tmp2(3) * cos_phi - mom_tmp2(2) * cos_theta * sin_phi - mom_tmp2(4) * sin_phi * sin_theta
         mom2(4) = mom_tmp2(4) * cos_theta - mom_tmp2(2) * sin_theta
      endif
      
  RETURN
  END SUBROUTINE


  
! careful, this is the inverse of Kallweit's boost
  SUBROUTINE boost_from_CMS_to_RefMom(RefMom,BoostMom,BoostMom2)
  implicit none
  real(8) :: BoostMom(1:4),RefMom(1:4)
  real(8), optional :: BoostMom2(1:4)
  real(8) :: EuklidSP,spacialFact,refMass

      refMass = dsqrt( RefMom(1)**2 - RefMom(2)**2 - RefMom(3)**2 - RefMom(4)**2 )

      EuklidSP = refMom(2)*boostMom(2) + refMom(3)*boostMom(3) + refMom(4)*boostMom(4)
      spacialFact = (EuklidSP/(refMom(1) + refMass) + boostMom(1))/refMass
      boostMom(1) = (refMom(1)*boostMom(1) + EuklidSP)/refMass
      boostMom(2) = boostMom(2) + refMom(2)*spacialFact
      boostMom(3) = boostMom(3) + refMom(3)*spacialFact
      boostMom(4) = boostMom(4) + refMom(4)*spacialFact  
 
      if( present(BoostMom2) ) then
          EuklidSP = refMom(2)*boostMom2(2) + refMom(3)*boostMom2(3) + refMom(4)*boostMom2(4)
          spacialFact = (EuklidSP/(refMom(1) + refMass) + boostMom2(1))/refMass
          boostMom2(1) = (refMom(1)*boostMom2(1) + EuklidSP)/refMass
          boostMom2(2) = boostMom2(2) + refMom(2)*spacialFact
          boostMom2(3) = boostMom2(3) + refMom(3)*spacialFact
          boostMom2(4) = boostMom2(4) + refMom(4)*spacialFact        
      endif
 
  RETURN
  END SUBROUTINE







! ----- auxilary functions ---------

  FUNCTION sqrt_lambda(x,y,z)
  implicit none
  real(8) :: sqrt_lambda,x,y,z
  
        sqrt_lambda = dsqrt(x**2 + y**2 + z**2 - 2d0*x*y - 2d0*x*z - 2d0*y*z)

  RETURN
  END FUNCTION



  FUNCTION g_s(s,m_sq,nu,smin,smax)! inverse Jacobian for s-channel mapping
  implicit none
  real(8) :: g_s,s,m_sq,nu,smin,smax
  
        if( nu.ne.1d0 ) then
          g_s =  (1d0-nu)/( (smax-m_sq)**(1d0-nu) - (smin-m_sq)**(1d0-nu) )/( s - m_sq )**nu
        else
          g_s = 1d0/( dlog(smax-m_sq) - dlog(smin-m_sq) )/( s - m_sq )
        endif
    
  RETURN
  END FUNCTION




  FUNCTION g_l(smin,smax)! inverse Jacobian for linear mapping
  implicit none
  real(8) :: g_l,smin,smax
  
      g_l =  1d0/(smax - smin)
          
  RETURN
  END FUNCTION
  
  

  FUNCTION g_s_BreitWigner(s,m_sq,gamma,nu,smin,smax)! inverse Jacobian for Breit-Wigner mapping
  implicit none
  real(8) :: g_s_BreitWigner,s,m_sq,gamma,nu,smin,smax,y1,y2,m

        m = dsqrt(m_sq)
        y1 = datan( (smin-m_sq)/m/gamma )
        y2 = datan( (smax-m_sq)/m/gamma )
        g_s_BreitWigner = m*gamma /(y2-y1) /( (s-m_sq)**2 + m_sq*gamma**2 )
        
  RETURN
  END FUNCTION



  FUNCTION g_d(sab,m1_sq,m2_sq)! inverse Jacobian for 1-->2 phase space
  implicit none
  real(8) :: g_d,sab,m1_sq,m2_sq
  
        g_d = 2d0*sab/sqrt_lambda(sab,m1_sq,m2_sq)/pi 
    
  RETURN
  END FUNCTION



  FUNCTION g_p(sab,pa_sq,pb_sq,t,m_sq,nu,tmin,tmax)! inverse Jacobian for 2-->2 phase space
  implicit none
  real(8) :: g_p,sab,pa_sq,pb_sq,t,m_sq,nu,tmin,tmax
  
        g_p = - 2d0/pi * sqrt_lambda(sab,pa_sq,pb_sq) * g_s(-t,-m_sq,nu,-tmin,-tmax)
    
  RETURN
  END FUNCTION



  FUNCTION h(RandomVar,m_sq,nu,smin,smax)! inv. mass mapping for s- and t-channel propagator
  implicit none
  real(8) :: h,RandomVar,m_sq,nu,smin,smax
  
        if( nu.ne.1d0 ) then
           h =   ( RandomVar*(smax-m_sq)**(1d0-nu) + (1d0-RandomVar)*(smin-m_sq)**(1d0-nu) )**(1d0/(1d0-nu))   +  m_sq
        else
           h = dexp( RandomVar*dlog(smax-m_sq) - (1d0-RandomVar)*dlog(smin-m_sq)  )    +  m_sq 
        endif

  RETURN
  END FUNCTION


  
  FUNCTION h_l(RandomVar,smin,smax)! linear inv. mass mapping
  implicit none
  real(8) :: h_l,RandomVar,smin,smax
  
         h_l = RandomVar*smax + (1d0-RandomVar)*smin

  RETURN
  END FUNCTION



  FUNCTION h_BreitWigner(RandomVar,m_sq,gamma,nu,smin,smax)! inv. mass mapping for BreitWigner propagator
  implicit none
  real(8) :: h_BreitWigner,RandomVar,m_sq,gamma,nu,smin,smax,m,y1,y2

        m  = dsqrt(m_sq)
        y1 = datan( (smin-m_sq)/m/gamma )
        y2 = datan( (smax-m_sq)/m/gamma )  
        
        h_BreitWigner = dsqrt(m_sq)*gamma * dtan(y1 + (y2-y1)*RandomVar) + m_sq

  RETURN
  END FUNCTION




  
  
  
! the k functions below are equal to the h functions above, only they return the invariant as argument and the g-jacobian as function return value

! this is h(RandomVar,m_sq,nu,smin,smax), i.e. h(..) for s-channel
  FUNCTION k_s(RandomVar,m_sq,nu,smin,smax,InvMass_sq)
  implicit none
  real(8) :: k_s,RandomVar,m_sq,nu,smin,smax,InvMass_sq
  
  
        if( nu.ne.1d0 ) then
           InvMass_sq =   ( RandomVar*(smax-m_sq)**(1d0-nu) + (1d0-RandomVar)*(smin-m_sq)**(1d0-nu) )**(1d0/(1d0-nu))  +  m_sq
           k_s =  ( (smax-m_sq)**(1d0-nu) - (smin-m_sq)**(1d0-nu) )*(( InvMass_sq - m_sq )**nu)/(1d0-nu) ! =  1/g_s
        else
           InvMass_sq = dexp( RandomVar*dlog(smax-m_sq) - (1d0-RandomVar)*dlog(smin-m_sq)  )  +  m_sq 
           k_s = ( dlog(smax-m_sq) - dlog(smin-m_sq) )*( InvMass_sq - m_sq ) ! =  1/g_s
        endif
        

  RETURN
  END FUNCTION

! this is -h(RandomVar,-m_sq,nu,-smin,-smax), i.e. h(..) for t-channel
  FUNCTION k_t(RandomVar,sab,pa_sq,pb_sq,m_sq,nu,smin,smax,InvMass_sq)
  implicit none
  real(8) :: k_t,RandomVar,m_sq,nu,smin,smax,InvMass_sq,sab,pa_sq,pb_sq
  
  
        if( nu.ne.1d0 ) then
           InvMass_sq =  - ( ( RandomVar*(-smax+m_sq)**(1d0-nu) + (1d0-RandomVar)*(-smin+m_sq)**(1d0-nu) )**(1d0/(1d0-nu))  -  m_sq   )
        else
           InvMass_sq = - ( dexp( RandomVar*dlog(-smax+m_sq) - (1d0-RandomVar)*dlog(-smin+m_sq)  )  -  m_sq )
        endif

        k_t =  -pi/2d0/sqrt_lambda(sab,pa_sq,pb_sq)/g_s(-InvMass_sq,-m_sq,nu,-smin,-smax)
                
  RETURN
  END FUNCTION


  
  FUNCTION k_l(RandomVar,smin,smax,InvMass_sq)
  implicit none
  real(8) :: k_l,RandomVar,smin,smax,InvMass_sq
  
         InvMass_sq = RandomVar*smax + (1d0-RandomVar)*smin
         k_l = (smax - smin) ! = 1/g_l(...)

  RETURN
  END FUNCTION



  FUNCTION k_BreitWigner(RandomVar,m_sq,gamma,nu,smin,smax,InvMass_sq)
  implicit none
  real(8) :: k_BreitWigner,RandomVar,m_sq,gamma,nu,smin,smax,m,y1,y2,InvMass_sq

        m  = dsqrt(m_sq)
        y1 = datan( (smin-m_sq)/m/gamma )
        y2 = datan( (smax-m_sq)/m/gamma )  
        InvMass_sq = dsqrt(m_sq)*gamma * dtan(y1 + (y2-y1)*RandomVar) + m_sq

        k_BreitWigner = (y2-y1)*( (InvMass_sq-m_sq)**2 + m_sq*gamma**2 )/(m*gamma)!  = 1/g_s_BreitWigner(...)
                
  RETURN
  END FUNCTION


  
  SUBROUTINE get_minmax_s(shat,mass1_sq,mass2_sq,sothers_min,minmax)
  implicit none
  real(8) :: shat,mass1_sq,mass2_sq,sothers_min,minmax(1:2)
     
     minmax(1) = (dsqrt(mass1_sq) + dsqrt(mass2_sq))**2   ! =min
     minmax(2) = (dsqrt(shat) - dsqrt(sothers_min))**2    ! =max
  
  RETURN
  END SUBROUTINE
  

  
  SUBROUTINE get_minmax_t(sab,pa_sq,pb_sq,ka_sq,kb_sq,minmax)
  implicit none
  real(8) :: sab,pa_sq,pb_sq,ka_sq,kb_sq,minmax(1:2),tmp
     
     tmp = ka_sq + pa_sq - (sab+ka_sq-kb_sq)*(sab+pa_sq-pb_sq)/(2d0*sab)
     minmax(1) = tmp - sqrt_lambda(sab,ka_sq,kb_sq)*sqrt_lambda(sab,pa_sq,pb_sq)/(2d0*sab)  ! =min
     minmax(2) = tmp + sqrt_lambda(sab,ka_sq,kb_sq)*sqrt_lambda(sab,pa_sq,pb_sq)/(2d0*sab)  ! =max
  
  RETURN
  END SUBROUTINE
  
  
  
  
!--------- private section ------------
 


! ! s-channel invariant mapping
!   SUBROUTINE remap_s_channel_propagator(this,InvMass_sq,Jac) 
!   implicit none
!   class(Propagator_S_channel) :: this
!   real(8) :: InvMass_sq,Jac
! 
!       if( this%Width.lt.1d-14 ) then
!            InvMass_sq = h(this%RandomVar, this%PropMass_sq, this%Power, this%InvMinMax_sq(1), this%InvMinMax_sq(2))
!            Jac = g_s(InvMass_sq, this%PropMass_sq, this%Power, this%InvMinMax_sq(1), this%InvMinMax_sq(2))
!       else
!            if( dabs(this%InvMinMax_sq(1) + this%InvMinMax_sq(2) -2d0*this%PropMass_sq).lt.1d-10 ) then! this is the narrow-width mapping (i.e. no integration)
!               InvMass_sq = this%PropMass_sq
!               Jac = pi/(dsqrt(this%PropMass_sq)*this%Width)
!            else ! this is the normal s-channel mapping     
!               InvMass_sq = h_BreitWigner(this%RandomVar, this%PropMass_sq, this%Width, this%Power, this%InvMinMax_sq(1), this%InvMinMax_sq(2))
!               Jac = g_s_BreitWigner(InvMass_sq, this%PropMass_sq, this%Width, this%Power, this%InvMinMax_sq(1), this%InvMinMax_sq(2))
!            endif
!       endif
! 
!   RETURN
!   END SUBROUTINE







! ! s-channel phase space
!   SUBROUTINE generate_s_channel_phasespace(this,Mom,Jac) 
!   implicit none
!   class(DecayPS_S_channel) :: this
!   real(8) :: Mom(1:4,1:2),E1,p1z,phi,theta,Jac,Mandelstam_S
! 
!         phi   = 2d0*pi * this%RandomVars(1)
!         theta = 1d0*pi * this%RandomVars(2)
!         Mandelstam_S = this%p0(1)**2 - this%p0(2)**2 - this%p0(3)**2 - this%p0(4)**2 
! 
!         E1 = ( Mandelstam_S + this%Mass1_sq - this%Mass2_sq )/2d0/dsqrt(Mandelstam_S)
!         p1z= sqrt_lambda(Mandelstam_S,this%Mass1_sq,this%Mass2_sq)/2d0/dsqrt(Mandelstam_S)
!         Mom(1:4,1) = (/ E1,0d0,0d0,p1z /)
!         call rotate3D_phi_theta(phi,theta,Mom(:,1))
!         
!         Mom(1:4,2) = this%p0(1:4) - Mom(1:4,1)
!         
!         call boost_from_CMS_to_RefMom(Mom(:,1),this%p0)
!         call boost_from_CMS_to_RefMom(Mom(:,2),this%p0)
! 
!         Jac = 1d0/g_d(Mandelstam_S, this%Mass1_sq,this%Mass2_sq)  *  dSin(pi*this%RandomVars(2))
! 
!   RETURN
!   END SUBROUTINE






! ! t-channel phase space
!   SUBROUTINE generate_t_channel_phasespace(this,Mom,Jac) 
!   implicit none
!   class(DecayPS_T_channel) :: this
!   real(8) :: Mom(1:4,1:2),p0(1:4)
!   real(8) :: phi,CosTheta,p0_sq,Jac,Mandelstam_S,Mandelstam_T,pa_sq,pb_sq,Eab,E1,p1z
!   real(8) :: phi_hat,CosTheta_hat,pa_hat(1:4),p_boost(1:4)
!   
!        p0(1:4) = this%pa(1:4) + this%pb(1:4)
!        pa_sq = this%pa(1)**2 - this%pa(2)**2 - this%pa(3)**2 - this%pa(4)**2
!        pb_sq = this%pb(1)**2 - this%pb(2)**2 - this%pb(3)**2 - this%pb(4)**2
!        Mandelstam_S = p0(1)*p0(1) -p0(2)*p0(2) -p0(3)*p0(3) -p0(4)*p0(4)
!        Mandelstam_T = -h(this%RandomVars(2), -this%PropMass_sq, this%Power, -this%InvMinMax_sq(1), -this%InvMinMax_sq(2))!  -m is correct ? 
! 
!        E1 = ( Mandelstam_S + this%Mass1_sq - this%Mass2_sq )/2d0/dsqrt(Mandelstam_S)
!        p1z= sqrt_lambda(Mandelstam_S,this%Mass1_sq,this%Mass2_sq)/2d0/dsqrt(Mandelstam_S)
!        phi   = 2d0*pi * this%RandomVars(1)
!        CosTheta = 1d0/sqrt_lambda(Mandelstam_S,this%Mass1_sq,this%Mass2_sq)/sqrt_lambda(Mandelstam_S,pa_sq,pb_sq)    &
!                 * ( (Mandelstam_S+this%Mass1_sq-this%Mass2_sq)*(Mandelstam_S+pa_sq-pb_sq) + 2d0*Mandelstam_S*(Mandelstam_T-this%Mass1_sq-pa_sq) )
!        Mom(1:4,1) = (/ E1,0d0,0d0,p1z /)
!        call rotate3D_phi_theta(phi,dacos(CosTheta),Mom(:,1))
! 
!        Eab = Mandelstam_S/dsqrt(Mandelstam_S)
!        Mom(1:4,2) = (/Eab,0d0,0d0,0d0/) - Mom(1:4,1)
! 
!        pa_hat(1:4) = this%pa(1:4)
!        p_boost(1:4) = this%pa(1:4)+this%pb(1:4)
!        call boost_from_CMS_to_RefMom(pa_hat,p_boost)
!        CosTheta_hat = pa_hat(4)/dsqrt( pa_hat(2)**2+pa_hat(3)**2+pa_hat(4)**2 )
!        phi_hat = datan( pa_hat(3)/pa_hat(2) )
!        if( pa_hat(2).lt.0d0 ) phi_hat = phi_hat + pi
!        
!        call rotate3D_phi_theta(phi_hat,dacos(CosTheta_hat),Mom(:,1))
!        call rotate3D_phi_theta(phi_hat,dacos(CosTheta_hat),Mom(:,2))
!        p_boost(1:4) = (/ this%pa(1)+this%pb(1), -(this%pa(2:4)+this%pb(2:4)) /)
!        call boost_from_CMS_to_RefMom(p_boost,Mom(:,1))
!        call boost_from_CMS_to_RefMom(p_boost,Mom(:,2))      
! 
!        Jac = 1d0/g_p( Mandelstam_S, pa_sq, pb_sq, Mandelstam_T, this%PropMass_sq, this%power, this%InvMinMax_sq(1), this%InvMinMax_sq(2))
! 
!   RETURN
!   END SUBROUTINE




  
  
END MODULE