!--YaofuZhou-----------------------------------------
! finite factor multiplied by LO cross section for
! "I bit" "K bit" and "P bit" from arXiv:0709.2881
module ModVHIKP
  use ModParameters
  use ModMisc
  use ModVHaux
  implicit none
  
public :: I_qq_vhg
public :: K_qq_vhg
public :: P_qq_vhg
public :: Kx_qq_vhg
public :: Px_qq_vhg
public :: K1_qq_vhg
public :: P1_qq_vhg
public :: K_gq_vhq
public :: P_gq_vhq

contains



real(8) function I_qq_vhg(shat)
  implicit none
  real(8), intent(in) :: shat
  I_qq_vhg = 10d0 - Pi**2 + 3d0*dlog(Mu_Fact**2/shat) + (dlog(Mu_Fact**2/shat))**2 &
           - 2d0*(dlog(alpha_dip))**2 + 3d0*(alpha_dip-1d0-dlog(alpha_dip))

  I_qq_vhg = I_qq_vhg-1d0!/2d0<--------fuck this last bug :)

  I_qq_vhg = I_qq_vhg *alphas*Cf/2d0/pi
  return
end function I_qq_vhg


!real(8) function K_qq_vhg(alphasx,me2_unpol,me2_unpolx,x)
real(8) function K_qq_vhg(me2_unpol,me2_unpolx,alphasx,x)
  implicit none
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: me2_unpol,me2_unpolx
  real(8), intent(in) :: x
  real(8) :: ThetaTerm

  K_qq_vhg = 0d0

  if(x.gt.(1d0-alpha_dip))then
    K_qq_vhg = K_qq_vhg + (me2_unpolx*alphasx - me2_unpol*alphas) &
                          *2d0/(1d0-x)*dlog((1d0-x)**2/x)
    ThetaTerm = dlog(2d0-x)
  else
    K_qq_vhg = K_qq_vhg + (me2_unpolx*alphasx - me2_unpol*alphas) &
                          *2d0/(1d0-x)*dlog((1d0-x)/x)
    ThetaTerm = dlog((2d0-x)/(1d0-x))
  endif

  K_qq_vhg = K_qq_vhg + me2_unpolx *alphasx &
                      *( -(1d0+x)*dlog((1d0-x)**2/x) + (1d0 - x) &
                       + (1d0+x**2)/(1d0-x) * dlog(min(1d0,alpha_dip/(1d0-x))) &
                       + 2d0/(1d0-x) * (dlog(2d0-x) -ThetaTerm) )

  K_qq_vhg = K_qq_vhg + me2_unpol *alphas &
                      *( 2d0/3d0*pi**2 - 5d0 + 2d0*(dlog(alpha_dip))**2 - 1.5d0*(alpha_dip-1d0-dlog(alpha_dip)) )

  K_qq_vhg = K_qq_vhg *Cf/2d0/pi

  return
end function K_qq_vhg






real(8) function P_qq_vhg(me2_unpol,me2_unpolx,alphasx,Mu_Fact_x,shat,x)
  implicit none
  real(8), intent(in) :: Mu_Fact_x
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: me2_unpol,me2_unpolx
  real(8), intent(in) :: shat
  real(8), intent(in) :: x

  P_qq_vhg = 0d0

  P_qq_vhg = -1d0 & !color
           * (me2_unpolx*alphasx*dlog(Mu_Fact_x**2/shat/x) &
            - me2_unpol *alphas *dlog(Mu_Fact**2  /shat  )) &
           * (1d0+x**2)/(1d0-x)

  P_qq_vhg = P_qq_vhg * Cf/2d0/pi

  return
end function P_qq_vhg










real(8) function Kx_qq_vhg(me2_unpolx,alphasx,x)
  implicit none
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: me2_unpolx
  real(8), intent(in) :: x
  real(8) :: ThetaTerm

  Kx_qq_vhg = 0d0

  if(x.gt.(1d0-alpha_dip))then
    Kx_qq_vhg = Kx_qq_vhg + (me2_unpolx*alphasx) &
                          *2d0/(1d0-x)*dlog((1d0-x)**2/x)
    ThetaTerm = dlog(2d0-x)
  else
    Kx_qq_vhg = Kx_qq_vhg + (me2_unpolx*alphasx) &
                          *2d0/(1d0-x)*dlog((1d0-x)/x)
    ThetaTerm = dlog((2d0-x)/(1d0-x))
  endif

  Kx_qq_vhg = Kx_qq_vhg + me2_unpolx *alphasx &
                      *( -(1d0+x)*dlog((1d0-x)**2/x) + (1d0 - x) &
                       + (1d0+x**2)/(1d0-x) * dlog(min(1d0,alpha_dip/(1d0-x))) &
                       + 2d0/(1d0-x) * (dlog(2d0-x) -ThetaTerm) )



  Kx_qq_vhg = Kx_qq_vhg *Cf/2d0/pi

  return
end function Kx_qq_vhg



real(8) function K1_qq_vhg(me2_unpol,x)
  implicit none
  real(8), intent(in) :: me2_unpol
  real(8), intent(in) :: x
  real(8) :: ThetaTerm

  K1_qq_vhg = 0d0

  if(x.gt.(1d0-alpha_dip))then
    K1_qq_vhg = K1_qq_vhg + ( - me2_unpol*alphas) &
                          *2d0/(1d0-x)*dlog((1d0-x)**2/x)
  else
    K1_qq_vhg = K1_qq_vhg + ( - me2_unpol*alphas) &
                          *2d0/(1d0-x)*dlog((1d0-x)/x)
  endif

  K1_qq_vhg = K1_qq_vhg + me2_unpol *alphas &
                      *( 2d0/3d0*pi**2 - 5d0 + 2d0*(dlog(alpha_dip))**2 - 1.5d0*(alpha_dip-1d0-dlog(alpha_dip)) )

  K1_qq_vhg = K1_qq_vhg *Cf/2d0/pi

  return
end function K1_qq_vhg


real(8) function Px_qq_vhg(me2_unpolx,alphasx,Mu_Fact_x,shat,x)
  implicit none
  real(8), intent(in) :: Mu_Fact_x
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: me2_unpolx
  real(8), intent(in) :: shat
  real(8), intent(in) :: x

  Px_qq_vhg = 0d0

  Px_qq_vhg = -me2_unpolx*alphasx*dlog(Mu_Fact_x**2/shat/x) * (1d0+x**2)/(1d0-x)

  Px_qq_vhg = Px_qq_vhg * Cf/2d0/pi

  return
end function Px_qq_vhg





real(8) function P1_qq_vhg(me2_unpol,shat,x)
  implicit none
  real(8), intent(in) :: me2_unpol
  real(8), intent(in) :: shat
  real(8), intent(in) :: x

  P1_qq_vhg = 0d0

  P1_qq_vhg = me2_unpol *alphas *dlog(Mu_Fact**2  /shat  )* (1d0+x**2)/(1d0-x)

  P1_qq_vhg = P1_qq_vhg * Cf/2d0/pi

  return
end function P1_qq_vhg











!real(8) function KP_qq_vhg(me2_unpol,me2_unpolx,shat,x)!Till's final expression
!  implicit none
!  real(8), intent(in) :: me2_unpol,me2_unpolx
!  real(8), intent(in) :: shat
!  real(8), intent(in) :: x
!
!  KP_qq_vhg = 0d0
!
!  KP_qq_vhg = KP_qq_vhg + me2_unpolx * ((1d0+x**2)/(1d0-x) * dlog( (1d0-x)**2 * shat / Mu_Fact**2 ) +1d0-x)
!
!  KP_qq_vhg = KP_qq_vhg - me2_unpol  * (2d0/(1d0-x) * dlog( (1d0-x)**2 * shat / Mu_Fact**2/x ) +5d0 + 1.5d0*dlog(Mu_Fact_x**2/shat) -2d0/3d0*pi**2)
!
!  KP_qq_vhg = KP_qq_vhg *alphas * Cf/2d0/pi
!
!  return
!end function KP_qq_vhg




real(8) function K_gq_vhq(alphasx,x)
  implicit none
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: x

  K_gq_vhq = 0d0

  K_gq_vhq = ( 2d0*(x**2) - 2d0*x +1d0 ) * ( dlog( (1d0-x)**2/x ) + dlog( min( 1d0, alpha_dip/(1d0-x) ) ) ) &
           + 2d0*x*(1d0-x) 

  K_gq_vhq = K_gq_vhq * alphasx * Tr / (2d0*pi)

  return
end function K_gq_vhq




real(8) function P_gq_vhq(alphasx,Mu_Fact_x,shat,x)
  implicit none
  real(8), intent(in) :: Mu_Fact_x
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: shat
  real(8), intent(in) :: x

  P_gq_vhq = (-1d0) * (2d0*(x**2) - 2d0*x +1d0) * dlog(Mu_Fact_x**2 / shat/x)

  P_gq_vhq = P_gq_vhq * Tr * alphasx / (2d0*pi)

  return
end function P_gq_vhq



!real(8) function theta_step(a,b)
!  implicit none
!  real(8), intent(in) :: a,b
!
!  if(a.gt.b)then
!    theta_step = 1d0
!  else
!    theta_step = 0d0
!  endif
!
!  return
!end function theta_step


end module ModVHIKP
!!--YaofuZhou-----------------------------------------
