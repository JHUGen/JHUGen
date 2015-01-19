include './variables.F90'

   integer, parameter :: ZZMode=00,ZgsMode=01,gsZMode=02,gsgsMode=03
   integer, parameter :: WWMode=10
   integer, parameter :: ggMode=20
   integer, parameter :: ZgMode=30,gsgMode=31

!-- parameters that define on-shell spin 0 coupling to SM fields, see note
   logical,  parameter :: generate_as = .false.! this cannot be changed
   complex(8),  parameter :: ahg1 = (1.0d0,0d0)! these parameters are not used
   complex(8),  parameter :: ahg2 = (0.0d0,0d0)
   complex(8),  parameter :: ahg3 = (0.0d0,0d0)  ! pseudoscalar
   complex(8),  parameter :: ahz1 = (1.0d0,0d0)
   complex(8),  parameter :: ahz2 = (0.0d0,0d0)  ! this coupling does not contribute for gamma+gamma final states
   complex(8),  parameter :: ahz3 = (0.0d0,0d0)  ! pseudoscalar

!-- parameters that define off-shell spin 0 coupling to SM fields, see note
   complex(8) :: ghg2
   complex(8) :: ghg3
   complex(8) :: ghg4   ! pseudoscalar
   complex(8) :: ghz1
   complex(8) :: ghz2
   complex(8) :: ghz3
   complex(8) :: ghz4  ! pseudoscalar

   complex(8) :: ghzgs2  
   complex(8) :: ghzgs3  
   complex(8) :: ghzgs4  
   complex(8) :: ghgsgs2 
   complex(8) :: ghgsgs3 
   complex(8) :: ghgsgs4 

   real(8),parameter :: MPhotonCutoff = 4d0*GeV


!---parameters that define spin 1 coupling to SM fields, see note
   complex(8) :: zprime_qq_left
   complex(8) :: zprime_qq_right
   complex(8) :: zprime_zz_1!  =1 for JP=1- vector
   complex(8) :: zprime_zz_2!  =1 for JP=1+ pseudovector


!-- parameters that define spin 2 coupling to SM fields, see note
! minimal coupling corresponds to a1 = b1 = b5 = 1 everything else 0
  complex(8) :: a1    ! g1  -- c.f. draft
  complex(8) :: a2    ! g2
  complex(8) :: a3    ! g3
  complex(8) :: a4    ! g4
  complex(8) :: a5    ! pseudoscalar, g8
  complex(8) :: graviton_qq_left! graviton coupling to quarks
  complex(8) :: graviton_qq_right

!-- see mod_Graviton
  logical,  parameter :: generate_bis = .true.! this cannot be changed
  logical,  parameter :: use_dynamic_MG = .true.

  complex(8) :: b1    !  all b' below are g's in the draft
  complex(8) :: b2
  complex(8) :: b3
  complex(8) :: b4
  complex(8) :: b5
  complex(8) :: b6
  complex(8) :: b7
  complex(8) :: b8
  complex(8) :: b9  ! this coupling does not contribute to gamma+gamma final states
  complex(8) :: b10  ! this coupling does not contribute to gamma+gamma final states


  complex(8),  parameter  :: c1 = (1.0d0,0d0)! these parameters are not used
  complex(8),  parameter  :: c2 = (0.0d0,0d0)
  complex(8),  parameter  :: c3 = (0.0d0,0d0)
  complex(8),  parameter  :: c41= (0.0d0,0d0)
  complex(8),  parameter  :: c42= (0.0d0,0d0)
  complex(8),  parameter  :: c5 = (0.0d0,0d0)
  complex(8),  parameter  :: c6 = (0.0d0,0d0)  ! this coupling does not contribute to gamma+gamma final states
  complex(8),  parameter  :: c7 = (0.0d0,0d0)  ! this coupling does not contribute to gamma+gamma final states



!  couplings for ttbar+H
!    complex(8), parameter :: kappa       = (1d0,0d0)
!    complex(8), parameter :: kappa_tilde = (0d0,0d0) 
!    complex(8), parameter :: couplHTT_right_dyn = m_top/vev/2d0 * ( kappa + (0d0,1d0)*kappa_tilde )
!    complex(8), parameter :: couplHTT_left_dyn  = m_top/vev/2d0 * ( kappa - (0d0,1d0)*kappa_tilde )





real(8),  parameter :: aR_lep =-2d0*sitW**2*(-1d0)
real(8),  parameter :: aL_lep =-2d0*sitW**2*(-1d0)-1d0
real(8),  parameter :: aR_neu =-2d0*sitW**2*(0d0)
real(8),  parameter :: aL_neu =-2d0*sitW**2*(0d0)+1d0

real(8),  parameter :: aR_QUp =-2d0*sitW**2*(2d0/3d0)
real(8),  parameter :: aL_QUp =-2d0*sitW**2*(2d0/3d0)+1d0
real(8),  parameter :: aR_QDn =-2d0*sitW**2*(-1d0/3d0)
real(8),  parameter :: aL_QDn =-2d0*sitW**2*(-1d0/3d0)-1d0

real(8),  parameter :: bL = dsqrt(2d0)*dsqrt(1d0-sitW**2)
real(8),  parameter :: bR = 0d0

real(8), parameter :: cR_lep = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0)
real(8), parameter :: cL_lep = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0)
real(8), parameter :: cR_neu = -2d0*sitW*dsqrt(1d0-sitW**2)*(0d0)
real(8), parameter :: cL_neu = -2d0*sitW*dsqrt(1d0-sitW**2)*(0d0)
real(8), parameter :: cR_QUp = -2d0*sitW*dsqrt(1d0-sitW**2)*(2d0/3d0)
real(8), parameter :: cL_QUp = -2d0*sitW*dsqrt(1d0-sitW**2)*(2d0/3d0)
real(8), parameter :: cR_QDn = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0/3d0)
real(8), parameter :: cL_QDn = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0/3d0)

real(8),  parameter :: fbGeV2=0.389379d12/(100d0**2)
real(8),  parameter :: SymmFac=1d0/2d0, SpinAvg=1d0/4d0, QuarkColAvg=1d0/3d0, GluonColAvg=1d0/8d0
integer,  target :: Up_  = 1
integer,  target :: Dn_  = 2
integer,  target :: Chm_ = 3
integer,  target :: Str_ = 4
integer,  target :: Top_ = 5
integer,  target :: Bot_ = 6
integer,  target :: ElP_ = 7
integer,  target :: MuP_ = 8
integer,  target :: TaP_ = 9
integer,  target :: Glu_ = 10
integer,  target :: Pho_ = 11
integer,  target :: Z0_  = 12
integer,  target :: Wp_  = 13
integer,  target :: NuE_ = 14
integer,  target :: NuM_ = 15
integer,  target :: NuT_ = 16
integer,  target :: Hig_ = 25

integer,  target :: AUp_  = -1
integer,  target :: ADn_  = -2
integer,  target :: AChm_ = -3
integer,  target :: AStr_ = -4
integer,  target :: ATop_ = -5
integer,  target :: ABot_ = -6
integer,  target :: ElM_  = -7
integer,  target :: MuM_  = -8
integer,  target :: TaM_  = -9
integer,  target :: Wm_   = -13
integer,  target :: ANuE_ = -14
integer,  target :: ANuM_ = -15
integer,  target :: ANuT_ = -16


real(8),  parameter :: pi =3.141592653589793238462643383279502884197d0
real(8),  parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
real(8),  parameter :: pisq = pi**2
real(8),  parameter :: one = 1.0d0, mone = -1.0d0
real(8),  parameter :: half  = 0.5d0,two = 2.0d0
real(8),  parameter :: zero  = 0.0d0
complex(8), parameter :: czero = 0.0d0
complex(8), parameter :: cone = 1.0d0
complex(8), parameter :: ci=(0.0d0,1.0d0)
complex(8), parameter :: ne=(0.0d0,1.0d0)


real(8) :: M_V,Ga_V
real(8), parameter :: tol = 0.0000001d0


    if( DecayMode1.le.3 ) then
       M_V = M_Z
       Ga_V= Ga_Z
    elseif( (DecayMode1.ge.4) .and. (DecayMode1.le.6) ) then
       M_V = M_W
       Ga_V= Ga_W    
    elseif( DecayMode1.eq.7 ) then
       M_V = 0d0
       Ga_V= 0d0    
    endif


