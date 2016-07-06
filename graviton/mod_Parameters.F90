MODULE ModParameters
implicit none
save

integer, public :: Collider, PDFSet,PChannel,Process,DecayMode
integer, public :: VegasIt1,VegasNc1,Collider_Energy
integer, public :: VegasIt1_default,VegasNc1_default
integer, public :: NumHistograms
logical, public :: unweighted
integer(8), public :: EvalCounter=0
integer(8), public :: RejeCounter=0
integer(8), public :: AccepCounter=0
integer(8), public :: AccepCounter_part(-5:5,-5:5)=0
character :: DataFile*(20)

logical, public, parameter :: seed_random = .true.

logical, public, parameter :: fix_channels_ratio = .false.
real(8), public, parameter :: channels_ratio_fix = 0.5d0    ! desired ratio of
                                                            ! N_qq/(N_qq+N_gg)


! we are using units of 100GeV, i.e. Lambda=10 is 1TeV
real(8), public, parameter :: m_Z     = 91.1876d0    /100d0  ! Z boson mass (PDG-2008)
real(8), public, parameter :: m_Grav  = 250d0      /100d0  !  X  mass (spin 0, spin 1, spin 2)
real(8), public, parameter :: Ga_Z    = 2.4952d0    /100d0  ! Z boson width(PDG-2008)! NOT USED; Z's are on-shell
real(8), public, parameter :: Ga_Grav = 10d0        /100d0  ! Graviton width (not used directly)
real(8), public, parameter :: Lambda  = 1000d0      /100d0  ! Lambda coupling enters in two places
                                                            ! overal scale for x-section and in power suppressed
                                                            ! operators/formfactors (former r).
real(8), public, parameter :: alpha_QED = 1d0/128.0d0       ! el.magn. coupling
real(8), public, parameter :: sitW = dsqrt(0.23119d0)       ! sin(Theta_Weinberg) (PDG-2008)
real(8), public, parameter :: Mu_Fact = m_Grav              ! pdf factorization scale
real(8), public, parameter :: LHC_Energy=14000d0    /100d0  ! LHC hadronic center of mass energy (LHC-expected)
real(8), public, parameter :: TEV_Energy=1960d0     /100d0  ! Tevatron hadronic center of mass energy

real(8), public, parameter :: Br_Z_up = 0.1657d0  ! branching fraction Ga(up)/Ga(hadronic)
real(8), public, parameter :: Br_Z_ch = 0.1657d0  ! branching fraction Ga(charm)/Ga(hadronic)
real(8), public, parameter :: Br_Z_dn = 0.2229d0  ! branching fraction Ga(down)/Ga(hadronic)
real(8), public, parameter :: Br_Z_st = 0.2229d0  ! branching fraction Ga(strange)/Ga(hadronic)
real(8), public, parameter :: Br_Z_bo = 1d0-Br_Z_up-Br_Z_ch-Br_Z_dn-Br_Z_st  ! branching fraction Ga(bottom)/Ga(hadronic)

!-- parameters that define spin 0 coupling to SM fields, see note
   complex(8), public, parameter :: ahg1 = (1.0d0,0d0)
   complex(8), public, parameter :: ahg2 = (0.0d0,0d0)
   complex(8), public, parameter :: ahg3 = (0.0d0,0d0)  ! pseudoscalar


   complex(8), public, parameter :: ahz1 = (1.0d0,0d0)
   complex(8), public, parameter :: ahz2 = (0.0d0,0d0)
   complex(8), public, parameter :: ahz3 = (0.0d0,0d0)  ! pseudoscalar

!---parameters that define spin 1 coupling to SM fields, see note

   complex(8), public, parameter :: zprime_qq_left  = (1.0d0,0d0)
   complex(8), public, parameter :: zprime_qq_right = (0.0d0,0d0)
   complex(8), public, parameter :: zprime_zz_v =  (1.0d0,0d0)
   complex(8), public, parameter :: zprime_zz_a =  (0.0d0,0d0)

!-- parameters that define spin 2 coupling to SM fields, see note
! minimal coupling corresponds to a1 = b1 = b5 = 1 everything else 0

  complex(8), public, parameter :: a1 = (1.0d0,0d0)    ! g1  -- c.f. draft
  complex(8), public, parameter :: a2 = (0.0d0,0d0)    ! g2
  complex(8), public, parameter :: a3 = (0.0d0,0d0)    ! g3
  complex(8), public, parameter :: a4 = (0.0d0,0d0)    ! g4
  complex(8), public, parameter :: a5 = (0.0d0,0d0)    ! pseudoscalar, g8

  logical, public, parameter :: generate_bis = .true.

  complex(8), public, parameter :: b1 = (1.0d0,0d0)    !  all b' below are g's in the draft
  complex(8), public, parameter :: b2 = (0.0d0,0d0)
  complex(8), public, parameter :: b3 = (0.0d0,0d0)
  complex(8), public, parameter :: b4 = (0.0d0,0d0)
  complex(8), public, parameter :: b5 = (1.0d0,0d0)
  complex(8), public, parameter :: b6 = (0.0d0,0d0)
  complex(8), public, parameter :: b7 = (0.0d0,0d0)
  complex(8), public, parameter :: b8 = (0.0d0,0d0)
  complex(8), public, parameter :: b9 = (0.0d0,0d0)
  complex(8), public, parameter :: b10 = (0.0d0,0d0)


  complex(8), public, parameter  :: c1 = (1.0d0,0d0)
  complex(8), public, parameter  :: c2 = (-0.5d0,0d0)
  complex(8), public, parameter  :: c3 = (0.0d0,0d0)
  complex(8), public, parameter  :: c4 = (-0.9d0,0d0)
  complex(8), public, parameter  :: c5 = (1d0,0d0)
  complex(8), public, parameter  :: c6 = (0d0,0d0)
  complex(8), public, parameter  :: c7 = (0d0,0d0)


! Z-f-fbar couplings:
!   g_R(f) = -sw/cw*Q(f)               = 1/2/sw/cw * aR,
!   g_L(f) = 1/sw/cw*T3(f) -sw/cw*Q(f) = 1/2/sw/cw * aL
! with
!   aR(f) = -2*sw**2*Q(f),
!   aL(f) = -2*sw**2*Q(f) + 2*T3(f).
! So that
real(8), public, parameter :: aR_lep =-2d0*sitW**2*(-1d0)
real(8), public, parameter :: aL_lep =-2d0*sitW**2*(-1d0)-1d0
real(8), public, parameter :: aR_QUp =-2d0*sitW**2*(2d0/3d0)
real(8), public, parameter :: aL_QUp =-2d0*sitW**2*(2d0/3d0)+1d0
real(8), public, parameter :: aR_QDn =-2d0*sitW**2*(-1d0/3d0)
real(8), public, parameter :: aL_QDn =-2d0*sitW**2*(-1d0/3d0)-1d0



real(8), public, parameter :: fbGeV2=0.389379d12/(100d0**2)
real(8), public, parameter :: SymmFac=1d0/2d0, SpinAvg=1d0/4d0, QuarkColAvg=1d0/3d0, GluonColAvg=1d0/8d0
integer, public, target :: Up_  = 1
integer, public, target :: Dn_  = 2
integer, public, target :: Chm_ = 3
integer, public, target :: Str_ = 4
integer, public, target :: Top_ = 5
integer, public, target :: Bot_ = 6
integer, public, target :: ElP_ = 7
integer, public, target :: MuP_ = 8
integer, public, target :: TaP_ = 9
integer, public, target :: Glu_ = 10
integer, public, target :: Pho_ = 11
integer, public, target :: Z0_  = 12
integer, public, target :: Wp_  = 13

integer, public, target :: AUp_  = -1
integer, public, target :: ADn_  = -2
integer, public, target :: AChm_ = -3
integer, public, target :: AStr_ = -4
integer, public, target :: ATop_ = -5
integer, public, target :: ABot_ = -6
integer, public, target :: ElM_  = -7
integer, public, target :: MuM_  = -8
integer, public, target :: TaM_  = -9
integer, public, target :: Wm_   = -13


real(8), public, parameter :: pi =3.141592653589793238462643383279502884197d0
real(8), public, parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
real(8), public, parameter :: pisq = pi**2
real(8), public, parameter :: one = 1.0d0, mone = -1.0d0
real(8), public, parameter :: half  = 0.5d0,two = 2.0d0
real(8), public, parameter :: zero  = 0.0d0
complex(8), parameter, public :: czero = 0.0d0
complex(8), parameter, public :: cone = 1.0d0
complex(8), parameter, public :: ci=(0.0d0,1.0d0)
complex(8), parameter, public :: ne=(0.0d0,1.0d0)



contains



FUNCTION convertLHE(Part)
implicit none
integer :: convertLHE
integer :: Part


  if(     Part.eq.0 ) then
      convertLHE = 0
  elseif(     Part.eq.Glu_ ) then
      convertLHE = 21
  elseif( Part.eq.ElM_ ) then
      convertLHE = 11
  elseif( Part.eq.ElP_ ) then
      convertLHE =-11
  elseif( Part.eq.MuM_ ) then
      convertLHE = 13
  elseif( Part.eq.MuP_ ) then
      convertLHE =-13
  elseif( Part.eq.Up_  ) then
      convertLHE = 2
  elseif( Part.eq.AUp_ ) then
      convertLHE =-2
  elseif( Part.eq.Dn_  ) then
      convertLHE = 1
  elseif( Part.eq.ADn_ ) then
      convertLHE =-1
  elseif( Part.eq.Chm_ ) then
      convertLHE = 4
  elseif( Part.eq.AChm_) then
      convertLHE =-4
  elseif( Part.eq.Str_ ) then
      convertLHE = 3
  elseif( Part.eq.AStr_) then
      convertLHE =-3
  elseif( Part.eq.Bot_ ) then
      convertLHE = 5
  elseif( Part.eq.ABot_) then
      convertLHE =-5
  elseif( Part.eq.Z0_) then
      convertLHE =23
  else
      print *, "LHE format not implemented for ",Part
      stop
  endif

END FUNCTION



END MODULE



