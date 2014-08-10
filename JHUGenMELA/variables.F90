real(8), parameter :: GeV=1d0   /100d0
real(8), parameter :: M_Z     = 91.1876d0 *GeV      ! Z boson mass (PDG-2011)
real(8), parameter :: Ga_Z    = 2.4952d0  *GeV      ! Z boson width(PDG-2011)
real(8), parameter :: M_W     = 80.399d0  *GeV      ! W boson mass (PDG-2011)
real(8), parameter :: m_tau = 1.8d0  *GeV           ! tau lepton mass
real(8), parameter :: Ga_W    = 2.085d0   *GeV      ! W boson width(PDG-2011)
real(8), parameter :: Lambda  = 1000d0    *GeV      ! Lambda coupling enters in two places
                                                            ! overal scale for x-section and in power suppressed
                                                            ! operators/formfactors (former r).


real(8), parameter :: Gf = 1.16639d-5/GeV**2        ! fermi constant
real(8), parameter :: vev = 1.0d0/dsqrt(Gf*dsqrt(2.0d0))
real(8), parameter :: gwsq = 4.0d0 * M_W**2/vev**2  ! weak constant squared
real(8), parameter :: alpha_QED = 1d0/128.0d0       ! el.magn. coupling
real(8), parameter :: alphas = 0.13229060d0         ! strong coupling
real(8), parameter :: sitW = dsqrt(0.23119d0)       ! sin(Theta_Weinberg) (PDG-2008)


integer, parameter :: DecayMode1 = 0
integer, parameter :: DecayMode2 = 0
logical, parameter :: includeInterference = .false.
logical, parameter :: OffShellReson = .false.
logical, parameter :: includeGammaStar = .false.
logical, parameter :: H_DK = .true.                 ! set to false so H in V* > VH does not decay


!-- momentum-dependent form factors
complex(8) :: ghg2_dyn
complex(8) :: ghg3_dyn
complex(8) :: ghg4_dyn

complex(8) :: ghz1_prime 
complex(8) :: ghz1_prime2
complex(8) :: ghz1_prime3
complex(8) :: ghz1_prime4
complex(8) :: ghz1_prime5
complex(8) :: ghz1_prime6
complex(8) :: ghz1_prime7

complex(8) :: ghz2_prime 
complex(8) :: ghz2_prime2
complex(8) :: ghz2_prime3
complex(8) :: ghz2_prime4
complex(8) :: ghz2_prime5
complex(8) :: ghz2_prime6
complex(8) :: ghz2_prime7

complex(8) :: ghz3_prime 
complex(8) :: ghz3_prime2
complex(8) :: ghz3_prime3
complex(8) :: ghz3_prime4
complex(8) :: ghz3_prime5
complex(8) :: ghz3_prime6
complex(8) :: ghz3_prime7

complex(8) :: ghz4_prime 
complex(8) :: ghz4_prime2
complex(8) :: ghz4_prime3
complex(8) :: ghz4_prime4
complex(8) :: ghz4_prime5
complex(8) :: ghz4_prime6
complex(8) :: ghz4_prime7

complex(8) :: ghzgs1_prime2

real(8), parameter :: lambda_z1 = 10000d0 * GeV
real(8), parameter :: lambda_z2 = 10000d0 * GeV
real(8), parameter :: lambda_z3 = 10000d0 * GeV
real(8), parameter :: lambda_z4 = 10000d0 * GeV
real(8), parameter :: lambda_z5 = 10000d0 * GeV
real(8), parameter :: lambda_Q  = 10000d0 * GeV

real(8), parameter :: scale_alpha_Z_uu = 1.04282d0 ! scaling factor of alpha (~partial width) for Z > u u~
real(8), parameter :: scale_alpha_Z_dd = 1.04282d0 ! scaling factor of alpha (~partial width) for Z > d d~
real(8), parameter :: scale_alpha_Z_ll = 1d0 ! scaling factor of alpha (~partial width) for Z > l+ l-
real(8), parameter :: scale_alpha_Z_nn = 1d0 ! scaling factor of alpha (~partial width) for Z > nu nu~
real(8), parameter :: scale_alpha_W_ud = 1.0993819d0 ! scaling factor of alpha (~partial width) for W > u d
real(8), parameter :: scale_alpha_W_cs = 1.0993819d0 ! scaling factor of alpha (~partial width) for W > c s
real(8), parameter :: scale_alpha_W_ln = 1d0 ! scaling factor of alpha (~partial width) for W > l nu