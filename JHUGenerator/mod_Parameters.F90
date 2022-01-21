MODULE ModParameters
implicit none
save
!
!
character(len=*),parameter :: JHUGen_Version="v7.5.1"
!
!
!=====================================================
!internal
integer, public, parameter :: dp = selected_real_kind(15)
real(8), public, parameter :: tol = 0.0000001d0
integer, public, parameter :: InvalidMode=-1,WWMode=00,ZZMode=01,ZgsMode=02,gsZMode=03,gsgsMode=04,ZgMode=05,gsgMode=06,ggMode=07
integer, public, parameter :: WWpMode=10,WpWMode=11,WpWpMode=12
integer, public, parameter :: ZZpMode=20,ZpZMode=21,ZpZpMode=22
integer, public, parameter :: gsZpMode=30,ZpgsMode=31,ZpgMode=32
integer, public :: Collider,PChannel,Process,DecayMode1,DecayMode2,TopDecays,TauDecays
integer, public :: VegasIt1,VegasNc0,VegasNc1,VegasNc2,PMZZEvals
real(8), public :: Collider_Energy
integer, public :: FacScheme,RenScheme
real(8), public :: MuFacMultiplier,MuRenMultiplier
integer, public :: VegasIt1_default,VegasNc0_default,VegasNc1_default,VegasNc2_default
integer, public :: NumHistograms
integer, public :: RequestNLeptons(1:2) = -1
integer, public :: RequestNJets(1:2) = -1
integer, public :: RequestOS(1:2) = -1
integer, public :: RequestOSSF(1:2) = -1
logical, public :: Unweighted,OffShellReson,OffShellV1,OffShellV2,ReadLHEFile,ConvertLHEFile,DoPrintPMZZ
logical, public :: ReadCSmax,GenerateEvents,CountTauAsAny,HasLeptonFilter, FoundHiggsMass, FoundHiggsWidth
integer, public :: WriteFailedEvents
logical, public :: FilesOpened = .false.
integer, public, parameter :: kRenFacScheme_default=0
integer, public, parameter :: kRenFacScheme_mhstar=1
integer, public, parameter :: kRenFacScheme_mjjhstar=2
integer, public, parameter :: kRenFacScheme_mjj_mhstar=3
integer, public, parameter :: kRenFacScheme_mj_mj_mhstar=4
integer, public, parameter :: kRenFacScheme_mjj=5
integer, public, parameter :: kRenFacScheme_mj_mj=6
integer, public, parameter :: kRenFacScheme_mjhstar=7
integer, public, parameter :: kRenFacScheme_mj_mhstar=8
integer, public, parameter :: kRenFacScheme_mj=9
integer, public, parameter :: kRenFacScheme_maxpTj=10
integer, public, parameter :: kRenFacScheme_minpTj=11
integer, public, parameter :: nRenFacSchemes=12
integer, public, parameter :: maxpart = 30
integer, public, parameter :: NMAXCHANNELS=200 ! Maximum 250 from vegas_common xi dimensions
integer(8), public :: EvalCounter=0
integer(8), public :: RejeCounter=0
integer(8), public :: AccepCounter=0
integer(8), public :: AlertCounter=0
integer(8), public :: AccepCounter_part(-6:6,-6:6)=0,RejeCounter_part(-6:6,-6:6)=0,RequEvents(-6:+6,-6:+6)=0
real(8), public :: CrossSecMax(-6:+6,-6:+6)=0,CrossSec(-6:+6,-6:+6)=0
integer(8), public :: RequEvents2(1:NMAXCHANNELS)=0
integer(8), public :: AccepCounter_part2(1:NMAXCHANNELS)=0
real(8), public :: CrossSec2(1:NMAXCHANNELS)=0
real(8), public :: CrossSecMax2(1:NMAXCHANNELS)=0
integer, public :: iPart_sel, jPart_sel, iChann_sel
real(8) :: time_start,time_end,time_int
logical, public :: warmup
character(len=500) :: DataFile
character(len=500) :: CSmaxFile
character(len=500) :: LogFile
character(len=500) :: LHEProdFile
! PDFset variables, present regardless of useLHAPDF value due to MELA
character(len=100) :: LHAPDFString
character(len=500) :: LHAPDF_DATA_PATH
integer, public :: LHAPDFMember, lenLHAPDFString ! lenLHAPDFString is needed in MELA
integer, public :: PDFSet
! End PDFset variables
#if useCollier==1
! COLLIER initialization variables
integer, public :: Collier_maxNLoopProps = -1
integer, public :: Collier_maxRank = -1
! End COLLIER initialization variables
#endif
logical, public :: includeInterference, ReweightInterference, writegit
real(8), public :: M_V,Ga_V, M_Vprime,Ga_Vprime, M_V_ps,Ga_V_ps, M_Z_ps,Ga_Z_ps, M_W_ps,Ga_W_ps
real(dp), public, parameter :: pi =3.141592653589793238462643383279502884197_dp
real(dp), public, parameter :: sqrt2 = 1.4142135623730950488016887242096980786_dp
real(dp), public, parameter :: gamma_0 = 0.5772156649015328606065120900824024310421_dp  !Euler–Mascheroni constant
real(dp), public, parameter :: pisq = pi**2
real(8), public, parameter :: one = 1.0d0, mone = -1.0d0
real(8), public, parameter :: half  = 0.5d0,two = 2.0d0
real(8), public, parameter :: zero  = 0d0
complex(8), parameter, public :: czero = (0d0,0d0)
complex(8), parameter, public :: cone = 1.0d0
complex(8), parameter, public :: ci=(0d0,1.0d0)
complex(8), parameter, public :: ne=(0d0,1.0d0)
real(8), public, parameter :: percent=1d0/100d0
real(8), public, parameter :: SymmFac=1d0/2d0, SpinAvg=1d0/4d0, QuarkColAvg=1d0/3d0, GluonColAvg=1d0/8d0
! Units
real(8), public, parameter :: GeV=1d0/100d0 ! we are using units of 100GeV, i.e. Lambda=10 is 1TeV
real(8), public, parameter :: TeV=1d+3*GeV
real(8), public, parameter :: MeV=1d-3*GeV
real(8), public, parameter :: keV=1d-6*GeV
real(8), public, parameter :: eV=1d-9*GeV
real(8), public, parameter :: meter=1d6
real(8), public, parameter :: milimeter=1d-3*meter
real(8), public, parameter :: ctauUnit=milimeter
real(8), public, parameter :: barn=1d-28*(meter**2)
real(8), public, parameter :: picobarn=1e-12*barn
real(8), public, parameter :: femtobarn=1e-15*barn
real(8), public, parameter :: attobarn=1e-18*barn
real(8), public, parameter :: second=1d15
real(8), public, parameter :: milisecond=1d-3*second
real(8), public, parameter :: SpeedOfLight=299792458*meter/second
real(8), public, parameter :: PlanckConstant=4.13567d-15*eV*second
real(8), public, parameter :: hbar=PlanckConstant/(2d0*pi)
real(8), public, parameter :: hbarc=hbar*SpeedOfLight
real(8), public, parameter :: hbarc2=hbarc**2
real(8), public, parameter :: hbarc2_fbGeV2 = hbarc2/(femtobarn*(GeV**2))
real(8), public, parameter :: hbarc2_pbGeV2 = hbarc2/(picobarn*(GeV**2))
real(8), public, parameter :: hbarc2XsecUnit = hbarc2_fbGeV2*(GeV**2) ! Units of xsec are reported in fb. Mutliplication by GeV**2 is to cancel out the units of ME*PS.
! real(8),public :: GlobalMax=-1d99
! real(8),public :: GlobalMin=+1d99
! integer,parameter :: NPart=200
! real(8),public :: PartitionMax(0:NPart,0:NPart)=-1d99
real(8),public :: minCS=1d10,maxCS=0d0,avgCS=0d0
integer, public :: Br_Z_ll_counter=0
integer, public :: Br_Z_inv_counter=0
integer, public :: Br_Z_uu_counter=0
integer, public :: Br_Z_dd_counter=0
integer, public :: Br_W_ll_counter=0
integer, public :: Br_W_ud_counter=0
integer, public :: Br_counter(1:5,1:5)=0
integer, public :: LeptInEvent(0:8) = 0
integer, public :: JetsInEvent(0:8) = 0
logical, public :: ReweightDecay = .false.
integer, public :: UserSeed = 0
integer, public  :: WidthScheme = 0   ! 1=running BW-width, 2=fixed BW-width (default), 3=Passarino's CPS
integer, public  :: WidthSchemeIn = 0   ! 1=running BW-width, 2=fixed BW-width (default), 3=Passarino's CPS
real(8), public :: mubarH = -999d0   !for CPS
real(8), public :: gabarH = -999d0   !for CPS
real(8), public :: maxInputmHstar = -999d0, minInputmHstar = 1d15, mHstarforphasespace
logical, public :: ReadPMZZ
character(len=500), public :: PMZZfile = "PMZZdistribution.out"
real(8), public :: PMZZ_mReso = -1d0
integer, public, parameter :: PMZZsize = 10000
real(8), public :: PMZZdistribution(1:PMZZsize,1:2)  !huge array, in normal cases will never get near the edge
integer, public :: PMZZminindex=-1, PMZZmaxindex=-1  !store the largest and smallest values currently used
complex(8), public :: PrintPMZZ   !real part is the minimum, imaginary part is the maximum
integer, public :: PrintPMZZIntervals
integer, public :: VBFoffsh_run=-1
logical, public :: FindCrossSectionWithWeights = .false.
real(8), public :: VprimeDecayLengthMassCutoffFactor = -1d0 ! x Ga_Vprime
real(8), public :: CrossSectionWithWeights = 0d0, CrossSectionWithWeightsErrorSquared = 0d0
!=====================================================



!=====================================================
!these fixed random seeds are NOT actually used to generate events
!the seed provided via the command line (or, if none is provided, one generated randomly using the system time)
! is put in place of the first seed that makes a difference in generating a random number
! then the seeds are used to generate the (compiler dependent) required number of seeds
! and THOSE are used for event generation
!note that even if the same seed is provided, results are compiler dependent
integer, public, parameter :: nmaxseeds = 33
integer, public, parameter :: DefaultSeeds(1:nmaxseeds) = (/847362834,470115596,392845769,191039475,372910496,192049687,695820194,218930493,902943834,471748302,123958674,390534012,938576849,386472918,938576483,891928354,593857698,938576432,948576849,192847564     ,218940493,902943854,471748304,123958670,390534013,938576846,386472917,938576484,891928359,593857698,938576434,948576848,192847565/)
integer, public            :: TheSeeds(1:nmaxseeds) = DefaultSeeds
!changing the default seeds is not advised, since then results from before the change
! will not be reproducible
!=====================================================

!=====================================================
!switches - should be set on the command line
logical, public :: fix_channels_ratio = .false.

real(8), public :: channels_ratio_fix = 0.25d0   ! desired ratio of  N_qq/(N_qq+N_gg)

logical, public :: writeWeightedLHE = .false.

logical, public :: RandomizeVVdecays = .true.    ! randomize DecayMode1 and DecayMode2 in H-->VV and TTBAR decays whenever appropriate

logical, public :: UseUnformattedRead = .false.  !Set this to true if the regular reading fails for whatever reason

logical, public :: HbbDecays =.false.                 ! default to false so H in V* > VH (Process = 51) does not decay to bbbar

logical, public :: H_DK =.false.                 ! default to false so H in V* > VH (Process = 50) does not decay to bbbar
!=====================================================

! new VH
character(len=2), public :: VH_PC = "lo"                ! VH partonic channel and mode selection, in development.
                                                        ! "ee" ( = e+ e- @LO)
                                                        ! "gg" ( = triangles + boxes of gg)
                                                        ! "qq" ( = q q~ @LO)
                                                        ! "lo" ( = q q~ @LO)
                                                        ! "tr" ( = triangles of gg)
                                                        ! "bo" ( = boxes of gg)
                                                        ! "in" ( = interference = 2*dble(box*dconjg(triangle)) of gg)
                                                        ! "qg" ( = real - dipoles, for g q/q~ > VH + q/q~, for development only)
                                                        ! "gq" ( = virtual + I + K + P, for g q/q~ > VH + q/q~, for development only)
                                                        ! "sb" ( = real - dipoles, for q q~ @NLO, for development only)
                                                        ! "sp" ( = virtual + I + K + P, for q q~ @NLO, for development only)
                                                        ! "nl" ( = full oneloop = q ~ @LO + NLO + gg + gq)
                                                        ! VH_PC overrides Pchannel.

real(8), public :: alpha_dip = 1d0 !extra non physical degree of freedom for dipoles. Vary to check indepedence (of alpha_dip).
! new VH
!=====================================================

!=====================================================
!cuts - should be set on the command line
real(8), public :: pTjetcut = -1d0*GeV                        ! jet min pt
real(8), public :: etajetcut = -1d0                           ! jet max |eta|
real(8), public :: detajetcut = -1d0                          ! min difference in eta between jets
real(8), public :: Rjet = -1d0                                ! jet deltaR
real(8), public :: mJJcut = 0d0*GeV                           ! minimum mJJ
real(8), public :: m4l_minmax(1:2) = (/ -1d0,-1d0 /)*GeV      ! min and max for m_4l. default is (-1,-1): m_4l ~ Higgs resonance (on-shell)
real(8), public :: m2l_minmax(1:2) = (/ 0d0,14000d0 /)*GeV   ! min and max for m_V in VH production;
real(8), public :: mVH_minmax(1:2) = (/ 0d0,14000d0 /)*GeV      ! min and max for m_VH in VH production;
logical, public :: includeGammaStar = .false.                 ! include offshell photons?
logical, public :: includeVprime = .false.
real(8), public :: MPhotonCutoff = -1d0*GeV                          ! minimum |mass_ll| for offshell photons when includeGammaStar = .true. or in VBF bkg
real(8), public :: pTlepcut = -1d0*GeV                         ! lepton min pt
real(8), public :: etalepcut = 999d0                           ! lepton max |eta|
real(8), public :: pTHcut = 0d0*GeV
logical, public :: JetsOppositeEta = .false.                   ! Ensures associated jets are in opposite hemispheres. Useful to obtain VBF/VBS topology from inclusive production
!=====================================================

!=====================================================
!constants
real(8), public            :: M_Top   = 173.2d0   *GeV      ! top quark mass
real(8), public            :: Ga_Top  = 2.0d0     *GeV      ! top quark width
real(8), public            :: M_Z     = 91.1876d0 *GeV      ! Z boson mass (PDG-2011)
real(8), public            :: Ga_Z    = 2.4952d0  *GeV      ! Z boson width(PDG-2011)
real(8), public            :: M_W     = 80.399d0  *GeV      ! W boson mass (PDG-2011)
real(8), public            :: Ga_W    = 2.085d0   *GeV      ! W boson width(PDG-2011)
real(8), public            :: M_Reso  = 125.0d0   *GeV      ! X resonance mass (spin 0, spin 1, spin 2, can be overwritten by command line argument)
real(8), public            :: Ga_Reso = 0.00407d0 *GeV      ! X resonance width
real(8), public            :: M_Reso2 = -1d0      *GeV      ! second resonance mass (spin 0 in off-shell MCFM, can be overwritten by command line argument)
real(8), public            :: Ga_Reso2= 0d0       *GeV      ! second resonance width

real(8), public            :: HiggsDecayLengthMM = 0d0      ! Higgs decay length in [mm]
real(8), public            :: VprimeDecayLengthMM = 0d0     ! Vprime decay length in [mm]

real(8), public            :: m_bot = 4.75d0       *GeV     ! bottom quark mass
real(8), public            :: m_charm = 1.275d0    *GeV     ! charm quark mass
real(8), public            :: m_el = 0.00051100d0  *GeV     ! electron mass
real(8), public            :: m_mu = 0.10566d0     *GeV     ! muon mass
real(8), public            :: m_tau = 1.7768d0     *GeV     ! tau mass
real(8), public            :: Ga_tau = 2.267d-12   *GeV     ! tau width
real(8), public            :: m_bot_4gen = 100000d0*GeV     ! b mass for fourth-generation loop in gg MCFM
real(8), public            :: m_top_4gen = 100000d0*GeV     ! t mass for fourth-generation loop in gg MCFM

real(8), public            :: Gf = 1.16639d-5/GeV**2        ! Fermi constant
real(8), public            :: alpha_QED = 1d0/128d0         ! el.magn. coupling
real(8), public            :: vev ! = 1.0d0/sqrt(Gf*sqrt(2.0d0))
real(8), public            :: gwsq ! = 4.0d0 * M_W**2/vev**2  ! weak constant squared
real(dp), public           :: esq ! = 4.0d0 * pi * alpha_QED  ! Fundamental charge

real(8), public            :: xw = 0.23119d0                ! sin**2(Theta_Weinberg) (PDG-2008)
real(8), public            :: sitW ! = dsqrt(xw)            ! sin(Theta_Weinberg) (PDG-2008)
real(8), public            :: twosc ! = sqrt(4.0_dp*xw*(1.0_dp-xw))
real(8), public, parameter :: LHC_Energy=13d0  *TeV      ! LHC hadronic center of mass energy
real(8), public, parameter :: TEV_Energy=1.96d0  *TeV       ! Tevatron hadronic center of mass energy
real(8), public, parameter :: ILC_Energy=250d0  *GeV        ! Linear collider center of mass energy
!command line: epPolarization, emPolarization
real(8), public            :: POL_A = 0d0                   ! e+ polarization. 0: no polarization, 100: helicity = 1, -100: helicity = -1
real(8), public            :: POL_B = 0d0                   ! e- polarization. 0: no polarization, 100: helicity = 1, -100: helicity = -1

! PDF and QCD scale variables, set in main::InitPDFNonConstVals if not a parameter
integer, public, parameter :: nQflavors_pdf = 5    ! Number of flavors enforced to the PDF, used in ModParameters::EvalAlphaS()
integer, public, parameter :: nloops_pdf = 1       ! alpha_s order
real(8), public            :: zmass_pdf            ! Z mass used in pdf toward the QCD scale, reset later in main per PDF if needed
real(8), public            :: Mu_Fact              ! pdf factorization scale (set to M_Reso in main.F90)
real(8), public            :: Mu_Ren               ! QCD renormalization (alpha_s) scale (set to M_Reso in main.F90)
real(dp), public           :: alphas               ! strong coupling per event, set to some reasonable value
real(dp), public           :: alphas_mz            ! strong coupling at M_Z, reset later in main per PDF
real(dp), public           :: gs                   ! = sqrt(alphas*4.0_dp*pi)


! CKM squared matrix entries
real(8), public            :: VCKM_ud! = 0.974285d0
real(8), public            :: VCKM_us! = 0.225290d0
real(8), public            :: VCKM_cs! = 0.9734244d0
real(8), public            :: VCKM_cd! =-0.225182d0
real(8), public            :: VCKM_tb! = 0.99912367d0
real(8), public            :: VCKM_ts! =-0.040920069d0
real(8), public            :: VCKM_cb! = dsqrt(1d0-VCKM_cd**2-VCKM_cs**2)
real(8), public            :: VCKM_ub! = dsqrt(1d0-VCKM_ud**2-VCKM_us**2)
real(8), public            :: VCKM_td! = dsqrt(1d0-VCKM_tb**2-VCKM_ts**2)


! absolute branching fraction (taken from PDG-2014)
real(8), public, parameter :: Br_Z_ll   = 10.10d0*percent                             ! leptonic Z branching
real(8), public, parameter :: Br_Z_hadr = 69.91d0*percent                             ! hadronic Z branching
real(8), public, parameter :: Br_Z_inv  = 100d0*percent - Br_Z_ll - Br_Z_hadr         ! invisible Z branching
real(8), public, parameter :: Br_Z_uu   = 11.6d0*percent                              ! up upbar Z branching
real(8), public, parameter :: Br_Z_cc   = 11.6d0*percent                              ! chm chmbar Z branching
real(8), public, parameter :: Br_Z_dd   = 15.6d0*percent                              ! dn dnbar Z branching
real(8), public, parameter :: Br_Z_ss   = 15.6d0*percent                              ! str strbar Z branching
real(8), public, parameter :: Br_Z_bb   = Br_Z_hadr - Br_Z_uu - Br_Z_dd - Br_Z_cc - Br_Z_ss
real(8), public, parameter :: Br_W_ll   = 32.72d0*percent                             ! leptonic W branching
real(8), public, parameter :: Br_W_hadr = 100d0*percent - Br_W_ll                     ! hadronic W branching


! derived branching fractions
real(8), public, parameter :: Br_Z_ee   = 1d0/3d0*Br_Z_ll                             ! electron Z branching
real(8), public, parameter :: Br_Z_mm   = 1d0/3d0*Br_Z_ll                             ! muon Z branching
real(8), public, parameter :: Br_Z_tt   = 1d0/3d0*Br_Z_ll                             ! tau Z branching
real(8), public, parameter :: Br_Z_nn   = 1d0/3d0*Br_Z_inv                            ! neutrino Z branching
real(8), public, parameter :: Br_W_en   = 1d0/3d0*Br_W_ll                             ! electron W branching
real(8), public, parameter :: Br_W_mn   = 1d0/3d0*Br_W_ll                             ! muon W branching
real(8), public, parameter :: Br_W_tn   = 1d0/3d0*Br_W_ll                             ! electron W branching
real(8), public, parameter :: Br_W_ud   = 1d0/2d0*Br_W_hadr                           ! u-d W branching
real(8), public, parameter :: Br_W_cs   = 1d0/2d0*Br_W_hadr                           ! c-s W branching

real(8), public, parameter :: Brlept_Z_ee = Br_Z_ee/Br_Z_ll                           ! Z branching fraction Ga(el)/Ga(leptonic)
real(8), public, parameter :: Brlept_Z_mm = Br_Z_mm/Br_Z_ll                           ! Z branching fraction Ga(mu)/Ga(leptonic)
real(8), public, parameter :: Brlept_Z_tt = Br_Z_tt/Br_Z_ll                           ! Z branching fraction Ga(tau)/Ga(leptonic)
real(8), public, parameter :: Brlept_Z_nn = Br_Z_nn/Br_Z_inv                          ! Z branching fraction Ga(neu)/Ga(invisible)
real(8), public, parameter :: Brlept_W_en = Br_W_en/Br_W_ll                           ! W branching fraction Ga(el)/Ga(leptonic)
real(8), public, parameter :: Brlept_W_mn = Br_W_mn/Br_W_ll                           ! W branching fraction Ga(mu)/Ga(leptonic)
real(8), public, parameter :: Brlept_W_tn = Br_W_tn/Br_W_ll                           ! W branching fraction Ga(tau)/Ga(leptonic)

real(8), public, parameter :: Brhadr_Z_uu = Br_Z_uu/Br_Z_hadr                         ! Z branching fraction Ga(up)/Ga(hadronic)
real(8), public, parameter :: Brhadr_Z_cc = Br_Z_cc/Br_Z_hadr                         ! Z branching fraction Ga(chm)/Ga(hadronic)
real(8), public, parameter :: Brhadr_Z_dd = Br_Z_dd/Br_Z_hadr                         ! Z branching fraction Ga(don)/Ga(hadronic)
real(8), public, parameter :: Brhadr_Z_ss = Br_Z_ss/Br_Z_hadr                         ! Z branching fraction Ga(str)/Ga(hadronic)
real(8), public, parameter :: Brhadr_Z_bb = Br_Z_bb/Br_Z_hadr                         ! Z branching fraction Ga(bot)/Ga(hadronic)
real(8), public, parameter :: Brhadr_W_ud = Br_W_ud/Br_W_hadr                         ! W branching fraction Ga(up)/Ga(hadronic)
real(8), public, parameter :: Brhadr_W_cs = Br_W_cs/Br_W_hadr                         ! W branching fraction Ga(chm)/Ga(hadronic)

real(8), public :: scale_alpha_Z_uu = 1.037560d0 ! scaling factor of alpha (~partial width) for Z > u u~, c c~
real(8), public :: scale_alpha_Z_dd = 1.037560d0 ! scaling factor of alpha (~partial width) for Z > d d~, s s~, b b~
real(8), public :: scale_alpha_Z_ll = 1d0        ! scaling factor of alpha (~partial width) for Z > l+ l-  (l=e,mu)
real(8), public :: scale_alpha_Z_tt = 1d0        ! scaling factor of alpha (~partial width) for Z > tau+ tau-
real(8), public :: scale_alpha_Z_nn = 1d0        ! scaling factor of alpha (~partial width) for Z > nu nu~
real(8), public :: scale_alpha_W_ud = 1.038200d0 ! scaling factor of alpha (~partial width) for W > u d
real(8), public :: scale_alpha_W_cs = 1.038200d0 ! scaling factor of alpha (~partial width) for W > c s
real(8), public :: scale_alpha_W_ln = 1d0        ! scaling factor of alpha (~partial width) for W > l nu (l=e,mu)
real(8), public :: scale_alpha_W_tn = 1d0        ! scaling factor of alpha (~partial width) for W > tau nu
! real(8), public, parameter :: scale_alpha_Z_ll = (1d0-(Br_Z_uu+Br_Z_cc)*scale_alpha_Z_uu-(Br_Z_dd+Br_Z_ss+Br_Z_bb)*scale_alpha_Z_dd)/(3d0*Br_Z_nn+3d0*Br_Z_ee) ! scaling factor of alpha (~partial width) for Z > l+ l- which restores the total width
! real(8), public, parameter :: scale_alpha_Z_nn = scale_alpha_Z_ll ! scaling factor of alpha (~partial width) for Z > nu nu~ which restores total width
! real(8), public, parameter :: scale_alpha_W_ln = (1d0-Br_W_ud*scale_alpha_W_ud-Br_W_cs*scale_alpha_W_cs)/(3d0*Br_W_en) ! scaling factor of alpha (~partial width) for W > l nu
!
! sum rules
! 1 = 3*Br_Z_nn + 3*Br_Z_ee + 2*Br_Z_uu + 3*Br_Z_dd
!
! sum rule with scaling factors
! 1 = 3*(Br_Z_nn*scale_alpha_Z_nn) + 3*(Br_Z_ee*scale_alpha_Z_ll) + 2*(Br_Z_uu*scale_alpha_Z_uu) + 3*(Br_Z_dd*scale_alpha_Z_dd)
!=====================================================


!=====================================================
!resonance couplings

! Lambda scale enters in two places
! overall scale for x-section and in power suppressed
! operators/formfactors (former r).
real(8), public, parameter :: Lambda  = 1000d0    *GeV
real(8), public, parameter :: Lambda2 = 1000d0    *GeV      ! for second resonance

!--------------------!
!-----! Spin-0 !-----!
!--------------------!
!-- parameters that define on-shell spin 0 coupling to SM fields, see note
   logical, public, parameter :: generate_as = .false. ! .true. uses ah* instead of gh*
   complex(8), public, parameter :: ahg1 = (1.0d0,0d0)
   complex(8), public, parameter :: ahg2 = (0d0,0d0)
   complex(8), public, parameter :: ahg3 = (0d0,0d0)  ! pseudoscalar
   complex(8), public, parameter :: ahz1 = (1.0d0,0d0)
   complex(8), public, parameter :: ahz2 = (0d0,0d0)  ! this coupling does not contribute for gamma+gamma final states
   complex(8), public, parameter :: ahz3 = (0d0,0d0)  ! pseudoscalar

!-- parameters that define off-shell spin 0 coupling to SM fields, see note
!-- Hgg couplings to gluons for point-like vertices
   complex(8), public :: ghg2 = (1.0d0,0d0)
   complex(8), public :: ghg3 = (0d0,0d0)
   complex(8), public :: ghg4 = (0d0,0d0)   ! pseudoscalar

!-- Hgg couplings to gluons for point-like vertices (4th generation vertex for MCFM)
   complex(8), public :: ghg2_4gen = (0d0,0d0)
   complex(8), public :: ghg3_4gen = (0d0,0d0)
   complex(8), public :: ghg4_4gen = (0d0,0d0)

!-- HVV' couplings to ZZ/ZA/AA and WW
   complex(8), public :: ghz1 = (2.0d0,0d0) ! SM=2 (MCFM => =1)
   complex(8), public :: ghz2 = (0d0,0d0)
   complex(8), public :: ghz3 = (0d0,0d0)
   complex(8), public :: ghz4 = (0d0,0d0)   ! pseudoscalar

!-- parameters that define q^2 dependent form factors
   complex(8), public :: ghz1_prime = (0d0,0d0)
   complex(8), public :: ghz1_prime2= (0d0,0d0)
   complex(8), public :: ghz1_prime3= (0d0,0d0)
   complex(8), public :: ghz1_prime4= (0d0,0d0)
   complex(8), public :: ghz1_prime5= (0d0,0d0)
   complex(8), public :: ghz1_prime6= (0d0,0d0)
   complex(8), public :: ghz1_prime7= (0d0,0d0)

   complex(8), public :: ghz2_prime = (0d0,0d0)
   complex(8), public :: ghz2_prime2= (0d0,0d0)
   complex(8), public :: ghz2_prime3= (0d0,0d0)
   complex(8), public :: ghz2_prime4= (0d0,0d0)
   complex(8), public :: ghz2_prime5= (0d0,0d0)
   complex(8), public :: ghz2_prime6= (0d0,0d0)
   complex(8), public :: ghz2_prime7= (0d0,0d0)

   complex(8), public :: ghz3_prime = (0d0,0d0)
   complex(8), public :: ghz3_prime2= (0d0,0d0)
   complex(8), public :: ghz3_prime3= (0d0,0d0)
   complex(8), public :: ghz3_prime4= (0d0,0d0)
   complex(8), public :: ghz3_prime5= (0d0,0d0)
   complex(8), public :: ghz3_prime6= (0d0,0d0)
   complex(8), public :: ghz3_prime7= (0d0,0d0)

   complex(8), public :: ghz4_prime = (0d0,0d0)
   complex(8), public :: ghz4_prime2= (0d0,0d0)
   complex(8), public :: ghz4_prime3= (0d0,0d0)
   complex(8), public :: ghz4_prime4= (0d0,0d0)
   complex(8), public :: ghz4_prime5= (0d0,0d0)
   complex(8), public :: ghz4_prime6= (0d0,0d0)
   complex(8), public :: ghz4_prime7= (0d0,0d0)

   complex(8), public :: ghzgs1_prime2= (0d0,0d0)
   complex(8), public :: ghzgs2  = (0d0,0d0)
   complex(8), public :: ghzgs3  = (0d0,0d0)
   complex(8), public :: ghzgs4  = (0d0,0d0)
   complex(8), public :: ghgsgs2 = (0d0,0d0)
   complex(8), public :: ghgsgs3 = (0d0,0d0)
   complex(8), public :: ghgsgs4 = (0d0,0d0)

   real(8),    public, parameter :: Lambda_z1 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_z2 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_z3 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_z4 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_zgs1 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_Q  = 10000d0*GeV

   integer,    public :: cz_q1sq = 0 ! Sign of q1,2,12**2 for the following Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
   integer,    public :: cz_q2sq = 0
   integer,    public :: cz_q12sq = 0
   ! These Lambdas all have a numerical value of 1d0
   real(8),    public :: Lambda_z11 = 100d0*GeV ! For Z1
   real(8),    public :: Lambda_z21 = 100d0*GeV
   real(8),    public :: Lambda_z31 = 100d0*GeV
   real(8),    public :: Lambda_z41 = 100d0*GeV
   real(8),    public :: Lambda_z12 = 100d0*GeV ! For Z2
   real(8),    public :: Lambda_z22 = 100d0*GeV
   real(8),    public :: Lambda_z32 = 100d0*GeV
   real(8),    public :: Lambda_z42 = 100d0*GeV
   real(8),    public :: Lambda_z10 = 100d0*GeV ! For the Higgs
   real(8),    public :: Lambda_z20 = 100d0*GeV
   real(8),    public :: Lambda_z30 = 100d0*GeV
   real(8),    public :: Lambda_z40 = 100d0*GeV

!-- extra HWW couplings for weak boson fusion when WW-spin-0 couplings are required to be different from ZZ-spin-0
!-- note: ZZ-spin-0 couplings are used in processes other than VBF, and WW is distinguished from ZZ only in case distinguish_HWWcouplings=.true.
   logical, public :: distinguish_HWWcouplings=.false.
   complex(8), public :: ghw1 = (0d0,0d0)
   complex(8), public :: ghw2 = (0d0,0d0)
   complex(8), public :: ghw3 = (0d0,0d0)
   complex(8), public :: ghw4 = (0d0,0d0)

!-- parameters that define q^2 dependent form factors in WBF WW-spin-0 case described above
   complex(8), public :: ghw1_prime = (0d0,0d0)
   complex(8), public :: ghw1_prime2= (0d0,0d0)
   complex(8), public :: ghw1_prime3= (0d0,0d0)
   complex(8), public :: ghw1_prime4= (0d0,0d0)
   complex(8), public :: ghw1_prime5= (0d0,0d0)
   complex(8), public :: ghw1_prime6= (0d0,0d0)
   complex(8), public :: ghw1_prime7= (0d0,0d0)

   complex(8), public :: ghw2_prime = (0d0,0d0)
   complex(8), public :: ghw2_prime2= (0d0,0d0)
   complex(8), public :: ghw2_prime3= (0d0,0d0)
   complex(8), public :: ghw2_prime4= (0d0,0d0)
   complex(8), public :: ghw2_prime5= (0d0,0d0)
   complex(8), public :: ghw2_prime6= (0d0,0d0)
   complex(8), public :: ghw2_prime7= (0d0,0d0)

   complex(8), public :: ghw3_prime = (0d0,0d0)
   complex(8), public :: ghw3_prime2= (0d0,0d0)
   complex(8), public :: ghw3_prime3= (0d0,0d0)
   complex(8), public :: ghw3_prime4= (0d0,0d0)
   complex(8), public :: ghw3_prime5= (0d0,0d0)
   complex(8), public :: ghw3_prime6= (0d0,0d0)
   complex(8), public :: ghw3_prime7= (0d0,0d0)

   complex(8), public :: ghw4_prime = (0d0,0d0)
   complex(8), public :: ghw4_prime2= (0d0,0d0)
   complex(8), public :: ghw4_prime3= (0d0,0d0)
   complex(8), public :: ghw4_prime4= (0d0,0d0)
   complex(8), public :: ghw4_prime5= (0d0,0d0)
   complex(8), public :: ghw4_prime6= (0d0,0d0)
   complex(8), public :: ghw4_prime7= (0d0,0d0)

   real(8),    public, parameter :: Lambda_w1 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_w2 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_w3 = 10000d0*GeV
   real(8),    public, parameter :: Lambda_w4 = 10000d0*GeV

   integer,    public :: cw_q1sq = 0 ! Sign of q1,2,12**2 for the following Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
   integer,    public :: cw_q2sq = 0
   integer,    public :: cw_q12sq = 0
   real(8),    public :: Lambda_w11 = 100d0*GeV ! For W+
   real(8),    public :: Lambda_w21 = 100d0*GeV
   real(8),    public :: Lambda_w31 = 100d0*GeV
   real(8),    public :: Lambda_w41 = 100d0*GeV
   real(8),    public :: Lambda_w12 = 100d0*GeV ! For W-
   real(8),    public :: Lambda_w22 = 100d0*GeV
   real(8),    public :: Lambda_w32 = 100d0*GeV
   real(8),    public :: Lambda_w42 = 100d0*GeV
   real(8),    public :: Lambda_w10 = 100d0*GeV ! For the Higgs
   real(8),    public :: Lambda_w20 = 100d0*GeV
   real(8),    public :: Lambda_w30 = 100d0*GeV
   real(8),    public :: Lambda_w40 = 100d0*GeV

!-- HVV contact terms
   complex(8), public :: ghzzp1 = (0d0,0d0)
   complex(8), public :: ghzzp2 = (0d0,0d0)
   complex(8), public :: ghzzp3 = (0d0,0d0)
   complex(8), public :: ghzzp4 = (0d0,0d0)

!-- parameters that define q^2 dependent form factors
   complex(8), public :: ghzzp1_prime = (0d0,0d0)
   complex(8), public :: ghzzp1_prime2= (0d0,0d0)
   complex(8), public :: ghzzp1_prime3= (0d0,0d0)
   complex(8), public :: ghzzp1_prime4= (0d0,0d0)
   complex(8), public :: ghzzp1_prime5= (0d0,0d0)
   complex(8), public :: ghzzp1_prime6= (0d0,0d0)
   complex(8), public :: ghzzp1_prime7= (0d0,0d0)

   complex(8), public :: ghzzp2_prime = (0d0,0d0)
   complex(8), public :: ghzzp2_prime2= (0d0,0d0)
   complex(8), public :: ghzzp2_prime3= (0d0,0d0)
   complex(8), public :: ghzzp2_prime4= (0d0,0d0)
   complex(8), public :: ghzzp2_prime5= (0d0,0d0)
   complex(8), public :: ghzzp2_prime6= (0d0,0d0)
   complex(8), public :: ghzzp2_prime7= (0d0,0d0)

   complex(8), public :: ghzzp3_prime = (0d0,0d0)
   complex(8), public :: ghzzp3_prime2= (0d0,0d0)
   complex(8), public :: ghzzp3_prime3= (0d0,0d0)
   complex(8), public :: ghzzp3_prime4= (0d0,0d0)
   complex(8), public :: ghzzp3_prime5= (0d0,0d0)
   complex(8), public :: ghzzp3_prime6= (0d0,0d0)
   complex(8), public :: ghzzp3_prime7= (0d0,0d0)

   complex(8), public :: ghzzp4_prime = (0d0,0d0)
   complex(8), public :: ghzzp4_prime2= (0d0,0d0)
   complex(8), public :: ghzzp4_prime3= (0d0,0d0)
   complex(8), public :: ghzzp4_prime4= (0d0,0d0)
   complex(8), public :: ghzzp4_prime5= (0d0,0d0)
   complex(8), public :: ghzzp4_prime6= (0d0,0d0)
   complex(8), public :: ghzzp4_prime7= (0d0,0d0)

!-- Zpgs
   complex(8), public :: ghzpgs1_prime2= (0d0,0d0)
   complex(8), public :: ghzpgs2  = (0d0,0d0)
   complex(8), public :: ghzpgs3  = (0d0,0d0)
   complex(8), public :: ghzpgs4  = (0d0,0d0)

!-- ZpZp
   complex(8), public :: ghzpzp1 = (0d0,0d0)
   complex(8), public :: ghzpzp2 = (0d0,0d0)
   complex(8), public :: ghzpzp3 = (0d0,0d0)
   complex(8), public :: ghzpzp4 = (0d0,0d0)

!-- parameters that define q^2 dependent form factors
   complex(8), public :: ghzpzp1_prime = (0d0,0d0)
   complex(8), public :: ghzpzp1_prime2= (0d0,0d0)
   complex(8), public :: ghzpzp1_prime3= (0d0,0d0)
   complex(8), public :: ghzpzp1_prime4= (0d0,0d0)
   complex(8), public :: ghzpzp1_prime5= (0d0,0d0)
   complex(8), public :: ghzpzp1_prime6= (0d0,0d0)
   complex(8), public :: ghzpzp1_prime7= (0d0,0d0)

   complex(8), public :: ghzpzp2_prime = (0d0,0d0)
   complex(8), public :: ghzpzp2_prime2= (0d0,0d0)
   complex(8), public :: ghzpzp2_prime3= (0d0,0d0)
   complex(8), public :: ghzpzp2_prime4= (0d0,0d0)
   complex(8), public :: ghzpzp2_prime5= (0d0,0d0)
   complex(8), public :: ghzpzp2_prime6= (0d0,0d0)
   complex(8), public :: ghzpzp2_prime7= (0d0,0d0)

   complex(8), public :: ghzpzp3_prime = (0d0,0d0)
   complex(8), public :: ghzpzp3_prime2= (0d0,0d0)
   complex(8), public :: ghzpzp3_prime3= (0d0,0d0)
   complex(8), public :: ghzpzp3_prime4= (0d0,0d0)
   complex(8), public :: ghzpzp3_prime5= (0d0,0d0)
   complex(8), public :: ghzpzp3_prime6= (0d0,0d0)
   complex(8), public :: ghzpzp3_prime7= (0d0,0d0)

   complex(8), public :: ghzpzp4_prime = (0d0,0d0)
   complex(8), public :: ghzpzp4_prime2= (0d0,0d0)
   complex(8), public :: ghzpzp4_prime3= (0d0,0d0)
   complex(8), public :: ghzpzp4_prime4= (0d0,0d0)
   complex(8), public :: ghzpzp4_prime5= (0d0,0d0)
   complex(8), public :: ghzpzp4_prime6= (0d0,0d0)
   complex(8), public :: ghzpzp4_prime7= (0d0,0d0)


   complex(8), public :: ezp_El_left  = (0d0,0d0)
   complex(8), public :: ezp_El_right  = (0d0,0d0)
   complex(8), public :: ezp_Mu_left  = (0d0,0d0)
   complex(8), public :: ezp_Mu_right  = (0d0,0d0)
   complex(8), public :: ezp_Ta_left  = (0d0,0d0)
   complex(8), public :: ezp_Ta_right  = (0d0,0d0)
   complex(8), public :: ezp_NuE_left  = (0d0,0d0)   !same for NuMu and NuTau
   complex(8), public :: ezp_NuE_right  = (0d0,0d0)  !same for NuMu and NuTau
   complex(8), public :: ezp_Up_left  = (0d0,0d0)
   complex(8), public :: ezp_Up_right  = (0d0,0d0)
   complex(8), public :: ezp_Chm_left  = (0d0,0d0)
   complex(8), public :: ezp_Chm_right  = (0d0,0d0)
   complex(8), public :: ezp_Dn_left  = (0d0,0d0)
   complex(8), public :: ezp_Dn_right  = (0d0,0d0)
   complex(8), public :: ezp_Str_left  = (0d0,0d0)
   complex(8), public :: ezp_Str_right  = (0d0,0d0)
   complex(8), public :: ezp_Bot_left  = (0d0,0d0)
   complex(8), public :: ezp_Bot_right  = (0d0,0d0)
   complex(8), public :: ezp_Top_left  = (0d0,0d0)
   complex(8), public :: ezp_Top_right  = (0d0,0d0)

   real(8), public :: M_Zprime = -1d0 ! <0: CT interaction, >=0: Heavy Zprime propagator
   real(8), public :: Ga_Zprime = 0d0
!--
!-- HVV contact terms
   complex(8), public :: ghwwp1 = (0d0,0d0)
   complex(8), public :: ghwwp2 = (0d0,0d0)
   complex(8), public :: ghwwp3 = (0d0,0d0)
   complex(8), public :: ghwwp4 = (0d0,0d0)

!-- parameters that define q^2 dependent form factors
   complex(8), public :: ghwwp1_prime = (0d0,0d0)
   complex(8), public :: ghwwp1_prime2= (0d0,0d0)
   complex(8), public :: ghwwp1_prime3= (0d0,0d0)
   complex(8), public :: ghwwp1_prime4= (0d0,0d0)
   complex(8), public :: ghwwp1_prime5= (0d0,0d0)
   complex(8), public :: ghwwp1_prime6= (0d0,0d0)
   complex(8), public :: ghwwp1_prime7= (0d0,0d0)

   complex(8), public :: ghwwp2_prime = (0d0,0d0)
   complex(8), public :: ghwwp2_prime2= (0d0,0d0)
   complex(8), public :: ghwwp2_prime3= (0d0,0d0)
   complex(8), public :: ghwwp2_prime4= (0d0,0d0)
   complex(8), public :: ghwwp2_prime5= (0d0,0d0)
   complex(8), public :: ghwwp2_prime6= (0d0,0d0)
   complex(8), public :: ghwwp2_prime7= (0d0,0d0)

   complex(8), public :: ghwwp3_prime = (0d0,0d0)
   complex(8), public :: ghwwp3_prime2= (0d0,0d0)
   complex(8), public :: ghwwp3_prime3= (0d0,0d0)
   complex(8), public :: ghwwp3_prime4= (0d0,0d0)
   complex(8), public :: ghwwp3_prime5= (0d0,0d0)
   complex(8), public :: ghwwp3_prime6= (0d0,0d0)
   complex(8), public :: ghwwp3_prime7= (0d0,0d0)

   complex(8), public :: ghwwp4_prime = (0d0,0d0)
   complex(8), public :: ghwwp4_prime2= (0d0,0d0)
   complex(8), public :: ghwwp4_prime3= (0d0,0d0)
   complex(8), public :: ghwwp4_prime4= (0d0,0d0)
   complex(8), public :: ghwwp4_prime5= (0d0,0d0)
   complex(8), public :: ghwwp4_prime6= (0d0,0d0)
   complex(8), public :: ghwwp4_prime7= (0d0,0d0)

   complex(8), public :: ghwpwp1 = (0d0,0d0)
   complex(8), public :: ghwpwp2 = (0d0,0d0)
   complex(8), public :: ghwpwp3 = (0d0,0d0)
   complex(8), public :: ghwpwp4 = (0d0,0d0)

!-- parameters that define q^2 dependent form factors
   complex(8), public :: ghwpwp1_prime = (0d0,0d0)
   complex(8), public :: ghwpwp1_prime2= (0d0,0d0)
   complex(8), public :: ghwpwp1_prime3= (0d0,0d0)
   complex(8), public :: ghwpwp1_prime4= (0d0,0d0)
   complex(8), public :: ghwpwp1_prime5= (0d0,0d0)
   complex(8), public :: ghwpwp1_prime6= (0d0,0d0)
   complex(8), public :: ghwpwp1_prime7= (0d0,0d0)

   complex(8), public :: ghwpwp2_prime = (0d0,0d0)
   complex(8), public :: ghwpwp2_prime2= (0d0,0d0)
   complex(8), public :: ghwpwp2_prime3= (0d0,0d0)
   complex(8), public :: ghwpwp2_prime4= (0d0,0d0)
   complex(8), public :: ghwpwp2_prime5= (0d0,0d0)
   complex(8), public :: ghwpwp2_prime6= (0d0,0d0)
   complex(8), public :: ghwpwp2_prime7= (0d0,0d0)

   complex(8), public :: ghwpwp3_prime = (0d0,0d0)
   complex(8), public :: ghwpwp3_prime2= (0d0,0d0)
   complex(8), public :: ghwpwp3_prime3= (0d0,0d0)
   complex(8), public :: ghwpwp3_prime4= (0d0,0d0)
   complex(8), public :: ghwpwp3_prime5= (0d0,0d0)
   complex(8), public :: ghwpwp3_prime6= (0d0,0d0)
   complex(8), public :: ghwpwp3_prime7= (0d0,0d0)

   complex(8), public :: ghwpwp4_prime = (0d0,0d0)
   complex(8), public :: ghwpwp4_prime2= (0d0,0d0)
   complex(8), public :: ghwpwp4_prime3= (0d0,0d0)
   complex(8), public :: ghwpwp4_prime4= (0d0,0d0)
   complex(8), public :: ghwpwp4_prime5= (0d0,0d0)
   complex(8), public :: ghwpwp4_prime6= (0d0,0d0)
   complex(8), public :: ghwpwp4_prime7= (0d0,0d0)

   complex(8), public :: ewp_El_left  = (0d0,0d0)
   complex(8), public :: ewp_El_right  = (0d0,0d0)
   complex(8), public :: ewp_Mu_left  = (0d0,0d0)
   complex(8), public :: ewp_Mu_right  = (0d0,0d0)
   complex(8), public :: ewp_Ta_left  = (0d0,0d0)
   complex(8), public :: ewp_Ta_right  = (0d0,0d0)
   complex(8), public :: ewp_Up_left  = (0d0,0d0)
   complex(8), public :: ewp_Up_right  = (0d0,0d0)
   complex(8), public :: ewp_Chm_left  = (0d0,0d0)
   complex(8), public :: ewp_Chm_right  = (0d0,0d0)
   complex(8), public :: ewp_Top_left  = (0d0,0d0)
   complex(8), public :: ewp_Top_right  = (0d0,0d0)

   real(8), public :: M_Wprime = -1d0 ! <0: CT interaction, >=0: Heavy Zprime propagator
   real(8), public :: Ga_Wprime = 0d0
!--




!-- second resonance (H2) couplings for MCFM interface
!-- Hgg couplings to gluons for point-like vertices
   complex(8), public :: gh2g2 = (0d0,0d0)
   complex(8), public :: gh2g3 = (0d0,0d0)
   complex(8), public :: gh2g4 = (0d0,0d0)   ! pseudoscalar

!-- Hgg couplings to gluons for point-like vertices (4th generation vertex for MCFM)
   complex(8), public :: gh2g2_4gen = (0d0,0d0)
   complex(8), public :: gh2g3_4gen = (0d0,0d0)
   complex(8), public :: gh2g4_4gen = (0d0,0d0)

!-- HVV' couplings to ZZ/ZA/AA and WW
   complex(8), public :: gh2z1 = (0.0d0,0d0)
   complex(8), public :: gh2z2 = (0.0d0,0d0)
   complex(8), public :: gh2z3 = (0.0d0,0d0)
   complex(8), public :: gh2z4 = (0.0d0,0d0)   ! pseudoscalar

   complex(8), public :: gh2zgs2  = (0.00d0,0d0)
   complex(8), public :: gh2zgs3  = (0.00d0,0d0)
   complex(8), public :: gh2zgs4  = (0.00d0,0d0)
   complex(8), public :: gh2gsgs2 = (0.00d0,0d0)
   complex(8), public :: gh2gsgs3 = (0.00d0,0d0)
   complex(8), public :: gh2gsgs4 = (0.00d0,0d0)

!-- parameters that define q^2 dependent form factors
   complex(8), public :: gh2z1_prime = (0.0d0,0d0)
   complex(8), public :: gh2z1_prime2= (0.0d0,0d0)
   complex(8), public :: gh2z1_prime3= (0.0d0,0d0)
   complex(8), public :: gh2z1_prime4= (0.0d0,0d0)
   complex(8), public :: gh2z1_prime5= (0.0d0,0d0)
   complex(8), public :: gh2z1_prime6= (0.0d0,0d0)
   complex(8), public :: gh2z1_prime7= (0.0d0,0d0)

   complex(8), public :: gh2z2_prime = (0.0d0,0d0)
   complex(8), public :: gh2z2_prime2= (0.0d0,0d0)
   complex(8), public :: gh2z2_prime3= (0.0d0,0d0)
   complex(8), public :: gh2z2_prime4= (0.0d0,0d0)
   complex(8), public :: gh2z2_prime5= (0.0d0,0d0)
   complex(8), public :: gh2z2_prime6= (0.0d0,0d0)
   complex(8), public :: gh2z2_prime7= (0.0d0,0d0)

   complex(8), public :: gh2z3_prime = (0.0d0,0d0)
   complex(8), public :: gh2z3_prime2= (0.0d0,0d0)
   complex(8), public :: gh2z3_prime3= (0.0d0,0d0)
   complex(8), public :: gh2z3_prime4= (0.0d0,0d0)
   complex(8), public :: gh2z3_prime5= (0.0d0,0d0)
   complex(8), public :: gh2z3_prime6= (0.0d0,0d0)
   complex(8), public :: gh2z3_prime7= (0.0d0,0d0)

   complex(8), public :: gh2z4_prime = (0.0d0,0d0)
   complex(8), public :: gh2z4_prime2= (0.0d0,0d0)
   complex(8), public :: gh2z4_prime3= (0.0d0,0d0)
   complex(8), public :: gh2z4_prime4= (0.0d0,0d0)
   complex(8), public :: gh2z4_prime5= (0.0d0,0d0)
   complex(8), public :: gh2z4_prime6= (0.0d0,0d0)
   complex(8), public :: gh2z4_prime7= (0.0d0,0d0)

   complex(8), public :: gh2zgs1_prime2= (0.0d0,0d0)

   real(8),    public, parameter :: Lambda2_z1 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_z2 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_z3 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_z4 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_zgs1 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_Q  = 10000d0*GeV

   integer,    public :: c2z_q1sq = 0 ! Sign of q1,2,12**2 for the following Lambda2's, set to 1 or -1 to get q**2-dependence from these form factor Lambda2s
   integer,    public :: c2z_q2sq = 0
   integer,    public :: c2z_q12sq = 0
   real(8),    public :: Lambda2_z11 = 100d0*GeV ! For Z1
   real(8),    public :: Lambda2_z21 = 100d0*GeV
   real(8),    public :: Lambda2_z31 = 100d0*GeV
   real(8),    public :: Lambda2_z41 = 100d0*GeV
   real(8),    public :: Lambda2_z12 = 100d0*GeV ! For Z2
   real(8),    public :: Lambda2_z22 = 100d0*GeV
   real(8),    public :: Lambda2_z32 = 100d0*GeV
   real(8),    public :: Lambda2_z42 = 100d0*GeV
   real(8),    public :: Lambda2_z10 = 100d0*GeV ! For the Higgs
   real(8),    public :: Lambda2_z20 = 100d0*GeV
   real(8),    public :: Lambda2_z30 = 100d0*GeV
   real(8),    public :: Lambda2_z40 = 100d0*GeV


!-- extra HWW couplings for weak boson fusion when WW-spin-0 couplings are required to be different from ZZ-spin-0
!-- note: ZZ-spin-0 couplings are used in processes other than VBF, and WW is distinguished from ZZ only in case distinguish_HWWcouplings=.true.
!    logical, public :: distinguish_HWWcouplings=.false.
   complex(8), public :: gh2w1 = (0.0d0,0d0)
   complex(8), public :: gh2w2 = (0.0d0,0d0)
   complex(8), public :: gh2w3 = (0.0d0,0d0)
   complex(8), public :: gh2w4 = (0.0d0,0d0)

!-- parameters that define q^2 dependent form factors in WBF WW-spin-0 case described above
   complex(8), public :: gh2w1_prime = (0.0d0,0d0)
   complex(8), public :: gh2w1_prime2= (0.0d0,0d0)
   complex(8), public :: gh2w1_prime3= (0.0d0,0d0)
   complex(8), public :: gh2w1_prime4= (0.0d0,0d0)
   complex(8), public :: gh2w1_prime5= (0.0d0,0d0)
   complex(8), public :: gh2w1_prime6= (0.0d0,0d0)
   complex(8), public :: gh2w1_prime7= (0.0d0,0d0)

   complex(8), public :: gh2w2_prime = (0.0d0,0d0)
   complex(8), public :: gh2w2_prime2= (0.0d0,0d0)
   complex(8), public :: gh2w2_prime3= (0.0d0,0d0)
   complex(8), public :: gh2w2_prime4= (0.0d0,0d0)
   complex(8), public :: gh2w2_prime5= (0.0d0,0d0)
   complex(8), public :: gh2w2_prime6= (0.0d0,0d0)
   complex(8), public :: gh2w2_prime7= (0.0d0,0d0)

   complex(8), public :: gh2w3_prime = (0.0d0,0d0)
   complex(8), public :: gh2w3_prime2= (0.0d0,0d0)
   complex(8), public :: gh2w3_prime3= (0.0d0,0d0)
   complex(8), public :: gh2w3_prime4= (0.0d0,0d0)
   complex(8), public :: gh2w3_prime5= (0.0d0,0d0)
   complex(8), public :: gh2w3_prime6= (0.0d0,0d0)
   complex(8), public :: gh2w3_prime7= (0.0d0,0d0)

   complex(8), public :: gh2w4_prime = (0.0d0,0d0)
   complex(8), public :: gh2w4_prime2= (0.0d0,0d0)
   complex(8), public :: gh2w4_prime3= (0.0d0,0d0)
   complex(8), public :: gh2w4_prime4= (0.0d0,0d0)
   complex(8), public :: gh2w4_prime5= (0.0d0,0d0)
   complex(8), public :: gh2w4_prime6= (0.0d0,0d0)
   complex(8), public :: gh2w4_prime7= (0.0d0,0d0)

   real(8),    public, parameter :: Lambda2_w1 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_w2 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_w3 = 10000d0*GeV
   real(8),    public, parameter :: Lambda2_w4 = 10000d0*GeV

   integer,    public :: c2w_q1sq = 0 ! Sign of q1,2,12**2 for the following Lambda2's, set to 1 or -1 to get q**2-dependence from these form factor Lambda2s
   integer,    public :: c2w_q2sq = 0
   integer,    public :: c2w_q12sq = 0
   real(8),    public :: Lambda2_w11 = 100d0*GeV ! For W+
   real(8),    public :: Lambda2_w21 = 100d0*GeV
   real(8),    public :: Lambda2_w31 = 100d0*GeV
   real(8),    public :: Lambda2_w41 = 100d0*GeV
   real(8),    public :: Lambda2_w12 = 100d0*GeV ! For W-
   real(8),    public :: Lambda2_w22 = 100d0*GeV
   real(8),    public :: Lambda2_w32 = 100d0*GeV
   real(8),    public :: Lambda2_w42 = 100d0*GeV
   real(8),    public :: Lambda2_w10 = 100d0*GeV ! For the Higgs
   real(8),    public :: Lambda2_w20 = 100d0*GeV
   real(8),    public :: Lambda2_w30 = 100d0*GeV
   real(8),    public :: Lambda2_w40 = 100d0*GeV


!-- Anomalous coupling parameters in triple and quartic gauge couplings
!     notation wrt paper: dV_N=d^{NWW}*d_1^N, dP_N=d^{NWW}*d_2^N, dM_N=d^{NWW}*d_3^N, dFour_N=d^{NWW}*d_4^N, for N=Z,\gamma
!                         dZZWpWm = d^{ZZWW}*sw^2/cw^2, dZAWpWm = d^{Z\gammaWW}*sw/cw, dAAWpWm = d^{\gamma\gammaWW}
   complex(8),    public :: dV_A= (1.0d0,0d0)
   complex(8),    public :: dP_A= (1.0d0,0d0)
   complex(8),    public :: dM_A= (1.0d0,0d0)
   complex(8),    public :: dV_Z= (1.0d0,0d0)
   complex(8),    public :: dP_Z= (1.0d0,0d0)
   complex(8),    public :: dM_Z= (1.0d0,0d0)
   complex(8),    public :: dAAWpWm= (1.0d0,0d0)
   complex(8),    public :: dZAWpWm= (1.0d0,0d0)
   complex(8),    public :: dZZWpWm= (1.0d0,0d0)
   complex(8),    public :: dFour_A= (0.0d0,0d0)
   complex(8),    public :: dFour_Z= (0.0d0,0d0)


!-- Hff couplings for ttbar+H and bbar+H
   complex(8), public :: kappa       = (1d0,0d0)
   complex(8), public :: kappa_tilde = (0d0,0d0)
!-- first resonance (H) couplings
   ! Usable in JHUGen-only processes
   complex(8), public :: kappa_top       = (1d0,0d0)
   complex(8), public :: kappa_tilde_top = (0d0,0d0)
   complex(8), public :: kappa_bot       = (1d0,0d0)
   complex(8), public :: kappa_tilde_bot = (0d0,0d0)
   ! Extra MCFM couplings to fourth-generation quarks
   complex(8), public :: kappa_4gen_top       = (0d0,0d0)
   complex(8), public :: kappa_tilde_4gen_top = (0d0,0d0)
   complex(8), public :: kappa_4gen_bot       = (0d0,0d0)
   complex(8), public :: kappa_tilde_4gen_bot = (0d0,0d0)

!-- seconds resonance (H2) couplings for MCFM interface
   complex(8), public :: kappa2_top       = (0d0,0d0)
   complex(8), public :: kappa2_tilde_top = (0d0,0d0)
   complex(8), public :: kappa2_bot       = (0d0,0d0)
   complex(8), public :: kappa2_tilde_bot = (0d0,0d0)

   complex(8), public :: kappa2_4gen_top       = (0d0,0d0)
   complex(8), public :: kappa2_tilde_4gen_top = (0d0,0d0)
   complex(8), public :: kappa2_4gen_bot       = (0d0,0d0)
   complex(8), public :: kappa2_tilde_4gen_bot = (0d0,0d0)


!--------------------!
!-----! Spin-1 !-----!
!--------------------!
!---parameters that define spin 1 coupling to SM fields, see note
   complex(8), public :: zprime_qq_left  = (1.0d0,0d0)
   complex(8), public :: zprime_qq_right = (1.0d0,0d0)
   complex(8), public :: zprime_zz_1 =  (0d0,0d0)!  =1 for JP=1- vector
   complex(8), public :: zprime_zz_2 =  (0d0,0d0)!  =1 for JP=1+ pseudovector


!--------------------!
!-----! Spin-2 !-----!
!--------------------!
!-- parameters that define spin 2 coupling to SM fields, see note
! minimal coupling corresponds to a1 = b1 = b5 = 1 everything else 0
  complex(8), public :: a1 = (0d0,0d0)    ! g1  -- c.f. draft
  complex(8), public :: a2 = (0d0,0d0)    ! g2
  complex(8), public :: a3 = (0d0,0d0)    ! g3
  complex(8), public :: a4 = (0d0,0d0)    ! g4
  complex(8), public :: a5 = (0d0,0d0)    ! pseudoscalar, g8
  complex(8), public :: graviton_qq_left  = (1.0d0,0d0)! graviton coupling to quarks
  complex(8), public :: graviton_qq_right = (1.0d0,0d0)

!-- see mod_Graviton for these two parameters
  logical, public, parameter :: generate_bis = .true.
  logical, public, parameter :: use_dynamic_MG = .true.

  !  all b' below are g's in the draft
  complex(8), public :: b1 = (0d0,0d0)
  complex(8), public :: b2 = (0d0,0d0)
  complex(8), public :: b3 = (0d0,0d0)
  complex(8), public :: b4 = (0d0,0d0)
  complex(8), public :: b5 = (0d0,0d0)  ! this coupling does not contribute to V+gamma final states
  complex(8), public :: b6 = (0d0,0d0)  ! this coupling does not contribute to V+gamma final states
  complex(8), public :: b7 = (0d0,0d0)  ! this coupling does not contribute to V+gamma final states
  complex(8), public :: b8 = (0d0,0d0)
  complex(8), public :: b9 = (0d0,0d0)  ! this coupling does not contribute to V+gamma final states
  complex(8), public :: b10 =(0d0,0d0)  ! this coupling does not contribute to V+gamma final states

  complex(8), public :: bzzp1 = (0d0,0d0)
  complex(8), public :: bzzp2 = (0d0,0d0)
  complex(8), public :: bzzp3 = (0d0,0d0)
  complex(8), public :: bzzp4 = (0d0,0d0)
  complex(8), public :: bzzp5 = (0d0,0d0)
  complex(8), public :: bzzp6 = (0d0,0d0)
  complex(8), public :: bzzp7 = (0d0,0d0)
  complex(8), public :: bzzp8 = (0d0,0d0)
  complex(8), public :: bzzp9 = (0d0,0d0)
  complex(8), public :: bzzp10 =(0d0,0d0)

  complex(8), public :: bzpzp1 = (0d0,0d0)
  complex(8), public :: bzpzp2 = (0d0,0d0)
  complex(8), public :: bzpzp3 = (0d0,0d0)
  complex(8), public :: bzpzp4 = (0d0,0d0)
  complex(8), public :: bzpzp5 = (0d0,0d0)
  complex(8), public :: bzpzp6 = (0d0,0d0)
  complex(8), public :: bzpzp7 = (0d0,0d0)
  complex(8), public :: bzpzp8 = (0d0,0d0)
  complex(8), public :: bzpzp9 = (0d0,0d0)
  complex(8), public :: bzpzp10 =(0d0,0d0)

  complex(8), public :: bzgs1 = (0d0,0d0)
  complex(8), public :: bzgs2 = (0d0,0d0)
  complex(8), public :: bzgs3 = (0d0,0d0)
  complex(8), public :: bzgs4 = (0d0,0d0)
  complex(8), public :: bzgs8 = (0d0,0d0)

  complex(8), public :: bgsgs1 = (0d0,0d0)
  complex(8), public :: bgsgs2 = (0d0,0d0)
  complex(8), public :: bgsgs3 = (0d0,0d0)
  complex(8), public :: bgsgs4 = (0d0,0d0)
  complex(8), public :: bgsgs8 = (0d0,0d0)

  complex(8), public :: bzpgs1 = (0d0,0d0)
  complex(8), public :: bzpgs2 = (0d0,0d0)
  complex(8), public :: bzpgs3 = (0d0,0d0)
  complex(8), public :: bzpgs4 = (0d0,0d0)
  complex(8), public :: bzpgs8 = (0d0,0d0)


  complex(8), public, parameter  :: c1 = (1.0d0,0d0)
  complex(8), public, parameter  :: c2 = (0d0,0d0)
  complex(8), public, parameter  :: c3 = (0d0,0d0)
  complex(8), public, parameter  :: c41= (0d0,0d0)
  complex(8), public, parameter  :: c42= (0d0,0d0)
  complex(8), public, parameter  :: c5 = (0d0,0d0)
  complex(8), public, parameter  :: c6 = (0d0,0d0) ! this coupling does not contribute to gamma+gamma final states
  complex(8), public, parameter  :: c7 = (0d0,0d0) ! this coupling does not contribute to gamma+gamma final states


! user Zff couplings
   real(8), public :: daz_Lep_left  = (0d0)
   real(8), public :: daz_Lep_right  = (0d0)
   real(8), public :: daz_Neu_left  = (0d0)
   real(8), public :: daz_Neu_right  = (0d0)
   real(8), public :: daz_QUp_left  = (0d0)
   real(8), public :: daz_QUp_right  = (0d0)
   real(8), public :: daz_QDn_left  = (0d0)
   real(8), public :: daz_QDn_right  = (0d0)

!=====================================================

!=====================================================
!internal
!---     B0_PDF=(11.-2.*NF/3.)/4./PI
real(dp), public, parameter :: B0_PDF(0:6) = (/ 0.8753521870054244D0,0.822300539308126D0,0.7692488916108274D0,0.716197243913529D0,0.6631455962162306D0,0.6100939485189321D0,0.5570423008216338D0 /)
real(dp), public, parameter :: nf = 5.0_dp
real(dp), public, parameter :: xn = 3.0_dp
real(dp), public, parameter :: Ca = 3.0_dp
real(dp), public, parameter :: Cf = 4.0_dp/3.0_dp
real(dp), public, parameter :: Tr = 0.5_dp
real(dp), public, parameter :: avegg = 1.0_dp/4.0_dp/64.0_dp
real(dp), public, parameter :: aveqg = 1.0_dp/4.0_dp/24.0_dp
real(dp), public, parameter :: aveqq = 1.0_dp/4.0_dp/9.0_dp

! Particle isospin and charges
real(8), public, parameter :: T3lL= -0.5d0
real(8), public, parameter :: T3lR=  0d0
real(8), public, parameter :: T3nL=  0.5d0
real(8), public, parameter :: T3nR=  0d0
real(8), public, parameter :: T3uL= 0.5d0
real(8), public, parameter :: T3uR= 0d0
real(8), public, parameter :: T3dL= -0.5d0
real(8), public, parameter :: T3dR= 0d0
real(8), public, parameter :: QlL = -1d0
real(8), public, parameter :: QlR = -1d0
real(8), public, parameter :: QnL =  0d0
real(8), public, parameter :: QnR =  0d0
real(8), public, parameter :: QuL = 2d0/3d0
real(8), public, parameter :: QuR = 2d0/3d0
real(8), public, parameter :: QdL = -1d0/3d0
real(8), public, parameter :: QdR = -1d0/3d0


! V-f-fbar couplings:
!   g_R(f) = -e*sw/cw*Q(f)                 = e/2/sw/cw * a(b,c)R,
!   g_L(f) = -e*sw/cw*Q(f) + e/sw/cw*T3(f) = e/2/sw/cw * a(b,c)L
! with
!   aR(f) = -2*sw**2*Q(f),
!   aL(f) = -2*sw**2*Q(f) + 2*T3(f).
! for V = Z-boson,
!   bR = 0
!   bL = dsqrt(2)*cw
! for V = W-boson,
! and
!   cR = -2*sw*cw*Q(f)
!   cL = -2*sw*cw*Q(f)
! for V = photon*
!
real(8), public            :: overallCouplVffsq ! Overall coupling squared that goes along with the ones below
real(8), public            :: aR_lep ! = -2d0*sitW**2*(-1d0)
real(8), public            :: aL_lep ! = -2d0*sitW**2*(-1d0)-1d0
real(8), public            :: aR_neu ! = -2d0*sitW**2*(0d0)
real(8), public            :: aL_neu ! = -2d0*sitW**2*(0d0)+1d0
real(8), public            :: aR_QUp ! = -2d0*sitW**2*(2d0/3d0)
real(8), public            :: aL_QUp ! = -2d0*sitW**2*(2d0/3d0)+1d0
real(8), public            :: aR_QDn ! = -2d0*sitW**2*(-1d0/3d0)
real(8), public            :: aL_QDn ! = -2d0*sitW**2*(-1d0/3d0)-1d0
real(8), public            :: bL ! = dsqrt(2d0)*dsqrt(1d0-sitW**2)
real(8), public            :: bR ! = 0d0
real(8), public            :: cR_lep ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0)
real(8), public            :: cL_lep ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0)
real(8), public            :: cR_neu ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(0d0)
real(8), public            :: cL_neu ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(0d0)
real(8), public            :: cR_QUp ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(2d0/3d0)
real(8), public            :: cL_QUp ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(2d0/3d0)
real(8), public            :: cR_QDn ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0/3d0)
real(8), public            :: cL_QDn ! = -2d0*sitW*dsqrt(1d0-sitW**2)*(-1d0/3d0)

! Coupling normalizations based on the common factor to be gwsq
real(dp), public           :: couplWffsq
real(dp), public           :: couplZffsq
real(dp), public           :: couplAZff
real(dp), public           :: couplAffsq
!----------------------------------------------------------------------------------

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
integer, public, target :: NuE_ = 14
integer, public, target :: NuM_ = 15
integer, public, target :: NuT_ = 16
integer, public, target :: Hig_ = 25
integer, public, target :: Zpr_ = 32
integer, public, target :: Zpr2_ = 33
integer, public, target :: Wppr_ = 34
integer, public, target :: Gra_ = 39

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
integer, public, target :: ANuE_ = -14
integer, public, target :: ANuM_ = -15
integer, public, target :: ANuT_ = -16
integer, public, target :: Wmpr_ = -34

integer, public, parameter :: Not_a_particle_  = -9000
real(8), public, parameter :: Mom_Not_a_particle(1:4) = (/0d0,0d0,0d0,0d0/)

integer, public, parameter :: pdfGlu_ = 0
integer, public, parameter :: pdfDn_ = 1
integer, public, parameter :: pdfUp_ = 2
integer, public, parameter :: pdfStr_ = 3
integer, public, parameter :: pdfChm_ = 4
integer, public, parameter :: pdfBot_ = 5
integer, public, parameter :: pdfTop_ = 6 ! Dummy
integer, public, parameter :: pdfADn_ = -1
integer, public, parameter :: pdfAUp_ = -2
integer, public, parameter :: pdfAStr_ = -3
integer, public, parameter :: pdfAChm_ = -4
integer, public, parameter :: pdfABot_ = -5
integer, public, parameter :: pdfATop_ = -6 ! Dummy

integer,parameter :: io_stdout=6
integer,parameter :: io_LHEOutFile=14
integer,parameter :: io_HistoFile=15
integer,parameter :: io_LHEInFile=16
integer,parameter :: io_LogFile=17
integer,parameter :: io_CSmaxFile=18
integer,parameter :: io_LHEOutFile2=19
integer,parameter :: io_LHEOutFile3=20
integer,parameter :: io_TmpFile=21   !to use for whatever purpose is needed, but close afterwards

integer, public :: DebugCounter(0:10) = 0
real(8), public :: debugvar(0:10) = 0d0


integer, public :: ijPartons(1:2)=0

!=====================================================


CONTAINS


function HVVSpinZeroDynamicCoupling (index,sWplus,sWminus,sWW,tryWWcoupl)
integer, intent(in) :: index
real(8), intent(in) :: sWplus, sWminus, sWW
real(8) :: sWplus_signed, sWminus_signed, sWW_signed, QsqCompoundFactor
logical,optional :: tryWWcoupl
complex(8) :: HVVSpinZeroDynamicCoupling
complex(8) :: vvcoupl(1:8)
real(8) :: lambda_v
real(8) :: lambda_v120(1:3)
logical :: forceZZcoupl
logical :: computeQsqCompundCoupl

   if(present(tryWWcoupl)) then
      forceZZcoupl = (.not.tryWWcoupl .or. .not.distinguish_HWWcouplings .or. (index.gt.4 .and. index.lt.12))
   else
      forceZZcoupl = .true.
   endif
   computeQsqCompundCoupl = .false.
   sWplus_signed=0d0
   sWminus_signed=0d0
   sWW_signed=0d0
   vvcoupl(:)=czero
   HVVSpinZeroDynamicCoupling=czero
   if( forceZZcoupl ) then
      if(cz_q1sq.ne.0) sWplus_signed=abs(sWplus)*dble(sign(1,cz_q1sq))
      if(cz_q2sq.ne.0) sWminus_signed=abs(sWminus)*dble(sign(1,cz_q2sq))
      if(cz_q12sq.ne.0) sWW_signed=abs(sWW)*dble(sign(1,cz_q12sq))
      if(cz_q1sq.ne.0 .or. cz_q2sq.ne.0 .or. cz_q12sq.ne.0) computeQsqCompundCoupl=.true.
      if(index.eq.1) then ! ZZ 1-4
         vvcoupl = (/ ghz1, ghz1_prime, ghz1_prime2, ghz1_prime3, ghz1_prime4, ghz1_prime5, ghz1_prime6, ghz1_prime7 /)
         lambda_v = Lambda_z1
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.2) then
         vvcoupl = (/ ghz2, ghz2_prime, ghz2_prime2, ghz2_prime3, ghz2_prime4, ghz2_prime5, ghz2_prime6, ghz2_prime7 /)
         lambda_v = Lambda_z2
         lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
      elseif(index.eq.3) then
         vvcoupl = (/ ghz3, ghz3_prime, ghz3_prime2, ghz3_prime3, ghz3_prime4, ghz3_prime5, ghz3_prime6, ghz3_prime7 /)
         lambda_v = Lambda_z3
         lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
      elseif(index.eq.4) then
         vvcoupl = (/ ghz4, ghz4_prime, ghz4_prime2, ghz4_prime3, ghz4_prime4, ghz4_prime5, ghz4_prime6, ghz4_prime7 /)
         lambda_v = Lambda_z4
         lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
      elseif(index.eq.5) then ! Zgs 1
         vvcoupl(3) = ghzgs1_prime2
         lambda_v = Lambda_zgs1
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.6) then ! Zgs 2-4
         vvcoupl(1) = ghzgs2
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
      elseif(index.eq.7) then
         vvcoupl(1) = ghzgs3
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
      elseif(index.eq.8) then
         vvcoupl(1) = ghzgs4
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
      elseif(index.eq.9) then ! gsgs 2-4
         vvcoupl(1) = ghgsgs2
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
      elseif(index.eq.10) then
         vvcoupl(1) = ghgsgs3
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
      elseif(index.eq.11) then
         vvcoupl(1) = ghgsgs4
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
      elseif(index.eq.12) then ! ZpZ 1-4
         vvcoupl = (/ ghzzp1, ghzzp1_prime, ghzzp1_prime2, ghzzp1_prime3, ghzzp1_prime4, ghzzp1_prime5, ghzzp1_prime6, ghzzp1_prime7 /)
         lambda_v = Lambda_z1
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.13) then
         vvcoupl = (/ ghzzp2, ghzzp2_prime, ghzzp2_prime2, ghzzp2_prime3, ghzzp2_prime4, ghzzp2_prime5, ghzzp2_prime6, ghzzp2_prime7 /)
         lambda_v = Lambda_z2
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.14) then
         vvcoupl = (/ ghzzp3, ghzzp3_prime, ghzzp3_prime2, ghzzp3_prime3, ghzzp3_prime4, ghzzp3_prime5, ghzzp3_prime6, ghzzp3_prime7 /)
         lambda_v = Lambda_z3
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.15) then
         vvcoupl = (/ ghzzp4, ghzzp4_prime, ghzzp4_prime2, ghzzp4_prime3, ghzzp4_prime4, ghzzp4_prime5, ghzzp4_prime6, ghzzp4_prime7 /)
         lambda_v = Lambda_z4
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.16) then ! ZpZp 1-4
         vvcoupl = (/ ghzpzp1, ghzpzp1_prime, ghzpzp1_prime2, ghzpzp1_prime3, ghzpzp1_prime4, ghzpzp1_prime5, ghzpzp1_prime6, ghzpzp1_prime7 /)
         lambda_v = Lambda_z1
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.17) then
         vvcoupl = (/ ghzpzp2, ghzpzp2_prime, ghzpzp2_prime2, ghzpzp2_prime3, ghzpzp2_prime4, ghzpzp2_prime5, ghzpzp2_prime6, ghzpzp2_prime7 /)
         lambda_v = Lambda_z2
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.18) then
         vvcoupl = (/ ghzpzp3, ghzpzp3_prime, ghzpzp3_prime2, ghzpzp3_prime3, ghzpzp3_prime4, ghzpzp3_prime5, ghzpzp3_prime6, ghzpzp3_prime7 /)
         lambda_v = Lambda_z3
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.19) then
         vvcoupl = (/ ghzpzp4, ghzpzp4_prime, ghzpzp4_prime2, ghzpzp4_prime3, ghzpzp4_prime4, ghzpzp4_prime5, ghzpzp4_prime6, ghzpzp4_prime7 /)
         lambda_v = Lambda_z4
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.20) then ! Zpgs 1
         vvcoupl(3) = ghzpgs1_prime2
         lambda_v = Lambda_zgs1
         lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
      elseif(index.eq.21) then ! Zpgs 2-4
         vvcoupl(1) = ghzpgs2
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
      elseif(index.eq.22) then
         vvcoupl(1) = ghzpgs3
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
      elseif(index.eq.23) then
         vvcoupl(1) = ghzpgs4
         lambda_v = 1d0 ! Not present
         lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
      endif
   else
      if(cw_q1sq.ne.0) sWplus_signed=abs(sWplus)*dble(sign(1,cw_q1sq))
      if(cw_q2sq.ne.0) sWminus_signed=abs(sWminus)*dble(sign(1,cw_q2sq))
      if(cw_q12sq.ne.0) sWW_signed=abs(sWW)*dble(sign(1,cw_q12sq))
      if(cw_q1sq.ne.0 .or. cw_q2sq.ne.0 .or. cw_q12sq.ne.0) computeQsqCompundCoupl=.true.
      if(index.eq.1) then
         vvcoupl = (/ ghw1, ghw1_prime, ghw1_prime2, ghw1_prime3, ghw1_prime4, ghw1_prime5, ghw1_prime6, ghw1_prime7 /)
         lambda_v = Lambda_w1
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.2) then
         vvcoupl = (/ ghw2, ghw2_prime, ghw2_prime2, ghw2_prime3, ghw2_prime4, ghw2_prime5, ghw2_prime6, ghw2_prime7 /)
         lambda_v = Lambda_w2
         lambda_v120 = (/ Lambda_w21, Lambda_w22, Lambda_w20 /)
      elseif(index.eq.3) then
         vvcoupl = (/ ghw3, ghw3_prime, ghw3_prime2, ghw3_prime3, ghw3_prime4, ghw3_prime5, ghw3_prime6, ghw3_prime7 /)
         lambda_v = Lambda_w3
         lambda_v120 = (/ Lambda_w31, Lambda_w32, Lambda_w30 /)
      elseif(index.eq.4) then
         vvcoupl = (/ ghw4, ghw4_prime, ghw4_prime2, ghw4_prime3, ghw4_prime4, ghw4_prime5, ghw4_prime6, ghw4_prime7 /)
         lambda_v = Lambda_w4
         lambda_v120 = (/ Lambda_w41, Lambda_w42, Lambda_w40 /)
      elseif(index.eq.12) then
         vvcoupl = (/ ghwwp1, ghwwp1_prime, ghwwp1_prime2, ghwwp1_prime3, ghwwp1_prime4, ghwwp1_prime5, ghwwp1_prime6, ghwwp1_prime7 /)
         lambda_v = Lambda_w1
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.13) then
         vvcoupl = (/ ghwwp2, ghwwp2_prime, ghwwp2_prime2, ghwwp2_prime3, ghwwp2_prime4, ghwwp2_prime5, ghwwp2_prime6, ghwwp2_prime7 /)
         lambda_v = Lambda_w2
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.14) then
         vvcoupl = (/ ghwwp3, ghwwp3_prime, ghwwp3_prime2, ghwwp3_prime3, ghwwp3_prime4, ghwwp3_prime5, ghwwp3_prime6, ghwwp3_prime7 /)
         lambda_v = Lambda_w3
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.15) then
         vvcoupl = (/ ghwwp4, ghwwp4_prime, ghwwp4_prime2, ghwwp4_prime3, ghwwp4_prime4, ghwwp4_prime5, ghwwp4_prime6, ghwwp4_prime7 /)
         lambda_v = Lambda_w4
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.16) then
         vvcoupl = (/ ghwpwp1, ghwpwp1_prime, ghwpwp1_prime2, ghwpwp1_prime3, ghwpwp1_prime4, ghwpwp1_prime5, ghwpwp1_prime6, ghwpwp1_prime7 /)
         lambda_v = Lambda_w1
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.17) then
         vvcoupl = (/ ghwpwp2, ghwpwp2_prime, ghwpwp2_prime2, ghwpwp2_prime3, ghwpwp2_prime4, ghwpwp2_prime5, ghwpwp2_prime6, ghwpwp2_prime7 /)
         lambda_v = Lambda_w2
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.18) then
         vvcoupl = (/ ghwpwp3, ghwpwp3_prime, ghwpwp3_prime2, ghwpwp3_prime3, ghwpwp3_prime4, ghwpwp3_prime5, ghwpwp3_prime6, ghwpwp3_prime7 /)
         lambda_v = Lambda_w3
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      elseif(index.eq.19) then
         vvcoupl = (/ ghwpwp4, ghwpwp4_prime, ghwpwp4_prime2, ghwpwp4_prime3, ghwpwp4_prime4, ghwpwp4_prime5, ghwpwp4_prime6, ghwpwp4_prime7 /)
         lambda_v = Lambda_w4
         lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
      endif
   endif

   if(vvcoupl(2).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(2) * lambda_v**4/(lambda_v**2 + abs(sWplus))/(lambda_v**2 + abs(sWminus))
   if(vvcoupl(3).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(3) * ( sWplus + sWminus )/lambda_v**2
   if(vvcoupl(4).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(4) * ( sWplus - sWminus )/lambda_v**2
   if(vvcoupl(5).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(5) * ( sWW )/Lambda_Q**2
   if(vvcoupl(6).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(6) * ( sWplus**2 + sWminus**2 )/lambda_v**4
   if(vvcoupl(7).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(7) * ( sWplus**2 - sWminus**2 )/lambda_v**4
   if(vvcoupl(8).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(8) * ( sWplus    * sWminus    )/lambda_v**4

   if(index.eq.1) then
      if(computeQsqCompundCoupl) then
         QsqCompoundFactor = (lambda_v120(1)**2 + sWplus_signed)*(lambda_v120(2)**2 + sWminus_signed)*(lambda_v120(3)**2 + sWW_signed)
         if(QsqCompoundFactor.ne.0d0) then
            QsqCompoundFactor = (lambda_v120(1)*lambda_v120(2)*lambda_v120(3))**2/QsqCompoundFactor
         endif
         HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling * QsqCompoundFactor
      endif
      if(vvcoupl(1).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(1)
   else
      if(vvcoupl(1).ne.czero) HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling + vvcoupl(1)
      if(computeQsqCompundCoupl) then
         QsqCompoundFactor = (lambda_v120(1)**2 + sWplus_signed)*(lambda_v120(2)**2 + sWminus_signed)*(lambda_v120(3)**2 + sWW_signed)
         if(QsqCompoundFactor.ne.0d0) then
            QsqCompoundFactor = (lambda_v120(1)*lambda_v120(2)*lambda_v120(3))**2/QsqCompoundFactor
         endif
         HVVSpinZeroDynamicCoupling = HVVSpinZeroDynamicCoupling * QsqCompoundFactor
      endif
   endif

return
end function


function VpffCoupling(jhuid, hel, useWp)
integer, intent(in) :: jhuid
integer, intent(in) :: hel
logical, intent(in) :: useWp
complex(8) :: VpffCoupling
integer :: absid

   VpffCoupling=czero
   if (abs(hel).ne.1) then
      return
   endif
   absid=abs(jhuid)

   if(useWp) then
      if (absid.eq.abs(Up_) .or. absid.eq.abs(Dn_)) then
         if(hel.eq.-1) then
            VpffCoupling=ewp_Up_left
         else
            VpffCoupling=ewp_Up_right
         endif
      elseif (absid.eq.abs(Chm_) .or. absid.eq.abs(Str_)) then
         if(hel.eq.-1) then
            VpffCoupling=ewp_Chm_left
         else
            VpffCoupling=ewp_Chm_right
         endif
      elseif (absid.eq.abs(Top_) .or. absid.eq.abs(Bot_)) then
         if(hel.eq.-1) then
            VpffCoupling=ewp_Top_left
         else
            VpffCoupling=ewp_Top_right
         endif
      elseif (absid.eq.abs(ElP_) .or. absid.eq.abs(NuE_)) then
         if(hel.eq.-1) then
            VpffCoupling=ewp_El_left
         else
            VpffCoupling=ewp_El_right
         endif
      elseif (absid.eq.abs(MuP_) .or. absid.eq.abs(NuM_)) then
         if(hel.eq.-1) then
            VpffCoupling=ewp_Mu_left
         else
            VpffCoupling=ewp_Mu_right
         endif
      elseif (absid.eq.abs(TaP_) .or. absid.eq.abs(NuT_)) then
         if(hel.eq.-1) then
            VpffCoupling=ewp_Ta_left
         else
            VpffCoupling=ewp_Ta_right
         endif
      endif
   else
      if (absid.eq.abs(Up_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Up_left
         else
            VpffCoupling=ezp_Up_right
         endif
      elseif (absid.eq.abs(Dn_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Dn_left
         else
            VpffCoupling=ezp_Dn_right
         endif
      elseif (absid.eq.abs(Chm_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Chm_left
         else
            VpffCoupling=ezp_Chm_right
         endif
      elseif (absid.eq.abs(Str_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Str_left
         else
            VpffCoupling=ezp_Str_right
         endif
      elseif (absid.eq.abs(Top_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Top_left
         else
            VpffCoupling=ezp_Top_right
         endif
      elseif (absid.eq.abs(Bot_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Bot_left
         else
            VpffCoupling=ezp_Bot_right
         endif
      elseif (absid.eq.abs(ElP_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_El_left
         else
            VpffCoupling=ezp_El_right
         endif
      elseif (absid.eq.abs(NuE_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_NuE_left
         else
            VpffCoupling=ezp_NuE_right
         endif
      elseif (absid.eq.abs(MuP_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Mu_left
         else
            VpffCoupling=ezp_Mu_right
         endif
      elseif (absid.eq.abs(NuM_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_NuE_left
         else
            VpffCoupling=ezp_NuE_right
         endif
      elseif (absid.eq.abs(TaP_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_Ta_left
         else
            VpffCoupling=ezp_Ta_right
         endif
      elseif (absid.eq.abs(NuT_)) then
         if(hel.eq.-1) then
            VpffCoupling=ezp_NuE_left
         else
            VpffCoupling=ezp_NuE_right
         endif
      endif
   endif

end function
function VpffCoupling_PDG(pdgid, hel, useWp)
integer, intent(in) :: pdgid
integer, intent(in) :: hel
logical, intent(in) :: useWp
complex(8) :: VpffCoupling_PDG
integer :: jhuid
   jhuid=convertLHEreverse(pdgid)
   VpffCoupling_PDG=VpffCoupling(jhuid,hel,useWp)
end function


FUNCTION ScaleFactor(id1in,id2in)
implicit none
real(8) :: ScaleFactor
integer :: id1, id2, id1in, id2in
id1 = abs(id1in)
id2 = abs(id2in)

! W->qq
if((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Up_)))then
  ScaleFactor = scale_alpha_W_ud
elseif((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Up_)))then
  ScaleFactor = scale_alpha_W_ud
elseif((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Up_)))then
  ScaleFactor = scale_alpha_W_ud
elseif((id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Chm_)))then
  ScaleFactor = scale_alpha_W_cs
elseif((id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Chm_)))then
  ScaleFactor = scale_alpha_W_cs
elseif((id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Chm_)))then
  ScaleFactor = scale_alpha_W_cs
! W-> td
elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Top_)))then
  ScaleFactor = 1d0
elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Top_)))then
  ScaleFactor = 1d0
elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Top_)))then
  ScaleFactor = 1d0
! W->lnu
elseif((abs(id1).eq.abs(convertLHE(NuT_))  .and.  abs(id2).eq.abs(convertLHE(TaP_)))  .or.  (abs(id1).eq.abs(convertLHE(TaP_))  .and.  abs(id2).eq.abs(convertLHE(NuT_))))then
  ScaleFactor = scale_alpha_W_tn
elseif((abs(id1).eq.abs(convertLHE(NuM_))  .and.  abs(id2).eq.abs(convertLHE(MuP_)))  .or.  (abs(id1).eq.abs(convertLHE(MuP_))  .and.  abs(id2).eq.abs(convertLHE(NuM_))))then
  ScaleFactor = scale_alpha_W_ln
elseif((abs(id1).eq.abs(convertLHE(NuE_))  .and.  abs(id2).eq.abs(convertLHE(ElP_)))  .or.  (abs(id1).eq.abs(convertLHE(ElP_))  .and.  abs(id2).eq.abs(convertLHE(NuE_))))then
  ScaleFactor = scale_alpha_W_ln
! Z->qq
elseif((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Up_))  .or.  (id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Chm_)))then
  ScaleFactor = scale_alpha_Z_uu
elseif((id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Bot_)))then
  ScaleFactor = scale_alpha_Z_dd
! Z-> ll, nunu
elseif(id1.eq.convertLHE(TaP_)  .and.  id2.eq.convertLHE(TaP_))then
  ScaleFactor = scale_alpha_Z_tt
elseif((id1.eq.convertLHE(MuP_)  .and.  id2.eq.convertLHE(MuP_))  .or.  (id1.eq.convertLHE(ElP_)  .and.  id2.eq.convertLHE(ElP_)))then
  ScaleFactor = scale_alpha_Z_ll
elseif((id1.eq.convertLHE(NuT_)  .and.  id2.eq.convertLHE(NuT_))  .or.  (id1.eq.convertLHE(NuM_)  .and.  id2.eq.convertLHE(NuM_))  .or.  (id1.eq.convertLHE(NuE_)  .and.  id2.eq.convertLHE(NuE_)))then
  ScaleFactor = scale_alpha_Z_nn
! Everything else
else
  ScaleFactor = 1d0
endif

return
END FUNCTION

FUNCTION CKMbare(id1in,id2in)
implicit none
real(8) :: CKMbare
integer :: id1, id2, id1in, id2in
id1 = abs(id1in)
id2 = abs(id2in)
if((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Up_)))then
  CKMbare= VCKM_ud
elseif((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Up_)))then
  CKMbare= VCKM_us
elseif((id1.eq.convertLHE(Up_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Up_)))then
  CKMbare= VCKM_ub
elseif((id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Chm_)))then
  CKMbare= VCKM_cd
elseif((id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Chm_)))then
  CKMbare= VCKM_cs
elseif((id1.eq.convertLHE(Chm_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Chm_)))then
  CKMbare= VCKM_cb
elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Dn_))  .or.  (id1.eq.convertLHE(Dn_)  .and.  id2.eq.convertLHE(Top_)))then
  CKMbare= VCKM_td
elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Str_))  .or.  (id1.eq.convertLHE(Str_)  .and.  id2.eq.convertLHE(Top_)))then
  CKMbare= VCKM_ts
elseif((id1.eq.convertLHE(Top_)  .and.  id2.eq.convertLHE(Bot_))  .or.  (id1.eq.convertLHE(Bot_)  .and.  id2.eq.convertLHE(Top_)))then
  CKMbare= VCKM_tb
else
  CKMbare= 1d0
endif

END FUNCTION

FUNCTION CKM(id1in,id2in)
implicit none
real(8) :: CKM
integer :: id1in, id2in
  CKM=CKMbare(id1in, id2in)*sqrt(ScaleFactor(id1in,id2in))
END FUNCTION

! for VH use
FUNCTION ZFFbare(id1in,id2in,h1,h2)
implicit none
real(8) :: ZFFbare
integer :: id1, id2, id1in, id2in
real(8) :: h1,h2

id1 = abs(id1in)
id2 = abs(id2in)

  if(id1.ne.id2)then
    print*,"Not a valid Z to fermion pair vertex!",id1in,id2in
    stop
  !Zll
  else if(id1.eq.convertLHE(ElM_).or.id1.eq.convertLHE(MuM_).or.id1.eq.convertLHE(TaM_))then
    if((id1in*h1).gt.0d0)then
      ZFFbare = aR_lep
    else
      ZFFbare = aL_lep
    endif
  !Zuu
  else if(id1.eq.convertLHE(Up_).or.id1.eq.convertLHE(Chm_).or.id1.eq.convertLHE(Top_))then
    if((id1in*h1).gt.0d0)then
      ZFFbare = aR_QUp
    else
      ZFFbare = aL_QUp
    endif
  !Zdd
  else if(id1.eq.convertLHE(Dn_).or.id1.eq.convertLHE(Str_).or.id1.eq.convertLHE(Bot_))then
    if((id1in*h1).gt.0d0)then
      ZFFbare = aR_QDn
    else
      ZFFbare = aL_QDn
    endif
  !Z nu nu
  else if(id1.eq.convertLHE(NuE_).or.id1.eq.convertLHE(NuM_).or.id1.eq.convertLHE(NuT_))then
    if((id1in*h1).lt.0d0)then
      ZFFbare = aL_neu
    else
      print*,"No right-handed neutrino here!",id1in,h1
      stop
    endif
  else
    print*,"Not a valid Z to fermion pair vertex!",id1in,id2in
    stop
  endif

return
END FUNCTION

FUNCTION ZFF(id1in,id2in,h1,h2)
implicit none
real(8) :: ZFF
integer :: id1in, id2in
real(8) :: h1,h2
  ZFF = ZFFbare(id1in,id2in,h1,h2) * sqrt(ScaleFactor(id1in,id2in))
return
END FUNCTION

FUNCTION AFF(id1in,id2in,h1,h2)
implicit none
real(8) :: AFF
integer :: id1, id2, id1in, id2in
real(8) :: h1,h2
logical, save :: printedneutrinowarning = .false.

id1 = abs(id1in)
id2 = abs(id2in)

  if(id1.ne.id2)then
    print*,"Not a valid Z to fermion pair vertex!",id1in,id2in
    stop
  !photon ll
  else if(id1.eq.convertLHE(ElM_).or.id1.eq.convertLHE(MuM_).or.id1.eq.convertLHE(TaM_))then
    if((id1in*h1).gt.0d0)then
      AFF = QlR
    else
      AFF = QlL
    endif
  !photon uu
  else if(id1.eq.convertLHE(Up_).or.id1.eq.convertLHE(Chm_).or.id1.eq.convertLHE(Top_))then
    if((id1in*h1).gt.0d0)then
      AFF = QuR
    else
      AFF = QuL
    endif
  !photon dd
  else if(id1.eq.convertLHE(Dn_).or.id1.eq.convertLHE(Str_).or.id1.eq.convertLHE(Bot_))then
    if((id1in*h1).gt.0d0)then
      AFF = QdR
    else
      AFF = QdL
    endif
  !photon nu nu
  else if(id1.eq.convertLHE(NuE_).or.id1.eq.convertLHE(NuM_).or.id1.eq.convertLHE(NuT_))then
    AFF = 0d0
    if (.not. printedneutrinowarning) then
      print*,"Warning, gamma > nu nu~ gives 0.",id1in,id2in
      printedneutrinowarning = .true.
    endif
  else
    print*,"Not a valid photon to fermion pair vertex!",id1in,id2in
    stop
  endif

return
END FUNCTION

! for VH use


FUNCTION convertLHEreverse(Part) ! PDG/PDF->JHU
implicit none
integer :: convertLHEreverse
integer :: Part

  if(     Part.eq.0 ) then      ! 0=Glu_ is not the official LHE convention
      convertLHEreverse = Glu_
  elseif( Part.eq.1 ) then
      convertLHEreverse = Dn_
  elseif( Part.eq.2 ) then
      convertLHEreverse = Up_
  elseif( Part.eq.3 ) then
      convertLHEreverse = Str_
  elseif( Part.eq.4 ) then
      convertLHEreverse = Chm_
  elseif( Part.eq.5 ) then
      convertLHEreverse = Bot_
  elseif( Part.eq.6 ) then
      convertLHEreverse = Top_
  elseif( Part.eq.-1 ) then
      convertLHEreverse = ADn_
  elseif( Part.eq.-2 ) then
      convertLHEreverse = AUp_
  elseif( Part.eq.-3 ) then
      convertLHEreverse = AStr_
  elseif( Part.eq.-4 ) then
      convertLHEreverse = AChm_
  elseif( Part.eq.-5 ) then
      convertLHEreverse = ABot_
  elseif( Part.eq.-6 ) then
      convertLHEreverse = ATop_
  elseif( Part.eq.21 ) then
      convertLHEreverse = Glu_
  elseif( Part.eq.11 ) then
      convertLHEreverse = ElM_
  elseif( Part.eq.22 ) then
      convertLHEreverse = Pho_
  elseif( Part.eq.23 ) then
      convertLHEreverse = Z0_
  elseif( Part.eq.24 ) then
      convertLHEreverse = Wp_
  elseif( Part.eq.-24 ) then
      convertLHEreverse = Wm_
  elseif( Part.eq.-11 ) then
      convertLHEreverse = ElP_
  elseif( Part.eq.13 ) then
      convertLHEreverse = MuM_
  elseif( Part.eq.-13 ) then
      convertLHEreverse = MuP_
  elseif( Part.eq.15 ) then
      convertLHEreverse = TaM_
  elseif( Part.eq.-15 ) then
      convertLHEreverse = TaP_
  elseif( Part.eq.12 ) then
      convertLHEreverse = NuE_
  elseif( Part.eq.-12) then
      convertLHEreverse = ANuE_
  elseif( Part.eq.14) then
      convertLHEreverse = NuM_
  elseif( Part.eq.-14) then
      convertLHEreverse = ANuM_
  elseif( Part.eq.16 ) then
      convertLHEreverse = NuT_
  elseif( Part.eq.-16) then
      convertLHEreverse = ANuT_
  elseif( Part.eq.+25) then
      convertLHEreverse = Hig_
  elseif( Part.eq.-25) then
      convertLHEreverse = Hig_
  elseif( Part.eq.32) then
      convertLHEreverse = Zpr_
  elseif( Part.eq.33) then
      convertLHEreverse = Zpr2_
  elseif( Part.eq.34) then
      convertLHEreverse = Wppr_
  elseif( Part.eq.-34) then
      convertLHEreverse = Wmpr_
  else
      print *, "MYLHE format not implemented for ",Part
      stop
  endif


END FUNCTION



FUNCTION convertLHE(Part) ! JHU->PDG/PDF
implicit none
integer :: convertLHE
integer :: Part


  if(     Part.eq.0 ) then
      convertLHE = 0
  elseif( Part.eq.Glu_ ) then
      convertLHE = 21
  elseif( Part.eq.ElM_ ) then
      convertLHE = 11
  elseif( Part.eq.ElP_ ) then
      convertLHE =-11
  elseif( Part.eq.MuM_ ) then
      convertLHE = 13
  elseif( Part.eq.MuP_ ) then
      convertLHE =-13
  elseif( Part.eq.TaM_ ) then
      convertLHE = 15
  elseif( Part.eq.TaP_ ) then
      convertLHE =-15
  elseif( Part.eq.NuE_ ) then
      convertLHE = 12
  elseif( Part.eq.ANuE_) then
      convertLHE =-12
  elseif( Part.eq.NuM_) then
      convertLHE = 14
  elseif( Part.eq.ANuM_) then
      convertLHE =-14
  elseif( Part.eq.NuT_ ) then
      convertLHE = 16
  elseif( Part.eq.ANuT_) then
      convertLHE =-16
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
  elseif( Part.eq.Top_ ) then
      convertLHE = 6
  elseif( Part.eq.ATop_) then
      convertLHE =-6
  elseif( Part.eq.Z0_) then
      convertLHE =23
  elseif( Part.eq.Wp_) then
      convertLHE =24
  elseif( Part.eq.Wm_) then
      convertLHE =-24
  elseif( Part.eq.Pho_) then
      convertLHE =22
  elseif( Part.eq.Hig_) then
      convertLHE =25
  elseif( Part.eq.Zpr_) then
      convertLHE =32
  elseif( Part.eq.Zpr2_) then
      convertLHE =33
  elseif( Part.eq.Wppr_) then
      convertLHE =34
  elseif( Part.eq.Wmpr_) then
      convertLHE =-34
  elseif( Part.eq.Gra_) then
      convertLHE =39
  elseif( Part.le.Not_a_particle_) then
      convertLHE = Not_a_particle_
  else
      print *, "LHE format not implemented for ",Part
      stop
  endif

END FUNCTION


FUNCTION convertToPartIndex(Part) ! JHU->PDF
implicit none
integer :: convertToPartIndex
integer :: Part


  if( Part.eq.Glu_ ) then
      convertToPartIndex = pdfGlu_
  elseif( Part.eq.Up_  ) then
      convertToPartIndex = pdfUp_
  elseif( Part.eq.AUp_ ) then
      convertToPartIndex = pdfAUp_
  elseif( Part.eq.Dn_  ) then
      convertToPartIndex = pdfDn_
  elseif( Part.eq.ADn_ ) then
      convertToPartIndex = pdfADn_
  elseif( Part.eq.Chm_ ) then
      convertToPartIndex = pdfChm_
  elseif( Part.eq.AChm_) then
      convertToPartIndex = pdfAChm_
  elseif( Part.eq.Str_ ) then
      convertToPartIndex = pdfStr_
  elseif( Part.eq.AStr_) then
      convertToPartIndex = pdfAStr_
  elseif( Part.eq.Bot_ ) then
      convertToPartIndex = pdfBot_
  elseif( Part.eq.ABot_) then
      convertToPartIndex = pdfABot_
  else
      print *, "Unsuccessful conversion to a parton ME array index from ",Part
      stop
  endif

END FUNCTION


FUNCTION convertFromPartIndex(Part) ! PDF->JHU
implicit none
integer :: convertFromPartIndex
integer :: Part


  if( Part.eq.pdfGlu_ ) then
      convertFromPartIndex = Glu_
  elseif( Part.eq.pdfUp_  ) then
      convertFromPartIndex = Up_
  elseif( Part.eq.pdfAUp_ ) then
      convertFromPartIndex = AUp_
  elseif( Part.eq.pdfDn_  ) then
      convertFromPartIndex = Dn_
  elseif( Part.eq.pdfADn_ ) then
      convertFromPartIndex = ADn_
  elseif( Part.eq.pdfChm_ ) then
      convertFromPartIndex = Chm_
  elseif( Part.eq.pdfAChm_) then
      convertFromPartIndex = AChm_
  elseif( Part.eq.pdfStr_ ) then
      convertFromPartIndex = Str_
  elseif( Part.eq.pdfAStr_) then
      convertFromPartIndex = AStr_
  elseif( Part.eq.pdfBot_ ) then
      convertFromPartIndex = Bot_
  elseif( Part.eq.pdfABot_) then
      convertFromPartIndex = ABot_
  elseif( Part.eq.pdfTop_ ) then
      convertFromPartIndex = Top_
  elseif( Part.eq.pdfATop_) then
      convertFromPartIndex = ATop_
  else
      print *, "Unsuccessful conversion to a parton id from the ME array index ",Part
      stop
  endif

END FUNCTION


subroutine SetMass(mass, ipart)
implicit none
real(8), intent(in) :: mass
integer, intent(in) :: ipart
integer :: Part

  Part=abs(ipart)
  if( Part.eq.abs(ElM_) ) then
      m_el = mass
  elseif( Part.eq.abs(MuM_) ) then
      m_mu = mass
  elseif( Part.eq.abs(TaM_) ) then
      m_tau = mass
  elseif( Part.eq.abs(Chm_) ) then
      m_charm = mass
  elseif( Part.eq.abs(Bot_) ) then
      m_bot = mass
  elseif( Part.eq.abs(Top_) ) then
      m_top = mass
  elseif( Part.eq.abs(Z0_) ) then
      M_Z = mass
  elseif( Part.eq.abs(Zpr_) ) then
      M_Zprime = mass
  elseif( Part.eq.abs(Wp_) ) then
      M_W = mass
  elseif( Part.eq.abs(Wppr_) ) then
      M_Wprime = mass
  elseif( Part.eq.abs(Hig_) ) then
      M_Reso = mass
  endif

end subroutine SetMass


subroutine SetDecayWidth(width, ipart)
implicit none
real(8), intent(in) :: width
integer, intent(in) :: ipart
integer :: Part

  Part=abs(ipart)
  if( Part.eq.abs(TaM_) ) then
      Ga_tau = width
  elseif( Part.eq.abs(Top_) ) then
      Ga_Top = width
  elseif( Part.eq.abs(Z0_) ) then
      Ga_Z = width
  elseif( Part.eq.abs(Zpr_) ) then
      Ga_Zprime = width
  elseif( Part.eq.abs(Wp_) ) then
      Ga_W = width
  elseif( Part.eq.abs(Wppr_) ) then
      Ga_Wprime = width
  elseif( Part.eq.abs(Hig_) ) then
      Ga_Reso = width
  endif

END subroutine SetDecayWidth

FUNCTION getMass(Part)
implicit none
real(8) :: getMass
integer :: Part


  if( Part.eq.Glu_ ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(ElM_) ) then
      getMass = m_el
  elseif( abs(Part).eq.abs(MuM_) ) then
      getMass = m_mu
  elseif( abs(Part).eq.abs(TaM_) ) then
      getMass = m_tau
  elseif( abs(Part).eq.abs(NuE_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(NuM_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(NuT_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(Up_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(Dn_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(Chm_) ) then
      getMass = m_charm
  elseif( abs(Part).eq.abs(Str_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(Bot_) ) then
      getMass = m_bot
  elseif( abs(Part).eq.abs(Top_) ) then
      getMass = m_top
  elseif( abs(Part).eq.abs(Z0_) ) then
      getMass = M_Z
  elseif( abs(Part).eq.abs(Zpr_) ) then
      getMass = M_Zprime
  elseif( abs(Part).eq.abs(Wp_) ) then
      getMass = M_W
  elseif( abs(Part).eq.abs(Wppr_) ) then
      getMass = M_Wprime
  elseif( abs(Part).eq.abs(Pho_) ) then
      getMass = 0d0
  elseif( abs(Part).eq.abs(Hig_) ) then
      getMass = M_Reso
  elseif( Part.eq.Not_a_particle_) then
      getMass = 0d0
  else
     print *, "Error in getMass",Part
     stop
  endif

END FUNCTION

FUNCTION getDecayWidth(Part)
implicit none
real(8) :: getDecayWidth
integer :: Part

  getDecayWidth = 0d0
  if( abs(Part).eq.abs(Top_) ) then
      getDecayWidth = Ga_Top
  elseif( abs(Part).eq.abs(Z0_) ) then
      getDecayWidth = Ga_Z
  elseif( abs(Part).eq.abs(Zpr_) ) then
      getDecayWidth = Ga_Zprime
  elseif( abs(Part).eq.abs(Wp_) ) then
      getDecayWidth = Ga_W
  elseif( abs(Part).eq.abs(Wppr_) ) then
      getDecayWidth = Ga_Wprime
  elseif( abs(Part).eq.abs(Hig_) ) then
      getDecayWidth = Ga_Reso
  elseif( abs(Part).eq.abs(TaM_) ) then
      getDecayWidth = Ga_tau
  endif

END FUNCTION


subroutine SetHiggsMass(jH,mass)
implicit none
real(8) :: mass
integer :: jH
   if (jH .eq. 1) then
      M_Reso=mass
   elseif (jH .eq.2) then
      M_Reso2=mass
   endif
END subroutine

subroutine SetHiggsDecayWidth(jH,width)
implicit none
real(8) :: width
integer :: jH
   if (jH .eq. 1) then
      Ga_Reso=width
   elseif (jH .eq.2) then
      Ga_Reso2=width
   endif
END subroutine

FUNCTION GetHiggsMass(jH)
implicit none
real(8) :: getHiggsMass
integer :: jH
   getHiggsMass=0d0
   if (jH .eq. 1) then
      getHiggsMass = M_Reso
   elseif (jH .eq.2) then
      getHiggsMass = M_Reso2
   endif
END FUNCTION

FUNCTION GetHiggsDecayWidth(jH)
implicit none
real(8) :: getHiggsDecayWidth
integer :: jH
   getHiggsDecayWidth=0d0
   if (jH .eq. 1) then
      getHiggsDecayWidth = Ga_Reso
   elseif (jH .eq.2) then
      getHiggsDecayWidth = Ga_Reso2
   endif
END FUNCTION


FUNCTION getParticle(Part)
implicit none
character(len=3) :: getParticle
integer :: Part


  if( Part.eq.Glu_ ) then
      getParticle = "glu"
  elseif( Part.eq.0 ) then
      getParticle = "glu"
  elseif( Part.eq.ElM_ ) then
      getParticle = "el-"
  elseif( Part.eq.ElP_ ) then
      getParticle = "el+"
  elseif( Part.eq.MuM_ ) then
      getParticle = "mu-"
  elseif( Part.eq.MuP_ ) then
      getParticle = "mu+"
  elseif( Part.eq.TaM_ ) then
      getParticle = "ta-"
  elseif( Part.eq.TaP_ ) then
      getParticle = "t+-"
  elseif( Part.eq.NuE_ ) then
      getParticle = "nuE"
  elseif( Part.eq.ANuE_ ) then
      getParticle = "AnE"
  elseif( Part.eq.NuM_ ) then
      getParticle = "nuM"
  elseif( Part.eq.ANuM_ ) then
      getParticle = "AnM"
  elseif( Part.eq.NuT_ ) then
      getParticle = "nuT"
  elseif( Part.eq.ANuT_ ) then
      getParticle = "AnT"
  elseif( Part.eq.Up_ ) then
      getParticle = " up"
  elseif( Part.eq.AUp_ ) then
      getParticle = "Aup"
  elseif( Part.eq.Dn_ ) then
      getParticle = " dn"
  elseif( Part.eq.ADn_ ) then
      getParticle = "Adn"
  elseif( Part.eq.Chm_ ) then
      getParticle = "chm"
  elseif( Part.eq.AChm_ ) then
      getParticle = "Achm"
  elseif( Part.eq.Str_ ) then
      getParticle = "str"
  elseif( Part.eq.AStr_ ) then
      getParticle = "Astr"
  elseif( Part.eq.Bot_ ) then
      getParticle = "bot"
  elseif( Part.eq.ABot_ ) then
      getParticle = "Abot"
  elseif( Part.eq.Top_ ) then
      getParticle = "top"
  elseif( Part.eq.ATop_ ) then
      getParticle = "Atop"
  elseif( Part.eq.Z0_ ) then
      getParticle = " Z0"
  elseif( Part.eq.Zpr_ ) then
      getParticle = " Zprime"
  elseif( Part.eq.Wp_ ) then
      getParticle = " W+"
  elseif( Part.eq.Wm_ ) then
      getParticle = " W-"
  elseif( Part.eq.Wppr_ ) then
      getParticle = " W+prime"
  elseif( Part.eq.Wmpr_ ) then
      getParticle = " W-prime"
  elseif( Part.eq.Pho_ ) then
      getParticle = "pho"
  elseif( Part.eq.Hig_ ) then
      getParticle = "Hig"
  else
     print *, "Error in getParticle",Part
     stop
  endif


END FUNCTION


FUNCTION getLHEParticle(PartLHE)
implicit none
character(len=3) :: getLHEParticle
integer :: PartLHE,Part


  Part = convertLHEreverse(PartLHE)
  getLHEParticle = getParticle(Part)

END FUNCTION




FUNCTION IsAZDecay(DKMode)
implicit none
logical :: IsAZDecay
integer :: DKMode


  if( DKMode.eq.0 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.1 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.2 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.3 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.8 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.9 ) then
     IsAZDecay = .true.
  elseif( DKMode.lt.0 ) then
     IsAZDecay = (            &
        DKMode.eq.-2*2   .or. & ! Z->ee
        DKMode.eq.-3*3   .or. & ! Z->mumu
        DKMode.eq.-5*5   .or. & ! Z->dd
        DKMode.eq.-7*7   .or. & ! Z->uu
        DKMode.eq.-11*11 .or. & ! Z->ss
        DKMode.eq.-13*13 .or. & ! Z->cc
        DKMode.eq.-17*17      & ! Z->bb
     )
  else
     IsAZDecay=.false.
  endif

END FUNCTION


FUNCTION IsAWDecay(DKMode)
implicit none
logical :: IsAWDecay
integer :: DKMode


  if( DKMode.eq.4 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.5 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.6 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.10 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.11 ) then
     IsAWDecay = .true.
  elseif( DKMode.lt.0 ) then
     IsAWDecay = (            &
        DKMode.eq.-2*1   .or. & ! W->enu
        DKMode.eq.-3*1   .or. & ! W->munu
        DKMode.eq.-5*7   .or. & ! W->du
        DKMode.eq.-5*13  .or. & ! W->dc
        DKMode.eq.-11*7  .or. & ! W->su
        DKMode.eq.-11*13 .or. & ! W->sc
        DKMode.eq.-17*7  .or. & ! W->bu
        DKMode.eq.-17*13      & ! W->bc
     )
  else
     IsAWDecay=.false.
  endif


END FUNCTION


FUNCTION IsAPhoton(DKMode)
implicit none
logical :: IsAPhoton
integer :: DKMode


  if( DKMode.eq.7 ) then
     IsAPhoton = .true.
  else
     IsAPhoton=.false.
  endif


END FUNCTION



FUNCTION IsDownTypeQuark(PartType)
implicit none
logical :: IsDownTypeQuark
integer :: PartType
   IsDownTypeQuark = ( abs(PartType).eq.abs(Dn_) .or. abs(PartType).eq.abs(Str_) .or. abs(PartType).eq.abs(Bot_) )
END FUNCTION
FUNCTION IsLHEDownTypeQuark(PartType)
implicit none
logical :: IsLHEDownTypeQuark
integer :: PartType
   IsLHEDownTypeQuark = ( abs(PartType).eq.1 .or. abs(PartType).eq.3 .or. abs(PartType).eq.5 )
END FUNCTION

FUNCTION IsUpTypeQuark(PartType)
implicit none
logical :: IsUpTypeQuark
integer :: PartType
   IsUpTypeQuark = ( abs(PartType).eq.abs(Up_) .or. abs(PartType).eq.abs(Chm_) .or. abs(PartType).eq.abs(Top_) )
END FUNCTION
FUNCTION IsLHEUpTypeQuark(PartType)
implicit none
logical :: IsLHEUpTypeQuark
integer :: PartType
   IsLHEUpTypeQuark = ( abs(PartType).eq.2 .or. abs(PartType).eq.4 .or. abs(PartType).eq.6 )
END FUNCTION

FUNCTION IsUpTypeLightQuark(PartType)
implicit none
logical :: IsUpTypeLightQuark
integer :: PartType
   IsUpTypeLightQuark = ( abs(PartType).eq.abs(Up_) .or. abs(PartType).eq.abs(Chm_))
END FUNCTION
FUNCTION IsLHEUpTypeLightQuark(PartType)
implicit none
logical :: IsLHEUpTypeLightQuark
integer :: PartType
   IsLHEUpTypeLightQuark = ( abs(PartType).eq.2 .or. abs(PartType).eq.4 )
END FUNCTION

FUNCTION IsAQuark(PartType)
implicit none
logical :: IsAQuark
integer :: PartType
   IsAQuark=IsUpTypeQuark(PartType) .or. IsDownTypeQuark(PartType)
END FUNCTION
FUNCTION IsALHEQuark(PartType)
implicit none
logical :: IsALHEQuark
integer :: PartType
   IsALHEQuark=IsLHEUpTypeQuark(PartType) .or. IsLHEDownTypeQuark(PartType)
END FUNCTION

FUNCTION IsALightQuark(PartType)
implicit none
logical :: IsALightQuark
integer :: PartType
   IsALightQuark=IsUpTypeLightQuark(PartType) .or. IsDownTypeQuark(PartType)
END FUNCTION
FUNCTION IsALHELightQuark(PartType)
implicit none
logical :: IsALHELightQuark
integer :: PartType
   IsALHELightQuark=IsLHEUpTypeLightQuark(PartType) .or. IsLHEDownTypeQuark(PartType)
END FUNCTION


FUNCTION IsANeutrino(PartType)
implicit none
logical :: IsANeutrino
integer :: PartType
   IsANeutrino = ( abs(PartType).eq.abs(NuE_) .or. abs(PartType).eq.abs(NuM_) .or. abs(PartType).eq.abs(NuT_) )
END FUNCTION
FUNCTION IsALHENeutrino(PartType)
implicit none
logical :: IsALHENeutrino
integer :: PartType
   IsALHENeutrino = ( abs(PartType).eq.12 .or. abs(PartType).eq.14 .or. abs(PartType).eq.16 )
END FUNCTION

FUNCTION IsALepton(PartType)! note that lepton means charged lepton here
implicit none
logical :: IsALepton
integer :: PartType
  IsALepton = ( abs(PartType).eq.abs(ElP_) .or. abs(PartType).eq.abs(MuP_) .or. abs(PartType).eq.abs(TaP_) )
END FUNCTION
FUNCTION IsALHELepton(PartType)! note that lepton means charged lepton here
implicit none
logical :: IsALHELepton
integer :: PartType
  IsALHELepton = ( abs(PartType).eq.11 .or. abs(PartType).eq.13 .or. abs(PartType).eq.15 )
END FUNCTION


FUNCTION IsAGluon(PartType)
implicit none
logical :: IsAGluon
integer :: PartType
  IsAGluon = (abs(PartType).eq.Glu_)
END FUNCTION

FUNCTION IsALHEGluon(PartType)
implicit none
logical :: IsALHEGluon
integer :: PartType
  IsALHEGluon = (abs(PartType).eq.convertLHE(Glu_))
END FUNCTION


FUNCTION IsAJet(PartType)
implicit none
logical :: IsAJet
integer :: PartType
  IsAJet = (IsALightQuark(PartType) .or. IsAGluon(PartType))
END FUNCTION

FUNCTION IsALHEJet(PartType)
implicit none
logical :: IsALHEJet
integer :: PartType
  IsALHEJet = (IsALHELightQuark(PartType) .or. IsALHEGluon(PartType))
END FUNCTION


FUNCTION IsABoson(PartType)
implicit none
logical :: IsABoson
integer :: PartType
  IsABoson = ( abs(PartType).eq.abs(Pho_) .or. abs(PartType).eq.abs(Z0_) .or. abs(PartType).eq.abs(Wp_) .or. abs(PartType).eq.abs(Hig_) )
END FUNCTION

FUNCTION IsALHEBoson(PartType)
implicit none
logical :: IsALHEBoson
integer :: PartType
  IsALHEBoson = ( abs(PartType).ge.22 .and. abs(PartType).le.25 )
END FUNCTION

function CoupledVertex(id,hel,useAHcoupl)
   implicit none
   integer, optional :: useAHcoupl
   integer, intent(in) :: id(1:2),hel
   integer :: testAHcoupl
   integer :: CoupledVertex

   testAHcoupl = 0
   if(present(useAHcoupl)) then
      testAHcoupl = useAHcoupl
   endif
   if( (&
   (id(1).eq.ElP_ .and. id(2).eq.NuE_) .or. (id(2).eq.ElP_ .and. id(1).eq.NuE_) .or. &
   (id(1).eq.MuP_ .and. id(2).eq.NuM_) .or. (id(2).eq.MuP_ .and. id(1).eq.NuM_) .or. &
   (id(1).eq.TaP_ .and. id(2).eq.NuT_) .or. (id(2).eq.TaP_ .and. id(1).eq.NuT_) .or. &
   (id(1).eq.Up_  .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Up_  .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_)) .or. &
   (id(1).eq.Chm_ .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Chm_ .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_)) .or. &
   (id(1).eq.Top_ .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Top_ .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_))      &
   ) .and. hel.lt.0) then
      CoupledVertex=Wp_
   elseif( (&
   (id(1).eq.ElM_ .and. id(2).eq.ANuE_) .or. (id(2).eq.ElM_ .and. id(1).eq.ANuE_) .or. &
   (id(1).eq.MuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.MuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.ANuT_) .or. (id(2).eq.TaM_ .and. id(1).eq.ANuT_) .or. &
   (id(1).eq.AUp_  .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.AUp_  .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_)) .or. &
   (id(1).eq.AChm_ .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.AChm_ .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_)) .or. &
   (id(1).eq.ATop_ .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.ATop_ .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_))      &
   ) .and. hel.lt.0) then
      CoupledVertex=Wm_
   elseif( (&
   (id(1).eq.ElM_ .and. id(2).eq.ElP_) .or. (id(2).eq.ElM_ .and. id(1).eq.ElP_) .or. &
   (id(1).eq.MuM_ .and. id(2).eq.MuP_) .or. (id(2).eq.MuM_ .and. id(1).eq.MuP_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.TaP_) .or. (id(2).eq.TaM_ .and. id(1).eq.TaP_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.AUp_) .or. (id(2).eq.Up_  .and. id(1).eq.AUp_) .or. &
   (id(1).eq.Dn_  .and. id(2).eq.ADn_) .or. (id(2).eq.Dn_  .and. id(1).eq.ADn_) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.AChm_) .or. (id(2).eq.Chm_ .and. id(1).eq.AChm_) .or. &
   (id(1).eq.Str_ .and. id(2).eq.AStr_) .or. (id(2).eq.Str_ .and. id(1).eq.Astr_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.ATop_) .or. (id(2).eq.Top_ .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.ABot_) .or. (id(2).eq.Bot_ .and. id(1).eq.ABot_)      &
   ) .and. hel.ne.0) then
      if(testAHcoupl.eq.1) then
         CoupledVertex=Pho_
      elseif(testAHcoupl.eq.2) then
         CoupledVertex=Hig_
      else
         CoupledVertex=Z0_
      endif
   elseif( (&
   (id(1).eq.NuE_ .and. id(2).eq.ANuE_) .or. (id(2).eq.NuE_ .and. id(1).eq.ANuE_) .or. &
   (id(1).eq.NuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.NuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.NuT_ .and. id(2).eq.ANuT_) .or. (id(2).eq.NuT_ .and. id(1).eq.ANuT_)      &
   ) .and. hel.lt.0) then
      CoupledVertex=Z0_ ! Only Z coupling to nuL-nubR
   else
      CoupledVertex=Not_a_particle_
   endif

   return
end function CoupledVertex

function CoupledVertexIsDiagonal(id,hel)
   implicit none
   integer, intent(in) :: id(1:2),hel
   logical :: CoupledVertexIsDiagonal

   if( (&
   (id(1).eq.ElP_ .and. id(2).eq.NuE_) .or. (id(2).eq.ElP_ .and. id(1).eq.NuE_) .or. & ! W+ -> l+ nu
   (id(1).eq.MuP_ .and. id(2).eq.NuM_) .or. (id(2).eq.MuP_ .and. id(1).eq.NuM_) .or. &
   (id(1).eq.TaP_ .and. id(2).eq.NuT_) .or. (id(2).eq.TaP_ .and. id(1).eq.NuT_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.ADn_) .or. (id(2).eq.Up_  .and. id(1).eq.ADn_) .or. & ! W+ -> qu qbd
   (id(1).eq.Chm_ .and. id(2).eq.AStr_) .or. (id(2).eq.Chm_ .and. id(1).eq.AStr_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.ABot_) .or. (id(2).eq.Top_ .and. id(1).eq.ABot_) .or. &
   (id(1).eq.ElM_ .and. id(2).eq.ANuE_) .or. (id(2).eq.ElM_ .and. id(1).eq.ANuE_) .or. & ! W- -> l- nub
   (id(1).eq.MuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.MuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.ANuT_) .or. (id(2).eq.TaM_ .and. id(1).eq.ANuT_) .or. &
   (id(1).eq.AUp_  .and. id(2).eq.Dn_ ) .or. (id(2).eq.AUp_  .and. id(1).eq.Dn_) .or. & ! W- -> qd qbu
   (id(1).eq.AChm_ .and. id(2).eq.Str_) .or. (id(2).eq.AChm_ .and. id(1).eq.Str_) .or. &
   (id(1).eq.ATop_ .and. id(2).eq.Bot_) .or. (id(2).eq.ATop_ .and. id(1).eq.Bot_) .or. &
   (id(1).eq.NuE_ .and. id(2).eq.ANuE_) .or. (id(2).eq.NuE_ .and. id(1).eq.ANuE_) .or. & ! Z -> nu nu
   (id(1).eq.NuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.NuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.NuT_ .and. id(2).eq.ANuT_) .or. (id(2).eq.NuT_ .and. id(1).eq.ANuT_)      &
   ) .and. hel.lt.0) then
      CoupledVertexIsDiagonal=.true.
   elseif( (&
   (id(1).eq.ElM_ .and. id(2).eq.ElP_) .or. (id(2).eq.ElM_ .and. id(1).eq.ElP_) .or. & ! Z -> ll
   (id(1).eq.MuM_ .and. id(2).eq.MuP_) .or. (id(2).eq.MuM_ .and. id(1).eq.MuP_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.TaP_) .or. (id(2).eq.TaM_ .and. id(1).eq.TaP_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.AUp_) .or. (id(2).eq.Up_  .and. id(1).eq.AUp_) .or. & ! Z -> qq
   (id(1).eq.Dn_  .and. id(2).eq.ADn_) .or. (id(2).eq.Dn_  .and. id(1).eq.ADn_) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.AChm_) .or. (id(2).eq.Chm_ .and. id(1).eq.AChm_) .or. &
   (id(1).eq.Str_ .and. id(2).eq.AStr_) .or. (id(2).eq.Str_ .and. id(1).eq.Astr_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.ATop_) .or. (id(2).eq.Top_ .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.ABot_) .or. (id(2).eq.Bot_ .and. id(1).eq.ABot_)      &
   ) .and. hel.ne.0) then
      CoupledVertexIsDiagonal=.true.
   else
      CoupledVertexIsDiagonal=.false.
   endif

   return
end function CoupledVertexIsDiagonal



! Counts *charged* leptons to be specific
FUNCTION CountLeptons( MY_IDUP )
implicit none
integer :: MY_IDUP(:),CountLeptons
integer :: i

   CountLeptons = 0
   do i = 1,size(MY_IDUP)
      if( IsALepton( MY_IDUP(i) ) ) CountLeptons=CountLeptons+1
   enddo


RETURN
END FUNCTION

! Counts light jets
FUNCTION CountJets( MY_IDUP )
implicit none
integer :: MY_IDUP(:),CountJets
integer :: i

   CountJets = 0
   do i = 1,size(MY_IDUP)
      if( IsAJet( MY_IDUP(i) ) ) CountJets=CountJets+1
   enddo


RETURN
END FUNCTION



FUNCTION SU2flip(Part)
implicit none
integer :: SU2flip
integer :: Part

  if( abs(Part).eq.Up_ ) then
      SU2flip = sign(1,Part)*Dn_
  elseif( abs(Part).eq.Dn_ ) then
      SU2flip = sign(1,Part)*Up_
  elseif( abs(Part).eq.Chm_ ) then
      SU2flip = sign(1,Part)*Str_
  elseif( abs(Part).eq.Str_ ) then
      SU2flip = sign(1,Part)*Chm_
  elseif( abs(Part).eq.Bot_ ) then
      SU2flip = sign(1,Part)*Top_
  elseif( abs(Part).eq.Top_ ) then
      SU2flip = sign(1,Part)*Bot_
  elseif( abs(Part).eq.ElP_ ) then
      SU2flip = sign(1,Part)*NuE_
  elseif( abs(Part).eq.MuP_ ) then
      SU2flip = sign(1,Part)*NuM_
  elseif( abs(Part).eq.TaP_ ) then
      SU2flip = sign(1,Part)*NuT_
  elseif( abs(Part).eq.NuE_ ) then
      SU2flip = sign(1,Part)*ElP_
  elseif( abs(Part).eq.NuM_ ) then
      SU2flip = sign(1,Part)*MuP_
  elseif( abs(Part).eq.NuT_ ) then
      SU2flip = sign(1,Part)*TaP_
  else
      print *, "Error: Invalid flavor in SU2flip ",Part
      stop
  endif

END FUNCTION


Function CouplToLHEWp(Part)
  implicit none
  logical :: CouplToLHEWp
  integer, intent(in) :: Part(1:2)

  CouplToLHEWp = &
     ((Part(1).eq.2 .and.Part(2).eq.-1) &
  .or.(Part(1).eq.2 .and.Part(2).eq.-3) &
  .or.(Part(1).eq.2 .and.Part(2).eq.-5) &
  .or.(Part(1).eq.4 .and.Part(2).eq.-1) &
  .or.(Part(1).eq.4 .and.Part(2).eq.-3) &
  .or.(Part(1).eq.4 .and.Part(2).eq.-5) &
  .or.(Part(1).eq.-1.and.Part(2).eq. 2) &
  .or.(Part(1).eq.-3.and.Part(2).eq. 2) &
  .or.(Part(1).eq.-5.and.Part(2).eq. 2) &
  .or.(Part(1).eq.-1.and.Part(2).eq. 4) &
  .or.(Part(1).eq.-3.and.Part(2).eq. 4) &
  .or.(Part(1).eq.-5.and.Part(2).eq. 4))

  return
end function CouplToLHEWp


Function CouplToLHEWm(Part)
  implicit none
  logical :: CouplToLHEWm
  integer, intent(in) :: Part(1:2)

  CouplToLHEWm = &
     ((Part(1).eq.-2.and.Part(2).eq. 1) &
  .or.(Part(1).eq.-2.and.Part(2).eq. 3) &
  .or.(Part(1).eq.-2.and.Part(2).eq. 5) &
  .or.(Part(1).eq.-4.and.Part(2).eq. 1) &
  .or.(Part(1).eq.-4.and.Part(2).eq. 3) &
  .or.(Part(1).eq.-4.and.Part(2).eq. 5) &
  .or.(Part(1).eq. 1.and.Part(2).eq.-2) &
  .or.(Part(1).eq. 3.and.Part(2).eq.-2) &
  .or.(Part(1).eq. 5.and.Part(2).eq.-2) &
  .or.(Part(1).eq. 1.and.Part(2).eq.-4) &
  .or.(Part(1).eq. 3.and.Part(2).eq.-4) &
  .or.(Part(1).eq. 5.and.Part(2).eq.-4))

  return
end function CouplToLHEWm





FUNCTION ChargeFlip(Part)
implicit none
integer :: ChargeFlip
integer :: Part

  if( (abs(Part).ge.1 .and. abs(Part).le.9)  .or. (abs(Part).ge.13 .and. abs(Part).le.16)) then! quarks, leptons, W's, neutrinos
      ChargeFlip = -Part
  elseif( (abs(Part).ge.10 .and. abs(Part).le.12) ) then!  glu,pho,Z0
      ChargeFlip = +Part
  elseif( Part.eq.25 .or. Part.eq.32 .or. Part.eq.39 ) then ! Higgs,Zprime,Graviton
      ChargeFlip = +Part
  else
      print *, "Error: Invalid flavor in ChargeFlip"
      stop
  endif

END FUNCTION


subroutine ComputeCKMElements(inVCKM_ud, inVCKM_us, inVCKM_cd, inVCKM_cs, inVCKM_ts, inVCKM_tb, inVCKM_ub, inVCKM_cb, inVCKM_td)
implicit none
real(8) :: inVCKM_ud
real(8) :: inVCKM_us
real(8) :: inVCKM_cd
real(8) :: inVCKM_cs
real(8) :: inVCKM_ts
real(8) :: inVCKM_tb
real(8), optional :: inVCKM_ub
real(8), optional :: inVCKM_cb
real(8), optional :: inVCKM_td
real(8) :: sumVsq(1:3),diffVsq

   VCKM_ud=inVCKM_ud
   VCKM_us=inVCKM_us
   VCKM_cd=inVCKM_cd
   VCKM_cs=inVCKM_cs
   VCKM_ts=inVCKM_ts
   VCKM_tb=inVCKM_tb

   if(present(inVCKM_ub)) then
      VCKM_ub = inVCKM_ub
   else
      diffVsq = 1d0-VCKM_ud**2-VCKM_us**2
      if (diffVsq.ge.0d0) then
         VCKM_ub = sqrt(diffVsq)
      else
         VCKM_ub = 0d0
      endif
   endif
   if(present(inVCKM_cb)) then
      VCKM_cb = inVCKM_cb
   else
      diffVsq = 1d0-VCKM_cd**2-VCKM_cs**2
      if (diffVsq.ge.0d0) then
         VCKM_cb = sqrt(diffVsq)
      else
         VCKM_cb = 0d0
      endif
   endif
   if(present(inVCKM_td)) then
      VCKM_td = inVCKM_td
   else
      diffVsq = 1d0-VCKM_td**2-VCKM_ts**2
      if (diffVsq.ge.0d0) then
         VCKM_tb = sqrt(diffVsq)
      else
         VCKM_tb = 0d0
      endif
   endif

   sumVsq(1) = VCKM_ub**2+VCKM_us**2+VCKM_ud**2
   VCKM_ud=VCKM_ud/sqrt(sumVsq(1))
   VCKM_us=VCKM_us/sqrt(sumVsq(1))
   VCKM_ub=VCKM_ub/sqrt(sumVsq(1))
   sumVsq(2) = VCKM_cb**2+VCKM_cs**2+VCKM_cd**2
   VCKM_cd=VCKM_cd/sqrt(sumVsq(2))
   VCKM_cs=VCKM_cs/sqrt(sumVsq(2))
   VCKM_cb=VCKM_cb/sqrt(sumVsq(2))
   sumVsq(3) = VCKM_tb**2+VCKM_ts**2+VCKM_td**2
   VCKM_td=VCKM_td/sqrt(sumVsq(3))
   VCKM_ts=VCKM_ts/sqrt(sumVsq(3))
   VCKM_tb=VCKM_tb/sqrt(sumVsq(3))

end subroutine ComputeCKMElements

subroutine SetDefaultCKM()
implicit none
   VCKM_ud = 0.974285d0
   VCKM_us = 0.225290d0
   VCKM_cs = 0.9734244d0
   VCKM_cd =-0.225182d0
   VCKM_tb = 0.99912367d0
   VCKM_ts =-0.040920069d0
   call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb)
end subroutine SetDefaultCKM

subroutine ComputeEWVariables()
implicit none

   ! Calculate fundamental couplings
   vev = 1.0d0/sqrt(Gf*sqrt(2.0d0))
   gwsq = 4.0d0 * M_W**2/vev**2
   sitW = sqrt(xw)
   twosc = sqrt(4d0*xw*(1d0-xw))

   esq = 4.0_dp*pi*alpha_QED
   overallCouplVffsq = esq/(twosc**2) ! ~= gwsq/4.0_dp/(1.0_dp-xw)
   ! Z couplings
   aR_lep =2d0*(T3lR-QlR*xw)
   aL_lep =2d0*(T3lL-QlL*xw)
   aR_neu =2d0*(T3nR-QnR*xw)
   aL_neu =2d0*(T3nL-QnL*xw)
   aR_QUp =2d0*(T3uR-QuR*xw)
   aL_QUp =2d0*(T3uL-QuL*xw)
   aR_QDn =2d0*(T3dR-QdR*xw)
   aL_QDn =2d0*(T3dL-QdL*xw)
   
   ! Zff couplings: adding user specified values (default 0.0) to SM values 
   aL_lep = aL_lep + daz_Lep_left 
   aR_lep = aR_lep + daz_Lep_right
   aL_neu = aL_neu + daz_Neu_left
   aR_neu = aR_neu + daz_Neu_right
   aL_QUp = aL_QUp + daz_Qup_left
   aR_Qup = aR_Qup + daz_Qup_right
   aL_QDn = aL_QDn + daz_QDn_left
   aR_QDn = aR_QDn + daz_QDn_right
   
   ! W couplings
   bL = sqrt(2d0*(1d0-xw))
   bR = 0d0
   ! A couplings
   cR_lep = -twosc*QlR
   cL_lep = -twosc*QlL
   cR_neu = -twosc*QnR
   cL_neu = -twosc*QnL
   cR_QUp = -twosc*QuR
   cL_QUp = -twosc*QuL
   cR_QDn = -twosc*QdR
   cL_QDn = -twosc*QdL

   ! Normalizations used in VH and VBF
   couplWffsq = gwsq/2.0_dp
   couplZffsq = gwsq/4.0_dp/(1.0_dp-xw)
   !couplZffsq = alpha_QED*4d0*pi/xw/4.0_dp/(1.0_dp-xw)
   couplAZff = -gwsq*sitW/2.0_dp/sqrt(1.0_dp-xw)
   couplAffsq = gwsq*xw


   

end subroutine ComputeEWVariables

subroutine ComputeQCDVariables()
implicit none
   gs = sqrt(alphas*4.0_dp*pi)
end subroutine ComputeQCDVariables



! QCD scale from MCFM
! Implementation into JHUGen by Ulascan Sarica, Dec. 2015
subroutine EvalAlphaS()
!   use ModParameters
   IMPLICIT NONE
#if useLHAPDF==1
!--- This is simply a wrapper to the LHAPDF implementation of the running coupling alphas, in the style of the native MCFM routine
   DOUBLE PRECISION alphasPDF
   REAL(DP) :: Q
      Q = Mu_Ren/GeV
      alphas=alphasPDF(Q)
#else
!     Evaluation of strong coupling constant alphas
!     Original Author: R.K. Ellis
!     q -- Scale at which alpha_s is to be evaluated
!     alphas_mz -- ModParameters value of alpha_s at the mass of the Z-boson
!     nloops_pdf -- ModParameters value of the number of loops (1,2, or 3) at which the beta function is evaluated to determine running.
!     If you somehow need a more complete implementation, check everything at or before commit 28472c5bfee128dde458fd4929b4d3ece9519ab8
   INTEGER, PARAMETER :: NF6=6
   INTEGER, PARAMETER :: NF5=5
   INTEGER, PARAMETER :: NF4=4
   INTEGER, PARAMETER :: NF3=3
   INTEGER, PARAMETER :: NF2=2
   INTEGER, PARAMETER :: NF1=1
   INTEGER, PARAMETER :: NF0=0

      IF (Mu_Ren .LE. 0d0) THEN
         WRITE(6,*) 'ModParameters::EvalAlphaS: Mu_Ren .le. 0, Mu_Ren (GeV) = ',(Mu_Ren*GeV)
         stop
      ENDIF
      IF (nQflavors_pdf .GT. NF6 .OR. nQflavors_pdf .LT. NF0) THEN
         WRITE(6,*) 'ModParameters::EvalAlphaS: nQflavors_pdf has to be between 0 and 6. nQflavors_pdf = ',nQflavors_pdf
         stop
      ENDIF
      IF (nQflavors_pdf .NE. NF5) THEN
         WRITE(6,*) 'ModParameters::EvalAlphaS: nQflavors_pdf invalid, nQflavors_pdf = ',nQflavors_pdf
         WRITE(6,*) 'ModParameters::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
         stop
      ENDIF
      IF (nloops_pdf .NE. 1) THEN
         WRITE(6,*) 'ModParameters::EvalAlphaS: nloops_pdf invalid, nloops_pdf = ',nloops_pdf
         WRITE(6,*) 'ModParameters::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
         stop
      ENDIF

      alphas=alphas_mz/(1.0_dp+alphas_mz*B0_PDF(nQflavors_pdf)*2.0_dp*dlog((Mu_Ren/zmass_pdf)))
#endif
      ! Calculate the derived couplings
      call ComputeQCDVariables()
   RETURN
end subroutine EvalAlphaS



!-----------------------------------------------------------------------------
!
      real(8) function massfrun(mf,scale)
!
!-----------------------------------------------------------------------------
!
!       This function returns the 'nloop' value of a MSbar fermion mass
!       at a given scale.
!
!       INPUT: mf    = MSbar mass of fermion at MSbar fermion mass scale
!              scale = scale at which the running mass is evaluated
!              asmz  = AS(MZ) : this is passed to alphas(scale,asmz,2)
!              nloop = # of loops in the evolutionC
!
!       COMMON BLOCKS: COMMON/QMASS/CMASS,BMASS,TMASS
!                      contains the MS-bar masses of the heavy quarks.
!
!       EXTERNAL:      double precision alphas(scale,asmz,2)
!
!-----------------------------------------------------------------------------
!
!      use ModParameters
      implicit none
!
!     ARGUMENTS
!
      real(8), intent(in) :: mf, scale
      real(8) scale_temp_ren
      !integer , intent(in) :: nloop
!
!     LOCAL
!
      real(8)  beta0, beta1,gamma0,gamma1
      real(8)  as,asmf,l2
      integer  nfrun

      scale_temp_ren=Mu_Ren
!
!     EXTERNAL
!
!      double precision  alphas
!      external          alphas
!
!     COMMON
!
!      real *8      cmass,bmass,tmass
!      COMMON/QMASS/CMASS,BMASS,TMASS
!
!     CONSTANTS
!
!      double precision  One, Two, Three, Pi
      !parameter( One = 1d0, Two = 2d0, Three = 3d0 )
      !parameter( Pi = 3.14159265358979323846d0)

      if ( mf.gt.m_top ) then
         nfrun = 6
      else
         nfrun = 5
      end if

      beta0 = ( 11d0 - 2d0/3d0 *nfrun )/4d0
      !beta1 = ( 102d0  - 38d0/3d0*nf )/16d0
      gamma0= 1d0
      !gamma1= ( 202d0/3d0  - 20d0/9d0*nf )/16d0
      !A1    = -beta1*gamma0/beta0**2+gamma1/beta0
      Mu_Ren=scale
      call EvalAlphaS()
      as=alphas

      Mu_Ren=mf
      call EvalAlphaS()
      asmf=alphas
      !l2    = (1d0+A1*as/Pi)/(one+A1*asmf/Pi)

      massfrun = mf * (as/asmf)**(gamma0/beta0)

      !if(nloop.eq.2) massfrun=massfrun*l2

      Mu_Ren=scale_temp_ren
      call EvalAlphaS()

      return
      end function massfrun









!========================================================================
!---- THESE ARE POLARIZATION ROUTINES
  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
  implicit none
    complex(dp), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------
    integer :: pol
    real(dp) :: p0,px,py,pz
    real(dp) :: pv,ct,st,cphi,sphi
    complex(dp) :: pol_mless(4)

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
!^^^END


    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0_dp-ct**2))

    if (st < tol) then
       cphi=1.0_dp
       sphi=0.0_dp
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0_dp) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif

    pol_mless(1)=czero
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless

  function pol_mless2(p,i,out)
  implicit none
    integer, intent(in) :: i
    complex(dp), intent(in) :: p(4)
    character(len=*), intent(in):: out
    complex(dp) :: pol_mless2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mless2 = pol_mless(p,i,outgoing=.true.)
    else
       pol_mless2 = pol_mless(p,i,outgoing=.false.)
    endif
  end function pol_mless2

  function pol_dk2mom(plepton,antilepton,i,outgoing)
  implicit none
    integer, intent(in) :: i
    integer :: j
    complex(dp), intent(in) :: plepton(1:4),antilepton(1:4)
    logical, intent(in),optional :: outgoing
    complex(dp) :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq


    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2

    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)
    !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
    !   This is an expression for (-i)/qsq* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_dk2mom(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_dk2mom(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_dk2mom(3)=-ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_dk2mom(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))


    do j=1,4
       pol_dk2mom(j)=pol_dk2mom(j)/qsq
    enddo

    ! -- do nothing in this case
    if (present(outgoing)) then
       !if (outgoing) pol_dk2mom = conjg(pol_dk2mom)
    endif

  end function pol_dk2mom

! -- massive vector polarization subroutine
   function pol_mass(p,i,outgoing)
   implicit none
   integer, intent(in) :: i
   integer :: pol
   complex(8), intent(in) :: p(4)
   logical, intent(in),optional :: outgoing
   complex(8) :: pol_mass(4)
   complex(8) :: msq,m
   real(8) :: p0,px,py,pz, pv,pvsq
   real(8) :: ct,st,cphi,sphi

      p0=dreal(p(1))
      px=dreal(p(2))
      py=dreal(p(3))
      pz=dreal(p(4))

      pv=px**2 + py**2 + pz**2
      m=p0**2-pv
      m=sqrt(m)
      pv=sqrt(pv)

      if(cdabs(pv/m).lt.1d-8) then
         if(i.eq.0) then
            pol_mass(1:3)=czero
            pol_mass( 4 )=cone
            return
         endif
         ct = 1d0; st=0d0
      else
         ct= pz/pv
         st= dsqrt(dabs(1.0d0-ct**2))
      endif


      if (st .lt. 1D-15) then
         cphi=1.0d0
         sphi=0d0
      else
         cphi= px/pv/st
         sphi= py/pv/st
      endif


!     i=0 is longitudinal polarization
!     the following ifstatement distinguishes between
!     positive and negative energies
      if ( p0 .gt. 0.0d0) then
         pol=i
      else
         pol=-i
      endif

      ! -- take complex conjugate for outgoing
      if (present(outgoing)) then
         if (outgoing) pol = -pol
      endif

      if(pol.eq.-1 .or. pol.eq.1) then
         pol_mass(1)=czero
         pol_mass(2)=(ct*cphi-pol*ci*sphi)/sqrt2
         pol_mass(3)=(ct*sphi+pol*ci*cphi)/sqrt2
         pol_mass(4)=-st/sqrt2
      else if(pol.eq.0) then
         pol_mass(1)= pv/m
         pol_mass(2)= p0/m/pv*px
         pol_mass(3)= p0/m/pv*py
         pol_mass(4)= p0/m/pv*pz
      else
         print *,"wrong helicity setting in pol_mass"
         stop
      endif

   end function pol_mass

  function pol_mass2(p,i,out)
   implicit none
    integer, intent(in) :: i
    complex(dp), intent(in) :: p(4)
    character(len=*), intent(in):: out
    complex(dp) :: pol_mass2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mass2 = pol_mass(p,i,outgoing=.true.)
    else
       pol_mass2 = pol_mass(p,i,outgoing=.false.)
    endif
  end function pol_mass2

!---- THESE ARE SPINOR ROUTINES
!     ubar spinor, massless
  function ubar0(p,i)
   implicit none
    complex(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    complex(dp) :: fc, fc2
    real(dp)    :: p0,px,py,pz,mass


    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    mass=dsqrt(dabs(p0**2-px**2-py**2-pz**2))
    if( mass.lt.1d-4 ) mass=0d0


    fc2 = p0 + pz
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then
       if (i.eq.1) then
          ubar0(1)=czero
          ubar0(2)=czero
          ubar0(3)=fc
          ubar0(4)=(px-ci*py)/fc
       elseif (i.eq.-1) then
          ubar0(1)=(px+ci*py)/fc
          ubar0(2)=-fc
          ubar0(3)=czero
          ubar0(4)=czero
       else
          stop 'ubar0: i out of range'
       endif
    else
       if (i.eq.1) then
          ubar0(1) = czero
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = sqrt(cone*two*p0)
       elseif (i.eq.-1) then
          ubar0(1) = sqrt(cone*(two*p0))
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = czero
       else
          stop 'ubar0: i out of range'
       endif
    endif


!       if (i.eq.1) then
!           ubar0(1)=dcmplx(mass,0d0)/fc
!           ubar0(2)=czero
!           ubar0(3)=fc
!           ubar0(4)=dcmplx(px,-py)/fc
!       elseif (i.eq.-1) then
!           ubar0(1)=dcmplx(px,py)/fc
!           ubar0(2)=-fc
!           ubar0(3)=czero
!           ubar0(4)=-dcmplx(mass,0d0)/fc
!        else
!           stop 'ubar0: i out of range'
!       endif



  end function ubar0

  ! -- v0  spinor, massless
  function v0(p,i)
   implicit none
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    complex(dp) :: v0(4)
    complex(dp) :: fc2, fc
    real(dp)    :: p0,px,py,pz,mass

    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    mass=dsqrt(dabs(p0**2-px**2-py**2-pz**2))
    if( mass.lt.1d-4 ) mass=0d0


    fc2 = p0 + pz
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then
       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=(px-ci*py)/fc
          v0(4)=-fc
       elseif (i.eq.-1) then
          v0(1)=fc
          v0(2)=(px+ci*py)/fc
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif
    else
       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=sqrt(cone*two*p0)
          v0(4)=czero
       elseif (i.eq.-1) then
          v0(1)=czero
          v0(2)=sqrt(cone*two*p0)
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif
    endif


!       if (i.eq.+1) then
!           v0(1)=czero
!           v0(2)=dcmplx(mass,0d0)/fc
!           v0(3)=dcmplx(px,-py)/fc
!           v0(4)=-fc
!       elseif (i.eq.-1) then
!           v0(1)=fc
!           v0(2)=dcmplx(px,py)/fc
!           v0(3)=dcmplx(-mass,0d0)/fc
!           v0(4)=czero
!        else
!           stop 'v0: i out of range'
!       endif



  end function v0

subroutine spinoru(p,za,zb,s)
!---Calculate spinor products
!---taken from MCFM & modified by R. Rontsch, May 2015
!---extended to deal with negative energies ie with all momenta outgoing
!---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
!---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      real(dp) :: p(:,:),two
      integer, parameter :: mxpart=14
      complex(dp):: c23(mxpart),f(mxpart),rt(mxpart),za(:,:),zb(:,:),czero,cone,ci
      real(dp)   :: s(:,:)
      integer i,j,N

      N = size(p,1)
!       if (size(p,1) .ne. N) then
!          call Error("spinorz: momentum mismatch",size(p,1))
!       endif
      two=2d0
      czero=dcmplx(0d0,0d0)
      cone=dcmplx(1d0,0d0)
      ci=dcmplx(0d0,1d0)


!---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czero
         zb(j,j)=za(j,j)

!-----positive energy case
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(p(j,4)+p(j,1))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
!-----negative energy case
            rt(j)=dsqrt(-p(j,4)-p(j,1))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=ci
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

         if (abs(s(i,j)).lt.1d-5) then
         zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
         else
         zb(i,j)=-dcmplx(s(i,j))/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

    end subroutine spinoru




!========================================================================

!========================================================================
! Common initialization functions that may be called multiple times if needed
! Check arXiv:1604.06792 for the parameters
subroutine InitCOLLIER(Nmax, Rmax)
#if useCollier==1
use COLLIER
implicit none
integer, intent(in) :: Nmax, Rmax
integer :: supNmax, supRmax
   supNmax = max(Nmax, Collier_maxNLoopProps)
   supRmax = max(Rmax, Collier_maxRank)
   if ((supNmax .gt. Collier_maxNLoopProps) .or. (supRmax .gt. Collier_maxRank)) then
      call Init_cll(supNmax,supRmax,'')
      call setMode_cll(1)
   endif
#else
implicit none
integer, intent(in) :: Nmax, Rmax
   return
#endif
end subroutine



END MODULE
