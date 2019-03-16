! If you use this program please cite Phys.Rev. D81 (2010) 075022; arXiv:1001.3396  [hep-ph],
!                                     Phys.Rev. D86 (2012) 095031; arXiv:1208.4018  [hep-ph],
!                                     Phys.Rev. D89 (2014) 035007; arXiv:1309.4819  [hep-ph],
!                                 and Phys.Rev. D94 (2016) 055023; arXiv:1606.03107 [hep-ph].
PROGRAM JHUGenerator
#if compiler==1
use ifport
#endif
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
use ModPMZZ
use ModHashCollection
implicit none
real(8) :: VG_Result,VG_Error

   call SetJHUGenDefaults()
   call GetCommandlineArgs()
   call InitProcessScaleSchemes()
   call InitPDFs()!
   call InitHisto()
   call InitParameters()
   call InitProcess()
   call InitVegas()
   call InitRandomSeed()
   call OpenFiles()
   call PrintLogo(io_stdout, "JHU Generator")
   call PrintLogo(io_LogFile, "JHU Generator")
   call WriteParameters(io_stdout)
   call WriteParameters(io_LogFile)
   if ( .not. ReadLHEFile .and. .not. ConvertLHEFile .and. .not.(CalculatesXsec(Process) .and. unweighted) ) then
      call InitOutput(1d0, 1d14)   !for VBF/HJJ the cross section is calculated, so use that in the <init> block
   endif
   call SetupHashes()
#if linkMELA==1
   call SetupMCFM(Process)
#endif
   write(io_stdout,*) " Running"
   if( ConvertLHEFile ) then
        call StartConvertLHE(VG_Result,VG_Error)
   elseif( ReadLHEFile ) then
        call StartReadLHE_NEW(VG_Result,VG_Error)
   elseif( DoPrintPMZZ ) then
        call PrintMZZdistribution()
   else
        if( Process.eq.0 .or. Process.eq.1  .or. Process.eq.2 .or. Process.eq.80 &
                         .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.66 .or. Process.eq.67 .or. Process.eq.68 .or. Process.eq.69 &
                         .or. Process.eq.90 .or. Process.eq.110 .or. Process.eq.111 .or. Process.eq.112 .or. Process.eq.113 .or. Process.eq.114 ) then
           call StartVegas_NEW(VG_Result,VG_Error)
        else
           call StartVegas(VG_Result,VG_Error)
        endif
   endif
   call WriteHisto(VG_Result,VG_Error,time_end-time_start)
   call FinalizeOutput()
   call CloseFiles()
   write(io_stdout,*) "Done"


END PROGRAM


! !       Scheme=kRenFacScheme_default: Defaults of each process with Mu_Fact==Mu_Ren
! !       Scheme=+-kRenFacScheme_mhstar: Scale ~ m_Hstar (or m_VV in bkg.); + for running, - for fixed
! !       Below, J stands for the particles immediately associated to the Higgs (e.g. ttH, VH, VBF, tqH)
! !       Scheme=+-kRenFacScheme_mjjhstar: Scale ~ m_JJHstar (or m_JJVV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mjj_mhstar: Scale ~ m_JJ + m_Hstar (or m_JJ+m_VV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mj_mj_mhstar: Scale ~ m_J + m_J + m_Hstar (or m_J+m_J+m_VV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mjj: Scale ~ m_JJ; + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mj_mj: Scale ~ m_J + m_J; + for running, - for fixed
! !       Scheme=+ kRenScheme_maxpTj: Scale ~ pT of the hardest associated jet
! !       Scheme=+ kRenScheme_minpTj: Scale ~ pT of the softest associated jet
! !       Below, J stands for either the top in JJ associated production (e.g. tqH), or the single jet (Hj)
! !       Scheme=+-kRenFacScheme_mjhstar: Scale ~ m_JHstar (or m_JVV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mj_mhstar: Scale ~ m_J + m_Hstar (or m_J+m_VV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mj: Scale ~ m_J for the heavy jet (ie. t or b); + for running, - for fixed
subroutine InitProcessScaleSchemes() ! If schemes are set to default, reset to the appropriate numbers
   use ModParameters
   use ModMisc
#if useCollier==1
   use COLLIER
#endif
   implicit none
   integer :: Nmax,rmax

      if( .not.                 &
         (                      &
            Process.eq. 0 .or.  &
            Process.eq. 1 .or.  &
            Process.eq. 2 .or.  &
            Process.eq.50 .or.  &
            Process.eq.51 .or.  &
            Process.eq.52 .or.  &
            Process.eq.60 .or.  &
            Process.eq.61 .or.  &
            Process.eq.62 .or.  &
            Process.eq.66 .or.  &
            Process.eq.67 .or.  &
            Process.eq.68 .or.  &
            Process.eq.69 .or.  &
            Process.eq.76 .or.  &
            Process.eq.77 .or.  &
            Process.eq.78 .or.  &
            Process.eq.80 .or.  &
            Process.eq.90 .or.  &
            Process.eq.110 .or. &
            Process.eq.111 .or. &
            Process.eq.112 .or. &
            Process.eq.113 .or. &
            Process.eq.114      &
         )                      &
      ) call Error("main::InitProcessScaleSchemes: Renormalization and factorization schemes are not implemented for process",Process)

      if(FacScheme.eq.kRenFacScheme_default) then
         if( &
         Process.eq. 0 .or. & !- ggH spin-0
         Process.eq. 1 .or. & !- ggH spin-1
         Process.eq. 2      & !- ggH spin-2
         ) then
            FacScheme = +kRenFacScheme_mhstar
            MuFacMultiplier = 0.5d0
         elseif( &
         Process.eq.60 .or. & !- HVBF without decays
         Process.eq.61 .or. & !- Hjj, gluon fusion
         Process.eq.62 .or. & !- Hj, gluon fusion
         Process.eq.50 .or. & !- VHiggs
         Process.eq.51 .or. & !- VH
         Process.eq.52      & !- HH
         ) then
            FacScheme = -kRenFacScheme_mhstar
            MuFacMultiplier = 1d0
         elseif( &
         Process.eq.66 .or. & !- VVHVV with decays
         Process.eq.67 .or. & !- VVVV with decays
         Process.eq.68 .or. & !- VVVV+VVHVV with decays
         Process.eq.69      & !- QCD JJ bkg with decays
         ) then
            FacScheme = +kRenFacScheme_mhstar
            MuFacMultiplier = 1d0
         elseif( &
         Process.eq.80 .or. & !- ttbar+H
         Process.eq.90      & !- bbbar+H
         ) then
            FacScheme = -kRenFacScheme_mj_mj_mhstar
            MuFacMultiplier = 0.5d0
         elseif( &
         Process.eq.110 .or. & !- t+H
         Process.eq.111 .or. & !- tb+H
         Process.eq.112 .or. & !- t+H s-channel
         Process.eq.113 .or. & !- tb+H s-channel
         Process.eq.114      & !- sum of all channels
         ) then
            FacScheme = -kRenFacScheme_mj_mhstar
            MuFacMultiplier = 0.25d0
         endif
      endif

      if(RenScheme.eq.kRenFacScheme_default) then
         if( &
         Process.eq. 0 .or. & !- ggH spin-0
         Process.eq. 1 .or. & !- ggH spin-1
         Process.eq. 2      & !- ggH spin-2
         ) then
            RenScheme = +kRenFacScheme_mhstar
            MuRenMultiplier = 0.5d0
         elseif( &
         Process.eq.60 .or. & !- HVBF without decays
         Process.eq.61 .or. & !- Hjj, gluon fusion
         Process.eq.62 .or. & !- Hj, gluon fusion
         Process.eq.50 .or. & !- VHiggs
         Process.eq.51 .or. & !- VH
         Process.eq.52      & !- HH
         ) then
            RenScheme = -kRenFacScheme_mhstar
            MuRenMultiplier = 1d0
         elseif( &
         Process.eq.66 .or. & !- VVHVV with decays
         Process.eq.67 .or. & !- VVVV with decays
         Process.eq.68 .or. & !- VVVV+VVHVV with decays
         Process.eq.69      & !- QCD JJ bkg with decays
         ) then
            RenScheme = +kRenFacScheme_mhstar
            MuRenMultiplier = 1d0
         elseif( &
         Process.eq.80 .or. & !- ttbar+H
         Process.eq.90      & !- bbbar+H
         ) then
            RenScheme = -kRenFacScheme_mj_mj_mhstar
            MuRenMultiplier = 0.5d0
         elseif( &
         Process.eq.110 .or. & !- t+H
         Process.eq.111 .or. & !- tb+H
         Process.eq.112 .or. & !- t+H s-channel
         Process.eq.113 .or. & !- tb+H s-channel
         Process.eq.114      & !- sum of all channels
         ) then
            RenScheme = -kRenFacScheme_mj_mhstar
            MuRenMultiplier = 0.25d0
         endif
      endif

      ! H+2j MEs
      if( &
         (                     &
            Process.eq.50 .or. &
            Process.eq.51 .or. &
            Process.eq.52 .or. &
            Process.eq.60 .or. &
            Process.eq.61 .or. &
            Process.eq.66 .or. &
            Process.eq.67 .or. &
            Process.eq.68 .or. &
            Process.eq.69 .or. &
            Process.eq.80 .or. &
            Process.eq.90      &
         ) .and. (             &
            (abs(FacScheme).eq.kRenFacScheme_mjhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj) .or. &
            (abs(RenScheme).eq.kRenFacScheme_mjhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj)      &
         )                     &
      ) call Error("ttH, bbH, HJJ, JJQCD, VBF and VH processes cannot distinguish the outgoing partons. Choose a different renormalization or factorization scheme.")

      if( &
         (                     &
            Process.eq.50 .or. &
            Process.eq.51 .or. &
            Process.eq.52 .or. &
            Process.eq.60 .or. &
            Process.eq.61 .or. &
            Process.eq.66 .or. &
            Process.eq.67 .or. &
            Process.eq.68 .or. &
            Process.eq.69      &
         ) .and. (             &
            (abs(FacScheme).eq.kRenFacScheme_mj_mj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj_mj) .or. &
            (abs(RenScheme).eq.kRenFacScheme_mj_mj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj_mj)      &
         )                     &
      ) call Error("HJJ, JJQCD, VBF and VH processes outgoing partons are mostly massless, and alpha_S at a scale ~0 GeV is very unstable. Choose a different renormalization or factorization scheme (e.g. kRenFacScheme_mhstar).")

      if( &
         (                     &
            Process.eq.50 .or. &
            Process.eq.51      &
         ) .and. (             &
            (abs(FacScheme).eq.kRenFacScheme_maxpTj) .or. (abs(RenScheme).eq.kRenFacScheme_maxpTj) .or.      &
            (abs(FacScheme).eq.kRenFacScheme_minpTj) .or. (abs(RenScheme).eq.kRenFacScheme_minpTj)           &
         ) .and. (             &
            DecayMode1.ne.1 .and. DecayMode2.ne.5 &
         )                     &
      ) call Error("For VH with V decay to something other than qq, can't base the scale on the pT of the a jet.  Choose a different renormalization or factorization scheme.")

      ! H+1j Me
      if( &
         (                     &
            Process.eq.62      &
         ) .and. (             &
            (FacScheme.eq.-kRenFacScheme_mjhstar) .or. (FacScheme.eq.-kRenFacScheme_mj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj) .or. &
            (RenScheme.eq.-kRenFacScheme_mjhstar) .or. (RenScheme.eq.-kRenFacScheme_mj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj) .or. &
            (abs(FacScheme).eq.kRenFacScheme_mjjhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mjj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mjj) .or. (abs(FacScheme).eq.kRenFacScheme_mj_mj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj_mj) .or. &
            (abs(RenScheme).eq.kRenFacScheme_mjjhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mjj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mjj) .or. (abs(RenScheme).eq.kRenFacScheme_mj_mj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj_mj)      &
         )                     &
      ) call Error("Invalid scheme for the HJ gluon fusion process. Choose a different renormalization or factorization scheme.")

      ! H+0j Me
      if( &
         (                     &
            Process.eq. 0 .or. & !- ggH spin-0
            Process.eq. 1 .or. & !- ggH spin-1
            Process.eq. 2      & !- ggH spin-2
         ) .and. (             &
            (abs(FacScheme).ne.kRenFacScheme_mhstar) .or. (abs(RenScheme).ne.kRenFacScheme_mhstar)      &
         )                     &
      ) call Error("Invalid scheme for the H+0J processes. Choose a different renormalization or factorization scheme.")


   return
end subroutine


SUBROUTINE SetJHUGenDefaults()
   use ModParameters
   implicit none

   call SetDefaultCKM()

   Collider=1
   VegasIt1=-1
   VegasNc0=-1
   VegasNc1=-1
   VegasNc2=-1
   Process = 0   ! select 0, 1 or 2 to represent the spin of the resonance
   PChannel=2
   DecayMode1=0  ! Z/W+
   DecayMode2=0  ! Z/W-
! !       DecayMode=0:  Z --> l+ l- (l=e,mu)
! !       DecayMode=1:  Z --> q qbar (q=u,d,c,s,b)
! !       DecayMode=2:  Z --> tau+ tau-
! !       DecayMode=3:  Z --> nu nubar (nu=nu_e,nu_mu,nu_tau)
! !       DecayMode=4:  W --> l nu_l (l=e,mu)
! !       DecayMode=5:  W --> q qbar' (q=u,c, qbar'=d,s)
! !       DecayMode=6:  W --> tau nu_tau
! !       DecayMode=7:  photon
! !       DecayMode=8:  Z --> l+ l- (l=e,mu,tau)
! !       DecayMode=9:  Z --> anything
! !       DecayMode=10: W --> l nu_l (l=e,mu,tau)
! !       DecayMode=11: W --> anything
   TopDecays=-1
   TauDecays=-1
   HbbDecays = .false.
   H_DK = .false.
   VH_PC="gg"
   alpha_dip=1d0
   Unweighted =.true.
   MuFacMultiplier = 1d0
   MuRenMultiplier = 1d0
   FacScheme = kRenFacScheme_default
   RenScheme = kRenFacScheme_default
   OffShellReson=.true.
   OffShellV1=.true.
   OffShellV2=.true.

   LHEProdFile=""
   ReadLHEFile=.false.
   ConvertLHEFile=.false.
   ReadCSmax=.false.
   ReadPMZZ = .false.
   PMZZFile="PMZZdistribution.out"
   PMZZEvals=-1
   DoPrintPMZZ = .false.
   PrintPMZZIntervals = 20
   GenerateEvents=.false.
   RequestNLeptons = -1
   RequestNJets = -1
   RequestOS=-1
   RequestOSSF=-1
   CountTauAsAny = .true.
   WriteFailedEvents=0
   VBFoffsh_run = -1

END SUBROUTINE SetJHUGenDefaults

SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
use ModMisc
use ModCommandLine
implicit none
character :: arg*(500)
integer :: NumArgs,NArg
logical :: help, DryRun, success, SetLastArgument, interfSet
logical :: SetRenScheme, SetMuRenMultiplier, SetFacScheme, SetMuFacMultiplier
logical :: SetAnomalousSpin0gg, Setghg2, SetAnomalousSpin0VV, Setghz1
logical :: SetZZcoupling, SetZZprimecoupling, SetZprimeZprimecoupling
logical :: SetZgammacoupling, Setgammagammacoupling, SetZprimegammacoupling
logical :: SetWWcoupling, SetWWprimecoupling, SetWprimeWprimecoupling
logical :: SetAnomalousSpin1qq, Setspin1qqleft, Setspin1qqright, SetSpin1VV
logical :: SetAnomalousSpin2gg, SetAnomalousSpin2qq, Setspin2qqleft, Setspin2qqright, SetSpin2VV
logical :: SetAnomalousHff, Setkappa
logical :: SetZprimeff, SetWprimeff, SetHZprime, SetHWprime
logical :: SetMZprime, SetGaZprime, SetMWprime, SetGaWprime
logical :: SetCKM,SetCKMub,SetCKMcb,SetCKMtd
logical :: SetpTjetcut, Setetajetcut, Setdetajetcut, SetdeltaRcut
logical :: SetpTlepcut, Setetalepcut, SetMPhotonCutoff
logical :: SetColliderEnergy
logical :: Setm2l_min,Setm2l_max,SetmVH_min,SetmVH_max,SetpTHcut
logical :: SetCSmaxFile, SetVBFoffsh_run
integer :: i
type(SaveValues) :: tosave, oldsavevalues

   help = .false.
   DryRun = .false.

#if useLHAPDF==1
   LHAPDFString = ""
   LHAPDFMember = 0
#else
   PDFSet=1      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40
#endif

   WidthScheme=-1
   WidthSchemeIn=-1

   interfSet = .false.
   SetMuFacMultiplier = .false.
   SetMuRenMultiplier = .false.
   SetFacScheme = .false.
   SetRenScheme = .false.

   SetCKM=.false.
   SetCKMub=.false.
   SetCKMcb=.false.
   SetCKMtd=.false.
   SetAnomalousSpin0gg=.false.
   Setghg2=.false.
   SetAnomalousSpin0VV=.false.
   Setghz1=.false.
   SetZZcoupling=.false.
   SetZZprimecoupling=.false.
   SetZprimeZprimecoupling=.false.
   SetZgammacoupling=.false.
   Setgammagammacoupling=.false.
   SetZprimegammacoupling=.false.
   SetWWcoupling=.false.
   SetWWprimecoupling=.false.
   SetWprimeWprimecoupling=.false.
   SetAnomalousSpin1qq=.false.
   Setspin1qqleft=.false.
   Setspin1qqright=.false.
   Setspin1VV=.false.
   SetAnomalousSpin2gg=.false.
   SetAnomalousSpin2qq=.false.
   Setspin2qqleft=.false.
   Setspin2qqright=.false.
   Setspin2VV=.false.
   SetAnomalousHff=.false.
   Setkappa=.false.

   SetHZprime=.false.
   SetZprimeff=.false.
   SetMZprime=.false.
   SetGaZprime=.false.
   SetHWprime=.false.
   SetWprimeff=.false.
   SetMWprime=.false.
   SetGaWprime=.false.

   SetpTjetcut=.false.
   Setetajetcut=.false.
   Setdetajetcut=.false.
   SetdeltaRcut=.false.

   SetpTlepcut=.false.
   Setetalepcut=.false.
   SetMPhotonCutoff=.false.

   SetColliderEnergy=.false.

   SetCSmaxFile=.false.
   SetVBFoffsh_run = .false.

   DataFile="./data/output"

   tosave = SaveValuesConstructor()

#if compiler==1
   NumArgs = NArgs()-1
#elif compiler==2
   NumArgs = COMMAND_ARGUMENT_COUNT()
#endif

   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    success = .false.
    call ReadCommandLineArgument(arg, "help", success, help)
    if( help ) then
        call PrintCommandLineArgs()
        stop 1
    endif
    !ReadCommandLineArgument is overloaded, it puts the value into the last argument
    ! by detecting the type.  It also sets success to .true. if the argument name (before =)
    ! is correct.
    call ReadCommandLineArgument(arg, "DryRun", success, DryRun)
    call ReadCommandLineArgument(arg, "Collider", success, Collider, tosave=tosave)
    call ReadCommandLineArgument(arg, "ColliderEnergy", success, Collider_Energy, multiply=1000*GeV, success2=SetColliderEnergy, tosave=tosave)
    call ReadCommandLineArgument(arg, "epPolarization", success, POL_A, tosave=tosave)
    call ReadCommandLineArgument(arg, "emPolarization", success, POL_B, tosave=tosave)
#if useLHAPDF==1
    call ReadCommandLineArgument(arg, "LHAPDF", success, LHAPDFString, tosave=tosave)
    call ReadCommandLineArgument(arg, "LHAPDFMem", success, LHAPDFMember, tosave=tosave)
#else
    call ReadCommandLineArgument(arg, "PDFSet", success, PDFSet, tosave=tosave)
#endif
    call ReadCommandLineArgument(arg, "MReso", success, M_Reso, multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "GaReso", success, Ga_Reso, multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "MReso2", success, M_Reso2, multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "GaReso2", success, Ga_Reso2, multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "ctauReso", success, HiggsDecayLengthMM, tosave=tosave)
    call ReadCommandLineArgument(arg, "VegasNc0", success, VegasNc0)
    call ReadCommandLineArgument(arg, "VegasNc1", success, VegasNc1)
    call ReadCommandLineArgument(arg, "VegasNc2", success, VegasNc2)
    call ReadCommandLineArgument(arg, "PChannel", success, PChannel, tosave=tosave)
    call ReadCommandLineArgument(arg, "VH_PC", success, VH_PC, tosave=tosave)
    call ReadCommandLineArgument(arg, "alpha_dip", success, alpha_dip, tosave=tosave)
    call ReadCommandLineArgument(arg, "DataFile", success, DataFile)
    call ReadCommandLineArgument(arg, "CSmaxFile", success, CSmaxFile, success2=SetCSmaxFile)
    call ReadCommandLineArgument(arg, "Process", success, Process, tosave=tosave)
    call ReadCommandLineArgument(arg, "DecayMode1", success, DecayMode1, tosave=tosave)
    call ReadCommandLineArgument(arg, "DecayMode2", success, DecayMode2, tosave=tosave)
    call ReadCommandLineArgument(arg, "FacScheme", success, FacScheme, success2=SetFacScheme, tosave=tosave)
    call ReadCommandLineArgument(arg, "RenScheme", success, RenScheme, success2=SetRenScheme, tosave=tosave)
    call ReadCommandLineArgument(arg, "MuFacMultiplier", success, MuFacMultiplier, success2=SetMuFacMultiplier, tosave=tosave)
    call ReadCommandLineArgument(arg, "MuRenMultiplier", success, MuRenMultiplier, success2=SetMuRenMultiplier, tosave=tosave)
    call ReadCommandLineArgument(arg, "TopDK", success, TopDecays, tosave=tosave)
    call ReadCommandLineArgument(arg, "TauDK", success, TauDecays, tosave=tosave)
    call ReadCommandLineArgument(arg, "HbbDK", success, H_DK, tosave=tosave)
    call ReadCommandLineArgument(arg, "HbbDK", success, HbbDecays, tosave=tosave)
    call ReadCommandLineArgument(arg, "ReweightDecay", success, ReweightDecay)
    call ReadCommandLineArgument(arg, "WidthScheme", success, WidthScheme, tosave=tosave)
    call ReadCommandLineArgument(arg, "WidthSchemeIn", success, WidthSchemeIn)

    call ReadCommandLineArgument(arg, "ReadPmHstar", success, ReadPMZZ)
    call ReadCommandLineArgument(arg, "PmHstarFile", success, PMZZfile)
    call ReadCommandLineArgument(arg, "PrintPmHstar", success, PrintPMZZ, multiplyreal=GeV, success2=DoPrintPMZZ) !undocumented, for internal testing.  PrintPMZZ is a complex(8), the real and imaginary parts are the minimum and maximum values to print
    call ReadCommandLineArgument(arg, "PrintPmHstarIntervals", success, PrintPMZZIntervals)                      !undocumented, for internal testing
    call ReadCommandLineArgument(arg, "PmHstarEvals", success, PMZZEvals)

    !same thing again for compatibility
    call ReadCommandLineArgument(arg, "ReadPMZZ", success, ReadPMZZ)   !undocumented, for compatibility
    call ReadCommandLineArgument(arg, "PMZZFile", success, PMZZfile)   !undocumented, for compatibility
    call ReadCommandLineArgument(arg, "PrintPMZZ", success, PrintPMZZ, multiplyreal=GeV, success2=DoPrintPMZZ)   !undocumented, for compatibility
    call ReadCommandLineArgument(arg, "PrintPMZZIntervals", success, PrintPMZZIntervals)                        !undocumented, for compatibility
    call ReadCommandLineArgument(arg, "PMZZEvals", success, PMZZEvals)   !undocumented, for compatibility

    call ReadCommandLineArgument(arg, "OffshellX", success, OffShellReson, tosave=tosave)
    call ReadCommandLineArgument(arg, "FilterNLept", success, RequestNLeptons)
    call ReadCommandLineArgument(arg, "FilterOSPairs", success, RequestOS)
    call ReadCommandLineArgument(arg, "FilterOSSFPairs", success, RequestOSSF)
    call ReadCommandLineArgument(arg, "CountTauAsAny", success, CountTauAsAny)
    call ReadCommandLineArgument(arg, "FilterNJets", success, RequestNJets)
    call ReadCommandLineArgument(arg, "Unweighted", success, Unweighted)
    call ReadCommandLineArgument(arg, "Interf", success, includeInterference, success2=interfSet, tosave=tosave)
    call ReadCommandLineArgument(arg, "ReweightInterf", success, reweightInterference, success2=interfSet, tosave=tosave)
    call ReadCommandLineArgument(arg, "ReadLHE", success, LHEProdFile, success2=ReadLHEFile)
    call ReadCommandLineArgument(arg, "ConvertLHE", success, LHEProdFile, success2=ConvertLHEFile)
    call ReadCommandLineArgument(arg, "ReadCSmax", success, ReadCSmax)
    call ReadCommandLineArgument(arg, "GenEvents", success, GenerateEvents, SetLastArgument)   !undocumented, not sure what this does
    if( SetLastArgument ) Unweighted = .false.
    call ReadCommandLineArgument(arg, "WriteFailedEvents", success, WriteFailedEvents)
    call ReadCommandLineArgument(arg, "Seed", success, UserSeed)
    call ReadCommandLineArgument(arg, "WriteGit", success, writegit) !undocumented, for testing purposes
    call ReadCommandLineArgument(arg, "RandomizeVVdecays", success, RandomizeVVdecays, tosave=tosave)
    call ReadCommandLineArgument(arg, "ChannelRatio", success, channels_ratio_fix, success2=fix_channels_ratio, tosave=tosave)
    call ReadCommandLineArgument(arg, "UnformattedRead", success, UseUnformattedRead)
    call ReadCommandLineArgument(arg, "WriteWeightedLHE", success, WriteWeightedLHE)
    call ReadCommandLineArgument(arg, "VBFoffsh_run", success, VBFoffsh_run, success2=setVBFoffsh_run)

    !anomalous couplings
    !If any anomalous couplings are set, the default ones have to be set explicitly to keep them on or turn them off
    !e.g. just setting ghz4=0.2982,0 is ambiguous if you mean to leave g1 on (so fa3=0.5)
    !                                                      or to turn it off (fa3=1 with a weird prefactor)
    !spin 0 gg couplings
    call ReadCommandLineArgument(arg, "ghg2", success, ghg2, success2=SetAnomalousSpin0gg, success3=Setghg2, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghg3", success, ghg3, success2=SetAnomalousSpin0gg, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghg4", success, ghg4, success2=SetAnomalousSpin0gg, tosave=tosave)

    !spin 0 ZZ couplings
    call ReadCommandLineArgument(arg, "ghz1", success, ghz1, success2=SetAnomalousSpin0VV, success3=Setghz1, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2", success, ghz2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3", success, ghz3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4", success, ghz4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    !spin 0 Zgamma couplings
    call ReadCommandLineArgument(arg, "ghzgs2", success, ghzgs2, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzgs3", success, ghzgs3, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzgs4", success, ghzgs4, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)

    !spin 0 gammagamma couplings
    call ReadCommandLineArgument(arg, "ghgsgs2", success, ghgsgs2, success2=SetAnomalousSpin0VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghgsgs3", success, ghgsgs3, success2=SetAnomalousSpin0VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghgsgs4", success, ghgsgs4, success2=SetAnomalousSpin0VV, success3=Setgammagammacoupling, tosave=tosave)

    !spin 0 ZZ momentum dependent couplings
    call ReadCommandLineArgument(arg, "ghz1_prime", success, ghz1_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz1_prime2", success, ghz1_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz1_prime3", success, ghz1_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz1_prime4", success, ghz1_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz1_prime5", success, ghz1_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz1_prime6", success, ghz1_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz1_prime7", success, ghz1_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghz2_prime", success, ghz2_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2_prime2", success, ghz2_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2_prime3", success, ghz2_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2_prime4", success, ghz2_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2_prime5", success, ghz2_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2_prime6", success, ghz2_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz2_prime7", success, ghz2_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghz3_prime", success, ghz3_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3_prime2", success, ghz3_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3_prime3", success, ghz3_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3_prime4", success, ghz3_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3_prime5", success, ghz3_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3_prime6", success, ghz3_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz3_prime7", success, ghz3_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghz4_prime", success, ghz4_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4_prime2", success, ghz4_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4_prime3", success, ghz4_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4_prime4", success, ghz4_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4_prime5", success, ghz4_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4_prime6", success, ghz4_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghz4_prime7", success, ghz4_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    !spin 0 Zgamma momentum dependent coupling
    call ReadCommandLineArgument(arg, "ghzgs1_prime2", success, ghzgs1_prime2, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)

    ! Sign of q1,2,12**2 for the Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
    call ReadCommandLineArgument(arg, "cz_q1sq", success, cz_q1sq, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "cz_q2sq", success, cz_q2sq, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "cz_q12sq", success, cz_q12sq, success2=SetAnomalousSpin0VV, tosave=tosave)
    ! Lambda's for q1,2,12**2 for the Lambda's
    call ReadCommandLineArgument(arg, "Lambda_z11", success, Lambda_z11, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z21", success, Lambda_z21, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z31", success, Lambda_z31, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z41", success, Lambda_z41, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z12", success, Lambda_z12, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z22", success, Lambda_z22, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z32", success, Lambda_z32, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z42", success, Lambda_z42, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z10", success, Lambda_z10, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z20", success, Lambda_z20, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z30", success, Lambda_z30, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_z40", success, Lambda_z40, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)

    !spin 0 WW couplings
    call ReadCommandLineArgument(arg, "ghw1", success, ghw1, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2", success, ghw2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3", success, ghw3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4", success, ghw4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghw1_prime", success, ghw1_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw1_prime2", success, ghw1_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw1_prime3", success, ghw1_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw1_prime4", success, ghw1_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw1_prime5", success, ghw1_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw1_prime6", success, ghw1_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw1_prime7", success, ghw1_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghw2_prime", success, ghw2_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2_prime2", success, ghw2_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2_prime3", success, ghw2_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2_prime4", success, ghw2_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2_prime5", success, ghw2_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2_prime6", success, ghw2_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw2_prime7", success, ghw2_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghw3_prime", success, ghw3_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3_prime2", success, ghw3_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3_prime3", success, ghw3_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3_prime4", success, ghw3_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3_prime5", success, ghw3_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3_prime6", success, ghw3_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw3_prime7", success, ghw3_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghw4_prime", success, ghw4_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4_prime2", success, ghw4_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4_prime3", success, ghw4_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4_prime4", success, ghw4_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4_prime5", success, ghw4_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4_prime6", success, ghw4_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghw4_prime7", success, ghw4_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    ! Sign of q1,2,12**2 for the Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
    call ReadCommandLineArgument(arg, "cw_q1sq", success, cw_q1sq, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "cw_q2sq", success, cw_q2sq, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "cw_q12sq", success, cw_q12sq, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    ! Lambda's for q1,2,12**2 for the Lambda's
    call ReadCommandLineArgument(arg, "Lambda_w11", success, Lambda_w11, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w21", success, Lambda_w21, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w31", success, Lambda_w31, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w41", success, Lambda_w41, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w12", success, Lambda_w12, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w22", success, Lambda_w22, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w32", success, Lambda_w32, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w42", success, Lambda_w42, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w10", success, Lambda_w10, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w20", success, Lambda_w20, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w30", success, Lambda_w30, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda_w40", success, Lambda_w40, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)

    !spin 1
    call ReadCommandLineArgument(arg, "zprime_qq_left", success, zprime_qq_left, success2=SetAnomalousSpin1qq, success3=Setspin1qqleft, tosave=tosave)
    call ReadCommandLineArgument(arg, "zprime_qq_right", success, zprime_qq_right, success2=SetAnomalousSpin1qq, success3=Setspin1qqright, tosave=tosave)
    call ReadCommandLineArgument(arg, "zprime_zz_1", success, zprime_zz_1, success2=SetSpin1VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "zprime_zz_2", success, zprime_zz_2, success2=SetSpin1VV, success3=SetZZcoupling, tosave=tosave)

    !spin 2
    call ReadCommandLineArgument(arg, "a1", success, a1, success2=SetAnomalousSpin2gg, tosave=tosave)
    call ReadCommandLineArgument(arg, "a2", success, a2, success2=SetAnomalousSpin2gg, tosave=tosave)
    call ReadCommandLineArgument(arg, "a3", success, a3, success2=SetAnomalousSpin2gg, tosave=tosave)
    call ReadCommandLineArgument(arg, "a4", success, a4, success2=SetAnomalousSpin2gg, tosave=tosave)
    call ReadCommandLineArgument(arg, "a5", success, a5, success2=SetAnomalousSpin2gg, tosave=tosave)
    call ReadCommandLineArgument(arg, "graviton_qq_left", success, graviton_qq_left, success2=SetAnomalousSpin2qq, success3=Setspin2qqleft, tosave=tosave)
    call ReadCommandLineArgument(arg, "graviton_qq_right", success, graviton_qq_right, success2=SetAnomalousSpin2qq, success3=Setspin2qqright, tosave=tosave)
    call ReadCommandLineArgument(arg, "b1", success, b1, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b2", success, b2, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b3", success, b3, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b4", success, b4, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b5", success, b5, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b6", success, b6, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b7", success, b7, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b8", success, b8, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b9", success, b9, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "b10", success, b10, success2=SetSpin2VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "bzgs1", success, bzgs1, success2=SetSpin2VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzgs2", success, bzgs2, success2=SetSpin2VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzgs3", success, bzgs3, success2=SetSpin2VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzgs4", success, bzgs4, success2=SetSpin2VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzgs8", success, bzgs8, success2=SetSpin2VV, success3=SetZgammacoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "bgsgs1", success, bgsgs1, success2=SetSpin2VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bgsgs2", success, bgsgs2, success2=SetSpin2VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bgsgs3", success, bgsgs3, success2=SetSpin2VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bgsgs4", success, bgsgs4, success2=SetSpin2VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bgsgs8", success, bgsgs8, success2=SetSpin2VV, success3=Setgammagammacoupling, tosave=tosave)

    !Hff couplings
    call ReadCommandLineArgument(arg, "kappa", success, kappa, success2=SetAnomalousHff, success3=Setkappa, tosave=tosave)
    call ReadCommandLineArgument(arg, "kappa_tilde", success, kappa_tilde, success2=SetAnomalousHff, tosave=tosave)

!   similar as above for the 2nd resonance in offshell VBF
    !spin 0 ZZ couplings
    call ReadCommandLineArgument(arg, "gh2z1", success, gh2z1, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2", success, gh2z2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3", success, gh2z3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4", success, gh2z4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    !spin 0 Zgamma couplings
    call ReadCommandLineArgument(arg, "gh2zgs2", success, gh2zgs2, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2zgs3", success, gh2zgs3, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2zgs4", success, gh2zgs4, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)

    !spin 0 gammagamma couplings
    call ReadCommandLineArgument(arg, "gh2gsgs2", success, gh2gsgs2, success2=SetAnomalousSpin0VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2gsgs3", success, gh2gsgs3, success2=SetAnomalousSpin0VV, success3=Setgammagammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2gsgs4", success, gh2gsgs4, success2=SetAnomalousSpin0VV, success3=Setgammagammacoupling, tosave=tosave)

    !spin 0 ZZ momentum dependent couplings
    call ReadCommandLineArgument(arg, "gh2z1_prime", success, gh2z1_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z1_prime2", success, gh2z1_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z1_prime3", success, gh2z1_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z1_prime4", success, gh2z1_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z1_prime5", success, gh2z1_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z1_prime6", success, gh2z1_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z1_prime7", success, gh2z1_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2z2_prime", success, gh2z2_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2_prime2", success, gh2z2_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2_prime3", success, gh2z2_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2_prime4", success, gh2z2_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2_prime5", success, gh2z2_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2_prime6", success, gh2z2_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z2_prime7", success, gh2z2_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2z3_prime", success, gh2z3_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3_prime2", success, gh2z3_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3_prime3", success, gh2z3_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3_prime4", success, gh2z3_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3_prime5", success, gh2z3_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3_prime6", success, gh2z3_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z3_prime7", success, gh2z3_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2z4_prime", success, gh2z4_prime, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4_prime2", success, gh2z4_prime2, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4_prime3", success, gh2z4_prime3, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4_prime4", success, gh2z4_prime4, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4_prime5", success, gh2z4_prime5, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4_prime6", success, gh2z4_prime6, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2z4_prime7", success, gh2z4_prime7, success2=SetAnomalousSpin0VV, success3=SetZZcoupling, tosave=tosave)

    !spin 0 Zgamma momentum dependent coupling
    call ReadCommandLineArgument(arg, "gh2zgs1_prime2", success, gh2zgs1_prime2, success2=SetAnomalousSpin0VV, success3=SetZgammacoupling, tosave=tosave)

    ! Sign of q1,2,12**2 for the Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
    call ReadCommandLineArgument(arg, "c2z_q1sq", success, c2z_q1sq, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "c2z_q2sq", success, c2z_q2sq, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "c2z_q12sq", success, c2z_q12sq, success2=SetAnomalousSpin0VV, tosave=tosave)
    ! Lambda's for q1,2,12**2 for the Lambda's
    call ReadCommandLineArgument(arg, "Lambda2_z11", success, Lambda2_z11, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z21", success, Lambda2_z21, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z31", success, Lambda2_z31, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z41", success, Lambda2_z41, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z12", success, Lambda2_z12, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z22", success, Lambda2_z22, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z32", success, Lambda2_z32, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z42", success, Lambda2_z42, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z10", success, Lambda2_z10, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z20", success, Lambda2_z20, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z30", success, Lambda2_z30, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_z40", success, Lambda2_z40, multiply=GeV, success2=SetAnomalousSpin0VV, tosave=tosave)

    !spin 0 WW couplings
    call ReadCommandLineArgument(arg, "gh2w1", success, gh2w1, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2", success, gh2w2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3", success, gh2w3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4", success, gh2w4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2w1_prime", success, gh2w1_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w1_prime2", success, gh2w1_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w1_prime3", success, gh2w1_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w1_prime4", success, gh2w1_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w1_prime5", success, gh2w1_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w1_prime6", success, gh2w1_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w1_prime7", success, gh2w1_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2w2_prime", success, gh2w2_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2_prime2", success, gh2w2_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2_prime3", success, gh2w2_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2_prime4", success, gh2w2_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2_prime5", success, gh2w2_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2_prime6", success, gh2w2_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w2_prime7", success, gh2w2_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2w3_prime", success, gh2w3_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3_prime2", success, gh2w3_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3_prime3", success, gh2w3_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3_prime4", success, gh2w3_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3_prime5", success, gh2w3_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3_prime6", success, gh2w3_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w3_prime7", success, gh2w3_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "gh2w4_prime", success, gh2w4_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4_prime2", success, gh2w4_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4_prime3", success, gh2w4_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4_prime4", success, gh2w4_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4_prime5", success, gh2w4_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4_prime6", success, gh2w4_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "gh2w4_prime7", success, gh2w4_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=SetWWcoupling, tosave=tosave)

    ! Sign of q1,2,12**2 for the Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
    call ReadCommandLineArgument(arg, "c2w_q1sq", success, c2w_q1sq, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "c2w_q2sq", success, c2w_q2sq, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "c2w_q12sq", success, c2w_q12sq, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    ! Lambda's for q1,2,12**2 for the Lambda's
    call ReadCommandLineArgument(arg, "Lambda2_w11", success, Lambda2_w11, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w21", success, Lambda2_w21, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w31", success, Lambda2_w31, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w41", success, Lambda2_w41, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w12", success, Lambda2_w12, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w22", success, Lambda2_w22, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w32", success, Lambda2_w32, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w42", success, Lambda2_w42, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w10", success, Lambda2_w10, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w20", success, Lambda2_w20, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w30", success, Lambda2_w30, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)
    call ReadCommandLineArgument(arg, "Lambda2_w40", success, Lambda2_w40, multiply=GeV, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, tosave=tosave)


    !contact interactions
    !spin 0 ZZp couplings
    call ReadCommandLineArgument(arg, "ghzzp1", success, ghzzp1, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2", success, ghzzp2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3", success, ghzzp3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4", success, ghzzp4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzzp1_prime", success, ghzzp1_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp1_prime2", success, ghzzp1_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp1_prime3", success, ghzzp1_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp1_prime4", success, ghzzp1_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp1_prime5", success, ghzzp1_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp1_prime6", success, ghzzp1_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp1_prime7", success, ghzzp1_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzzp2_prime", success, ghzzp2_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2_prime2", success, ghzzp2_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2_prime3", success, ghzzp2_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2_prime4", success, ghzzp2_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2_prime5", success, ghzzp2_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2_prime6", success, ghzzp2_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp2_prime7", success, ghzzp2_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzzp3_prime", success, ghzzp3_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3_prime2", success, ghzzp3_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3_prime3", success, ghzzp3_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3_prime4", success, ghzzp3_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3_prime5", success, ghzzp3_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3_prime6", success, ghzzp3_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp3_prime7", success, ghzzp3_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzzp4_prime", success, ghzzp4_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4_prime2", success, ghzzp4_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4_prime3", success, ghzzp4_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4_prime4", success, ghzzp4_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4_prime5", success, ghzzp4_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4_prime6", success, ghzzp4_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzzp4_prime7", success, ghzzp4_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)

    !spin 0 Zpgamma couplings
    call ReadCommandLineArgument(arg, "ghzpgs1_prime2", success, ghzpgs1_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpgs2", success, ghzpgs2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpgs3", success, ghzpgs3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpgs4", success, ghzpgs4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)

    !spin 0 ZpZp couplings
    call ReadCommandLineArgument(arg, "ghzpzp1", success, ghzpzp1, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2", success, ghzpzp2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3", success, ghzpzp3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4", success, ghzpzp4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzpzp1_prime", success, ghzpzp1_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp1_prime2", success, ghzpzp1_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp1_prime3", success, ghzpzp1_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp1_prime4", success, ghzpzp1_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp1_prime5", success, ghzpzp1_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp1_prime6", success, ghzpzp1_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp1_prime7", success, ghzpzp1_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzpzp2_prime", success, ghzpzp2_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2_prime2", success, ghzpzp2_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2_prime3", success, ghzpzp2_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2_prime4", success, ghzpzp2_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2_prime5", success, ghzpzp2_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2_prime6", success, ghzpzp2_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp2_prime7", success, ghzpzp2_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzpzp3_prime", success, ghzpzp3_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3_prime2", success, ghzpzp3_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3_prime3", success, ghzpzp3_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3_prime4", success, ghzpzp3_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3_prime5", success, ghzpzp3_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3_prime6", success, ghzpzp3_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp3_prime7", success, ghzpzp3_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghzpzp4_prime", success, ghzpzp4_prime, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4_prime2", success, ghzpzp4_prime2, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4_prime3", success, ghzpzp4_prime3, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4_prime4", success, ghzpzp4_prime4, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4_prime5", success, ghzpzp4_prime5, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4_prime6", success, ghzpzp4_prime6, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghzpzp4_prime7", success, ghzpzp4_prime7, success2=SetAnomalousSpin0VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "bzzp1", success, bzzp1, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp2", success, bzzp2, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp3", success, bzzp3, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp4", success, bzzp4, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp5", success, bzzp5, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp6", success, bzzp6, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp7", success, bzzp7, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp8", success, bzzp8, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp9", success, bzzp9, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzzp10", success, bzzp10, success2=SetSpin2VV, success3=includeVprime, success4=SetZZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "bzpzp1", success, bzpzp1, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp2", success, bzpzp2, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp3", success, bzpzp3, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp4", success, bzpzp4, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp5", success, bzpzp5, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp6", success, bzpzp6, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp7", success, bzpzp7, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp8", success, bzpzp8, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp9", success, bzpzp9, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpzp10", success, bzpzp10, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimeZprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "bzpgs1", success, bzpgs1, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpgs2", success, bzpgs2, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpgs3", success, bzpgs3, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpgs4", success, bzpgs4, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "bzpgs8", success, bzpgs8, success2=SetSpin2VV, success3=includeVprime, success4=SetZprimegammacoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ezp_El_left", success, ezp_El_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_El_right", success, ezp_El_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Mu_left", success, ezp_Mu_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Mu_right", success, ezp_Mu_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Ta_left", success, ezp_Ta_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Ta_right", success, ezp_Ta_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Up_left", success, ezp_Up_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Up_right", success, ezp_Up_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Chm_left", success, ezp_Chm_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Chm_right", success, ezp_Chm_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Top_left", success, ezp_Top_left, success2=SetZprimeff, tosave=tosave)    !undocumented because it's useless
    call ReadCommandLineArgument(arg, "ezp_Top_right", success, ezp_Top_right, success2=SetZprimeff, tosave=tosave)  !undocumented because it's useless
    call ReadCommandLineArgument(arg, "ezp_Dn_left", success, ezp_Dn_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Dn_right", success, ezp_Dn_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Str_left", success, ezp_Str_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Str_right", success, ezp_Str_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Bot_left", success, ezp_Bot_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_Bot_right", success, ezp_Bot_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_NuE_left", success, ezp_NuE_left, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ezp_NuE_right", success, ezp_NuE_right, success2=SetZprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "MZprime", success, M_Zprime, multiply=GeV, success2=SetMZprime, tosave=tosave)
    call ReadCommandLineArgument(arg, "GaZprime", success, Ga_Zprime, multiply=GeV, success2=SetGaZprime, tosave=tosave)

    !spin 0 WWp couplings
    call ReadCommandLineArgument(arg, "ghwwp1", success, ghwwp1, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2", success, ghwwp2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3", success, ghwwp3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4", success, ghwwp4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwwp1_prime", success, ghwwp1_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp1_prime2", success, ghwwp1_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp1_prime3", success, ghwwp1_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp1_prime4", success, ghwwp1_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp1_prime5", success, ghwwp1_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp1_prime6", success, ghwwp1_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp1_prime7", success, ghwwp1_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwwp2_prime", success, ghwwp2_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2_prime2", success, ghwwp2_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2_prime3", success, ghwwp2_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2_prime4", success, ghwwp2_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2_prime5", success, ghwwp2_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2_prime6", success, ghwwp2_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp2_prime7", success, ghwwp2_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwwp3_prime", success, ghwwp3_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3_prime2", success, ghwwp3_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3_prime3", success, ghwwp3_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3_prime4", success, ghwwp3_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3_prime5", success, ghwwp3_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3_prime6", success, ghwwp3_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp3_prime7", success, ghwwp3_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwwp4_prime", success, ghwwp4_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4_prime2", success, ghwwp4_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4_prime3", success, ghwwp4_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4_prime4", success, ghwwp4_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4_prime5", success, ghwwp4_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4_prime6", success, ghwwp4_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwwp4_prime7", success, ghwwp4_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWWprimecoupling, tosave=tosave)

    !spin 0 WpWp couplings
    call ReadCommandLineArgument(arg, "ghwpwp1", success, ghwpwp1, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2", success, ghwpwp2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3", success, ghwpwp3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4", success, ghwpwp4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwpwp1_prime", success, ghwpwp1_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp1_prime2", success, ghwpwp1_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp1_prime3", success, ghwpwp1_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp1_prime4", success, ghwpwp1_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp1_prime5", success, ghwpwp1_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp1_prime6", success, ghwpwp1_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp1_prime7", success, ghwpwp1_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwpwp2_prime", success, ghwpwp2_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2_prime2", success, ghwpwp2_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2_prime3", success, ghwpwp2_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2_prime4", success, ghwpwp2_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2_prime5", success, ghwpwp2_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2_prime6", success, ghwpwp2_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp2_prime7", success, ghwpwp2_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwpwp3_prime", success, ghwpwp3_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3_prime2", success, ghwpwp3_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3_prime3", success, ghwpwp3_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3_prime4", success, ghwpwp3_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3_prime5", success, ghwpwp3_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3_prime6", success, ghwpwp3_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp3_prime7", success, ghwpwp3_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ghwpwp4_prime", success, ghwpwp4_prime, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4_prime2", success, ghwpwp4_prime2, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4_prime3", success, ghwpwp4_prime3, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4_prime4", success, ghwpwp4_prime4, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4_prime5", success, ghwpwp4_prime5, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4_prime6", success, ghwpwp4_prime6, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)
    call ReadCommandLineArgument(arg, "ghwpwp4_prime7", success, ghwpwp4_prime7, success2=distinguish_HWWcouplings, success3=SetAnomalousSpin0VV, success4=includeVprime, success5=SetWprimeWprimecoupling, tosave=tosave)

    call ReadCommandLineArgument(arg, "ewp_El_left", success, ewp_El_left, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_El_right", success, ewp_El_right, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Mu_left", success, ewp_Mu_left, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Mu_right", success, ewp_Mu_right, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Ta_left", success, ewp_Ta_left, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Ta_right", success, ewp_Ta_right, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Up_left", success, ewp_Up_left, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Up_right", success, ewp_Up_right, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Chm_left", success, ewp_Chm_left, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Chm_right", success, ewp_Chm_right, success2=SetWprimeff, tosave=tosave)
    call ReadCommandLineArgument(arg, "ewp_Top_left", success, ewp_Top_left, success2=SetWprimeff, tosave=tosave)    !undocumented because it's useless (until contact terms are included in tH)
    call ReadCommandLineArgument(arg, "ewp_Top_right", success, ewp_Top_right, success2=SetWprimeff, tosave=tosave)  !undocumented because it's useless (until contact terms are included in tH)
    call ReadCommandLineArgument(arg, "MWprime", success, M_Wprime, multiply=GeV, success2=SetMWprime, tosave=tosave)
    call ReadCommandLineArgument(arg, "GaWprime", success, Ga_Wprime, multiply=GeV, success2=SetGaWprime, tosave=tosave)


    ! CKM elements
    call ReadCommandLineArgument(arg, "Vud", success, VCKM_ud, success2=SetCKM, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vus", success, VCKM_us, success2=SetCKM, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vcd", success, VCKM_cd, success2=SetCKM, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vcs", success, VCKM_cs, success2=SetCKM, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vts", success, VCKM_ts, success2=SetCKM, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vtb", success, VCKM_tb, success2=SetCKM, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vub", success, VCKM_ub, success2=SetCKM, success3=SetCKMub, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vcb", success, VCKM_cb, success2=SetCKM, success3=SetCKMcb, tosave=tosave)
    call ReadCommandLineArgument(arg, "Vtd", success, VCKM_td, success2=SetCKM, success3=SetCKMtd, tosave=tosave)


    !cuts
    call ReadCommandLineArgument(arg, "pTjetcut", success, pTjetcut, multiply=GeV, success2=SetpTjetcut, tosave=tosave)
    call ReadCommandLineArgument(arg, "etajetcut", success, etajetcut, success2=Setetajetcut, tosave=tosave)
    call ReadCommandLineArgument(arg, "detajetcut", success, detajetcut, success2=Setdetajetcut, tosave=tosave)
    call ReadCommandLineArgument(arg, "deltaRcut", success, Rjet, success2=SetdeltaRcut, tosave=tosave)
    call ReadCommandLineArgument(arg, "mJJcut", success, mJJcut, multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "m4l_min", success, m4l_minmax(1), multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "m4l_max", success, m4l_minmax(2), multiply=GeV, tosave=tosave)
    call ReadCommandLineArgument(arg, "m2l_min", success, m2l_minmax(1), multiply=GeV, success2=Setm2l_min, tosave=tosave)
    call ReadCommandLineArgument(arg, "m2l_max", success, m2l_minmax(2), multiply=GeV, success2=Setm2l_max, tosave=tosave)
    call ReadCommandLineArgument(arg, "mVH_min", success, mVH_minmax(1), multiply=GeV, success2=SetmVH_min, tosave=tosave)
    call ReadCommandLineArgument(arg, "mVH_max", success, mVH_minmax(2), multiply=GeV, success2=SetmVH_max, tosave=tosave)
    call ReadCommandLineArgument(arg, "MPhotonCutoff", success, MPhotonCutoff, multiply=GeV, success2=SetMPhotonCutoff, tosave=tosave)
    call ReadCommandLineArgument(arg, "pTlepcut", success, pTlepcut, multiply=GeV, success2=SetpTlepcut, tosave=tosave)
    call ReadCommandLineArgument(arg, "etalepcut", success, etalepcut, success2=Setetalepcut, tosave=tosave)
    call ReadCommandLineArgument(arg, "JetsOppositeEta", success, JetsOppositeEta, tosave=tosave)
    call ReadCommandLineArgument(arg, "pTHcut", success, pTHcut, success2=SetpTHcut, tosave=tosave)

    if( .not.success ) then
        call Error("Unknown command line argument: " // trim(arg))
    endif
   enddo

   SetHZprime = SetZZprimecoupling .or. SetZprimeZprimecoupling .or. SetZprimegammacoupling
   SetHWprime = SetWWprimecoupling .or. SetWprimeWprimecoupling

    !================================
    !Command line argument processing
    !================================

    !output file

    i = len(trim(DataFile))
    if( i.gt.3 ) then
       if( DataFile(i-3:i).eq.".lhe" ) then
          DataFile = DataFile(1:i-4)
       endif
    endif
    call system('mkdir -p ./data')! -p is suppressing error messages if directory already exists

    if( VBFoffsh_run.eq.1 ) DataFile = trim(DataFile)//"_1"
    if( VBFoffsh_run.eq.2 ) DataFile = trim(DataFile)//"_2"
    if( VBFoffsh_run.eq.3 ) DataFile = trim(DataFile)//"_3"
    if( VBFoffsh_run.eq.4 ) DataFile = trim(DataFile)//"_4"
    if( VBFoffsh_run.eq.5 ) DataFile = trim(DataFile)//"_5"
    if( VBFoffsh_run.eq.6 ) DataFile = trim(DataFile)//"_6"
    if( VBFoffsh_run.eq.7 ) DataFile = trim(DataFile)//"_7"
    if( VBFoffsh_run.eq.8 ) DataFile = trim(DataFile)//"_8"

    if( SetCSmaxFile ) then
      i = len(trim(CSmaxFile))
      if( CSmaxFile(i-3:i).eq.".lhe" ) then
          CSmaxFile = CSmaxFile(1:i-4)
      endif
      if( VBFoffsh_run.eq.1 ) CSmaxFile = trim(CSmaxFile)//"_1"
      if( VBFoffsh_run.eq.2 ) CSmaxFile = trim(CSmaxFile)//"_2"
      if( VBFoffsh_run.eq.3 ) CSmaxFile = trim(CSmaxFile)//"_3"
      if( VBFoffsh_run.eq.4 ) CSmaxFile = trim(CSmaxFile)//"_4"
      if( VBFoffsh_run.eq.5 ) CSmaxFile = trim(CSmaxFile)//"_5"
      if( VBFoffsh_run.eq.6 ) CSmaxFile = trim(CSmaxFile)//"_6"
      if( VBFoffsh_run.eq.7 ) CSmaxFile = trim(CSmaxFile)//"_7"
      if( VBFoffsh_run.eq.8 ) CSmaxFile = trim(CSmaxFile)//"_8"
    else
      CSmaxFile = DataFile
    endif

    if (ReadCSmax) then
      call oldsavevalues%ReadFromFile(trim(CSmaxFile)//"_commandlineinfo.txt")
      call CompareSaveValues(tosave, oldsavevalues)
    else
      call tosave%WriteToFile(trim(CSmaxFile)//"_commandlineinfo.txt")
    endif

    !VBF offshell
    if( .not. (Process.ge.66 .and. Process.le.69) .and. SetVBFoffsh_run ) then
      call Error("VBFoffsh_run is only for VBF offshell")
    endif

    !PChannel
    if (Process.eq.0) PChannel = 0   !only gluons
!Yaofu Zhou ggZH
    !if (Process.eq.1 .or. Process.eq.51 .or. Process.eq.60 .or. Process.eq.66 .or. Process.eq.67 .or. Process.eq.68) PChannel = 1   !only quarks
    if (Process.eq.1 .or. Process.eq.60 .or. Process.eq.66 .or. Process.eq.67 .or. Process.eq.68) PChannel = 1   !only quarks

    !LHAPDF

#if useLHAPDF==1
    call GET_ENVIRONMENT_VARIABLE("LHAPDF_DATA_PATH", LHAPDF_DATA_PATH)
    if( LHAPDFString.eq."" ) then
       if (.not. (ReadLHEFile .or. ConvertLHEFile)) then
         print *, "Need to specify pdf file name in command line argument LHAPDF"
         stop 1
       endif
    endif
#endif

    !Collider and energy
     if (.not.SetColliderEnergy) then
      IF( COLLIDER.EQ.1) THEN
        Collider_Energy  = LHC_Energy
      ELSEIF( COLLIDER.EQ.2 ) THEN
        Collider_Energy  = TEV_Energy
      ELSEIF( COLLIDER.EQ.0 ) THEN
        Collider_Energy  = ILC_Energy
      ENDIF
    endif

    !Renormalization/factorization schemes

    if( &
      (abs(FacScheme) .ge. nRenFacSchemes) .or. (abs(RenScheme) .ge. nRenFacSchemes)     .or. &
      (FacScheme .eq. -kRenFacScheme_maxpTj) .or. (RenScheme .eq. -kRenFacScheme_maxpTj) .or. &
      (FacScheme .eq. -kRenFacScheme_minpTj) .or. (RenScheme .eq. -kRenFacScheme_minpTj) .or. &
      (MuFacMultiplier.le.0d0) .or. (MuRenMultiplier.le.0d0)                                  &
    ) call Error("The renormalization or factorization scheme is invalid, or the scale multiplier to either is not positive.")
    if( SetFacScheme.neqv.SetMuFacMultiplier ) call Error("If you want to set the factorization scale, please set both the scheme (FacScheme) and the multiplier (MuFacMultiplier)")
    if( SetRenScheme.neqv.SetMuRenMultiplier ) call Error("If you want to set the renormalization scale, please set both the scheme (RenScheme) and the multiplier (MuRenMultiplier)")

    !DecayModes

    if( ConvertLHEFile ) then
       DecayMode2 = DecayMode1
    endif

    if( Process.eq.50 ) then
        DecayMode2=DecayMode1
        if( Collider.eq.2 ) call Error("Collider 2 not available for VHiggs")
        if( (IsAWDecay(DecayMode1) ) .and. (Collider.ne.1) ) call Error("WHiggs with Collider 1 only")
    endif

    if( Process.eq.51 ) then
        DecayMode2=DecayMode1
        if( Collider.eq.2 ) call Error("Collider 2 not available for VH")
        if( (IsAWDecay(DecayMode1) ) .and. (Collider.ne.1) ) call Error("WH with Collider 1 only")
    endif

    if( Process.ge.110 .and. Process.le.114 ) DecayMode2 = DecayMode1

    if( (TopDecays.ne.0 .and. TopDecays.ne.1) .and. (Process.eq.80 .or. (Process.ge.110 .and. Process.le.114)) ) call Error("Specify TopDK=0,1")
    if( (TopDecays.eq.1) .and. .not. IsAWDecay(DecayMode1) ) call Error("Invalid DecayMode1 for top decays")
    if( (TopDecays.eq.1) .and. .not. IsAWDecay(DecayMode2) ) call Error("Invalid DecayMode2 for top decays")

    if( (TauDecays.eq.1) .and. .not. IsAWDecay(DecayMode1) ) call Error("Invalid DecayMode1 for tau decays")
    if( (TauDecays.eq.1) .and. .not. IsAWDecay(DecayMode2) ) call Error("Invalid DecayMode2 for tau decays")

    ! Checks on DecayMode1 and 2
    ! First check!
    if( IsAPhoton(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed. Swapping the decay modes."
       DecayMode1 = DecayMode1 + DecayMode2
       DecayMode2 = DecayMode1 - DecayMode2
       DecayMode1 = DecayMode1 - DecayMode2
    endif

    if( IsAZDecay(DecayMode1) .or. (Process.eq.51.and.IsAPhoton(DecayMode1)) ) then
       M_V = M_Z
       Ga_V= Ga_Z
    elseif( IsAWDecay(DecayMode1) ) then
       M_V = M_W
       Ga_V= Ga_W
    elseif( IsAPhoton(DecayMode1) ) then
       M_V = 0d0
       Ga_V= 0d0
    endif

    if( IsAZDecay(DecayMode1) .or. (Process.eq.50.and.IsAPhoton(DecayMode1)) ) then
       M_V = M_Z
       Ga_V= Ga_Z
       M_Vprime = M_Zprime
       Ga_Vprime = Ga_Zprime
    elseif( IsAWDecay(DecayMode1) ) then
       M_V = M_W
       Ga_V= Ga_W
       M_Vprime = M_Wprime
       Ga_Vprime = Ga_Wprime
    elseif( IsAPhoton(DecayMode1) ) then
       M_V = 0d0
       Ga_V= 0d0
       M_Vprime = -1d0
       Ga_Vprime = 0d0
    endif

    M_V_ps = M_V
    Ga_V_ps = Ga_V
    M_Z_ps = M_Z
    Ga_Z_ps = Ga_Z
    M_W_ps = M_W
    Ga_W_ps = Ga_W

    if(      .not.SetZZcoupling .and. .not.SetZgammacoupling .and. .not.SetZZprimecoupling &
       .and. (Process.ne.0 .or. ghz1.eq.0d0) &    !for Process=0 you have to have explicitly turned off the SM coupling
       .and. (SetZprimeZprimecoupling .or. SetZprimegammacoupling) &
       .and. Process.le.2 .and. (SetMZprime.and.IsAZDecay(DecayMode1) .or. SetMWprime.and.IsAWDecay(DecayMode1))) then
       !need more complicated logic here if this is done for VBF
       M_V_ps = M_Vprime
       Ga_V_ps = Ga_Vprime
       M_Z_ps = M_Zprime
       Ga_Z_ps = Ga_Zprime
       M_W_ps = M_Wprime
       Ga_W_ps = Ga_Wprime
    endif

    !ReadLHE and ConvertLHE
    !MUST HAPPEN BEFORE DETERMINING INTERFERENCE
    !so that the mass and width can be read
    !AND BEFORE DEALING WITH WIDTHSCHEMES
    !so that WidthSchemeIn can be read

    if( ReadLHEFile .and. Process.ne.0  ) then
        print *, "ReadLHE option is only allowed for spin-0 resonances"
        stop 1
    endif
    if( ConvertLHEFile .and. Process.ne.0  ) then
        print *, "ConvertLHE option is only allowed for spin-0 resonances"
        stop 1
    endif
    if( ReadLHEFile .and. .not. Unweighted ) then
        print *, "ReadLHE option is only allowed for generating unweighted events"
        stop 1
    endif

    if (ReadLHEFile .or. ConvertLHEFile) then
       call OpenFiles()
       call ReadMassWidth()
    endif

    !interference, photon couplings, ...
    if( (DecayMode1.eq.DecayMode2 .and. IsAZDecay(DecayMode1)) .or.  &
        (DecayMode1.eq.9) .or. (DecayMode2.eq.9)               .or.  &
        (DecayMode1.eq.8  .and. DecayMode2.eq.0)               .or.  &
        (DecayMode1.eq.8  .and. DecayMode2.eq.2)               .or.  &
        (DecayMode1.eq.0  .and. DecayMode2.eq.8)               .or.  &
        (DecayMode1.eq.2  .and. DecayMode2.eq.8)               ) then !  allow interference
            if( .not.interfSet ) then!  set default interference switch
                if (Process.ge.66 .and. Process.le.69) then
                    if (Unweighted) then
                        includeInterference = .false.
                        reweightInterference = .true.
                    else if (m4l_minmax(1) .gt. 2d0*M_Z) then
                        includeInterference = .false.
                        reweightInterference = .false.
                    else
                        includeInterference = .true.
                        reweightInterference = .false.
                    endif
                elseif( M_Reso.gt.2d0*M_Z ) then
                    includeInterference = .false.
                    reweightInterference = .false.
                else
                    includeInterference = .true.
                    reweightInterference = .false.
                endif
            endif
    else
        includeInterference = .false.   ! no interference if decay mode does not allow 4 same flavor leptons
        reweightInterference = .false.
    endif

    if (reweightInterference .and. includeInterference) call Error("Can't set both Interf and ReweightInterf")
    if (reweightInterference .and. .not.unweighted) call Error("Can't reweight interference for weighted events, try setting Interf=1 instead")
    if (reweightInterference .and. ReadLHEFile) call Error("Interference reweighting is not implemented for ReadLHE mode")
    if (reweightInterference .and. .not. (Process.ge.66 .and. Process.le.69)) then
      print *, "Interference reweighting is not implemented for process ", Process
      stop 1
    endif

    !decay mode checks
    if( (IsAZDecay(DecayMode1) .and. IsAZDecay(DecayMode2) .and. (Process.eq.0 .or. Process.eq.2) .and. TauDecays.lt.0) .or. (Process.ge.50 .and. Process.le.51 .and. IsAZDecay(DecayMode1)) .or. Process.eq.60 .or. Process.eq.66 ) then
        includeGammaStar = (SetZgammacoupling .or. Setgammagammacoupling .or. SetZprimegammacoupling)
    elseif( (IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) .and. (Process.eq.0 .or. Process.eq.2) .and. TauDecays.lt.0) .or. (Process.ge.50 .and. Process.le.51 .and. IsAPhoton(DecayMode1)) ) then
        includeGammaStar = Setgammagammacoupling
    else if (Process.eq.67 .or. Process.eq.68 .or. Process.eq.69) then
        includeGammaStar = .true. ! Not really gamma*, but rather gamma* or gluon, set to true to manipulate phasespace generation
    endif

    if( (DecayMode1.ge.12) .or. (DecayMode2.ge.12) .or. (DecayMode1.lt..0) .or. (DecayMode2.lt.0) ) then
       print *," DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif

    if( &
        ( IsAWDecay(DecayMode1).and.(IsAZDecay(DecayMode2).or.IsAPhoton(DecayMode2)) ) .or. &
        ( IsAWDecay(DecayMode2).and.(IsAZDecay(DecayMode1).or.IsAPhoton(DecayMode1)) )     &
        ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif


    if( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) .and. .not.SetZgammacoupling .and. .not.Setgammagammacoupling .and. .not.SetZprimegammacoupling .and. (Process.eq.0 .or. Process.eq.2)) then
        if (Process.eq.0) then
          print *, "To decay the resonance to Z+photon, you need to set one of the HZgamma (ghzgs*) or Hgammagamma (ghgsgs*) couplings."
        else if (Process.eq.2) then
          print *, "To decay the resonance to Z+photon, you need to set one of the GZgamma (bzgs*) or Ggammagamma (bgsgs*) couplings."
        end if
        stop 1
    endif

    if( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) .and. .not.Setgammagammacoupling .and. (Process.eq.0 .or. Process.eq.2)) then
        if (Process.eq.0) then
          print *, "To decay the resonance to photon+photon, you need to set one of the Hgammagamma (ghgsgs*) couplings."
        else if (Process.eq.2) then
          print *, "To decay the resonance to photon+photon, you need to set one of the Ggammagamma (bgsgs*) couplings."
        end if
        stop 1
    endif

    if( Process.ge.50 .and. Process.le.51 .and. IsAPhoton(DecayMode1) .and. .not.SetZgammacoupling .and. .not.Setgammagammacoupling ) then
        print *, "To produce gammaH, you need to set one of the HZgamma (ghzgs*) or Hgammagamma(ghgsgs*) couplings."
        stop 1
    endif

    !cut checks
    if(.not.SetpTjetcut) then
        if(Process.eq.50 .or. Process.eq.51) then
            pTjetcut = 0d0*GeV
        else
            pTjetcut = 15d0*GeV
        endif
    endif
    if(.not.Setetajetcut) then
        if(Process.ge.66 .and. Process.le.69) then
            etajetcut = 4d0
        else
            etajetcut = infinity()
        endif
    endif
    if(.not.Setdetajetcut) then
        if(Process.ge.66 .and. Process.le.69) then
            detajetcut = 2d0
        else
            detajetcut = 0d0
        endif
    endif
    if(.not.SetdeltaRcut) then
        if(Process.eq.50 .or. Process.eq.51) then
            Rjet = 0d0
        else
            Rjet = 0.3d0
        endif
    endif
    if(.not.SetpTlepcut) then
        if(Process.ge.66 .and. Process.le.69) then
            pTlepcut = 3d0*GeV
        else
            pTlepcut = 0d0*GeV
        endif
    endif
    if(.not.Setetalepcut) then
        if(Process.ge.66 .and. Process.le.69) then
            etalepcut = 2.7d0
        else
            etalepcut = infinity()
        endif
    endif
    if(.not.SetMPhotonCutoff) then
        if(includeGammaStar) then
            MPhotonCutoff = 4d0*GeV
        elseif(Process.ge.66.and.Process.le.69) then
            MPhotonCutoff = 2.5d0*GeV
        else
            MPhotonCutoff = 0d0
        endif
    endif
    if(.not.Setm2l_min)then
      if(Process.eq.51)then
         m2l_minmax(1) = 0d0 * GeV
      else
         m2l_minmax(1) = 0d0 * GeV
      endif
    endif
    if(.not.Setm2l_max)then
      if(Process.eq.51)then
         m2l_minmax(2) = infinity()
      else
         m2l_minmax(2) = infinity()
      endif
    endif
    if(.not.SetmVH_min)then
      if(Process.eq.51)then
         mVH_minmax(1) = 0d0 * GeV
      else
         mVH_minmax(1) = 0d0 * GeV
      endif
    endif
    if(.not.SetmVH_max)then
      if(Process.eq.51)then
         mVH_minmax(2) = infinity()
      else
         mVH_minmax(2) = infinity()
      endif
    endif
    if(.not.SetpTHcut)then
      if(Process.eq.51)then
         pTHcut = 0d0 * GeV
      else
         pTHcut = 0d0 * GeV
      endif
    endif
    if((Process.eq.60 .or. Process.eq.66 .or. Process.eq.67 .or. Process.eq.68 .or. Process.eq.69) .and. includeGammaStar .and. pTjetcut.le.0d0) then
       print *, " Process=",Process," with offshell photons requires a non-zero pT cut. Current setting cut ",pTjetcut/GeV," GeV is not allowed."
       stop 1
    endif
    if((Process.eq.61 .or. Process.eq.62) .and. pTjetcut.le.0d0) then
       print *, " Process=",Process," requires a non-zero pT cut. Current setting cut ",pTjetcut/GeV," GeV is not allowed."
       stop 1
    endif

    !---------------------------------------!
    !           Set OffShellV1/V2           !
    !---------------------------------------!
    if( IsAPhoton(DecayMode2) .and. IsAZDecay(DecayMode1) ) then
       print *,"Z is offshell and photon is on-shell in Z+photon production."
       print *,"Randomization of the order of writing of the decay products to the LHE file is disabled."
       RandomizeVVdecays = .false.
    elseif( IsAPhoton(DecayMode2) .and. IsAPhoton(DecayMode1) .and. Process.le.2 ) then
       print *,"Both photons are on-shell in photon+photon production."
    endif
    OffShellV1=.not.IsAPhoton(DecayMode1)
    OffShellV2=.not.IsAPhoton(DecayMode2)


    !lepton filter

    if( RequestOS.lt.RequestOSSF ) then
        RequestOS = RequestOSSF
    endif
    if( RequestNLeptons .lt. 2*RequestOS ) then
        RequestNLeptons = 2*RequestOS
    endif

    !WidthScheme and reweighting
    if( (WidthScheme.eq.0 .or. WidthSchemeIn.eq.0) .and. .not.DoPrintPMZZ ) then
        print *, "WidthScheme=0 removes the propagator entirely!  This generally only makes sense in PrintPMZZ mode."
        print *, "If you really want to use it anyway, remove this error in main.F90 and recompile."
        stop 1
    endif
    if( .not.ReadLHEFile .and. .not.DoPrintPMZZ ) then
        if( ReweightDecay .or. WidthSchemeIn.gt.0 ) then
            call Error("ReweightDecay and WidthSchemeIn only make sense in ReadLHE mode")
        endif
        if( WidthScheme.lt.0 ) then
            WidthScheme = 2
        endif
        WidthSchemeIn = WidthScheme
    else
        if( WidthSchemeIn.eq.0 .and. ReweightDecay ) then
            print *, "Can't ReweightDecay for WidthSchemeIn=0"
            stop 1
        endif
        if( WidthScheme.lt.0 .and. WidthSchemeIn.lt.0 ) then
            if( ReweightDecay ) then
                print *, "If you want to reweight the decay, you need to specify a width scheme to correct"
                print *, " for the VV branching fraction/matrix element."
                stop 1
            endif
            WidthScheme = 2
            WidthSchemeIn = 2
        elseif( WidthScheme.lt.0 .and. WidthSchemeIn.ge.0 ) then
            WidthScheme = WidthSchemeIn
        elseif( WidthScheme.ge.0 .and. WidthSchemeIn.lt.0 ) then
            WidthSchemeIn = WidthScheme
        else !both > 0
            !nothing
        endif
    endif

    if( PMZZEvals.lt.0 ) then
        !more evals for a lower mass Higgs because the integration converges slower there
        PMZZEvals = int(200000 * (1 + 24*dexp((1d0-m_Reso/(125d0*GeV))*4d0)))
        if( PMZZEvals.gt.10000000 ) PMZZEvals = 10000000
    endif

    !WriteFailedEvents

    if( WriteFailedEvents.lt.0 .or. WriteFailedEvents.gt.2 ) then
        call Error("WriteFailedEvents can only be 0, 1, or 2.  Please see the manual.")
    endif


    !---------------------------------------!
    !            Check couplings            !
    !---------------------------------------!

    if( Process.eq.1 ) then
      if( SetAnomalousSpin0gg .or. SetAnomalousSpin0VV .or. SetAnomalousHff ) then
        call Error("There is no point setting spin 0 couplings for spin 1 production")
      endif
      if( SetAnomalousSpin2gg .or. SetAnomalousSpin2qq .or. SetSpin2VV ) then
        call Error("There is no point setting spin 2 couplings for spin 1 production")
      endif
    else if( Process.eq.2 ) then
      if( SetAnomalousSpin0gg .or. SetAnomalousSpin0VV .or. SetAnomalousHff ) then
        call Error("There is no point setting spin 0 couplings for spin 2 production")
      endif
      if( SetAnomalousSpin1qq .or. SetSpin1VV ) then
        call Error("There is no point setting spin 1 couplings for spin 2 production")
      endif
    else !spin 0
      if( SetAnomalousSpin1qq .or. SetSpin1VV ) then
        call Error("There is no point setting spin 1 couplings for spin 0 production")
      endif
      if( SetAnomalousSpin2gg .or. SetAnomalousSpin2qq .or. SetSpin2VV ) then
        call Error("There is no point setting spin 2 couplings for spin 0 production")
      endif
    endif

    ! Spin-0 (incomplete)
    if( SetAnomalousSpin0gg .and. .not.Setghg2 ) then
        call Error("If you set an anomalous spin 0 gg coupling, you need to explicitly set ghg2 as well. This coupling is initialized to a non-zero value.")
    endif
    if( .not.Setghz1 .and. OffShellV1 .and. OffShellV2 .and. SetAnomalousSpin0VV) then
        call Error("If you set an anomalous spin 0 VV coupling, you need to explicitly set ghz1 as well. This couplings is initialized to a non-zero value.")
    endif
    if( SetAnomalousHff .and. .not.Setkappa ) then
        call Error("If you set an anomalous Hff coupling, you need to explicitly set kappa as well. This coupling is initialized to a non-zero value.")
    endif

    ! Contact terms
    if (Process.le.2 .or. Process.eq.50) then
        if (IsAZDecay(DecayMode1) .or. (Process.eq.50 .and. IsAPhoton(DecayMode1))) then
            if ((SetHZprime .and. .not.SetZprimeff) .or. (.not.SetHZprime .and. SetZprimeff)) then
                call Error("To use Z' contact terms, you have to set both HVZ' and Z'ff couplings")
            endif
            if ((SetMZprime.or.SetGaZprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of Z' doesn't do anything if you don't set HVZ' couplings")
            endif
            if (SetMWprime .or. SetGaWprime) then
                call Error("Don't set the W' mass and width in ZZ decay")
            endif
        elseif (IsAWDecay(DecayMode1)) then
            if ((SetHZprime .and. .not.SetWprimeff) .or. (.not.SetHZprime .and. SetWprimeff)) then
                call Error("To use W' contact terms, you have to set both HZZ'/HZ'Z' (which are used for HWW'/HW'W') and W'ff couplings")
            endif
            if ((SetMWprime.or.SetGaWprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of W' doesn't do anything if you don't set HZZ'/HZ'Z' couplings (which are used for HWW'/HW'W')")
            endif
            if (SetMZprime .or. SetGaZprime) then
                call Error("Don't set the Z' mass and width in WW decay")
            endif
        endif
    endif

    if( (Process.eq.50 .or. Process.eq.60) .and. SetZprimegammacoupling ) then
        call Error("Z'gamma couplings are not implemented for VBF or VH")
        !If you implement them and remove this error, also edit the Vprimekwargs function
        !in MELA/test/testME_more.py to not remove the Z'gamma couplings for process = 50 or 60
    endif

    ! Contact terms
    if (Process.le.2 .or. Process.eq.50) then
        if (IsAZDecay(DecayMode1) .or. (Process.eq.50 .and. IsAPhoton(DecayMode1))) then
            if ((SetHZprime .and. .not.SetZprimeff) .or. (.not.SetHZprime .and. SetZprimeff)) then
                call Error("To use Z' contact terms, you have to set both HVZ' and Z'ff couplings")
            endif
            if ((SetMZprime.or.SetGaZprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of Z' doesn't do anything if you don't set HVZ' couplings")
            endif
            if (SetMWprime .or. SetGaWprime) then
                call Error("Don't set the W' mass and width in ZZ decay")
            endif
        elseif (IsAWDecay(DecayMode1)) then
            if ((SetHZprime .and. .not.SetWprimeff) .or. (.not.SetHZprime .and. SetWprimeff)) then
                call Error("To use W' contact terms, you have to set both HZZ'/HZ'Z' (which are used for HWW'/HW'W') and W'ff couplings")
            endif
            if ((SetMWprime.or.SetGaWprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of W' doesn't do anything if you don't set HZZ'/HZ'Z' couplings (which are used for HWW'/HW'W')")
            endif
            if (SetMZprime .or. SetGaZprime) then
                call Error("Don't set the Z' mass and width in WW decay")
            endif
        endif
    endif
    if (Process.eq.60 .or. (Process.ge.66 .and. Process.lt.69)) then
        if ((SetHZprime .and. .not.SetZprimeff) .or. (.not.SetHZprime .and. SetZprimeff)) then
            call Error("To use Z' contact terms, you have to set both HVZ' and Z'ff couplings")
        endif

        if (distinguish_HWWcouplings) then
            if ((SetHWprime .and. .not.SetWprimeff) .or. (.not.SetHWprime .and. SetWprimeff)) then
                call Error("To use W' contact terms, you have to set both HVW' and W'ff couplings")
            endif
            if ((SetMZprime.or.SetGaZprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of Z' doesn't do anything if you don't set HVZ' couplings")
            endif
            if ((SetMWprime.or.SetGaWprime) .and. .not.SetHWprime) then
                call Error("Setting the mass and width of W' doesn't do anything if you don't set HVW' couplings")
            endif
        else
            if (SetHZprime .and. .not.SetWprimeff) then
                call Error("HZZ'/HZ'Z' couplings are also used for HWW'/HW'W', so if you set them you also need to set W'ff couplings (possibly to 0).")
            endif
            if (SetWprimeff .and. .not.SetHZprime) then
                call Error("If you set W'ff couplings, and you don't distinguish HZZ and HWW couplings, then you also have to set HZZ'/HZ'Z' couplings, which are also used for HWW'/HW'W'")
            endif
            if ((SetMZprime.or.SetGaZprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of Z' doesn't do anything if you don't set HVZ' couplings")
            endif
            if ((SetMWprime.or.SetGaWprime) .and. .not.SetHZprime) then
                call Error("Setting the mass and width of W' doesn't do anything if you don't set HZZ'/HZ'Z' couplings (which are used for HWW'/HW'W')")
            endif
        endif
    endif


    if( SetMZprime .and. .not.SetGaZprime ) then
        call Error("If you set the mass of Zprime, you also have to set the width! GaZprime=...")
    endif
    if( SetGaZprime .and. .not.SetMZprime ) then
        call Error("If you set the width of Zprime, you also have to set the mass! MZprime=...")
    endif
    if( SetMWprime .and. .not.SetGaWprime ) then
        call Error("If you set the mass of Wprime, you also have to set the width! GaWprime=...")
    endif
    if( SetGaWprime .and. .not.SetMWprime ) then
        call Error("If you set the width of Wprime, you also have to set the mass! MWprime=...")
    endif

    if( (Process.eq.50 .or. Process.eq.60) .and. SetZprimegammacoupling ) then
        call Error("Z'gamma couplings are not implemented for VBF or VH")
        !If you implement them and remove this error, also edit the Vprimekwargs function
        !in MELA/test/testME_more.py to not remove the Z'gamma couplings for process = 50 or 60
    endif

    ! Spin-1
    if( Process.eq.1) then
       if( SetAnomalousSpin1qq .and. .not.(Setspin1qqleft.and.Setspin1qqright) ) then
           call Error("If you set an anomalous spin 1 qq coupling, you need to set both zprime_qq_left and zprime_qq_right.")
       endif
       if (zprime_zz_1.eq.czero .and. zprime_zz_2.eq.czero) then
          call Error("For spin 1 production, you need to explicitly set zprime_zz_1 (1-) or zprime_zz_2 (1+) non-zero.")
       endif
       if (zprime_qq_left.eq.czero .and. zprime_qq_right.eq.czero) then
          call Error("For spin 1 production, zprime_qq_left and zprime_qq_right cannot both be 0.")
       endif
    endif

    ! Spin-2
    if( Process.eq.2) then
       if( SetAnomalousSpin2qq .and. .not.(Setspin2qqleft.and.Setspin2qqright) ) then
           call Error("If you set an anomalous spin 2 qq coupling, you need to explicitly set both graviton_qq_left and graviton_qq_right.")
       endif
       if ((graviton_qq_left.eq.czero .and. graviton_qq_right.eq.czero) .and. PChannel.ne.0) then
          call Error("In spin 2 qq production, the couplings cannot all be zero. You can use PChannel=0 for gg-only production.")
       endif
       if ((a1.eq.czero .and. a2.eq.czero .and. a3.eq.czero .and. a4.eq.czero .and. a5.eq.czero) .and. PChannel.ne.1) then
          call Error("In spin 2 gg production, the couplings cannot all be zero. You can use PChannel=1 for qq-only production, or explicitly set at least one of the gg couplings a1, a2, a3, a4, or a5 non-zero.")
       endif
       if (.not. SetSpin2VV) then
          call Error("For spin 2 production you need to explicitly set one of the spin 2 couplings to ZZ, Zgamma, gammagamma, ZZ', Z'Z', or Z'gamma non-zero.")
       endif
    endif

    ! VBF
    if( distinguish_HWWcouplings .and. Process.ne.60 .and. Process.ne.66 .and. Process.ne.68 ) then
        call Error("The separate HWW couplings are only used for VBF.  For H->WW decay or WH production, please set ghz* instead.")
    endif

    !Special for offshell VBF - default ghz1 is 1 to match the MCFM convention
    if( .not.Setghz1 .and. Process.ge.66 .and. Process.le.68 ) then
        !note this implies .not. SetAnomalousSpin0VV because of earlier errors
        ghz1=(1d0,0d0)
    endif


   call ComputeEWVariables()
   if(.not.SetCKMub .and. .not.SetCKMcb .and. .not.SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb)
   else if(.not.SetCKMub .and. .not.SetCKMcb .and. SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_td=VCKM_td)
   else if(.not.SetCKMub .and. SetCKMcb .and. .not.SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_cb=VCKM_cb)
   else if(SetCKMub .and. .not.SetCKMcb .and. .not.SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_ub=VCKM_ub)
   else if(SetCKMub .and. SetCKMcb .and. .not.SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_ub=VCKM_ub, inVCKM_cb=VCKM_cb)
   else if(SetCKMub .and. .not.SetCKMcb .and. SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_ub=VCKM_ub, inVCKM_td=VCKM_td)
   else if(.not.SetCKMub .and. SetCKMcb .and. SetCKMtd) then
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_cb=VCKM_cb, inVCKM_td=VCKM_td)
   else
      call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb, inVCKM_ub=VCKM_ub, inVCKM_cb=VCKM_cb, inVCKM_td=VCKM_td)
   endif

   if (DryRun) then
     print *, "Running JHUGen with these options should work."
     stop
   endif

return
END SUBROUTINE






SUBROUTINE InitVegas()
use ModKinematics
implicit none
include "vegas_common.f"
include 'maxwt.f'


  xl(1:mxdim) = 0d0
  xu(1:mxdim) = 1d0
  acc = -1d0
  nprn = 1
  readin=.false.
  writeout=.false.
  stopvegas=.false.
  evtgen=.false.

return
END SUBROUTINE



SUBROUTINE InitPDFValues()
   use ModParameters
   use ModKinematics
   implicit none

#if useLHAPDF==0
   IF (alphas_mz .LE. 0d0) THEN
      WRITE(6,*) 'alphas_mz .le. 0:',alphas_mz
      WRITE(6,*) 'continuing with alphas_mz=0.118'
      alphas_mz=0.118d0
   ENDIF
#endif

   Mu_Fact = M_Reso ! Set pdf scale to resonance mass by default, later changed as necessary in the EvalWeighted/EvalUnweighted subroutines
   Mu_Ren = M_Reso ! Set renorm. scale to resonance mass by default, later changed as necessary in the EvalWeighted/EvalUnweighted subroutines
   if (.not.ReadLHEFile .and. .not.ConvertLHEFile) call EvalAlphaS() ! Set alphas at default Mu_Ren. Notice ModParameters::ComputeQCDVariables is automatically called!
   return
END SUBROUTINE


SUBROUTINE InitPDFs()

#if useLHAPDF==1

   use ModParameters
   implicit none
   DOUBLE PRECISION alphasPDF
   character(len=5) :: LHAPDFversionnumber

   if (.not.ReadLHEFile .and. .not.ConvertLHEFile) then
     call LHAPDFversion(LHAPDFversionnumber)
     if (LHAPDFversionnumber .ge. "6.2.1") then
       call InitPDFSetByName(trim(LHAPDFString)) ! Let LHAPDF handle everything
     else
       call InitPDFSet(trim(LHAPDFString)) ! Let LHAPDF handle everything
     endif
     call InitPDF(LHAPDFMember)

     alphas_mz=alphasPDF(zmass_pdf)
     ! Dummy initialization, just in case. These values are not used.
     !nloops_pdf = 1
     zmass_pdf = M_Z
   endif

#else

   use ModParameters
   use ModKinematics
   implicit none
   character :: pdftable*(100)

   if (.not.ReadLHEFile .and. .not.ConvertLHEFile) then

     zmass_pdf = M_Z ! Take zmass_pdf=M_Z in pdfs that do not specify this value

     if( PDFSet.eq.1 ) then ! CTEQ6L1
        call SetCtq6(4)  ! 4    CTEQ6L1  Leading Order cteq6l1.tbl

        alphas_mz=0.130d0
        !nloops_pdf=1
     elseif( PDFSet.eq.3 ) then  ! NNPDF 3.0 LO with a_s=0.13
        pdftable(:)="./pdfs/NNPDF30_lo_as_0130.LHgrid"
        call NNPDFDriver(pdftable)
        call NNinitPDF(0)

        alphas_mz=0.130d0
        !nloops_pdf=1
        zmass_pdf=91.199996948242188d0*GeV
     elseif( (PDFSet.eq.2) .or. (PDFSet.ge.201 .and. PDFSet.le.240) ) then ! MSTW2008 and variations
        alphas_mz=0.13939d0
        !nloops_pdf=1
     else ! Everything else
        write(6,*) "main.F90::InitPDFs: PDFSet",PDFSet,"QCD parameters are unknown. Please double-check! Stopping JHUGen..."
        stop
        ! Could also have used these instead of the stop statement, but why introduce arbitrary number?
        !alphas_mz = 0.13229060d0
        !nloops_pdf = 1
     endif

   endif

#endif

     call InitPDFValues() ! Call this only once
   return
END SUBROUTINE




SUBROUTINE InitParameters
use ModParameters
implicit none

! rescale V branchings to preserve the correct branching proportions in partial decays
if( (DecayMode1.eq.8 .and. DecayMode2.eq.9) .or.  &
    (DecayMode1.eq.9 .and. DecayMode2.eq.8) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 2d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 2d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 2d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 1d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 1d0


elseif( (DecayMode1.eq.0 .and. DecayMode2.eq.9) .or.  &
        (DecayMode1.eq.9 .and. DecayMode2.eq.0) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 2d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 2d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 2d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 1d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 2d0


elseif( (DecayMode1.eq.1 .and. DecayMode2.eq.9) .or.  &
        (DecayMode1.eq.9 .and. DecayMode2.eq.1) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 1d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 1d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 2d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 2d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 2d0


elseif( (DecayMode1.eq.2 .and. DecayMode2.eq.9) .or.  &
        (DecayMode1.eq.9 .and. DecayMode2.eq.2) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 2d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 2d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 2d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 2d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 1d0


elseif( (DecayMode1.eq.3 .and. DecayMode2.eq.9) .or.  &
        (DecayMode1.eq.9 .and. DecayMode2.eq.3) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 2d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 2d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 1d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 2d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 2d0

elseif( (DecayMode1.eq.8 .and. DecayMode2.eq.0) .or.  &
        (DecayMode1.eq.0 .and. DecayMode2.eq.8) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 1d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 1d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 1d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 1d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 2d0

elseif( (DecayMode1.eq.8 .and. DecayMode2.eq.2) .or.  &
        (DecayMode1.eq.2 .and. DecayMode2.eq.8) ) then
        scale_alpha_Z_uu = scale_alpha_Z_uu * 1d0
        scale_alpha_Z_dd = scale_alpha_Z_dd * 1d0
        scale_alpha_Z_nn = scale_alpha_Z_nn * 1d0
        scale_alpha_Z_ll = scale_alpha_Z_ll * 2d0
        scale_alpha_Z_tt = scale_alpha_Z_tt * 1d0

elseif( (DecayMode1.eq.4 .and. DecayMode2.eq.11) .or.  &
        (DecayMode1.eq.11.and. DecayMode2.eq.4) ) then
        scale_alpha_W_ud = scale_alpha_W_ud * 2d0
        scale_alpha_W_cs = scale_alpha_W_cs * 2d0
        scale_alpha_W_ln = scale_alpha_W_ln * 1d0
        scale_alpha_W_tn = scale_alpha_W_tn * 2d0


elseif( (DecayMode1.eq.5 .and. DecayMode2.eq.11) .or.  &
        (DecayMode1.eq.11.and. DecayMode2.eq.5) ) then
        scale_alpha_W_ud = scale_alpha_W_ud * 1d0
        scale_alpha_W_cs = scale_alpha_W_cs * 1d0
        scale_alpha_W_ln = scale_alpha_W_ln * 2d0
        scale_alpha_W_tn = scale_alpha_W_tn * 2d0

elseif( (DecayMode1.eq.6 .and. DecayMode2.eq.11) .or.  &
        (DecayMode1.eq.11.and. DecayMode2.eq.6) ) then
        scale_alpha_W_ud = scale_alpha_W_ud * 2d0
        scale_alpha_W_cs = scale_alpha_W_cs * 2d0
        scale_alpha_W_ln = scale_alpha_W_ln * 2d0
        scale_alpha_W_tn = scale_alpha_W_tn * 1d0

elseif( (DecayMode1.eq.10.and. DecayMode2.eq.11) .or.  &
        (DecayMode1.eq.11.and. DecayMode2.eq.10) ) then
        scale_alpha_W_ud = scale_alpha_W_ud * 2d0
        scale_alpha_W_cs = scale_alpha_W_cs * 2d0
        scale_alpha_W_ln = scale_alpha_W_ln * 1d0
        scale_alpha_W_tn = scale_alpha_W_tn * 1d0

elseif( (DecayMode1.eq.4 .and. DecayMode2.eq.10) .or.  &
        (DecayMode1.eq.10.and. DecayMode2.eq.4) ) then
        scale_alpha_W_ud = scale_alpha_W_ud * 1d0
        scale_alpha_W_cs = scale_alpha_W_cs * 1d0
        scale_alpha_W_ln = scale_alpha_W_ln * 1d0
        scale_alpha_W_tn = scale_alpha_W_tn * 2d0

elseif( (DecayMode1.eq.6 .and. DecayMode2.eq.10) .or.  &
        (DecayMode1.eq.10.and. DecayMode2.eq.6) ) then
        scale_alpha_W_ud = scale_alpha_W_ud * 1d0
        scale_alpha_W_cs = scale_alpha_W_cs * 1d0
        scale_alpha_W_ln = scale_alpha_W_ln * 2d0
        scale_alpha_W_tn = scale_alpha_W_tn * 1d0
endif



END SUBROUTINE



SUBROUTINE InitProcess()
use ModParameters
use ModMisc
use ModCrossSection
use ModCrossSection_VH
use ModTTBHiggs
#if useCollier==1
use ModCrossSection_HH
use Collier
#endif
implicit none
include "vegas_common.f"

! NOTE: NDim for weighted vegas run is not minimal

      NDim = 0
      NDim = NDim + 7    ! PS integration
      if( .not.ReadLHEFile ) NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! offzchannel
      if( includeInterference ) NDim = NDim + 1    ! for interference
      if( OffShellV1 .or. OffShellV2 ) then
          NDim = NDim + 2    ! integration over Z's invariant mass
      else
          NDim = NDim + 3    ! integration over mz1, mz2, mg for adjusted kinematics but with on-shell matrix elements
      endif
      VegasIt1_default = 5
      if( .not.ReadLHEFile ) then
          VegasNc0_default = 10000000
          VegasNc1_default = 50000000
          VegasNc2_default = 10000
      else
          VegasNc0_default = 10000000
          VegasNc1_default = 1000000000! this should be "infinity" such that the loop runs until the whole input LHE file is read
          VegasNc2_default = 10000
      endif

      !- gg-->spin0-->4f
      if(Process.eq.0 ) then
        NDim = 11
      endif
      !- gg-->spin1-->4f
      if(Process.eq.1 ) then
        NDim = 12
      endif
      !- gg-->spin2-->4f
      if(Process.eq.2 ) then
        NDim = 12
      endif

      !- HVBF
      if(Process.eq.60) then
         NDim = 5        ! for phase space
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! for pdf sampling
         NDim = NDim + 1 ! for PS sampling
         VegasIt1_default = 15
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif

      !- ppJJ with decays
      if(Process.ge.66 .and. Process.le.69) then
         NDim = 5
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 8
         NDim = NDim + 1
         NDim = NDim + 1
         NDim = NDim + 1

         VegasIt1_default = 15
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif

      !- Hjj, gluon fusion
      if(Process.eq.61) then
         NDim = 5        ! phase space
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! for pdf sampling
         NDim = NDim + 1 ! for PS sampling
         VegasIt1_default = 15
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- Hj, gluon fusion
      if(Process.eq.62) then
         NDim = 5+1 !1 for color flow ramdomization in gg > Hg
         NDim = NDim + 2 ! sHat integration
         if( unweighted ) NDim = NDim + 1  ! random number which decides if event is accepted

         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- VHiggs
      if(Process.eq.50) then
         NDim = 17
         NDim = NDim + 2 ! sHat integration
         if( unweighted ) NDim = NDim + 1  ! random number which decides if event is accepted

         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- VH
      if(Process.eq.51) then
#if useCollier==1
         call Init_cll(4,3,trim(DataFile)//"_collierfiles")
         call InitCacheSystem_cll(1,4)
         call InitEvent_cll
         call setMode_cll(1)!1. use COLI branch; 2. use DD branch; 3. use both branches and compare.
#else
         print *, "Need to link COLLIER for this process."
         print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
         print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
         print *, "specified in the makefile."
         stop 1
#endif
! if Collier is used
         NDim = 19
         if( unweighted ) NDim = NDim + 3  ! partonic channel and acceptance
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5): flavor in Z/W decay mode
! yrnd(6:7): phi_4 and cos(theta_4) in the CM frame of Z*(3)
! yrnd(8:9): phi_6 and cos(theta_6) in the CM frame of decay product of Z(4)
! yrnd(10:11): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
! yRnd(12): inv_mass(4)
! yRnd(13): inv_mass(5)
! yRnd(14:16): real emission phase space
! yRnd(17): NLO integration for + distribution for PDF renormalization
! yrnd(18:19): PDF mapping
! yrnd(20): flavor of j in gq > WH+j
! yrnd(21): partonic channel in unweighted events
! yRnd(22): accept or reject in unweighted events

         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- HH
      if(Process.eq.52) then
#if useCollier==1
         call Init_cll(4,4,trim(DataFile)//"_collierfiles")
         call InitCacheSystem_cll(1,4)
         call InitEvent_cll
         call setMode_cll(1)!1. use COLI branch; 2. use DD branch; 3. use both branches and compare.
#else
         print *, "Need to link COLLIER for this process."
         print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
         print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
         print *, "specified in the makefile."
         stop 1
#endif
! if Collier is used
         NDim = 15
         if( unweighted ) NDim = NDim + 2  ! partonic channel and acceptance
! yrnd(1:2): helicity of parton 1:2
! yrnd(3): helicities of parton 6,7
! yrnd(4): helicities of parton 8,9
! yrnd(5:6): phi_4 and cos(theta_4) in the CM frame of 1+2(3)
! yrnd(7:8): phi_6 and cos(theta_6) in the CM frame of decay product of H(4)
! yrnd(9:10): phi_8 and cos(theta_8) in the CM frame of decay product of H(5)
! yRnd(11): inv_mass(4)
! yRnd(12): inv_mass(5)
! yRnd(13): swap momenta in PS for stability
! yrnd(14:15): PDF mapping
! yrnd(16): partonic channel in unweighted events
! yRnd(17): accept or reject in unweighted events
         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- ttbar+H
      if(Process.eq.80) then
         call InitProcess_TTBH()
         NDim = 5                                  ! phase space
         if( TopDecays.ne.0 ) NDim = Ndim + 8      ! phase space
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1  ! MC sampling for gg and qqb channel

         VegasIt1_default = 5
         VegasNc0_default =  100000
         VegasNc1_default =  500000
         VegasNc2_default =    1000
      endif
      !- bbbar+H
      if(Process.eq.90) then
         TopDecays = 0
         m_Top = m_Bot

         call InitProcess_TTBH()
         NDim = 12
         NDim = NDim + 2 ! sHat integration
         if( unweighted ) NDim = NDim + 1  ! random number which decides if event is accepted

         VegasIt1_default = 5
         VegasNc0_default =  100000
         VegasNc1_default =  500000
         VegasNc2_default =    1000
      endif
     ! RR added -- t+H
      if(Process.eq.110) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! select partonic channel

         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif
     ! RR added -- tb+H
      if(Process.eq.111) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! select partonic channel

         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif

     ! RR added -- t+H s-schannel
      if(Process.eq.112) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! select partonic channel

         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif
     ! RR added -- tb+H s-channel
      if(Process.eq.113) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! select partonic channel

         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif
     ! sum of all TH channels
      if(Process.eq.114) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 1 ! select partonic channel

         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif

END SUBROUTINE


#if linkMELA==1
subroutine SetupMCFM(Process)
use ModMCFMWrapper
use ModMisc
implicit none
integer, intent(in) :: Process

   if( &
      (Process.ge.66 .and. Process.le.69) &
      ) then
      call MCFM_firsttime()
      if( &
         (Process.ge.66 .and. Process.le.69) &
         ) then
         call Setup_MCFM_qqVVqq_firsttime(Process)
      !else if ... add more processes here
      else
         call Error("Don't know how to set up this process for MCFM")
      endif
   endif

end subroutine
#endif




SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModCrossSection
use ModCrossSection_VH
#if useCollier==1
use ModCrossSection_HH
#endif
use ModCrossSection_HJJ
use ModCrossSection_TTBH
use ModCrossSection_BBBH
use ModCrossSection_TH
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22)
real(8) :: dum, RES(-5:5,-5:5),VG2(-5:5,-5:5)
integer :: i, i1, j1,PChannel_aux, PChannel_aux1,NHisto
include 'csmaxvalue.f'
integer :: flav1,flav2,StatusPercent,LastStatusPercent=0
integer, dimension(:), allocatable :: gen_seed


if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 .and.  (unweighted) ) then
      VegasNc1 = VegasNc1_default
      VegasNc2 = VegasNc2_default
endif
if( VegasNc1.eq.-1 .and.  .not. (unweighted) ) VegasNc1 = VegasNc1_default
if( VegasNc2.eq.-1 .and.  .not. (unweighted) ) VegasNc2 = VegasNc2_default



   warmup = .false.
   itmx = VegasIt1
   ncall= VegasNc1

   PChannel_aux = PChannel

if ( (unweighted.eqv..false.) .or. (GenerateEvents.eqv..true.) ) then  !----------------------- weighted events


    if( (GenerateEvents.eqv..true.) ) then
        itmx = 1
        ncall= VegasNc1
        call vegas(EvalOnlyPS,VG_Result,VG_Error,VG_Chi2)
        return
    endif

    ! WARM-UP RUN
    itmx = VegasIt1
    ncall= VegasNc1
    warmup = .true.

    if (Process.eq.60 .or. Process.eq.61) then
      call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.62) then
      call vegas(EvalWeighted_HJ,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.50) then
      call vegas(EvalWeighted_VHiggs,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.51 .or. Process.eq.52) then
#if useCollier==1
      if (Process.eq.51) then
        call vegas(EvalWeighted_VH,VG_Result,VG_Error,VG_Chi2)
      elseif (Process.eq.52) then
        call vegas(EvalWeighted_HH,VG_Result,VG_Error,VG_Chi2)
      endif
#else
      print *, "Need to link COLLIER for this process."
      print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
      print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
      print *, "specified in the makefile."
      stop 1
#endif
    elseif (Process.eq.80) then
      call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.90) then
      call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.110 .or. Process.eq.111 .or. Process.eq.112 .or. Process.eq.113)  then
      call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
   else
      call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events
!       call vegas(EvalWeighted_tautau,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events
    endif

    !DATA RUN
    call ClearHisto()
    warmup = .false.
    EvalCounter=0
    RejeCounter=0
    AccepCounter=0
    Br_counter(:,:) = 0
    AlertCounter=0

    avgcs = 0d0

    itmx = 1!!!!!!!!!!!!!!!!!!!!
    ncall= VegasNc2
    if (process.eq.60 .or. process.eq.61) then
      call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    elseif (process.eq.62 .or. process.eq.61) then
      call vegas1(EvalWeighted_HJ,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.50) then
      call vegas1(EvalWeighted_VHiggs,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.51 .or. Process.eq.52) then
#if useCollier==1
      if (Process.eq.51) then
        call vegas1(EvalWeighted_VH,VG_Result,VG_Error,VG_Chi2)
      elseif (Process.eq.52) then
        call vegas1(EvalWeighted_HH,VG_Result,VG_Error,VG_Chi2)
      endif
#else
      print *, "Need to link COLLIER for this process."
      print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
      print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
      print *, "specified in the makefile."
      stop 1
#endif
    elseif (Process.eq.80) then
      call vegas1(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.90) then
      call vegas1(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.110 .or. Process.eq.111 .or. Process.eq.112 .or. Process.eq.113) then
      call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
    else
      call vegas1(EvalWeighted,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events
!       call vegas1(EvalWeighted_tautau,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events
    endif



elseif(unweighted.eqv..true.) then  !----------------------- unweighted events

    VG = zero
    VG2= zero
    csmax = zero

    if( .not. ReadCSmax ) then
        print *, " finding maximal weight with ",VegasNc0," evaluations"
        warmup = .true.
        do i=1,VegasNc0
            call random_number(yRnd)
            if (Process.eq.60 .or. Process.eq.61) then
                RES = 0d0
                dum = EvalUnWeighted_HJJ(yRnd,.false.,(/-99,-99/),RES)
                VG = VG + RES
            elseif (Process.eq.62) then
                RES = 0d0
                dum = EvalUnWeighted_HJ(yRnd,.false.,RES)
                VG = VG + RES
            elseif (Process.eq.50) then
                RES = 0d0
                dum = EvalUnWeighted_VHiggs(yRnd,.false.,RES)
                VG = VG + RES
            elseif (Process.eq.51 .or. Process.eq.52) then
#if useCollier==1
                if (Process.eq.51) then
                    RES = 0d0
                    dum = EvalUnWeighted_VH(yRnd,.false.,RES)
                    VG = VG + RES
                elseif (Process.eq.52) then
                    RES = 0d0
                    dum = EvalUnWeighted_HH(yRnd,.false.,RES)
                    VG = VG + RES
                endif
#else
                print *, "Need to link COLLIER for this process."
                print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
                print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
                print *, "specified in the makefile."
                stop 1
#endif
            else
                if (PChannel_aux.eq.0.or.PChannel_aux.eq.2) then
                    PChannel= 0
                    RES = 0d0
                    dum = EvalUnWeighted(yRnd,.false.,RES)
                    VG(0,0) = VG(0,0) + RES(0,0)
                endif
                if(PChannel_aux.eq.1.or.PChannel_aux.eq.2) then
                    PChannel = 1
                    RES = 0d0
                    dum = EvalUnWeighted(yRnd,.false.,RES)
                    VG = VG + RES
                endif
            endif
            PChannel = PChannel_aux
        enddo

        open(unit=io_CSmaxFile,file=trim(CSmaxFile)//'_CSmax.bin',form='unformatted',status='replace')
        WRITE(io_CSmaxFile) CSMAX
        close(io_CSmaxFile)
    else
        open(unit=io_CSmaxFile,file=trim(CSmaxFile)//'_CSmax.bin',form='unformatted')
        READ(io_CSmaxFile) CSMAX
        close(io_CSmaxFile)
    endif

     VG = VG/dble(VegasNc0)
     csmax   = 1.5d0*csmax    !  adjustment factors, can be choosen  separately channel/by/channel



!        print *, " gg/qqb ratio = ", VG(0,0)/(VG(+1,-1) + VG(+2,-2) + VG(+3,-3) + VG(+4,-4) + VG(+5,-5)   &
!                                             +VG(-1,+1) + VG(-2,+2) + VG(-3,+3) + VG(-4,+4) + VG(-5,+5))

!------------------adj_par fixes by how much the quark-induced channels need to be adjusted

    if (PChannel.eq.2.and.fix_channels_ratio) then
        adj_par = VG(0,0)/(VG(-5,5)+VG(-4,4)+VG(-3,3)+VG(-2,2)+VG(-1,1)  &
                + VG(1,-1)+VG(2,-2)+VG(3,-3)+VG(4,-4)+VG(5,-5))*channels_ratio_fix/(one-channels_ratio_fix)
    else
        adj_par = 1d0
    endif


!--- rescale the gluon induced channel
    csmax(0,0) = csmax(0,0)/adj_par


!------------------------------- set counts to zero for actual evaluation

    EvalCounter = 0
    AccepCounter = 0
    RejeCounter = 0
    AlertCounter = 0
    AccepCounter_part = 0

    call cpu_time(time_start)
    warmup = .true.
    if( VegasNc1.ne.-1 ) then
        print *, " generating events with ",VegasNc1," tries"
        do i=1,VegasNc1
            call random_number(yRnd)
            if (Process.eq.60 .or. Process.eq.61) then
                      dum = EvalUnWeighted_HJJ(yRnd,.true.,(/-99,-99/),RES)! RES is a dummy here
            elseif (Process.eq.62) then
                      dum = EvalUnWeighted_HJ(yRnd,.true.,RES)! RES is a dummy here
            elseif (Process.eq.50) then
                      dum = EvalUnWeighted_VHiggs(yRnd,.true.,RES)! RES is a dummy here
            elseif (Process.eq.51 .or. Process.eq.52) then
#if useCollier==1
                if (Process.eq.51) then
                          dum = EvalUnWeighted_VH(yRnd,.true.,RES)! RES is a dummy here
                elseif (Process.eq.52) then
                          dum = EvalUnWeighted_HH(yRnd,.true.,RES)! RES is a dummy here
                endif
#else
                print *, "Need to link COLLIER for this process."
                print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
                print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
                print *, "specified in the makefile."
                stop 1
#endif
            else
                dum = EvalUnWeighted(yRnd,.true.,RES)! RES is a dummy here
            endif
        enddo

    elseif( VegasNc2.ne.-1 ) then
        AccepCounter = 0
        print *, " generating ",VegasNc2," events"
        do while( AccepCounter.lt.VegasNc2 )
              call random_number(yRnd)
              if (Process.eq.60 .or. Process.eq.61) then
                dum = EvalUnWeighted_HJJ(yRnd,.true.,(/-99,-99/),RES)! RES is a dummy here
              elseif (Process.eq.62) then
                dum = EvalUnWeighted_HJ(yRnd,.true.,RES)! RES is a dummy here
              elseif (Process.eq.50) then
                  dum = EvalUnWeighted_VHiggs(yRnd,.true.,RES)! RES is a dummy here
              elseif (Process.eq.51 .or. Process.eq.52) then
#if useCollier==1
                  if (Process.eq.51) then
                      dum = EvalUnWeighted_VH(yRnd,.true.,RES)! RES is a dummy here
                  elseif (Process.eq.52) then
                      dum = EvalUnWeighted_HH(yRnd,.true.,RES)! RES is a dummy here
                  endif
#else
                  print *, "Need to link COLLIER for this process."
                  print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
                  print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
                  print *, "specified in the makefile."
                  stop 1
#endif
              else
                  dum = EvalUnWeighted(yRnd,.true.,RES)! RES is a dummy here
              endif
              StatusPercent = int(100d0*AccepCounter/dble(VegasNc2))
              if( mod(StatusPercent,10).eq.0 .and. LastStatusPercent.ne.StatusPercent ) then
                 write(io_stdout,"(X,I3,A)") int(100d0*AccepCounter/dble(VegasNc2)),"% "
!                 write(io_stdout,"(X,I3,A)",advance='no') int(100d0*AccepCounter/dble(VegasNc2)),"% "
!                 flush(io_stdout)
                 LastStatusPercent = StatusPercent
              endif
        enddo




    else
        print *, "ERROR: VegasNc1 and VegasNc2 must not be set at the same time"
        stop
    endif
    call cpu_time(time_end)
    print *, ""



    print *, " Evaluation Counter: ",EvalCounter
    print *, " Acceptance Counter: ",AccepCounter
    print *, " Rejection  Counter: ",RejeCounter
    do i1=-5,5
      print *, " Acceptance  Counter_part: ", i1, AccepCounter_part(i1,-i1)
    enddo
    print *, " Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        print *, "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        print *, "       Increase CSMAX in main.F90."
    endif
    write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)


  endif! unweighted


    if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) then
          write(*,*) ""
          write(*,"(A)") "                 el              mu             tau             neu              jet"
          write(*,"(A,5F16.4)") " el ",dble(Br_counter(1,1:5))/dble(AccepCounter)
          write(*,"(A,5F16.4)") " mu ",dble(Br_counter(2,1:5))/dble(AccepCounter)
          write(*,"(A,5F16.4)") " tau",dble(Br_counter(3,1:5))/dble(AccepCounter)
          write(*,"(A,5F16.4)") " neu",dble(Br_counter(4,1:5))/dble(AccepCounter)
          write(*,"(A,5F16.4)") " jet",dble(Br_counter(5,1:5))/dble(AccepCounter)
          write(*,*) ""
          write(*,"(A,5F16.3)") "llll: ",(dble(Br_counter(1,1))+dble(Br_counter(2,2))+dble(Br_counter(3,3)))/dble(AccepCounter)
          write(*,"(A,5F16.3)") "llLL: ",(dble(Br_counter(1,2))+dble(Br_counter(1,3))+  &
                                          dble(Br_counter(2,1))+dble(Br_counter(2,3))+  &
                                          dble(Br_counter(3,1))+dble(Br_counter(3,2)))/dble(AccepCounter)
          write(*,"(A,5F16.3)") "2l2q: ",(dble(Br_counter(1,5))+dble(Br_counter(2,5))+dble(Br_counter(3,5))+  &
                                          dble(Br_counter(5,1))+dble(Br_counter(5,2))+dble(Br_counter(5,3)) )/dble(AccepCounter)

          write(*,"(A,5F16.3)") "4l/2q2l: ",(dble(Br_counter(1,2))+dble(Br_counter(1,3))+ dble(Br_counter(1,1))+dble(Br_counter(2,2))+dble(Br_counter(3,3)) &
                                        + dble(Br_counter(2,1))+dble(Br_counter(2,3))+   &
                                          dble(Br_counter(3,1))+dble(Br_counter(3,2)))/  &
                                          (dble(Br_counter(1,5))+dble(Br_counter(2,5))+dble(Br_counter(3,5))+  &
                                          dble(Br_counter(5,1))+dble(Br_counter(5,2))+dble(Br_counter(5,3)) )

      !     print *, alpha_QED/12d0*M_Z * (   (aR_lep+aL_lep)**2 + (aR_lep-aL_lep)**2        &
      !                                      +(aR_lep+aL_lep)**2 + (aR_lep-aL_lep)**2        &
      !                                      +(aR_lep+aL_lep)**2 + (aR_lep-aL_lep)**2        &
      !                                      +(aR_neu+aL_neu)**2 + (aR_neu-aL_neu)**2        &
      !                                      +(aR_neu+aL_neu)**2 + (aR_neu-aL_neu)**2        &
      !                                      +(aR_neu+aL_neu)**2 + (aR_neu-aL_neu)**2        &
      !                                      +((aR_Qup+aL_Qup)**2 + (aR_Qup-aL_Qup)**2 )*3d0 *1.0366d0 &
      !                                      +((aR_Qdn+aL_Qdn)**2 + (aR_Qdn-aL_Qdn)**2 )*3d0 *1.0366d0 &
      !                                      +((aR_Qup+aL_Qup)**2 + (aR_Qup-aL_Qup)**2 )*3d0 *1.0366d0 &
      !                                      +((aR_Qdn+aL_Qdn)**2 + (aR_Qdn-aL_Qdn)**2 )*3d0 *1.0366d0 &
      !                                      +((aR_Qdn+aL_Qdn)**2 + (aR_Qdn-aL_Qdn)**2 )*3d0 *1.0366d0 *(1d0-3d0/2d0/(M_Z/2d0/m_bot)**2) &
      !                                  )/4d0/(one-sitW**2)/sitW**2

    endif



return
END SUBROUTINE








SUBROUTINE StartVegas_NEW(VG_Result,VG_Error)
use ModCrossSection
use ModCrossSection_VH
#if useCollier==1
use ModCrossSection_HH
#endif
use ModCrossSection_BBBH
use ModCrossSection_HJJ
use ModCrossSection_TH
use ModCrossSection_TTBH
use modHiggsJJ
use modTHiggs
use ModKinematics
use ModMisc
use ModParameters
implicit none
include "vegas_common.f"
include 'maxwt.f'
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),calls1,calls2,calls_rescale
real(8) :: dum, RES(-5:5,-5:5),ResFrac(-5:5,-5:5),TotalXSec
integer :: i, j, i1, j1,PChannel_aux, PChannel_aux1,NHisto
integer, pointer :: ijSel(:,:)
include 'csmaxvalue.f'
integer :: flav1,flav2,StatusPercent,MissingEvents,MaxEvts,imax
integer :: VegasSeed,PreviousSum,ios,NumPartonicChannels
integer :: VBFoffsh_Hash_Size,VBFoffsh_run_size
integer, parameter :: VBFoffsh_run_maxsize=NMAXPSPARTITIONS
character :: ProcessStr*(3)
logical :: UseBetaVersion=.false.
real(8) :: VG_Result_in(1:NMAXPSPARTITIONS),VG_Error_in(1:NMAXPSPARTITIONS),calls1_in(1:NMAXPSPARTITIONS),calls2_in(1:NMAXPSPARTITIONS)
real(8) :: CrossSec2_in(1:VBFoffsh_run_maxsize,1:NMAXCHANNELS),CrossSecMax2_in(1:VBFoffsh_run_maxsize,NMAXCHANNELS),CrossSectionWithWeights_in(1:VBFoffsh_run_maxsize),CrossSectionWithWeightsErrorSquared_in(1:VBFoffsh_run_maxsize)

    VBFoffsh_Hash_Size = 0
    VBFoffsh_run_size = 0
    VG_Result = -13d0
    VG_Error  = -13d0

    if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
    if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
    if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 .and.  (unweighted) ) then
          VegasNc1 = VegasNc1_default
          VegasNc2 = VegasNc2_default
    endif
    if( VegasNc1.eq.-1 .and.  .not. (unweighted) ) VegasNc1 = VegasNc1_default
    if( VegasNc2.eq.-1 .and.  .not. (unweighted) ) VegasNc2 = VegasNc2_default



    if(Process.lt.10) then
      write(ProcessStr,"(I1)") Process
      ProcessStr="0"//trim(ProcessStr)
    elseif(Process.lt.100) then
      write(ProcessStr,"(I2)") Process
    else
      write(ProcessStr,"(I3)") Process
    endif

! nprn = -1


    call cpu_time(time_start)
    warmup = .false.
    itmx = VegasIt1
    ncall= VegasNc1
    PChannel_aux = PChannel


    outgridfile=trim(CSmaxFile)//'.grid'
    ingridfile=trim(outgridfile)


    if( Process.gt.69 .and. Process.le.79 ) call Error("Missing ChannelHash for the process") !can change 79 if we add more processes in between

if ( (unweighted.eqv..false.) .or. (GenerateEvents.eqv..true.) ) then  !----------------------- weighted events


    if( (GenerateEvents.eqv..true.) ) then
        itmx = 1
        ncall= VegasNc1
        call vegas(EvalOnlyPS,VG_Result,VG_Error,VG_Chi2)
        return
    endif



    ! WARM-UP RUN
    if( .not. ReadCSmax ) then
      readin=.false.
      writeout=.true.

      itmx = VegasIt1
      ncall= VegasNc1
      warmup = .true.
      if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.80 ) call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.90 ) call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)

      if( Process.eq.60 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
      if( Process.ge.66 .and. Process.le.69 ) call vegas(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.61 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)

      if( Process.eq.110) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.111) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.112) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.113) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
      if( Process.eq.114) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
    endif

    !DATA RUN
    writeout=.false.
    if( ReadCSmax ) then
        readin=.true.
    else
        readin=.false.
    endif

    call ClearHisto()
    warmup = .false.
    Br_counter(:,:) = 0
    EvalCounter=0
    RejeCounter=0
    AccepCounter=0
    AlertCounter=0
    avgcs = 0d0
    itmx = 2
    ncall= VegasNc2
    if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) call vegas1(EvalWeighted,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.80 ) call vegas1(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.90 ) call vegas1(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)

    if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)

    if( Process.eq.110) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.111) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.112) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.113) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.114) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)




elseif(unweighted.eqv..true.) then  !----------------------- unweighted events

UseBetaVersion = CalculatesXsec(Process)


if( UseBetaVersion ) then
! !-------------------new stuff -------------------

IF( .NOT. (Process.ge.66 .and. Process.le.69) ) THEN! special treatment for offshell VBF

    write(io_stdout,"(2X,A)")  "Scanning the integrand"
    warmup = .true.
    itmx = 10
    ncall= VegasNc0


    if( ReadCSmax ) then
        readin=.true.
        writeout=.false.
        ingridfile = trim(CSmaxFile)//'_step2.grid'
        open(unit=io_TmpFile,file=trim(CSmaxFile)//'_gridinfo.txt',form='formatted',status='old')
        read(io_TmpFile,fmt=*) calls1
        read(io_TmpFile,fmt=*) CrossSec
        read(io_TmpFile,fmt=*) CrossSecMax
        read(io_TmpFile,fmt=*) VG_Result
        read(io_TmpFile,fmt=*) VG_Error
        close(unit=io_TmpFile)
    else
        readin=.false.
        writeout=.true.
        if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.60 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.61 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
!        if( Process.ge.66 .and. Process.le.69 ) call vegas(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.110) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.111) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.112) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.113) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.114) call vegas(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)

        itmx = 2
        writeout=.true.
        outgridfile = trim(CSmaxFile)//'_step2.grid'
        if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) call vegas1(EvalWeighted,VG_Result,VG_Error,VG_Chi2)
!         if( Process.eq.80 ) call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2) ! adjust to LHE format
!         if( Process.eq.90 ) call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
!        if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.110) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.111) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.112) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.113) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.114) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)


        CrossSecMax(:,:) = 0d0
        CrossSec(:,:) = 0d0
        print *, "resetting CrossSecMax(:,:)"
        itmx = 1
        if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) call vegas1(EvalWeighted,VG_Result,VG_Error,VG_Chi2)
!         if( Process.eq.80 ) call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2) ! adjust to LHE format
!         if( Process.eq.90 ) call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
!        if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.110) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.111) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.112) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.113) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.114) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        writeout=.false.
        ingridfile=trim(outgridfile)

        call vegas_get_calls(calls1)
        CrossSec(:,:) = CrossSec(:,:)/dble(itmx)

        open(unit=io_TmpFile,file=trim(CSmaxFile)//'_gridinfo.txt',form='formatted',status='replace')
        write(io_TmpFile,fmt=*) calls1
        write(io_TmpFile,fmt=*) CrossSec
        write(io_TmpFile,fmt=*) CrossSecMax
        write(io_TmpFile,fmt=*) VG_Result
        write(io_TmpFile,fmt=*) VG_Error
        close(unit=io_TmpFile)
    endif

    write(io_stdout,"(A)")  ""
    write(io_stdout,*) "Total xsec: ",VG_Result, " +/-",VG_Error, " fb    vs.",sum(CrossSec(:,:))
    call InitOutput(VG_Result, VG_Error)


    RequEvents(:,:) = 0
    if (VegasNc2.ne.-1) then
      call HouseOfRepresentatives(CrossSec(-5:5,-5:5), RequEvents(-5:5,-5:5), VegasNc2)
    endif

    if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2) then
       call getRef_PPXchannelHash(ijSel)
    elseif( Process.eq.60 ) then
       call getRef_VBFchannelHash(ijSel)
    elseif( Process.eq.61 ) then
       call getRef_HJJchannelHash(ijSel)
    elseif( Process.ge.110 .and. Process .le.114 ) then
       call getRef_THchannelHash(ijSel)
    else
       call getRef_GENchannelHash(ijSel)
    endif

    do i=1,121
         i1 = ijSel(i,1)
         j1 = ijSel(i,2)
         if( RequEvents(i1,j1).gt.0 .and. ijSel(i,3).ge.0 ) write(io_stdout,"(1X,I3,A,I3,I3,I4,A,3X,F8.3,I9)") i," Fractional partonic xsec ",i1,j1,ijSel(i,3)," "//getLHEParticle(i1)//" "//getLHEParticle(j1)//" ",CrossSec(i1,j1)/VG_Result,RequEvents(i1,j1)
    enddo
    write(io_stdout,"(2X,A,F8.3,I9)") "Sum        partonic xsec   x   x    ",sum(CrossSec(:,:))/VG_Result,sum(RequEvents(:,:))






    write(io_stdout,"(A)")  ""
    write(io_stdout,"(2X,A)")  "Event generation"

    call ClearHisto()
    warmup = .false.
    evtgen = .true.
    itmx = 1
!     nprn = 0
    EvalCounter = 0
    Br_counter(:,:) = 0
    RejeCounter = 0
    AlertCounter = 0
    AccepCounter_part(:,:) = 0
    StatusPercent = 0d0

    CrossSecMax(:,:) = 1.0d0 * CrossSecMax(:,:)    !  adjustment factor
    call cpu_time(time_start)



    itmx=1
    ncall= VegasNc0    !/10     dmax inside vegas needs to be adapted for this. or at least thisdmax inside mod_Crosssection
    nprn =-1 !0
    writeout=.false.

    call vegas_get_calls(calls2)
    calls_rescale = calls1/calls2
    CrossSecMax(:,:) = CrossSecMax(:,:) * calls_rescale
    print *, "Rescale CrossSecMax by ",calls_rescale

    PreviousSum = 0
    if( sum(RequEvents(:,:)).eq.0 ) StatusPercent = 100d0
    do while( StatusPercent.lt.100d0  )
!     do while( AccepCounter_part(iPart_sel,jPart_sel).lt.RequEvents(iPart_sel,jPart_sel) )
        call cpu_time(time_start)
        readin=.true.  ! this prevents adapting the grid during this while-loop

        if( Process.eq.0 .or. Process.eq.1 .or. Process.eq.2 ) call vegas1(EvalWeighted,VG_Result,VG_Error,VG_Chi2)
!         if( Process.eq.80 ) call vegas1(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)! adjust to LHE format
    !     if( Process.eq.90 ) call vegas1(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
 !       if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.110) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.111) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.112) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.113) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.114) call vegas1(EvalWeighted2_TH,VG_Result,VG_Error,VG_Chi2)


!         call system('clear')
        write(io_stdout,*) ""
        do i1=-5,5
        do j1=-5,5
            if( RequEvents(i1,j1).gt.0 .and. AccepCounter_part(i1,j1).lt.RequEvents(i1,j1)  ) then
!                write(io_stdout,"(1X,A,I4,I4,2X,I7,2X,I7,2X,F8.3,1X,A)") "Generated events ", i1,j1,(AccepCounter_part(i1,j1)),(RequEvents(i1,j1)),dble(AccepCounter_part(i1,j1))/dble(RequEvents(i1,j1))*100d0,"%"
!                call PrintStatusBar2(int(dble(AccepCounter_part(i1,j1))/(dble(RequEvents(i1,j1)))*100),"channel "//getLHEParticle(i1)//" "//getLHEParticle(j1)//" ")
               write(*,"(I3,I3,I7,I7,F16.6,E16.3)") i1,j1, &
               AccepCounter_part(i1,j1), RequEvents(i1,j1), &
               dble(AccepCounter_part(i1,j1))/dble(RejeCounter_part(i1,j1))*100d0, &
               CrossSecMax(i1,j1)
            endif
        enddo
        enddo
        write(*,"(A,I10,I10,F16.6)") "PS gen eff. ",DebugCounter(10),DebugCounter(9),dble(DebugCounter(10))/dble(DebugCounter(9))*100d0
        DebugCounter(9:10)=0
        StatusPercent = int(100d0*dble(sum(AccepCounter_part(:,:)))/dble(sum(RequEvents(:,:)))  )
        print *, "StatusPercent=",StatusPercent, "  Events=",sum(AccepCounter_part(:,:))
        call cpu_time(time_end)
        write(io_stdout,*)  "Event generation rate (events,events/sec)",sum(AccepCounter_part(:,:))-PreviousSum,dble(sum(AccepCounter_part(:,:))-PreviousSum)/(time_end-time_start+1d-10)
        PreviousSum = sum(AccepCounter_part(:,:))
    enddo

! enddo



    print *, " Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter+1d-10) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90."
        write(io_stdout, *) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_stdout, *) "       Increase CSMAX in main.F90."
    endif
    write(io_stdout,*)  " event generation rate (events/sec)",dble(sum(AccepCounter_part(:,:)))/(time_end-time_start+1d-10)




ELSEIF( Process.ge.66 .and. Process.le.69 ) THEN! special treatment for offshell VBF

    if (Process.ge.66 .and. Process.lt.69) then
       VBFoffsh_Hash_Size = Hash_MCFM_qqVVqq_Size
       VBFoffsh_run_size = VBFoffsh_run_qqVVqq_size
    else
       VBFoffsh_Hash_Size = Hash_MCFM_qqVVqqStrong_Size
       VBFoffsh_run_size = VBFoffsh_run_qqVVqqStrong_size
    endif

    !write(6,*) "VBFoffsh_Hash_Size, VBFoffsh_run_size", VBFoffsh_Hash_Size, VBFoffsh_run_size

    if( .not.(VBFoffsh_run.ge.1 .and. VBFoffsh_run.le.VBFoffsh_run_size) ) call Error("Please specify VBFoffsh_run")
    print *, "VBFoffsh_run = ",VBFoffsh_run

    write(io_stdout,"(2X,A)")  "Scanning the integrand"
    warmup = .true.

    if( ReadCSmax ) then
        i=len(trim(CSmaxFile))
        open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'1_gridinfo.txt',form='formatted',status='old',iostat=ios)
        read(io_TmpFile,fmt=*) calls1_in(1)
        read(io_TmpFile,fmt=*) CrossSec2_in(1,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) CrossSecMax2_in(1,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) VG_Result_in(1)
        read(io_TmpFile,fmt=*) VG_Error_in(1)
        read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(1), CrossSectionWithWeightsErrorSquared_in(1)
        close(unit=io_TmpFile)
        if( ios.eq.0 ) print *, "read ",trim(CSmaxFile(1:i-1))//'1_gridinfo.txt'

        open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'2_gridinfo.txt',form='formatted',status='old',iostat=ios)
        read(io_TmpFile,fmt=*) calls1_in(2)
        read(io_TmpFile,fmt=*) CrossSec2_in(2,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) CrossSecMax2_in(2,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) VG_Result_in(2)
        read(io_TmpFile,fmt=*) VG_Error_in(2)
        read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(2), CrossSectionWithWeightsErrorSquared_in(2)
        close(unit=io_TmpFile)
        if( ios.eq.0 ) print *, "read ",trim(CSmaxFile(1:i-1))//'2_gridinfo.txt'

        open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'3_gridinfo.txt',form='formatted',status='old',iostat=ios)
        read(io_TmpFile,fmt=*) calls1_in(3)
        read(io_TmpFile,fmt=*) CrossSec2_in(3,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) CrossSecMax2_in(3,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) VG_Result_in(3)
        read(io_TmpFile,fmt=*) VG_Error_in(3)
        read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(3), CrossSectionWithWeightsErrorSquared_in(3)
        close(unit=io_TmpFile)
        if( ios.eq.0 ) print *, "read ",trim(CSmaxFile(1:i-1))//'3_gridinfo.txt'

        open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'4_gridinfo.txt',form='formatted',status='old',iostat=ios)
        read(io_TmpFile,fmt=*) calls1_in(4)
        read(io_TmpFile,fmt=*) CrossSec2_in(4,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) CrossSecMax2_in(4,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) VG_Result_in(4)
        read(io_TmpFile,fmt=*) VG_Error_in(4)
        read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(4), CrossSectionWithWeightsErrorSquared_in(4)
        close(unit=io_TmpFile)
        if( ios.eq.0 ) print *, "read ",trim(CSmaxFile(1:i-1))//'4_gridinfo.txt'

        open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'5_gridinfo.txt',form='formatted',status='old',iostat=ios)
        read(io_TmpFile,fmt=*) calls1_in(5)
        read(io_TmpFile,fmt=*) CrossSec2_in(5,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) CrossSecMax2_in(5,1:VBFoffsh_Hash_Size)
        read(io_TmpFile,fmt=*) VG_Result_in(5)
        read(io_TmpFile,fmt=*) VG_Error_in(5)
        read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(5), CrossSectionWithWeightsErrorSquared_in(5)
        close(unit=io_TmpFile)
        if( ios.eq.0 ) print *, "read ",trim(CSmaxFile(1:i-1))//'5_gridinfo.txt'

        if (Process .eq. 69) then
           open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'6_gridinfo.txt',form='formatted',status='old',iostat=ios)
           read(io_TmpFile,fmt=*) calls1_in(6)
           read(io_TmpFile,fmt=*) CrossSec2_in(6,1:VBFoffsh_Hash_Size)
           read(io_TmpFile,fmt=*) CrossSecMax2_in(6,1:VBFoffsh_Hash_Size)
           read(io_TmpFile,fmt=*) VG_Result_in(6)
           read(io_TmpFile,fmt=*) VG_Error_in(6)
           read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(6), CrossSectionWithWeightsErrorSquared_in(6)
           close(unit=io_TmpFile)
           if( ios.eq.0 ) print *, "read ",trim(DataFile(1:i-1))//'6_gridinfo.txt'

           open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'7_gridinfo.txt',form='formatted',status='old',iostat=ios)
           read(io_TmpFile,fmt=*) calls1_in(7)
           read(io_TmpFile,fmt=*) CrossSec2_in(7,1:VBFoffsh_Hash_Size)
           read(io_TmpFile,fmt=*) CrossSecMax2_in(7,1:VBFoffsh_Hash_Size)
           read(io_TmpFile,fmt=*) VG_Result_in(7)
           read(io_TmpFile,fmt=*) VG_Error_in(7)
           read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(7), CrossSectionWithWeightsErrorSquared_in(7)
           close(unit=io_TmpFile)
           if( ios.eq.0 ) print *, "read ",trim(DataFile(1:i-1))//'7_gridinfo.txt'

           open(unit=io_TmpFile,file=trim(CSmaxFile(1:i-1))//'8_gridinfo.txt',form='formatted',status='old',iostat=ios)
           read(io_TmpFile,fmt=*) calls1_in(8)
           read(io_TmpFile,fmt=*) CrossSec2_in(8,1:VBFoffsh_Hash_Size)
           read(io_TmpFile,fmt=*) CrossSecMax2_in(8,1:VBFoffsh_Hash_Size)
           read(io_TmpFile,fmt=*) VG_Result_in(8)
           read(io_TmpFile,fmt=*) VG_Error_in(8)
           read(io_TmpFile,fmt=*) CrossSectionWithWeights_in(8), CrossSectionWithWeightsErrorSquared_in(8)
           close(unit=io_TmpFile)
           if( ios.eq.0 ) print *, "read ",trim(DataFile(1:i-1))//'8_gridinfo.txt'
        endif

        !write(6,*) "calls1_in:",calls1_in
        do j=2,VBFoffsh_run_size
           if( calls1_in(1).ne.calls1_in(j) ) call Error("Mismatch in calls1")
        enddo
        calls1 = calls1_in(1)

        CrossSec2(:) = -1d0
        do i=1,VBFoffsh_Hash_Size
            do j=1,VBFoffsh_run_size
               if( CrossSec2_in(j,i).ne.0d0 .and. CrossSec2(i).eq.-1d0 ) CrossSec2(i) = CrossSec2_in(j,i)
            enddo
            if( CrossSec2(i).lt.0d0 ) then
                print *, "WARNING: CrossSec2(i) has not been filled",i,CrossSec2_in(1:VBFoffsh_run_size,i)
                CrossSec2(i) = 0d0
            !else
            !    print *, "Filled CrossSec2(i) with",i,CrossSec2_in(1:VBFoffsh_run_size,i)
            endif
        enddo
        if (VBFoffsh_Hash_Size .lt. NMAXCHANNELS) CrossSec2(VBFoffsh_Hash_Size+1:NMAXCHANNELS)=0d0


        CrossSecMax2(:) = -1d0
        do i=1,VBFoffsh_Hash_Size
            CrossSecMax2(i) = CrossSecMax2_in( max(1,VBFoffsh_run),i)
        enddo
        if (VBFoffsh_Hash_Size .lt. NMAXCHANNELS) CrossSecMax2(VBFoffsh_Hash_Size+1:NMAXCHANNELS)=0d0

        VG_Result = 0d0; VG_Error = 0d0
        CrossSectionWithWeights = 0d0; CrossSectionWithWeightsErrorSquared = 0d0
        do j=1,VBFoffsh_run_size
           VG_Result = VG_Result_in(j) + VG_Result
           VG_Error  = VG_Error_in(j)**2 + VG_Error
           CrossSectionWithWeights = CrossSectionWithWeights_in(j) + CrossSectionWithWeights
           CrossSectionWithWeightsErrorSquared = CrossSectionWithWeightsErrorSquared_in(j) + CrossSectionWithWeightsErrorSquared
        enddo
        VG_Error = dsqrt(VG_Error)

        !write(6,*) "CrossSec2:",CrossSec2
        !write(6,*) "CrossSecMax2:",CrossSecMax2
        !write(6,*) "CrossSectionWithWeights,CrossSectionWithWeightsErrorSquared:",CrossSectionWithWeights,CrossSectionWithWeightsErrorSquared
        !write(6,*) "VG_Result,VG_Error:",VG_Result,VG_Error
        !pause


    else
        itmx = 10
        ncall= VegasNc0
        readin=.false.
        writeout=.true.
        if( Process.ge.66 .and. Process.le.69 ) call vegas(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)

        itmx = 2
        writeout=.true.
        outgridfile = trim(CSmaxFile)//'_step2.grid'
        if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)


        CrossSecMax2(:) = 0d0
        CrossSec2(:) = 0d0
        print *, "resetting CrossSecMax2(:)"
        itmx = 1
        call ClearHisto()
        FindCrossSectionWithWeights = .true.
        CrossSectionWithWeights = 0d0
        CrossSectionWithWeightsErrorSquared = 0d0
        if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        FindCrossSectionWithWeights = .false.
        writeout=.false.
        ingridfile=trim(outgridfile)

        call vegas_get_calls(calls1)
        CrossSec2(:) = CrossSec2(:)/dble(itmx)

        open(unit=io_TmpFile,file=trim(CSmaxFile)//'_gridinfo.txt',form='formatted',status='replace')
        write(io_TmpFile,fmt=*) calls1
        write(io_TmpFile,fmt=*) CrossSec2(1:VBFoffsh_Hash_Size)
        write(io_TmpFile,fmt=*) CrossSecMax2(1:VBFoffsh_Hash_Size)
        write(io_TmpFile,fmt=*) VG_Result
        write(io_TmpFile,fmt=*) VG_Error
        write(io_TmpFile,fmt=*) CrossSectionWithWeights, CrossSectionWithWeightsErrorSquared
        close(unit=io_TmpFile)
    endif

    write(io_stdout,"(A)")  ""
    write(io_stdout,*) "Total unweighted xsec (used by Vegas): ", VG_Result, " +/-", VG_Error, " fb    vs.",sum(CrossSec2(:))
    write(io_stdout,*) "Total xsec with weights (use for physics): ", CrossSectionWithWeights, " +/-", sqrt(CrossSectionWithWeightsErrorSquared)
    call InitOutput(CrossSectionWithWeights, sqrt(CrossSectionWithWeightsErrorSquared))

    RequEvents2(:) = 0
    call HouseOfRepresentatives2(CrossSec2, RequEvents2, VegasNc2)

    if (Process.eq.69) then
       call getRef_MCFM_qqVVqqStrong_Hash(ijSel)
    else
       call getRef_MCFM_qqVVqq_Hash(ijSel)
    endif
    RequEvents(:,:) = 0d0
    do i=1,VBFoffsh_Hash_Size
         i1 = convertToPartIndex(ijSel(i,1))
         j1 = convertToPartIndex(ijSel(i,2))
         if( RequEvents2(i).ne.0 ) write(*,"(I4,I4,I4,F18.8,I10,I10)") i, i1,j1,CrossSec2(i),RequEvents2(i)
    enddo


    write(io_stdout,"(2X,A,F18.3,I10)") "Sum        partonic xsec:",sum(CrossSec2(:))/VG_Result,sum(RequEvents2(:))


    if( .not. ReadCSmax ) then! For ReadCSmax=.false. the program ends here
       if( VBFoffsh_run.ge.1 .and. VBFoffsh_run.le.VBFoffsh_run_size ) print *, "WARNING: These are not the final number of events because the total cross section needs to be assembled from VBFoffsh_run=1,2,3,4,5 (,6,7,8 for process=69)"
       RETURN
    endif

    if (Process .eq. 69) then
       ! Hashes go as 1-25, 26-50, 51-100, 101-140, 141-150, 151-160, 161-170, 171-175
       if( VBFoffsh_run.eq.1 ) then ! removing the requested events for the wrong VBFoffsh_run
            RequEvents2(26:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 25
       elseif( VBFoffsh_run.eq.2 ) then
            RequEvents2(1:25) = 0
            RequEvents2(51:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 25
       elseif( VBFoffsh_run.eq.3 ) then
            RequEvents2(1:50) = 0
            RequEvents2(101:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 50
       elseif( VBFoffsh_run.eq.4 ) then
            RequEvents2(1:100) = 0
            RequEvents2(141:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 40
       elseif( VBFoffsh_run.eq.5 ) then
            RequEvents2(1:140) = 0
            RequEvents2(151:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 10
       elseif( VBFoffsh_run.eq.6 ) then
            RequEvents2(1:150) = 0
            RequEvents2(161:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 10
       elseif( VBFoffsh_run.eq.7 ) then
            RequEvents2(1:160) = 0
            RequEvents2(171:Hash_MCFM_qqVVqqStrong_Size) = 0
            NumPartonicChannels = 10
       elseif( VBFoffsh_run.eq.8 ) then
            RequEvents2(1:170) = 0
            NumPartonicChannels = 5
       endif
    else
       if( VBFoffsh_run.eq.1 ) then ! removing the requested events for the wrong VBFoffsh_run
            RequEvents2(3:Hash_MCFM_qqVVqq_Size) = 0
            NumPartonicChannels = 2
       elseif( VBFoffsh_run.eq.2 ) then
            RequEvents2(1:2) = 0
            RequEvents2(10:Hash_MCFM_qqVVqq_Size) = 0
            NumPartonicChannels = 7
       elseif( VBFoffsh_run.eq.3 ) then
            RequEvents2(1:9) = 0
            RequEvents2(41:Hash_MCFM_qqVVqq_Size) = 0
            NumPartonicChannels = 31
       elseif( VBFoffsh_run.eq.4 ) then
            RequEvents2(1:40) = 0
            RequEvents2(104:Hash_MCFM_qqVVqq_Size) = 0
            NumPartonicChannels = 63
       elseif( VBFoffsh_run.eq.5 ) then
            RequEvents2(1:103) = 0
            NumPartonicChannels = 61
       endif
    endif
    ingridfile = trim(CSmaxFile)//'_step2.grid'

    !write(6,*) "NumPartonicChannels | RequEvents2",NumPartonicChannels," | ",RequEvents2


!     print *, "New sorted hash for VBFoffsh_run=",VBFoffsh_run
!     do i=1,VBFoffsh_Hash_Size
!          i1 = convertToPartIndex(ijSel(i,1))
!          j1 = convertToPartIndex(ijSel(i,2))
!          if( RequEvents2(i).ne.0 ) write(*,"(I4,I4,I4,I4,F18.8,I10,I10)") i,SortedHash(i), i1,j1,CrossSec2(i)/VG_Result,RequEvents2(i)
! !         write(*,"(I4,I4,I4,I4,F18.8,I10,I10)") i, SortedHash(i), i1,j1,CrossSec2(i)/VG_Result,RequEvents2(i)
!     enddo


!--------------------------------------------------------------------
    write(io_stdout,"(A)")  ""
    write(io_stdout,"(2X,A)")  "Event generation"

    call ClearHisto()
    warmup = .false.
    evtgen = .true.
!     nprn = 0
    EvalCounter = 0
    Br_counter(:,:) = 0
    RejeCounter = 0
    AlertCounter = 0
    AccepCounter_part2(:) = 0
    StatusPercent = 0d0

    CrossSecMax2(:) = 1.5d0 * CrossSecMax2(:)    !  adjustment factor
    call cpu_time(time_start)



    itmx=1
    ncall= VegasNc0    !/10     dmax inside vegas needs to be adapted for this. or at least thisdmax inside mod_Crosssection
    nprn =-1 !0
    writeout=.false.

    call vegas_get_calls(calls2)
    calls_rescale = calls1/calls2
    CrossSecMax2(:) = CrossSecMax2(:) * calls_rescale
    print *, "Rescale CrossSecMax2 by ",calls_rescale

    PreviousSum = 0
    if( sum(RequEvents2(:)).eq.0 ) StatusPercent = 100d0
    do while( StatusPercent.lt.100d0  )
        call cpu_time(time_start)
        readin=.true.  ! this prevents adapting the grid during this while-loop

        if( Process.ge.66 .and. Process.le.69 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        write(io_stdout,*) ""

        do i1=1,VBFoffsh_Hash_Size
            if( RequEvents2(i1).gt.0 .and. AccepCounter_part2(i1).lt.RequEvents2(i1)  ) then
               write(*,"(I3,I7,I7,E16.3)") i1, &
               AccepCounter_part2(i1), RequEvents2(i1), &
               CrossSecMax2(i1)
            endif
        enddo

        StatusPercent = int(100d0*dble(sum(AccepCounter_part2(:)))/dble(sum(RequEvents2(:)))  )
        print *, "StatusPercent=",StatusPercent, "  Events=",sum(AccepCounter_part2(:))
        call cpu_time(time_end)
        write(io_stdout,*)  "Event generation rate (events,events/sec)",sum(AccepCounter_part2(:))-PreviousSum,dble(sum(AccepCounter_part2(:))-PreviousSum)/(time_end-time_start+1d-10)
        PreviousSum = sum(AccepCounter_part2(:))
    enddo





    print *, " Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter+1d-10) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90."
        write(io_stdout, *) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_stdout, *) "       Increase CSMAX in main.F90."
    endif
    write(io_stdout,*)  " event generation rate (events/sec)",dble(sum(AccepCounter_part(:,:)))/(time_end-time_start+1d-10)




    !pause


ENDIF


else! not beta version




!-------------------old stuff -------------------
    VG(:,:) = zero
    CSmax(:,:) = zero

    if( .not. ReadCSmax ) then
        print *, " finding maximal weight with ",VegasNc0," evaluations"
        do i=1,VegasNc0
                call random_number(yRnd)
                RES(:,:) = 0d0
                if( Process.eq.80 ) then
                    dum = EvalUnWeighted_TTBH(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.90 ) then
                    dum = EvalUnWeighted_BBBH(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.60 ) then
                    dum = EvalUnWeighted_HJJ(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.ge.66 .and. Process.le.69 ) then
                    dum = EvalUnWeighted_HJJ_fulldecay(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.61 ) then
                    dum = EvalUnWeighted_HJJ(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.110 ) then
                    dum = EvalUnWeighted_TH(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.111 ) then
                    dum = EvalUnWeighted_TH(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.112 ) then
                    dum = EvalUnWeighted_TH(yRnd,.false.,(/-99,-99/),RES)
                elseif( Process.eq.113 ) then
                    dum = EvalUnWeighted_TH(yRnd,.false.,(/-99,-99/),RES)
                endif
                VG(:,:) = VG(:,:) + RES(:,:)
                PChannel = PChannel_aux
        enddo
        open(unit=io_CSmaxFile,file=trim(CSmaxFile)//'_CSmax.bin',form='unformatted',status='replace')
        WRITE(io_CSmaxFile) CSMAX,VG
        close(io_CSmaxFile)
    else
        open(unit=io_CSmaxFile,file=trim(CSmaxFile)//'_CSmax.bin',form='unformatted')
        READ(io_CSmaxFile) CSMAX,VG
        close(io_CSmaxFile)
    endif

   CSmax(:,:)   = 1.5d0 * CSmax(:,:)    !  adjustment factor

   VG(:,:) = VG(:,:)/dble(VegasNc0)
   TotalXSec = sum(  VG(:,:) )
   print *, ""
   write(io_stdout,"(1X,A,F10.3)") "Total xsec: ",TotalXSec


    RequEvents(:,:) = 0
    if (VegasNc2.ne.-1) then
      call HouseOfRepresentatives(VG(:,:), RequEvents(-5:5,-5:5), VegasNc2)
    endif

    do i1=-5,5
    do j1=-5,5
         if( VG(i1,j1).gt.1d-9 ) write(io_stdout,"(1X,A,I4,I4,F8.3,I9)") "Fractional partonic xsec ", i1,j1,VG(i1,j1)/TotalXSec,RequEvents(i1,j1)
    enddo
    enddo





!------------------------------- set counts to zero for actual evaluation

    EvalCounter = 0
    AccepCounter = 0
    RejeCounter = 0
    AlertCounter = 0
    AccepCounter_part(:,:) = 0

    call cpu_time(time_start)


    do i1=-5,5!! idea: instead of these 2 do-loop introduce randomized loop
    do j1=-5,5
    if(RequEvents(i1,j1).ne.0) then
          if( (PChannel.eq.0) .and. (abs(i1)+abs(j1).ne.0) ) cycle
          if( (PChannel.eq.1) .and. i1*j1.eq.0 ) cycle

          write(io_stdout,*) ""
          write(io_stdout,*) ""
          write(io_stdout,"(X,A,I8,A,I4,I4,A)",advance='no') "generating ",RequEvents(i1,j1)," events for channel",i1,j1,":  "
          flush(io_stdout)

          do while( AccepCounter_part(i1,j1)  .lt. RequEvents(i1,j1) )
              call random_number(yRnd)
              if( Process.eq.80 ) then
                  dum = EvalUnWeighted_TTBH(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.90 ) then
                  dum = EvalUnWeighted_BBBH(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.60 ) then
                  dum = EvalUnWeighted_HJJ(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.ge.66 .and. Process.le.69 ) then
                  dum = EvalUnWeighted_HJJ_fulldecay(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.61 ) then
                  dum = EvalUnWeighted_HJJ(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.110 ) then
                  dum = EvalUnWeighted_TH(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.111 ) then
                  dum = EvalUnWeighted_TH(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.112 ) then
                  dum = EvalUnWeighted_TH(yRnd,.true.,(/i1,j1/),RES)
              elseif( Process.eq.113 ) then
                  dum = EvalUnWeighted_TH(yRnd,.true.,(/i1,j1/),RES)
              endif
              StatusPercent = int(100d0*(AccepCounter_part(i1,j1))  /  dble(RequEvents(i1,j1))  )
              call PrintStatusBar( StatusPercent )
          enddo

    endif
    enddo
    enddo




    call cpu_time(time_end)
    print *, ""
    print *, ""
    print *, " Evaluation Counter: ",EvalCounter
    print *, " Acceptance Counter: ",AccepCounter
    do i1=-5,+5
    do j1=-5,+5
      if( AccepCounter_part(i1,j1).ne.0 ) print *, " Acceptance  Counter_part: ", i1,j1, AccepCounter_part(i1,j1)
    enddo
    enddo
    print *, " Alert  Counter: ",AlertCounter
    print *, " gg/qqb ratio = ", dble(AccepCounter_part(0,0))/dble(sum(AccepCounter_part(:,:))-AccepCounter_part(0,0)+1d-10)
    if( dble(AlertCounter)/dble(AccepCounter+1d-10) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90."
        write(io_stdout, *) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_stdout, *) "       Increase CSMAX in main.F90."
    endif
    write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start+1d-10)


endif! beta version

  endif! unweighted




return
END SUBROUTINE




SUBROUTINE ReadMassWidth()
use ModParameters
use ModMisc
implicit none
logical :: FirstEvent,InMadgraphMassBlock
character(len=160) :: FirstLines
integer :: i, j, stat
character(len=150) :: InputFmt0,InputFmt1
integer :: nline, tries
integer :: EventNumPart, OtherInts(1:5), LHE_IDUP(1:maxpart)
real(8) :: Mass(1:maxpart), OtherReal8s(1:4)
character(len=160) :: OtherLines
character(len=160) :: EventLine(0:maxpart)

!    search for line with first event
     FirstEvent = .false.
     FoundHiggsMass=.false.
     FoundHiggsWidth=.false.
     InMadgraphMassBlock=.false.
     do while ( .not.FirstEvent )
        read(io_LHEInFile,fmt="(A160)",IOSTAT=stat,END=99) FirstLines

        !Read the Higgs mass
        !JHUGen
        i = Index(FirstLines, "Resonance: spin=")
        if ( i .ne. 0) then
            do while (FirstLines(i:i+4) .ne. "mass=")
                i = i+1
            enddo
            i = i+5 !after the =
            do while (FirstLines(i:i).eq." ")
                i = i+1
            enddo
            j = i
            do while (FirstLines(j:j).ne." ")
                j = j+1
            enddo
            read(FirstLines(i:j),*) M_Reso
            M_Reso = M_Reso*GeV
            FoundHiggsMass=.true.
        endif
        !POWHEG
        if( FirstLines(1:5).eq."hmass" ) then
               i = 6
               do while(FirstLines(i:i).eq." ")
                   i = i+1
               enddo
               do while(FirstLines(i:i).ne." ")
                   i = i+1
               enddo
               read(FirstLines(6:i),*) M_Reso
               M_Reso = M_Reso*GeV!  convert to units of 100GeV
               FoundHiggsMass=.true.
        endif
        !Madgraph
        if (Index(Capitalize(FirstLines),"BLOCK MASS").ne.0 .and. .not.FoundHiggsMass) then
               InMadgraphMassBlock=.true.
        elseif (Index(Capitalize(FirstLines),"BLOCK").ne.0) then
               InMadgraphMassBlock=.false.
        endif
        if (InMadgraphMassBlock) then
               i = Index(FirstLines, " 25 ")
               if (i.ne.0) then
                   j = Index(FirstLines(i+4:len(FirstLines)), " ") + i+4 - 1
                   if (j.eq.0) then
                       j = len(FirstLines)
                   endif
                   read(FirstLines(i+4:j),*) M_Reso
                   M_Reso = M_Reso*GeV
                   FoundHiggsMass=.true.
               endif
        endif

        !Read the Higgs width
        !JHUGen
        i = Index(FirstLines, "Resonance: spin=")
        if ( i .ne. 0) then
            do while (FirstLines(i:i+5) .ne. "width=")
                i = i+1
            enddo
            i = i+6 !after the =
            do while (FirstLines(i:i).eq." ")
                i = i+1
            enddo
            j = i
            do while (FirstLines(j:j).ne." ")
                j = j+1
            enddo
            read(FirstLines(i:j),*) Ga_Reso
            Ga_Reso = Ga_Reso*GeV
            FoundHiggsWidth=.true.
        endif
        !POWHEG
        if( FirstLines(1:6).eq."hwidth" ) then
               i = 7
               do while(FirstLines(i:i).eq." ")
                   i = i+1
               enddo
               do while(FirstLines(i:i).ne." ")
                   i = i+1
               enddo
               read(FirstLines(7:i),*) Ga_Reso
               Ga_Reso = Ga_Reso*GeV!  convert to units of 100GeV
               FoundHiggsWidth=.true.
        endif
        !Madgraph
        if (Index(Capitalize(FirstLines),"DECAY ").ne.0) then
               i = Index(FirstLines, " 25 ")
               if (i.ne.0) then
                   j = Index(FirstLines(i+4:len(FirstLines)), " ") + i+4 - 1
                   if (j.eq.0) then
                       j = len(FirstLines)
                   endif
                   read(FirstLines(i+4:j),*) Ga_Reso
                   Ga_Reso = Ga_Reso*GeV
                   FoundHiggsWidth=.true.
               endif
        endif

        !Read the generated mass shape for reweighting, POWHEG only
        if (FirstLines(1:7).eq."bwshape") then
            i = 8
            do while(FirstLines(i:i).eq." ")
                i = i+1
            enddo
            do while(FirstLines(i:i).ne." ")
                i = i+1
            enddo
            read(FirstLines(8:i),*) j
            if( WidthSchemeIn.gt.0 .and. j.ne.WidthSchemeIn ) then
                print *, "WidthSchemeIn is specified to be ", WidthSchemeIn, " but the LHE file says:"
                print *, FirstLines
                stop 1
            elseif( WidthScheme.gt.0 .and. WidthSchemeIn.le.0 .and. j.ne.WidthScheme ) then
                print *, "WidthScheme is specified to be ", WidthScheme, " but the LHE file says:"
                print *, FirstLines
                print "(A,I1,A,I1,A,I1,A,I1)", "If you want to reweight the propagator from ", j, " to ", WidthScheme, " please specify WidthSchemeIn=", j, " WidthScheme=", WidthScheme
                stop 1
            endif
            if( WidthScheme.lt.0 ) read(FirstLines(8:i),*) WidthScheme
            if( WidthSchemeIn.lt.0 ) read(FirstLines(8:i),*) WidthSchemeIn
        endif

        if( Index(FirstLines, "<event").ne.0 ) FirstEvent=.true.
     enddo
99   continue

     if(.not.ReadLHEFile .or. (m_Reso.lt.2*m_V .and. .not.ReweightDecay)) then
         call ReopenInFile()
         maxInputmHstar = m_Reso
         minInputmHstar = m_Reso
         mHstarforphasespace = m_Reso
         return
     endif

     print *, "Finding range of mHstar in the LHE input file..."


     InputFmt0 = ""
     InputFmt1 = ""

     do while ( .true. )
         read(io_LHEInFile,"(A)") EventLine(0)
         if (UseUnformattedRead) then
             read(EventLine(0),*) EventNumPart
         else
             if (InputFmt0.eq."") then
                 InputFmt0 = FindInputFmt0(EventLine(0))
             endif
             read(EventLine(0),InputFmt0) EventNumPart, OtherInts(1), OtherReal8s(1:4)
         endif
!        read event lines
         do nline=1,EventNumPart
            read(io_LHEInFile,fmt="(A160)") EventLine(nline)
         enddo
         if( EventNumPart.lt.3 .or. EventNumPart.gt.maxpart ) then
            call Error("Number of particles in LHE input exceeds allowed limit",EventNumPart)
         endif
         do nline=1,EventNumPart
            if (UseUnformattedRead) then
                read(EventLine(nline),*) LHE_IDUP(nline),OtherInts(1:5),OtherReal8s(1:4),Mass(nline)
            else
                if (InputFmt1.eq."") then
                    InputFmt1 = FindInputFmt1(EventLine(nline))
                endif
                read(EventLine(nline),InputFmt1) LHE_IDUP(nline),OtherInts(1:5),OtherReal8s(1:4),Mass(nline)
            endif
            if( abs(LHE_IDUP(nline)).eq.25 ) then!   select the Higgs (ID=25, h0) ! Ulascan: Should be safer to have 'if( abs(LHE_IDUP(nline)).eq.abs(convertLHE(Hig_)) )'
                  Mass(nline) = Mass(nline)*GeV            !  convert to units of 100GeV
                  if( Mass(nline).gt.maxInputmHstar ) maxInputmHstar = Mass(nline)
                  if( Mass(nline).lt.minInputmHstar ) minInputmHstar = Mass(nline)
                  exit
            endif
         enddo
!        read optional lines
         tries = 0
         do while (.true.)
              tries = tries +1
              read(io_LHEInFile,fmt="(A160)",IOSTAT=stat,END=98) OtherLines(1:160)
              if(OtherLines(1:30).eq."</LesHouchesEvents>") then
                  goto 98
              elseif( Index(OtherLines,"<event").ne.0 ) then
                  exit
              elseif( tries.gt.10000000 ) then
                  write(io_LHEOutFile,"(A)") "</event>"
                  print *, "ERROR: cannot find </event>"
                  exit
              endif
         enddo
     enddo
98   continue
     if(maxinputmHstar.lt.0) then !no events in the file
         mininputmHstar = m_Reso
         maxinputmHstar = m_Reso
     endif
     print *, "... and it's ", minInputmHstar/GeV, "  -  ", maxInputmHstar/GeV, " GeV"

     if(m_Reso.lt.2*m_V) then
         mHstarforphasespace = m_Reso
     else
         mHstarforphasespace = maxinputmHstar
     endif

     call ReopenInFile()

return
END SUBROUTINE

SUBROUTINE InitReadLHE(BeginEventLine)
use ModParameters
use ModMisc
implicit none
logical :: FirstEvent, WroteHeader, IsThisTheFirstInitLine
character(len=160) :: FirstLines
integer :: stat, InitFieldIndex, i, IDWTUPIndex(1:2), IDWTUP
character(len=4) :: IDWTUPFmt
character(len=100), intent(out) :: BeginEventLine

     write(io_LHEOutFile ,'(A)') '<LesHouchesEvents version="1.0">'
     FirstEvent = .false.
     WroteHeader = .false.
     IsThisTheFirstInitLine = .false.
     do while ( .not.FirstEvent )
        read(io_LHEInFile,fmt="(A160)",IOSTAT=stat,END=99) FirstLines
        if ( FirstLines(1:4).eq."<!--" .and. .not.WroteHeader ) then
            call InitOutput(1d0, 1d14)
            WroteHeader = .true.
        endif
        if (index(FirstLines,"<MG").ne.0 .and. .not.WroteHeader) then  !Sometimes MadGraph doesn't have a comment at the beginning
            call InitOutput(1d0, 1d14)                                 !In that case put the JHUGen header before the MadGraph
            write(io_LHEOutFile, "(A)") "-->"                          ! proc card, etc.
            WroteHeader = .true.                                       !and put the Higgs mass/width in a separate comment
        endif
        if( IsThisTheFirstInitLine ) then
          !See https://arxiv.org/pdf/hep-ph/0109068v1.pdf under IDWTUP
          !For some processes, POWHEG writes all events with a weight of 1 (or +/-1)
          !and writes IDWTUP as 3 (-3).  When Pythia sees IDWTUP=3 (-3), it assumes
          !all weights are 1 (+/-1) and ignores the weights (abs(weights)) written
          !in the LHE file. If we want to modify those weights, we have to set
          !IDWTUP to 4 (-4).
          if( ReweightDecay .or. WidthScheme.ne.WidthSchemeIn ) then
            InitFieldIndex = 1
            i = 1
            do while( InitFieldIndex.le.9 )
              do while( FirstLines(i:i).eq." " )
                i=i+1
              enddo
              if( InitFieldIndex.eq.9 ) IDWTUPIndex(1) = i
              do while( FirstLines(i:i).ne." " )
                i=i+1
              enddo
              if( InitFieldIndex.eq.9 ) IDWTUPIndex(2) = i-1
              InitFieldIndex=InitFieldIndex+1
            enddo
            read(FirstLines(IDWTUPIndex(1):IDWTUPIndex(2)), fmt=*) IDWTUP
            if( abs(IDWTUP).eq.3 ) IDWTUP = IDWTUP * 4 / 3

            write(IDWTUPFmt, fmt="(A,I1,A)") "(I", IDWTUPIndex(2)-IDWTUPIndex(1)+1, ")"

            write(FirstLines(IDWTUPIndex(1):IDWTUPIndex(2)), fmt=IDWTUPFmt) IDWTUP
          endif
          IsThisTheFirstInitLine = .false.
        endif
        if (Index(FirstLines,"<init>").ne.0 ) then
            if( .not.WroteHeader ) then !If not now, when?
                call InitOutput(1d0, 1d14)
                WroteHeader = .true.
            endif
            IsThisTheFirstInitLine = .true.
        endif

        if( Index(FirstLines, "<event").ne.0 ) then
            FirstEvent=.true.
            BeginEventLine = trim(FirstLines)
        else
            if( Index(FirstLines,"<LesHouchesEvents").ne.0 .or. Index(FirstLines,"<!--").ne.0 ) then
            else
              write(io_LHEOutFile,"(A)") trim(firstlines)
            endif
        endif

    enddo

    if( .not. FoundHiggsMass ) then
       write(io_stdout,"(2X,A,1F7.2)")  "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
       write(io_LogFile,"(2X,A,1F7.2)") "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
    else
       write(io_stdout,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
       write(io_LogFile,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
    endif
    if( .not. FoundHiggsWidth ) then
       write(io_stdout,"(2X,A,1F10.5)")  "ERROR: Higgs width could not be read from LHE input file. Assuming default value",Ga_Reso*100d0
       write(io_LogFile,"(2X,A,1F10.5)") "ERROR: Higgs width could not be read from LHE input file. Assuming default value",Ga_Reso*100d0
    else
       write(io_stdout,"(2X,A,1F10.5,A)") "A Higgs width of ",Ga_Reso*100d0," GeV was determined from the LHE input file."
       write(io_LogFile,"(2X,A,1F10.5,A)") "A Higgs width of ",Ga_Reso*100d0," GeV was determined from the LHE input file."
    endif
    write(io_stdout,"(A)") ""
    write(io_LogFile,"(A)") ""

99  continue

return
END SUBROUTINE


SUBROUTINE StartReadLHE_NEW(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
use ModPMZZ
implicit none
include 'csmaxvalue.f'
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),Res,EMcheck(1:4),DecayWeight,DecayWidth,DecayWidth0
real(8) :: HiggsDK_Mom(1:4,1:13),Ehat
real(8) :: MomExt(1:4,1:maxpart),MomHiggs(1:4),Mass(1:maxpart),pH2sq
integer :: tries, nParticle,  ICOLUP(1:2,1:7+maxpart),LHE_IntExt(1:7+maxpart),HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13)
character(len=150) :: InputFmt0,InputFmt1
integer :: nline,intDummy,Nevent
integer :: LHE_IDUP(1:maxpart),LHE_ICOLUP(1:2,1:maxpart),LHE_MOTHUP(1:2,1:maxpart)
integer :: EventNumPart, EventProcessId
real(8) :: WeightScaleAqedAqcd(1:4)
character(len=160) :: OtherLines
character(len=160) :: EventLine(0:maxpart)
integer :: n,stat,iHiggs,VegasSeed,AccepLastPrinted
character(len=100) :: BeginEventLine
logical :: Empty


if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2
AccepLastPrinted = 0

call InitReadLHE(BeginEventLine)

     if( ReweightDecay ) then
         print *, " finding P_decay(m4l) distribution with ", PMZZEvals, " calls per point"
         call InitMZZdistribution()
     endif
     DecayWidth0 = GetMZZProbability(M_Reso,-1d0,.true.)

     print *, " finding maximal weight for mZZ=", maxInputmHstar/GeV, " GeV with ",VegasNc0," points"
     VG = zero
     CSmax = zero
     EHat = mHstarforphasespace! for m > 2mV: max mHstar found in the LHE file, which should give the max value of the integrand
                               ! for m < 2mV: m_Reso
     if( TauDecays.lt.0 ) then
         do tries=1,VegasNc0
             call random_number(yRnd)
             DecayWeight = EvalUnWeighted_DecayToVV(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP(1:9),HiggsDK_ICOLUP)
         enddo

     else
         do tries=1,VegasNc0
             call random_number(yRnd)
             DecayWeight = EvalUnWeighted_DecayToTauTau(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,4:13),HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13))
         enddo
     endif
     csmax(0,0)   = 1.5d0*csmax(0,0)    !  safety buffer

     InputFmt0 = ""
     InputFmt1 = ""

     print *, " generating events"
     EvalCounter = 0
     AccepCounter = 0
     RejeCounter = 0
     AccepCounter_part = 0
     call cpu_time(time_start)
     NEvent=0
     do while ( .true. )
         NEvent=NEvent + 1
         LeptInEvent(:) = 0
         JetsInEvent(:) = 0
         read(io_LHEInFile,"(A)") EventLine(0)
         if (UseUnformattedRead) then
             read(EventLine(0),*) EventNumPart, EventProcessId, WeightScaleAqedAqcd
         else
             if (InputFmt0.eq."") then
                 InputFmt0 = FindInputFmt0(EventLine(0))
             endif
             read(EventLine(0),InputFmt0) EventNumPart, EventProcessId, WeightScaleAqedAqcd
         endif
!        read event lines
         do nline=1,EventNumPart
            read(io_LHEInFile,fmt="(A160)") EventLine(nline)
         enddo
         if( EventNumPart.lt.3 .or. EventNumPart.gt.maxpart ) then
            call Error("Number of particles in LHE input exceeds allowed limit",EventNumPart)
         endif
         do nline=1,EventNumPart
            if (UseUnformattedRead) then
                read(EventLine(nline),*) LHE_IDUP(nline),LHE_IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline)
            else
                if (InputFmt1.eq."") then
                    InputFmt1 = FindInputFmt1(EventLine(nline))
                endif
                read(EventLine(nline),InputFmt1) LHE_IDUP(nline),LHE_IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline)
            endif
            MomExt(1:4,nline) = MomExt(1:4,nline)*GeV!  convert to units of 100GeV
            Mass(nline) = Mass(nline)*GeV            !  convert to units of 100GeV
            if( abs(LHE_IDUP(nline)).eq.25 ) then!   select the Higgs (ID=25, h0)
                  MomHiggs(1:4) = MomExt(1:4,nline)
                  pH2sq = dsqrt(abs(MomHiggs(1:4).dot.MomHiggs(1:4)))
                  iHiggs = nline
            endif
            if( IsALHELepton(LHE_IDUP(nline)) .and. LHE_IntExt(nline).eq.1 ) then
                  LeptInEvent(0) = LeptInEvent(0) + 1
                  LeptInEvent( LeptInEvent(0) ) = LHE_IDUP(nline)
            endif
            if( IsALHEJet(LHE_IDUP(nline)) .and. LHE_IntExt(nline).eq.1 ) then
                  JetsInEvent(0) = JetsInEvent(0) + 1
                  JetsInEvent( JetsInEvent(0) ) = LHE_IDUP(nline)
            endif
         enddo


!         accept/reject sampling for H->VV decay contribution
          EHat = pH2sq
          DecayWeight = 0d0

          DecayWidth = GetMZZProbability(EHat,-1d0,.true.)!  could also be used to determine csmax for this particular event to improve efficiency (-->future work)
          WeightScaleAqedAqcd(1) = WeightScaleAqedAqcd(1) * DecayWidth/DecayWidth0 ! This line does nothing unless either WidthScheme.ne.WidthSchemeIn or ReweightDecay

          if( TauDecays.lt.0 ) then
                do tries=1,5000000
                    call random_number(yRnd)
                    DecayWeight = EvalUnWeighted_DecayToVV(yRnd,.true.,EHat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP(1:9),HiggsDK_ICOLUP(1:2,1:9))
                    if( Res.ne.0d0 ) exit
                enddo
          else
                do tries=1,5000000
                    call random_number(yRnd)
                    DecayWeight =  EvalUnWeighted_DecayToTauTau(yRnd,.true.,EHat,Res,HiggsDK_Mom(1:4,1:13),HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13))
                    if( Res.ne.0d0 ) exit
                enddo
          endif

          Empty = .false.
          if( Res.le.0 .and. WriteFailedEvents.ne.0 ) then
              if( WriteFailedEvents.eq.1 ) then
                  WeightScaleAqedAqcd(1) = 0d0! events that were not accepted after 50 Mio. tries are assigned weight zero
                  Res = 1d0
              elseif( WriteFailedEvents.eq.2 ) then
                  WeightScaleAqedAqcd(1) = 0d0
                  Res = 1d0
                  Empty = .true.
              endif
          endif

          if( Res.gt.0d0 ) then ! decay event was accepted
             if( TauDecays.lt.0 ) then!  H->VV->4f
                call boost(HiggsDK_Mom(1:4,6),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,7),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,8),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,9),MomHiggs(1:4),pH2sq)
                HiggsDK_Mom(1:4,4) = HiggsDK_Mom(1:4,6) + HiggsDK_Mom(1:4,7)
                HiggsDK_Mom(1:4,5) = HiggsDK_Mom(1:4,8) + HiggsDK_Mom(1:4,9)
                HiggsDK_IDUP(4) = convertLHE(HiggsDK_IDUP(4))
                HiggsDK_IDUP(5) = convertLHE(HiggsDK_IDUP(5))
                HiggsDK_IDUP(6) = convertLHE(HiggsDK_IDUP(6))
                HiggsDK_IDUP(7) = convertLHE(HiggsDK_IDUP(7))
                HiggsDK_IDUP(8) = convertLHE(HiggsDK_IDUP(8))
                HiggsDK_IDUP(9) = convertLHE(HiggsDK_IDUP(9))
                if (UseUnformattedRead) then
                    call WriteOutEvent_NEW(EventNumPart,LHE_IDUP,LHE_IntExt,LHE_MOTHUP,LHE_ICOLUP,MomExt,HiggsDK_Mom(1:4,4:9),Mass,iHiggs,HiggsDK_IDUP,HiggsDK_ICOLUP,EventProcessId,EventWeight=WeightScaleAqedAqcd(1),EventScaleAqedAqcd=WeightScaleAqedAqcd(2:4),BeginEventLine=BeginEventLine,Empty=Empty)
                else
                    call WriteOutEvent_NEW(EventNumPart,LHE_IDUP,LHE_IntExt,LHE_MOTHUP,LHE_ICOLUP,MomExt,HiggsDK_Mom(1:4,4:9),Mass,iHiggs,HiggsDK_IDUP,HiggsDK_ICOLUP,EventProcessId,EventWeight=WeightScaleAqedAqcd(1),EventScaleAqedAqcd=WeightScaleAqedAqcd(2:4),BeginEventLine=BeginEventLine,InputFmt0=InputFmt0,Empty=Empty)
                endif
             else! H->tautau
                call boost(HiggsDK_Mom(1:4,4),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,5),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,6),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,7),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,8),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,9),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,10),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,11),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,12),MomHiggs(1:4),pH2sq)
                call boost(HiggsDK_Mom(1:4,13),MomHiggs(1:4),pH2sq)

                HiggsDK_IDUP(4) = convertLHE(HiggsDK_IDUP(4))
                HiggsDK_IDUP(5) = convertLHE(HiggsDK_IDUP(5))
                HiggsDK_IDUP(6) = convertLHE(HiggsDK_IDUP(6))
                HiggsDK_IDUP(7) = convertLHE(HiggsDK_IDUP(7))
                HiggsDK_IDUP(8) = convertLHE(HiggsDK_IDUP(8))
                HiggsDK_IDUP(9) = convertLHE(HiggsDK_IDUP(9))
                HiggsDK_IDUP(10)= convertLHE(HiggsDK_IDUP(10))
                HiggsDK_IDUP(11)= convertLHE(HiggsDK_IDUP(11))
                HiggsDK_IDUP(12)= convertLHE(HiggsDK_IDUP(12))
                HiggsDK_IDUP(13)= convertLHE(HiggsDK_IDUP(13))
                if (UseUnformattedRead) then
                    call WriteOutEvent_HFF(EventNumPart,LHE_IDUP,LHE_IntExt,LHE_MOTHUP,LHE_ICOLUP,MomExt,HiggsDK_Mom(1:4,1:13),Mass,iHiggs,HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13),EventProcessId,EventWeight=WeightScaleAqedAqcd(1),EventScaleAqedAqcd=WeightScaleAqedAqcd(2:4),BeginEventLine=BeginEventLine,Empty=Empty)
                else
                    call WriteOutEvent_HFF(EventNumPart,LHE_IDUP,LHE_IntExt,LHE_MOTHUP,LHE_ICOLUP,MomExt,HiggsDK_Mom(1:4,1:13),Mass,iHiggs,HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13),EventProcessId,EventWeight=WeightScaleAqedAqcd(1),EventScaleAqedAqcd=WeightScaleAqedAqcd(2:4),BeginEventLine=BeginEventLine,InputFmt0=InputFmt0,Empty=Empty)
                endif
             endif

             if( mod(NEvent,5000).eq.0 .and. AccepCounter.ne.AccepLastPrinted) then
                  call cpu_time(time_int)
                  write(io_stdout,*)  NEvent," events processed (",time_int-time_start, ") seconds"
                  write(io_LogFile,*) NEvent," events processed (",time_int-time_start, ") seconds"
                  AccepLastPrinted = AccepCounter
             endif

          elseif( Res.eq.0d0 ) then ! decay event was not accepted after ncall evaluations, read next production event
             print *, "Rejected event after ",tries-1," evaluations"
             AlertCounter = AlertCounter + 1
          endif




!        read optional lines
         tries = 0
         do while (.true.)
              tries = tries +1
              read(io_LHEInFile,fmt="(A160)",IOSTAT=stat,END=99) OtherLines(1:160)
              if(OtherLines(1:30).eq."</LesHouchesEvents>") then
                  if( RequestNLeptons.gt.0 .or. RequestNJets.gt.0 ) then
                    write(io_LHEOutFile,"(A)") "<!-- Filter information:"
                    write(io_LHEOutFile,"(A,I8)") "     events processed:  ", NEvent
                    write(io_LHEOutFile,"(A,I8)") "     events accepted:   ", AccepCounter
                    write(io_LHEOutFile,"(A,1F6.2,A)") "     filter efficiency: ", dble(AccepCounter)/dble(NEvent)*100d0,"% -->"
                  endif
                  goto 99
              elseif( (Index(OtherLines(1:7),"</event").ne.0) .and. Res.gt.0d0 ) then
                  write(io_LHEOutFile,"(A)") trim(OtherLines)
              elseif( Index(OtherLines,"<event").ne.0 ) then
                  BeginEventLine = trim(OtherLines)
                  exit
              elseif( Res.gt.0d0 ) then !if there are "#" comments
                  write(io_LHEOutFile,fmt="(A)") trim(OtherLines)
              elseif( tries.gt.10000000 ) then
                  write(io_LHEOutFile,"(A)") "</event>"
                  print *, "ERROR: cannot find </event>"
                  exit
              endif
         enddo

     enddo
99   continue
     call cpu_time(time_end)




    write(io_stdout,*) ""
    write(io_stdout,*) "Evaluation Counter: ",EvalCounter
    write(io_stdout,*) "Acceptance Counter: ",AccepCounter
    write(io_stdout,*) "Rejection  Counter: ",RejeCounter
    write(io_stdout,*) "Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        write(io_stdout,*) "ALERT: The number of rejected events exceeds 1%."
        write(io_stdout,*) "       Increase CSMAX in main.F90 or VegasNc0."
    endif
   write(io_stdout,*)  "Event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)
   if( RequestNLeptons.gt.0 .or. RequestNJets.gt.0 ) write(io_stdout,"(A,1F6.2,A)") " Filter efficiency:",dble(AccepCounter)/dble(NEvent)*100d0," %"


    write(io_LogFile,*) ""
    write(io_LogFile,*) "Evaluation Counter: ",EvalCounter
    write(io_LogFile,*) "Acceptance Counter: ",AccepCounter
    write(io_LogFile,*) "Rejection  Counter: ",RejeCounter
    write(io_LogFile,*) "Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90 or VegasNc0."
    endif
   write(io_LogFile,*)  "Event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)
   if( RequestNLeptons.gt.0 .or. RequestNJets.gt.0 ) write(io_LogFile,"(A,1F6.2,A)") " Filter efficiency:",dble(AccepCounter)/dble(NEvent)*100d0," %"



return
END SUBROUTINE




SUBROUTINE StartConvertLHE(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
use ModMisc
implicit none
include 'csmaxvalue.f'
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),Res,dum,EMcheck(1:4),xRnd
real(8) :: AcceptedEvent(1:4,1:maxpart),Ehat,pH2sq
real(8) :: MomExt(1:4,1:maxpart),MomShift(1:4,1:maxpart),MomHiggs(1:4),MomParton(1:4,1:maxpart),Mass(1:maxpart),Lifetime(1:maxpart),Spin(1:maxpart)
integer :: tries, nParticle, MY_IDUP(1:7+maxpart), ICOLUP(1:2,1:7+maxpart),IntExt(1:7+maxpart),convertparent
integer :: nline,intDummy,Nevent
integer :: LHE_IDUP(1:maxpart+3),   LHE_ICOLUP(1:2,1:maxpart+3),   LHE_MOTHUP(1:2,1:maxpart+3)
integer :: LHE_IDUP_Part(1:maxpart),LHE_ICOLUP_Part(1:2,1:maxpart),LHE_MOTHUP_Part(1:2,1:maxpart+3)
integer :: EventNumPart,nparton,EventProcessId
real(8) :: WeightScaleAqedAqcd(1:4)
character(len=120) :: PDFLine
character(len=160) :: EventLine(0:maxpart+3)
integer :: VegasSeed,i,stat,DecayParticles(1:2)
integer, dimension(:), allocatable :: gen_seed
character(len=*),parameter :: DefaultFmt0 = "I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7"
character(len=*),parameter :: Fmt1 = "6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0"
character(len=150) :: IndentedFmt0, IndentedFmt1, InputFmt0, InputFmt1
character(len=100) :: BeginEventLine
integer :: indent



if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2

Lifetime(1:maxpart) = 0d0
Spin(1:maxpart) = 0d0

call InitReadLHE(BeginEventLine)

     InputFmt0 = ""
     InputFmt1 = ""

     print *, " converting events"
     call cpu_time(time_start)
     NEvent=0
     do while ( .true. )
         NEvent=NEvent + 1
         read(io_LHEInFile,"(A)") EventLine(0)
         if (UseUnformattedRead) then
             read(EventLine(0),*) EventNumPart, EventProcessId, WeightScaleAqedAqcd
         else
             if (InputFmt0.eq."") then
                 InputFmt0 = FindInputFmt0(EventLine(0))
             endif
             read(EventLine(0),InputFmt0) EventNumPart, EventProcessId, WeightScaleAqedAqcd
         endif

!        read event lines
         do nline=1,EventNumPart
            read(io_LHEInFile,"(A)") EventLine(nline)
         enddo
         if( EventNumPart.lt.3 .or. EventNumPart.gt.maxpart ) then
            call Error("Number of particles in LHE input exceeds allowed limit",EventNumPart)
         endif

 !       convert event lines into variables assuming that the Higgs resonance has ID 25
         nparton = 0
         do nline=1,EventNumPart
            if (UseUnformattedRead) then
                read(EventLine(nline),*) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline),Lifetime(nline),Spin(nline)
            else
                if (InputFmt1.eq."") then
                    InputFmt1 = FindInputFmt1(EventLine(nline))
                endif
                read(EventLine(nline),InputFmt1) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline),Lifetime(nline),Spin(nline)
            endif
            Spin(1:maxpart) = 0d0  !discard spin information
                                   !would not make sense if converting from a generator that has spin

            if( IntExt(nline).eq.2 .and. (LHE_IDUP(nline).eq.convertLHE(Z0_) .or. LHE_IDUP(nline).eq.convertLHE(Wp_) .or. LHE_IDUP(nline).eq.convertLHE(Wm_)) ) then
               convertparent = nline
            endif
! print *, nline
! print *, MomExt(1:4,nline)
! print *, get_MInv(MomExt(1:4,nline))
         enddo! nline
! pause




! print *, get_MInv2(MomExt(1:4,5)),get_MInv2(MomExt(1:4,6))
! call ShiftMass(MomExt(1:4,5),MomExt(1:4,6),0d0,0d0,MomExt(1:4,1),MomExt(1:4,2))
! print *, get_MInv2(MomExt(1:4,1)),get_MInv2(MomExt(1:4,2))
!
! print *, MomExt(1:4,5)+MomExt(1:4,6) - MomExt(1:4,1)-MomExt(1:4,2)
! print *, MomExt(1:4,5) - MomExt(1:4,1)
! print *, MomExt(1:4,6) - MomExt(1:4,2)
! pause




! MARKUS: NOTE converts only Z-->Z and W-->W, and not Z-->W;
!              converts only to quarks, 2leptons, 3leptons, anything
!              i.e. DecayMode1=0,4,8,10, 1,5, 9,11

         call random_number(xRnd)
         MomShift(:,:) = MomExt(:,:)
         i=1
         do nline=1,EventNumPart
              if( LHE_MOTHUP(1,nline).eq.convertparent .and. LHE_MOTHUP(2,nline).eq.convertparent ) then! found a decay particle
               if( DecayMode1.eq.0 .and. LHE_IDUP(convertparent).eq.convertLHE(Z0_) ) then! convert Z decay products to 2 leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( ZLepBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -ZLepBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.4 .and. LHE_IDUP(convertparent).eq.convertLHE(Wp_) ) then! convert W+ decay products to 2 leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( WLepBranching(xRnd) )
                         LHE_IDUP(nline) = -LHE_IDUP(nline)! converts LepM to LepP
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( - SU2flip(WLepBranching(xRnd)) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.4 .and. LHE_IDUP(convertparent).eq.convertLHE(Wm_) ) then! convert W- decay products to 2 leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( Su2flip(WLepBranching(xRnd)) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -WLepBranching(xRnd) )
                         LHE_IDUP(nline) = -LHE_IDUP(nline)! converts -LepM to +LepM
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.8 .and. LHE_IDUP(convertparent).eq.convertLHE(Z0_) ) then! convert Z decay products to 3 leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( ZLepPlusTauBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -ZLepPlusTauBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.10 .and. LHE_IDUP(convertparent).eq.convertLHE(Wp_) ) then! convert W+ decay products to 3 leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( WLepPlusTauBranching(xRnd) )
                         LHE_IDUP(nline) = -LHE_IDUP(nline)! converts LepM to LepP
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( - SU2flip(WLepPlusTauBranching(xRnd)) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.10 .and. LHE_IDUP(convertparent).eq.convertLHE(Wm_) ) then! convert W- decay products to 3 leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( Su2flip(WLepPlusTauBranching(xRnd)) )
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -WLepPlusTauBranching(xRnd) )
                         LHE_IDUP(nline) = -LHE_IDUP(nline)! converts -LepM to +LepM
                         LHE_ICOLUP(1:2,nline) = (/0,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.1 .and. LHE_IDUP(convertparent).eq.convertLHE(Z0_) ) then! convert Z decay products to quarks
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( ZQuaBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/805,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -ZQuaBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/0,805/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.5 .and. LHE_IDUP(convertparent).eq.convertLHE(Wp_) ) then! convert W+ decay products to quarks
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( WQuaUpBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/805,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( - SU2flip(WQuaUpBranching(xRnd)) )
                         LHE_ICOLUP(1:2,nline) = (/0,805/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.5 .and. LHE_IDUP(convertparent).eq.convertLHE(Wm_) ) then! convert W- decay products to quarks
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( Su2flip(WQuaUpBranching(xRnd)) )
                         LHE_ICOLUP(1:2,nline) = (/805,0/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -WQuaUpBranching(xRnd) )
                         LHE_ICOLUP(1:2,nline) = (/0,805/)
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.9 .and. LHE_IDUP(convertparent).eq.convertLHE(Z0_) ) then! convert Z decay products to quarks and leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( ZAnyBranching(xRnd) )
                         if( IsAQuark(convertLHEreverse(LHE_IDUP(nline))) ) then
                             LHE_ICOLUP(1:2,nline) = (/805,0/)
                         else
                             LHE_ICOLUP(1:2,nline) = (/0,0/)
                         endif
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  else
                         LHE_IDUP(nline) = convertLHE( -ZAnyBranching(xRnd) )
                         if( IsAQuark(convertLHEreverse(LHE_IDUP(nline))) ) then
                             LHE_ICOLUP(1:2,nline) = (/0,805/)
                         else
                             LHE_ICOLUP(1:2,nline) = (/0,0/)
                         endif
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

! branching counter
! if( abs(convertLHEreverse(LHE_IDUP(nline))).ge.7 .and. abs(convertLHEreverse(LHE_IDUP(nline))).le.9 ) Br_Z_ll_counter=Br_Z_ll_counter+1
! if( abs(convertLHEreverse(LHE_IDUP(nline))).ge.14 .and. abs(convertLHEreverse(LHE_IDUP(nline))).le.16 ) Br_Z_inv_counter=Br_Z_inv_counter+1
! if( abs(convertLHEreverse(LHE_IDUP(nline))).eq.Up_ .or.abs(convertLHEreverse(LHE_IDUP(nline))).eq.Chm_ ) Br_Z_uu_counter=Br_Z_uu_counter+1
! if( abs(convertLHEreverse(LHE_IDUP(nline))).eq.Dn_ .or. abs(convertLHEreverse(LHE_IDUP(nline))).eq.Str_ .or. abs(convertLHEreverse(LHE_IDUP(nline))).eq.Bot_) Br_Z_dd_counter=Br_Z_dd_counter+1
! EvalCounter=EvalCounter+1

               elseif( DecayMode1.eq.11 .and. LHE_IDUP(convertparent).eq.convertLHE(Wp_) ) then! convert W+ decay products to quarks and leptons
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( WAnyBranching(xRnd) )
                         if( IsAQuark(convertLHEreverse(LHE_IDUP(nline))) ) then
                              LHE_ICOLUP(1:2,nline) = (/805,0/)
                         else
                              LHE_ICOLUP(1:2,nline) = (/0,0/)
                              LHE_IDUP(nline) = -LHE_IDUP(nline)! converts LepM to LepP
                         endif
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
!                          print *, "here 2",  -WAnyBranching(xRnd), - SU2flip(WAnyBranching(xRnd))
!                          pause
                  else
                         LHE_IDUP(nline) = convertLHE( - SU2flip(WAnyBranching(xRnd)) )
                         if( IsAQuark(convertLHEreverse(LHE_IDUP(nline))) ) then
                              LHE_ICOLUP(1:2,nline) = (/0,805/)
                         else
                              LHE_ICOLUP(1:2,nline) = (/0,0/)
                         endif
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               elseif( DecayMode1.eq.11 .and. LHE_IDUP(convertparent).eq.convertLHE(Wm_) ) then! convert W- decay products to quarks
                  if( LHE_IDUP(nline).gt.0 ) then
                         LHE_IDUP(nline) = convertLHE( SU2flip(WAnyBranching(xRnd)) )
                         if( IsAQuark(convertLHEreverse(LHE_IDUP(nline))) ) then
                              LHE_ICOLUP(1:2,nline) = (/805,0/)
                         else
                              LHE_ICOLUP(1:2,nline) = (/0,0/)
                         endif
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
!                          print *, "here", ( SU2flip(WAnyBranching(xRnd)) ),-( -WAnyBranching(xRnd) )
!                          pause
                  else
                         LHE_IDUP(nline) = convertLHE( -WAnyBranching(xRnd) )
                         if( IsAQuark(convertLHEreverse(LHE_IDUP(nline))) ) then
                              LHE_ICOLUP(1:2,nline) = (/0,805/)
                         else
                              LHE_ICOLUP(1:2,nline) = (/0,0/)
                              LHE_IDUP(nline) = -LHE_IDUP(nline)! converts -LepM to +LepM
                         endif
                         Mass(nline) = getMass( convertLHEreverse(LHE_IDUP(nline)) ) /GeV
                         DecayParticles(i) = nline; i=i+1;
                  endif

               else
                  call Error("Invalid DecayMode1 in StartConvertLHE")
               endif! DecayMode

              endif
         enddo! nline

         call ShiftMass(MomExt(1:4,DecayParticles(1)),MomExt(1:4,DecayParticles(2)),         &
                        getMass( convertLHEreverse(LHE_IDUP(DecayParticles(1))) )/GeV, getMass( convertLHEreverse(LHE_IDUP(DecayParticles(2))) )/GeV,                                                &
                        MomShift(1:4,DecayParticles(1)),MomShift(1:4,DecayParticles(2)))

!          print *, getMass( convertLHEreverse(LHE_IDUP(DecayParticles(1))) ), getMass( convertLHEreverse(LHE_IDUP(DecayParticles(2))) )
!          print *, get_MInv(MomShift(1:4,DecayParticles(1)))
!          print *, get_MInv(MomShift(1:4,DecayParticles(2)))
!          print *, MomShift(1:4,DecayParticles(1))
!          print *, MomShift(1:4,DecayParticles(2))
!          pause

         write(io_LHEOutFile,"(A)") trim(BeginEventLine)
         indent = 0
         do while (BeginEventLine(indent+1:indent+1).eq." ")
             indent = indent+1
         enddo
         if (.not.UseUnformattedRead) then
             IndentedFmt0=InputFmt0
         else
             if (indent.eq.0) then
                 write(IndentedFmt0, "(A,A,A)") "(", DefaultFmt0, ")"
             else
                 write(IndentedFmt0, "(A,I1,A,A,A)") "(", indent, "X,", DefaultFmt0, ")"
             endif
         endif
         if (indent.eq.0) then
             write(IndentedFmt1, "(A,A,A)") "(", Fmt1, ")"
         else
             write(IndentedFmt1, "(A,I1,A,A,A)") "(", indent, "X,", Fmt1, ")"
         endif
         write(io_LHEOutFile,fmt=IndentedFmt0) EventNumPart,EventProcessId,WeightScaleAqedAqcd!  read number of particle from the first line after <event> and other info
         do nline=1,EventNumPart
            write(io_LHEOutFile,fmt=IndentedFmt1) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomShift(2,nline),MomShift(3,nline),MomShift(4,nline),MomShift(1,nline),Mass(nline),Lifetime(nline),Spin(nline)
         enddo

!        read optional lines
         do while (.true.)
              read(io_LHEInFile,fmt="(A120)",IOSTAT=stat,END=99) PDFLine(1:120)
              if(PDFLine(1:30).eq."</LesHouchesEvents>") then
                  goto 99
              elseif( Index(PDFLine,"</event").ne.0 ) then
                  write(io_LHEOutFile,"(A)") "</event>"
              elseif( Index(PDFLine,"<event").ne.0 ) then
                  exit
              else
                  write(io_LHEOutFile,fmt="(A)") trim(PDFLine)
              endif
         enddo
     enddo

99   continue
     call cpu_time(time_end)

RETURN
END SUBROUTINE



SUBROUTINE OpenFiles()
use ModParameters
implicit none
integer :: i

   if( .not.FilesOpened ) then
       FilesOpened = .true.

       print *, ""
       if( unweighted ) then
          open(unit=io_LHEOutFile,file=trim(DataFile)//'.lhe',form='formatted',access= 'sequential',status='replace')        ! LHE event file
          open(unit=io_HistoFile, file=trim(DataFile)//'.dat',form='formatted',access= 'sequential',status='replace')        ! histogram file
       else
          if( (writeWeightedLHE) .or. (GenerateEvents) ) open(unit=io_LHEOutFile,file=trim(DataFile)//'.lhe',form='formatted',access= 'sequential',status='replace')        ! LHE event file
          open(unit=io_HistoFile,file=trim(DataFile)//'.dat',form='formatted',access= 'sequential',status='replace')         ! histogram file
       endif
       open(unit=io_LogFile,file=trim(DataFile)//'.log',form='formatted',access= 'sequential',status='replace')              ! log file


       if( ReadLHEFile .or. ConvertLHEFile ) then
          open(unit=io_LHEInFile,file=trim(LHEProdFile),form='formatted',access= 'sequential',status='old')                  ! LHE input file
       endif
   endif

return
END SUBROUTINE


SUBROUTINE ReopenInFile()
use ModParameters
implicit none

    if ( ReadLHEFile .or. ConvertLHEFile ) then
        close(io_LHEInFile)
        open(unit=io_LHEInFile,file=trim(LHEProdFile),form='formatted',access= 'sequential',status='old')
    endif

return
END SUBROUTINE


SUBROUTINE CloseFiles()
use ModParameters
implicit none

   close(io_LHEOutFile)
   close(io_HistoFile)
   if( ReadLHEFile .or. ConvertLHEFile ) close(io_LHEInFile)
   close(io_LogFile)


return
END SUBROUTINE




SUBROUTINE ClearHisto()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: NHisto


  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE


SUBROUTINE InitHisto()
use modParameters
implicit none

  if( Process.eq.60 .or. (Process.ge.66 .and. Process.le.69) ) then
     call InitHisto_HVBF()
  elseif( Process.eq.61) then
     call InitHisto_HJJ()
  elseif( Process.eq.62) then
     call InitHisto_HJ()
  elseif (Process.eq.50) then
     call InitHisto_VHiggs()
  elseif (Process.eq.51) then
     call InitHisto_VH()
  elseif (Process.eq.52) then
     call InitHisto_HH()
  elseif (Process.eq.80) then
     call InitHisto_TTBH()
  elseif (Process.eq.90) then
     call InitHisto_BBBH()
  elseif (Process.eq.110 .or. Process .eq. 111 .or. Process.eq.112 .or. Process .eq. 113 .or. Process .eq. 114) then
     call InitHisto_TH()
  else

     if( TauDecays.eq.0 .or. TauDecays.eq.1 ) then
         call InitHisto_Htautau()
     else
         call InitHisto_HZZ()
     endif

!      call InitHisto_Htoptop()
  endif

RETURN
END SUBROUTINE




SUBROUTINE InitHisto_HZZ()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

         it_sav = 1
          NumHistograms = 18
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(LepP)"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "pT(LepM)"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "y(LepP)"
          Histo(3)%NBins  = 32
          Histo(3)%BinSize= 0.125
          Histo(3)%LowVal =-2d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y(LepM)"
          Histo(4)%NBins  = 32
          Histo(4)%BinSize= 0.125
          Histo(4)%LowVal =-2d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "cos(theta*)"   ! scattering angle of Z in resonance rest frame
          Histo(5)%NBins  = 200
          Histo(5)%BinSize= 0.01d0
          Histo(5)%LowVal =-1d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "Phi1"  ! angle between plane of beam-scatterin axis and the lepton plane of Z1 in the resonance rest frame
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 0.078539d0 *2d0
          Histo(6)%LowVal =-3.14159d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "cos(theta1)"    ! angle between lepton momentum in the rest frame of its own Z(call it Z1) and
          Histo(7)%NBins  = 40                  ! the direction of the other Z (call it Z2) in the rest frame of Z1
          Histo(7)%BinSize= 0.05d0
          Histo(7)%LowVal =-1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "Phi"  ! angle between lepton planes in resonance rest frame
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.078539d0*2d0
          Histo(8)%LowVal = -3.14159d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "MV1 invariant mass"
          Histo(9)%NBins  = 200
          Histo(9)%BinSize= 0.5d0*GeV
          Histo(9)%LowVal = 00d0*GeV
          Histo(9)%SetScale= 1d0/GeV

          Histo(10)%Info   = "MV2 invariant mass"
          Histo(10)%NBins  = 200
          Histo(10)%BinSize= 0.5d0*GeV
          Histo(10)%LowVal = 00d0*GeV
          Histo(10)%SetScale= 1d0/GeV

          Histo(11)%Info   = "MG resonance invariant mass"
          Histo(11)%NBins  = 200
          Histo(11)%BinSize= 1d0*GeV
          Histo(11)%LowVal = 50d0*GeV
!           Histo(11)%LowVal = 650*GeV !50d0*GeV
          Histo(11)%SetScale= 1d0/GeV

          Histo(12)%Info   = "mLep12 invariant mass"
          Histo(12)%NBins  = 200
          Histo(12)%BinSize= 0.5d0*GeV
          Histo(12)%LowVal = 00d0*GeV
          Histo(12)%SetScale= 1d0/GeV

          Histo(13)%Info   = "mLep34 invariant mass"
          Histo(13)%NBins  = 200
          Histo(13)%BinSize= 0.5d0*GeV
          Histo(13)%LowVal = 00d0*GeV
          Histo(13)%SetScale= 1d0/GeV

          Histo(14)%Info   = "mLep14 invariant mass"
          Histo(14)%NBins  = 200
          Histo(14)%BinSize= 0.5d0*GeV
          Histo(14)%LowVal = 00d0*GeV
          Histo(14)%SetScale= 1d0/GeV

          Histo(15)%Info   = "mLep23 invariant mass"
          Histo(15)%NBins  = 200
          Histo(15)%BinSize= 0.5d0*GeV
          Histo(15)%LowVal = 00d0*GeV
          Histo(15)%SetScale= 1d0/GeV

!           Histo(12)%Info   = "total ee cross section"
!           Histo(12)%NBins  = 1
!           Histo(12)%BinSize= 1d0
!           Histo(12)%LowVal =  0d0
!           Histo(12)%SetScale= 1d0
!
!           Histo(13)%Info   = "total mm cross section"
!           Histo(13)%NBins  = 1
!           Histo(13)%BinSize= 1d0
!           Histo(13)%LowVal = 0d0
!           Histo(13)%SetScale= 1d0
!
!           Histo(14)%Info   = "total tt cross section"
!           Histo(14)%NBins  = 1
!           Histo(14)%BinSize= 1d0
!           Histo(14)%LowVal = 0d010
!           Histo(14)%SetScale= 1d0
!
!           Histo(15)%Info   = "total em cross section"
!           Histo(15)%NBins  = 1
!           Histo(15)%BinSize= 1d0
!           Histo(15)%LowVal = 0d0
!           Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "total et cross section"
          Histo(16)%NBins  = 1
          Histo(16)%BinSize= 1d0
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "total mt cross section"
          Histo(17)%NBins  = 1
          Histo(17)%BinSize= 1d0
          Histo(17)%LowVal = 0d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "total cross section"
          Histo(18)%NBins  = 1
          Histo(18)%BinSize= 1d0
          Histo(18)%LowVal = 0d0
          Histo(18)%SetScale= 1d0


  do NHisto=1,NumHistograms
!       if( .not.allocated(Histo(NHisto)%Value) ) then
!         allocate( Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
!         if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!       endif
!       if( .not.allocated(Histo(NHisto)%Value2) ) then
!         allocate( Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
!         if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!       endif
!       if( .not.allocated(Histo(NHisto)%Hits) ) then
!         allocate( Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
!         if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!       endif
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE





SUBROUTINE InitHisto_HJ()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_j1"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV


  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo

RETURN
END SUBROUTINE




SUBROUTINE InitHisto_HJJ()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m_jj"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "dphi_jj"
          Histo(2)%NBins  = 25
          Histo(2)%BinSize= 0.125d0
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT(H)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "y(j1)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal = -5d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "y(j2)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 0.2d0
          Histo(5)%LowVal = -5d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "dy(y1-j2)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal = -5d0
          Histo(6)%SetScale= 1d0


  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE
SUBROUTINE InitHisto_TTBH()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 10
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_top"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "pT_H"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "mt"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.4d0*GeV
          Histo(3)%LowVal = 160d0*GeV
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "mtbar"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.4d0*GeV
          Histo(4)%LowVal = 160d0*GeV
          Histo(4)%SetScale= 1d0/GeV

          Histo(5)%Info   = "mWp"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 0.4d0*GeV
          Histo(5)%LowVal = 70d0*GeV
          Histo(5)%SetScale= 1d0/GeV

          Histo(6)%Info   = "mWm"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.4d0*GeV
          Histo(6)%LowVal = 70d0*GeV
          Histo(6)%SetScale= 1d0/GeV

          Histo(7)%Info   = "pT_b"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 1d0/GeV

          Histo(8)%Info   = "pT_l"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 1d0/GeV

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 10d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 1d0/GeV

          Histo(10)%Info   = "D_0minus"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 0.02
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 1d0


  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo

RETURN
END SUBROUTINE



SUBROUTINE InitHisto_BBBH()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 2
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_bot"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "pT_H"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0/GeV


  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo

RETURN
END SUBROUTINE




SUBROUTINE InitHisto_TH()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 7
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_top"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "y(top)"
          Histo(2)%NBins  = 60
          Histo(2)%BinSize= 0.1d0
          Histo(2)%LowVal = -3d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT(jet)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "y(jet)"
          Histo(4)%NBins  = 60
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal = -6d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT(Higgs)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal =  0d0*GeV
          Histo(5)%SetScale= 1d0/GeV

          Histo(6)%Info   = "y(Higgs)"
          Histo(6)%NBins  = 100
          Histo(6)%BinSize= 0.1d0
          Histo(6)%LowVal = -5d0
          Histo(6)%SetScale= 1d0
          Histo(7)%Info   = "D_0minus"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.02
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 1d0



  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo

RETURN
END SUBROUTINE

SUBROUTINE InitHisto_Htautau()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto


          it_sav = 1
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m_tauP"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 0.05d0*GeV
          Histo(1)%LowVal = 0d0*GeV
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "m_tauM"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.05d0*GeV
          Histo(2)%LowVal = 0d0*GeV
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "m_Wp"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.05d0*GeV
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "m_Wm"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.05d0*GeV
          Histo(4)%LowVal = 0d0*GeV
          Histo(4)%SetScale= 1d0/GeV

          Histo(5)%Info   = "y_tau+"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 0.25d0
          Histo(5)%LowVal = -5d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "y_tau-"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 0.25d0
          Histo(6)%LowVal = -5d0
          Histo(6)%SetScale= 1d0


          do NHisto=1,NumHistograms
              Histo(NHisto)%Value(:) = 0d0
              Histo(NHisto)%Value2(:)= 0d0
              Histo(NHisto)%Hits(:)  = 0
          enddo

RETURN
END SUBROUTINE





SUBROUTINE InitHisto_Htoptop()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto


          it_sav = 1
          NumHistograms = 4
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m_top"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 0.25d0*GeV
          Histo(1)%LowVal = 166d0*GeV
          Histo(1)%SetScale= 1d0/GeV


          Histo(2)%Info   = "m_topbar"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.25d0*GeV
          Histo(2)%LowVal = 166d0*GeV
          Histo(2)%SetScale= 1d0/GeV


          Histo(3)%Info   = "m_Wp"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0*GeV
          Histo(3)%LowVal = 75d0*GeV
          Histo(3)%SetScale= 1d0/GeV


          Histo(4)%Info   = "m_Wm"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0*GeV
          Histo(4)%LowVal = 75d0*GeV
          Histo(4)%SetScale= 1d0/GeV


          do NHisto=1,NumHistograms
              Histo(NHisto)%Value(:) = 0d0
              Histo(NHisto)%Value2(:)= 0d0
              Histo(NHisto)%Hits(:)  = 0
          enddo

RETURN
END SUBROUTINE




SUBROUTINE InitHisto_HVBF()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms =  10
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m_jj"
          Histo(1)%NBins  = 100
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "m_jj"
          Histo(2)%NBins  = 100
          Histo(2)%BinSize= 2d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "pT(H)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "pT(j1)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 1d0/GeV

          Histo(5)%Info   = "dy(j1,j2)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 0.2d0
          Histo(5)%LowVal = -5d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "mz1"
          Histo(6)%NBins  = 100
          Histo(6)%BinSize= 2d0*GeV
          Histo(6)%LowVal = 0d0*GeV
          Histo(6)%SetScale= 1d0/GeV

          Histo(7)%Info   = "mz2"
          Histo(7)%NBins  = 100
          Histo(7)%BinSize= 2d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 1d0/GeV

          Histo(8)%Info   = "m_4l"
          Histo(8)%NBins  = 100
          Histo(8)%BinSize= 40d0*GeV
          Histo(8)%LowVal = 100d0*GeV
          Histo(8)%SetScale= 1d0/GeV

          Histo(9)%Info   = "m_4l"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 1d0*GeV
          Histo(9)%LowVal = 85d0*GeV
          Histo(9)%SetScale= 1d0/GeV

          Histo(10)%Info   = "Phi1"  ! angle between plane of beam-scatterin axis and the lepton plane of Z1 in the resonance rest frame
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.078539d0 *2d0
          Histo(10)%LowVal =-3.14159d0
          Histo(10)%SetScale= 1d0



  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE






SUBROUTINE InitHisto_VHiggs()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 10
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m(jj)"
          Histo(1)%NBins  = 80
          Histo(1)%BinSize= 20d0*GeV/80d0
          Histo(1)%LowVal = 115d0*GeV
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "m(ll)"
          Histo(2)%NBins  = 80
          Histo(2)%BinSize= 20d0*GeV/80d0
          Histo(2)%LowVal = 75d0*GeV
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "pt(V)"
          Histo(3)%NBins  = 80
          Histo(3)%BinSize= 300d0*GeV/80d0
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "pt(H)"
          Histo(4)%NBins  = 100
          Histo(4)%BinSize= 5d0*GeV
          Histo(4)%LowVal = 0d0*GeV
          Histo(4)%SetScale= 1d0/GeV

          Histo(5)%Info   = "m(V*)"   ! scattering angle of Z in resonance rest frame
          Histo(5)%NBins  = 80
          Histo(5)%BinSize= 300d0*GeV/80d0
          Histo(5)%LowVal = 200d0*GeV
          Histo(5)%SetScale= 1d0/GeV

          Histo(6)%Info   = "costheta1"
          Histo(6)%NBins  = 80
          Histo(6)%BinSize= 2d0/80d0
          Histo(6)%LowVal = -1d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "costheta2"
          Histo(7)%NBins  = 80
          Histo(7)%BinSize= 2d0/80d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "phistar1"
          Histo(8)%NBins  = 80
          Histo(8)%BinSize= 6.4d0/80d0
          Histo(8)%LowVal = -3.2d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "phi"
          Histo(9)%NBins  = 80
          Histo(9)%BinSize= 6.4d0/80d0
          Histo(9)%LowVal = -3.2d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "EHat"
          Histo(10)%NBins  = 80
          Histo(10)%BinSize= 10d0*GeV
          Histo(10)%LowVal = 200d0*GeV
          Histo(10)%SetScale= 1d0/GeV

  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE




SUBROUTINE InitHisto_VH()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 10
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m(jj)"
          Histo(1)%NBins  = 80
          Histo(1)%BinSize= 20d0*GeV/80d0
          Histo(1)%LowVal = 115d0*GeV
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "m(ll)"
          Histo(2)%NBins  = 80
          Histo(2)%BinSize= 20d0*GeV/80d0
          Histo(2)%LowVal = 75d0*GeV
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "pt(V)"
          Histo(3)%NBins  = 80
          Histo(3)%BinSize= 300d0*GeV/80d0
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "pt(H)"
          Histo(4)%NBins  = 100
          Histo(4)%BinSize= 5d0*GeV
          Histo(4)%LowVal = 0d0*GeV
          Histo(4)%SetScale= 1d0/GeV

          Histo(5)%Info   = "m(V*)"   ! scattering angle of Z in resonance rest frame
          Histo(5)%NBins  = 80
          Histo(5)%BinSize= 1000d0*GeV/80d0
          Histo(5)%LowVal = 200d0*GeV
          Histo(5)%SetScale= 1d0/GeV

          Histo(6)%Info   = "costheta1"
          Histo(6)%NBins  = 80
          Histo(6)%BinSize= 2d0/80d0
          Histo(6)%LowVal = -1d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "costheta2"
          Histo(7)%NBins  = 80
          Histo(7)%BinSize= 2d0/80d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "phistar1"
          Histo(8)%NBins  = 80
          Histo(8)%BinSize= 6.4d0/80d0
          Histo(8)%LowVal = -3.2d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "phi"
          Histo(9)%NBins  = 80
          Histo(9)%BinSize= 6.4d0/80d0
          Histo(9)%LowVal = -3.2d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "EHat"
          Histo(10)%NBins  = 80
          Histo(10)%BinSize= 10d0*GeV
          Histo(10)%LowVal = 200d0*GeV
          Histo(10)%SetScale= 1d0/GeV

  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE








SUBROUTINE InitHisto_HH()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 10
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "m(H1)"
          Histo(1)%NBins  = 80
          Histo(1)%BinSize= 20d0*GeV/80d0
          Histo(1)%LowVal = 115d0*GeV
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "m(H2)"
          Histo(2)%NBins  = 80
          Histo(2)%BinSize= 20d0*GeV/80d0
          Histo(2)%LowVal = 75d0*GeV
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "pt(H1)"
          Histo(3)%NBins  = 80
          Histo(3)%BinSize= 300d0*GeV/80d0
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 1d0/GeV

          Histo(4)%Info   = "pt(H2)"
          Histo(4)%NBins  = 80
          Histo(4)%BinSize= 300d0*GeV/80d0
          Histo(4)%LowVal = 0d0*GeV
          Histo(4)%SetScale= 1d0/GeV

          Histo(5)%Info   = "m(HH)"
          Histo(5)%NBins  = 80
          Histo(5)%BinSize= 1000d0*GeV/80d0
          Histo(5)%LowVal = 200d0*GeV
          Histo(5)%SetScale= 1d0/GeV

          Histo(6)%Info   = "costheta*"
          Histo(6)%NBins  = 80
          Histo(6)%BinSize= 2d0/80d0
          Histo(6)%LowVal = -1d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "costheta1"
          Histo(7)%NBins  = 80
          Histo(7)%BinSize= 2d0/80d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "costheta2"
          Histo(8)%NBins  = 80
          Histo(8)%BinSize= 2d0/80d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "phi1"
          Histo(9)%NBins  = 80
          Histo(9)%BinSize= 6.4d0/80d0
          Histo(9)%LowVal = -3.2d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "phi"
          Histo(10)%NBins  = 80
          Histo(10)%BinSize= 6.4d0/80d0
          Histo(10)%LowVal = -3.2d0
          Histo(10)%SetScale= 1d0

  do NHisto=1,NumHistograms
      Histo(NHisto)%Value(:) = 0d0
      Histo(NHisto)%Value2(:)= 0d0
      Histo(NHisto)%Hits(:)  = 0
  enddo


RETURN
END SUBROUTINE









SUBROUTINE InitOutput(CrossSection, CrossSectionError)
use ModParameters
implicit none

integer :: exitstatus
integer :: incoming1, incoming2
real(8) :: beamenergy1, beamenergy2
integer :: pdfgup1, pdfgup2, pdfsup1, pdfsup2  !pdfgup is outdated, set to 0.  pdfsup = LHAGLUE code
integer :: weightscheme, nprocesses
real(8) :: CrossSection, CrossSectionError
character(len=500) :: ReadLines
integer :: stat

    if( Collider.eq.0 ) then
        incoming1 = -11
        incoming2 = 11
    else
        incoming1=2212
        incoming2=2212
    endif
    beamenergy1 = Collider_Energy/GeV / 2d0
    beamenergy2 = Collider_Energy/GeV / 2d0
    pdfgup1 = 0
    pdfgup2 = 0
    if( Collider.eq.0 .or. ReadLHEFile .or. ConvertLHEFile ) then
        pdfsup1 = 0
        pdfsup2 = 0
    else
#if useLHAPDF==1
        !temporary fix, this seems to be the only way
        !email from Andy Buckley:
        !  Hi,
        !  The Fortran interface was originally written in the era when PDFs could *only* be accessed by number rather than name. Accordingly it doesn't provide a way to look up the ID code, because it expects that you need to have known it in advance!
        !  We can add an ID-lookup function in the "new" Fortran interface for the next version, but I think for now this information isn't available other than through the C++ and Python APIs.
        !  Andy
        open(unit=io_TmpFile,file=trim(LHAPDF_DATA_PATH)//"/"//LHAPDFString,form='formatted',access= 'sequential',status='old')
        pdfsup1 = -999
        do while( pdfsup1.lt.0 )
            read(io_TmpFile,fmt="(A160)",IOSTAT=stat,END=99) ReadLines
            if( ReadLines(1:8).eq."SetIndex" ) then
                read(ReadLines(10:500),*) pdfsup1
                pdfsup1 = pdfsup1 + LHAPDFMember
                exit
            endif
        enddo
99      CONTINUE
        close(io_TmpFile)
        if( pdfsup1.lt.0 ) then
            print *, "Error: couldn't get the PDF set index.  Writing 0."
            pdfsup1 = 0
        endif
        pdfsup2 = pdfsup1
#else
        if( PDFSet.eq.1 ) then
            pdfsup1 = 10042
            pdfsup2 = 10042
        elseif( PDFSet.eq.2 ) then
            pdfsup1 = 21000
            pdfsup2 = 21000
        elseif( PDFSet.ge.201.and.PDFSet.le.240 ) then !21041-21080
            pdfsup1 = 21040 + PDFSet - 200
            pdfsup2 = 21040 + PDFSet - 200
        elseif( PDFSet.eq.3 ) then
            pdfsup1 = 263000
            pdfsup2 = 263000
        endif
#endif
    endif
    if( ReweightInterference .or. .not.unweighted ) then
        weightscheme = 4
    else
        weightscheme = 3
    endif
    nprocesses = 1

    if( (unweighted) .or. ( (.not.unweighted) .and. (writeWeightedLHE) )  ) then
        if ( .not. ReadLHEFile .and. .not. ConvertLHEFile ) then
            write(io_LHEOutFile ,'(A)') '<LesHouchesEvents version="1.0">'
        endif
        write(io_LHEOutFile ,'(A)') '<!--'
        write(io_LHEOutFile ,'(A,A,A)') 'Output from the JHUGenerator ',trim(JHUGen_Version),' described in arXiv:1001.3396 [hep-ph], arXiv:1208.4018 [hep-ph], arXiv:1309.4819 [hep-ph], arXiv:1606.03107 [hep-ph]'

        if( writegit ) then
            write(io_LHEOutFile,'(A)') ""
            write(io_LHEOutFile,'(A)') "Current git commit:"
            close(io_LHEOutFile)
            call system("git rev-parse HEAD >> "//trim(DataFile)//".lhe")
            call system("! git diff --name-only | grep .", exitstatus)
            if( exitstatus.ne.0 ) then
                print *, "can't write git commit to the lhe file with a dirty working tree!"
                stop 1
            endif
            open(unit=io_LHEOutFile,file=trim(DataFile)//'.lhe',form='formatted',access='append',status='old')
            write(io_LHEOutFile,'(A)') ""
        endif
        call WriteParameters(io_LHEOutFile)

        if( ReadLHEFile .or. ConvertLHEFile ) then
            write(io_LHEOutFile ,'(A)') ''
        else
            write(io_LHEOutFile ,'(A)') '-->'
            write(io_LHEOutFile ,'(A)') '<init>'
            write(io_LHEOutFile ,'(I4,X,I4,F24.16,X,F24.16,X,I1,X,I1,X,I7,X,I7,X,I2,X,I2)') incoming1, incoming2, beamenergy1, beamenergy2, pdfgup1, pdfgup2, pdfsup1, pdfsup2, weightscheme, nprocesses
! in order of appearance:  (see also http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017)
! (*) incoming particle1 (2212=proton), incoming particle2,
! (*) energies of colliding particles,
! (*) out-dated pdf information for colliding particles (supposed to be 0),
! (*) pdf code of LHAGLUE for colliding particles (10042=CTEQ6Ll, MSTW2008=21000,21041-21080)
! (*) weighting strategy (3=unweighted events, 4=weighted events,  otherwise=see LHE manuals)
! (*) number of process types to be accepted (default=1, otherwise=see manual)
            write(io_LHEOutFile ,'(1PE14.7,1X,1PE14.7,1X,1PE14.7,1X,I3)') CrossSection, CrossSectionError, 1.00000000000E-00, Process
! in order of appearance:
! (*) total cross section in pb
! (*) stat. error in the total cross section in pb
! (*) maximum weight
! (*) list of all user process ID's
            write(io_LHEOutFile ,'(A)') '</init>'
        endif

    endif


END SUBROUTINE




SUBROUTINE FinalizeOutput
use ModParameters
implicit none

    if( (unweighted) .or. ( (.not.unweighted) .and. (writeWeightedLHE) )  ) then
        write(io_LHEOutFile ,'(A)') '</LesHouchesEvents>'

        write(io_LogFile,'(A)') ''
        write(io_LogFile,'(A,I9)') 'Accepted Events: ',AccepCounter
        write(io_LogFile,'(A)') 'Branchings:'
        write(io_LogFile,'(12X,A,F6.3,A,F6.3,A)') 'Br_Z_ll_counter  ',dble(Br_Z_ll_counter)/AccepCounter *100d0," %    +/-",dsqrt(dble(Br_Z_ll_counter))/AccepCounter *100d0," %"
        write(io_LogFile,'(12X,A,F6.3,A,F6.3,A)') 'Br_Z_inv_counter ',dble(Br_Z_inv_counter)/AccepCounter*100d0," %    +/-",dsqrt(dble(Br_Z_inv_counter))/AccepCounter *100d0," %"
        write(io_LogFile,'(12X,A,F6.3,A,F6.3,A)') 'Br_Z_uu_counter  ',dble(Br_Z_uu_counter)/AccepCounter *100d0," %    +/-",dsqrt(dble(Br_Z_uu_counter))/AccepCounter *100d0," %"
        write(io_LogFile,'(12X,A,F6.3,A,F6.3,A)') 'Br_Z_dd_counter  ',dble(Br_Z_dd_counter)/AccepCounter *100d0," %    +/-",dsqrt(dble(Br_Z_dd_counter))/AccepCounter *100d0," %"
        write(io_LogFile,'(12X,A,F6.3,A,F6.3,A)') 'Br_W_ll_counter  ',dble(Br_W_ll_counter)/AccepCounter *100d0," %    +/-",dsqrt(dble(Br_W_ll_counter))/AccepCounter *100d0," %"
        write(io_LogFile,'(12X,A,F6.3,A,F6.3,A)') 'Br_W_ud_counter  ',dble(Br_W_ud_counter)/AccepCounter *100d0," %    +/-",dsqrt(dble(Br_W_ud_counter))/AccepCounter *100d0," %"
    endif

    if( any(DebugCounter.ne.0) ) then

       print *, "DebugCounter(:)=",DebugCounter(:)

    endif

END SUBROUTINE





SUBROUTINE WriteParameters(TheUnit)
use ModParameters
implicit none
integer :: TheUnit
character :: arg*(500)

    call Get_Command(arg)
    write(TheUnit,"(3X,A)") ""
    write(TheUnit,"(3X,A,A)") "Command line: ",trim(arg)
    write(TheUnit,"(3X,A)") ""

    write(TheUnit,"(3X,A)") "Input Parameter:"
    if( Collider.eq.0 ) write(TheUnit,"(4X,A,1F8.2)") "Collider: e+ e-, sqrt(s)=",Collider_Energy*100d0
    if( Collider.eq.1 ) write(TheUnit,"(4X,A,1F8.2)") "Collider: P-P, sqrt(s)=",Collider_Energy*100d0
    if( Collider.eq.2 ) write(TheUnit,"(4X,A,1F8.2)") "Collider: P-Pbar, sqrt(s)=",Collider_Energy*100d0
    if( Process.eq.0 ) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.1 ) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=1, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.2 ) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=2, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.60) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.61) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.62) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( (Process.eq.66 .or. Process.eq.68) .and. M_Reso.ge.0d0  ) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "1st Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( (Process.eq.66 .or. Process.eq.68) .and. M_Reso2.ge.0d0 ) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "2nd Resonance: spin=0, mass=",M_Reso2*100d0," width=",Ga_Reso2*100d0
    if( Process.eq.50) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.51) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.80) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.90) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.110) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.111) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.112) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.113) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.114) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( ReadLHEFile )    write(TheUnit,"(4X,A)") "           (This is ReadLHEFile mode. Resonance mass/width are read from LHE input parameters.)"
    if( ConvertLHEFile ) write(TheUnit,"(4X,A)") "           (This is ConvertLHEFile mode. Resonance mass/width are read from LHE input parameters.)"
    if( HiggsDecayLengthMM.ne.0d0 ) write(TheUnit,"(4X,A,F10.5,A)") "           ctau=", HiggsDecayLengthMM, " mm"
    if( &
         (.not.ReadLHEFile .and. (Process.le.2 .or. Process.eq.50 .or. Process.eq.60 .or. (Process.ge.66 .and. Process.le.69) .or. ((TopDecays.eq.1).and.Process.eq.80) .or. (Process.ge.110 .and. Process.le.113))) &
    .or. (ReadLHEFile .and. TauDecays.ne.0) &
    .or. ConvertLHEFile ) &
    then
        if( .not.ReadLHEFile .and. (ConvertLHEFile .or. Process.eq.50 .or. (Process.ge.110 .and. Process.le.114)) ) then
            write(TheUnit,"(4X,A,I2,2X,A,I2)") "DecayMode1:",DecayMode1
        else if( ReadLHEFile .or. Process.le.2 .or. Process .eq. 80 ) then
            write(TheUnit,"(4X,A,I2,2X,A,I2)") "DecayMode1:",DecayMode1, "DecayMode2:",DecayMode2
        endif
        if( Process.eq.60 .or. (Process.ge.66 .and. Process.le.69) .or. IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z boson: mass=",M_Z*100d0,", width=",Ga_Z*100d0
        if( Process.eq.60 .or. (Process.ge.66 .and. Process.le.69) .or. IsAWDecay(DecayMode1) .or. IsAWDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W boson: mass=",M_W*100d0,", width=",Ga_W*100d0
    endif
    if( Process.eq.80 .or. Process.eq.110 .or. Process.eq.111 .or.Process.eq.112 .or. Process.eq.113 .or. Process.eq.114) write(TheUnit,"(4X,A,F8.4,A,F6.4)") "Top quark mass=",m_top*100d0,", width=",Ga_top*100d0
    if( Process.eq.80 .or. Process.eq.110 .or. Process.eq.111 .or. Process.eq.112 .or. Process.eq.113 .or. Process.eq.114) write(TheUnit,"(4X,A,I2)") "Top quark decay=",TOPDECAYS
    if( Process.eq.90 ) write(TheUnit,"(4X,A,F8.4,A,F6.4)") "Bottom quark mass=",m_top*100d0
    if( Process.eq.50 .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.62 .or. (Process.ge.66 .and. Process.le.69) .or. Process.eq.90 .or. &
       ((Process.eq.80 .or. (Process.ge.110 .and. Process.le.114)) .and. m_Top.lt.10d0*GeV) ) then
        write(TheUnit,"(4X,A)") "Jet cuts:"
        write(TheUnit,"(12X,A,F8.2,A)") "pT >= ", pTjetcut/GeV, " GeV"
        if( Process.ge.66 .and. Process.le.69 ) then
            write(TheUnit,"(9X,A,F8.2)") "|eta| <= ", etajetcut
            write(TheUnit,"(4X,A,F8.2)") "|Deltaeta| >= ", detajetcut
            if (JetsOppositeEta) write(TheUnit,"(5X,A,F8.2)") "eta1*eta2 <= ", 0d0
        endif
        if( Process.eq.50 .or. Process.eq.60 .or. Process.eq.61 .or. (Process.ge.66 .and. Process.le.69) .or. Process.eq.80 .or. Process.eq.90) then
            write(TheUnit,"(8X,A,F8.2)") "DeltaR >= ", Rjet
            write(TheUnit,"(11X,A,F8.2,A)") "mJJ >= ", mJJcut/GeV, " GeV"
        endif
    endif
    if( Process.ge.66 .and. Process.le.69 ) then
        write(TheUnit,"(4X,A)") "4l cuts:"
        if( Process.ge.66 .and. Process.le.69 ) then
            write(TheUnit,"(12X,A,F8.2,A)") "pT >= ", pTlepcut/GeV, " GeV"
            write(TheUnit,"(9X,A,F8.2)") "|eta| <= ", etalepcut
        endif
        write(TheUnit,"(11X,A,F8.2,A)") "mll >= ", MPhotonCutoff/GeV, " GeV"
        if( Process.ge.66 .and. Process.le.69 ) then
            write(TheUnit,"(F10.2,A,F10.2,A)") m4l_minmax(1)/GeV, " GeV <= m4l <= ", m4l_minmax(2)/GeV, " GeV"
        endif
    endif
    if( (Process.eq.0 .or. Process.eq.2 .or. Process.eq.50) .and. includeGammaStar ) then
        write(TheUnit,"(6X,A,F8.2,A)") "m(gammastar) >= ", MPhotonCutoff/GeV, " GeV"
    endif
    if( (ReadLHEFile) .and. (RequestNLeptons.gt.0) ) then
        if ( RequestOS .le. 0 ) then
            if ( RequestNLeptons .eq. 1 ) then
                write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," lepton."
            else
                write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons."
            endif
        elseif ( RequestOSSF .le. 0 ) then
            if ( RequestOS*2 .eq. RequestNLeptons ) then
                if ( RequestOS .eq. 1 ) then
                    write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons in an OS pair."
                else
                    write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons in OS pairs."
                endif
            else
                if ( RequestOS .eq. 1 ) then
                    write(TheUnit,"(4X,A,I2,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons including ", RequestOS, " OS pair."
                else
                    write(TheUnit,"(4X,A,I2,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons including ", RequestOS, " OS pairs."
                endif
            endif
        elseif ( RequestOSSF*2 .eq. RequestNLeptons ) then
            if ( RequestOSSF .eq. 1 ) then
                write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons in an OSSF pair."
            else
                write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons in OSSF pairs."
            endif
        elseif ( RequestOSSF .eq. RequestOS ) then
            if ( RequestOSSF .eq. 1 ) then
                write(TheUnit,"(4X,A,I2,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons including ", RequestOSSF, " OSSF pair."
            else
                write(TheUnit,"(4X,A,I2,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons including ", RequestOSSF, " OSSF pairs."
            endif
        elseif ( RequestOS*2 .eq. RequestNLeptons ) then
            write(TheUnit,"(4X,A,I2,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons in OS pairs, ", RequestOSSF, " of the pairs OSSF."
        else ! this will never happen
            write(TheUnit,"(4X,A,I2,A,I2,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons including ", RequestOS, " OS pairs, ", RequestOSSF, " of them OSSF."
        endif
    endif
    if( CountTauAsAny .and. RequestOSSF.gt.0 ) then
        write(TheUnit,"(8X,A)") "(counting tau in place of e or mu of the same sign, if necessary)"
    endif
    if( (ReadLHEFile) .and. (RequestNJets.gt.0) ) then
        if ( RequestNJets .eq. 1 ) then
            write(TheUnit,"(4X,A,I2,A)") "Jet filter activated. Requesting ",RequestNJets," quark/gluon."
        else
            write(TheUnit,"(4X,A,I2,A)") "Jet filter activated. Requesting ",RequestNJets," quarks/gluons."
        endif
    endif
    write(TheUnit,"(4X,A,20I11)") "Random seed: ",UserSeed
    write(TheUnit,"(4X,A)") "To reproduce results using this seed, JHUGen should be compiled with the same compiler:"
#if compiler==1
    if( modulo(__INTEL_COMPILER, 10) == 0 ) then !last digit is 0, e.g. 1110.  ifort --version writes this as 11.1, not 11.10
        write(TheUnit,"(6X,A,I2,A,I1,A,I8)") "ifort (IFORT) ", __INTEL_COMPILER/100, ".", modulo(__INTEL_COMPILER, 100)/10, " ", __INTEL_COMPILER_BUILD_DATE
    else
        write(TheUnit,"(6X,A,I2,A,I2,A,I8)") "ifort (IFORT) ", __INTEL_COMPILER/100, ".", modulo(__INTEL_COMPILER, 100), " ", __INTEL_COMPILER_BUILD_DATE
    endif
#elif compiler==2
    write(TheUnit,"(6X,A,A)") "GNU Fortran (GCC) ", __VERSION__
#endif

    if( .not. (ReadLHEFile .or. ConvertLHEFile) ) then
        write(TheUnit,"(4X,A)") ""
        if( Process.le.2 ) write(TheUnit,"(4X,A,L,L,L)") "OffXVV: ",OffShellReson,OffShellV1,OffShellV2
        write(TheUnit,"(4X,A,I3)") "WidthScheme: ",WidthScheme
        if( Process.le.2 .or. Process.eq.80 .or. Process.eq.90 ) write(TheUnit,"(4X,A,I1)") "PChannel: ",PChannel
#if useLHAPDF==1
        write(TheUnit,"(4X,A,A,A,I3)") "LHAPDF set ",trim(LHAPDFString), ", member=",LHAPDFMember
#else
        write(TheUnit,"(4X,A,I3)") "PDFSet: ",PDFSet
#endif
        write(TheUnit,"(4X,A,L)") "Unweighted: ",Unweighted
    endif
    if( WidthScheme.ne.WidthSchemeIn ) then
        write(TheUnit,"(4X,A,I1,A,I1)") "Reweighting propagator from WidthScheme ", WidthSchemeIn, " to WidthScheme ", WidthScheme
    endif
    if( ReweightDecay ) then
        write(TheUnit,"(4X,A,I1)") "Reweighting events using the decay matrix element, using input WidthScheme ", WidthSchemeIn
    endif
    if(Process.ge.66 .and. Process.le.69 .and. ReweightInterference) then
      write(TheUnit, "(4X,A)") "Interference: included through event weights"
    elseif( Process.le.2 .or. (Process.ge.66 .and. Process.le.69) .or. ReadLHEFile ) then
      write(TheUnit,"(4X,A,L)") "Interference: ",includeInterference
    endif

    if( &
        ( (Process.le.2 .or. ReadLHEFile) .and. (IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2)) ) .or. &
        Process.eq.60 .or. (Process.ge.66 .and. Process.le.68)                                            &
      ) write(TheUnit,"(4X,A,L)") "Intermediate off-shell photons: ",includeGammaStar
    if( Process.eq.69 ) write(TheUnit,"(4X,A,L)") "Intermediate off-shell gluons: ",includeGammaStar

    write(TheUnit,"(4X,A)") ""
    if( (Process.eq.0 .and. TauDecays.lt.0) .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.62 .or. Process.eq.66 .or. Process.eq.68 .or. Process.eq.50 .or. (Process.eq.51 .and. VH_PC.ne."bo") ) then
        write(TheUnit,"(4X,A)") "spin-0-VV couplings: "
        write(TheUnit,"(6X,A,L)") "generate_as=",generate_as
        if( generate_as ) then
            write(TheUnit,"(6X,A,2E16.8,A1)") "ahg1=",ahg1,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ahg2=",ahg2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ahg3=",ahg3,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ahz1=",ahz1,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ahz2=",ahz2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ahz3=",ahz3,"i"
        else
            if( cdabs(ghg2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghg2=",ghg2,"i"
            if( cdabs(ghg3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghg3=",ghg3,"i"
            if( cdabs(ghg4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghg4=",ghg4,"i"
            if( cdabs(ghz1 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghz1=",ghz1,"i"
            if( cdabs(ghz2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghz2=",ghz2,"i"
            if( cdabs(ghz3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghz3=",ghz3,"i"
            if( cdabs(ghz4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghz4=",ghz4,"i"
            if( cdabs(ghz1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime= ",ghz1_prime ,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime2=",ghz1_prime2,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime3=",ghz1_prime3,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime4=",ghz1_prime4,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz1_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime5=",ghz1_prime5,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz1_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime6=",ghz1_prime6,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz1_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime7=",ghz1_prime7,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghz2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime= ",ghz2_prime ,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime2=",ghz2_prime2,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime3=",ghz2_prime3,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime4=",ghz2_prime4,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz2_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime5=",ghz2_prime5,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz2_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime6=",ghz2_prime6,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz2_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime7=",ghz2_prime7,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghz3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime= ",ghz3_prime ,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime2=",ghz3_prime2,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime3=",ghz3_prime3,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime4=",ghz3_prime4,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz3_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime5=",ghz3_prime5,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz3_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime6=",ghz3_prime6,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz3_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime7=",ghz3_prime7,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghz4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime= ",ghz4_prime ,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghz4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime2=",ghz4_prime2,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghz4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime3=",ghz4_prime3,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghz4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime4=",ghz4_prime4,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghz4_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime5=",ghz4_prime5,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghz4_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime6=",ghz4_prime6,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghz4_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime7=",ghz4_prime7,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzgs2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs2=",ghzgs2,"i"
            if( cdabs(ghzgs3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs3=",ghzgs3,"i"
            if( cdabs(ghzgs4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs4=",ghzgs4,"i"
            if( cdabs(ghzgs1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,2X,A,1PE12.4)") "ghzgs1_prime2=",ghzgs1_prime2,"i,","Lambda_zgs1=",Lambda_zgs1/GeV
            if( cdabs(ghgsgs2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs2=",ghgsgs2,"i"
            if( cdabs(ghgsgs3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs3=",ghgsgs3,"i"
            if( cdabs(ghgsgs4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs4=",ghgsgs4,"i"
            if( cdabs(ghzzp1 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzzp1=",ghzzp1,"i"
            if( cdabs(ghzzp2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzzp2=",ghzzp2,"i"
            if( cdabs(ghzzp3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzzp3=",ghzzp3,"i"
            if( cdabs(ghzzp4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzzp4=",ghzzp4,"i"
            if( cdabs(ghzzp1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime= ",ghzzp1_prime ,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime2=",ghzzp1_prime2,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime3=",ghzzp1_prime3,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime4=",ghzzp1_prime4,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp1_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime5=",ghzzp1_prime5,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp1_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime6=",ghzzp1_prime6,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp1_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp1_prime7=",ghzzp1_prime7,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzzp2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime= ",ghzzp2_prime ,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime2=",ghzzp2_prime2,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime3=",ghzzp2_prime3,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime4=",ghzzp2_prime4,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp2_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime5=",ghzzp2_prime5,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp2_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime6=",ghzzp2_prime6,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp2_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp2_prime7=",ghzzp2_prime7,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzzp3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime= ",ghzzp3_prime ,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime2=",ghzzp3_prime2,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime3=",ghzzp3_prime3,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime4=",ghzzp3_prime4,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp3_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime5=",ghzzp3_prime5,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp3_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime6=",ghzzp3_prime6,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp3_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp3_prime7=",ghzzp3_prime7,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzzp4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime= ",ghzzp4_prime ,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzzp4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime2=",ghzzp4_prime2,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzzp4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime3=",ghzzp4_prime3,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzzp4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime4=",ghzzp4_prime4,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzzp4_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime5=",ghzzp4_prime5,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzzp4_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime6=",ghzzp4_prime6,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzzp4_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzzp4_prime7=",ghzzp4_prime7,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpgs2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpgs2=",ghzpgs2,"i"
            if( cdabs(ghzpgs3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpgs3=",ghzpgs3,"i"
            if( cdabs(ghzpgs4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpgs4=",ghzpgs4,"i"
            if( cdabs(ghzpgs1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,2X,A,1PE12.4)") "ghzpgs1_prime2=",ghzpgs1_prime2,"i,","Lambda_zgs1=",Lambda_zgs1/GeV
            if( cdabs(ghzpzp1 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpzp1=",ghzpzp1,"i"
            if( cdabs(ghzpzp2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpzp2=",ghzpzp2,"i"
            if( cdabs(ghzpzp3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpzp3=",ghzpzp3,"i"
            if( cdabs(ghzpzp4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghzpzp4=",ghzpzp4,"i"
            if( cdabs(ghzpzp1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime= ",ghzpzp1_prime ,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime2=",ghzpzp1_prime2,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime3=",ghzpzp1_prime3,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime4=",ghzpzp1_prime4,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp1_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime5=",ghzpzp1_prime5,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp1_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime6=",ghzpzp1_prime6,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp1_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp1_prime7=",ghzpzp1_prime7,"i,","Lambda_z1=",Lambda_z1/GeV
            if( cdabs(ghzpzp2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime= ",ghzpzp2_prime ,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime2=",ghzpzp2_prime2,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime3=",ghzpzp2_prime3,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime4=",ghzpzp2_prime4,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp2_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime5=",ghzpzp2_prime5,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp2_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime6=",ghzpzp2_prime6,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp2_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp2_prime7=",ghzpzp2_prime7,"i,","Lambda_z2=",Lambda_z2/GeV
            if( cdabs(ghzpzp3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime= ",ghzpzp3_prime ,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime2=",ghzpzp3_prime2,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime3=",ghzpzp3_prime3,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime4=",ghzpzp3_prime4,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp3_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime5=",ghzpzp3_prime5,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp3_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime6=",ghzpzp3_prime6,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp3_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp3_prime7=",ghzpzp3_prime7,"i,","Lambda_z3=",Lambda_z3/GeV
            if( cdabs(ghzpzp4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime= ",ghzpzp4_prime ,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpzp4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime2=",ghzpzp4_prime2,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpzp4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime3=",ghzpzp4_prime3,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpzp4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime4=",ghzpzp4_prime4,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpzp4_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime5=",ghzpzp4_prime5,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpzp4_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime6=",ghzpzp4_prime6,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cdabs(ghzpzp4_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghzpzp4_prime7=",ghzpzp4_prime7,"i,","Lambda_z4=",Lambda_z4/GeV
            if( cz_q1sq.ne.0) then
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z11= ",Lambda_z11/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z21= ",Lambda_z21/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z31= ",Lambda_z31/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z41= ",Lambda_z41/GeV
            endif
            if( cz_q2sq.ne.0) then
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z12= ",Lambda_z12/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z22= ",Lambda_z22/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z32= ",Lambda_z32/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z42= ",Lambda_z42/GeV
            endif
            if( cz_q12sq.ne.0) then
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z10= ",Lambda_z10/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z20= ",Lambda_z20/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z30= ",Lambda_z30/GeV
               write(TheUnit,"(6X,A,1PE12.4)") "Lambda_z40= ",Lambda_z40/GeV
            endif

            if(distinguish_HWWcouplings) then
               if( cdabs(ghw1 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghw1=",ghw1,"i"
               if( cdabs(ghw2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghw2=",ghw2,"i"
               if( cdabs(ghw3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghw3=",ghw3,"i"
               if( cdabs(ghw4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghw4=",ghw4,"i"
               if( cdabs(ghw1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime= ",ghw1_prime ,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime2=",ghw1_prime2,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime3=",ghw1_prime3,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime4=",ghw1_prime4,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw1_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime5=",ghw1_prime5,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw1_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime6=",ghw1_prime6,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw1_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw1_prime7=",ghw1_prime7,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghw2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime= ",ghw2_prime ,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime2=",ghw2_prime2,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime3=",ghw2_prime3,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime4=",ghw2_prime4,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw2_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime5=",ghw2_prime5,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw2_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime6=",ghw2_prime6,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw2_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw2_prime7=",ghw2_prime7,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghw3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime= ",ghw3_prime ,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime2=",ghw3_prime2,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime3=",ghw3_prime3,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime4=",ghw3_prime4,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw3_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime5=",ghw3_prime5,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw3_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime6=",ghw3_prime6,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw3_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw3_prime7=",ghw3_prime7,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghw4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime= ",ghw4_prime ,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghw4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime2=",ghw4_prime2,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghw4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime3=",ghw4_prime3,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghw4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime4=",ghw4_prime4,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghw4_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime5=",ghw4_prime5,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghw4_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime6=",ghw4_prime6,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghw4_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghw4_prime7=",ghw4_prime7,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp1 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwwp1=",ghwwp1,"i"
               if( cdabs(ghwwp2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwwp2=",ghwwp2,"i"
               if( cdabs(ghwwp3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwwp3=",ghwwp3,"i"
               if( cdabs(ghwwp4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwwp4=",ghwwp4,"i"
               if( cdabs(ghwwp1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime= ",ghwwp1_prime ,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime2=",ghwwp1_prime2,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime3=",ghwwp1_prime3,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime4=",ghwwp1_prime4,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp1_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime5=",ghwwp1_prime5,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp1_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime6=",ghwwp1_prime6,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp1_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp1_prime7=",ghwwp1_prime7,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwwp2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime= ",ghwwp2_prime ,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime2=",ghwwp2_prime2,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime3=",ghwwp2_prime3,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime4=",ghwwp2_prime4,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp2_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime5=",ghwwp2_prime5,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp2_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime6=",ghwwp2_prime6,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp2_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp2_prime7=",ghwwp2_prime7,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwwp3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime= ",ghwwp3_prime ,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime2=",ghwwp3_prime2,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime3=",ghwwp3_prime3,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime4=",ghwwp3_prime4,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp3_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime5=",ghwwp3_prime5,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp3_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime6=",ghwwp3_prime6,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp3_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp3_prime7=",ghwwp3_prime7,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwwp4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime= ",ghwwp4_prime ,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime2=",ghwwp4_prime2,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime3=",ghwwp4_prime3,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime4=",ghwwp4_prime4,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp4_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime5=",ghwwp4_prime5,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp4_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime6=",ghwwp4_prime6,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwwp4_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwwp4_prime7=",ghwwp4_prime7,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp1 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwpwp1=",ghwpwp1,"i"
               if( cdabs(ghwpwp2 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwpwp2=",ghwpwp2,"i"
               if( cdabs(ghwpwp3 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwpwp3=",ghwpwp3,"i"
               if( cdabs(ghwpwp4 ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ghwpwp4=",ghwpwp4,"i"
               if( cdabs(ghwpwp1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime= ",ghwpwp1_prime ,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime2=",ghwpwp1_prime2,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime3=",ghwpwp1_prime3,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime4=",ghwpwp1_prime4,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp1_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime5=",ghwpwp1_prime5,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp1_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime6=",ghwpwp1_prime6,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp1_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp1_prime7=",ghwpwp1_prime7,"i,","Lambda_w1=",Lambda_w1/GeV
               if( cdabs(ghwpwp2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime= ",ghwpwp2_prime ,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime2=",ghwpwp2_prime2,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime3=",ghwpwp2_prime3,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime4=",ghwpwp2_prime4,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp2_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime5=",ghwpwp2_prime5,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp2_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime6=",ghwpwp2_prime6,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp2_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp2_prime7=",ghwpwp2_prime7,"i,","Lambda_w2=",Lambda_w2/GeV
               if( cdabs(ghwpwp3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime= ",ghwpwp3_prime ,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime2=",ghwpwp3_prime2,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime3=",ghwpwp3_prime3,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime4=",ghwpwp3_prime4,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp3_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime5=",ghwpwp3_prime5,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp3_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime6=",ghwpwp3_prime6,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp3_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp3_prime7=",ghwpwp3_prime7,"i,","Lambda_w3=",Lambda_w3/GeV
               if( cdabs(ghwpwp4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime= ",ghwpwp4_prime ,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime2=",ghwpwp4_prime2,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime3=",ghwpwp4_prime3,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime4=",ghwpwp4_prime4,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp4_prime5).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime5=",ghwpwp4_prime5,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp4_prime6).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime6=",ghwpwp4_prime6,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cdabs(ghwpwp4_prime7).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghwpwp4_prime7=",ghwpwp4_prime7,"i,","Lambda_w4=",Lambda_w4/GeV
               if( cw_q1sq.ne.0) then
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w11= ",Lambda_w11/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w21= ",Lambda_w21/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w31= ",Lambda_w31/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w41= ",Lambda_w41/GeV
               endif
               if( cw_q2sq.ne.0) then
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w12= ",Lambda_w12/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w22= ",Lambda_w22/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w32= ",Lambda_w32/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w42= ",Lambda_w42/GeV
               endif
               if( cw_q12sq.ne.0) then
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w10= ",Lambda_w10/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w20= ",Lambda_w20/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w30= ",Lambda_w30/GeV
                  write(TheUnit,"(6X,A,1PE12.4)") "Lambda_w40= ",Lambda_w40/GeV
               endif
            endif
        endif
        if(includeVprime) then
            if( Process.eq.60 .or. (Process.ge.66 .and. Process.le.69) .or. IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2) ) then
                if(M_Zprime.gt.0d0) then
                  write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z' boson: mass=",M_Zprime*100d0,", width=",Ga_Zprime*100d0
                else
                  write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z' boson: heavy mass limit (contact interaction)"
                endif
                if( cdabs(ezp_El_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_El_left=  ",ezp_El_left  ,"i"
                if( cdabs(ezp_El_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_El_right= ",ezp_El_right ,"i"
                if( cdabs(ezp_Mu_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Mu_left=  ",ezp_Mu_left  ,"i"
                if( cdabs(ezp_Mu_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Mu_right= ",ezp_Mu_right ,"i"
                if( cdabs(ezp_Ta_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Ta_left=  ",ezp_Ta_left  ,"i"
                if( cdabs(ezp_Ta_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Ta_right= ",ezp_Ta_right ,"i"
                if( cdabs(ezp_NuE_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_NuE_left= ",ezp_NuE_left ,"i"
                if( cdabs(ezp_NuE_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_NuE_right=",ezp_NuE_right,"i"
                if( cdabs(ezp_Up_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Up_left=  ",ezp_Up_left  ,"i"
                if( cdabs(ezp_Up_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Up_right= ",ezp_Up_right ,"i"
                if( cdabs(ezp_Dn_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Dn_left=  ",ezp_Dn_left  ,"i"
                if( cdabs(ezp_Dn_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Dn_right= ",ezp_Dn_right ,"i"
                if( cdabs(ezp_Str_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Str_left= ",ezp_Str_left ,"i"
                if( cdabs(ezp_Str_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Str_right=",ezp_Str_right,"i"
                if( cdabs(ezp_Chm_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Chm_left= ",ezp_Chm_left ,"i"
                if( cdabs(ezp_Chm_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Chm_right=",ezp_Chm_right,"i"
                if( cdabs(ezp_Bot_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Bot_left= ",ezp_Bot_left ,"i"
                if( cdabs(ezp_Bot_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Bot_right=",ezp_Bot_right,"i"
                if( cdabs(ezp_Top_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Top_left= ",ezp_Top_left ,"i"
                if( cdabs(ezp_Top_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ezp_Top_right=",ezp_Top_right,"i"
            endif
            if( Process.eq.60 .or. (Process.ge.66 .and. Process.le.69) .or. IsAWDecay(DecayMode1) .or. IsAWDecay(DecayMode2) ) then
                if(M_Wprime.gt.0d0) then
                  write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W' boson: mass=",M_Wprime*100d0,", width=",Ga_Wprime*100d0
                else
                  write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W' boson: mass=heavy (contact interaction)"
                endif
                if( cdabs(ewp_El_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_El_left=  ",ewp_El_left  ,"i"
                if( cdabs(ewp_El_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_El_right= ",ewp_El_right ,"i"
                if( cdabs(ewp_Mu_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Mu_left=  ",ewp_Mu_left  ,"i"
                if( cdabs(ewp_Mu_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Mu_right= ",ewp_Mu_right ,"i"
                if( cdabs(ewp_Ta_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Ta_left=  ",ewp_Ta_left  ,"i"
                if( cdabs(ewp_Ta_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Ta_right= ",ewp_Ta_right ,"i"
                if( cdabs(ewp_Up_left  ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Up_left=  ",ewp_Up_left  ,"i"
                if( cdabs(ewp_Up_right ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Up_right= ",ewp_Up_right ,"i"
                if( cdabs(ewp_Chm_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Chm_left= ",ewp_Chm_left ,"i"
                if( cdabs(ewp_Chm_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Chm_right=",ewp_Chm_right,"i"
                if( cdabs(ewp_Top_left ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Top_left= ",ewp_Top_left ,"i"
                if( cdabs(ewp_Top_right).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "ewp_Top_right=",ewp_Top_right,"i"
            endif
        endif
    endif
    if( (Process.eq.0 .and. TauDecays.ge.0) .or. Process.eq.80 .or. Process.eq.90 .or. (Process.eq.51 .and. VH_PC.ne."tr" .and. VH_PC.ne."ee" .and. VH_PC.ne."qq") ) then
        write(TheUnit,"(4X,A)") "spin-0-ff couplings: "
        if( cdabs(kappa ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "kappa=",kappa,"i"
        if( cdabs(kappa_tilde ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A1)") "kappa_tilde=",kappa_tilde,"i"
    endif
    if( Process.eq.1 ) then
        write(TheUnit,"(4X,A)") "spin-1-VV couplings: "
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_qq_left =",zprime_qq_left,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_qq_right=",zprime_qq_right,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_zz_1=",zprime_zz_1,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_zz_2=",zprime_zz_2,"i"
    endif
    if( Process.eq.2 ) then
        write(TheUnit,"(4X,A)") "spin-2-VV couplings: "
        write(TheUnit,"(6X,A,L)") "generate_bis=",generate_bis
        write(TheUnit,"(6X,A,L)") "use_dynamic_MG=",use_dynamic_MG
        write(TheUnit,"(6X,A,2E16.8,A1)") "a1 =",a1,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a2 =",a2,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a3 =",a3,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a4 =",a4,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a5 =",a5,"i"
        if( generate_bis ) then
            if (b1.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b1 =",b1,"i"
            if (b2.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b2 =",b2,"i"
            if (b3.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b3 =",b3,"i"
            if (b4.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b4 =",b4,"i"
            if (b5.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b5 =",b5,"i"
            if (b6.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b6 =",b6,"i"
            if (b7.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b7 =",b7,"i"
            if (b8.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b8 =",b8,"i"
            if (b9.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b9 =",b9,"i"
            if (b10.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "b10=",b10,"i"

            if (bzgs1.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzgs1 =",bzgs1,"i"
            if (bzgs2.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzgs2 =",bzgs2,"i"
            if (bzgs3.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzgs3 =",bzgs3,"i"
            if (bzgs4.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzgs4 =",bzgs4,"i"
            if (bzgs8.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzgs8 =",bzgs8,"i"

            if (bgsgs1.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bgsgs1 =",bgsgs1,"i"
            if (bgsgs2.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bgsgs2 =",bgsgs2,"i"
            if (bgsgs3.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bgsgs3 =",bgsgs3,"i"
            if (bgsgs4.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bgsgs4 =",bgsgs4,"i"
            if (bgsgs8.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bgsgs8 =",bgsgs8,"i"

            if (includeVprime) then
                if (IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2)) then
                    if(M_Zprime.gt.0d0) then
                      write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z' boson: mass=",M_Zprime*100d0,", width=",Ga_Zprime*100d0
                    else
                      write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z' boson: heavy mass limit (contact interaction)"
                    endif
                elseif (IsAWDecay(DecayMode1) .or. IsAWDecay(DecayMode2)) then
                    if(M_Wprime.gt.0d0) then
                      write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W' boson: mass=",M_Wprime*100d0,", width=",Ga_Wprime*100d0
                    else
                      write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W' boson: heavy mass limit (contact interaction)"
                    endif
                endif
            endif

            if (bzzp1.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp1 =",bzzp1,"i"
            if (bzzp2.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp2 =",bzzp2,"i"
            if (bzzp3.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp3 =",bzzp3,"i"
            if (bzzp4.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp4 =",bzzp4,"i"
            if (bzzp5.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp5 =",bzzp5,"i"
            if (bzzp6.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp6 =",bzzp6,"i"
            if (bzzp7.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp7 =",bzzp7,"i"
            if (bzzp8.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp8 =",bzzp8,"i"
            if (bzzp9.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp9 =",bzzp9,"i"
            if (bzzp10.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzzp10=",bzzp10,"i"

            if (bzpzp1.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp1 =",bzpzp1,"i"
            if (bzpzp2.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp2 =",bzpzp2,"i"
            if (bzpzp3.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp3 =",bzpzp3,"i"
            if (bzpzp4.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp4 =",bzpzp4,"i"
            if (bzpzp5.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp5 =",bzpzp5,"i"
            if (bzpzp6.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp6 =",bzpzp6,"i"
            if (bzpzp7.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp7 =",bzpzp7,"i"
            if (bzpzp8.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp8 =",bzpzp8,"i"
            if (bzpzp9.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp9 =",bzpzp9,"i"
            if (bzpzp10.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpzp10=",bzpzp10,"i"

            if (bzpgs1.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpgs1 =",bzpgs1,"i"
            if (bzpgs2.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpgs2 =",bzpgs2,"i"
            if (bzpgs3.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpgs3 =",bzpgs3,"i"
            if (bzpgs4.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpgs4 =",bzpgs4,"i"
            if (bzpgs8.ne.0) write(TheUnit,"(6X,A,2E16.8,A1)") "bzpgs8 =",bzpgs8,"i"
        else
            write(TheUnit,"(6X,A,2E16.8,A1)") "c1 =",c1,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c2 =",c2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c3 =",c3,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c41=",c41,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c42=",c42,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c5 =",c5,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c6 =",c6,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c7 =",c7,"i"
        endif
    endif
    write(TheUnit,"(6X,A,2F8.1,A1)") "Lambda=",Lambda*100d0


    write(TheUnit,"(4X,A)") ""
    write(TheUnit,"(4X,A)") "LHE output: "//trim(DataFile)//'.lhe'
    write(TheUnit,"(4X,A)") "Histogram output: "//trim(DataFile)//'.dat'
    write(TheUnit,"(4X,A)") "Log file: "//trim(DataFile)//'.log'
    if( ReadLHEFile .or. ConvertLHEFile ) write(TheUnit,"(4X,A)") "LHE input: "//trim(LHEProdFile)
    write(TheUnit,"(4X,A)") ""


END SUBROUTINE






SUBROUTINE WriteHisto(VG_Result,VG_Error,RunTime)
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
integer :: NBin,Hits,NHisto
real(8) :: BinSize,LowVal,BinVal,Value,Error,Integral
real(8),parameter :: ToGeV=1d2, ToPb=1d-3
real(8) :: VG_Result,VG_Error,RunTime


  write(io_HistoFile,"(A,2X,1F9.2,A)") "#"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "# column #1=histogram, #2=x-axis bin, #3=binned cross-section, #4=error, #5=number of events "
  write(io_HistoFile,"(A,2X,1F9.2,A)") "# this can be easily plotted in gnuplot with the commands:"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "#    NHi=3"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "#    plot 'output.dat' using ($2):($1==NHi ? $3 : 1/0) w st"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "# where NHi selects the histogram"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "#"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "#"
  write(io_HistoFile,"(A,2X,1F9.2,A)") "# run time =",RunTime/60d0,"min"
  write(io_HistoFile,"(A,2X,1PE20.10,2X,1PE20.5)") "# TotCS =",VG_Result,VG_Error
  write(io_HistoFile,"(A)") "#"
  do NHisto=1,NumHistograms
      write(io_HistoFile,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      Integral = 0d0
      BinSize = Histo(NHisto)%BinSize * Histo(NHisto)%SetScale
      LowVal  = Histo(NHisto)%LowVal  * Histo(NHisto)%SetScale
      do NBin=1, Histo(NHisto)%NBins
          BinVal = (LowVal+(NBin-1)*BinSize)
          Hits   = Histo(NHisto)%Hits(NBin)
          if( unweighted ) then
              it = 1
              Value  = Histo(NHisto)%Value(NBin)/BinSize / dble(EvalCounter)
              Integral = Integral + Histo(NHisto)%Value(NBin)
              if( Hits.gt.0 ) then
                  Error  = 1d0/dsqrt(dble(Hits))
              else
                  Error = 0d0
              endif
          else
              Value  = Histo(NHisto)%Value(NBin)/BinSize/it
              Integral = Integral + Histo(NHisto)%Value(NBin)/it
              Error  = 1d0/(BinSize)/it * dsqrt( Histo(NHisto)%Value2(NBin) - 1d0/it/ncall*Histo(NHisto)%Value(NBin)**2 )
          endif
          write(io_HistoFile,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
      enddo
      Integral = Integral + (Histo(NHisto)%Value(0)+Histo(NHisto)%Value(Histo(NHisto)%NBins+1))/it
      write(io_HistoFile,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(io_HistoFile,"(A)") "#"
  enddo

return
END SUBROUTINE



SUBROUTINE InitRandomSeed()
use modParameters
use modMisc
#if compiler==1
use ifport   !needed for getpid()
#endif
implicit none
integer, dimension(:), allocatable :: gen_seed
integer :: n,i,sclock
real :: tmp_real, tmp_real2
logical :: finished


    if (UserSeed.eq.0) then
#if compiler==1
        call random_seed()
#elif compiler==2
        call random_seed(size=n)
        allocate(gen_seed(n))
        call system_clock(count=sclock)
        gen_seed = sclock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = gen_seed)
        deallocate(gen_seed)
#endif
        call random_number(tmp_real)
        UserSeed = floor(tmp_real * 1000000000.)
        UserSeed = UserSeed + getpid()
    endif

    call random_seed(size=n)
    if( n.gt.nmaxseeds ) then
        print *, "Your compiler wants ", n, " seeds, but nmaxseeds=", nmaxseeds
        print *, "Try changing nmaxseeds and adding to DefaultSeeds (both in mod_Parameters)"
        stop 1
    endif

    call random_seed(put = DefaultSeeds(1:n))
    call random_number(tmp_real)
    do i=1,n   !Sometimes not all of the seeds requested matter.  Put UserSeed in place of the first one that does
        TheSeeds = DefaultSeeds
        TheSeeds(i) = UserSeed
        call random_seed(put = TheSeeds(1:n))
        call random_number(tmp_real2)
        if( abs(tmp_real2-tmp_real).ge.0.001 ) then
            exit
        endif
    enddo
    if( abs(tmp_real2-tmp_real).lt.0.001 ) then
        print *, "Fatal error: the random seed doesn't seem to make a difference!?!?"
        stop 1
    endif

    allocate(gen_seed(n))
    do i=1,n
        call random_number(tmp_real)
        gen_seed(i) = floor(tmp_real * 1000000000.)
    enddo
    call random_seed(put = gen_seed)
    deallocate(gen_seed)

return
END SUBROUTINE


SUBROUTINE PrintStatusBar(StatusPercent)
use modParameters
implicit none
integer :: StatusPercent
integer,save :: LastPercent=0

     if( mod(StatusPercent,10).eq.0 .and. LastPercent.ne.StatusPercent ) then
        write(io_stdout,"(X,I3,A)",advance='no') StatusPercent,"% "
        flush(io_stdout)
        LastPercent = StatusPercent
     endif

return
END SUBROUTINE


SUBROUTINE PrintStatusBar2(StatusPercent,ThePreface)
implicit none
integer :: StatusPercent
character(len=*):: ThePreface

      if( StatusPercent.lt.4d0*1   ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |------------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*2   ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |*-----------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*3 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-*----------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*4 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |--*---------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*5 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |---*--------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*6 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |----*-------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*7 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-----*------------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*8 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |------*-----------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*9 ) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-------*----------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*10) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |--------*---------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*11) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |---------*--------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*12) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |----------*-------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*13) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-----------*------------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*14) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |------------*-----------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*15) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-------------*----------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*16) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |--------------*---------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*17) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |---------------*--------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*18) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |----------------*-------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*19) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-----------------*------| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*20) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |------------------*-----| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*21) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-------------------*----| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*22) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |--------------------*---| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*23) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |---------------------*--| (",StatusPercent,"%)"
      elseif( StatusPercent.lt.4d0*24) then
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |----------------------*-| (",StatusPercent,"%)"
      else
          write(*,"(2X,A,A,I4,A)") trim(ThePreface)," |-----------------------*| (",StatusPercent,"%)"
      endif

return
END SUBROUTINE


SUBROUTINE PrintCommandLineArgs()
use modParameters
implicit none

        print *, ""
        print *, " help:                Print all command line options"
        print *, " DryRun:              Check that the command line is valid, then exit"
        print *, " Process configuration:"
        print *, "   Collider:          1=LHC (default), 2=Tevatron, 0=e+e-"
        print *, "   ColliderEnergy:    in TeV.  default is 13 TeV for LHC, 1.96 TeV for Tevatron,"
        print *, "                      250 GeV for e+e-"
        print *, "   Process:           0=spin-0, 1=spin-1, 2=spin-2 resonance,"
        print *, "                      50=qq/ee->VH, 51=gg->ZH,"
        print *, "                      60=weakVBF, 61=pp->Hjj, 62=pp->Hj,"
        print *, "                      66=VVHVV offshell, 67=VVVVbkg, 68=VVHVV+VVVV,"
        print *, "                      80=ttH, 90=bbH,"
        print *, "                      110=t+H t channel, 111=tbar+H t channel,"
        print *, "                      112=t+H s channel, 113=tbar+H s channel"
        print *, "                      114=t/tbar+H t/s channels"
        print *, "   DecayMode1:        decay mode for vector boson 1 (Z/W/gamma)"
        print *, "   DecayMode2:        decay mode for vector boson 2 (Z/W/gamma)"
        print *, "                        0=Z->2l,  1=Z->2q, 2=Z->2tau, 3=Z->2nu,"
        print *, "                        4=W->lnu, 5=W->2q, 6=W->taunu,"
        print *, "                        7=gamma, 8=Z->2l+2tau,"
        print *, "                        9=Z->anything, 10=W->lnu+taunu, 11=W->anything"
        print *, "   Interf:            0=neglect interference for 4f final states,"
        print *, "                      1=include interference"
        print *, "   ReweightInterf:    if true, include interference as LHE event weights for"
        print *, "                      offshell VBF events"
        print *, "   RandomizeVVdecays: Randomizes the order of DecayMode1 and DecayMode2,"
        print *, "                      per event (default true)"
        print *, "                      For a WW decay, turning this off will mean"
        print *, "                      DecayMode1 is W+ and DecayMode2 is W-"
        print *, "   PChannel:          0=g+g, 1=q+qb, 2=both"
        print *, "   ChannelRatio:      ratio of qqb / (qqb + gg), for Process=2 PChannel=2"
        print *, "                       default is to allow this ratio to come from the couplings"
        print *, "   PDFSet:            1=CTEQ6L1(2001),  2=MSTW(2008),"
        print *, "                      2xx=MSTW with eigenvector set xx=01..40,"
        print *, "                      3=NNPDF3.0LO"
        print *, "                      (only valid if not interfaced with LHAPDF)"
        print *, "   LHAPDF:            info file to use if interfaced with LHAPDF"
        print *, "                      (example: NNPDF30_lo_as_0130/NNPDF30_lo_as_0130.info)"
        print *, "   LHAPDFMem:         member number in LHAPDF set"
        print *, "   epPolarization:    Polarization of e+ for e+e- collider"
        print *, "   emPolarization:    Polarization of e- for e+e- collider"
        print *, "                        0:      no polarization"
        print *, "                        +/-100: helicity=+/-1"
        print *, "   TopDK:             For ttH or t+H, 0=leave top quarks as stable, 1=decay top quarks"
        print *, "   TauDK:             In ReadLHE mode, specify this option as either 0 or 1"
        print *, "                      to decay H->tautau.  If it is 0, the taus are written as"
        print *, "                      stable; if it is 1, they decay to Wnu, with the W's decaying"
        print *, "                      according to DecayModes1,2."
        print *, "   HbbDK:             For VH production, decay H->bb"
        print *, "   VH_PC:             VH partonic channel and mode selection"
        print *, "                      ee ( = e+ e- @LO)"
        print *, "                      gg ( = triangles + boxes of gg)"
        print *, "                      qq ( = q q~ @LO)"
        print *, "                      lo ( = q q~ @LO)"
        print *, "                      tr ( = triangles of gg)"
        print *, "                      bo ( = boxes of gg)"
        print *, "                      in ( = interference = 2*dble(box*dconjg(triangle)) of gg)"
        print *, "                      qg ( = real - dipoles, for g q/q~ > VH + q/q~, for development only)"
        print *, "                      gq ( = K + P, for g q/q~ > VH + q/q~, for development only)"
        print *, "                      sb ( = real - dipoles, for q q~ @NLO, for development only)"
        print *, "                      sp ( = virtual + I + K + P, for q q~ @NLO, for development only)"
        print *, "                      nl ( = NLO = q q~ @LO + NLO + gq)"
        print *, "                      VH_PC overrides Pchannel."
        print *, "   alpha_dip          extra non-physical degree of freedom for Process=51 & VH_PC=nl, defaulted at 1."
        print *, "                      Vary to check indepedence (of alpha_dip)."
        print *, "   VBFoffsh_run:      For VBF offshell production, set this to a number from 1-5"
        print *, "                      for each of the 5 jobs.  See manual for more details."
        print *, " Resonance parameters:"
        print *, "   MReso:             resonance mass in GeV (default=125.00)"
        print *, "   GaReso:            resonance width in GeV (default=0.00407)"
        print *, "   ctauReso:          resonance decay length in mm (default=0)"
        print *, "   OffshellX:         Whether to allow resonance (X) to go offshell"
        print *, "                      in processes 0, 1 or 2"
        print *, "   MReso2:            2nd resonance mass in GeV in offshell VBF"
        print *, "   GaReso2:           2nd resonance width in GeV in offshell VBF"
        print *, " EW coupling parameters:"
        print *, "   Vud:               CKM element for W-ud couplings"
        print *, "   Vus:               CKM element for W-us couplings"
        print *, "   Vub:               CKM element for W-ub couplings"
        print *, "   Vcd:               CKM element for W-cd couplings"
        print *, "   Vcs:               CKM element for W-cs couplings"
        print *, "   Vcb:               CKM element for W-cb couplings"
        print *, "   Vtd:               CKM element for W-td couplings"
        print *, "   Vts:               CKM element for W-ts couplings"
        print *, "   Vtb:               CKM element for W-tb couplings"
        print *, " Cuts:"
        print *, "   pTjetcut:          Minimum pT for jets in GeV (default: 15)"
        print *, "   deltaRcut:         Minimum deltaR for jets (default: 0.3)"
        print *, "   mJJcut:            Minimum dijet mass in GeV (default: 0)"
        print *, "   MPhotonCutoff:     Minimum mass for offshell photons in GeV, when included (default: 4)"
        print *, "   etajetcut:         Maximum |eta| for jets in offshell VBF (default: 4)"
        print *, "   detajetcut:        Minimum deltaeta between jets in offshell VBF (default: 2)"
        print *, "   JetsOppositeEta:   Require sgn(eta) to be opposite for the two jets in offshell VBF"
        print *, "                      (default: true)"
        print *, "   pTlepcut:          Minimum pT for leptons in offshell VBF, in GeV (default: 3)"
        print *, "   etalepcut:         Maximum |eta| for leptons in offshell VBF (default: 2.7)"
        print *, "   m4l_min, m4l_max:  Minimum and maximum four-lepton mass in offshell VBF"
        print *, "   m2l_min:    Minimum invariant mass of V (onshell) in new VH (\texttt{Process=51}) (default: 0)"
        print *, "   m2l_max:   Maximum invariant mass of V (onshell) in newVH (\texttt{Process=51}) (default: infinity)"
        print *, "   mVH_min:   Minimum invariant mass of VH in new VH (\texttt{Process=51}) (default: 0)"
        print *, "   mVH_max:   Maximum invariant mass of VH in new VH (\texttt{Process=51}) (default: infinity)"
        print *, " Renormalization and factorization scales:"
        print *, "   FacScheme:         PDF factorization scale scheme"
        print *, "   MuFacMultiplier:   Multiplier for the factorization scale chosen by FacScheme"
        print *, "   RenScheme:         QCD renormalization scale scheme"
        print *, "   MuRenMultiplier:   Multiplier for the renormalization scale chosen by RenScheme"
        print *, " Lepton and jet filter:"
        print *, "   FilterNLept:       For decay mode, reject events that have less than FilterNLept leptons"
        print *, "   FilterOSPairs:     For decay mode, reject events that have less than FilterOSPairs pairs of"
        print *, "                      sign leptons of any flavor."
        print *, "   FilterOSSFPairs:   For decay mode, reject events that have less than FilterOSSFPairs pairs of"
        print *, "                      opposite-sign-same-flavor leptons."
        print *, "   CountTauAsAny:     For FilterOSSFPairs, taus can stand in place of electrons or muons"
        print *, "                      of the same charge."
        print *, "   FilterNJets:       For decay mode, reject events that have less than FilterNJets quarks"
        print *, "                      and/or gluons"
        print *, "   WriteFailedEvents: Write events that fail in the LHE file, but with a weight of 0"
        print *, "                      (off by default)"
        print *, " Higgs propagator and decay width:"
        print *, "   WidthScheme:       Higgs width scheme: 1 for running width, 2 for fixed width (default),"
        print *, "                      and 3 for the CPS"
        print *, "   WidthSchemeIn:     For decay mode, reweight from one propagator to another by setting"
        print *, "                      WidthScheme and WidthSchemeIn to different values"
        print *, "   ReweightDecay:     For decay mode, reweight input decay by the decay probability"
        print *, "   PmHstarEvals:      For ReweightDecay, number of evaluations per mass point (default: 200000)"
        print *, "   ReadPmHstar:       For ReweightDecay, read the decay probability distribution from a file"
        print *, "   PmHstarFile:       File to write and read the decay probability distribution"
        print *, " Statistics options:"
        print *, "   VegasNc0:          number of evaluations for integrand scan"
        print *, "   VegasNc1:          number of evaluations for accept-reject sampling"
        print *, "   VegasNc2:          number of events for accept-reject sampling"
        print *, "   ReadCSmax:         Read the results of the grid generation step from a file"
        print *, "   CSmaxFile:         File to use for reading (if ReadCSmax is set) or writing (otherwise)"
        print *, "                      the results of the grid generation step.  Depending on the process,"
        print *, "                      suffixes are appended to this base name. (default: DataFile without .lhe)"
        print *, "   Seed:              Random seed for event generation"
        print *, " I/O options:"
        print *, "   Unweighted:        0=weighted events, 1=unweighted events"
        print *, "   WriteWeightedLHE:  For Unweighted=0, write weighted events to an LHE file"
        print *, "                      (note that the output could be huge)"
        print *, "   DataFile:          LHE output file"
        print *, "   ReadLHE:           LHE input file from external file (only spin-0)"
        print *, "   ConvertLHE:        Convert decay of the V from VH production."
        print *, "                      Use DecayMode1 to specify the decay."
        print *, "                      (should be a Z or W mode, depending on the input file)"
        print *, "   UnformattedRead:   Turn this on if the normal, faster reading fails"
        print *, " Couplings:"
        print *, "   See manual for the full list"
        print *, ""

END SUBROUTINE

