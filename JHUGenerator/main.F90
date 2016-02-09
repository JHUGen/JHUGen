! If you use this program please cite Phys.Rev. D81 (2010) 075022; arXiv:1001.3396 [hep-ph],
!                                     Phys.Rev. D86 (2012) 095031; arXiv:1208.4018 [hep-ph],
!                                 and Phys.Rev. D89 (2014) 035007; arXiv:1309.4819 [hep-ph].
PROGRAM JHUGenerator
use ModParameters
use ModKinematics
use ModCrossSection
#if compiler==1
use ifport
#endif
implicit none
real(8) :: VG_Result,VG_Error

   call GetCommandlineArgs()
   call InitProcessScaleSchemes()
   call InitPDFs()!  
   call InitHisto()
   call InitParameters()
   call InitProcess()
   call InitVegas()
   call InitRandomSeed()
   call OpenFiles()
   call PrintLogo(io_stdout)
   call PrintLogo(io_LogFile)
   call WriteParameters(io_stdout)
   call WriteParameters(io_LogFile)
   if ( .not. ReadLHEFile .and. .not. ConvertLHEFile .and. .not.((Process.eq.60 .or. Process.eq.61) .and. unweighted) ) then
      call InitOutput(1d0, 1d14)   !for VBF/HJJ the cross section is calculated, so use that in the <init> block
   endif
   write(io_stdout,*) " Running"
   if( ConvertLHEFile ) then
        call StartConvertLHE(VG_Result,VG_Error)
   elseif( ReadLHEFile ) then
        call StartReadLHE_NEW(VG_Result,VG_Error)
   elseif( CalcPMZZ ) then
        call GetMZZdistribution()
   else
        if( Process.eq.80 .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.66 .or. Process.eq.90 .or. &
            Process.eq.110 .or. Process.eq.111 .or. Process.eq.112 .or. Process.eq.113 ) then
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
! !       Below, J stands for either the top in JJ associated production (e.g. tqH), or the single jet (Hj)
! !       Scheme=+-kRenFacScheme_mjhstar: Scale ~ m_JHstar (or m_JVV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mj_mhstar: Scale ~ m_J + m_Hstar (or m_J+m_VV in bkg.); + for running, - for fixed
! !       Scheme=+-kRenFacScheme_mj: Scale ~ m_J for the heavy jet (ie. t or b); + for running, - for fixed
subroutine InitProcessScaleSchemes() ! If schemes are set to default, reset to the appropriate numbers
   use ModParameters
   use ModMisc
   implicit none

      if( .not.                 &
         (                      &
            Process.eq. 0 .or.  &
            Process.eq. 1 .or.  &
            Process.eq. 2 .or.  &
            Process.eq.50 .or.  &
            Process.eq.60 .or.  &
            Process.eq.61 .or.  &
            Process.eq.62 .or.  &
            Process.eq.66 .or.  &
            Process.eq.80 .or.  &
            Process.eq.90 .or.  &
            Process.eq.110 .or. &
            Process.eq.111 .or. &
            Process.eq.112 .or. &
            Process.eq.113      &
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
         Process.eq.50 & !- VHiggs
         ) then
            FacScheme = -kRenFacScheme_mhstar
            MuFacMultiplier = 1d0
         elseif( &
         Process.eq.66      & !- HVBF with decays
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
         Process.eq.113      & !- tb+H s-channel
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
         Process.eq.50 & !- VHiggs
         ) then
            RenScheme = -kRenFacScheme_mhstar
            MuRenMultiplier = 1d0
         elseif( &
         Process.eq.66     & !- HVBF with decays
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
         Process.eq.113      & !- tb+H s-channel
         ) then
            RenScheme = -kRenFacScheme_mj_mhstar
            MuRenMultiplier = 0.25d0
         endif
      endif

      ! H+2j MEs
      if( &
         (                     &
            Process.eq.50 .or. &
            Process.eq.60 .or. &
            Process.eq.61 .or. &
            Process.eq.66 .or. &
            Process.eq.80 .or. &
            Process.eq.90      &
         ) .and. (             &
            (abs(FacScheme).eq.kRenFacScheme_mjhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj) .or. &
            (abs(RenScheme).eq.kRenFacScheme_mjhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj)      &
         )                     &
      ) call Error("ttH, bbH, HJJ, VBF and VH processes cannot distinguish the outgoing partons. Choose a different renormalization or factorization scheme.")

      if( &
         (                     &
            Process.eq.50 .or. &
            Process.eq.60 .or. &
            Process.eq.61 .or. &
            Process.eq.66      &
         ) .and. (             &
            (abs(FacScheme).eq.kRenFacScheme_mj_mj_mhstar) .or. (abs(FacScheme).eq.kRenFacScheme_mj_mj) .or. &
            (abs(RenScheme).eq.kRenFacScheme_mj_mj_mhstar) .or. (abs(RenScheme).eq.kRenFacScheme_mj_mj)      &
         )                     &
      ) call Error("HJJ, VBF and VH processes outgoing partons are mostly massless, and alpha_S at a scale ~0 GeV is very unstable. Choose a different renormalization or factorization scheme (e.g. kRenFacScheme_mhstar).")

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



SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
use ModMisc
implicit none
character :: arg*(500)
integer :: NumArgs,NArg,OffShell_XVV
logical :: help, success, SetLastArgument, interfSet
logical :: SetAnomalousSpin0gg, Setghg2, SetAnomalousSpin0ZZ, Setghz1
logical :: SetZgammacoupling, Setgammagammacoupling
logical :: SetAnomalousSpin1qq, Setspin1qqleft, Setspin1qqright, SetAnomalousSpin1ZZ, Set1minus
logical :: SetAnomalousSpin2gg, Seta2, SetAnomalousSpin2qq, Setspin2qqleft, Setspin2qqright,SetAnomalousSpin2ZZ, Setb2
logical :: SetAnomalousHff, Setkappa

   help = .false.

   Collider=1
   VegasIt1=-1
#if useLHAPDF==1
   LHAPDFString = ""
   LHAPDFMember = 0
#else
   PDFSet=1      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40
#endif
   VegasNc0=-1
   VegasNc1=-1
   VegasNc2=-1
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

   WidthScheme=0
   WidthSchemeIn=0
   TopDecays=-1
   TauDecays=-1
   Process = 0   ! select 0, 1 or 2 to represent the spin of the resonance
   Unweighted =.true.
   OffShell_XVV=011! 000: X,V1,V2 on-shell; 010: X,V2 on-shell, V1 off-shell; and so on
   LHEProdFile=""
   ReadLHEFile=.false.
   ConvertLHEFile=.false.
   ReadCSmax=.false.
   CalcPMZZ = .false.
   GenerateEvents=.false.
   RequestNLeptons = -1
   RequestOS=-1
   RequestOSSF=-1
   CountTauAsAny = .true.
   interfSet = .false.
   WriteFailedEvents=0

   MuFacMultiplier = 1d0
   MuRenMultiplier = 1d0
   FacScheme = kRenFacScheme_default
   RenScheme = kRenFacScheme_default

   SetAnomalousSpin0gg=.false.
   Setghg2=.false.
   SetAnomalousSpin0ZZ=.false.
   Setghz1=.false.
   SetZgammacoupling=.false.
   Setgammagammacoupling=.false.
   SetAnomalousSpin1qq=.false.
   Setspin1qqleft=.false.
   Setspin1qqright=.false.
   SetAnomalousSpin1ZZ=.false.
   Set1minus=.false.
   SetAnomalousSpin2gg=.false.
   Seta2=.false.
   SetAnomalousSpin2qq=.false.
   Setspin2qqleft=.false.
   Setspin2qqright=.false.
   SetAnomalousSpin2ZZ=.false.
   Setb2=.false.
   SetAnomalousHff=.false.
   Setkappa=.false.

   DataFile="./data/output"


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
    call ReadCommandLineArgument(arg, "Collider", success, Collider)
#if useLHAPDF==1
    call ReadCommandLineArgument(arg, "LHAPDF", success, LHAPDFString)
    call ReadCommandLineArgument(arg, "LHAPDFMem", success, LHAPDFMember)
#else
    call ReadCommandLineArgument(arg, "PDFSet", success, PDFSet)
#endif
    call ReadCommandLineArgument(arg, "MReso", success, M_Reso, SetLastArgument)
    if( SetLastArgument ) M_Reso = M_Reso*GeV
    call ReadCommandLineArgument(arg, "GaReso", success, Ga_Reso, SetLastArgument)
    if( SetLastArgument ) Ga_Reso = Ga_Reso*GeV
    call ReadCommandLineArgument(arg, "VegasNc0", success, VegasNc0)
    call ReadCommandLineArgument(arg, "VegasNc1", success, VegasNc1)
    call ReadCommandLineArgument(arg, "VegasNc2", success, VegasNc2)
    call ReadCommandLineArgument(arg, "PChannel", success, PChannel)
    call ReadCommandLineArgument(arg, "DataFile", success, DataFile)
    call ReadCommandLineArgument(arg, "Process", success, Process)
    call ReadCommandLineArgument(arg, "DecayMode1", success, DecayMode1)
    call ReadCommandLineArgument(arg, "DecayMode2", success, DecayMode2)
    call ReadCommandLineArgument(arg, "FacScheme", success, FacScheme)
    call ReadCommandLineArgument(arg, "RenScheme", success, RenScheme)
    call ReadCommandLineArgument(arg, "MuFacMultiplier", success, MuFacMultiplier)
    call ReadCommandLineArgument(arg, "MuRenMultiplier", success, MuRenMultiplier)
    call ReadCommandLineArgument(arg, "TopDK", success, TopDecays)
    call ReadCommandLineArgument(arg, "TauDK", success, TauDecays)
    call ReadCommandLineArgument(arg, "ReweightDecay", success, ReweightDecay)
    call ReadCommandLineArgument(arg, "WidthScheme", success, WidthScheme)
    call ReadCommandLineArgument(arg, "WidthSchemeIn", success, WidthSchemeIn)
    call ReadCommandLineArgument(arg, "OffXVV", success, OffShell_XVV)
    call ReadCommandLineArgument(arg, "FilterNLept", success, RequestNLeptons)
    call ReadCommandLineArgument(arg, "FilterOSPairs", success, RequestOS)
    call ReadCommandLineArgument(arg, "FilterOSSFPairs", success, RequestOSSF)
    call ReadCommandLineArgument(arg, "CountTauAsAny", success, CountTauAsAny)
    call ReadCommandLineArgument(arg, "Unweighted", success, Unweighted)
    call ReadCommandLineArgument(arg, "Interf", success, includeInterference, success2=interfSet)
    call ReadCommandLineArgument(arg, "ReadLHE", success, LHEProdFile, success2=ReadLHEFile)
    call ReadCommandLineArgument(arg, "ConvertLHE", success, LHEProdFile, success2=ConvertLHEFile)
    call ReadCommandLineArgument(arg, "ReadCSmax", success, ReadCSmax)
    call ReadCommandLineArgument(arg, "GenEvents", success, GenerateEvents, SetLastArgument)
    if( SetLastArgument ) Unweighted = .false.
    call ReadCommandLineArgument(arg, "CalcPMZZ", success, CalcPMZZ)
    call ReadCommandLineArgument(arg, "WriteFailedEvents", success, WriteFailedEvents)
    call ReadCommandLineArgument(arg, "Seed", success, UserSeed)
    call ReadCommandLineArgument(arg, "WriteGit", success, writegit) !for testing purposes

    !anomalous couplings
    !If any anomalous couplings are set, the default ones have to be set explicitly to keep them on or turn them off
    !e.g. just setting ghz4=0.2982,0 is ambiguous if you mean to leave g1 on (so fa3=0.5)
    !                                                      or to turn it off (fa3=1 with a weird prefactor)
    !spin 0 gg couplings
    call ReadCommandLineArgument(arg, "ghg2", success, ghg2, success2=SetAnomalousSpin0gg, success3=Setghg2)
    call ReadCommandLineArgument(arg, "ghg3", success, ghg3, success2=SetAnomalousSpin0gg)
    call ReadCommandLineArgument(arg, "ghg4", success, ghg4, success2=SetAnomalousSpin0gg)

    !spin 0 ZZ couplings
    call ReadCommandLineArgument(arg, "ghz1", success, ghz1, success2=SetAnomalousSpin0ZZ, success3=Setghz1)
    call ReadCommandLineArgument(arg, "ghz2", success, ghz2, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3", success, ghz3, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4", success, ghz4, success2=SetAnomalousSpin0ZZ)

    !spin 0 Zgamma couplings
    call ReadCommandLineArgument(arg, "ghzgs2", success, ghzgs2, success2=SetZgammacoupling)
    call ReadCommandLineArgument(arg, "ghzgs3", success, ghzgs3, success2=SetZgammacoupling)
    call ReadCommandLineArgument(arg, "ghzgs4", success, ghzgs4, success2=SetZgammacoupling)

    !spin 0 gammagamma couplings
    call ReadCommandLineArgument(arg, "ghgsgs2", success, ghgsgs2, success2=Setgammagammacoupling)
    call ReadCommandLineArgument(arg, "ghgsgs3", success, ghgsgs3, success2=Setgammagammacoupling)
    call ReadCommandLineArgument(arg, "ghgsgs4", success, ghgsgs4, success2=Setgammagammacoupling)

    !spin 0 ZZ momentum dependent couplings
    call ReadCommandLineArgument(arg, "ghz1_prime", success, ghz1_prime, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz1_prime2", success, ghz1_prime2, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz1_prime3", success, ghz1_prime3, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz1_prime4", success, ghz1_prime4, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz1_prime5", success, ghz1_prime5, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz1_prime6", success, ghz1_prime6, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz1_prime7", success, ghz1_prime7, success2=SetAnomalousSpin0ZZ)

    call ReadCommandLineArgument(arg, "ghz2_prime", success, ghz2_prime, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz2_prime2", success, ghz2_prime2, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz2_prime3", success, ghz2_prime3, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz2_prime4", success, ghz2_prime4, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz2_prime5", success, ghz2_prime5, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz2_prime6", success, ghz2_prime6, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz2_prime7", success, ghz2_prime7, success2=SetAnomalousSpin0ZZ)

    call ReadCommandLineArgument(arg, "ghz3_prime", success, ghz3_prime, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3_prime2", success, ghz3_prime2, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3_prime3", success, ghz3_prime3, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3_prime4", success, ghz3_prime4, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3_prime5", success, ghz3_prime5, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3_prime6", success, ghz3_prime6, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz3_prime7", success, ghz3_prime7, success2=SetAnomalousSpin0ZZ)

    call ReadCommandLineArgument(arg, "ghz4_prime", success, ghz4_prime, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4_prime2", success, ghz4_prime2, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4_prime3", success, ghz4_prime3, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4_prime4", success, ghz4_prime4, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4_prime5", success, ghz4_prime5, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4_prime6", success, ghz4_prime6, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "ghz4_prime7", success, ghz4_prime7, success2=SetAnomalousSpin0ZZ)

    !spin 0 Zgamma momentum dependent coupling
    call ReadCommandLineArgument(arg, "ghzgs1_prime2", success, ghzgs1_prime2, success2=SetZgammacoupling)

    ! Sign of q1,2,12**2 for the Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
    call ReadCommandLineArgument(arg, "cz_q1sq", success, cz_q1sq, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "cz_q2sq", success, cz_q1sq, success2=SetAnomalousSpin0ZZ)
    call ReadCommandLineArgument(arg, "cz_q12sq", success, cz_q1sq, success2=SetAnomalousSpin0ZZ)

    !spin 0 WW couplings
    call ReadCommandLineArgument(arg, "ghw1", success, ghw1, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2", success, ghw2, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3", success, ghw3, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4", success, ghw4, success2=distinguish_HWWcouplings)

    call ReadCommandLineArgument(arg, "ghw1_prime", success, ghw1_prime, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw1_prime2", success, ghw1_prime2, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw1_prime3", success, ghw1_prime3, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw1_prime4", success, ghw1_prime4, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw1_prime5", success, ghw1_prime5, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw1_prime6", success, ghw1_prime6, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw1_prime7", success, ghw1_prime7, success2=distinguish_HWWcouplings)

    call ReadCommandLineArgument(arg, "ghw2_prime", success, ghw2_prime, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2_prime2", success, ghw2_prime2, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2_prime3", success, ghw2_prime3, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2_prime4", success, ghw2_prime4, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2_prime5", success, ghw2_prime5, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2_prime6", success, ghw2_prime6, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw2_prime7", success, ghw2_prime7, success2=distinguish_HWWcouplings)

    call ReadCommandLineArgument(arg, "ghw3_prime", success, ghw3_prime, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3_prime2", success, ghw3_prime2, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3_prime3", success, ghw3_prime3, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3_prime4", success, ghw3_prime4, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3_prime5", success, ghw3_prime5, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3_prime6", success, ghw3_prime6, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw3_prime7", success, ghw3_prime7, success2=distinguish_HWWcouplings)

    call ReadCommandLineArgument(arg, "ghw4_prime", success, ghw4_prime, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4_prime2", success, ghw4_prime2, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4_prime3", success, ghw4_prime3, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4_prime4", success, ghw4_prime4, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4_prime5", success, ghw4_prime5, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4_prime6", success, ghw4_prime6, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "ghw4_prime7", success, ghw4_prime7, success2=distinguish_HWWcouplings)

    ! Sign of q1,2,12**2 for the Lambda's, set to 1 or -1 to get q**2-dependence from these form factor Lambdas
    call ReadCommandLineArgument(arg, "cw_q1sq", success, cw_q1sq, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "cw_q2sq", success, cw_q1sq, success2=distinguish_HWWcouplings)
    call ReadCommandLineArgument(arg, "cw_q12sq", success, cw_q1sq, success2=distinguish_HWWcouplings)

    !spin 1
    call ReadCommandLineArgument(arg, "zprime_qq_left", success, zprime_qq_left, success2=SetAnomalousSpin1qq, success3=Setspin1qqleft)
    call ReadCommandLineArgument(arg, "zprime_qq_right", success, zprime_qq_right, success2=SetAnomalousSpin1qq, success3=Setspin1qqright)
    call ReadCommandLineArgument(arg, "zprime_zz_1", success, zprime_zz_1, success2=SetAnomalousSpin1ZZ, success3=Set1minus)
    call ReadCommandLineArgument(arg, "zprime_zz_2", success, zprime_zz_2, success2=SetAnomalousSpin1ZZ)

    !spin 2
    call ReadCommandLineArgument(arg, "a1", success, a1, success2=SetAnomalousSpin2gg)
    call ReadCommandLineArgument(arg, "a2", success, a2, success2=SetAnomalousSpin2gg, success3=Seta2)
    call ReadCommandLineArgument(arg, "a3", success, a3, success2=SetAnomalousSpin2gg)
    call ReadCommandLineArgument(arg, "a4", success, a4, success2=SetAnomalousSpin2gg)
    call ReadCommandLineArgument(arg, "a5", success, a5, success2=SetAnomalousSpin2gg)
    call ReadCommandLineArgument(arg, "graviton_qq_left", success, graviton_qq_left, success2=SetAnomalousSpin2qq, success3=Setspin2qqleft)
    call ReadCommandLineArgument(arg, "graviton_qq_right", success, graviton_qq_right, success2=SetAnomalousSpin2qq, success3=Setspin2qqright)
    call ReadCommandLineArgument(arg, "b1", success, b1, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b2", success, b2, success2=SetAnomalousSpin2ZZ, success3=Setb2)
    call ReadCommandLineArgument(arg, "b3", success, b3, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b4", success, b4, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b5", success, b5, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b6", success, b6, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b7", success, b7, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b8", success, b8, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b9", success, b9, success2=SetAnomalousSpin2ZZ)
    call ReadCommandLineArgument(arg, "b10", success, b10, success2=SetAnomalousSpin2ZZ)

    !Hff couplings
    call ReadCommandLineArgument(arg, "kappa", success, kappa, success2=SetAnomalousHff, success3=Setkappa)
    call ReadCommandLineArgument(arg, "kappa_tilde", success, kappa_tilde, success2=SetAnomalousHff)

    !cuts
    call ReadCommandLineArgument(arg, "pTjetcut", success, pTjetcut, SetLastArgument)
    if( SetLastArgument ) pTjetcut = pTjetcut*GeV
    call ReadCommandLineArgument(arg, "deltaRcut", success, Rjet)
    call ReadCommandLineArgument(arg, "mJJcut", success, mJJcut, SetLastArgument)
    if( SetLastArgument ) mJJcut = mJJcut*GeV
    call ReadCommandLineArgument(arg, "VBF_m4l_min", success, m4l_minmax(1), SetLastArgument)
    if( SetLastArgument ) m4l_minmax(1) = m4l_minmax(1)*GeV
    call ReadCommandLineArgument(arg, "VBF_m4l_max", success, m4l_minmax(2), SetLastArgument)
    if( SetLastArgument ) m4l_minmax(2) = m4l_minmax(2)*GeV
    call ReadCommandLineArgument(arg, "MPhotonCutoff", success, MPhotonCutoff, SetLastArgument)
    if( SetLastArgument ) MPhotonCutoff = MPhotonCutoff*GeV

    if( .not.success ) then
        call Error("Unknown command line argument: " // trim(arg))
    endif
   enddo

    !================================
    !Command line argument processing
    !================================

    !PChannel

    if (Process.eq.0) PChannel = 0   !only gluons
    if (Process.eq.1 .or. Process.eq.50 .or. Process.eq.60 .or. Process.eq.66) PChannel = 1   !only quarks

    !OffXVV

    if(OffShell_XVV.ge.100) then
        OffShell_XVV=OffShell_XVV-100
        OffShellReson=.true.
    else
        OffShellReson=.false.
    endif
    if(OffShell_XVV.ge.10) then
        OffShell_XVV=OffShell_XVV-10
        OffShellV1=.true.
    else
        OffShellV1=.false.
    endif
    if(OffShell_XVV.ge.1) then
        OffShell_XVV=OffShell_XVV-1
        OffShellV2=.true.
    else
        OffShellV2=.false.
    endif

    !LHAPDF

#if useLHAPDF==1
    if( LHAPDFString.eq."" ) then
       print *, "Need to specify pdf file name in command line argument LHAPDF"
       stop 1
    endif
    call GET_ENVIRONMENT_VARIABLE("LHAPDF_DATA_PATH", LHAPDF_DATA_PATH)
#endif    
    
    !Renormalization/factorization schemes

    if((abs(FacScheme) .ge. nRenFacSchemes) .or. (abs(RenScheme) .ge. nRenFacSchemes) .or. (MuFacMultiplier.le.0d0) .or. (MuRenMultiplier.le.0d0)) call Error("The renormalization or factorization scheme is invalid, or the scale multiplier to either is not positive.")

    !ReadLHE and ConvertLHE
    !MUST HAPPEN BEFORE DETERMINING INTERFERENCE
    !so that the mass can be read

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

    !DecayModes, interference, photon couplings, ...

    if( ConvertLHEFile ) then
       DecayMode2 = DecayMode1
    endif 

    if( Process.eq.50 ) then
        DecayMode2=DecayMode1
        if( Collider.eq.2 ) call Error("Collider 2 not available for VH")
        if( (IsAZDecay(DecayMode1).eqv..false.) .and. (Collider.ne.1) ) call Error("WH with Collider 1 only")
    endif

    if( Process.ge.110 .and. Process.le.113 ) DecayMode2 = DecayMode1
    
    if( (TopDecays.ne.0 .and. TopDecays.ne.1) .and. (Process.eq.80 .or. (Process.ge.110 .and. Process.le.113)) ) call Error("Specify TopDK=0,1")
    if( (TopDecays.eq.1) .and. .not. IsAWDecay(DecayMode1) ) call Error("Invalid DecayMode1 for top decays")
    if( (TopDecays.eq.1) .and. .not. IsAWDecay(DecayMode2) ) call Error("Invalid DecayMode2 for top decays")

    if( (TauDecays.eq.1) .and. .not. IsAWDecay(DecayMode1) ) call Error("Invalid DecayMode1 for tau decays")
    if( (TauDecays.eq.1) .and. .not. IsAWDecay(DecayMode2) ) call Error("Invalid DecayMode2 for tau decays")

    if( IsAZDecay(DecayMode1) ) then
       M_V = M_Z
       Ga_V= Ga_Z
    elseif( IsAWDecay(DecayMode1) ) then
       M_V = M_W
       Ga_V= Ga_W    
    elseif( IsAPhoton(DecayMode1) ) then
       M_V = 0d0
       Ga_V= 0d0    
    endif

    if( (DecayMode1.eq.DecayMode2 .and. IsAZDecay(DecayMode1)) .or.  &
        (DecayMode1.eq.9) .or. (DecayMode2.eq.9)               .or.  &
        (DecayMode1.eq.8  .and. DecayMode2.eq.0)               .or.  &
        (DecayMode1.eq.8  .and. DecayMode2.eq.2)               .or.  &
        (DecayMode1.eq.0  .and. DecayMode2.eq.8)               .or.  &
        (DecayMode1.eq.2  .and. DecayMode2.eq.8)               ) then !  allow interference
            if( .not.interfSet ) then!  set default interference switch
                if( M_Reso.gt.2d0*M_Z ) then
                    includeInterference = .false.
                else
                    includeInterference = .true.
                endif
            endif
    else
        includeInterference = .false.   ! no interference if decay mode does not allow 4 same flavor leptons
    endif

    if( IsAZDecay(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
        includeGammaStar = (SetZgammacoupling .or. Setgammagammacoupling)
    elseif( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
        includeGammaStar = Setgammagammacoupling
    endif

    if( IsAZDecay(DecayMode1) .and. IsAWDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif

!     if( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
!        print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
!        stop
!     endif

    if( IsAWDecay(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif

    if( (IsAWDecay(DecayMode1) .and. IsAPhoton(DecayMode2)) .or. (IsAPhoton(DecayMode1) .and. IsAWDecay(DecayMode2)) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif

    if( IsAPhoton(DecayMode1) .and. (IsAZDecay(DecayMode2) .or. IsAWDecay(DecayMode2)) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       print *, " Please try swapping the decay modes."
       stop 1
    endif

    if( IsAPhoton(DecayMode2) .and. (IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1)) ) then! require OffXVV=010 for Z+photon
       if( .not. ((.not.OffShellReson) .and. OffShellV1 .and. (.not.OffShellV2)) ) then
          print *, "OffXVV has to be 010 for Z+photon production"
          stop 1
       endif
    endif

    if( IsAPhoton(DecayMode2) .and. (.not.IsAPhoton(DecayMode1) .and. .not. IsAZDecay(DecayMode1)) ) then
       print *, "DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif

    if( IsAPhoton(DecayMode1) .and. (OffShellV1 .or. OffShellV2) ) then
       print *, "Photons have to be on-shell."
       stop 1
    endif

    if( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) .and. .not.SetZgammacoupling .and. .not.Setgammagammacoupling ) then
        print *, "To decay to Zgamma, you need to set one of the HZgamma (ghzgs*) or Hgammagamma (ghgsgs*) couplings."
        stop 1
    endif

    if( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) .and. .not.Setgammagammacoupling ) then
        print *, "To decay to gammagamma, you need to set one of the Hgammagamma couplings (ghgsgs*)."
        stop 1
    endif

    if( (DecayMode1.ge.12) .or. (DecayMode2.ge.12) .or. (DecayMode1.lt..0) .or. (DecayMode2.lt.0) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop 1
    endif

    !lepton filter

    if( RequestOS.lt.RequestOSSF ) then
        RequestOS = RequestOSSF
    endif
    if( RequestNLeptons .lt. 2*RequestOS ) then
        RequestNLeptons = 2*RequestOS
    endif

    !WidthScheme and reweighting
    if( CalcPMZZ .and. (WidthScheme.le.0 .or. WidthSchemeIn.le.0) ) then
        ReweightDecay=.true.
    endif
    if( .not.ReadLHEFile .and. .not.ConvertLHEFile .and. .not.CalcPMZZ ) then
        if( ReweightDecay .or. WidthSchemeIn.gt.0 ) then
            call Error("ReweightDecay and WidthSchemeIn only make sense in ReadLHE mode")
        endif
        if( WidthScheme.le.0 ) then
            WidthScheme = 2
        endif
    else
        if( WidthScheme.le.0 .and. WidthSchemeIn.le.0 ) then
            if( ReweightDecay.and..not.CalcPMZZ ) then
                print *, "If you want to reweight the decay, you need to specify a width scheme to correct"
                print *, " for the VV branching fraction/matrix element."
                stop 1
            endif
            WidthScheme = 2
            WidthSchemeIn = 2
        elseif( WidthScheme.le.0 .and. WidthSchemeIn.gt.0 ) then
            WidthScheme = WidthSchemeIn
        elseif( WidthScheme.gt.0 .and. WidthSchemeIn.le.0 ) then
            WidthSchemeIn = WidthScheme
        else !both > 0
            !nothing
        endif
    endif

    !WriteFailedEvents

    if( WriteFailedEvents.lt.0 .or. WriteFailedEvents.gt.2 ) then
        call Error("WriteFailedEvents can only be 0, 1, or 2.  Please see the manual.")
    endif

    !couplings

    if( SetAnomalousSpin0gg .and. .not.Setghg2 ) then
        call Error("If you set an anomalous spin 0 gg coupling, you need to explicitly set ghg2 as well")
    endif
    if( SetAnomalousSpin0ZZ .and. .not.Setghz1 ) then
        call Error("If you set an anomalous spin 0 ZZ coupling, you need to explicitly set ghz1 as well")
    endif
    if( SetAnomalousSpin1qq .and. .not.(Setspin1qqleft.and.Setspin1qqright) ) then
        call Error("If you set an anomalous spin 1 qq coupling, you need to set both zprime_qq_left and zprime_qq_right")
    endif
    if( SetAnomalousSpin1ZZ .and. .not.Set1minus ) then
        call Error("If you set an anomalous spin 1 ZZ coupling, you need to explicitly set zprime_zz_1 as well")
    endif
    if( SetAnomalousSpin2gg .and. .not.Seta2 ) then
        call Error("If you set an anomalous spin 2 gg coupling, you need to explicitly set a2 as well")
    endif
    if( SetAnomalousSpin2qq .and. .not.(Setspin2qqleft.and.Setspin2qqright) ) then
        call Error("If you set an anomalous spin 2 qq coupling, you need to explicitly set both graviton_qq_left and graviton_qq_right")
    endif
    if( SetAnomalousSpin2ZZ .and. .not.Setb2 ) then
        call Error("If you set an anomalous spin 2 ZZ coupling, you need to explicitly set b2 as well")
    endif
    if( SetAnomalousHff .and. .not.Setkappa ) then
        call Error("If you set an anomalous Hff coupling, you need to explicitly set kappa as well")
    endif
    if( distinguish_HWWcouplings .and. Process.ne.60 .and. Process.ne.66 ) then
        call Error("The separate HWW couplings are only used for VBF.  For H->WW decay or WH production, please set ghz* instead.")
    endif

return
END SUBROUTINE






SUBROUTINE InitVegas()
use ModKinematics
implicit none
include "vegas_common.f"

  xl(1:mxdim) = 0d0
  xu(1:mxdim) = 1d0
  acc = -1d0
  nprn = 1
  readin=.false.
  writeout=.false.
  stopvegas=.false.

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
   call EvalAlphaS() ! Set alphas at default Mu_Ren. Notice ModParameters::ComputeQCDVariables is automatically called!
   return
END SUBROUTINE


SUBROUTINE InitPDFs()

#if useLHAPDF==1

   use ModParameters
   implicit none
   DOUBLE PRECISION alphasPDF

     call InitPDFset(trim(LHAPDFString)) ! Let LHAPDF handle everything
     call InitPDF(LHAPDFMember)

     alphas_mz=alphasPDF(zmass_pdf)
     ! Dummy initialization, just in case. These values are not used.
     !nloops_pdf = 1
     zmass_pdf = M_Z

#else

   use ModParameters
   use ModKinematics
   implicit none
   character :: pdftable*(100)

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
        write(6,*),"main.F90::InitPDFs: PDFSet",PDFSet,"QCD parameters are unknown. Please double-check! Stopping JHUGen..."
        stop
        ! Could also have used these instead of the stop statement, but why introduce arbitrary number?
        !alphas_mz = 0.13229060d0
        !nloops_pdf = 1
     endif

#endif

     call InitPDFValues() ! Call this only once
   return
END SUBROUTINE




SUBROUTINE InitParameters
use ModParameters
implicit none

IF( COLLIDER.EQ.1) THEN 
  Collider_Energy  = LHC_Energy
ELSEIF( COLLIDER.EQ.2 ) THEN
  Collider_Energy  = TEV_Energy
ELSEIF( COLLIDER.EQ.0 ) THEN
  Collider_Energy  = ILC_Energy
ENDIF


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
use ModTTBHiggs
implicit none
include "vegas_common.f"

! NOTE: NDim for weighted vegas run is not minimal

      NDim = 0
      NDim = NDim + 6    ! PS integration
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
      !- HVBF with decays
      if(Process.eq.66) then
         NDim = 5
         NDim = NDim + 2 ! sHat integration
         NDim = NDim + 8
         NDim = NDim + 1
         NDim = NDim + 1
         if( unweighted ) NDim = NDim + 1  ! random number which decides if event is accepted
         
         VegasIt1_default = 5
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
         if( unweighted ) NDim = NDim + 1  ! random number which decides if event is accepted
         
         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif
     ! RR added -- tb+H
      if(Process.eq.111) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif

     ! RR added -- t+H s-schannel
      if(Process.eq.112) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif
     ! RR added -- tb+H s-channel
      if(Process.eq.113) then
         NDim = 9
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default =  500000
         VegasNc1_default =  500000
         VegasNc2_default =  500000
      endif

END SUBROUTINE






SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModCrossSection
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
    AlertCounter=0

    avgcs = 0d0

    itmx = 1
    ncall= VegasNc2
    if (process.eq.60 .or. process.eq.61) then 
      call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    elseif (process.eq.62 .or. process.eq.61) then 
      call vegas1(EvalWeighted_HJ,VG_Result,VG_Error,VG_Chi2)
    elseif (Process.eq.50) then
      call vegas1(EvalWeighted_VHiggs,VG_Result,VG_Error,VG_Chi2)
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

        open(unit=io_CSmaxFile,file='CSmax.bin',form='unformatted',status='replace')
        WRITE(io_CSmaxFile) CSMAX
        close(io_CSmaxFile)
    else
        open(unit=io_CSmaxFile,file='CSmax.bin',form='unformatted')
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
               
   
  endif! unweighted
  
  


return
END SUBROUTINE








SUBROUTINE StartVegas_NEW(VG_Result,VG_Error)
use ModCrossSection
use ModCrossSection_HJJ
use ModCrossSection_TTBH
use ModCrossSection_BBBH
use ModCrossSection_TH
use ModKinematics
use ModParameters
use modHiggsJJ
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),calls1,calls2,calls_rescale
real(8) :: dum, RES(-5:5,-5:5),ResFrac(-5:5,-5:5),TotalXSec
integer :: i, i1, j1,PChannel_aux, PChannel_aux1,NHisto,ijSel(1:121,1:3)
include 'csmaxvalue.f'
integer :: flav1,flav2,StatusPercent,MissingEvents,MaxEvts,imax
integer :: VegasSeed
character :: ProcessStr*(3)
logical :: UseBetaVersion=.false.

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


    call cpu_time(time_start)    
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
    if( Process.eq.80 ) call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.90 ) call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
    
    if( Process.eq.60 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.66 ) call vegas(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.61 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)

    if( Process.eq.110) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.111) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)

    if( Process.eq.112) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.113) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)

    !DATA RUN
    call ClearHisto()   
    warmup = .false.
    EvalCounter=0
    RejeCounter=0
    AccepCounter=0
    AlertCounter=0
    avgcs = 0d0
    itmx = 2
    ncall= VegasNc2
    if( Process.eq.80 ) call vegas1(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.90 ) call vegas1(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)

    if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.66 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)

    if( Process.eq.110) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.111) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.112) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.113) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)




elseif(unweighted.eqv..true.) then  !----------------------- unweighted events

if( Process.eq.60 ) UseBetaVersion=.true.
if( Process.eq.61 ) UseBetaVersion=.true.


if( UseBetaVersion ) then 
! !-------------------new stuff -------------------
    write(io_stdout,"(1X,A)")  "Scanning the integrand"
    warmup = .true.
    itmx = 5
    ncall= VegasNc0
    outgridfile="vegas_"//trim(ProcessStr)//".grid"  
    ingridfile=trim(outgridfile)
    
    if( ReadCSmax ) then
        readin=.true.
        writeout=.false.
        itmx = 3 
    else
        readin=.false.
        writeout=.true.
        if( Process.eq.60 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.61 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
        itmx = 3
    endif
      
    
    CrossSecMax(:,:) = 0d0
    CrossSec(:,:) = 0d0
        
!     if( Process.eq.80 ) call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2) ! adjust to LHE format
!     if( Process.eq.90 ) call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
!     if( Process.eq.66 ) call vegas(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
!     if( Process.eq.110) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
!     if( Process.eq.111) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)    
!     if( Process.eq.112) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)    
!     if( Process.eq.113) call vegas(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)    


    call vegas_get_calls(calls1)
    CrossSec(:,:) = CrossSec(:,:)/dble(itmx)    
    write(io_stdout,"(A)")  ""
    write(io_stdout,"(2X,A,F10.3,A,F10.3,A,F10.3)") "Total xsec: ",VG_Result, " +/-",VG_Error, " fb    vs.",sum(CrossSec(:,:))
    call InitOutput(VG_Result, VG_Error)

    RequEvents(:,:)=0
    do i1=-5,5
    do j1=-5,5
        RequEvents(i1,j1) = RequEvents(i1,j1) + nint( CrossSec(i1,j1)/VG_Result * VegasNc2 )
    enddo
    enddo



    if( Process.eq.60 ) then
       call get_VBFchannelHash(ijSel)
    else
!        call get_GENchannelHash(ijSel)
       call get_HJJchannelHash(ijSel)
    endif
! do i=1,121
!          i1 = ijSel(i,1)
!          j1 = ijSel(i,2)
!          if( j1.gt.i1 ) cycle
!          write(io_stdout,"(1X,I3,A,I3,I3,A,3X,F8.3,I9)") i," Fractional partonic xsec ",i1,j1," "//getLHEParticle(i1)//" "//getLHEParticle(j1)//" ",CrossSec(i1,j1)/VG_Result,RequEvents(i1,j1) 
! enddo
! pause
    
    do i=1,121
         i1 = ijSel(i,1)
         j1 = ijSel(i,2)
         if( RequEvents(i1,j1).gt.0 .and. ijSel(i,3).eq.1 ) write(io_stdout,"(1X,I3,A,I3,I3,A,3X,F8.3,I9)") i," Fractional partonic xsec ",i1,j1," "//getLHEParticle(i1)//" "//getLHEParticle(j1)//" ",CrossSec(i1,j1)/VG_Result,RequEvents(i1,j1) 
    enddo
    write(io_stdout,"(2X,A,F8.3,I9)") "Sum        partonic xsec   x   x    ",sum(CrossSec(:,:))/VG_Result,sum(RequEvents(:,:))
  
  
  
  
!   add some events that got lost due to rounding errors
!   distribute them according to the partonic cross section fractions and finally add the last pieces to the largest partonic contribution
    MissingEvents = VegasNc2 - sum(RequEvents(:,:))
    if( MissingEvents.ne.0 ) then
!         print *, "MISSING EVENTS",MissingEvents
        MaxEvts = -10000
        do i=1,121
            i1 = ijSel(i,1)
            j1 = ijSel(i,2)
            RequEvents(i1,j1) = RequEvents(i1,j1) + nint( CrossSec(i1,j1)/VG_Result * MissingEvents )
            if( RequEvents(i1,j1).gt.MaxEvts ) then
              MaxEvts = RequEvents(i1,j1)
              imax=i
            endif
!             print *, "adding",i1,j1,nint( CrossSec(i1,j1)/VG_Result * MissingEvents )
        enddo       
        MissingEvents = VegasNc2 - sum(RequEvents(:,:))
!         print *, "MISSING EVENTS",MissingEvents
        i1 = ijSel(imax,1)
        j1 = ijSel(imax,2)
        RequEvents(i1,j1) = RequEvents(i1,j1) + MissingEvents
        write(*,"(2X,A,I9)") "Adjusting number of events. New event count=",sum(RequEvents(:,:))
    endif
    
      
      
    
    write(io_stdout,"(A)")  ""
    write(io_stdout,"(1X,A)")  "Event generation"
    call ClearHisto()   
    warmup = .false.
    itmx = 1
!     nprn = 0  
    EvalCounter = 0
    RejeCounter = 0
    AlertCounter = 0
    AccepCounter_part(:,:) = 0 
    StatusPercent = 0d0
    
    CrossSecMax(:,:) = 1.0d0 * CrossSecMax(:,:)    !  adjustment factor
    call cpu_time(time_start)    
    

    itmx=200000
    ncall= 1000000
    call vegas_get_calls(calls2)
    calls_rescale = calls1/calls2
    CrossSecMax(:,:) = CrossSecMax(:,:) * calls_rescale    

    
!     do while( StatusPercent.lt.100d0  )
!         if( Process.eq.80 ) call vegas1(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)! adjust to LHE format
    !     if( Process.eq.90 ) call vegas1(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    !     if( Process.eq.66 ) call vegas1(EvalWeighted_HJJ_fulldecay,VG_Result,VG_Error,VG_Chi2)
        if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
!         if( Process.eq.110) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)
!         if( Process.eq.111) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)      
!         if( Process.eq.112) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)      
!         if( Process.eq.113) call vegas1(EvalWeighted_TH,VG_Result,VG_Error,VG_Chi2)      
        call system('clear')
        write(io_stdout,*) ""
        do i1=-5,5
        do j1=-5,5
            if( RequEvents(i1,j1).gt.0 ) then 
!                write(io_stdout,"(1X,A,I4,I4,2X,I7,2X,I7,2X,F8.3,1X,A)") "Generated events ", i1,j1,(AccepCounter_part(i1,j1)),(RequEvents(i1,j1)),dble(AccepCounter_part(i1,j1))/dble(RequEvents(i1,j1))*100d0,"%"
               call PrintStatusBar2(int(dble(AccepCounter_part(i1,j1))/(dble(RequEvents(i1,j1)))*100),"channel "//getLHEParticle(i1)//" "//getLHEParticle(j1)//" ")
            endif   
        enddo
        enddo  
        StatusPercent = int(100d0*dble(sum(AccepCounter_part(:,:)))/dble(sum(RequEvents(:,:)))  )   
!     enddo
    call cpu_time(time_end)  
    
    print *, " Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter+1d-10) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90."
        write(io_stdout, *) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_stdout, *) "       Increase CSMAX in main.F90."
    endif
    write(io_stdout,*)  " event generation rate (events/sec)",dble(sum(AccepCounter_part(:,:)))/(time_end-time_start+1d-10)

    


else! beta version




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
                elseif( Process.eq.66 ) then
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
        open(unit=io_CSmaxFile,file='CSmax.bin',form='unformatted',status='replace')
        WRITE(io_CSmaxFile) CSMAX,VG
        close(io_CSmaxFile)
    else
        open(unit=io_CSmaxFile,file='CSmax.bin',form='unformatted')
        READ(io_CSmaxFile) CSMAX,VG
        close(io_CSmaxFile)
    endif
   
   CSmax(:,:)   = 1.5d0 * CSmax(:,:)    !  adjustment factor

   VG(:,:) = VG(:,:)/dble(VegasNc0)
   TotalXSec = sum(  VG(:,:) )
   print *, ""    
   write(io_stdout,"(1X,A,F10.3)") "Total xsec: ",TotalXSec


    RequEvents(:,:)=0
    do i1=-5,5
    do j1=-5,5
         RequEvents(i1,j1) = RequEvents(i1,j1) + int( VG(i1,j1)/TotalXSec * VegasNc2 )
    enddo
    enddo
    
   
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
              elseif( Process.eq.66 ) then
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

!    search for line with first event
     FirstEvent = .false.
     FoundHiggsMass=.false.
     FoundHiggsWidth=.false.
     InMadgraphMassBlock=.false.
     do while ( .not.FirstEvent )
        read(16,fmt="(A160)",IOSTAT=stat,END=99) FirstLines

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
        if (ReweightDecay .and. FirstLines(1:7).eq."bwshape") then
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
            elseif( WidthScheme.gt.0 .and. WidthSchemeIn.lt.0 .and. .not.ReweightDecay .and. j.ne.WidthScheme ) then
                print *, "WidthScheme is specified to be ", WidthSchemeIn, " but the LHE file says:"
                print *, FirstLines
                print *, "If you want to reweight the propagator from ", WidthSchemeIn, " to ", WidthScheme, " please specify ReweightDecay=1"
                stop 1
            endif
            if( WidthScheme.le.0 ) read(FirstLines(8:i),*) WidthScheme
            if( WidthSchemeIn.le.0 ) read(FirstLines(8:i),*) WidthSchemeIn
        endif

        if( Index(FirstLines, "<event").ne.0 ) FirstEvent=.true.
     enddo
99   continue

     call ReopenInFile()

return
END SUBROUTINE

SUBROUTINE InitReadLHE(BeginEventLine)
use ModParameters
use ModMisc
implicit none
logical :: FirstEvent, WroteHeader
character(len=160) :: FirstLines
integer :: stat
character(len=100), intent(out) :: BeginEventLine

     write(io_LHEOutFile ,'(A)') '<LesHouchesEvents version="1.0">'
     FirstEvent = .false.
     WroteHeader = .false.
     do while ( .not.FirstEvent )
        read(16,fmt="(A160)",IOSTAT=stat,END=99) FirstLines
        if ( FirstLines(1:4).eq."<!--" .and. .not.WroteHeader ) then
            call InitOutput(1d0, 1d14)
            WroteHeader = .true.
        endif
        if (index(FirstLines,"<MG").ne.0 .and. .not.WroteHeader) then  !Sometimes MadGraph doesn't have a comment at the beginning
            call InitOutput(1d0, 1d14)                                 !In that case put the JHUGen header before the MadGraph
            write(io_LHEOutFile, "(A)") "-->"                          ! proc card, etc.
            WroteHeader = .true.                                       !and put the Higgs mass/width in a separate comment
        endif
        if (Index(FirstLines,"<init>").ne.0 .and. .not.WroteHeader ) then !If not now, when?
            call InitOutput(1d0, 1d14)
            WroteHeader = .true.
        endif

        if( Index(FirstLines, "<event").ne.0 ) then
            FirstEvent=.true.
            BeginEventLine = trim(FirstLines)
        else
            if( importExternal_LHEinit ) then
                if( Index(FirstLines,"<LesHouchesEvents").ne.0 .or. Index(FirstLines,"<!--").ne.0 ) then
                else
                  write(io_LHEOutFile,"(A)") trim(firstlines)
                endif
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
use ModParameters
use ModMisc
implicit none
include 'csmaxvalue.f'
integer,parameter :: maxpart=30!=max.part particles in LHE file; this parameter should match the one in WriteOutEvent of mod_Kinematics
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),Res,EMcheck(1:4),DecayWeight,DecayWidth,DecayWidth0
real(8) :: HiggsDK_Mom(1:4,1:13),Ehat,GetMZZProbability
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
integer,parameter :: PMZZcalls = 200000
integer,parameter :: PMZZ0calls = 1000000
logical :: Empty


if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2
AccepLastPrinted = 0

call InitReadLHE(BeginEventLine)

     if( ReweightDecay ) print *, " finding P_decay(m_Reso) with ",PMZZ0calls," points" !otherwise it happens instantaneously, so no need to print
     DecayWidth0 = GetMZZProbability(m_Reso,PMZZ0calls)

     print *, " finding maximal weight for mZZ=mReso with ",VegasNc0," points"
     VG = zero
     CSmax = zero
     EHat = M_Reso! fixing Ehat to M_Reso which should determine the max. of the integrand
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
     csmax(0,0)   = 1.5d0*csmax(0,0)    !  savety buffer

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
         read(16,"(A)") EventLine(0)
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
            read(16,fmt="(A160)") EventLine(nline)
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
            if( IsALHELepton(LHE_IDUP(nline)) ) then
                  LeptInEvent(0) = LeptInEvent(0) + 1
                  LeptInEvent( LeptInEvent(0) ) = LHE_IDUP(nline)
            endif
         enddo
            
         
!         accept/reject sampling for H->VV decay contribution
          EHat = pH2sq
          DecayWeight = 0d0
          
          DecayWidth = GetMZZProbability(EHat,PMZZcalls)!  could also be used to determine csmax for this particular event to improve efficiency (-->future work)
          WeightScaleAqedAqcd(1) = WeightScaleAqedAqcd(1) * DecayWidth/DecayWidth0

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
              read(16,fmt="(A160)",IOSTAT=stat,END=99) OtherLines(1:160)
              if(OtherLines(1:30).eq."</LesHouchesEvents>") then
                  if( RequestNLeptons.gt.0 ) then
                    write(io_LHEOutFile,"(A)") "<!-- Lepton filter information:"
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
        write(io_stdout,*) "       Increase CSMAX in main.F90 or VegasNc1."
    endif
   write(io_stdout,*)  "Event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)
   if( RequestNLeptons.gt.0 ) write(io_stdout,"(A,1F6.2,A)") " Lepton filter efficiency:",dble(AccepCounter)/dble(NEvent)*100d0," %"

   
    write(io_LogFile,*) ""
    write(io_LogFile,*) "Evaluation Counter: ",EvalCounter
    write(io_LogFile,*) "Acceptance Counter: ",AccepCounter
    write(io_LogFile,*) "Rejection  Counter: ",RejeCounter
    write(io_LogFile,*) "Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90 or VegasNc1."
    endif
   write(io_LogFile,*)  "Event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)
   if( RequestNLeptons.gt.0 ) write(io_LogFile,"(A,1F6.2,A)") " Lepton filter efficiency:",dble(AccepCounter)/dble(NEvent)*100d0," %"



return
END SUBROUTINE




SUBROUTINE StartConvertLHE(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
use ModMisc
implicit none
include 'csmaxvalue.f'
integer,parameter :: maxpart=15!=max.partons; this parameter should match the one in WriteOutEvent of mod_Kinematics
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
         read(16,"(A)") EventLine(0)
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
            read(16,"(A)") EventLine(nline)
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
              read(16,fmt="(A120)",IOSTAT=stat,END=99) PDFLine(1:120)
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





SUBROUTINE GetMZZdistribution()
use ModCrossSection
use ModParameters
implicit none
real(8) :: DecayWeight,DecayWidth,DecayWidth0
real(8) :: Ehat,GetMZZProbability
integer :: nscan
integer,parameter :: nmax=20, Ncalls=200000
real(8),parameter :: ScanRange=120d0*GeV



  DecayWidth0 = GetMZZProbability(M_Reso,Ncalls)
  
  do nscan=-nmax,+nmax,1
     
     EHat = M_Reso+ScanRange*nscan/dble(nmax)
     DecayWidth = GetMZZProbability(EHat,Ncalls)
     write(*,"(1F10.5,1PE16.9)") EHat*100d0,DecayWidth/DecayWidth0
     
  enddo

RETURN
END SUBROUTINE



FUNCTION CalcMZZProbability(EHat,Ncalls)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
implicit none
real(8) :: DecayWeight,yRnd(1:22),Res,CalcMZZProbability
real(8) :: HiggsDK_Mom(1:4,1:13),Ehat
integer :: HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13)
integer :: evals,Ncalls

     CalcMZZProbability = 0d0
     do evals=1,Ncalls
         call random_number(yRnd)
         if( TauDecays.lt.0 ) then
             DecayWeight = EvalUnWeighted_DecayToVV(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP(1:9),HiggsDK_ICOLUP)
         else
             DecayWeight = EvalUnWeighted_DecayToTauTau(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,4:13),HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13))
         endif
         CalcMZZProbability = CalcMZZProbability + DecayWeight
     enddo
     CalcMZZProbability = CalcMZZProbability/dble(evals)

RETURN
END FUNCTION



FUNCTION GetMZZProbability(EHat,Ncalls)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
use ModYRdata
implicit none
!Turning on FastVersion is not advised, even though it it faster
!For the MZZProbability, it uses EHat*Gamma*BR, with Gamma and BR obtained from the yellow report
!1) with anomalous couplings, these values are no longer correct
!2) at high mass, which is where this function is most useful, the YR branching ratio integrates over
!     the wide Higgs mass shape, making the value less correct for this purpose
logical, parameter :: FastVersion = .false.
real(8) :: EHat, BigGamma, CalcMZZProbability, BranchingRatio, GetMZZProbability
integer :: Ncalls

     GetMZZProbability = 1d0

     if( ReweightDecay ) then
         call HTO_gridHt(EHat/GeV,BigGamma)
         call YR_GetBranchingFraction(EHat/GeV, BranchingRatio)

         if( FastVersion ) then
             GetMZZProbability = GetMZZProbability * EHat * BigGamma*BranchingRatio
         else
             GetMZZProbability = GetMZZProbability * CalcMZZProbability(EHat, Ncalls)
         endif

         if( WidthSchemeIn.eq.3 ) then
             !need to divide by overcompensated factor of m4l*BigGamma from POWHEG
             !if( FastVersion ) this is a bit redundant, but BigGamma is basically a lookup table
             !  so it doesn't take much time
             GetMZZProbability = GetMZZProbability / (EHat*BigGamma)
         elseif( WidthSchemeIn.eq.1 ) then
             !divide by overcompensated m4l*gamma_running in the numerator
             GetMZZProbability = GetMZZProbability / (EHat**2*Ga_Reso/M_Reso)
         elseif( WidthSchemeIn.eq.2 ) then
             !numerator is a constant M_Reso*Ga_Reso, do nothing
         else
             call Error("Invalid WidthSchemeIn!")
         endif
     endif

     if( WidthScheme.ne.WidthSchemeIn ) then
         !reweight the propagator
         GetMZZProbability = GetMZZProbability * GetBWPropagator(EHat**2, WidthScheme) / GetBWPropagator(EHat**2, WidthSchemeIn)
     endif


RETURN
END FUNCTION



SUBROUTINE OpenFiles()
use ModParameters
implicit none
integer :: i

   if( .not.FilesOpened ) then
       FilesOpened = .true.
       call system('mkdir -p ./data')! -p is suppressing error messages if directory already exists

       i = len(trim(DataFile))
       if( DataFile(i-3:i).eq.".lhe" ) then
           !print *, DataFile
           DataFile = DataFile(1:i-4)
           !print *, DataFile
       endif

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

  if( Process.eq.60 .or. Process.eq.66 ) then
     call InitHisto_HVBF()
  elseif( Process.eq.61) then
     call InitHisto_HJJ()
  elseif( Process.eq.62) then
     call InitHisto_HJ()
  elseif (Process.eq.50) then
     call InitHisto_VHiggs()
  elseif (Process.eq.80) then
     call InitHisto_TTBH()
  elseif (Process.eq.90) then
     call InitHisto_BBBH()
  elseif (Process.eq.110 .or. Process .eq. 111 .or. Process.eq.112 .or. Process .eq. 113) then
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
          NumHistograms = 8
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
          
          Histo(6)%Info   = "y(j1)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal = -5d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "y(j2)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal = -5d0
          Histo(7)%SetScale= 1d0
          
          Histo(8)%Info   = "m_4l"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal = 120d0*GeV
          Histo(8)%SetScale= 1d0/GeV
! 
!           Histo(7)%Info   = "weights"
!           Histo(7)%NBins  = 50
!           Histo(7)%BinSize= 0.25d0
!           Histo(7)%LowVal = -18d0
!           Histo(7)%SetScale= 1d0



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
          NumHistograms = 9
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
          Histo(4)%NBins  = 80
          Histo(4)%BinSize= 300d0*GeV/80d0
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
    if( Collider.eq.0 ) then
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
        open(unit=io_LHEInFile,file=trim(LHAPDF_DATA_PATH)//"/"//LHAPDFString,form='formatted',access= 'sequential',status='old')
        pdfsup1 = -999
        do while( pdfsup1.lt.0 )
            read(16,fmt="(A160)",IOSTAT=stat,END=99) ReadLines
            if( ReadLines(1:8).eq."SetIndex" ) then
                read(ReadLines(10:500),*) pdfsup1
                pdfsup1 = pdfsup1 + LHAPDFMember
                exit
            endif
        enddo
99      CONTINUE
        if( pdfsup1.lt.0 ) then
            print *, "Error: couldn't get the PDF set index.  Writing 0."
            pdfsup2 = 0
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
    if( .not.unweighted ) then
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
        write(io_LHEOutFile ,'(A,A6,A)') 'Output from the JHUGenerator ',trim(JHUGen_Version),' described in arXiv:1001.3396 [hep-ph],arXiv:1208.4018 [hep-ph],arXiv:1309.4819 [hep-ph]'

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

        if( (ReadLHEFile .or. ConvertLHEFile) .and. (importExternal_LHEinit) ) then
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
    if( Process.eq.66) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.50) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.80) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.90) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.110) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.111) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.112) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.113) write(TheUnit,"(4X,A,F7.2,A,F10.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( ReadLHEFile )    write(TheUnit,"(4X,A)") "           (This is ReadLHEFile mode. Resonance mass/width are read from LHE input parameters.)"
    if( ConvertLHEFile ) write(TheUnit,"(4X,A)") "           (This is ConvertLHEFile mode. Resonance mass/width are read from LHE input parameters.)"
    if( &
         (.not.ReadLHEFile .and. (Process.le.2 .or. Process.eq.50 .or. Process.eq.60 .or. Process.eq.66 .or. ((TopDecays.eq.1).and.Process.eq.80) .or. (Process.ge.110 .and. Process.le.113))) &
    .or. (ReadLHEFile .and. TauDecays.ne.0) &
    .or. ConvertLHEFile ) &
    then
        if( .not.ReadLHEFile .and. (ConvertLHEFile .or. Process.eq.50 .or. (Process.ge.110 .and. Process.le.113)) ) then
            write(TheUnit,"(4X,A,I2,2X,A,I2)") "DecayMode1:",DecayMode1
        else if( ReadLHEFile .or. Process.le.2 .or. Process .eq. 80 ) then
            write(TheUnit,"(4X,A,I2,2X,A,I2)") "DecayMode1:",DecayMode1, "DecayMode2:",DecayMode2
        endif
        if( Process.eq.60 .or. Process.eq.66 .or. IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z-boson: mass=",M_Z*100d0,", width=",Ga_Z*100d0
        if( Process.eq.60 .or. Process.eq.66 .or. IsAWDecay(DecayMode1) .or. IsAWDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W-boson: mass=",M_W*100d0,", width=",Ga_W*100d0
    endif
    if( Process.eq.80 .or. Process.eq.110 .or. Process.eq.111 .or.Process.eq.112 .or. Process.eq.113 ) write(TheUnit,"(4X,A,F8.4,A,F6.4)") "Top quark mass=",m_top*100d0,", width=",Ga_top*100d0
    if( Process.eq.80 .or. Process.eq.110 .or. Process.eq.111 .or. Process.eq.112 .or. Process.eq.113) write(TheUnit,"(4X,A,I2)") "Top quark decay=",TOPDECAYS
    if( Process.eq.90 ) write(TheUnit,"(4X,A,F8.4,A,F6.4)") "Bottom quark mass=",m_top*100d0
    if( Process.eq.60 .or. Process.eq.61 .or. Process.eq.62 .or. Process.eq.90 .or. &
       ((Process.eq.80 .or. (Process.ge.110 .and. Process.le.113)) .and. m_Top.lt.10d0*GeV) ) then
        write(TheUnit,"(4X,A)") "Jet cuts:"
        write(TheUnit,"(12X,A,F8.2,A)") "pT >= ", pTjetcut/GeV, " GeV"
        if( Process.eq.60 .or. Process.eq.61 .or. Process.eq.80 .or. Process.eq.90) then
            write(TheUnit,"(8X,A,F8.2)") "DeltaR >= ", Rjet
            write(TheUnit,"(11X,A,F8.2,A)") "mJJ >= ", mJJcut/GeV, " GeV"
        endif
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
    write(TheUnit,"(4X,A,20I11)") "Random seed: ",UserSeed

    if( .not. (ReadLHEFile .or. ConvertLHEFile) ) then
        write(TheUnit,"(4X,A)") ""
        if( Process.le.2 ) write(TheUnit,"(4X,A,L,L,L)") "OffXVV: ",OffShellReson,OffShellV1,OffShellV2
        write(TheUnit,"(4X,A,I3)") "WidthScheme: ",WidthScheme
        if( Process.le.2 .or. Process.eq.80 .or. Process.eq.90 ) write(TheUnit,"(4X,A,I1)") "PChannel: ",PChannel
#if useLHAPDF==1
        write(TheUnit,"(4X,A,A,A,I3)") "LHAPDF set ",trim(LHAPDFString), ", member=",LHAPDFMember
#else
        write(TheUnit,"(4X,A,I1)") "PDFSet: ",PDFSet
#endif
        write(TheUnit,"(4X,A,L)") "Unweighted: ",Unweighted
    endif
    if( ReweightDecay ) then
        write(TheUnit,"(4X,A,I3,A,I3)") "Reweighting from WidthScheme ", WidthSchemeIn, " to WidthScheme ", WidthScheme
    endif
    if( Process.le.2 .or. ReadLHEFile ) write(TheUnit,"(4X,A,L)") "Interference: ",includeInterference
    if( (Process.le.2 .or. ReadLHEFile) .and. (IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2))  ) write(TheUnit,"(4X,A,L)") "Intermediate off-shell photons: ",includeGammaStar

    write(TheUnit,"(4X,A)") ""
    if( Process.eq.0 .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.62 .or. Process.eq.66 .or. Process.eq.50 ) then
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
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghg2=",ghg2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghg3=",ghg3,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghg4=",ghg4,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghz1=",ghz1,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghz2=",ghz2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghz3=",ghz3,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "ghz4=",ghz4,"i"
            if( includeGammaStar .or. IsAPhoton(DecayMode2) ) then
                if( includeGammaStar .or. IsAZDecay(DecayMode1) ) then
                    write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs2=",ghzgs2,"i"
                    write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs3=",ghzgs3,"i"
                    write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs4=",ghzgs4,"i"
                endif
                if( includeGammaStar .or. IsAPhoton(DecayMode1) ) then
                    write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs2=",ghgsgs2,"i"
                    write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs3=",ghgsgs3,"i"
                    write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs4=",ghgsgs4,"i"
                endif
                if( includeGammaStar) then
                    write(TheUnit,"(6X,A,F8.2,A)") "m(gammastar) >= ", MPhotonCutoff/GeV, " GeV"
                endif
            endif
            if( cdabs(ghz1_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime= ",ghz1_prime ,"i,","Lambda_z1=",Lambda_z1*100d0
            if( cdabs(ghz1_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime2=",ghz1_prime2,"i,","Lambda_z1=",Lambda_z1*100d0
            if( cdabs(ghz1_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime3=",ghz1_prime3,"i,","Lambda_z1=",Lambda_z1*100d0
            if( cdabs(ghz1_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz1_prime4=",ghz1_prime4,"i,","Lambda_z1=",Lambda_z1*100d0
            if( cdabs(ghz2_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime= ",ghz2_prime ,"i,","Lambda_z2=",Lambda_z2*100d0
            if( cdabs(ghz2_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime2=",ghz2_prime2,"i,","Lambda_z2=",Lambda_z2*100d0
            if( cdabs(ghz2_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime3=",ghz2_prime3,"i,","Lambda_z2=",Lambda_z2*100d0
            if( cdabs(ghz2_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz2_prime4=",ghz2_prime4,"i,","Lambda_z2=",Lambda_z2*100d0
            if( cdabs(ghz3_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime= ",ghz3_prime ,"i,","Lambda_z3=",Lambda_z3*100d0
            if( cdabs(ghz3_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime2=",ghz3_prime2,"i,","Lambda_z3=",Lambda_z3*100d0
            if( cdabs(ghz3_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime3=",ghz3_prime3,"i,","Lambda_z3=",Lambda_z3*100d0
            if( cdabs(ghz3_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz3_prime4=",ghz3_prime4,"i,","Lambda_z3=",Lambda_z3*100d0
            if( cdabs(ghz4_prime ).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime= ",ghz4_prime ,"i,","Lambda_z4=",Lambda_z4*100d0
            if( cdabs(ghz4_prime2).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime2=",ghz4_prime2,"i,","Lambda_z4=",Lambda_z4*100d0
            if( cdabs(ghz4_prime3).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime3=",ghz4_prime3,"i,","Lambda_z4=",Lambda_z4*100d0
            if( cdabs(ghz4_prime4).ne.0d0 ) write(TheUnit,"(6X,A,2E16.8,A2,4X,A,1PE12.4)") "ghz4_prime4=",ghz4_prime4,"i,","Lambda_z4=",Lambda_z4*100d0
        endif
    elseif( Process.eq.1 ) then
        write(TheUnit,"(4X,A)") "spin-1-VV couplings: "
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_qq_left =",zprime_qq_left,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_qq_right=",zprime_qq_right,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_zz_1=",zprime_zz_1,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "zprime_zz_2=",zprime_zz_2,"i"
    elseif( Process.eq.2 ) then
        write(TheUnit,"(4X,A)") "spin-2-VV couplings: "
        write(TheUnit,"(6X,A,L)") "generate_bis=",generate_bis
        write(TheUnit,"(6X,A,L)") "use_dynamic_MG=",use_dynamic_MG
        write(TheUnit,"(6X,A,2E16.8,A1)") "a1=",a1,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a2=",a2,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a3=",a3,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a4=",a4,"i"
        write(TheUnit,"(6X,A,2E16.8,A1)") "a5=",a5,"i"
        if( generate_bis ) then 
            write(TheUnit,"(6X,A,2E16.8,A1)") "b1 =",b1,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b2= ",b2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b3 =",b3,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b4 =",b4,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b5 =",b5,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b6 =",b6,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b7 =",b7,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b8 =",b8,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b9 =",b9,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "b10=",b10,"i"
        else
            write(TheUnit,"(6X,A,2E16.8,A1)") "c1 =",c1,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c2 =",c2,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c3 =",c3,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c41 =",c41,"i"
            write(TheUnit,"(6X,A,2E16.8,A1)") "c42 =",c42,"i"
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

        write(io_stdout,*) ""
        write(io_stdout,"(2X,A)") "Command line arguments:"
        write(io_stdout,"(4X,A)") "Collider:   1=LHC, 2=Tevatron, 0=e+e-"
        write(io_stdout,"(4X,A)") "Process: 0=spin-0, 1=spin-1, 2=spin-2 resonance"
        write(io_stdout,"(4X,A)") "         50=pp/ee->VH, 60=weakVBF, 61=pp->Hjj, 62=pp->Hj"
!         write(io_stdout,"(4X,A)") "         50=pp/ee->VH, 60/66=weakVBF (without/with decay+SM bkg), 61=pp->Hjj, 62=pp->Hj"
        write(io_stdout,"(4X,A)") "         80=pp->ttbar+H, 90=pp->bbbar+H"
        write(io_stdout,"(4X,A)") "         110=pp->t+H (t-channel), 111=pp->tbar+H (t-ch.), 112=pp->t+H (s-ch.), 113=pp->tbar+H (s-ch.)"
        write(io_stdout,"(4X,A)") "MReso:      resonance mass (default=125.60), format: yyy.xx"
        write(io_stdout,"(4X,A)") "DecayMode1: decay mode for vector boson 1 (Z/W+/gamma)"
        write(io_stdout,"(4X,A)") "DecayMode2: decay mode for vector boson 2 (Z/W-/gamma)"
        write(io_stdout,"(4X,A)") "              0=Z->2l,  1=Z->2q, 2=Z->2tau, 3=Z->3nu,"
        write(io_stdout,"(4X,A)") "              4=W->lnu, 5=W->2q, 6=W->taunu,"
        write(io_stdout,"(4X,A)") "              7=gamma, 8=Z->2l+2tau,"
        write(io_stdout,"(4X,A)") "              9=Z->anything, 10=W->lnu+taunu, 11=W->anything"
        write(io_stdout,"(4X,A)") "              20=tau->nu + W(select with TauDK)"
        write(io_stdout,"(4X,A)") "              30=top->b bbar"
        write(io_stdout,"(4X,A)") "              40=top->b + W(select with TopDK)"
        write(io_stdout,"(4X,A)") "TopDK:      decay mode for tops in ttbar+H or H->ttbar, 0=stable, 1=decaying (use DecayMode1/2 = 4,5 for W+/W-"
        write(io_stdout,"(4X,A)") "TauDK:      decay mode for taus in H->tautau, 0=stable, 1=decaying (use DecayMode1/2 = 4,5 for W+/W-"
        write(io_stdout,"(4X,A)") "BotDK:      decay mode for bottom quarks in H->bbar, 0=deactivated, 1=activated"
        write(io_stdout,"(4X,A)") "PChannel:   0=g+g, 1=q+qb, 2=both"
        write(io_stdout,"(4X,A)") "OffXVV:     off-shell option for resonance(X),or vector bosons(VV)"
        write(io_stdout,"(4X,A)") "WidthScheme:1=running width, 2=fixed width (default), 3=complex pole scheme"
#if useLHAPDF==1
        write(io_stdout,"(4X,A)") "LHAPDF:     name of the LHA PDF file, e.g. NNPDF30_lo_as_0130/NNPDF30_lo_as_0130.info"
        write(io_stdout,"(4X,A)") "LHAPDFMem:  member PDF number, default=0 (best fit)"
#else
        write(io_stdout,"(4X,A)") "PDFSet:     1=CTEQ6L1(default), 2=MSTW2008LO,  2xx=MSTW with eigenvector set xx=01..40), 3=NNPDF3.0LO"
#endif        
        write(io_stdout,"(4X,A)") "VegasNc0:   number of evaluations for integrand scan"
        write(io_stdout,"(4X,A)") "VegasNc1:   number of evaluations for accept-reject sampling"
        write(io_stdout,"(4X,A)") "VegasNc2:   number of events for accept-reject sampling"
        write(io_stdout,"(4X,A)") "Unweighted: 0=weighted events, 1=unweighted events"
        write(io_stdout,"(4X,A)") "Interf:     0=neglect interference for 4f final states, 1=include interference"
        write(io_stdout,"(4X,A)") "DataFile:   LHE output file"
        write(io_stdout,"(4X,A)") "ReadLHE:    LHE input file from external file (only spin-0)"
        write(io_stdout,"(4X,A)") "ConvertLHE: LHE input file from external file (only spin-0)"
        write(io_stdout,*) ""

END SUBROUTINE




SUBROUTINE PrintLogo(TheUnit)
use modParameters
implicit none
integer :: TheUnit

    write(TheUnit,"(A90)") " "
    write(TheUnit,"(A90)") " ***************************************************************************************"
    write(TheUnit,"(A50,A6,A34)") " *                                JHU Generator ",trim(JHUGen_Version),"             *"
    write(TheUnit,"(A90)") " ***************************************************************************************"
    write(TheUnit,"(A90)") " *                                                                                     *"
    write(TheUnit,"(A90)") " *   Spin and parity determination of single-produced resonances at hadron colliders   *"
    write(TheUnit,"(A90)") " *                                                                                     *"
    write(TheUnit,"(A90)") " *          I. Anderson, S. Bolognesi, F. Caola, Y. Gao, A. Gritsan, Z. Guo,           *"
    write(TheUnit,"(A90)") " *        C. Martin, K. Melnikov, R.Rontsch, H. Roskes, U. Sarica, M. Schulze,         *" 
    write(TheUnit,"(A90)") " *                   N. Tran, A. Whitbeck, M. Xiao, C. You, Y. Zhou                    *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D81 (2010) 075022;  arXiv:1001.3396 [hep-ph],              *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D86 (2012) 095031;  arXiv:1208.4018 [hep-ph],              *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D89 (2014) 035007;  arXiv:1309.4819 [hep-ph].              *"
    write(TheUnit,"(A90)") " *                                                                                     *"
    write(TheUnit,"(A90)") " ***************************************************************************************"
    write(TheUnit,"(A90)") " "
return
END SUBROUTINE



