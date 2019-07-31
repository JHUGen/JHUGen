MODULE ModJHUGen
   use ModParameters
   use ModJHUGenMELA
   use ModKinematics
#if compiler==1
   use ifport
#endif
   implicit none

   public :: InitFirstTime,ResetPDFs

   CONTAINS



SUBROUTINE InitFirstTime(pdftable,pdfstrlength,pdfmember,collider_sqrts)
   use ModParameters
   use ModKinematics
   use ModMisc
   implicit none
   integer :: pdfstrlength
   character(len=pdfstrlength) pdftable
   integer :: pdfmember
   real(8), intent(in) :: collider_sqrts

   call SetJHUGenDefaults()

   includeInterference=.true.
   includeGammaStar=.true.
   MPhotonCutoff=4d0*GeV
   WidthScheme=0
   PDFSet=3      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40
   LHAPDFString = pdftable
   LHAPDFMember = pdfmember
   lenLHAPDFString = pdfstrlength

!---------------------------
   call PrintLogo(io_stdout, "JHUGen MELA")
!---------------------------
   call ResetMubarHGabarH()
!---------------------------
   call ComputeEWVariables()
   call ComputeCKMElements(VCKM_ud, VCKM_us, VCKM_cd, VCKM_cs, VCKM_ts, VCKM_tb)
   print *, "JHUGen CKM initialization"
   print *, "Vud = ",VCKM_ud
   print *, "Vus = ",VCKM_us
   print *, "Vub = ",VCKM_ub
   print *, "Vcd = ",VCKM_cd
   print *, "Vcs = ",VCKM_cs
   print *, "Vcb = ",VCKM_cb
   print *, "Vtd = ",VCKM_td
   print *, "Vts = ",VCKM_ts
   print *, "Vtb = ",VCKM_tb
!---------------------------

#if useLHAPDF==1
   if( LHAPDFString.eq."" ) then
      print *, "Need to specify pdf file for LHAPDF"
      stop
   endif
#else
   if( LHAPDFString.eq."" .and. PDFSet.eq.3) then
      print *, "Need to specify pdf file for NNPDF"
      stop
   endif
#endif

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

   call InitPDFs()
   IF( COLLIDER.EQ.1) THEN
      Collider_Energy  = collider_sqrts
   ELSE
      print *,"Collider not implemented."
      stop
   ENDIF

#if useCollier==1
   Collier_maxNLoopProps = -1
   Collier_maxRank = -1
#endif
   call InitCOLLIER(4,3) ! Arguments for ggZH

   return

END SUBROUTINE

SUBROUTINE SetJHUGenDefaults()
   use ModParameters
   use ModHashCollection
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
   H_DK = .false.
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
   RequestOS=-1
   RequestOSSF=-1
   CountTauAsAny = .true.
   WriteFailedEvents=0

   call SetupHashes()

END SUBROUTINE SetJHUGenDefaults

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

   use ModParameters
   use ModKinematics
   implicit none

#if useLHAPDF==1

   implicit none
   DOUBLE PRECISION alphasPDF

     call InitPDFSetByName(trim(LHAPDFString,lenLHAPDFString)) ! Let LHAPDF handle everything
     call InitPDF(LHAPDFMember)

     alphas_mz=alphasPDF(zmass_pdf)
     ! Dummy initialization, just in case. These values are not used.
     !nloops_pdf = 1
     zmass_pdf = M_Z

#else

     zmass_pdf = M_Z ! Take zmass_pdf=M_Z in pdfs that do not specify this value

     if( PDFSet.eq.1 ) then ! CTEQ6L1
        call SetCtq6(4)  ! 4    CTEQ6L1  Leading Order cteq6l1.tbl

        alphas_mz=0.130d0
        !nloops_pdf=1
     elseif( PDFSet.eq.3 ) then  ! NNPDF 3.0 LO with a_s=0.13
        call NNPDFDriver(LHAPDFString,lenLHAPDFString)
        call NNinitPDF(LHAPDFMember)

        alphas_mz=0.130d0
        !nloops_pdf=1
        zmass_pdf=91.199996948242188d0*GeV
     elseif( (PDFSet.eq.2) .or. (PDFSet.ge.201 .and. PDFSet.le.240) ) then ! MSTW2008 and variations
        alphas_mz=0.13939d0
        !nloops_pdf=1
     else ! Everything else
        write(6,*) "mod_JHUGen.F90::InitPDFs: PDFSet",PDFSet,"QCD parameters are unknown. Please double-check! Stopping JHUGen..."
        stop
        ! Could also have used these instead of the stop statement, but why introduce arbitrary number?
        !alphas_mz = 0.13229060d0
        !nloops_pdf = 1
     endif

#endif

     call InitPDFValues() ! Call this only once
   return
END SUBROUTINE

subroutine ResetPDFs(pdftable,pdfstrlength,pdfmember,pdfid)
implicit none
integer :: pdfstrlength
character(len=pdfstrlength) pdftable
integer :: pdfmember
integer, intent(in) :: pdfid
   PDFSet=pdfid      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40    3: NNPDF3.0
   LHAPDFString = pdftable
   LHAPDFMember = pdfmember
   lenLHAPDFString = pdfstrlength
#if useLHAPDF==1
   if( LHAPDFString.eq."" ) then
      print *, "Need to specify pdf file for LHAPDF"
      stop
   endif
#else
   if( LHAPDFString.eq."" .and. PDFSet.eq.3) then
      print *, "Need to specify pdf file for NNPDF"
      stop
   endif
#endif
   call InitPDFs()
end subroutine


END MODULE ModJHUGen
