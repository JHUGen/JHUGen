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
logical,parameter :: useBetaVersion=.false.! this should be set to .false.


   call GetCommandlineArgs()
   call InitPDFs()!  
   call InitHisto()
   call InitParameters()
   call InitProcess()
   call InitVegas()
   call InitRandomSeed(TheSeeds)
   call OpenFiles()
   call PrintLogo(io_stdout)
   call PrintLogo(io_LogFile)
   call WriteParameters(io_stdout)
   call WriteParameters(io_LogFile)
   call InitOutput()
   write(io_stdout,*) " Running"
   if( .not.useBetaVersion .and.   ConvertLHEFile ) then
        call StartConvertLHE(VG_Result,VG_Error)
   elseif( .not.useBetaVersion .and. .not.ReadLHEFile ) then
        if( Process.eq.80 .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.90 ) then
           call StartVegas_NEW(VG_Result,VG_Error)
        else
           call StartVegas(VG_Result,VG_Error)
        endif
   elseif( useBetaVersion .and. .not.ReadLHEFile ) then
        call StartVegas_BETA(VG_Result,VG_Error)
   elseif( .not.useBetaVersion .and.      ReadLHEFile ) then
!         call StartReadLHE(VG_Result,VG_Error)
        call StartReadLHE_NEW(VG_Result,VG_Error)
   endif
   call WriteHisto(VG_Result,VG_Error,time_end-time_start)
   call FinalizeOutput()
   call CloseFiles()
   write(io_stdout,*) "Done"

END PROGRAM





SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
use ModMisc
implicit none
character :: arg*(500)
integer :: NumArgs,NArg,OffShell_XVV,iargument,CountArg,iinterf

   Collider=1
   VegasIt1=-1
   PDFSet=1      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40
   VegasNc0=-1
   VegasNc1=-1
   VegasNc2=-1
   PChannel=2
   DecayMode1=0  ! Z/W+
   DecayMode2=0  ! Z/W-
   TopDecays=1
   Process = 0   ! select 0, 1 or 2 to represent the spin of the resonance
   Unweighted =.true.
   OffShell_XVV=011! 000: X,V1,V2 on-shell; 010: X,V2 on-shell, V1 off-shell; and so on
   LHEProdFile=""
   ReadLHEFile=.false.
   ConvertLHEFile=.false.
   ReadCSmax=.false.
   GenerateEvents=.false.
   RequestNLeptons = -1
   RequestOSSF=.true.
   LHAPDFString = ""
   LHAPDFMember = 0
   iinterf = -1

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

   DataFile="./data/output"


#if compiler==1
   NumArgs = NArgs()-1
#elif compiler==2
   NumArgs = COMMAND_ARGUMENT_COUNT()
#endif

   CountArg = 0
   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    if( arg(1:4).eq."help" .or. arg(1:4).eq."Help" .or. arg(1:4).eq."HELP") then
        call PrintCommandLineArgs()
        stop
    elseif( arg(1:9).eq."Collider=" ) then
        read(arg(10:10),*) Collider
        CountArg = CountArg + 1
    elseif( arg(1:7).eq."PDFSet=" ) then
        read(arg(8:10),*) PDFSet
        CountArg = CountArg + 1
    elseif( arg(1:7).eq."LHAPDF=" ) then
        read(arg(8:107),"(A)") LHAPDFString
        CountArg = CountArg + 1
    elseif( arg(1:10).eq."LHAPDFMem=" ) then
        read(arg(11:13),*) LHAPDFMember
        CountArg = CountArg + 1
    elseif( arg(1:6).eq."MReso=" ) then
        read(arg(7:12),*) M_Reso
        M_Reso = M_Reso*GeV
        CountArg = CountArg + 1
    elseif( arg(1:9).eq."VegasNc0=" ) then
        read(arg(10:17),*) VegasNc0
        CountArg = CountArg + 1
    elseif( arg(1:9).eq."VegasNc1=" ) then
        read(arg(10:17),*) VegasNc1
        CountArg = CountArg + 1
    elseif( arg(1:9).eq."VegasNc2=" ) then
        read(arg(10:17),*) VegasNc2
        CountArg = CountArg + 1
    elseif( arg(1:9).eq."PChannel=" ) then
        read(arg(10:11),*) PChannel
        CountArg = CountArg + 1
    elseif( arg(1:9).eq."DataFile=" ) then
        read(arg(10:109),"(A)") DataFile
        CountArg = CountArg + 1
    elseif( arg(1:8).eq."Process=" ) then
        read(arg(9:10),*) Process
        CountArg = CountArg + 1
    elseif( arg(1:11).eq."DecayMode1=" ) then
        read(arg(12:13),*) DecayMode1
        CountArg = CountArg + 1
    elseif( arg(1:11).eq."DecayMode2=" ) then
        read(arg(12:13),*) DecayMode2
        CountArg = CountArg + 1
    elseif( arg(1:6).eq."TopDK=" ) then
        read(arg(7:7),*) TopDecays
        CountArg = CountArg + 1
    elseif( arg(1:7) .eq."OffXVV=" ) then
        read(arg(8:10),*) OffShell_XVV
        CountArg = CountArg + 1
    elseif( arg(1:12).eq."FilterNLept=" ) then
        read(arg(13:13),*) RequestNLeptons
        CountArg = CountArg + 1
    elseif( arg(1:11) .eq."FilterOSSF=" ) then
        read(arg(12:12),*) iargument
        if( iargument.eq.1 ) then
            RequestOSSF = .true.
        else
            RequestOSSF = .false.
        endif        
        CountArg = CountArg + 1
    elseif( arg(1:11) .eq."Unweighted=" ) then
        read(arg(12:12),*) iargument
        if( iargument.eq.0 ) then
            Unweighted = .false.
        else
            Unweighted = .true.
        endif
        CountArg = CountArg + 1
    elseif( arg(1:7) .eq."Interf=" ) then
        read(arg(8:8),*) iinterf
        CountArg = CountArg + 1
    elseif( arg(1:8) .eq."ReadLHE=" ) then
        read(arg(9:500),"(A)") LHEProdFile
        ReadLHEFile=.true.
        CountArg = CountArg + 1
    elseif( arg(1:11) .eq."ConvertLHE=" ) then
        read(arg(12:500),"(A)") LHEProdFile
        ConvertLHEFile=.true.
        CountArg = CountArg + 1
    elseif( arg(1:9) .eq."ReadCSmax" ) then
        ReadCSmax=.true.
        CountArg = CountArg + 1
    elseif( arg(1:9) .eq."GenEvents" ) then
        GenerateEvents=.true.
        Unweighted=.false.
        CountArg = CountArg + 1
    endif
   enddo

    Mu_Fact = M_Reso! setting pdf scale to resonance mass

    if( CountArg.ne.NumArgs ) then
        print *, "unknown command line argument"
        stop
    endif

    if (Process.eq.0) PChannel = 0   !only gluons
    if (Process.eq.1 .or. Process.eq.50 .or. Process.eq.60) PChannel = 1   !only quarks

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

#if useLHAPDF==1
    if( LHAPDFString.eq."" ) then
       print *, "Need to specify pdf file name in command line argument LHAPDF"
       stop    
    endif
#endif    
    
!     if( ((OffShellV1).or.(OffShellV2).or.(OffShellReson)) ) then
!         print *, "off shell Z/W's only allowed for spin 0,2 resonance"
! !         stop
!     endif

    if( ConvertLHEFile ) then
       DecayMode2 = DecayMode1
    endif 

    if( Process.eq.50 ) then
        DecayMode2=DecayMode1
        if( Collider.eq.2 ) then
          print *, "Collider 2 not available for VH"
          stop
        endif
    endif
    if( (IsAZDecay(DecayMode1).eqv..false.) .and. (Collider.ne.1) ) then
      print *, "WH with Collider 1 only"
      stop
    endif

    
    if( (TopDecays.ne.0) .and. (Process.eq.80) ) then! TTBH
       if( TopDecays.ne.1 ) call Error("TopDecays=2,3,4 are no longer supported. Use DecayMode1/2.")
       if( .not. IsAWDecay(DecayMode1) ) call Error("Invalid DecayMode1 for top decays")
       if( .not. IsAWDecay(DecayMode2) ) call Error("Invalid DecayMode2 for top decays")
!        if( DecayMode1.eq.4 .and. DecayMode2.eq.4 ) then
          TopDecays = 1
!        elseif( DecayMode1.eq.5 .and. DecayMode2.eq.5 ) then
!           TopDecays = 2
!        elseif( DecayMode1.eq.5 .and. DecayMode2.eq.4 ) then 
!           TopDecays = 3
!        elseif( DecayMode1.eq.4 .and. DecayMode2.eq.6 ) then 
!           TopDecays = 4
!        else
!           call Error("Tau decay modes not yet supported in top decays")
!        endif
    endif

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


    if( ((DecayMode1.eq.0) .and. (DecayMode2.eq.0)) .or.  &
        ((DecayMode1.eq.2) .and. (DecayMode2.eq.2)) .or.  &
        ((DecayMode1.eq.8) .and. (DecayMode2.eq.8)) .or.  &
        ((DecayMode1.eq.9) .and. (DecayMode2.eq.9)) .or.  &
        ((DecayMode1.eq.0) .and. (DecayMode2.eq.8)) .or.  &
        ((DecayMode1.eq.0) .and. (DecayMode2.eq.9)) .or.  &
        ((DecayMode1.eq.2) .and. (DecayMode2.eq.8)) .or.  &
        ((DecayMode1.eq.2) .and. (DecayMode2.eq.9)) .or.  &
        ((DecayMode1.eq.8) .and. (DecayMode2.eq.9)) .or.  &
        ((DecayMode2.eq.0) .and. (DecayMode1.eq.0)) .or.  &
        ((DecayMode2.eq.0) .and. (DecayMode1.eq.8)) .or.  &
        ((DecayMode2.eq.0) .and. (DecayMode1.eq.9)) .or.  &
        ((DecayMode2.eq.2) .and. (DecayMode1.eq.8)) .or.  &
        ((DecayMode2.eq.2) .and. (DecayMode1.eq.9)) .or.  &
        ((DecayMode2.eq.8) .and. (DecayMode1.eq.9))       ) then !  allow interference
            if( iinterf.eq.-1 ) then!  set default interference switch
                if( M_Reso.gt.2d0*M_Z ) then
                    includeInterference = .false.
                else
                    includeInterference = .true.
                endif
            elseif( iinterf.eq.0 ) then! overwrite default and use user selection
                includeInterference = .false.
            else
                includeInterference = .true.
            endif
    else
        includeInterference = .false.   ! no interference if decay mode does not allow 4 same flavor leptons
    endif

    if( IsAZDecay(DecayMode1) .and. IsAWDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

!     if( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
!        print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
!        stop
!     endif

    if( IsAWDecay(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAWDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAPhoton(DecayMode1) .and. (IsAZDecay(DecayMode2) .or. IsAWDecay(DecayMode2)) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       print *, " Please try swapping the decay modes."
       stop
    endif

    if( IsAPhoton(DecayMode2) .and. (IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1)) ) then! require OffXVV=010 for Z+photon
       if( .not. ((.not.OffShellReson) .and. OffShellV1 .and. (.not.OffShellV2)) ) then
          print *, "OffXVV has to be 010 for Z+photon production"
          stop
       endif
    endif

    if( IsAPhoton(DecayMode2) .and. (.not.IsAPhoton(DecayMode1) .and. .not. IsAZDecay(DecayMode1)) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAPhoton(DecayMode1) .and. (OffShellV1 .or. OffShellV2) ) then
       print *, " Photons have to be on-shell."
       stop
    endif

    if( (DecayMode1.ge.12) .or. (DecayMode2.ge.12) .or. (DecayMode1.lt..0) .or. (DecayMode2.lt.0) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( ReadLHEFile .and. Process.ne.0  ) then
        print *, "ReadLHE option is only allowed for spin-0 resonances"
        stop
    endif
    if( ConvertLHEFile .and. Process.ne.0  ) then
        print *, "ConvertLHE option is only allowed for spin-0 resonances"
        stop
    endif
    if( ReadLHEFile .and. .not. Unweighted ) then
        print *, "ReadLHE option is only allowed for generating unweighted events"
        stop
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

return
END SUBROUTINE







SUBROUTINE InitPDFs()
use ModParameters
implicit none
character :: pdftable*(100)

#if useLHAPDF==1
     call InitPDFset(trim(LHAPDFString))
     call InitPDF(LHAPDFMember)  
#else
     if( PDFSet.eq.1 ) then
        call SetCtq6(4)  ! 4    CTEQ6L1  Leading Order cteq6l1.tbl
     elseif( PDFSet.eq.3 ) then
!         pdftable(:)="./pdfs/NNPDF23_lo_as_0130.LHgrid"
        pdftable(:)="./pdfs/NNPDF30_lo_as_0130.LHgrid"
        call NNPDFDriver(pdftable)
        call NNinitPDF(0)
     endif
#endif


     
     
     
     
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
         NDim = 5
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- Hjj, gluon fusion
      if(Process.eq.61) then
         NDim = 5
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- Hj, gluon fusion
      if(Process.eq.62) then
         NDim = 5+1 !1 for color flow ramdomization in gg > Hg
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- VHiggs
      if(Process.eq.50) then
         NDim = 17
         NDim = NDim + 2 ! sHat integration
         VegasIt1_default = 5
         VegasNc0_default = 10000000
         VegasNc1_default = 500000
         VegasNc2_default = 10000
      endif
      !- ttbar+H
      if(Process.eq.80) then
         call InitProcess_TTBH()
         NDim = 12
         NDim = NDim + 2 ! sHat integration
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
         VegasIt1_default = 5
         VegasNc0_default =  100000
         VegasNc1_default =  500000
         VegasNc2_default =    1000
      endif

      if( unweighted ) then
          NDim = NDim + 1  ! random number which decides if event is accepted
      endif
      NDim = NDim + 1  ! MC sampling for gg and qqb channel

END SUBROUTINE








SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModCrossSection
use ModCrossSection_HJJ
use ModCrossSection_TTBH
use ModCrossSection_BBBH
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
    else
      call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events
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
    else
      call vegas1(EvalWeighted,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events
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
    print *, ""

    call cpu_time(time_end)


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
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22)
real(8) :: dum, RES(-5:5,-5:5),ResFrac(-5:5,-5:5),TotalXSec
integer :: i, i1, j1,PChannel_aux, PChannel_aux1,NHisto,RequEvents(-5:+5,-5:+5)
include 'csmaxvalue.f'
integer :: flav1,flav2,StatusPercent
integer :: VegasSeed


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
    if( Process.eq.80 ) call vegas(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.90 ) call vegas(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.60 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.61 ) call vegas(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)


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
    if( Process.eq.80 ) call vegas1(EvalWeighted_TTBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.90 ) call vegas1(EvalWeighted_BBBH,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.60 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)
    if( Process.eq.61 ) call vegas1(EvalWeighted_HJJ,VG_Result,VG_Error,VG_Chi2)





elseif(unweighted.eqv..true.) then  !----------------------- unweighted events

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
                elseif( Process.eq.61 ) then
                    dum = EvalUnWeighted_HJJ(yRnd,.false.,(/-99,-99/),RES)
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
   write(io_stdout,"(1X,A,F8.3)") "Total xsec",TotalXSec


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
              elseif( Process.eq.61 ) then
                  dum = EvalUnWeighted_HJJ(yRnd,.true.,(/i1,j1/),RES)
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
    print *, " gg/qqb ratio = ", dble(AccepCounter_part(0,0))/dble(sum(AccepCounter_part(:,:))-AccepCounter_part(0,0))
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90."
        write(io_stdout, *) "ALERT: The number of rejected events with too small CSMAX exceeds 1%."
        write(io_stdout, *) "       Increase CSMAX in main.F90."
    endif
    write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)

                     
  endif! unweighted
  
  


return
END SUBROUTINE





SUBROUTINE StartVegas_BETA(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22)
real(8) :: dum, RES(-5:5,-5:5)
integer :: i, i1, PChannel_aux, PChannel_aux1,n1,n2,n3,n4
include 'csmaxvalue.f'
integer :: VegasSeed


if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc1.eq.-1 ) VegasNc1 = 5000



   warmup = .false.
   itmx = 3  !VegasIt1        ! overwrite this for development phase
   ncall= 500000 !VegasNc1   ! overwrite this for development phase

   if( (Process.eq.2 .and. PChannel.ne.0) .or. (Process.eq.1) ) then
        print *, " ERROR:"
        print *, " You are running the beta version for generating events."
        print *, " This version does not yet support qqb initial states."
        print *, " Please set  useBetaVersion=.false.  in main.F90."
        print *, ""
        stop
   endif


   PChannel_aux = PChannel

if (unweighted.eqv..false.) then  !----------------------- weighted events

  call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events



elseif(unweighted.eqv..true.) then  !----------------------- unweighted events


    print *, " scanning the integrand"

    call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)

!print *, globalMin
!print *, globalMax
!pause
do n1=1,NPart-1
do n2=1,NPart-1
!do n3=1,NPart-1
!do n4=1,NPart-1
     if( PartitionMax(n1,n2).eq.-1d99 ) then
          print *, " unfilled partition: ",n1,n2
         PartitionMax(n1,n2) = globalMax
     endif
!      if( PartitionMax(n1,n2,n3,n4).eq.-1d99 ) then
! !          print *, "unfilled partitions",n1,n2,n3,n4
!          PartitionMax(n1,n2,n3,n4) = globalMax
!      endif
!enddo
!enddo
enddo
enddo

!pause
   call ClearHisto()
   
   print *, " generating events"
   call cpu_time(time_start)
   do while( AccepCounter.lt.VegasNc1 )! generate a fixed number of VegasNc1 events
!   do i=1,ncall
         call random_number(yRnd)
! print *, "generated event no.",AccepCounter
         dum = EvalUnWeighted_BETA(yRnd)
   enddo
   call cpu_time(time_end)

   print *, "globalMin",globalMin
   print *, "globalMax",globalMax
   print *, "accepted",AccepCounter
   print *, "rejected",RejeCounter
   print *, "efficiency",dble(AccepCounter)/dble(RejeCounter)
   print *, "time (sec)",time_end-time_start
   print *, "rate (events/sec)",dble(AccepCounter)/(time_end-time_start)

   write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)

endif



return
END SUBROUTINE






! SUBROUTINE StartReadLHE(VG_Result,VG_Error)
! use ModCrossSection
! use ModKinematics
! use ModParameters
! use ModMisc
! implicit none
! include 'csmaxvalue.f'
! integer,parameter :: maxpart=15!=max.partons; this parameter should match the one in WriteOutEvent of mod_Kinematics
! real(8) :: VG_Result,VG_Error,VG_Chi2
! real(8) :: yRnd(1:22),Res,dum,EMcheck(1:4)
! real(8) :: AcceptedEvent(1:4,1:maxpart),Ehat
! real(8) :: MomExt(1:4,1:maxpart),MomHiggs(1:4),MomParton(1:4,1:maxpart),Mass(1:maxpart),pH2sq
! integer :: tries, nParticle, MY_IDUP(1:7+maxpart), ICOLUP(1:2,1:7+maxpart),IntExt
! character(len=*),parameter :: POWHEG_Fmt0 = "(5X,I2,A160)"
! character(len=*),parameter :: POWHEG_Fmt1 = "(5X,I3,4X,I2,4X,I2,4X,I2,2X,I4,2X,I4,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9)"
! character(len=*),parameter :: JHUGen_Fmt0 = "(I2,A160)"
! character(len=*),parameter :: JHUGen_Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
! character(len=*),parameter :: JHUGen_old_Fmt0 = "(2X,I2,A160)"
! character(len=*),parameter :: JHUGen_old_Fmt1 = "(I3,X,I2,X,I2,X,I2,X,I3,X,I3,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7)"
! character(len=*),parameter :: MadGra_Fmt0 = "(I2,A160)"
! character(len=*),parameter :: MadGra_Fmt1 = "(7X,I3,2X,I3,3X,I2,3X,I2,3X,I3,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1F3.0,X,1F3.0)"
! character(len=150) :: InputFmt0,InputFmt1
! logical :: FirstEvent,M_ResoSet
! integer :: nline,intDummy,Nevent
! integer :: LHE_IDUP(1:maxpart+3),   LHE_ICOLUP(1:2,1:maxpart+3),   LHE_MOTHUP(1:2,1:maxpart+3)
! integer :: LHE_IDUP_Part(1:maxpart),LHE_ICOLUP_Part(1:2,1:maxpart),LHE_MOTHUP_Part(1:2,1:maxpart+3)
! integer :: EventNumPart,nparton
! character(len=160) :: FirstLines,EventInfoLine,PDFLine
! character(len=160) :: EventLine(1:maxpart+3)
! integer :: VegasSeed,stat,n
! integer,parameter :: InputLHEFormat = 1  !  1=POWHEG, 2=JHUGen (old format), 3=JHUGen (new format), 4=MadGraph
! 
! 
! if(InputLHEFormat.eq.1) then
!   InputFmt0 = trim(POWHEG_Fmt0)
!   InputFmt1 = trim(POWHEG_Fmt1)
! elseif(InputLHEFormat.eq.4) then
!   InputFmt0 = trim(MadGra_Fmt0)
!   InputFmt1 = trim(MadGra_Fmt1)
! elseif(InputLHEFormat.eq.2) then
!   InputFmt0 = trim(JHUGen_old_Fmt0)
!   InputFmt1 = trim(JHUGen_old_Fmt1)
! else
!   InputFmt0 = trim(JHUGen_Fmt0)
!   InputFmt1 = trim(JHUGen_Fmt1)
! endif
! 
! 
! 
! if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
! if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
! if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
! if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2
! 
! 
!     call InitRandomSeed(TheSeed)
! 
! !    search for line with first event
!      FirstEvent = .false.
!      M_ResoSet=.false.
!      do while ( .not.FirstEvent )
!         read(16,fmt="(A160)",IOSTAT=stat,END=99) FirstLines
!         if( FirstLines(1:5).eq."hmass" ) then 
!                read(FirstLines(6:13),fmt="(F7.0)") M_Reso
!                M_Reso = M_Reso*GeV!  convert to units of 100GeV
!                M_ResoSet=.true.
!         endif
!         if( FirstLines(1:7).eq."<event>" ) then 
!                FirstEvent=.true.
!         else
!             if( importExternal_LHEinit ) then
!                 if( FirstLines(1:17).eq."<LesHouchesEvents" .or. FirstLines(1:4).eq."<!--" ) then
!                 else
!                   write(io_LHEOutFile,"(A)") trim(firstlines)
!                 endif
!             endif
!         endif
!      enddo
!      if( .not. M_ResoSet ) then
!         write(io_stdout,"(2X,A,1F7.2)")  "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
!         write(io_LogFile,"(2X,A,1F7.2)") "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
!      else
!         write(io_stdout,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
!         write(io_LogFile,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
!      endif
!      write(io_stdout,"(A)") ""
!      write(io_LogFile,"(A)") ""
! 
! 
! 
!       print *, " finding maximal weight with ",VegasNc0," points"
!       VG = zero
!       CSmax = zero
!       Ehat = M_Reso! fixing Ehat to M_Reso which should determine the max. of the integrand
!       do tries=1,VegasNc0
!           call random_number(yRnd)
!           dum = EvalUnWeighted_withoutProduction(yRnd,.false.,EHat,Res,AcceptedEvent,MY_IDUP(1:9),ICOLUP(1:2,1:9))
!       enddo
!       csmax(0,0)   = 1.5d0*csmax(0,0)    !  savety buffer
! 
! 
! 
! 
! 
!      print *, " generating events"
!      EvalCounter = 0
!      AccepCounter = 0
!      RejeCounter = 0
!      AccepCounter_part = 0
!      call cpu_time(time_start)
!      NEvent=0
!      do while ( .true. ) 
!          NEvent=NEvent + 1
!          read(16,fmt=InputFmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
! !        read event lines
!          do nline=1,EventNumPart
!             read(16,fmt="(A160)") EventLine(nline)
!          enddo
!          if( EventNumPart.lt.3 .or. EventNumPart.gt.maxpart ) then
!             call Error("Number of particles in LHE input exceeds allowed limit",EventNumPart)
!          endif
! 
!  !       convert event lines into variables assuming that the Higgs resonance has ID 25
!          nparton = 0
!          do nline=1,EventNumPart
!             read(EventLine(nline),fmt=InputFmt1) LHE_IDUP(nline),IntExt,LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline)
!             if( IntExt.ne.-1 .and. IntExt.ne.+1 ) cycle ! remove internal particles
!             MomExt(1:4,nline) = MomExt(1:4,nline)*GeV!  convert to units of 100GeV
!             Mass(nline) = Mass(nline)*GeV            !  convert to units of 100GeV
!             if( abs(LHE_IDUP(nline)).eq.25 ) then!   select the Higgs (ID=25, h0)
!                   MomHiggs(1:4) = MomExt(1:4,nline)
!                   pH2sq = dsqrt(abs(MomHiggs(1:4).dot.MomHiggs(1:4)))
! !             elseif( (abs(LHE_IDUP(nline)).eq.21) .or. (abs(LHE_IDUP(nline)).ge.1 .and. abs(LHE_IDUP(nline)).le.5) ) then!   select the gluons (ID=21) or quarks (ID=1,.,5)
!             else
!                   nparton = nparton + 1
!                   MomParton(1:4,nparton) = MomExt(1:4,nline)
!                   LHE_IDUP_Part(nparton) = LHE_IDUP(nline)
!                   LHE_ICOLUP_Part(1:2,nparton) = LHE_ICOLUP(1:2,nline)
!                   LHE_MOTHUP_Part(1:2,nparton) = LHE_MOTHUP(1:2,nline)
! !              else
! !                   call Error("Unknown particle in LHE input file",LHE_IDUP(nline))
!             endif
!          enddo
! !        read optional pdf line
!          read(16,fmt="(A160)",IOSTAT=stat,END=99) PDFLine(1:160)
! !          if( .not. PDFLine(1:4).eq."#pdf") then
!          if( .not. (PDFLine(1:4).eq."#pdf" .or. PDFLine(1:5).eq."#rwgt")) then         
!              PDFLine(:)=""
!              backspace(16)
!          endif
! 
! 
! ! !        reject event if M_Reso and pH2sq deviate by more than 20 GeV 
! !          if( abs( M_Reso - pH2sq) .gt. 20d0*GeV ) then
! !               write(io_stdout,"(2X,A,2F10.4)")  "WARNING: Higgs mass and momentum squared deviate by more than 20 GeV!",pH2sq*100d0,M_Reso*100d0
! !               write(io_LogFile,"(2X,A,2F10.4)") "WARNING: Higgs mass and momentum squared deviate by more than 20 GeV!",pH2sq*100d0,M_Reso*100d0
! !               cycle
! !          endif
! 
!          EMcheck(1:4) = MomParton(1:4,1) + MomParton(1:4,2) - MomHiggs(1:4)
!          do n=3,nparton
!              EMcheck(1:4) = EMcheck(1:4) - MomParton(1:4,n)
!          enddo
!          if( any(abs(EMcheck(1:4)).gt.1d0*GeV) ) then
!               print *, EMcheck(1:4)
!               call Error("energy momentum violation while reading LHE production momenta.",nparton)
!          endif
! 
! !         accept/reject sampling for H->VV decay contribution
!           EHat = pH2sq
!           do tries=1,5000000
!               call random_number(yRnd)
!               dum = EvalUnWeighted_withoutProduction(yRnd,.true.,Ehat,RES,AcceptedEvent,MY_IDUP(1:9),ICOLUP(1:2,1:9))
!               if( Res.ne.0d0 ) exit
!           enddo
!           if( Res.ne.0d0 ) then ! decay event was accepted
!              call boost(AcceptedEvent(1:4,1),MomHiggs(1:4),pH2sq)
!              call boost(AcceptedEvent(1:4,2),MomHiggs(1:4),pH2sq)
!              call boost(AcceptedEvent(1:4,3),MomHiggs(1:4),pH2sq)
!              call boost(AcceptedEvent(1:4,4),MomHiggs(1:4),pH2sq)
! 
!              MY_IDUP(1) = convertLHEreverse(LHE_IDUP_Part(1))! overwrite inital ID with external ones
!              MY_IDUP(2) = convertLHEreverse(LHE_IDUP_Part(2))
!              ICOLUP(1:2,1:2) = LHE_ICOLUP_Part(1:2,1:2)!  overwrite inital colors with external ones
!              do n=3,nparton
!                 MY_IDUP(7+n)= convertLHEreverse(LHE_IDUP_Part(n))
!                 ICOLUP(1:2,7+n) = LHE_ICOLUP_Part(1:2,n)!  overwrite final colors with external ones
!              enddo
!              call WriteOutEvent((/MomParton(1:4,1),MomParton(1:4,2),AcceptedEvent(1:4,1),AcceptedEvent(1:4,2),AcceptedEvent(1:4,3),AcceptedEvent(1:4,4)/), &
!                                 MY_IDUP(1:7+nparton),ICOLUP(1:2,1:7+nparton),MomFSPartons=MomParton(1:4,3:nparton),EventInfoLine=EventInfoLine,PDFLine=PDFLine,MOTHUP_Parton=LHE_MOTHUP_Part)
!              if( mod(AccepCounter,5000).eq.0 ) then
!                   call cpu_time(time_int)
!                   write(io_stdout,*)  NEvent," events accepted (",time_int-time_start, ") seconds"
!                   write(io_LogFile,*) NEvent," events accepted (",time_int-time_start, ") seconds"
!              endif
!           else! decay event was not accepted after ncall evaluations, read next production event
!              print *, "rejected event after ",tries-1," evaluations"
!              AlertCounter = AlertCounter + 1 
!           endif
! 
! !        skip event lines 
! !          read(16,fmt="(A7)",IOSTAT=stat,END=99) FirstLines! skip <\event>
! ! !          if( stat.lt.0 ) exit
! !          read(16,fmt="(A30)",IOSTAT=stat,END=99) FirstLines!   skip <event> or </LesHouchesEvents>
! !          if( FirstLines(1:30).eq."</LesHouchesEvents>" ) exit
! !          if( NEvent.eq. VegasNc1 ) exit
! 
! 
! !        read optional lines
!          FirstEvent = .true.
!          do while (.true.) 
!               read(16,fmt="(A160)",IOSTAT=stat,END=99) EventInfoLine(1:160)
!               if(EventInfoLine(1:30).eq."</LesHouchesEvents>") then
!                   goto 99
!               elseif( EventInfoLine(1:8).eq."<event>" ) then
!                   exit
!               else!if there are "#" comments
!                   if( FirstEvent ) then 
!                       backspace(io_LHEOutFile)! remove "</event>" from WriteOutEvent
!                       FirstEvent = .false.
!                   endif
!                   write(io_LHEOutFile,fmt="(A)") trim(EventInfoLine)
!               endif
!          enddo         
! 
!          
!      enddo
! 99   continue
!      call cpu_time(time_end)
! 
! 
! 
! 
!     write(io_stdout,*) ""
!     write(io_stdout,*) "Evaluation Counter: ",EvalCounter
!     write(io_stdout,*) "Acceptance Counter: ",AccepCounter
!     write(io_stdout,*) "Rejection  Counter: ",RejeCounter
!     write(io_stdout,*) " Alert  Counter: ",AlertCounter
!     if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
!         write(io_stdout,*) "ALERT: The number of rejected events exceeds 1%."
!         write(io_stdout,*) "       Increase CSMAX in main.F90 or VegasNc1."
!     endif
!    write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)
! 
!     write(io_LogFile,*) ""
!     write(io_LogFile,*) "Evaluation Counter: ",EvalCounter
!     write(io_LogFile,*) "Acceptance Counter: ",AccepCounter
!     write(io_LogFile,*) "Rejection  Counter: ",RejeCounter
!     write(io_LogFile,*) " Alert  Counter: ",AlertCounter
!     if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
!         write(io_LogFile,*) "ALERT: The number of rejected events exceeds 1%."
!         write(io_LogFile,*) "       Increase CSMAX in main.F90 or VegasNc1."
!     endif
!    write(io_LogFile,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)
! 
! 
! 
! return
! END SUBROUTINE





SUBROUTINE StartReadLHE_NEW(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
use ModMisc
implicit none
include 'csmaxvalue.f'
integer,parameter :: maxpart=30!=max.part particles in LHE file; this parameter should match the one in WriteOutEvent of mod_Kinematics
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),Res,dum,EMcheck(1:4)
real(8) :: HiggsDK_Mom(1:4,4:9),Ehat
real(8) :: MomExt(1:4,1:maxpart),MomHiggs(1:4),Mass(1:maxpart),pH2sq
integer :: tries, nParticle, HiggsDK_IDUP(1:9), ICOLUP(1:2,1:7+maxpart),LHE_IntExt(1:7+maxpart),HiggsDK_ICOLUP(1:2,1:9)
character(len=*),parameter :: POWHEG_Fmt0 = "(5X,I2,A160)"
character(len=*),parameter :: POWHEG_Fmt1 = "(5X,I3,4X,I2,4X,I2,4X,I2,2X,I4,2X,I4,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9)"
character(len=*),parameter :: JHUGen_Fmt0 = "(I2,A160)"
character(len=*),parameter :: JHUGen_Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
character(len=*),parameter :: JHUGen_old_Fmt0 = "(2X,I2,A160)"
character(len=*),parameter :: JHUGen_old_Fmt1 = "(I3,X,I2,X,I2,X,I2,X,I3,X,I3,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7)"
character(len=*),parameter :: MadGra_Fmt0 = "(I2,A160)"
character(len=*),parameter :: MadGra_Fmt1 = "(7X,I3,2X,I3,3X,I2,3X,I2,3X,I3,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1F3.0,X,1F3.0)"
character(len=150) :: InputFmt0,InputFmt1
logical :: FirstEvent,M_ResoSet,Ga_ResoSet
integer :: nline,intDummy,Nevent
integer :: LHE_IDUP(1:maxpart),LHE_ICOLUP(1:2,1:maxpart),LHE_MOTHUP(1:2,1:maxpart)
integer :: EventNumPart
character(len=160) :: FirstLines,EventInfoLine,OtherLines
character(len=160) :: EventLine(1:maxpart)
integer :: n,stat,iHiggs,VegasSeed
integer,parameter :: InputLHEFormat = 1  !  1=POWHEG, 2=JHUGen (old format), 3=JHUGen (new format), 4=MadGraph



if(InputLHEFormat.eq.1) then
  InputFmt0 = trim(POWHEG_Fmt0)
  InputFmt1 = trim(POWHEG_Fmt1)
elseif(InputLHEFormat.eq.4) then
  InputFmt0 = trim(MadGra_Fmt0)
  InputFmt1 = trim(MadGra_Fmt1)
elseif(InputLHEFormat.eq.2) then
  InputFmt0 = trim(JHUGen_old_Fmt0)
  InputFmt1 = trim(JHUGen_old_Fmt1)
else
  InputFmt0 = trim(JHUGen_Fmt0)
  InputFmt1 = trim(JHUGen_Fmt1)
endif



if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2

!    search for line with first event
     FirstEvent = .false.
     M_ResoSet=.false.
     Ga_ResoSet=.false.
     do while ( .not.FirstEvent )
        read(16,fmt="(A160)",IOSTAT=stat,END=99) FirstLines
        if( FirstLines(1:5).eq."hmass" ) then 
               read(FirstLines(6:16),*) M_Reso
               M_Reso = M_Reso*GeV!  convert to units of 100GeV
               M_ResoSet=.true.
        endif
        if( FirstLines(1:6).eq."hwidth" ) then 
               read(FirstLines(7:16),*) Ga_Reso
               Ga_Reso = Ga_Reso*GeV!  convert to units of 100GeV
               Ga_ResoSet=.true.
        endif
        if( FirstLines(1:3).eq."-->" ) then
            write(io_LHEOutFile ,"(A)") ""
            write(io_LHEOutFile ,"(A)") "JHUGen Resonance parameters used for event generation:"
            write(io_LHEOutFile ,"(A,F6.1,A)") "hmass  ",M_Reso*100d0,"        ! Higgs boson mass"
            write(io_LHEOutFile ,"(A,F10.5,A)") "hwidth",Ga_Reso*100d0,"      ! Higgs boson width"
            write(io_LHEOutFile ,"(A)") ""
        endif
        if( FirstLines(1:7).eq."<event>" ) then 
               FirstEvent=.true.
        else
            if( importExternal_LHEinit ) then
                if( FirstLines(1:17).eq."<LesHouchesEvents" .or. FirstLines(1:4).eq."<!--" ) then
                else
                  write(io_LHEOutFile,"(A)") trim(firstlines)
                endif
            endif
        endif
     enddo
     if( .not. M_ResoSet ) then
        write(io_stdout,"(2X,A,1F7.2)")  "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
        write(io_LogFile,"(2X,A,1F7.2)") "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
     else
        write(io_stdout,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
        write(io_LogFile,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
     endif
     if( .not. Ga_ResoSet ) then
        write(io_stdout,"(2X,A,1F7.2)")  "ERROR: Higgs width could not be read from LHE input file. Assuming default value",Ga_Reso*100d0
        write(io_LogFile,"(2X,A,1F7.2)") "ERROR: Higgs width could not be read from LHE input file. Assuming default value",Ga_Reso*100d0
     else
        write(io_stdout,"(2X,A,1F9.5,A)") "A Higgs width of ",Ga_Reso*100d0," GeV was determined from the LHE input file."
        write(io_LogFile,"(2X,A,1F9.5,A)") "A Higgs width of ",Ga_Reso*100d0," GeV was determined from the LHE input file."
     endif
     write(io_stdout,"(A)") ""
     write(io_LogFile,"(A)") ""


      print *, " finding maximal weight with ",VegasNc0," points"
      VG = zero
      CSmax = zero
      Ehat = M_Reso! fixing Ehat to M_Reso which should determine the max. of the integrand
      do tries=1,VegasNc0
          call random_number(yRnd)
          dum = EvalUnWeighted_withoutProduction(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP,HiggsDK_ICOLUP)
      enddo
      csmax(0,0)   = 1.5d0*csmax(0,0)    !  savety buffer



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
         read(16,fmt=InputFmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
!        read event lines
         do nline=1,EventNumPart
            read(16,fmt="(A160)") EventLine(nline)
         enddo
         if( EventNumPart.lt.3 .or. EventNumPart.gt.maxpart ) then
            call Error("Number of particles in LHE input exceeds allowed limit",EventNumPart)
         endif
         do nline=1,EventNumPart
            read(EventLine(nline),fmt=InputFmt1) LHE_IDUP(nline),LHE_IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline)
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
          do tries=1,5000000
              call random_number(yRnd)
              dum = EvalUnWeighted_withoutProduction(yRnd,.true.,Ehat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP,HiggsDK_ICOLUP)
              if( Res.ne.0d0 ) exit
          enddo
          if( Res.gt.0d0 ) then ! decay event was accepted
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
             
             call WriteOutEvent_NEW(EventNumPart,LHE_IDUP,LHE_IntExt,LHE_MOTHUP,LHE_ICOLUP,MomExt,HiggsDK_Mom,Mass,iHiggs,HiggsDK_IDUP,HiggsDK_ICOLUP,EventInfoLine)

             if( mod(AccepCounter,5000).eq.0 ) then
                  call cpu_time(time_int)
                  write(io_stdout,*)  NEvent," events accepted (",time_int-time_start, ") seconds"
                  write(io_LogFile,*) NEvent," events accepted (",time_int-time_start, ") seconds"
             endif

          elseif( Res.eq.0d0 ) then ! decay event was not accepted after ncall evaluations, read next production event
             print *, "Rejected event after ",tries-1," evaluations"
             AlertCounter = AlertCounter + 1 
          endif

!        read optional lines
         FirstEvent = .true.
         tries = 0
         do while (.true.) 
              tries = tries +1 
              read(16,fmt="(A160)",IOSTAT=stat,END=99) OtherLines(1:160)
              if(OtherLines(1:30).eq."</LesHouchesEvents>") then
                  if( RequestNLeptons.gt.0 ) write(io_LHEOutFile,"(A,1F6.2,A)") "<!-- Lepton filter efficiency:",dble(AccepCounter)/dble(NEvent)*100d0," % -->"
                  goto 99
              elseif( OtherLines(1:8).eq."</event>" .and. Res.gt.0d0 ) then
                  write(io_LHEOutFile,"(A)") "</event>"
              elseif( OtherLines(1:8).eq."<event>" ) then
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
real(8) :: MomExt(1:4,1:maxpart),MomShift(1:4,1:maxpart),MomHiggs(1:4),MomParton(1:4,1:maxpart),Mass(1:maxpart),Spin(1:maxpart),Lifetime(1:maxpart)
integer :: tries, nParticle, MY_IDUP(1:7+maxpart), ICOLUP(1:2,1:7+maxpart),IntExt(1:7+maxpart),convertparent
character(len=*),parameter :: POWHEG_Fmt0 = "(5X,I2,A160)"
character(len=*),parameter :: POWHEG_Fmt1 = "(5X,I3,4X,I2,4X,I2,4X,I2,2X,I4,2X,I4,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE12.5,1X,1PE10.3)"
character(len=*),parameter :: JHUGen_Fmt0 = "(I2,A120)"
character(len=*),parameter :: JHUGen_Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
character(len=*),parameter :: JHUGen_old_Fmt0 = "(2X,I2,A160)"
character(len=*),parameter :: JHUGen_old_Fmt1 = "(I3,X,I2,X,I2,X,I2,X,I3,X,I3,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7)"
character(len=*),parameter :: MadGra_Fmt0 = "(I2,A120)"
character(len=*),parameter :: MadGra_Fmt1 = "(7X,I3,2X,I3,3X,I2,3X,I2,3X,I3,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1F3.0,X,1F3.0)"
character(len=150) :: InputFmt0,InputFmt1
logical :: FirstEvent,M_ResoSet,Ga_ResoSet
integer :: nline,intDummy,Nevent
integer :: LHE_IDUP(1:maxpart+3),   LHE_ICOLUP(1:2,1:maxpart+3),   LHE_MOTHUP(1:2,1:maxpart+3)
integer :: LHE_IDUP_Part(1:maxpart),LHE_ICOLUP_Part(1:2,1:maxpart),LHE_MOTHUP_Part(1:2,1:maxpart+3)
integer :: EventNumPart,nparton
character(len=160) :: FirstLines
character(len=120) :: EventInfoLine,PDFLine
character(len=160) :: EventLine(1:maxpart+3)
integer :: VegasSeed,i,stat,DecayParticles(1:2)
integer, dimension(:), allocatable :: gen_seed
integer,parameter :: InputLHEFormat = 1  !  1=POWHEG, 2=JHUGen (old format), 3=JHUGen (new format), 4=MadGraph


if(InputLHEFormat.eq.1) then
  InputFmt0 = trim(POWHEG_Fmt0)
  InputFmt1 = trim(POWHEG_Fmt1)
elseif(InputLHEFormat.eq.4) then
  InputFmt0 = trim(MadGra_Fmt0)
  InputFmt1 = trim(MadGra_Fmt1)
elseif(InputLHEFormat.eq.2) then
  InputFmt0 = trim(JHUGen_old_Fmt0)
  InputFmt1 = trim(JHUGen_old_Fmt1)
else
  InputFmt0 = trim(JHUGen_Fmt0)
  InputFmt1 = trim(JHUGen_Fmt1)
endif



if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2


!    search for line with first event
     FirstEvent = .false.
     M_ResoSet=.false.
     do while ( .not.FirstEvent )
        read(16,fmt="(A160)",IOSTAT=stat,END=99) FirstLines
        if( FirstLines(1:5).eq."hmass" ) then 
               read(FirstLines(6:13),fmt="(F7.1)") M_Reso
               M_Reso = M_Reso*GeV!  convert to units of 100GeV
               M_ResoSet=.true.
        endif
        if( FirstLines(1:6).eq."hwidth" ) then 
               read(FirstLines(7:17),fmt="(F9.3)") Ga_Reso
               Ga_Reso = Ga_Reso*GeV!  convert to units of 100GeV
               Ga_ResoSet=.true.
        endif       
        if( FirstLines(1:7).eq."<event>" ) then 
               FirstEvent=.true.
        else
            if( importExternal_LHEinit ) then
                if( FirstLines(1:17).eq."<LesHouchesEvents" .or. FirstLines(1:4).eq."<!--" ) then
                else
                  write(io_LHEOutFile,"(A)") trim(firstlines)
                endif
            endif
        endif
     enddo
     if( .not. M_ResoSet ) then
        write(io_stdout,"(2X,A,1F7.2)")  "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
        write(io_LogFile,"(2X,A,1F7.2)") "ERROR: Higgs mass could not be read from LHE input file. Assuming default value",M_Reso*100d0
     else
        write(io_stdout,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
        write(io_LogFile,"(2X,A,1F7.2,A)") "A Higgs mass of ",M_Reso*100d0," GeV was determined from the LHE input file."
     endif
     if( .not. Ga_ResoSet ) then
        write(io_stdout,"(2X,A,1F7.2)")  "ERROR: Higgs width could not be read from LHE input file. Assuming default value",Ga_Reso*100d0
        write(io_LogFile,"(2X,A,1F7.2)") "ERROR: Higgs width could not be read from LHE input file. Assuming default value",Ga_Reso*100d0
     else
        write(io_stdout,"(2X,A,1F9.5,A)") "A Higgs width of ",Ga_Reso*100d0," GeV was determined from the LHE input file."
        write(io_LogFile,"(2X,A,1F9.5,A)") "A Higgs width of ",Ga_Reso*100d0," GeV was determined from the LHE input file."
     endif
     write(io_stdout,"(A)") ""
     write(io_LogFile,"(A)") ""


     print *, " converting events"
     call cpu_time(time_start)
     NEvent=0
     do while ( .true. ) 
         NEvent=NEvent + 1
         read(16,fmt=InputFmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
!        read event lines
         do nline=1,EventNumPart
            read(16,fmt="(A160)") EventLine(nline)
         enddo
         if( EventNumPart.lt.3 .or. EventNumPart.gt.maxpart ) then
            call Error("Number of particles in LHE input exceeds allowed limit",EventNumPart)
         endif

 !       convert event lines into variables assuming that the Higgs resonance has ID 25
         nparton = 0
         do nline=1,EventNumPart
            read(EventLine(nline),fmt=InputFmt1) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline),Spin(nline),Lifetime(nline)
            if( abs(LHE_IDUP(nline)).eq.25 ) then!   select the Higgs (ID=25, h0)
                  MomHiggs(1:4) = MomExt(1:4,nline)
                  pH2sq = dsqrt(abs(MomHiggs(1:4).dot.MomHiggs(1:4)))
            else
                  nparton = nparton + 1
                  MomParton(1:4,nparton) = MomExt(1:4,nline)
                  LHE_IDUP_Part(nparton) = LHE_IDUP(nline)
                  LHE_ICOLUP_Part(1:2,nparton) = LHE_ICOLUP(1:2,nline)
                  LHE_MOTHUP_Part(1:2,nparton) = LHE_MOTHUP(1:2,nline)
            endif
            
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
         
         write(io_LHEOutFile,"(A)") "<event>"
         write(io_LHEOutFile,fmt=InputFmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
         do nline=1,EventNumPart
            write(io_LHEOutFile,fmt=InputFmt1) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomShift(2,nline),MomShift(3,nline),MomShift(4,nline),MomShift(1,nline),Mass(nline),Spin(nline),Lifetime(nline)
         enddo
         
! !        read optional lines
!          read(16,fmt="(A160)",IOSTAT=stat,END=99) PDFLine(1:160)
!          if( PDFLine(1:1).ne."#" ) then
!              PDFLine(:)=""
!              backspace(16)
!          else
!              write(io_LHEOutFile,fmt="(A)") trim(PDFLine)
!          endif
!          write(io_LHEOutFile,"(A)") "</event>"         
!          
! !        skip event lines 
!          read(16,fmt="(A7)",IOSTAT=stat,END=99) FirstLines! skip <\event>
!          read(16,fmt="(A30)",IOSTAT=stat,END=99) FirstLines!   skip <event> or </LesHouchesEvents>
!          if( FirstLines(1:30).eq."</LesHouchesEvents>" ) exit


!        read optional lines
         do while (.true.) 
              read(16,fmt="(A120)",IOSTAT=stat,END=99) PDFLine(1:120)
              if(PDFLine(1:30).eq."</LesHouchesEvents>") then
                  goto 99
              elseif( PDFLine(1:8).eq."</event>" ) then
                  write(io_LHEOutFile,"(A)") "</event>"
              elseif( PDFLine(1:8).eq."<event>" ) then
                  exit
              else
                  write(io_LHEOutFile,fmt="(A)") trim(PDFLine)
              endif
         enddo
     enddo
     
99   continue
     call cpu_time(time_end)

return
END SUBROUTINE








! SUBROUTINE ScrambleLHE(Infile,Outfile)
! use ModParameters
! use ModMisc
! implicit none
! character(len=*),parameter :: POWHEG_Fmt0 = "(5X,I2,A160)"
! character(len=*),parameter :: POWHEG_Fmt1 = "(5X,I3,4X,I2,4X,I2,4X,I2,2X,I4,2X,I4,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9)"
! character(len=*),parameter :: JHUGen_Fmt0 = "(I2,A160)"
! character(len=*),parameter :: JHUGen_Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
! character(len=*),parameter :: JHUGen_old_Fmt0 = "(2X,I2,A160)"
! character(len=*),parameter :: JHUGen_old_Fmt1 = "(I3,X,I2,X,I2,X,I2,X,I3,X,I3,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7)"
! character(len=*),parameter :: MadGra_Fmt0 = "(I2,A160)"
! character(len=*),parameter :: MadGra_Fmt1 = "(7X,I3,2X,I3,3X,I2,3X,I2,3X,I3,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1F3.0,X,1F3.0)"
! character(len=150) :: InputFmt0,InputFmt1
! logical :: FirstEvent
! character(len=160) :: HeaderLines,EventInfoLine,OtherLines
! character(len=160) :: EventLine(1:maxpart)
! integer,parameter :: InputLHEFormat = 1  !  1=POWHEG, 2=JHUGen (old format), 3=JHUGen (new format), 4=MadGraph
! 
! 
! if(InputLHEFormat.eq.1) then
!   InputFmt0 = trim(POWHEG_Fmt0)
!   InputFmt1 = trim(POWHEG_Fmt1)
! elseif(InputLHEFormat.eq.4) then
!   InputFmt0 = trim(MadGra_Fmt0)
!   InputFmt1 = trim(MadGra_Fmt1)
! elseif(InputLHEFormat.eq.2) then
!   InputFmt0 = trim(JHUGen_old_Fmt0)
!   InputFmt1 = trim(JHUGen_old_Fmt1)
! else
!   InputFmt0 = trim(JHUGen_Fmt0)
!   InputFmt1 = trim(JHUGen_Fmt1)
! endif
! 
! 
!      open(unit=io_LHEOutFile,file=trim(DataFile)//'.lhe',form='formatted',access= 'sequential') ! in1 (part1)
!      open(unit=io_LHEOutFile2,file=trim(DataFile)//'.lhe',form='formatted',access= 'sequential')! in2 (part2)
!      open(unit=io_LHEOutFile3,file=trim(DataFile)//'.new.lhe',form='formatted',access= 'sequential',status='replace')! out
! 
! 
! !    search for line with first event in gg-file
!      FirstEvent = .false.
!      do while ( .not.FirstEvent )
!         read(io_LHEOutFile,fmt="(A160)",IOSTAT=stat,END=99) HeaderLines
!         if( HeaderLines(1:7).eq."<event>" ) then 
!                FirstEvent=.true.
!         endif
!      enddo
!      read(io_LHEOutFile,fmt=InputFmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
!      do nline=1,EventNumPart!  read event lines
!         read(io_LHEOutFile,fmt="(A160)") EventLine(nline)
!      enddo
!      do nline=1,EventNumPart    
!         read(EventLine(nline),fmt=InputFmt1) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline),Spin(nline),Lifetime(nline)
!      enddo
!      if( LHE_IDUP(1).eq.21 .and. LHE_IDUP(2).eq.21 ) then!  gg initial state
!      endif
! 
!      
!      
!      
!      
!      
!      
!      
! !    search for line with first event in gg-file
!      FirstEvent = .false.
!      do while ( .not.FirstEvent )
!         read(io_LHEOutFile2,fmt="(A160)",IOSTAT=stat,END=99) HeaderLines
!         if( HeaderLines(1:7).eq."<event>" ) then 
!                FirstEvent=.true.
!         endif
!      enddo
!      read(io_LHEOutFile2,fmt=InputFmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
!      do nline=1,EventNumPart!  read event lines
!         read(io_LHEOutFile2,fmt="(A160)") EventLine(nline)
!      enddo
!      do nline=1,EventNumPart    
!         read(EventLine(nline),fmt=InputFmt1) LHE_IDUP(nline),IntExt(nline),LHE_MOTHUP(1,nline),LHE_MOTHUP(2,nline),LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline),Spin(nline),Lifetime(nline)
!      enddo
!      if( IsAQuark(LHE_IDUP(1)) .and. IsAQuark(LHE_IDUP(2)) ) then!  qq initial state
!      endif
!      read(io_LHEOutFile2,fmt="(A160)",IOSTAT=stat,END=99) HeaderLines! read </event> 
!      
!      
!      
!      
! 
! return
! END SUBROUTINE
SUBROUTINE OpenFiles()
use ModParameters
implicit none
logical :: dirresult


   call system('mkdir -p ./data')! -p is suppressing error messages if directory already exists

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

   
return
END SUBROUTINE




SUBROUTINE CloseFiles()
use ModParameters
implicit none

   close(io_LHEOutFile)
   close(io_HistoFile)
   if( ReadLHEFile ) close(io_LHEInFile)
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

  if( Process.eq.60 .or. Process.eq.61 ) then
     call InitHisto_HVBF()
  elseif( Process.eq.62) then
     call InitHisto_HJ()
  elseif (Process.eq.50) then
     call InitHisto_VHiggs()
  elseif (Process.eq.80) then
     call InitHisto_TTBH()
  elseif (Process.eq.90) then
     call InitHisto_BBBH()
  else
     call InitHisto_HZZ()
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
          Histo(11)%BinSize= 10d0*GeV
          Histo(11)%LowVal = 0d0*GeV
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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




SUBROUTINE InitHisto_TTBH()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

          it_sav = 1
          NumHistograms = 3
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_top"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "pT_H"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "D_0minus"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.02
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 1d0


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


SUBROUTINE InitHisto_HVBF()
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

!           Histo(3)%Info   = "log(MEsq)"
!           Histo(3)%NBins  = 11
!           Histo(3)%BinSize= 1d0
!           Histo(3)%LowVal = -10d0
!           Histo(3)%SetScale= 1d0
! 
!           Histo(4)%Info   = "log(MEsq*pdf)"
!           Histo(4)%NBins  = 11
!           Histo(4)%BinSize= 1d0
!           Histo(4)%LowVal = -10d0
!           Histo(4)%SetScale= 1d0

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












SUBROUTINE InitOutput
use ModParameters
implicit none

    if( (unweighted) .or. ( (.not.unweighted) .and. (writeWeightedLHE) )  ) then 
        write(io_LHEOutFile ,'(A)') '<LesHouchesEvents version="1.0">'
        write(io_LHEOutFile ,'(A)') '<!--'
        write(io_LHEOutFile ,'(A,A6,A)') 'Output from the JHUGenerator ',trim(JHUGen_Version),' described in arXiv:1001.3396 [hep-ph],arXiv:1208.4018 [hep-ph],arXiv:1309.4819 [hep-ph]'

        if( .not. ReadLHEFile ) then
            write(io_LHEOutFile ,'(A)') ''
            write(io_LHEOutFile ,'(A,F5.1,A)') 'hmass   ',M_Reso*100d0,'     ! Higgs boson mass'
            write(io_LHEOutFile ,'(A,F8.5,A)') 'hwidth ',Ga_Reso*100d0,'   ! Higgs boson width'
            write(io_LHEOutFile ,'(A)') ''
        endif

        call WriteParameters(io_LHEOutFile)

        if( (ReadLHEFile .or. ConvertLHEFile) .and. (importExternal_LHEinit) ) then
            write(io_LHEOutFile ,'(A)') ''
        else
            write(io_LHEOutFile ,'(A)') '-->'
            write(io_LHEOutFile ,'(A)') '<init>'
            if( (.not.unweighted) .and. (writeWeightedLHE) ) then
              if(Collider.eq.0)then
                write(io_LHEOutFile ,'(A,2F24.16,A)') ' -11   11',(Collider_Energy*50d0),(Collider_Energy*50d0),' 0 0     0     0 4  1' 
              else
                write(io_LHEOutFile ,'(A,2F24.16,A)') '2212 2212',(Collider_Energy*50d0),(Collider_Energy*50d0),' 0 0 10042 10042 4  1' 
              endif
            else
              if(Collider.eq.0)then
                write(io_LHEOutFile ,'(A,2F24.16,A)') ' -11   11',(Collider_Energy*50d0),(Collider_Energy*50d0),' 0 0     0     0 3  1' 
              else
                write(io_LHEOutFile ,'(A,2F24.16,A)') '2212 2212',(Collider_Energy*50d0),(Collider_Energy*50d0),' 0 0 10042 10042 3  1' 
              endif
            endif
! in order of appearance:  (see also http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017)
! (*) incoming particle1 (2212=proton), incoming particle2, 
! (*) energies of colliding particles, 
! (*) out-dated pdf information for colliding particles (supposed to be 0), 
! (*) pdf code of LHAGLUE for colliding particles (10042=CTEQ6Ll, MSTW2008=21000,21041-21080)    
! (*) weighting strategy (3=unweighted events, 4=weighted events,  otherwise=see LHE manuals)
! (*) number of process types to be accepted (default=1, otherwise=see manual)
! 
            write(io_LHEOutFile ,'(A)') '0.43538820803E-02  0.72559367904E-05  1.00000000000E-00 100'
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
    if( Process.eq.0 ) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.1 ) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=1, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.2 ) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=2, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.60) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.61) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.50) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.80) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.90) write(TheUnit,"(4X,A,F7.2,A,F7.5)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( ReadLHEFile )    write(TheUnit,"(4X,A)") "           (This is ReadLHEFile mode. Resonance mass/width might be overwritten by LHE input parameters. See below.)"
    if( ConvertLHEFile ) write(TheUnit,"(4X,A)") "           (This is ConvertLHEFile mode. Resonance mass/width might be overwritten by LHE input parameters. See below.)"
    write(TheUnit,"(4X,A,I2,2X,A,I2)") "DecayMode1:",DecayMode1, "DecayMode2:",DecayMode2
    if( IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z-boson: mass=",M_Z*100d0,", width=",Ga_Z*100d0
    if( IsAWDecay(DecayMode1) .or. IsAWDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W-boson: mass=",M_W*100d0,", width=",Ga_W*100d0
    if( Process.eq.80 ) write(TheUnit,"(4X,A,F8.4,A,F6.4)") "Top quark mass=",m_top*100d0,", width=",Ga_top*100d0
!     if( Process.eq.80 ) write(TheUnit,"(4X,A,I2)") "Top quark decay=",TOPDECAYS
    if( Process.eq.90 ) write(TheUnit,"(4X,A,F8.4,A,F6.4)") "Bottom quark mass=",m_top*100d0
    if( (ReadLHEFile) .and. (RequestNLeptons.gt.0) .and. (RequestOSSF) ) write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," OSSF leptons."
    if( (ReadLHEFile) .and. (RequestNLeptons.gt.0) .and. .not. (RequestOSSF)) write(TheUnit,"(4X,A,I2,A)") "Lepton filter activated. Requesting ",RequestNLeptons," leptons."
    write(TheUnit,"(4X,A,20I11)") "Random seeds: ",TheSeeds(1:TheSeeds(0))

    if( .not. (ReadLHEFile .or. ConvertLHEFile) ) then
        write(TheUnit,"(4X,A)") ""
        write(TheUnit,"(4X,A,L,L,L)") "OffXVV: ",OffShellReson,OffShellV1,OffShellV2
        write(TheUnit,"(4X,A,I1)") "PChannel: ",PChannel
#if useLHAPDF==1
        write(TheUnit,"(4X,A,A,A,I3)") "LHAPDF set ",trim(LHAPDFString), ", member=",LHAPDFMember
#else
        write(TheUnit,"(4X,A,I1)") "PDFSet: ",PDFSet
#endif
        write(TheUnit,"(4X,A,L)") "Unweighted: ",Unweighted
    endif
    write(TheUnit,"(4X,A,L)") "Interference: ",includeInterference
    if( (IsAZDecay(DecayMode1)) .or. (IsAZDecay(DecayMode2))  ) write(TheUnit,"(4X,A,L)") "Intermediate off-shell photons: ",includeGammaStar

    if( .not. seed_random ) write(TheUnit,"(4X,A)") "NOTE: seed_random==FALSE (switched off)"

    write(TheUnit,"(4X,A)") ""
    if( Process.eq.0 .or. Process.eq.60 .or. Process.eq.61 .or. Process.eq.62.or. Process.eq.50 ) then
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
            if( (includeGammaStar .and. ((IsAZDecay(DecayMode1)) .or. (IsAZDecay(DecayMode2)))) ) then
                write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs2=",ghzgs2,"i"
                write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs3=",ghzgs3,"i"
                write(TheUnit,"(6X,A,2E16.8,A1)") "ghzgs4=",ghzgs4,"i"
                write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs2=",ghgsgs2,"i"
                write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs3=",ghgsgs3,"i"
                write(TheUnit,"(6X,A,2E16.8,A1)") "ghgsgs4=",ghgsgs4,"i"
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
              Error  = 1d0/dsqrt(dble(Hits))
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



! if seed_random=true:  "Seeds" returns the seeds chosen according to system clock
! is seed_random=false: "Seeds" is used as input for initializing the random number generator
SUBROUTINE InitRandomSeed(Seeds)
use modParameters
use modMisc
implicit none
integer :: Seeds(0:20)
integer, dimension(:), allocatable :: gen_seed
integer :: n,i,sclock,SeedSize


    if (seed_random) then 
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
        call random_seed(size=SeedSize)
        Seeds(0) = SeedSize
        call random_seed(get=Seeds(1:SeedSize))
    else        
        call random_seed(size=n)
        if( n.ne.Seeds(0) ) call Error("Number of input seeds does not match random_seed(size=n)",n)
        call random_seed(put = Seeds(1:n))
    endif

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




SUBROUTINE PrintCommandLineArgs()
use modParameters
implicit none

        write(io_stdout,*) ""
        write(io_stdout,"(2X,A)") "Command line arguments:"
        write(io_stdout,"(4X,A)") "Collider:   1=LHC, 2=Tevatron, 0=e+e-"
        write(io_stdout,"(4X,A)") "Process: 0=spin-0, 1=spin-1, 2=spin-2 resonance, 50=pp/ee->VH, 60=weakVBF, 61=pp->Hjj, 62=pp->Hj, 80=pp->ttbar+H, 90=pp->bbbar+H"
        write(io_stdout,"(4X,A)") "MReso:      resonance mass (default=125.60), format: yyy.xx"
        write(io_stdout,"(4X,A)") "DecayMode1: decay mode for vector boson 1 (Z/W+/gamma)"
        write(io_stdout,"(4X,A)") "DecayMode2: decay mode for vector boson 2 (Z/W-/gamma)"
        write(io_stdout,"(4X,A)") "              0=Z->2l,  1=Z->2q, 2=Z->2tau, 3=Z->3nu,"
        write(io_stdout,"(4X,A)") "              4=W->lnu, 5=W->2q, 6=W->taunu,"
        write(io_stdout,"(4X,A)") "              7=gamma, 8=Z->2l+2tau,"
        write(io_stdout,"(4X,A)") "              9=Z->anything, 10=W->lnu+taunu, 11=W->anything"
        write(io_stdout,"(4X,A)") "TopDK:      decay mode for tops in ttbar+H, 0=stable, 1=decaying (use DecayMode1/2 = 4,5 for W+/W-"
        write(io_stdout,"(4X,A)") "PChannel:   0=g+g, 1=q+qb, 2=both"
        write(io_stdout,"(4X,A)") "OffXVV:     off-shell option for resonance(X),or vector bosons(VV)"
        write(io_stdout,"(4X,A)") "PDFSet:     1=CTEQ6L1(default), 2=MSTW2008LO,  2xx=MSTW with eigenvector set xx=01..40), 3=NNPDF3.0LO"
#if useLHAPDF==1
        write(io_stdout,"(4X,A)") "LHAPDF:     name of the LHA PDF file, e.g. NNPDF30_lo_as_0130/NNPDF30_lo_as_0130.info"
        write(io_stdout,"(4X,A)") "LHAPDFMem:  member PDF number, default=0 (best fit)"
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
    write(TheUnit,"(A90)") " *          I. Anderson, S. Bolognesi, F. Caola, Y. Gao, A. Gritsan, C. Martin,        *"
    write(TheUnit,"(A90)") " *                Z. Guo, K. Melnikov, H. Roskes, U. Sarica, M. Schulze,               *" 
    write(TheUnit,"(A90)") " *                   N. Tran, A. Whitbeck, M. Xiao, C. You, Y. Zhou                    *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D81 (2010) 075022;  arXiv:1001.3396 [hep-ph],              *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D86 (2012) 095031;  arXiv:1208.4018 [hep-ph],              *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D89 (2014) 035007;  arXiv:1309.4819 [hep-ph].              *"
    write(TheUnit,"(A90)") " *                                                                                     *"
    write(TheUnit,"(A90)") " ***************************************************************************************"
    write(TheUnit,"(A90)") " "
return
END SUBROUTINE



