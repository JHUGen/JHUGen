! If you use this program please cite Phys.Rev. D81 (2010) 075022; arXiv:1001.3396 [hep-ph].
!                                 and Phys.Rev. D86 (2012) 095031; arXiv:1208.4018 [hep-ph]
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
   call InitPDFs()
   call InitHisto()
   call InitParameters()
   call InitProcess()
   call InitVegas()
   call OpenFiles()
   call PrintLogo(io_stdout)
   call PrintLogo(io_LogFile)
   call WriteParameters(io_stdout)
   call WriteParameters(io_LogFile)
   call InitOutput()
   write(io_stdout,*) " Running"
   if( .not.useBetaVersion .and. .not.ReadLHEFile ) call StartVegas(VG_Result,VG_Error)
   if( useBetaVersion      .and. .not.ReadLHEFile ) call StartVegas_BETA(VG_Result,VG_Error)
   if( .not.useBetaVersion .and.      ReadLHEFile ) call StartReadLHE(VG_Result,VG_Error)
   call WriteHisto(VG_Result,VG_Error,time_end-time_start)
   call FinalizeOutput()
   call CloseFiles()
   write(io_stdout,*) " Done"


END PROGRAM





SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
implicit none
character :: arg*(120)
integer :: NumArgs,NArg,OffShell_XVV,iunwgt,CountArg

   Collider=1
   VegasIt1=-1
   PDFSet=1      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40
   VegasNc0=-1
   VegasNc1=-1
   VegasNc2=-1
   PChannel=2
   DecayMode1=0  ! Z/W+
   DecayMode2=0  ! Z/W-
   Process = 0   ! select 0, 1 or 2 to represent the spin of the resonance
   Unweighted =.true.
   OffShell_XVV=000! 000: X,V1,V2 on-shell; 010: X,V2 on-shell, V1 off-shell; and so on
   LHEProdFile=""
   ReadLHEFile=.false.

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
! !       DecayMode=11: W --> anyting

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
    elseif( arg(1:7) .eq."OffXVV=" ) then
        read(arg(8:10),*) OffShell_XVV
        CountArg = CountArg + 1
    elseif( arg(1:11) .eq."Unweighted=" ) then
        read(arg(12:12),*) iunwgt
        CountArg = CountArg + 1
        if( iunwgt.eq.0 ) then
            Unweighted = .false.
        else
            Unweighted = .true.
        endif
    elseif( arg(1:8) .eq."ReadLHE=" ) then
        read(arg(9:108),"(A)") LHEProdFile
        ReadLHEFile=.true.
        CountArg = CountArg + 1
    endif
   enddo


    Mu_Fact = M_Reso! setting pdf scale to resonance mass

   if( CountArg.ne.NumArgs ) then
        print *, "unknown command line argument"
        stop
   endif

    if (Process.eq.0) PChannel = 0   !only gluons
    if (Process.eq.1) PChannel = 1   !only quarks

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

!     if( ((OffShellV1).or.(OffShellV2).or.(OffShellReson)) ) then
!         print *, "off shell Z/W's only allowed for spin 0,2 resonance"
! !         stop
!     endif
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

    if( IsAZDecay(DecayMode1) .and. IsAWDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAWDecay(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAWDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
       print *, " DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( IsAPhoton(DecayMode1) .and. .not.IsAPhoton(DecayMode2) ) then
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

    call SetCtq6(4)  ! 4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl

return
END SUBROUTINE




SUBROUTINE InitParameters
use ModParameters
implicit none

IF( COLLIDER.EQ.1 ) THEN
  Collider_Energy  = LHC_Energy
ELSEIF( COLLIDER.EQ.2 ) THEN
  Collider_Energy  = TEV_Energy
ENDIF

END SUBROUTINE



SUBROUTINE InitProcess()
use ModParameters
use ModMisc
use ModCrossSection
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


      if( unweighted ) then
          NDim = NDim + 1  ! random number which decides if event is accepted
      endif
      NDim = NDim + 1  ! MC sampling for gg and qqb channel

END SUBROUTINE








SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22)
real(8) :: dum, RES(-5:5,-5:5)
logical :: warmup
integer :: i, i1, PChannel_aux, PChannel_aux1
include 'csmaxvalue.f'
integer :: n,clock
integer, dimension(:), allocatable :: gfort_seed


if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1) then 
      VegasNc1 = VegasNc1_default
      VegasNc2 = VegasNc2_default
endif


   warmup = .false.
   itmx = VegasIt1
   ncall= VegasNc1

   PChannel_aux = PChannel

if (unweighted.eqv..false.) then  !----------------------- weighted events

  call vegas(EvalWeighted,VG_Result,VG_Error,VG_Chi2)    ! usual call of vegas for weighted events



elseif(unweighted.eqv..true.) then  !----------------------- unweighted events

    VG = zero
    csmax = zero
    if (seed_random) then 
#if compiler==1
        call random_seed()
#elif compiler==2
        call random_seed(size=n)
        allocate(gfort_seed(n))
        call system_clock(count=clock)
        gfort_seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = gfort_seed)
        deallocate(gfort_seed)        
#endif
    endif


    print *, " finding maximal weight with ",VegasNc0," evaluations"
    do i=1,VegasNc0
        call random_number(yRnd)
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

        PChannel = PChannel_aux
    enddo



    csmax   = 1.5d0*csmax    !  adjustment factors, can be choosen  separately channel/by/channel
!      do i1=-5,5
!          print *, i1, csmax(i1,-i1), VG(i1,-i1)
!      enddo
!      pause

!------------------adj_par fixes by how much the quark-induced channels need to be adjusted

    if (PChannel.eq.2.and.fix_channels_ratio) then
        adj_par = VG(0,0)/(VG(-5,5)+VG(-4,4)+VG(-3,3)+VG(-2,2)+VG(-1,1)  &
                + VG(1,-1)+VG(2,-2)+VG(3,-3)+VG(4,-4)+VG(5,-5))*channels_ratio_fix/(one-channels_ratio_fix)
    else
        adj_par = 100000
    endif

!--- rescale the gluon induced channel
    csmax(0,0) = csmax(0,0)/adj_par


!------------------------------- set counts to zero for actual evaluation

    EvalCounter = 0
    AccepCounter = 0
    RejeCounter = 0
    AlertCounter = 0
    AccepCounter_part = 0

    if (seed_random) then 
#if compiler==1
        call random_seed()
#elif compiler==2
        call random_seed(size=n)
        allocate(gfort_seed(n))
        call system_clock(count=clock)
        gfort_seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = gfort_seed)
        deallocate(gfort_seed)        
#endif
    endif

    call cpu_time(time_start)
    if( VegasNc1.ne.-1 ) then
        print *, " generating events with ",VegasNc1," tries"
        do i=1,VegasNc1
            call random_number(yRnd)
            dum = EvalUnWeighted(yRnd,.true.,RES)! RES is a dummy here
        enddo

    elseif( VegasNc2.ne.-1 ) then
        print *, " generating ",VegasNc2," events"
        do while( AccepCounter.lt.VegasNc2 )
              call random_number(yRnd)
              dum = EvalUnWeighted(yRnd,.true.,RES)! RES is a dummy here
        enddo

    else
        print *, "ERROR: VegasNc1 and VegasNc2 must not be set at the same time"
        stop
    endif

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

  endif! unweighted

  write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)

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
logical :: warmup
integer :: i, i1, PChannel_aux, PChannel_aux1,n1,n2,n3,n4
include 'csmaxvalue.f'
integer :: n,clock
integer, dimension(:), allocatable :: gfort_seed

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

    if (seed_random) then 
#if compiler==1
        call random_seed()
#elif compiler==2
        call random_seed(size=n)
        allocate(gfort_seed(n))
        call system_clock(count=clock)
        gfort_seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = gfort_seed)
        deallocate(gfort_seed)        
#endif
    endif

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

endif

   write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)


return
END SUBROUTINE






SUBROUTINE StartReadLHE(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
use ModMisc
implicit none
include 'csmaxvalue.f'
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22),Res,dum,EMcheck(1:4)
real(8) :: AcceptedEvent(1:4,1:4),Ehat
real(8) :: MomExt(1:4,1:4),MomHiggs(1:4),MomGlu(1:4,1:3),Mass(1:4),pH2sq
integer :: tries, nParticle, MY_IDUP(1:10), ICOLUP(1:2,1:10)
character(len=*),parameter :: POWHEG_Fmt0 = "(6X,I2,A160)"
character(len=*),parameter :: POWHEG_Fmt1 = "(6X,I2,4X,I3,4X,I3,3X,I3,1X,I3,3X,I3,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9,1X,1PE16.9)"
logical :: FirstEvent,M_ResoSet
integer :: nline,LHE_IDUP(1:4),LHE_ICOLUP(1:2,1:4),LHE_ICOLUP_Glu(1:2,1:3),LHE_IDUP_Glu(1:3),intDummy,Nevent
integer :: EventNumPart,nglu
character(len=160) :: FirstLines,EventInfoLine
character(len=160) :: EventLine(1:4)
integer :: n,clock,i,stat
integer, dimension(:), allocatable :: gfort_seed


if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc0.eq.-1 ) VegasNc0 = VegasNc0_default
if( VegasNc1.eq.-1 .and. VegasNc2.eq.-1 ) VegasNc1 = VegasNc1_default
if( VegasNc1.eq.-1 .and. .not.VegasNc2.eq.-1 ) VegasNc1 = VegasNc2

    if (seed_random) then 
#if compiler==1
        call random_seed()
#elif compiler==2
        call random_seed(size=n)
        allocate(gfort_seed(n))
        call system_clock(count=clock)
        gfort_seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = gfort_seed)
        deallocate(gfort_seed)        
#endif
    endif



!    search for line with first event
     FirstEvent = .false.
     M_ResoSet=.false.
     do while ( .not.FirstEvent )
        read(16,fmt="(A160)",IOSTAT=stat,END=99) FirstLines
        if( FirstLines(1:5).eq."hmass" ) then 
               read(FirstLines(6:13),fmt="(F7.0)") M_Reso
               M_Reso = M_Reso*GeV!  convert to units of 100GeV
               M_ResoSet=.true.
        endif
        if( FirstLines(1:7).eq."<event>" ) then 
               FirstEvent=.true.
        else
            if( importPOWHEG_LHEinit ) then 
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
     write(io_stdout,"(A)") ""
     write(io_LogFile,"(A)") ""




      print *, " finding maximal weight with ",VegasNc0," points"
      VG = zero
      CSmax = zero
      Ehat = M_Reso! fixing Ehat to M_Reso which should determine the max. of the integrand
      do tries=1,VegasNc0
          call random_number(yRnd)
          dum = EvalUnWeighted_withoutProduction(yRnd,.false.,EHat,Res,AcceptedEvent,MY_IDUP(1:9),ICOLUP(1:2,1:9))
      enddo
      csmax(0,0)   = 1.5d0*csmax(0,0)    !  savety factor





     print *, " generating events"
     EvalCounter = 0
     AccepCounter = 0
     RejeCounter = 0
     AccepCounter_part = 0
     call cpu_time(time_start)
     NEvent=0
     do while ( .true. ) 
         NEvent=NEvent + 1

         read(16,fmt=POWHEG_Fmt0) EventNumPart,EventInfoLine!  read number of particle from the first line after <event> and other info
!        read event lines
         do nline=1,EventNumPart
            read(16,fmt="(A160)") EventLine(nline)
         enddo
         if( EventNumPart.ne.3 .and. EventNumPart.ne.4 ) then
            print *, "ERROR: number of particles in LHE input file is neither 3 nor 4",EventNumPart
            stop
         endif

 !       convert event lines into variables assuming that Higgs is always in 3rd place
         nglu = 0 
         do nline=1,EventNumPart
            read(EventLine(nline),fmt=POWHEG_Fmt1) LHE_IDUP(nline),intDummy,intDummy,intDummy,LHE_ICOLUP(1,nline),LHE_ICOLUP(2,nline),MomExt(2,nline),MomExt(3,nline),MomExt(4,nline),MomExt(1,nline),Mass(nline)
            MomExt(1:4,nline) = MomExt(1:4,nline)*GeV!  convert to units of 100GeV
            Mass(nline) = Mass(nline)*GeV            !  convert to units of 100GeV

            if( abs(LHE_IDUP(nline)).eq.25 ) then!   select the Higgs (ID=25, h0)
                  MomHiggs(1:4) = MomExt(1:4,nline)
                  pH2sq = dsqrt(abs(MomHiggs(1:4).dot.MomHiggs(1:4)))
            elseif( abs(LHE_IDUP(nline)).eq.21 .or. (abs(LHE_IDUP(nline)).ge.1 .and. abs(LHE_IDUP(nline)).le.5) ) then!   select the gluons (ID=21) or quarks (ID=1,.,5)
                  nglu = nglu + 1
                  MomGlu(1:4,nglu) = MomExt(1:4,nline)
                  LHE_IDUP_Glu(nglu) = LHE_IDUP(nline)
                  LHE_ICOLUP_Glu(1:2,nglu) = LHE_ICOLUP(1:2,nline)
             else
                  print *, "ERROR: unknown particle in LHE input file",LHE_IDUP(nline)
                  stop
            endif
         enddo

! !        reject event if M_Reso and pH2sq deviate by more than 20 GeV 
!          if( abs( M_Reso - pH2sq) .gt. 20d0*GeV ) then
!               write(io_stdout,"(2X,A,2F10.4)")  "WARNING: Higgs mass and momentum squared deviate by more than 20 GeV!",pH2sq*100d0,M_Reso*100d0
!               write(io_LogFile,"(2X,A,2F10.4)") "WARNING: Higgs mass and momentum squared deviate by more than 20 GeV!",pH2sq*100d0,M_Reso*100d0
!               cycle
!          endif

         EMcheck(1:4) = MomGlu(1:4,1) + MomGlu(1:4,2) - MomHiggs(1:4)
         if( nglu.eq.3 ) EMcheck(1:4) = EMcheck(1:4) - MomGlu(1:4,3)
         if( any(abs(EMcheck(1:4)).gt.1d0*GeV) ) then
              print *, "ERROR: energy momentum violation while reading LHE production momenta.",EMcheck(1:4)
              stop
         endif


!         accept/reject sampling for H->VV decay contribution
          EHat = pH2sq
          do tries=1,5000000
              call random_number(yRnd)
              dum = EvalUnWeighted_withoutProduction(yRnd,.true.,Ehat,RES,AcceptedEvent,MY_IDUP(1:9),ICOLUP(1:2,1:9))
              if( Res.ne.0d0 ) exit
          enddo
          if( Res.ne.0d0 ) then ! decay event was accepted
             call boost(AcceptedEvent(1:4,1),MomHiggs(1:4),pH2sq)
             call boost(AcceptedEvent(1:4,2),MomHiggs(1:4),pH2sq)
             call boost(AcceptedEvent(1:4,3),MomHiggs(1:4),pH2sq)
             call boost(AcceptedEvent(1:4,4),MomHiggs(1:4),pH2sq)
             if(EventNumPart.eq.3) then
                    MY_IDUP(1) = convertLHEreverse(LHE_IDUP_Glu(1))
                    MY_IDUP(2) = convertLHEreverse(LHE_IDUP_Glu(2))
                    ICOLUP(1:2,1:2) = LHE_ICOLUP_Glu(1:2,1:2)!  overwrite inital colors with POWHEG ones
                    call WriteOutEvent((/MomGlu(1:4,1),MomGlu(1:4,2),AcceptedEvent(1:4,1),AcceptedEvent(1:4,2),AcceptedEvent(1:4,3),AcceptedEvent(1:4,4)/),MY_IDUP(1:9),ICOLUP(1:2,1:9),EventInfoLine=EventInfoLine)
             elseif(EventNumPart.eq.4) then
                    MY_IDUP(1) = convertLHEreverse(LHE_IDUP_Glu(1))!  overwrite inital ID's with POWHEG ones
                    MY_IDUP(2) = convertLHEreverse(LHE_IDUP_Glu(2))
                    MY_IDUP(10)= convertLHEreverse(LHE_IDUP_Glu(3))
                    ICOLUP(1:2,1:2) = LHE_ICOLUP_Glu(1:2,1:2)!  overwrite inital colors with POWHEG ones
                    ICOLUP(1:2,10)  = LHE_ICOLUP_Glu(1:2,3)  !  set third parton color
                    call WriteOutEvent((/MomGlu(1:4,1),MomGlu(1:4,2),AcceptedEvent(1:4,1),AcceptedEvent(1:4,2),AcceptedEvent(1:4,3),AcceptedEvent(1:4,4)/),MY_IDUP(1:10),ICOLUP(1:2,1:10),MomRealGlu=MomGlu(1:4,3),EventInfoLine=EventInfoLine)
             endif
             if( mod(AccepCounter,5000).eq.0 ) then
                  call cpu_time(time_int)
                  write(io_stdout,*)  NEvent," events accepted (",time_int-time_start, ") seconds"
                  write(io_LogFile,*) NEvent," events accepted (",time_int-time_start, ") seconds"
             endif
          else! decay event was not accepted after ncall evaluations, read next production event
             print *, "rejected event after ",tries-1," evaluations"
             AlertCounter = AlertCounter + 1 
          endif


!        skip event lines 
         read(16,fmt="(A7)",IOSTAT=stat,END=99) FirstLines! skip <\event>
!          if( stat.lt.0 ) exit
         read(16,fmt="(A30)",IOSTAT=stat,END=99) FirstLines!   skip <event> or </LesHouchesEvents>
         if( FirstLines(1:30).eq."</LesHouchesEvents>" ) exit
         if( NEvent.eq. VegasNc1 ) exit

     enddo
99   continue
     call cpu_time(time_end)




    write(io_stdout,*) ""
    write(io_stdout,*) "Evaluation Counter: ",EvalCounter
    write(io_stdout,*) "Acceptance Counter: ",AccepCounter
    write(io_stdout,*) "Rejection  Counter: ",RejeCounter
    write(io_stdout,*) " Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        write(io_stdout,*) "ALERT: The number of rejected events exceeds 1%."
        write(io_stdout,*) "       Increase CSMAX in main.F90 or VegasNc1."
    endif
   write(io_stdout,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)

    write(io_LogFile,*) ""
    write(io_LogFile,*) "Evaluation Counter: ",EvalCounter
    write(io_LogFile,*) "Acceptance Counter: ",AccepCounter
    write(io_LogFile,*) "Rejection  Counter: ",RejeCounter
    write(io_LogFile,*) " Alert  Counter: ",AlertCounter
    if( dble(AlertCounter)/dble(AccepCounter) .gt. 1d0*percent ) then
        write(io_LogFile,*) "ALERT: The number of rejected events exceeds 1%."
        write(io_LogFile,*) "       Increase CSMAX in main.F90 or VegasNc1."
    endif
   write(io_LogFile,*)  " event generation rate (events/sec)",dble(AccepCounter)/(time_end-time_start)



return
END SUBROUTINE




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
      open(unit=io_HistoFile,file=trim(DataFile)//'.dat',form='formatted',access= 'sequential',status='replace')         ! histogram file
   endif
   open(unit=io_LogFile,file=trim(DataFile)//'.log',form='formatted',access= 'sequential',status='replace')              ! log file


   if( ReadLHEFile ) then 
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

          Histo(12)%Info   = "total ee cross section"
          Histo(12)%NBins  = 1
          Histo(12)%BinSize= 1d0
          Histo(12)%LowVal =  0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "total mm cross section"
          Histo(13)%NBins  = 1
          Histo(13)%BinSize= 1d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "total tt cross section"
          Histo(14)%NBins  = 1
          Histo(14)%BinSize= 1d0
          Histo(14)%LowVal = 0d010
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "total em cross section"
          Histo(15)%NBins  = 1
          Histo(15)%BinSize= 1d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0

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





SUBROUTINE InitOutput
use ModParameters
implicit none

    if( unweighted ) then 
        write(io_LHEOutFile ,'(A)') '<LesHouchesEvents version="1.0">'
        write(io_LHEOutFile ,'(A)') '<!--'
        write(io_LHEOutFile ,'(A,A6,A)') 'Output from the JHUGenerator ',trim(JHUGen_Version),' described in arXiv:1001.3396 [hep-ph],arXiv:1208.4018 [hep-ph]'

        write(io_LHEOutFile ,'(A)') ''

        call WriteParameters(io_LHEOutFile)

        if( ReadLHEFile .and. importPOWHEG_LHEinit ) then
            write(io_LHEOutFile ,'(A)') '------------------------------------------------------'
        else
            write(io_LHEOutFile ,'(A)') '-->'
            write(io_LHEOutFile ,'(A)') '<init>'
            write(io_LHEOutFile ,'(A,2F24.16,A)') '2212 2212',(Collider_Energy*50d0),(Collider_Energy*50d0),' 0 0 10042 10042 3  1' 
! in order of appearance:  (see also http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017)
! (*) incoming particle1 (2212=proton), incoming particle2, 
! (*) energies of colliding particles, 
! (*) out-dated pdf information for colliding particles (supposed to be 0), 
! (*) pdf code of LHAGLUE for colliding particles (10042=CTEQ6Ll, MSTW2008=21000,21041-21080)    
! (*) weighting strategy (3=accept all weights, otherwise=see LHE manuals)
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

    if( unweighted ) then 
        write(io_LHEOutFile ,'(A)') '</LesHouchesEvents>'
    endif
  
END SUBROUTINE




SUBROUTINE WriteParameters(TheUnit)
use ModParameters
implicit none
integer :: TheUnit


    write(TheUnit,"(3X,A)") "Input Parameter:"
    if( Collider.eq.1 ) write(TheUnit,"(4X,A,1F8.2)") "Collider: P-P, sqrt(s)=",Collider_Energy*100d0
    if( Collider.eq.2 ) write(TheUnit,"(4X,A,1F8.2)") "Collider: P-Pbar, sqrt(s)=",Collider_Energy*100d0
    if( Process.eq.0 ) write(TheUnit,"(4X,A,F7.2,A,F6.3)") "Resonance: spin=0, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.1 ) write(TheUnit,"(4X,A,F7.2,A,F6.3)") "Resonance: spin=1, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( Process.eq.2 ) write(TheUnit,"(4X,A,F7.2,A,F6.3)") "Resonance: spin=2, mass=",M_Reso*100d0," width=",Ga_Reso*100d0
    if( ReadLHEFile ) write(TheUnit,"(4X,A)") "           (This is ReadLHEFile mode. Resonance mass is read from LHE input file.)"
    write(TheUnit,"(4X,A,I2,2X,A,I2)") "DecayMode1:",DecayMode1, "DecayMode2:",DecayMode2
    if( IsAZDecay(DecayMode1) .or. IsAZDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "Z-boson: mass=",M_Z*100d0,", width=",Ga_Z*100d0
    if( IsAWDecay(DecayMode1) .or. IsAWDecay(DecayMode2) ) write(TheUnit,"(4X,A,F6.3,A,F6.4)") "W-boson: mass=",M_W*100d0,", width=",Ga_W*100d0


    if( .not.ReadLHEFile ) then
        write(TheUnit,"(4X,A)") ""
        write(TheUnit,"(4X,A,L,L,L)") "OffXVV: ",OffShellReson,OffShellV1,OffShellV2
        write(TheUnit,"(4X,A,I1)") "PChannel: ",PChannel
        write(TheUnit,"(4X,A,I1)") "PDFSet: ",PDFSet
        write(TheUnit,"(4X,A,L)") "Unweighted: ",Unweighted
    endif
    write(TheUnit,"(4X,A,L)") "Interference: ",includeInterference


    write(TheUnit,"(4X,A)") ""
    if( Process.eq.0 ) then
        write(TheUnit,"(4X,A)") "spin-0-VV couplings: "
        write(TheUnit,"(6X,A,L)") "generate_as=",generate_as
        if( generate_as ) then 
            write(TheUnit,"(6X,A,2F7.4,A1)") "ahg1=",ahg1,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ahg2=",ahg2,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ahg3=",ahg3,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ahz1=",ahz1,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ahz2=",ahz2,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ahz3=",ahz3,"i"
        else
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghg2=",ghg2,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghg3=",ghg3,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghg4=",ghg4,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghz1=",ghz1,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghz2=",ghz2,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghz3=",ghz3,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "ghz4=",ghz4,"i"
        endif
    elseif( Process.eq.1 ) then
        write(TheUnit,"(4X,A)") "spin-1-VV couplings: "
        write(TheUnit,"(6X,A,2F7.4,A1)") "zprime_qq_left =",zprime_qq_left,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "zprime_qq_right=",zprime_qq_right,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "zprime_zz_1=",zprime_zz_1,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "zprime_zz_2=",zprime_zz_2,"i"
    elseif( Process.eq.2 ) then
        write(TheUnit,"(4X,A)") "spin-2-VV couplings: "
        write(TheUnit,"(6X,A,L)") "generate_bis=",generate_bis
        write(TheUnit,"(6X,A,L)") "use_dynamic_MG=",use_dynamic_MG
        write(TheUnit,"(6X,A,2F7.4,A1)") "a1=",a1,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "a2=",a2,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "a3=",a3,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "a4=",a4,"i"
        write(TheUnit,"(6X,A,2F7.4,A1)") "a5=",a5,"i"
        if( generate_bis ) then 
            write(TheUnit,"(6X,A,2F7.4,A1)") "b1 =",b1,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b2= ",b2,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b3 =",b3,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b4 =",b4,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b5 =",b5,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b6 =",b6,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b7 =",b7,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b8 =",b8,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b9 =",b9,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "b10=",b10,"i"
        else
            write(TheUnit,"(6X,A,2F7.4,A1)") "c1 =",c1,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c2 =",c2,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c3 =",c3,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c41 =",c41,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c42 =",c42,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c5 =",c5,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c6 =",c6,"i"
            write(TheUnit,"(6X,A,2F7.4,A1)") "c7 =",c7,"i"
        endif
    endif
    write(TheUnit,"(6X,A,2F8.1,A1)") "Lambda=",Lambda*100d0


    write(TheUnit,"(4X,A)") ""
    write(TheUnit,"(4X,A)") "LHE output: "//trim(DataFile)//'.lhe'
    write(TheUnit,"(4X,A)") "Histogram output: "//trim(DataFile)//'.dat'
    write(TheUnit,"(4X,A)") "Log file: "//trim(DataFile)//'.log'
    if( ReadLHEFile ) write(TheUnit,"(4X,A)") "LHE input: "//trim(LHEProdFile)
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




SUBROUTINE PrintCommandLineArgs()
use modParameters
implicit none

        write(io_stdout,*) ""
        write(io_stdout,"(2X,A)") "Command line arguments:"
        write(io_stdout,"(4X,A)") "Collider:   1=LHC, 2=Tevatron"
        write(io_stdout,"(4X,A)") "Process:    0=spin-0, 1=spin-1, 2=spin-2 resonance"
        write(io_stdout,"(4X,A)") "MReso:      resonance mass (default=126.00), format: yyy.xx"
        write(io_stdout,"(4X,A)") "DecayMode1: decay mode for vector boson 1 (Z/W+/gamma)"
        write(io_stdout,"(4X,A)") "DecayMode2: decay mode for vector boson 2 (Z/W-/gamma)"
        write(io_stdout,"(4X,A)") "              0=Z->2l,  1=Z->2q, 2=Z->2tau, 3=Z->2nu,"
        write(io_stdout,"(4X,A)") "              4=W->lnu, 5=W->2q, 6=W->taunu,"
        write(io_stdout,"(4X,A)") "              7=gamma, 8=Z->2l+2tau,"
        write(io_stdout,"(4X,A)") "              9=Z->anything, 10=W->lnu+taunu, 11=W->anything"
        write(io_stdout,"(4X,A)") "PChannel:   0=g+g, 1=q+qb, 2=both"
        write(io_stdout,"(4X,A)") "OffXVV:     off-shell option for resonance(X),or vector bosons(VV)"
        write(io_stdout,"(4X,A)") "PDFSet:     1=CTEQ6L1(2001), 2=MSTW(2008),  2xx=MSTW with eigenvector set xx=01..40)"
        write(io_stdout,"(4X,A)") "VegasNc0:   number of evaluations for integrand scan"
        write(io_stdout,"(4X,A)") "VegasNc1:   number of evaluations for accept-reject sampling"
        write(io_stdout,"(4X,A)") "Unweighted: 0=weighted events, 1=unweighted events"
        write(io_stdout,"(4X,A)") "DataFile:   LHE output file"
        write(io_stdout,"(4X,A)") "ReadLHE:    LHE input file from POWHEG (only spin-0)"
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
    write(TheUnit,"(A90)") " *                      S. Bolognesi, Y. Gao, A. Gritsan, Z. Guo,                      *"
    write(TheUnit,"(A90)") " *                    K. Melnikov, M. Schulze, N. Tran, A. Whitbeck                    *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D81 (2010) 075022;  arXiv:1001.3396 [hep-ph]               *"
    write(TheUnit,"(A90)") " *                Phys.Rev. D86 (2012) 095031;  arXiv:1208.4018 [hep-ph]               *"
    write(TheUnit,"(A90)") " *                                                                                     *"
    write(TheUnit,"(A90)") " ***************************************************************************************"
    write(TheUnit,"(A90)") " "
return
END SUBROUTINE


