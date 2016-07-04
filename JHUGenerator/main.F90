PROGRAM JHUGenerator
use ModParameters
use ModKinematics
use ModCrossSection
#if compiler==1
use ifport
#endif
implicit none
real(8) :: time_start,time_end
real(8) :: VG_Result,VG_Error

   call GetCommandlineArgs()
   call InitPDFs()
   call InitHisto()
   call InitParameters()
   call InitProcess()
   call InitVegas()
   call OpenFiles()
   call InitOutput()
   print *, "Running"
   call cpu_time(time_start)
       call StartVegas(VG_Result,VG_Error)
   call cpu_time(time_end)
   call WriteHisto(VG_Result,VG_Error,time_end-time_start)
   call CloseFiles()
   print *, "Done"

END PROGRAM





SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
implicit none
character :: arg*(31)
integer :: NumArgs,NArg,OffShell_XVV

   Collider=1
   VegasIt1=-1
   PDFSet=1      ! 1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40
   VegasNc1=-1
   PChannel=2
   DecayMode1=0
   DecayMode2=0
   Process = 0   ! select 0, 1 or 2 to represent the spin of the resonance
   Unweighted =.true.
   OffShell_XVV=000! 000: X,V1,V2 on-shell; 010: X,V2 on-shell, V1 off-shell; and so on

! !       DecayMode=0: Z --> l+ l- (l=e,mu)
! !       DecayMode=1: Z --> q qbar (q=u,d,c,s,b)
! !       DecayMode=2: Z --> tau+ tau-
! !       DecayMode=3: Z --> nu nubar (nu=nu_e,nu_mu,nu_tau)
! !       DecayMode=4: W --> l nu_l (l=e,mu)
! !       DecayMode=5: W --> q qbar' (q=u,c, qbar'=d,s)
! !       DecayMode=6: W --> tau nu_tau
! !       DecayMode=7: photon

   DataFile="./data/output"

#if compiler==1
   NumArgs = NArgs()-1
#elif compiler==2
   NumArgs = COMMAND_ARGUMENT_COUNT()
#endif
   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    if( arg(1:9).eq."Collider=" ) then
        read(arg(10:10),*) Collider
    elseif( arg(1:7).eq."PDFSet=" ) then
        read(arg(8:10),*) PDFSet
    elseif( arg(1:9).eq."VegasNc1=" ) then
        read(arg(10:17),*) VegasNc1
    elseif( arg(1:9).eq."PChannel=" ) then
        read(arg(10:11),*) PChannel
    elseif( arg(1:9).eq."DataFile=" ) then
        read(arg(10:31),*) DataFile
    elseif( arg(1:8).eq."Process=" ) then
        read(arg(9:10),*) Process
    elseif( arg(1:11).eq."DecayMode1=" ) then
        read(arg(12:13),*) DecayMode1
    elseif( arg(1:11).eq."DecayMode2=" ) then
        read(arg(12:13),*) DecayMode2
    elseif( arg(1:7) .eq."OffXVV=" ) then
        read(arg(8:10),*) OffShell_XVV
    endif
   enddo

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

    if( (DecayMode1.le.3) .and. (DecayMode2.ge.4) ) then
       print *, "1 DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( (DecayMode1.ge.4) .and. (DecayMode2.le.3) ) then
       print *, "2 DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( (DecayMode1.eq.7) .and. (DecayMode2.ne.7) ) then
       print *, "3 DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
       stop
    endif

    if( (DecayMode1.eq.7) .and. (OffShellV1 .or. OffShellV2) ) then
       print *, "4 Photons have to be on-shell."
       stop
    endif

    if( (DecayMode1.ge.8) .or. (DecayMode2.ge.8) .or. (DecayMode1.lt..0) .or. (DecayMode2.lt.0) ) then
       print *, "5 DecayMode1=",DecayMode1," and DecayMode2=",DecayMode2," are not allowed."
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! added by Nhan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE InitOutput
use ModParameters
implicit none

        write(14 ,'(A)') '<LesHouchesEvents version="1.0">'
        write(14 ,'(A)') '<!--'
        write(14 ,'(A)') 'Output from the generator given in http://arxiv.org/abs/1001.3396'
        write(14 ,'(A)') '-->'
        write(14 ,'(A)') '<init>'
        write(14 ,'(A,2F24.16,A)') '2212 2212',(Collider_Energy*50d0),(Collider_Energy*50d0),' 0 0 10042 10042 3  1'
        write(14 ,'(A)') '0.43538820803E-02  0.72559367904E-05  0.87076000000E-07 100'
        write(14 ,'(A)') '</init>'

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
      NDim = NDim + 2    ! shat integration
      NDim = NDim + 1    ! offzchannel
      if( OffShellV1 .or. OffShellV2 ) then
          NDim = NDim + 2    ! integration over Z's invariant mass
      else
          NDim = NDim + 3    ! integration over mz1, mz2, mg for adjusted kinematics but with on-shell matrix elements
      endif
      VegasIt1_default = 1
      VegasNc1_default = 5000000

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
integer :: n,clock
integer, dimension(:), allocatable :: gfort_seed
include 'csmaxvalue.f'

if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc1.eq.-1 ) VegasNc1 = VegasNc1_default

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
    print *, "finding maximal weight"
    do i=1,ncall
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

!    print *, PChannel
!    pause

    csmax   = 5d0*csmax    !  adjustment factors, can be choosen  separately channel/by/channel
!      do i1=-5,5
!          print *, i1, csmax(i1,-i1), VG(i1,-i1)
!      enddo
!      pause

!------------------adj_par fixes by how much the quark-induced channels need to be adjusted

    if (PChannel.eq.2.and.fix_channels_ratio) then
        adj_par = VG(0,0)/(VG(-5,5)+VG(-4,4)+VG(-3,3)+VG(-2,2)+VG(-1,1)  &
                + VG(1,-1)+VG(2,-2)+VG(3,-3)+VG(4,-4)+VG(5,-5))*channels_ratio_fix/(one-channels_ratio_fix)
    else
        adj_par = one
    endif

!--- rescale the gluon induced channel
    csmax(0,0) = csmax(0,0)/adj_par


!------------------------------- set counts to zero for actual evaluation

    EvalCounter = 0
    AccepCounter = 0
    RejeCounter = 0
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
        
    print *, "generating events"
    do i=1,ncall
        call random_number(yRnd)
        dum = EvalUnWeighted(yRnd,.true.,RES)! RES is a dummy here
    enddo

    print *, "Evaluation Counter: ",EvalCounter
    print *, "Acceptance  Counter: ",AccepCounter
    print *, "Rejection  Counter: ",RejeCounter
    do i1=-5,5
      print *, "Acceptance  Counter_part: ", i1, AccepCounter_part(i1,-i1)
    enddo

  endif! unweighted


return
END SUBROUTINE





SUBROUTINE OpenFiles()
use ModParameters
implicit none

   print *, ""
   print *, "Data file:         "//trim(DataFile)//'.lhe'
   print *, "Vegas status file: "//trim(DataFile)//'.status'
   open(unit=14,file=trim(DataFile)//'.lhe',form='formatted',access= 'sequential',status='replace')            ! LHE event file
   open(unit=15,file=trim(DataFile)//'.dat',form='formatted',access= 'sequential',status='replace')         ! histogram file

return
END SUBROUTINE




SUBROUTINE CloseFiles()
implicit none

   close(14)
   close(15)

return
END SUBROUTINE




SUBROUTINE InitHisto()
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

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0/GeV

          Histo(2)%Info   = "pT_LepM"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 1d0/GeV

          Histo(3)%Info   = "cos(Psi_LepPZ)"    ! angle between lepton momentum in the rest frame of its own Z(call it Z1) and
          Histo(3)%NBins  = 40                  ! the direction of the other Z (call it Z2) in the rest frame of Z1
          Histo(3)%BinSize= 0.05d0
          Histo(3)%LowVal =-1d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "Psi_LepPlanes"  ! angle between lepton planes
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.078539d0
          Histo(4)%LowVal =0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "cos(ThetaZ)"   ! scattering angle of Z in graviton rest frame
          Histo(5)%NBins  = 200
          Histo(5)%BinSize= 0.01d0
          Histo(5)%LowVal =-1d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "MV1 invariant mass"
          Histo(6)%NBins  = 200
          Histo(6)%BinSize= 0.5d0*GeV
          Histo(6)%LowVal = 10d0*GeV
          Histo(6)%SetScale= 1d0/GeV

          Histo(7)%Info   = "MV2 invariant mass"
          Histo(7)%NBins  = 200
          Histo(7)%BinSize= 0.5d0*GeV
          Histo(7)%LowVal = 10d0*GeV
          Histo(7)%SetScale= 1d0/GeV

          Histo(8)%Info   = "MG resonance invariant mass"
          Histo(8)%NBins  = 200
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 1d0/GeV


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





SUBROUTINE WriteHisto(VG_Result,VG_Error,RunTime)
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
integer :: NBin,Hits,NHisto
real(8) :: BinSize,LowVal,BinVal,Value,Error,Integral
real(8),parameter :: ToGeV=1d2, ToPb=1d-3
real(8) :: VG_Result,VG_Error,RunTime


  write(* ,"(A,2X,1F9.2,A)") "# run time =",RunTime/60d0,"min"
  write(15,"(A,2X,1F9.2,A)") "# run time =",RunTime/60d0,"min"
  write(* ,"(A,2X,1PE20.10,2X,1PE20.5)") "# TotCS =",VG_Result,VG_Error
  write(15,"(A,2X,1PE20.10,2X,1PE20.5)") "# TotCS =",VG_Result,VG_Error
  write(* ,"(A)") "#"
  write(15,"(A)") "#"
  do NHisto=1,NumHistograms
      write(* ,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      write(15,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
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
          write(* ,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
          write(15,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
      enddo
      Integral = Integral + (Histo(NHisto)%Value(0)+Histo(NHisto)%Value(Histo(NHisto)%NBins+1))/it
      write(*,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(15,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(*,"(A)") "#"
      write(15,"(A)") "#"
  enddo

  write(14, '(A,X,I9,X,A)') '<!-- Number of events:', AccepCounter,' -->'
  write(14 ,'(A)') '</LesHouchesEvents>'

return
END SUBROUTINE
