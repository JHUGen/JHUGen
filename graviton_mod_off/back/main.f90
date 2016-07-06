PROGRAM Graviton
use ModParameters
use ModKinematics
use ModCrossSection
use ifport
implicit none
real(8) :: time_start,time_end
real(8) :: VG_Result,VG_Error
include 'csmaxvalue.f'

    csmax_gg = 0d0
    csmax_qq = 0d0

   call GetCommandlineArgs()
   call InitPDFs()
   call InitHisto()
   call InitParameters()
   call InitProcess()
   call InitVegas()
   call OpenFiles()
   print *, "Running";
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
integer :: NumArgs,NArg

   Collider=1
   VegasIt1=-1
   VegasNc1=-1
   PChannel=2
   Process = 2   ! here, select 0, 1 or 2 to represent the spin of the resonance
   Unweighted = .false.
   DataFile="./data/test"

   NumArgs = NArgs()-1
   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    if( arg(1:9).eq."Collider=" ) then
        read(arg(10:10),*) Collider
!     elseif( arg(1:9).eq."VegasIt1=" ) then
!         read(arg(10:11),*) VegasIt1
    elseif( arg(1:9).eq."VegasNc1=" ) then
        read(arg(10:17),*) VegasNc1
    elseif( arg(1:9).eq."PChannel=" ) then
        read(arg(10:11),*) PChannel
    elseif( arg(1:9).eq."DataFile=" ) then
        read(arg(10:31),*) DataFile
    endif
   enddo

    if (Process.eq.0) PChannel = 0   !only gluons 
    if (Process.eq.1) PChannel = 1   !only quarks

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


      NDim = 0
      NDim = NDim + 6    ! PS integration
      NDim = NDim + 2    ! shat integration
      VegasIt1_default = 1
      VegasNc1_default = 1000000
!      VegasNc1_default = 30000000
      unweighted = .true.

      if( unweighted ) then
          NDim = NDim + 1  ! unweighting
      endif
      NDim = NDim + 1  ! MC sampling for gg and qqb channel

END SUBROUTINE




SUBROUTINE InitHisto()
use ModMisc
use ModKinematics
use ModParameters
implicit none
integer :: AllocStatus,NHisto

  it_sav = 1

          NumHistograms = 5
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0/100d0
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_LepM"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0/100d0
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

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

  do NHisto=1,NumHistograms
      if( .not.allocated(Histo(NHisto)%Value) ) then
        allocate( Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo(NHisto)%Value2) ) then
        allocate( Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo(NHisto)%Hits) ) then
        allocate( Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1) = 0d0
      Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1)= 0d0
      Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1)  = 0
  enddo


RETURN
END SUBROUTINE






SUBROUTINE StartVegas(VG_Result,VG_Error)
use ModCrossSection
use ModKinematics
use ModParameters
implicit none
include "vegas_common.f"
real(8) :: VG_Result,VG_Error,VG_Chi2
real(8) :: yRnd(1:22)
real(8) :: dum
logical :: warmup
integer :: i, PChannel_aux, PChannel_aux1
include 'csmaxvalue.f'

if( VegasIt1.eq.-1 ) VegasIt1 = VegasIt1_default
if( VegasNc1.eq.-1 ) VegasNc1 = VegasNc1_default

   warmup = .false.
   itmx = VegasIt1
   ncall= VegasNc1

   PChannel_aux = PChannel

   if (unweighted.eq..false.) then  !----------------------- for unweighted events  
       
  call vegas(EvalCS,VG_Result,VG_Error,VG_Chi2)

   elseif(unweighted.eq..true.) then 

      VG_gg = zero 
      VG_qq = zero 
      csmax_gg = zero 
      csmax_qq = zero

 if (PChannel_aux.eq.0.or.PChannel_aux.eq.2) then 
    PChannel= 0
    call vegas(EvalCS,VG_Result,VG_Error,VG_Chi2)
    VG_gg = VG_Result
    csmax_gg = 5d0*csmax_gg   ! adjustment factor 
 endif    

if(PChannel_aux.eq.1.or.PChannel_aux.eq.2) then 
    PChannel = 1 
    call vegas(EvalCS,VG_Result,VG_Error,VG_Chi2)
    VG_qq = VG_Result
   csmax_qq = 10d0*csmax_qq    ! adjustment factor 
endif 

  PChannel = PChannel_aux


  if (PChannel.eq.2.and.fix_channels_ratio) then 
  adj_par = VG_gg/VG_qq*channels_ratio_fix/(one-channels_ratio_fix)
     else 
  adj_par = one 
  endif 


!------------------------------- set counts to zero for actual evaluation

      EvalCounter = 0d0
      AccepCounter = 0d0
      AccepCounter_q = 0d0
      AccepCounter_g = 0d0


      if (seed_random) call random_seed()

      do i=1,ncall
 
        call random_number(yRnd)


      dum = EvalCS_LO_ppllll(yRnd)
 
!  call vegas(EvalCS_LO_ppllll,VG_Result,VG_Error,VG_Chi2)
 

      enddo 

      print *, "Evaluation Counter: ",EvalCounter
      print *, "Acceptance  Counter: ",AccepCounter
      print *, "Acceptance  Counter_q: ",AccepCounter_q
      print *, "Acceptance  Counter_g: ",AccepCounter_g


  endif

return
END SUBROUTINE





SUBROUTINE OpenFiles()
use ModParameters
implicit none

   print *, ""
   print *, "Data file:         "//trim(DataFile)//'.dat'
   print *, "Vegas status file: "//trim(DataFile)//'.status'
   open(unit=14,file=trim(DataFile)//'.dat',form='formatted',access= 'sequential',status='replace')            ! Histogram file
   open(unit=15,file=trim(DataFile)//'.status',form='formatted',access= 'sequential',status='replace')         ! Vegas status file

return
END SUBROUTINE




SUBROUTINE CloseFiles()
implicit none

   close(14)
   close(15)

return
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
  write(14,"(A,2X,1F9.2,A)") "# run time =",RunTime/60d0,"min"
  write(* ,"(A,2X,1PE20.10,2X,1PE20.5)") "# TotCS =",VG_Result,VG_Error
  write(14,"(A,2X,1PE20.10,2X,1PE20.5)") "# TotCS =",VG_Result,VG_Error
  write(* ,"(A)") ""
  write(14,"(A)") ""
  do NHisto=1,NumHistograms
      write(* ,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      write(14,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      Integral = 0d0
      BinSize = Histo(NHisto)%BinSize * Histo(NHisto)%SetScale
      LowVal  = Histo(NHisto)%LowVal  * Histo(NHisto)%SetScale
      do NBin=1, Histo(NHisto)%NBins
          BinVal = (LowVal+(NBin-1)*BinSize)
          Hits   = Histo(NHisto)%Hits(NBin)
          if( unweighted ) then
              Value  = Histo(NHisto)%Value(NBin)/BinSize / dble(EvalCounter)
              Integral = Integral + Histo(NHisto)%Value(NBin)
              Error  = 1d0/dsqrt(dble(Hits))
          else
              Value  = Histo(NHisto)%Value(NBin)/BinSize/it
              Integral = Integral + Histo(NHisto)%Value(NBin)/it
              Error  = 1d0/(BinSize)/it * dsqrt( Histo(NHisto)%Value2(NBin) - 1d0/it/ncall*Histo(NHisto)%Value(NBin)**2 )
          endif
          write(* ,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
          write(14,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"
      enddo
      Integral = Integral + (Histo(NHisto)%Value(0)+Histo(NHisto)%Value(Histo(NHisto)%NBins+1))/it
      write(*,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(14,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(*,*) ""
      write(14,*) ""
  enddo

return
END SUBROUTINE
