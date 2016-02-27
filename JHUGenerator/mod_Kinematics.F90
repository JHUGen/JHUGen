MODULE ModKinematics
implicit none
save


type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8) :: Value(0:300)
    real(8) :: Value2(0:300)
    integer :: Hits(0:300)
    character :: Info*(50)
end type

integer,public :: it_sav
type(Histogram),allocatable :: Histo(:)

contains




SUBROUTINE WriteOutEvent(Mom,MY_IDUP,ICOLUP,MomFSPartons,EventWeight,EventInfoLine,PDFLine,MOTHUP_Parton)
use ModParameters
use modMisc
implicit none
real(8) :: Mom(1:4,1:6)
integer,parameter :: maxpart=15! this parameter should match the one in main.F90
real(8),optional :: MomFSPartons(:,:)
real(8),optional :: EventWeight
character(len=160),optional :: EventInfoLine,PDFLine
integer,optional :: MOTHUP_Parton(:,:)
real(8) :: Spin, Lifetime,s34,s36,s45,s56,smallestInv
real(8) :: XFV(1:4), Z1FV(1:4), Z2FV(1:4)
real(8) :: MomDummy(1:4,1:4+maxpart)
real(8) :: Part1Mass,Part2Mass,XMass,V1Mass,V2Mass,L11Mass,L12Mass,L21Mass,L22Mass,tmp,PartonMass(1:4+maxpart)
integer :: a,b,c,NumFSPartons
integer :: MY_IDUP(:),ICOLUP(:,:)
integer :: LHE_IDUP(1:7+maxpart),i,ISTUP(1:7+maxpart),MOTHUP(1:2,1:7+maxpart)
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,mZ1,mZ2,HiggsDKLength
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"



! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA



do i=1,9
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
enddo

! NUP changes for gamma gamma final state
if ( LHE_IDUP(4).eq.22 .and. LHE_IDUP(5).eq.22 ) then! photon+photon FS
    NUP=5
elseif ( LHE_IDUP(4).ne.22 .and. LHE_IDUP(5).eq.22 ) then! Z+photon FS
    NUP=7
else! ZZ or WW FS
    NUP=9
endif



IDPRUP=Process
if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif

SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas

ISTUP(1) = -1
ISTUP(2) = -1
ISTUP(3) = 2
ISTUP(4) = 2
ISTUP(5) = 2
ISTUP(6) = 1
ISTUP(7) = 1
ISTUP(8) = 1
ISTUP(9) = 1

if ( LHE_IDUP(4).eq.22 .and. LHE_IDUP(5).eq.22 ) then! photon+photon FS
    ISTUP(4) = 1
    ISTUP(5) = 1
elseif ( LHE_IDUP(4).ne.22 .and. LHE_IDUP(5).eq.22 ) then! Z+photon FS
    ISTUP(4) = 2
    ISTUP(5) = 1
    ISTUP(6) = 1
    ISTUP(7) = 1
endif



MOTHUP(1,1) = 0
MOTHUP(2,1) = 0
MOTHUP(1,2) = 0
MOTHUP(2,2) = 0
MOTHUP(1,3) = 1
MOTHUP(2,3) = 2
MOTHUP(1,4) = 3
MOTHUP(2,4) = 3
MOTHUP(1,5) = 3
MOTHUP(2,5) = 3
MOTHUP(1,6) = 4
MOTHUP(2,6) = 4
MOTHUP(1,7) = 4
MOTHUP(2,7) = 4
MOTHUP(1,8) = 5
MOTHUP(2,8) = 5
MOTHUP(1,9) = 5
MOTHUP(2,9) = 5




if( present(MOTHUP_Parton) ) then
   MOTHUP(1:2,1) = MOTHUP_Parton(1:2,1)
   MOTHUP(1:2,2) = MOTHUP_Parton(1:2,2)
endif


NumFSPartons=0
if( present(MomFSPartons) ) then ! add additional FS partons
    NumFSPartons = size(MomFSPartons,2)
    do a=1,NumFSPartons
       NUP=NUP+1
       LHE_IDUP(9+a) = convertLHE( MY_IDUP(9+a  ) )
       ISTUP(9+a) = 1
       if( present(MOTHUP_Parton) ) then
          MOTHUP(1:2,9+a) = MOTHUP_Parton(1:2,2+a)
       else
          MOTHUP(1:2,9+a) = (/ mod(a,2)+1,mod(a,2)+1 /) !  = (1,1), (2,2), (1,1), ...
       endif
       MomDummy(1,6+a) = 100.0d0*MomFSPartons(1,a)
       MomDummy(2,6+a) = 100.0d0*MomFSPartons(2,a)
       MomDummy(3,6+a) = 100.0d0*MomFSPartons(3,a)
       MomDummy(4,6+a) = 100.0d0*MomFSPartons(4,a)
    enddo
endif


LHE_IDUP(3) = 25
if( Process.eq.1 ) LHE_IDUP(3) = 32
if( Process.eq.2 ) LHE_IDUP(3) = 39
Lifetime = 0.0d0
Spin = 0.1d0
call getHiggsDecayLength(HiggsDKLength)


do a=1,6
    MomDummy(1,a) = 100.0d0*Mom(1,a)
    MomDummy(2,a) = 100.0d0*Mom(2,a)
    MomDummy(3,a) = 100.0d0*Mom(3,a)
    MomDummy(4,a) = 100.0d0*Mom(4,a)
enddo

do b=1,4! V boson momenta
    Z1FV(b) = MomDummy(b,3)+MomDummy(b,4)
    Z2FV(b) = MomDummy(b,5)+MomDummy(b,6)
enddo

do c=1,4! X resonance momentum
    XFV(c) = Z1FV(c) + Z2FV(c)
enddo


! check energy-momentum conservation when we don't use "adjusted kinematics"
if( (OffShellV1).or.(OffShellV2).or.(IsAPhoton(DecayMode2)) ) then
    do c=1,4! loop over components of 4-momentum
          tmp=Mom(c,1)+Mom(c,2)-Mom(c,3)-Mom(c,4)-Mom(c,5)-Mom(c,6)
          do a=1,NumFSPartons
              tmp=tmp-MomFSPartons(c,a)
          enddo
          if( abs(tmp)/Mom(1,1).gt.1d-5 ) print *, "Error: energy-momentum violation!",c,abs(tmp)/Mom(1,1)
    enddo
endif



!  associte lepton pairs to MOTHUP
if( (IsAZDecay(DecayMode1)).and.(IsAZDecay(DecayMode2)) .and. abs(LHE_IDUP(7)).eq.abs(LHE_IDUP(9)) ) then 
     s34 = Get_MInv( Mom(1:4,3)+Mom(1:4,4) )
     s56 = Get_MInv( Mom(1:4,5)+Mom(1:4,6) )
     s36 = Get_MInv( Mom(1:4,3)+Mom(1:4,6) )
     s45 = Get_MInv( Mom(1:4,4)+Mom(1:4,5) )
     smallestInv = minloc((/dabs(s34-M_V),dabs(s56-M_V),dabs(s36-M_V),dabs(s45-M_V)/),1)        
     if( smallestInv.eq.3 .or. smallestInv.eq.4 ) then
        call swapi(MOTHUP(1,6),MOTHUP(1,8))
        call swapi(MOTHUP(2,6),MOTHUP(2,8))
        call swapi(ICOLUP(1,6),ICOLUP(1,8))
        call swapi(ICOLUP(2,6),ICOLUP(2,8))
        Z1FV(1:4) = MomDummy(1:4,3)+MomDummy(1:4,6)
        Z2FV(1:4) = MomDummy(1:4,5)+MomDummy(1:4,4)
     endif
endif





! calculating and checking masses
tmp = (MomDummy(1:4,1)).dot.(MomDummy(1:4,1))
if( tmp.lt. -1d-3 ) print *, "Error 1: large negative mass!",tmp
Part1Mass = dSQRT(dabs(tmp))

tmp = (MomDummy(1:4,2)).dot.(MomDummy(1:4,2))
if( tmp.lt. -1d-3 ) print *, "Error 2: large negative mass!",tmp
Part2Mass = dSQRT(dabs(tmp))

tmp = (XFV(1:4)).dot.(XFV(1:4))
if( tmp.lt. -1d-3 ) print *, "Error 3: large negative mass!",tmp
XMass = dSQRT(dabs(tmp))

tmp = (Z1FV(1:4)).dot.(Z1FV(1:4))
if( tmp.lt. -1d-3 ) print *, "Error 4: large negative mass!",tmp
V1Mass = dSQRT(dabs(tmp))
if( V1Mass.lt.1d-5 ) V1Mass=0d0

tmp = (Z2FV(1:4)).dot.(Z2FV(1:4))
if( tmp.lt. -1d-3 ) print *, "Error 5: large negative mass!",tmp
V2Mass = dSQRT(dabs(tmp))
if( V2Mass.lt.1d-5 ) V2Mass=0d0

tmp = (MomDummy(1:4,3)).dot.(MomDummy(1:4,3))
if( tmp.lt. -1d-3 ) print *, "Error 6: large negative mass!",tmp
L12Mass = dSQRT(dABS(tmp))
if( L12Mass.lt.1d-6 ) L12Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,3) = MomDummy(1,3) + 1d-7

tmp = (MomDummy(1:4,4)).dot.(MomDummy(1:4,4))
if( tmp.lt. -1d-3 ) print *, "Error 7: large negative mass!",tmp
L11Mass = dSQRT(dABS(tmp))
if( L11Mass.lt.1d-6 ) L11Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,4) = MomDummy(1,4) + 1d-7

tmp = (MomDummy(1:4,5)).dot.(MomDummy(1:4,5))
if( tmp.lt. -1d-3 ) print *, "Error 8: large negative mass!",tmp
L22Mass = dSQRT(dABS(tmp))
if( L22Mass.lt.1d-6 ) L22Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,5) = MomDummy(1,5) + 1d-7

tmp = (MomDummy(1:4,6)).dot.(MomDummy(1:4,6))
if( tmp.lt. -1d-3 ) print *, "Error 9: large negative mass!",tmp
L21Mass = dSQRT(dABS(tmp))
if( L21Mass.lt.1d-6 ) L21Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,6) = MomDummy(1,6) + 1d-7

do a=1,NumFSPartons
    tmp = (MomDummy(1:4,6+a)).dot.(MomDummy(1:4,6+a))
    PartonMass(a) = dSQRT(dabs(tmp))
enddo






write(io_LHEOutFile,"(A)") "<event>"
if( ReadLHEFile .and. importExternal_LHEinit .and. present(EventInfoLine) ) then
   write(io_LHEOutFile,"(I2,X,A)") NUP,trim(EventInfoLine)
else
   write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
!  in order of appearance:
!  (*) number of particles in the event
!  (*) process ID (user defined)
!  (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
!  (*) pdf factorization scale in GeV
!  (*) alpha_QED coupling for this event
!  (*) alpha_s coupling for this event
endif

! parton_a
i=1
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,1),MomDummy(1,1),Part1Mass,Lifetime,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,2),MomDummy(1,2),Part2Mass,Lifetime,Spin

! X
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),XFV(2:4),XFV(1),XMass,HiggsDKLength,Spin

! V1
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Z1FV(2:4),Z1FV(1),V1Mass,Lifetime,Spin

! V2
i=5
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Z2FV(2:4),Z2FV(1),V2Mass,Lifetime,Spin

! decay product 1 (V1): l-, nu or q
i=7
if (LHE_IDUP(i).gt.-9000) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,3),MomDummy(1,3),L12Mass,Lifetime,Spin
endif

! decay product 2 (V1): l+, nubar or qbar
i=6
if (LHE_IDUP(i).gt.-9000) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,4),MomDummy(1,4),L11Mass,Lifetime,Spin
endif

! decay product 1 (V2): l-, nu or q
i=9
if (LHE_IDUP(i).gt.-9000) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,5),MomDummy(1,5),L22Mass,Lifetime,Spin
endif

! decay product 2 (V2): l+, nubar or qbar
i=8
if (LHE_IDUP(i).gt.-9000) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,6),MomDummy(1,6),L21Mass,Lifetime,Spin
endif

! additional F.S. partons
do a=1,NumFSPartons
    write(io_LHEOutFile,fmt1) LHE_IDUP(9+a),ISTUP(9+a), MOTHUP(1,9+a),MOTHUP(2,9+a), ICOLUP(1,9+a),ICOLUP(2,9+a),MomDummy(2:4,6+a),MomDummy(1,6+a),PartonMass(a),Lifetime,Spin
enddo

if( present(PDFLine) ) then
  if( PDFLine.ne."" ) write(io_LHEOutFile,"(A)") trim(PDFLine)
endif

write(io_LHEOutFile,"(A)") "</event>"


! print * ,"check ", LHE_IDUP(6),MomDummy(1:4,4)
! print * ,"check ", LHE_IDUP(7),MomDummy(1:4,3)
! print * ,"check ", LHE_IDUP(8),MomDummy(1:4,6)
! print * ,"check ", LHE_IDUP(9),MomDummy(1:4,5)

! print * ,"check ", LHE_IDUP(6),L11Mass
! print * ,"check ", LHE_IDUP(7),L12Mass
! print * ,"check ", LHE_IDUP(8),L21Mass
! print * ,"check ", LHE_IDUP(9),L22Mass
! pause

RETURN
END SUBROUTINE




SUBROUTINE WriteOutEvent_NEW(NUP,IDUP,ISTUP,MOTHUP,ICOLUP,Mom,HiggsDK_Mom,Mass,iHiggs,HiggsDK_IDUP,HiggsDK_ICOLUP,EventProcessId,EventWeight,EventScaleAqedAqcd,BeginEventLine,InputFmt0,Empty)
use ModParameters
use modMisc
implicit none
real(8) :: Mom(:,:),HiggsDK_Mom(:,:),Mass(:)
! real(8),optional :: MomFSPartons(:,:)
integer,optional :: EventProcessId
real(8),optional :: EventWeight
real(8),optional :: EventScaleAqedAqcd(1:3)
character(len=*),optional :: BeginEventLine
character(len=*),optional :: InputFmt0
! integer,optional :: MOTHUP_Parton(:,:)
real(8) :: Spin, Lifetime, s34,s56,s36,s45,smallestInv
integer :: IDUP(:),ISTUP(:),MOTHUP(:,:),ICOLUP(:,:)
integer :: HiggsDK_IDUP(:),HiggsDK_ICOLUP(:,:),HiggsDK_ISTUP(4:9),HiggsDK_MOTHUP(1:2,4:9)
integer,parameter :: maxpart=30
integer :: i,iHiggs
integer :: NUP,NUP_NEW,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: DefaultFmt0 = "I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7"
character(len=*),parameter :: Fmt1 = "6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0"
integer :: indent
character(len=150) :: IndentedFmt0, IndentedFmt1
logical,optional :: Empty
logical :: IsEmpty


!   For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
!   The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


!     assignments:
!     HiggsDK_Mom(1:4,1) --> HiggsDK_IDUP(4)! V1
!     HiggsDK_Mom(1:4,2) --> HiggsDK_IDUP(5)! V2
!     HiggsDK_Mom(1:4,3) --> HiggsDK_IDUP(7)! l- 
!     HiggsDK_Mom(1:4,4) --> HiggsDK_IDUP(6)! l+
!     HiggsDK_Mom(1:4,5) --> HiggsDK_IDUP(9)! l-
!     HiggsDK_Mom(1:4,6) --> HiggsDK_IDUP(8)! l+

    if( present(Empty) ) then
        IsEmpty = Empty
    else
        IsEmpty = .false.
    endif

    ! NUP changes for gamma gamma final state
    if( IsEmpty ) then
        NUP_NEW = -NUP
    elseif ( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then! photon+photon FS
        NUP_NEW = 2
    elseif ( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then! Z+photon FS
        NUP_NEW = 4
    else! ZZ or WW FS
        NUP_NEW = 6
    endif

    if( present(EventProcessId) .and. importExternal_LHEinit) then
        IDPRUP=EventProcessId
    else
        IDPRUP=Process
    endif
    if( present(EventWeight) ) then
        XWGTUP=EventWeight
    else
        XWGTUP=1.0d0
    endif
    if( present(EventScaleAqedAqcd)  .and. importExternal_LHEinit) then
        SCALUP=EventScaleAqedAqcd(1)
        AQEDUP=EventScaleAqedAqcd(2)
        AQCDUP=EventScaleAqedAqcd(3)
    else
        SCALUP=Mu_Fact * 100d0
        AQEDUP=alpha_QED
        AQCDUP=alphas
    endif
    ISTUP(iHiggs) = 2
    if ( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then! photon+photon FS
        HiggsDK_ISTUP(4:9) = (/1,1,0,0,0,0/)
        HiggsDK_MOTHUP(1:2,4) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,5) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,6:9) = 0
    elseif ( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then! Z+photon FS
        HiggsDK_ISTUP(4:9) = (/2,1,1,1,0,0/)
        HiggsDK_MOTHUP(1:2,4) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,5) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,6) = (/1,1/) + NUP
        HiggsDK_MOTHUP(1:2,7) = (/1,1/) + NUP
        HiggsDK_MOTHUP(1:2,8:9) = 0        
    else    
        HiggsDK_ISTUP(4:9) = (/2,2,1,1,1,1/)
        HiggsDK_MOTHUP(1:2,4) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,5) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,6) = (/1,1/) + NUP
        HiggsDK_MOTHUP(1:2,7) = (/1,1/) + NUP
        HiggsDK_MOTHUP(1:2,8) = (/2,2/) + NUP
        HiggsDK_MOTHUP(1:2,9) = (/2,2/) + NUP
    endif

    Lifetime = 0.0d0
    Spin = 0.1d0
    call getHiggsDecayLength(HiggsDKLength)

    !  associte lepton pairs to MOTHUP
    if( (IsAZDecay(DecayMode1)).and.(IsAZDecay(DecayMode2)) .and. abs(HiggsDK_IDUP(7)).eq.abs(HiggsDK_IDUP(9)) ) then 
        s34 = Get_MInv( HiggsDK_Mom(1:4,3)+HiggsDK_Mom(1:4,4) )
        s56 = Get_MInv( HiggsDK_Mom(1:4,5)+HiggsDK_Mom(1:4,6) )
        s36 = Get_MInv( HiggsDK_Mom(1:4,3)+HiggsDK_Mom(1:4,6) )
        s45 = Get_MInv( HiggsDK_Mom(1:4,4)+HiggsDK_Mom(1:4,5) )
        smallestInv = minloc((/dabs(s34-M_V),dabs(s56-M_V),dabs(s36-M_V),dabs(s45-M_V)/),1)        
        if( smallestInv.eq.3 .or. smallestInv.eq.4 ) then
            call swapi(HiggsDK_MOTHUP(1,6),HiggsDK_MOTHUP(1,8))
            call swapi(HiggsDK_MOTHUP(2,6),HiggsDK_MOTHUP(2,8))
            call swapi(HiggsDK_ICOLUP(1,6),HiggsDK_ICOLUP(1,8))
            call swapi(HiggsDK_ICOLUP(2,6),HiggsDK_ICOLUP(2,8))
            HiggsDK_Mom(1:4,1) = HiggsDK_Mom(1:4,3)+ HiggsDK_Mom(1:4,6)
            HiggsDK_Mom(1:4,2) = HiggsDK_Mom(1:4,4)+ HiggsDK_Mom(1:4,5)
        endif
    endif

    if (present(BeginEventLine)) then
        write(io_LHEOutFile, "(A)") trim(BeginEventLine)
        indent = 0
        do while (BeginEventLine(indent+1:indent+1).eq." ")
            indent = indent+1
        end do
    else
        write(io_LHEOutFile,"(A)") "<event>"
        indent = 0
    endif
    if (present(InputFmt0)) then
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
    write(io_LHEOutFile,IndentedFmt0) NUP+NUP_NEW,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    !  in order of appearance:
    !  (*) number of particles in the event
    !  (*) process ID (user defined)
    !  (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
    !  (*) pdf factorization scale in GeV
    !  (*) alpha_QED coupling for this event
    !  (*) alpha_s coupling for this event

    
    if( .not. IsEmpty ) then
!       write out existing particles
        do i = 1, NUP
            if( i.eq.iHiggs ) then 
               write(io_LHEOutFile,IndentedFmt1) IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),  &
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,HiggsDKLength, Spin           
            else
               write(io_LHEOutFile,IndentedFmt1) IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),  &
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,Lifetime, Spin   
            endif                          
        enddo


!       write new intermediate particles and Higgs decay products
        call swap_mom(HiggsDK_Mom(1:4,3),HiggsDK_Mom(1:4,4))! swap to account for flipped asignments
        call swap_mom(HiggsDK_Mom(1:4,5),HiggsDK_Mom(1:4,6))! swap to account for flipped asignments
        do i = 4,4 + (NUP_NEW-1)
            write(io_LHEOutFile,IndentedFmt1) HiggsDK_IDUP(i),HiggsDK_ISTUP(i), HiggsDK_MOTHUP(1,i),HiggsDK_MOTHUP(2,i), HiggsDK_ICOLUP(1,i),HiggsDK_ICOLUP(2,i),  &
                                              HiggsDK_Mom(2:4,i-3)/GeV,HiggsDK_Mom(1,i-3)/GeV, get_MInv(HiggsDK_Mom(1:4,i-3))/GeV, Lifetime, Spin   
        enddo
    endif

RETURN
END SUBROUTINE




SUBROUTINE WriteOutEvent_HFF(NUP,IDUP,ISTUP,MOTHUP,ICOLUP,Mom,HiggsDK_Mom,Mass,iHiggs,HiggsDK_IDUP,HiggsDK_ICOLUP,EventProcessId,EventWeight,EventScaleAqedAqcd,BeginEventLine,InputFmt0,Empty)
use ModParameters
use modMisc
implicit none
real(8) :: Mom(:,:),HiggsDK_Mom(:,:),Mass(:)
integer,optional :: EventProcessId
real(8),optional :: EventWeight
real(8),optional :: EventScaleAqedAqcd(1:3)
character(len=*),optional :: BeginEventLine
character(len=*),optional :: InputFmt0
real(8) :: Spin, Lifetime, s34,s56,s36,s45,smallestInv
integer :: IDUP(:),ISTUP(:),MOTHUP(:,:),ICOLUP(:,:)
integer :: HiggsDK_IDUP(:),HiggsDK_ICOLUP(:,:),HiggsDK_ISTUP(4:13),HiggsDK_MOTHUP(1:2,4:13)
integer,parameter :: maxpart=30
integer :: i,iHiggs
integer :: NUP,NUP_NEW,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: DefaultFmt0 = "I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7"
character(len=*),parameter :: Fmt1 = "6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0"
integer :: indent
character(len=150) :: IndentedFmt0, IndentedFmt1
integer, parameter :: inLeft=1, inRight=2, Hig=3, tauP=4, tauM=5, Wp=6, Wm=7,   nu=8, nubar_tau=9, lepP=10,   lepM=11, nu_tau=12, nubar=13
logical,optional :: Empty
logical :: IsEmpty


!   For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
!   The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


    if( present(Empty) ) then
        IsEmpty = Empty
    else
        IsEmpty = .false.
    endif

    if( IsEmpty ) then
        NUP_NEW = -NUP
    elseif( TauDecays.eq.0 ) then
        NUP_NEW = 2
    else
        NUP_NEW = 10
    endif


    if( present(EventProcessId) .and. importExternal_LHEinit) then
        IDPRUP=EventProcessId
    else
        IDPRUP=Process
    endif
    if( present(EventWeight) ) then
        XWGTUP=EventWeight
    else
        XWGTUP=1.0d0
    endif
    if( present(EventScaleAqedAqcd)  .and. importExternal_LHEinit) then
        SCALUP=EventScaleAqedAqcd(1)
        AQEDUP=EventScaleAqedAqcd(2)
        AQCDUP=EventScaleAqedAqcd(3)
    else
        SCALUP=Mu_Fact * 100d0
        AQEDUP=alpha_QED
        AQCDUP=alphas
    endif
    ISTUP(iHiggs) = 2
    if ( TauDecays.eq.0 ) then
        HiggsDK_ISTUP(4:5) = (/1,1/)
        HiggsDK_MOTHUP(1:2,4) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,5) = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,6:11) = 0
    else    
        HiggsDK_ISTUP(4:13) = (/2,2,2,2,1,1,1,1,1,1/)
        HiggsDK_MOTHUP(1:2,tauP)     = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,tauM)     = (/iHiggs,iHiggs/)
        HiggsDK_MOTHUP(1:2,Wp)       = (/1,1/) + NUP
        HiggsDK_MOTHUP(1:2,nubar_tau)= (/1,1/) + NUP
        HiggsDK_MOTHUP(1:2,Wm)       = (/2,2/) + NUP
        HiggsDK_MOTHUP(1:2,nu_tau)   = (/2,2/) + NUP
        HiggsDK_MOTHUP(1:2,nu)       = (/3,3/) + NUP
        HiggsDK_MOTHUP(1:2,lepP)     = (/3,3/) + NUP
        HiggsDK_MOTHUP(1:2,lepM)     = (/4,4/) + NUP
        HiggsDK_MOTHUP(1:2,nubar)    = (/4,4/) + NUP
    endif

    Lifetime = 0.0d0
    Spin = 0.1d0
    call getHiggsDecayLength(HiggsDKLength)


    if (present(BeginEventLine)) then
        write(io_LHEOutFile, "(A)") trim(BeginEventLine)
        indent = 0
        do while (BeginEventLine(indent+1:indent+1).eq." ")
            indent = indent+1
        end do
    else
        write(io_LHEOutFile,"(A)") "<event>"
        indent = 0
    endif
    if (present(InputFmt0)) then
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
    write(io_LHEOutFile,IndentedFmt0) NUP+NUP_NEW,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    !  in order of appearance:
    !  (*) number of particles in the event
    !  (*) process ID (user defined)
    !  (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
    !  (*) pdf factorization scale in GeV
    !  (*) alpha_QED coupling for this event
    !  (*) alpha_s coupling for this event


    if( .not. IsEmpty ) then
!       write out existing particles
        do i = 1, NUP
            if( i.eq.iHiggs ) then 
               write(io_LHEOutFile,IndentedFmt1) IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),  &
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,HiggsDKLength, Spin           
            else
               write(io_LHEOutFile,IndentedFmt1) IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),  &
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,Lifetime, Spin   
            endif                          
        enddo

!       write new intermediate particles and Higgs decay products
        do i = 4,4 + (NUP_NEW-1)
            write(io_LHEOutFile,IndentedFmt1) HiggsDK_IDUP(i),HiggsDK_ISTUP(i), HiggsDK_MOTHUP(1,i),HiggsDK_MOTHUP(2,i), HiggsDK_ICOLUP(1,i),HiggsDK_ICOLUP(2,i),  &
                                              HiggsDK_Mom(2:4,i)/GeV,HiggsDK_Mom(1,i)/GeV, get_MInv(HiggsDK_Mom(1:4,i))/GeV, Lifetime, Spin   
        enddo
    endif

RETURN
END SUBROUTINE


SUBROUTINE ShiftMass(p1,p2,m1,m2,p1hat,p2hat)
use ModMisc
implicit none
real(8),intent(in) :: p1(1:4),p2(1:4)
real(8) :: m1,m2,p1hat(1:4),p2hat(1:4)
real(8) :: xi,eta,a,b,c,p1sq,p2sq,p1p2

  p1sq = p1(1:4).dot.p1(1:4)
  p2sq = p2(1:4).dot.p2(1:4)
  p1p2 = p1(1:4).dot.p2(1:4)
  
  a = ( p1sq*p2(1:4) - p2sq*p1(1:4) + p1p2*(p2(1:4)-p1(1:4)) ).dot.( p1sq*p2(1:4) - p2sq*p1(1:4) + p1p2*(p2(1:4)-p1(1:4)) )
  b = ( p1sq+p2sq+2d0*p1p2+m2**2-m1**2 ) * ( p1p2**2 - p1sq*p2sq )
  c = 0.25d0*( p1sq+p2sq+2d0*p1p2+m2**2-m1**2 )**2*p1sq - (p1sq+p1p2)**2*m2**2
  eta = 1d0/2d0/a * ( -b - dsqrt( dabs(b**2 -4d0*a*c) ) )
  xi = ( p1sq+p2sq+2d0*p1p2 + m2**2 - m1**2 - 2d0*eta*(p2sq+p1p2) )/2d0/( p1sq + p1p2 )

  p2hat(1:4) = xi*p1(1:4) + eta*p2(1:4)
  p1hat(1:4) = (1d0-xi)*p1(1:4) + (1d0-eta)*p2(1:4)


! if( dabs( (p1hat.dot.p1hat)-m1**2 )/m1**2.gt.1d-3 ) then
!     print *, "1",p1hat.dot.p1hat , m1**2
!     print *, p1
!     print *, p1hat
!     print *, a,b,c,eta,xi,p1sq + p1p2
!     pause
! endif
! if( dabs( (p2hat.dot.p2hat)-m2**2 )/m2**2.gt.1d-3 ) then
!     print *, "2",p2hat.dot.p2hat , m2**2
!     print *, p2
!     print *, p2hat
!     print *, a,b,c,eta,xi,p1sq + p1p2
!     pause
! endif


RETURN
END SUBROUTINE




SUBROUTINE WriteOutEvent_HJ(Mom,MY_IDUP,ICOLUP,MomRealGlu,EventWeight,EventInfoLine)
use ModParameters
implicit none
real(8) :: Mom(1:4,1:4)
real(8),optional :: MomRealGlu(1:4)
character(len=160),optional :: EventInfoLine
real(8),optional :: EventWeight
real(8) :: Spin, Lifetime
! real(8) :: XFV(1:4), Z1FV(1:4), Z2FV(1:4)
real(8) :: MomDummy(1:4,1:5)
real(8) :: Part1Mass,Part2Mass,Part3Mass,Part4Mass,XMass
integer :: a,b,c
integer :: MY_IDUP(:),LHE_IDUP(1:10),i,ISTUP(1:10),MOTHUP(1:2,1:10),ICOLUP(:,:)
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"



! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


do i=1,4
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
enddo


IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas

ISTUP(1) = -1
ISTUP(2) = -1
ISTUP(3) = 1
ISTUP(4) = 1


MOTHUP(1,1) = 0
MOTHUP(2,1) = 0
MOTHUP(1,2) = 0
MOTHUP(2,2) = 0
MOTHUP(1,3) = 1
MOTHUP(2,3) = 2
MOTHUP(1,4) = 1
MOTHUP(2,4) = 2


NUP=4

if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif

Lifetime = 0.0d0
Spin = 0.1d0
call getHiggsDecayLength(HiggsDKLength)

do a=1,4
    MomDummy(1,a) = 100.0d0*Mom(1,a)
    MomDummy(2,a) = 100.0d0*Mom(2,a)
    MomDummy(3,a) = 100.0d0*Mom(3,a)
    MomDummy(4,a) = 100.0d0*Mom(4,a)
enddo

! do b=1,4
!     Z1FV(b) = MomDummy(b,1)-MomDummy(b,3)
!     Z2FV(b) = MomDummy(b,2)+MomDummy(b,4)
! enddo


Part1Mass = 0d0
Part2Mass = 0d0
Part4Mass = 0d0
XMass = M_Reso* 100d0



write(io_LHEOutFile,"(A)") "<event>"
if( .not. ReadLHEFile ) write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
! in order of appearance:
! (*) number of particles in the event
! (*) process ID (user defined)
! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
! (*) pdf factorization scale in GeV
! (*) alpha_QED coupling for this event 
! (*) alpha_s coupling for this event



! parton_a
i=1
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part1Mass,Lifetime,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part2Mass,Lifetime,Spin

! H
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),XMass,HiggsDKLength,Spin

! j1
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part4Mass,Lifetime,Spin


write(io_LHEOutFile,"(A)") "</event>"



END SUBROUTINE


SUBROUTINE WriteOutEvent_HJJ_fulldecay(Mom,MY_IDUP,ICOLUP,EventWeight)
use ModParameters
use ModMisc
implicit none
integer,parameter :: NUP=10
real(8) :: Mom(1:4,1:NUP)
real(8),optional :: EventWeight
integer :: MY_IDUP(1:NUP),ICOLUP(1:2,1:NUP),LHE_IDUP(1:NUP),ISTUP(1:NUP),MOTHUP(1:2,1:NUP)
integer :: IDPRUP,i
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,MomDummy(1:4,1:NUP),TheMass
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10


! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA

IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas


MOTHUP(1:2,inTop) = (/0,0/);          ISTUP(inTop) = -1
MOTHUP(1:2,inBot) = (/0,0/);          ISTUP(inBot) = -1

MOTHUP(1:2,outTop)= (/inTop,inTop/);  ISTUP(outTop)= +1
MOTHUP(1:2,outBot)= (/inBot,inBot/);  ISTUP(outBot)= +1

MOTHUP(1:2,V1)    = (/inTop,inBot/);  ISTUP(V1)    = +2
MOTHUP(1:2,V2)    = (/inTop,inBot/);  ISTUP(V2)    = +2

MOTHUP(1:2,Lep1P) = (/V1,V1/);        ISTUP(Lep1P) = +1
MOTHUP(1:2,Lep1M) = (/V1,V1/);        ISTUP(Lep1M) = +1
MOTHUP(1:2,Lep2P) = (/V2,V2/);        ISTUP(Lep2P) = +1
MOTHUP(1:2,Lep2M) = (/V2,V2/);        ISTUP(Lep2M) = +1


   
if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif
Lifetime = 0.0d0
Spin     = 0.1d0


do i=1,NUP
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
    MomDummy(1,i) = 100.0d0*Mom(1,i)
    MomDummy(2,i) = 100.0d0*Mom(2,i)
    MomDummy(3,i) = 100.0d0*Mom(3,i)
    MomDummy(4,i) = 100.0d0*Mom(4,i)
enddo



write(io_LHEOutFile,"(A)") "<event>"
write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
! in order of appearance:
! (*) number of particles in the event
! (*) process ID (user defined)
! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
! (*) pdf factorization scale in GeV
! (*) alpha_QED coupling for this event 
! (*) alpha_s coupling for this event


do i=1,NUP
     TheMass = get_Minv(MomDummy(:,i))
     if( i.le.4  ) TheMass = 0.0d0  ! setting quark masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime,Spin
enddo

write(io_LHEOutFile,"(A)") "</event>"


RETURN
END SUBROUTINE





SUBROUTINE WriteOutEvent_TTBH(Mom,MY_IDUP,ICOLUP,EventWeight)
use ModParameters
use ModMisc
implicit none
real(8) :: Mom(1:4,1:13)
real(8),optional :: EventWeight
integer :: MY_IDUP(1:13),ICOLUP(1:2,1:13),LHE_IDUP(1:13),ISTUP(1:13),MOTHUP(1:2,1:13)
integer :: NUP,IDPRUP,i
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,MomDummy(1:4,1:13),TheMass
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13


! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA



IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas



MOTHUP(1:2,inLeft) = (/0,0/);             ISTUP(inLeft) = -1
MOTHUP(1:2,inRight)= (/0,0/);             ISTUP(inRight)= -1

MOTHUP(1:2,Hbos)   = (/inLeft,inRight/);  ISTUP(Hbos)   = +1
MOTHUP(1:2,tbar)   = (/inLeft,inRight/);  ISTUP(tbar)   = +2
MOTHUP(1:2,t)      = (/inLeft,inRight/);  ISTUP(t)      = +2

MOTHUP(1:2,bbar)   = (/tbar,tbar/);       ISTUP(bbar)   = +1
MOTHUP(1:2,Wm)     = (/tbar,tbar/);       ISTUP(Wm)     = +2
MOTHUP(1:2,lepM)   = (/Wm,Wm/);           ISTUP(lepM)   = +1
MOTHUP(1:2,nubar)  = (/Wm,Wm/);           ISTUP(nubar)  = +1

MOTHUP(1:2,b)      = (/t,t/);             ISTUP(b)      = +1
MOTHUP(1:2,Wp)     = (/t,t/);             ISTUP(Wp)     = +2
MOTHUP(1:2,lepP)   = (/Wp,Wp/);           ISTUP(lepP)   = +1
MOTHUP(1:2,nu)     = (/Wp,Wp/);           ISTUP(nu)     = +1


if( TopDecays.eq.0 ) then
   NUP = 5
   ISTUP(tbar)   = +1
   ISTUP(t)      = +1 
else
   NUP=13
endif

if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif
Lifetime = 0.0d0
Spin     = 0.1d0


do i=1,5
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
    MomDummy(1,i) = 100.0d0*Mom(1,i)
    MomDummy(2,i) = 100.0d0*Mom(2,i)
    MomDummy(3,i) = 100.0d0*Mom(3,i)
    MomDummy(4,i) = 100.0d0*Mom(4,i)
enddo

if( TopDecays.ne.0 ) then
      ! introduce b-quark mass for LHE output 
!       call ShiftMass(Mom(1:4,b),   Mom(1:4,Wp),m_bot,M_W,  MomDummy(1:4,b),   MomDummy(1:4,Wp) )
!       call ShiftMass(Mom(1:4,bbar),Mom(1:4,Wm),m_bot,M_W,  MomDummy(1:4,bbar),MomDummy(1:4,Wm) )
      MomDummy(1:4,b)   = Mom(1:4,b)
      MomDummy(1:4,bbar)= Mom(1:4,bbar)
      MomDummy(1:4,Wp)  = Mom(1:4,Wp)
      MomDummy(1:4,Wm)  = Mom(1:4,Wm)

      ! introduce lepton/quark masses for LHE output  
      call ShiftMass(Mom(1:4,LepP),Mom(1:4,Wp)-Mom(1:4,LepP), GetMass(MY_IDUP(LepP)),0d0,  MomDummy(1:4,LepP),MomDummy(1:4,Nu) )
      call ShiftMass(Mom(1:4,LepM),Mom(1:4,Wm)-Mom(1:4,LepM), GetMass(MY_IDUP(LepM)),0d0,  MomDummy(1:4,LepM),MomDummy(1:4,Nubar) )

      do i=6,13
          LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
          MomDummy(1,i) = 100.0d0*MomDummy(1,i)
          MomDummy(2,i) = 100.0d0*MomDummy(2,i)
          MomDummy(3,i) = 100.0d0*MomDummy(3,i)
          MomDummy(4,i) = 100.0d0*MomDummy(4,i)
      enddo
endif


write(io_LHEOutFile,"(A)") "<event>"
write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
! in order of appearance:
! (*) number of particles in the event
! (*) process ID (user defined)
! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
! (*) pdf factorization scale in GeV
! (*) alpha_QED coupling for this event 
! (*) alpha_s coupling for this event

do i=1,NUP
!      TheMass = GetMass( MY_IDUP(i) )*100d0
     TheMass = get_Minv(MomDummy(:,i))
     if( i.le.2  ) TheMass = 0.0d0  ! setting initial parton masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime,Spin
enddo
write(io_LHEOutFile,"(A)") "</event>"


RETURN
END SUBROUTINE





SUBROUTINE WriteOutEvent_BBBH(Mom,MY_IDUP,ICOLUP,EventWeight)
use ModParameters
implicit none
real(8) :: Mom(1:4,1:11)
real(8),optional :: EventWeight
integer :: MY_IDUP(1:11),ICOLUP(1:2,1:11),LHE_IDUP(1:5),ISTUP(1:5),MOTHUP(1:2,1:5)
integer :: NUP,IDPRUP,i
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,MomDummy(1:4,1:13)
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
integer, parameter :: bbar=4,b=5,Hbos=3,inLeft=1,inRight=2


! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


do i=1,5
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
enddo


IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas

ISTUP(1:5) = (/-1,-1,1,1,1/)


MOTHUP(1:2,inLeft)  = (/0,0/)
MOTHUP(1:2,inRight) = (/0,0/)
MOTHUP(1:2,Hbos)    = (/1,2/)
MOTHUP(1:2,bbar)    = (/1,2/)
MOTHUP(1:2,b)       = (/1,2/)


NUP = 5

   
if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif

Lifetime = 0.0d0
Spin     = 0.1d0

do i=1,11
    MomDummy(1,i) = 100.0d0*Mom(1,i)
    MomDummy(2,i) = 100.0d0*Mom(2,i)
    MomDummy(3,i) = 100.0d0*Mom(3,i)
    MomDummy(4,i) = 100.0d0*Mom(4,i)
enddo



write(io_LHEOutFile,"(A)") "<event>"
write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
! in order of appearance:
! (*) number of particles in the event
! (*) process ID (user defined)
! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
! (*) pdf factorization scale in GeV
! (*) alpha_QED coupling for this event 
! (*) alpha_s coupling for this event



! parton_a
i=1
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),0d0,Lifetime,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),0d0,Lifetime,Spin

! H
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),M_Reso*100d0,Lifetime,Spin

! bb
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),m_top*100d0,Lifetime,Spin

! b
i=5
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),m_top*100d0,Lifetime,Spin



write(io_LHEOutFile,"(A)") "</event>"




END SUBROUTINE







SUBROUTINE WriteOutEvent_TH(Mom,MY_IDUP,ICOLUP,EventWeight)
use ModParameters
use ModMisc
implicit none
real(8) :: Mom(1:4,1:9)
real(8),optional :: EventWeight
integer :: MY_IDUP(1:9),ICOLUP(1:2,1:9),LHE_IDUP(1:9),ISTUP(1:9),MOTHUP(1:2,1:9)
integer :: NUP,IDPRUP,i
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,MomDummy(1:4,1:13),TheMass
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9


! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA



IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas



MOTHUP(1:2,inLeft) = (/0,0/);             ISTUP(inLeft) = -1
MOTHUP(1:2,inRight)= (/0,0/);             ISTUP(inRight)= -1

MOTHUP(1:2,Hbos)   = (/inLeft,inRight/);  ISTUP(Hbos)   = +1
MOTHUP(1:2,t)      = (/inLeft,inRight/);  ISTUP(t)      = +2
MOTHUP(1:2,qout)   = (/inLeft,inRight/);  ISTUP(qout)   = +1

MOTHUP(1:2,b)      = (/t,t/);             ISTUP(b)      = +1
MOTHUP(1:2,W)      = (/t,t/);             ISTUP(W)      = +2
MOTHUP(1:2,lep)    = (/W,W/);             ISTUP(lep)    = +1
MOTHUP(1:2,nu)     = (/W,W/);             ISTUP(nu)     = +1


if( TopDecays.eq.0 ) then
   NUP = 5
   ISTUP(t)      = +1 
else
   NUP=9
endif

if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif
Lifetime = 0.0d0
Spin     = 0.1d0


do i=1,5
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
    MomDummy(1,i) = 100.0d0*Mom(1,i)
    MomDummy(2,i) = 100.0d0*Mom(2,i)
    MomDummy(3,i) = 100.0d0*Mom(3,i)
    MomDummy(4,i) = 100.0d0*Mom(4,i)
enddo

if( TopDecays.ne.0 ) then
      ! introduce b-quark mass for LHE output 
!       call ShiftMass(Mom(1:4,b),   Mom(1:4,W),m_bot,M_W,  MomDummy(1:4,b),   MomDummy(1:4,W) )
      MomDummy(1:4,b) = Mom(1:4,b)
      MomDummy(1:4,W) = Mom(1:4,W)

      ! introduce lepton/quark masses for LHE output  
      call ShiftMass(Mom(1:4,Lep),Mom(1:4,W)-Mom(1:4,Lep), GetMass(MY_IDUP(Lep)),0d0,  MomDummy(1:4,Lep),MomDummy(1:4,Nu) )

      do i=6,9
          LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
          MomDummy(1,i) = 100.0d0*MomDummy(1,i)
          MomDummy(2,i) = 100.0d0*MomDummy(2,i)
          MomDummy(3,i) = 100.0d0*MomDummy(3,i)
          MomDummy(4,i) = 100.0d0*MomDummy(4,i)
      enddo
endif



write(io_LHEOutFile,"(A)") "<event>"
write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
! in order of appearance:
! (*) number of particles in the event
! (*) process ID (user defined)
! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
! (*) pdf factorization scale in GeV
! (*) alpha_QED coupling for this event 
! (*) alpha_s coupling for this event

do i=1,NUP
!      TheMass = GetMass( MY_IDUP(i) )*100d0
     TheMass = get_Minv(MomDummy(:,i))
     if( i.le.2  ) TheMass = 0.0d0  ! setting initial parton masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime,Spin
enddo
write(io_LHEOutFile,"(A)") "</event>"


RETURN
END SUBROUTINE






SUBROUTINE WriteOutEvent_HVBF(Mom,MY_IDUP,ICOLUP,MomRealGlu,EventWeight,EventInfoLine)
use ModParameters
implicit none
real(8) :: Mom(1:4,1:5)
real(8),optional :: MomRealGlu(1:4)
character(len=170),optional :: EventInfoLine
real(8),optional :: EventWeight
real(8) :: Spin, Lifetime
! real(8) :: XFV(1:4), Z1FV(1:4), Z2FV(1:4)
real(8) :: MomDummy(1:4,1:5)
real(8) :: Part1Mass,Part2Mass,Part3Mass,Part4Mass,XMass
integer :: a,b,c
integer :: MY_IDUP(:),LHE_IDUP(1:10),i,ISTUP(1:10),MOTHUP(1:2,1:10),ICOLUP(:,:)
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"



! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


do i=1,5
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
enddo


IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas

ISTUP(1) = -1
ISTUP(2) = -1
ISTUP(3) = 1
ISTUP(4) = 1
ISTUP(5) = 1


! this is not correct: intermediate VV's need to be added
! MOTHUP(1,1) = 0
! MOTHUP(2,1) = 0
! MOTHUP(1,2) = 0
! MOTHUP(2,2) = 0
! MOTHUP(1,3) = 1
! MOTHUP(2,3) = 1
! MOTHUP(1,4) = 2
! MOTHUP(2,4) = 2
! MOTHUP(1,5) = 1
! MOTHUP(2,5) = 2

MOTHUP(1,1) = 0
MOTHUP(2,1) = 0
MOTHUP(1,2) = 0
MOTHUP(2,2) = 0
MOTHUP(1,3) = 1
MOTHUP(2,3) = 2
MOTHUP(1,4) = 1
MOTHUP(2,4) = 2
MOTHUP(1,5) = 1
MOTHUP(2,5) = 2


NUP=5

if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif

Lifetime = 0.0d0
Spin = 0.1d0
call getHiggsDecayLength(HiggsDKLength)

do a=1,5
    MomDummy(1,a) = 100.0d0*Mom(1,a)
    MomDummy(2,a) = 100.0d0*Mom(2,a)
    MomDummy(3,a) = 100.0d0*Mom(3,a)
    MomDummy(4,a) = 100.0d0*Mom(4,a)
enddo

! do b=1,4
!     Z1FV(b) = MomDummy(b,1)-MomDummy(b,3)
!     Z2FV(b) = MomDummy(b,2)+MomDummy(b,4)
! enddo


Part1Mass = 0d0
Part2Mass = 0d0
Part3Mass = 0d0
Part4Mass = 0d0
XMass = M_Reso* 100d0



write(io_LHEOutFile,"(A)") "<event>"
if( .not. ReadLHEFile ) write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
! in order of appearance:
! (*) number of particles in the event
! (*) process ID (user defined)
! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
! (*) pdf factorization scale in GeV
! (*) alpha_QED coupling for this event 
! (*) alpha_s coupling for this event



! parton_a
i=1
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part1Mass,Lifetime,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part2Mass,Lifetime,Spin

! j1
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part3Mass,Lifetime,Spin

! j2
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part4Mass,Lifetime,Spin

! H
i=5
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),XMass,HiggsDKLength,Spin


write(io_LHEOutFile,"(A)") "</event>"

! if( dabs(MomDummy(2,1)+MomDummy(2,2)+MomDummy(2,3)+MomDummy(2,4)+MomDummy(2,5)).gt. 1d-10 ) then
!     print *, "checker",MomDummy(2,1)+MomDummy(2,2)+MomDummy(2,3)+MomDummy(2,4)+MomDummy(2,5)
!     pause
! endif


END SUBROUTINE




SUBROUTINE WriteOutEvent_VHiggs(id,helicity,MomExt,inv_mass,EventWeight)
use ModParameters
implicit none
double precision, intent(in) :: MomExt(1:4,1:9), inv_mass(9)
!double precision, intent(in) :: beam_h(2)
double precision helicity(9)
real(8) :: MomDummy(1:4,1:9), MassDummy(9)
integer :: id(9)
integer :: ICOLUP(4,2)
real(8) :: EventWeight
real(8) :: Spin
integer :: i
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,HiggsDKLength
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"



MomDummy = MomExt*1d2
MassDummy = inv_mass*1d2

IDPRUP=Process
SCALUP=Mu_Fact * 100d0
AQEDUP=alpha_QED
AQCDUP=alphas

call getHiggsDecayLength(HiggsDKLength)


if(H_DK.eqv..false.)then
    NUP=6
else
    NUP=8
endif

    XWGTUP=EventWeight

write(io_LHEOutFile,"(A)") "<event>"
if( .not. ReadLHEFile ) write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP


!if((COLLIDER.ne.0) .and. (unweighted.eqv..false.))then
!    beam_id(1)=2212
!    beam_id(2)=2212
!endif
helicity(3:5)=0d0
Spin = 0.1d0

if(COLLIDER.eq.0)then
  ICOLUP=0
else
  if(id(1).gt.0)then
    ICOLUP(1,1)=503
    ICOLUP(1,2)=0
    ICOLUP(2,1)=0
    ICOLUP(2,2)=503
  else
    ICOLUP(1,1)=0
    ICOLUP(1,2)=503
    ICOLUP(2,1)=503
    ICOLUP(2,2)=0
  endif
endif

if((id(6).eq.convertLHE(Up_)).or.(id(6).eq.convertLHE(Dn_)).or.(id(6).eq.convertLHE(Str_)).or.(id(6).eq.convertLHE(Chm_)).or.(id(6).eq.convertLHE(Bot_)))then
    ICOLUP(3,1)=502
    ICOLUP(3,2)=0
    ICOLUP(4,1)=0
    ICOLUP(4,2)=502
elseif((id(6).eq.convertLHE(AUp_)).or.(id(6).eq.convertLHE(ADn_)).or.(id(6).eq.convertLHE(AStr_)).or.(id(6).eq.convertLHE(AChm_)).or.(id(6).eq.convertLHE(ABot_)))then
    ICOLUP(3,1)=0
    ICOLUP(3,2)=502
    ICOLUP(4,1)=502
    ICOLUP(4,2)=0
else
    ICOLUP(3:4,1:2)=0
endif
!!print *, helicity!!!!!!!!!!
do i=1,2
    write(io_LHEOutFile,fmt1) id(i), -1,0,0,ICOLUP(i,1),ICOLUP(i,2),MomDummy(2:4,i), MomDummy(1,i), 0.0d0, 0.0d0, Spin
enddo

write(io_LHEOutFile,fmt1) id(4), 2,1,2,0,0,MomDummy(2:4,4), MomDummy(1,4), MassDummy(4), 0d0, Spin

if(H_DK.eqv..true.)then
  write(io_LHEOutFile,fmt1) id(5), 2,1,2,0,0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength, Spin
else
  write(io_LHEOutFile,fmt1) id(5), 1,1,2,0,0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength, Spin
endif

write(io_LHEOutFile,fmt1) id(6), 1,3,3,ICOLUP(3,1),ICOLUP(3,2),MomDummy(2:4,6), MomDummy(1,6), MassDummy(6), 0.0d0, Spin

write(io_LHEOutFile,fmt1) id(7), 1,3,3,ICOLUP(4,1),ICOLUP(4,2),MomDummy(2:4,7), MomDummy(1,7), MassDummy(7), 0.0d0, Spin

if(H_DK.eqv..true.)then
write(io_LHEOutFile,fmt1) id(8), 1,4,4,501,0,MomDummy(2:4,8), MomDummy(1,8), MassDummy(8), 0.0d0, Spin

write(io_LHEOutFile,fmt1) id(9), 1,4,4,0,501,MomDummy(2:4,9), MomDummy(1,9), MassDummy(9), 0.0d0, Spin
endif
write(io_LHEOutFile,"(A)") "</event>"


END SUBROUTINE







SUBROUTINE EvalPhasespace_VDecay(VMom,MV,ML1,ML2,xRndPS,MomDK,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: VMom(1:4),MomChk(1:4,1:3)
real(8) :: MomDK(1:4,1:2)
real(8) :: xRndPS(1:2),MV,ML1,ML2
integer,parameter :: N2=2
real(8),parameter :: PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


   if( MV.ne.0d0 ) then
      call genps(2,MV,xRndPS(1:2),(/ML1,ML2/),MomDK(1:4,1:2),PSWgt2)
!     boost all guys to the V boson frame:
      call boost(MomDK(1:4,1),VMom(1:4),MV)
      call boost(MomDK(1:4,2),VMom(1:4),MV)
      PSWgt = PSWgt2*PiWgt2
   else! this is for photons
      MomDK(1:4,1) = VMom(1:4)
      MomDK(1:4,2) = 0d0
   endif


RETURN
END SUBROUTINE


SUBROUTINE EvalPhasespace_VBF(EHat,M_H,xRndPS,Mom,PSWgt)
  use modParameters
  use modMisc
  implicit none
  real(8), intent(in) :: Ehat, M_H, xRndPS(5)
  real(8), intent(out) :: Mom(1:4,1:5), PSWgt
  real(8) :: slocal, mhsq, emax, x7, x8, x9, x10, x13
  real(8) :: E6, cos2, sin2, phi2, cos6, sin6, phi6, qhsq, PSWup, PSWdn, PSWdc
  real(8) :: qh(4), n6(4)

  Mom = zero
  PSWgt = zero

  x7 = xRndPS(1)
  x8 = xRndPS(2)
  x9 = xRndPS(3)
  x10 = xRndPS(4)
  x13 = xRndPS(5)

  slocal = Ehat**2
  mhsq = m_H**2

  cos2 = one-two*x7
  sin2 = dsqrt(dabs(one-cos2**2))
  phi2 = two*pi * x8
  cos6 = -one + two*x9
  sin6 = dsqrt(dabs(one-cos6**2))
  phi6 = two*pi * x10
  qhsq = mhsq + (slocal-mhsq)*x13

  emax = (slocal-qhsq)/(two*Ehat)
  
  Mom(:,1) = Ehat/two * (/one,zero,zero,one/)
  Mom(:,2) = Ehat/two * (/one,zero,zero,-one/)

  Mom(:,3) = emax * (/one,sin2*dcos(phi2),sin2*dsin(phi2),cos2/)
  
  qh = Mom(:,1)+Mom(:,2)-Mom(:,3)
  n6 = (/one,sin6*dcos(phi6),sin6*dsin(phi6),cos6/)
  E6 = (qhsq-mhsq)/two/MinkowskyProduct(qh,n6)

  Mom(:,4) = E6 * n6(:)

  Mom(:,5) = Mom(:,1)+Mom(:,2)-Mom(:,3)-Mom(:,4) !-- Higgs momentum

  !-- PS weight below
  PSWup = one/8.0d0/pi * (one-qhsq/slocal)
  PSWdn = E6/4.0d0/pi/MinkowskyProduct(qh,n6)
  PSWdc = one !-- no decay for the time being

  PSWgt = PSWup * PSWdn * PSWdc * (slocal-mhsq)/(two*pi)
  
  RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_VBF_NEW(EHat,xRndPS,Mom,PSWgt)
use modParameters
use modMisc
implicit none
real(8) :: EHat,xRndPS(:),Mom(:,:),PSWgt
real(8) :: s1,s2,y1,y2,phi1,phi2,Th1,Th2
real(8) :: E1,cos_Th1,sin_Th1,sin_phi1,cos_phi1
real(8) :: E2,cos_Th2,sin_Th2,sin_phi2,cos_phi2
real(8) :: x1,x2,x3
integer,parameter :: N2=3
real(8),parameter :: PiWgt = (2d0*Pi)**(4-N2*3)

! xRndPS(1:5) = (/ 0.9d0, 0.9d0, 0.5d0, 0.4d0, 0.78d0 /)


!   s1   = EHat*0.5d0 * ( -1d0 + xRndPS(1) )
!   s2   = EHat*0.5d0 * ( -1d0 + xRndPS(2) )
  s1   = -EHat**2 * (1d0-M_Reso**2/EHat**2) * ( xRndPS(1) )
  s2   = -EHat**2 * (1d0-M_Reso**2/EHat**2) * ( xRndPS(2) )
  
  y1   = -10d0 + xRndPS(3)*(20d0)
  y2   = -10d0 + xRndPS(4)*(20d0) 
  phi1 = 2d0*Pi * xRndPS(5)
  
  E1 = -s1/Ehat * dexp(y1) * dcosh(y1)
  E2 = -s2/Ehat * dexp(-y2) * dcosh(-y2)
   
  Th1 = 2d0*datan( exp(-y1) )
  Th2 = 2d0*datan( exp(+y2) )
  
  x1 = 2d0*E1*E2*dsin(phi1)*dsin(Th1)*dsin(Th2)
  x2 = 2d0*E1*E2*dcos(phi1)*dsin(Th1)*dsin(Th2)
  x3 = 2d0*E1*E2 - (E1+E2)*Ehat - M_Reso**2 + s1 + s2 - E2*Ehat*dcos(Th2) + E1*dcos(Th1)*(Ehat + 2*E2*dcos(Th2))
 
if( x1**4 + x1**2*x2**2 - x1**2*x3**2 .lt. 0d0 ) then
 PSWgt = 0d0
 return
endif
 
 
!  print *, ""
!  print *, "Ehat",Ehat*100d0
!  print *, "E1",e1*100d0
!  print *, "E2",e2*100d0
!  print *, "-sqrt(s1)",-dsqrt(dabs(s1*100d0**2))
!  print *, "-sqrt(s2)",-dsqrt(dabs(s2*100d0**2))
!  print *, "y1",y1
!  print *, "y2",y2
!  
!  print *, "Th1,cos(Th1)",Th1,dcos(th1)
!  print *, "Th2,cos(Th2)",Th2,dcos(th2)
!  print *, "phi1,cos(phi1)",phi1,dcos(phi1)
!   
!  print *, "x1",x1
!  print *, "x2",x2
!  print *, "x3",x3
!  print *, x1**4 + x1**2*x2**2 - x1**2*x3**2
!  print *, (x1**2 + x2**2)    
!  
!   phi2 = dATan2( (-x3 + (x2**2*x3)/(x1**2 + x2**2) + x2*dsqrt(-(x1**2*(-x1**2 - x2**2 + x3**2)))/(x1**2 + x2**2))/x1,   &
!                  (-(x2*x3) - dsqrt(x1**4 + x1**2*x2**2 - x1**2*x3**2))/(x1**2 + x2**2)                                  &
!                )
!    print *, "sol1",phi2
!    
!    
  phi2 = dATan2( (-x3 + (x2**2*x3)/(x1**2 + x2**2) - x2*dsqrt(-(x1**2*(-x1**2 - x2**2 + x3**2)))/(x1**2 + x2**2))/x1,   &
                 (-(x2*x3) + dsqrt(x1**4 + x1**2*x2**2 - x1**2*x3**2))/(x1**2 + x2**2)                                  &
               )
!    print *, "sol2",phi2
   
!    print *, "chekcer",   -(x1*x3)/(x1**2+x2**2)  +  dsqrt((x1**2 * x3**2)/(x1**2+x2**2)**2/4d0 +(x2**2-x3**2)/(x1**2+x2**2) ) 
!    print *, "chekcer",   -(x1*x3)/(x1**2+x2**2)  -  dsqrt((x1**2 * x3**2)/(x1**2+x2**2)**2/4d0 +(x2**2-x3**2)/(x1**2+x2**2) ) 
!    print *, "chekcer",   (x1**2 * x3**2)/(x1**2+x2**2)**2/4d0 + (x2**2-x3**2)/(x1**2+x2**2)    
!    pause
 
  
  sin_Th1 = dsin(Th1)
  cos_Th1 = dcos(Th1)
  sin_Th2 = dsin(Th2)
  cos_Th2 = dcos(Th2)
  sin_phi1= dsin(phi1)
  cos_phi1= dcos(phi1)
  sin_phi2= dsin(phi2)
  cos_phi2= dcos(phi2)
    
  Mom(1:4,1) = EHat*0.5d0 *(/1d0,0d0,0d0,+1d0/)
  Mom(1:4,2) = EHat*0.5d0 *(/1d0,0d0,0d0,-1d0/)
  Mom(1:4,3) = E1 * (/ 1d0, sin_Th1*sin_phi1, sin_Th1*cos_phi1, cos_Th1 /)
  Mom(1:4,4) = E2 * (/ 1d0, sin_Th2*sin_phi2, sin_Th2*cos_phi2, cos_Th2 /)
  Mom(1:4,5) = Mom(1:4,1)+Mom(1:4,2) - Mom(1:4,3)- Mom(1:4,4)- Mom(1:4,5) 

  PSWgt = 1d0/dcosh(y1)**2/dcosh(y2)**2 * PiWgt * (EHat**2 * (1d0-M_Reso**2/EHat**2))**2 * (20d0)**2 * 2d0*pi  & 
          / dabs( x1*dcos(phi2) - x2*dsin(phi2) )
  
RETURN
END SUBROUTINE




SUBROUTINE EvalPhasespace_2to2(EHat,Masses,xRndPS,Mom,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: EHat,Masses(1:2)
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
integer,parameter :: N2=2
real(8),parameter :: PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)

   call genps(2,Ehat,xRndPS(1:2),Masses,Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

return
END SUBROUTINE


SUBROUTINE EvalPhasespace_2(EHat,Masses,xRndPS,Mom,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: EHat,Masses(1:2)
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:2),xRndPS(1:2)
integer,parameter :: N2=2
real(8),parameter :: PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)

   call genps(2,Ehat,xRndPS(1:2),Masses,Mom(1:4,1:2),PSWgt)
   PSWgt = PSWgt*PiWgt2

return
END SUBROUTINE




subroutine AdjustKinematics(eta1,eta2,MomExt,MomDK,xgr,xz2,xz1,MomExt_f,MomDK_f)
 use ModParameters
 use ModMisc
#if compiler==1
 use ifport
#endif
 implicit none
 real(8) :: eta1,eta2,MomExt(1:4,1:4),MomDK(1:4,1:4),xgr,xz1,xz2
 real(8) :: MomG(1:4),MomBoost(1:4),MomZ1(1:4),MomZ2(1:4)
 real(8) :: Moml1_1(1:4),Moml1_2(1:4),Moml2_1(1:4),Moml2_2(1:4)
 real(8) :: EZ,  EZ1, EZ2, pz, MZ1, MZ2,  xmax, xx1, xx2
 real(8) :: pz12, MomExt_f(1:4,1:4), MomDK_f(1:4,1:4)
 real(8) :: MZ, MG, MomG_f(1:4)

!---- summarize all the momenta / before the invariant mass adjustments

    MomZ1(1:4) = MomExt(1:4,3)
    MomZ2(1:4) = MomExt(1:4,4)

    MomG = MomZ1+MomZ2

    Moml1_1(1:4) = MomDK(1:4,1)
    Moml1_2(1:4) = MomDK(1:4,2)

    Moml2_1(1:4) = MomDK(1:4,3)
    Moml2_2(1:4) = MomDK(1:4,4)

!----- begin with the invariant mass adjustment starting from the ``Gravition''

    MomExt_f(:,:) = 0d0

    xmax = atan((Collider_Energy**2-M_Reso**2)/M_Reso/Ga_Reso) + atan(M_Reso/Ga_Reso)
    MG=dsqrt(M_Reso**2 + M_Reso*Ga_Reso*tan(xgr*xmax-atan(M_Reso/Ga_Reso)))

    if ((MG/M_Reso*eta1.gt.1d0).or.(MG/M_Reso*eta2.gt.1d0)) return

    MomExt_f(:,1) =MG/M_Reso*MomExt(:,1)
    MomExt_f(:,2) =MG/M_Reso*MomExt(:,2)

    MomG_f(1:4) = MomExt_f(1:4,1)+MomExt_f(1:4,2)


!--------------------------------------------------------

    MomBoost(1) = MomG(1)
    MomBoost(2:4) = -MomG(2:4)

    call boost(MomZ1(1:4),MomBoost(1:4),M_Reso)
    call boost(MomZ2(1:4),MomBoost(1:4),M_Reso)


 !-- energies and momenta of the Z's in Gr rest frame
     EZ = M_Reso/two
     pz = sqrt(EZ**2 - M_V**2)

 !--- generate two more random numbers, to get invariant masses of the two Z's

    xmax = atan((mG**2-M_V**2)/M_V/Ga_V) + atan(M_V/Ga_V)
    MZ1=dsqrt(M_V**2 + M_V*Ga_V*tan(xz2*xmax-atan(M_V/Ga_V)))
    MZ2=dsqrt(M_V**2 + M_V*Ga_V*tan(xz1*xmax-atan(M_V/Ga_V)))


    if (mG.lt.(MZ1 + MZ2))  return                ! reject events which can not happen

    EZ1 = (mG**2 + MZ1**2 - MZ2**2)/2d0/mG
    EZ2 = (mG**2 + MZ2**2 - MZ1**2)/2d0/mG
    pz12 = sqrt(EZ1**2 - MZ1**2)

 !-- calculating the proper Z-four vectors in the Grav rest frame

    MomZ1(1) = MomZ1(1)*EZ1/Ez
    MomZ1(2:4) = MomZ1(2:4)*pz12/pz

    MomZ2(1) = MomZ2(1)*EZ2/Ez
    MomZ2(2:4)= MomZ2(2:4)*pz12/pz

 !-- boost the new vectors into the lab frame

    call boost(MomZ1(1:4),MomG_f(1:4),mG)
    call boost(MomZ2(1:4),MomG_f(1:4),mG)


 !-- now, do the same exercise with the leptons

 !--- first Z
    MomBoost(1)   =  MomExt(1,3)
    MomBoost(2:4) = -MomExt(2:4,3)

    call boost(Moml1_1(1:4),MomBoost,M_V)
    call boost(Moml1_2(1:4),MomBoost,M_V)

    Moml1_1 = MZ1/M_V*Moml1_1
    Moml1_2 = MZ1/M_V*Moml1_2

    call boost(Moml1_1(1:4),MomZ1(1:4),mZ1)
    call boost(Moml1_2(1:4),MomZ1(1:4),mZ1)

 !-- second Z
    MomBoost(1)   =  MomExt(1,4)
    MomBoost(2:4) = -MomExt(2:4,4)

    call boost(Moml2_1(1:4),MomBoost,M_V)
    call boost(Moml2_2(1:4),MomBoost,M_V)

    Moml2_1 = MZ2/M_V*Moml2_1
    Moml2_2 = MZ2/M_V*Moml2_2

    call boost(Moml2_1(1:4),MomZ2(1:4),mZ2)
    call boost(Moml2_2(1:4),MomZ2(1:4),mZ2)


 !--- now the full collection of all 4-vectors ; with off-shell kinematics

    MomExt_f(1:4,3) = MomZ1(1:4)
    MomExt_f(1:4,4) = MomZ2(1:4)

    MomDK_f(1:4,1) = Moml1_1(1:4)
    MomDK_f(1:4,2) = Moml1_2(1:4)
    MomDK_f(1:4,3) = Moml2_1(1:4)
    MomDK_f(1:4,4) = Moml2_2(1:4)


end subroutine AdjustKinematics





SUBROUTINE Kinematics(NumPart,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:),MomDK(:,:), mZ1, mZ2, MReso, ml12,ml34,ml14,ml32
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NumPart,NBin(:)
real(8) :: pT_lepM,pT_lepP,y_lepM,y_lepP,MomFerm(1:4),MomZ2(1:4),MomReso(1:4),CosTheta1,Phi,Phi1,signPhi,signPhi1
real(8) :: CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,CosThetaZ,CosThetaStar
real(8),parameter :: Rsep_ll=0.2d0

!  W^+/Z(MomExt(:,3))    -->    e^+(MomDK(:,2)) + nu_e(MomDK(:,1))       /   e^+(MomDK(:,2)) + e^-(MomDK(:,1))
!  W^-/Z(MomExt(:,4))    -->    e^-(MomDK(:,3)) + nubar_e(MomDK(:,4))    /   e^-(MomDK(:,3)) + e^+(MomDK(:,4))

!   MomDK(:,i): i=1 fermion
!   MomDK(:,i): i=2 anti-fermion
!   MomDK(:,i): i=3 fermion'
!   MomDK(:,i): i=4 anti-fermion'

      applyPSCut = .false.




!--- compute the invariant mass of the resonance
      MReso = sqrt(abs(   2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))  ))



! !     associte lepton pairs 
!       mZ1 = Get_MInv( MomDK(1:4,1)+MomDK(1:4,2) )
!       mZ2 = Get_MInv( MomDK(1:4,1)+MomDK(1:4,4) )
!       if( dabs(mZ1-m_V) .lt. dabs(mZ2-m_V)  ) then
!           MomLept(1:4,1) = MomDK(1:4,1)
!           MomLept(1:4,2) = MomDK(1:4,2)
!           MomLept(1:4,3) = MomDK(1:4,3)
!           MomLept(1:4,4) = MomDK(1:4,4)
!       else
!           MomLept(1:4,1) = MomDK(1:4,1)
!           MomLept(1:4,2) = MomDK(1:4,4)
!           MomLept(1:4,3) = MomDK(1:4,3)
!           MomLept(1:4,4) = MomDK(1:4,2)
!       endif
      
!     MC truth for lepton pairs 
      MomLept(1:4,1:4) = MomDK(1:4,1:4)


!--- compute the invariant mass of the two vector bosons, assuming Z1->l1+l2 and Z2->l3+l4
      mZ1 = sqrt(abs( 2d0*(MomLept(1:4,1).dot.MomLept(1:4,2))  ))
      mZ2 = sqrt(abs( 2d0*(MomLept(1:4,3).dot.MomLept(1:4,4))  ))

      ml12 = mZ1
      ml34 = mZ2
      ml14 = sqrt(abs( 2d0*(MomLept(1:4,1).dot.MomLept(1:4,4))  ))
      ml32 = sqrt(abs( 2d0*(MomLept(1:4,3).dot.MomLept(1:4,2))  ))


!   compute pT of leptons in the lab frame
      pT_lepP = get_PT(MomLept(1:4,2))
      pT_lepM = get_PT(MomLept(1:4,3))

      y_lepP = get_eta(MomLept(1:4,2))
      y_lepM = get_eta(MomLept(1:4,3))



!       if( includeGammaStar ) then ! cut out collinear singularity when intermediate photons are present
!           if( get_R(MomLept(1:4,1),MomLept(1:4,2)).lt.Rsep_ll ) then
!               applyPSCut = .true.
!               return
!           endif
!           if( get_R(MomLept(1:4,1),MomLept(1:4,3)).lt.Rsep_ll ) then
!               applyPSCut = .true.
!               return
!           endif
!           if( get_R(MomLept(1:4,1),MomLept(1:4,4)).lt.Rsep_ll ) then
!               applyPSCut = .true.
!               return
!           endif
!           if( get_R(MomLept(1:4,2),MomLept(1:4,3)).lt.Rsep_ll ) then
!               applyPSCut = .true.
!               return
!           endif
!           if( get_R(MomLept(1:4,2),MomLept(1:4,4)).lt.Rsep_ll ) then
!               applyPSCut = .true.
!               return
!           endif
!           if( get_R(MomLept(1:4,3),MomLept(1:4,4)).lt.Rsep_ll ) then
!               applyPSCut = .true.
!               return
!           endif
!       endif


! construct cos(theta1): angle between direction of fermion from Z1 and negative direction of opposite Z in Z1 rest frame

! this is not fully correct: first, all momenta should be boosted into the resonance rest frame

      MomBoost(1)   = +MomExt(1,3)
      MomBoost(2:4) = -MomExt(2:4,3)

      MomFerm(1:4)  = MomLept(1:4,1)
      call boost(MomFerm(1:4),MomBoost(1:4),mZ1)! boost fermion from Z1 into Z1 rest frame

      MomZ2(1) = MomExt(1,4)
      MomZ2(2:4) = -MomExt(2:4,4)
      call boost(MomZ2(1:4),MomBoost(1:4),mZ1)! boost -Z2 into Z1 rest frame

      CosTheta1 = Get_CosAlpha( MomFerm(1:4),MomZ2(1:4) )






! construct Phi: angle between lepton planes in resonance rest frame
!                equivalent to angle between normal vectors of lepton planes in resonance rest frame

      MomReso(1:4)= MomExt(1:4,3) + MomExt(1:4,4)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
      MomLeptX(1:4,1:4) = MomLept(1:4,1:4)
      ScatteringAxis(1:4) = MomExt(1:4,3)
      call boost(MomLeptX(1:4,1),MomBoost(1:4),MReso)! boost all leptons into the resonance frame
      call boost(MomLeptX(1:4,2),MomBoost(1:4),MReso)
      call boost(MomLeptX(1:4,3),MomBoost(1:4),MReso)
      call boost(MomLeptX(1:4,4),MomBoost(1:4),MReso)
      call boost(ScatteringAxis(1:4),MomBoost(1:4),MReso)


!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt( MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 )! normalize
      
      MomLeptPlane2(2:4) = (MomLeptX(2:4,3)).cross.(MomLeptX(2:4,4))! orthogonal vector to lepton plane
      MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt( MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 )! normalize

!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
      signPhi = sign(1d0,  (dummy(2)*ScatteringAxis(2)+dummy(3)*ScatteringAxis(3)+dummy(4)*ScatteringAxis(4))  )! use q1

      Phi = signPhi * acos(-1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))

!print *, "phi",phi
!pause

! phi(ll)
! Phi = acos((MomLept(2,2)*MomLept(2,3) + MomLept(3,2)*MomLept(3,3))/dsqrt(MomLept(2,2)**2+MomLept(3,2)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2) )




! construct cos(theta*):    scattering angle of Z's in graviton rest frame

      MomReso(1:4)= MomExt(1:4,3) + MomExt(1:4,4)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
      MomZ2(1:4) = MomExt(1:4,4)
!       if( mz1.lt.mz2 ) then 
!           MomZ2(1:4) = MomExt(1:4,3)
!       else
!           MomZ2(1:4) = MomExt(1:4,4)          
!       endif
      call boost(MomZ2(1:4),MomBoost(1:4),MReso)
      CosThetaStar = Get_CosTheta( MomZ2(1:4) )

! print *, MomReso(1:4)
!       MomZ2(1:4) = MomReso(1:4)
!       call boost(MomZ2(1:4),MomBoost(1:4),MReso)
! print *, MomZ2(1:4)
! pause



! construct Phi1:  angle between beam-scattering plane and the lepton plane of Z1 in the resonance rest frame

      BeamAxis(1:4) = (/1d0,0d0,0d0,1d0/)!  energy components are dummies here and will not be used
      ScatteringAxis(1:4) = MomExt(1:4,3)
      MomReso(1:4)= MomExt(1:4,3) + MomExt(1:4,4)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
!       call boost(BeamAxis(1:4),MomBoost(1:4),MReso)
      call boost(ScatteringAxis(1:4),MomBoost(1:4),MReso)


      MomBeamScatterPlane(2:4) = (BeamAxis(2:4)).cross.(ScatteringAxis(2:4))! orthogonal vector to beam-scattering plane
      MomBeamScatterPlane(2:4) = MomBeamScatterPlane(2:4)/dsqrt( MomBeamScatterPlane(2)**2+MomBeamScatterPlane(3)**2+MomBeamScatterPlane(4)**2 ) 

!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomBeamScatterPlane(2:4))
      signPhi1 = sign(1d0,  (dummy(2)*ScatteringAxis(2)+dummy(3)*ScatteringAxis(3)+dummy(4)*ScatteringAxis(4))  )! use q1

      Phi1 = signPhi1 * acos(MomLeptPlane1(2)*MomBeamScatterPlane(2) + MomLeptPlane1(3)*MomBeamScatterPlane(3) + MomLeptPlane1(4)*MomBeamScatterPlane(4))


!     binning
      NBin(1)  = WhichBin(1,pT_lepP)
      NBin(2)  = WhichBin(2,pT_lepM)
      NBin(3)  = WhichBin(3,y_lepP)
      NBin(4)  = WhichBin(4,y_lepM)
      NBin(5)  = WhichBin(5,CosThetaStar)
      NBin(6)  = WhichBin(6,Phi1)
      NBin(7)  = WhichBin(7,CosTheta1)
      NBin(8)  = WhichBin(8,Phi)
if( mz1.gt.mz2 ) then 
      NBin(9)  = WhichBin(9,mZ2)
      NBin(10) = WhichBin(10,mZ1)
else
      NBin(9)  = WhichBin(9,mZ1)
      NBin(10) = WhichBin(10,mZ2)
endif
      NBin(11) = WhichBin(11,mReso)
      NBin(12) = WhichBin(12,ml12)
      NBin(13) = WhichBin(13,ml34)
      NBin(14) = WhichBin(14,ml14)
      NBin(15) = WhichBin(15,ml32)

      NBin(16:18) = 1!  this is for the l+/l- total cross sections

RETURN
END SUBROUTINE





SUBROUTINE Kinematics_HVBF(NumPart,MomExt,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:), mZ1, mZ2, MReso
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NumPart,NBin(:)
real(8) :: m_jj,y_j1,y_j2,dphi_jj,dy_j1j2,pT_jl,pT_j1,pT_j2,pT_H,m_4l
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, Higgs=5


       applyPSCut = .false.
 
       m_jj = get_MInv( MomExt(1:4,outTop)+MomExt(1:4,outBot) )
       y_j1 = get_eta(MomExt(1:4,outTop))
       y_j2 = get_eta(MomExt(1:4,outBot))
       pT_H = get_PT(MomExt(1:4,Higgs))
       pT_j1= get_PT(MomExt(1:4,outTop))
       pT_j2= get_PT(MomExt(1:4,outBot))
       pT_jl = max(pT_j1,pT_j2)
       dy_j1j2 = y_j1 - y_j2
       m_4l = get_MInv(MomExt(1:4,Higgs))

!        if( abs(y_j1).lt.1d0 .or. abs(y_j2).lt.1d0 .or. y_j1*y_j2.gt.0d0 ) then
!           applyPSCut=.true.
!           return
!        endif

       if(  pT_j1.lt.pTjetcut .or. pT_j2.lt.pTjetcut .or. m_jj.lt.mJJcut)  then
          applyPSCut=.true.
          return
       endif


       dphi_jj = abs( Get_PHI(MomExt(1:4,3)) - Get_PHI(MomExt(1:4,4)) )
       if( dphi_jj.gt.Pi ) dphi_jj=2d0*Pi-dphi_jj


 !     binning
       NBin(:) = 1
       NBin(1)  = WhichBin(1,m_jj)
       NBin(2)  = WhichBin(2,dphi_jj)
       NBin(3)  = WhichBin(3,pT_H)
       NBin(4)  = WhichBin(4,pT_jl)
       NBin(5)  = WhichBin(5,dy_j1j2)
       NBin(6)  = WhichBin(6,y_j1)
       NBin(7)  = WhichBin(7,y_j2)
       NBin(8)  = WhichBin(8,m_4l)


RETURN
END SUBROUTINE





SUBROUTINE Kinematics_HVBF_fulldecay(MomExt,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:), mZ1, mZ2, MReso
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NBin(:)
real(8) :: m_jj,y_j1,y_j2,dphi_jj,dy_j1j2,pT_jl,pT_j1,pT_j2,pT_H,m_4l
real(8) :: pT_l1,pT_l2,pT_l3,pT_l4,y_l1,y_l2,y_l3,y_l4
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10


       applyPSCut = .false.
 
       m_jj = get_MInv( MomExt(1:4,outTop)+MomExt(1:4,outBot) )
       m_4l = get_MInv( MomExt(1:4,Lep1P)+MomExt(1:4,Lep1M)+MomExt(1:4,Lep2P)+MomExt(1:4,Lep2M) )
       y_j1 = get_eta(MomExt(1:4,outTop))
       y_j2 = get_eta(MomExt(1:4,outBot))
       pT_H = get_PT(MomExt(1:4,V1)+MomExt(1:4,V2))
       pT_j1= get_PT(MomExt(1:4,outTop))
       pT_j2= get_PT(MomExt(1:4,outBot))
       pT_l1= get_PT(MomExt(1:4,Lep1P))
       pT_l2= get_PT(MomExt(1:4,Lep1M))
       pT_l3= get_PT(MomExt(1:4,Lep2P))
       pT_l4= get_PT(MomExt(1:4,Lep2M))
       y_l1= get_eta(MomExt(1:4,Lep1P))
       y_l2= get_eta(MomExt(1:4,Lep1M))
       y_l3= get_eta(MomExt(1:4,Lep2P))
       y_l4= get_eta(MomExt(1:4,Lep2M))
       pT_jl = max(pT_j1,pT_j2)
       dy_j1j2 = y_j1 - y_j2
       
       mZ1 = get_MInv(MomExt(1:4,Lep1P)+MomExt(1:4,Lep1M))
       mZ2 = get_MInv(MomExt(1:4,Lep2P)+MomExt(1:4,Lep2M))


       if( mZ1.lt.2.5d0*GeV ) then
          applyPSCut=.true.
          return
       endif
       if( mZ2.lt.2.5d0*GeV ) then
          applyPSCut=.true.
          return
       endif
       if( m_4l.lt.70d0*GeV ) then
          applyPSCut=.true.
          return
       endif
       
       if( pT_l1.lt.3d0*GeV .or. pT_l2.lt.3d0*GeV .or. pT_l3.lt.3d0*GeV .or. pT_l4.lt.3d0*GeV ) then
          applyPSCut=.true.
          return
       endif
       
       if( dabs(y_l1).gt.2.7d0 .or. dabs(y_l2).gt.2.7d0 .or. dabs(y_l3).gt.2.7d0 .or. dabs(y_l4).gt.2.7d0 ) then
          applyPSCut=.true.
          return
       endif
       

!      VERY loose VBF cuts
!        if( m_jj.lt.400d0*GeV ) then
!           applyPSCut=.true.
!           return
!        endif
!        if( dabs(m_4l-m_reso).gt.20d0*GeV ) then
!           applyPSCut=.true.
!           return
!        endif

       if( dabs(m_4l).lt.70d0*GeV ) then
          applyPSCut=.true.
          return
       endif

!        if( abs(y_j1).gt.5d0 .or. abs(y_j2).gt.5d0 ) then
!           applyPSCut=.true.
!           return
!        endif
! 
!        if( abs(y_j1-y_j2).lt.2.0d0 .or. y_j1*y_j2.gt.0d0 ) then
!           applyPSCut=.true.
!           return
!        endif

!        if(  pT_j1.lt.pTjetcut .or. pT_j2.lt.pTjetcut .or. m_jj.lt.mJJcut )  then
!           applyPSCut=.true.
!           return
!        endif

       dphi_jj = abs( Get_PHI(MomExt(1:4,3)) - Get_PHI(MomExt(1:4,4)) )
       if( dphi_jj.gt.Pi ) dphi_jj=2d0*Pi-dphi_jj


 !     binning
       NBin(:) = 1
       NBin(1)  = WhichBin(1,m_jj)
       NBin(2)  = WhichBin(2,dphi_jj)
       NBin(3)  = WhichBin(3,pT_H)
       NBin(4)  = WhichBin(4,pT_jl)
       NBin(5)  = WhichBin(5,dy_j1j2)
       NBin(6)  = WhichBin(6,y_j1)
       NBin(7)  = WhichBin(7,y_j2)
       NBin(8)  = WhichBin(8,m_4l)


RETURN
END SUBROUTINE







SUBROUTINE Kinematics_HJ(NumPart,MomExt,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:)!,mZ1, mZ2, MReso
!real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
!real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NBin(:),NumPart!NumHadr,JetList(1:2),NJEt,PartList(1:2)
!real(8) :: pT_lepM,pT_lepP,y_lepM,y_lepP,MomFerm(1:4),MomZ2(1:4),MomReso(1:4),CosTheta1,Phi,Phi1,signPhi,signPhi1
!real(8) :: CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,CosThetaZ,CosThetaStar
real(8) :: MomHadr(1:4,1:2),MomJet(1:4,1:2),pT_j1!,pT_j2,pT_3, pT_4,m_jj,dphi_jj



applyPSCut = .false.


!    NJet=0
    MomHadr(1:4,1)=MomExt(1:4,4)
!    NumHadr=1
!    PartList=(/1,2/)
!    JetList(1:2) = (/1,2/)
    MomJet(1:4,1) = MomHadr(1:4,1)
!    call JetAlgo_kt(Rjet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
!    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order

!    if( Njet.ne.2 ) then
!      applyPSCut=.true.
!      return
!    endif

    pT_j1 = get_PT(MomJet(1:4,1))

    if( pT_j1.lt.ptjetcut) then
    !print *, "cut",pT_j1
      applyPSCut=.true.
      return
    endif

 !     binning
       NBin(1)  = WhichBin(1,pT_j1)


RETURN
END SUBROUTINE









SUBROUTINE Kinematics_HJJ(NumPart,MomExt,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:),mZ1, mZ2, MReso
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NumPart,NBin(:),NumHadr,JetList(1:2),NJEt,PartList(1:2)
real(8) :: pT_lepM,pT_lepP,y_lepM,y_lepP,MomFerm(1:4),MomZ2(1:4),MomReso(1:4),CosTheta1,Phi,Phi1,signPhi,signPhi1
real(8) :: CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,CosThetaZ,CosThetaStar,y_j1,y_j2,dy_j1j2
real(8) :: pT_3, pT_4,MomHadr(1:4,1:2),MomJet(1:4,1:2),pT_j1,pT_j2,m_jj,dphi_jj,pT_H



applyPSCut = .false.


    NJet=0
    MomHadr(1:4,1:2)=MomExt(1:4,3:4)
    NumHadr=2
    PartList=(/1,2/)
    JetList(1:2) = (/1,2/)
    MomJet(1:4,1:2) = MomHadr(1:4,1:2)
    call JetAlgo_kt(Rjet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order


    if( Njet.ne.2 ) then
      applyPSCut=.true.
      return
    endif

    pT_j1 = get_PT(momjet(1:4,1))
    pT_j2 = get_PT(momjet(1:4,2))

    pT_H = get_PT(MomExt(1:4,5))
    
    y_j1 = get_ETA(momjet(1:4,1))
    y_j2 = get_ETA(momjet(1:4,2))

   dy_j1j2 = y_j1 - y_j2


  
    m_jj = get_MInv( MomJet(1:4,1)+MomJet(1:4,2) )
  
    if( pT_j1.lt.ptjetcut .or. pT_j2.lt.ptjetcut .or. m_jj.lt.mJJcut) then
      applyPSCut=.true.
      return
    endif

    dphi_jj = abs( Get_PHI(MomJet(1:4,1)) - Get_PHI(MomJet(1:4,2)) )
    if( dphi_jj.gt.Pi ) dphi_jj=2d0*Pi-dphi_jj


 !     binning
       NBin(:) = 1
       NBin(1)  = WhichBin(1,m_jj)
       NBin(2)  = WhichBin(2,dphi_jj)
       NBin(3)  = WhichBin(3,pT_H)
       NBin(4)  = WhichBin(4,y_j1)
       NBin(5)  = WhichBin(5,y_j2)
       NBin(6)  = WhichBin(6,dy_j1j2)

       
RETURN
END SUBROUTINE



SUBROUTINE Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut)
use ModMisc
use ModParameters
implicit none

logical :: applyPSCut
integer :: NumPart,NBin(:),id(:)
real(8) :: m_jj,y_j1,y_j2,dphi_jj, m_ll, pt_V, pt_H, m_Vstar, costheta1, costheta2, phistar1, phi
double precision MomBoost(1:4), MomFerm(1:4), inv_mass(1:9), MomLeptX(1:4,1:4), ScatteringAxis(1:4), MomReso(1:4)
double precision MomLeptPlane1(2:4), MomLeptPlane2(2:4), dummy(2:4), signPhi
double precision, intent(in) :: MomExt(1:4,1:9) !,beam_momentum(2,4),four_momentum(7,4)


     applyPSCut = .false.
     m_jj = inv_mass(5)
     m_ll = inv_mass(4)

     if(inv_mass(4).le.getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))then
      applyPSCut=.true.
     endif
     if(inv_mass(5).le.getMass(convertLHEreverse(id(8)))+getMass(convertLHEreverse(id(9))))then
      applyPSCut=.true.
     endif

     pt_V = get_PT(MomExt(1:4,6)+MomExt(1:4,7))
     pt_H = get_PT(MomExt(1:4,8)+MomExt(1:4,9))
     !inv_mass(4,5,6,7)
     m_Vstar = get_MInv(MomExt(1:4,6)+MomExt(1:4,7)+MomExt(1:4,8)+MomExt(1:4,9))
!     y_j1 = get_eta(MomExt(1:4,3))
!     y_j2 = get_eta(MomExt(1:4,4))

!        if( abs(y_j1).lt.1d0 .or. abs(y_j2).lt.1d0 .or. y_j1*y_j2.gt.0d0 ) then
!           applyPSCut=.true.
!           return
!        endif

!     dphi_jj = abs( Get_PHI(MomExt(1:4,3)) - Get_PHI(MomExt(1:4,4)) )
!     if( dphi_jj.gt.Pi ) dphi_jj=2d0*Pi-dphi_jj

!costheta2 - Z decay angle
      MomBoost(1)   = +MomExt(1,4)
      MomBoost(2:4) = -MomExt(2:4,4)

      MomFerm(1:4)  = MomExt(1:4,6)
      call boost(MomFerm(1:4),MomBoost(1:4),inv_mass(4))! boost fermion from Z1 into Z1 rest frame

      costheta2 = Get_CosAlpha( MomFerm(1:4),MomExt(1:4,4) )

     
!costheta1 - production angle
! no boost for now
      MomBoost(1)   = +MomExt(1,3)
      MomBoost(2:4) = -MomExt(2:4,3)

      MomFerm(1:4)  = MomExt(1:4,4)
      !call boost(MomFerm(1:4),MomBoost(1:4),inv_mass(1))! boost fermion from Z1 into Z1 rest frame
      call LORENTZ(MomFerm(1:4),MomBoost(1:4))

      !if(beam_momentum(1,1).gt.beam_momentum(2,1))then
      !  costheta1 = Get_CosAlpha( MomFerm(1:4),beam_momentum(1,1:4) )
      !else
      !  costheta1 = Get_CosAlpha( MomFerm(1:4),beam_momentum(2,1:4) )
      !endif

      ScatteringAxis(1:4)=(/1,0,0,1/)
      costheta1 = Get_CosAlpha( MomFerm(1:4), ScatteringAxis(1:4) )

      !if(beam_momentum(1,4).lt.0d0)then
       ! print *, '!!!'
      !endif


     !costheta1 = -(four_momentum(2,2)*beam_momentum(1,2)+four_momentum(2,3)*beam_momentum(1,3)+four_momentum(2,4)*beam_momentum(1,4)) &
     !         / dsqrt((four_momentum(2,2)**2+four_momentum(2,3)**2+four_momentum(2,4)**2) &
     !               *(beam_momentum(1,2)**2+beam_momentum(1,3)**2+beam_momentum(1,4)**2))


!phistar1
      MomReso(1:4)= MomExt(1:4,5)!MomExt(1:4,3) + MomExt(1:4,4)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
      MomLeptX(1:4,1) = MomExt(1:4,6)
      MomLeptX(1:4,2) = MomExt(1:4,7)
      MomLeptX(1:4,3) = MomExt(1:4,8)
      MomLeptX(1:4,4) = MomExt(1:4,9)
      ScatteringAxis(1:4) = MomExt(1:4,4)
      call boost(MomLeptX(1:4,1),MomBoost(1:4),inv_mass(5))! boost all leptons into the resonance frame
      call boost(MomLeptX(1:4,2),MomBoost(1:4),inv_mass(5))
      call boost(MomLeptX(1:4,3),MomBoost(1:4),inv_mass(5))
      call boost(MomLeptX(1:4,4),MomBoost(1:4),inv_mass(5))
      call boost(ScatteringAxis(1:4),MomBoost(1:4),inv_mass(5))
!test
!print *, MomLeptX(1:4,3)+MomLeptX(1:4,4)
!print *, four_momentum(2,1:4)
!print *, four_momentum(3,1:4)
!pause


!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt(dabs(MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15) )! normalize
      
      MomLeptPlane2(2:4) = (MomLeptX(2:4,3)).cross.(MomLeptX(2:4,4))! orthogonal vector to lepton plane
      MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt(dabs(MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15 ))! normalize

!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
      signPhi = sign(1d0,  (dummy(2)*ScatteringAxis(2)+dummy(3)*ScatteringAxis(3)+dummy(4)*ScatteringAxis(4))  )! use q1
    !test
    !signPhi = 1d0
      Phistar1 = signPhi * acos(-1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))


!phi
      !MomReso(1:4)= MomExt(1:4,5)!MomExt(1:4,8) + MomExt(1:4,9)
      !MomBoost(1)   = +MomReso(1)
      !MomBoost(2:4) = -MomReso(2:4)
      MomLeptX(1:4,1) = MomExt(1:4,6)
      MomLeptX(1:4,2) = MomExt(1:4,7)
      !MomLeptX(1:4,3) = MomExt(1:4,1)
      !MomLeptX(1:4,4) = MomExt(1:4,2)

      !if(beam_momentum(1,1).gt.beam_momentum(2,1))then
        !MomLeptX(1:4,3) = beam_momentum(1,1:4)
      !else
      !  MomLeptX(1:4,3) = beam_momentum(2,1:4)
      !endif


      !ScatteringAxis(1:4) = MomExt(1:4,5)
      ScatteringAxis(1:4)=(/1,0,0,1/)
      !call boost(MomLeptX(1:4,1),MomBoost(1:4),inv_mass(5))! boost all leptons into the resonance frame
      !call boost(MomLeptX(1:4,2),MomBoost(1:4),inv_mass(5))
      !call boost(MomLeptX(1:4,3),MomBoost(1:4),inv_mass(5))
      !call boost(MomLeptX(1:4,4),MomBoost(1:4),inv_mass(5))
      !call boost(ScatteringAxis(1:4),MomBoost(1:4),inv_mass(5))
!test
!print *, MomLeptX(1:4,3)+MomLeptX(1:4,4)
!print *, four_momentum(2,1:4)
!print *, four_momentum(3,1:4)
!pause


!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt( MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 )! normalize
      
      MomLeptPlane2(2:4) = (ScatteringAxis(2:4)).cross.(MomExt(2:4,4))! orthogonal vector to lepton plane
      MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt( MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 )! normalize
      !MomLeptPlane2(2:4) = (/0d0,1d0,0d0/)

!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
      signPhi = -sign(1d0,  (dummy(2)*ScatteringAxis(2)+dummy(3)*ScatteringAxis(3)+dummy(4)*ScatteringAxis(4))  )! use q1
    !test
    !signPhi = 1d0
      phi = signPhi * acos(1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))


!     binning
     NBin(1)  = WhichBin(1,m_jj)
     NBin(2)  = WhichBin(2,m_ll)
     NBin(3)  = WhichBin(3,pt_V)
     NBin(4)  = WhichBin(4,pt_H)
     NBin(5)  = WhichBin(5,m_Vstar)
     Nbin(6)  = WhichBin(6,costheta1)
     Nbin(7)  = WhichBin(7,costheta2)
     Nbin(8)  = WhichBin(8,phistar1)
     Nbin(9)  = WhichBin(9,phi)

RETURN
END SUBROUTINE Kinematics_VHiggs




SUBROUTINE Kinematics_TTBH(Mom,applyPSCut,NBin)
use ModParameters
use ModMisc
! use modTTBH
implicit none
real(8) :: Mom(1:4,1:13),MomMELA(1:4,1:13)
logical :: applyPSCut
integer :: NBin(:)
real(8) :: pT_t,pT_H,pT_tbar,MatElSq_H0,MatElSq_H1,D_0minus
real(8) :: mt,mtbar,mttbar,mWp,mWm,pT_b,pT_l,pT_miss
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13
logical,save :: FirstTime=.true.


    applyPSCut = .false.

    pT_t = get_PT(Mom(1:4,t))
    pT_tbar = get_PT(Mom(1:4,tbar))
    pT_H = get_PT(Mom(1:4,Hbos))
    pT_b = get_PT(Mom(1:4,b))
    pT_l = get_PT(Mom(1:4,LepP))
    pT_miss = get_PT(Mom(1:4,nu)+Mom(1:4,nubar))
    mt = get_MInv(Mom(1:4,t))
    mtbar = get_MInv(Mom(1:4,tbar))
    mttbar = get_MInv(Mom(1:4,t)+Mom(1:4,tbar))
    mWp = get_MInv(Mom(1:4,Wp))
    mWm = get_MInv(Mom(1:4,Wm))
    
    if( m_Top.lt.10d0*GeV  .and. (pT_t.lt.pTjetcut .or. pT_tbar.lt.pTjetcut .or. mttbar.lt.mJJcut) ) applyPSCut=.true.

    
!     if( FirstTime ) then
! !       call NNPDFDriver("./pdfs/NNPDF30_lo_as_0130.LHgrid",33)
! !       call NNinitPDF(0)
!       call InitProcess_TTBH(m_Reso,m_top)
!       FirstTime = .false.
!     endif
!     MomMELA(1:4,1) = -(/         65d0,           0.0000000000000000d0, 0.0000000000000000d0,      65d0           /)
!     MomMELA(1:4,2) = -(/         65d0,           0.0000000000000000d0, 0.0000000000000000d0,     -65d0           /)  
!     MomMELA(1:4,3:11) = Mom(1:4,3:11)  remap here because of new W bosons
!     MomMELA(1:4,12:13) = 0d0
!     
!     call EvalXSec_PP_TTBH(MomMELA(1:4,1:13),(/(1d0,0d0),(0d0,0d0)/),TopDecays,2,MatElSq_H0)
!     call EvalXSec_PP_TTBH(MomMELA(1:4,1:13),(/(0d0,0d0),(1d0,0d0)/),TopDecays,2,MatElSq_H1)
!     
!     D_0minus = MatElSq_H0/(MatElSq_H0 + 2d0*MatElSq_H1 )

    D_0minus=0d0

!   binning
    NBin(1)  = WhichBin(1,pT_t)
    NBin(2)  = WhichBin(2,pT_H)
    NBin(3)  = WhichBin(3,mt)
    NBin(4)  = WhichBin(4,mtbar)
    NBin(5)  = WhichBin(5,mWp)
    NBin(6)  = WhichBin(6,mWm)
    NBin(7)  = WhichBin(7,pT_b)
    NBin(8)  = WhichBin(8,pT_l)
    NBin(9)  = WhichBin(9,pT_miss)
    NBin(10) = WhichBin(10,D_0minus)

    
    
RETURN
END SUBROUTINE



SUBROUTINE Kinematics_BBBH(Mom,applyPSCut,NBin)
use ModParameters
use ModMisc
implicit none
real(8) :: Mom(1:4,1:11)
logical :: applyPSCut
integer :: NBin(:)
real(8) :: pT_b,pT_H,pT_bbar,mbbbar
integer, parameter :: bbar=4,b=5,Hbos=3,inLeft=1,inRight=2


    applyPSCut = .false.

    pT_b = get_PT(Mom(1:4,b))
    pT_bbar = get_PT(Mom(1:4,bbar))
    pT_H = get_PT(Mom(1:4,Hbos))
    mbbbar = get_MInv(Mom(1:4,b)+Mom(1:4,bbar))
    
    if( pT_b.lt.pTjetcut .or. pT_bbar.lt.pTjetcut .or. mbbbar.lt.mJJcut) applyPSCut=.true.
    
!   binning
    NBin(1)  = WhichBin(1,pT_b)
    NBin(2)  = WhichBin(2,pT_H)
    
    
RETURN
END SUBROUTINE





SUBROUTINE Kinematics_TH(Mom,applyPSCut,NBin)
use ModParameters
use ModMisc
! use ModTH
implicit none
real(8) :: Mom(1:4,1:9),MomMELA(1:4,1:9)
logical :: applyPSCut
integer :: NBin(:)
real(8) :: pT_Top,pT_Higgs,pT_j,eta_j,eta_top,eta_Higgs,y_Higgs,y_top,y_j,MatElSq_H0,MatElSq_H1,D_0minus
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4,ljet=5,bbar=6, lepP=8,nu=9
logical,save :: FirstTime=.true.


    applyPSCut = .false.


       pT_Top  = get_PT(Mom(1:4,t))
       pT_Higgs= get_PT(Mom(1:4,Hbos))
       pT_j= get_PT(Mom(1:4,ljet))

       y_top=get_eta(Mom(1:4,t))
       y_Higgs=get_eta(Mom(1:4,Hbos))
       y_j=get_eta(Mom(1:4,ljet))

    
    if( m_Top.lt.10d0*GeV  .and. (pT_top.lt.pTjetcut) ) applyPSCut=.true.
    
   
   
D_0minus = 0d0
! !     if( FirstTime ) then
! ! !       call NNPDFDriver("./pdfs/NNPDF30_lo_as_0130.LHgrid",33)
! ! !       call NNinitPDF(0)
! !       call InitProcess_TTBH(m_Reso,m_top)
! !       FirstTime = .false.
! !     endif
!     MomMELA(1:4,1) = -(/         65d0,           0.0000000000000000d0, 0.0000000000000000d0,      65d0           /)
!     MomMELA(1:4,2) = -(/         65d0,           0.0000000000000000d0, 0.0000000000000000d0,     -65d0           /)  
!     MomMELA(1:4,3:9) = Mom(1:4,3:9)
!     
!     call EvalXSec_PP_TH(MomMELA,(/(1d0,0d0),(0d0,0d0)/),TopDecays,MatElSq_H0)
!     call EvalXSec_PP_TH(MomMELA,(/(0d0,0d0),(1d0,0d0)/),TopDecays,MatElSq_H1)
!     
!     MatElSq_H0 = MatElSq_H0 !/40d0
!     MatElSq_H1 = MatElSq_H1 !/147d0
        
!     D_0minus = MatElSq_H0/(MatElSq_H0 + 1d0*MatElSq_H1 )
   
   
    
!   binning
                                                                             
       NBin(1) = WhichBin(1,pT_Top)
       NBin(2) = WhichBin(2,y_Top)
       NBin(3) = WhichBin(3,pT_j)
       NBin(4) = WhichBin(4,y_j)
       NBin(5) = WhichBin(5,pT_Higgs)
       NBin(6) = WhichBin(6,y_Higgs)
       NBin(7) = WhichBin(7,D_0minus)
    
    
    
RETURN
END SUBROUTINE





SUBROUTINE Kinematics_Htautau(Mom,applyPSCut,NBin)
use ModParameters
use ModMisc
implicit none
real(8) :: Mom(:,:)
logical :: applyPSCut
integer :: NBin(:)
real(8) :: m_tauP,m_tauM,m_Wp,m_Wm,y_tauM,y_tauP
integer, parameter :: inLeft=1, inRight=2, Hig=3, tauP=4, tauM=5, Wp=6, Wm=7,   nu=8, nubar_tau=9, lepP=10,   lepM=11, nu_tau=12, nubar=13


    applyPSCut = .false.
    
    m_tauP = get_MInv(Mom(1:4,tauP))
    m_tauM = get_MInv(Mom(1:4,tauM))
    m_Wp   = get_MInv(Mom(1:4,Wp))
    m_Wm   = get_MInv(Mom(1:4,Wm))
    y_tauP = get_ETA(Mom(1:4,tauP))  
    y_tauM = get_ETA(Mom(1:4,tauM))  
    
    
    
!   binning
    NBin(:)=1
    NBin(1) = WhichBin(1,m_tauP)
    NBin(2) = WhichBin(2,m_tauM)
    NBin(3) = WhichBin(3,m_Wp)
    NBin(4) = WhichBin(4,m_Wm)
    NBin(5) = WhichBin(5,y_tauP)
    NBin(6) = WhichBin(6,y_tauM)
    
RETURN
END SUBROUTINE
FUNCTION ZLepBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZLepBranching


  if( xRnd .le. Brlept_Z_ee/(100d0*percent-Brlept_Z_tt) ) then
      ZLepBranching = ElM_
  elseif(xRnd .le. (Brlept_Z_ee+Brlept_Z_mm)/(100d0*percent-Brlept_Z_tt) ) then
      ZLepBranching = MuM_
  else
      print *, "error ",xRnd
      stop
  endif

!print *, "checker 2",(Brlept_Z_ee)/(100d0*percent-Brlept_Z_tt)
!print *, "checker 2",(Brlept_Z_ee+Brlept_Z_mm)/(100d0*percent-Brlept_Z_tt)

RETURN
END FUNCTION




FUNCTION ZLepPlusTauBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZLepPlusTauBranching

  if( xRnd .le. Brlept_Z_ee ) then
      ZLepPlusTauBranching = ElM_
  elseif(xRnd .le. Brlept_Z_ee+Brlept_Z_mm ) then
      ZLepPlusTauBranching = MuM_
  elseif(xRnd .le. Brlept_Z_ee+Brlept_Z_mm+Brlept_Z_tt ) then
      ZLepPlusTauBranching = TaM_
  else
      print *, "error ",xRnd
      stop
  endif

! print *, "checker 3",Brlept_Z_ee
! print *, "checker 3",Brlept_Z_ee+Brlept_Z_mm
! print *, "checker 3",Brlept_Z_ee+Brlept_Z_mm+Brlept_Z_tt

RETURN
END FUNCTION



FUNCTION ZNuBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZNuBranching

  if( xRnd .le. Brlept_Z_nn ) then
      ZNuBranching = NuE_
  elseif(xRnd .le. Brlept_Z_nn+Brlept_Z_nn ) then
      ZNuBranching = NuM_
  elseif(xRnd .le. Brlept_Z_nn+Brlept_Z_nn+Brlept_Z_nn ) then
      ZNuBranching = NuT_
  else
      print *, "error ",xRnd
      stop
  endif

!print *, "checker 4",Brlept_Z_nn
!print *, "checker 4",Brlept_Z_nn+Brlept_Z_nn
!print *, "checker 4",Brlept_Z_nn+Brlept_Z_nn+Brlept_Z_nn

RETURN
END FUNCTION





FUNCTION ZQuaBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZQuaBranching_flat
real(8),parameter :: Ncol=3d0
real(8),parameter :: xxxx=1d0/15d0
real(8),parameter :: yyyy=Ncol*xxxx

  if( xRnd .le. yyyy ) then
      ZQuaBranching_flat = Up_
  elseif(xRnd .le. yyyy+yyyy) then
      ZQuaBranching_flat = Chm_
  elseif(xRnd .le. yyyy+yyyy+yyyy) then
      ZQuaBranching_flat = Dn_
  elseif(xRnd .le. yyyy+yyyy+yyyy+yyyy) then
      ZQuaBranching_flat = Str_
  elseif(xRnd .le. yyyy+yyyy+yyyy+yyyy+yyyy) then
      ZQuaBranching_flat = Bot_
  else
      print *, "error ",xRnd
      stop
  endif

!print *, "checker 1",Brhadr_Z_uu,Brhadr_Z_dd
!print *, "checker 1",Brhadr_Z_uu+Brhadr_Z_cc+Brhadr_Z_dd+Brhadr_Z_ss+Brhadr_Z_bb

RETURN
END FUNCTION



FUNCTION ZQuaBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZQuaBranching
real(8),parameter :: Ncol=3d0
real(8),parameter :: xxxx=1d0/15d0
real(8),parameter :: yyyy=Ncol*xxxx

  if( xRnd .le. Brhadr_Z_uu ) then
      ZQuaBranching = Up_
  elseif(xRnd .le. Brhadr_Z_uu+Brhadr_Z_cc) then
      ZQuaBranching = Chm_
  elseif(xRnd .le. Brhadr_Z_uu+Brhadr_Z_cc+Brhadr_Z_dd) then
      ZQuaBranching = Dn_
  elseif(xRnd .le. Brhadr_Z_uu+Brhadr_Z_cc+Brhadr_Z_dd+Brhadr_Z_ss) then
      ZQuaBranching = Str_
  elseif(xRnd .le. Brhadr_Z_uu+Brhadr_Z_cc+Brhadr_Z_dd+Brhadr_Z_ss+Brhadr_Z_bb) then
      ZQuaBranching = Bot_
  else
      print *, "error ",xRnd
      stop
  endif
  

RETURN
END FUNCTION








FUNCTION ZAnyBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZAnyBranching_flat
real(8),parameter :: Ncol=3d0
real(8),parameter :: xx=1d0/21d0
real(8),parameter :: yy=Ncol*xx

! real(8),parameter :: xx=1d0/11d0
! real(8),parameter :: yy=xx


  if( xRnd .le. yy ) then
      ZAnyBranching_flat = Up_
  elseif(xRnd .le. 2*yy ) then
      ZAnyBranching_flat = Chm_
  elseif(xRnd .le. 2*yy+yy ) then
      ZAnyBranching_flat = Dn_
  elseif(xRnd .le. 2*yy+2*yy ) then
      ZAnyBranching_flat = Str_
  elseif(xRnd .le. 2*yy+3*yy ) then
      ZAnyBranching_flat = Bot_
  elseif(xRnd .le. yy*(2+3) + xx ) then
      ZAnyBranching_flat = ElM_
  elseif(xRnd .le. yy*(2+3) + xx*2 ) then
      ZAnyBranching_flat = MuM_
  elseif(xRnd .le. yy*(2+3) + xx*3 ) then
      ZAnyBranching_flat = TaM_
  elseif(xRnd .le. yy*(2+3)+xx*3 + xx ) then
      ZAnyBranching_flat = NuE_
  elseif(xRnd .le. yy*(2+3)+xx*3 + xx*2 ) then
      ZAnyBranching_flat = NuM_
  elseif(xRnd .le. yy*(2+3)+xx*3 + xx*3 ) then
      ZAnyBranching_flat = NuT_
  else
      print *, "error ",xRnd
      stop
  endif

  
 
  

!   if( xRnd .le. yy*scale_alpha_Z_uu ) then
!       ZAnyBranching = Up_
!   elseif(xRnd .le. 2*yy*scale_alpha_Z_uu ) then
!       ZAnyBranching = Chm_
!   elseif(xRnd .le. 2*yy*scale_alpha_Z_uu+yy*scale_alpha_Z_dd ) then
!       ZAnyBranching = Dn_
!   elseif(xRnd .le. 2*yy*scale_alpha_Z_uu+2*yy*scale_alpha_Z_dd ) then
!       ZAnyBranching = Str_
!   elseif(xRnd .le. 2*yy*scale_alpha_Z_uu+3*yy*scale_alpha_Z_dd ) then
!       ZAnyBranching = Bot_
!   elseif(xRnd .le. yy*(2*scale_alpha_Z_uu+3*scale_alpha_Z_dd) + xx*scale_alpha_Z_ll ) then
!       ZAnyBranching = ElM_
!   elseif(xRnd .le. yy*(2*scale_alpha_Z_uu+3*scale_alpha_Z_dd) + xx*2*scale_alpha_Z_ll ) then
!       ZAnyBranching = MuM_
!   elseif(xRnd .le. yy*(2*scale_alpha_Z_uu+3*scale_alpha_Z_dd) + xx*3*scale_alpha_Z_ll ) then
!       ZAnyBranching = TaM_
!   elseif(xRnd .le. yy*(2*scale_alpha_Z_uu+3*scale_alpha_Z_dd)+xx*3*scale_alpha_Z_ll + xx*scale_alpha_Z_nn ) then
!       ZAnyBranching = NuE_
!   elseif(xRnd .le. yy*(2*scale_alpha_Z_uu+3*scale_alpha_Z_dd)+xx*3*scale_alpha_Z_ll + xx*2*scale_alpha_Z_nn ) then
!       ZAnyBranching = NuM_
!   elseif(xRnd .le. yy*(2*scale_alpha_Z_uu+3*scale_alpha_Z_dd)+xx*3*scale_alpha_Z_ll + xx*3*scale_alpha_Z_nn ) then
!       ZAnyBranching = NuT_
!   else
!       print *, "error ",xRnd
!       stop
!   endif



RETURN
END FUNCTION




FUNCTION ZAnyBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZAnyBranching


  if( xRnd .le. Br_Z_uu ) then
      ZAnyBranching = Up_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc) then
      ZAnyBranching = Chm_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd) then
      ZAnyBranching = Dn_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss) then
      ZAnyBranching = Str_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb) then
      ZAnyBranching = Bot_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee) then
      ZAnyBranching = ElM_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm) then
      ZAnyBranching = MuM_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt) then
      ZAnyBranching = TaM_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt+Br_Z_nn) then
      ZAnyBranching = NuE_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt+Br_Z_nn+Br_Z_nn) then
      ZAnyBranching = NuM_
  elseif(xRnd .le. Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt+Br_Z_nn+Br_Z_nn+Br_Z_nn) then
      ZAnyBranching = NuT_
  else
      print *, "error ",xRnd
      stop
  endif


!   if( xRnd .le. scale_alpha_Z_uu*Br_Z_uu ) then
!       ZAnyBranching = Up_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc) then
!       ZAnyBranching = Chm_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd) then
!       ZAnyBranching = Dn_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss) then
!       ZAnyBranching = Str_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb) then
!       ZAnyBranching = Bot_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee) then
!       ZAnyBranching = ElM_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee+scale_alpha_Z_ll*Br_Z_mm) then
!       ZAnyBranching = MuM_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee+scale_alpha_Z_ll*Br_Z_mm+scale_alpha_Z_ll*Br_Z_tt) then
!       ZAnyBranching = TaM_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee+scale_alpha_Z_ll*Br_Z_mm+scale_alpha_Z_ll*Br_Z_tt+scale_alpha_Z_nn*Br_Z_nn) then
!       ZAnyBranching = NuE_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee+scale_alpha_Z_ll*Br_Z_mm+scale_alpha_Z_ll*Br_Z_tt+scale_alpha_Z_nn*Br_Z_nn+scale_alpha_Z_nn*Br_Z_nn) then
!       ZAnyBranching = NuM_
!   elseif(xRnd .le. scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee+scale_alpha_Z_ll*Br_Z_mm+scale_alpha_Z_ll*Br_Z_tt+scale_alpha_Z_nn*Br_Z_nn+scale_alpha_Z_nn*Br_Z_nn+scale_alpha_Z_nn*Br_Z_nn) then
!       ZAnyBranching = NuT_
!   else
!       print *, "error ",xRnd
!       stop
!   endif

   
! print *, "checker ",Br_Z_uu+Br_Z_cc+Br_Z_dd+Br_Z_ss+Br_Z_bb+Br_Z_ee+Br_Z_mm+Br_Z_tt+Br_Z_nn+Br_Z_nn+Br_Z_nn
! print *, "checker ",scale_alpha_Z_uu*Br_Z_uu+scale_alpha_Z_uu*Br_Z_cc+scale_alpha_Z_dd*Br_Z_dd+scale_alpha_Z_dd*Br_Z_ss+scale_alpha_Z_dd*Br_Z_bb+scale_alpha_Z_ll*Br_Z_ee+scale_alpha_Z_ll*Br_Z_mm+scale_alpha_Z_ll*Br_Z_tt+scale_alpha_Z_nn*Br_Z_nn+scale_alpha_Z_nn*Br_Z_nn+scale_alpha_Z_nn*Br_Z_nn
! pause
  
  
RETURN
END FUNCTION






FUNCTION WLepBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WLepBranching

  if( xRnd .le. Brlept_W_en /(100d0*percent-Brlept_W_tn) ) then
      WLepBranching = ElM_
  elseif(xRnd .le. (Brlept_W_en+Brlept_W_mn)/(100d0*percent-Brlept_W_tn) ) then
      WLepBranching = MuM_
  else
      print *, "error ",xRnd
      stop
  endif

!print *, "checker 6",Brlept_W_en /(100d0*percent-Brlept_W_tn)
!print *, "checker 6",(Brlept_W_en+Brlept_W_mn)/(100d0*percent-Brlept_W_tn)

RETURN
END FUNCTION



FUNCTION WLepPlusTauBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WLepPlusTauBranching

  if( xRnd .le. Brlept_W_en ) then
      WLepPlusTauBranching = ElM_
  elseif(xRnd .le. Brlept_W_en+Brlept_W_mn ) then
      WLepPlusTauBranching = MuM_
  elseif(xRnd .le. Brlept_W_en+Brlept_W_mn+Brlept_W_tn ) then
      WLepPlusTauBranching = TaM_
  else
      print *, "error ",xRnd
      stop
  endif

! print *, "checker 7",Brlept_W_en
! print *, "checker 7",Brlept_W_en+Brlept_W_mn
! print *, "checker 7",Brlept_W_en+Brlept_W_mn+Brlept_W_tn

RETURN
END FUNCTION




FUNCTION WQuaUpBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WQuaUpBranching

  if( xRnd .le. Brhadr_W_ud ) then
      WQuaUpBranching = Up_
  elseif(xRnd .le. Brhadr_W_ud+Brhadr_W_cs ) then
      WQuaUpBranching = Chm_
  else
      print *, "error ",xRnd
      stop
  endif


!print *, "checker 8",Brhadr_W_ud
!print *, "checker 8",Brhadr_W_ud+Brhadr_W_cs

RETURN
END FUNCTION




FUNCTION WAnyBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WAnyBranching_flat
real(8),parameter :: Ncol=3d0
real(8),parameter :: xx=1d0/9d0
real(8),parameter :: yy=Ncol*xx


  if( xRnd .le. yy ) then
      WAnyBranching_flat = Up_
  elseif(xRnd .le. 2*yy ) then
      WAnyBranching_flat = Chm_
  elseif(xRnd .le. 2*yy+xx ) then
      WAnyBranching_flat = ElM_
  elseif(xRnd .le. 2*yy+2*xx ) then
      WAnyBranching_flat = MuM_
  elseif(xRnd .le. 2*yy+3*xx ) then
      WAnyBranching_flat = TaM_
  else
      print *, "error ",xRnd
      stop
  endif


RETURN
END FUNCTION





FUNCTION WAnyBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WAnyBranching


  if( xRnd .le. Br_W_ud ) then
      WAnyBranching = Up_
  elseif(xRnd .le. Br_W_ud+Br_W_cs ) then
      WAnyBranching = Chm_
  elseif(xRnd .le. Br_W_ud+Br_W_cs+Br_W_en ) then
      WAnyBranching = ElM_
  elseif(xRnd .le. Br_W_ud+Br_W_cs+Br_W_en+Br_W_mn ) then
      WAnyBranching = MuM_
  elseif(xRnd .le. Br_W_ud+Br_W_cs+Br_W_en+Br_W_mn+Br_W_tn ) then
      WAnyBranching = TaM_
  else
      print *, "error ",xRnd
      stop
  endif


!   if( xRnd .le. scale_alpha_W_ud*Br_W_ud ) then
!       WAnyBranching = Up_
!   elseif(xRnd .le. scale_alpha_W_ud*Br_W_ud+scale_alpha_W_cs*Br_W_cs ) then
!       WAnyBranching = Chm_
!   elseif(xRnd .le. scale_alpha_W_ud*Br_W_ud+scale_alpha_W_cs*Br_W_cs+scale_alpha_W_ln*Br_W_en ) then
!       WAnyBranching = ElM_
!   elseif(xRnd .le. scale_alpha_W_ud*Br_W_ud+scale_alpha_W_cs*Br_W_cs+scale_alpha_W_ln*Br_W_en+scale_alpha_W_ln*Br_W_mn ) then
!       WAnyBranching = MuM_
!   elseif(xRnd .le. scale_alpha_W_ud*Br_W_ud+scale_alpha_W_cs*Br_W_cs+scale_alpha_W_ln*Br_W_en+scale_alpha_W_ln*Br_W_mn+scale_alpha_W_ln*Br_W_tn ) then
!       WAnyBranching = TaM_
!   else
!       print *, "error ",xRnd
!       stop
!   endif

! print *, "checker 2",Br_W_ud+Br_W_cs+Br_W_en+Br_W_mn+Br_W_tn
! print *, "checker 2",scale_alpha_W_ud*Br_W_ud+scale_alpha_W_cs*Br_W_cs+scale_alpha_W_ln*Br_W_en+scale_alpha_W_ln*Br_W_mn+scale_alpha_W_ln*Br_W_tn
! pause
  
RETURN
END FUNCTION





SUBROUTINE VVBranchings(MY_IDUP,ICOLUP,ColorBase)
use ModParameters
implicit none
integer :: MY_IDUP(4:9),ICOLUP(1:2,6:9),DKFlavor,ICOLUP_Base
integer, optional ::ColorBase
real(8) :: DKRnd

!    particle associations:
!    
!    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor

   if (present(ColorBase)) then
       ICOLUP_BASE = ColorBase
   else
       ICOLUP_BASE = 800
   endif

   ICOLUP(:,:) = 0
   if( DecayMode1.eq.0 ) then! Z1->2l
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZLepBranching( DKRnd )!= ElM or MuM
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
   elseif( DecayMode1.eq.1 ) then! Z1->2q
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZQuaBranching_flat( DKRnd )!= Up,Dn,Chm,Str,Bot
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
        ICOLUP(1:2,6) = (/            0,ICOLUP_BASE+3/)
        ICOLUP(1:2,7) = (/ICOLUP_BASE+3,            0/)
   elseif( DecayMode1.eq.2 ) then! Z1->2tau
        MY_IDUP(4) = Z0_
        MY_IDUP(6) = TaP_
        MY_IDUP(7) = TaM_
   elseif( DecayMode1.eq.3 ) then! Z1->2nu
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZNuBranching( DKRnd )!= NuE,NuM,NuT
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
   elseif( DecayMode1.eq.4 ) then! W1(+)->lnu
        call random_number(DKRnd)
        MY_IDUP(4) = Wp_
        DKFlavor = WLepBranching( DKRnd )!= ElM or MuM
        MY_IDUP(6) = +abs(DKFlavor)     ! lepton(+)
        MY_IDUP(7) = +abs(DKFlavor)+7   ! neutrino        
   elseif( DecayMode1.eq.5 ) then! W1(+)->2q
        call random_number(DKRnd)
        MY_IDUP(4) = Wp_
        DKFlavor = WQuaUpBranching( DKRnd )!= Up,Chm
!         MY_IDUP(6) = -abs(DKFlavor)-1  ! anti-dn flavor
!         MY_IDUP(7) = +abs(DKFlavor)    ! up flavor
        MY_IDUP(7) = +abs(DKFlavor)           ! up flavor
        MY_IDUP(6) = GetCKMPartner(MY_IDUP(7))! anti-dn flavor         
        ICOLUP(1:2,6) = (/            0,ICOLUP_BASE+3/)
        ICOLUP(1:2,7) = (/ICOLUP_BASE+3,            0/)
   elseif( DecayMode1.eq.6 ) then! W1(+)->taunu
        MY_IDUP(4) = Wp_
        MY_IDUP(6) = TaP_
        MY_IDUP(7) = NuT_
   elseif( DecayMode1.eq.7 ) then! photon
        MY_IDUP(4) = Pho_
        MY_IDUP(6) = -9999
        MY_IDUP(7) = -9999
   elseif( DecayMode1.eq.8 ) then! Z1->2l+2tau
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZLepPlusTauBranching( DKRnd )!= ElM or MuM or TaM
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
   elseif( DecayMode1.eq.9 ) then! Z1-> anything
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZAnyBranching_flat( DKRnd )
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
        if(IsAQuark(DKFlavor)) then
           ICOLUP(1:2,6) = (/            0,ICOLUP_BASE+3/)
           ICOLUP(1:2,7) = (/ICOLUP_BASE+3,            0/)
        endif
   elseif( DecayMode1.eq.10 ) then! W1(+)->l+tau  +nu
        call random_number(DKRnd)
        MY_IDUP(4) = Wp_
        DKFlavor = WLepPlusTauBranching( DKRnd )!= ElM or MuM or TaM
        MY_IDUP(6) = +abs(DKFlavor)     ! lepton(+)
        MY_IDUP(7) = +abs(DKFlavor)+7   ! neutrino
   elseif( DecayMode1.eq.11 ) then! W1(+)-> anything
        call random_number(DKRnd)
        MY_IDUP(4) = Wp_
        DKFlavor = WAnyBranching_flat( DKRnd )
        if(IsAQuark(DKFlavor)) then
!            MY_IDUP(6) = -abs(DKFlavor)-1  ! anti-dn flavor  
!            MY_IDUP(7) = +abs(DKFlavor)    ! up flavor
           MY_IDUP(7) = +abs(DKFlavor)           ! up flavor
           MY_IDUP(6) = GetCKMPartner(MY_IDUP(7))! anti-dn flavor  
           ICOLUP(1:2,6) = (/            0,ICOLUP_BASE+3/)
           ICOLUP(1:2,7) = (/ICOLUP_BASE+3,            0/)
        else
           MY_IDUP(6) = +abs(DKFlavor)     ! lepton(+)
           MY_IDUP(7) = +abs(DKFlavor)+7   ! neutrino
        endif
   endif


   if( DecayMode2.eq.0 ) then! Z2->2l (sample over el,mu)
        call random_number(DKRnd)
        MY_IDUP(5) = Z0_
        DKFlavor = ZLepBranching( DKRnd )!= ElM or MuM
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
   elseif( DecayMode2.eq.1 ) then! Z2->2q
        call random_number(DKRnd)
        MY_IDUP(5) = Z0_
        DKFlavor = ZQuaBranching_flat( DKRnd )!= Up,Dn,Chm,Str,Bot
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
        ICOLUP(1:2,8) = (/            0,ICOLUP_BASE+4/)
        ICOLUP(1:2,9) = (/ICOLUP_BASE+4,            0/)
   elseif( DecayMode2.eq.2 ) then! Z2->2tau
        MY_IDUP(5) = Z0_
        MY_IDUP(8) = TaP_
        MY_IDUP(9) = TaM_
   elseif( DecayMode2.eq.3 ) then! Z2->2nu
        call random_number(DKRnd)
        MY_IDUP(5) = Z0_
        DKFlavor = ZNuBranching( DKRnd )!= NuE,NuM,NuT
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
   elseif( DecayMode2.eq.4 ) then! W2(-)->lnu
        call random_number(DKRnd)
        MY_IDUP(5) = Wm_
        DKFlavor = WLepBranching( DKRnd )!= ElM or MuM
        MY_IDUP(8) = -abs(DKFlavor)-7   ! anti-neutrino
        MY_IDUP(9) = -abs(DKFlavor)     ! lepton(-)
   elseif( DecayMode2.eq.5 ) then! W2(-)->2q (sample over u,d,s,c)
        call random_number(DKRnd)
        MY_IDUP(5) = Wm_
        DKFlavor = WQuaUpBranching( DKRnd )!= Up,Chm
!         MY_IDUP(8) = -abs(DKFlavor)    ! anti-up flavor
!         MY_IDUP(9) = +abs(DKFlavor)+1  ! dn flavor
        MY_IDUP(8) = -abs(DKFlavor)           ! up flavor
        MY_IDUP(9) = GetCKMPartner(MY_IDUP(8))! dn flavor
        ICOLUP(1:2,8) = (/            0,ICOLUP_BASE+4/)
        ICOLUP(1:2,9) = (/ICOLUP_BASE+4,            0/)
   elseif( DecayMode2.eq.6 ) then! W2(-)->taunu
        MY_IDUP(5) = Wm_
        MY_IDUP(8) = ANuT_
        MY_IDUP(9) = TaM_
   elseif( DecayMode2.eq.7 ) then! photon
        MY_IDUP(5) = Pho_
        MY_IDUP(8) = -9999
        MY_IDUP(9) = -9999
   elseif( DecayMode2.eq.8 ) then! Z2->2l+2tau
        call random_number(DKRnd)
        MY_IDUP(5) = Z0_
        DKFlavor = ZLepPlusTauBranching( DKRnd )!= ElM or MuM or TaM
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
   elseif( DecayMode2.eq.9 ) then! Z2-> anything
        call random_number(DKRnd)
        MY_IDUP(5) = Z0_
        DKFlavor = ZAnyBranching_flat( DKRnd )
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
        if(IsAQuark(DKFlavor)) then
           ICOLUP(1:2,8) = (/            0,ICOLUP_BASE+4/)
           ICOLUP(1:2,9) = (/ICOLUP_BASE+4,            0/)
        endif
   elseif( DecayMode2.eq.10 ) then! W2(-)->l+tau + nu
        call random_number(DKRnd)
        MY_IDUP(5) = Wm_
        DKFlavor = WLepPlusTauBranching( DKRnd )!= ElM or MuM or TaM
        MY_IDUP(8) = -abs(DKFlavor)-7   ! anti-neutrino
        MY_IDUP(9) = -abs(DKFlavor)     ! lepton(-)
   elseif( DecayMode2.eq.11 ) then! W2(-)-> anything
        call random_number(DKRnd)
        MY_IDUP(5) = Wm_
        DKFlavor = WAnyBranching_flat( DKRnd )
        if(IsAQuark(DKFlavor)) then
!            MY_IDUP(8) = -abs(DKFlavor)    ! anti-up flavor
!            MY_IDUP(9) = +abs(DKFlavor)+1  ! dn flavor
           MY_IDUP(8) = -abs(DKFlavor)           ! up flavor
           MY_IDUP(9) = GetCKMPartner(MY_IDUP(8))! dn flavor
           ICOLUP(1:2,8) = (/            0,ICOLUP_BASE+4/)
           ICOLUP(1:2,9) = (/ICOLUP_BASE+4,            0/)
        else
           MY_IDUP(8) = -abs(DKFlavor)-7   ! anti-neutrino
           MY_IDUP(9) = -abs(DKFlavor)     ! lepton(-)
        endif
   endif


RETURN
END SUBROUTINE





FUNCTION GetCKMPartner( Flavor )
use modMisc
use modParameters
implicit none
integer :: Flavor,GetCKMPartner
real(8) :: FlavorRnd,sumCKM,Vsq(1:3)

    call random_number(FlavorRnd)

    
    if( abs(Flavor).eq.abs(Up_) ) then
        Vsq(1) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Dn_)) ))**2
        Vsq(2) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Str_)) ))**2
        Vsq(3) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Bot_)) ))**2
        Vsq(:) = Vsq(:)/scale_alpha_W_ud

        sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
        FlavorRnd = FlavorRnd*sumCKM

        if( FlavorRnd.le.Vsq(1) ) then!  u-->d
           GetCKMPartner = -sign(1,Flavor) * abs(Dn_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(1)) ) then!  u-->s
           GetCKMPartner = -sign(1,Flavor) * abs(Str_)
        else!  u-->b
           GetCKMPartner = -sign(1,Flavor) * abs(Bot_)
        endif
        
    elseif( abs(Flavor).eq.abs(Chm_) ) then
        Vsq(1) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Dn_)) ))**2
        Vsq(2) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Str_)) ))**2
        Vsq(3) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Bot_)) ))**2    
        Vsq(:) = Vsq(:)/scale_alpha_W_cs
        
        sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
        FlavorRnd = FlavorRnd*sumCKM    

        if( FlavorRnd.le.Vsq(2) ) then!  c-->s
           GetCKMPartner = -sign(1,Flavor) * abs(Str_)     
        elseif( FlavorRnd.le.(Vsq(1)+Vsq(2)) ) then!  c-->d
           GetCKMPartner = -sign(1,Flavor) * abs(Dn_)
        else!  c-->b
           GetCKMPartner = -sign(1,Flavor) * abs(Bot_)
        endif

    elseif( abs(Flavor).eq.abs(Top_) ) then
        Vsq(1) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Dn_)) ))**2
        Vsq(2) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Str_)) ))**2
        Vsq(3) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Bot_)) ))**2

        sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
        FlavorRnd = FlavorRnd*sumCKM    
        
        if( FlavorRnd.le.Vsq(3) ) then!  t-->b
           GetCKMPartner = -sign(1,Flavor) * abs(Bot_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(3)) ) then!  t-->s
           GetCKMPartner = -sign(1,Flavor) * abs(Str_)
        else!  t-->d
           GetCKMPartner = -sign(1,Flavor) * abs(Dn_)
        endif


    elseif( abs(Flavor).eq.abs(Dn_) ) then
        Vsq(1) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Up_)) ))**2
        Vsq(2) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Chm_)) ))**2
        Vsq(3) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Top_)) ))**2
        Vsq(:) = Vsq(:)/scale_alpha_W_ud

        sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
        FlavorRnd = FlavorRnd*sumCKM

        if( FlavorRnd.le.Vsq(1) ) then!  d-->u
           GetCKMPartner = -sign(1,Flavor) * abs(Up_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(1)) ) then!  d-->c
           GetCKMPartner = -sign(1,Flavor) * abs(Chm_)
        else!  d-->t
           GetCKMPartner = -sign(1,Flavor) * abs(Top_)
        endif
        
    elseif( abs(Flavor).eq.abs(Str_) ) then
        Vsq(1) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Up_)) ))**2
        Vsq(2) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Chm_)) ))**2
        Vsq(3) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Top_)) ))**2    
        Vsq(:) = Vsq(:)/scale_alpha_W_cs
        
        sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
        FlavorRnd = FlavorRnd*sumCKM    

        if( FlavorRnd.le.Vsq(2) ) then!  s-->c
           GetCKMPartner = -sign(1,Flavor) * abs(Chm_)     
        elseif( FlavorRnd.le.(Vsq(1)+Vsq(2)) ) then!  s-->u
           GetCKMPartner = -sign(1,Flavor) * abs(Up_)
        else!  s-->t
           GetCKMPartner = -sign(1,Flavor) * abs(Top_)
        endif

    elseif( abs(Flavor).eq.abs(Bot_) ) then
        Vsq(1) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Up_)) ))**2
        Vsq(2) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Chm_)) ))**2
        Vsq(3) = (CKM( convertLHEreverse(abs(Flavor)), convertLHEreverse(abs(Top_)) ))**2

        sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
        FlavorRnd = FlavorRnd*sumCKM    
        
        if( FlavorRnd.le.Vsq(3) ) then!  b-->t
           GetCKMPartner = -sign(1,Flavor) * abs(Top_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(3)) ) then!  b-->c
           GetCKMPartner = -sign(1,Flavor) * abs(Chm_)
        else!  b -->u
           GetCKMPartner = -sign(1,Flavor) * abs(Up_)
        endif
    
    else
        call Error("Dn flavor conversion not yet implemented")
    endif

RETURN
END FUNCTION


FUNCTION GetCKMPartner_flat( Flavor )
use modMisc
use modParameters
implicit none
integer :: Flavor,GetCKMPartner_flat
real(8) :: FlavorRnd,sumCKM,Vsq(1:3)

    call random_number(FlavorRnd)
    Vsq(:) = 1d0
    sumCKM = Vsq(1)+Vsq(2)+Vsq(3)
    FlavorRnd = FlavorRnd*sumCKM
    
    if( abs(Flavor).eq.abs(Up_) ) then
        if( FlavorRnd.le.Vsq(1) ) then!  u-->d
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Dn_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(1)) ) then!  u-->s
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Str_)
        else!  u-->b
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Bot_)
        endif
        
    elseif( abs(Flavor).eq.abs(Chm_) ) then
        if( FlavorRnd.le.Vsq(2) ) then!  c-->s
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Str_)     
        elseif( FlavorRnd.le.(Vsq(1)+Vsq(2)) ) then!  c-->d
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Dn_)
        else!  c-->b
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Bot_)
        endif

    elseif( abs(Flavor).eq.abs(Top_) ) then
        if( FlavorRnd.le.Vsq(3) ) then!  t-->b
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Bot_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(3)) ) then!  t-->s
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Str_)
        else!  t-->d
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Dn_)
        endif


    elseif( abs(Flavor).eq.abs(Dn_) ) then
        if( FlavorRnd.le.Vsq(1) ) then!  d-->u
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Up_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(1)) ) then!  d-->c
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Chm_)
        else!  d-->t
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Top_)
        endif
        
    elseif( abs(Flavor).eq.abs(Str_) ) then
        if( FlavorRnd.le.Vsq(2) ) then!  s-->c
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Chm_)     
        elseif( FlavorRnd.le.(Vsq(1)+Vsq(2)) ) then!  s-->u
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Up_)
        else!  s-->t
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Top_)
        endif

    elseif( abs(Flavor).eq.abs(Bot_) ) then
        if( FlavorRnd.le.Vsq(3) ) then!  b-->t
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Top_)
        elseif( FlavorRnd.le.(Vsq(2)+Vsq(3)) ) then!  b-->c
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Chm_)
        else!  b -->u
           GetCKMPartner_flat = -sign(1,Flavor) * abs(Up_)
        endif
    
    else
        call Error("Dn flavor conversion not yet implemented")
    endif

RETURN
END FUNCTION



SUBROUTINE BranchingCounter(MY_IDUP)
use ModParameters
implicit none
integer :: MY_IDUP(6:9)


    if( IsAZDecay(DecayMode1) ) then
        if(  all(MY_IDUP(6:7)-(/ElP_,ElM_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/MuP_,MuM_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/TaP_,TaM_/).eq.0) &
          )  Br_Z_ll_counter=Br_Z_ll_counter+1 

        if(  all(MY_IDUP(6:7)-(/ANuE_,NuE_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/ANuM_,NuM_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/ANuT_,NuT_/).eq.0) &
          )  Br_Z_inv_counter=Br_Z_inv_counter+1 

        if(  all(MY_IDUP(6:7)-(/ADn_,Dn_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/AStr_,Str_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/ABot_,Bot_/).eq.0) &
          )  Br_Z_dd_counter=Br_Z_dd_counter+1 

        if(  all(MY_IDUP(6:7)-(/AUp_,Up_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/AChm_,Chm_/).eq.0) &
          )  Br_Z_uu_counter=Br_Z_uu_counter+1 

    elseif( IsAWDecay(DecayMode1) ) then
        if(  all(MY_IDUP(6:7)-(/ElP_,NuE_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/MuP_,NuM_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/TaP_,NuT_/).eq.0) &
          )  Br_W_ll_counter=Br_W_ll_counter+1 

        if(  all(MY_IDUP(6:7)-(/ADn_,Up_/).eq.0) &
        .or. all(MY_IDUP(6:7)-(/AStr_,Chm_/).eq.0) &
          )  Br_W_ud_counter=Br_W_ud_counter+1 

    elseif( IsAPhoton(DecayMode1) ) then
    endif

RETURN
END SUBROUTINE





FUNCTION WhichBin(NHisto,Value)
implicit none
integer :: WhichBin,NHisto
real(8) :: Value

   WhichBin = (Value-Histo(NHisto)%LowVal)/Histo(NHisto)%BinSize + 1
   if( WhichBin.lt.0 ) then
      WhichBin = 0
   elseif( WhichBin.gt.Histo(NHisto)%NBins ) then
      WhichBin = Histo(NHisto)%NBins+1
   endif

RETURN
END FUNCTION





! FUNCTION WhichXBin(NHisto,XValue)
! use ModParameters
! implicit none
! integer :: WhichXBin,NHisto
! real(8) :: XValue
! integer :: i
! include "vegas_common.f"
! 
!     whichxbin = int( xValue*NPart )!  uniform distribution
! 
! 
! !    do i=1,50!                         distribution according to vegas grid
! !       if( XValue .lt. xi(i,NHisto) ) then
! !           WhichXBin=i
! !           return
! !       endif
! !    enddo
! RETURN
! END FUNCTION







SUBROUTINE IntoHisto(NHisto,NBin,Value)
implicit none
integer :: NHisto,NBin
real(8) :: Value

     Histo(NHisto)%Value(NBin) = Histo(NHisto)%Value(NBin)  + Value
     Histo(NHisto)%Value2(NBin)= Histo(NHisto)%Value2(NBin) + Value**2
     Histo(NHisto)%Hits(NBin)  = Histo(NHisto)%Hits(NBin)+1

RETURN
END SUBROUTINE







SUBROUTINE boost2Lab(x1,x2,NumPart,Mom)
implicit none
real(8) Mom(1:4,1:NumPart)
real(8) x1,x2
real(8) gamma,betagamma,MomTmp1,MomTmp4
integer :: i,NumPart

  gamma     = (x1+x2)/2d0/dsqrt(x1*x2)
  betagamma = (x2-x1)/2d0/dsqrt(x1*x2)

  do i=1,NumPart
      MomTmp1=Mom(1,i)
      MomTmp4=Mom(4,i)
      Mom(1,i)= gamma*MomTmp1 - betagamma*MomTmp4
      Mom(4,i)= gamma*MomTmp4 - betagamma*MomTmp1
  enddo

RETURN
END SUBROUTINE



SUBROUTINE getHiggsDecayLength(ctau)
use ModParameters
implicit none
real(8) :: x,xp,xpp,Len0,propa,ctau,ctau0
integer :: loop


     ctau  = 0d0
     ctau0 = HiggsDecayLengthMM
     if( ctau0.lt.1d-16 ) RETURN

     do loop=1,4000000!  4Mio. tries otherwise return zero
          call random_number(x)
          xp = 10*x*ctau0             ! scan between 0..10*ctau0

          propa = dexp( -xp/(ctau0) ) ! the max of propa is 1.0
          call random_number(xpp)
          if( xpp.lt.propa ) then!   accept
                ctau = xp
                RETURN
          endif
     enddo


RETURN
END SUBROUTINE




SUBROUTINE SmearExternal(xRnd,Mass,Width,MinEnergy,MaxEnergy,invMass,Jacobi)
use ModParameters
real(8) :: xRnd,Mass,Width,invMass,Jacobi,MinEnergy,MaxEnergy
real(8) :: r,rmin,rmax,BW


IF( Width.lt.(1d-6)*GeV ) THEN
   invMass = Mass
   Jacobi  = 1d0
ELSE
   rmin=1d0/(Width*Mass) * datan( (MinEnergy**2-Mass**2)/(Width*Mass)  )    
   rmax=1d0/(Width*Mass) * datan( (MaxEnergy**2-Mass**2)/(Width*Mass)  )
   r = xRnd*(rmax-rmin) + rmin
   
   invMass = dsqrt(dabs( Mass*Width * dtan( Mass*Width*r )  +  Mass**2 ))
!    BW = (invMass**2 - Mass**2)**2 + Mass**2 * Width**2            ! Breit-Wiegner propagator 
   BW = (Mass**2 * Width**2) * ( dtan( Mass*Width*r )**2 + 1d0 )    ! equivalent to above

   Jacobi  = (rmax-rmin) * BW   /(2d0*Pi)                           ! this Jacobian has [GeV^2]; it becomes [GeV^-2] when hitting the sq. propagator 
ENDIF                                                               ! factor 2*Pi comes from integr. measure


return
END SUBROUTINE




SUBROUTINE PDFMapping(MapType,yRnd,eta1,eta2,Ehat,sHatJacobi)
use ModParameters
use ModMisc
implicit none
integer :: MapType
real(8) :: yRnd(1:2),eta1,eta2,EHat,sHatJacobi,tau,nPotMap,z,sbar,fmax
real(8) :: etamin, Ymax, Y, Ymin, MThresh

  MThresh = M_Reso

  if( MapType.eq.1 ) then!  no mapping
      eta1 = yRnd(1)
      eta2 = yRnd(2)
      sHatJacobi = 1d0
  elseif( MapType.eq.2 ) then!  exponential mapping
      tau = (MThresh/Collider_Energy)**2
      eta1 = tau**yRnd(1)
      eta2 = tau**( (1d0-yRnd(1))*yRnd(2) )
      sHatJacobi = dlog(tau)**2*(1d0-yRnd(1))*eta1*eta2
  elseif( MapType.eq.3 ) then!  linear mapping
      tau = (MThresh/Collider_Energy)**2
      eta1 = (1d0-tau)*yRnd(1) + tau
      eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
      sHatJacobi = (1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)
  elseif( MapType.eq.4 ) then!  MCFM mapping
      tau = dexp(dlog(((MThresh/Collider_Energy)**2))*yRnd(1))
      eta1 = dsqrt(tau)*dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
      eta2 = dsqrt(tau)/dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
      sHatJacobi = dlog(((MThresh/Collider_Energy)**2))*tau*dlog(tau)
  elseif( MapType.eq.5 ) then!  nPotMap mapping
!       nPotMap = 0.5d0
!       tau = (MThresh/Collider_Energy)**2
!       yRnd(1) = yRnd(1)**nPotMap
!       yRnd(2) = yRnd(2)**nPotMap
!       eta1 = (1d0-tau) * yRnd(1) + tau
!       eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
!       sHatJacobi=(1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*nPotMap**2*((yRnd(1)*yRnd(2))**(1d0/nPotMap))**(nPotMap-1d0)
  elseif( MapType.eq.10 ) then!  Breit-Wigner mapping
      fmax = 1d0/M_Reso/Ga_Reso * ( datan((Collider_Energy**2-M_Reso**2)/M_Reso/Ga_Reso) - datan(-M_Reso/Ga_Reso) )
      sbar = M_Reso*Ga_Reso * dtan(fmax*yRnd(1)*M_Reso*Ga_Reso - datan(M_Reso/Ga_Reso) ) + M_Reso**2
      z = sbar/Collider_Energy**2
      eta1 = z + (1d0-z)*yRnd(2)
      eta2 = z/eta1
      sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (sbar - M_Reso**2)**2 + M_Reso**2*Ga_Reso**2 )

  elseif (MapType.eq.11) then ! delta-function map

     etamin = M_Reso**2/Collider_Energy**2
     eta2 = etamin + (1d0-etamin)*yRnd(2)
     eta1 = etamin/eta2
     fmax = 0.5d0*pi/M_Reso**3/Ga_Reso
     sHatJacobi = fmax*etamin/eta2* ( (Collider_Energy**2*eta1*eta2 - M_Reso**2)**2 + M_Reso**2*Ga_Reso**2 )

  elseif (MapType.eq.12) then ! delta-function map / new
     Ymax = dlog(Collider_Energy/M_Reso)
     Y = -Ymax + 2d0*Ymax*yRnd(2)
     eta1 = M_Reso/Collider_Energy*exp(Y)
     eta2 = M_Reso/Collider_Energy*exp(-Y)     
     fmax = 0.5d0*pi/M_Reso**3/Ga_Reso*2d0*Ymax
     sHatJacobi = fmax*(M_Reso**2*Ga_Reso**2 )
     
! ! ! !      
!       fmax = 1d0/M_Reso/Ga_Reso * ( datan((Collider_Energy**2-M_Reso**2)/M_Reso/Ga_Reso) - datan(-M_Reso/Ga_Reso) )
!       sbar = M_Reso**2
!       z = sbar/Collider_Energy**2
!       eta1 = z + (1d0-z)*yRnd(2)
!       eta2 = z/eta1
!       sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (sbar - M_Reso**2)**2 + M_Reso**2*Ga_Reso**2 )     
! ! ! ! ! !      


  elseif( MapType.eq.13 ) then!  Breit-Wigner mapping with M = M_Z + M_h
      fmax = 1d0/(M_Reso+M_Z)/Ga_Z * ( datan((Collider_Energy**2-(M_Reso+M_Z)**2)/(M_Reso+M_Z)/Ga_Z) - datan(-(M_Reso+M_Z)/Ga_Z) )
      sbar = (M_Reso+M_Z)*Ga_Z * dtan(fmax*yRnd(1)*(M_Reso+M_Z)*Ga_Z - datan((M_Reso+M_Z)/Ga_Z) ) + (M_Reso+M_Z)**2
      z = sbar/Collider_Energy**2
      eta1 = z + (1d0-z)*yRnd(2)
      eta2 = z/eta1
      sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (sbar - (M_Reso+M_Z)**2)**2 + (M_Reso+M_Z)**2*Ga_Z**2 )
      
  elseif( MapType.eq.14 ) then! Z-Higgs associate production   
      etamin = (M_Reso+M_Z-2d0*Ga_Z)/Collider_Energy
      z = 1d0/(1d0-etamin)
      Ymin = ((z-1d0)/(z-yRnd(1)))**2
      Ymax = 1d0 / Ymin
      Ymin = 0.5d0*dlog(Ymin)
      Ymax = 0.5d0*dlog(Ymax)
      Y = Ymin + yrnd(2)*(Ymax-Ymin)
      eta1 = (z-1d0)/(z-yRnd(1))*dexp(Y)
      eta2 = (z-1d0)/(z-yRnd(1))*dexp(-Y)
      sHatJacobi = 2d0*(z-1d0)**2/(z-yRnd(1))**3*(Ymax-Ymin)
      
  elseif( MapType.eq.15 ) then! W-Higgs associate production   
      etamin = (M_Reso+M_W-2d0*Ga_W)/Collider_Energy
      z = 1d0/(1d0-etamin)
      Ymin = ((z-1d0)/(z-yRnd(1)))**2
      Ymax = 1d0 / Ymin
      Ymin = 0.5d0*dlog(Ymin)
      Ymax = 0.5d0*dlog(Ymax)
      Y = Ymin + yrnd(2)*(Ymax-Ymin)
      eta1 = (z-1d0)/(z-yRnd(1))*dexp(Y)
      eta2 = (z-1d0)/(z-yRnd(1))*dexp(-Y)
      sHatJacobi = 2d0*(z-1d0)**2/(z-yRnd(1))**3*(Ymax-Ymin)
      
  elseif( MapType.eq.16 ) then! H+j   
      etamin = (M_Reso-3d0*Ga_Reso)/Collider_Energy
      z = 1d0/(1d0-etamin)
      Ymin = ((z-1d0)/(z-yRnd(1)))**2
      Ymax = 1d0 / Ymin
      Ymin = 0.5d0*dlog(Ymin)
      Ymax = 0.5d0*dlog(Ymax)
      Y = Ymin + yrnd(2)*(Ymax-Ymin)
      eta1 = (z-1d0)/(z-yRnd(1))*dexp(Y)
      eta2 = (z-1d0)/(z-yRnd(1))*dexp(-Y)
      sHatJacobi = 2d0*(z-1d0)**2/(z-yRnd(1))**3*(Ymax-Ymin)
  else
      call Error("PDF mapping not available")
  endif

  EHat = Collider_Energy*dsqrt(eta1*eta2)

RETURN
END SUBROUTINE



FUNCTION ReweightBWPropagator(shat)! shat is the resonance inv. mass squared
use modParameters
implicit none
real(8) :: ReweightBWPropagator,shat
real(8) :: BreitWigner,BreitWigner_Run,muH,gaH,mubarH,gabarH


    ReweightBWPropagator = 1d0
    
    if( WidthScheme.eq.1) then! running width
        BreitWigner = 1d0/( (shat-M_Reso**2)**2 + (M_Reso*Ga_Reso)**2 )
        BreitWigner_Run =  1d0/( (shat-M_Reso**2)**2 + (shat*Ga_Reso/M_Reso)**2 )
       
    elseif( WidthScheme.eq.2) then! Passarino'S CPS
        BreitWigner = 1d0/( (shat-M_Reso**2)**2 + (M_Reso*Ga_Reso)**2 )

        call CALL_HTO(dsqrt(dabs(shat))*100d0,m_top*100d0,gabarH,mubarH)
        if( IsNaN(gabarH) .or. IsNaN(mubarH) ) then
          print *, "Passarino's CALL_HTO returned a NaN"
          print *, "gabarH,mubarH,Ehat)",gabarH,mubarH,dsqrt(dabs(shat))*100d0
          print *, "returning weight 1.0"          
          ReweightBWPropagator = 1d0
          RETURN
        endif       
        mubarH = mubarH/100d0
        gabarH = gabarH/100d0

        muH = dsqrt( mubarH**2/(1d0+(gabarH/mubarH)**2) )
        gaH = muH/mubarH*gabarH
        
        BreitWigner_Run = 1d0 /( (shat-muH**2)**2 + (muH*gaH)**2 )
    endif

    ReweightBWPropagator = BreitWigner_Run/BreitWigner
    

RETURN
END FUNCTION



RECURSIVE SUBROUTINE JetAlgo_kt(Rsep_jet,PartonList,MomParton,NJet,JetList,MomJet)  ! initial call must have NJet=0 and MomJet(1:4,:) = MomPartons(1:4,:)
use ModMisc
use ModParameters
implicit none
integer :: PartonList(:), JetList(:)
real(8) :: MomParton(:,:)
integer :: NJet
real(8) :: MomJet(:,:)
real(8) :: Rsep_jet
integer :: NParton,i,j,k,dii_minp,dij_minp,ij(1:120,1:2)  ! max. partons=8, (8-1)! =5040  ! max. partons=6, (6-1)! =120
real(8) :: dii(1:6),dij(1:120),Rij,dii_min,dij_min!,eta(1:5),phi(1:5)


NParton = size(PartonList)
if( NParton.eq.0 ) then
    return
elseif( NParton.eq.1 ) then
!   print *, "HARD JET", PartonList(NParton)
   NJet = NJet +1
   JetList(NJet) = PartonList(NParton)
   return
endif

!generate dii, dij
do i=1,NParton
   dii(i) = get_PT2(MomParton(1:4,PartonList(i)))**(-1d0)
enddo

k=0
do i=1,NParton-1
do j=i+1,NParton
   k = k+1
   Rij = get_R( MomParton(1:4,PartonList(i)), MomParton(1:4,PartonList(j)) )
   dij(k) = dmin1(dii(i),dii(j)) * (Rij/Rsep_jet)**2
   ij(k,1)=i
   ij(k,2)=j
enddo
enddo
!!print *, dii(1:NParton)
!!print *, dij(1:k)


! find minima
dii_min=dii(1)
dii_minp=1
do i=2,NParton
  if( dii(i).lt.dii_min ) then
    dii_min = dii(i)
    dii_minp= i
  endif
enddo

dij_min=dij(1)
dij_minp=1
do i=2,k
  if( dij(i).lt.dij_min ) then
    dij_min = dij(i)
    dij_minp= i
  endif
enddo


if( dii_min.lt.dij_min ) then ! no recombination
   NJet = NJet +1
   k=0
!   print *, "HARD JET", PartonList(dii_minp)
   JetList(NJet) = PartonList(dii_minp)
   MomJet(1:4,PartonList(dii_minp)) = MomParton(1:4,PartonList(dii_minp))
   do i=1,NParton!  remove momenta dii_minp from parton list
      if( i.eq.dii_minp ) cycle
      k=k+1
      PartonList(k) = PartonList(i)
   enddo
   call JetAlgo_kt(Rsep_jet,PartonList(1:k),MomParton(:,:),NJet,JetList(:),MomJet(:,:))
else ! recombination of dij(dij_min)
!    if( RecombPrescr.eq.1 ) then
!         MomJet(1:4,PartonList(ij(dij_minp,1))) = EllisSoperComb(MomJet(1:4,PartonList(ij(dij_minp,1))),MomJet(1:4,PartonList(ij(dij_minp,2))))   ! Ellis-Soper combination
!    else
        MomJet(1:4,PartonList(ij(dij_minp,1))) = MomJet(1:4,PartonList(ij(dij_minp,1))) + MomJet(1:4,PartonList(ij(dij_minp,2)))                 ! four vector addition
!    endif
   MomJet(1:4,PartonList(ij(dij_minp,2))) = 0d0
   k=0
!   print *, "RECOMB.",PartonList(ij(dij_minp,1)),PartonList(ij(dij_minp,2))
   do i=1,NParton!  remove momenta j from parton list
      if( i.eq.ij(dij_minp,2) ) then
         cycle
      endif
      k=k+1
      PartonList(k) = PartonList(i)
   enddo
   call JetAlgo_kt(Rsep_jet,PartonList(1:k),MomJet(:,:),NJet,JetList(:),MomJet(:,:))
endif

return
END SUBROUTINE





SUBROUTINE EvalPhaseSpace_VH(yRnd,MomExt,inv_mass,mass,PSWgt)
!use modMisc
implicit none

      real(8), intent(in) :: yRnd(1:20),mass(9,2)
      real(8) :: phi, beta, gamma
      real(8) :: temp_vector(4), temp_boost(4)
      integer :: i
      real(8), parameter :: Pi = 3.14159265358979323846d0
      real(8), parameter :: Twopi = 2d0 * Pi
      real(8) :: MomExt(1:4,1:9)
      real(8) :: MomDummy(1:4)
      !double precision four_momentum(7,4)
      real(8), intent(out) :: PSWgt,inv_mass(9)
! 1=E, 2,3,4=p_x,y,z
!psg_mass(mass, width)
      real(8) :: cm_abs3p(9)
      real(8) :: cm_sin_theta(9), cm_cos_theta(9)
      real(8) :: cm_sin_phi(9), cm_cos_phi(9)
!use Cauchy distribution for Breit-Wigner distribution for the invariant mass of 2 and 3?
!      logical, parameter :: breit_wigner = .true.
      real(8) :: jacobian4, jacobian5

!3333333333
!invariant mass of 3
      inv_mass(3) = dsqrt((MomExt(1,3)+MomExt(4,3))*(MomExt(1,3)-MomExt(4,3)))
!generating invariant mass of 4 and 5
!if using uniform distribution
!      if(breit_wigner.eqv..false.)then
!        inv_mass(4) = yRnd(12)
!        inv_mass(5) = yRnd(13)*(1d0-inv_mass(4))
!        jacobian4 = (1d0-inv_mass(4))
!        inv_mass(4) = inv_mass(4) * inv_mass(3)
!        jacobian4 = inv_mass(3)
!        inv_mass(5) = inv_mass(5) * inv_mass(3)
!        jacobian5 = jacobian5 * inv_mass(3)
!      else
!if using Breit-Wigner distribution
!
        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), inv_mass(3)**2, jacobian4)))
        inv_mass(5) = dsqrt(dabs(bw_sq(yRnd(13),mass(5,1), mass(5,2), (inv_mass(3)-inv_mass(4))**2, jacobian5)))
!print *, bw_sq(yRnd(13),mass(5,1), mass(5,2), (inv_mass(3)-inv_mass(4))**2, jacobian5), inv_mass(4:5)
!      endif


!444444444444
!energy of 4 in the CM frame of 3
      MomExt(1,4)=(inv_mass(3)**2+(inv_mass(4)+inv_mass(5))*(inv_mass(4)-inv_mass(5)))/2d0/inv_mass(3)
!|3-momentum| of 4 in the CM frame of 3
      cm_abs3p(4) = dsqrt((MomExt(1,4)+inv_mass(4)) * (MomExt(1,4)-inv_mass(4)))
!generating cos(theta_4) and phi_4 in the CM frame of 3
      cm_cos_theta(4) = yRnd(6)
      cm_cos_theta(4) = cm_cos_theta(4)*2d0-1d0
      cm_sin_theta(4) = dsqrt((1d0+cm_cos_theta(4))  *(1d0-cm_cos_theta(4)))
      phi = yRnd(7)
      phi=Twopi*phi
      cm_cos_phi(4) = dcos(phi)
      cm_sin_phi(4) = dsin(phi)
!3-momentum of 4 in the CM frame of 3
      MomExt(2,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_cos_phi(4)
      MomExt(3,4)=cm_abs3p(4)*cm_sin_theta(4)*cm_sin_phi(4)
      MomExt(4,4)=cm_abs3p(4)*cm_cos_theta(4)
!boost the 4-momentum of 4 to the lab frame
!x and y components stay the same
!z component
      beta = MomExt(4,3)/MomExt(1,3)
      gamma = 1d0/dsqrt((1d0+beta)*(1d0-beta))
      MomDummy(1:4)=MomExt(1:4,4)
      MomExt(4,4)=(MomDummy(4)+MomDummy(1)*beta) *gamma
!energy 
      MomExt(1,4)=(MomDummy(1)+MomDummy(4)*beta) *gamma
!555555555555555555
      MomExt(1:4,5) = MomExt(1:4,3) - MomExt(1:4,4)

!666666666666666666
!invariant mass of 6
      inv_mass(6)=0d0
!energy of 6 in the CM frame of 4
      MomExt(1,6)=inv_mass(4)/2d0
!|3-momentum| of 6 in the CM frame of 4
      cm_abs3p(6)=MomExt(1,6)
!generating cos(theta_6) and phi_6 in the CM frame of 4
!z-axis is along the boost of 2
      cm_cos_theta(6) = yRnd(8)
      cm_cos_theta(6) = cm_cos_theta(6)*2d0-1d0
      cm_sin_theta(6) = dsqrt((1d0+cm_cos_theta(6)) *(1d0-cm_cos_theta(6)))
      cm_cos_phi(6) = yRnd(9)
      phi=Twopi*cm_cos_phi(6)
      cm_cos_phi(6) = dcos(phi)
      cm_sin_phi(6) = dsin(phi)
!3-momentum of 6 in the CM frame of 4
      MomExt(2,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_cos_phi(6)
      MomExt(3,6)=cm_abs3p(6)*cm_sin_theta(6)*cm_sin_phi(6)
      MomExt(4,6)=cm_abs3p(6)*cm_cos_theta(6)
!boost the 4-momentum of 6 to the lab frame
      temp_vector = MomExt(1:4,6)
      temp_boost = MomExt(1:4,4)

      call LORENTZ(temp_vector, temp_boost)

      MomExt(1:4,6) = temp_vector
!7777777777777777777777
!invariant mass of 7
      inv_mass(7)=0d0
!4-momentum of 7 (lab frame) by energy-momentum conservation
      MomExt(1:4,7)=MomExt(1:4,4)-MomExt(1:4,6)
!8888888888888888888888
!invariant mass of 8
      inv_mass(8)=0d0
!energy of 8 in the CM frame of 5
      MomExt(1,8)=inv_mass(5)/2d0
!|3-momentum| of 8 in the CM frame of 5
      cm_abs3p(8)=MomExt(1,8)
!generating cos(theta_8) and phi_8 in the CM frame of 5
!z-axis is along the boost of 5
      cm_cos_theta(8) = yRnd(10)
      cm_cos_theta(8) = cm_cos_theta(8)*2d0-1d0
      cm_sin_theta(8) = dsqrt((1d0+cm_cos_theta(8)) *(1d0-cm_cos_theta(8)))
      cm_cos_phi(8) = yRnd(11)
      phi=Twopi*cm_cos_phi(8)
      cm_cos_phi(8) = dcos(phi)
      cm_sin_phi(8) = dsin(phi)
!3-momentum of 8 in the CM frame of 5
!x and y components are not necessary, yet
      MomExt(2,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_cos_phi(8)
      MomExt(3,8)=cm_abs3p(8)*cm_sin_theta(8)*cm_sin_phi(8)
      MomExt(4,8)=cm_abs3p(8)*cm_cos_theta(8)
!boost the 4-momentum of 6 to the lab frame
      temp_vector = MomExt(1:4,8)
      temp_boost = MomExt(1:4,5)

      call LORENTZ(temp_vector, temp_boost)

      MomExt(1:4,8) = temp_vector
!9999999999999999999999
!invariant mass of 9
      inv_mass(9)=0d0     
!4-momentum of 9 (lab frame) by energy-momentum conservation
      MomExt(1:4,9)=MomExt(1:4,5)-MomExt(1:4,8)

      PSWgt = jacobian4*jacobian5*cm_abs3p(4)/(4d0*pi)/inv_mass(3)/(8d0*pi)**2/(2*pi)**2
      !print *,  "()",inv_mass, jacobian4, jacobian5, cm_abs3p(4), PSWgt, "()"

!do i=4,7
!print *, dsqrt(dabs(four_momentum(i,:).dot.four_momentum(i,:)))
!enddo
!pause

RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_2to3M(EHat,Mass,xRndPS,Mom,PSWgt,Width)
use ModParameters
implicit none
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,Mass(1:3)
real(8), optional :: Width(1:3)
real(8) :: xRndPS(1:5),xRndWidth(1:3),BW_Jacobi(1:3)
real(8) :: Mom(1:4,1:5),TmpMom(1:4),BW_Mass(1:3)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)

  if( present(Width) ) then
      call random_number(xRndWidth)
      
      call SmearExternal(xRndWidth(1),Mass(1),Width(1),Mass(1)-6d0*Width(1),Mass(1)+6d0*Width(1),BW_Mass(1),BW_Jacobi(1))
      call SmearExternal(xRndWidth(2),Mass(2),Width(2),Mass(2)-6d0*Width(2),Mass(2)+6d0*Width(2),BW_Mass(2),BW_Jacobi(2))
      call SmearExternal(xRndWidth(3),Mass(3),Width(3),Mass(3)-6d0*Width(3),Mass(3)+6d0*Width(3),BW_Mass(3),BW_Jacobi(3))
      
      call genps(3,Ehat,xRndPS(1:5),BW_Mass,Mom(1:4,3:5),PSWgt)
      PSWgt = PSWgt*PiWgtPr                             &
              * BW_Jacobi(1)*BW_Jacobi(2)*BW_Jacobi(3)  &  ! maybe add factors for higgs  
              * (2d0*Ga_Top*m_Top)**2                   &  ! remove narrow-width prefactor 
              * (Ga_Top**2 * m_Top**2)**2                  ! and replace by on-shell propagator to restore correct scaling                   
                   
  else
      call genps(3,Ehat,xRndPS(1:5),Mass,Mom(1:4,3:5),PSWgt)
      PSWgt = PSWgt*PiWgtPr
  endif


!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,Mass/),Mom(1:4,3:5),PSWgt)
!   TmpMom(1:4) = Mom(1:4,3)
!   Mom(1:4,3)  = Mom(1:4,5)
!   Mom(1:4,5)  = TmpMom(1:4)

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE TTbar_OnShellProjection(MomIn,MomOut)
use modParameters
use modMisc
implicit none
real(8) :: MomIn(:,:),MomOut(:,:),MomTmp(1:4)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13


    call ShiftMass(MomIn(1:4,tbar),MomIn(1:4,t),m_top,m_top,MomOut(1:4,tbar),MomOut(1:4,t))! project top momenta on-shell
    
    MomTmp(1:4) = MomOut(1:4,t) - MomIn(1:4,Wp)
    call ShiftMass(MomTmp,MomIn(1:4,Wp),m_Bot,m_W,MomOut(1:4,b),MomOut(1:4,Wp))! project b and W+ on-shell
    
    MomTmp(1:4) = MomOut(1:4,tbar) - MomIn(1:4,Wm)
    call ShiftMass(MomTmp,MomIn(1:4,Wm),m_Bot,m_W,MomOut(1:4,bbar),MomOut(1:4,Wm))! project bbar and W- on-shell
    
    MomTmp(1:4) = MomOut(1:4,Wp) - MomIn(1:4,lepP)                                ! project W+ decay products on-shell
    MomOut(1:4,nu)   = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP)) * MomIn(1:4,lepP)
    MomOut(1:4,lepP) = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP))) * MomIn(1:4,lepP)
    
    MomTmp(1:4) = MomOut(1:4,Wm) - MomIn(1:4,lepM)                                ! project W- decay products on-shell
    MomOut(1:4,nubar) = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM)) * MomIn(1:4,lepM)
    MomOut(1:4,lepM)  = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM))) * MomIn(1:4,lepM)

return
END SUBROUTINE




SUBROUTINE TTbar_OffShellProjection(MomIn,MomOut,Jacobian)
use modParameters
use modMisc
implicit none
real(8) :: MomIn(:,:),MomOut(:,:),MomTmp(1:4),Jacobian
real(8) :: xRndWidth(2:5),BW_Mass(2:5),BW_Jacobi(2:5)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13


    call random_number(xRndWidth)
      
    call SmearExternal(xRndWidth(2),m_top,Ga_Top,m_top-6d0*Ga_Top,m_top+6d0*Ga_Top,BW_Mass(2),BW_Jacobi(2))
    call SmearExternal(xRndWidth(3),m_top,Ga_Top,m_top-6d0*Ga_Top,m_top+6d0*Ga_Top,BW_Mass(3),BW_Jacobi(3))
    call SmearExternal(xRndWidth(4),m_W,Ga_W,m_W-6d0*Ga_W,m_W+6d0*Ga_W,BW_Mass(4),BW_Jacobi(4))
    call SmearExternal(xRndWidth(5),m_W,Ga_W,m_W-6d0*Ga_W,m_W+6d0*Ga_W,BW_Mass(5),BW_Jacobi(5))
    Jacobian = BW_Jacobi(2) * BW_Jacobi(3) * BW_Jacobi(4) * BW_Jacobi(5)    

! print *, "smeared mt",(BW_Mass(2:3)-m_top)*100d0
! print *, "smeared mw",(BW_Mass(4:5)-m_w)*100d0

    call ShiftMass(MomIn(1:4,tbar),MomIn(1:4,t),BW_Mass(2),BW_Mass(3),MomOut(1:4,tbar),MomOut(1:4,t))
    
    MomTmp(1:4) = MomOut(1:4,t) - MomIn(1:4,Wp)
    call ShiftMass(MomTmp,MomIn(1:4,Wp),m_Bot,BW_Mass(4),MomOut(1:4,b),MomOut(1:4,Wp))! project b and W+ on-shell
    
    MomTmp(1:4) = MomOut(1:4,tbar) - MomIn(1:4,Wm)
    call ShiftMass(MomTmp,MomIn(1:4,Wm),m_Bot,BW_Mass(5),MomOut(1:4,bbar),MomOut(1:4,Wm))! project bbar and W- on-shell
    
    MomTmp(1:4) = MomOut(1:4,Wp) - MomIn(1:4,lepP)                                ! project W+ decay products on-shell
    MomOut(1:4,nu)   = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP)) * MomIn(1:4,lepP)
    MomOut(1:4,lepP) = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP))) * MomIn(1:4,lepP)
    
    MomTmp(1:4) = MomOut(1:4,Wm) - MomIn(1:4,lepM)                                ! project W- decay products on-shell
    MomOut(1:4,nubar) = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM)) * MomIn(1:4,lepM)
    MomOut(1:4,lepM)  = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM))) * MomIn(1:4,lepM)

return
END SUBROUTINE




SUBROUTINE Top_OffShellProjection(MomIn,MomOut,Jacobian)
use modParameters
use modMisc
implicit none
real(8) :: MomIn(:,:),MomOut(:,:),MomTmp(1:4),Jacobian
real(8) :: xRndWidth(1:2),BW_Mass(1:2),BW_Jacobi(1:2)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9


    call random_number(xRndWidth)
    
    call SmearExternal(xRndWidth(1),m_top,Ga_Top,m_top-6*Ga_Top,m_top+6*Ga_Top,BW_Mass(1),BW_Jacobi(1))
    call SmearExternal(xRndWidth(2),m_W,Ga_W,m_W-6d0*Ga_W,m_W+6d0*Ga_W,BW_Mass(2),BW_Jacobi(2))
    Jacobian = BW_Jacobi(1) * BW_Jacobi(2)

! print *, "smeared mt",BW_Mass(1)
! print *, "smeared mw",BW_Mass(2)

    call ShiftMass(MomIn(1:4,t),MomIn(1:4,qout),BW_Mass(1),0d0,MomOut(1:4,t),MomOut(1:4,qout))
    
    MomTmp(1:4) = MomOut(1:4,t) - MomIn(1:4,W)
    call ShiftMass(MomTmp,MomIn(1:4,W),m_Bot,BW_Mass(2),MomOut(1:4,b),MomOut(1:4,W))! project b and W on-shell
    
    MomTmp(1:4) = MomOut(1:4,W) - MomIn(1:4,lep)                                ! project W decay products on-shell
    MomOut(1:4,nu)   = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lep)) * MomIn(1:4,lep)
    MomOut(1:4,lep) = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lep))) * MomIn(1:4,lep)

return
END SUBROUTINE




SUBROUTINE EvalPhasespace_2to3ArbMass(EHat,Mass,xRndPS,Mom,PSWgt)
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,Mass(1:3)
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)                                                                                                              
! integer :: NPart,i                                                                                                              
! real(8) :: vel,parx,theta ! for checks                                                              
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


  call genps(3,Ehat,xRndPS(1:5),Mass,Mom(1:4,3:5),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,Mass/),Mom(1:4,3:5),PSWgt)         
!   TmpMom(1:4) = Mom(1:4,3)                                                                                         
!   Mom(1:4,3)  = Mom(1:4,5)                                                                                             
!   Mom(1:4,5)  = TmpMom(1:4)                                                                                            

!  particles on the beam axis:                                                                             
                                    
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE EvalPhasespace_TopDecay(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)

!     MomDK(1:4,i): i= 1:bottom, 2:lepton, 3:neutrino
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      call genps(2,m_W,xRndPS(3:4),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      PSWgt = PSWgt2*PiWgt2 * PSWgt3*PiWgt2

RETURN
END SUBROUTINE
SUBROUTINE EvalPhasespace_tautau(xRnd,pHiggs,MY_IDUP,Mom,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xRnd(:), pHiggs(:), Mom(:,:)
integer :: MY_IDUP(1:2)
real(8) :: Jac,Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9
real(8) :: Minvsq_tau1,Minvsq_tau2,Minvsq_Wp,Minvsq_Wm
real(8),parameter :: m_nu = 0d0, m_lep=0d0
integer, parameter :: inLeft=1, inRight=2, Hig=3, tauP=4, tauM=5, Wp=6, Wm=7,   nu=8, nubar_tau=9, lepP=10,   lepM=11, nu_tau=12, nubar=13

  
   
   
!  H-->tau tau (NWA)
   Jac1 = s_channel_prop_decay(pHiggs,(/m_tau,ga_tau,m_tau,m_tau/),(/m_tau,ga_tau,m_tau,m_tau/),xRnd(1:2),Mom(:,tauP),Mom(:,tauM)) 

!  tau-->W nu (BW)
   Jac2 = s_channel_prop_decay(Mom(:,tauP),(/m_W,ga_W,0d0,m_tau/),(/m_nu,0d0,0d0,0d0/),xRnd(3:5),Mom(:,Wp),Mom(:,nubar_tau)) 
   Jac3 = s_channel_prop_decay(Mom(:,tauM),(/m_W,ga_W,0d0,m_tau/),(/m_nu,0d0,0d0,0d0/),xRnd(6:8),Mom(:,Wm),Mom(:,nu_tau)) 

!  W-->l nu (ONSH)
   Jac4 = s_channel_prop_decay(Mom(:,Wp),(/m_Lep,0d0,0d0,0d0/),(/m_nu,0d0,0d0,0d0/),xRnd( 9:10),Mom(:,LepP),Mom(:,nu)) 
   Jac5 = s_channel_prop_decay(Mom(:,Wm),(/m_Lep,0d0,0d0,0d0/),(/m_nu,0d0,0d0,0d0/),xRnd(11:12),Mom(:,LepM),Mom(:,nubar)) 
   
   Jac = Jac1*Jac2*Jac4*Jac3*Jac5 * PSNorm6
   
   
!    print *, "OS checker",dsqrt(pHiggs.dot.pHiggs )
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,tauP).dot.Mom(1:4,tauP) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,tauM).dot.Mom(1:4,tauM) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,Wp).dot.Mom(1:4,Wp) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,Wm).dot.Mom(1:4,Wm) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,nu).dot.Mom(1:4,nu) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,nu_tau).dot.Mom(1:4,nu_tau) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,nubar_tau).dot.Mom(1:4,nubar_tau) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,nubar).dot.Mom(1:4,nubar) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,lepP).dot.Mom(1:4,lepP) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,lepM).dot.Mom(1:4,lepM) ))
!    print *, "----------"
!    print *, "Mom.cons. ",pHiggs(1:4)-Mom(1:4,tauP)-Mom(1:4,tauM)
!    print *, "Mom.cons. ",Mom(1:4,tauP) - Mom(1:4,nubar_tau)-Mom(1:4,Wp)
!    print *, "Mom.cons. ",Mom(1:4,tauM) - Mom(1:4,nu_tau)-Mom(1:4,Wm)
!    print *, "Mom.cons. ",Mom(1:4,Wp) - Mom(1:4,nu)-Mom(1:4,lepP)
!    print *, "Mom.cons. ",Mom(1:4,Wm) - Mom(1:4,nubar)-Mom(1:4,lepM)   
!    pause
   

RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_VBF_NEW2(xchannel,xRnd,Energy,Mom,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xchannel,xRnd(:), Energy, Mom(:,:)
integer :: iChannel
real(8) :: Jac,Jac1,Jac2,Jac3,Mom_ij_Dummy(1:4),s35,s45,minmax(1:2)
integer, parameter :: NumChannels=4
integer, parameter :: inLeft=1, inRight=2, qup=3, qdn=4, Higgs=5


   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

  
   iChannel = int(xchannel * NumChannels)+1

IF( iChannel.EQ.1 ) THEN   
   
   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                           !  int d(s35)  
!  equival.: call get_minmax_s(Energy**2,0d0,M_Reso**2,0d0,minmax); Jac1 = k_l(xRnd(5),minmax(1),minmax(2),s35)
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),M_W**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,4))              !  1+2 --> (35)+4
   Jac3 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2)-Mom(1:4,4),M_W**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,3),Mom(1:4,5))    !  1+(24) --> 3+5
   Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                        !  combine   
   
ELSEIF( iChannel.EQ.2 ) THEN   

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                           !  int d(s35)  
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),M_Z**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,4))              !  1+2 --> (35)+4 
   Jac3 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2)-Mom(1:4,4),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,3),Mom(1:4,5))    !  1+(24) --> 3+5
   Jac  = Jac1*Jac2*Jac3 * PSNorm3                                                                                       !  combine  
   
ELSEIF( iChannel.EQ.3 ) THEN   

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                           !  int d(s45)  
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),M_W**2,s45,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,3))              !  1+2 --> (45)+3 
   Jac3 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2)-Mom(1:4,3),M_W**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,4),Mom(1:4,5))    !  1+(23) --> 4+5
   Jac  = Jac1*Jac2*Jac3 * PSNorm3                                                                                       !  combine  
   
ELSEIF( iChannel.EQ.4 ) THEN   

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                           !  int d(s45)  
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),M_Z**2,s45,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,3))              !  1+2 --> (45)+3 
   Jac3 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2)-Mom(1:4,3),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,4),Mom(1:4,5))    !  1+(23) --> 4+5
   Jac  = Jac1*Jac2*Jac3 * PSNorm3                                                                                       !  combine  
   
ELSE
   print *, xchannel , NumChannels
   call Error("PS channel not available in EvalPhasespace_VBF_NEW2",ichannel)
ENDIF 
      
      
      if( IsNaN(Jac) ) then! THIS SHOULD BE REMOVED AFTER DEBUGGING
         Jac = 0d0
         print *, "ERROR in EvalPhasespace_VBF_NEW2, NaN Jac",Energy,xchannel,xRnd
         print *, "ERROR in Channel",ichannel
      endif

      
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,3).dot.Mom(1:4,3) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,4).dot.Mom(1:4,4) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,5).dot.Mom(1:4,5) ))
!    print *, "----------"
!    print *, "Mom.cons. ",Mom(1:4,1)+Mom(1:4,2)-Mom(1:4,3)-Mom(1:4,4)-Mom(1:4,5)
!    pause
   
   
   
RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_VBF_deterministic(xchannel,xRnd,Energy,iSel,jSel,rSel,sSel,Mom,Jac) ! TESTING PURPOSES
use ModParameters
use ModPhasespace
use ModMisc
implicit none
integer, intent(in) :: iSel,jSel,rSel,sSel
real(8) :: xchannel, xRnd(:), Energy, Mom(:,:)
integer :: jz1, jz2, jw1, jw2 ! incoming partons
integer :: kz1, kz2, kw1, kw2 ! outgoing partons
logical :: ZZ_fusion,WW_fusion
real(8) :: Jac,Jac1,Jac2,Jac3,Mom_ij_Dummy(1:4),s35,s45

   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

	ZZ_fusion=.false.
	WW_fusion=.false.

	jz1 = 1
	jz2 = 2
	jw1 = jz1
	jw2 = jz2

	if( &
	(iSel.eq.rSel .and. jSel.eq.sSel) &
	.or. &
	(iSel.eq.sSel .and. jSel.eq.rSel) &
	) then
	   ZZ_fusion=.true.
	   if( iSel.eq.sSel .and. jSel.eq.rSel ) then
	      kz1 = 4
	      kz2 = 3
	   else
	      kz1 = 3
		  kz2 = 4
	   endif
	endif

	if( &
	( (sign(iSel,rSel).eq.iSel .and. sign(jSel,sSel).eq.jSel) .and. (abs(iSel-rSel).eq.1 .or. abs(iSel-rSel).eq.3 .or. abs(iSel-rSel).eq.5) .and. (abs(jSel-sSel).eq.1 .or. abs(jSel-sSel).eq.3 .or. abs(jSel-sSel).eq.5) ) &
	.or. &
	( (sign(iSel,sSel).eq.iSel .and. sign(jSel,rSel).eq.jSel) .and. (abs(iSel-sSel).eq.1 .or. abs(iSel-sSel).eq.3 .or. abs(iSel-sSel).eq.5) .and. (abs(jSel-rSel).eq.1 .or. abs(jSel-rSel).eq.3 .or. abs(jSel-rSel).eq.5) ) &
	) then
	   WW_fusion=.true.
	   ! W_is W_jr fusion
	   if( (sign(iSel,sSel).eq.iSel .and. sign(jSel,rSel).eq.jSel) .and. (abs(iSel-sSel).eq.1 .or. abs(iSel-sSel).eq.3 .or. abs(iSel-sSel).eq.5) .and. (abs(jSel-rSel).eq.1 .or. abs(jSel-rSel).eq.3 .or. abs(jSel-rSel).eq.5) ) then
	      kw1 = 4
	      kw2 = 3
	   else	! W_ir W_js fusion
	      kw1 = 3
	      kw2 = 4	      
	   endif
	endif

   if(WW_fusion .and. .not.ZZ_fusion) then
      if(xchannel .lt. 0.5) then
         Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                                   !  int d(s35)  
         Jac2 = t_channel_prop_decay(Mom(1:4,jw1),Mom(1:4,jw2),M_W**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,kw2))                !  1+2 --> (35)+4
         Jac3 = t_channel_prop_decay(Mom(1:4,jw1),Mom(1:4,jw2)-Mom(1:4,kw2),M_W**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kw1),Mom(1:4,5))    !  1+(24) --> 3+5
         Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
      else
         Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                                   !  int d(s45)  
         Jac2 = t_channel_prop_decay(Mom(1:4,jw1),Mom(1:4,jw2),M_W**2,0d0,s45,xRnd(2:3),Mom(1:4,kw1),Mom_ij_Dummy(1:4))                !  1+2 --> 3+(45) 
         Jac3 = t_channel_prop_decay(Mom(1:4,jw1)-Mom(1:4,kw1),Mom(1:4,jw2),M_W**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kw2),Mom(1:4,5))    !  (13)+2 --> 4+5
         Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
      endif
   elseif(ZZ_fusion .and. .not.WW_fusion) then
      if(iSel .eq.jSel) then
         if(xchannel .lt. 0.25) then
            Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                                   !  int d(s35)  
            Jac2 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2),M_Z**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,kz2))                !  1+2 --> (35)+4
            Jac3 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2)-Mom(1:4,kz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz1),Mom(1:4,5))    !  1+(24) --> 3+5
            Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
         elseif(xchannel .lt. 0.5) then
            Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                                   !  int d(s45)  
            Jac2 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2),M_Z**2,0d0,s45,xRnd(2:3),Mom(1:4,kz1),Mom_ij_Dummy(1:4))                !  1+2 --> 3+(45) 
            Jac3 = t_channel_prop_decay(Mom(1:4,jz1)-Mom(1:4,kz1),Mom(1:4,jz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz2),Mom(1:4,5))    !  (13)+2 --> 4+5
            Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
         elseif(xchannel .lt. 0.75) then
            Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                                   !  int d(s35)  
            Jac2 = t_channel_prop_decay(Mom(1:4,jz2),Mom(1:4,jz1),M_Z**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,kz2))                !  2+1 --> (35)+4
            Jac3 = t_channel_prop_decay(Mom(1:4,jz2),Mom(1:4,jz1)-Mom(1:4,kz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz1),Mom(1:4,5))    !  2+(14) --> 3+5
            Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
         else
            Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                                   !  int d(s45)  
            Jac2 = t_channel_prop_decay(Mom(1:4,jz2),Mom(1:4,jz1),M_Z**2,0d0,s45,xRnd(2:3),Mom(1:4,kz1),Mom_ij_Dummy(1:4))                !  2+1 --> 3+(45) 
            Jac3 = t_channel_prop_decay(Mom(1:4,jz2)-Mom(1:4,kz1),Mom(1:4,jz1),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz2),Mom(1:4,5))    !  (23)+1 --> 4+5
            Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
         endif
	  else
         if(xchannel .lt. 0.5) then
            Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                                   !  int d(s35)  
            Jac2 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2),M_Z**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,kz2))                !  1+2 --> (35)+4
            Jac3 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2)-Mom(1:4,kz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz1),Mom(1:4,5))    !  1+(24) --> 3+5
            Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
         else
            Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                                   !  int d(s45)  
            Jac2 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2),M_Z**2,0d0,s45,xRnd(2:3),Mom(1:4,kz1),Mom_ij_Dummy(1:4))                !  1+2 --> 3+(45) 
            Jac3 = t_channel_prop_decay(Mom(1:4,jz1)-Mom(1:4,kz1),Mom(1:4,jz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz2),Mom(1:4,5))    !  (13)+2 --> 4+5
            Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
         endif
	  endif
   elseif(ZZ_fusion .and. WW_fusion) then
      if(xchannel .lt. 0.25) then
         Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                                   !  int d(s35)  
         Jac2 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2),M_Z**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,kz2))                !  1+2 --> (35)+4
         Jac3 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2)-Mom(1:4,kz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz1),Mom(1:4,5))    !  1+(24) --> 3+5
         Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
      elseif(xchannel .lt. 0.5) then
         Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                                   !  int d(s45)  
         Jac2 = t_channel_prop_decay(Mom(1:4,jz1),Mom(1:4,jz2),M_Z**2,0d0,s45,xRnd(2:3),Mom(1:4,kz1),Mom_ij_Dummy(1:4))                !  1+2 --> 3+(45) 
         Jac3 = t_channel_prop_decay(Mom(1:4,jz1)-Mom(1:4,kz1),Mom(1:4,jz2),M_Z**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kz2),Mom(1:4,5))    !  (13)+2 --> 4+5
         Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
      elseif(xchannel .lt. 0.75) then
         Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                                   !  int d(s35)  
         Jac2 = t_channel_prop_decay(Mom(1:4,jw1),Mom(1:4,jw2),M_W**2,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,kw2))                !  1+2 --> (35)+4
         Jac3 = t_channel_prop_decay(Mom(1:4,jw1),Mom(1:4,jw2)-Mom(1:4,kw2),M_W**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kw1),Mom(1:4,5))    !  1+(24) --> 3+5
         Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
      else
         Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                                   !  int d(s45)  
         Jac2 = t_channel_prop_decay(Mom(1:4,jw1),Mom(1:4,jw2),M_W**2,0d0,s45,xRnd(2:3),Mom(1:4,kw1),Mom_ij_Dummy(1:4))                !  1+2 --> 3+(45) 
         Jac3 = t_channel_prop_decay(Mom(1:4,jw1)-Mom(1:4,kw1),Mom(1:4,jw2),M_W**2,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,kw2),Mom(1:4,5))    !  (13)+2 --> 4+5
         Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                                !  combine   
      endif

   endif

      if( IsNaN(Jac) ) then! THIS SHOULD BE REMOVED AFTER DEBUGGING
         Jac = 0d0
         print *, "ERROR in EvalPhasespace_VBF_deterministic, NaN Jac",Energy,xchannel,xRnd
      endif
   
RETURN
END SUBROUTINE




SUBROUTINE EvalPhasespace_VBF_H4f(xchannel,xRnd,Energy,Mom,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xchannel,xRnd(:), Energy, Mom(:,:)
integer :: iChannel
real(8) :: Jac,Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9
real(8) :: s3H,s4H,s56,s78,s910,Mom_Dummy(1:4),xRndOffShellZ
real(8), parameter :: RescaleWidth=10d0
integer, parameter :: NumChannels=4
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10


   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)
   
   iChannel = int(xchannel * NumChannels -1d-10)+1
!    print *, "PS channel ",iChannel
   
   
   
IF( iChannel.EQ.1 ) THEN! 34 + WW-->H-->ZZ

!  masses
   if( m4l_minmax(1).lt.0d0 ) then
      Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
   else
      Jac1 = k_l(xRnd(1),m4l_minmax(1)**2,min(Energy**2,m4l_minmax(2)**2),s56)                                            !  int d(s56)    = linear mapping
   endif
   Jac2 = k_l(xRnd(2),s56,Energy**2,s3H)                                                                                          !  int d(s3H)     
   Jac3 = s_channel_propagator(M_Z**2,Ga_Z,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
   Jac4 = s_channel_propagator(M_Z**2,Ga_Z,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2

!  splittings
   Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_W**2,s3H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,4))                                !  1+2 --> (3H)+4
   Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,4),M_W**2,0d0,s56,xRnd(7:8),Mom(:,3),Mom_Dummy(1:4))                       !  1+(24) --> 3+H
   Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6       
   Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8       
   Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10      
      
      
ELSEIF( iChannel.EQ.3 ) THEN! 34 + ZZ-->H-->ZZ

!  masses
   if( m4l_minmax(1).lt.0d0 ) then
      Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
   else
      Jac1 = k_l(xRnd(1),m4l_minmax(1)**2,min(Energy**2,m4l_minmax(2)**2),s56)                                            !  int d(s56)    = linear mapping
   endif    
   Jac2 = k_l(xRnd(2),s56,Energy**2,s3H)                                                                                          !  int d(s3H)     
   Jac3 = s_channel_propagator(M_Z**2,Ga_Z,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
   Jac4 = s_channel_propagator(M_Z**2,Ga_Z,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2
     
!  splittings
   Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_Z**2,s3H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,4))                                !  1+2 --> (3H)+4
   Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,4),M_Z**2,0d0,s56,xRnd(7:8),Mom(:,3),Mom_Dummy(1:4))                       !  1+(24) --> 3+H
   Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6       
   Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8       
   Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10      
   

ELSEIF( iChannel.EQ.2 ) THEN! 43 + WW-->H-->ZZ   

!  masses
   if( m4l_minmax(1).lt.0d0 ) then
      Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
   else
      Jac1 = k_l(xRnd(1),m4l_minmax(1)**2,min(Energy**2,m4l_minmax(2)**2),s56)                                            !  int d(s56)    = linear mapping
   endif
   Jac2 = k_l(xRnd(2),s56,Energy**2,s4H)                                                                                          !  int d(s4H)     
   Jac3 = s_channel_propagator(M_Z**2,Ga_Z,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
   Jac4 = s_channel_propagator(M_Z**2,Ga_Z,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2
     
!  splittings
   Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_W**2,s4H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,3))                                !  1+2 --> (4H)+3
   Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,3),M_W**2,0d0,s56,xRnd(7:8),Mom(:,4),Mom_Dummy(1:4))                       !  1+(23) --> 4+H
   Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6       
   Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8       
   Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10      
   

ELSEIF( iChannel.EQ.4 ) THEN! 43 + ZZ-->H-->ZZ   

!  masses
   if( m4l_minmax(1).lt.0d0 ) then
      Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
   else
      Jac1 = k_l(xRnd(1),m4l_minmax(1)**2,min(Energy**2,m4l_minmax(2)**2),s56)                                            !  int d(s56)    = linear mapping
   endif  
   Jac2 = k_l(xRnd(2),s56,Energy**2,s4H)                                                                                          !  int d(s4H)     
   Jac3 = s_channel_propagator(M_Z**2,Ga_Z,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
   Jac4 = s_channel_propagator(M_Z**2,Ga_Z,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2

!  splittings
   Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_Z**2,s4H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,3))                                !  1+2 --> (4H)+3
   Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,3),M_Z**2,0d0,s56,xRnd(7:8),Mom(:,4),Mom_Dummy(1:4))                       !  1+(23) --> 4+H
   Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6       
   Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8       
   Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10      
   
ELSE
   call Error("PS channel not available in EvalPhasespace_VBF_H4f",ichannel)
ENDIF 

   
!    call random_number(xRndOffShellZ)   ! switching this off for test purposes
!    if( xRndOffShellZ.gt.0.5d0 ) then
!         call swap_mom(Mom(:,5),Mom(:,6))
!         call swap_mom(Mom(:,7),Mom(:,9))
!         call swap_mom(Mom(:,8),Mom(:,10))
!    endif
   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6*Jac7*Jac8*Jac9 * PSNorm6                                                                   !  combine   




   if( isNan(jac) ) then
      print *, "EvalPhasespace_VBF_H4f NaN"
      print *, Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9,ichannel
      Jac = 0d0
   endif

!    print *, "OS checker", dsqrt( dabs(Mom(1:4,3).dot.Mom(1:4,3) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,4).dot.Mom(1:4,4) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,7).dot.Mom(1:4,7) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,8).dot.Mom(1:4,8) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,9).dot.Mom(1:4,9) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,10).dot.Mom(1:4,10) ))
!    print *, "----------"
!    print *, "Mom.cons. ",Mom(1:4,1)+Mom(1:4,2)-Mom(1:4,3)-Mom(1:4,4)-Mom(1:4,7)-Mom(1:4,8)-Mom(1:4,9)-Mom(1:4,10)
!    print *, "Mom.cons. ",Mom_Dummy(1:4)-Mom(1:4,5)-Mom(1:4,6)
!    print *, "Mom.cons. ",Mom(1:4,5)-Mom(1:4,7)-Mom(1:4,8)
!    print *, "Mom.cons. ",Mom(1:4,6)-Mom(1:4,9)-Mom(1:4,10)
!    print *, "----------"
!    print *, "Inv.mass  ",get_MInv(Mom_Dummy(1:4))*100d0
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,5))*100d0
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,6))*100d0
!    pause
   
   
   
RETURN
END SUBROUTINE




SUBROUTINE EvalPhasespace_H4f(xchannel,xRnd,Energy,Mom,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xchannel,xRnd(1:8), Energy, Mom(4,8)
integer :: iChannel
real(8) :: Jac,Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9
real(8) :: s3H,s4H,s56,s78,s910,Mom_Dummy(1:4)
real(8), parameter :: RescaleWidth=10d0
integer :: NumChannels=   4


   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)
   
   if( .not. includeInterference ) NumChannels=2
   iChannel = int(xchannel * NumChannels -1d-10)+1
!    print *, "PS channel ",iChannel

! if( iChannel.eq.2 ) iChannel=1! turn this on only when mX > 2mV
! if( iChannel.eq.4 ) iChannel=3
   
   
   
IF( iChannel.EQ.1 ) THEN


!  masses
   Jac1=1d0; s56=Energy**2   
! print *, "entering Jac2"
   Jac2 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(1),s78)                                                                   !  int d(s78)    = Z1
! print *, "entering Jac3"
   Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)                                           !  int d(s910) = Z2
   
   
!    print *, "x",xrnd(1:2)
!    print *, "s",s56,s78,s910;pause
   
!  splittings
   Mom_Dummy(1:4) = (/Energy,0d0,0d0,0d0/)
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,3),Mom(:,4))                                                   !  H --> 5+6       
   Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,5),Mom(:,6))                                                         !  5 --> 7+8       
   Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,7),Mom(:,8))                                                        !  6 --> 9+10      
   
   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6 * PSNorm4  !*NumChannels                                                         !  combine   


ELSEIF( iChannel.EQ.2 ) THEN

!  masses
   Jac1=1d0; s56=Energy**2   
   Jac2 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(1),s78)                                                                   !  int d(s78)    = Z1
   Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)                                           !  int d(s910) = Z2

!  splittings
   Mom_Dummy(1:4) = (/Energy,0d0,0d0,0d0/)
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,4),Mom(:,3))                                                   !  H --> 5+6       
   Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,5),Mom(:,6))                                                         !  5 --> 7+8       
   Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,7),Mom(:,8))                                                        !  6 --> 9+10      
      
   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6 * PSNorm4  !*NumChannels                                                         !  combine   

   
   
ELSEIF( iChannel.EQ.3 ) THEN

!  masses
   Jac1=1d0; s56=Energy**2   
   Jac2 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(1),s78)                                                                   !  int d(s78)    = Z1
   Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)                                           !  int d(s910) = Z2

!  splittings
   Mom_Dummy(1:4) = (/Energy,0d0,0d0,0d0/)
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,3),Mom(:,4))                                                   !  H --> 5+6       
   Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,7),Mom(:,6))                                                         !  5 --> 7+8       
   Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,5),Mom(:,8))                                                        !  6 --> 9+10      
      
   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6 * PSNorm4  !*NumChannels                                                         !  combine   


ELSEIF( iChannel.EQ.4 ) THEN

!  masses
   Jac1=1d0; s56=Energy**2   
   Jac2 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(1),s78)                                                                   !  int d(s78)    = Z1
   Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)                                           !  int d(s910) = Z2

!  splittings
   Mom_Dummy(1:4) = (/Energy,0d0,0d0,0d0/)
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,4),Mom(:,3))                                                   !  H --> 5+6       
   Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,7),Mom(:,6))                                                         !  5 --> 7+8       
   Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,5),Mom(:,8))                                                        !  6 --> 9+10      
      
   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6 * PSNorm4  !*NumChannels                                                         !  combine   
   
 
ENDIF

! print *, energy,dsqrt(s56);pause


   if( isNan(jac) ) then
      print *, "EvalPhasespace_H4f NaN"
      print *, Energy
      print *, s56,s78,s910
      print *, Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,ichannel
      print *, xchannel,xRnd
!       pause
      Jac = 0d0
   endif

!    print *, "iChannel",ichannel
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,3).dot.Mom(1:4,3) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,4).dot.Mom(1:4,4) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,7).dot.Mom(1:4,7) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,8).dot.Mom(1:4,8) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,9).dot.Mom(1:4,9) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,10).dot.Mom(1:4,10) ))
!    print *, "----------"
!    print *, "Mom.cons. ",Mom(1:4,1)+Mom(1:4,2)-Mom(1:4,5)-Mom(1:4,6)-Mom(1:4,7)-Mom(1:4,8)
!    print *, "Mom.cons. ",Mom_Dummy(1:4)-Mom(1:4,3)-Mom(1:4,4)
!    print *, "Mom.cons. ",Mom(1:4,3)-Mom(1:4,5)-Mom(1:4,6)
!    print *, "Mom.cons. ",Mom(1:4,4)-Mom(1:4,7)-Mom(1:4,8)
!    print *, "----------"
!    print *, "Inv.mass  ",get_MInv(Mom_Dummy(1:4))*100d0
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,3)+Mom(1:4,4))*100d0
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,3))*100d0
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,4))*100d0
!    pause
   
   
   
RETURN
END SUBROUTINE




! H-->gaga or H-->Zga phase space
SUBROUTINE EvalPhasespace_HVga(xRnd,Energy,Mom,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xRnd(:), Energy, s56,s78, Mom(:,:)
real(8) :: Jac,Jac1,Jac2,Jac3,Mom_Dummy(1:4)


   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)
      
   Mom_Dummy(1:4) = (/Energy,0d0,0d0,0d0/)
   if( .not.IsAPhoton(DecayMode1) ) then
      s56=Energy**2   
      Jac1 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(1),s78)
      Jac2 = s_channel_decay(Mom_Dummy(1:4),s78,0d0,xRnd(2:3),Mom(:,3),Mom(:,7)) 
      Jac3 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(4:5),Mom(:,5),Mom(:,6))
      Jac = Jac1 * Jac2 * Jac3 * PSNorm3
   else
      Jac1 = s_channel_decay(Mom_Dummy(1:4),0d0,0d0,xRnd(1:2),Mom(:,5),Mom(:,7))
      Jac = Jac1 * PSNorm2
   endif
   
   if( isNan(jac) ) then
      print *, "EvalPhasespace_Hgaga NaN"
      print *, Energy,Jac
      print *, xRnd
!       pause
      Jac = 0d0
   endif

!    print *, "OS checker", dsqrt( dabs(Mom(1:4,3).dot.Mom(1:4,3) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,4).dot.Mom(1:4,4) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,5).dot.Mom(1:4,5) ))
!    print *, "OS checker", dsqrt( dabs(Mom(1:4,6).dot.Mom(1:4,6) ))
!    print *, "----------"
!    if( .not.IsAPhoton(DecayMode1) ) then   
!       print *, "Mom.cons. ",Mom(1:4,1)+Mom(1:4,2)-Mom(1:4,4) -Mom(1:4,5)-Mom(1:4,6)
!       print *, "Mom.cons. ",Mom(1:4,3)-Mom(1:4,5)-Mom(1:4,6)
!    else   
!       print *, "Mom.cons. ",Mom(1:4,1)+Mom(1:4,2)-Mom(1:4,3) -Mom(1:4,4)
!    endif
!    print *, "----------"
!    pause
   
   
   
RETURN
END SUBROUTINE



! Breit-Wigner mass^2
function bw_sq(x, m, ga, smax, jacobian)
implicit none
real(8) :: bw_sq
real(8), intent(in) :: m, ga, smax,x
real(8) :: xmin, xmax,xprime
real(8), intent(out) :: jacobian

xmin=-datan(m/ga)/m/ga
xmax=-datan((-smax+m**2)/ga/m)/ga/m
xprime=x*(xmax-xmin)+xmin
bw_sq=m**2+dtan(xprime*ga*m)*ga*m
jacobian=(ga*m)**2 * (1d0+dtan(ga*m*xprime)**2) * (xmax-xmin)

return
end function bw_sq




!LORENTZ.F
!VERSION 20130123
!
!A subroutine that performs a general boost to a four vector
!(vector) based on another four vector (boost). The primed and
!unprimed frames have their axes in parallel to one another.
!Rotation is not performed by this subroutine.

      subroutine LORENTZ(vector, boost)

      implicit none

      double precision vector(4), boost(4) 
      double precision lambdaMtrx(4,4), vector_copy(4)
      double precision beta(2:4), beta_sq, gamma
      integer i,j
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

!      double precision KRONECKER_DELTA
!      external KRONECKER_DELTA

      do i=2,4
        beta(i) = boost(i)/boost(1)
      enddo

      beta_sq = beta(2)**2+beta(3)**2+beta(4)**2

  if(beta_sq.ge.epsilon)then

      gamma = 1d0/dsqrt(1d0-beta_sq)

      lambdaMtrx(1,1) = gamma

      do i=2,4
        lambdaMtrx(1,i) = gamma*beta(i)
        lambdaMtrx(i,1) = lambdaMtrx(1,i)
      enddo

      do i=2,4
      do j=2,4
        lambdaMtrx(i,j) = (gamma-1d0)*beta(i)*beta(j)/beta_sq + KRONECKER_DELTA(i,j)
      enddo
      enddo


!apply boost to vector1
      vector_copy = vector
      vector = 0d0
      do i=1,4
      do j=1,4
        vector(i) = vector(i) + lambdaMtrx(i,j)*vector_copy(j)
      enddo
      enddo
  endif

      return
      END subroutine LORENTZ



! KRONECKER_DELTA.F
!
! KRONECKER_DELTA(i,j)
! A function that returns 1 if i=j, and 0 otherwise.
      double precision function KRONECKER_DELTA(i,j)
      integer i,j
      if(i.eq.j)then
        KRONECKER_DELTA = 1d0
      else
        KRONECKER_DELTA = 0d0
      endif

      return
      end function KRONECKER_DELTA





subroutine SetRunningScales(p,id) ! p in JHU-GeV, id in JHUGen conventions
use ModParameters
use ModMisc
implicit none
real(dp), intent(in) :: p(1:4,4:6) ! No need to run the second index from 3 to 7: pH, pJ1, pJ2
integer, intent(in) :: id(4:7) ! id_JJH/id_JJVV, id_J1, id_J2, id_JJ (if applicable)  
real(8) :: polemass(3:7) ! mJJH, mH, mJ1, mJ2, mJJ (if applicable)
real(8) :: pJJHstar(4),pHstar(4),pJ(4,2),pJJ(4),pJHstar(4)
integer idx,ip

   pHstar(:) = 0d0
   pJJ(:) = 0d0
   polemass(3) = getMass(id(4)) ! Pole mass of the JJH system
   polemass(4) = M_Reso
   do idx=4,6
      if(idx.eq.4) then
         do ip=1,4
            pHstar(ip) = pHstar(ip) + p(ip,idx)
         enddo
      else
         polemass(idx) = getMass(id(idx))
         do ip=1,4
            pJJ(ip) = pJJ(ip) + p(ip,idx)
         enddo
      endif
   enddo
   polemass(7) = getMass(id(7)) ! Pole mass of the JJ system

   pJJHstar = pJJ + pHstar
   if(polemass(5).lt.polemass(6)) then
      pJ(:,1)=p(:,5)
      pJ(:,2)=p(:,6)
   else
      pJ(:,1)=p(:,6)
      pJ(:,2)=p(:,5)
      call swapr(polemass(5),polemass(6)) ! will use polemass(5) as the greater mass below
   endif
   pJHstar(1:4) = pJ(1:4,1) + pHstar(1:4)

   ! Determine the appropriate factorization scale for the chosen scheme from pole and invariant masses
   if(FacScheme .eq. kRenFacScheme_mhstar) then
      Mu_Fact = Get_MInv(pHstar(1:4))

   elseif(FacScheme .eq. -kRenFacScheme_mhstar) then
      Mu_Fact = polemass(4)

   elseif(FacScheme .eq. kRenFacScheme_mjjhstar) then
      Mu_Fact = Get_MInv(pJJHstar(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mjjhstar) then
      Mu_Fact = polemass(3)

   elseif(FacScheme .eq. kRenFacScheme_mjj_mhstar) then
      Mu_Fact = Get_MInv(pJJ(1:4))+Get_MInv(pHstar(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mjj_mhstar) then
      Mu_Fact = polemass(4)+polemass(7)
   elseif(FacScheme .eq. kRenFacScheme_mj_mj_mhstar) then
      Mu_Fact = Get_MInv(pJ(1:4,1))+Get_MInv(pJ(1:4,2))+Get_MInv(pHstar(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mj_mj_mhstar) then
      Mu_Fact = polemass(4)+polemass(5)+polemass(6)

   elseif(FacScheme .eq. kRenFacScheme_mjj) then
      Mu_Fact = Get_MInv(pJJ(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mjj) then
      Mu_Fact = polemass(7)
   elseif(FacScheme .eq. kRenFacScheme_mj_mj) then
      Mu_Fact = Get_MInv(pJJ(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mj_mj) then
      Mu_Fact = polemass(5)+polemass(6)

   elseif(FacScheme .eq. kRenFacScheme_mjhstar) then
      Mu_Fact = Get_MInv(pJHstar(1:4))
   elseif(FacScheme .eq. kRenFacScheme_mj_mhstar) then
      Mu_Fact = Get_MInv(pJ(1:4,1))+Get_MInv(pHstar(1:4))
   elseif((FacScheme .eq. -kRenFacScheme_mjhstar) .or. (FacScheme .eq. -kRenFacScheme_mj_mhstar)) then
      Mu_Fact = polemass(4)+polemass(5)
   elseif(FacScheme .eq. kRenFacScheme_mj) then
      Mu_Fact = Get_MInv(pJ(1:4,1))
   elseif(FacScheme .eq. -kRenFacScheme_mj) then
      Mu_Fact = polemass(5)
   endif

   ! Do the same for the renormalization scale
   if(RenScheme .eq. kRenFacScheme_mhstar) then
      Mu_Ren = Get_MInv(pHstar(1:4))

   elseif(RenScheme .eq. -kRenFacScheme_mhstar) then
      Mu_Ren = polemass(4)

   elseif(RenScheme .eq. kRenFacScheme_mjjhstar) then
      Mu_Ren = Get_MInv(pJJHstar(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mjjhstar) then
      Mu_Ren = polemass(3)

   elseif(RenScheme .eq. kRenFacScheme_mjj_mhstar) then
      Mu_Ren = Get_MInv(pJJ(1:4))+Get_MInv(pHstar(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mjj_mhstar) then
      Mu_Ren = polemass(4)+polemass(7)
   elseif(RenScheme .eq. kRenFacScheme_mj_mj_mhstar) then
      Mu_Ren = Get_MInv(pJ(1:4,1))+Get_MInv(pJ(1:4,2))+Get_MInv(pHstar(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mj_mj_mhstar) then
      Mu_Ren = polemass(4)+polemass(5)+polemass(6)

   elseif(RenScheme .eq. kRenFacScheme_mjj) then
      Mu_Ren = Get_MInv(pJJ(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mjj) then
      Mu_Ren = polemass(7)
   elseif(RenScheme .eq. kRenFacScheme_mj_mj) then
      Mu_Ren = Get_MInv(pJJ(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mj_mj) then
      Mu_Ren = polemass(5)+polemass(6)

   elseif(RenScheme .eq. kRenFacScheme_mjhstar) then
      Mu_Ren = Get_MInv(pJHstar(1:4))
   elseif(RenScheme .eq. kRenFacScheme_mj_mhstar) then
      Mu_Ren = Get_MInv(pJ(1:4,1))+Get_MInv(pHstar(1:4))
   elseif((RenScheme .eq. -kRenFacScheme_mjhstar) .or. (RenScheme .eq. -kRenFacScheme_mj_mhstar)) then
      Mu_Ren = polemass(4)+polemass(5)
   elseif(RenScheme .eq. kRenFacScheme_mj) then
      Mu_Ren = Get_MInv(pJ(1:4,1))
   elseif(RenScheme .eq. -kRenFacScheme_mj) then
      Mu_Ren = polemass(5)
   endif

   ! Never ever allow the scales to go negative
   Mu_Fact = abs(Mu_Fact) * MuFacMultiplier
   Mu_Ren = abs(Mu_Ren) * MuRenMultiplier


return
end subroutine SetRunningScales




SUBROUTINE setPDFs(x1,x2,pdf)
use ModParameters
implicit none
real(8) :: x1,x2,PDFScale
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),chm(1:2),bot(1:2),glu(1:2),phot(1:2),sbar(1:2),cbar(1:2),bbar(1:2)
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2),NNpdf(1:2,-6:7)

        PDFScale=Mu_Fact*100d0
        pdf(:,:) = 0d0
        
#if useLHAPDF==1
        call evolvePDF(x1,PDFScale,NNpdf(1,-6:7))
        call evolvePDF(x2,PDFScale,NNpdf(2,-6:7))
            NNpdf(1,-6:7) = NNpdf(1,-6:7)/x1
            NNpdf(2,-6:7) = NNpdf(2,-6:7)/x2
            
            pdf(Up_,1)   = NNpdf(1,+2)         * swPDF_u
            pdf(AUp_,1)  = NNpdf(1,-2)         * swPDF_u
            pdf(Dn_,1)   = NNpdf(1,+1)         * swPDF_d
            pdf(ADn_,1)  = NNpdf(1,-1)         * swPDF_d
            pdf(Chm_,1)  = NNpdf(1,+4)         * swPDF_c
            pdf(AChm_,1) = NNpdf(1,-4)         * swPDF_c
            pdf(Str_,1)  = NNpdf(1,+3)         * swPDF_s
            pdf(AStr_,1) = NNpdf(1,-3)         * swPDF_s
            pdf(Bot_,1)  = NNpdf(1,+5)         * swPDF_b
            pdf(ABot_,1) = NNpdf(1,-5)         * swPDF_b
            pdf(0,1)     = NNpdf(1,+0)         * swPDF_g            
            
            pdf(Up_,2)   = NNpdf(2,+2)         * swPDF_u
            pdf(AUp_,2)  = NNpdf(2,-2)         * swPDF_u
            pdf(Dn_,2)   = NNpdf(2,+1)         * swPDF_d
            pdf(ADn_,2)  = NNpdf(2,-1)         * swPDF_d
            pdf(Chm_,2)  = NNpdf(2,+4)         * swPDF_c
            pdf(AChm_,2) = NNpdf(2,-4)         * swPDF_c
            pdf(Str_,2)  = NNpdf(2,+3)         * swPDF_s
            pdf(AStr_,2) = NNpdf(2,-3)         * swPDF_s
            pdf(Bot_,2)  = NNpdf(2,+5)         * swPDF_b
            pdf(ABot_,2) = NNpdf(2,-5)         * swPDF_b
            pdf(0,2)     = NNpdf(2,+0)         * swPDF_g
            
            pdf(:,:) = dabs(pdf(:,:))
            
            RETURN
            
#else
        if( PDFSet.eq.1 ) then
            call cteq6(x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
            call cteq6(x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
        elseif( PDFSet.eq.2 ) then
            call GetAllPDFs("pdfs/mstw2008lo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot(1))
            str(1)= (str(1)+sbar(1))/2d0
            chm(1)= (chm(1)+cbar(1))/2d0
            bot(1)= (bot(1)+bbar(1))/2d0
            upv(1)=upv(1)/x1
            dnv(1)=dnv(1)/x1
            usea(1)=usea(1)/x1
            dsea(1)=dsea(1)/x1
            str(1)=str(1)/x1
            chm(1)=chm(1)/x1
            bot(1)=bot(1)/x1
            glu(1)=glu(1)/x1
            phot(1)=phot(1)/x1

            call GetAllPDFs("pdfs/mstw2008lo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot(2))
            str(2)= (str(2)+sbar(2))/2d0
            chm(2)= (chm(2)+cbar(2))/2d0
            bot(2)= (bot(2)+bbar(2))/2d0
            upv(2)=upv(2)/x2
            dnv(2)=dnv(2)/x2
            usea(2)=usea(2)/x2
            dsea(2)=dsea(2)/x2
            str(2)=str(2)/x2
            chm(2)=chm(2)/x2
            bot(2)=bot(2)/x2
            glu(2)=glu(2)/x2
            phot(2)=phot(2)/x2

        elseif( PDFSet.ge.201 .and. PDFSet.le.240) then
            call GetAllPDFs("pdfs/mstw2008lo.90cl",PDFSet-200,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot(1))
            str(1)= (str(1)+sbar(1))/2d0
            chm(1)= (chm(1)+cbar(1))/2d0
            bot(1)= (bot(1)+bbar(1))/2d0
            upv(1)=upv(1)/x1
            dnv(1)=dnv(1)/x1
            usea(1)=usea(1)/x1
            dsea(1)=dsea(1)/x1
            str(1)=str(1)/x1
            chm(1)=chm(1)/x1
            bot(1)=bot(1)/x1
            glu(1)=glu(1)/x1
            phot(1)=phot(1)/x1

            call GetAllPDFs("pdfs/mstw2008lo.90cl",PDFSet-200,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot(2))
            str(2)= (str(2)+sbar(2))/2d0
            chm(2)= (chm(2)+cbar(2))/2d0
            bot(2)= (bot(2)+bbar(2))/2d0
            upv(2)=upv(2)/x2
            dnv(2)=dnv(2)/x2
            usea(2)=usea(2)/x2
            dsea(2)=dsea(2)/x2
            str(2)=str(2)/x2
            chm(2)=chm(2)/x2
            bot(2)=bot(2)/x2
            glu(2)=glu(2)/x2
            phot(2)=phot(2)/x2
        elseif( PDFSet.eq.3 ) then

            call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
            call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))
            NNpdf(1,-6:7) = NNpdf(1,-6:7)/x1
            NNpdf(2,-6:7) = NNpdf(2,-6:7)/x2
            
    !       PROTON CONTENT
            pdf(Up_,1)   = NNpdf(1,+2)         * swPDF_u
            pdf(AUp_,1)  = NNpdf(1,-2)         * swPDF_u
            pdf(Dn_,1)   = NNpdf(1,+1)         * swPDF_d
            pdf(ADn_,1)  = NNpdf(1,-1)         * swPDF_d
            pdf(Chm_,1)  = NNpdf(1,+4)         * swPDF_c
            pdf(AChm_,1) = NNpdf(1,-4)         * swPDF_c
            pdf(Str_,1)  = NNpdf(1,+3)         * swPDF_s
            pdf(AStr_,1) = NNpdf(1,-3)         * swPDF_s
            pdf(Bot_,1)  = NNpdf(1,+5)         * swPDF_b
            pdf(ABot_,1) = NNpdf(1,-5)         * swPDF_b
            pdf(0,1)     = NNpdf(1,+0)         * swPDF_g            
            
            pdf(Up_,2)   = NNpdf(2,+2)         * swPDF_u
            pdf(AUp_,2)  = NNpdf(2,-2)         * swPDF_u
            pdf(Dn_,2)   = NNpdf(2,+1)         * swPDF_d
            pdf(ADn_,2)  = NNpdf(2,-1)         * swPDF_d
            pdf(Chm_,2)  = NNpdf(2,+4)         * swPDF_c
            pdf(AChm_,2) = NNpdf(2,-4)         * swPDF_c
            pdf(Str_,2)  = NNpdf(2,+3)         * swPDF_s
            pdf(AStr_,2) = NNpdf(2,-3)         * swPDF_s
            pdf(Bot_,2)  = NNpdf(2,+5)         * swPDF_b
            pdf(ABot_,2) = NNpdf(2,-5)         * swPDF_b
            pdf(0,2)     = NNpdf(2,+0)         * swPDF_g            

            pdf(:,:) = dabs(pdf(:,:))            
            RETURN
        else
            print *, "PDFSet",PDFSet,"not available!"
            stop
        endif
#endif
        
IF( COLLIDER.EQ.1 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u
        pdf(AUp_,1)  = usea(1)             * swPDF_u
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d
        pdf(ADn_,1)  = dsea(1)             * swPDF_d
        pdf(Chm_,1)  = chm(1)              * swPDF_c
        pdf(AChm_,1) = chm(1)              * swPDF_c
        pdf(Str_,1)  = str(1)              * swPDF_s
        pdf(AStr_,1) = str(1)              * swPDF_s
        pdf(Bot_,1)  = bot(1)              * swPDF_b
        pdf(ABot_,1) = bot(1)              * swPDF_b
        pdf(0,1)     = glu(1)              * swPDF_g

!       PROTON CONTENT
        pdf(Up_,2)   = (upv(2) + usea(2))  * swPDF_u
        pdf(AUp_,2)  = usea(2)             * swPDF_u
        pdf(Dn_,2)   = (dnv(2) + dsea(2))  * swPDF_d
        pdf(ADn_,2)  = dsea(2)             * swPDF_d
        pdf(Chm_,2)  = chm(2)              * swPDF_c
        pdf(AChm_,2) = chm(2)              * swPDF_c
        pdf(Str_,2)  = str(2)              * swPDF_s
        pdf(AStr_,2) = str(2)              * swPDF_s
        pdf(Bot_,2)  = bot(2)              * swPDF_b
        pdf(ABot_,2) = bot(2)              * swPDF_b
        pdf(0,2)     = glu(2)              * swPDF_g

ELSEIF( COLLIDER.EQ.2 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u
        pdf(AUp_,1)  = usea(1)             * swPDF_u
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d
        pdf(ADn_,1)  = dsea(1)             * swPDF_d
        pdf(Chm_,1)  = chm(1)              * swPDF_c
        pdf(AChm_,1) = chm(1)              * swPDF_c
        pdf(Str_,1)  = str(1)              * swPDF_s
        pdf(AStr_,1) = str(1)              * swPDF_s
        pdf(Bot_,1)  = bot(1)              * swPDF_b
        pdf(ABot_,1) = bot(1)              * swPDF_b
        pdf(0,1)     = glu(1)              * swPDF_g

!       ANTI-PROTON CONTENT
        pdf(Up_,2)   = usea(2)             * swPDF_u
        pdf(AUp_,2)  = (upv(2)+usea(2))    * swPDF_u
        pdf(Dn_,2)   = dsea(2)             * swPDF_d
        pdf(ADn_,2)  = (dnv(2) + dsea(2))  * swPDF_d
        pdf(Chm_,2)  = chm(2)              * swPDF_c
        pdf(AChm_,2) = chm(2)              * swPDF_c
        pdf(Str_,2)  = str(2)              * swPDF_s
        pdf(AStr_,2) = str(2)              * swPDF_s
        pdf(Bot_,2)  = bot(2)              * swPDF_b
        pdf(ABot_,2) = bot(2)              * swPDF_b
        pdf(0,2)     = glu(2)              * swPDF_g

ENDIF

pdf(:,:) = dabs(pdf(:,:))


RETURN
END SUBROUTINE




SUBROUTINE CTEQ6(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
implicit none
double precision X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU
double precision Q,xsave,qsave,Ctq6Pdf,D,U

         Q=SCALE
         xsave=X
         qsave=Q
         U =         Ctq6Pdf(1,X,Q)
         D =         Ctq6Pdf(2,X,Q)
         USEA =      Ctq6Pdf(-1,X,Q)
         DSEA =      Ctq6Pdf(-2,X,Q)
         STR =       Ctq6Pdf(3,X,Q)
         CHM =       Ctq6Pdf(4,X,Q)
         BOT =       Ctq6Pdf(5,X,Q)
         GLU  =      Ctq6Pdf(0,X,Q)
         UPV=U-USEA
         DNV=D-DSEA
         X=xsave
         Q=qsave
RETURN
END SUBROUTINE





! QCD scale from MCFM
! Implementation into JHUGen by Ulascan Sarica, Dec. 2015
subroutine EvalAlphaS()
   use ModParameters
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
   
      IF (Mu_Ren .LE. 0d0) THEN 
         WRITE(6,*) 'ModKinematics::EvalAlphaS: Mu_Ren .le. 0, Mu_Ren (GeV) = ',(Mu_Ren*GeV)
         stop
      ENDIF
      IF (nQflavors_pdf .NE. NF5) THEN 
         WRITE(6,*) 'ModKinematics::EvalAlphaS: nQflavors_pdf invalid, nQflavors_pdf = ',nQflavors_pdf
         WRITE(6,*) 'ModKinematics::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
         stop
      ENDIF
      IF (nloops_pdf .NE. 1) THEN 
         WRITE(6,*) 'ModKinematics::EvalAlphaS: nloops_pdf invalid, nloops_pdf = ',nloops_pdf
         WRITE(6,*) 'ModKinematics::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
         stop
      ENDIF

      alphas=alphas_mz/(1_dp+alphas_mz*B0_PDF(NF5)*2_dp*dlog((Mu_Ren/zmass_pdf)))
#endif
      ! Calculate the derived couplings
      call ComputeQCDVariables()
   RETURN
end subroutine EvalAlphaS




END MODULE
