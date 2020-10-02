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

!public :: massfrun

contains




SUBROUTINE WriteOutEvent(Mom,MY_IDUP,ICOLUP,MomFSPartons,EventWeight,EventInfoLine,PDFLine,MOTHUP_Parton)
use ModParameters
use modMisc
implicit none
real(8) :: Mom(1:4,1:6)
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
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,mZ1,mZ2,HiggsDKLength,VprimeDKLength(1:2)
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"



! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA



do i=1,9
    LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
enddo

SCALUP=Mu_Fact/GeV
AQEDUP=alpha_QED
AQCDUP=alphas

ISTUP(1:2) = -1 ! Mother status
ISTUP(3:5) = 2 ! Intermediate particle status
ISTUP(6:9) = 1 ! Final hard process particle status
if ( LHE_IDUP(4).eq.22 .and. LHE_IDUP(5).eq.22 ) then! photon+photon FS
    NUP=5
    ISTUP(4:5) = 1
elseif ( LHE_IDUP(4).ne.22 .and. LHE_IDUP(5).eq.22 ) then! Z+photon FS
    NUP=7
    ISTUP(4) = 2
    ISTUP(5:7) = 1
else! ZZ or WW FS
    NUP=9
endif



IDPRUP=Process
if( present(EventWeight) ) then
    XWGTUP=EventWeight
else
    XWGTUP=1.0d0
endif



MOTHUP(1:2,1) = 0
MOTHUP(1:2,2) = 0
MOTHUP(1,3) = 1;MOTHUP(2,3) = 2
MOTHUP(1:2,4) = 3
MOTHUP(1:2,5) = 3
MOTHUP(1:2,6) = 4
MOTHUP(1:2,7) = 4
MOTHUP(1:2,8) = 5
MOTHUP(1:2,9) = 5




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
       MomDummy(1,6+a) = MomFSPartons(1,a)/GeV
       MomDummy(2,6+a) = MomFSPartons(2,a)/GeV
       MomDummy(3,6+a) = MomFSPartons(3,a)/GeV
       MomDummy(4,6+a) = MomFSPartons(4,a)/GeV
    enddo
endif


LHE_IDUP(3) = 25
if( Process.eq.1 ) LHE_IDUP(3) = 32
if( Process.eq.2 ) LHE_IDUP(3) = 39
VprimeDKLength(:) = 0d0
Lifetime = 0.0d0
Spin = 0.1d0
call getHiggsDecayLength(HiggsDKLength)


do a=1,6
    MomDummy(1,a) = Mom(1,a)/GeV
    MomDummy(2,a) = Mom(2,a)/GeV
    MomDummy(3,a) = Mom(3,a)/GeV
    MomDummy(4,a) = Mom(4,a)/GeV
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



!  associate lepton pairs to MOTHUP
if( (IsAZDecay(DecayMode1)).and.(IsAZDecay(DecayMode2)) .and. abs(LHE_IDUP(7)).eq.abs(LHE_IDUP(9)) ) then
     s34 = Get_MInv( Mom(1:4,3)+Mom(1:4,4) )
     s56 = Get_MInv( Mom(1:4,5)+Mom(1:4,6) )
     s36 = Get_MInv( Mom(1:4,3)+Mom(1:4,6) )
     s45 = Get_MInv( Mom(1:4,4)+Mom(1:4,5) )
     smallestInv = minloc((/dabs(s34-M_V),dabs(s56-M_V),dabs(s36-M_V),dabs(s45-M_V)/),1)
     if (includeVprime .and. M_Vprime .ge. 0d0 .and. ( &
         VprimeDecayLengthMassCutoffFactor.le.0d0 .or. (&
            abs(s36-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime .and. abs(s45-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime &
            ) .or. (&
            abs(s34-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime .and. abs(s56-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime &
            ) &
         ) &
      ) then
        smallestInv = minloc((/dabs(s34-M_Vprime),dabs(s56-M_Vprime),dabs(s36-M_Vprime),dabs(s45-M_Vprime)/),1)
     endif
     if( smallestInv.eq.3 .or. smallestInv.eq.4 ) then
        call swapi(MOTHUP(1,7),MOTHUP(1,9))
        call swapi(MOTHUP(2,7),MOTHUP(2,9))
        call swapi(ICOLUP(1,7),ICOLUP(1,9))
        call swapi(ICOLUP(2,7),ICOLUP(2,9))
        Z1FV(1:4) = MomDummy(1:4,3)+MomDummy(1:4,6)
        Z2FV(1:4) = MomDummy(1:4,5)+MomDummy(1:4,4)
     endif
endif





! calculating and checking masses
tmp = Get_MInv2(MomDummy(1:4,1))
if( tmp.lt. -1d-3 ) print *, "Error 1: large negative mass!",tmp
Part1Mass = dSQRT(dabs(tmp))

tmp = Get_MInv2(MomDummy(1:4,2))
if( tmp.lt. -1d-3 ) print *, "Error 2: large negative mass!",tmp
Part2Mass = dSQRT(dabs(tmp))

tmp = Get_MInv2(XFV(1:4))
if( tmp.lt. -1d-3 ) print *, "Error 3: large negative mass!",tmp
XMass = dSQRT(dabs(tmp))

tmp = Get_MInv2(Z1FV(1:4))
if( tmp.lt. -1d-3 ) print *, "Error 4: large negative mass!",tmp
V1Mass = dSQRT(dabs(tmp))
if( V1Mass.lt.1d-5 ) V1Mass=0d0

tmp = Get_MInv2(Z2FV(1:4))
if( tmp.lt. -1d-3 ) print *, "Error 5: large negative mass!",tmp
V2Mass = dSQRT(dabs(tmp))
if( V2Mass.lt.1d-5 ) V2Mass=0d0

tmp = Get_MInv2(MomDummy(1:4,3))
if( tmp.lt. -1d-3 ) print *, "Error 6: large negative mass!",tmp
L11Mass = dSQRT(dABS(tmp))
if( L11Mass.lt.1d-6 ) L11Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,3) = MomDummy(1,3) + 1d-7

tmp = Get_MInv2(MomDummy(1:4,4))
if( tmp.lt. -1d-3 ) print *, "Error 7: large negative mass!",tmp
L12Mass = dSQRT(dABS(tmp))
if( L12Mass.lt.1d-6 ) L12Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,4) = MomDummy(1,4) + 1d-7

tmp = Get_MInv2(MomDummy(1:4,5))
if( tmp.lt. -1d-3 ) print *, "Error 8: large negative mass!",tmp
L21Mass = dSQRT(dABS(tmp))
if( L21Mass.lt.1d-6 ) L21Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,5) = MomDummy(1,5) + 1d-7

tmp = Get_MInv2(MomDummy(1:4,6))
if( tmp.lt. -1d-3 ) print *, "Error 9: large negative mass!",tmp
L22Mass = dSQRT(dABS(tmp))
if( L22Mass.lt.1d-6 ) L22Mass=0d0
if( tmp.lt.0d0 ) MomDummy(1,6) = MomDummy(1,6) + 1d-7

do a=1,NumFSPartons
    tmp = Get_MInv2(MomDummy(1:4,6+a))
    if( tmp.lt. -1d-3 ) print *, "Error: large negative mass",tmp,"for particle",6+a
    PartonMass(a) = dSQRT(dabs(tmp))
enddo






write(io_LHEOutFile,"(A)") "<event>"
if( ReadLHEFile .and. present(EventInfoLine) ) then
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
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,1),MomDummy(1,1),Part1Mass,Lifetime/ctauUnit,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,2),MomDummy(1,2),Part2Mass,Lifetime/ctauUnit,Spin

! X
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),XFV(2:4),XFV(1),XMass,HiggsDKLength/ctauUnit,Spin

! V1
i=4
if (includeVprime .and. M_Vprime .ge. 0d0 .and. (VprimeDecayLengthMassCutoffFactor.le.0d0 .or. abs(V1Mass*GeV-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime) .and. LHE_IDUP(6) .ne. Not_a_particle_ .and. LHE_IDUP(7) .ne. Not_a_particle_) then
  call getVprimeDecayLength(VprimeDKLength(1))
endif
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Z1FV(2:4),Z1FV(1),V1Mass,VprimeDKLength(1)/ctauUnit,Spin

! V2
i=5
if (includeVprime .and. M_Vprime .ge. 0d0 .and. (VprimeDecayLengthMassCutoffFactor.le.0d0 .or. abs(V2Mass*GeV-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime) .and. LHE_IDUP(8) .ne. Not_a_particle_ .and. LHE_IDUP(9) .ne. Not_a_particle_) then
   call getVprimeDecayLength(VprimeDKLength(2))
endif
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Z2FV(2:4),Z2FV(1),V2Mass,VprimeDKLength(2)/ctauUnit,Spin

! decay product 1 (V1): l-, nu or q
i=7
if (LHE_IDUP(i).gt.Not_a_particle_) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,4),MomDummy(1,4),L12Mass,Lifetime/ctauUnit,Spin
endif

! decay product 2 (V1): l+, nubar or qbar
i=6
if (LHE_IDUP(i).gt.Not_a_particle_) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,3),MomDummy(1,3),L11Mass,Lifetime/ctauUnit,Spin
endif

! decay product 1 (V2): l-, nu or q
i=9
if (LHE_IDUP(i).gt.Not_a_particle_) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,6),MomDummy(1,6),L22Mass,Lifetime/ctauUnit,Spin
endif

! decay product 2 (V2): l+, nubar or qbar
i=8
if (LHE_IDUP(i).gt.Not_a_particle_) then
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,5),MomDummy(1,5),L21Mass,Lifetime/ctauUnit,Spin
endif

! additional F.S. partons
do a=1,NumFSPartons
    write(io_LHEOutFile,fmt1) LHE_IDUP(9+a),ISTUP(9+a), MOTHUP(1,9+a),MOTHUP(2,9+a), ICOLUP(1,9+a),ICOLUP(2,9+a),MomDummy(2:4,6+a),MomDummy(1,6+a),PartonMass(a),Lifetime/ctauUnit,Spin
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
real(8) :: Spin, s34,s56,s36,s45,smallestInv
integer :: IDUP(:),ISTUP(:),MOTHUP(:,:),ICOLUP(:,:)
integer :: HiggsDK_IDUP(:),HiggsDK_ICOLUP(:,:),HiggsDK_ISTUP(4:9),HiggsDK_MOTHUP(1:2,4:9)
integer :: i,iHiggs
integer :: NUP,NUP_NEW,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime(3:9)
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

    if( present(EventProcessId) ) then
        IDPRUP=EventProcessId
    else
        IDPRUP=Process
    endif
    if( present(EventWeight) ) then
        XWGTUP=EventWeight
    else
        XWGTUP=1.0d0
    endif
    if( present(EventScaleAqedAqcd) ) then
        SCALUP=EventScaleAqedAqcd(1)
        AQEDUP=EventScaleAqedAqcd(2)
        AQCDUP=EventScaleAqedAqcd(3)
    else
        SCALUP=Mu_Fact/GeV
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

    Lifetime(:) = 0d0
    Spin = 0.1d0
    call getHiggsDecayLength(Lifetime(3))

    !  associate lepton pairs to MOTHUP
    if( (IsAZDecay(DecayMode1)).and.(IsAZDecay(DecayMode2)) .and. abs(HiggsDK_IDUP(7)).eq.abs(HiggsDK_IDUP(9)) ) then
        s34 = Get_MInv( HiggsDK_Mom(1:4,3)+HiggsDK_Mom(1:4,4) )
        s56 = Get_MInv( HiggsDK_Mom(1:4,5)+HiggsDK_Mom(1:4,6) )
        s36 = Get_MInv( HiggsDK_Mom(1:4,3)+HiggsDK_Mom(1:4,6) )
        s45 = Get_MInv( HiggsDK_Mom(1:4,4)+HiggsDK_Mom(1:4,5) )
        smallestInv = minloc((/dabs(s34-M_V),dabs(s56-M_V),dabs(s36-M_V),dabs(s45-M_V)/),1)
        if (includeVprime .and. M_Vprime .ge. 0d0 .and. ( &
            VprimeDecayLengthMassCutoffFactor.le.0d0 .or. (&
               abs(s36-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime .and. abs(s45-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime &
               ) .or. (&
               abs(s34-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime .and. abs(s56-M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime &
               ) &
            ) &
         ) then
           smallestInv = minloc((/dabs(s34-M_Vprime),dabs(s56-M_Vprime),dabs(s36-M_Vprime),dabs(s45-M_Vprime)/),1)
        endif
        if( smallestInv.eq.3 .or. smallestInv.eq.4 ) then
            call swapi(HiggsDK_MOTHUP(1,7),HiggsDK_MOTHUP(1,9))
            call swapi(HiggsDK_MOTHUP(2,7),HiggsDK_MOTHUP(2,9))
            call swapi(HiggsDK_ICOLUP(1,7),HiggsDK_ICOLUP(1,9))
            call swapi(HiggsDK_ICOLUP(2,7),HiggsDK_ICOLUP(2,9))
            HiggsDK_Mom(1:4,1) = HiggsDK_Mom(1:4,3)+ HiggsDK_Mom(1:4,6)
            HiggsDK_Mom(1:4,2) = HiggsDK_Mom(1:4,4)+ HiggsDK_Mom(1:4,5)
        endif
    endif

    if (includeVprime .and. M_Vprime .ge. 0d0 .and. (VprimeDecayLengthMassCutoffFactor.le.0d0 .or. abs(Get_MInv(HiggsDK_Mom(1:4,1)) - M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime) .and. .not.IsAPhoton(DecayMode1)) then
      call getVprimeDecayLength(Lifetime(4))
    endif
    if (includeVprime .and. M_Vprime .ge. 0d0 .and. (VprimeDecayLengthMassCutoffFactor.le.0d0 .or. abs(Get_MInv(HiggsDK_Mom(1:4,2)) - M_Vprime).lt.VprimeDecayLengthMassCutoffFactor*Ga_Vprime) .and. .not.IsAPhoton(DecayMode2)) then
      call getVprimeDecayLength(Lifetime(5))
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
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,Lifetime(3)/ctauUnit, Spin
            else
               write(io_LHEOutFile,IndentedFmt1) IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),  &
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,0d0, Spin
            endif
        enddo


!       write new intermediate particles and Higgs decay products
        !call swap_mom(HiggsDK_Mom(1:4,3),HiggsDK_Mom(1:4,4))! swap to account for flipped asignments
        !call swap_mom(HiggsDK_Mom(1:4,5),HiggsDK_Mom(1:4,6))! swap to account for flipped asignments
        do i = 4,4 + (NUP_NEW-1)
            write(io_LHEOutFile,IndentedFmt1) HiggsDK_IDUP(i),HiggsDK_ISTUP(i), HiggsDK_MOTHUP(1,i),HiggsDK_MOTHUP(2,i), HiggsDK_ICOLUP(1,i),HiggsDK_ICOLUP(2,i),  &
                                              HiggsDK_Mom(2:4,i-3)/GeV,HiggsDK_Mom(1,i-3)/GeV, get_MInv(HiggsDK_Mom(1:4,i-3))/GeV, Lifetime(i)/ctauUnit, Spin
        enddo
    endif

RETURN
END SUBROUTINE


SUBROUTINE WriteOutEvent_gg4f_fullproddec(Mom,MY_IDUP,ICOLUP,EventWeight)
use ModParameters
use ModMisc
implicit none
integer,parameter :: inTop=1, inBot=2, V1=3, V2=4, Lep1P=5, Lep1M=6, Lep2P=7, Lep2M=8, NUP=8
real(8) :: Mom(1:4,1:NUP),xRnd,s34,s36,s45,s56
real(8),optional :: EventWeight
integer :: MY_IDUP(1:NUP),ICOLUP(1:2,1:NUP),LHE_IDUP(1:NUP),ISTUP(1:NUP),MOTHUP(1:2,1:NUP)
integer :: IDPRUP,i,smallestInv,j
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,MomDummy(1:4,1:NUP),TheMass,mjj
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
integer, parameter :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)
! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


      IDPRUP=Process
      SCALUP=Mu_Fact/GeV
      AQEDUP=alpha_QED
      AQCDUP=alphas

      ISTUP(inTop:inBot) = -1
      ISTUP(V1:V2) = 2
      ISTUP(V2+1:NUP) = +1

      ICOLUP(1:2,V1:V2) = 0
      if( IsAGluon(MY_IDUP(1)) .and. IsAGluon(MY_IDUP(2)) ) then! gg->VV
         ICOLUP(1:2,1) = (/501,502/)
         ICOLUP(1:2,2) = (/502,501/)
      elseif (MY_IDUP(1).gt.0 .and. MY_IDUP(2).lt.0) then  ! qqb->VV
         ICOLUP(1:2,1) = (/501,000/)
         ICOLUP(1:2,2) = (/000,501/)
      elseif (MY_IDUP(2).gt.0 .and. MY_IDUP(1).lt.0) then  ! qbq->VV
         ICOLUP(1:2,2) = (/501,000/)
         ICOLUP(1:2,1) = (/000,501/)
      else
        print *, "Colors for the gg4f process cannot be resolved."
        print *, MY_IDUP(inTop:inBot)
        stop 1
      endif


      MOTHUP(1:2,inTop:inBot) = 0
      MOTHUP(1:2,V1)=(/inTop,inBot/)
      MOTHUP(1:2,V2)=(/inTop,inBot/)
      MOTHUP(1:2,Lep1P)= (/V1,V1/)
      MOTHUP(1:2,Lep1M)= (/V1,V1/)
      MOTHUP(1:2,Lep2P)= (/V2,V2/)
      MOTHUP(1:2,Lep2M)= (/V2,V2/)


      if( present(EventWeight) ) then
          XWGTUP=EventWeight
      else
          XWGTUP=1d0
      endif
      Lifetime = 0.0d0
      Spin     = 0.1d0


      do i=1,NUP
          LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
          MomDummy(1,i) = Mom(1,i)/GeV
          MomDummy(2,i) = Mom(2,i)/GeV
          MomDummy(3,i) = Mom(3,i)/GeV
          MomDummy(4,i) = Mom(4,i)/GeV
      enddo


!  associate lepton pairs to MOTHUP
      if( MY_IDUP(Lep1P).eq.MY_IDUP(Lep2P) ) then
          s34 = Get_MInv( Mom(1:4,Lep1P)+Mom(1:4,Lep1M) )
          s56 = Get_MInv( Mom(1:4,Lep2P)+Mom(1:4,Lep2M) )
          s36 = Get_MInv( Mom(1:4,Lep1P)+Mom(1:4,Lep2M) )
          s45 = Get_MInv( Mom(1:4,Lep2P)+Mom(1:4,Lep1M) )
          smallestInv = minloc((/dabs(s34-M_V),dabs(s56-M_V),dabs(s36-M_V),dabs(s45-M_V)/),1)
          if( smallestInv.eq.3 .or. smallestInv.eq.4 ) then
              call swap(MOTHUP(1,Lep1P),MOTHUP(1,Lep2P))
              call swap(MOTHUP(2,Lep1P),MOTHUP(2,Lep2P))
              call swap(ICOLUP(1,Lep1P),ICOLUP(1,Lep2P))
              call swap(ICOLUP(2,Lep1P),ICOLUP(2,Lep2P))
              call swap(MomDummy(1:4,Lep1P),MomDummy(1:4,Lep2P))
          endif
      endif
      MomDummy(1:4,V1) = MomDummy(1:4,Lep1P)+MomDummy(1:4,Lep1M)
      MomDummy(1:4,V2) = MomDummy(1:4,Lep2P)+MomDummy(1:4,Lep2M)


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
     if( i.le.inBot  ) TheMass = 0d0  ! setting incoming quark/gluon masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime/ctauUnit,Spin
enddo

write(io_LHEOutFile,"(A)") "</event>"


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


    if( present(EventProcessId) ) then
        IDPRUP=EventProcessId
    else
        IDPRUP=Process
    endif
    if( present(EventWeight) ) then
        XWGTUP=EventWeight
    else
        XWGTUP=1.0d0
    endif
    if( present(EventScaleAqedAqcd) ) then
        SCALUP=EventScaleAqedAqcd(1)
        AQEDUP=EventScaleAqedAqcd(2)
        AQCDUP=EventScaleAqedAqcd(3)
    else
        SCALUP=Mu_Fact/GeV
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
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,HiggsDKLength/ctauUnit, Spin
            else
               write(io_LHEOutFile,IndentedFmt1) IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),  &
                                                 Mom(2:4,i)/GeV,Mom(1,i)/GeV, Mass(i)/GeV,Lifetime/ctauUnit, Spin
            endif
        enddo

!       write new intermediate particles and Higgs decay products
        do i = 4,4 + (NUP_NEW-1)
            write(io_LHEOutFile,IndentedFmt1) HiggsDK_IDUP(i),HiggsDK_ISTUP(i), HiggsDK_MOTHUP(1,i),HiggsDK_MOTHUP(2,i), HiggsDK_ICOLUP(1,i),HiggsDK_ICOLUP(2,i),  &
                                              HiggsDK_Mom(2:4,i)/GeV,HiggsDK_Mom(1,i)/GeV, get_MInv(HiggsDK_Mom(1:4,i))/GeV, Lifetime/ctauUnit, Spin
        enddo
    endif

RETURN
END SUBROUTINE


SUBROUTINE ShiftMass(p1,p2,m1,m2,p1hat,p2hat,MassWeight)
use ModMisc
implicit none
real(8),intent(in) :: p1(1:4),p2(1:4)
real(8) :: m1,m2,p1hat(1:4),p2hat(1:4)
real(8),optional :: MassWeight
real(8) :: xi,eta,a,b,c,p1sq,p2sq,p1p2
real(8) :: p1hatsq, p2hatsq, p12hatsq

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

  if( present(MassWeight) ) then
     p1hatsq = p1hat(1:4).dot.p1hat(1:4)
     p2hatsq = p2hat(1:4).dot.p2hat(1:4)
     p12hatsq = (p1hat(1:4)+p2hat(1:4)).dot.(p1hat(1:4)+p2hat(1:4)) ! Should be p1p2*2d0, but better to avoid un-anticipated uses

     ! Below is the same as 2d0/pi/g_d(p12hatsq,p1hatsq,p2hatsq)
     MassWeight = (p12hatsq**2+p1hatsq**2+p2hatsq**2-2d0*(p1hatsq*p2hatsq+p1hatsq*p12hatsq+p12hatsq*p2hatsq)) ! Writing this way instead of get_MInv should avoid the issue of - vs + invariant masses
     if(MassWeight.ge.0d0 .and. p12hatsq.ne.0d0) then
        MassWeight = sqrt(MassWeight/(p12hatsq**2))
     else
        MassWeight = 0d0
     endif
  endif

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
SCALUP=Mu_Fact/GeV
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
    MomDummy(1,a) = Mom(1,a)/GeV
    MomDummy(2,a) = Mom(2,a)/GeV
    MomDummy(3,a) = Mom(3,a)/GeV
    MomDummy(4,a) = Mom(4,a)/GeV
enddo

! do b=1,4
!     Z1FV(b) = MomDummy(b,1)-MomDummy(b,3)
!     Z2FV(b) = MomDummy(b,2)+MomDummy(b,4)
! enddo


Part1Mass = 0d0
Part2Mass = 0d0
Part4Mass = 0d0
XMass = M_Reso/GeV



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
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part1Mass,Lifetime/ctauUnit,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part2Mass,Lifetime/ctauUnit,Spin

! H
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),XMass,HiggsDKLength/ctauUnit,Spin

! j1
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part4Mass,Lifetime/ctauUnit,Spin


write(io_LHEOutFile,"(A)") "</event>"



END SUBROUTINE


SUBROUTINE WriteOutEvent_HJJ_fulldecay(Mom,MY_IDUP,ICOLUP,do78,EventWeight)
use ModParameters
use ModMisc
implicit none
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10, NUP=10
real(8) :: Mom(1:4,1:NUP),xRnd,s34,s36,s45,s56
real(8),optional :: EventWeight
integer :: MY_IDUP(1:NUP),ICOLUP(1:2,1:NUP),LHE_IDUP(1:NUP),ISTUP(1:NUP),MOTHUP(1:2,1:NUP)
integer :: IDPRUP,i,smallestInv,j
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,MomDummy(1:4,1:NUP),TheMass,mjj
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.0)"
integer, parameter :: LHA2M_ID(-6:6)  = (/-5,-6,-3,-4,-1,-2,10,2,1,4,3,6,5/)
logical :: do78, canbeVBF, canbeVH, isVHlike
! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


      IDPRUP=Process
      SCALUP=Mu_Fact/GeV
      AQEDUP=alpha_QED
      AQCDUP=alphas

      ISTUP(1:2) = -1
      MOTHUP(1:2,1:2) = 0

      if (Process .eq. 69) then
         if( IsAGluon(MY_IDUP(1)) .and. IsAGluon(MY_IDUP(2)) ) then! gg->gg/qqb
            ICOLUP(1:2,1) = (/501,502/)
            ICOLUP(1:2,2) = (/503,501/)
            if (IsAGluon(MY_IDUP(3))) then ! gg->gg
               call random_number(xRnd)
               if (xRnd.gt.0.5) then
                  ICOLUP(1:2,3) = (/504,502/)
                  ICOLUP(1:2,4) = (/503,504/)
               else
                  ICOLUP(1:2,4) = (/504,502/)
                  ICOLUP(1:2,3) = (/503,504/)
               endif
            elseif (MY_IDUP(3).gt.0) then  ! gg->qqb
               ICOLUP(1:2,3) = (/503,000/)
               ICOLUP(1:2,4) = (/000,502/)
            else                           ! gg->qbq
               ICOLUP(1:2,4) = (/503,000/)
               ICOLUP(1:2,3) = (/000,502/)
            endif
         elseif (IsAGluon(MY_IDUP(3)) .and. IsAGluon(MY_IDUP(4))) then ! qqb->gg
            ICOLUP(1:2,3) = (/501,502/)
            ICOLUP(1:2,4) = (/503,501/)
            if (MY_IDUP(1).gt.0) then  ! qqb->gg
               ICOLUP(1:2,1) = (/503,000/)
               ICOLUP(1:2,2) = (/000,502/)
            else                           ! qbq->gg
               ICOLUP(1:2,2) = (/503,000/)
               ICOLUP(1:2,1) = (/000,502/)
            endif
         elseif( (IsAQuark(MY_IDUP(1)) .and. IsAGluon(MY_IDUP(2))) .or. (IsAQuark(MY_IDUP(2)) .and. IsAGluon(MY_IDUP(1))) ) then! qg/gq->qg/gq or qbg/gqb->qbg/gqb
            ICOLUP(1:2,1) = (/501,000/)
            ICOLUP(1:2,2) = (/502,501/)
            ICOLUP(1:2,3) = (/503,000/)
            ICOLUP(1:2,4) = (/502,503/)

            if ( IsAGluon(MY_IDUP(1)) ) then
              do j=1,2; call swap(ICOLUP(j,1),ICOLUP(j,2)); enddo
              if (MY_IDUP(2).lt.0) then
                 call swap(ICOLUP(1,1),ICOLUP(2,1))
                 call swap(ICOLUP(1,2),ICOLUP(2,2))
              endif
            else
               if (MY_IDUP(1).lt.0) then
                 call swap(ICOLUP(1,1),ICOLUP(2,1))
                 call swap(ICOLUP(1,2),ICOLUP(2,2))
               endif
            endif

            if ( IsAGluon(MY_IDUP(3)) ) then
              do j=1,2; call swap(ICOLUP(j,3),ICOLUP(j,4)); enddo
              if (MY_IDUP(4).lt.0) then
                 call swap(ICOLUP(1,3),ICOLUP(2,3))
                 call swap(ICOLUP(1,4),ICOLUP(2,4))
              endif
            else
              if (MY_IDUP(3).lt.0) then
                 call swap(ICOLUP(1,3),ICOLUP(2,3))
                 call swap(ICOLUP(1,4),ICOLUP(2,4))
              endif
            endif
         ! qq->qq
         elseif (abs(MY_IDUP(1)).eq.abs(MY_IDUP(2))) then
            ICOLUP(1:2,1) = (/501,000/)
            ICOLUP(1:2,2) = (/501,000/)
            ICOLUP(1:2,3) = (/502,000/)
            ICOLUP(1:2,4) = (/502,000/)
            if (MY_IDUP(1).eq.MY_IDUP(3) .and. MY_IDUP(1).eq.MY_IDUP(4)) then
               call random_number(xRnd)
               if (xRnd.lt.0.5d0) then
                  call swap(ICOLUP(1,2),ICOLUP(1,3))
                  call swap(ICOLUP(2,2),ICOLUP(2,3))
               else
                  call swap(ICOLUP(1,2),ICOLUP(1,4))
                  call swap(ICOLUP(2,2),ICOLUP(2,4))
               endif
            else if (MY_IDUP(1).eq.MY_IDUP(3)) then
               call random_number(xRnd)
               if (xRnd.lt.0.5d0) then
                  call swap(ICOLUP(1,2),ICOLUP(1,3))
                  call swap(ICOLUP(2,2),ICOLUP(2,3))
               ! else leave colors alone
               endif
            else if (MY_IDUP(1).eq.MY_IDUP(4)) then
               call random_number(xRnd)
               if (xRnd.lt.0.5d0) then
                  call swap(ICOLUP(1,2),ICOLUP(1,4))
                  call swap(ICOLUP(2,2),ICOLUP(2,4))
               ! else leave colors alone
               endif
            endif
            if (MY_IDUP(1).lt.0) then
               call swap(ICOLUP(1,1),ICOLUP(2,1))
            endif
            if (MY_IDUP(2).lt.0) then
               call swap(ICOLUP(1,2),ICOLUP(2,2))
            endif
            if (MY_IDUP(3).lt.0) then
               call swap(ICOLUP(1,3),ICOLUP(2,3))
            endif
            if (MY_IDUP(4).lt.0) then
               call swap(ICOLUP(1,4),ICOLUP(2,4))
            endif
         elseif (abs(MY_IDUP(1)).eq.abs(MY_IDUP(3))) then
            ICOLUP(1:2,1) = (/501,000/)
            ICOLUP(1:2,2) = (/502,000/)
            ICOLUP(1:2,3) = (/501,000/)
            ICOLUP(1:2,4) = (/502,000/)
            if (MY_IDUP(1).lt.0) then
               call swap(ICOLUP(1,1),ICOLUP(2,1))
            endif
            if (MY_IDUP(2).lt.0) then
               call swap(ICOLUP(1,2),ICOLUP(2,2))
            endif
            if (MY_IDUP(3).lt.0) then
               call swap(ICOLUP(1,3),ICOLUP(2,3))
            endif
            if (MY_IDUP(4).lt.0) then
               call swap(ICOLUP(1,4),ICOLUP(2,4))
            endif
         elseif (abs(MY_IDUP(1)).eq.abs(MY_IDUP(4))) then
            ICOLUP(1:2,1) = (/501,000/)
            ICOLUP(1:2,2) = (/502,000/)
            ICOLUP(1:2,3) = (/502,000/)
            ICOLUP(1:2,4) = (/501,000/)
            if (MY_IDUP(1).lt.0) then
               call swap(ICOLUP(1,1),ICOLUP(2,1))
            endif
            if (MY_IDUP(2).lt.0) then
               call swap(ICOLUP(1,2),ICOLUP(2,2))
            endif
            if (MY_IDUP(3).lt.0) then
               call swap(ICOLUP(1,3),ICOLUP(2,3))
            endif
            if (MY_IDUP(4).lt.0) then
               call swap(ICOLUP(1,4),ICOLUP(2,4))
            endif
         else
           print *, "Color for this Process 69 configuration cannot be resolved"
           print *, MY_IDUP(1:4)
           stop 1
         endif

      else
         canbeVH = .false.
         canbeVBF = .false.

         !try VBF first because there are more options there, but overwrite with VH if that's better

         if( MY_IDUP(1).gt.0 ) then ! quark
           ICOLUP(1:2,1) = (/501,000/)
         else! anti-quark
           ICOLUP(1:2,1) = (/000,501/)
         endif
         if( MY_IDUP(2).gt.0 ) then! quark
           ICOLUP(1:2,2) = (/502,000/)
         else! anti-quark
           ICOLUP(1:2,2) = (/000,502/)
         endif

         !WW is typically larger xsec, so choose that if possible, then ZZ
         if (.not.IsAQuark(MY_IDUP(3))) then
           ICOLUP(1:2,3) = 0
           ICOLUP(1:2,4) = 0
           canbeVH = .true.
         elseif(       abs(CoupledVertex((/MY_IDUP(1),-MY_IDUP(3)/), -1)) .eq. abs(Wp_) &
             .and. abs(CoupledVertex((/MY_IDUP(2),-MY_IDUP(4)/), -1)) .eq. abs(Wp_)) then
           ICOLUP(1:2,3) = ICOLUP(1:2,1)
           ICOLUP(1:2,4) = ICOLUP(1:2,2)
           canbeVBF = .true.
         elseif(       abs(CoupledVertex((/MY_IDUP(1),-MY_IDUP(4)/), -1)) .eq. abs(Wp_) &
                 .and. abs(CoupledVertex((/MY_IDUP(2),-MY_IDUP(3)/), -1)) .eq. abs(Wp_)) then
           ICOLUP(1:2,3) = ICOLUP(1:2,2)
           ICOLUP(1:2,4) = ICOLUP(1:2,1)
           canbeVBF = .true.
         elseif(       abs(CoupledVertex((/MY_IDUP(1),-MY_IDUP(3)/), -1)) .eq. abs(Z0_) &
                 .and. abs(CoupledVertex((/MY_IDUP(2),-MY_IDUP(4)/), -1)) .eq. abs(Z0_)) then
           ICOLUP(1:2,3) = ICOLUP(1:2,1)
           ICOLUP(1:2,4) = ICOLUP(1:2,2)
           canbeVBF = .true.
         elseif(       abs(CoupledVertex((/MY_IDUP(1),-MY_IDUP(4)/), -1)) .eq. abs(Z0_) &
                 .and. abs(CoupledVertex((/MY_IDUP(2),-MY_IDUP(3)/), -1)) .eq. abs(Z0_)) then
           ICOLUP(1:2,3) = ICOLUP(1:2,2)
           ICOLUP(1:2,4) = ICOLUP(1:2,1)
           canbeVBF = .true.
         endif

         if( CoupledVertex(MY_IDUP(1:2), -1) .eq. CoupledVertex(MY_IDUP(3:4), -1) .and. CoupledVertex(MY_IDUP(3:4), -1) .ne. Not_a_particle_) then
           canbeVH = .true.
         endif

         if( .not. canbeVH .and. .not. canbeVBF) then
           print *, "Event doesn't make sense"
           print *, MY_IDUP(1:4)
           stop 1
         endif

         if( canbeVH .and. canbeVBF ) then
           mjj = Get_MInv( Mom(1:4,3)+Mom(1:4,4) )
           if( abs(CoupledVertex(MY_IDUP(1:2), -1)).eq.abs(Wp_) ) then
             isVHlike = ( mjj .gt. M_W-2d0*Ga_W .and. mjj .lt. M_W+2d0*Ga_W )
           else
             isVHlike = ( mjj .gt. M_Z-2d0*Ga_Z .and. mjj .lt. M_Z+2d0*Ga_Z )
           endif
         elseif( canbeVH ) then
           isVHlike = .true.
         endif

         if( isVHlike ) then
           if( MY_IDUP(1).gt.0 ) then ! quark
             ICOLUP(1:2,1) = (/501,000/)
             ICOLUP(1:2,2) = (/000,501/)
           else! anti-quark
             ICOLUP(1:2,1) = (/000,501/)
             ICOLUP(1:2,2) = (/501,000/)
           endif
           if (IsAQuark(MY_IDUP(3))) then
              if( MY_IDUP(3).gt.0 ) then ! quark
                ICOLUP(1:2,3) = (/502,000/)
                ICOLUP(1:2,4) = (/000,502/)
              else! anti-quark
                ICOLUP(1:2,3) = (/000,502/)
                ICOLUP(1:2,4) = (/502,000/)
              endif
           endif
         endif

         if( MY_IDUP(3).eq.Glu_) then ! 10 = LHA2M_ID(0) for the case where MY_IDUP(3:4) are not set (i.e. deterministicME=.false.). 10 doesn't mean gluon here
            call Error("This shouldn't be able to happen anymore, assert false")
            call random_number(xRnd)
            if( xRnd.lt.0.5d0 ) then! ZZ fusion
                   MY_IDUP(3:4)= (/MY_IDUP(1),MY_IDUP(2)/)
            else! WW fusion
                   MY_IDUP(3) = -GetCKMPartner( MY_IDUP(1) )
                   MY_IDUP(4) = -GetCKMPartner( MY_IDUP(2) )
                   if( abs(MY_IDUP(3)).eq.Top_ ) MY_IDUP(3) = MY_IDUP(1)
                   if( abs(MY_IDUP(4)).eq.Top_ ) MY_IDUP(4) = MY_IDUP(2)
            endif
         endif

      endif

      ISTUP(3:10) = +1
      MOTHUP(1:2,3)= (/1,2/)
      MOTHUP(1:2,4)= (/1,2/)


      ISTUP(5:6) = 2
      ICOLUP(1:2,5:6) = 0
      MOTHUP(1:2,5)=(/1,2/)
      MOTHUP(1:2,6)=(/1,2/)

      MOTHUP(1:2,7 )= (/5,5/)
      MOTHUP(1:2,8 )= (/5,5/)
      MOTHUP(1:2,9 )= (/6,6/)
      MOTHUP(1:2,10)= (/6,6/)


      if( present(EventWeight) ) then
          XWGTUP=EventWeight
      else
          XWGTUP=1.0d0
      endif
      Lifetime = 0.0d0
      Spin     = 0.1d0


      do i=1,NUP
          LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
          MomDummy(1,i) = Mom(1,i)/GeV
          MomDummy(2,i) = Mom(2,i)/GeV
          MomDummy(3,i) = Mom(3,i)/GeV
          MomDummy(4,i) = Mom(4,i)/GeV
      enddo


!  associate lepton pairs to MOTHUP
      if( MY_IDUP(7).eq.MY_IDUP(9) ) then
          s34 = Get_MInv( Mom(1:4,7)+Mom(1:4,8) )
          s56 = Get_MInv( Mom(1:4,9)+Mom(1:4,10) )
          s36 = Get_MInv( Mom(1:4,7)+Mom(1:4,10) )
          s45 = Get_MInv( Mom(1:4,8)+Mom(1:4,9) )
          smallestInv = minloc((/dabs(s34-M_V),dabs(s56-M_V),dabs(s36-M_V),dabs(s45-M_V)/),1)
          if( smallestInv.eq.3 .or. smallestInv.eq.4 ) then
              call swapi(MOTHUP(1,7),MOTHUP(1,9))
              call swapi(MOTHUP(2,7),MOTHUP(2,9))
              call swapi(ICOLUP(1,7),ICOLUP(1,9))
              call swapi(ICOLUP(2,7),ICOLUP(2,9))
              MomDummy(1:4,5) = MomDummy(1:4,9)+MomDummy(1:4,8)
              MomDummy(1:4,6) = MomDummy(1:4,7)+MomDummy(1:4,10)
          else
              !The Z's are not always lined up this way already, for reasons I don't understand.
              MomDummy(1:4,5) = MomDummy(1:4,7)+MomDummy(1:4,8)
              MomDummy(1:4,6) = MomDummy(1:4,9)+MomDummy(1:4,10)
          endif
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
     TheMass = get_Minv(MomDummy(:,i))
     if( i.le.inbot .or. (.not. do78 .and. i.le.outBot)  ) TheMass = 0.0d0  ! setting quark masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime/ctauUnit,Spin
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
SCALUP=Mu_Fact/GeV
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
    MomDummy(1,i) = Mom(1,i)/GeV
    MomDummy(2,i) = Mom(2,i)/GeV
    MomDummy(3,i) = Mom(3,i)/GeV
    MomDummy(4,i) = Mom(4,i)/GeV
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
          MomDummy(1,i) = MomDummy(1,i)/GeV
          MomDummy(2,i) = MomDummy(2,i)/GeV
          MomDummy(3,i) = MomDummy(3,i)/GeV
          MomDummy(4,i) = MomDummy(4,i)/GeV
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
     TheMass = get_Minv(MomDummy(:,i))
     if( i.le.2  ) TheMass = 0d0  ! setting initial parton masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime/ctauUnit,Spin
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
SCALUP=Mu_Fact/GeV
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
    MomDummy(1,i) = Mom(1,i)/GeV
    MomDummy(2,i) = Mom(2,i)/GeV
    MomDummy(3,i) = Mom(3,i)/GeV
    MomDummy(4,i) = Mom(4,i)/GeV
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
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),0d0,Lifetime/ctauUnit,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),0d0,Lifetime/ctauUnit,Spin

! H
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),M_Reso/GeV,Lifetime/ctauUnit,Spin

! bb
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),m_top/GeV,Lifetime/ctauUnit,Spin

! b
i=5
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),m_top/GeV,Lifetime/ctauUnit,Spin



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
SCALUP=Mu_Fact/GeV
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
    MomDummy(1,i) = Mom(1,i)/GeV
    MomDummy(2,i) = Mom(2,i)/GeV
    MomDummy(3,i) = Mom(3,i)/GeV
    MomDummy(4,i) = Mom(4,i)/GeV
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
          MomDummy(1,i) = MomDummy(1,i)/GeV
          MomDummy(2,i) = MomDummy(2,i)/GeV
          MomDummy(3,i) = MomDummy(3,i)/GeV
          MomDummy(4,i) = MomDummy(4,i)/GeV
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
     TheMass = get_Minv(MomDummy(:,i))
     if( i.le.2  ) TheMass = 0.0d0  ! setting initial parton masses to zero
     write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),TheMass,Lifetime/ctauUnit,Spin
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
SCALUP=Mu_Fact/GeV
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
    MomDummy(1,a) = Mom(1,a)/GeV
    MomDummy(2,a) = Mom(2,a)/GeV
    MomDummy(3,a) = Mom(3,a)/GeV
    MomDummy(4,a) = Mom(4,a)/GeV
enddo

! do b=1,4
!     Z1FV(b) = MomDummy(b,1)-MomDummy(b,3)
!     Z2FV(b) = MomDummy(b,2)+MomDummy(b,4)
! enddo


Part1Mass = 0d0
Part2Mass = 0d0
Part3Mass = 0d0
Part4Mass = 0d0
XMass = M_Reso/GeV



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
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part1Mass,Lifetime/ctauUnit,Spin

! parton_b
i=2
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part2Mass,Lifetime/ctauUnit,Spin

! j1
i=3
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part3Mass,Lifetime/ctauUnit,Spin

! j2
i=4
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),Part4Mass,Lifetime/ctauUnit,Spin

! H
i=5
write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,i),MomDummy(1,i),XMass,HiggsDKLength/ctauUnit,Spin


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



MomDummy = MomExt/GeV
MassDummy = inv_mass/GeV

IDPRUP=Process
SCALUP=Mu_Fact/GeV
AQEDUP=alpha_QED
AQCDUP=alphas

call getHiggsDecayLength(HiggsDKLength)

NUP = 4
if(H_DK) NUP=NUP+2
if(.not.IsAPhoton(DecayMode1)) NUP=NUP+2

    XWGTUP=EventWeight

write(io_LHEOutFile,"(A)") "<event>"
if( .not. ReadLHEFile ) write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

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

do i=1,2
    write(io_LHEOutFile,fmt1) id(i), -1,0,0,ICOLUP(i,1),ICOLUP(i,2),MomDummy(2:4,i), MomDummy(1,i), 0.0d0, 0.0d0, Spin
enddo

if(IsAPhoton(DecayMode1)) then
  write(io_LHEOutFile,fmt1) id(4), 1,1,2,0,0,MomDummy(2:4,4), MomDummy(1,4), MassDummy(4), 0d0, Spin
else
  write(io_LHEOutFile,fmt1) id(4), 2,1,2,0,0,MomDummy(2:4,4), MomDummy(1,4), MassDummy(4), 0d0, Spin
endif

if(H_DK.eqv..true.)then
  write(io_LHEOutFile,fmt1) id(5), 2,1,2,0,0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength/ctauUnit, Spin
else
  write(io_LHEOutFile,fmt1) id(5), 1,1,2,0,0,MomDummy(2:4,5), MomDummy(1,5), MassDummy(5), HiggsDKLength/ctauUnit, Spin
endif

if(.not.IsAPhoton(DecayMode1)) then
  write(io_LHEOutFile,fmt1) id(6), 1,3,3,ICOLUP(3,1),ICOLUP(3,2),MomDummy(2:4,6), MomDummy(1,6), MassDummy(6), 0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(7), 1,3,3,ICOLUP(4,1),ICOLUP(4,2),MomDummy(2:4,7), MomDummy(1,7), MassDummy(7), 0.0d0, Spin
endif

if(H_DK.eqv..true.)then
  write(io_LHEOutFile,fmt1) id(8), 1,4,4,501,0,MomDummy(2:4,8), MomDummy(1,8), MassDummy(8), 0.0d0, Spin
  write(io_LHEOutFile,fmt1) id(9), 1,4,4,0,501,MomDummy(2:4,9), MomDummy(1,9), MassDummy(9), 0.0d0, Spin
endif

write(io_LHEOutFile,"(A)") "</event>"


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



SUBROUTINE Kinematics_gg4f_fullproddec(MomExt,ids,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer,parameter :: inTop=1, inBot=2, V1=3, V2=4, Lep1P=5, Lep1M=6, Lep2P=7, Lep2M=8, NUP=8
real(8), intent(in) :: MomExt(1:4,1:NUP)
integer, intent(in) :: ids(1:NUP-2) ! excludes V1, V2
logical, intent(out) :: applyPSCut
integer, intent(inout) :: NBin(:)
integer :: ipart,jpart,lepcount,jetcount,id_ij,idpho_ij
real(8) :: m_4l,pT_part,eta_part,eta_jpart,dphi_ij,deta_ij,dR_ij,m_ij,pT_ij
logical, parameter :: doPrintFailReason=.false.

      applyPSCut = .false.
      lepcount=0
      jetcount=0

      if (doPrintFailReason) then
         write(6,*) "Begin debugging of Kinematics_gg4f_fullproddec"
         do jpart=1,6
            write(6,*) "id(",jpart,")=",convertLHE(ids(jpart))
         enddo
         do jpart=1,NUP
            write(6,*) "MomExt(",jpart,")=",MomExt(:,jpart)
         enddo
      endif

      if ( any(RequestNLeptons.ge.0) ) then
         lepcount = CountLeptons(ids(V2+1-2:NUP-2))
         applyPSCut = ( (RequestNLeptons(1).ge.0 .and. lepcount .lt. RequestNLeptons(1)) .or. (RequestNLeptons(2).ge.0 .and. lepcount .gt. RequestNLeptons(2)) )
         if (applyPSCut) then
            if (doPrintFailReason) write(6,*) "Failed nlep cutoff (nleps=",lepcount,")"
            return
         endif
      endif
      if ( any(RequestNJets.ge.0) ) then
         jetcount = CountJets(ids(V2+1-2:NUP-2))
         applyPSCut = ( (RequestNJets(1).ge.0 .and. jetcount .lt. RequestNJets(1)) .or. (RequestNJets(2).ge.0 .and. jetcount .gt. RequestNJets(2)) )
         if (applyPSCut) then
            if (doPrintFailReason) write(6,*) "Failed njet cutoff (njets=",jetcount,")"
            return
         endif
      endif

      m_4l = get_MInv( MomExt(1:4,Lep1P)+MomExt(1:4,Lep1M)+MomExt(1:4,Lep2P)+MomExt(1:4,Lep2M) )

      do ipart=V2+1,NUP
         pT_part = get_PT(MomExt(1:4,ipart))
         eta_part = get_Eta(MomExt(1:4,ipart))
         if (IsALepton(ids(ipart-2))) then
            if (pT_part.lt.pTlepcut) then
               if (doPrintFailReason) write(6,*) "Failed pTlep cutoff. pT=",pT_part,"<",pTlepcut
               applyPSCut=.true.
               return
            endif
            if (abs(eta_part).gt.etalepcut) then
               if (doPrintFailReason) write(6,*) "Failed etalep cutoff. eta=",eta_part,">",etalepcut
               applyPSCut=.true.
               return
            endif
         endif
         if (IsAJet(ids(ipart-2))) then
            if (pT_part.lt.pTjetcut) then
               if (doPrintFailReason) write(6,*) "Failed pTjet cutoff. pT=",pT_part,"<",pTjetcut
               applyPSCut=.true.
               return
            endif
            if (abs(eta_part).gt.etalepcut) then
               if (doPrintFailReason) write(6,*) "Failed etajet cutoff. eta=",eta_part,">",etajetcut
               applyPSCut=.true.
               return
            endif
         endif

         do jpart=ipart+1,NUP
            dphi_ij = abs( Get_PHI(MomExt(1:4,ipart)) - Get_PHI(MomExt(1:4,jpart)) )
            if( dphi_ij.gt.Pi ) then
               dphi_ij=2d0*Pi-dphi_ij
            endif
            eta_jpart = get_Eta(MomExt(1:4,jpart))
            deta_ij = abs(eta_part - eta_jpart)
            dR_ij = sqrt(deta_ij**2 + dphi_ij**2)
            m_ij = get_MInv(MomExt(1:4,ipart)+MomExt(1:4,jpart))
            pT_ij = get_PT(MomExt(1:4,ipart)+MomExt(1:4,jpart))

            idpho_ij = CoupledVertex((/ids(ipart-2),ids(jpart-2)/),-1,useAHcoupl=1)
            id_ij = CoupledVertex((/ids(ipart-2),ids(jpart-2)/),-1)
            if (idpho_ij.eq.Pho_) then
               if (m_ij.lt.MPhotonCutoff) then
                  if (doPrintFailReason) write(6,*) "Failed mphoton cutoff. m(",ipart,jpart,")=",m_ij,"<",MPhotonCutoff
                  applyPSCut=.true.
                  return
               endif
            endif

            if ( &
               id_ij.eq.Z0_ .and. (                       &
                  (ipart.eq.Lep1P .and. jpart.eq.Lep1M) .or. &
                  (ipart.eq.Lep2P .and. jpart.eq.Lep2M) .or. &
                  includeInterference .and. ( &
                     (ipart.eq.Lep1P .and. jpart.eq.Lep2M) .or. &
                     (ipart.eq.Lep2P .and. jpart.eq.Lep1M)      &
                  ) &
               ) &
               ) then
               if (pT_ij.lt.0.1d0*GeV) then
                  if (doPrintFailReason) write(6,*) "MCFM does not allow p_T^Z<0.1 GeV. pT(",ipart,jpart,")=",pT_ij,"<0.1 GeV"
                  applyPSCut=.true.
                  return
               endif
            endif
            if ( &
               abs(id_ij).eq.abs(Wp_) .and. ( &
                  (ipart.eq.Lep1P .and. jpart.eq.Lep1M) .or. &
                  (ipart.eq.Lep2P .and. jpart.eq.Lep2M) .or. &
                  includeInterference .and. ( &
                     (ipart.eq.Lep1P .and. jpart.eq.Lep2M) .or. &
                     (ipart.eq.Lep2P .and. jpart.eq.Lep1M)      &
                  ) &
               ) &
               ) then
               if (pT_ij.lt.0.05d0*GeV) then
                  if (doPrintFailReason) write(6,*) "MCFM does not allow p_T^W<0.05 GeV. pT(",ipart,jpart,")=",pT_ij,"<0.05 GeV"
                  applyPSCut=.true.
                  return
               endif
            endif

            if( IsAJet(ids(ipart-2)) .and. IsAJet(ids(jpart-2)) ) then
               if (deta_ij.lt.detajetcut .or. (JetsOppositeEta .and. eta_part*eta_jpart.gt.0d0)) then
                  if (doPrintFailReason) write(6,*) "Failed detajet cutoff. deta=|",eta_part,"-",eta_jpart,"|<",detajetcut
                  applyPSCut=.true.
                  return
               endif
               if( m_ij.lt.mJJcut .or. dR_ij.lt.Rjet )  then
                  if (doPrintFailReason) write(6,*) "Failed mjj cutoff. mjj=",m_ij,"<",mJJcut,"or dRjj=",dR_ij,"<",Rjet
                  applyPSCut=.true.
                  return
               endif
            endif
         enddo
      enddo

!     binning
      NBin(:) = 1


   if (doPrintFailReason) write(6,*) "End debugging of Kinematics_gg4f_fullproddec"

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





SUBROUTINE Kinematics_HVBF_fulldecay(MomExt,ids,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:), mZ1, mZ2, MReso, mZ1alt, mZ2alt
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NBin(:),ids(1:8),jpart,lepcount,jetcount
real(8) :: m_jj,y_j1,y_j2,dphi_jj,dy_j1j2,pT_jl,pT_j1,pT_j2,pT_H,m_4l,dR_j1j2
real(8) :: pT_l1,pT_l2,pT_l3,pT_l4,y_l1,y_l2,y_l3,y_l4
real(8) :: Phi1,signPhi1,MomReso(1:4)
logical, parameter :: doPrintFailReason=.false.
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=8, Lep1M=7, Lep2P=10, Lep2M=9


       applyPSCut = .false.
       lepcount=0
       jetcount=0

       if ( any(RequestNLeptons.ge.0) ) then
          lepcount = CountLeptons(ids(3:8))
          applyPSCut = ( (RequestNLeptons(1).ge.0 .and. lepcount .lt. RequestNLeptons(1)) .or. (RequestNLeptons(2).ge.0 .and. lepcount .gt. RequestNLeptons(2)) )
          if (applyPSCut) then
             if (doPrintFailReason) write(6,*) "Failed nlep cutoff (nleps=",lepcount,")"
             return
          endif
       endif
       if ( any(RequestNJets.ge.0) ) then
          jetcount = CountJets(ids(3:8))
          applyPSCut = ( (RequestNJets(1).ge.0 .and. jetcount .lt. RequestNJets(1)) .or. (RequestNJets(2).ge.0 .and. jetcount .gt. RequestNJets(2)) )
          if (applyPSCut) then
             if (doPrintFailReason) write(6,*) "Failed njet cutoff (njets=",jetcount,")"
             return
          endif
       endif

       m_jj = get_MInv( MomExt(1:4,outTop)+MomExt(1:4,outBot) )
       m_4l = get_MInv( MomExt(1:4,Lep1P)+MomExt(1:4,Lep1M)+MomExt(1:4,Lep2P)+MomExt(1:4,Lep2M) )
       y_j1 = get_eta(MomExt(1:4,outTop))
       y_j2 = get_eta(MomExt(1:4,outBot))
       pT_H = get_PT(MomExt(1:4,V1)+MomExt(1:4,V2))
       pT_j1= get_PT(MomExt(1:4,outTop))
       pT_j2= get_PT(MomExt(1:4,outBot))
       pT_l1= get_PT(MomExt(1:4,Lep1M))
       pT_l2= get_PT(MomExt(1:4,Lep1P))
       pT_l3= get_PT(MomExt(1:4,Lep2M))
       pT_l4= get_PT(MomExt(1:4,Lep2P))
       y_l1= get_eta(MomExt(1:4,Lep1M))
       y_l2= get_eta(MomExt(1:4,Lep1P))
       y_l3= get_eta(MomExt(1:4,Lep2M))
       y_l4= get_eta(MomExt(1:4,Lep2P))
       pT_jl = max(pT_j1,pT_j2)
       dy_j1j2 = y_j1 - y_j2
       dR_j1j2 = get_R(MomExt(1:4,outTop), MomExt(1:4,outBot))

       mZ1 = get_MInv(MomExt(1:4,Lep1M)+MomExt(1:4,Lep1P))
       mZ2 = get_MInv(MomExt(1:4,Lep2M)+MomExt(1:4,Lep2P))
       mZ1alt = get_MInv(MomExt(1:4,Lep1M)+MomExt(1:4,Lep2P))
       mZ2alt = get_MInv(MomExt(1:4,Lep2M)+MomExt(1:4,Lep1P))

       dphi_jj = abs( Get_PHI(MomExt(1:4,3)) - Get_PHI(MomExt(1:4,4)) )
       if( dphi_jj.gt.Pi ) dphi_jj=2d0*Pi-dphi_jj


       if (doPrintFailReason) then
         do jpart=1,10
            write(6,*) "MomExt(",jpart,")=",MomExt(:,jpart)
         enddo
       endif

       ! Off-shell photon cut-offs
       if ( m_jj.lt.MPhotonCutoff .and. .not.IsANeutrino(ids(7)) .and. CoupledVertex(ids(7),ids(8)).eq.Z0_) then
          if (doPrintFailReason) write(6,*) "Failed mphoton cutoff. mll=",m_jj,"<",MPhotonCutoff
          applyPSCut=.true.
          return
       endif
       if( mZ1.lt.MPhotonCutoff .and. .not.IsANeutrino(ids(3)) .and. CoupledVertex(ids(3),ids(4)).eq.Z0_ ) then
          if (doPrintFailReason) write(6,*) "Failed mphoton cutoff. mZ1=",mZ1,"<",MPhotonCutoff
          applyPSCut=.true.
          return
       endif
       if( mZ2.lt.MPhotonCutoff .and. .not.IsANeutrino(ids(5)) .and. CoupledVertex(ids(5),ids(6)).eq.Z0_ ) then
          if (doPrintFailReason) write(6,*) "Failed mphoton cutoff. mZ2=",mZ2,"<",MPhotonCutoff
          applyPSCut=.true.
          return
       endif

       if(ids(3).eq.ids(5) .and. ids(4).eq.ids(6)) then
          if( mZ1alt.lt.MPhotonCutoff .and. .not.IsANeutrino(ids(3)) .and. CoupledVertex(ids(3),ids(6)).eq.Z0_ ) then
             if (doPrintFailReason) write(6,*) "Failed mphoton cutoff. mZ1alt=",mZ1alt,"<",MPhotonCutoff
             applyPSCut=.true.
             return
          endif
          if( mZ2alt.lt.MPhotonCutoff .and. .not.IsANeutrino(ids(5)) .and. CoupledVertex(ids(5),ids(4)).eq.Z0_ ) then
             if (doPrintFailReason) write(6,*) "Failed mphoton cutoff. mZ2alt=",mZ2alt,"<",MPhotonCutoff
             applyPSCut=.true.
             return
          endif
       endif

       if( &
         (pT_l1.lt.pTlepcut .and. IsALepton(ids(3))) .or. &
         (pT_l2.lt.pTlepcut .and. IsALepton(ids(4))) .or. &
         (pT_l3.lt.pTlepcut .and. IsALepton(ids(5))) .or. &
         (pT_l4.lt.pTlepcut .and. IsALepton(ids(6))) .or. &
         (pT_j1.lt.pTlepcut .and. IsALepton(ids(7))) .or. &
         (pT_j2.lt.pTlepcut .and. IsALepton(ids(8)))      &
         ) then
          if (doPrintFailReason) write(6,*) "Failed pTlep cutoff. pTls=",pT_l1,pT_l2,pT_l3,pT_l4,pT_j1,pT_j2,"<",pTlepcut
          applyPSCut=.true.
          return
       endif
       if( &
         (pT_l1.lt.pTjetcut .and. IsAJet(ids(3))) .or. &
         (pT_l2.lt.pTjetcut .and. IsAJet(ids(4))) .or. &
         (pT_l3.lt.pTjetcut .and. IsAJet(ids(5))) .or. &
         (pT_l4.lt.pTjetcut .and. IsAJet(ids(6))) .or. &
         (pT_j1.lt.pTjetcut .and. IsAJet(ids(7))) .or. &
         (pT_j2.lt.pTjetcut .and. IsAJet(ids(8)))      &
         ) then
          if (doPrintFailReason) write(6,*) "Failed pTjet cutoff. pTls=",pT_l1,pT_l2,pT_l3,pT_l4,pT_j1,pT_j2,"<",pTjetcut
          applyPSCut=.true.
          return
       endif

       if( &
         (abs(y_l1).gt.etalepcut .and. IsALepton(ids(3))) .or. &
         (abs(y_l2).gt.etalepcut .and. IsALepton(ids(4))) .or. &
         (abs(y_l3).gt.etalepcut .and. IsALepton(ids(5))) .or. &
         (abs(y_l4).gt.etalepcut .and. IsALepton(ids(6))) .or. &
         (abs(y_j1).gt.etalepcut .and. IsALepton(ids(7))) .or. &
         (abs(y_j2).gt.etalepcut .and. IsALepton(ids(8)))      &
         ) then
          if (doPrintFailReason) write(6,*) "Failed ylep cutoff. ys=",y_l1,y_l2,y_l3,y_l4,y_j1,y_j2,">",etalepcut
          applyPSCut=.true.
          return
       endif
       if( &
         (abs(y_l1).gt.etajetcut .and. IsAJet(ids(3))) .or. &
         (abs(y_l2).gt.etajetcut .and. IsAJet(ids(4))) .or. &
         (abs(y_l3).gt.etajetcut .and. IsAJet(ids(5))) .or. &
         (abs(y_l4).gt.etajetcut .and. IsAJet(ids(6))) .or. &
         (abs(y_j1).gt.etajetcut .and. IsAJet(ids(7))) .or. &
         (abs(y_j2).gt.etajetcut .and. IsAJet(ids(8)))      &
         ) then
          if (doPrintFailReason) write(6,*) "Failed yjet cutoff. ys=",y_l1,y_l2,y_l3,y_l4,y_j1,y_j2,">",etajetcut
          applyPSCut=.true.
          return
       endif

        if( IsAJet(ids(7)) .and. (abs(y_j1-y_j2).lt.detajetcut .or. (JetsOppositeEta .and. y_j1*y_j2.gt.0d0)) ) then
           if (doPrintFailReason) write(6,*) "Failed detajet cutoff. deta=",y_j1,"-",y_j2,"<",detajetcut
           applyPSCut=.true.
           return
        endif

        if( IsAJet(ids(7)) .and. (m_jj.lt.mJJcut .or. dR_j1j2.lt.Rjet) )  then
           if (doPrintFailReason) write(6,*) "Failed mjj cutoff. mjj=",m_jj,"<",mJJcut,"or dRjj=",dR_j1j2,"<",Rjet
           applyPSCut=.true.
           return
        endif


! construct Phi1:  angle between beam-scattering plane and the lepton plane of Z1 in the resonance rest frame
      MReso=m_4l

      MomReso(1:4)= MomExt(1:4,V1) + MomExt(1:4,V2)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)

!!      if( dabs( get_MInv(MomExt(1:4,V1)) - get_MInv(MomLept(1:4,1)+MomLept(1:4,2))) .lt.     )

      MomLeptX(1:4,1) = MomExt(1:4,Lep1P)
      MomLeptX(1:4,2) = MomExt(1:4,Lep1M)
      call boost(MomLeptX(1:4,1),MomBoost(1:4),MReso)! boost all leptons into the resonance frame
      call boost(MomLeptX(1:4,2),MomBoost(1:4),MReso)

!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt( MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 )! normalize




      BeamAxis(1:4) = (/1d0,0d0,0d0,1d0/)!  energy components are dummies here and will not be used
      ScatteringAxis(1:4) = MomExt(1:4,V1)  ! Z1
      MomReso(1:4)= MomExt(1:4,V1) + MomExt(1:4,V2)   !Z1+Z2
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
       NBin(:) = 1
       NBin(1)  = WhichBin(1,m_jj)
       NBin(2)  = WhichBin(2,m_jj)
       NBin(3)  = WhichBin(3,pT_H)
       NBin(4)  = WhichBin(4,pT_jl)
       NBin(5)  = WhichBin(5,dy_j1j2)
       NBin(6)  = WhichBin(6,mZ1)
       NBin(7)  = WhichBin(7,mZ2)
       NBin(8)  = WhichBin(8,m_4l)
       NBin(9)  = WhichBin(9,m_4l)
       NBin(10) = WhichBin(10,Phi1)


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




SUBROUTINE Kinematics_VHiggs(id,MomExt,inv_mass,NBin,applyPSCut,useAonshell)
use ModMisc
use ModParameters
implicit none

logical, optional :: useAonshell
logical :: applyPSCut
integer :: NumPart,NBin(:),id(:)
real(8) :: m_jj,y_j1,y_j2,dphi_jj, m_ll, m_inv_V, pt_V, pt_H, pt1, pt2, deltaR, m_Vstar, m_inv_H, costheta1, costheta2, phistar1, phi
double precision MomBoost(1:4), MomFerm(1:4), inv_mass(1:9), MomLeptX(1:4,1:4), ScatteringAxis(1:4), MomReso(1:4), MomZ(1:4), MomFerm1(1:4), MomFerm2(1:4)
double precision MomLeptPlane1(2:4), MomLeptPlane2(2:4), dummy(2:4), signPhi,EHat
double precision, intent(in) :: MomExt(1:4,1:9)! 1:in 2:in 3:V* 4:V 5:H 6,7: q(Z)q(Z) 8,9: q(h)q(H)
logical :: hasAonshell
integer :: i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     hasAonshell = .false.
     if(present(useAonshell)) then
        hasAonshell=useAonshell
     endif

     applyPSCut = .false.
     m_jj = get_MInv(MomExt(1:4,5))
     m_ll = get_MInv(MomExt(1:4,4))
     m_Vstar = get_MInv(MomExt(1:4,3))
     m_inv_H = get_MInv(MomExt(1:4,5))
     m_inv_V = m_ll

     pt_H = get_PT(MomExt(1:4,5))
     pt_V = get_PT(MomExt(1:4,4))

     EHat = MomExt(1,1)+MomExt(1,2)

!     do i=6,7
!        if(IsALHELepton(id(i)))then
!            if(get_PT(MomExt(:,i)).lt.pTlepcut)applyPSCut=.true.
!            !if(get_PT(MomExt(:,i)).lt.5d0*GeV)applyPSCut=.true.
!            if(dabs(get_eta(MomExt(:,i))).gt.etalepcut)applyPSCut=.true.
!        endif
!    enddo

      if(.not.hasAonshell) then
         if(m_ll.le.getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))then
            applyPSCut=.true.
         endif
         if(includeGammaStar .and. .not.IsAWDecay(DecayMode1) .and. (m_ll.lt.MPhotonCutoff .or. m_Vstar.lt.MPhotonCutoff))then
            applyPSCut=.true.
         endif
         if(IsAQuark(convertLHEreverse(id(6)))) then
            pt1 = get_PT(MomExt(1:4,6))
            pt2 = get_PT(MomExt(1:4,7))
            deltaR = get_R(MomExt(1:4,6), MomExt(1:4,7))
            if(m_ll.lt.mJJcut .or. pt1.lt.ptjetcut .or. pt2.lt.ptjetcut .or. deltaR.lt.Rjet) then
               applyPSCut=.true.
            endif
         endif
      else
         if(includeGammaStar .and. m_Vstar.lt.MPhotonCutoff)then
            applyPSCut=.true.
         endif
      endif
      if(H_DK) then
         if(m_jj.le.(getMass(convertLHEreverse(id(8)))+getMass(convertLHEreverse(id(9)))))then
            applyPSCut=.true.
         endif
         if(IsAQuark(convertLHEreverse(id(8)))) then
            pt1 = get_PT(MomExt(1:4,8))
            pt2 = get_PT(MomExt(1:4,9))
            deltaR = get_R(MomExt(1:4,8), MomExt(1:4,9))
            if(m_jj.lt.mJJcut .or. pt1.lt.ptjetcut .or. pt2.lt.ptjetcut .or. deltaR.lt.Rjet) then
               applyPSCut=.true.
            endif
         endif
      else
         if(dabs(m_jj-M_Reso).gt.10d0*Ga_Reso) then
            applyPSCut=.true.
         endif
      endif

!-- FIND PLOTTING BINS
!costheta1 - production angle
!      MomBoost(1)   = -MomExt(1,3)
!      MomBoost(2:4) = +MomExt(2:4,3)
!      if(id(2).lt.0) then
!         MomFerm(1:4)  = -MomExt(1:4,2)
!      else
!         MomFerm(1:4)  = -MomExt(1:4,1)
!      endif
!      call boost(MomFerm(1:4),MomBoost(1:4),inv_mass(3))! boost fermion from Z1 into Z1 rest frame
!      costheta1 = Get_CosAlpha( MomFerm(1:4), -MomExt(1:4,3) )
!
!!      if(id(1).gt.0)then
!!        costheta1=0.5d0
!!      else
!!        costheta1=-0.5d0
!!      endif


      MomZ = MomExt(:,4)
      MomBoost(1)   = +MomExt(1,3)
      MomBoost(2:4) = -MomExt(2:4,3)
      if( (MomExt(2,3)**2+MomExt(3,3)**2+MomExt(4,3)**2).gt.1d-15 )then
        !if(MomExt(4,1)*dble(id(1)).gt.0d0)then
            ScatteringAxis = MomExt(:,1) + MomExt(:,2)
            !ScatteringAxis(1:4) = (/0d0,0d0,0d0,1d0/)
        !else
            !ScatteringAxis =-MomExt(:,1) - MomExt(:,2)
        !endif
        call boost(MomZ(1:4),MomBoost(1:4),m_Vstar)! boost Z into Z* rest frame
      else
        !if(MomExt(4,1)*dble(id(1)).gt.0d0)then
            ScatteringAxis(1:4) = (/0d0,0d0,0d0,1d0/)
        !else
            !ScatteringAxis(1:4) = (/0d0,0d0,0d0,-1d0/)
        !endif
      endif
      costheta1 = Get_CosAlpha( ScatteringAxis, MomZ )



!phistar1
      MomReso(1:4)  = -MomExt(1:4,5)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
      if(id(2).lt.0) then
         MomLeptX(1:4,1) = -MomExt(1:4,2)
         MomLeptX(1:4,2) = -MomExt(1:4,1)
      else
         MomLeptX(1:4,1) = -MomExt(1:4,1)
         MomLeptX(1:4,2) = -MomExt(1:4,2)
      endif
      if(id(8).gt.0) then
         MomLeptX(1:4,3) = MomExt(1:4,8)
         MomLeptX(1:4,4) = MomExt(1:4,9)
      else
         MomLeptX(1:4,3) = MomExt(1:4,9)
         MomLeptX(1:4,4) = MomExt(1:4,8)
      endif
      call boost(MomLeptX(1:4,1),MomBoost(1:4),inv_mass(5))! boost all leptons into the resonance frame
      call boost(MomLeptX(1:4,2),MomBoost(1:4),inv_mass(5))
      call boost(MomLeptX(1:4,3),MomBoost(1:4),inv_mass(5))
      call boost(MomLeptX(1:4,4),MomBoost(1:4),inv_mass(5))
      ScatteringAxis(1:4) = MomLeptX(1:4,3)-MomLeptX(1:4,4)
      ScatteringAxis(1:4) = ScatteringAxis(1:4) / dsqrt(dabs(ScatteringAxis(2)**2+ScatteringAxis(3)**2+ScatteringAxis(4)**2 +1d-15 ))
!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt(dabs(MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15) )! normalize
      MomLeptPlane2(2:4) = (ScatteringAxis(2:4)).cross.(MomLeptX(2:4,1)+MomLeptX(2:4,2))
      MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt(dabs(MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15 ))! normalize
!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
      signPhi = sign(1d0,  (dummy(2)*(MomLeptX(2,1)+MomLeptX(2,2))+dummy(3)*(MomLeptX(3,1)+MomLeptX(3,2))+dummy(4)*(MomLeptX(4,1)+MomLeptX(4,2)))  )! use q1
      Phistar1 = signPhi * acos(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4))


      costheta2=0d0
      phi=0d0
      if(.not.hasAonshell) then

!costheta2 - Z2 decay angle
         MomZ(1:4) = MomExt(1:4,3)
         MomBoost(1)   = +MomExt(1,4)
         MomBoost(2:4) = -MomExt(2:4,4)

         if(id(6).gt.0) then
            MomFerm1(1:4)  = MomExt(1:4,6)
            MomFerm2(1:4)  = MomExt(1:4,7)
         else
            MomFerm1(1:4)  = MomExt(1:4,7)
            MomFerm2(1:4)  = MomExt(1:4,6)
         endif

         call boost(MomFerm1(1:4),MomBoost(1:4),m_inv_V)! boost fermion from Z2 into Z2 rest frame
         call boost(MomZ(1:4),MomBoost(1:4),m_inv_V)
         costheta2 = Get_CosAlpha( MomFerm1(1:4),MomZ(1:4) )

!phi
         MomReso(1:4)  = -MomExt(1:4,5)
         MomBoost(1)   = +MomReso(1)
         MomBoost(2:4) = -MomReso(2:4)
         if(id(2).lt.0) then
            MomLeptX(1:4,1) = -MomExt(1:4,2)
            MomLeptX(1:4,2) = -MomExt(1:4,1)
         else
            MomLeptX(1:4,1) = -MomExt(1:4,1)
            MomLeptX(1:4,2) = -MomExt(1:4,2)
         endif
         if(id(6).gt.0) then
            MomLeptX(1:4,3) = MomExt(1:4,6)
            MomLeptX(1:4,4) = MomExt(1:4,7)
         else
            MomLeptX(1:4,3) = MomExt(1:4,7)
            MomLeptX(1:4,4) = MomExt(1:4,6)
         endif

         call boost(MomLeptX(1:4,1),MomBoost(1:4),m_inv_H)! boost all fermions into the resonance frame
         call boost(MomLeptX(1:4,2),MomBoost(1:4),m_inv_H)
         call boost(MomLeptX(1:4,3),MomBoost(1:4),m_inv_H)
         call boost(MomLeptX(1:4,4),MomBoost(1:4),m_inv_H)
!     orthogonal vectors defined as p(fermion) x p(antifermion)
         MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
         MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt( MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15  )! normalize
         MomLeptPlane2(2:4) = (MomLeptX(2:4,3)).cross.(MomLeptX(2:4,4))! orthogonal vector to lepton plane
         MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt( MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15  )! normalize
!     get the sign
         dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
         signPhi = sign(1d0,  (dummy(2)*(MomLeptX(2,1)+MomLeptX(2,2))+dummy(3)*(MomLeptX(3,1)+MomLeptX(3,2))+dummy(4)*(MomLeptX(4,1)+MomLeptX(4,2)))  )! use q1
         !phi = signPhi * acos(-1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))
         phi = signPhi * acos(1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))

      endif

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
     Nbin(10) = WhichBin(10,EHat)


RETURN
END SUBROUTINE Kinematics_VHiggs







SUBROUTINE Kinematics_eeVH(id,MomExt,NBin,applyPSCut,HDecays,PhoOnshell)
use ModMisc
use ModParameters
use ModVHLO
#if useCollier==1
use ModVHgg
#endif
implicit none

real(8) :: helicity(1:9)
integer, intent(in) :: id(:)
logical, intent(in) :: HDecays
logical, intent(in), optional :: PhoOnshell
logical :: applyPSCut
integer :: NBin(:)
integer :: i,j,k,l,idME2(1:9)
real(8) :: m_jj,y_j1,y_j2,dphi_jj, m_ll, pt_V, pt_H, pt1, pt2, deltaR, m_Vstar, m_inv_H, m_inv_V, costheta1, costheta2, phistar1, phi, eta_Higgs
double precision MomBoost(1:4), MomFerm(1:4), MomLeptX(1:4,1:4), ScatteringAxis(1:4), MomReso(1:4), MomZ(1:4), MomFerm1(1:4), MomFerm2(1:4), MomVStar(1:4)
double precision MomLeptPlane1(2:4), MomLeptPlane2(2:4), dummy(2:4), signPhi
double precision, intent(in) :: MomExt(1:4,1:10)
logical :: hasAonshell

     hasAonshell = .false.
      if(present(PhoOnshell)) then
         hasAonshell=PhoOnshell
      endif

     applyPSCut = .false.
     m_jj = get_MInv(MomExt(1:4,5))
     m_ll = get_MInv(MomExt(1:4,4))
!print*, "m_ll", m_ll
     m_Vstar = get_MInv(MomExt(1:4,3))
     m_inv_H = m_jj
     m_inv_V = m_ll
!print*, "m_inv_V", m_inv_V
     pt_H = get_PT(MomExt(1:4,5))
     pt_V = get_PT(MomExt(1:4,4))


     if(m_Vstar.gt.mVH_minmax(2))applyPSCut=.true.
     if(m_Vstar.lt.mVH_minmax(1))applyPSCut=.true.

     if(pt_H.lt.pTHcut)applyPSCut=.true.

     if(HDecays)then
        do i=8,9
            if(get_PT(MomExt(:,i)).lt.pTjetcut)applyPSCut=.true.
            if(dabs(get_eta(MomExt(:,i))).gt.etajetcut)applyPSCut=.true.
        enddo
        if(m_jj.lt.mJJcut)applyPSCut=.true.
        if(deltaR.lt.Rjet)applyPSCut=.true.
!     else
!        if(dabs(m_jj-M_Reso).gt.10d0*Ga_Reso)then
!            applyPSCut=.true.
!        endif
     endif

     do i=6,7
        if(IsALHELepton(id(i)))then
            if(get_PT(MomExt(:,i)).lt.pTlepcut)applyPSCut=.true.
            !if(get_PT(MomExt(:,i)).lt.5d0*GeV)applyPSCut=.true.
            if(dabs(get_eta(MomExt(:,i))).gt.etalepcut)applyPSCut=.true.
        endif
     enddo

     if(IsALHEQuark(id(6)))then
        do i=6,7
            if(get_PT(MomExt(:,i)).lt.pTjetcut)applyPSCut=.true.
            if(etajetcut.ge.0d0)then
                if(dabs(get_eta(MomExt(:,i))).gt.etajetcut)applyPSCut=.true.
            endif
        enddo
        if(m_ll.lt.mJJcut)applyPSCut=.true.
        if(get_R(MomExt(1:4,6), MomExt(1:4,7)).lt.Rjet)applyPSCut=.true.
     endif

     if(m_ll.lt.m2l_minmax(1))applyPSCut=.true.
     if(m_ll.gt.m2l_minmax(2))applyPSCut=.true.

!     if(.not.hasAonshell) then
!        if(m_ll.le.getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))applyPSCut=.true.
!        if(includeGammaStar .and. .not.IsAWDecay(DecayMode1) .and. (m_ll.lt.MPhotonCutoff .or. m_Vstar.lt.MPhotonCutoff))applyPSCut=.true.
!     else
!        if(includeGammaStar .and. m_Vstar.lt.MPhotonCutoff)applyPSCut=.true.
!     endif

!-- FIND PLOTTING BINS
!costheta1 - production angle
      MomZ = MomExt(:,4)
      !MomBoost(1)   = +MomExt(1,3)
      !MomBoost(2:4) = -MomExt(2:4,3)
      !if( (MomExt(2,3)**2+MomExt(3,3)**2+MomExt(4,3)**2).gt.1d-8 )then
        !ScatteringAxis = MomExt(:,1) + MomExt(:,2)
      !  ScatteringAxis = MomExt(:,3)
      !  call boost(MomZ(1:4),MomBoost(1:4),m_Vstar)! boost Z into Z* rest frame
      !else
        ScatteringAxis(1:4) = (/0d0,0d0,0d0,1d0/)
      !endif
      costheta1 = Get_CosAlpha( ScatteringAxis, MomZ )


!costheta2 - Z2 decay angle
         MomZ(1:4) = MomExt(1:4,3)
         MomBoost(1)   = +MomExt(1,4)
         MomBoost(2:4) = -MomExt(2:4,4)

         if(id(6).gt.0) then
            MomFerm1(1:4)  = MomExt(1:4,6)
            MomFerm2(1:4)  = MomExt(1:4,7)
         else
            MomFerm1(1:4)  = MomExt(1:4,7)
            MomFerm2(1:4)  = MomExt(1:4,6)
         endif

         call boost(MomFerm1(1:4),MomBoost(1:4),m_inv_V)! boost fermion from Z2 into Z2 rest frame
         call boost(MomZ(1:4),MomBoost(1:4),m_inv_V)
         costheta2 = Get_CosAlpha( MomFerm1(1:4),MomZ(1:4) )



!phistar1, which is really phi1 in arXiv:1309.4819 [hep-ph] FIG.1 middle
      MomReso(1:4)  = -MomExt(1:4,5)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
      !if(id(2).lt.0) then
      !   MomLeptX(1:4,1) = -MomExt(1:4,2)
      !   MomLeptX(1:4,2) = -MomExt(1:4,1)
      !else
         MomLeptX(1:4,1) = -MomExt(1:4,1)
         MomLeptX(1:4,2) = -MomExt(1:4,2)
      !endif
      !if(id(8).gt.0) then
         MomLeptX(1:4,3) = MomExt(1:4,8)
         MomLeptX(1:4,4) = MomExt(1:4,9)
      !else
      !   MomLeptX(1:4,3) = MomExt(1:4,9)
      !   MomLeptX(1:4,4) = MomExt(1:4,8)
      !endif
      call boost(MomLeptX(1:4,1),MomBoost(1:4),m_inv_H)! boost all fermions into the resonance frame
      call boost(MomLeptX(1:4,2),MomBoost(1:4),m_inv_H)
      call boost(MomLeptX(1:4,3),MomBoost(1:4),m_inv_H)
      call boost(MomLeptX(1:4,4),MomBoost(1:4),m_inv_H)
      ScatteringAxis(1:4) = MomLeptX(1:4,3)-MomLeptX(1:4,4)
      ScatteringAxis(1:4) = ScatteringAxis(1:4) / dsqrt(dabs(ScatteringAxis(2)**2+ScatteringAxis(3)**2+ScatteringAxis(4)**2 +1d-15 ))
!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt(dabs(MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15) )! normalize
      MomLeptPlane2(2:4) = (ScatteringAxis(2:4)).cross.(MomLeptX(2:4,1)+MomLeptX(2:4,2))
      MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt(dabs(MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15 ))! normalize
!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
      signPhi = sign(1d0,  (dummy(2)*(MomLeptX(2,1)+MomLeptX(2,2))+dummy(3)*(MomLeptX(3,1)+MomLeptX(3,2))+dummy(4)*(MomLeptX(4,1)+MomLeptX(4,2)))  )! use q1
      Phistar1 = signPhi * acos(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4))


!phi
         MomReso(1:4)  = -MomExt(1:4,5)
         MomBoost(1)   = +MomReso(1)
         MomBoost(2:4) = -MomReso(2:4)
         !if(id(2).lt.0) then
         !   MomLeptX(1:4,1) = -MomExt(1:4,2)
         !   MomLeptX(1:4,2) = -MomExt(1:4,1)
         !else
            MomLeptX(1:4,1) = -MomExt(1:4,1)
            MomLeptX(1:4,2) = -MomExt(1:4,2)
         !endif
         if(id(6).gt.0) then
            MomLeptX(1:4,3) = MomExt(1:4,6)
            MomLeptX(1:4,4) = MomExt(1:4,7)
         else
            MomLeptX(1:4,3) = MomExt(1:4,7)
            MomLeptX(1:4,4) = MomExt(1:4,6)
         endif

         call boost(MomLeptX(1:4,1),MomBoost(1:4),m_inv_H)! boost all fermions into the resonance frame
         call boost(MomLeptX(1:4,2),MomBoost(1:4),m_inv_H)
         call boost(MomLeptX(1:4,3),MomBoost(1:4),m_inv_H)
         call boost(MomLeptX(1:4,4),MomBoost(1:4),m_inv_H)
!     orthogonal vectors defined as p(fermion) x p(antifermion)
         MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
         MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt( MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15  )! normalize
         MomLeptPlane2(2:4) = (MomLeptX(2:4,3)).cross.(MomLeptX(2:4,4))! orthogonal vector to lepton plane
         MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt( MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15  )! normalize
!     get the sign
         dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
         signPhi = sign(1d0,  (dummy(2)*(MomLeptX(2,1)+MomLeptX(2,2))+dummy(3)*(MomLeptX(3,1)+MomLeptX(3,2))+dummy(4)*(MomLeptX(4,1)+MomLeptX(4,2)))  )! use q1
         !phi = signPhi * acos(-1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))
         phi = signPhi * acos(1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))

!      endif

! eta(H)
     eta_Higgs = get_eta(MomExt(1:4,5))


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
     Nbin(10)  = WhichBin(10,eta_Higgs)

RETURN
END SUBROUTINE Kinematics_eeVH










SUBROUTINE Kinematics_ppVH(id,Mom,NBin,applyPSCut,HDecays,PhoOnshell)
use ModMisc
use ModParameters
use ModVHLO
#if useCollier==1
use ModVHgg
#endif
implicit none

real(8) :: helicity(1:9)
integer, intent(in) :: id(:)
logical, intent(in) :: HDecays
logical, intent(in), optional :: PhoOnshell
logical :: applyPSCut
integer :: NBin(:)
integer :: i,j,k,l,idME2(1:9)
real(8) :: m_jj,y_j1,y_j2,dphi_jj, m_ll, pt_V, pt_H, pt1, pt2, deltaR, m_Vstar, m_inv_H, m_inv_V, costheta1, costheta2, phistar1, phi, eta_Higgs
real(8) :: MomBoost(1:4), MomFerm(1:4), MomLeptX(1:4,1:4), ScatteringAxis(1:4), MomReso(1:4), MomZ(1:4), MomFerm1(1:4), MomFerm2(1:4), MomVStar(1:4)
real(8) :: MomLeptPlane1(2:4), MomLeptPlane2(2:4), dummy(2:4), signPhi
real(8), intent(in) :: Mom(1:4,1:10)
real(8) :: MomExt(1:4,1:10)
logical :: hasAonshell
real(8) :: x1,x2,momx1x2(1:4,1:3),momab(1:4,1:2),momb(1:4),ehat
real(8) :: MomBoostT(1:4),m_Vstar_L,m_Vstar_T
!real(8) :: eta_Vstar

     MomExt=Mom!keep original copy of Mom.

     hasAonshell = .false.
      if(present(PhoOnshell)) then
         hasAonshell=PhoOnshell
      endif

     applyPSCut = .false.
     m_jj = get_MInv(MomExt(1:4,5))
     m_ll = get_MInv(MomExt(1:4,4))
!print*, "m_ll", m_ll
     m_Vstar = get_MInv(MomExt(1:4,3))
     m_inv_H = m_jj
     m_inv_V = m_ll
!print*, "m_inv_V", m_inv_V
     pt_H = get_PT(MomExt(1:4,5))
     pt_V = get_PT(MomExt(1:4,4))

     if(m_Vstar.gt.mVH_minmax(2))applyPSCut=.true.
     if(m_Vstar.lt.mVH_minmax(1))applyPSCut=.true.

     if(pt_H.lt.pTHcut)applyPSCut=.true.

     if(HDecays)then
        do i=8,9
            if(get_PT(MomExt(:,i)).lt.pTjetcut)applyPSCut=.true.
            if(dabs(get_eta(MomExt(:,i))).gt.etajetcut)applyPSCut=.true.
        enddo
        if(m_jj.lt.mJJcut)applyPSCut=.true.
        if(deltaR.lt.Rjet)applyPSCut=.true.
!     else
!        if(dabs(m_jj-M_Reso).gt.10d0*Ga_Reso)then
!            applyPSCut=.true.
!        endif
     endif

     do i=6,7
        if(IsALHELepton(id(i)))then
            if(get_PT(MomExt(:,i)).lt.pTlepcut)applyPSCut=.true.
            !if(get_PT(MomExt(:,i)).lt.5d0*GeV)applyPSCut=.true.
            if(dabs(get_eta(MomExt(:,i))).gt.etalepcut)applyPSCut=.true.
        endif
     enddo

     if(IsALHEQuark(id(6)))then
        do i=6,7
            if(get_PT(MomExt(:,i)).lt.pTjetcut)applyPSCut=.true.
            if(etajetcut.ge.0d0)then
                if(dabs(get_eta(MomExt(:,i))).gt.etajetcut)applyPSCut=.true.
            endif
        enddo
        if(m_ll.lt.mJJcut)applyPSCut=.true.
        if(get_R(MomExt(1:4,6), MomExt(1:4,7)).lt.Rjet)applyPSCut=.true.
     endif

     if(m_ll.lt.m2l_minmax(1))applyPSCut=.true.
     if(m_ll.gt.m2l_minmax(2))applyPSCut=.true.

!     if(.not.hasAonshell) then
!        if(m_ll.le.getMass(convertLHEreverse(id(6)))+getMass(convertLHEreverse(id(7))))applyPSCut=.true.
!        if(includeGammaStar .and. .not.IsAWDecay(DecayMode1) .and. (m_ll.lt.MPhotonCutoff .or. m_Vstar.lt.MPhotonCutoff))applyPSCut=.true.
!     else
!        if(includeGammaStar .and. m_Vstar.lt.MPhotonCutoff)applyPSCut=.true.
!     endif

      m_Vstar_T=dsqrt(mom(1,3)**2-mom(2,3)**2-mom(3,3)**2)

      if(Get_PT2(momext(:,3)).gt.1d-8)then
        MomBoostT(1:4)=(/mom(1,3),-mom(2,3),-mom(3,3),0d0/)!momentum/vector used to boost away the transverse moment of VH
        call boost(momext(1:4,3),MomBoostT(1:4),m_Vstar_T)
      endif

!print*,"---------"
!print*,mom(1:4,3)
!print*,momext(1:4,3)

      m_Vstar_L=get_Minv(momext(1:4,3))

      call x1x2(1,momext(1:4,3),x1,x2)

!print*,m_Vstar_L

      momext(:,1)=0.5d0*m_Vstar_L * (/+1d0,0d0,0d0,+1d0/)
      momext(:,2)=0.5d0*m_Vstar_L * (/+1d0,0d0,0d0,-1d0/)
      call boost2lab(x1,x2,2,momext(1:4,1:2))

!print*,momext(1:4,1),"1"
!print*,momext(4,3),x1-x2
!print*,momext(1:4,1)+momext(1:4,2)-momext(1:4,3)

      if(Get_PT2(mom(:,3)).gt.1d-8)then

        MomBoostT(1:4)=(/mom(1,3),-mom(2,3),-mom(3,3),0d0/)!momentum/vector used to boost away the transverse moment of VH

        if((.not.HbbDecays).and.hasAonshell)then   ! if H, V
          call boost(momext(1:4,4),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,5),MomBoostT(1:4),m_Vstar_T)
        elseif(HbbDecays.and.(.not.hasAonshell))then ! if H > bb, V > ff
          call boost(momext(1:4,6),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,7),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,8),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,9),MomBoostT(1:4),m_Vstar_T)
          momext(:,4)=momext(:,6)+momext(:,7)
          momext(:,5)=momext(:,8)+momext(:,9)
        elseif((.not.HbbDecays).and.(.not.hasAonshell))then  ! if H, V > ff
          call boost(momext(1:4,5),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,6),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,7),MomBoostT(1:4),m_Vstar_T)
          momext(:,4)=momext(:,6)+momext(:,7)
        else                                       ! if H > bb, V
          call boost(momext(1:4,4),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,8),MomBoostT(1:4),m_Vstar_T)
          call boost(momext(1:4,9),MomBoostT(1:4),m_Vstar_T)
          momext(:,5)=momext(:,8)+momext(:,9)
        endif

      endif

!print*,-momext(:,6)-momext(:,7)-momext(:,8)-momext(:,9)
!print*,momext(:,1)+momext(:,2)-momext(:,6)-momext(:,7)-momext(:,8)-momext(:,9)

!-- FIND PLOTTING BINS
!costheta1 - production angle
      MomZ = MomExt(:,4)
      MomBoost(1)   = +MomExt(1,3)
      MomBoost(2:4) = -MomExt(2:4,3)
      if( (MomExt(2,3)**2+MomExt(3,3)**2+MomExt(4,3)**2).gt.1d-8 )then
        !ScatteringAxis = MomExt(:,1) + MomExt(:,2)
        ScatteringAxis = MomExt(:,3)
        call boost(MomZ(1:4),MomBoost(1:4),m_Vstar)! boost Z into Z* rest frame
      else
        ScatteringAxis(1:4) = (/0d0,0d0,0d0,1d0/)
      endif
      costheta1 = Get_CosAlpha( ScatteringAxis, MomZ )


!costheta2 - Z2 decay angle
         MomZ(1:4) = MomExt(1:4,3)
         MomBoost(1)   = +MomExt(1,4)
         MomBoost(2:4) = -MomExt(2:4,4)

         if(id(6).gt.0) then
            MomFerm1(1:4)  = MomExt(1:4,6)
            MomFerm2(1:4)  = MomExt(1:4,7)
         else
            MomFerm1(1:4)  = MomExt(1:4,7)
            MomFerm2(1:4)  = MomExt(1:4,6)
         endif

         call boost(MomFerm1(1:4),MomBoost(1:4),m_inv_V)! boost fermion from Z2 into Z2 rest frame
         call boost(MomZ(1:4),MomBoost(1:4),m_inv_V)
         costheta2 = Get_CosAlpha( MomFerm1(1:4),MomZ(1:4) )



!phistar1, which is really phi1 in arXiv:1309.4819 [hep-ph] FIG.1 middle
      MomReso(1:4)  = -MomExt(1:4,5)
      MomBoost(1)   = +MomReso(1)
      MomBoost(2:4) = -MomReso(2:4)
      !if(id(2).lt.0) then
      !   MomLeptX(1:4,1) = -MomExt(1:4,2)
      !   MomLeptX(1:4,2) = -MomExt(1:4,1)
      !else
         MomLeptX(1:4,1) = -MomExt(1:4,1)
         MomLeptX(1:4,2) = -MomExt(1:4,2)
      !endif
      !if(id(8).gt.0) then
         MomLeptX(1:4,3) = MomExt(1:4,8)
         MomLeptX(1:4,4) = MomExt(1:4,9)
      !else
      !   MomLeptX(1:4,3) = MomExt(1:4,9)
      !   MomLeptX(1:4,4) = MomExt(1:4,8)
      !endif
      call boost(MomLeptX(1:4,1),MomBoost(1:4),m_inv_H)! boost all fermions into the resonance frame
      call boost(MomLeptX(1:4,2),MomBoost(1:4),m_inv_H)
      call boost(MomLeptX(1:4,3),MomBoost(1:4),m_inv_H)
      call boost(MomLeptX(1:4,4),MomBoost(1:4),m_inv_H)
      ScatteringAxis(1:4) = MomLeptX(1:4,3)-MomLeptX(1:4,4)
      ScatteringAxis(1:4) = ScatteringAxis(1:4) / dsqrt(dabs(ScatteringAxis(2)**2+ScatteringAxis(3)**2+ScatteringAxis(4)**2 +1d-15 ))
!     orthogonal vectors defined as p(fermion) x p(antifermion)
      MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
      MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt(dabs(MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15) )! normalize
      MomLeptPlane2(2:4) = (ScatteringAxis(2:4)).cross.(MomLeptX(2:4,1)+MomLeptX(2:4,2))
      MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt(dabs(MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15 ))! normalize
!     get the sign
      dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
      signPhi = sign(1d0,  (dummy(2)*(MomLeptX(2,1)+MomLeptX(2,2))+dummy(3)*(MomLeptX(3,1)+MomLeptX(3,2))+dummy(4)*(MomLeptX(4,1)+MomLeptX(4,2)))  )! use q1
      Phistar1 = signPhi * acos(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4))


!phi
         MomReso(1:4)  = -MomExt(1:4,5)
         MomBoost(1)   = +MomReso(1)
         MomBoost(2:4) = -MomReso(2:4)
         !if(id(2).lt.0) then
         !   MomLeptX(1:4,1) = -MomExt(1:4,2)
         !   MomLeptX(1:4,2) = -MomExt(1:4,1)
         !else
            MomLeptX(1:4,1) = -MomExt(1:4,1)
            MomLeptX(1:4,2) = -MomExt(1:4,2)
         !endif
         if(id(6).gt.0) then
            MomLeptX(1:4,3) = MomExt(1:4,6)
            MomLeptX(1:4,4) = MomExt(1:4,7)
         else
            MomLeptX(1:4,3) = MomExt(1:4,7)
            MomLeptX(1:4,4) = MomExt(1:4,6)
         endif

         call boost(MomLeptX(1:4,1),MomBoost(1:4),m_inv_H)! boost all fermions into the resonance frame
         call boost(MomLeptX(1:4,2),MomBoost(1:4),m_inv_H)
         call boost(MomLeptX(1:4,3),MomBoost(1:4),m_inv_H)
         call boost(MomLeptX(1:4,4),MomBoost(1:4),m_inv_H)
!     orthogonal vectors defined as p(fermion) x p(antifermion)
         MomLeptPlane1(2:4) = (MomLeptX(2:4,1)).cross.(MomLeptX(2:4,2))! orthogonal vector to lepton plane
         MomLeptPlane1(2:4) = MomLeptPlane1(2:4)/dsqrt( MomLeptPlane1(2)**2+MomLeptPlane1(3)**2+MomLeptPlane1(4)**2 +1d-15  )! normalize
         MomLeptPlane2(2:4) = (MomLeptX(2:4,3)).cross.(MomLeptX(2:4,4))! orthogonal vector to lepton plane
         MomLeptPlane2(2:4) = MomLeptPlane2(2:4)/dsqrt( MomLeptPlane2(2)**2+MomLeptPlane2(3)**2+MomLeptPlane2(4)**2 +1d-15  )! normalize
!     get the sign
         dummy(2:4) = (MomLeptPlane1(2:4)).cross.(MomLeptPlane2(2:4))
         signPhi = sign(1d0,  (dummy(2)*(MomLeptX(2,1)+MomLeptX(2,2))+dummy(3)*(MomLeptX(3,1)+MomLeptX(3,2))+dummy(4)*(MomLeptX(4,1)+MomLeptX(4,2)))  )! use q1
         !phi = signPhi * acos(-1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))
         phi = signPhi * acos(1d0*(MomLeptPlane1(2)*MomLeptPlane2(2) + MomLeptPlane1(3)*MomLeptPlane2(3) + MomLeptPlane1(4)*MomLeptPlane2(4)))

!      endif

! eta(H)
     eta_Higgs = get_eta(Mom(1:4,5))

! eta(V*)
!     eta_Vstar = get_eta(Mom(1:4,3))


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
     Nbin(10)  = WhichBin(10,eta_Higgs)
!     Nbin(11)  = WhichBin(11,eta_Vstar)

RETURN
END SUBROUTINE Kinematics_ppVH









subroutine x1x2(NumPart,p_final,x1,x2)
use modParameters
implicit none
integer, intent(in) :: NumPart
real(8), intent(in) :: p_final(1:4,1:NumPart)
real(8), intent(out) :: x1,x2
real(8) :: Etot, pztot
integer :: i
!print*, p_final
Etot=0d0
pztot=0d0
do i=1,NumPart
    Etot=Etot+p_final(1,i)
    pztot=pztot+p_final(4,i)
enddo

x1 = (Etot+Pztot)/Collider_Energy
x2 = (Etot-Pztot)/Collider_Energy
!print*,x1,x2

end subroutine x1x2











SUBROUTINE Kinematics_HH(id,Mom,NBin,applyPSCut)
  use ModMisc
  use ModParameters
  implicit none

  logical :: applyPSCut
  double precision, intent(in) :: Mom(1:4,1:9)
  integer :: NBin(:),id(:)
  real(8) :: m_H1,m_H2,pt_H1,pt_H2,m_HH,costhetastar,costheta1,costheta2,phi1,phi
  real(8) :: deltaR1,deltaR2
  real(8) :: Mom_temp(1:4,4:9),Mom_hist(1:4,1:9)
  real(8) :: MomH1(1:4),MomH1oppo(1:4),MomH2oppo(1:4),MomBoost(1:4),MomFerm(1:4),Momb(1:4,1:4),ScatteringAxis(1:4),BeamAxis(1:4)
  real(8) :: MombbPlane1(2:4),MombbPlane2(2:4),dummy(2:4),MomBeamScatterPlane(2:4)
  real(8) :: signPhi,signPhi1
  integer :: i

  applyPSCut = .false.
  Mom_hist(1:4,1:9)=Mom(1:4,1:9)
  m_HH = get_MInv(Mom_hist(1:4,3))
  m_H1 = get_MInv(Mom_hist(1:4,4))
  m_H2 = get_MInv(Mom_hist(1:4,5))
  if(m_H2.gt.m_H1)then
    Mom_temp(1:4,4:9)=Mom(1:4,4:9)
    Mom_hist(:,4)=Mom_temp(:,5)
    Mom_hist(:,5)=Mom_temp(:,4)
    Mom_hist(:,6)=Mom_temp(:,8)
    Mom_hist(:,7)=Mom_temp(:,9)
    Mom_hist(:,8)=Mom_temp(:,6)
    Mom_hist(:,9)=Mom_temp(:,7)
  endif
  pt_H1 = get_PT(Mom_hist(1:4,4))
  pt_H2 = get_PT(Mom_hist(1:4,5))

  if(m_H1.lt.mJJcut .or. m_H2.lt.mJJcut)applyPSCut=.true.

  deltaR1 = get_R(Mom_hist(1:4,6), Mom_hist(1:4,7))
  deltaR2 = get_R(Mom_hist(1:4,8), Mom_hist(1:4,9))
  if(deltaR1.lt.Rjet .or. deltaR2.lt.Rjet)applyPSCut=.true.

  do i=6,9
    if(get_PT(Mom_hist(:,i)).lt.pTjetcut)applyPSCut=.true.
    if(get_eta(Mom_hist(:,i)).gt.etajetcut)applyPSCut=.true.
  enddo

!-- FIND PLOTTING BINS
!costheta* - production angle
  MomBoost(1)   = -Mom_hist(1,3)
  MomBoost(2:4) = +Mom_hist(2:4,3)
  MomH1(1:4)=Mom_hist(:,4)
  call boost(MomH1(1:4),MomBoost(1:4),m_HH)! boost H1 from H1 into HH rest frame
  costhetastar = Get_CosTheta( MomH1(1:4) )

!costheta1 - H1 > b6 b~7 decay angle
  MomBoost(1)   = +Mom_hist(1,3)
  MomBoost(2:4) = -Mom_hist(2:4,3)
  MomFerm(1:4)  = Mom_hist(1:4,6)
  call boost(MomFerm(1:4),MomBoost(1:4),m_H1)! boost b from H1 into H1 rest frame

  MomH1oppo(1) = Mom_hist(1,5)
  MomH1oppo(2:4) = -Mom_hist(2:4,5)
  call boost(MomH1oppo(1:4),MomBoost(1:4),m_H1)! boost -H2 into H1 rest frame
  CosTheta1 = Get_CosAlpha( MomFerm(1:4),MomH1oppo(1:4) )

!costheta1 - H2 > b8 b~9 decay angle
  MomBoost(1)   = +Mom_hist(1,3)
  MomBoost(2:4) = -Mom_hist(2:4,3)
  MomFerm(1:4)  = Mom_hist(1:4,8)
  call boost(MomFerm(1:4),MomBoost(1:4),m_H2)! boost b from H2 into H2 rest frame

  MomH2oppo(1) = Mom_hist(1,4)
  MomH2oppo(2:4) = -Mom_hist(2:4,4)
  call boost(MomH2oppo(1:4),MomBoost(1:4),m_H2)! boost -H1 into H2 rest frame
  CosTheta2 = Get_CosAlpha( MomFerm(1:4),MomH2oppo(1:4) )

!phi
  MomBoost(1)   = +Mom_hist(1,3)
  MomBoost(2:4) = -Mom_hist(2:4,3)
  Momb(1:4,1:4) = Mom_hist(1:4,6:9)
  ScatteringAxis(1:4) = Mom_hist(1:4,3)
  call boost(Momb(1:4,1),MomBoost(1:4),m_HH)! boost all leptons into the 1+2 frame
  call boost(Momb(1:4,2),MomBoost(1:4),m_HH)
  call boost(Momb(1:4,3),MomBoost(1:4),m_HH)
  call boost(Momb(1:4,4),MomBoost(1:4),m_HH)
  call boost(ScatteringAxis(1:4),MomBoost(1:4),m_HH)
!     orthogonal vectors defined as p(b) x p(b~)
  MombbPlane1(2:4) = (Momb(2:4,1)).cross.(Momb(2:4,2))! orthogonal vector to bb~ plane
  MombbPlane1(2:4) = MombbPlane1(2:4)/dsqrt( MombbPlane1(2)**2+MombbPlane1(3)**2+MombbPlane1(4)**2 )! normalize
  MombbPlane2(2:4) = (Momb(2:4,3)).cross.(Momb(2:4,4))! orthogonal vector to bb~ plane
  MombbPlane2(2:4) = MombbPlane2(2:4)/dsqrt( MombbPlane2(2)**2+MombbPlane2(3)**2+MombbPlane2(4)**2 )! normalize
!     get the sign
  dummy(2:4) = (MombbPlane1(2:4)).cross.(MombbPlane2(2:4))
  signPhi = sign(1d0,  (dummy(2)*ScatteringAxis(2)+dummy(3)*ScatteringAxis(3)+dummy(4)*ScatteringAxis(4))  )! use q1
  Phi = signPhi * acos(-1d0*(MombbPlane1(2)*MombbPlane2(2) + MombbPlane1(3)*MombbPlane2(3) + MombbPlane1(4)*MombbPlane2(4)))

! construct Phi1:  angle between beam-scattering plane and the lepton plane of H1 in the resonance rest frame

  BeamAxis(1:4) = (/1d0,0d0,0d0,1d0/)!  energy components are dummies here and will not be used
  ScatteringAxis(1:4) = Mom_hist(1:4,4)
  MomBoost(1)   = +Mom_hist(1,3)
  MomBoost(2:4) = -Mom_hist(2:4,3)
  call boost(ScatteringAxis(1:4),MomBoost(1:4),m_HH)

  MomBeamScatterPlane(2:4) = (BeamAxis(2:4)).cross.(ScatteringAxis(2:4))! orthogonal vector to beam-scattering plane
  MomBeamScatterPlane(2:4) = MomBeamScatterPlane(2:4)/dsqrt( MomBeamScatterPlane(2)**2+MomBeamScatterPlane(3)**2+MomBeamScatterPlane(4)**2 )
!     get the sign
  dummy(2:4) = (MombbPlane1(2:4)).cross.(MomBeamScatterPlane(2:4))
  signPhi1 = sign(1d0,  (dummy(2)*ScatteringAxis(2)+dummy(3)*ScatteringAxis(3)+dummy(4)*ScatteringAxis(4))  )! use q1
  Phi1 = signPhi1 * acos(MombbPlane1(2)*MomBeamScatterPlane(2) + MombbPlane1(3)*MomBeamScatterPlane(3) + MombbPlane1(4)*MomBeamScatterPlane(4))

!     binning
     NBin(1)  = WhichBin(1,m_H1)
     NBin(2)  = WhichBin(2,m_H2)
     NBin(3)  = WhichBin(3,pt_H1)
     NBin(4)  = WhichBin(4,pt_H2)
     NBin(5)  = WhichBin(5,m_HH)
     Nbin(6)  = WhichBin(6,costhetastar)
     Nbin(7)  = WhichBin(7,costheta1)
     Nbin(8)  = WhichBin(8,costheta2)
     Nbin(9)  = WhichBin(9,phi1)
     Nbin(10)  = WhichBin(10,phi)

RETURN
END SUBROUTINE Kinematics_HH




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

SUBROUTINE Kinematics_TWH(Mom,applyPSCut,NBin)
use ModParameters
use ModMisc
! use modTTBH
implicit none
real(8) :: Mom(1:4,1:11),MomMELA(1:4,1:11)
logical :: applyPSCut
integer :: NBin(:)
real(8) :: pT_t,pT_H,pT_Wm,MatElSq_H0,MatElSq_H1,D_0minus
real(8) :: mt,mWm,mtWm,mWp,pT_b,pT_l,pT_lm,pT_miss,y_top,y_Higgs,y_Wm
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4,Wm=5,  b=6,Wp=7,lepP=8,nu=9,  lepM=10,nubar=11
logical,save :: FirstTime=.true.


    applyPSCut = .false.

    pT_t = get_PT(Mom(1:4,t))
    pT_Wm = get_PT(Mom(1:4,Wm))
    pT_H = get_PT(Mom(1:4,Hbos))
    pT_b = get_PT(Mom(1:4,b))
    pT_l = get_PT(Mom(1:4,LepP))
    pT_lm = get_PT(Mom(1:4,LepM))
    pT_miss = get_PT(Mom(1:4,nu)+Mom(1:4,nubar))
    mt = get_MInv(Mom(1:4,t))
    mWm = get_MInv(Mom(1:4,Wm))
    mtWm = get_MInv(Mom(1:4,t)+Mom(1:4,Wm))
    mWp = get_MInv(Mom(1:4,Wp))
    mWm = get_MInv(Mom(1:4,Wm))
    
    y_top=get_eta(Mom(1:4,t))
    y_Higgs=get_eta(Mom(1:4,Hbos))
    y_Wm=get_eta(Mom(1:4,Wm))

    if( m_Top.lt.10d0*GeV  .and. (pT_t.lt.pTjetcut) ) applyPSCut=.true.


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
    NBin(4)  = WhichBin(4,mWm)
    NBin(5)  = WhichBin(5,mWp)
    NBin(6)  = WhichBin(6,pT_b)
    NBin(7)  = WhichBin(7,pT_l)    
    NBin(8)  = WhichBin(8,pT_lm)
    NBin(9)  = WhichBin(9,pT_miss)
    NBin(10) = WhichBin(10,y_top)
    NBin(11)  = WhichBin(11,y_Wm)
    NBin(12)  = WhichBin(12,y_higgs)
    NBin(13) = WhichBin(13,D_0minus)



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

FUNCTION ZLepBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZLepBranching_flat


  if( xRnd .le. 0.5d0 ) then
      ZLepBranching_flat = ElM_
  else
      ZLepBranching_flat = MuM_
  endif

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

FUNCTION ZLepPlusTauBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZLepPlusTauBranching_flat

  if( xRnd .le. (1d0/3d0) ) then
      ZLepPlusTauBranching_flat = ElM_
  elseif(xRnd .le. (2d0/3d0) ) then
      ZLepPlusTauBranching_flat = MuM_
  else
      ZLepPlusTauBranching_flat = TaM_
  endif

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

FUNCTION ZNuBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZNuBranching_flat

  if( xRnd .le. (1d0/3d0) ) then
      ZNuBranching_flat = NuE_
  elseif(xRnd .le. (2d0/3d0) ) then
      ZNuBranching_flat = NuM_
  else
      ZNuBranching_flat = NuT_
  endif

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
  elseif(xRnd .le. 2d0*yy ) then
      ZAnyBranching_flat = Chm_
  elseif(xRnd .le. 3d0*yy ) then
      ZAnyBranching_flat = Dn_
  elseif(xRnd .le. 4d0*yy ) then
      ZAnyBranching_flat = Str_
  elseif(xRnd .le. 5d0*yy ) then
      ZAnyBranching_flat = Bot_
  elseif(xRnd .le. 5d0*yy + xx ) then
      ZAnyBranching_flat = ElM_
  elseif(xRnd .le. 5d0*yy + 2d0*xx ) then
      ZAnyBranching_flat = MuM_
  elseif(xRnd .le. 5d0*yy + 3d0*xx ) then
      ZAnyBranching_flat = TaM_
  elseif(xRnd .le. 5d0*yy + 4d0*xx ) then
      ZAnyBranching_flat = NuE_
  elseif(xRnd .le. 5d0*yy + 5d0*xx ) then
      ZAnyBranching_flat = NuM_
  elseif(xRnd .le. 5d0*yy + 6d0*xx ) then
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

FUNCTION WLepBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WLepBranching_flat

  if( xRnd .le. 0.5d0 ) then
      WLepBranching_flat = ElM_
  else
      WLepBranching_flat = MuM_
  endif

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

FUNCTION WLepPlusTauBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WLepPlusTauBranching_flat

  if( xRnd .le. (1d0/3d0) ) then
      WLepPlusTauBranching_flat = ElM_
  elseif(xRnd .le. (2d0/3d0) ) then
      WLepPlusTauBranching_flat = MuM_
  else
      WLepPlusTauBranching_flat = TaM_
  endif

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

FUNCTION WQuaUpBranching_flat(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WQuaUpBranching_flat

  if( xRnd .le. 0.5d0 ) then
      WQuaUpBranching_flat = Up_
  else
      WQuaUpBranching_flat = Chm_
  endif

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


  if( xRnd .le.        yy ) then
      WAnyBranching_flat = Up_
  elseif(xRnd .le. 2d0*yy ) then
      WAnyBranching_flat = Chm_
  elseif(xRnd .le. 2d0*yy +     xx ) then
      WAnyBranching_flat = ElM_
  elseif(xRnd .le. 2d0*yy + 2d0*xx ) then
      WAnyBranching_flat = MuM_
  elseif(xRnd .le. 2d0*yy + 3d0*xx ) then
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





! VBranching and VVBranchings also take into account color multiplicity, so it should be used only when color is not taken into account inside the ME!
! This is why some final states count 3 times instead of just once.
! Rule of thumb for counting final states in this function is e,mu,ta,nus=1, d,u,s,c,b=3
! When a CKM partner is called, since sum(VCKM**@)=1 along a row, counting should not be unaffected.
SUBROUTINE VBranching(DecayMode,MY_IDUP,ICOLUP,CombWeight,ColorBase)
use ModParameters
use ModMisc
implicit none
integer, intent(in) :: DecayMode
integer :: MY_IDUP(1:3),ICOLUP(1:2,1:2),DKFlavor,ICOLUP_Base
integer, optional ::ColorBase
real(8), intent(out) :: CombWeight
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
   CombWeight = 1d0
   if( DecayMode.eq.0 ) then! Z1->2l
        call random_number(DKRnd)
        MY_IDUP(1) = Z0_
        DKFlavor = ZLepBranching_flat( DKRnd )!= ElM or MuM
        MY_IDUP(2) =-DKFlavor
        MY_IDUP(3) =+DKFlavor
        CombWeight = CombWeight*2d0
   elseif( DecayMode.eq.1 ) then! Z1->2q
        call random_number(DKRnd)
        MY_IDUP(1) = Z0_
        DKFlavor = ZQuaBranching_flat( DKRnd )!= Up,Dn,Chm,Str,Bot
        MY_IDUP(2) =-DKFlavor
        MY_IDUP(3) =+DKFlavor
        ICOLUP(1:2,1) = (/            0,ICOLUP_BASE+3/)
        ICOLUP(1:2,2) = (/ICOLUP_BASE+3,            0/)
        CombWeight = CombWeight*5d0*3d0
   elseif( DecayMode.eq.2 ) then! Z1->2tau
        MY_IDUP(1) = Z0_
        MY_IDUP(2) = TaP_
        MY_IDUP(3) = TaM_
   elseif( DecayMode.eq.3 ) then! Z1->2nu
        call random_number(DKRnd)
        MY_IDUP(1) = Z0_
        DKFlavor = ZNuBranching_flat( DKRnd )!= NuE,NuM,NuT
        MY_IDUP(2) =-DKFlavor
        MY_IDUP(3) =+DKFlavor
        CombWeight = CombWeight*3d0
   elseif( DecayMode.eq.4 ) then! W1(+)->lnu
        call random_number(DKRnd)
        MY_IDUP(1) = Wp_
        DKFlavor = WLepBranching_flat( DKRnd )!= ElM or MuM
        MY_IDUP(2) = +abs(DKFlavor)     ! lepton(+)
        MY_IDUP(3) = +abs(DKFlavor)+7   ! neutrino
        CombWeight = CombWeight*2d0
   elseif( DecayMode.eq.5 ) then! W1(+)->2q
        call random_number(DKRnd)
        MY_IDUP(1) = Wp_
        DKFlavor = WQuaUpBranching_flat( DKRnd )!= Up,Chm
        MY_IDUP(3) = +abs(DKFlavor)           ! up flavor
        MY_IDUP(2) = GetCKMPartner(MY_IDUP(3))! anti-dn flavor
        ICOLUP(1:2,1) = (/            0,ICOLUP_BASE+3/)
        ICOLUP(1:2,2) = (/ICOLUP_BASE+3,            0/)
        CombWeight = CombWeight*2d0*3d0
   elseif( DecayMode.eq.6 ) then! W1(+)->taunu
        MY_IDUP(1) = Wp_
        MY_IDUP(2) = TaP_
        MY_IDUP(3) = NuT_
   elseif( DecayMode.eq.7 ) then! photon
        MY_IDUP(1) = Pho_
        MY_IDUP(2) = Not_a_particle_
        MY_IDUP(3) = Not_a_particle_
   elseif( DecayMode.eq.8 ) then! Z1->2l+2tau
        call random_number(DKRnd)
        MY_IDUP(1) = Z0_
        DKFlavor = ZLepPlusTauBranching_flat( DKRnd )!= ElM or MuM or TaM
        MY_IDUP(2) =-DKFlavor
        MY_IDUP(3) =+DKFlavor
        CombWeight = CombWeight*3d0
   elseif( DecayMode.eq.9 ) then! Z1-> anything
        call random_number(DKRnd)
        MY_IDUP(1) = Z0_
        DKFlavor = ZAnyBranching_flat( DKRnd )
        MY_IDUP(2) =-DKFlavor
        MY_IDUP(3) =+DKFlavor
        if(IsAQuark(DKFlavor)) then
           ICOLUP(1:2,1) = (/            0,ICOLUP_BASE+3/)
           ICOLUP(1:2,2) = (/ICOLUP_BASE+3,            0/)
        endif
        CombWeight = CombWeight*21d0!*(6d0 + 5d0*3d0)
   elseif( DecayMode.eq.10 ) then! W1(+)->l+tau  +nu
        call random_number(DKRnd)
        MY_IDUP(1) = Wp_
        DKFlavor = WLepPlusTauBranching_flat( DKRnd )!= ElM or MuM or TaM
        MY_IDUP(2) = +abs(DKFlavor)     ! lepton(+)
        MY_IDUP(3) = +abs(DKFlavor)+7   ! neutrino
        CombWeight = CombWeight*3d0
   elseif( DecayMode.eq.11 ) then! W1(+)-> anything
        call random_number(DKRnd)
        MY_IDUP(1) = Wp_
        DKFlavor = WAnyBranching_flat( DKRnd )
        if(IsAQuark(DKFlavor)) then
           MY_IDUP(3) = +abs(DKFlavor)           ! up flavor
           MY_IDUP(2) = GetCKMPartner(MY_IDUP(3))! anti-dn flavor
           ICOLUP(1:2,1) = (/            0,ICOLUP_BASE+3/)
           ICOLUP(1:2,2) = (/ICOLUP_BASE+3,            0/)
        else
           MY_IDUP(2) = +abs(DKFlavor)     ! lepton(+)
           MY_IDUP(3) = +abs(DKFlavor)+7   ! neutrino
        endif
        CombWeight = CombWeight*9d0!*(3d0 + 2d0*3d0)
   ! Exclusive Z1 decay modes
   elseif( &
      DecayMode.eq.-2*2   .or. & ! Z->ee
      DecayMode.eq.-3*3   .or. & ! Z->mumu
      DecayMode.eq.-5*5   .or. & ! Z->dd
      DecayMode.eq.-7*7   .or. & ! Z->uu
      DecayMode.eq.-11*11 .or. & ! Z->ss
      DecayMode.eq.-13*13 .or. & ! Z->cc
      DecayMode.eq.-17*17      & ! Z->bb
      ) then
        MY_IDUP(1) = Z0_
        DKFlavor = (                      &
            ElM_*LogicalToInteger(DecayMode.eq.-2*2)   + & ! Z->ee
            MuM_*LogicalToInteger(DecayMode.eq.-3*3)   + & ! Z->mumu
            Dn_*LogicalToInteger(DecayMode.eq.-5*5)    + & ! Z->dd
            Up_*LogicalToInteger(DecayMode.eq.-7*7)    + & ! Z->uu
            Str_*LogicalToInteger(DecayMode.eq.-11*11) + & ! Z->ss
            Chm_*LogicalToInteger(DecayMode.eq.-13*13) + & ! Z->cc
            Bot_*LogicalToInteger(DecayMode.eq.-17*17)   & ! Z->bb
        )
        MY_IDUP(2) =-DKFlavor
        MY_IDUP(3) =+DKFlavor
        if(IsAQuark(DKFlavor)) then
           ICOLUP(1:2,1) = (/            0,ICOLUP_BASE+3/)
           ICOLUP(1:2,2) = (/ICOLUP_BASE+3,            0/)
        endif
   ! Exclusive W+ decay modes
   elseif( &
      DecayMode.eq.-2*1   .or. & ! W->enu
      DecayMode.eq.-3*1   .or. & ! W->munu
      DecayMode.eq.-5*7   .or. & ! W->du
      DecayMode.eq.-5*13  .or. & ! W->dc
      DecayMode.eq.-11*7  .or. & ! W->su
      DecayMode.eq.-11*13 .or. & ! W->sc
      DecayMode.eq.-17*7  .or. & ! W->bu
      DecayMode.eq.-17*13      & ! W->bc
      ) then
        MY_IDUP(1) = Wp_
        MY_IDUP(2) = (                    &
           ElP_*LogicalToInteger(DecayMode.eq.-2*1)    + & ! W->enu
           MuP_*LogicalToInteger(DecayMode.eq.-3*1)    + & ! W->munu
           ADn_*LogicalToInteger(DecayMode.eq.-5*7)    + & ! W->du
           ADn_*LogicalToInteger(DecayMode.eq.-5*13)   + & ! W->dc
           AStr_*LogicalToInteger(DecayMode.eq.-11*7)  + & ! W->su
           AStr_*LogicalToInteger(DecayMode.eq.-11*13) + & ! W->sc
           ABot_*LogicalToInteger(DecayMode.eq.-17*7)  + & ! W->bu
           ABot_*LogicalToInteger(DecayMode.eq.-17*13)   & ! W->bc
        )
        MY_IDUP(3) = (                   &
           NuE_*LogicalToInteger(DecayMode.eq.-2*1)   + & ! W->enu
           NuM_*LogicalToInteger(DecayMode.eq.-3*1)   + & ! W->munu
           Up_*LogicalToInteger(DecayMode.eq.-5*7)    + & ! W->du
           Chm_*LogicalToInteger(DecayMode.eq.-5*13)  + & ! W->dc
           Up_*LogicalToInteger(DecayMode.eq.-11*7)   + & ! W->su
           Chm_*LogicalToInteger(DecayMode.eq.-11*13) + & ! W->sc
           Up_*LogicalToInteger(DecayMode.eq.-17*7)   + & ! W->bu
           Chm_*LogicalToInteger(DecayMode.eq.-17*13)   & ! W->bc
        )
        if(IsAQuark(DKFlavor)) then
           ICOLUP(1:2,1) = (/            0,ICOLUP_BASE+3/)
           ICOLUP(1:2,2) = (/ICOLUP_BASE+3,            0/)
        endif
   endif


RETURN
END SUBROUTINE

SUBROUTINE VVBranchings(MY_IDUP,ICOLUP,CombWeight,ColorBase)
use ModParameters
use ModMisc
implicit none
integer :: MY_IDUP(4:9),ICOLUP(1:2,6:9),ICOLUP_Base
integer, optional ::ColorBase
real(8), intent(out) :: CombWeight
integer :: tmp_idup(1:3),tmp_icolup(1:2,1:2)
real(8) :: tmp_CombWeight
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
   CombWeight = 1d0

   call VBranching(DecayMode1, tmp_idup, tmp_icolup, tmp_CombWeight, ICOLUP_BASE)
   MY_IDUP(4) = tmp_idup(1)
   MY_IDUP(6) = tmp_idup(2)
   MY_IDUP(7) = tmp_idup(3)
   ICOLUP(1:2,6) = tmp_icolup(1:2,1)
   ICOLUP(1:2,7) = tmp_icolup(1:2,2)
   CombWeight = CombWeight * tmp_CombWeight

   if( Process.lt.110 .or. Process.gt.114) then
   ICOLUP_BASE = ICOLUP_BASE+1
   call VBranching(DecayMode2, tmp_idup, tmp_icolup, tmp_CombWeight, ICOLUP_BASE)
   MY_IDUP(5) = tmp_idup(1)
   MY_IDUP(8) = tmp_idup(2)
   MY_IDUP(9) = tmp_idup(3)
   ICOLUP(1:2,8) = tmp_icolup(1:2,1)
   ICOLUP(1:2,9) = tmp_icolup(1:2,2)
   if (isAWDecay(DecayMode2)) then
      MY_IDUP(5) = -MY_IDUP(5)
      call swap(MY_IDUP(8),MY_IDUP(9))
      MY_IDUP(8) = -MY_IDUP(8)
      MY_IDUP(9) = -MY_IDUP(9)
   endif
   CombWeight = CombWeight * tmp_CombWeight
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
real(8),intent(out) :: ctau
   call getSimpleDecayLength(HiggsDecayLengthMM,ctau)
END SUBROUTINE

SUBROUTINE getVprimeDecayLength(ctau)
use ModParameters
implicit none
real(8),intent(out) :: ctau
   call getSimpleDecayLength(VprimeDecayLengthMM,ctau)
END SUBROUTINE

subroutine getSimpleDecayLength(ctau0,ctau)
implicit none
real(8), intent(in) :: ctau0
real(8), intent(out) :: ctau
real(8) :: x,xp,xpp,Len0,propa
integer :: loop

     ctau  = 0d0
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
end subroutine


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




SUBROUTINE PDFMapping(MapType,yRnd,eta1,eta2,Ehat,sHatJacobi,EhatMin)
use ModParameters
use ModMisc
implicit none
integer :: MapType
real(8) :: yRnd(1:2),eta1,eta2,EHat,sHatJacobi,tau,nPotMap,z,sbar,fmax
real(8) :: etamin, Ymax, Y, Ymin, MThresh
real(8), optional :: EhatMin

  if( present(EhatMin) ) then
    MThresh = EhatMin
  else
    MThresh = M_Reso
  endif

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

     etamin = (M_Reso/Collider_Energy)**2
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

  elseif( MapType.eq.14 .or. MapType.eq.15 .or. MapType.eq.16 .or. MapType.eq.17 ) then! V-Higgs associated production, or H+j (16), or HH (17)
      if (MapType.eq.17) then
         etamin = (M_Reso+M_Reso-10d0*Ga_Reso)/Collider_Energy
      else
         etamin = (M_Reso-10d0*Ga_Reso)/Collider_Energy
         if (includeGammaStar .and. MPhotonCutoff.gt.0d0) then
            etamin = max(etamin, MPhotonCutoff/Collider_Energy)
         endif
      endif
      etamin = max(0d0, etamin)
      z = 1d0/(1d0-etamin)
      if (isnan(z) .or. z.le.yRnd(1)) then
         sHatJacobi=0d0
         eta1=0d0
         eta2=0d0
      else
         Ymin = ((z-1d0)/(z-yRnd(1)))**2
         Ymax = 1d0 / Ymin
         Ymin = 0.5d0*dlog(Ymin)
         Ymax = 0.5d0*dlog(Ymax)
         Y = Ymin + yrnd(2)*(Ymax-Ymin)
         eta1 = (z-1d0)/(z-yRnd(1))*dexp(Y)
         eta2 = (z-1d0)/(z-yRnd(1))*dexp(-Y)
         sHatJacobi = 2d0*(z-1d0)**2/(z-yRnd(1))**3*(Ymax-Ymin)
      endif

  else
      call Error("PDF mapping not available")
  endif

  EHat = Collider_Energy*dsqrt(eta1*eta2)

RETURN
END SUBROUTINE


FUNCTION GetBWPropagator(sHat, scheme)
use modMisc
use modParameters
implicit none
real(8) :: GetBWPropagator,sHat
real(8) :: mhb, ghb, BigGamma
integer :: scheme

    if( scheme.eq.1 ) then! running width
        GetBWPropagator =  1d0/( (sHat-M_Reso**2)**2 + (sHat*Ga_Reso/M_Reso)**2 )
    elseif( scheme.eq.2 ) then! fixed width
        GetBWPropagator = 1d0/( (sHat-M_Reso**2)**2 + (M_Reso*Ga_Reso)**2 )
    elseif( scheme.eq.3 ) then! Passarino's CPS
        if( mubarH.lt.0d0 .or. gabarH.lt.0d0 ) then
          call CALL_HTO(M_Reso/GeV, m_top/GeV, mhb, ghb)
          if( IsNaN(mubarH).or.IsNaN(gabarH) ) then
            print *, "Passarino's CALL_HTO returned a NaN"
            print *, "(gabarH,Ehat)",gabarH,dsqrt(dabs(sHat))/GeV
            stop 1
            RETURN
          endif
          mhb = mhb*GeV
          ghb = ghb*GeV

          mubarH = sqrt(mhb**2/(1d0+(ghb/mhb)**2))
          gabarH = mubarH/mhb*ghb
        endif

        GetBWPropagator = 1d0/( (sHat-mubarH**2)**2 + (mubarH*gabarH)**2 )

        !call HTO_gridHt(dsqrt(dabs(sHat))/GeV,BigGamma)
        !BigGamma = BigGamma*GeV

        !print *, dsqrt(dabs(sHat))/GeV, gabarH/GeV, BigGamma/GeV
    elseif( scheme.eq.0 ) then  !remove the propagator completely
        GetBWPropagator = 1d0
    else
        print *, "Invalid scheme: ", scheme
        stop 1
    endif


RETURN
END FUNCTION

FUNCTION ReweightBWPropagator(sHat)! sHat is the resonance inv. mass squared
use modMisc
use modParameters
implicit none
real(8) :: ReweightBWPropagator,sHat
real(8) :: BreitWigner,BreitWigner_Run


     ReweightBWPropagator = 1d0
     BreitWigner = GetBWPropagator(sHat, 2)
     BreitWigner_Run = GetBWPropagator(sHat, WidthScheme)

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

! print *, "smeared mt",(BW_Mass(2:3)-m_top)/GeV
! print *, "smeared mw",(BW_Mass(4:5)-m_w)/GeV

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


SUBROUTINE TW_OffShellProjection(MomIn,MomOut,Jacobian)
use modParameters
use modMisc
implicit none
real(8) :: MomIn(:,:),MomOut(:,:),MomTmp(1:4),Jacobian
real(8) :: xRndWidth(2:5),BW_Mass(2:5),BW_Jacobi(2:5)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4,Wm=5,  b=6,Wp=7,lepP=8,nu=9, lepM=10,nubar=11


    call random_number(xRndWidth)

    call SmearExternal(xRndWidth(3),m_top,Ga_Top,m_top-6d0*Ga_Top,m_top+6d0*Ga_Top,BW_Mass(3),BW_Jacobi(3)) !top
    call SmearExternal(xRndWidth(4),m_W,Ga_W,m_W-6d0*Ga_W,m_W+6d0*Ga_W,BW_Mass(4),BW_Jacobi(4)) !wp from top-decay
    call SmearExternal(xRndWidth(5),m_W,Ga_W,m_W-6d0*Ga_W,m_W+6d0*Ga_W,BW_Mass(5),BW_Jacobi(5)) !wm
    Jacobian = BW_Jacobi(3) * BW_Jacobi(4) * BW_Jacobi(5)

! print *, "smeared mt",(BW_Mass(3)-m_top)/GeV
! print *, "smeared mw",(BW_Mass(4:5)-m_w)/GeV

    call ShiftMass(MomIn(1:4,t),MomIn(1:4,Wm),BW_Mass(3),BW_Mass(5),MomOut(1:4,t),MomOut(1:4,Wm))

    MomTmp(1:4) = MomOut(1:4,t) - MomIn(1:4,Wp)
    call ShiftMass(MomTmp,MomIn(1:4,Wp),m_Bot,BW_Mass(4),MomOut(1:4,b),MomOut(1:4,Wp))! project b and W+ on-shell


    MomTmp(1:4) = MomOut(1:4,Wp) - MomIn(1:4,lepP)                                ! project W+ decay products on-shell
    MomOut(1:4,nu)   = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP)) * MomIn(1:4,lepP)
    MomOut(1:4,lepP) = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP))) * MomIn(1:4,lepP)

    MomTmp(1:4) = MomOut(1:4,Wm) - MomIn(1:4,lepM)                                ! project W- decay products on-shell
    MomOut(1:4,nubar) = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM)) * MomIn(1:4,lepM)
    MomOut(1:4,lepM)  = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM))) * MomIn(1:4,lepM)

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




SUBROUTINE EvalPhaseSpace_VHiggs(yRnd,MomExt,inv_mass,mass,PSWgt,useAonshell)
use ModParameters
implicit none
      logical, intent(in) :: useAonshell
      !logical, intent(in) :: HDecays
      !logical, intent(in), optional :: PhoOnshell
      !logical, optional :: ZAinterference
      real(8), intent(in) :: yRnd(1:20),mass(9,2)
      real(8) :: phi, beta, gamma
      real(8) :: temp_vector(4), temp_boost(4)
      integer :: i
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
      logical :: hasInterference
      !if(present(PhoOnshell)) then
      !   hasAonshell=PhoOnshell
      !endif
      !hasInterference = .false.
      !if(present(ZAinterference)) then
      !   hasInterference=ZAinterference
      !endif

      MomExt(:,4:9)=0d0
      inv_mass(4:9)=0d0

!333333333
      inv_mass(3) = dsqrt((MomExt(1,3)+MomExt(4,3))*(MomExt(1,3)-MomExt(4,3)))

!555555555
      if(H_DK)then
        inv_mass(5) = dsqrt(dabs(bw_sq(yRnd(13),mass(5,1), mass(5,2), inv_mass(3)**2, jacobian5)))
        jacobian5 = jacobian5 /16d0 /Pi**2 != (ds5/2pi)*(1/8pi)
      else
        inv_mass(5) = mass(5,1)
        jacobian5=1d0
      endif

!4444444444
      if(useAonshell) then
        inv_mass(4) = mass(4,1)
!print*,inv_mass(4),mass(4,1)
        jacobian4=1d0
      elseif(.not.hasInterference)then
        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
        jacobian4 = jacobian4 /16d0 /Pi**2 != (ds4/2pi)*(1/8pi)
      else
        inv_mass(4) = dsqrt(dabs(bw_sq(yRnd(12),mass(4,1), mass(4,2), (inv_mass(3)-inv_mass(5))**2, jacobian4)))
        !inv_mass(4) = dsqrt(one_over_s_sq(yRnd(12), MPhotonCutoff**2, (inv_mass(3)-inv_mass(5))**2, jacobian4))
        jacobian4 = jacobian4 /16d0 /Pi**2
      endif

!444444444444
!energy of 4 in the CM frame of 3
      MomExt(1,4)=(inv_mass(3)**2+(inv_mass(4)+inv_mass(5))*(inv_mass(4)-inv_mass(5)))/2d0/inv_mass(3)

      if(MomExt(1,4).lt.0d0)then
        MomExt=0d0
        PSWgt=0d0
      endif

!|3-momentum| of 4 in the CM frame of 3
      cm_abs3p(4) = dsqrt((MomExt(1,4)+inv_mass(4)) * (MomExt(1,4)-inv_mass(4)))
!generating cos(theta_4) and phi_4 in the CM frame of 3
      cm_cos_theta(4) = yRnd(7)
      cm_cos_theta(4) = cm_cos_theta(4)*2d0-1d0
      cm_sin_theta(4) = dsqrt((1d0+cm_cos_theta(4))  *(1d0-cm_cos_theta(4)))
      phi = yRnd(6)
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
      if(.not.useAonshell) then
!invariant mass of 6
         inv_mass(6)=0d0
!energy of 6 in the CM frame of 4
         MomExt(1,6)=inv_mass(4)/2d0
!|3-momentum| of 6 in the CM frame of 4
         cm_abs3p(6)=MomExt(1,6)
!generating cos(theta_6) and phi_6 in the CM frame of 4
!z-axis is along the boost of 2
         cm_cos_theta(6) = yRnd(9)
         cm_cos_theta(6) = cm_cos_theta(6)*2d0-1d0
         cm_sin_theta(6) = dsqrt((1d0+cm_cos_theta(6)) *(1d0-cm_cos_theta(6)))
         cm_cos_phi(6) = yRnd(8)
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
      else
         MomExt(1:4,6)=MomExt(1:4,4)
      endif

!7777777777777777777777
!invariant mass of 7
      inv_mass(7)=0d0
!4-momentum of 7 (lab frame) by energy-momentum conservation
      MomExt(1:4,7)=MomExt(1:4,4)-MomExt(1:4,6)

!8888888888888888888888
      if(H_DK)then
!invariant mass of 8
        inv_mass(8)=0d0
!energy of 8 in the CM frame of 5
        MomExt(1,8)=inv_mass(5)/2d0
!|3-momentum| of 8 in the CM frame of 5
        cm_abs3p(8)=MomExt(1,8)
!generating cos(theta_8) and phi_8 in the CM frame of 5
!z-axis is along the boost of 5
        cm_cos_theta(8) = yRnd(11)
        cm_cos_theta(8) = cm_cos_theta(8)*2d0-1d0
        cm_sin_theta(8) = dsqrt((1d0+cm_cos_theta(8)) *(1d0-cm_cos_theta(8)))
        cm_cos_phi(8) = yRnd(10)
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
!4-momentum of 9 (lab frame) by energy-momentum conservation
        MomExt(1:4,9)=MomExt(1:4,5)-MomExt(1:4,8)
      endif

      PSWgt = jacobian4*jacobian5 * cm_abs3p(4)/(4d0*pi)/inv_mass(3)
      !if(isnan(PSWgt).or.(PSWgt.eq.0d0)) print *,  "()",inv_mass, jacobian4, jacobian5, PSWgt, "()"

!do i=4,7
!print *, dsqrt(dabs(four_momentum(i,:).dot.four_momentum(i,:)))
!enddo
!pause

      if(isnan(PSWgt).or.PSWgt.eq.0d0)then
        MomExt=0d0
        PSWgt=0d0
      endif


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
!  print *, "Ehat",Ehat/GeV
!  print *, "E1",e1/GeV
!  print *, "E2",e2/GeV
!  print *, "-sqrt(s1)",-dsqrt(dabs(s1/GeV**2))
!  print *, "-sqrt(s2)",-dsqrt(dabs(s2/GeV**2))
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

SUBROUTINE EvalPhasespace_VBF_or_HJJ(xchannel,xRnd,Energy,Mom,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xchannel,xRnd(:), Energy, Mom(:,:)
real(8) :: Jac,Jac1,Jac2,Jac3,Mom_ij_Dummy(1:4),s35,s45,minmax(1:2)
real(8) :: tMassSq

!tMassSq=(M_Reso*0.5d0)**2
tMassSq=(M_Reso)**2

Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

if ( xchannel.lt.0.25 ) then

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                           !  int d(s35)
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),tMassSq,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,4))              !  1+2 --> (35)+4
   Jac3 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2)-Mom(1:4,4),tMassSq,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,3),Mom(1:4,5))    !  1+(24) --> 3+5
   Jac = Jac1*Jac2*Jac3 * PSNorm3                                                                                        !  combine

elseif ( xchannel.lt.0.5 ) then

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s35)                                                                           !  int d(s35)
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),tMassSq,s35,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,4))              !  1+2 --> (35)+4
   Jac3 = t_channel_prop_decay(Mom(1:4,1)-Mom(1:4,4),Mom(1:4,2),tMassSq,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,3),Mom(1:4,5))    !  (14)+2 --> 3+5
   Jac  = Jac1*Jac2*Jac3 * PSNorm3                                                                                       !  combine

elseif ( xchannel.lt.0.75 ) then

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                           !  int d(s45)
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),tMassSq,s45,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,3))              !  1+2 --> (45)+3
   Jac3 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2)-Mom(1:4,3),tMassSq,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,4),Mom(1:4,5))    !  1+(23) --> 4+5
   Jac  = Jac1*Jac2*Jac3 * PSNorm3                                                                                       !  combine

else

   Jac1 = k_l(xRnd(1),M_Reso**2,Energy**2,s45)                                                                           !  int d(s45)
   Jac2 = t_channel_prop_decay(Mom(1:4,1),Mom(1:4,2),tMassSq,s45,0d0,xRnd(2:3),Mom_ij_Dummy(1:4),Mom(1:4,3))              !  1+2 --> (45)+3
   Jac3 = t_channel_prop_decay(Mom(1:4,1)-Mom(1:4,3),Mom(1:4,2),tMassSq,0d0,M_Reso**2,xRnd(4:5),Mom(1:4,4),Mom(1:4,5))    !  (13)+2 --> 4+5
   Jac  = Jac1*Jac2*Jac3 * PSNorm3                                                                                       !  combine

endif

if( isnan(Jac) ) then! THIS SHOULD BE REMOVED AFTER DEBUGGING
   Jac = 0d0
   print *, "ERROR in EvalPhasespace_VBF_or_HJJ: NaN Jac",Energy,xchannel,xRnd
endif

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

SUBROUTINE EvalPhasespace_VBF_H4f(xchannel,xRnd,Energy,Mom,Jac,ids,swap34_56,do78,id12_78)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xchannel,xRnd(:), Energy, Mom(:,:)
integer :: EqualLeptons,ids(1:8)
integer :: iChannel
real(8) :: Jac,Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9
real(8) :: s3H,s4H,s56,s78,s910,s34,s35,s46,Mom_Dummy(1:4),Mom_Dummy2(1:4),xRndLeptInterf,Emin,Emax
real(8) :: BWmass_ps, BWwidth_ps
real(8) :: s1min, s2min
integer :: NumChannels, it_chan, ch_ctr
integer :: id12, id78, id17, id28, id18, id27, id12_78
logical :: swap34_56, do78, isZH, isWH, isVBF
integer,parameter :: inTop=1, inBot=2, outTop=3, outBot=4, V1=5, V2=6, Lep1P=7, Lep1M=8, Lep2P=9, Lep2M=10
logical,parameter :: includeNewBWPSinEW = .true.

   isZH = .false.
   isWH = .false.
   isVBF = .false.
   swap34_56 = .false.
   do78 = .false.
   id12_78 = Not_a_particle_
   BWmass_ps=-1d0
   BWwidth_ps=-1d0

   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

   EMin = min(Energy,m4l_minmax(1))
   EMax = min(Energy,m4l_minmax(2))

   s1min = (max(MPhotonCutoff,0d0))**2
   s2min = (max(MPhotonCutoff,0d0))**2

   if( EMin.lt.0d0 .or. EMin.gt.EMax ) call Error("m4l_minmax is not set correctly")


   if ( Process .ne. 69) then
      id12=CoupledVertex((/-ids(1),-ids(2)/),-1)
      id78=CoupledVertex(ids(7:8),-1)
      id17=CoupledVertex((/-ids(1),ids(7)/),-1)
      id28=CoupledVertex((/-ids(2),ids(8)/),-1)
      id18=CoupledVertex((/-ids(1),ids(8)/),-1)
      id27=CoupledVertex((/-ids(2),ids(7)/),-1)

      isZH = (id12.eq.Z0_ .and. id78.eq.id12)
      isWH = (abs(id12).eq.abs(Wp_) .and. id78.eq.-id12 .and. CoupledVertexIsDiagonal(ids(1:2),-1)) ! Must require W from diagonal elements of CKM as in the ME
      isVBF = ( &
         ( (id17.eq.Z0_ .and. id28.eq.id17) ) .or. &
         ( (abs(id17).eq.abs(Wp_) .and. id28.eq.-id17) .and. CoupledVertexIsDiagonal((/-ids(1),ids(7)/),-1) ) .or. &
         ( (id18.eq.Z0_ .and. id27.eq.id18) ) .or. &
         ( (abs(id18).eq.abs(Wp_) .and. id27.eq.-id18) .and. CoupledVertexIsDiagonal((/-ids(1),ids(8)/),-1) ) &
         )

      if (((isWH .and. .not.isZH) .or. (isZH .and. .not.isWH)) .and. .not.isVBF) id12_78 = id78

      NumChannels = 0
      if (.not.(isZH .or. isWH .or. isVBF)) call Error("EW off-shell process is not ZVV, WVV or VBF/VBS")
      if (isZH) NumChannels = NumChannels+1
      if (isWH) NumChannels = NumChannels+1
      if (isVBF) NumChannels = NumChannels+1
      iChannel = int(xchannel * NumChannels -1d-10)+1

      !write(6,*) "ids=",ids
      !write(6,*) "Various ids=",id12,id78,id17,id28,id18,id27
      !write(6,*) "isZH?",isZH
      !write(6,*) "isWH?",isWH
      !write(6,*) "isVBF?",isVBF
      !write(6,*) "iChannel = ",iChannel,"/",NumChannels,"(xchannel: ",xchannel,")"
      !pause

      if (includeNewBWPSinEW) then
         if (Process.eq.67 .or. Process.eq.71) then
            ! This is almost flat but not quite
            BWmass_ps = M_Z
            BWwidth_ps = 4d0*M_V_ps-BWmass_ps ! Since 2*M_V needs to be covered
            if (Emin.lt.(BWmass_ps-10d0*BWwidth_ps) .or. Emax.gt.(BWmass_ps+10d0*BWwidth_ps)) then
               BWmass_ps = -1d0
               BWwidth_ps = -1d0
            endif
         else if (Process.eq.66 .or. Process.eq.70) then
            if( M_Reso.ge.0d0 .and. M_Reso2.ge.0d0 ) then ! Both resonances are present
               BWmass_ps = max(M_Reso,M_Reso2)
               BWwidth_ps = max(max(abs(M_Reso-M_Reso2),Ga_Reso), Ga_Reso2) ! Cover the full mass difference
            else if( M_Reso.ge.0d0 ) then
               BWmass_ps = M_Reso
               BWwidth_ps = Ga_Reso
            else if( M_Reso2.ge.0d0 ) then
               BWmass_ps = M_Reso2
               BWwidth_ps = Ga_Reso2
            endif
            if (Emin.lt.(BWmass_ps-10d0*BWwidth_ps) .or. Emax.gt.(BWmass_ps+10d0*BWwidth_ps)) then
               BWmass_ps = -1d0
               BWwidth_ps = -1d0
            endif
         else if (Process.eq.68 .or. Process.eq.72) then
            if( M_Reso.ge.0d0 .and. M_Reso2.ge.0d0 ) then ! Both resonances are present
               BWmass_ps = max(M_Reso,M_Reso2)
               BWwidth_ps = max(max(abs(M_Reso-M_Reso2),Ga_Reso), Ga_Reso2) ! Cover the full mass difference
            else if( M_Reso.ge.0d0 ) then
               BWmass_ps = M_Reso
               BWwidth_ps = Ga_Reso
            else if( M_Reso2.ge.0d0 ) then
               BWmass_ps = M_Reso2
               BWwidth_ps = Ga_Reso2
            endif
            if ( &
               Emin.lt.(BWmass_ps-10d0*BWwidth_ps) .or. Emax.gt.(BWmass_ps+10d0*BWwidth_ps) &
               .or. (Emin.lt.M_Z .and. abs(BWmass_ps-M_Z).gt.10d0*BWwidth_ps) & ! Check distance from M_Z as well for BSI
               ) then
               BWmass_ps = -1d0
               BWwidth_ps = -1d0
            endif
         endif
      else
         if( .not.(Emin.gt.M_Reso+2d0*Ga_Reso .or. Emax.lt.M_Reso-2d0*Ga_Reso .or. Process.eq.67 .or. Process.eq.71) ) then
            BWmass_ps = M_Reso
            BWwidth_ps = Ga_Reso
         endif
      endif
      if( BWmass_ps.lt.0d0 .or. BWwidth_ps.lt.0d0 ) then ! Create flat 4f mass
         Jac1 = k_l(xRnd(1),Emin**2,Emax**2,s56)
      else ! Create a BW 4f mass
         Jac1 = k_BreitWigner(xRnd(1),BWmass_ps**2,BWwidth_ps,Emin**2,Emax**2,s56)
      endif
      Jac3 = s_channel_propagator(M_V_ps**2,Ga_V_ps,s1min,s56,xRnd(3),s78) ! m1
      Jac4 = s_channel_propagator(M_V_ps**2,Ga_V_ps,s2min,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910) ! m2
      !if (s2min .ge.(dsqrt(s56)-dsqrt(s78))**2) then
      !   write(6,*) "Invalid s910 generated: s56, s78, s910, s2min, comp = ",s56,s78,s910,s2min,(dsqrt(s56)-dsqrt(s78))**2
      !endif

      ch_ctr=1
      do it_chan=1,4
         if (it_chan .eq. 1) then ! ZH PS
            if (.not. isZH) cycle
            if (iChannel .ne. ch_ctr) then
               ch_ctr = ch_ctr+1
               cycle
            endif

            do78 = .true.
            Jac2 = s_channel_propagator(M_Z**2,Ga_Z*5,mJJcut**2,(Energy-dsqrt(s56))**2,xRnd(2),s34) ! Associated mZ for Z->ff
            Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:))
            Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4))
            Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))

            !write(6,*) "Did ZH PS"

            exit
         elseif (it_chan .eq. 2) then ! WH PS
            if (.not. isWH) cycle
            if (iChannel .ne. ch_ctr) then
               ch_ctr = ch_ctr+1
               cycle
            endif

            do78 = .true.
            Jac2 = s_channel_propagator(M_W**2,Ga_W*5,mJJcut**2,(Energy-dsqrt(s56))**2,xRnd(2),s34) ! Associated mW for W->ff
            Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:))
            Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4))
            Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))

            !write(6,*) "Did WH PS"

            exit
         elseif (it_chan .eq. 3) then ! VBF unswapped PS
            if (.not. isVBF) cycle
            if (iChannel .ne. ch_ctr) then
               ch_ctr = ch_ctr+1
               cycle
            endif

            do78 = .true.
            Jac2 = k_l(xRnd(2),mJJcut**2,(Energy-dsqrt(s56))**2,s34)
            Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:))
            Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4))
            Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))

            !write(6,*) "Did VBF PS"

            exit
         endif

      enddo
   else
      NumChannels=1

      BWmass_ps = M_Z
      BWwidth_ps = 4d0*M_V_ps-BWmass_ps ! Since 2*M_V needs to be covered
      if (Emin.lt.(BWmass_ps-10d0*BWwidth_ps) .or. Emax.gt.(BWmass_ps+10d0*BWwidth_ps)) then
         BWmass_ps = -1d0
         BWwidth_ps = -1d0
      endif
      if( &
         BWmass_ps.lt.0d0 .or. BWwidth_ps.lt.0d0 .or. &
         Emin.gt.BWmass_ps+10d0*BWwidth_ps .or. Emax.lt.BWmass_ps-10d0*BWwidth_ps &
         ) then ! Create flat 4f mass
         Jac1 = k_l(xRnd(1),Emin**2,Emax**2,s56)
      else ! Create a BW 4f mass
         Jac1 = k_BreitWigner(xRnd(1),BWmass_ps**2,BWwidth_ps,Emin**2,Emax**2,s56)
      endif

      Jac3 = s_channel_propagator(M_V_ps**2,Ga_V_ps,0d0,s56,xRnd(3),s78) ! m1
      Jac4 = s_channel_propagator(M_V_ps**2,Ga_V_ps,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910) ! m2

      do78 = .true.
      Jac2 = k_l(xRnd(2),mJJcut**2,(Energy-dsqrt(s56))**2,s34)
      Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:))
      Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4))
      Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))
   endif

   if( includeInterference .and. ids(3).eq.ids(5) .and. ids(4).eq.ids(6) ) then
      call random_number(xRndLeptInterf)
      if( xRndLeptInterf.gt.0.5d0 ) then ! Swapped config.
         Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,9),Mom(:,8))       !   Z --> ffbar
         Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,7),Mom(:,10))      !   Z --> ffbar
         swap34_56 = .true.
      else ! Normal config
         Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))       !   Z --> ffbar
         Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))      !   Z --> ffbar
      endif
   else ! Normal config
      Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))          !   Z --> ffbar
      Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))         !   Z --> ffbar
   endif



! IF( iChannel.EQ.1 ) THEN! 34 + WW-->H-->ZZ
!
! !  masses
!    if( Emin.lt.0d0 ) then
!       Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                   !  int d(s56)    = Higgs resonance
!    else
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)                                            !  int d(s56)    = linear mapping
!    endif
!    Jac2 = k_l(xRnd(2),s56,Energy**2,s3H)                                                                                          !  int d(s3H)
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2
!
! !  splittings
!    Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_W**2,s3H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,4))                                !  1+2 --> (3H)+4
!    Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,4),M_W**2,0d0,s56,xRnd(7:8),Mom(:,3),Mom_Dummy(1:4))                       !  1+(24) --> 3+H
!    Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6
!    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8
!    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10
!
!
!
! ELSEIF( iChannel.EQ.2 ) THEN! 43 + WW-->H-->ZZ
!
! !  masses
!    if( Emin.lt.0d0 ) then
!       Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
!    else
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)                                            !  int d(s56)    = linear mapping
!    endif
!    Jac2 = k_l(xRnd(2),s56,Energy**2,s4H)                                                                                          !  int d(s4H)
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2
!
! !  splittings
!    Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_W**2,s4H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,3))                                !  1+2 --> (4H)+3
!    Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,3),M_W**2,0d0,s56,xRnd(7:8),Mom(:,4),Mom_Dummy(1:4))                       !  1+(23) --> 4+H
!    Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6
!    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8
!    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10
!
!
!
! ELSEIF( iChannel.EQ.3 ) THEN! 34 + ZZ-->H-->ZZ
!
! !  masses
!    if( Emin.lt.0d0 ) then
!       Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
!    else
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)                                            !  int d(s56)    = linear mapping
!    endif
!    Jac2 = k_l(xRnd(2),s56,Energy**2,s3H)                                                                                          !  int d(s3H)
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2
!
! !  splittings
!    Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_Z**2,s3H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,4))                                !  1+2 --> (3H)+4
!    Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,4),M_Z**2,0d0,s56,xRnd(7:8),Mom(:,3),Mom_Dummy(1:4))                       !  1+(24) --> 3+H
!    Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6
!    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8
!    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10
!
!
! ELSEIF( iChannel.EQ.4 ) THEN! 43 + ZZ-->H-->ZZ
!
! !  masses
!    if( Emin.lt.0d0 ) then
!       Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)                                                    !  int d(s56)    = Higgs resonance
!    else
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)                                            !  int d(s56)    = linear mapping
!    endif
!    Jac2 = k_l(xRnd(2),s56,Energy**2,s4H)                                                                                          !  int d(s4H)
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)                                                                   !  int d(s78)    = Z1
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)                                           !  int d(s910) = Z2
!
! !  splittings
!    Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_Z**2,s4H,0d0,xRnd(5:6),Mom_Dummy(1:4),Mom(:,3))                                !  1+2 --> (4H)+3
!    Jac6 = t_channel_prop_decay(Mom(:,1),Mom(:,2)-Mom(:,3),M_Z**2,0d0,s56,xRnd(7:8),Mom(:,4),Mom_Dummy(1:4))                       !  1+(23) --> 4+H
!    Jac7 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))                                                   !  H --> 5+6
!    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))                                                         !  5 --> 7+8
!    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))                                                        !  6 --> 9+10
!
!
!
! ! ELSEIF( iChannel.EQ.5 ) THEN
! !
! !
! !
! ! !  Z boson emissions from either of the two quark lines with t-channel Z
! !
! !    Jac1 = s_channel_propagator(M_Z**2,Ga_Z,0d0,Energy**2,xRnd(1),s78)
! !    Jac2 = s_channel_propagator(M_Z**2,Ga_Z,0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)
! !
! !    Jac3 = s_channel_propagator(0d0,0d0,s78,(Energy-dsqrt(s910))**2,xRnd(3),s35)
! !    Jac4 = s_channel_propagator(0d0,0d0,s910,(Energy-dsqrt(s35))**2,xRnd(4),s46)
! !    Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_Z**2,s35,s46,xRnd(5:6),Mom_Dummy(1:4),Mom_Dummy2(1:4))
! !
! !    Jac6 = s_channel_decay(Mom_Dummy(1:4),0d0,s78,xRnd(7:8),Mom(:,3),Mom(:,5))
! !    Jac7 = s_channel_decay(Mom_Dummy2(1:4),0d0,s910,xRnd(9:10),Mom(:,4),Mom(:,6))
! !
! !    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))
! !    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))
! !
! !
! ! ! -------------------------------
! !
! !
! ! ELSEIF( iChannel.EQ.6 ) THEN
! !
! ! !  W boson emissions from either of the two quark lines with t-channel W   (??? why??)
! !
! !    Jac1 = s_channel_propagator(M_W**2,Ga_W,0d0,Energy**2,xRnd(1),s78)
! !    Jac2 = s_channel_propagator(M_W**2,Ga_W,0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)
! !
! !    Jac3 = s_channel_propagator(0d0,0d0,s78,(Energy-dsqrt(s910))**2,xRnd(3),s35)
! !    Jac4 = s_channel_propagator(0d0,0d0,s910,(Energy-dsqrt(s35))**2,xRnd(4),s46)
! !    Jac5 = t_channel_prop_decay(Mom(:,1),Mom(:,2),M_W**2,s35,s46,xRnd(5:6),Mom_Dummy(1:4),Mom_Dummy2(1:4))
! !
! !    Jac6 = s_channel_decay(Mom_Dummy(1:4),0d0,s78,xRnd(7:8),Mom(:,3),Mom(:,5))
! !    Jac7 = s_channel_decay(Mom_Dummy2(1:4),0d0,s910,xRnd(9:10),Mom(:,4),Mom(:,6))
! !
! !    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))
! !    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))
! !
!
!
!
!
!
! ! -------------------------------
!
!
!
! ELSEIF( iChannel.EQ.5 ) THEN
!
!
! !  ZH phase space with flat Z-->jj propagator
!
!
! !  masses
!    if( Emin.gt.M_Reso+2*Ga_Reso .or. Emax.lt.M_Reso-2*Ga_Reso .or. Process.eq.67 ) then
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)
!    else
! !       Jac1 = k_BreitWigner_Quadr(xRnd(1),M_Reso**2,Ga_Reso,Emin**2,Emax**2,s56)
!      Jac1 = k_BreitWigner(xRnd(1),M_Reso**2,Ga_Reso,Emin**2,Emax**2,s56)
!    endif
!
!
! !    Jac2 = s_channel_propagator(M_Z**2,Ga_Z,0d0,(Energy-dsqrt(s56))**2,xRnd(2),s34)
!    Jac2 = k_l(xRnd(2),mJJcut**2,(Energy-dsqrt(s56))**2,s34)    ! s34 = mjj^2, hence the minumum is set to mJJcut
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)
!
! !  splittings
!    Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:)) !   Z* --> Z+H
!    Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4)) !  Z --> 34
!    Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))  !   H --> ZZ
!
!    if( includeInterference .and. EqualLeptons.eq.0 ) then!   EqualLeptons=0 means equal leptons
!       call random_number(xRndLeptInterf)
!       if( xRndLeptInterf.gt.0.5d0 ) then
!           Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,9),Mom(:,8))       !   Z --> ffbar
!           Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,7),Mom(:,10))      !   Z --> ffbar
!       else
!           Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))       !   Z --> ffbar
!           Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))      !   Z --> ffbar
!       endif
! !       Jac8 = Jac8 * 2d0
!    else
!       Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))       !   Z --> ffbar
!       Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))      !   Z --> ffbar
!    endif
!
!
! ELSEIF( iChannel.EQ.6 ) THEN
!
!
! !  ZH phase space with BW Z-->jj propagator
!
!
! !  masses
!    if( Emin.lt.0d0 ) then
!       Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)
!    else
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)
!    endif
!    Jac2 = s_channel_propagator(M_Z**2,Ga_Z,mJJcut**2,(Energy-dsqrt(s56))**2,xRnd(2),s34) ! s34 = mjj^2, hence the minumum is set to mJJcut
! !    Jac2 = k_l(xRnd(2),0d0,(Energy-dsqrt(s56))**2,s34)
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)
!
!
! !  splittings
!    Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:)) !   Z* --> Z+H
!    Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4)) !  Z --> 34
!    Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))  !   H --> ZZ
!    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))       !   Z --> ffbar
!    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))      !   Z --> ffbar
!
!
!
!
!
! ELSEIF( iChannel.EQ.7 ) THEN
!
!
! !  WH phase space with BW W-->jj propagator
!
!
! !  masses
!    if( Emin.lt.0d0 ) then
!       Jac1 = s_channel_propagator(M_Reso**2,RescaleWidth*Ga_Reso,0d0,Energy**2,xRnd(1),s56)
!    else
!       Jac1 = k_l(xRnd(1),Emin**2,min(Energy**2,Emax**2),s56)
!    endif
!    Jac2 = s_channel_propagator(M_W**2,Ga_W,mJJcut**2,(Energy-dsqrt(s56))**2,xRnd(2),s34) ! s34 = mjj^2, hence the minumum is set to mJJcut
! !    Jac2 = k_l(xRnd(2),0d0,(Energy-dsqrt(s56))**2,s34)
!    Jac3 = s_channel_propagator(M_V**2,Ga_V,0d0,s56,xRnd(3),s78)
!    Jac4 = s_channel_propagator(M_V**2,Ga_V,0d0,(dsqrt(s56)-dsqrt(s78))**2,xRnd(4),s910)
!
!
! !  splittings
!    Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s34,s56,xRnd(5:6),Mom_Dummy(:),Mom_Dummy2(:)) !   Z* --> Z+H
!    Jac6 = s_channel_decay(Mom_Dummy(:),0d0,0d0,xRnd(7:8),Mom(:,3),Mom(:,4)) !  Z --> 34
!    Jac7 = s_channel_decay(Mom_Dummy2(:),s78,s910,xRnd(9:10),Mom(:,5),Mom(:,6))  !   H --> ZZ
!    Jac8 = s_channel_decay(Mom(:,5),0d0,0d0,xRnd(11:12),Mom(:,7),Mom(:,8))       !   Z --> ffbar
!    Jac9 = s_channel_decay(Mom(:,6),0d0,0d0,xRnd(13:14),Mom(:,9),Mom(:,10))      !   Z --> ffbar
!
!
!
!
! ELSEIF( iChannel.EQ.999 ) THEN
!
!    call genps(6,Energy,xRnd(1:14),(/0d0,0d0,0d0,0d0,0d0,0d0/),Mom(1:4,3:8),Jac)
!    Mom(1:4,10)= Mom(1:4,8)
!    Mom(1:4,9) = Mom(1:4,7)
!    Mom(1:4,8) = Mom(1:4,6)
!    Mom(1:4,7) = Mom(1:4,5)
!    Mom(1:4,5) = Mom(1:4,7)+Mom(1:4,8)
!    Mom(1:4,6) = Mom(1:4,9)+Mom(1:4,10)
!    Jac = Jac * (2d0*Pi)**(4-(6)*3) * (4d0*Pi)**((6)-1)
!    RETURN
!
! ELSE
!    call Error("PS channel not available in EvalPhasespace_VBF_H4f",ichannel)
! ENDIF


   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6*Jac7*Jac8*Jac9 * PSNorm6

!    Jac = Jac * NumChannels

   if( isNan(jac) ) then
      print *, "EvalPhasespace_VBF_H4f NaN"
      print *, Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9,ichannel
      if( isNan(jac) ) Jac = 0d0

      write(6,*) "ids=",ids
      write(6,*) "Various ids=",id12,id78,id17,id28,id18,id27
      write(6,*) "isZH?",isZH
      write(6,*) "isWH?",isWH
      write(6,*) "isVBF?",isVBF
      write(6,*) "iChannel = ",iChannel,"/",NumChannels,"(xchannel: ",xchannel,")"
      write(6,*) "m56:",sqrt(s56)/GeV
      write(6,*) "m78:",sqrt(s78)/GeV
      write(6,*) "m910:",sqrt(s910)/GeV
      write(6,*) "m34:",sqrt(s34)/GeV
      write(6,*) "Last case:",ch_ctr

      call Exit(1)
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
!    print *, "Inv.mass  ",get_MInv(Mom_Dummy(1:4))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,5))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,6))/GeV
!    pause



RETURN
END SUBROUTINE

SUBROUTINE EvalPhasespace_gg4f(xRnd,eta1,eta2,Energy,Mom,Jac,ids,swap34_56)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
integer,parameter :: inTop=1, inBot=2, V1=3, V2=4, Lep1P=5, Lep1M=6, Lep2P=7, Lep2M=8, NUP=8
real(8) :: xRnd(:), eta1, eta2, Energy, Mom(:,:)
integer :: ids(1:6)
real(8) :: Jac,Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9
real(8) :: sysY,s34,s56,s78,xRndLeptInterf,Emin,Emax
real(8) :: BWmass_ps, BWwidth_ps
real(8) :: s1min, s2min
logical :: swap34_56

   swap34_56 = .false.

   EMin = 0d0
   EMax = Collider_Energy
   Jac = 0d0
   if (m4l_minmax(1).ge.0d0) then
      EMin = m4l_minmax(1)
   endif
   if (m4l_minmax(2).ge.0d0) then
      EMax = min(m4l_minmax(2),Emax)
   endif
   if( EMin.gt.EMax ) call Error("m4l_minmax is not set correctly")

   if( &
      includeInterference .and. &
      (CoupledVertex((/ids(3),ids(6)/),-1).ne.Not_a_particle_ .and. CoupledVertex((/ids(5),ids(4)/),-1).ne.Not_a_particle_) &
      ) then
      call random_number(xRndLeptInterf)
      if( xRndLeptInterf.gt.0.5d0 ) then ! Swap 46
         swap34_56 = .true.
      endif
   endif
   if (.not. swap34_56) then
      if (CoupledVertex(ids(3:4),-1,useAHcoupl=1).eq.Pho_) then
         s1min = (max(MPhotonCutoff,getMass(ids(3))+getMass(ids(4))))**2
      else
         s1min = getMass(ids(3))+getMass(ids(4))
      endif
      if (CoupledVertex(ids(5:6),-1,useAHcoupl=1).eq.Pho_) then
         s2min = (max(MPhotonCutoff,getMass(ids(5))+getMass(ids(6))))**2
      else
         s2min = getMass(ids(5))+getMass(ids(6))
      endif
   else
      if (CoupledVertex((/ids(3),ids(6)/),-1,useAHcoupl=1).eq.Pho_) then
         s1min = (max(MPhotonCutoff,getMass(ids(3))+getMass(ids(6))))**2
      else
         s1min = getMass(ids(3))+getMass(ids(6))
      endif
      if (CoupledVertex((/ids(5),ids(4)/),-1,useAHcoupl=1).eq.Pho_) then
         s2min = (max(MPhotonCutoff,getMass(ids(5))+getMass(ids(4))))**2
      else
         s2min = getMass(ids(5))+getMass(ids(4))
      endif
   endif
   Emin = max(Emin,sqrt(max(s1min,s2min)))
   ! s1,2max=s34

   ! Find s34 = Energy**2
   BWmass_ps=-1d0
   BWwidth_ps=-1d0
   if (Process.eq.73) then
      if( M_Reso.ge.0d0 .and. M_Reso2.ge.0d0 ) then ! Both resonances are present
         BWmass_ps = max(M_Reso,M_Reso2)
         BWwidth_ps = max(max(abs(M_Reso-M_Reso2),Ga_Reso), Ga_Reso2) ! Cover the full mass difference
         if (Emin.gt.max(M_Reso+10d0*Ga_Reso,M_Reso2+10d0*Ga_Reso2) .or. Emax.lt.min(M_Reso-10d0*Ga_Reso,M_Reso2-10d0*Ga_Reso2)) then
            BWmass_ps = -1d0
            BWwidth_ps = -1d0
         endif
      else if( M_Reso.ge.0d0 ) then
         BWmass_ps = M_Reso
         BWwidth_ps = Ga_Reso
         if (Emin.gt.(BWmass_ps+10d0*BWwidth_ps) .or. Emax.lt.(BWmass_ps-10d0*BWwidth_ps)) then
            BWmass_ps = -1d0
            BWwidth_ps = -1d0
         endif
      else if( M_Reso2.ge.0d0 ) then
         BWmass_ps = M_Reso2
         BWwidth_ps = Ga_Reso2
         if (Emin.gt.(BWmass_ps+10d0*BWwidth_ps) .or. Emax.lt.(BWmass_ps-10d0*BWwidth_ps)) then
            BWmass_ps = -1d0
            BWwidth_ps = -1d0
         endif
      endif
   else if (Process.eq.75) then
      if( M_Reso.ge.0d0 .and. M_Reso2.ge.0d0 ) then ! Both resonances are present
         BWmass_ps = max(M_Reso,M_Reso2)
         BWwidth_ps = max(abs(M_Reso-M_Reso2), abs(M_Reso-2d0*M_V_ps), abs(M_Reso2-2d0*M_V_ps), Ga_Reso, Ga_Reso2) ! Cover the full mass difference
         if (Emin.gt.max(M_Reso+10d0*Ga_Reso,M_Reso2+10d0*Ga_Reso2) .or. Emax.lt.min(M_Reso-10d0*Ga_Reso,M_Reso2-10d0*Ga_Reso2)) then
            BWmass_ps = -1d0
            BWwidth_ps = -1d0
         endif
      else if( M_Reso.ge.0d0 ) then
         BWmass_ps = M_Reso
         BWwidth_ps = max(abs(M_Reso-2d0*M_V_ps), Ga_Reso)
         if (Emin.gt.(BWmass_ps+10d0*BWwidth_ps) .or. Emax.lt.(BWmass_ps-10d0*BWwidth_ps)) then
            BWmass_ps = -1d0
            BWwidth_ps = -1d0
         endif
      else if( M_Reso2.ge.0d0 ) then
         BWmass_ps = M_Reso2
         BWwidth_ps = max(abs(M_Reso2-2d0*M_V_ps), Ga_Reso2)
         if (Emin.gt.(BWmass_ps+10d0*BWwidth_ps) .or. Emax.lt.(BWmass_ps-10d0*BWwidth_ps)) then
            BWmass_ps = -1d0
            BWwidth_ps = -1d0
         endif
      endif
   else if (Process.eq.74) then
      BWmass_ps = 2d0*M_V_ps
      BWwidth_ps = max(2d0*M_V_ps, Ga_V_ps)
      if (Emin.gt.(BWmass_ps+2d0*BWwidth_ps) .or. Emax.lt.(BWmass_ps-10d0*Ga_V_ps)) then
         BWmass_ps = -1d0
         BWwidth_ps = -1d0
      endif
   endif
   if( BWmass_ps.lt.0d0 .or. BWwidth_ps.lt.0d0 ) then ! Create flat 4f mass
      Jac1 = k_l(xRnd(1),Emin**2,Emax**2,s34)
   else ! Create a BW 4f mass
      Jac1 = k_BreitWigner(xRnd(1),BWmass_ps**2,BWwidth_ps,Emin**2,Emax**2,s34)
   endif
   Jac1 = Jac1 / (Collider_Energy**2) ! This is because Jac1*Jac2 should be divided by s
   if (s34.le.0d0) then
      return
   endif
   Energy = sqrt(s34)

   ! Find y and set eta1, eta2
   Jac2 = rapidity_tan_map(xRnd(2),sysY,ywidthset=sqrt(2d0))
   eta1 = Energy/Collider_Energy*exp(sysY)  ! x1
   eta2 = Energy/Collider_Energy*exp(-sysY) ! x2
   if (eta1.ge.1d0 .or. eta2.ge.1d0) then
      return
   endif

   ! Begin four-momenta
   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

   Jac3 = s_channel_propagator(M_V_ps**2,Ga_V_ps,s1min,s34,xRnd(3),s56) ! Find s56=m1**2
   Jac4 = s_channel_propagator(M_V_ps**2,Ga_V_ps,s2min,(dsqrt(s34)-dsqrt(s56))**2,xRnd(4),s78) ! Find s78=m2**2
   Jac5 = s_channel_decay((/Energy,0d0,0d0,0d0/),s56,s78,xRnd(5:6),Mom(:,V1),Mom(:,V2)) ! Decay pVV to pV1, pV2 in CoM
   if( swap34_56 ) then ! Swapped config.
      Jac6 = s_channel_decay(Mom(:,V1),0d0,0d0,xRnd(7:8),Mom(:,Lep1P),Mom(:,Lep2M)) ! Decay pV1
      Jac7 = s_channel_decay(Mom(:,V2),0d0,0d0,xRnd(9:10),Mom(:,Lep2P),Mom(:,Lep1M)) ! Decay pV2
   else ! Normal config
      Jac6 = s_channel_decay(Mom(:,V1),0d0,0d0,xRnd(7:8),Mom(:,Lep1P),Mom(:,Lep1M)) ! Decay pV1
      Jac7 = s_channel_decay(Mom(:,V2),0d0,0d0,xRnd(9:10),Mom(:,Lep2P),Mom(:,Lep2M)) ! Decay pV2
   endif

   Jac = Jac1*Jac2*Jac3*Jac4*Jac5*Jac6*Jac7 * PSNorm4

   if( isNan(jac) ) then
      print *, "EvalPhasespace_gg4f NaN"
      print *, Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7
      if( isNan(jac) ) Jac = 0d0

      write(6,*) "ids=",ids
      write(6,*) "m34:",sqrt(s34)/GeV
      write(6,*) "m56:",sqrt(s56)/GeV
      write(6,*) "m78:",sqrt(s78)/GeV

      call Exit(1)
   endif

RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_H4f(xchannel,xRnd,Energy,Mom,id,Jac)
use ModParameters
use ModPhasespace
use ModMisc
implicit none
real(8) :: xchannel,xRnd(1:8), Energy, Mom(4,8)
integer, intent(in) :: id(1:4) ! 4f ids
integer :: iChannel
real(8) :: Jac,Jac1,Jac2,Jac3,Jac4,Jac5,Jac6,Jac7,Jac8,Jac9
real(8) :: s3H,s4H,s56,s78,s910,Mom_Dummy(1:4),xRndOffShellZ,vectormass(1:2,1:2)
integer :: NumChannels

   NumChannels=4 ! IMPORTANT: Have to set it here rather than declaration!
   iChannel = 1
   Jac1=1d0;Jac2=1d0;Jac3=1d0;Jac4=1d0;Jac5=1d0;Jac6=1d0;Jac7=1d0;Jac8=1d0;Jac9=1d0;Jac=1d0;
   s56=Energy**2

   vectormass(1,1)=M_V_ps  ! Test is already done in main.90
   vectormass(1,2)=Ga_V_ps ! Test is already done in main.90
   if( IsAZDecay(DecayMode2) ) then
      vectormass(2,1)=M_Z_ps
      vectormass(2,2)=Ga_Z_ps
   elseif( IsAWDecay(DecayMode2) ) then
      vectormass(2,1)=M_W_ps
      vectormass(2,2)=Ga_W_ps
   elseif( IsAPhoton(DecayMode2) ) then
      vectormass(2,1)=0d0
      vectormass(2,2)=0d0
   endif

   Mom(1:4,1) = 0.5d0*Energy * (/+1d0,0d0,0d0,+1d0/)
   Mom(1:4,2) = 0.5d0*Energy * (/+1d0,0d0,0d0,-1d0/)

   ! If does not have 4f-interference, then use the first two channels only
   if(.not.includeInterference .or. (IsAWDecay(DecayMode1).and.IsAWDecay(DecayMode2)) .or. .not.(id(1).eq.id(3) .and. id(2).eq.id(4)) ) then
      NumChannels=2
   endif

   ! Ordering becomes important if the on-shellness of V1 and V2 are not specified to be the same.
   ! If one V is off-shell and the other is on-shell, or if exactly one decay mode is a photon, use the first channel.
   if((OffShellV1.neqv.OffShellV2) .or. ((IsAPhoton(DecayMode1) .or. IsAPhoton(DecayMode2)) .and. DecayMode1.ne.DecayMode2 )) then
      !print *,"Passes this if-statement!"
      NumChannels=1
   endif
   if(NumChannels.gt.1) then
      iChannel = int(xchannel * NumChannels -1d-10)+1
   endif

   !print *, "PS channel / NumChannels ",iChannel,NumChannels

   ! masses
   if(OffShellV1) then
      if(OffShellV2) then ! OffXVV=X11
         Jac2 = s_channel_propagator((vectormass(1,1))**2,vectormass(1,2),0d0,s56,xRnd(1),s78)                                          !  int d(s78)    = Z1
         Jac3 = s_channel_propagator((vectormass(2,1))**2,vectormass(2,2),0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)                  !  int d(s910)   = Z2
      else                ! OffXVV=X10
         s910 = (vectormass(2,1))**2                                                                                                    !  int d(s910)   = Z2
         Jac2 = s_channel_propagator((vectormass(1,1))**2,vectormass(1,2),0d0,(Energy-dsqrt(s910))**2,xRnd(1),s78)                  !  int d(s78)    = Z1
      endif
   else
      s78 = (vectormass(1,1))**2                                                                                                        !  int d(s78)    = Z1
      if(OffShellV2) then ! OffXVV=X01
         Jac3 = s_channel_propagator((vectormass(2,1))**2,vectormass(2,2),0d0,(Energy-dsqrt(s78))**2,xRnd(2),s910)                  !  int d(s910)   = Z2
      else                ! OffXVV=X00
         s910 = (vectormass(2,1))**2                                                                                                    !  int d(s910)   = Z2
         if( s910.gt.((dsqrt(s56)-dsqrt(s910))**2) ) Jac3=0d0
      endif
   endif

   !print *, "x",xrnd(1:2)
   !print *, "s",s56,s78,s910
   !pause

   Jac = Jac1*Jac2*Jac3
   if(Jac.eq.0d0) return ! Checkpoint for on/off-shellness

   Mom_Dummy(1:4) = (/Energy,0d0,0d0,0d0/)

!  splittings
IF( iChannel.EQ.1 ) THEN
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,3),Mom(:,4))                                                   !  H --> 5+6
   if(.not.IsAPhoton(DecayMode1)) then
      Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,5),Mom(:,6))                                                       !  5 --> 7+8
   else
      Mom(:,5) = Mom(:,3)
      Mom(:,6) = 0d0
   endif
   if(.not.IsAPhoton(DecayMode2)) then
      Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,7),Mom(:,8))                                                       !  6 --> 9+10
   else
      Mom(:,7) = Mom(:,4)
      Mom(:,8) = 0d0
   endif
ELSEIF( iChannel.EQ.2 ) THEN
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,4),Mom(:,3))                                                   !  H --> 6+5
   if(.not.IsAPhoton(DecayMode2)) then
      Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,5),Mom(:,6))                                                       !  5 --> 7+8
   else
      Mom(:,5) = Mom(:,3)
      Mom(:,6) = 0d0
   endif
   if(.not.IsAPhoton(DecayMode1)) then
      Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,7),Mom(:,8))                                                       !  6 --> 9+10
   else
      Mom(:,7) = Mom(:,4)
      Mom(:,8) = 0d0
   endif

ELSEIF( iChannel.EQ.3 ) THEN                                                                                                     !  Interference scenarios
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,3),Mom(:,4))                                                   !  H --> 5+6
   Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,7),Mom(:,6))                                                          !  5 --> 9+8
   Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,5),Mom(:,8))                                                          !  6 --> 7+10
ELSEIF( iChannel.EQ.4 ) THEN
   Jac4 = s_channel_decay(Mom_Dummy(1:4),s78,s910,xRnd(3:4),Mom(:,4),Mom(:,3))                                                   !  H --> 6+5
   Jac5 = s_channel_decay(Mom(:,3),0d0,0d0,xRnd(5:6),Mom(:,7),Mom(:,6))                                                          !  5 --> 9+8
   Jac6 = s_channel_decay(Mom(:,4),0d0,0d0,xRnd(7:8),Mom(:,5),Mom(:,8))                                                          !  6 --> 7+10
ENDIF

   Jac = Jac*Jac4*Jac5*Jac6 * PSNorm4  !*NumChannels                                                                                 !  combine

   !print *, energy,dsqrt(s56)
   !print *,"Generated momenta: "
   !print *,Mom(:,3)
   !print *,Mom(:,4)
   !print *,Mom(:,5)
   !print *,Mom(:,6)
   !print *,Mom(:,7)
   !print *,Mom(:,8)
   !pause


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
!    print *, "Inv.mass  ",get_MInv(Mom_Dummy(1:4))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,3)+Mom(1:4,4))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,3))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,4))/GeV
!    pause



RETURN
END SUBROUTINE

SUBROUTINE EvalPhasespace_H4f_OLD(xchannel,xRnd,Energy,Mom,Jac)
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
!    print *, "Inv.mass  ",get_MInv(Mom_Dummy(1:4))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,3)+Mom(1:4,4))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,3))/GeV
!    print *, "Inv.mass  ",get_MInv(Mom(1:4,4))/GeV
!    pause



RETURN
END SUBROUTINE

SUBROUTINE EvalPhasespace_HVga(xRnd,Energy,Mom,Jac) ! H-->gaga or H-->Zga phase space
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







subroutine SetRunningScales(p,id) ! p in JHU-GeV, id in JHUGen conventions
use ModParameters
use ModMisc
implicit none
real(dp), intent(in) :: p(1:4,4:6) ! No need to run the second index from 4 to 7: pH, pJ1, pJ2
integer, intent(in) :: id(4:7) ! id_JJH/id_JJVV, id_J1, id_J2, id_JJ (if applicable)
real(8) :: polemass(3:7) ! mJJH, mH, mJ1, mJ2, mJJ (if applicable)
real(8) :: pJJHstar(4),pHstar(4),pJ(4,2),pJJ(4),pJHstar(4),pTjet(5:6),maxpTjet,minpTjet
integer idx,ip

   pHstar(:) = 0d0
   pJJ(:) = 0d0
   pTjet(:) = 0d0
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
         pTjet(idx) = get_PT(p(1:4,idx))
      endif
   enddo
   maxpTjet = maxval(pTjet)
   minpTjet = minval(pTjet, mask=.not.all(p(1:4,5:6).eq.0d0, 1))
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
   elseif(FacScheme .eq. kRenFacScheme_maxpTj) then
      Mu_Fact = maxpTjet
   elseif(FacScheme .eq. kRenFacScheme_minpTj) then
      Mu_Fact = minpTjet
   else
      call Error("This should never be able to happen.")
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
   elseif(RenScheme .eq. kRenFacScheme_maxpTj) then
      Mu_Ren = maxpTjet
   elseif(RenScheme .eq. kRenFacScheme_minpTj) then
      Mu_Ren = minpTjet
   else
      call Error("This should never be able to happen.")
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

        PDFScale=Mu_Fact/GeV
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



!
!
!! QCD scale from MCFM
!! Implementation into JHUGen by Ulascan Sarica, Dec. 2015
!subroutine EvalAlphaS()
!   use ModParameters
!   IMPLICIT NONE
!#if useLHAPDF==1
!!--- This is simply a wrapper to the LHAPDF implementation of the running coupling alphas, in the style of the native MCFM routine
!   DOUBLE PRECISION alphasPDF
!   REAL(DP) :: Q
!      Q = Mu_Ren/GeV
!      alphas=alphasPDF(Q)
!#else
!!     Evaluation of strong coupling constant alphas
!!     Original Author: R.K. Ellis
!!     q -- Scale at which alpha_s is to be evaluated
!!     alphas_mz -- ModParameters value of alpha_s at the mass of the Z-boson
!!     nloops_pdf -- ModParameters value of the number of loops (1,2, or 3) at which the beta function is evaluated to determine running.
!!     If you somehow need a more complete implementation, check everything at or before commit 28472c5bfee128dde458fd4929b4d3ece9519ab8
!   INTEGER, PARAMETER :: NF6=6
!   INTEGER, PARAMETER :: NF5=5
!   INTEGER, PARAMETER :: NF4=4
!   INTEGER, PARAMETER :: NF3=3
!   INTEGER, PARAMETER :: NF2=2
!   INTEGER, PARAMETER :: NF1=1
!
!      IF (Mu_Ren .LE. 0d0) THEN
!         WRITE(6,*) 'ModKinematics::EvalAlphaS: Mu_Ren .le. 0, Mu_Ren (GeV) = ',(Mu_Ren*GeV)
!         stop
!      ENDIF
!      IF (nQflavors_pdf .NE. NF5) THEN
!         WRITE(6,*) 'ModKinematics::EvalAlphaS: nQflavors_pdf invalid, nQflavors_pdf = ',nQflavors_pdf
!         WRITE(6,*) 'ModKinematics::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
!         stop
!      ENDIF
!      IF (nloops_pdf .NE. 1) THEN
!         WRITE(6,*) 'ModKinematics::EvalAlphaS: nloops_pdf invalid, nloops_pdf = ',nloops_pdf
!         WRITE(6,*) 'ModKinematics::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
!         stop
!      ENDIF
!
!      alphas=alphas_mz/(1.0_dp+alphas_mz*B0_PDF(NF5)*2.0_dp*dlog((Mu_Ren/zmass_pdf)))
!#endif
!      ! Calculate the derived couplings
!      call ComputeQCDVariables()
!   RETURN
!end subroutine EvalAlphaS
!
!
!
!!-----------------------------------------------------------------------------
!!
!      real(8) function massfrun(mf,scale)
!!
!!-----------------------------------------------------------------------------
!!
!!       This function returns the 'nloop' value of a MSbar fermion mass
!!       at a given scale.
!!
!!       INPUT: mf    = MSbar mass of fermion at MSbar fermion mass scale
!!              scale = scale at which the running mass is evaluated
!!              asmz  = AS(MZ) : this is passed to alphas(scale,asmz,2)
!!              nloop = # of loops in the evolutionC
!!
!!       COMMON BLOCKS: COMMON/QMASS/CMASS,BMASS,TMASS
!!                      contains the MS-bar masses of the heavy quarks.
!!
!!       EXTERNAL:      double precision alphas(scale,asmz,2)
!!
!!-----------------------------------------------------------------------------
!!
!      use ModParameters
!      implicit none
!!
!!     ARGUMENTS
!!
!      real(8), intent(in) :: mf, scale
!      real(8) scale_temp_ren
!      !integer , intent(in) :: nloop
!!
!!     LOCAL
!!
!      real(8)  beta0, beta1,gamma0,gamma1
!      real(8)  as,asmf,l2
!      integer  nfrun
!
!      scale_temp_ren=Mu_Ren
!!
!!     EXTERNAL
!!
!!      double precision  alphas
!!      external          alphas
!!
!!     COMMON
!!
!!      real *8      cmass,bmass,tmass
!!      COMMON/QMASS/CMASS,BMASS,TMASS
!!
!!     CONSTANTS
!!
!!      double precision  One, Two, Three, Pi
!      !parameter( One = 1d0, Two = 2d0, Three = 3d0 )
!      !parameter( Pi = 3.14159265358979323846d0)
!
!      if ( mf.gt.m_top ) then
!         nfrun = 6
!      else
!         nfrun = 5
!      end if
!
!      beta0 = ( 11d0 - 2d0/3d0 *nfrun )/4d0
!      !beta1 = ( 102d0  - 38d0/3d0*nf )/16d0
!      gamma0= 1d0
!      !gamma1= ( 202d0/3d0  - 20d0/9d0*nf )/16d0
!      !A1    = -beta1*gamma0/beta0**2+gamma1/beta0
!      Mu_Ren=scale
!      call EvalAlphaS()
!      as=alphas
!
!      Mu_Ren=mf
!      call EvalAlphaS()
!      asmf=alphas
!      !l2    = (1d0+A1*as/Pi)/(one+A1*asmf/Pi)
!
!      massfrun = mf * (as/asmf)**(gamma0/beta0)
!
!      !if(nloop.eq.2) massfrun=massfrun*l2
!
!      Mu_Ren=scale_temp_ren
!      call EvalAlphaS()
!
!      return
!      end function massfrun
!





END MODULE
