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



SUBROUTINE WriteOutEvent(Mom,MY_IDUP,ICOLUP)
use ModParameters
implicit none
real(8) :: Mom(1:4,1:6)
real(8) :: Spin, Lifetime
real(8) :: XFV(1:4), Z1FV(1:4), Z2FV(1:4)
real(8) :: MomDummy(1:4,1:6)
real(8) :: Part1Mass,Part2Mass,XMass,V1Mass,V2Mass,L11Mass,L12Mass,L21Mass,L22Mass,tmp
integer :: a,b,c
integer :: MY_IDUP(1:9),LHE_IDUP(1:9),i,ISTUP(1:9),MOTHUP(1:2,1:9),ICOLUP(1:2,1:9)
integer :: NUP,IDPRUP
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP
real(8) :: ntRnd
character(len=*),parameter :: fmt1 = "(I3,X,I2,X,I2,X,I2,X,I3,X,I3,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7,X,1PE14.7)"

do i=1,9
!        print *, "my_idup(i) ",MY_IDUP(i)
LHE_IDUP(i) = convertLHE( MY_IDUP(i) )
!        print *, "LHE_IDUP(i) ",LHE_IDUP(i)
enddo

! NUP changes for gamma gamma final state
if (LHE_IDUP(4).eq.22) then
NUP=5
else
NUP=9
endif

IDPRUP=100
XWGTUP=1.
SCALUP=1000.
AQEDUP=alpha_QED
AQCDUP=0d0

ISTUP(1) = - 1
ISTUP(2) =  -1
ISTUP(3) = 2
ISTUP(4) = 2
ISTUP(5) = 2
ISTUP(6) = 1
ISTUP(7) = 1
ISTUP(8) = 1
ISTUP(9) = 1

if (LHE_IDUP(4).eq.22) then
ISTUP(4) = 1
ISTUP(5) = 1
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

! Added by Nhan
LHE_IDUP(3) = 39 ! X particle
Lifetime = 0.0
Spin = 1.0

do a=1,6
MomDummy(1,a) = 100.0d0*Mom(1,a)
MomDummy(2,a) = 100.0d0*Mom(2,a)
MomDummy(3,a) = 100.0d0*Mom(3,a)
MomDummy(4,a) = 100.0d0*Mom(4,a)
enddo

do b=1,4
Z1FV(b) = MomDummy(b,3)+MomDummy(b,4)
Z2FV(b) = MomDummy(b,5)+MomDummy(b,6)
enddo

do c=1,4
XFV(c) = Z1FV(c) + Z2FV(c)
enddo

tmp = MomDummy(1,1)*MomDummy(1,1)-MomDummy(2,1)*MomDummy(2,1)-MomDummy(3,1)*MomDummy(3,1)-MomDummy(4,1)*MomDummy(4,1)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
Part1Mass = dSQRT(dabs(tmp))

tmp = MomDummy(1,2)*MomDummy(1,2)-MomDummy(2,2)*MomDummy(2,2)-MomDummy(3,2)*MomDummy(3,2)-MomDummy(4,2)*MomDummy(4,2)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
Part2Mass = dSQRT(dabs(tmp))

tmp = XFV(1)*XFV(1)-XFV(2)*XFV(2)-XFV(3)*XFV(3)-XFV(4)*XFV(4)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
XMass = dSQRT(dabs(tmp))

tmp = Z1FV(1)*Z1FV(1)-Z1FV(2)*Z1FV(2)-Z1FV(3)*Z1FV(3)-Z1FV(4)*Z1FV(4)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
V1Mass = dSQRT(dabs(tmp))
if( V1Mass.lt.1d-5 ) then
V1Mass=0d0
endif

tmp = Z2FV(1)*Z2FV(1)-Z2FV(2)*Z2FV(2)-Z2FV(3)*Z2FV(3)-Z2FV(4)*Z2FV(4)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
V2Mass = dSQRT(dabs(tmp))
if( V2Mass.lt.1d-5 ) then
V2Mass=0d0
endif

tmp = MomDummy(1,3)*MomDummy(1,3)-MomDummy(2,3)*MomDummy(2,3)-MomDummy(3,3)*MomDummy(3,3)-MomDummy(4,3)*MomDummy(4,3)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
L12Mass = dSQRT(dABS(tmp))
if( L12Mass.lt.1d-5 ) then
L12Mass=0d0
endif
if( tmp.lt.0d0 ) then
MomDummy(1,3) = MomDummy(1,3) + 1d-7
endif

tmp = MomDummy(1,4)*MomDummy(1,4)-MomDummy(2,4)*MomDummy(2,4)-MomDummy(3,4)*MomDummy(3,4)-MomDummy(4,4)*MomDummy(4,4)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
L11Mass = dSQRT(dABS(tmp))
if( L11Mass.lt.1d-5 ) then
L11Mass=0d0
endif
if( tmp.lt.0d0 ) then
MomDummy(1,4) = MomDummy(1,4) + 1d-7
endif

tmp = MomDummy(1,5)*MomDummy(1,5)-MomDummy(2,5)*MomDummy(2,5)-MomDummy(3,5)*MomDummy(3,5)-MomDummy(4,5)*MomDummy(4,5)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
L22Mass = dSQRT(dABS(tmp))
if( L22Mass.lt.1d-5 ) then
L22Mass=0d0
endif
if( tmp.lt.0d0 ) then
MomDummy(1,5) = MomDummy(1,5) + 1d-7
endif

tmp = MomDummy(1,6)*MomDummy(1,6)-MomDummy(2,6)*MomDummy(2,6)-MomDummy(3,6)*MomDummy(3,6)-MomDummy(4,6)*MomDummy(4,6)
if( tmp.lt. -1d-3 ) print *, "Error: large negative mass!"
L21Mass = dSQRT(dABS(tmp))
if( L21Mass.lt.1d-5 ) then
L21Mass=0d0
endif
if( tmp.lt.0d0 ) then
MomDummy(1,6) = MomDummy(1,6) + 1d-7
endif


write(14,"(A)") "<event>"
write(14,"(I1,X,I3,X,1PE13.7,X,1PE13.7,X,1PE13.7,X,1PE13.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP


! parton_a
i=1
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,1),MomDummy(1,1),Part1Mass,Lifetime,Spin

! parton_b
i=2
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,2),MomDummy(1,2),Part2Mass,Lifetime,Spin

! X
i=3
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),XFV(2:4),XFV(1),XMass,Lifetime,Spin

! V1
i=4
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Z1FV(2:4),Z1FV(1),V1Mass,Lifetime,Spin

! V2
i=5
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Z2FV(2:4),Z2FV(1),V2Mass,Lifetime,Spin


! decay product 1 (V1): l-, nu or q
i=7
if (LHE_IDUP(i).gt.-9000) then
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,3),MomDummy(1,3),L12Mass,Lifetime,Spin
endif

! decay product 2 (V1): l+, nubar or qbar
i=6
if (LHE_IDUP(i).gt.-9000) then
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,4),MomDummy(1,4),L11Mass,Lifetime,Spin
endif

! decay product 1 (V2): l-, nu or q
i=9
if (LHE_IDUP(i).gt.-9000) then
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,5),MomDummy(1,5),L22Mass,Lifetime,Spin
endif

! decay product 2 (V2): l+, nubar or qbar
i=8
if (LHE_IDUP(i).gt.-9000) then
write(14,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),MomDummy(2:4,6),MomDummy(1,6),L21Mass,Lifetime,Spin
endif

write(14,"(A)") "</event>"

! print * ,"check ", LHE_IDUP(6),MomDummy(1:4,4)
! print * ,"check ", LHE_IDUP(7),MomDummy(1:4,3)
! print * ,"check ", LHE_IDUP(8),MomDummy(1:4,6)
! print * ,"check ", LHE_IDUP(9),MomDummy(1:4,5)
! pause

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


      call genps(2,MV,xRndPS(1:2),(/ML1,ML2/),MomDK(1:4,1:2),PSWgt2)

!     boost all guys to the V boson frame:
      call boost(MomDK(1:4,1),VMom(1:4),MV)
      call boost(MomDK(1:4,2),VMom(1:4),MV)
      PSWgt = PSWgt2*PiWgt2

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

    if ((MG/M_Reso*eta1.gt.1).or.(MG/M_Reso*eta2.gt.1)) return

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
real(8) :: MomExt(:,:),MomDK(:,:), mZ1, mZ2, MReso
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),BeamAxis(1:4),ScatteringAxis(1:4),dummy(1:4)
real(8) :: MomLept(1:4,1:4),MomLeptX(1:4,1:4),MomLeptPlane1(2:4),MomLeptPlane2(2:4),MomBeamScatterPlane(2:4)
logical :: applyPSCut
integer :: NumPart,NBin(:)
real(8) :: pT_lepM,pT_lepP,y_lepM,y_lepP,MomFerm(1:4),MomZ2(1:4),MomReso(1:4),CosTheta1,Phi,Phi1,signPhi,signPhi1
real(8) :: CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,CosThetaZ,CosThetaStar


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





!   compute pT of leptons in the lab frame
      pT_lepP = get_PT(MomLept(1:4,2))
      pT_lepM = get_PT(MomLept(1:4,3))

      y_lepP = get_eta(MomLept(1:4,2))
      y_lepM = get_eta(MomLept(1:4,3))






! construct cos(theta1): angle between direction of fermion from Z1 and negative direction of opposite Z in Z1 rest frame

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

return
END SUBROUTINE


FUNCTION ZQuaBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZQuaBranching

  if( xRnd .le. Br_Z_up ) then
      ZQuaBranching = Up_
  elseif(xRnd .le. Br_Z_up+Br_Z_ch) then
      ZQuaBranching = Chm_
  elseif(xRnd .le. Br_Z_up+Br_Z_ch+Br_Z_dn) then
      ZQuaBranching = Dn_
  elseif(xRnd .le. Br_Z_up+Br_Z_ch+Br_Z_dn+Br_Z_st) then
      ZQuaBranching = Str_
  elseif(xRnd .le. Br_Z_up+Br_Z_ch+Br_Z_dn+Br_Z_st+Br_Z_bo) then
      ZQuaBranching = Bot_
  else
      print *, "error ",xRnd
      stop
  endif

RETURN
END FUNCTION



FUNCTION ZLepBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZLepBranching

  if( xRnd .le. 0.5d0 ) then
      ZLepBranching = ElM_
  elseif(xRnd .le. 1.0d0) then
      ZLepBranching = MuM_
  else
      print *, "error ",xRnd
      stop
  endif

RETURN
END FUNCTION



FUNCTION ZNuBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: ZNuBranching

  if( xRnd .le. 1d0/3d0 ) then
      ZNuBranching = NuE_
  elseif(xRnd .le. 2.0d0/3d0) then
      ZNuBranching = NuM_
  elseif(xRnd .le. 1.0d0) then
      ZNuBranching = NuT_
  else
      print *, "error ",xRnd
      stop
  endif

RETURN
END FUNCTION





FUNCTION WLepBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WLepBranching

  if( xRnd .le. 0.5d0 ) then
      WLepBranching = ElM_
  elseif(xRnd .le. 1.0d0) then
      WLepBranching = MuM_
  else
      print *, "error ",xRnd
      stop
  endif

RETURN
END FUNCTION


FUNCTION WQuaUpBranching(xRnd)
use ModParameters
implicit none
real(8) :: xRnd
integer :: WQuaUpBranching

  if( xRnd .le. 0.5d0 ) then
      WQuaUpBranching = Up_
  elseif(xRnd .le. 1d0 ) then
      WQuaUpBranching = Chm_
  else
      print *, "error ",xRnd
      stop
  endif

RETURN
END FUNCTION




SUBROUTINE VVBranchings(MY_IDUP,ICOLUP)
use ModParameters
implicit none
integer :: MY_IDUP(4:9),ICOLUP(1:2,6:9),DKFlavor
real(8) :: DKRnd

!    particle associations:
!    
!    IDUP(6)  -->  MomDK(:,2)  -->     v-spinor
!    IDUP(7)  -->  MomDK(:,1)  -->  ubar-spinor
!    IDUP(8)  -->  MomDK(:,4)  -->     v-spinor
!    IDUP(9)  -->  MomDK(:,3)  -->  ubar-spinor
!


   if( DecayMode1.eq.0 ) then! Z1->2l
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZLepBranching( DKRnd )!= ElM or MuM
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
   elseif( DecayMode1.eq.1 ) then! Z1->2q
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZQuaBranching( DKRnd )!= Up,Dn,Chm,Str,Bot
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
        ICOLUP(1:2,6) = (/0,503/)
        ICOLUP(1:2,7) = (/503,0/)
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
        MY_IDUP(6) = -abs(DKFlavor)-1  ! anti-dn flavor
        MY_IDUP(7) = +abs(DKFlavor)    ! up flavor
        ICOLUP(1:2,6) = (/0,503/)
        ICOLUP(1:2,7) = (/503,0/)
   elseif( DecayMode1.eq.6 ) then! W1(+)->taunu
        MY_IDUP(4) = Wp_
        MY_IDUP(6) = TaP_
        MY_IDUP(7) = NuT_
   elseif( DecayMode1.eq.7 ) then! photon
        MY_IDUP(4) = Pho_
        MY_IDUP(6) = -9999
        MY_IDUP(7) = -9999
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
        DKFlavor = ZQuaBranching( DKRnd )!= Up,Dn,Chm,Str,Bot
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
        ICOLUP(1:2,8) = (/0,504/)
        ICOLUP(1:2,9) = (/504,0/)
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
        MY_IDUP(8) = -abs(DKFlavor)    ! anti-up flavor
        MY_IDUP(9) = +abs(DKFlavor)+1  ! dn flavor
        ICOLUP(1:2,8) = (/0,504/)
        ICOLUP(1:2,9) = (/504,0/)
   elseif( DecayMode2.eq.6 ) then! W2(-)->taunu
        MY_IDUP(5) = Wm_
        MY_IDUP(8) = ANuT_
        MY_IDUP(9) = TaM_
   elseif( DecayMode2.eq.7 ) then! photon
        MY_IDUP(5) = Pho_
        MY_IDUP(8) = -9999
        MY_IDUP(9) = -9999
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





FUNCTION WhichXBin(NHisto,XValue)
use ModParameters
implicit none
integer :: WhichXBin,NHisto
real(8) :: XValue
integer :: i
include "vegas_common.f"

    whichxbin = int( xValue*NPart )!  uniform distribution


!    do i=1,50!                         distribution according to vegas grid
!       if( XValue .lt. xi(i,NHisto) ) then
!           WhichXBin=i
!           return
!       endif
!    enddo
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




SUBROUTINE PDFMapping(MapType,yRnd,eta1,eta2,Ehat,sHatJacobi)
use ModParameters
use ModMisc
implicit none
integer :: MapType
real(8) :: yRnd(1:2),eta1,eta2,EHat,sHatJacobi,tau,nPotMap,z,sbar,fmax
real(8) :: etamin, Ymax, Y

  if( MapType.eq.1 ) then!  no mapping
      eta1 = yRnd(1)
      eta2 = yRnd(2)
      sHatJacobi = 1d0
  elseif( MapType.eq.2 ) then!  exponential mapping
!       tau = (2d0*m_Top/Collider_Energy)**2
!       eta1 = tau**yRnd(1)
!       eta2 = tau**( (1d0-yRnd(1))*yRnd(2) )
!       sHatJacobi = dlog(tau)**2*(1d0-yRnd(1))*eta1*eta2
  elseif( MapType.eq.3 ) then!  linear mapping
!       tau = (2d0*m_Top/Collider_Energy)**2
!       eta1 = (1d0-tau)*yRnd(1) + tau
!       eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
!       sHatJacobi = (1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)
  elseif( MapType.eq.4 ) then!  MCFM mapping
!       tau = dexp(dlog(((2d0*m_Top/Collider_Energy)**2))*yRnd(1))
!       eta1 = dsqrt(tau)*dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
!       eta2 = dsqrt(tau)/dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
!       sHatJacobi = dlog(((2d0*m_Top/Collider_Energy)**2))*tau*dlog(tau)
  elseif( MapType.eq.5 ) then!  nPotMap mapping
!       nPotMap = 0.5d0
!       tau = (2d0*m_Top/Collider_Energy)**2
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
     Ymax = log(Collider_Energy/M_Reso)
     Y = -Ymax + 2d0*Ymax*yRnd(2)
     eta1 = M_Reso/Collider_Energy*exp(Y)
     eta2 = M_Reso/Collider_Energy*exp(-Y)
     fmax = 0.5d0*pi/M_Reso**3/Ga_Reso*2d0*Ymax
     sHatJacobi = fmax*(M_Reso**2*Ga_Reso**2 )
  else
      call Error("PDF mapping not available")
  endif

  EHat = Collider_Energy*dsqrt(eta1*eta2)

RETURN
END SUBROUTINE







SUBROUTINE setPDFs(x1,x2,MuFac,pdf)
use ModParameters
implicit none
real(8) :: x1,x2,PDFScale,MuFac
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),chm(1:2),bot(1:2),glu(1:2),phot(1:2),sbar(1:2),cbar(1:2),bbar(1:2)
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2)

        PDFScale=MuFac*100d0
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
        else
            print *, "PDFSet",PDFSet,"not available!"
            stop
        endif

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



END MODULE
