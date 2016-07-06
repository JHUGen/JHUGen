MODULE ModKinematics
implicit none
save


type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8),allocatable :: Value(:)
    real(8),allocatable :: Value2(:)
    integer,allocatable :: Hits(:)
    character :: Info*(50)
end type

integer,public :: it_sav
type(Histogram),allocatable :: Histo(:)

contains


SUBROUTINE WriteOutEvent(eta1,eta2,Mom)
implicit none
real(8) :: eta1,eta2,Mom(1:4,1:6)

   write(14 ,"(1X,1F24.16,1X,1F24.16,1X,24F24.16)") eta1,eta2,Mom(1:4,1:6)

END SUBROUTINE



SUBROUTINE EvalPhasespace_ZDecay(ZMom,xRndPS,MomDK,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: ZMom(1:4),MomChk(1:4,1:3)
real(8) :: MomDK(1:4,1:2)
real(8) :: xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


!     MomDK(1:4,i): i= 1:l-, 2:l+
      call genps(2,m_Z,xRndPS(1:2),(/0d0,0d0/),MomDK(1:4,1:2),PSWgt2)
!     boost all guys to the Z frame:
      call boost(MomDK(1:4,1),ZMom(1:4),m_Z)
      call boost(MomDK(1:4,2),ZMom(1:4),m_Z)
      PSWgt = PSWgt2*PiWgt2

RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_2to2Z(EHat,xRndPS,Mom,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_Z,m_Z/),Mom(1:4,3:4),PSWgt)
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







SUBROUTINE Kinematics(NumPart,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
real(8) :: MomExt(:,:),MomDK(:,:)
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),MomZ(1:4),MomG(1:4)
real(8) :: MomLept(1:4,1:4)
logical :: applyPSCut
integer :: NumPart,NBin(:)
real(8) :: pT_lepM,pT_lepP
real(8) :: CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,CosThetaZ


      applyPSCut = .false.



!  eval kinematic variables
!     angle
      MomBoost(1)   =+MomExt(1,3)
      MomBoost(2:4) =-MomExt(2:4,3)

      MomLepP(1:4)  = MomDK(1:4,1)
      call boost(MomLepP(1:4),MomBoost(1:4),m_Z)

      MomZ(1:4) = MomExt(1:4,4)
      call boost(MomZ(1:4),MomBoost(1:4),m_Z)

      CosPhi_LepPZ = (MomLepP(2)*MomZ(2)+MomLepP(3)*MomZ(3)  &
      +MomLepP(4)*MomZ(4))/MomLepP(1)/dsqrt(MomZ(1)**2-m_Z**2)

!     pT of leptons
      pT_lepP = get_PT(MomDK(1:4,1))
      pT_lepM = get_PT(MomDK(1:4,2))

!     angle between lepton planes
!      MomLepP(2:4) = MomDK(2:4,1).cross.MomDK(2:4,2)
!      MomLepP(2:4) = MomLepP(2:4)/dsqrt( MomLepP(2)**2+MomLepP(3)**2+MomLepP(4)**2 )
!      MomLepM(2:4) = MomDK(2:4,3).cross.MomDK(2:4,4)
!      MomLepM(2:4) = MomLepM(2:4)/dsqrt( MomLepM(2)**2+MomLepM(3)**2+MomLepM(4)**2 )

      MomG(1:4)= MomExt(1:4,3) + MomExt(1:4,4)
      MomBoost(1)   =+MomG(1)
      MomBoost(2:4) =-MomG(2:4)
      MomLept = MomDK 
      call boost(MomLept(1:4,1),MomBoost(1:4),m_Grav)
      call boost(MomLept(1:4,2),MomBoost(1:4),m_Grav)
      call boost(MomLept(1:4,3),MomBoost(1:4),m_Grav)
      call boost(MomLept(1:4,4),MomBoost(1:4),m_Grav)

      MomLepP(2:4) = MomLept(2:4,1).cross.MomLept(2:4,2)
      MomLepP(2:4) = MomLepP(2:4)/dsqrt( MomLepP(2)**2+MomLepP(3)**2+MomLepP(4)**2 )
      MomLepM(2:4) = MomLept(2:4,3).cross.MomLept(2:4,4)
      MomLepM(2:4) = MomLepM(2:4)/dsqrt( MomLepM(2)**2+MomLepM(3)**2+MomLepM(4)**2 )

      CosPhi_LepPlanes = acos(MomLepP(2)*MomLepM(2)+MomLepP(3)*MomLepM(3)+MomLepP(4)*MomLepM(4))

!     scattering angle of Z in graviton rest frame
      MomG(1:4)= MomExt(1:4,3) + MomExt(1:4,4)
      MomBoost(1)   =+MomG(1)
      MomBoost(2:4) =-MomG(2:4)
      MomZ(1:4) = MomExt(1:4,4)
      call boost(MomZ(1:4),MomBoost(1:4),m_Grav)
      CosThetaZ = MomZ(4)/dsqrt(MomZ(2)**2+MomZ(3)**2+MomZ(4)**2)


!     binning
      NBin(1) = WhichBin(1,pT_lepP)
      NBin(2) = WhichBin(2,pT_lepM)
      NBin(3) = WhichBin(3,CosPhi_LepPZ)
      NBin(4) = WhichBin(4,CosPhi_LepPlanes)
      NBin(5) = WhichBin(5,CosThetaZ)


return
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




SUBROUTINE PDFMapping(MapType,yRnd,eta1,eta2,Ehat,sHatJacobi)
use ModParameters
use ModMisc
implicit none
integer :: MapType
real(8) :: yRnd(1:2),eta1,eta2,EHat,sHatJacobi,tau,nPotMap,z,sbar,fmax
real(8) :: etamin

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
      fmax = 1d0/m_Grav/Ga_Grav * ( datan((Collider_Energy**2-m_Grav**2)/m_Grav/Ga_Grav) - datan(-m_Grav/Ga_Grav) )
      sbar = m_Grav*Ga_Grav * dtan(fmax*yRnd(1)*m_Grav*Ga_Grav - atan(m_Grav/Ga_Grav) ) + m_Grav**2
      z = sbar/Collider_Energy**2
      eta1 = z + (1d0-z)*yRnd(2)
      eta2 = z/eta1
      sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (Collider_Energy**2*eta1*eta2 - m_Grav**2)**2 + m_Grav**2*Ga_Grav**2 )
  elseif (MapType.eq.11) then ! delta-function map

     etamin = m_Grav**2/Collider_Energy**2
     eta2 = etamin + (1d0-etamin)*yRnd(2)
     eta1 = etamin/eta2
     fmax = 0.5d0*pi/m_Grav**3/Ga_Grav
     sHatJacobi = fmax*etamin/eta2* ( (Collider_Energy**2*eta1*eta2 - m_Grav**2)**2 + m_Grav**2*Ga_Grav**2 )
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
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),chm(1:2),bot(1:2),glu(1:2)
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2)

        PDFScale=MuFac*100d0
        call cteq6(x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
        call cteq6(x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))

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
