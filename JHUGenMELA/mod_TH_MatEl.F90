MODULE modTH
IMPLICIT NONE

 public :: InitProcess_TH      
 public :: EvalXSec_PP_TH,EvalXSec_PP_TBH         ! hadronic
 public :: EvalAmp_QB_TH, EvalAmp_QbarBbar_TbarH  ! t-channel
 public :: EvalAmp_QQB_THBbar,EvalAmp_QQB_TbarHB  ! s-channel
 
 private

 real(8), parameter :: GeV=1d0/100d0 ! we are using units of 100GeV, i.e. Lambda=10 is 1TeV 
 real(8), parameter :: Gf = 1.16639d-5/GeV**2        ! Fermi constant
 real(8), parameter :: vev = 1.0d0/sqrt(Gf*sqrt(2.0d0))

 real(8) :: M_Reso=125d0*GeV
 

 CONTAINS 


 

SUBROUTINE InitProcess_TH(m_Higgs)
implicit none
real(8) :: m_Higgs


  M_Reso = m_Higgs
  
RETURN
END SUBROUTINE
  
  
SUBROUTINE EvalXSec_PP_TH(Mom,TTBHcoupl,TopDecays,Channel,Res)
implicit none
real(8) :: Mom(1:4,1:9),Res
complex(8) :: TTBHcoupl(1:2)
integer :: Channel! 0=s+t channel, 1=t channel, 2=s channel
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: eta1,eta2,Etot,Pztot,MatElSq_GG,MatElSq_QQB,MatElSq_QBQ
real(8) :: x1,x2,PDFScale,Collider_Energy,E_CMS
real(8) :: NNpdf(1:2,-6:7),m_ferm,LO_Res_Unpol(-6:6,-6:6)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'
      

      
      Collider_Energy = -(Mom(1,inLeft)+Mom(1,inRight))
      if( TopDecays.eq.0 ) then
          Etot = Mom(1,Hbos)+Mom(1,t)+Mom(1,qout)
          Pztot= Mom(4,Hbos)+Mom(4,t)+Mom(4,qout)
      else
          Etot = Mom(1,Hbos) + Mom(1,qout) + Mom(1,b)+Mom(1,lep)+Mom(1,nu)
          Pztot= Mom(4,Hbos) + Mom(4,qout) + Mom(4,b)+Mom(4,lep)+Mom(4,nu)
      endif
      x1 = (Etot+Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      x2 = (Etot-Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      E_CMS = dsqrt(x1*x2)*Collider_Energy
      PDFScale = 0.25d0*( m_Top + m_Reso ) * 100d0

      Mom(1:4,inLeft)  = x1 * Mom(1:4,inLeft)
      Mom(1:4,inRight) = x2 * Mom(1:4,inRight)
      
      call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
      call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))
      
      Res = 0d0
      if( Channel.eq.1 .or. Channel.eq.0 ) then
          call EvalAmp_QB_TH(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)      
          Res =   LO_Res_Unpol(Up_,Bot_)   * ( NNpdf(1,+2)*NNpdf(2,+5)  +  NNpdf(1,+4)*NNpdf(2,+5) )   &
                + LO_Res_Unpol(Bot_,Up_)   * ( NNpdf(1,+5)*NNpdf(2,+2)  +  NNpdf(1,+5)*NNpdf(2,+4) )   &
                + LO_Res_Unpol(ADn_,Bot_)  * ( NNpdf(1,-1)*NNpdf(2,+5)  +  NNpdf(1,-3)*NNpdf(2,+5) )   &
                + LO_Res_Unpol(Bot_,ADn_)  * ( NNpdf(1,+5)*NNpdf(2,-1)  +  NNpdf(1,+5)*NNpdf(2,-3) )
      endif
      if( Channel.eq.2 .or. Channel.eq.0 ) then
          call EvalAmp_QQB_THBBAR(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
          Res = Res &
                + LO_Res_Unpol(Up_,ADn_)   *  NNpdf(1,+2)*NNpdf(2,-1)     &
                + LO_Res_Unpol(ADn_,Up_)   *  NNpdf(1,-1)*NNpdf(2,+2)     &
                + LO_Res_Unpol(Chm_,AStr_) *  NNpdf(1,+4)*NNpdf(2,-3)     &
                + LO_Res_Unpol(AStr_,Chm_) *  NNpdf(1,-3)*NNpdf(2,+4)  
      endif
      Res = Res/x1/x2/(2d0*E_CMS**2)
            
!     restore incoming momenta (in all-outgoing convention)
      Mom(1,1:2) = -0.5d0*Collider_Energy
      Mom(4,1)   = -0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
      Mom(4,2)   = +0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
  
RETURN
END SUBROUTINE



      
SUBROUTINE EvalXSec_PP_TbH(Mom,TTBHcoupl,TopDecays,Channel,Res)
implicit none
real(8) :: Mom(1:4,1:13),Res
complex(8) :: TTBHcoupl(1:2)
integer :: Channel! 0=s+t channel, 1=t channel, 2=s channel
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: eta1,eta2,Etot,Pztot,MatElSq_GG,MatElSq_QQB,MatElSq_QBQ
real(8) :: x1,x2,PDFScale,Collider_Energy,E_CMS
real(8) :: NNpdf(1:2,-6:7),m_ferm,LO_Res_Unpol(-6:6,-6:6)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'
      

      
      Collider_Energy = -(Mom(1,inLeft)+Mom(1,inRight))
      if( TopDecays.eq.0 ) then
          Etot = Mom(1,Hbos)+Mom(1,t)+Mom(1,qout)
          Pztot= Mom(4,Hbos)+Mom(4,t)+Mom(4,qout)
      else
          Etot = Mom(1,Hbos) + Mom(1,qout) + Mom(1,b)+Mom(1,lep)+Mom(1,nu)
          Pztot= Mom(4,Hbos) + Mom(4,qout) + Mom(4,b)+Mom(4,lep)+Mom(4,nu)
      endif
      x1 = (Etot+Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      x2 = (Etot-Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      E_CMS = dsqrt(x1*x2)*Collider_Energy
      PDFScale = 0.25d0*( m_Top + m_Reso ) * 100d0

      Mom(1:4,inLeft)  = x1 * Mom(1:4,inLeft)
      Mom(1:4,inRight) = x2 * Mom(1:4,inRight)
      
      call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
      call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))

      Res = 0d0
      if( Channel.eq.1 .or. Channel.eq.0 ) then
          call EvalAmp_QbarBbar_TbarH(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
          Res = + LO_Res_Unpol(Dn_,ABot_)  * ( NNpdf(Dn_,1)*NNpdf(ABot_,2)  + NNpdf(Str_,1)*NNpdf(ABot_,2) )    &
                + LO_Res_Unpol(ABot_,Dn_)  * ( NNpdf(ABot_,1)*NNpdf(Dn_,2)  + NNpdf(ABot_,1)*NNpdf(Str_,2) )    &
                + LO_Res_Unpol(AUp_,ABot_) * ( NNpdf(AUp_,1)*NNpdf(ABot_,2) + NNpdf(AChm_,1)*NNpdf(ABot_,2))    &
                + LO_Res_Unpol(ABot_,AUp_) * ( NNpdf(ABot_,1)*NNpdf(AUp_,2) + NNpdf(ABot_,1)*NNpdf(AChm_,2))
      endif
      if( Channel.eq.2 .or. Channel.eq.0 ) then
          call EvalAmp_QQB_TBARHB(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
          Res = Res &
                + LO_Res_Unpol(AUp_,Dn_)   *  NNpdf(1,-2)*NNpdf(2,+1)     &
                + LO_Res_Unpol(Dn_,AUp_)   *  NNpdf(1,+1)*NNpdf(2,-2)     &
                + LO_Res_Unpol(AChm_,Str_) *  NNpdf(1,-4)*NNpdf(2,+3)     &
                + LO_Res_Unpol(Str_,AChm_) *  NNpdf(1,+3)*NNpdf(2,-4)  
      endif
      Res = Res/x1/x2/(2d0*E_CMS**2)

!     restore incoming momenta (in all-outgoing convention)
      Mom(1,1:2) = -0.5d0*Collider_Energy
      Mom(4,1)   = -0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
      Mom(4,2)   = +0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
  
RETURN
END SUBROUTINE



      
 
 
! t-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing
SUBROUTINE EvalAmp_QB_TH(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'



! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,qout)*MomExt(1,t)-MomExt(2,qout)*MomExt(2,t)-MomExt(3,qout)*MomExt(3,t)-MomExt(4,qout)*MomExt(4,t)
      p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
      p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)
      MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
      MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
      MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
      MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
      MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,qout)/2d0/p4Dp5
      MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
      MomExtFlat(7,1:4)=MomExt(1:4,qout)

    ! use different flattened momenta for top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
!          ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,b)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)
      ENDIF


! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp(:,:,:)=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call tdecay(5,6,8,9,10,za,zb,decay_amp)
      ENDIF
            
      call ubhtdamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Up_,Bot_,1:2))
      call ubhtdamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Bot_,Up_,1:2))        
      call ubhtdamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ADn_,Bot_,1:2))
      call ubhtdamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Bot_,ADn_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Up_,Bot_)  = cdabs(LOAmp(Up_,Bot_,1))**2  + cdabs(LOAmp(Up_,Bot_,2))**2
      LO_Res_Unpol(Bot_,Up_)  = cdabs(LOAmp(Bot_,Up_,1))**2  + cdabs(LOAmp(Bot_,Up_,2))**2
      LO_Res_Unpol(ADn_,Bot_) = cdabs(LOAmp(ADn_,Bot_,1))**2 + cdabs(LOAmp(ADn_,Bot_,2))**2
      LO_Res_Unpol(Bot_,ADn_) = cdabs(LOAmp(Bot_,ADn_,1))**2 + cdabs(LOAmp(Bot_,ADn_,2))**2
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
! t-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing   
SUBROUTINE EvalAmp_QbarBbar_TbarH(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'



! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,qout)*MomExt(1,t)-MomExt(2,qout)*MomExt(2,t)-MomExt(3,qout)*MomExt(3,t)-MomExt(4,qout)*MomExt(4,t)
   p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
   p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)
   MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
   MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
   MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
   MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
   MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,qout)/2d0/p4Dp5
   MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
   MomExtFlat(7,1:4)=MomExt(1:4,qout)

    ! use different flattened momenta for anti-top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,b)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)

      ENDIF



! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call atdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
      call dbbarhtbaruamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Dn_,ABot_,1:2))
      call dbbarhtbaruamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ABot_,Dn_,1:2))
      call dbbarhtbaruamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,TTBHcoupl,LOAmp(AUp_,ABot_,1:2))
      call dbbarhtbaruamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ABot_,AUp_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Dn_,ABot_)  = cdabs(LOAmp(Dn_,ABot_,1))**2  + cdabs(LOAmp(Dn_,ABot_,2))**2
      LO_Res_Unpol(ABot_,Dn_)  = cdabs(LOAmp(ABot_,Dn_,1))**2  + cdabs(LOAmp(ABot_,Dn_,2))**2
      LO_Res_Unpol(AUp_,ABot_) = cdabs(LOAmp(AUp_,ABot_,1))**2 + cdabs(LOAmp(AUp_,ABot_,2))**2
      LO_Res_Unpol(ABot_,AUp_) = cdabs(LOAmp(ABot_,AUp_,1))**2 + cdabs(LOAmp(ABot_,AUp_,2))**2
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
   

! s-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing   
SUBROUTINE EvalAmp_QQB_THBBAR(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, bout=5, bdk=6,W=7,lep=8,nu=9
include 'includeVars.F90'


! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,bout)*MomExt(1,t)-MomExt(2,bout)*MomExt(2,t)-MomExt(3,bout)*MomExt(3,t)-MomExt(4,bout)*MomExt(4,t)
      p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
      p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)

      MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
      MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
      MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
      MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
      MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,bout)/2d0/p4Dp5
      MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
      MomExtFlat(7,1:4)=MomExt(1:4,bout)

    ! use different flattened momenta for top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
!          ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,bdk)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)
      ENDIF


! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp(:,:,:)=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call tdecay(5,6,8,9,10,za,zb,decay_amp)
      ENDIF
      call udbar_htbbaramp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Up_,ADn_,1:2))
      call udbar_htbbaramp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ADn_,Up_,1:2))       
      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Up_,ADn_)  = cdabs(LOAmp(Up_,ADn_,1))**2  + cdabs(LOAmp(Up_,ADn_,2))**2
      LO_Res_Unpol(ADn_,Up_)  = cdabs(LOAmp(ADn_,Up_,1))**2  + cdabs(LOAmp(ADn_,Up_,2))**2
      LO_Res_Unpol(Chm_,AStr_)  = LO_Res_Unpol(Up_,ADn_) 
      LO_Res_Unpol(AStr_,Chm_)  = LO_Res_Unpol(ADn_,Up_) 
           
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2

   
RETURN
END SUBROUTINE EvalAmp_QQB_THBBAR




! s-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing
SUBROUTINE EvalAmp_QQB_TBARHB(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, bout=5, bdk=6,W=7,lep=8,nu=9
include 'includeVars.F90'


! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,bout)*MomExt(1,t)-MomExt(2,bout)*MomExt(2,t)-MomExt(3,bout)*MomExt(3,t)-MomExt(4,bout)*MomExt(4,t)
   p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
   p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)
   MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
   MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
   MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
   MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
   MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,bout)/2d0/p4Dp5
   MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
   MomExtFlat(7,1:4)=MomExt(1:4,bout)

    ! use different flattened momenta for anti-top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,bdk)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)

      ENDIF



! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call atdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
      call ubard_Htbarbamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(AUp_,Dn_,1:2))
      call ubard_Htbarbamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Dn_,AUp_,1:2))


      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Dn_,AUp_)  = cdabs(LOAmp(Dn_,AUp_,1))**2  + cdabs(LOAmp(Dn_,AUp_,2))**2
      LO_Res_Unpol(AUp_,Dn_)  = cdabs(LOAmp(AUp_,Dn_,1))**2  + cdabs(LOAmp(AUp_,Dn_,2))**2
      LO_Res_Unpol(Str_,AChm_)  = LO_Res_Unpol(Dn_,AUp_) 
      LO_Res_Unpol(AChm_,Str_)  = LO_Res_Unpol(AUp_,Dn_)       
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
   
   
   
! t-channel   
      SUBROUTINE ubhtdamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production u(p1)+b(p2)->H(p3)+t(p4)+d(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2),TTBHcoupl(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s24,s34,s15,mt,mw
        include 'includeVars.F90'

        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top

! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))        
        
        ampw(1) = 1/( - mw**2 + s24)/( - mw**2 + s15)*za(p5,e4)*zb(p1,&
     & p2)*mw**2 + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)/(zb(k4,e4)&
     & )*za(p5,e3)*zb(p2,k4)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + s24)&
     & /( - mw**2 + s15)/(zb(k4,e4))*za(p5,k3)*zb(p2,k4)*zb(k3,p1)*&
     & mt**2
        
        ampw(2) = 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)*za(p5,e3)*&
     & zb(p2,e4)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + &
     & s15)*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*mt + 1/( - mw**2 + s24)/( - &
     & mw**2 + s15)/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*mw**2*mt
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p5,e4&
     & )*zb(p1,p2)*mt*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15&
     & )/(zb(k4,e4))*za(p5,p1)*zb(p1,p2)*zb(p1,k4)*mt*vev*KL - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s15)/(zb(k4,e4))*za(p5,p2)*zb(p1,&
     & p2)*zb(p2,k4)*mt*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p5,p1&
     & )*zb(p1,p2)*zb(p1,e4)*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s15)*za(p5,p2)*zb(p1,p2)*zb(p2,e4)*vev*KL - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s15)/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*&
     & mt**2*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      RETURN
      END SUBROUTINE
        
    
! t-channel    
      SUBROUTINE dbbarhtbaruamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production d(p1)+bbar(p2)->H(p3)+tbar(p4)+u(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2),TTBHcoupl(1:2)
        real(8)    :: s24,s34,s15,mt,mw
        include 'includeVars.F90'

        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top

! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))
        
        ampw(1) = 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)*za(p2,e4)*&
     & za(p5,e3)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + &
     & s15)*za(p2,e4)*za(p5,k3)*zb(k3,p1)*mt - 1/( - mw**2 + s24)/( - &
     & mw**2 + s15)/(zb(k4,e4))*za(p2,p5)*zb(p1,k4)*mw**2*mt
        
        ampw(2) =  - 1/( - mw**2 + s24)/( - mw**2 + s15)*za(p2,p5)*zb(&
     & p1,e4)*mw**2 + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)/(za(k4,&
     & e4))*za(p2,k4)*za(p5,e3)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + &
     & s24)/( - mw**2 + s15)/(za(k4,e4))*za(p2,k4)*za(p5,k3)*zb(k3,p1)*&
     & mt**2
        
        ampt(1) = 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p2,p5)*&
     & za(e4,p2)*zb(p2,p1)*vev*KR + 1d0/2d0/( - mt**2 + s34)/( - mw**2&
     &  + s15)*za(p2,p5)*za(e4,p5)*zb(p5,p1)*vev*KR + 1d0/2d0/( - mt**2&
     &  + s34)/( - mw**2 + s15)/(zb(k4,e4))*za(p2,p5)*zb(p1,k4)*mt**2*&
     & vev*KL
        
        ampt(2) = 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p2,p5)*&
     & zb(p1,e4)*mt*vev*KL + 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)&
     & /(za(k4,e4))*za(p2,p5)*za(k4,p2)*zb(p2,p1)*mt*vev*KR + 1d0/2d0/(&
     &  - mt**2 + s34)/( - mw**2 + s15)/(za(k4,e4))*za(p2,p5)*za(k4,p5)&
     & *zb(p5,p1)*mt*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)

      RETURN        
      END SUBROUTINE 
    


    
! s-channel    
      SUBROUTINE udbar_htbbaramp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production u(p1)+dbar(p2)->H(p3)+t(p4)+bbar(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2),TTBHcoupl(1:2)
        real(8)    :: s12,s34,s45,mt,mw
        include 'includeVars.F90'    


        s45=s(k4,p5)+s(e4,p5)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s12=s(p1,p2)
        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))

        
        ampw(1) = 1/( - mw**2 + s45)/( - mw**2 + s12)*za(p2,e4)*zb(p1,&
     & p5)*mw**2 + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)/(zb(k4,e4)&
     & )*za(p2,e3)*zb(p5,k4)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + s45)&
     & /( - mw**2 + s12)/(zb(k4,e4))*za(p2,k3)*zb(p5,k4)*zb(k3,p1)*&
     & mt**2
        
        ampw(2) = 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)*za(p2,e3)*&
     & zb(p5,e4)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + &
     & s12)*za(p2,k3)*zb(p5,e4)*zb(k3,p1)*mt + 1/( - mw**2 + s45)/( - &
     & mw**2 + s12)/(za(k4,e4))*za(p2,k4)*zb(p1,p5)*mw**2*mt
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p2,e4&
     & )*zb(p1,p5)*mt*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12&
     & )/(zb(k4,e4))*za(p2,p1)*zb(p1,p5)*zb(p1,k4)*mt*vev*KL - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s12)/(zb(k4,e4))*za(p2,p5)*zb(p1,&
     & p5)*zb(p5,k4)*mt*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p2,p1&
     & )*zb(p1,p5)*zb(p1,e4)*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s12)*za(p2,p5)*zb(p1,p5)*zb(p5,e4)*vev*KL - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s12)/(za(k4,e4))*za(p2,k4)*zb(p1,p5)*&
     & mt**2*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      RETURN  
      END SUBROUTINE
        


! s-channel    
      SUBROUTINE ubard_Htbarbamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production ubar(p1)+d(p2)->H(p3)+tbar(p4)+b(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2),TTBHcoupl(1:2)
        real(8)    :: s45,s34,s12,mt,mw
        include 'includeVars.F90'    


        s45=s(k4,p5)+s(e4,p5)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s12=s(p1,p2)

        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))

       
        ampw(1) = 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)*za(p1,e3)*&
     & za(p5,e4)*zb(e3,p2)*mt + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + &
     & s12)*za(p1,k3)*za(p5,e4)*zb(k3,p2)*mt + 1/( - mw**2 + s45)/( - &
     & mw**2 + s12)/(zb(k4,e4))*za(p1,p5)*zb(p2,k4)*mw**2*mt
        
        ampw(2) = 1/( - mw**2 + s45)/( - mw**2 + s12)*za(p1,p5)*zb(p2,&
     & e4)*mw**2 + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)/(za(k4,e4)&
     & )*za(p1,e3)*za(p5,k4)*zb(e3,p2)*mt**2 + 1d0/2d0/( - mw**2 + s45)&
     & /( - mw**2 + s12)/(za(k4,e4))*za(p1,k3)*za(p5,k4)*zb(k3,p2)*&
     & mt**2
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p1,p5&
     & )*za(e4,p1)*zb(p1,p2)*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s12)*za(p1,p5)*za(e4,p5)*zb(p5,p2)*vev*KR - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s12)/(zb(k4,e4))*za(p1,p5)*zb(p2,k4)*&
     & mt**2*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p1,p5&
     & )*zb(p2,e4)*mt*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12&
     & )/(za(k4,e4))*za(p1,p5)*za(k4,p1)*zb(p1,p2)*mt*vev*KR - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s12)/(za(k4,e4))*za(p1,p5)*za(k4,&
     & p5)*zb(p5,p2)*mt*vev*KR                

        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      RETURN  
      END SUBROUTINE 
    
    
    
    
    
    
    
    SUBROUTINE TDECAY(k4,e4,b,ep,nu,za,zb,dkamp)
! top decay routine, taken from MCFM, see hep-ph:/1204.1513
       implicit none
       integer :: k4,e4,b,ep,nu
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp
       include 'includeVars.F90'
 
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = za(b,nu)*zb(ep,e4)
       dkamp(2) = m_top * za(b,nu)*zb(ep,k4)/za(e4,k4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * (4d0*dsqrt(2d0)*m_W**2*GF)        
   
     END SUBROUTINE
   
    
    SUBROUTINE ATDECAY(k4,e4,bbar,em,nubar,za,zb,dkamp)
! anti-top decay routine, taken from MCFM, see hep-ph:/1204.1513
       implicit none
       integer :: k4,e4,bbar,em,nubar
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp
       include 'includeVars.F90'
  
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = -m_top * zb(bbar,nubar)*za(em,k4)/zb(e4,k4)
       dkamp(2) = -zb(bbar,nubar)*za(em,e4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * (4d0*dsqrt(2d0)*m_W**2*GF)
   
   
     END SUBROUTINE 
   

 


    SUBROUTINE convert_to_MCFM(p,pout)
      implicit none
! converts from (E,px,py,pz) to (px,py,pz,E)
      real(8) :: p(1:4),tmp(1:4)
      real(8), optional :: pout(1:4)

      if( present(pout) ) then
          pout(1)=p(2)  
          pout(2)=p(3)  
          pout(3)=p(4) 
          pout(4)=p(1)  
      else
          tmp(1)=p(1)
          tmp(2)=p(2)
          tmp(3)=p(3)
          tmp(4)=p(4)

          p(1)=tmp(2)  
          p(2)=tmp(3) 
          p(3)=tmp(4)  
          p(4)=tmp(1)  
      endif  
      
    END SUBROUTINE




subroutine spinoru(N,p,za,zb,s)
!---Calculate spinor products      
!---taken from MCFM & modified by R. Rontsch, May 2015
!---extended to deal with negative energies ie with all momenta outgoing                                                                
!---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,                                                                                  
!---za(i,j)*zb(j,i)=s(i,j)                      
      implicit none
      real(8) :: p(:,:),two
      integer, parameter :: mxpart=14
      complex(8):: c23(N),f(N),rt(N),za(:,:),zb(:,:),czero,cone,ci
      real(8)   :: s(:,:)
      integer i,j,N
      
      if (size(p,1) .ne. N) then
         print *, "spinorz: momentum mismatch"
         stop
      endif
      two=2d0
      czero=dcmplx(0d0,0d0)
      cone=dcmplx(1d0,0d0)
      ci=dcmplx(0d0,1d0)
      

!---if one of the vectors happens to be zero this routine fails.                                                                                                                
      do j=1,N
         za(j,j)=czero
         zb(j,j)=za(j,j)

!-----positive energy case                                                                                                                                                      
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(p(j,4)+p(j,1))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
!-----negative energy case                                                                                                                                                      
            rt(j)=dsqrt(-p(j,4)-p(j,1))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=ci
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

         if (abs(s(i,j)).lt.1d-5) then
         zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
         else
         zb(i,j)=-dcmplx(s(i,j))/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

    end subroutine spinoru

    
    
    
END MODULE

