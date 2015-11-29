MODULE ModTHiggs
use ModParameters
implicit none


public :: EvalAmp_QB_TH, EvalAmp_QbarBbar_TH     ! t-channel
public :: EvalAmp_QQB_THBBAR,EvalAmp_QQB_TBARHB  ! s-channel
private

 CONTAINS 


 
   
SUBROUTINE EvalAmp_QB_TH(MomExt,LO_Res_Unpol)
use ModMisc
use modParameters
implicit none
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9



! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,qout)*MomExt(1,t)-MomExt(2,qout)*MomExt(2,t)-MomExt(3,qout)*MomExt(3,t)-MomExt(4,qout)*MomExt(4,t)
      p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
      p2Dp3=MomExt(1,inRight)*MomExt(1,Hbos)-MomExt(2,inRight)*MomExt(2,Hbos)-MomExt(3,inRight)*MomExt(3,Hbos)-MomExt(4,inRight)*MomExt(4,Hbos)
      MomExtFlat(1,1:4)=MomExt(1:4,inLeft)
      MomExtFlat(2,1:4)=MomExt(1:4,inRight)
      MomExtFlat(3,1:4)=m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
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
      
      call ubhtdamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Up_,Bot_,1:2))
      call ubhtdamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Bot_,Up_,1:2))        
      call ubhtdamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,LOAmp(ADn_,Bot_,1:2))
      call ubhtdamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,LOAmp(Bot_,ADn_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci
      LO_Res_Unpol(Up_,Bot_)  = cdabs(LOAmp(Up_,Bot_,1))**2  + cdabs(LOAmp(Up_,Bot_,2))**2
      
! print *, "raoul",LO_Res_Unpol(Up_,Bot_)
! call myAmp((/MomExt(1:4,4),-MomExt(1:4,1),MomExt(1:4,5),-MomExt(1:4,2)/),ghz1,ghz2,ghz3,kappa,kappa_tilde,LO_Res_Unpol(Chm_,Bot_))
! print *, "marku",LO_Res_Unpol(Chm_,Bot_)
! print *, "ratio",LO_Res_Unpol(Up_,Bot_)/LO_Res_Unpol(Chm_,Bot_)
! pause
      
      LO_Res_Unpol(Bot_,Up_)  = cdabs(LOAmp(Bot_,Up_,1))**2  + cdabs(LOAmp(Bot_,Up_,2))**2
      LO_Res_Unpol(Chm_,Bot_) = LO_Res_Unpol(Up_,Bot_)
      LO_Res_Unpol(Bot_,Chm_) = LO_Res_Unpol(Bot_,Up_)
      LO_Res_Unpol(ADn_,Bot_) = cdabs(LOAmp(ADn_,Bot_,1))**2 + cdabs(LOAmp(ADn_,Bot_,2))**2
      LO_Res_Unpol(Bot_,ADn_) = cdabs(LOAmp(Bot_,ADn_,1))**2 + cdabs(LOAmp(Bot_,ADn_,2))**2
      LO_Res_Unpol(AStr_,Bot_)= LO_Res_Unpol(ADn_,Bot_)
      LO_Res_Unpol(Bot_,AStr_)= LO_Res_Unpol(Bot_,ADn_)
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
   

SUBROUTINE EvalAmp_QbarBbar_TH(MomExt,LO_Res_Unpol)
use ModMisc
use modParameters
implicit none
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9



! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,qout)*MomExt(1,t)-MomExt(2,qout)*MomExt(2,t)-MomExt(3,qout)*MomExt(3,t)-MomExt(4,qout)*MomExt(4,t)
   p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
   p2Dp3=MomExt(1,inRight)*MomExt(1,Hbos)-MomExt(2,inRight)*MomExt(2,Hbos)-MomExt(3,inRight)*MomExt(3,Hbos)-MomExt(4,inRight)*MomExt(4,Hbos)
   MomExtFlat(1,1:4)=MomExt(1:4,inLeft)
   MomExtFlat(2,1:4)=MomExt(1:4,inRight)
   MomExtFlat(3,1:4)=m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
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
      call dbbarhtbaruamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Dn_,ABot_,1:2))
      call dbbarhtbaruamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(ABot_,Dn_,1:2))
      call dbbarhtbaruamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,LOAmp(AUp_,ABot_,1:2))
      call dbbarhtbaruamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,LOAmp(ABot_,AUp_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Dn_,ABot_)  = cdabs(LOAmp(Dn_,ABot_,1))**2  + cdabs(LOAmp(Dn_,ABot_,2))**2
      LO_Res_Unpol(ABot_,Dn_)  = cdabs(LOAmp(ABot_,Dn_,1))**2  + cdabs(LOAmp(ABot_,Dn_,2))**2
      LO_Res_Unpol(Str_,ABot_) = LO_Res_Unpol(Dn_,ABot_)
      LO_Res_Unpol(ABot_,Str_) = LO_Res_Unpol(ABot_,Dn_)
      
      LO_Res_Unpol(AUp_,ABot_) = cdabs(LOAmp(AUp_,ABot_,1))**2 + cdabs(LOAmp(AUp_,ABot_,2))**2
      LO_Res_Unpol(ABot_,AUp_) = cdabs(LOAmp(ABot_,AUp_,1))**2 + cdabs(LOAmp(ABot_,AUp_,2))**2
      LO_Res_Unpol(AChm_,ABot_)= LO_Res_Unpol(AUp_,ABot_)
      LO_Res_Unpol(ABot_,AChm_)= LO_Res_Unpol(ABot_,AUp_)
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
 
   
SUBROUTINE EvalAmp_QQB_THBBAR(MomExt,LO_Res_Unpol)
use ModMisc
use modParameters
implicit none
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, bout=5, bdk=6,W=7,lep=8,nu=9


! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,bout)*MomExt(1,t)-MomExt(2,bout)*MomExt(2,t)-MomExt(3,bout)*MomExt(3,t)-MomExt(4,bout)*MomExt(4,t)
      p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
      p2Dp3=MomExt(1,inRight)*MomExt(1,Hbos)-MomExt(2,inRight)*MomExt(2,Hbos)-MomExt(3,inRight)*MomExt(3,Hbos)-MomExt(4,inRight)*MomExt(4,Hbos)

      MomExtFlat(1,1:4)=MomExt(1:4,inLeft)
      MomExtFlat(2,1:4)=MomExt(1:4,inRight)
      MomExtFlat(3,1:4)=m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
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
      call udbar_htbbaramp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Up_,ADn_,1:2))
      call udbar_htbbaramp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(ADn_,Up_,1:2))       
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


SUBROUTINE EvalAmp_QQB_TBARHB(MomExt,LO_Res_Unpol)
use ModMisc
use modParameters
implicit none
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, bout=5, bdk=6,W=7,lep=8,nu=9



! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,bout)*MomExt(1,t)-MomExt(2,bout)*MomExt(2,t)-MomExt(3,bout)*MomExt(3,t)-MomExt(4,bout)*MomExt(4,t)
   p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
   p2Dp3=MomExt(1,inRight)*MomExt(1,Hbos)-MomExt(2,inRight)*MomExt(2,Hbos)-MomExt(3,inRight)*MomExt(3,Hbos)-MomExt(4,inRight)*MomExt(4,Hbos)
   MomExtFlat(1,1:4)=MomExt(1:4,inLeft)
   MomExtFlat(2,1:4)=MomExt(1:4,inRight)
   MomExtFlat(3,1:4)=m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
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
      call ubard_Htbarbamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(AUp_,Dn_,1:2))
      call ubard_Htbarbamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Dn_,AUp_,1:2))


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
   
      
      
   
   
   
      SUBROUTINE ubhtdamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production u(p1)+b(p2)->H(p3)+t(p4)+d(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        use ModParameters
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s24,s34,s15,mt,mw

        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)        
                
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
        
      END SUBROUTINE
        
    
      SUBROUTINE dbbarhtbaruamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production d(p1)+bbar(p2)->H(p3)+tbar(p4)+u(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        use ModParameters
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s24,s34,s15,mt,mw
    
        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)
        
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
        
      END SUBROUTINE 
      
      
      
      subroutine udbar_htbbaramp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production u(p1)+dbar(p2)->H(p3)+t(p4)+bbar(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        use ModParameters
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s12,s34,s45,mt,mw
    
        s45=s(k4,p5)+s(e4,p5)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s12=s(p1,p2)
        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)

        
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
        
      end subroutine udbar_htbbaramp
        
    
      subroutine ubard_Htbarbamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production ubar(p1)+d(p2)->H(p3)+tbar(p4)+b(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        use ModParameters
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:),KL,KR
        complex(8) :: amp(2),ampw(2),ampt(2)
        real(8)    :: s45,s34,s12,mt,mw
    
        s45=s(k4,p5)+s(e4,p5)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s12=s(p1,p2)

        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)

       
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
        
      end subroutine ubard_Htbarbamp
    
    
    
    
    SUBROUTINE TDECAY(k4,e4,b,ep,nu,za,zb,dkamp)
! top decay routine, taken from MCFM, see hep-ph:/1204.1513
       use ModParameters
       implicit none
       integer :: k4,e4,b,ep,nu
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp
     
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = za(b,nu)*zb(ep,e4)
       dkamp(2) = m_top * za(b,nu)*zb(ep,k4)/za(e4,k4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * g2_weak
   
   
     END SUBROUTINE
   
    
    SUBROUTINE ATDECAY(k4,e4,bbar,em,nubar,za,zb,dkamp)
! anti-top decay routine, taken from MCFM, see hep-ph:/1204.1513
       use ModParameters
       implicit none
       integer :: k4,e4,bbar,em,nubar
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       real(8),parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
       complex(8) :: WProp
     
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = -m_top * zb(bbar,nubar)*za(em,k4)/zb(e4,k4)
       dkamp(2) = -zb(bbar,nubar)*za(em,e4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * g2_weak
   
   
     END SUBROUTINE 
   


     SUBROUTINE myAmp(p,VVHg1,VVHg2,VVHg3,TTHg1,TTHg2,res)
     implicit none
     real(8) :: p(1:4,1:4),s12,s34,s234,IV1,IV2,sprod(1:4,1:4),res
     complex(8) :: VVHg1,VVHg2,VVHg3,TTHg1,TTHg2
     complex(8) :: HelAmp(1:2),Props1,Props2
     complex(8) :: za(1:4,1:4),zb(1:4,1:4)
     complex(8),parameter :: cI=(0d0,1d0)
     
        call my_spinoru(4,(/p(1:4,1),p(1:4,2),p(1:4,3),p(1:4,4)/),za,zb,sprod)
     
        s12 = m_Top**2 + sprod(1,2)
        s34 = sprod(3,4)
        Props1 =  cI/(s12-M_W**2+cI*M_W*Ga_W) * cI/(s34-M_W**2+cI*M_W*Ga_W)           ! all propagator denominators for the HVV diagram
        s234 = sprod(2,3)+sprod(2,4)+sprod(3,4)
        Props2 =  cI/(s234-M_Top**2+cI*M_Top*Ga_Top) * cI/(s34-M_W**2+cI*M_W*Ga_W)    ! all propagator denominators for the Yukawa diagram
     
        IV1 = 1d0/dsqrt(2d0)/sitW
        IV2 = 1d0/dsqrt(2d0)/sitW
     
        HelAmp(1) =  2*Props2*(TTHg1 - cI*TTHg2)*IV1*IV2*zb(2,4)*(za(2,3)*zb(2,1) - za(3,4)*zb(4,1)) +      &
       (m_top*Props1*IV1*IV2*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1))*(VVHg2*za(3,1)**2*zb(2,1)*zb(3,1)**2 + VVHg2*za(2,3)*zb(2,1)*zb(2,3)*(-2*m_top**2 + za(3,1)*zb(3,1)) +       &
            za(3,4)*(M_W**2*(VVHg2 - cI*VVHg3)*zb(2,4)*zb(3,1) - m_top**2*VVHg2*zb(2,1)*zb(3,4)) +       &
            za(3,1)*zb(2,1)*zb(3,1)*(-VVHg1 - 2*m_top**2*VVHg2 + M_W**2*VVHg2 + 2*VVHg2*za(2,1)*zb(2,1) + VVHg2*za(2,4)*zb(2,4) +       &
               VVHg2*za(4,1)*zb(4,1))))/(M_W**2*za(3,1)*zb(3,1))    

         HelAmp(2) =   (-2*m_top*Props2*IV1*IV2*zb(2,4)*((-TTHg1 + cI*TTHg2)*za(2,3)*zb(2,3) + (TTHg1 + cI*TTHg2)*za(3,1)*zb(3,1) +       &
            (-TTHg1 + cI*TTHg2)*za(3,4)*zb(3,4)))/zb(3,1) +       &
            (Props1*IV1*IV2*(m_top**2*VVHg2*za(2,3)**2*zb(2,3)**2*zb(2,4)*(-2*m_top**2 + za(3,1)*zb(3,1)) +       &
            za(2,3)*(zb(2,4)*(za(3,1)*zb(3,1)*(m_top**2*zb(2,3)*(-VVHg1 - 2*m_top**2*VVHg2 + M_W**2*VVHg2 + 2*VVHg2*za(2,1)*zb(2,1) + VVHg2*za(2,4)*zb(2,4)) +       &
                     VVHg2*((m_top - M_W)*(m_top + M_W)*za(3,1)*zb(2,3) - M_W**2*za(4,1)*zb(2,4))*zb(3,1)) - m_top**4*VVHg2*za(3,4)*zb(2,3)*zb(3,4)) +       &
               m_top**2*VVHg2*za(3,1)*zb(2,3)*(za(4,1)*zb(2,4)*zb(3,1) + zb(2,3)*(-2*m_top**2 + za(3,1)*zb(3,1)))*zb(4,1)) +       &
            za(3,1)*(M_W**2*zb(2,4)*zb(3,1)*(2*VVHg1*za(3,1)*zb(3,1) + cI*VVHg3*za(3,1)**2*zb(3,1)**2 +       &
                  cI*VVHg3*za(3,4)*(-(za(2,1)*zb(2,4)*zb(3,1)) + m_top**2*zb(3,4))) -       &
               (za(3,1)*zb(3,1)*(m_top**2*zb(2,3)*(VVHg1 + 2*m_top**2*VVHg2 - M_W**2*VVHg2 - 2*VVHg2*za(2,1)*zb(2,1) - VVHg2*za(2,4)*zb(2,4)) +       &
                     ((-(m_top**2*VVHg2) + M_W**2*(VVHg2 + cI*VVHg3))*za(3,1)*zb(2,3) + M_W**2*VVHg2*za(4,1)*zb(2,4))*zb(3,1)) +       &
                  m_top**4*VVHg2*za(3,4)*zb(2,3)*zb(3,4))*zb(4,1) + m_top**2*VVHg2*za(3,1)*za(4,1)*zb(2,3)*zb(3,1)*zb(4,1)**2)))/(M_W**2*za(3,1)*zb(3,1)**2)
     
          Res =  HelAmp(1)*dconjg(HelAmp(1)) + HelAmp(2)*dconjg(HelAmp(2)) 
     
     END SUBROUTINE 

  subroutine my_spinoru(n,p,za,zb,s)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: p(4,n)
    complex(8), intent(out) :: za(n,n), zb(n,n)
    real(8), intent(out) :: s(n,n)
    integer :: i,j
    complex(8) :: c23(n), f(n)
    real(8) :: rt(n)
      
    !---if one of the vectors happens to be zero this routine fails.
    do j=1,N
       za(j,j)=czero
       zb(j,j)=za(j,j)

       !-----positive energy case
       if (p(1,j) .gt. zero) then
          rt(j)=dsqrt(dabs(p(2,j)+p(1,j)))
          c23(j)=dcmplx(p(4,j),-p(3,j))
          f(j)=(one,zero)
       else
       !-----negative energy case
          rt(j)=dsqrt(dabs(-p(1,j)-p(2,j)))
          c23(j)=dcmplx(-p(4,j),p(3,j))
          f(j)=ci
       endif
    enddo

    do i=2,N
  
     do j=1,i-1
          s(i,j)=two*(p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)-p(4,i)*p(4,j))
          za(i,j)=f(i)*f(j)  * ( c23(i)*dcmplx(rt(j)/(rt(i)+1d-16))-c23(j)*dcmplx(rt(i)/(rt(j)+1d-16)) )
          
          if (dabs(s(i,j)).lt.1d-5) then
             zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
          else
             zb(i,j)=-dcmplx(s(i,j))/(za(i,j)+1d-16)
          endif
          
          za(j,i)=-za(i,j)
          zb(j,i)=-zb(i,j)
          s(j,i)=s(i,j)
          
       enddo

    enddo

    return
    
  end subroutine 


END MODULE 