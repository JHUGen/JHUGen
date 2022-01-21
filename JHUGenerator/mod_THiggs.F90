MODULE ModTHiggs
use ModParameters
implicit none


public :: EvalXsec_pp_TH, EvalXsec_pp_THW
public :: EvalAmp_QB_TH, EvalAmp_QbarBbar_TH     ! t-channel
public :: EvalAmp_QQB_THBBAR,EvalAmp_QQB_TBARHB  ! s-channel
public :: EvalAmp_GB_TWMH,EvalAmp_GBB_TBWPH      ! tw-channel
private

 CONTAINS


subroutine EvalXsec_pp_TH(Mom,Channel,Res)
implicit none
real(8), intent(in) :: Mom(1:4,1:9)
integer, intent(in) :: Channel! 0=s+t channel, 1=t channel, 2=s channel
real(8), intent(out) :: Res(-5:5,-5:5)
real(8) :: MatElSq_GG,MatElSq_QQB,MatElSq_QBQ
real(8) :: LO_Res_Unpol(-6:6,-6:6)
integer :: i,j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9


   Res(:,:) = 0d0
   if( Channel.eq.1 .or. Channel.eq.0 ) then
      LO_Res_Unpol(:,:)=0d0
      call EvalAmp_QB_TH(Mom,LO_Res_Unpol) ! Computes u-b, b-u, db-b, b-db initial states
      do i=-6,6; do j=-6,6
         if(abs(i).ne.abs(Top_) .and. abs(j).ne.abs(Top_)) then
            Res(i,j) = Res(i,j) + LO_Res_Unpol(convertToPartIndex(i),convertToPartIndex(j))
         endif
      enddo; enddo

      LO_Res_Unpol(:,:)=0d0
      call EvalAmp_QbarBbar_TH(Mom,LO_Res_Unpol) ! Computes ub-bb, bb-ub, d-bb, bb-d initial states
      do i=-6,6; do j=-6,6
         if(abs(i).ne.abs(Top_) .and. abs(j).ne.abs(Top_)) then
            Res(i,j) = Res(i,j) + LO_Res_Unpol(convertToPartIndex(i),convertToPartIndex(j))
         endif
      enddo; enddo
   endif
   if( Channel.eq.2 .or. Channel.eq.0 ) then
      LO_Res_Unpol(:,:)=0d0
      call EvalAmp_QQB_THBBAR(Mom,LO_Res_Unpol) ! Computes u-db, db-u, c-sb, sb-c initial states
      do i=-6,6; do j=-6,6
         if(abs(i).ne.abs(Top_) .and. abs(j).ne.abs(Top_)) then
            Res(i,j) = Res(i,j) + LO_Res_Unpol(convertToPartIndex(i),convertToPartIndex(j))
         endif
      enddo; enddo

      LO_Res_Unpol(:,:)=0d0
      call EvalAmp_QQB_TBARHB(Mom,LO_Res_Unpol) ! Computes ub-d, d-ub, cb-s, s-cb initial states
      do i=-6,6; do j=-6,6
         if(abs(i).ne.abs(Top_) .and. abs(j).ne.abs(Top_)) then
            Res(i,j) = Res(i,j) + LO_Res_Unpol(convertToPartIndex(i),convertToPartIndex(j))
         endif
      enddo; enddo
   endif
   return
end subroutine

subroutine EvalXsec_pp_THW(Mom,Res)
implicit none
real(8), intent(in) :: Mom(1:4,1:11) !additional Mom(1:4,10:11) are lep and nu from W-decay in tw channel
real(8), intent(out) :: Res(-5:5,-5:5)
real(8) :: LO_Res_Unpol(-6:6,-6:6)
integer :: i,j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9, lepW=10, nuW=11

      LO_Res_Unpol(:,:)=0d0
      call EvalAmp_GB_TWMH(Mom,LO_Res_Unpol) ! Computes g-b, b-g initial states
      do i=-6,6; do j=-6,6
         if(abs(i).ne.abs(Top_) .and. abs(j).ne.abs(Top_)) then
            Res(i,j) = Res(i,j) + LO_Res_Unpol(convertToPartIndex(i),convertToPartIndex(j))
         endif
      enddo; enddo

      LO_Res_Unpol(:,:)=0d0
      call EvalAmp_GBB_TBWPH(Mom,LO_Res_Unpol) ! Computes g-bb, bb-g initial states
      do i=-6,6; do j=-6,6
         if(abs(i).ne.abs(Top_) .and. abs(j).ne.abs(Top_)) then
            Res(i,j) = Res(i,j) + LO_Res_Unpol(convertToPartIndex(i),convertToPartIndex(j))
         endif
      enddo; enddo
      
   return
end subroutine


SUBROUTINE EvalAmp_QB_TH(MomExt,LO_Res_Unpol)
use ModMisc
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

!          call spinoru(7,MomExtFlat,za,zb,s)
         call spinoru(MomExtFlat,za,zb,s)

      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

!          call spinoru(10,MomExtFlatDK,za,zb,s)
         call spinoru(MomExtFlatDK,za,zb,s)
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
      LO_Res_Unpol(Bot_,Up_)  = cdabs(LOAmp(Bot_,Up_,1))**2  + cdabs(LOAmp(Bot_,Up_,2))**2


! print *, "coupl",ghz1,ghz2,ghz3,kappa,kappa_tilde
! call myAmp((/MomExt(1:4,t),-MomExt(1:4,inLeft),MomExt(1:4,qout),-MomExt(1:4,inRight)/),   M_W**2/vev*ghz1,    +2d0/vev*ghz2,    +2d0/vev*ghz3,  &
!                                                                                           kappa*(-m_top/vev),   kappa_tilde*(-m_top/vev)    ,LO_Res_Unpol(Bot_,Chm_))
! print *, "marku",LO_Res_Unpol(Bot_,Chm_)
! print *, "raoul",LO_Res_Unpol(Bot_,Up_)
! print *, "ratio",LO_Res_Unpol(Bot_,Up_)/LO_Res_Unpol(Bot_,Chm_)
! pause


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

!          call spinoru(7,MomExtFlat,za,zb,s)
         call spinoru(MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

!          call spinoru(10,MomExtFlatDK,za,zb,s)
         call spinoru(MomExtFlatDK,za,zb,s)
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

!          call spinoru(7,MomExtFlat,za,zb,s)
         call spinoru(MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

!          call spinoru(10,MomExtFlatDK,za,zb,s)
         call spinoru(MomExtFlatDK,za,zb,s)
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

!          call spinoru(7,MomExtFlat,za,zb,s)
         call spinoru(MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

!          call spinoru(10,MomExtFlatDK,za,zb,s)
         call spinoru(MomExtFlatDK,za,zb,s)
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

!amplitude squared for partonic reaction 
!
! g(iL) + b(iR) -> t(t) + H(H) + Wm(w) 
!              (-> b(bt) + Wp(wt) + Wm(w)
!               -> b(bt) + f1~(lept) + f'1(nut) + f'2(lepW) + f2~(nuW))
!
!with the fermion-antifermion pairs f'1/f1~ and f'2/f2~ 
!from leptonic or hadronic decays of the Wp/Wm   
!
SUBROUTINE EvalAmp_GB_TWMH(MomExt,LO_Res_Unpol)
use ModMisc
use ModParameters
use ModTopDecay

implicit none
real(8) :: MomExt(1:4,1:11),LO_Res_UnPol(-6:6,-6:6)
complex(8) :: LOAmpab=0d0, LOAmpba=0d0

complex(8) :: UBbD(4),UBbD2(4),Eg(4),Eg2(4),UBtD(4),ECw(4)
integer, parameter :: iL=1,iR=2, H=3, t=4, w=5, bt=6,Wt=7,lept=8,nut=9, lepW=10,nuW=11
integer :: i,j,bhel,thel,gpol,wpol,theld=-1,wpold=-1

      
      LO_Res_UnPol=0d0
      
      do bhel=-1,1,2
        call ubarSpi_Dirac(dcmplx(MomExt(:,iR)),0d0,bhel,UBbD) 
        call ubarSpi_Dirac(dcmplx(MomExt(:,iL)),0d0,bhel,UBbD2)  
        
        
        do gpol=-1,1,2  

          Eg=pol_mless(dcmplx(MomExt(:,iL)),gpol,.false.)
          Eg2=pol_mless(dcmplx(MomExt(:,iR)),gpol,.false.)
                          
                                
          if( TOPDECAYS .EQ. 0 ) then 
            theld=+1
            wpold=+1
          endif     
          do wpol=-1,wpold,1
            if( TOPDECAYS .EQ. 0 ) then 
              ECw=pol_mass(dcmplx(MomExt(:,w)),wpol,.true.)                
            elseif( TOPDECAYS .NE. 0 ) then
              call WDecay(Wm_,(/MomExt(1:4,lepW),MomExt(1:4,nuW)/),ECw)
              ECw=ECw*(0d0,-1d0)/sqrt2      
            endif          
            do thel=-1,theld,2   
               
              if( TOPDECAYS .EQ. 0 ) then 
                call ubarSpi_Dirac(dcmplx(MomExt(:,t)),M_Top,thel,UBtD) 
              elseif( TOPDECAYS .NE. 0 ) then
                call TopDecay(Top_,(/MomExt(1:4,bt),MomExt(1:4,lept),MomExt(1:4,nut)/),UBtD)
                UBtD=UBtD*(0d0,-1d0)/sqrt2
              endif  
                  


              call gb_twmHamp(MomExt(1:4,iL),Eg,MomExt(1:4,iR),UBbD,MomExt(1:4,t),UBtD,MomExt(1:4,w),ECw,MomExt(1:4,H),LOAmpab)
                        
              LO_Res_Unpol(0,6) = LO_Res_Unpol(0,6) + cdabs(LOAmpab)**2

                
                
              call gb_twmHamp(MomExt(1:4,iR),Eg2,MomExt(1:4,iL),UBbD2,MomExt(1:4,t),UBtD,MomExt(1:4,w),ECw,MomExt(1:4,H),LOAmpba)


              LO_Res_Unpol(6,0) = LO_Res_Unpol(6,0) + cdabs(LOAmpba)**2
                
            enddo
          enddo 
        enddo
        

      enddo

        LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:)*1d0/2d0*1d0/2d0*1d0/8d0*(4d0/3d0) !* ColFac * SpinAvg * QuarkColAvg**2
      

RETURN
END SUBROUTINE EvalAmp_GB_TWMH


!amplitude squared for partonic reaction 
!
! g(iL) + b~(iR) -> t~(t) + H(H) + Wp(w) 
!               (-> b~(bt) + Wm(wt) + Wp(w)
!                -> b~(bt) + f'1(lept) + f1~(nut) + f2~(lepW) + f2(nuW))
!
!with the fermion-antifermion pairs f'1/f1~ and f'2/f2~ 
!from leptonic or hadronic decays of the Wm/Wp   
!
SUBROUTINE EvalAmp_GBB_TBWPH(MomExt,LO_Res_Unpol)
use ModMisc
use ModParameters
use ModTopDecay

implicit none
real(8) :: MomExt(1:4,1:11),LO_Res_UnPol(-6:6,-6:6)
complex(8) :: LOAmpab=0d0, LOAmpba=0d0

complex(8) :: VbD(4),Eg(4),VbD2(4),Eg2(4),VtD(4),ECw(4)
integer, parameter :: iL=1,iR=2, H=3, t=4, w=5, bt=6,Wt=7,lept=8,nut=9, lepW=10,nuW=11
integer :: i,j,bhel,thel,gpol,wpol,theld=-1,wpold=-1


      LO_Res_UnPol=0d0
      
      do bhel=-1,1,2
        call vSpi_Dirac(dcmplx(MomExt(:,iR)),0d0,bhel,VbD)       
        call vSpi_Dirac(dcmplx(MomExt(:,iL)),0d0,bhel,VbD2)  
                                                                           
        do gpol=-1,1,2  

          Eg=pol_mless(dcmplx(MomExt(:,iL)),gpol,.false.)
          Eg2=pol_mless(dcmplx(MomExt(:,iR)),gpol,.false.)
                        
                                 
          if( TOPDECAYS .EQ. 0 ) then  
            theld=+1
            wpold=+1
          endif    
          do wpol=-1,wpold,1         

            if( TOPDECAYS .EQ. 0 ) then  
              ECw=pol_mass(dcmplx(MomExt(:,w)),wpol,.true.) 
            elseif( TOPDECAYS .NE. 0 ) then
              call WDecay(Wp_,(/MomExt(1:4,lepW),MomExt(1:4,nuW)/),ECw)               
              ECw=ECw*(0d0,-1d0)/sqrt2       
            endif  
                      
            do thel=-1,theld,2   
    
              if( TOPDECAYS .EQ. 0 ) then  
                call vSpi_Dirac(dcmplx(MomExt(:,t)),M_Top,thel,VtD) 
              elseif( TOPDECAYS .NE. 0 ) then
                call TopDecay(ATop_,(/MomExt(1:4,bt),MomExt(1:4,lept),MomExt(1:4,nut)/),VtD)
                VtD=VtD*(0d0,-1d0)/sqrt2  
              endif  

              call gbb_tbwpHamp(MomExt(1:4,iL),Eg,MomExt(1:4,iR),VbD,MomExt(1:4,t),VtD,MomExt(1:4,w),ECw,MomExt(1:4,H),LOampab)
              
              LO_Res_Unpol(0,-6) = LO_Res_Unpol(0,-6) + cdabs(LOAmpab)**2
                
                
              call gbb_tbwpHamp(MomExt(1:4,iR),Eg2,MomExt(1:4,iL),VbD2,MomExt(1:4,t),VtD,MomExt(1:4,w),ECw,MomExt(1:4,H),LOampba)

              LO_Res_Unpol(-6,0) = LO_Res_Unpol(-6,0) + cdabs(LOAmpba)**2
                
            enddo
          enddo 
        enddo
      enddo

      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:)*1d0/2d0*1d0/2d0*1d0/8d0*(4d0/3d0) !* ColFac * SpinAvg * QuarkColAvg**2

RETURN
END SUBROUTINE EvalAmp_GBB_TBWPH



!     normalization in FORM +2*i_/vev*(om1*mw^2*d_(ro,si)+om2*p3(si)*p3(ro)+om3*e_(ro,si,al,be)*p24(al)*p15(be));
      subroutine ubhtdamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),om1,om2,om3,KL,KR
        real(8)    :: s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,s125,mt,mw

        s12=s(p1,p2)
        s13=s(p1,e3)+s(p1,k3)+s(e3,k3)
        s14=s(p1,k4)+s(p1,e4)+s(e4,k4)
        s15=s(p1,p5)
        s23=s(p2,e3)+s(p2,k3)+s(e3,k3)
        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s25=s(p2,p5)
        s45=s(p5,k4)+s(p5,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s35=s(e3,p5)+s(k3,p5)+s(e3,k3)
        s125=s(p1,p2)+s(p1,p5)+s(p2,p5)
        mw=M_W
        mt=m_Top

        KL=-mt/vev*(kappa-(0d0,1d0)*kappa_tilde)
        KR=-mt/vev*(kappa+(0d0,1d0)*kappa_tilde)

         om1 = ghz1/2d0  ! this should be proper a1,a2,a3
         om2 = ghz2
         om3 = ghz4    * (-1d0) ! correct buggy minus sign in Raoul's code


        ampw(1) =  + om3 * (  - 1d0/2d0*za(p5,p2)*za(e4,p1)*zb(p1,p2)*&
     &    zb(p2,p1)*ci - 1d0/2d0*za(p5,p2)*za(e4,p5)*zb(p2,p1)*zb(p5,p2&
     &    )*ci + za(p5,k4)*za(e4,p1)*zb(p1,p2)*zb(p1,k4)*ci + 1d0/2d0*&
     &    za(p5,k4)*za(e4,p1)*zb(p1,p2)*zb(k4,p1)*ci + za(p5,k4)*za(e4,&
     &    p5)*zb(p1,p2)*zb(p5,k4)*ci + 1d0/2d0*za(p5,k4)*za(e4,p5)*zb(&
     &    p5,p2)*zb(k4,p1)*ci + za(p5,e4)*za(e4,p1)*zb(p1,p2)*zb(p1,e4)&
     &    *ci + 1d0/2d0*za(p5,e4)*za(e4,p1)*zb(p1,p2)*zb(e4,p1)*ci + &
     &    za(p5,e4)*za(e4,p5)*zb(p1,p2)*zb(p5,e4)*ci + 1d0/2d0*za(p5,e4&
     &    )*za(e4,p5)*zb(p5,p2)*zb(e4,p1)*ci + za(p5,e4)*zb(p1,p2)*&
     &    mt**2*ci - 1d0/2d0*za(p5,e4)*zb(p1,p2)*s12*ci - 1d0/2d0*za(p5&
     &    ,e4)*zb(p1,p2)*s14*ci - 1d0/2d0*za(p5,e4)*zb(p1,p2)*s45*ci - &
     &    1d0/2d0*za(p5,e4)*zb(p1,p2)*s25*ci )
      ampw(1) = ampw(1) + om2 * (  - 1d0/2d0*za(p5,e3)*za(e4,e3)*zb(e3,&
     &    p1)*zb(e3,p2) - 1d0/2d0*za(p5,e3)*za(e4,k3)*zb(e3,p1)*zb(k3,&
     &    p2) - 1d0/2d0*za(p5,k3)*za(e4,e3)*zb(e3,p2)*zb(k3,p1) - 1d0/2d0&
      &    *za(p5,k3)*za(e4,k3)*zb(k3,p1)*zb(k3,p2) - 1d0/2d0/(zb(k4,e4&
     &    ))*za(p5,e3)*zb(p2,k4)*zb(e3,p1)*mw**(-2)*mt**4 + 1d0/2d0/(&
     &    zb(k4,e4))*za(p5,e3)*zb(p2,k4)*zb(e3,p1)*s24*mw**(-2)*mt**2&
     &     + 1d0/4d0/(zb(k4,e4))*za(p5,e3)*zb(p2,k4)*zb(e3,p1)*s12*&
     &    mw**(-2)*mt**2 + 1d0/4d0/(zb(k4,e4))*za(p5,e3)*zb(p2,k4)*zb(&
     &    e3,p1)*s14*mw**(-2)*mt**2 + 1d0/4d0/(zb(k4,e4))*za(p5,e3)*zb(&
     &    p2,k4)*zb(e3,p1)*s45*mw**(-2)*mt**2 + 1d0/4d0/(zb(k4,e4))*za(&
     &    p5,e3)*zb(p2,k4)*zb(e3,p1)*s25*mw**(-2)*mt**2 - 1d0/2d0/(zb(&
     &    k4,e4))*za(p5,k3)*zb(p2,k4)*zb(k3,p1)*mw**(-2)*mt**4 + 1d0/2d0&
     &    /(zb(k4,e4))*za(p5,k3)*zb(p2,k4)*zb(k3,p1)*s24*mw**(-2)*mt**2&
     &     + 1d0/4d0/(zb(k4,e4))*za(p5,k3)*zb(p2,k4)*zb(k3,p1)*s12*&
     &    mw**(-2)*mt**2 + 1d0/4d0/(zb(k4,e4))*za(p5,k3)*zb(p2,k4)*zb(&
     &    k3,p1)*s14*mw**(-2)*mt**2 )
      ampw(1) = ampw(1) + om2 * ( 1d0/4d0/(zb(k4,e4))*za(p5,k3)*zb(p2,&
     &    k4)*zb(k3,p1)*s45*mw**(-2)*mt**2 + 1d0/4d0/(zb(k4,e4))*za(p5,&
     &    k3)*zb(p2,k4)*zb(k3,p1)*s25*mw**(-2)*mt**2 )
      ampw(1) = ampw(1) + om1 * ( za(p5,e4)*zb(p1,p2)*mw**2 - 1d0/2d0/(&
     &    zb(k4,e4))*za(p5,p2)*zb(p2,p1)*zb(p2,k4)*mt**2 - 1d0/2d0/(zb(&
     &    k4,e4))*za(p5,k4)*zb(p2,k4)*zb(k4,p1)*mt**2 - 1d0/2d0/(zb(k4,&
     &    e4))*za(p5,e4)*zb(p2,k4)*zb(e4,p1)*mt**2 )

        ampw(2) =  + om3 * (  - 1d0/2d0/(za(k4,e4))*za(p5,p2)*za(k4,p1)&
     &    *zb(p1,p2)*zb(p2,p1)*mt*ci - 1d0/2d0/(za(k4,e4))*za(p5,p2)*&
     &    za(k4,p5)*zb(p2,p1)*zb(p5,p2)*mt*ci - 1d0/2d0/(za(k4,e4))*za(&
     &    p5,k4)*za(k4,p1)*zb(p1,p2)*zb(k4,p1)*mt*ci - 1/(za(k4,e4))*&
     &    za(p5,k4)*za(k4,p5)*zb(p1,p2)*zb(k4,p5)*mt*ci + 1d0/2d0/(za(&
     &    k4,e4))*za(p5,k4)*za(k4,p5)*zb(p5,p2)*zb(k4,p1)*mt*ci + 1/(&
     &    za(k4,e4))*za(p5,k4)*zb(p1,p2)*mt**3*ci - 1d0/2d0/(za(k4,e4))&
     &    *za(p5,k4)*zb(p1,p2)*s12*mt*ci - 1d0/2d0/(za(k4,e4))*za(p5,k4&
     &    )*zb(p1,p2)*s14*mt*ci - 1d0/2d0/(za(k4,e4))*za(p5,k4)*zb(p1,&
     &    p2)*s45*mt*ci - 1d0/2d0/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*s25*&
     &    mt*ci - 1d0/2d0/(za(k4,e4))*za(p5,e4)*za(k4,p1)*zb(p1,p2)*zb(&
     &    e4,p1)*mt*ci - 1/(za(k4,e4))*za(p5,e4)*za(k4,p5)*zb(p1,p2)*&
     &    zb(e4,p5)*mt*ci + 1d0/2d0/(za(k4,e4))*za(p5,e4)*za(k4,p5)*zb(&
     &    p5,p2)*zb(e4,p1)*mt*ci )
      ampw(2) = ampw(2) + om2 * (  - 1d0/2d0*za(p5,e3)*zb(p2,e4)*zb(e3,&
     &    p1)*mw**(-2)*mt**3 + 1d0/2d0*za(p5,e3)*zb(p2,e4)*zb(e3,p1)*&
     &    s24*mw**(-2)*mt + 1d0/4d0*za(p5,e3)*zb(p2,e4)*zb(e3,p1)*s12*&
     &    mw**(-2)*mt + 1d0/4d0*za(p5,e3)*zb(p2,e4)*zb(e3,p1)*s14*&
     &    mw**(-2)*mt + 1d0/4d0*za(p5,e3)*zb(p2,e4)*zb(e3,p1)*s45*&
     &    mw**(-2)*mt + 1d0/4d0*za(p5,e3)*zb(p2,e4)*zb(e3,p1)*s25*&
     &    mw**(-2)*mt - 1d0/2d0*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*mw**(-2)*&
     &    mt**3 + 1d0/2d0*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*s24*mw**(-2)*mt&
     &     + 1d0/4d0*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*s12*mw**(-2)*mt + 1d0&
      &   /4d0*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*s14*mw**(-2)*mt + 1d0/4d0*&
     &    za(p5,k3)*zb(p2,e4)*zb(k3,p1)*s45*mw**(-2)*mt + 1d0/4d0*za(p5&
     &    ,k3)*zb(p2,e4)*zb(k3,p1)*s25*mw**(-2)*mt - 1d0/2d0/(za(k4,e4)&
     &    )*za(p5,e3)*za(k4,e3)*zb(e3,p1)*zb(e3,p2)*mt - 1d0/2d0/(za(k4&
     &    ,e4))*za(p5,e3)*za(k4,k3)*zb(e3,p1)*zb(k3,p2)*mt - 1d0/2d0/(&
     &    za(k4,e4))*za(p5,k3)*za(k4,e3)*zb(e3,p2)*zb(k3,p1)*mt - 1d0/2d0&
      &    /(za(k4,e4))*za(p5,k3)*za(k4,k3)*zb(k3,p1)*zb(k3,p2)*mt )
      ampw(2) = ampw(2) + om1 * (  - 1d0/2d0*za(p5,p2)*zb(p2,p1)*zb(p2,&
     &    e4)*mt - 1d0/2d0*za(p5,k4)*zb(p2,e4)*zb(k4,p1)*mt - 1d0/2d0*&
     &    za(p5,e4)*zb(p2,e4)*zb(e4,p1)*mt + 1/(za(k4,e4))*za(p5,k4)*&
     &    zb(p1,p2)*mw**2*mt )

        ampt(1) =  - 1d0/2d0*za(p5,e4)*zb(p1,p2)*mt*vev*KR - 1d0/2d0/(&
     &    zb(k4,e4))*za(p5,p1)*zb(p1,p2)*zb(p1,k4)*mt*vev*KL + 1d0/2d0&
     &    /(zb(k4,e4))*za(p5,p2)*zb(p2,p1)*zb(p2,k4)*mt*vev*KL

        ampt(2) =  - 1d0/2d0*za(p5,p1)*zb(p1,p2)*zb(p1,e4)*vev*KL + 1d0/&
     &    2d0*za(p5,p2)*zb(p2,p1)*zb(p2,e4)*vev*KL - 1d0/2d0/(za(k4,e4)&
     &    )*za(p5,k4)*zb(p1,p2)*mt**2*vev*KR

! -- include propagators
        ampt = ampt / (s125-mt**2+cI*M_Top*Ga_Top)/(s15-mw**2+cI*M_W*Ga_W)
        ampw = ampw / (s24-mw**2+cI*M_W*Ga_W)/(s15-mw**2+cI*M_W*Ga_W)

        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
      end subroutine ubhtdamp



      SUBROUTINE dbbarhtbaruamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production d(p1)+bbar(p2)->H(p3)+tbar(p4)+u(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),KL,KR
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

        amp(1)=(ampw(1)*(ghz1/2d0)  +ampt(1))*mdecay(1)
        amp(2)=(ampw(2)*(ghz1/2d0)  +ampt(2))*mdecay(2)

      END SUBROUTINE



      subroutine udbar_htbbaramp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production u(p1)+dbar(p2)->H(p3)+t(p4)+bbar(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),KL,KR
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

        amp(1)=(ampw(1)*(ghz1/2d0)  +ampt(1))*mdecay(1)
        amp(2)=(ampw(2)*(ghz1/2d0)  +ampt(2))*mdecay(2)

      end subroutine udbar_htbbaramp


      subroutine ubard_Htbarbamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,amp)
! amplitude for production ubar(p1)+d(p2)->H(p3)+tbar(p4)+b(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),KL,KR
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

        amp(1)=(ampw(1)*(ghz1/2d0)  +ampt(1))*mdecay(1)
        amp(2)=(ampw(2)*(ghz1/2d0)  +ampt(2))*mdecay(2)

      end subroutine ubard_Htbarbamp
      
      subroutine gb_twmHamp(Pg,Eg,Pb,UBb,Pt,UBt,Pw,ECw,Ph,amp)
! amplitude for production g(Pg)+b(Pb)->t(Pt)+wm(Pw)+h(Ph)
! allowing for scalar & pseudoscalar couplings of Higgs to top (Spinors and Pol-Vectors in Dirac representation!)
        implicit none
        real(8)    :: Pg(4),Pb(4),Pt(4),Pw(4),Ph(4)
        complex(8) :: UBb(4),Eg(4),UBt(4),ECw(4)
        complex(8) :: amp,ampts,ampws,amptt1,amptt2,ampwt
        real(8)    :: mt,mw,v
        complex(8)    :: kap,kapt,a1WW,a2WW,a4WW   


        mw=M_W
        mt=m_Top
        v=vev
        
        kap=kappa
        kapt=kappa_tilde
        a1WW=ghz1
        a2WW=ghz2
        a4WW=ghz4
        

        ampts = -((gs*mw*(-(DCONJG(UBb(2))*(Eg(1)*((Pb(1) + Pg(1))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   Eg(4)*(-((Pb(1) + Pg(1))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Pb(1) + Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))) -&
     &   DCONJG(UBb(4))*(Eg(1)*((Pb(1) + Pg(1))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   Eg(4)*(-((Pb(1) + Pg(1))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Pb(1) + Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))))) -&
     &   DCONJG(UBb(1))*(-((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(-((Pb(1) + Pg(1))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   + Eg(1)*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + Eg(4)*(-((Pb(2) + Pg(2)&
     &   + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))))) -&
     &   DCONJG(UBb(3))*(-((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(-((Pb(1) + Pg(1))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   + Eg(1)*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + Eg(4)*(-((Pb(2) + Pg(2)&
     &   + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))))))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) + Pg(2))**2 +&
     &   (Pb(3) + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mt**2 - (Ph(1) + Pt(1))**2 + (Ph(2) +&
     &   Pt(2))**2 + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)))

        ampws = (dcmplx(0d0,-1d0)*gs*mw*(DCONJG(UBb(1))*(Eg(1)*((Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) -&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)))) +&
     &   DCONJG(UBb(3))*(Eg(1)*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) +&
     &   Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1))&
     &   + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) -&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)))) +&
     &   DCONJG(UBb(2))*((Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1)&
     &   + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(1)*(-((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4))) + (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) + (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)))) +&
     &   DCONJG(UBb(4))*((Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1)&
     &   + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(1)*(-((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4))) + (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) + (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*UBt(4))))))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) + Pg(2))**2 +&
     &   (Pb(3) + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mw**2 - (Ph(1) + Pw(1))**2 + (Ph(2) +&
     &   Pw(2))**2 + (Ph(3) + Pw(3))**2 + (Ph(4) + Pw(4))**2))

        amptt1 = (dcmplx(0d0,-1d0)*gs*mt*mw*(DCONJG(UBb(1))*((ECw(1) + ECw(4))*((Pb(4) -&
     &   Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(1) + ECw(4))*((mt - Pb(1) +&
     &   Pw(1))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) + Pw(1))*(-(kapt*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kapt*(Pg(4) - Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) -&
     &   Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3)&
     &   - Eg(4)*UBt(4)) - dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4)))))) + DCONJG(UBb(3))*((ECw(1) +&
     &   ECw(4))*((Pb(4) - Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(1) + ECw(4))*((mt - Pb(1) +&
     &   Pw(1))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) + Pw(1))*(-(kapt*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kapt*(Pg(4) - Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) -&
     &   Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3)&
     &   - Eg(4)*UBt(4)) - dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4)))))) + DCONJG(UBb(2))*((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(4) - Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt +&
     &   Pb(1) - Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(1) - ECw(4))*((Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(1) -&
     &   ECw(4))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt +&
     &   Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) +&
     &   Pw(1))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4)))))) + DCONJG(UBb(4))*((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(4) - Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt +&
     &   Pb(1) - Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(1) - ECw(4))*((Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(1) -&
     &   ECw(4))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt +&
     &   Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) +&
     &   Pw(1))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))))))/(sqrt2*v**2*(mt**2 - (Pg(1) -&
     &   Pt(1))**2 + (Pg(2) - Pt(2))**2 + (Pg(3) - Pt(3))**2 + (Pg(4) - Pt(4))**2)*(mt**2 -&
     &   (Pb(1) - Pw(1))**2 + (Pb(2) - Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))
     
     amptt2 = (gs*mw*(-(DCONJG(UBb(2))*((ECw(1) - ECw(4))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2)&
     &   - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) - ECw(4))*((mt -&
     &   Pb(1) + Pw(1))*(mt*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) +&
     &   mt*(Ph(4) + Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))))))&
     &   - DCONJG(UBb(4))*((ECw(1) - ECw(4))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2)&
     &   - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) - ECw(4))*((mt -&
     &   Pb(1) + Pw(1))*(mt*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) +&
     &   mt*(Ph(4) + Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   - DCONJG(UBb(1))*((ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) +&
     &   ECw(4))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) - (Pb(4) - Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3)&
     &   + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(1) + ECw(4))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2)&
     &   + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) +&
     &   (mt - Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   - DCONJG(UBb(3))*((ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) +&
     &   ECw(4))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) - (Pb(4) - Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3)&
     &   + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(1) + ECw(4))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2)&
     &   + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) +&
     &   (mt - Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))))))/(sqrt2*v**2*(mt**2 - (Ph(1) + Pt(1))**2 + (Ph(2) + Pt(2))**2&
     &   + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)*(mt**2 - (Pb(1) - Pw(1))**2 + (Pb(2) -&
     &   Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))
     
        ampwt = (gs*mw*(dcmplx(0d0,-1d0)*DCONJG(UBb(1))*((a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) +&
     &   (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3))&
     &   - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) -&
     &   ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*(Eg(2)*Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) +&
     &   Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) - Eg(2)*Pt(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) -&
     &   Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) - Eg(4)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) - (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) +&
     &   Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1))&
     &   + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) +&
     &   Pw(4)))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) -&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) +&
     &   Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4)))) - dcmplx(0d0,1d0)*DCONJG(UBb(3))*((a1WW*mw**2*(ECw(2)&
     &   + dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) +&
     &   (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3))&
     &   - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) -&
     &   ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*(Eg(2)*Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) +&
     &   Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) - Eg(2)*Pt(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) -&
     &   Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) - Eg(4)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) - (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) +&
     &   Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1))&
     &   + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) +&
     &   Pw(4)))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) -&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) +&
     &   Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4)))) - DCONJG(UBb(4))*(dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(1)&
     &   - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) +&
     &   (a1WW*(mw**2*(dcmplx(0d0,1d0)*ECw(2) + ECw(3)) + (dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) + Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4))) +&
     &   dcmplx(0d0,2d0)*(a4WW*(ECw(3)*Ph(4)*Pw(1) - dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) -&
     &   ECw(1)*Ph(4)*Pw(3) + ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) +&
     &   Ph(1)*(dcmplx(0d0,1d0)*Pw(2) + Pw(3))) - ECw(3)*Ph(1)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4))))))*(Eg(2)*Pg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) + Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) -&
     &   Eg(2)*Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) - Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) -&
     &   Eg(4)*Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) - dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pw(2)&
     &   + Pw(3)))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) + ECw(2)*(Ph(2) +&
     &   Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) -&
     &   dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4)))) + DCONJG(UBb(2))*(dcmplx(0d0,-1d0)*(a1WW*mw**2*(ECw(1)&
     &   - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) -&
     &   dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pw(2) + Pw(3)))*(ECw(3)*Ph(3) + ECw(4)*Ph(4)&
     &   - ECw(1)*(Ph(1) + Pw(1)) + ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) +&
     &   2*a4WW*(ECw(3)*Ph(4)*Pw(1) - dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(Eg(2)*Pg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) + Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) -&
     &   Eg(2)*Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) - Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) -&
     &   Eg(4)*Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) + (a1WW*(mw**2*(dcmplx(0d0,1d0)*ECw(2) + ECw(3)) +&
     &   (dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pw(2) + Pw(3))*(ECw(3)*Ph(3) +&
     &   ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) + ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) +&
     &   ECw(4)*Pw(4))) + dcmplx(0d0,2d0)*(a4WW*(ECw(3)*Ph(4)*Pw(1) - dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2)&
     &   - ECw(1)*Ph(4)*Pw(3) + ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) +&
     &   Ph(1)*(dcmplx(0d0,1d0)*Pw(2) + Pw(3))) - ECw(3)*Ph(1)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4))))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) +&
     &   dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))))))/(sqrt2*v**2*(mt**2 - Pg(1)**2 + Pg(2)**2 +&
     &   Pg(3)**2 + Pg(4)**2 + 2*Pg(1)*Pt(1) - Pt(1)**2 - 2*Pg(2)*Pt(2) + Pt(2)**2 -&
     &   2*Pg(3)*Pt(3) + Pt(3)**2 - 2*Pg(4)*Pt(4) + Pt(4)**2)*(mw**2 - Ph(1)**2 + Ph(2)**2 +&
     &   Ph(3)**2 + Ph(4)**2 - 2*Ph(1)*Pw(1) - Pw(1)**2 + 2*Ph(2)*Pw(2) + Pw(2)**2 +&
     &   2*Ph(3)*Pw(3) + Pw(3)**2 + 2*Ph(4)*Pw(4) + Pw(4)**2))

         amp = ampts + ampws + amptt1 + amptt2 + ampwt

      end subroutine gb_twmHamp      


      subroutine gbb_tbwpHamp(Pg,Eg,Pb,Vb,Pt,Vt,Pw,ECw,Ph,amp)
! amplitude for production g(Pg)+b~(Pb)->t~(Pt)+wp(Pw)+h(Ph)
! allowing for scalar & pseudoscalar couplings of Higgs to top (Spinors and Pol-Vectors in Dirac representation!)
        implicit none
        real(8)    :: Pg(4),Pb(4),Pt(4),Pw(4),Ph(4)
        complex(8) :: Vb(4),Eg(4),Vt(4),ECw(4)
        complex(8) :: amp,ampts,ampws,amptt1,amptt2,ampwt
        real(8)    :: mt,mw,v
        complex(8)    :: kap,kapt,a1WW,a2WW,a4WW   


        mw=M_W
        mt=m_Top
        v=vev
                
        kap=kappa
        kapt=kappa_tilde
        a1WW=ghz1
        a2WW=ghz2
        a4WW=ghz4
        

        ampts = -((gs*mt*mw*(dcmplx(0d0,1d0)*kapt*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(1) - kap*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(1) - kap*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) - ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(2) + dcmplx(0d0,1d0)*kapt*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) + ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(2) + kap*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(3) - dcmplx(0d0,1d0)*kapt*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(3) - dcmplx(0d0,1d0)*kapt*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) - ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(4) + kap*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) + ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(4)))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) +&
     &   Pg(2))**2 + (Pb(3) + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mt**2 - (Ph(1) + Pt(1))**2 +&
     &   (Ph(2) + Pt(2))**2 + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)))

        ampws = -((gs*mw*(-((dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   ((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) - (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) + ECw(4))&
     &   - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1)&
     &   + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))*Vt(1)) -&
     &   (((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) +&
     &   Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) - (-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - (DCONJG(Vb(4))*Eg(1) -&
     &   DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))*Vt(2) +&
     &   (dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) - (-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - (DCONJG(Vb(4))*Eg(1) -&
     &   DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   ((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) - (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) + ECw(4))&
     &   - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1)&
     &   + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))*Vt(3) +&
     &   (((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) +&
     &   Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) - (-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - (DCONJG(Vb(4))*Eg(1) -&
     &   DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4))))))*Vt(4)))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) + Pg(2))**2 + (Pb(3)&
     &   + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mw**2 - (Ph(1) + Pw(1))**2 + (Ph(2) + Pw(2))**2 +&
     &   (Ph(3) + Pw(3))**2 + (Ph(4) + Pw(4))**2)))

        amptt1 = (gs*mt*mw*(-((dcmplx(0d0,1d0)*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((mt - Pg(1)&
     &   + Pt(1))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kapt*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))))) - Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   dcmplx(0d0,1d0)*(mt - Pg(1) + Pt(1))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + dcmplx(0d0,1d0)*Eg(1)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Pg(4) - Pt(4))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(1)) + (Eg(4)*(dcmplx(0d0,1d0)*(mt - Pg(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) -&
     &   Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   dcmplx(0d0,1d0)*(mt - Pg(1) + Pt(1))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + dcmplx(0d0,1d0)*Eg(1)*((Pg(4) - Pt(4))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) - kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))))&
     &   + dcmplx(0d0,1d0)*(mt + Pg(1) - Pt(1))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))))&
     &   + (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(2) - (Eg(1)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   dcmplx(0d0,1d0)*(mt - Pg(1) + Pt(1))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + dcmplx(0d0,1d0)*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(2) - dcmplx(0d0,1d0)*Pg(3) -&
     &   Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   - dcmplx(0d0,1d0)*Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Pg(4) - Pt(4))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(3) + dcmplx(0d0,1d0)*(Eg(1)*((mt - Pg(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kapt*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))))) + Eg(4)*((Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(2) - dcmplx(0d0,1d0)*Pg(3) -&
     &   Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Pg(4) - Pt(4))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(4)))/(sqrt2*v**2*(mt**2 - (Pg(1) - Pt(1))**2 + (Pg(2) - Pt(2))**2 +&
     &   (Pg(3) - Pt(3))**2 + (Pg(4) - Pt(4))**2)*(mt**2 - (Pb(1) - Pw(1))**2 + (Pb(2) -&
     &   Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))

        amptt2 = (gs*mt*mw*(kapt*((Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   Eg(1)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(1) + dcmplx(0d0,1d0)*kap*((mt - Ph(1) - Pt(1))*((Eg(2)&
     &   + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) + (Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))*Vt(1) - kapt*((Ph(4) +&
     &   Pt(4))*(Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) + (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3)))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*(Eg(4)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4)))))*Vt(2) - dcmplx(0d0,1d0)*kap*((mt - Ph(1) - Pt(1))*(Eg(1)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(4) + Pt(4))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3)))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(2) - dcmplx(0d0,1d0)*kap*((Ph(2) +&
     &   dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   Eg(1)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(3) - kapt*((mt - Ph(1) - Pt(1))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) + (Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))*Vt(3) +&
     &   dcmplx(0d0,1d0)*kap*((Ph(4) + Pt(4))*(Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3)))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*(Eg(4)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4)))))*Vt(4) + kapt*((mt - Ph(1) - Pt(1))*(Eg(1)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(4) + Pt(4))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3)))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(4)))/(sqrt2*v**2*(mt**2 - (Ph(1) + Pt(1))**2 +&
     &   (Ph(2) + Pt(2))**2 + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)*(mt**2 - (Pb(1) -&
     &   Pw(1))**2 + (Pb(2) - Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))
     
        ampwt = (dcmplx(0d0,-1d0)*gs*mw*((-((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))))) +&
     &   Eg(1)*(-((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) + (mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) -&
     &   Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))))*Vt(1) -&
     &   (Eg(1)*(-((mt + Pg(1) - Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) - Eg(4)*((mt -&
     &   Pg(1) + Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))))*Vt(2) +&
     &   ((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(-((mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) +&
     &   Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) +&
     &   Eg(1)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))))*Vt(3) -&
     &   (Eg(4)*(-((mt + Pg(1) - Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) - Eg(1)*((mt -&
     &   Pg(1) + Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(-((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) + (mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) +&
     &   Pw(4))))))))*Vt(4)))/(sqrt2*v**2*(mt**2 - (Pg(1) - Pt(1))**2 + (Pg(2) - Pt(2))**2&
     &   + (Pg(3) - Pt(3))**2 + (Pg(4) - Pt(4))**2)*(mw**2 - (Ph(1) + Pw(1))**2 + (Ph(2) +&
     &   Pw(2))**2 + (Ph(3) + Pw(3))**2 + (Ph(4) + Pw(4))**2))

         amp = ampts + ampws + amptt1 + amptt2 + ampwt

      end subroutine gbb_tbwpHamp
      



    SUBROUTINE TDECAY(k4,e4,b,ep,nu,za,zb,dkamp)
! top decay routine, taken from MCFM, see hep-ph:/1204.1513
       implicit none
       integer :: k4,e4,b,ep,nu
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp


   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = za(b,nu)*zb(ep,e4)
       dkamp(2) = m_top * za(b,nu)*zb(ep,k4)/za(e4,k4)

       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top)
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W

       dkamp = dkamp * WProp * NWAFactor_Top * gwsq


     END SUBROUTINE


    SUBROUTINE ATDECAY(k4,e4,bbar,em,nubar,za,zb,dkamp)
! anti-top decay routine, taken from MCFM, see hep-ph:/1204.1513
       implicit none
       integer :: k4,e4,bbar,em,nubar
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp

   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = -m_top * zb(bbar,nubar)*za(em,k4)/zb(e4,k4)
       dkamp(2) = -zb(bbar,nubar)*za(em,e4)

       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top)
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W

       dkamp = dkamp * WProp * NWAFactor_Top * gwsq


     END SUBROUTINE



     SUBROUTINE myAmp(p,VVHg1,VVHg2,VVHg3,TTHg1,TTHg2,res)! 1:t_out 2:b_out 3:q_out 4:qbar_out
     implicit none
     real(8) :: p(1:4,1:4),pt_light(1:4),res
     real(8) :: s12,s34,s234,IV1,IV2,sprod(1:4,1:4)
     complex(8) :: VVHg1,VVHg2,VVHg3,TTHg1,TTHg2
     complex(8) :: HelAmp(1:2),Props1,Props2
     complex(8) :: za(1:4,1:4),zb(1:4,1:4)
     complex(8),parameter :: cI=(0d0,1d0)

        pt_light(1:4) = p(1:4,1) - m_top**2/(2d0*(p(1,3)*p(1,1)-p(2,3)*p(2,1)-p(3,3)*p(3,1)-p(4,3)*p(4,1)))*p(1:4,3)
        call my_spinoru(4,(/pt_light(1:4),p(1:4,2),p(1:4,3),p(1:4,4)/),za,zb,sprod)
        sprod(1,2) = 2d0*(p(1,2)*p(1,1)-p(2,2)*p(2,1)-p(3,2)*p(3,1)-p(4,2)*p(4,1))!  overwriting sprod(1,2) because 1=not top momentum above

        s12 = m_Top**2 + sprod(1,2)
        s34 = sprod(3,4)
        Props1 =  cI/(s12-M_W**2+cI*M_W*Ga_W) * cI/(s34-M_W**2+cI*M_W*Ga_W)           ! all propagator denominators for the HVV diagram
        s234 = sprod(2,3)+sprod(2,4)+sprod(3,4)
        Props2 =  cI/(s234-M_Top**2+cI*M_Top*Ga_Top) * cI/(s34-M_W**2+cI*M_W*Ga_W)    ! all propagator denominators for the Yukawa diagram

        IV1 = dsqrt(gwsq)/dSqrt(2d0)
        IV2 = dsqrt(gwsq)/dSqrt(2d0)

        ! assuming tbar*( TTHg1 + i*TTHg2*gamma^5 )*t
        HelAmp(1) =  (IV1*IV2*(-(m_top**3*Props1*VVHg2*zb(2,1)*(2*za(2,3)*zb(2,3) + 2*za(3,1)*zb(3,1) + za(3,4)*zb(3,4))*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1))) +   &
           4*M_W**2*Props2*(TTHg1 - cI*TTHg2)*za(3,1)*zb(2,4)*zb(3,1)*(za(2,3)*zb(2,1) - za(3,4)*zb(4,1)) +   &
           m_top*Props1*zb(3,1)*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1))*   &
            (-2*VVHg1*za(3,1)*zb(2,1) + 2*M_W**2*(VVHg2*za(3,1)*zb(2,1) + (VVHg2 - cI*VVHg3)*za(3,4)*zb(2,4)) +    &
              VVHg2*za(3,1)*zb(2,1)*(2*za(2,1)*zb(2,1) + za(2,3)*zb(2,3) + za(2,4)*zb(2,4) + za(3,1)*zb(3,1) + za(4,1)*zb(4,1)))))/(2*M_W**2*za(3,1)*zb(3,1))

         HelAmp(2) =  (IV1*IV2*(-4*m_top*M_W**2*Props2*za(3,1)*zb(2,4)*zb(3,1)*   &
            (-(TTHg1*za(2,3)*zb(2,3)) + cI*TTHg2*za(2,3)*zb(2,3) + TTHg1*za(3,1)*zb(3,1) + cI*TTHg2*za(3,1)*zb(3,1) + (-TTHg1 + cI*TTHg2)*za(3,4)*zb(3,4)) -    &
           m_top**4*Props1*VVHg2*zb(2,3)*(2*za(2,3)*zb(2,3) + 2*za(3,1)*zb(3,1) + za(3,4)*zb(3,4))*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1)) +     &
           m_top**2*Props1*za(3,1)*zb(3,1)*(-2*VVHg1*zb(2,3)*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1)) +     &
              VVHg2*zb(2,3)*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1))*(2*za(2,1)*zb(2,1) + za(2,3)*zb(2,3) + za(2,4)*zb(2,4) + za(3,1)*zb(3,1) + za(4,1)*zb(4,1)) +   &
              2*M_W**2*(cI*VVHg3*za(3,4)*zb(2,4)*zb(3,4) + VVHg2*zb(2,3)*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1)))) -   &
           2*M_W**2*Props1*za(3,1)*zb(3,1)**2*(-2*VVHg1*za(3,1)*zb(2,4) + VVHg2*(za(3,1)*zb(2,3) + za(4,1)*zb(2,4))*(za(2,3)*zb(2,4) + za(3,1)*zb(4,1)) +   &
              cI*VVHg3*(za(2,1)*za(3,4)*zb(2,4)**2 + za(3,1)**2*(-(zb(2,4)*zb(3,1)) + zb(2,3)*zb(4,1))))))/(2*M_W**2*za(3,1)*zb(3,1)**2)

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
