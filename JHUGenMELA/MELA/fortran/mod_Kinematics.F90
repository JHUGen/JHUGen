MODULE ModKinematics
implicit none
save

public :: GetBWPropagator,ReweightBWPropagator,SetRunningScales,setPDFs,EvalAlphaS

CONTAINS


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
        DKFlavor = ZLepBranching_flat( DKRnd )!= ElM or MuM
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
        DKFlavor = ZNuBranching_flat( DKRnd )!= NuE,NuM,NuT
        MY_IDUP(6) =-DKFlavor
        MY_IDUP(7) =+DKFlavor
   elseif( DecayMode1.eq.4 ) then! W1(+)->lnu
        call random_number(DKRnd)
        MY_IDUP(4) = Wp_
        DKFlavor = WLepBranching_flat( DKRnd )!= ElM or MuM
        MY_IDUP(6) = +abs(DKFlavor)     ! lepton(+)
        MY_IDUP(7) = +abs(DKFlavor)+7   ! neutrino
   elseif( DecayMode1.eq.5 ) then! W1(+)->2q
        call random_number(DKRnd)
        MY_IDUP(4) = Wp_
        DKFlavor = WQuaUpBranching_flat( DKRnd )!= Up,Chm
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
        MY_IDUP(6) = Not_a_particle_
        MY_IDUP(7) = Not_a_particle_
   elseif( DecayMode1.eq.8 ) then! Z1->2l+2tau
        call random_number(DKRnd)
        MY_IDUP(4) = Z0_
        DKFlavor = ZLepPlusTauBranching_flat( DKRnd )!= ElM or MuM or TaM
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
        DKFlavor = WLepPlusTauBranching_flat( DKRnd )!= ElM or MuM or TaM
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
        DKFlavor = ZLepBranching_flat( DKRnd )!= ElM or MuM
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
        DKFlavor = ZNuBranching_flat( DKRnd )!= NuE,NuM,NuT
        MY_IDUP(8) =-DKFlavor
        MY_IDUP(9) =+DKFlavor
   elseif( DecayMode2.eq.4 ) then! W2(-)->lnu
        call random_number(DKRnd)
        MY_IDUP(5) = Wm_
        DKFlavor = WLepBranching_flat( DKRnd )!= ElM or MuM
        MY_IDUP(8) = -abs(DKFlavor)-7   ! anti-neutrino
        MY_IDUP(9) = -abs(DKFlavor)     ! lepton(-)
   elseif( DecayMode2.eq.5 ) then! W2(-)->2q (sample over u,d,s,c)
        call random_number(DKRnd)
        MY_IDUP(5) = Wm_
        DKFlavor = WQuaUpBranching_flat( DKRnd )!= Up,Chm
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
        MY_IDUP(8) = Not_a_particle_
        MY_IDUP(9) = Not_a_particle_
   elseif( DecayMode2.eq.8 ) then! Z2->2l+2tau
        call random_number(DKRnd)
        MY_IDUP(5) = Z0_
        DKFlavor = ZLepPlusTauBranching_flat( DKRnd )!= ElM or MuM or TaM
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
        DKFlavor = WLepPlusTauBranching_flat( DKRnd )!= ElM or MuM or TaM
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

      alphas=alphas_mz/(1.0_dp+alphas_mz*B0_PDF(NF5)*2.0_dp*dlog((Mu_Ren/zmass_pdf)))
#endif
      ! Calculate the derived couplings
      call ComputeQCDVariables()
   RETURN
end subroutine EvalAlphaS




END MODULE
