!--YaofuZhou-----------------------------------------
! finite factor multiplied by LO cross section for
! interference of virtual and LO, plus "I bit" "K bit" and "P bit" from arXiv:0709.2881
module ModVHvirtual
  use ModParameters
  use ModMisc
  use ModVHaux
  implicit none
  
public :: fac_VH_virtual

contains

real(8) function fac_VH_virtual(shat)
!  use Collier
  implicit none
  real(8), intent(in) :: shat
  complex(8) :: B0_0, B0_shat, C0
!  real(8) :: ReB0_0, ReB0_shat, ReC0

!  call SetMuUV2_cll(Mu_Ren**2)
!  call SetMuIR2_cll(Mu_Ren**2)

!  call SetMuUV2_cll(9d0)
!  call SetMuIR2_cll(9d0)
!  shat=16d0

!  call B0_cll(B0_0,dcmplx(0d0),dcmplx(0d0),dcmplx(0d0))
!  call B0_cll(B0_shat,dcmplx(shat),dcmplx(0d0),dcmplx(0d0))
!  call C0_cll(C0,dcmplx(0d0),dcmplx(shat),dcmplx(0d0),dcmplx(0d0),dcmplx(0d0),dcmplx(0d0))

!print*,B0_shat !COLLIER
!print*,2d0-cdlog( -dcmplx(shat)/dcmplx(Mu_Ren**2) ) !Denner + COLLIER eps-normalization
!print*,2d0-dlog( shat/Mu_Ren**2 ) ! real part of Denner + COLLIER eps-normalization
!print*,B0_0
!print*,B0_shat
!print*,C0
!print*,(-pi**2+3d0*(cdlog(-dcmplx(Mu_Ren**2)/dcmplx(shat)))**2)/6/shat !imaginary part has a negative sign
!print*,real((cdlog(-dcmplx(Mu_Ren**2)/dcmplx(shat)))**2,8)
!print*,(dlog(Mu_Ren**2/shat))**2-pisq
!print*,"----------"
!pause

!  B0_0 = B0_0 * ci * pisq / (2d0*pi)**4
!  B0_shat = B0_shat * ci * pisq / (2d0*pi)**4
!  C0 = C0 * ci * pisq / (2d0*pi)**4

!  ReB0_0 = real(B0_0,8)
!  ReB0_shat = real(B0_shat,8)
!  ReC0 = real(C0,8)

!  fac_VH_virtual = !- real(B0_shat,8)*3d0/4d0 &
                    !- real( 2d0-cdlog( -dcmplx(shat)/dcmplx(Mu_Ren**2) ) ,8 )*3d0/4d0 & !- real(B0_shat,8)*3d0/4d0
  fac_VH_virtual = - ( 2d0-dlog( shat/Mu_Ren**2 ) )*3d0/4d0 & !- real(B0_shat,8)*3d0/4d0
                   - 0.5d0 &
                   - ((dlog(Mu_Ren**2/shat))**2-pisq)/4d0
                    !- real((cdlog(-dcmplx(Mu_Ren**2)/dcmplx(shat)))**2,8) /4d0
      
  fac_VH_virtual = 2d0 * fac_VH_virtual * Cf * alphas / pi !for 2Re{(iM_virt + iM_c)^* iM_LO}

!if(isNan(fac_VH_virtual))then
!print*,B0_0, B0_shat, C0, shat
!pause
!endif

  return
END function fac_VH_virtual

end module ModVHvirtual
!!--YaofuZhou-----------------------------------------
