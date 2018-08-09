!--YaofuZhou-----------------------------------------
! finite factor multiplied by LO cross section for
! interference of virtual and LO, plus "I bit" "K bit" and "P bit" from arXiv:0709.2881
module ModVHvirtual
  use ModParameters
  use ModMisc
  use ModVHaux
  implicit none
  
public :: fac_ZH_virtual
public :: IPK_ZH
private :: theta_step

contains

real(8) function fac_ZH_virtual(s12)
  use Collier
  implicit none
  real(8), intent(in) :: s12
  complex(8) :: B0_0, B0_s12, C0

  call SetMuUV2_cll(Mu_Ren**2)
  call SetMuIR2_cll(Mu_Ren**2)

!  call SetMuUV2_cll(9d0)
!  call SetMuIR2_cll(9d0)
!  s12=16d0

  call B0_cll(B0_0,dcmplx(0d0),dcmplx(0d0),dcmplx(0d0))
  call B0_cll(B0_s12,dcmplx(s12),dcmplx(0d0),dcmplx(0d0))
  call C0_cll(C0,dcmplx(0d0),dcmplx(s12),dcmplx(0d0),dcmplx(0d0),dcmplx(0d0),dcmplx(0d0))
!  call C0_cll(C0,dcmplx(0d0),dcmplx(16d0),dcmplx(0d0),dcmplx(0d0),dcmplx(0d0),dcmplx(0d0))
!print*,"B0_0 = ",B0_0
!print*,"B0_s12 = ",B0_s12
!print*,"C0 = ",C0
  fac_ZH_virtual = dreal(B0_0)*2d0 &
                 - dreal(B0_s12)*3d0/2d0 &
                 - dreal(C0)*s12 &
                 - 1d0 &
                 + 1d0/2d0*(dlog(4d0*pi*Mu_Ren**2/s12))**2 &
                 - 1d0/2d0*(2d0*gamma_0-3d0)*dlog(4d0*pi*Mu_Ren**2/s12) &
                 + 1d0/12d0*pisq &
                 + 1d0/2d0*gamma_0**2 &
                 - 3d0/2d0*gamma_0 &
                 + (dlog(alpha_dip))**2
      
  fac_ZH_virtual = fac_ZH_virtual * Cf * alphas / pi

if(isNan(fac_ZH_virtual))then
print*,B0_0, B0_s12, C0, s12
pause
endif

  return
END function fac_ZH_virtual



real(8) function IPK_ZH(Mu_Fact_x,alphasx,xsection,xsectionx,shat,shatx,pdf,pdfx,PreFac,PreFacx,x)
  implicit none
  real(8), intent(in) :: Mu_Fact_x
  real(8), intent(in) :: alphasx
  real(8), intent(in) :: xsection,xsectionx
  real(8), intent(in) :: shat, shatx
  real(8), intent(in) :: pdf(1:2), pdfx(1:2)
  real(8), intent(in) :: PreFac,PreFacx
  real(8), intent(in) :: x

  IPK_ZH = 0d0

  IPK_ZH = IPK_ZH + xsectionx * pdfx(1) * pdfx(2) * PreFacx * &
                   ( (1d0-x) &
                   - (1d0+x)*dlog((1d0-x)**2/x/alpha_dip) &
                   + 2d0/(1d0-x)*dlog(2d0-x) * (1d0 - theta_step(x,1d0-alpha_dip)) &
                   + (1d0+x**2)/(1d0-x)*dlog(min(1d0,alpha_dip/(1d0-x))) &
                   - 2d0/(1d0-x)*dlog((2d0-x)/(1d0-x)) * theta_step(1d0-alpha_dip,x) )

  IPK_ZH = IPK_ZH + (xsectionx * pdfx(1) * pdfx(2) * PreFacx &
                   - xsection  * pdf(1)  * pdf(2)  * PreFac) &
                  * 2d0/(1d0-x) * ( dlog((1d0-x)/x) + dlog(1d0-x)*theta_step(x,1d0-alpha_dip) )

  IPK_ZH = IPK_ZH - (xsectionx * pdfx(1) * pdfx(2) * PreFacx * dlog(4d0*pi*Mu_Fact_x**2/shatx) &
                   - xsection  * pdf(1)  * pdf(2)  * PreFac  * dlog(4d0*pi*Mu_Fact_x**2/shat) ) &
                  * (1d0+x**2)/(1d0-x)

  IPK_ZH = IPK_ZH * alphasx*Cf/2d0/pi

!if(IPK_ZH.ge.10000000d0)  then
!  print*,"x",xsectionx, pdfx(1), pdfx(2), x 
!  print*,"1",xsection ,pdf(1)  , pdf(2) 
!  print*,IPK_ZH,dsqrt(shatx)-M_Z-M_reso,dsqrt(shat),dsqrt(shatx)
!endif
  return
end function IPK_ZH


real(8) function theta_step(a,b)
  implicit none
  real(8), intent(in) :: a,b

  if(a.gt.b)then
    theta_step = 1d0
  else
    theta_step = 0d0
  endif

  return
end function theta_step

end module ModVHvirtual
!!--YaofuZhou-----------------------------------------
