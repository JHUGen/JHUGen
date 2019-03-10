!--YaofuZhou-----------------------------------------
module ModffbHBox11_12mpm
  use ModParameters
  use ModMisc
  use COLLIER
  implicit none
  public :: ffbHBox11_12mpm

contains

subroutine ffbHBox11_12mpm(MomExt,Spaa,Spbb,sprod,Box11_12)
  implicit none
  real(8), intent(in) :: MomExt(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  complex(8), intent(out) :: Box11_12
  complex(8) :: MomInv(1:6),masses2(0:3)
  complex(8) :: Dcoeff(0:1,0:3,0:3,0:3),Dcoeffuv(0:1,0:3,0:3,0:3)
  real(8) :: Derr(0:3)
  complex(8) :: D0,D1,D2,D3
  complex(8) :: D00,D11,D22,D33,D12,D21,D13,D31,D23,D32
  complex(8) :: D001,D002,D003,D111,D222,D333,D112,D113,D221,D223,D331,D332,D123
  complex(8) :: cA,cV
  real(8) :: mloop
  integer rank

  rank = 3

  
  mloop = m_top

  MomInv(1) = sprod(3,4)
  MomInv(2) = 0d0
  MomInv(3) = MomExt(:,5).dot.MomExt(:,5)
  MomInv(4) = 0d0
  MomInv(5) = (MomExt(:,1)-MomExt(:,6)-MomExt(:,7)).dot.(MomExt(:,1)-MomExt(:,6)-MomExt(:,7))
  MomInv(6) = (MomExt(:,2)-MomExt(:,6)-MomExt(:,7)).dot.(MomExt(:,2)-MomExt(:,6)-MomExt(:,7))

  masses2(0:3) = (/m_top**2,m_top**2,m_top**2,m_top**2/)
	
  call SetMuUV2_cll(Mu_Ren**2)
  call SetMuIR2_cll(Mu_Ren**2)

  call D_cll(Dcoeff,Dcoeffuv,MomInv(1:6),masses2(0:3),rank,Derr(0:3))

  D0 = Dcoeff(0,0,0,0)
  D1 = Dcoeff(0,1,0,0)
  D2 = Dcoeff(0,0,1,0)
  D3 = Dcoeff(0,0,0,1)
  D00 = Dcoeff(1,0,0,0)
  D11 = Dcoeff(0,2,0,0)
  D22 = Dcoeff(0,0,2,0)
  D33 = Dcoeff(0,0,0,2)
  D12 = Dcoeff(0,1,1,0)
  D21 = Dcoeff(0,1,1,0)
  D13 = Dcoeff(0,1,0,1)
  D31 = Dcoeff(0,1,0,1)
  D23 = Dcoeff(0,0,1,1)
  D32 = Dcoeff(0,0,1,1)
  D001 = Dcoeff(1,1,0,0)
  D002 = Dcoeff(1,0,1,0)
  D003 = Dcoeff(1,0,0,1)
  D111 = Dcoeff(0,3,0,0)
  D222 = Dcoeff(0,0,3,0)
  D333 = Dcoeff(0,0,0,3)
  D112 = Dcoeff(0,2,1,0)
  D113 = Dcoeff(0,2,0,1)
  D221 = Dcoeff(0,1,2,0)
  D223 = Dcoeff(0,0,2,1)
  D331 = Dcoeff(0,1,0,2)
  D332 = Dcoeff(0,0,1,2)
  D123 = Dcoeff(0,1,1,1)

      Box11_12 = &
     &  - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2&
     & )*Spaa(2,3)**2*Spbb(1,2)*Spbb(1,3)*Spbb(1,4)*D2*mloop - 32d0/(&
     & dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,&
     & 3)**2*Spbb(1,2)*Spbb(1,3)*Spbb(1,4)*D22*mloop - 32d0/(dsqrt(2d0))&
     & /(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,3)**2*&
     & Spbb(1,2)*Spbb(1,3)*Spbb(1,4)*D12*mloop - 32d0/(dsqrt(2d0))/(Spaa(&
     & 1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)*&
     & Spbb(1,3)*Spbb(1,4)*D13*mloop - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)*Spbb(1,3&
     & )*Spbb(1,4)*D23*mloop - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))&
     & /(Spbb(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)**&
     & 2*D2*mloop - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))&
     & *Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)**2*D22*mloop&
     &  - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)&
     & *Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)**2*D12*mloop - 32d0/(&
     & dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,&
     & 3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)**2*D13*mloop
      Box11_12 = Box11_12 - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(&
     & 1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)**2*D23*mloop - 16d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(&
     & 2,3)*Spaa(3,4)*Spbb(1,3)*Spbb(1,4)**2*D22*mloop - 16d0/(dsqrt(2d0)&
     & )/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(3&
     & ,4)*Spbb(1,3)*Spbb(1,4)**2*D12*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2&
     & ))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,2)*Spaa(2,4)*Spaa(3,4)*Spbb(1,&
     & 4)**3*D22*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(&
     & 1,2))*Spaa(1,2)*Spaa(2,4)*Spaa(3,4)*Spbb(1,4)**3*D12*mloop + 16d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,3)*Spaa(&
     & 2,3)*Spaa(2,4)*Spbb(1,3)*Spbb(1,4)**2*D1*mloop + 16d0/(dsqrt(2d0))&
     & /(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2,&
     & 4)*Spbb(1,3)*Spbb(1,4)**2*D2*mloop + 32d0/(dsqrt(2d0))/(Spaa(1,2))&
     & /(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2,4)*Spbb(1,3)&
     & *Spbb(1,4)**2*D11*mloop + 48d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0)&
     & )/(Spbb(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2,4)*Spbb(1,3)*Spbb(1,4)&
     & **2*D22*mloop
      Box11_12 = Box11_12 + 80d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(&
     & 1,3)*Spaa(2,3)*Spaa(2,4)*Spbb(1,3)*Spbb(1,4)**2*D12*mloop + 16d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,3)*Spaa(&
     & 2,4)**2*Spbb(1,4)**3*D1*mloop + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(1,3)*Spaa(2,4)**2*Spbb(1,4)**3*D2*&
     & mloop + 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*&
     & Spaa(1,3)*Spaa(2,4)**2*Spbb(1,4)**3*D11*mloop + 48d0/(dsqrt(2d0))&
     & /(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,3)*Spaa(2,4)**2*&
     & Spbb(1,4)**3*D22*mloop + 80d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))&
     & /(Spbb(1,2))*Spaa(1,3)*Spaa(2,4)**2*Spbb(1,4)**3*D12*mloop - 16d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,4)*Spaa(&
     & 2,3)**2*Spbb(1,3)*Spbb(1,4)**2*D1*mloop - 16d0/(dsqrt(2d0))/(Spaa(&
     & 1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,4)*Spaa(2,3)**2*Spbb(1,3)*&
     & Spbb(1,4)**2*D2*mloop - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))&
     & /(Spbb(1,2))*Spaa(1,4)*Spaa(2,3)**2*Spbb(1,3)*Spbb(1,4)**2*D11*&
     & mloop
      Box11_12 = Box11_12 - 48d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(&
     & 1,4)*Spaa(2,3)**2*Spbb(1,3)*Spbb(1,4)**2*D22*mloop - 80d0/(dsqrt(&
     & 2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,4)*Spaa(2,3)**2*&
     & Spbb(1,3)*Spbb(1,4)**2*D12*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4)**&
     & 3*D1*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))&
     & *Spaa(1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4)**3*D2*mloop - 32d0/(&
     & dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,4)*Spaa(2,&
     & 3)*Spaa(2,4)*Spbb(1,4)**3*D11*mloop - 48d0/(dsqrt(2d0))/(Spaa(1,2)&
     & )/(dsqrt(2d0))/(Spbb(1,2))*Spaa(1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4&
     & )**3*D22*mloop - 80d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1&
     & ,2))*Spaa(1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4)**3*D12*mloop - 16d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*&
     & Spbb(1,2)*Spbb(1,3)*Spbb(3,4)*D0*mloop - 32d0/(dsqrt(2d0))/(Spaa(1&
     & ,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,2)*Spbb(1,3)*&
     & Spbb(3,4)*D1*mloop
      Box11_12 = Box11_12 - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(&
     & 2,3)**3*Spbb(1,2)*Spbb(1,3)*Spbb(3,4)*D2*mloop - 32d0/(dsqrt(2d0))&
     & /(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,2)*&
     & Spbb(1,3)*Spbb(3,4)*D11*mloop - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,2)*Spbb(1,3)*Spbb(3,4&
     & )*D22*mloop - 64d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2)&
     & )*Spaa(2,3)**3*Spbb(1,2)*Spbb(1,3)*Spbb(3,4)*D12*mloop - 16d0/(&
     & dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*&
     & Spbb(1,2)*Spbb(1,3)*Spbb(3,4)*D13*mloop - 16d0/(dsqrt(2d0))/(Spaa(&
     & 1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,2)*Spbb(1,3)*&
     & Spbb(3,4)*D23*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(&
     & Spbb(1,2))*Spaa(2,3)**3*Spbb(1,3)**2*Spbb(2,4)*D1*mloop - 16d0/(&
     & dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*&
     & Spbb(1,3)**2*Spbb(2,4)*D2*mloop + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,3)**2*Spbb(2,4)*D13*&
     & mloop
      Box11_12 = Box11_12 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(&
     & 2,3)**3*Spbb(1,3)**2*Spbb(2,4)*D23*mloop + 16d0/(dsqrt(2d0))/(&
     & Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,3)*Spbb(1&
     & ,4)*Spbb(2,3)*D1*mloop + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))&
     & /(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,3)*Spbb(1,4)*Spbb(2,3)*D2*mloop&
     &  - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)&
     & **3*Spbb(1,3)*Spbb(1,4)*Spbb(2,3)*D13*mloop - 16d0/(dsqrt(2d0))/(&
     & Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**3*Spbb(1,3)*Spbb(1&
     & ,4)*Spbb(2,3)*D23*mloop - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))&
     & /(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4&
     & )*D0*mloop - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))&
     & *Spaa(2,3)**2*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D1*mloop&
     &  - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)&
     & **2*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D2*mloop - 32d0/(&
     & dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*&
     & Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D11*mloop
      Box11_12 = Box11_12 - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(&
     & 2,3)**2*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D22*mloop - 64d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*&
     & Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D12*mloop - 16d0/(dsqrt(&
     & 2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*&
     & Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D13*mloop - 16d0/(dsqrt(2d0))/(Spaa(&
     & 1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*Spbb(1,2)*&
     & Spbb(1,4)*Spbb(3,4)*D23*mloop - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*Spbb(1,3)*Spbb(1,4&
     & )*Spbb(2,4)*D0*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(&
     & Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*Spbb(1,3)*Spbb(1,4)*Spbb(2,4)*&
     & D1*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*&
     & Spaa(2,3)**2*Spaa(2,4)*Spbb(1,3)*Spbb(1,4)*Spbb(2,4)*D2*mloop + &
     & 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2&
     & *Spaa(2,4)*Spbb(1,3)*Spbb(1,4)*Spbb(2,4)*D13*mloop + 16d0/(dsqrt(&
     & 2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*&
     & Spbb(1,3)*Spbb(1,4)*Spbb(2,4)*D23*mloop
      Box11_12 = Box11_12 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2&
     & ,3)**2*Spaa(2,4)*Spbb(1,4)**2*Spbb(2,3)*D0*mloop + 16d0/(dsqrt(2d0&
     & ))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*&
     & Spbb(1,4)**2*Spbb(2,3)*D1*mloop + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(&
     & dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*Spaa(2,4)*Spbb(1,4)**2*Spbb(&
     & 2,3)*D2*mloop - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,&
     & 2))*Spaa(2,3)**2*Spaa(2,4)*Spbb(1,4)**2*Spbb(2,3)*D13*mloop - 16d0&
     & /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spbb(1,2))*Spaa(2,3)**2*&
     & Spaa(2,4)*Spbb(1,4)**2*Spbb(2,3)*D23*mloop
  
  return
end subroutine ffbHBox11_12mpm

end module ModffbHBox11_12mpm
!!--YaofuZhou-----------------------------------------
