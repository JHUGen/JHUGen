!--YaofuZhou-----------------------------------------
module ModggboxHH3pp
  use ModParameters
  use ModMisc
  use COLLIER
  implicit none
  public :: ggboxHH3pp

contains

subroutine ggboxHH3pp(mom,Spaa,Spbb,sprod,ggboxHH3)
  implicit none
  real(8), intent(in) :: mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  complex(8), intent(out) :: ggboxHH3
  complex(8) :: MomInv(1:6),masses2(0:3)
  complex(8) :: Dcoeff(0:2,0:4,0:4,0:4),Dcoeffuv(0:2,0:4,0:4,0:4)
  real(8) :: Derr(0:4)
  complex(8) :: D0,D1,D2,D3
  complex(8) :: D00,D11,D22,D33,D12,D21,D13,D31,D23,D32
  complex(8) :: D001,D002,D003,D111,D222,D333,D112,D113,D221,D223,D331,D332,D123
  complex(8) :: D0000,D0011,D0012,D0013,D0022,D0023,D0033,D1111,D1112,D1113,D1122,D1123,D1133,D1222,D1223,D1233,D1333,D2222,D2223,D2233,D2333,D3333
  real(8) :: DimST,mloop
  integer rank

  rank = 4

  DimST = 4d0
  mloop = m_top

  MomInv(1) = 0d0
  MomInv(2) = mom(:,5).dot.mom(:,5)
  MomInv(3) = sprod(3,4)
  MomInv(4) = 0d0
  MomInv(5) = (mom(:,2)-mom(:,6)-mom(:,7)).dot.(mom(:,2)-mom(:,6)-mom(:,7))
  MomInv(6) = sprod(1,2)

  masses2(0:3) = (/m_top**2,m_top**2,m_top**2,m_top**2/)
     
  call SetMuUV2_cll(Mu_Ren**2)
  call SetMuIR2_cll(Mu_Ren**2)

  call D_cll(Dcoeff,Dcoeffuv,MomInv(1:6),masses2(0:3),rank,Derr(0:4))

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
  D0000 = Dcoeff(2,0,0,0)
  D0011 = Dcoeff(1,2,0,0)
  D0012 = Dcoeff(1,1,1,0)
  D0013 = Dcoeff(1,1,0,1)
  D0022 = Dcoeff(1,0,2,0)
  D0023 = Dcoeff(1,0,1,1)
  D0033 = Dcoeff(1,0,0,2)
  D1111 = Dcoeff(0,4,0,0)
  D1112 = Dcoeff(0,3,1,0)
  D1113 = Dcoeff(0,3,0,1)
  D1122 = Dcoeff(0,2,2,0)
  D1123 = Dcoeff(0,2,1,1)
  D1133 = Dcoeff(0,2,0,2)
  D1222 = Dcoeff(0,1,3,0)
  D1223 = Dcoeff(0,1,2,1)
  D1233 = Dcoeff(0,1,1,2)
  D1333 = Dcoeff(0,1,0,3)
  D2222 = Dcoeff(0,0,4,0)
  D2223 = Dcoeff(0,0,3,1)
  D2233 = Dcoeff(0,0,2,2)
  D2333 = Dcoeff(0,0,1,3)
  D3333 = Dcoeff(0,0,0,4)

      ggboxHH3 = &
     8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**3* &
     Spbb(1,2)**3*D123 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**3*Spbb(1,2)**3*D112 + 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**3*Spbb(1,2)**3* &
     D113 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1 &
     ,2)**3*Spbb(1,2)**3*D221 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0) &
     )/(Spaa(1,2))*Spaa(1,2)**3*Spbb(1,2)**3*D331 + 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**3*Spbb(1,2)**3* &
     D1122 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)**3*Spbb(1,2)**3*D1123 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**3*Spbb(1,2)**3*D1133 - 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(1,3)*Spbb(1,2)**2*Spbb(1,3)*D12 - 4d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(1,3)*Spbb(1,2)**2* &
     Spbb(1,3)*D13
      ggboxHH3 = ggboxHH3 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(1,3)*Spbb(1,2)**2*Spbb(1,3)*D123 -  &
     12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2 &
     *Spaa(1,3)*Spbb(1,2)**2*Spbb(1,3)*D112 - 4d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(1,3)*Spbb(1,2)**2* &
     Spbb(1,3)*D113 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,2)**2*Spaa(1,3)*Spbb(1,2)**2*Spbb(1,3)*D221 - 16d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(1,3)*Spbb(1,2)**2*Spbb(1,3)*D1122 - 16d0/(dsqrt(2d0))/(Spaa(1 &
     ,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(1,3)*Spbb(1,2)**2 &
     *Spbb(1,3)*D1123 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)**2*Spaa(1,4)*Spbb(1,2)**2*Spbb(1,4)*D12 - 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(1,4)*Spbb(1,2)**2*Spbb(1,4)*D13 - 8d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(1,4)*Spbb(1,2)**2* &
     Spbb(1,4)*D123
      ggboxHH3 = ggboxHH3 - 12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(1,4)*Spbb(1,2)**2*Spbb(1,4)*D112 -  &
     4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(1,4)*Spbb(1,2)**2*Spbb(1,4)*D113 - 8d0/(dsqrt(2d0))/(Spaa(1,2 &
     ))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(1,4)*Spbb(1,2)**2* &
     Spbb(1,4)*D221 - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)**2*Spaa(1,4)*Spbb(1,2)**2*Spbb(1,4)*D1122 - 16d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(1,4)*Spbb(1,2)**2*Spbb(1,4)*D1123 + 8d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2* &
     Spbb(2,3)*D22 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2 &
     ))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D12 + 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D13 + 8d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2* &
     Spbb(2,3)*D23
      ggboxHH3 = ggboxHH3 + 24d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D123 +  &
     4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D222 + 20d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2* &
     Spbb(2,3)*D221 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D223 + 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D331 + 4d0/(dsqrt(2d0))/(Spaa(1,2 &
     ))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2* &
     Spbb(2,3)*D332 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D1222 + 32d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,3)*Spbb(1,2)**2*Spbb(2,3)*D1223 + 16d0/(dsqrt(2d0))/(Spaa(1 &
     ,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,3)*Spbb(1,2)**2 &
     *Spbb(2,3)*D1233
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D22 + 4d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D12 + 4d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2* &
     Spbb(2,4)*D13 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2 &
     ))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D23 + 24d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D123 + 4d0/(dsqrt(2d0))/(Spaa(1,2 &
     ))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2* &
     Spbb(2,4)*D222 + 20d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D221 + 8d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D223 + 4d0/(dsqrt(2d0))/(Spaa(1,2 &
     ))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2* &
     Spbb(2,4)*D331
      ggboxHH3 = ggboxHH3 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D332 +  &
     16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2 &
     *Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D1222 + 32d0/(dsqrt(2d0))/(Spaa( &
     1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)** &
     2*Spbb(2,4)*D1223 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(2,4)*Spbb(1,2)**2*Spbb(2,4)*D1233 &
      + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2) &
     **2*Spaa(3,4)*Spbb(1,2)**2*Spbb(3,4)*D2 - 8d0/(dsqrt(2d0))/(Spaa(1 &
     ,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(3,4)*Spbb(1,2)**2 &
     *Spbb(3,4)*D22 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,2)**2*Spaa(3,4)*Spbb(1,2)**2*Spbb(3,4)*D123 - 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spaa(3,4)*Spbb(1,2)**2*Spbb(3,4)*D222 - 16d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(3,4)*Spbb(1,2)**2* &
     Spbb(3,4)*D221
      ggboxHH3 = ggboxHH3 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spaa(3,4)*Spbb(1,2)**2*Spbb(3,4)*D223 -  &
     16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2 &
     *Spaa(3,4)*Spbb(1,2)**2*Spbb(3,4)*D1222 - 16d0/(dsqrt(2d0))/(Spaa( &
     1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spaa(3,4)*Spbb(1,2)** &
     2*Spbb(3,4)*D1223 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spbb(1,2)**2*D0*mloop**2 - 8d0/(dsqrt(2d0) &
     )/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spbb(1,2)**2* &
     D1*mloop**2 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2)) &
     *Spaa(1,2)**2*Spbb(1,2)**2*D2*mloop**2 - 4d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spbb(1,2)**2*D3* &
     mloop**2 + 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)**2*Spbb(1,2)**2*D00 + 32d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spbb(1,2)**2*D001 + 24d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2* &
     Spbb(1,2)**2*D002
      ggboxHH3 = ggboxHH3 + 24d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)**2*Spbb(1,2)**2*D003 + 32d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)**2*Spbb(1,2)**2* &
     D0012 + 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)**2*Spbb(1,2)**2*D0013 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)**2*Spbb(1,2)*Spbb(1,3 &
     )**2*D12 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(1,3)**2*Spbb(1,2)*Spbb(1,3)**2*D112 + 4d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)**2* &
     Spbb(1,2)*Spbb(1,3)**2*D221 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)**2*Spbb(1,2)*Spbb(1,3)**2* &
     D1122 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(1,3)*Spaa(1,4)*Spbb(1,2)*Spbb(1,3)*Spbb(1,4)*D12 + 8d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa( &
     1,3)*Spaa(1,4)*Spbb(1,2)*Spbb(1,3)*Spbb(1,4)*D112 + 8d0/(dsqrt(2d0 &
     ))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa( &
     1,4)*Spbb(1,2)*Spbb(1,3)*Spbb(1,4)*D221
      ggboxHH3 = ggboxHH3 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(1,4)*Spbb(1,2)*Spbb(1,3)* &
     Spbb(1,4)*D1122 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,3) &
     *D222 - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(1,3)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,3)*D221 &
      - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)* &
     Spaa(1,3)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,3)*D223 - 24d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     3)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,3)*D1222 - 24d0/(dsqrt(2d0 &
     ))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa( &
     2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,3)*D1223 - 4d0/(dsqrt(2d0))/(Spaa( &
     1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)* &
     Spbb(1,2)*Spbb(1,3)*Spbb(2,4)*D2 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)* &
     Spbb(1,3)*Spbb(2,4)*D22
      ggboxHH3 = ggboxHH3 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,3)* &
     Spbb(2,4)*D12 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2 &
     ))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)* &
     D23 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1, &
     2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)*D123 - 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)*D222 - 16d0/(dsqrt(2d0) &
     )/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2 &
     ,4)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)*D221 - 4d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1 &
     ,2)*Spbb(1,3)*Spbb(2,4)*D223 - 16d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)* &
     Spbb(1,3)*Spbb(2,4)*D1222 - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1, &
     3)*Spbb(2,4)*D1223
      ggboxHH3 = ggboxHH3 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)* &
     Spbb(2,3)*D2 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2) &
     )*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)* &
     D22 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1, &
     2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D12 + 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     3)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D23 + 8d0/(dsqrt(2d0)) &
     /(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2, &
     4)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D123 - 8d0/(dsqrt(2d0))/(Spaa(1,2 &
     ))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1, &
     2)*Spbb(1,4)*Spbb(2,3)*D1222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(2,4)*Spbb(1,2)* &
     Spbb(1,4)*Spbb(2,3)*D1223 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0 &
     ))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(3,4)*Spbb(1,2)*Spbb(1,3) &
     *Spbb(3,4)*D22
      ggboxHH3 = ggboxHH3 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spaa(3,4)*Spbb(1,2)*Spbb(1,3)* &
     Spbb(3,4)*D222 + 12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)*Spaa(1,3)*Spaa(3,4)*Spbb(1,2)*Spbb(1,3)*Spbb(3,4) &
     *D221 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(1,3)*Spaa(3,4)*Spbb(1,2)*Spbb(1,3)*Spbb(3,4)* &
     D1222 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(1,3)*Spbb(1,2)*Spbb(1,3)*D0*mloop**2 + 4d0/(dsqrt(2d0)) &
     /(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spbb(1, &
     2)*Spbb(1,3)*D1*mloop**2 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0) &
     )/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spbb(1,2)*Spbb(1,3)*D2* &
     mloop**2 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(1,3)*Spbb(1,2)*Spbb(1,3)*D00 - 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spbb(1,2) &
     *Spbb(1,3)*D001 - 24d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spaa(1,3)*Spbb(1,2)*Spbb(1,3)*D002
      ggboxHH3 = ggboxHH3 - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,3)*Spbb(1,2)*Spbb(1,3)*D0012 + 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     4)**2*Spbb(1,2)*Spbb(1,4)**2*D12 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)**2*Spbb(1,2)*Spbb(1,4 &
     )**2*D112 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(1,4)**2*Spbb(1,2)*Spbb(1,4)**2*D221 + 8d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)**2* &
     Spbb(1,2)*Spbb(1,4)**2*D1122 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)* &
     Spbb(1,3)*Spbb(2,4)*D2 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0)) &
     /(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)* &
     Spbb(2,4)*D22 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2 &
     ))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)* &
     D12 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1, &
     2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)*D23
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)* &
     Spbb(2,4)*D123 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)* &
     D1222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,3)*Spbb(2,4)*D1223 - 4d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa( &
     1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D2 - 4d0/(dsqrt(2d0)) &
     /(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2, &
     3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D22 - 4d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2 &
     )*Spbb(1,4)*Spbb(2,3)*D12 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0 &
     ))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4) &
     *Spbb(2,3)*D23 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)* &
     D123
      ggboxHH3 = ggboxHH3 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4)* &
     Spbb(2,3)*D222 - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3) &
     *D221 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D223 - 16d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa( &
     1,4)*Spaa(2,3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D1222 - 16d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)* &
     Spaa(2,3)*Spbb(1,2)*Spbb(1,4)*Spbb(2,3)*D1223 - 4d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,4) &
     *Spbb(1,2)*Spbb(1,4)*Spbb(2,4)*D222 - 16d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,4)*Spbb(1,2 &
     )*Spbb(1,4)*Spbb(2,4)*D221 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,4)*Spbb(1,2)*Spbb(1, &
     4)*Spbb(2,4)*D223
      ggboxHH3 = ggboxHH3 - 24d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)* &
     Spbb(2,4)*D1222 - 24d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(2,4)*Spbb(1,2)*Spbb(1,4)*Spbb(2,4 &
     )*D1223 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(1,4)*Spaa(3,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D22 &
      + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)* &
     Spaa(1,4)*Spaa(3,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D222 + 12d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     4)*Spaa(3,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D221 + 16d0/(dsqrt(2d0) &
     )/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spaa(3 &
     ,4)*Spbb(1,2)*Spbb(1,4)*Spbb(3,4)*D1222 + 4d0/(dsqrt(2d0))/(Spaa(1 &
     ,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spbb(1,2)*Spbb( &
     1,4)*D0*mloop**2 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spaa(1,4)*Spbb(1,2)*Spbb(1,4)*D1*mloop**2 + 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     4)*Spbb(1,2)*Spbb(1,4)*D2*mloop**2
      ggboxHH3 = ggboxHH3 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spbb(1,2)*Spbb(1,4)*D00 - 8d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1, &
     4)*Spbb(1,2)*Spbb(1,4)*D001 - 24d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(1,4)*Spbb(1,2)*Spbb(1,4)* &
     D002 - 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(1,4)*Spbb(1,2)*Spbb(1,4)*D0012 + 4d0/(dsqrt(2d0))/(Spaa( &
     1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)* &
     Spbb(2,3)**2*D22 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)*Spbb(2,3)**2*D23 + 12d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2, &
     3)**2*Spbb(1,2)*Spbb(2,3)**2*D222 + 16d0/(dsqrt(2d0))/(Spaa(1,2)) &
     /(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)*Spbb(2 &
     ,3)**2*D223 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2)) &
     *Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)*Spbb(2,3)**2*D332 + 8d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)**2* &
     Spbb(1,2)*Spbb(2,3)**2*D2222
      ggboxHH3 = ggboxHH3 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(2,3)**2*Spbb(1,2)*Spbb(2,3)**2*D2223 &
      + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)* &
     Spaa(2,3)**2*Spbb(1,2)*Spbb(2,3)**2*D2233 + 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(2,4) &
     *Spbb(1,2)*Spbb(2,3)*Spbb(2,4)*D22 + 8d0/(dsqrt(2d0))/(Spaa(1,2)) &
     /(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2) &
     *Spbb(2,3)*Spbb(2,4)*D23 + 24d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0 &
     ))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(2,3) &
     *Spbb(2,4)*D222 + 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(2,3)*Spbb(2,4 &
     )*D223 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(2,3)*Spbb(2,4)*D332 &
      + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2) &
     *Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(2,3)*Spbb(2,4)*D2222 + 32d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2, &
     3)*Spaa(2,4)*Spbb(1,2)*Spbb(2,3)*Spbb(2,4)*D2223
      ggboxHH3 = ggboxHH3 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(2,4)*Spbb(1,2)*Spbb(2,3)* &
     Spbb(2,4)*D2233 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)*Spaa(2,3)*Spaa(3,4)*Spbb(1,2)*Spbb(2,3)*Spbb(3,4) &
     *D22 - 20d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(2,3)*Spaa(3,4)*Spbb(1,2)*Spbb(2,3)*Spbb(3,4)*D222 - 12d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa( &
     2,3)*Spaa(3,4)*Spbb(1,2)*Spbb(2,3)*Spbb(3,4)*D223 - 16d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)* &
     Spaa(3,4)*Spbb(1,2)*Spbb(2,3)*Spbb(3,4)*D2222 - 16d0/(dsqrt(2d0)) &
     /(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spaa(3, &
     4)*Spbb(1,2)*Spbb(2,3)*Spbb(3,4)*D2223 - 4d0/(dsqrt(2d0))/(Spaa(1, &
     2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spbb(1,2)*Spbb(2 &
     ,3)*D0*mloop**2 - 12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spaa(2,3)*Spbb(1,2)*Spbb(2,3)*D2*mloop**2 - 4d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2, &
     3)*Spbb(1,2)*Spbb(2,3)*D3*mloop**2
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spbb(1,2)*Spbb(2,3)*D00 + 40d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2, &
     3)*Spbb(1,2)*Spbb(2,3)*D002 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spbb(1,2)*Spbb(2,3)*D003 +  &
     32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)* &
     Spaa(2,3)*Spbb(1,2)*Spbb(2,3)*D0022 + 32d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,3)*Spbb(1,2)*Spbb(2,3 &
     )*D0023 + 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(2,4)**2*Spbb(1,2)*Spbb(2,4)**2*D22 + 4d0/(dsqrt(2d0 &
     ))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)**2* &
     Spbb(1,2)*Spbb(2,4)**2*D23 + 12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)**2*Spbb(1,2)*Spbb(2,4)**2* &
     D222 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(2,4)**2*Spbb(1,2)*Spbb(2,4)**2*D223 + 4d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)**2*Spbb(1 &
     ,2)*Spbb(2,4)**2*D332
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(2,4)**2*Spbb(1,2)*Spbb(2,4)**2*D2222 &
      + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2) &
     *Spaa(2,4)**2*Spbb(1,2)*Spbb(2,4)**2*D2223 + 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)**2*Spbb(1 &
     ,2)*Spbb(2,4)**2*D2233 - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0)) &
     /(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)*Spaa(3,4)*Spbb(1,2)*Spbb(2,4)* &
     Spbb(3,4)*D22 - 20d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,2)*Spaa(2,4)*Spaa(3,4)*Spbb(1,2)*Spbb(2,4)*Spbb(3,4)* &
     D222 - 12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,2)*Spaa(2,4)*Spaa(3,4)*Spbb(1,2)*Spbb(2,4)*Spbb(3,4)*D223 - 16d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa( &
     2,4)*Spaa(3,4)*Spbb(1,2)*Spbb(2,4)*Spbb(3,4)*D2222 - 16d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)* &
     Spaa(3,4)*Spbb(1,2)*Spbb(2,4)*Spbb(3,4)*D2223 - 4d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)*Spbb(1,2) &
     *Spbb(2,4)*D0*mloop**2
      ggboxHH3 = ggboxHH3 - 12d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(2,4)*Spbb(1,2)*Spbb(2,4)*D2*mloop**2 &
      - 4d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)* &
     Spaa(2,4)*Spbb(1,2)*Spbb(2,4)*D3*mloop**2 + 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)*Spbb(1,2) &
     *Spbb(2,4)*D00 + 40d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1 &
     ,2))*Spaa(1,2)*Spaa(2,4)*Spbb(1,2)*Spbb(2,4)*D002 + 8d0/(dsqrt(2d0 &
     ))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2,4)*Spbb( &
     1,2)*Spbb(2,4)*D003 + 32d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(2,4)*Spbb(1,2)*Spbb(2,4)*D0022 + 32d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(2, &
     4)*Spbb(1,2)*Spbb(2,4)*D0023 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(3,4)**2*Spbb(1,2)*Spbb(3,4 &
     )**2*D222 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,2)*Spaa(3,4)**2*Spbb(1,2)*Spbb(3,4)**2*D2222 + 8d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(3,4)* &
     Spbb(1,2)*Spbb(3,4)*D2*mloop**2
      ggboxHH3 = ggboxHH3 - 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,2)*Spaa(3,4)*Spbb(1,2)*Spbb(3,4)*D002 - 32d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spaa(3, &
     4)*Spbb(1,2)*Spbb(3,4)*D0022 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spbb(1,2)*D0*mloop**4 - 48d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,2)*Spbb(1, &
     2)*D00*mloop**2 + 96d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa( &
     1,2))*Spaa(1,2)*Spbb(1,2)*D0000 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)**2*Spaa(2,3)*Spbb(1,3)**2*Spbb( &
     2,3)*D1222 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))* &
     Spaa(1,3)**2*Spaa(2,4)*Spbb(1,3)*Spbb(1,4)*Spbb(2,3)*D1222 + 8d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa( &
     1,4)*Spaa(2,3)*Spbb(1,3)**2*Spbb(2,4)*D1222 + 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(1,4)*Spaa(2,3) &
     *Spbb(1,3)*Spbb(1,4)*Spbb(2,3)*D1222 + 8d0/(dsqrt(2d0))/(Spaa(1,2) &
     )/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(1,4)*Spaa(2,4)*Spbb(1,3 &
     )*Spbb(1,4)*Spbb(2,4)*D1222
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,3)*Spaa(1,4)*Spaa(2,4)*Spbb(1,4)**2*Spbb(2,3)* &
     D1222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,3)*Spaa(2,3)**2*Spbb(1,3)*Spbb(2,3)**2*D222 - 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)**2*Spbb(1 &
     ,3)*Spbb(2,3)**2*D2222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0)) &
     /(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)**2*Spbb(1,3)*Spbb(2,3)**2*D2223 &
      - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)* &
     Spaa(2,3)*Spaa(2,4)*Spbb(1,3)*Spbb(2,3)*Spbb(2,4)*D222 - 8d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2, &
     3)*Spaa(2,4)*Spbb(1,3)*Spbb(2,3)*Spbb(2,4)*D2222 - 8d0/(dsqrt(2d0) &
     )/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2 &
     ,4)*Spbb(1,3)*Spbb(2,3)*Spbb(2,4)*D2223 - 8d0/(dsqrt(2d0))/(Spaa(1 &
     ,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2,4)*Spbb( &
     1,4)*Spbb(2,3)**2*D222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0)) &
     /(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,3)** &
     2*D2222
      ggboxHH3 = ggboxHH3 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,3)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,3)**2* &
     D2223 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,3)*Spaa(2,3)*Spaa(3,4)*Spbb(1,3)*Spbb(2,3)*Spbb(3,4)*D222 + 8d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa( &
     2,3)*Spaa(3,4)*Spbb(1,3)*Spbb(2,3)*Spbb(3,4)*D2222 - 8d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)* &
     Spbb(1,3)*Spbb(2,3)*D2*mloop**2 - 24d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)*Spbb(1,3)*Spbb(2,3)* &
     D22*mloop**2 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2 &
     ))*Spaa(1,3)*Spaa(2,3)*Spbb(1,3)*Spbb(2,3)*D002 - 16d0/(dsqrt(2d0) &
     )/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,3)*Spbb(1 &
     ,3)*Spbb(2,3)*D0022 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,3)*Spaa(2,4)**2*Spbb(1,4)*Spbb(2,3)*Spbb(2,4)* &
     D222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1 &
     ,3)*Spaa(2,4)**2*Spbb(1,4)*Spbb(2,3)*Spbb(2,4)*D2222
      ggboxHH3 = ggboxHH3 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,3)*Spaa(2,4)**2*Spbb(1,4)*Spbb(2,3)*Spbb(2,4)* &
     D2223 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,3)*Spaa(2,4)*Spaa(3,4)*Spbb(1,4)*Spbb(2,3)*Spbb(3,4)*D222 + 8d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa( &
     2,4)*Spaa(3,4)*Spbb(1,4)*Spbb(2,3)*Spbb(3,4)*D2222 - 8d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,4)* &
     Spbb(1,4)*Spbb(2,3)*D2*mloop**2 - 24d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,3)* &
     D22*mloop**2 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2 &
     ))*Spaa(1,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,3)*D002 - 16d0/(dsqrt(2d0) &
     )/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,3)*Spaa(2,4)*Spbb(1 &
     ,4)*Spbb(2,3)*D0022 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,4)**2*Spaa(2,3)*Spbb(1,3)*Spbb(1,4)*Spbb(2,4)* &
     D1222 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,4)**2*Spaa(2,4)*Spbb(1,4)**2*Spbb(2,4)*D1222
      ggboxHH3 = ggboxHH3 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,4)*Spaa(2,3)**2*Spbb(1,3)*Spbb(2,3)*Spbb(2,4)* &
     D222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1 &
     ,4)*Spaa(2,3)**2*Spbb(1,3)*Spbb(2,3)*Spbb(2,4)*D2222 - 8d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2, &
     3)**2*Spbb(1,3)*Spbb(2,3)*Spbb(2,4)*D2223 - 8d0/(dsqrt(2d0))/( &
     Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,3)*Spaa(2,4) &
     *Spbb(1,3)*Spbb(2,4)**2*D222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/( &
     dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,3)* &
     Spbb(2,4)**2*D2222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,3)*Spbb(2,4)**2* &
     D2223 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,4)*Spaa(2,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,3)*Spbb(2,4)*D222 - 8d0 &
     /(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa( &
     2,3)*Spaa(2,4)*Spbb(1,4)*Spbb(2,3)*Spbb(2,4)*D2222 - 8d0/(dsqrt( &
     2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,3)* &
     Spaa(2,4)*Spbb(1,4)*Spbb(2,3)*Spbb(2,4)*D2223
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,4)*Spaa(2,3)*Spaa(3,4)*Spbb(1,3)*Spbb(2,4)* &
     Spbb(3,4)*D222 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,4)*Spaa(2,3)*Spaa(3,4)*Spbb(1,3)*Spbb(2,4)*Spbb(3,4)* &
     D2222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,4)*Spaa(2,3)*Spbb(1,3)*Spbb(2,4)*D2*mloop**2 - 24d0/(dsqrt(2d0)) &
     /(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,3)*Spbb(1, &
     3)*Spbb(2,4)*D22*mloop**2 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,3)*Spbb(1,3)*Spbb(2,4)*D002 -  &
     16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)* &
     Spaa(2,3)*Spbb(1,3)*Spbb(2,4)*D0022 - 8d0/(dsqrt(2d0))/(Spaa(1,2)) &
     /(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,4)**2*Spbb(1,4)*Spbb(2 &
     ,4)**2*D222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2)) &
     *Spaa(1,4)*Spaa(2,4)**2*Spbb(1,4)*Spbb(2,4)**2*D2222 - 8d0/( &
     dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2, &
     4)**2*Spbb(1,4)*Spbb(2,4)**2*D2223
      ggboxHH3 = ggboxHH3 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/( &
     Spaa(1,2))*Spaa(1,4)*Spaa(2,4)*Spaa(3,4)*Spbb(1,4)*Spbb(2,4)* &
     Spbb(3,4)*D222 + 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1, &
     2))*Spaa(1,4)*Spaa(2,4)*Spaa(3,4)*Spbb(1,4)*Spbb(2,4)*Spbb(3,4)* &
     D2222 - 8d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa( &
     1,4)*Spaa(2,4)*Spbb(1,4)*Spbb(2,4)*D2*mloop**2 - 24d0/(dsqrt(2d0)) &
     /(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,4)*Spbb(1, &
     4)*Spbb(2,4)*D22*mloop**2 + 16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt( &
     2d0))/(Spaa(1,2))*Spaa(1,4)*Spaa(2,4)*Spbb(1,4)*Spbb(2,4)*D002 -  &
     16d0/(dsqrt(2d0))/(Spaa(1,2))/(dsqrt(2d0))/(Spaa(1,2))*Spaa(1,4)* &
     Spaa(2,4)*Spbb(1,4)*Spbb(2,4)*D0022
  
  return
end subroutine ggboxHH3pp

end module ModggboxHH3pp
!!--YaofuZhou-----------------------------------------