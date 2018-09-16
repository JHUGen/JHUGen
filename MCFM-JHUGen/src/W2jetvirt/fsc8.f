      function Fsc8(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Fsc8
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1,Lsm1_2me,Lnrat
      real(dp):: t
      Fsc8=
     .-((za(j1,j5)*za(j2,j5)*((Lnrat(-s(j2,j3),-s(j5,j6))*za(j1,j2))/
     .(za(j1,j4)*za(j2,j3))+Lnrat(-s(j1,j2),-s(j5,j6))/za(j3,j4)))/
     .(za(j1,j3)*za(j1,j4)*za(j5,j6)))+
     .(Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))*za(j1,j5)**2*
     .za(j2,j3))/
     .(za(j1,j3)**2*za(j1,j4)*za(j3,j4)*za(j5,j6))+
     .(Lsm1_2me(t(j1,j2,j3),t(j2,j3,j4),s(j2,j3),s(j5,j6))*za(j1,j5)**2*
     .za(j2,j4)**2)/(za(j1,j4)**3*za(j2,j3)*za(j3,j4)*za(j5,j6))+
     .(Lnrat(-t(j1,j2,j3),-s(j5,j6))*za(j2,j5)**2)/
     .(za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))+
     .(Lnrat(-t(j2,j3,j4),-s(j5,j6))*za(j2,j5)**2)/
     .(2._dp*za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))+
     .(Lsm1_2me(t(j1,j2,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*za(j2,j3)*
     .za(j4,j5)**2)/(za(j1,j4)*za(j3,j4)**3*za(j5,j6))-
     .(L0(-t(j1,j2,j3),-s(j2,j3))*za(j1,j2)**2*
     .(za(j1,j5)*za(j2,j5)*zb(j1,j2)+za(j1,j5)*za(j3,j5)*zb(j1,j3)))/
     .(s(j2,j3)*za(j1,j3)*za(j1,j4)**2*za(j2,j3)*za(j5,j6))+
     .(L1(-t(j2,j3,j4),-s(j5,j6))*za(j1,j2)**2*za(j5,j6)*zb(j1,j6)**2)/
     .(2._dp*s(j5,j6)**2*za(j1,j4)*za(j2,j3)*za(j3,j4))+
     .(L0(-t(j2,j3,j4),-s(j5,j6))*((za(j1,j5)**2*za(j2,j4)*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4)))/
     .(za(j1,j4)**2*za(j2,j3)*za(j3,j4)*za(j5,j6))+
     .(za(j1,j2)*za(j2,j5)*zb(j1,j6))/(za(j1,j4)*za(j2,j3)*za(j3,j4))))/
     .s(j5,j6)+(L0(-t(j1,j2,j3),-s(j1,j2))*za(j2,j3)*
     .(-(za(j1,j5)*za(j3,j5)*zb(j1,j3))-za(j2,j5)*za(j3,j5)*zb(j2,j3)))/
     .(s(j1,j2)*za(j1,j3)*za(j3,j4)**2*za(j5,j6))-
     .(L0(-t(j1,j2,j4),-s(j1,j2))*za(j2,j4)*za(j4,j5)*
     .(-(za(j1,j5)*zb(j1,j4))-za(j2,j5)*zb(j2,j4)))/
     .(s(j1,j2)*za(j1,j4)*za(j3,j4)**2*za(j5,j6))+
     .(L0(-t(j2,j3,j4),-s(j2,j3))*(za(j1,j5)*za(j2,j4)+za(j1,j4)*za(j2,j
     .5))*
     .za(j4,j5)*zb(j3,j4))/(s(j2,j3)*za(j1,j4)**2*za(j3,j4)*za(j5,j6))-
     .(L0(-t(j1,j2,j4),-s(j5,j6))*za(j2,j4)*za(j3,j5)*zb(j3,j6))/
     .(s(j5,j6)*za(j1,j4)*za(j3,j4)**2)+
     .(L1(-t(j1,j2,j4),-s(j5,j6))*za(j2,j3)*za(j5,j6)*zb(j3,j6)**2)/
     .(s(j5,j6)**2*za(j1,j4)*za(j3,j4))+
     .(L1(-t(j1,j2,j3),-s(j5,j6))*za(j2,j4)**2*za(j5,j6)*zb(j4,j6)**2)/
     .(s(j5,j6)**2*za(j1,j4)*za(j2,j3)*za(j3,j4))+
     .(L0(-t(j1,j2,j3),-s(j5,j6))*(-((za(j2,j5)*za(j4,j5)*
     .(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4)))/
     .(za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6)))-
     .(za(j2,j4)*(-(za(j1,j5)*za(j2,j4))+2._dp*za(j1,j4)*za(j2,j5))*
     .zb(j4,j6))/(za(j1,j4)**2*za(j2,j3)*za(j3,j4))+
     .(za(j2,j4)*za(j4,j5)*zb(j4,j6))/(za(j1,j4)*za(j3,j4)**2)))/s(j5,j6
     .)+
     .(-(za(j2,j5)**2/(za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6)))+
     .(za(j2,j5)*za(j4,j5)*zb(j3,j4))/
     .(za(j1,j4)*za(j3,j4)*za(j5,j6)*t(j2,j3,j4))-
     .(za(j1,j5)*za(j2,j5)*zb(j1,j3)*zb(j3,j4))/
     .(za(j1,j4)*za(j5,j6)*t(j1,j3,j4)*t(j2,j3,j4))-
     .(zb(j3,j4)*(-(za(j1,j4)*zb(j1,j6))-za(j3,j4)*zb(j3,j6))*
     .(za(j2,j3)*zb(j3,j6)+za(j2,j4)*zb(j4,j6)))/
     .(za(j1,j4)*za(j3,j4)*zb(j5,j6)*t(j1,j3,j4)*t(j2,j3,j4)))/2-
     .(L1(-s(j2,j3),-t(j2,j3,j4))*za(j2,j3)*za(j4,j5)**2*zb(j3,j4)**2)/
     .(2._dp*za(j1,j4)*za(j3,j4)*za(j5,j6)*t(j2,j3,j4)**2)


      return
      end
