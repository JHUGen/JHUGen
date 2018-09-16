      function Fcc_qpgpgmqm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Fcc_qpgpgmqm

      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: L0,Lsm1,Lsm1_2mht,I3m
      real(dp):: t

      Fcc_qpgpgmqm=
     .(I3m(s(j1,j2),s(j3,j4),s(j5,j6))*zb(j1,j2)*
     .(-(s(j1,j2)*s(j3,j4)*za(j4,j5)**2)+
     .(-(za(j1,j4)*(za(j3,j5)*zb(j1,j3)+za(j4,j5)*zb(j1,j4)))-
     .za(j2,j4)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))**2))/
     .(2._dp*za(j1,j2)*za(j5,j6)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))
     .*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j3,j4))+
     .(I3m(s(j3,j4),s(j1,j2),s(j5,j6))*za(j3,j4)*
     .(-(s(j1,j2)*s(j3,j4)*zb(j1,j6)**2)+
     .(zb(j1,j3)*(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))+
     .zb(j1,j4)*(-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6)))**2))/
     .(2._dp*za(j1,j2)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j3,j4)*zb(j5,j6))-
     .(2._dp*za(j1,j3)*(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))*
     .((L0(-s(j2,j3),-t(j1,j2,j3))*
     .(za(j1,j2)*zb(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j1,j2)*zb(j3,j6))
     .)/t(j1,j2,j3)**2+(L0(-s(j5,j6),-t(j1,j2,j3))*za(j3,j4)*zb(j4,j6))/
     .(za(j2,j3)*t(j1,j2,j3))))/
     .(za(j1,j2)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))+
     .Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))*
     .((za(j4,j5)**2*zb(j1,j2)**3)/
     .(za(j5,j6)*zb(j1,j3)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .zb(j2,j3)*t(j1,j2,j3))+
     .(za(j1,j3)*(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))**2)/
     .(za(j1,j2)*za(j2,j3)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .zb(j5,j6)*t(j1,j2,j3)))+
     .Lsm1_2mht(s(j3,j4),t(j1,j2,j3),s(j5,j6),s(j1,j2))*
     .((za(j4,j5)**2*zb(j1,j2)**2*(-(za(j1,j4)*zb(j1,j2))+
     .za(j3,j4)*zb(j2,j3)))/
     .(za(j5,j6)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j2,j3)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*t(j1,j2,j3))+
     .(za(j1,j3)*(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))**2)/
     .(za(j1,j2)*za(j2,j3)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .zb(j5,j6)*t(j1,j2,j3)))+
     .(2._dp*zb(j2,j4)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))*
     .((L0(-s(j2,j3),-t(j2,j3,j4))*
     .(za(j2,j5)*za(j3,j4)*zb(j2,j4)+za(j3,j4)*za(j3,j5)*zb(j3,j4))
     .)/t(j2,j3,j4)**2+(L0(-s(j5,j6),-t(j2,j3,j4))*za(j1,j5)*zb(j1,j2))/
     .(zb(j2,j3)*t(j2,j3,j4))))/
     .(za(j5,j6)*zb(j3,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))+
     .Lsm1(-s(j3,j4),-t(j2,j3,j4),-s(j2,j3),-t(j2,j3,j4))*
     .((zb(j2,j4)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))**2)/
     .(za(j5,j6)*zb(j2,j3)*zb(j3,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*t(j2,j3,j4))+
     .(za(j3,j4)**3*zb(j1,j6)**2)/
     .(za(j2,j3)*za(j2,j4)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .zb(j5,j6)*t(j2,j3,j4)))+
     .Lsm1_2mht(s(j1,j2),t(j2,j3,j4),s(j5,j6),s(j3,j4))*
     .((zb(j2,j4)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))**2)/
     .(za(j5,j6)*zb(j2,j3)*zb(j3,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*t(j2,j3,j4))+
     .(za(j3,j4)**2*(za(j2,j3)*zb(j1,j2)-za(j3,j4)*zb(j1,j4))*zb(j1,j6)*
     .*2)/
     .(za(j2,j3)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j5,j6)*t(j2,j3,j4)
     .))

      return
      end

