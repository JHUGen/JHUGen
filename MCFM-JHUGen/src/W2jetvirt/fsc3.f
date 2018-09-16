      function Fsc3(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Fsc3
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L1
      real(dp):: t
      Fsc3=
     .-((za(j4,j5)**2*((L1(-s(j5,j6),-t(j2,j3,j4))*
     .(za(j1,j5)*za(j2,j4)*zb(j1,j2)+
     .za(j1,j5)*za(j3,j4)*zb(j1,j3))**2)/
     .(2._dp*za(j4,j5)**2*t(j2,j3,j4)**2)-
     .(L1(-s(j3,j4),-t(j2,j3,j4))*za(j2,j5)**2*za(j3,j4)**2*zb(j2,j3)**2
     .)/
     .(2._dp*za(j4,j5)**2*t(j2,j3,j4)**2)))/
     .(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6)))+
     .(-((za(j4,j5)*(za(j1,j3)*zb(j1,j3)*zb(j1,j6)+
     .za(j2,j3)*zb(j1,j3)*zb(j2,j6)))/t(j1,j2,j3))+
     .(zb(j1,j6)*(-(za(j2,j5)*za(j3,j4)*zb(j2,j3))+
     .za(j3,j4)*za(j4,j5)*zb(j3,j4)))/t(j2,j3,j4))/
     .(3._dp*s(j5,j6)*za(j2,j3)**2)+(-((za(j2,j5)*za(j4,j5)*zb(j2,j3))/
     .(za(j1,j2)*za(j2,j3)*za(j5,j6)*t(j2,j3,j4)))+
     .(-(za(j2,j4)*zb(j2,j6))-za(j3,j4)*zb(j3,j6))**2/
     .(za(j1,j2)*za(j2,j3)*za(j3,j4)*zb(j5,j6)*t(j2,j3,j4))-
     .(za(j4,j5)*zb(j2,j3)*(-(za(j2,j5)*zb(j2,j3))+za(j4,j5)*zb(j3,j4)))
     ./
     .(za(j1,j2)*za(j5,j6)*t(j1,j2,j3)*t(j2,j3,j4))+
     .(za(j2,j4)*zb(j1,j6)*zb(j2,j3)*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6)))/
     .(za(j1,j2)*za(j2,j3)*zb(j5,j6)*t(j1,j2,j3)*t(j2,j3,j4)))/2
      return
      end


