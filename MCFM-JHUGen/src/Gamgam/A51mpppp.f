      function A51mpppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A51mpppp
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5
C----Eq.(4) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(-i)*(-1)
C--- to give (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      A51mpppp=-((s(j2,j3)+s(j3,j4)+s(j4,j5))*zb(j2,j5)**2
     & -zb(j2,j4)*za(j4,j3)*zb(j3,j5)*zb(j2,j5)
     & -zb(j1,j2)*zb(j1,j5)/(za(j1,j2)*za(j1,j5))
     & *(za(j1,j2)**2*za(j1,j3)**2*zb(j2,j3)/za(j2,j3)
     &  +za(j1,j3)**2*za(j1,j4)**2*zb(j3,j4)/za(j3,j4)
     &  +za(j1,j4)**2*za(j1,j5)**2*zb(j4,j5)/za(j4,j5)))
     & /(3._dp*zb(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))
      return
      end
