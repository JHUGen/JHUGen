      function Fcc_qpgpgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Fcc_qpgpgpqm

      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: L0,Lsm1,Lsm1_2me
      real(dp):: t

      Fcc_qpgpgpqm=
     .-(((-Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))-
     .Lsm1(-s(j2,j3),-t(j2,j3,j4),-s(j3,j4),-t(j2,j3,j4))-
     .Lsm1_2me(t(j1,j2,j3),t(j2,j3,j4),s(j2,j3),s(j5,j6)))*za(j4,j5)**2)
     ./
     .(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6)))-
     .(2._dp*L0(-t(j2,j3,j4),-s(j3,j4))*za(j2,j5)*za(j4,j5)*zb(j2,j3))/
     .(s(j3,j4)*za(j1,j2)*za(j2,j3)*za(j5,j6))-
     .(2._dp*L0(-s(j5,j6),-t(j2,j3,j4))*za(j1,j5)*za(j4,j5)*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3)))/
     .(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6)*t(j2,j3,j4))

      return
      end

