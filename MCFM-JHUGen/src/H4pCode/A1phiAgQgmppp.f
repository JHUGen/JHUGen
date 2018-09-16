c---Note: These routines are not used within this code
c---      to compute the A4;3 amplitude, instead it is
c---      calculated directly
      function A1phiagqgmpppL(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiagqgmpppL

C     implementation of arXiv:0906.0008v1, Eq. 4.27
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      real(dp):: s123
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)

      A1phiagqgmpppL=0.5_dp*(
     & za(j1,j3)*zab2(j1,j3,j4,j2)/(za(j2,j3)*za(j3,j4)*za(j4,j1))
     .+za(j1,j3)**2*zb(j3,j4)/(za(j1,j2)*za(j2,j3)*za(j3,j4))
     &                    )
     &             -2._dp*zab2(j1,j2,j3,j4)**2/(za(j1,j2)*za(j2,j3)*s123)

      return
      end

      function A1phiagqgmpppR(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiagqgmpppR

C     implementation of arXiv:0906.0008v1 reflection relation (4.2.3)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiagqgmpppL

      A1phiagqgmpppR=A1phiagqgmpppL(j1,j4,j3,j2,za,zb)

      return
      end

      function A1phiagqgmpppf(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiagqgmpppf

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
C     implementation of arXiv:0906.0008v1, Eq. 4.28
      integer:: j1,j2,j3,j4

      A1phiagqgmpppf=czip

      return
      end

