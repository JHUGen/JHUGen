      function A0phiggggpmmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiggggpmmm

C----Expresssion of Eq. (3.8) of hep-th/0411092v2
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      real(dp):: s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      A0phiggggpmmm=
     .+(zab2(j3,j2,j4,j1)*za(j2,j4))**2/(s3(j1,j2,j4)*s(j1,j2)*s(j1,j4))
     .+(zab2(j4,j2,j3,j1)*za(j2,j3))**2/(s3(j1,j2,j3)*s(j1,j2)*s(j2,j3))
     .+(zab2(j2,j3,j4,j1)*za(j3,j4))**2/(s3(j1,j3,j4)*s(j1,j4)*s(j3,j4))
     .-za(j2,j4)/(za(j1,j2)*zb(j2,j3)*zb(j3,j4)*za(j4,j1))
     .*(-s(j2,j3)*zab2(j2,j3,j4,j1)/zb(j4,j1)
     &  -s(j3,j4)*zab2(j4,j2,j3,j1)/zb(j1,j2)
     .-s3(j2,j3,j4)*za(j2,j4))
      return
      end

