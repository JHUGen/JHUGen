      function A1phiAQggmpppL(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiAQggmpppL

C     implementation of arXiv:0906.0008v1, Eq. 4.16
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      real(dp):: s3,mhsq
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      A1phiAQggmpppL=
     & 0.5_dp*za(j1,j2)*zab2(j1,j3,j4,j2)/(za(j2,j3)*za(j3,j4)*za(j4,j1))
     & +0.5_dp*za(j1,j3)*zb(j3,j4)/(za(j2,j3)*za(j3,j4))
     & +2._dp*zab2(j1,j3,j4,j2)**2/(za(j3,j4)*za(j4,j1)*zab2(j3,j1,j4,j2))
     & -2._dp*zab2(j1,j2,j3,j4)**2*zab2(j2,j1,j3,j4)
     & /(za(j1,j2)*za(j2,j3)*s3(j1,j2,j3)*zab2(j3,j1,j2,j4))
     & -2._dp*zb(j2,j4)**3*mhsq**2
     & /(zb(j1,j2)*s3(j4,j1,j2)*zab2(j3,j1,j2,j4)*zab2(j3,j1,j4,j2))
     & -1._dp/3._dp*za(j1,j3)*zb(j3,j4)*za(j4,j1)/(za(j1,j2)*za(j3,j4)**2)
      return
      end

      function A1phiAQggmpppR(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiAQggmpppR

C     implementation of arXiv:0906.0008v1, Eq. 4.17
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      A1phiAQggmpppR=
     & -0.5_dp*zab2(j1,j2,j3,j4)/(za(j2,j3)*za(j3,j4))
     & -0.5_dp*za(j1,j2)*zb(j2,j3)*za(j3,j1)
     &  /(za(j2,j3)*za(j3,j4)*za(j4,j1))
      return
      end

      function A1phiAQggmpppF(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiAQggmpppF

C     implementation of arXiv:0906.0008v1, Eq. 4.18
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A1phiAQggmpppF=
     & +1._dp/3._dp*za(j1,j3)*zb(j3,j4)*za(j4,j1)/(za(j1,j2)*za(j3,j4)**2)
      return
      end
