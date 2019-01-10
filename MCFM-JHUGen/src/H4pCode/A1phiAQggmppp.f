      double complex function A1phiAQggmpppL(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.16
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double complex zab2
      double precision s3,mhsq
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      A1phiAQggmpppL=
     . 0.5d0*za(j1,j2)*zab2(j1,j3,j4,j2)/(za(j2,j3)*za(j3,j4)*za(j4,j1))
     . +0.5d0*za(j1,j3)*zb(j3,j4)/(za(j2,j3)*za(j3,j4))
     . +2d0*zab2(j1,j3,j4,j2)**2/(za(j3,j4)*za(j4,j1)*zab2(j3,j1,j4,j2))
     . -2d0*zab2(j1,j2,j3,j4)**2*zab2(j2,j1,j3,j4)
     . /(za(j1,j2)*za(j2,j3)*s3(j1,j2,j3)*zab2(j3,j1,j2,j4))
     . -2d0*zb(j2,j4)**3*mhsq**2
     . /(zb(j1,j2)*s3(j4,j1,j2)*zab2(j3,j1,j2,j4)*zab2(j3,j1,j4,j2))
     . -1d0/3d0*za(j1,j3)*zb(j3,j4)*za(j4,j1)/(za(j1,j2)*za(j3,j4)**2)
      return
      end

      double complex function A1phiAQggmpppR(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.17
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex zab2
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      A1phiAQggmpppR=
     . -0.5d0*zab2(j1,j2,j3,j4)/(za(j2,j3)*za(j3,j4))
     . -0.5d0*za(j1,j2)*zb(j2,j3)*za(j3,j1)
     .  /(za(j2,j3)*za(j3,j4)*za(j4,j1))
      return
      end

      double complex function A1phiAQggmpppF(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.18
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      A1phiAQggmpppF=
     . +1d0/3d0*za(j1,j3)*zb(j3,j4)*za(j4,j1)/(za(j1,j2)*za(j3,j4)**2)
      return
      end
