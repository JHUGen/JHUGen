c---Note: These routines are not used within this code
c---      to compute the A4;3 amplitude, instead it is
c---      calculated directly
      function A1phiAgQgmmppL(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiAgQgmmppL

C     implementation of arXiv:0906.0008v1, Eq. 4.32
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2,V3L,lnrat,L0,L1,Lsm1DS,Lsm1_2me,
     & l341,l34,l41,A0phiAgQgmmpp,A0phidAgQgmmpp
      real(dp):: s123,s234,s412,s341,mhsq
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      s412=s(j4,j1)+s(j4,j2)+s(j1,j2)
      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s341=s(j3,j4)+s(j3,j1)+s(j4,j1)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      l341=lnrat(musq,-s341)
      l34=lnrat(musq,-s(j3,j4))
      l41=lnrat(musq,-s(j4,j1))

      V3L=
     & -epinv**2-epinv*l34-0.5_dp*l34**2
     & -epinv**2-epinv*l41-0.5_dp*l41**2
     & -3._dp/2._dp*(epinv+l341)-7._dp/2._dp-deltar/2._dp
     & -Lsm1_2me(s412,s123,s(j1,j2),mhsq)
     & -Lsm1_2me(s123,s234,s(j2,j3),mhsq)

      A1phiAgQgmmppL=
     & -0.5_dp*za(j2,j4)**2*zab2(j1,j2,j3,j4)**2/(za(j3,j4)*za(j1,j4))
     &  *L1(-s234,-s(j2,j3))/s(j2,j3)**2
     & +0.5_dp*za(j2,j3)**2*zb(j3,j4)**2*za(j1,j4)/za(j3,j4)
     &  *L1(-s341,-s(j4,j1))/s(j4,j1)**2
     & +2._dp*za(j1,j2)*za(j2,j3)*zb(j3,j4)/za(j3,j4)
     &  *(L0(-s234,-s(j2,j3))/s(j2,j3)+L0(-s341,-s(j4,j1))/s(j4,j1))
     & +0.5_dp*A0phiAgQgmmpp(j1,j2,j3,j4,za,zb)*lnrat(-s234,-s(j2,j3))
     & -0.5_dp*zab2(j2,j1,j3,j4)**2/(s341*zb(j1,j4)*za(j3,j4))
     & +0.5_dp*s123*zb(j1,j3)*zb(j3,j4)
     &  /(zb(j1,j2)*zb(j2,j3)*zab2(j4,j2,j3,j1))
     & -0.5_dp*za(j1,j2)**2*za(j2,j4)*(s(j1,j4)+s(j2,j4)+s(j3,j4))
     &  /(za(j1,j4)*za(j2,j3)*za(j3,j4)*zab2(j4,j1,j2,j3))
     & -0.5_dp*(s123**2*za(j2,j4)**2*zb(j3,j4))
     &  /(s(j2,j3)*za(j3,j4)*zab2(j4,j2,j3,j1)*zab2(j4,j1,j2,j3))
     & -2._dp*A0phidAgQgmmpp(j1,j2,j3,j4,za,zb)

      A1phiAgQgmmppL=A1phiAgQgmmppL
     & +A0phiAgQgmmpp(j1,j2,j3,j4,za,zb)*(
     &   V3L
     &  -Lsm1DS(s(j3,j4),s(j4,j1),s341)
     &  -Lsm1DS(s(j1,j2),s(j2,j3),s123))

      return
      end

      function A1phiAgQgmmppR(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiAgQgmmppR

C     implementation of arXiv:0906.0008v1 reflection relation (4.2.3)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiAgQgmppmL

      A1phiAgQgmmppR=A1phiAgQgmppmL(j1,j4,j3,j2,za,zb)

      return
      end

      function A1phiAgQgmmppf(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiAgQgmmppf

C     implementation of arXiv:0906.0008v1, Eq. 4.34
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4

      A1phiAgQgmmppf=czip

      return
      end

