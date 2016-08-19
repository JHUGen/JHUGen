c---Note: These routines are not used within this code
c---      to compute the A4;3 amplitude, instead it is
c---      calculated directly
      double complex function A1phiAgQgmmppL(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.32
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer j1,j2,j3,j4
      double complex zab2,V3L,lnrat,L0,L1,Lsm1DS,Lsm1_2me,
     . l341,l34,l41,A0phiAgQgmmpp,A0phidAgQgmmpp
      double precision s123,s234,s412,s341,mhsq
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
     . -epinv**2-epinv*l34-0.5d0*l34**2
     . -epinv**2-epinv*l41-0.5d0*l41**2
     . -3d0/2d0*(epinv+l341)-7d0/2d0-deltar/2d0
     . -Lsm1_2me(s412,s123,s(j1,j2),mhsq)
     . -Lsm1_2me(s123,s234,s(j2,j3),mhsq)
      
      A1phiAgQgmmppL=
     . -0.5d0*za(j2,j4)**2*zab2(j1,j2,j3,j4)**2/(za(j3,j4)*za(j1,j4))
     .  *L1(-s234,-s(j2,j3))/s(j2,j3)**2
     . +0.5d0*za(j2,j3)**2*zb(j3,j4)**2*za(j1,j4)/za(j3,j4)
     .  *L1(-s341,-s(j4,j1))/s(j4,j1)**2
     . +2d0*za(j1,j2)*za(j2,j3)*zb(j3,j4)/za(j3,j4)
     .  *(L0(-s234,-s(j2,j3))/s(j2,j3)+L0(-s341,-s(j4,j1))/s(j4,j1))
     . +0.5d0*A0phiAgQgmmpp(j1,j2,j3,j4,za,zb)*lnrat(-s234,-s(j2,j3))
     . -0.5d0*zab2(j2,j1,j3,j4)**2/(s341*zb(j1,j4)*za(j3,j4))
     . +0.5d0*s123*zb(j1,j3)*zb(j3,j4)
     .  /(zb(j1,j2)*zb(j2,j3)*zab2(j4,j2,j3,j1))
     . -0.5d0*za(j1,j2)**2*za(j2,j4)*(s(j1,j4)+s(j2,j4)+s(j3,j4))
     .  /(za(j1,j4)*za(j2,j3)*za(j3,j4)*zab2(j4,j1,j2,j3))
     . -0.5d0*(s123**2*za(j2,j4)**2*zb(j3,j4))
     .  /(s(j2,j3)*za(j3,j4)*zab2(j4,j2,j3,j1)*zab2(j4,j1,j2,j3))
     . -2d0*A0phidAgQgmmpp(j1,j2,j3,j4,za,zb)
      
      A1phiAgQgmmppL=A1phiAgQgmmppL
     . +A0phiAgQgmmpp(j1,j2,j3,j4,za,zb)*(
     .   V3L
     .  -Lsm1DS(s(j3,j4),s(j4,j1),s341)
     .  -Lsm1DS(s(j1,j2),s(j2,j3),s123))
      
      return
      end
      
      double complex function A1phiAgQgmmppR(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1 reflection relation (4.2.3)
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A1phiAgQgmppmL
 
      A1phiAgQgmmppR=A1phiAgQgmppmL(j1,j4,j3,j2,za,zb)
      
      return
      end
      
      double complex function A1phiAgQgmmppf(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.34
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
 
      A1phiAgQgmmppf=czip
      
      return
      end
      
