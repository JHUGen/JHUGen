c---Note: These routines are not used within this code
c---      to compute the A4;3 amplitude, instead it is
c---      calculated directly
      double complex function A1phiAgQgmppmL(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.27
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer j1,j2,j3,j4
      double complex zab2,V2L,l123,l34,l41,lnrat,L0,L1,Lsm1DS,
     . Lsm1_2me,A0phiAgQgmppm
      double precision s123,s234,s412,s341,mhsq
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      s412=s(j4,j1)+s(j4,j2)+s(j1,j2)
      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s341=s(j3,j4)+s(j3,j1)+s(j4,j1)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)      
      l123=lnrat(musq,-s123)
      l34=lnrat(musq,-s(j3,j4))
      l41=lnrat(musq,-s(j4,j1))

      V2L=
     . -epinv**2-epinv*l34-0.5d0*l34**2
     . -epinv**2-epinv*l41-0.5d0*l41**2
     . +13d0/6d0*(epinv+l123)
     . +119d0/18d0-deltar/6d0
     . -Lsm1_2me(s412,s123,s(j1,j2),mhsq)
     . -Lsm1_2me(s123,s234,s(j2,j3),mhsq)

      A1phiAgQgmppmL=
     . +0.5d0*za(j1,j3)**2*zab2(j4,j1,j2,j3)**2/(za(j1,j2)*za(j2,j3))
     .  *L1(-s123,-s(j1,j2))/s(j1,j2)**2
     . -0.5d0*za(j1,j2)*za(j3,j4)**2*zb(j2,j3)**2/za(j2,j3)
     .  *L1(-s234,-s(j3,j4))/s(j3,j4)**2
     . -2d0*za(j1,j4)*zb(j2,j3)*za(j3,j4)/za(j2,j3)
     .  *(L0(-s123,-s(j1,j2))/s(j1,j2)+L0(-s234,-s(j3,j4))/s(j3,j4))
     . +0.5d0*A0phiAgQgmppm(j1,j2,j3,j4,za,zb)*lnrat(-s123,-s(j1,j2))
     . -0.5d0*s341*zb(j2,j3)*zb(j1,j3)
     .  /(zb(j3,j4)*zb(j1,j4)*zab2(j2,j3,j4,j1))
     . -0.5d0*zab2(j4,j1,j3,j2)**2/(s123*zb(j1,j2)*za(j2,j3))
     . -za(j1,j4)*zb(j2,j3)*za(j3,j4)/(s(j1,j2)*za(j2,j3))
     . +0.5d0*za(j1,j4)**2*(s(j1,j3)+s(j2,j3))
     .  /(s(j1,j2)*za(j1,j2)*za(j2,j3))
     . -0.5d0*zb(j2,j3)*za(j3,j4)*zab2(j2,j1,j4,j3)
     .  /(za(j2,j3)*zb(j3,j4)*zab2(j2,j3,j4,j1))
            
      A1phiAgQgmppmL=A1phiAgQgmppmL
     . +A0phiAgQgmppm(j1,j2,j3,j4,za,zb)*(
     .   V2L
     .  -Lsm1DS(s(j4,j1),s(j1,j2),s412)
     .  -Lsm1DS(s(j2,j3),s(j3,j4),s234))
      
      return
      end
      
      double complex function A1phiAgQgmppmR(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1 reflection relation (4.2.3)
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A1phiAgQgmmppL
 
      A1phiAgQgmppmR=A1phiAgQgmmppL(j1,j4,j3,j2,za,zb)
      
      return
      end
      
      double complex function A1phiAgQgmppmf(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.28
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double complex A0phiAgQgmppm,lnrat,l123
      
      l123=lnrat(musq,-(s(j1,j2)+s(j1,j3)+s(j2,j3)))
 
      A1phiAgQgmppmf=A0phiAgQgmppm(j1,j2,j3,j4,za,zb)*(
     . -2d0/3d0*(epinv+l123)-10d0/9d0)
      
      return
      end
      
