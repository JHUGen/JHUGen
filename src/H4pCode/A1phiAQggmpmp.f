      double complex function A1phiAQggmpmpL(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.24
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer j1,j2,j3,j4
      double complex zab2,V1L,A0phiAQggmpmp,L2,L1,L0,lnrat,sum,
     . A0phidAQggmpmp
      double complex l23,l34,l41,l12,Lsm1_2me,Lsm1DS
      double precision s3,mhsq,s123,s234,s341,s412
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      
      s123=s3(j1,j2,j3)
      s234=s3(j2,j3,j4)
      s341=s3(j3,j4,j1)
      s412=s3(j4,j1,j2)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))
      l34=lnrat(musq,-s(j3,j4))
      l41=lnrat(musq,-s(j4,j1))

      V1L=
     . -epinv**2-epinv*l23-0.5d0*l23**2
     . -epinv**2-epinv*l34-0.5d0*l34**2
     . -epinv**2-epinv*l41-0.5d0*l41**2
     . +13d0/6d0*(epinv+l12)
     . +119d0/18d0-deltar/6d0
     . -Lsm1_2me(s123,s234,s(j2,j3),mhsq) 
     . -Lsm1_2me(s341,s412,s(j4,j1),mhsq) 
     . -Lsm1_2me(s412,s123,s(j1,j2),mhsq) 

      A1phiAQggmpmpL=A0phiAQggmpmp(j1,j2,j3,j4,za,zb)
     . *(V1L-13d0/6d0*lnrat(-s412,-s(j1,j2))
     .      -Lsm1DS(s(j3,j4),s(j4,j1),s341)
     .      -Lsm1DS(s(j1,j2),s(j2,j3),s123))

      sum=
     . +za(j1,j4)**2*za(j2,j3)**3/(za(j1,j2)*za(j3,j4)*za(j2,j4)**3)
     .  *(Lsm1DS(s(j2,j3),s(j3,j4),s234)+Lsm1DS(s(j4,j1),s(j1,j2),s412))
     . +2d0/3d0*za(j1,j2)**2*za(j3,j4)**2*zb(j2,j4)**3/za(j1,j4)
     .  *L2(-s412,-s(j1,j2))/s(j1,j2)**3
     . -0.5d0*za(j1,j2)*za(j2,j3)*za(j3,j4)*zb(j2,j4)**2/za(j2,j4)
     .  *L1(-s234,-s(j3,j4))/s(j3,j4)**2
     . +(0.5d0*za(j1,j4)*zab2(j3,j1,j2,j4)**2/za(j2,j4)
     .  -1d0/3d0*za(j1,j3)*za(j1,j4)*zab2(j3,j1,j2,j4)**2
     .          /(za(j1,j2)*za(j3,j4))
     .  -2d0/3d0*za(j1,j3)*za(j1,j2)*za(j3,j4)*zb(j2,j4)**2/za(j1,j4))
     .  *L1(-s412,-s(j1,j2))/s(j1,j2)**2
     . -(za(j1,j2)*zb(j2,j4)*zab2(j3,j1,j4,j2)**2/(za(j1,j4)*zb(j1,j2))
     .  +0.5d0*za(j1,j4)*za(j2,j3)**2*zb(j2,j4)**2/za(j2,j4))
     .  *L1(-s412,-s(j4,j1))/s(j4,j1)**2
     . -(zb(j2,j4)*za(j3,j4)*zab2(j1,j2,j3,j4)**2/(za(j1,j4)*zb(j3,j4))
     .  +0.5d0*za(j1,j4)*za(j2,j3)**2*zb(j2,j4)**2/za(j2,j4))
     .  *L1(-s234,-s(j2,j3))/s(j2,j3)**2
     . +(3d0*za(j1,j3)**2*zab2(j3,j1,j2,j4)/(za(j1,j2)*za(j3,j4))
     .  +2d0*za(j1,j3)*zab2(j3,j1,j2,j4)**2
     .      /(za(j1,j2)*za(j3,j4)*zb(j1,j4))
     .  +1d0/3d0*za(j1,j3)**2*zb(j2,j4)/za(j1,j4)
     .  -zab2(j3,j1,j2,j4)**2/(zb(j1,j4)*za(j2,j4)))
     .  *L0(-s412,-s(j1,j2))/s(j1,j2)
     . +3d0*za(j2,j3)*za(j1,j3)*zb(j2,j4)/za(j2,j4)
     .  *(L0(-s234,-s(j2,j3))/s(j2,j3)+L0(-s412,-s(j4,j1))/s(j4,j1))
     . +za(j2,j3)*zb(j2,j4)*(za(j1,j2)*za(j3,j4)
     .                      +2d0*za(j1,j4)*za(j2,j3))/za(j2,j4)**2
     .  *L0(-s234,-s(j3,j4))/s(j3,j4)
     . -(1d0/3d0*za(j1,j3)**3/(za(j1,j2)*za(j3,j4)*za(j1,j4))
     .  +0.5d0*za(j2,j3)*za(j1,j3)**2/(za(j1,j2)*za(j2,j4)*za(j3,j4))
     .  +za(j2,j3)**2*zb(j2,j4)/(za(j2,j4)**2*zb(j1,j4))
     .  +2d0*za(j2,j3)**3*za(j1,j4)*zb(j2,j4)
     .      /(za(j2,j4)**2*za(j3,j4)*za(j1,j2)*zb(j1,j4)))
     .  *lnrat(-s412,-s(j1,j2))     
     . +za(j1,j2)**2*za(j3,j4)*zb(j2,j4)
     .  /(za(j2,j4)**2*za(j1,j4)*zb(j3,j4))*lnrat(-s234,-s(j2,j3))
     . +za(j3,j4)**2*za(j1,j2)*zb(j2,j4)
     .  /(za(j2,j4)**2*za(j1,j4)*zb(j1,j2))*lnrat(-s412,-s(j4,j1))
     . -5d0/6d0*za(j1,j3)**2*zb(j2,j4)/(s(j1,j2)*za(j1,j4))
     . -1d0/3d0*za(j1,j3)**2*zab2(j3,j1,j2,j4)
     .  /(s(j1,j2)*za(j1,j2)*za(j3,j4))
     . -1d0/3d0*za(j1,j3)*zb(j2,j4)
     .  *(2d0*za(j3,j4)*zb(j4,j2)+za(j3,j1)*zb(j1,j2))
     .  /(s412*za(j1,j4)*zb(j1,j2))
     . +0.5d0*(
     .  zb(j2,j4)*zab2(j3,j1,j2,j4)*zab2(j3,j2,j4,j1)
     .  /(s412*zb(j1,j4)*zb(j1,j2)*za(j2,j4))
     . -za(j1,j3)**2*zb(j1,j4)/(s(j1,j2)*za(j2,j4))
     . -za(j1,j3)*zb(j2,j4)*za(j3,j4)/(za(j1,j4)*zb(j1,j2)*za(j2,j4))
     . +za(j1,j2)*zb(j2,j4)**2/(zb(j2,j3)*zb(j3,j4)*za(j2,j4)) )
     . +za(j1,j3)*zb(j2,j4)*zab2(j1,j2,j3,j4)
     .  /(s(j2,j3)*za(j1,j4)*zb(j3,j4))
     . -za(j1,j3)*zb(j2,j4)*zab2(j3,j1,j4,j2)
     .  /(s(j4,j1)*za(j1,j4)*zb(j1,j2))
     . -zb(j2,j4)**2*za(j3,j4)*za(j2,j3)
     .  /(s(j4,j1)*za(j2,j4)*zb(j1,j2))
     . -2d0*A0phidAQggmpmp(j1,j2,j3,j4,za,zb)
     
      A1phiAQggmpmpL=A1phiAQggmpmpL+sum

      return
      end

      double complex function A1phiAQggmpmpR(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.25
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      include 'deltar.f'
      integer j1,j2,j3,j4
      double complex zab2,VR,A0phiAQggmpmp,Lsm1_2me,lnrat,l12,L0,L1,
     & Lsm1DS
      double precision s3,mhsq,s341,s234,s412
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      l12=lnrat(musq,-s(j1,j2))
      s341=s3(j3,j4,j1)
      s234=s3(j2,j3,j4)
      s412=s3(j4,j1,j2)

      VR=-epinv**2-epinv*l12-0.5d0*l12**2
     .   -3d0/2d0*(epinv+l12)-7d0/2d0-deltar/2d0
     .   -Lsm1_2me(s234,s341,s(j3,j4),mhsq)
      
      A1phiAQggmpmpR=
     . +za(j1,j2)**2*za(j3,j4)**2/(za(j1,j4)*za(j2,j4)**3)
     . *(Lsm1DS(s(j2,j3),s(j3,j4),s234)+Lsm1DS(s(j4,j1),s(j1,j2),s412))
     . -0.5d0*za(j1,j2)**2*za(j3,j4)**2*zb(j2,j4)**2
     .       /(za(j1,j4)*za(j2,j4))*L1(-s412,-s(j1,j2))/s(j1,j2)**2
     . +zb(j2,j4)*za(j3,j4)*zab2(j1,j2,j3,j4)**2/(za(j1,j4)*zb(j3,j4))
     .       *L1(-s234,-s(j2,j3))/s(j2,j3)**2
     . -(za(j1,j2)*za(j3,j4)*zb(j2,j4)**3/zb(j2,j3)
     .  +0.5d0*za(j2,j3)**3*zab2(j1,j3,j4,j2)**2
     .       /(za(j1,j2)*za(j3,j4)*za(j2,j4)))
     .  *L1(-s234,-s(j3,j4))/s(j3,j4)**2
     . -0.5d0*za(j1,j4)*za(j2,j3)**2*zb(j2,j4)**2/za(j2,j4)
     .       *(L1(-s412,-s(j4,j1))/s(j4,j1)**2
     .        -L1(-s234,-s(j2,j3))/s(j2,j3)**2)
     . -za(j1,j2)**2*za(j3,j4)**2*zb(j2,j4)/(za(j1,j4)*za(j2,j4)**2)
     .       *L0(-s412,-s(j1,j2))/s(j1,j2)
     . +za(j2,j3)*zb(j2,j4)*(
     .  2d0*za(j1,j2)*za(j3,j4)+za(j1,j4)*za(j2,j3))/(za(j2,j4)**2)
     .       *L0(-s412,-s(j4,j1))/s(j4,j1)
     . -za(j1,j2)**2*za(j3,j4)*zb(j2,j4)
     .  /(za(j2,j4)**2*za(j1,j4)*zb(j3,j4))*lnrat(-s234,-s(j2,j3))
     . +(za(j3,j4)*za(j1,j2)*zb(j2,j4)/(za(j2,j4)**2*zb(j2,j3))
     .  +0.5d0*za(j2,j3)*za(j1,j3)**2/(za(j1,j2)*za(j2,j4)*za(j3,j4)))
     .       *lnrat(-s234,-s(j3,j4))
     . -0.5d0*zb(j2,j4)*zab2(j3,j1,j2,j4)*zab2(j3,j1,j4,j2)
     .  /(s(j4,j1)*s412*zb(j1,j2))
     . -0.5d0*zb(j2,j4)**2*za(j3,j4)*za(j2,j3)
     .  /(s(j4,j1)*za(j2,j4)*zb(j1,j2))
     . -0.5d0*(za(j1,j3)*za(j2,j3)**2*zab2(j1,j3,j4,j2))
     .  /(s(j3,j4)*za(j3,j4)*za(j1,j2)*za(j2,j4))
     . +0.5d0*za(j2,j3)*zb(j2,j4)*zab2(j1,j3,j4,j2)*(s(j2,j3)+s(j3,j4))
     .  /(s(j3,j4)*s234*zb(j2,j3)*za(j2,j4))
     . +0.5d0*zb(j2,j4)**2*zab2(j1,j2,j3,j4)/(s234*zb(j2,j3)*zb(j3,j4))
     . -za(j1,j2)*zb(j2,j4)*za(j3,j4)*zab2(j1,j2,j3,j4)
     .  /(s(j2,j3)*za(j1,j4)*zb(j3,j4)*za(j2,j4))
     . +za(j1,j3)*zb(j2,j4)/(zb(j2,j3)*za(j2,j4))
     
      A1phiAQggmpmpR=A1phiAQggmpmpR
     .             +A0phiAQggmpmp(j1,j2,j3,j4,za,zb)*VR
      
      return
      end

      double complex function A1phiAQggmpmpF(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.26
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'epinv.f'
      integer j1,j2,j3,j4
      double complex l12,zab2,L2,lnrat,A0phiAQggmpmp
      double precision s412
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      s412=s(j4,j1)+s(j4,j2)+s(j1,j2)

      l12=lnrat(musq,-s(j1,j2))

      A1phiAQggmpmpF=A0phiAQggmpmp(j1,j2,j3,j4,za,zb)
     . *(-2d0/3d0*(epinv+l12)-10d0/9d0)
     .-1d0/3d0*(za(j1,j4)**2*zab2(j3,j1,j2,j4)**3/(za(j1,j2)*za(j3,j4))
     .         +za(j1,j2)**2*zb(j2,j4)**3*za(j3,j4)**2/za(j1,j4))
     .        *L2(-s412,-s(j1,j2))/s(j1,j2)**3
     .+1d0/3d0*A0phiAQggmpmp(j1,j2,j3,j4,za,zb)*lnrat(-s412,-s(j1,j2))
     .+0.5d0*za(j1,j3)**2*zb(j2,j4)/(s(j1,j2)*za(j1,j4))
     .+1d0/6d0*za(j1,j3)*(za(j3,j1)*zb(j1,j2)*zab2(j3,j1,j4,j2)
     .                   -(za(j3,j4)*zb(j4,j2))**2)
     .                  /(zb(j1,j2)*za(j3,j4)*za(j1,j4)*s412)
     .+1d0/6d0*za(j1,j3)**3*s412
     .        /(s(j1,j2)*za(j1,j2)*za(j3,j4)*za(j1,j4))
      
      return
      end
