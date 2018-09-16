      double complex function A1phiAQggmppmL(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.19
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer j1,j2,j3,j4
      double complex zab2,V1L,A0phiAQggmppm,L2,L1,L0,lnrat,sum,
     . A0phidAQggmppm
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
      A1phiAQggmppmL=A0phiAQggmppm(j1,j2,j3,j4,za,zb)
     . *(V1L-Lsm1DS(s(j2,j3),s(j3,j4),s234)
     .      -Lsm1DS(s(j4,j1),s(j1,j2),s3(j4,j1,j2)))
      sum=
     . +za(j1,j4)**3/(za(j1,j2)*za(j3,j4)*za(j1,j3))
     . *(Lsm1DS(s(j1,j2),s(j2,j3),s123)
     .  +Lsm1DS(s(j3,j4),s(j4,j1),s3(j3,j4,j1)))
     . +(4d0/3d0*za(j1,j3)**2*zab2(j4,j1,j2,j3)**3/(za(j1,j2)*za(j3,j4))
     . -za(j1,j2)*zb(j2,j3)**2*za(j3,j4)*zab2(j4,j1,j2,j3)
     . -1d0/3d0*za(j1,j2)**2*zb(j2,j3)**3*za(j3,j4)**2/za(j1,j3))
     . *L2(-s123,-s(j1,j2))/s(j1,j2)**3
     . +(0.5d0*za(j1,j3)**2*za(j2,j4)*zab2(j4,j1,j2,j3)**2
     . /(za(j1,j2)*za(j2,j3)*za(j3,j4))
     . +za(j1,j3)*za(j1,j4)*zab2(j4,j1,j2,j3)**2/(za(j1,j2)*za(j3,j4))
     . +0.5d0*za(j1,j2)*za(j3,j4)*zb(j2,j3)**2*za(j1,j4)/za(j1,j3))
     . *L1(-s123,-s(j1,j2))/s(j1,j2)**2
     . -0.5d0*za(j1,j2)*za(j3,j4)*za(j2,j4)*zb(j2,j3)**2
     . /za(j2,j3)*L1(-s234,-s(j3,j4))/s(j3,j4)**2
     . +za(j1,j4)**2*zab2(j4,j1,j2,j3)/(za(j1,j2)*za(j3,j4))
     . *L0(-s123,-s(j1,j2))/s(j1,j2)

     . -2d0*za(j1,j4)*za(j2,j4)*zb(j2,j3)/za(j2,j3)
     . *(L0(-s123,-s(j1,j2))/s(j1,j2)
     .  +L0(-s234,-s(j3,j4))/s(j3,j4))

c--- check: looks like there should be a factor of
c---        (-im) in front of A0 here
     . -5d0/6d0*(2d0*A0phiAQggmppm(j1,j2,j3,j4,za,zb)
     . +za(j1,j4)**3/(za(j1,j2)*za(j3,j4)*za(j1,j3)))
     . *lnrat(-s123,-s(j1,j2))

     . +5d0/6d0*za(j1,j4)**2*zab2(j4,j1,j2,j3)
     . /(s(j1,j2)*za(j1,j2)*za(j3,j4))
     . -1d0/6d0*za(j1,j4)**2*zb(j2,j3)*za(j3,j4)
     . /(za(j2,j3)*za(j1,j3)*zab2(j4,j1,j3,j2))
     . +2d0/3d0*za(j1,j4)*zb(j2,j3)*za(j3,j4)
     . /(zb(j1,j2)*za(j2,j3)*za(j1,j3))

     . -2d0/3d0*za(j1,j4)*za(j2,j4)*zab2(j4,j1,j3,j2)
     . /(s(j1,j2)*za(j2,j3)*za(j3,j4))
     . +1d0/3d0*zb(j2,j3)*zab2(j4,j1,j3,j2)*zab2(j4,j2,j3,j1)
     . /(s123*zb(j1,j2)**2*za(j2,j3))

     . -1d0/6d0*zab2(j4,j2,j3,j1)
     . *(za(j4,j1)*zb(j1,j2)+2d0*za(j4,j3)*zb(j3,j2))
     . *(2d0*za(j4,j1)*zb(j1,j2)+za(j4,j3)*zb(j3,j2))**2
     . /(s123*zb(j1,j2)**2*za(j2,j3)*za(j3,j4)*zab2(j4,j1,j3,j2))

     . +0.5d0*zb(j2,j3)*zab2(j2,j1,j4,j3)
     . /(zb(j1,j4)*za(j2,j3)*zb(j3,j4))

     . +0.5d0*zab2(j4,j1,j3,j2)*zab2(j4,j1,j2,j3)
     . /(s123*za(j2,j3)*zb(j1,j2))

     . -0.5d0*za(j1,j4)**2*zb(j1,j3)/(s(j1,j2)*za(j2,j3))

     . -2d0*A0phidAQggmppm(j1,j2,j3,j4,za,zb)
     
      A1phiAQggmppmL=A1phiAQggmppmL+sum

      return
      end

      double complex function A1phiAQggmppmR(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.21
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      include 'deltar.f'
      integer j1,j2,j3,j4
      double complex zab2,VR,A0phiAQggmppm,Lsm1_2me,lnrat,l12,L0,L1,
     . Lsm1DS
      double precision s3,mhsq,s123,s234,s341
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      l12=lnrat(musq,-s(j1,j2))
      s123=s3(j1,j2,j3)
      s234=s3(j2,j3,j4)
      s341=s3(j3,j4,j1)

      VR=-epinv**2-epinv*l12-0.5d0*l12**2
     .   -3d0/2d0*(epinv+l12)-7d0/2d0-deltar/2d0
     .   -Lsm1_2me(s3(j2,j3,j4),s3(j3,j4,j1),s(j3,j4),mhsq)
      
      A1phiAQggmppmR=
     . za(j1,j4)**2/(za(j2,j3)*za(j1,j3))*
     .  (Lsm1DS(s(j1,j2),s(j2,j3),s123)+Lsm1DS(s(j3,j4),s(j4,j1),s341))
     . -0.5d0*za(j1,j2)**2*zb(j2,j3)**2*za(j3,j4)**2
     .  /(za(j2,j3)*za(j1,j3))*L1(-s123,-s(j1,j2))/(s(j1,j2)**2)
     . +0.5d0*za(j2,j4)**3*zab2(j1,j3,j4,j2)**2
     .  /(za(j1,j2)*za(j2,j3)*za(j3,j4))
     .  *L1(-s234,-s(j3,j4))/(s(j3,j4)**2)
     . -2d0*za(j1,j2)*za(j3,j4)*za(j1,j4)*zb(j2,j3)
     .  /(za(j2,j3)*za(j1,j3))*L0(-s123,-s(j1,j2))/s(j1,j2)
     . -2d0*za(j1,j4)*za(j2,j4)*zb(j2,j3)/za(j2,j3)
     .  *L0(-s234,-s(j3,j4))/s(j3,j4)
     . -3d0/2d0*za(j1,j4)**2/(za(j2,j3)*za(j1,j3))
     .  *lnrat(-s123,-s(j1,j2))
     . +1d0/2d0*A0phiAQggmppm(j1,j2,j3,j4,za,zb)*lnrat(-s234,-s(j3,j4))
     . +0.5d0*(
     .  za(j1,j4)*zb(j2,j3)*za(j3,j4)/(zb(j1,j2)*za(j2,j3)*za(j1,j3))
     . +zb(j2,j3)*zb(j1,j3)*zab2(j2,j1,j4,j3)
     .  /(zb(j3,j4)*zb(j1,j4)*zab2(j2,j3,j4,j1)) 
     . -zab2(j4,j1,j3,j2)*zab2(j4,j1,j2,j3)
     .  /(s123*za(j2,j3)*zb(j1,j2))
     . +za(j1,j4)**2*za(j2,j4)**2*(s(j2,j1)+s(j2,j3)+s(j2,j4))
     .  /(za(j1,j2)*za(j2,j3)*za(j3,j4)**2*zab2(j2,j1,j4,j3))
     . -s341**2*zb(j2,j3)*za(j2,j4)**3
     .  /(s(j3,j4)*za(j2,j3)*za(j3,j4)
     .    *zab2(j2,j1,j4,j3)*zab2(j2,j3,j4,j1))
     .        )
      
      A1phiAQggmppmR=A1phiAQggmppmR
     .             +A0phiAQggmppm(j1,j2,j3,j4,za,zb)*VR
      
      return
      end

      double complex function A1phiAQggmppmF(j1,j2,j3,j4,za,zb)
      implicit none
C     implementation of arXiv:0906.0008v1, Eq. 4.23
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'epinv.f'
      integer j1,j2,j3,j4
      double complex l12,zab2,L2,lnrat,A0phiAQggmppm
      double precision s3
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)

      l12=lnrat(musq,-s(j1,j2))

      A1phiAQggmppmF=A0phiAQggmppm(j1,j2,j3,j4,za,zb)
     . *(-2d0/3d0*(epinv+l12)-10d0/9d0)
     .+1d0/3d0*(za(j1,j2)**2*zb(j2,j3)**3*za(j3,j4)**2/za(j1,j3)
     .         -zab2(j4,j1,j2,j3)**3*za(j1,j3)**2/(za(j1,j2)*za(j3,j4)))
     .        *L2(-s3(j1,j2,j3),-s(j1,j2))/s(j1,j2)**3
     .-1d0/3d0*(za(j1,j4)**2/(za(j2,j3)*za(j1,j3))
     .         +za(j1,j4)**2*za(j2,j4)/(za(j1,j2)*za(j2,j3)*za(j3,j4)))
     .        *lnrat(-s3(j1,j2,j3),-s(j1,j2))
     .-0.5d0*za(j1,j4)**2*zb(j2,j3)/(s(j1,j2)*za(j1,j3))
     .+1d0/6d0*za(j1,j4)*(za(j4,j1)*zb(j1,j2)*zab2(j4,j1,j3,j2)
     .                   -(za(j4,j3)*zb(j3,j2))**2)
     .                  /(zb(j1,j2)*za(j3,j4)*za(j1,j3)*s3(j1,j2,j3))
     .+1d0/6d0*za(j1,j4)**3*s3(j1,j2,j3)
     .        /(s(j1,j2)*za(j1,j2)*za(j3,j4)*za(j1,j3))
      
      return
      end
