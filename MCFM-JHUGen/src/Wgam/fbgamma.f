      double complex function fbgamma(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5
      double complex L0,L1,Lsm1
      fbgamma=
     . (za(i1,i2)*za(i3,i5))**2/(za(i3,i4)*za(i1,i5)*za(i2,i5)**3)
     . *Lsm1(-s(i1,i2),-s(i3,i4),-s(i1,i5),-s(i3,i4))

     . -0.5d0*(za(i2,i3)*zb(i2,i5))**2*za(i1,i5)
     . /(za(i3,i4)*za(i2,i5)*s(i1,i5)**2)*L1(-s(i3,i4),-s(i1,i5))

     . +za(i2,i3)*zb(i2,i5)
     . *(za(i1,i2)*za(i3,i5)+za(i1,i3)*za(i2,i5))
     .  /(za(i3,i4)*za(i2,i5)**2*s(i1,i5))*L0(-s(i3,i4),-s(i1,i5))

     . -0.5d0*(za(i1,i2)*za(i3,i5))**2*zb(i2,i5)
     . /(za(i3,i4)*za(i1,i5)*za(i2,i5)**2*s(i1,i2))
     . *(s(i1,i5)/s(i1,i2)*L1(-s(i3,i4),-s(i1,i2))
     . +3d0*L0(-s(i3,i4),-s(i1,i2)))

     . -0.5d0*za(i2,i3)*zb(i2,i5)*za(i3,i5)
     . /(zb(i1,i2)*za(i3,i4)*za(i2,i5)**2)

     . -0.5d0*zb(i1,i4)*zb(i2,i5)*zb(i4,i5)
     . /(zb(i1,i2)*zb(i3,i4)*zb(i1,i5)*za(i2,i5))
      return
      end
