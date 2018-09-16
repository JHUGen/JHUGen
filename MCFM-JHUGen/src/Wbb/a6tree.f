      double complex function a6treea(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c---  DKS Eq. 2.8 (multiplied by a factor of (-i))
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      integer i1,i2,i3,i4
      double complex z2
      double precision t
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      a6treea=za(j1,j3)*zb(j2,j5)*z2(j6,j2,j5,j4)
     . /(s(j3,j4)*s(j5,j6)*t(j1,j3,j4))

      return 
      end

      double complex function a6treeb(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c---  DKS Eq. 2.9 (multiplied by a factor of (-i))
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      integer i1,i2,i3,i4
      double complex z2
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      a6treeb=(za(j1,j3)*zb(j2,j5)*z2(j6,j2,j5,j4)
     .        +zb(j2,j4)*za(j1,j6)*z2(j3,j1,j6,j5)
     . )/(s(j1,j2)*s(j3,j4)*s(j5,j6))

      return 
      end

      double complex function a6trees(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      double complex a6treea
      a6trees=a6treea(j1,j2,j3,j4,j5,j6,za,zb)
     .       +a6treea(j1,j2,j6,j5,j4,j3,za,zb)
      return 
      end
