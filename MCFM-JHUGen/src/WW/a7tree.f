      function A7trees(j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A7trees

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
c---  DKS Eq. 3.15
      integer:: j1,j2,j3,j4,j5,j6,j7
      complex(dp):: a7treea

      a7trees=a7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)
     &       +a7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)

      return
      end

      function A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A7treea
C---This function taken from Eq. 2.22 of DKS
c---Multiplied by a factor of (-i)
C---in the natural way this function refers to the process
c---u(2)+dbar(1)-->nu(3)+e^+(3)+b(5)+bbar(6)+g(7)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: i1,i2,i3,i4
      complex(dp):: z2
      real(dp):: t134,t256
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      t256=s(j2,j5)+s(j2,j6)+s(j5,j6)
      A7treea=za(j1,j3)
     & *(za(j1,j3)*zb(j3,j4)*zb(j2,j5)*z2(j6,j2,j5,j7)/t256
     & +z2(j6,j1,j3,j4)*z2(j1,j2,j7,j5)/za(j7,j2))
     & /(za(j1,j7)*s(j3,j4)*s(j5,j6)*t134)
      return
      end


c      function A7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)
c      implicit none
c      include 'types.f'
c      complex(dp):: A7treeb
C---This function taken from Eq. 2.23 of DKS
c---Multiplied by a factor of (-i)
C---in the natural way this function refers to the process
c---u(2)+dbar(1)-->nu(3)+e^+(3)+b(5)+bbar(6)+g(7)
c
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'sprods_com.f'
c      include 'zprods_decl.f'
c      integer:: j1,j2,j3,j4,j5,j6,j7
c      integer:: i1,i2,i3,i4
c      complex(dp):: z2
c      real(dp):: t127
c      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
c      t127=s(j1,j2)+s(j1,j7)+s(j2,j7)
c      A7treeb=(-za(j3,j6)*zb(j4,j5)
c     & *(z2(j1,j5,j6,j2)*za(j2,j1)+z2(j1,j5,j6,j7)*za(j7,j1))
c     & +za(j1,j3)*z2(j1,j2,j7,j4)*z2(j6,j3,j4,j5)
c     & -za(j1,j6)*z2(j1,j2,j7,j5)*z2(j3,j5,j6,j4))
c     & /(za(j1,j7)*za(j7,j2)*s(j3,j4)*s(j5,j6)*t127)
c      return
c      end
