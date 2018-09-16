      function atreepm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: atreepm

c---Atreepm is the amplitude for
c---q-(-p4)+Q-(-p2)+l-(-p5) ---> q+(p1)+Q+(p3)+l+(p6)
c---All outgoing particles are right-handed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c-----Note all tree amplitudes divided by im wrt to BDK paper

      atreepm=
     & -zb(j1,j3)*za(j5,j4)*(za(j2,j1)*zb(j1,j6)+za(j2,j3)*zb(j3,j6))
     & /(s(j2,j3)*s(j5,j6)*t(j1,j2,j3))
     & -za(j2,j4)*zb(j6,j1)*(za(j5,j2)*zb(j2,j3)+za(j5,j4)*zb(j4,j3))
     & /(s(j2,j3)*s(j5,j6)*t(j2,j3,j4))

      return
      end


