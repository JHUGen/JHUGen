      function atreesl(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: atreesl

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c-----Note all tree amplitudes divided by im wrt to BDK paper
      atreesl=-((za(j4,j5)*zb(j1,j3)
     & *(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     & (s(j1,j2)*s(j5,j6)*t(j1,j2,j3)))
     & +(za(j2,j4)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*zb(j3,j6))/
     & (s(j1,j2)*s(j5,j6)*t(j1,j2,j4))
      return
      end



