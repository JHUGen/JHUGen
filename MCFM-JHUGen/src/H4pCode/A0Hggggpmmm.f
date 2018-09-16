      function A0Hggggpmmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hggggpmmm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiggggpmmm
      A0Hggggpmmm=A0phiggggpmmm(j1,j2,j3,j4,za,zb)
c     &           +A0phiggggmppp(j1,j2,j3,j4,zb,za) ! This term is zero
      return
      end
