      function A0Hggggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hggggmpmp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiggggmpmp
      A0Hggggmpmp=A0phiggggmpmp(j1,j2,j3,j4,za,zb)
     &           +A0phiggggmpmp(j4,j1,j2,j3,zb,za)
      return
      end
