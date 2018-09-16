c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      function A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hqarbmppm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiqarbmppm

      A0Hqarbmppm=A0phiqarbmppm(j1,j2,j3,j4,za,zb)
     &           +A0phiqarbmppm(j2,j1,j4,j3,zb,za)
      
      return
      end

      function A0Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0Hqarbmpmp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiqarbmpmp

      A0Hqarbmpmp=A0phiqarbmpmp(j1,j2,j3,j4,za,zb)
     &           +A0phiqarbmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
