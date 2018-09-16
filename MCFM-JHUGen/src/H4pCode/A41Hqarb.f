c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      function A41Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41Hqarbmppm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A41phiqarbmppm

      A41Hqarbmppm=A41phiqarbmppm(j1,j2,j3,j4,za,zb)
     &            +A41phiqarbmppm(j2,j1,j4,j3,zb,za)
      
      return
      end

      function A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41Hqarbmpmp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A41phiqarbmpmp

      A41Hqarbmpmp=A41phiqarbmpmp(j1,j2,j3,j4,za,zb)
     &            +A41phiqarbmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
