c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      double complex function A41Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A41phiqarbmppm

      A41Hqarbmppm=A41phiqarbmppm(j1,j2,j3,j4,za,zb)
     .            +A41phiqarbmppm(j2,j1,j4,j3,zb,za)
      
      return
      end

      double complex function A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A41phiqarbmpmp

      A41Hqarbmpmp=A41phiqarbmpmp(j1,j2,j3,j4,za,zb)
     .            +A41phiqarbmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
