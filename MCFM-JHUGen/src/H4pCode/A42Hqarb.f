c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      double complex function A42Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A42phiqarbmppm

      A42Hqarbmppm=A42phiqarbmppm(j1,j2,j3,j4,za,zb)
     .            +A42phiqarbmppm(j2,j1,j4,j3,zb,za)
      
      return
      end

      double complex function A42Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A42phiqarbmpmp

      A42Hqarbmpmp=A42phiqarbmpmp(j1,j2,j3,j4,za,zb)
     .            +A42phiqarbmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
