c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      double complex function A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiqarbmppm

      A0Hqarbmppm=A0phiqarbmppm(j1,j2,j3,j4,za,zb)
     .           +A0phiqarbmppm(j2,j1,j4,j3,zb,za)
      
      return
      end

      double complex function A0Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiqarbmpmp

      A0Hqarbmpmp=A0phiqarbmpmp(j1,j2,j3,j4,za,zb)
     .           +A0phiqarbmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
