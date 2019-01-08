      double complex function Aslcphiqarbmpmp(j1,j2,j3,j4,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
C---arXIv:09060008v1 Eq.(4.14)
      double complex Aslcphiqarbmppm
      integer j1,j2,j3,j4
      Aslcphiqarbmpmp=-Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

