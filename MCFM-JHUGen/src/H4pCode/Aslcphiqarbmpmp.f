      function Aslcphiqarbmpmp(j1,j2,j3,j4,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Aslcphiqarbmpmp
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
C---arXIv:09060008v1 Eq.(4.14)
      complex(dp):: Aslcphiqarbmppm
      integer:: j1,j2,j3,j4
      Aslcphiqarbmpmp=-Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

