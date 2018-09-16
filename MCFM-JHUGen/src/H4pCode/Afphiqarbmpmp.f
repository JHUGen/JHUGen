      function Afphiqarbmpmp(j1,j2,j3,j4,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Afphiqarbmpmp
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
C---arXIv:09060008v1 Eq.(4.15)
      integer:: j1,j2,j3,j4
      complex(dp):: Afphiqarbmppm
      Afphiqarbmpmp=-Afphiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

