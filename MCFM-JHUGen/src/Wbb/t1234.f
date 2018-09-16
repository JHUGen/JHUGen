      function t2(j1,j2,j3,j4)
      implicit none
      include 'types.f'
      complex(dp):: t2

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: j1,j2,j3,j4
      t2=zb(j1,j2)*za(j4,j2)+zb(j1,j3)*za(j4,j3)
      return
      end

