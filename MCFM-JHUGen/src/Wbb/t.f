      function t(j1,j2,j3)
      implicit none
      include 'types.f'
      real(dp):: t

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      integer:: j1,j2,j3
      t=s(j1,j2)+s(j2,j3)+s(j1,j3)
      end
