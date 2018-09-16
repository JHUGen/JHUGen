      function dot(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: dot
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i,j
      real(dp):: p(mxpart,4)
      dot=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      return
      end

