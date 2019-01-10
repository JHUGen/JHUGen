      double precision function dot(p,i,j)
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4)
      dot=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      return
      end

