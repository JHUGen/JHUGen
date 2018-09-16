      function F33m(p1sq,p2sq,p3sq)
      implicit none
      include 'types.f'
      complex(dp):: F33m

c      include 'scale.f'
      real(dp):: p1sq,p2sq,p3sq
c      complex(dp):: qlI3
      complex(dp):: I3m

c--- NOTE: checked on 8/30/09 that qlI3 == -I3m
c---       and F33m is defined to be (-1)*(scalar integral)
c      F33m=-qlI3(p1sq,p2sq,p3sq,0._dp,0._dp,0._dp,musq,0)
      F33m=I3m(p1sq,p2sq,p3sq)

      return
      end

