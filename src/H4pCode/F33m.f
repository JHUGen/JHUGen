      double complex function F33m(p1sq,p2sq,p3sq)
      implicit none
c      include 'scale.f'
      double precision p1sq,p2sq,p3sq
c      double complex qlI3
      double complex I3m

c--- NOTE: checked on 8/30/09 that qlI3 == -I3m
c---       and F33m is defined to be (-1)*(scalar integral)
c      F33m=-qlI3(p1sq,p2sq,p3sq,0d0,0d0,0d0,musq,0)
      F33m=I3m(p1sq,p2sq,p3sq)

      return
      end

