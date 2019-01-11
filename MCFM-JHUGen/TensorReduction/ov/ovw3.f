      double precision function ovw3(n1,n2,p1,p2,p3,del3in)
      implicit none
      include 'ovupdown.f'
      double precision p1(4),p2(4),p3(4),n1up(4),n2up(4),del3in,del4
      integer n1,n2
      n1up(1:4)=up(1:4,n1)
      n2up(1:4)=up(1:4,n2)
      ovw3=del4(n1up,p1,p2,p3,n2up,p1,p2,p3)/del3in
      return
      end
