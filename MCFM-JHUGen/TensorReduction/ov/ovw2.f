      double precision function ovw2(n1,n2,p1,p2,del2in)
      implicit none
      include 'ovupdown.f'
      double precision p1(4),p2(4),n1up(4),n2up(4),del2in,del3
      integer n1,n2
      n1up(1:4)=up(1:4,n1)
      n2up(1:4)=up(1:4,n2)
      ovw2=del3(n1up,p1,p2,n2up,p1,p2)/del2in
      return
      end
