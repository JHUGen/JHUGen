      subroutine ovswitch(p1,F0,F1,F2,F3,F4,F1a,F2a,F3a,F4a)
C----from F0,F1,F2,F3,F4 calculate F1a,F2a,F3a,F4a 
C----which has first denominator
C----shifted by p1
      implicit none
      include 'TRydef.f'
      include 'TRmaxindex.f'
      integer n1,n2,n3,n4,ep
      double complex F0(-2:0),F1(y1max,-2:0),F2(y2max,-2:0),
     & F3(y3max,-2:0),F4(y4max,-2:0)
      double complex F1a(y1max,-2:0),F2a(y2max,-2:0),
     & F3a(y3max,-2:0),F4a(y4max,-2:0)
      double precision p1(4)

      do ep=-2,0
      do n1=1,4
      if (maxeindex .lt. 2) cycle
      F1a(n1,ep)=F1(n1,ep)-p1(n1)*F0(ep)
      do n2=n1,4
      if (maxeindex .lt. 3) cycle
      F2a(y2(n1,n2),ep)=F2(y2(n1,n2),ep)
     & -p1(n1)*F1(n2,ep)-p1(n2)*F1(n1,ep)
     & +p1(n1)*p1(n2)*F0(ep)
      do n3=n2,4
      if (maxeindex .lt. 4) cycle
      F3a(y3(n1,n2,n3),ep)=F3(y3(n1,n2,n3),ep)
     & -p1(n1)*F2(y2(n2,n3),ep)
     & -p1(n2)*F2(y2(n3,n1),ep)
     & -p1(n3)*F2(y2(n1,n2),ep)

     & +p1(n1)*p1(n2)*F1(n3,ep)
     & +p1(n2)*p1(n3)*F1(n1,ep)
     & +p1(n3)*p1(n1)*F1(n2,ep)
  
     & -p1(n1)*p1(n2)*p1(n3)*F0(ep)

      do n4=n3,4
      if (maxeindex .lt. 5) cycle
      F4a(y4(n1,n2,n3,n4),ep)=F4(y4(n1,n2,n3,n4),ep)
     & -p1(n1)*F3(y3(n2,n3,n4),ep)
     & -p1(n2)*F3(y3(n3,n4,n1),ep)
     & -p1(n3)*F3(y3(n4,n1,n2),ep)
     & -p1(n4)*F3(y3(n1,n2,n3),ep)

     & +p1(n1)*p1(n2)*F2(y2(n3,n4),ep)
     & +p1(n1)*p1(n3)*F2(y2(n2,n4),ep)
     & +p1(n1)*p1(n4)*F2(y2(n2,n3),ep)
     & +p1(n2)*p1(n3)*F2(y2(n1,n4),ep)
     & +p1(n2)*p1(n4)*F2(y2(n1,n3),ep)
     & +p1(n3)*p1(n4)*F2(y2(n1,n2),ep)

     & -p1(n1)*p1(n2)*p1(n3)*F1(n4,ep)
     & -p1(n2)*p1(n3)*p1(n4)*F1(n1,ep)
     & -p1(n3)*p1(n4)*p1(n1)*F1(n2,ep)
     & -p1(n4)*p1(n1)*p1(n2)*F1(n3,ep)

     & +p1(n1)*p1(n2)*p1(n3)*p1(n4)*F0(ep)

      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
