      subroutine pvswitch2(q1,F0,F1,F2,F2a)
C----from F0,F1,F2 calculate F2a which has all denominators
C----shifted by q1
      implicit none
      include 'TRydef.f'
      integer n1,n2,ep
      double complex F0(-2:0),F1(y1max,-2:0),F2(y2max,-2:0),
     . F2a(y2max,-2:0)
      double precision q1(4)
      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      F2a(y2(n1,n2),ep)=F2(y2(n1,n2),ep)
     . -q1(n1)*F1(n2,ep)
     . -q1(n2)*F1(n1,ep)
     . +q1(n1)*q1(n2)*F0(ep)
      enddo
      enddo
      enddo
      return
      end
