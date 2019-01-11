      subroutine pvswitch1(q1,F0,F1,F1a)
C----from F0,F1 calculate F1a which has all denominators
C----shifted by q1
      implicit none
      include 'TRydef.f'
      integer n1,ep
      double complex F0(-2:0),F1(y1max,-2:0),F1a(y1max,-2:0)
      double precision q1(4)
      do ep=-2,0
      do n1=1,4
      F1a(n1,ep)=F1(n1,ep)-q1(n1)*F0(ep)
      enddo
      enddo
      return
      end
