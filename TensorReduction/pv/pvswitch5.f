      subroutine pvswitch5(q1,F0,F1,F2,F3,F4,F5,F5a)
      implicit none
      include 'TRydef.f'
      integer n1,n2,n3,n4,n5,ep
      double complex F0(-2:0),F1(y1max,-2:0),F2(y2max,-2:0),
     . F3(y3max,-2:0),F4(y4max,-2:0),F5(y5max,-2:0),
     . F5a(y5max,-2:0)
      double precision q1(4)
      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4

      F5a(y5(n1,n2,n3,n4,n5),ep)=
     . +F5(y5(n1,n2,n3,n4,n5),ep)
     . -q1(n1)*F4(y4(n2,n3,n4,n5),ep)
     . -q1(n2)*F4(y4(n1,n3,n4,n5),ep)
     . -q1(n3)*F4(y4(n1,n2,n4,n5),ep)
     . -q1(n4)*F4(y4(n1,n2,n3,n5),ep)
     . -q1(n5)*F4(y4(n1,n2,n3,n4),ep)

     . +q1(n1)*q1(n2)*F3(y3(n3,n4,n5),ep)
     . +q1(n1)*q1(n3)*F3(y3(n2,n4,n5),ep)
     . +q1(n1)*q1(n4)*F3(y3(n2,n3,n5),ep)
     . +q1(n1)*q1(n5)*F3(y3(n2,n3,n4),ep)
     . +q1(n2)*q1(n3)*F3(y3(n1,n4,n5),ep)
     . +q1(n2)*q1(n4)*F3(y3(n1,n3,n5),ep)
     . +q1(n2)*q1(n5)*F3(y3(n1,n3,n4),ep)
     . +q1(n3)*q1(n4)*F3(y3(n1,n2,n5),ep)
     . +q1(n3)*q1(n5)*F3(y3(n1,n2,n4),ep)
     . +q1(n4)*q1(n5)*F3(y3(n1,n2,n3),ep)

     . -q1(n1)*q1(n2)*q1(n3)*F2(y2(n4,n5),ep)
     . -q1(n1)*q1(n2)*q1(n4)*F2(y2(n3,n5),ep)
     . -q1(n1)*q1(n2)*q1(n5)*F2(y2(n3,n4),ep)
     . -q1(n1)*q1(n3)*q1(n4)*F2(y2(n2,n5),ep)
     . -q1(n1)*q1(n3)*q1(n5)*F2(y2(n2,n4),ep)
     . -q1(n1)*q1(n4)*q1(n5)*F2(y2(n2,n3),ep)
     . -q1(n2)*q1(n3)*q1(n4)*F2(y2(n1,n5),ep)
     . -q1(n2)*q1(n3)*q1(n5)*F2(y2(n1,n4),ep)
     . -q1(n2)*q1(n4)*q1(n5)*F2(y2(n1,n3),ep)
     . -q1(n3)*q1(n4)*q1(n5)*F2(y2(n1,n2),ep)

     . +q1(n1)*q1(n2)*q1(n3)*q1(n4)*F1(n5,ep)
     . +q1(n1)*q1(n2)*q1(n3)*q1(n5)*F1(n4,ep)
     . +q1(n1)*q1(n2)*q1(n4)*q1(n5)*F1(n3,ep)
     . +q1(n1)*q1(n3)*q1(n4)*q1(n5)*F1(n2,ep)
     . +q1(n2)*q1(n3)*q1(n4)*q1(n5)*F1(n1,ep)
 
     . -q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)*F0(ep)


      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end

