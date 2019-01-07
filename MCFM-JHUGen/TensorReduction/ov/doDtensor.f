c--- This is just a wrapping routine that calls ovDtensor or pvDtensor
      subroutine doDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     . FD0,FD1,FD2,FD3,FD4,FD5,FD6)
      implicit none
C     q1,q2,q3 are the loop offset momenta
C     m1s,m2s,m3s,m4s are the squares of the masses in the propagators
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRtensorcontrol.f'
      double complex FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     . FD3(y3max,-2:0),FD4(y4max,-2:0),FD5(y5max,-2:0),FD6(y6max,-2:0)
      double precision p2(4),p3(4),q1(4),q2(4),q3(4),m1s,m2s,m3s,m4s
      
      if (doovred) then
        p3(:)=q3(:)-q2(:)
        p2(:)=q2(:)-q1(:)
        call ovDtensor(q1,p2,p3,m1s,m2s,m3s,m4s,FD0,FD1,FD2,FD3,FD4)
        FD5=czip
        FD6=czip
      endif
      
      if (dopvred) then
        call pvDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     &   FD0,FD1,FD2,FD3,FD4,FD5,FD6)
      endif
      
      return
      end

