c--- This is just a wrapping routine that calls ovDtensor or pvDtensor
      subroutine doEtensor(q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     & FE0,FE1,FE2,FE3,FE4,FE5,FE6)
      implicit none
C     q1,q2,q3 are the loop offset momenta
C     m1s,m2s,m3s,m4s are the squares of the masses in the propagators
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRtensorcontrol.f'
      double complex FE0(-2:0),FE1(y1max,-2:0),FE2(y2max,-2:0),
     & FE3(y3max,-2:0),FE4(y4max,-2:0),FE5(y5max,-2:0),FE6(y6max,-2:0)
      double precision p2(4),p3(4),p4(4),q1(4),q2(4),q3(4),q4(4),
     & m1s,m2s,m3s,m4s,m5s
      
      if (doovred) then
        p4(:)=q4(:)-q3(:)
        p3(:)=q3(:)-q2(:)
        p2(:)=q2(:)-q1(:)
        call ovEtensor(q1,p2,p3,p4,m1s,m2s,m3s,m4s,m5s,
     &                 FE0,FE1,FE2,FE3,FE4,FE5)
      endif
      
      if (dopvred) then
        call pvEtensor(q1,p2,p3,p4,m1s,m2s,m3s,m4s,m5s,
     &                 FE0,FE1,FE2,FE3,FE4,FE5)
      endif
      
      FE6=czip

      return
      end

