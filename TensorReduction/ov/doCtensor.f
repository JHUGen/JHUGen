c--- This is just a wrapping routine that calls ovCtensor or pvCtensor
      subroutine doCtensor(q1,q2,m1s,m2s,m3s,
     . FC0,FC1,FC2,FC3,FC4,FC5,FC6)
      implicit none
C     q1,q2 are the momenta in the propagators
C     m1s,m2s,m3s are the squares of the masses in the propagators
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRtensorcontrol.f'
      double complex FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     . FC3(y3max,-2:0),FC4(y4max,-2:0),FC5(y5max,-2:0),FC6(y6max,-2:0),
     . C00(-2:0),tau3(4,-2:0)
      double precision q1(4),q2(4),p2(4),m1s,m2s,m3s

      if (doovred) then
        p2(:)=q2(:)-q1(:)
        call ovCtensor(q1,p2,m1s,m2s,m3s,FC0,FC1,FC2,FC3,C00,tau3)
        FC4=czip
        FC5=czip
        FC6=czip
      endif
      
      if (dopvred) then
        call pvCtensor(q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6)
      endif
      
      return
      end

