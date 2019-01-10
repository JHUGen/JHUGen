c--- This is just a wrapping routine that calls ovBtensor or pvBtensor
      subroutine doBtensor(q1,m1s,m2s,FB0,FB1,FB2,FB3,FB4,FB5,FB6)
      implicit none
C     q1 is the momentum in the loop = p1 the external momenta
C     m1s,m2s are the squares of the internal masses
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRtensorcontrol.f'
      double complex FB0(-2:0),FB1(y1max,-2:0),
     . FB2(y2max,-2:0),FB3(y3max,-2:0),FB4(y4max,-2:0),FB5(y5max,-2:0)
     . ,FB6(y6max,-2:0),B00(-2:0)
      double precision q1(4),m1s,m2s
      
      if (doovred) then
        call ovBtensor(q1,m1s,m2s,FB0,FB1,FB2,B00)      
        FB3=czip
        FB4=czip
        FB5=czip
        FB6=czip
      endif
      
      if (dopvred) then
        call pvBtensor(q1,m1s,m2s,FB0,FB1,FB2,FB3,FB4,FB5,FB6)
      endif
      
      return
      end

