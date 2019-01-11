      integer y1(4),y2(4,4),y3(4,4,4),y4(4,4,4,4),
     . y5(4,4,4,4,4),y6(4,4,4,4,4,4),y7(4,4,4,4,4,4,4)
      integer y1max,y2max,y3max,y4max,y5max,y6max,y7max
      parameter(y1max=4,y2max=10,y3max=20,y4max=35,y5max=56,y6max=84,
     .          y7max=120)
      common/ydef/y1,y2,y3,y4,y5,y6,y7
!$omp threadprivate(/ydef/)
