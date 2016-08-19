      integer z1(2),z2(2,2),z3(2,2,2),z4(2,2,2,2),z5(2,2,2,2,2),
     . z6(2,2,2,2,2,2),z7(2,2,2,2,2,2,2)
      double precision delta(2,2)
      common/Carraydef/z1,z2,z3,z4,z5,z6,z7,delta
      integer z1max,z2max,z3max,z4max,z5max,z6max,z7max
      parameter(z1max=2,z2max=3,z3max=4,z4max=5,z5max=6,z6max=7,
     .          z7max=8)
!$omp threadprivate(/Carraydef/)
