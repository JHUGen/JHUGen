      integer z1(3),z2(3,3),z3(3,3,3),z4(3,3,3,3),z5(3,3,3,3,3),
     . z6(3,3,3,3,3,3),z7(3,3,3,3,3,3,3)
      double precision delta(3,3)
      common/Darraydef/delta,z1,z2,z3,z4,z5,z6,z7
      integer z1max,z2max,z3max,z4max,z5max,z6max,z7max
      parameter(z1max=3,z2max=6,z3max=10,z4max=15,z5max=21,z6max=28,
     .          z7max=36)
!$omp threadprivate(/Darraydef/)
