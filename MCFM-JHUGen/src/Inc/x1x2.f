! include 'x1x2.f' with proper openmp pragma
      double precision xx(2)
      common/x1x2/xx
!$omp threadprivate(/x1x2/)







