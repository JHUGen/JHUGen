c--- NOTE: this version avoids conflicts with Looptools
      real(dp):: b0
      common/QCDb0/b0
!$omp threadprivate(/QCDb0/)
