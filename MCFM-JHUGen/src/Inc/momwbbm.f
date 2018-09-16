      real(dp):: mom(mxpart,4),bp,bm
      common/momwbbm/mom,bp,bm
!$omp threadprivate(/momwbbm/)

