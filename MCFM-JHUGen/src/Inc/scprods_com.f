      real(dp):: sc(mxpart,mxpart)
      common/scprods/sc
!$omp threadprivate(/scprods/)
