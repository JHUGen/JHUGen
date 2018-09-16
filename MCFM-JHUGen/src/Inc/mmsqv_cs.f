      real(dp):: mmsqv_cs(0:2,2,2)
      common/mmsqv_cs/mmsqv_cs
!$omp threadprivate(/mmsqv_cs/)
