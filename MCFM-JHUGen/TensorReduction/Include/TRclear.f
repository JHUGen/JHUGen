      logical clear(1:5)
      common/clear/clear
!$omp threadprivate(/clear/)
