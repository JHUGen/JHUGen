      logical doovred,dopvred
      integer TRtensorcontrol
      common/TRtensorcontrol/doovred,dopvred,TRtensorcontrol
!$omp threadprivate(/TRtensorcontrol/)
