      real(dp):: ptjetmin,etajetmin,etajetmax,ptbjetmin,etabjetmax
      common/jetcuts/ptjetmin,etajetmin,etajetmax,ptbjetmin,etabjetmax
!$omp threadprivate(/jetcuts/)
