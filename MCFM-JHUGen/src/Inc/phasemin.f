      real(dp):: taumin
      common/taumin/taumin
      include 'xmin.f'
!$omp threadprivate(/taumin/)

