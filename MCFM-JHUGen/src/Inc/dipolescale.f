c--- contains the dynamic scale for each dipole
      real(dp):: dipscale(0:maxd)
      common/dipolescale/dipscale
!$omp threadprivate(/dipolescale/)
