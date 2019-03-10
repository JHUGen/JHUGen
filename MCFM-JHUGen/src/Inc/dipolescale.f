c--- contains the dynamic scale for each dipole
      double precision dipscale(0:maxd)
      common/dipolescale/dipscale
!$omp threadprivate(/dipolescale/)
