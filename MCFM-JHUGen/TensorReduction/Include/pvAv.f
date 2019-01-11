      double complex Av(Naa*Namax,-2:0)
      common/AAv/Av
!$omp threadprivate(/AAv/)
