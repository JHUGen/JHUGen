      double complex Epv(Nee*Nemax,-2:0)
      common/Epv/Epv
!$omp threadprivate(/Epv/)
