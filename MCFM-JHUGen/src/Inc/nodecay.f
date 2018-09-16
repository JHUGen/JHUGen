************************************************************************
* nodecay should be set to TRUE if the Higgs decay is not included     *
************************************************************************
      logical:: nodecay
      common/nodecay/nodecay
!$omp threadprivate(/nodecay/)
