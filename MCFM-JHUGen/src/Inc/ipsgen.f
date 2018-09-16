      logical:: doipsgen
      integer:: ipsgen,maxipsgen
      common/ipsgen/doipsgen,ipsgen,maxipsgen
!$omp threadprivate(/ipsgen/)
