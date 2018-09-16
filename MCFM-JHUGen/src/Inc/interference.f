      logical:: interference,bw34_56
      common/interference/interference,bw34_56
      real(dp):: vsymfact
      common/vsymfact/vsymfact
!$omp threadprivate(/interference/,/vsymfact/)
