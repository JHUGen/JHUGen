      real(dp):: msq_cs(0:2,-nf:nf,-nf:nf)
      common/msq_cs/msq_cs
!$omp threadprivate(/msq_cs/)
