      real(dp):: msqv_cs(0:2,-nf:nf,-nf:nf)
      common/msqv_cs/msqv_cs
!$omp threadprivate(/msqv_cs/)
