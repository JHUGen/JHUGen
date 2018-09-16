      real(dp):: xdelg1_z,xdelg1_g,xlambda_g,xlambda_z,xdelk_g,
     & xdelk_z
      common/xanomcoup/xdelg1_z,xdelg1_g,xlambda_g,xlambda_z,xdelk_g,
     & xdelk_z
!$omp threadprivate(/xanomcoup/)
