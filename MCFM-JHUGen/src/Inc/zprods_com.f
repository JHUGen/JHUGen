      double complex za(mxpart,mxpart),zb(mxpart,mxpart)
      common/zprods/za,zb
!$omp threadprivate(/zprods/)
