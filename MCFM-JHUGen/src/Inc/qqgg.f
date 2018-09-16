      integer:: qq,qg,gq,gg
      parameter (qq=1,qg=2,gq=3,gg=4)
      logical:: qqproc,qgproc,gqproc,ggproc
      common/dipproc/qqproc,qgproc,gqproc,ggproc
!$omp threadprivate(/dipproc/)
