      integer n2,n3 
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
!$omp threadprivate(/breit/)
