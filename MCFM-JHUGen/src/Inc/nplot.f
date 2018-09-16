      integer:: nplot,maxhisto
      parameter(maxhisto=100,nplot=4*maxhisto)
      character*3 linlog(nplot)
      character*8 titlearray(nplot)
      common/topd/titlearray,linlog
      integer:: nextnplot
      common/plotindex/nextnplot
!$omp threadprivate(/plotindex/)
