      integer:: ibin(maxd+1,maxhisto),icont(maxhisto)
      real(dp):: xbinwgt(maxd+1,maxhisto),binmin(maxhisto),
     & bindel(maxhisto)
      common/smarthisto/ibin,icont,xbinwgt,binmin,bindel
!$omp threadprivate(/smarthisto/)
