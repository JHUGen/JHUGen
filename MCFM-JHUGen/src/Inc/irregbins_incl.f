c--- Variable declarations and common blocks for producing
c--- histograms with irregular bin widths
      logical:: irregbin(maxhisto)
      integer:: nirreg,irregbin_ptr(maxhisto)
c--- Maximum 10 irregular histograms, with maximum 30 bins in each
      real(dp):: irregbinedges(10,30)
      common/irregbinmcfm/irregbinedges,irregbin,nirreg,irregbin_ptr

