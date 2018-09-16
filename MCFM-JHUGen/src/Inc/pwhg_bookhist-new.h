c -*- Fortran -*-
      integer nhist,maxbins
      parameter (nhist=400,maxbins=500)
      character * 50 stringhist(nhist)
      real * 8 xhistarr(maxbins+1,nhist),yhistarr(0:maxbins+1,nhist)
      integer nhits(0:maxbins+1,nhist),nbins(nhist),jhist,ient1(nhist)
      real * 8 yhistarr1(0:maxbins+1,nhist),
     1       errhistarr1(0:maxbins+1,nhist),
     1       yhistarr2(0:maxbins+1,nhist),
     3       errhistarr2(0:maxbins+1,nhist)
      common/histnew/xhistarr,yhistarr,
     1       yhistarr1,errhistarr1,
     3       yhistarr2,errhistarr2,
     5       nhits,nbins,jhist,ient1,
     6       stringhist
      save /histnew/

