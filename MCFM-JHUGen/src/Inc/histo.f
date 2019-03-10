      include 'nplot.f'
      integer maxnbin
      parameter(maxnbin=200)
      DOUBLE PRECISION HIST(nplot,maxnbin),XHIS(nplot,maxnbin),
     & HDEL(nplot),
     & HMIN(nplot),HMAX(nplot),HAVG(nplot),HINT(nplot),HSIG(nplot)
      COMMON/HISTOR/HIST,XHIS,HDEL,HMIN,HMAX,HAVG,HINT,HSIG

      CHARACTER TITLE*100,BOOK*3
      COMMON/HISTOC/BOOK(nplot),TITLE(nplot)                         
      INTEGER NBIN(nplot),IHIS(nplot,maxnbin),IUSCORE(nplot),
     & IOSCORE(nplot),IENT(nplot),NHIST
      COMMON/HISTOI/NBIN,IHIS,IUSCORE,IOSCORE,IENT,NHIST
