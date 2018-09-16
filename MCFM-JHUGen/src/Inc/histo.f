      include 'nplot.f'
      integer maxnbin
      parameter(maxnbin=200)
      real(dp)::HIST(nplot,maxnbin),XHIS(nplot,maxnbin),
     & HDEL(nplot),
     & HMIN(nplot),HMAX(nplot),HAVG(nplot),HINT(nplot),HSIG(nplot)
      COMMON/HISTOR/HIST,XHIS,HDEL,HMIN,HMAX,HAVG,HINT,HSIG

      CHARACTER TITLE*100,BOOK*3
      COMMON/HISTOC/BOOK(nplot),TITLE(nplot)
      INTEGER NBIN(nplot),NHIST
      integer(kind=8) IHIS(nplot,maxnbin),IUSCORE(nplot),
     & IOSCORE(nplot),IENT(nplot)
      COMMON/HISTOI/IHIS,IUSCORE,IOSCORE,IENT,NBIN,NHIST

      integer, parameter:: tagbook=1, tagplot=2

