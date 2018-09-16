      integer:: nplot,maxhisto
      parameter(maxhisto=100)
      integer:: maxnbin
      parameter(maxnbin=200)
      real(dp):: TMPHIST(maxhisto,maxnbin),
     & TMPXHIS(maxhisto,maxnbin),TMPHDEL(maxhisto),TMPHMIN(maxhisto),
     & TMPHMAX(maxhisto),TMPHAVG(maxhisto),TMPHINT(maxhisto),
     & TMPHSIG(maxhisto)
      COMMON/TMPHISTOR/TMPHIST,TMPXHIS,TMPHDEL,TMPHMIN,TMPHMAX,
     & TMPHAVG,TMPHINT,TMPHSIG

      CHARACTER TMPTITLE*100,TMPBOOK*3
      COMMON/TMPHISTOC/TMPBOOK(maxhisto),TMPTITLE(maxhisto)
      integer:: TMPNBIN(maxhisto),TMPIHIS(maxhisto,maxnbin),
     & TMPIUSCORE(maxhisto),TMPIOSCORE(maxhisto),TMPIENT(maxhisto),
     & TMPNHIST
      COMMON/TMPHISTOI/TMPNBIN,TMPIHIS,TMPIUSCORE,TMPIOSCORE,
     & TMPIENT,TMPNHIST
!$omp threadprivate(/TMPHISTOR/,/TMPHISTOC/,/TMPHISTOI/)
