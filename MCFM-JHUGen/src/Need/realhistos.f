      subroutine zerorealhistos
      implicit none
      include 'types.f'
c--- zero out all entries in the temporary histograms used for
c--- binning the weights in the real contribution
      
      include 'nplot.f'
      integer:: nplotmax,j
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(/nplotmax/)
      
      do j=1,nplotmax
      call tmpmzero(j)
      enddo
      
      return
      end
      

      subroutine addrealhistos(wgt)
      implicit none
      include 'types.f'
c--- add temporay histograms to the cumulative ones
       
c      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT integer:: (I-N)
      integer:: I,L
      include 'histo.f'
      real(dp):: wgt
      integer:: nplotmax
      logical:: added
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(/nplotmax/)

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

c--- loop over plots
      do I=1,nplotmax

      added=.false.

      DO L=1,NBIN(I)

c--- add weights
      HIST(I,L)=HIST(I,L) + TMPHIST(I,L)
c--- add errors
      HIST(maxhisto+I,L)=HIST(maxhisto+I,L)
     & + TMPHIST(I,L)**2*HDEL(I)
c--- the extra factor of HDEL(I) is to account for the normalization
c--- by the bin width (c.f. MFILL in mbook.f)
      
c--- count entries
      if (TMPIHIS(I,L) .GT. 0) then
        IHIS(I,L)=IHIS(I,L)+1
        IHIS(maxhisto+I,L)=IHIS(maxhisto+I,L)+1
        added=.true.
      endif

      ENDDO

      added=.true.

c--- if any bin has been filled, increment histogram counter
      if (added) then
        IENT(I)=IENT(I)+1
        IENT(maxhisto+I)=IENT(maxhisto+I)+1
c--- otherwise, increment out of bounds counters if necessary
      else
        if     (TMPIUSCORE(maxhisto+I) .GT. 0) then
          IUSCORE(I)=IUSCORE(I)+1      
          IUSCORE(maxhisto+I)=IUSCORE(maxhisto+I)+1
        elseif (TMPIOSCORE(maxhisto+I) .GT. 0) then
          IOSCORE(I)=IOSCORE(I)+1      
          IOSCORE(maxhisto+I)=IOSCORE(maxhisto+I)+1
        endif
      endif
      
      
      enddo
      
      return
      end
      
