c--- Modification of mbook.f with different SUBROUTINE TMPnames
c--- "TMPxxxxx" and threadprivate histogram common blocks

**********************************************************************
C    SIMPLE HISTOGRAMMING PACKAGE --  SIMPLIFIED VERSION OF HBOOK
C    BY Michelangelo Mangano    NOVEMBER 1988
C    LAST REVISED NOVEMBER 9, 1988  
c     (minor modifications by I Hinchliffe   1 May, 89)
C**********************************************************************
C
C Fills up to 100 histograms with up to 100 bins. 
C Gives a data file (to be specified in the calling program by assigning 
C a file name to unit 98) and a topdrawer file (to be specified in the 
C calling program by assigning a file name to unit 99).
C
C INITIALIZATION:
C Call once INIHIST; this just resets a few counters and logical::s
C Call MBOOK(N,'TITLE',DEL,XMIN,XMAX) for each histogram to be booked.
C N (an integer::) is the label of the histogram;
C 'TITLE' is the name of the histogram (no more then 100 characters);
C DEL (real*8) is the bin size;
C XMIN (real*8) is the lower limit of the first bin;
C XMAX (real*8)is the upper limit of the last  bin
C Example:
C      call mbook(2,'pt distribution',1.,10.,70.)
C This call initializes histogram number 2, called 'pt distribution';
C The bin size will be 1. (possibly GeV, if that's what you want), the
C first bin being  10.<x<11. and the last one being 69.<x<70.
C
C FILLING:
C When it's time, call MFILL(N,X,Y); this will add Y (real*8) to the bin 
C in which X (real*8) happens to be, within histogram N. 
C
C--------------------------------------------------------------------------
      SUBROUTINE TMPMBOOK(N,TIT,DEL,XMIN,XMAX)
      
c      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT integer:: (I-N)
      include 'types.f'
      integer:: N,I,NNBIN
      real(dp):: DEL,XMIN,XMAX
      CHARACTER*(*) TIT
      include 'tmphisto.f'
      TMPNHIST=MAX(N,TMPNHIST)
      TMPTITLE(N)=TIT
      TMPBOOK(N)='YES'
      TMPHDEL(N)=DEL
      TMPHMIN(N)=XMIN
      TMPHMAX(N)=XMAX+1d-8
      NNBIN=INT((XMAX+1d-8-XMIN)/DEL)
      IF (NNBIN .GT. 100) THEN
      WRITE(6,*) XMAX,XMIN,DEL,NNBIN,' BIN SIZE TOO LARGE'
      DEL=(XMAX-XMIN)/99.d0
      NNBIN=INT((XMAX-XMIN)/DEL)
      ENDIF
      TMPNBIN(N)=NNBIN
      TMPIENT(N)=0
      TMPIUSCORE(N)=0
      TMPIOSCORE(N)=0
      TMPHAVG(N)=0.d0
      TMPHINT(N)=0.d0
      DO 1 I=1,TMPNBIN(N)
      TMPXHIS(N,I)=TMPHMIN(N)+TMPHDEL(N)*(real(I,dp)-0.5d0)
      TMPIHIS(N,I)=0
   1  TMPHIST(N,I)=0.d0
      END

      SUBROUTINE TMPMFILL(N,X,Y)
      
c      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT integer:: (I-N)
      include 'types.f'
      integer:: N,I
      real(dp):: X,Y
      include 'tmphisto.f'
      I=INT((X-TMPHMIN(N))/TMPHDEL(N)+1)
      IF(I.GT.0.AND.I.LE.TMPNBIN(N))  THEN
      TMPIENT(N)=TMPIENT(N)+1
      TMPIHIS(N,I)=TMPIHIS(N,I)+1
      TMPHIST(N,I)=TMPHIST(N,I)+Y/TMPhdel(n)
c     we are renormalising the weights by the bin width
      ELSEIF(I.LE.0) THEN
      TMPIUSCORE(N)=TMPIUSCORE(N)+1
      ELSEIF(I.GT.TMPNBIN(N)) THEN
      TMPIOSCORE(N)=TMPIOSCORE(N)+1
      ENDIF
      END

     
      SUBROUTINE TMPMZERO(N)
      
c      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT integer:: (I-N)
      include 'types.f'
      integer:: N,I
      include 'tmphisto.f'
      TMPBOOK(N)='RES'
      TMPIENT(N)=0
      TMPIUSCORE(N)=0
      TMPIOSCORE(N)=0
      TMPHAVG(N)=0.d0
      TMPHINT(N)=0.d0
      DO 1 I=1,TMPNBIN(N)
      TMPIHIS(N,I)=0
   1  TMPHIST(N,I)=0.d0
      END

