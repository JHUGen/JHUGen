C**********************************************************************
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
C Call once INIHIST; this just resets a few counters and logicals
C Call MBOOK(N,'TITLE',DEL,XMIN,XMAX) for each histogram to be booked.
C N (an integer) is the label of the histogram;
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
C PLAYING AROUND:
C At the end of the day you may want to sum, divide, cancel, etc.etc.
C various histograms (bin by bin). Then you call MOPERA(I,'O',J,K,X,Y). 
C The 1-character string O can take the following values:
C +  : sums       X*(hist I) with Y*(hist J) and puts the result in hist K;
C -  : subtracts  X*(hist I) with Y*(hist J) and puts the result in hist K;
C *  : multiplies X*(hist I) with Y*(hist J) and puts the result in hist K;
C /  : divides    X*(hist I) with Y*(hist J) and puts the result in hist K;
C F  : multiplies hist I by the factor X, and puts the result in hist K;
C R  : takes the square root of  hist  I, and puts the result in hist K;if
C      the value at a given bin is less than or equal to 0., puts 0. in K
C S  : takes the square      of  hist  I, and puts the result in hist K;
C L  : takes the log_10 of  hist  I, and puts the result in hist K; if the
C      value at a given bin is less than or equal to 0., puts 0. in K
C M  : statistical analysis; if I contains the weights (let's say WGT),
C      J contains variable times weight (F*WGT) and K contains the
C      variable squared times the weight (F**2*WGT), then, after using 'M',
C      J will contain the average value of the variable <F> and K will 
C      contain the sigma of the average: sigma=sqrt(<F**2>-<F>**2).
C      If WGT=1. for all the entries, then it is enough to put I=J, and
C      it is not necessary to book a hist with the weights.
C V  : estimates errors for vegas evaluation of differential distributions.
C      Fill I with the values of
C      the functions do integrate times the Vegas weight (fun*wgt); fill
C      J with fun**2*wgt; then K will contain an estimate of the error
C      of the integration. Putting X=1/(#of iterations) performs the 
C      average over the iterations, and gives the right normalization to 
C      the differential distribution, I, and to the errors, K. J stays the same.
C U  : same as V, but I also remains the same (i.e. only K is filled)
C
C FINAL ACCOUNTING:
C Now we can finalize our histograms; MFINAL(N) will calculate the integral
C of the histogram N, the mean value of the X variable and its RMS.
C If we now want to renormalize the hist's, we can call MNORM(N,X), which
C will normalize the integral to X  -- CAUTION: do not call MNORM before
C MFINAL, it will blow up.
C
C OUTPUT:
C To get a .dat file containing the values of the histograms, together with
C some information (like integral, mean values, etc.etc.) call MPRINT(N),
C for each hist N that you want in the .dat file. Before the call to MPRINT
C you want to open unit 98 and give it a name:
C     OPEN(UNIT=98,NAME='NAME.DAT',STATUS='NEW')
C If you want a topdrawer file with a plot of the hist values, call 
C MTOP(N,M,'X','Y','SCALE'). The points of the plot will be taken from histogram
C N, the error bars from histogram M. 'SCALE', character*(*), determines
C the scale for y, logarithmic or linear (SCALE=LOG,LIN).
C If you do not want error bars, keep
C a histogram of zeros, or just call a hist that had not been booked.
C X will appear as a 'bottom title', and Y will appear as a 'left title'.
C The top title is by default the name of the histogram itself.
C A little box below the plot will contain some information on the plot
C itself. Before calling MTOP,                     
C     OPEN(UNIT=99,NAME='NAME.TOP',STATUS='NEW')
c Empty histograms are not put out by MTOP.
C--------------------------------------------------------------------------

      SUBROUTINE MBOOK(N,TIT,DEL,XMIN,XMAX)
      implicit none
      include 'types.f'
      real(dp) DEL,XMIN,XMAX
      integer N,I,NNBIN
      CHARACTER*(*) TIT
      include 'histo.f'
      NHIST=MAX(N,NHIST)
      TITLE(N)=TIT                     
      BOOK(N)='YES'
      HDEL(N)=DEL
      HMIN(N)=XMIN
      HMAX(N)=XMAX+1d-8
      NNBIN=INT((XMAX+1d-8-XMIN)/DEL)
      IF (NNBIN .GT. 100) THEN
      WRITE(6,*) XMAX,XMIN,DEL,NNBIN,' BIN SIZE TOO LARGE'
      DEL=(XMAX-XMIN)/99._dp
      NNBIN=INT((XMAX-XMIN)/DEL)
      ENDIF
      NBIN(N)=NNBIN
      IENT(N)=0
      IUSCORE(N)=0
      IOSCORE(N)=0
      HAVG(N)=0._dp
      HINT(N)=0._dp
      DO 1 I=1,NBIN(N)
      XHIS(N,I)=HMIN(N)+HDEL(N)*(real(I,dp)-0.5_dp)
      IHIS(N,I)=0
   1  HIST(N,I)=0._dp
      END

      SUBROUTINE MFILL(N,X,Y)
      implicit none
      include 'types.f'
      integer N,I,K,L
      real(dp) X,Y,Yadd
      include 'histo.f'
      I=INT((X-HMIN(N))/HDEL(N)+1)
      IF(I.GT.0.AND.I.LE.NBIN(N))  THEN
!$omp atomic
      IENT(N)=IENT(N)+1
!$omp atomic
      IHIS(N,I)=IHIS(N,I)+1
      Yadd=Y/hdel(n)
!$omp atomic
      HIST(N,I)=HIST(N,I)+Yadd
c     we are renormalising the weights by the bin width
      ELSEIF(I.LE.0) THEN
!$omp atomic
      IUSCORE(N)=IUSCORE(N)+1
      ELSEIF(I.GT.NBIN(N)) THEN
!$omp atomic
      IOSCORE(N)=IOSCORE(N)+1
      ENDIF
      END

      SUBROUTINE MOPERA(I,OPER,J,K,X,Y)
      implicit none
      include 'types.f'
      integer I,J,K,L
      real(dp) XAVG,XSQAVG,XNORM,X,Y
      CHARACTER OPER*1
      include 'histo.f'
      IF(NBIN(I).NE.NBIN(J).AND.(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.
     &'*'.OR.OPER.EQ.'/'.OR.OPER.EQ.'M')) GO TO 10
      DO L=1,NBIN(I)
      IF(OPER.EQ.'+') THEN
      HIST(K,L)=X*HIST(I,L) + Y*HIST(J,L)
      ELSEIF(OPER.EQ.'-') THEN
      HIST(K,L)=X*HIST(I,L) - Y*HIST(J,L)
      ELSEIF(OPER.EQ.'*') THEN
      HIST(K,L)=X*HIST(I,L) * Y*HIST(J,L)
      ELSEIF(OPER.EQ.'/') THEN
        IF(Y.EQ.0._dp.OR.HIST(J,L).EQ.0._dp) THEN
          HIST(K,L)=0._dp
          ELSE
          HIST(K,L)=X*HIST(I,L) / (Y*HIST(J,L))
        ENDIF
      ELSEIF(OPER.EQ.'F') THEN
      HIST(K,L)=X*HIST(I,L)
      ELSEIF(OPER.EQ.'R') THEN
        IF(HIST(I,L).GT.0._dp) THEN
        HIST(K,L)=X*SQRT(HIST(I,L))
        ELSE
        HIST(K,L)=0._dp
        ENDIF
      ELSEIF(OPER.EQ.'S') THEN
      HIST(K,L)=X*HIST(I,L)**2
      ELSEIF(OPER.EQ.'l') THEN
        IF(HIST(I,L).EQ.0._dp.OR.J.EQ.0) THEN
             HIST(K,L)=0._dp
             ELSE
             HIST(K,L)=X*LOG10(Y*HIST(I,L))
        ENDIF
      ELSEIF(OPER.EQ.'M') THEN
        IF(I.NE.J) XNORM=HIST(I,L)
        IF(I.EQ.J) XNORM=real(IHIS(J,L),dp)
        IF(XNORM.NE.0._dp) THEN
        XAVG=HIST(J,L)/XNORM
        HIST(K,L)=sqrt(ABS(-XAVG**2+HIST(K,L)/XNORM)/real(IHIS(I,L),dp))
        HIST(J,L)=XAVG 
        ELSE
        HIST(K,L)=0._dp
        HIST(J,L)=0._dp                           
        ENDIF
      ELSEIF(OPER.EQ.'V') THEN                 
        XAVG=HIST(I,L)*X
        XSQAVG=HIST(J,L)*X
c--- need extra factor to account for renormalization by bin width
        XSQAVG=XSQAVG/hdel(i)
        XNORM=(IENT(I)+IUSCORE(I)+IOSCORE(I))*X
        IF(XNORM.NE.0._dp) THEN
        HIST(K,L)=sqrt(X*ABS(XSQAVG-XAVG**2/XNORM))
        HIST(I,L)=XAVG
        ELSE
        HIST(K,L)=0._dp
        ENDIF
      ELSEIF(OPER.EQ.'U') THEN ! same as 'V', but write errors only   
        XAVG=HIST(I,L)*X
        XSQAVG=HIST(J,L)*X
c--- need extra factor to account for renormalization by bin width
        XSQAVG=XSQAVG/hdel(i)
        XNORM=(IENT(I)+IUSCORE(I)+IOSCORE(I))*X
        IF(XNORM.NE.0._dp) THEN
        HIST(K,L)=sqrt(X*ABS(XSQAVG-XAVG**2/XNORM))
c        HIST(I,L)=XAVG ! removed from 'V'
        ELSE
        HIST(K,L)=0._dp
        ENDIF
      ELSE
      WRITE(98,5) OPER
   5  FORMAT(' ****** OPERATION ="',A1,'" UNKNOWN ********'/)
      RETURN
      ENDIF
      END DO
      RETURN
  10  WRITE(98,20) I,J
  20  FORMAT(' ****** INCOMPATIBLE OPERATION HIST ',I3,' &',I3,
     &                                                   '*******'/)
      END
     
      SUBROUTINE MZERO(N)
      implicit none
      include 'types.f'
      integer N,I
      include 'histo.f'
      BOOK(N)='RES'
      IENT(N)=0
      IUSCORE(N)=0
      IOSCORE(N)=0
      HAVG(N)=0._dp
      HINT(N)=0._dp
      DO 1 I=1,NBIN(N)
      IHIS(N,I)=0
   1  HIST(N,I)=0._dp
      END

      SUBROUTINE MRESET(N)
      implicit none
      include 'types.f'
      integer N
      include 'histo.f'
      BOOK(N)='RES'
      END

      SUBROUTINE MFINAL(N)
      implicit none
      include 'types.f'
      real(dp) AVG,XIN
      integer N,J
      include 'histo.f'
      IF(BOOK(N).NE.'YES') RETURN
      AVG=0._dp
      XIN=0._dp                                
c      SIG=0._dp
      DO 1, J=1,NBIN(N)
      AVG=AVG+HIST(N,J)*XHIS(N,J)
   1  XIN=XIN+HIST(N,J)
      IF(XIN.EQ.0._dp) GO TO 10
      HAVG(N)=AVG/XIN
c      DO 2, J=1,NBIN(N)
c   2  SIG=HIST(N,J)*(XHIS(N,J)-HAVG(N))**2+SIG
c      IF(SIG.GE.0.)HSIG(N)=SQRT(SIG/XIN)
      HINT(N)=XIN*hdel(n)
      RETURN
  10  BOOK(N)=' NO'
      END               

      SUBROUTINE MNORM(N,X)
      implicit none
      include 'types.f'
      real(dp) X
      integer N,I
      include 'histo.f'
      IF(BOOK(N).NE.'YES')RETURN
      DO 1, I=1,NBIN(N)
    1 HIST(N,I)=HIST(N,I)/HINT(N)*X
      HINT(N)=X
      END

      SUBROUTINE MPRINT(N,L)
      implicit none
      include 'types.f'
      integer N,L,I
      include 'histo.f'
      logical scaleplots                  
      real(dp) scalefac
      common/scaleplots/scalefac,scaleplots
c      DATA INI/0/
c      IF(INI.EQ.0) THEN
c      CALL IDATE(IMON,IDAY,IYEAR)
c      CALL TIME(CTIME)
c      INI=1
c      ENDIF
      IF(BOOK(N).NE.'YES') then
      write(98,21) n
      RETURN
      end if
c      WRITE(98,7) N,IYEAR,IMON,IDAY,CTIME(1:5)
      WRITE(98,8) N
      WRITE(98,*) TITLE(N)
      DO I=1,NBIN(N)
      if (scaleplots) then
      WRITE(98,10) XHIS(N,I),scalefac*HIST(N,I),HIST(L,I)
      else
      WRITE(98,10) XHIS(N,I),HIST(N,I),HIST(L,I)
      endif
      ENDDO
      WRITE(98,15) HAVG(N),HSIG(N),HINT(N)
      WRITE(98,20) IENT(N),IUSCORE(N),IOSCORE(N)
c    7 FORMAT(4X,'HIST = ',I3,'   19',I2,'-',I2,'-',I2,1X,A5/)
    8 FORMAT(4X,'HIST = ',I3)
   10 FORMAT(3(3x,G13.6))
   15 FORMAT(/' AVG =',E10.3,4X,' RMS =',E10.3,' INTEGRAL =',E10.3,/)
   20 FORMAT('ENTRIES=',I10,1X,'U`FLOW=',I10,1X,'O`FLOW=',I10,//)
   21 FORMAT(' HISTOGRAM ',I3,' IS EMPTY')
      END

      SUBROUTINE MTOP(N,M,BTIT,LTIT,SCALE)
      implicit none
      include 'types.f'
      integer N,M,J
      CHARACTER*(*) LTIT,BTIT,SCALE
      include 'histo.f'
c--- added these variables to scale plots at intermediate steps
      logical scaleplots                  
      double precision scalefac
      common/scaleplots/scalefac,scaleplots
c      DATA INI/0/
c      IF(INI.EQ.0) THEN
c      CALL IDATE(IMON,IDAY,IYEAR)
c      CALL TIME(CTIME)
c      INI=1
c      ENDIF
      
      IF(BOOK(N).NE.'YES') RETURN
c      WRITE(99,100) TITLE(N),BTIT,LTIT,SCALE,HMIN(N),HMAX(N)
      WRITE(99,101) TITLE(N),TITLE(N),TITLE(N),SCALE,HMIN(N),HMAX(N)
c  100 FORMAT( /1x,                               
c     &' SET WINDOW Y 2.5 TO 7.'/,1X,
c     &' SET WINDOW X 2.5 TO 10.'/,1X,
c     &' SET FONT DUPLEX '/1X,
c     &' SET SYMBOL 5O SIZE 1.8'/,1X,
c     &' TITLE TOP ','"',A50,'"',/1X,
c     &' TITLE BOTTOM ','"',A50,'"',/1X,
c     &' TITLE LEFT ','"',A50,'"',/1X,
c     &' SET SCALE Y ',A5,/1X,
c     &' (SET TICKS TOP OFF)   '/1x,     
c     &' SET LIMITS X ',F12.5,' ',F12.5,/1X,
c     &' SET ORDER X Y DY ')
  100 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A50,'"',/1X,
     &' TITLE BOTTOM ','"',A50,'"',/1X,
     &' TITLE LEFT ','"',A50,'"',/1X,
     &' SET SCALE Y ',A5,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F12.5,' ',F12.5,/1X,
     &' SET ORDER X Y DY')
  101 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP SIZE=3','"',A,/1X,
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE LEFT ','"dS/d',A,' [fb]"',/1X,
     &' CASE       ','" G"',/1X,
     &' SET SCALE Y ',A5,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F12.5,' ',F12.5,/1X,
     &' SET ORDER X Y DY')
      DO 1 J=1,NBIN(N)
c      IF(HIST(N,J).EQ.0.) GO TO 1
      if (scaleplots) then
      WRITE(99,'(3X,G13.6,2(2X,G13.4))')  
     & XHIS(N,J),scalefac*HIST(N,J),HIST(M,J)
      else
      WRITE(99,'(3X,G13.6,2(2X,G13.4))')  
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
      endif
    1 CONTINUE
      WRITE(99,200)
  200 FORMAT('   PLOT')
      if (scaleplots) then
      WRITE(99,300) scalefac*HINT(N),scalefac*HAVG(N),scalefac*HSIG(N),
     & IENT(N),IUSCORE(N),IOSCORE(N)
      else
      WRITE(99,300) HINT(N),HAVG(N),HSIG(N),IENT(N),IUSCORE(N)
     &   ,IOSCORE(N)
      endif
c  300 FORMAT( /1x,                               
c     &' BOX 7. 0.75 SIZE 9. 1.5'/,1X,
c     &' SET WINDOW Y 0. TO 2.'/,1X,
c     &' SET TITLE SIZE -1.5'/1X,
c     &' SET FONT DUPLEX '/1X,
c     &' TITLE 2.8 1.2 "INTEGRAL =',E10.3,'   AVERAGE =',E10.3,
c     &             '   RMS =',E10.3,'"',/1X,
c     &' TITLE 2.8 0.8 "Entries =',I10,4x,'Underflow =',I10,4X
c     &                                 ,'Overflow =',I10,'"',/1X,
c     &' SET TITLE SIZE -2')
  300 FORMAT( /1x,                               
     &' BOX 7. 0.75 SIZE 9. 1.5'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' TITLE 2.8 1.2 "INTGRL =',E12.5,'   AVGE =',E12.5,
     &             '   RMS =',E12.5,'"',/1X,
     &' TITLE 2.8 0.8 "Entries =',I9,2x,'U`flow =',I9,2X
     &                                 ,'O`flow =',I9,'"',/1X,
     &' SET TITLE SIZE -2')
      WRITE(99,400)
  400 FORMAT('   NEW PLOT')
      END


C*******************************************************************
C     END OF THE HISTOGRAMMING PACKAGE
C*******************************************************************


c--F  Add gnuplot output
      
      SUBROUTINE MGNUPLOT(N,M,BTIT,LTIT,SCALE)
      implicit none
      include 'types.f'
      integer N,M,idlength,LEN_TRIM,J,istring
      CHARACTER*(*) LTIT,BTIT,SCALE
      include 'histo.f'
c--- added these variables to scale plots at intermediate steps
      logical scaleplots                  
      double precision scalefac
      common/scaleplots/scalefac,scaleplots

      istring=LEN_TRIM(TITLE(N))

      IF(BOOK(N).NE.'YES') RETURN
      WRITE(97,101) TITLE(N)(1:istring), TITLE(N)(1:istring),
     & TITLE(N)(1:istring), HMIN(N),HMAX(N)
  101 FORMAT( /1x,                               
     &' set title ','"',A,' distribution" font "Helvetica, 20"',/1X,
     &' set xlabel ','"',A,'" font "Helvetica, 20"',/1X,
     &' set ylabel ','"d{/Symbol s}/d',A,
     &' [fb]" font "Helvetica, 20"',/1X,
     &' set xrange [ ',F10.5,':',F10.5,']')
      if(SCALE .eq. 'log') then
         write(97,*) ' set logscale y'
      endif
      write(97,*) ' plot "-" with histeps'

      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      if (scaleplots) then
      WRITE(97,'(3(2X,G13.6))')  
     & XHIS(N,J),scalefac*HIST(N,J),HIST(M,J)
      else
      WRITE(97,'(3X,G13.6,2(2X,G13.6))')  
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
      endif
    1 CONTINUE
      WRITE(97,200)
  200 FORMAT(' e')
      END



c-- Root output (works in versions <=5)
      SUBROUTINE MROOTPLOT5(N,M,BTIT,LTIT)
      implicit none
      include 'types.f'
      integer N,M,idlength,LEN_TRIM,J,istring
      CHARACTER*(*) LTIT,BTIT
      CHARACTER*5 histoid
      include 'histo.f'
c--- added these variables to scale plots at intermediate steps
      logical scaleplots                  
      double precision scalefac
      common/scaleplots/scalefac,scaleplots

      istring=LEN_TRIM(TITLE(N))

      if (N.lt.10) then
         write(histoid, '(A2,I1,A2)') 'id', N, '  '
      elseif ((N.ge.10).and.(N.lt.100)) then
         write(histoid, '(A2,I2,A1)') 'id', N, ' '
      elseif (N.ge.100) then
         write(histoid, '(A2,I3)') 'id', N
      endif

      idlength=LEN_TRIM(histoid)

      IF(BOOK(N).NE.'YES') RETURN

      WRITE(96,131) histoid(1:idlength), TITLE(N)(1:istring), NBIN(N), 
     & HMIN(N), HMAX(N)
 131  FORMAT ( /1X,
     & ' mcfmhisto -> cd();', /1X,
     & ' TH1F *hist = new TH1F( "', A, '", "', A, '", ',
     & I0, ', ', F10.5, ', ', F10.5, ');')

      WRITE(96, 132) histoid(1:idlength), TITLE(N)(1:istring),
     & histoid(1:idlength), TITLE(N)(1:istring)
 132  FORMAT ( /1X, 
     & ' ', A, ' -> GetXaxis() -> SetTitle("', A, '");', /1X,
     & ' ', A, ' -> GetYaxis() -> SetTitle(" d#sigma/d', A, 
     & ' [fb]");', /1X)

      WRITE (96,*) ' ', histoid(1:idlength), ' -> GetYaxis() -> ',
     & 'SetTitleOffset(1.2);'

      WRITE(96,*) ' ', histoid(1:idlength), ' -> SetStats(false);'

      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      if (scaleplots) then
      WRITE(96,'(3(2X,G13.6))')  
     & XHIS(N,J),scalefac*HIST(N,J),HIST(M,J)
      else
C         write(96,*) ' ', histoid(1:idlength), ' -> Fill(', 
C     &        XHIS(N,J), ', ', HIST(N,J), ');'
         write(96,*) '  int xbin = ', histoid(1:idlength),'->FindBin(',
     & XHIS(n,j),');' 
         write(96,*) ' ', histoid(1:idlength), ' -> SetBinContent(', 
     &       ' xbin', ', ', HIST(N,J), ');'
         write(96,*) ' ', histoid(1:idlength), ' -> SetBinError(', 
     &       ' xbin', ', ', HIST(M,J), ');'
      endif
    1 CONTINUE

      WRITE (96, *) ' histos -> Add(hist); '
      WRITE (96, *) ''
      WRITE (96, *) ''

      END


c-- Root output (works in versions >=5)
      SUBROUTINE MROOTPLOT(N,M,BTIT,LTIT)
      implicit none
      include 'types.f'
      integer N,M,J,idlength,LEN_TRIM,istring
      CHARACTER*(*) LTIT,BTIT
      CHARACTER*5 histoid
      include 'histo.f'
c--- added these variables to scale plots at intermediate steps
      logical scaleplots                  
      double precision scalefac
      common/scaleplots/scalefac,scaleplots

      istring=LEN_TRIM(TITLE(N))

      if (N.lt.10) then
         write(histoid, '(A2,I1,A2)') 'id', N, '  '
      elseif ((N.ge.10).and.(N.lt.100)) then
         write(histoid, '(A2,I2,A1)') 'id', N, ' '
      elseif (N.ge.100) then
         write(histoid, '(A2,I3)') 'id', N
      endif

      idlength=LEN_TRIM(histoid)

      IF(BOOK(N).NE.'YES') RETURN

      WRITE(96,141) histoid(3:idlength), histoid(1:idlength),
     &  TITLE(N)(1:istring), NBIN(N), HMIN(N), HMAX(N)
 141  FORMAT ( /1X,
     & ' mcfmhisto -> cd();', /1X,
     & ' TH1F *hist', A, ' = new TH1F( "', A, '", "', A, '", ',
     & I0, ', ', F12.4, ', ', F12.4, ');')

      WRITE(96, 142) histoid(1:idlength), TITLE(N)(1:istring),
     & histoid(1:idlength), TITLE(N)(1:istring)
 142  FORMAT ( /1X, 
     & ' ', A, ' -> GetXaxis() -> SetTitle("', A, '");', /1X,
     & ' ', A, ' -> GetYaxis() -> SetTitle(" d#sigma/d', A, 
     & ' [fb]");', /1X)

      WRITE (96,*) ' ', histoid(1:idlength), ' -> GetYaxis() -> ',
     & 'SetTitleOffset(1.2);'

      WRITE(96,*) ' ', histoid(1:idlength), ' -> SetStats(false);'

      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      if (scaleplots) then
      WRITE(96,'(3(2X,G13.6))')  
     & XHIS(N,J),scalefac*HIST(N,J),HIST(M,J)
      else
C         write(96,*) ' ', histoid(1:idlength), ' -> Fill(', 
C     &        XHIS(N,J), ', ', HIST(N,J), ');'
         write(96,*) '  xbin = ', histoid(1:idlength),'->FindBin(',
     & XHIS(n,j),');' 
         write(96,*) ' ', histoid(1:idlength), ' -> SetBinContent(', 
     &       ' xbin', ', ', HIST(N,J), ');'
         write(96,*) ' ', histoid(1:idlength), ' -> SetBinError(', 
     &       ' xbin', ', ', HIST(M,J), ');'
      endif
    1 CONTINUE

      WRITE (96, 143) histoid(3:idlength)
 143  FORMAT ('//  hist', A, ' -> Draw("hist");')
      WRITE (96, *) ''
      WRITE (96, *) ''

      END

