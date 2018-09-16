       SUBROUTINE MRSEBH(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C***************************************************************C
C                                                               C
C                                                               C
C     NEW VERSIONS !!!! JANUARY 1990  (AS DESCRIBED IN          C
C     "PARTON DISTRIBUTIONS ... " P.N. HARRIMAN, A.D. MARTIN,   C
C     R.G. ROBERTS AND W.J. STIRLING PREPRINT DTP-90-04 )       C
C                                                               C
C         ********* DEBUGGED APRIL 1990********                 C 
C                                                               C
C         ****** NOW DOWN TO X=10^-5 **********                 C
C                                                               C
C  MODE 1 CORRESPONDS TO  HARRIMAN,                             C
C  MARTIN, ROBERTS, STIRLING (EMC FIT)    WITH LAMBDA= 100 MEV  C
C                                                               C
C  MODE 2  CORRESPONDS TO HARRIMAN,                             C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C    THE FIRST NUMBER IN THE "E" GRID IS  .01969                C
C    THE FIRST NUMBER IN THE "B" GRID IS  .03058                C
C                                                               C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE LIGHT QUARK SEA I.E. UBAR(X)=DBAR(X)=               C
C      SEA/X FOR A PROTON.  IF IN DOUBT, CHECK THE              C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV)          C
C                                                               C
C                         -*-                                   C
C                                                               C
C     (THE RANGE OF APPLICABILITY IS CURRENTLY:                 C
C     10**-5 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO WJS@UK.AC.DUR.HEP             C
C                                                               C
C                                                               C
C***************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      IF(MODE.EQ.1) CALL STRC27(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.2) CALL STRC28(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      RETURN
      END
C
      SUBROUTINE STRC27(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU) 

C :::::::::::: HARRIMAN MARTIN ROBERTS STIRLING (E) sea=ubar=dbar=2sbar at Q0**2

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(7,nx,19),G(7),XX(nx),N0(7)
      character*72 filename,checkpath
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0,5/
      DATA INIT/0/
      save XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT
!$omp threadprivate(XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT)

      xsave=x  ! don't let x be altered if it's out of range!!

      IF(INIT.NE.0) GOTO 10
      filename=checkpath('Pdfdata/mrse.dat')
      open(unit=27,file=filename,status='old')
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(27,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M)
C 1=UV 2=DV 3=GLUE 4=(UBAR+DBAR)/2 5=CBAR 7=BBAR 6=SBAR
         DO 25 I=1,7
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,ntenth-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,6
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  50  FORMAT(7F10.5)
      DO 40 I=1,7
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
      close(unit=27)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,7
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4) ! THIS SEA IS (UBAR+DBAR)/2
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)

      x=xsave  !restore x

      RETURN
      END

      SUBROUTINE STRC28(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU) 

C :::::::::::: HARRIMAN MARTIN ROBERTS STIRLING (B) sea=ubar=dbar=2sbar at Q0**2

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(7,nx,19),G(7),XX(nx),N0(7)
      character*72 filename,checkpath
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0,5/
      DATA INIT/0/
      save XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT
!$omp threadprivate(XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT)


      xsave=x

      IF(INIT.NE.0) GOTO 10
      filename=checkpath('Pdfdata/mrsb.dat')
      open(unit=28,file=filename,status='old')
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(28,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M)
C 1=UV 2=DV 3=GLUE 4=(UBAR+DBAR)/2 5=CBAR 7=BBAR 6=SBAR
         DO 25 I=1,7
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,ntenth-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,6
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  50  FORMAT(7F10.5)
      DO 40 I=1,7
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
      close(unit=28)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,7
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4) ! THIS SEA IS (UBAR+DBAR)/2
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)

      x=xsave

      RETURN
      END
