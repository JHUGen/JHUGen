       SUBROUTINE MRSEB(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C***************************************************************C
C                                                               C
C     This is a package for the new MRS(A',G) parton            C
C     distributions. The minimum Q^2  value is 5 GeV^2 and the  C
C     x range is, as before 10^-5 < x < 1. MSbar factorization  C
C     is used. The package reads 2 grids, which are in separate C
C     files (A'=for020.dat/ftn20, G=for021.dat/ftn21).          C  
C     Note that x times the parton distribution is returned,    C
C     Q is the scale in GeV,                                    C
C     and Lambda(MSbar,nf=4) = 231/255 MeV for A'/G.            C
C                                                               C
C       MODE=20 for MRS(A')                                      C
C       MODE=21 for MRS(G)                                       C
C                                                               C
C         The reference is :                                    C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         Phys. Lett. B354 (1995) 155-162                       C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE 20 GRID IS 0.00341            C
C         THE FIRST NUMBER IN THE 21 GRID IS 0.00269            C
C                                                               C
C***************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      Q2=SCALE**2
      IF(Q2.LT.5D0.OR.Q2.GT.1310720.D0)    PRINT 99
      IF(X.LT.1D-5.OR.X.GT.1D0)            PRINT 98
      IF(MODE.EQ.20) 
     .   CALL STRC20(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.21) 
     .   CALL STRC21(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
  99  FORMAT('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  FORMAT('  WARNING:   X  VALUE IS OUT OF RANGE   ')
      RETURN
      END
C
      SUBROUTINE STRC20(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C     THIS IS THE NEW  "Aprime" FIT -- Feb 1995 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
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
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/
      save XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT
!$omp threadprivate(XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT)
 
 
      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      filename=checkpath('Pdfdata/mrs95ap.dat')
      OPEN(UNIT=20,file=filename,status='unknown')
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(20,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
      close(unit=20)
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
      DO 60 I=1,8
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
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C
      SUBROUTINE STRC21(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)

C     THIS IS THE NEW  "G" FIT -- Feb 1995 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
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
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/
      save XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT
!$omp threadprivate(XX,F,XMIN,XMAX,QSQMIN,QSQMAX,N0,INIT)
 
 
      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      filename=checkpath('Pdfdata/mrs95_g.dat')
      OPEN(UNIT=21,file=filename,status='unknown')
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(21,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
      close(unit=21)
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
      DO 60 I=1,8
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
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C
