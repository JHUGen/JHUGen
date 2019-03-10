      DOUBLE PRECISION FUNCTION ALPHAS(Q,AMZ,NLOOP)
c     Evaluation of strong coupling constant alpha_S
c     Author: R.K. Ellis

c     q -- scale at which alpha_s is to be evaluated
c     amz -- value of alpha_s at the mass of the Z-boson
c     nloop -- the number of loops (1,2, or 3) at which beta 
c     function is evaluated to determine running.
c     the values of the cmass and the bmass should be set
c     in common block qmass.

      IMPLICIT NONE
c--- added by JMC: use consistent value of MZ from common block
      include 'masses.f'
      DOUBLE PRECISION Q,T,AMZ,AMZ0,AMB,AMC,BMASS,CMASS,AS_OUT
      INTEGER NLOOP,NLOOP0,NF3,NF4,NF5
c      PARAMETER(ZMASS=91.188D0)
      PARAMETER(NF5=5,NF4=4,NF3=3)
      COMMON/QMASS/CMASS,BMASS
      SAVE AMZ0,NLOOP0,AMB,AMC
      DATA AMZ0,NLOOP0/0D0,0/

          

      IF (Q .LE. 0D0) THEN 
         WRITE(6,*) 'q .le. 0 in alphas'
         WRITE(6,*) 'q= ',Q
	 stop
c         PAUSE
      ENDIF
      IF (AMZ .LE. 0D0) THEN 
         WRITE(6,*) 'amz .le. 0 in alphas',AMZ
         WRITE(6,*) 'continuing with amz=0.117'
c         PAUSE
         AMZ=0.117D0
      ENDIF
      IF (CMASS .LE. 0.3D0) THEN 
         WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
         WRITE(6,*) 'continuing with cmass=1.5GeV'
c         PAUSE
         CMASS=1.5D0
      ENDIF
      IF (BMASS .LE. 0D0) THEN 
         WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
         WRITE(6,*) 'continuing with bmass=5.0GeV'
c         PAUSE
         BMASS=5D0
      ENDIF
c--- 3-flavour running only
      if     (cmass .gt. 999d0) then
        T=2D0*DLOG(Q/ZMASS)
        CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF3)
	ALPHAS=AS_OUT
	RETURN
c--- 4-flavour running only
      elseif (bmass .gt. 999d0) then
        T=2D0*DLOG(Q/ZMASS)
        CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF4)
	ALPHAS=AS_OUT
	RETURN
      endif
c--- establish value of coupling at b- and c-mass and save
      IF ((AMZ .NE. AMZ0) .OR. (NLOOP .NE. NLOOP0)) THEN
         AMZ0=AMZ
         NLOOP0=NLOOP
         T=2D0*DLOG(BMASS/ZMASS)
         CALL NEWTON1(T,AMZ,AMB,NLOOP,NF5)
         T=2D0*DLOG(CMASS/BMASS)
         CALL NEWTON1(T,AMB,AMC,NLOOP,NF4)
      ENDIF
c--- evaluate strong coupling at scale q
      IF (Q  .LT. BMASS) THEN
           IF (Q  .LT. CMASS) THEN
             T=2D0*DLOG(Q/CMASS)
             CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3)
           ELSE
             T=2D0*DLOG(Q/BMASS)
             CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4)
           ENDIF
      ELSE
      T=2D0*DLOG(Q/ZMASS)
      CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF5)
      ENDIF
      ALPHAS=AS_OUT
      RETURN
      END


      SUBROUTINE DIFF(Q,AMZ,NLOOP)
      IMPLICIT NONE
      DOUBLE PRECISION BETA(3:5),B0(3:5),C1(3:5),C2(3:5)
      INTEGER NLOOP,J
      DOUBLE PRECISION Q,QP,QM,AMZ,CMASS,BMASS,X1,X2,X3,EP,DIFF1,ALPHAS
      COMMON/QMASS/CMASS,BMASS
C---     B0=(11.-2.*F/3.)/4./PI
      DATA B0/0.716197243913527D0,0.66314559621623D0,0.61009394851893D0/
C---     C1=(102.D0-38.D0/3.D0*F)/4.D0/PI/(11.D0-2.D0/3.D0*F)
      DATA C1/.565884242104515D0,0.49019722472304D0,0.40134724779695D0/
C---     C2=(2857.D0/2.D0-5033*F/18.D0+325*F**2/54)
C---     /16.D0/PI**2/(11.D0-2.D0/3.D0*F)
      DATA C2/0.453013579178645D0,0.30879037953664D0,0.14942733137107D0/
C---     DEL=SQRT(4*C2-C1**2)

      X1=ALPHAS(Q,AMZ,1)
      X2=ALPHAS(Q,AMZ,2)
      X3=ALPHAS(Q,AMZ,3)
      J=3
      IF (Q .GT. CMASS) J=4
      IF (Q .GT. BMASS) J=5
      EP=.001D0
      QP=Q*(1D0+EP)
      QM=Q*(1D0-EP)
      IF (NLOOP .EQ.1) THEN 
      BETA(J)=-B0(J)*X1**2
      DIFF1=(ALPHAS(QP,AMZ,1)-ALPHAS(QM,AMZ,1))/4d0/EP/BETA(J)
      ENDIF
      IF (NLOOP .EQ.2) THEN 
      BETA(J)=-B0(J)*X2**2*(1D0+C1(J)*X2)
      DIFF1=(ALPHAS(QP,AMZ,2)-ALPHAS(QM,AMZ,2))/4d0/EP/BETA(J)
      ENDIF
      IF (NLOOP .EQ.3) THEN 
      BETA(J)=-B0(J)*X3**2*(1D0+C1(J)*X3+C2(J)*X3**2)
      DIFF1=(ALPHAS(QP,AMZ,3)-ALPHAS(QM,AMZ,3))/4d0/EP/BETA(J)
      ENDIF
      WRITE(6,*) Q,DIFF1,NLOOP
      RETURN
      END

C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      DOUBLE PRECISION B0(3:5),C1(3:5),C2(3:5),DEL(3:5)
C      PARAMETER(PI=3.1415926535898D0)
C      NLOOP=2
C      AMZ=0.113D0
C      DO N=3,5
C      F=DFLOAT(N)
C      B0(N)=(11.D0-2.D0*F/3.D0)/4.D0/PI
C      C1(N)=(102.D0-38.D0/3.D0*F)/4.D0/PI/(11.D0-2.D0/3.D0*F)
C      C2(N)=(2857.D0/2.D0-5033*F/18.D0+325D0*F**2/54D0)
C     &   /16D0/PI**2/(11.D0-2D0/3D0*F)
C      DEL(N)=SQRT(4D0*C2(N)-C1(N)**2)
C      ENDDO
C      OPEN(UNIT=67,FILE='TEMP.DAT')
C      WRITE(67,*) B0
C      WRITE(67,*) C1
C      WRITE(67,*) C2
C      WRITE(67,*) DEL
C      DO N=1,100
C      Q=DFLOAT(N)+0.1
C      WRITE(6,*)
C      CALL DIFF(Q,AMZ,1)
C      CALL DIFF(Q,AMZ,2)
C      CALL DIFF(Q,AMZ,3)
C      ENDDO
C      STOP
C      END


