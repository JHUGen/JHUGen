      FUNCTION ALPHAS(Q,AMZ,NLOOPS)
c     Evaluation of strong coupling constant alpha_S
c     Author: R.K. Ellis

c     q -- scale at which alpha_s is to be evaluated
c     amz -- value of alpha_s at the mass of the Z-boson
c     nloop -- the number of loops (1,2, or 3) at which beta 
c     function is evaluated to determine running.
c     the values of the cmass and the bmass should be set
c     in common block qmass.
      IMPLICIT NONE
      include 'types.f'
      include 'constants.f'
c--- added by JMC: use consistent value of MZ from common block
      include 'masses.f'
      real(dp)::ALPHAS
      real(dp)::Q,T,AMZ,AMZ0,AMB,AMC,BMASS,CMASS,AS_OUT
      INTEGER::NLOOPS,NLOOP0,NF3,NF4,NF5
      PARAMETER(NF5=5,NF4=4,NF3=3)
      COMMON/QMASS/CMASS,BMASS
      SAVE AMZ0,NLOOP0,AMB,AMC
      DATA AMZ0,NLOOP0/0._dp,0/

          

      IF (Q .LE. 0._dp) THEN 
   !      WRITE(6,*) 'q .le. 0 in alphas'
   !      WRITE(6,*) 'q= ',Q
	 stop
c         PAUSE
      ENDIF
      IF (AMZ .LE. 0._dp) THEN 
   !      WRITE(6,*) 'amz .le. 0 in alphas',AMZ
   !      WRITE(6,*) 'continuing with amz=0.117'
c         PAUSE
         AMZ=0.117_dp
      ENDIF
      IF (CMASS .LE. 0.3_dp) THEN 
   !      WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
   !      WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
   !      WRITE(6,*) 'continuing with cmass=1.5GeV'
c         PAUSE
         CMASS=1.5_dp
      ENDIF
      IF (BMASS .LE. 0._dp) THEN 
   !      WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
   !      WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
   !      WRITE(6,*) 'continuing with bmass=5.0GeV'
c         PAUSE
         BMASS=5._dp
      ENDIF
c--- 3-flavour running only
      if     (cmass .gt. 999._dp) then
        T=2._dp*log(Q/ZMASS)
        CALL NEWTON1(T,AMZ,AS_OUT,NLOOPS,NF3)
	ALPHAS=AS_OUT
	RETURN
c--- 4-flavour running only
      elseif (bmass .gt. 999._dp) then
        T=2._dp*log(Q/ZMASS)
        CALL NEWTON1(T,AMZ,AS_OUT,NLOOPS,NF4)
	ALPHAS=AS_OUT
	RETURN
      endif
c--- establish value of coupling at b- and c-mass and save
      IF ((AMZ .NE. AMZ0) .OR. (NLOOPS .NE. NLOOP0)) THEN
         AMZ0=AMZ
         NLOOP0=NLOOPS
         T=2._dp*log(BMASS/ZMASS)
         CALL NEWTON1(T,AMZ,AMB,NLOOPS,NF5)
         T=2._dp*log(CMASS/BMASS)
         CALL NEWTON1(T,AMB,AMC,NLOOPS,NF4)
      ENDIF
c--- evaluate strong coupling at scale q
      IF (Q  .LT. BMASS) THEN
           IF (Q  .LT. CMASS) THEN
             T=2._dp*log(Q/CMASS)
             CALL NEWTON1(T,AMC,AS_OUT,NLOOPS,NF3)
           ELSE
             T=2._dp*log(Q/BMASS)
             CALL NEWTON1(T,AMB,AS_OUT,NLOOPS,NF4)
           ENDIF
      ELSE
      T=2._dp*log(Q/ZMASS)
      CALL NEWTON1(T,AMZ,AS_OUT,NLOOPS,NF5)
      ENDIF
      ALPHAS=AS_OUT
      RETURN
      END


C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      real(dp):: B0(3:5),C1(3:5),C2(3:5),DEL(3:5)
C      PARAMETER(PI=3.1415926535898D0)
C      NLOOP=2
C      AMZ=0.113D0
C      DO N=3,5
C      F=real(N,dp)
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
C      Q=real(N,dp)+0.1
C      WRITE(6,*)
C      CALL DIFF(Q,AMZ,1)
C      CALL DIFF(Q,AMZ,2)
C      CALL DIFF(Q,AMZ,3)
C      ENDDO
C      STOP
C      END


