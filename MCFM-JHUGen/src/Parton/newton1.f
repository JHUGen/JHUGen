      SUBROUTINE NEWTON1(T,A_IN,A_OUT,NLOOP,NF)
C     Author: R.K. Ellis

c---  calculate a_out using nloop beta-function evolution 
c---  with nf flavours, given starting value as-in
c---  given as_in and logarithmic separation between 
c---  input scale and output scale t.
c---  Evolution is performed using Newton's method,
c---  with a precision given by tol.

      IMPLICIT NONE
      include 'types.f'
      INTEGER NLOOP,NF
      real(dp):: T,A_IN,A_OUT,AS,TOL,F2,F3,F,FP,DELTA
      real(dp):: B0(0:6),C1(0:6),C2(0:6),DEL(0:6)
      PARAMETER(TOL=5e-4_dp)
C---     B0=(11.-2.*NF/3.)/4./PI
      DATA B0/
     &0.8753521870054244_dp,0.822300539308126_dp,0.7692488916108274_dp,
     &0.716197243913529_dp,0.6631455962162306_dp,0.6100939485189321_dp,
     &0.5570423008216338_dp/
C---   C1=(102._dp-38._dp/3._dp*NF)/4._dp/PI/(11._dp-2._dp/3._dp*NF)
      DATA C1/
     &0.7379001906987874_dp,0.6879600765907734_dp,0.631131670881654_dp,
     &0.5658842421045168_dp,0.4901972247230377_dp,0.4013472477969535_dp,
     &0.2955734657420913_dp/
C---   C2=(2857._dp/2._dp-5033*NF/18._dp+325*NF**2/54)
C---   /16._dp/PI**2/(11._dp-2._dp/3._dp*NF)
      DATA C2/
     &0.8223710842788609_dp,0.7077616059424726_dp,0.5852293127502415_dp,
     &0.4530135791786467_dp,0.3087903795366415_dp,0.1494273313710745_dp,
     &-0.02940123632478559_dp/
C---     DEL=SQRT(4*C2-C1**2)  (DEL(6) imaginary, set equal to zero
      DATA DEL/
     & 1.656800424215946_dp,1.535499057891964_dp,1.393768296744871_dp,
     & 1.221404659092302_dp,0.9974307991136014_dp,0.660779624511916_dp,
     & 0._dp/

      DATA F,FP/0._dp,1._dp/
      F2(AS)=1._dp/AS+C1(NF)*LOG((C1(NF)*AS)/(1._dp+C1(NF)*AS))
      F3(AS)=1._dp/AS+0.5D0*C1(NF)
     & *LOG((C2(NF)*AS**2)/(1._dp+C1(NF)*AS+C2(NF)*AS**2))
     & -(C1(NF)**2-2D0*C2(NF))/DEL(NF)
     & *ATAN((2D0*C2(NF)*AS+C1(NF))/DEL(NF))

      IF ((NF .lt. 0) .or. (NF .gt. 6) 
     & .or. ((NF.eq.6) .and. (NLOOP.gt.2))
     & .or. (NLOOP.lt.1)) then
          write(6,*) 'unimplemented value of NF/NLOOP in newton1.f'
          write(6,*) 'NF,NLOOP=',NF,NLOOP
          STOP 
      ENDIF           

      A_OUT=A_IN/(1._dp+A_IN*B0(NF)*T)
      IF (NLOOP .EQ. 1) RETURN
      A_OUT=A_IN/(1._dp+B0(NF)*A_IN*T
     &     +C1(NF)*A_IN*LOG(1._dp+A_IN*B0(NF)*T))
      IF (A_OUT .LT. 0D0) AS=0.3D0
 30   AS=A_OUT

      IF (NLOOP .EQ. 2) THEN
          F=B0(NF)*T+F2(A_IN)-F2(AS)
          FP=1._dp/(AS**2*(1._dp+C1(NF)*AS))
      ELSEIF (NLOOP .EQ. 3) THEN
          F=B0(NF)*T+F3(A_IN)-F3(AS)
          FP=1._dp/(AS**2*(1._dp+C1(NF)*AS+C2(NF)*AS**2))
      ELSE
          WRITE(6,*) 'Unimplemented value of NLOOP in newton1'
          stop
      ENDIF
      A_OUT=AS-F/FP
      DELTA=ABS(F/FP/AS)
      IF (DELTA .GT. TOL) GO TO 30
      RETURN
      END

