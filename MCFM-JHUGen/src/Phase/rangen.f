C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RANGEN(N,R)
      IMPLICIT NONE
C---RANDOM NUMBER GENERATOR
C   USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C   RETURNS A VECTOR OF N RANDOM VALUES
C   IF (N.EQ.0) THE FIRST TWO VALUES IN R SET THE SEEDS
C   IF (N.LT.0) PRINT THE CURRENT VALUES OF THE SEEDS
      include 'types.f'
      real(dp):: R(*)
      integer:: N,I,ISEED(2),K,IZ
      DATA ISEED/12345,678900/
      IF (N.LT.0) WRITE (*,'(I10,A,I10,I11)') -N-1,', ISEED=',ISEED
      IF (N.GT.0) THEN
        DO I=1,N
          K=ISEED(1)/53668
          ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
          IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
          K=ISEED(2)/52774
          ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
          IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
          IZ=ISEED(1)-ISEED(2)
          IF (IZ.LT.1) IZ=IZ+2147483562
          R(I)=real(IZ,dp)*4.656613E-10_dp
        ENDDO
      ELSEIF (N.EQ.0) THEN
        ISEED(1)=NINT(R(1))
        ISEED(2)=NINT(R(2))
      ENDIF
      END
C-----------------------------------------------------------------------
