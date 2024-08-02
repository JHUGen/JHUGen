C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(-1,1)*P(-1,2) + P(-1,1)*P(-1,3) + P(-1,2)*P(-1,3)
C     
      SUBROUTINE SSS2_0(S1, S2, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      REAL*8 P3(0:3)
      COMPLEX*16 S1(*)
      COMPLEX*16 S2(*)
      COMPLEX*16 S3(*)
      COMPLEX*16 TMP18
      COMPLEX*16 TMP29
      COMPLEX*16 TMP30
      COMPLEX*16 VERTEX
      P1(0) = DBLE(S1(1))
      P1(1) = DBLE(S1(2))
      P1(2) = DIMAG(S1(2))
      P1(3) = DIMAG(S1(1))
      P2(0) = DBLE(S2(1))
      P2(1) = DBLE(S2(2))
      P2(2) = DIMAG(S2(2))
      P2(3) = DIMAG(S2(1))
      P3(0) = DBLE(S3(1))
      P3(1) = DBLE(S3(2))
      P3(2) = DIMAG(S3(2))
      P3(3) = DIMAG(S3(1))
      TMP18 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP29 = (P1(0)*P3(0)-P1(1)*P3(1)-P1(2)*P3(2)-P1(3)*P3(3))
      TMP30 = (P2(0)*P3(0)-P2(1)*P3(1)-P2(2)*P3(2)-P2(3)*P3(3))
      VERTEX = COUP*(-S1(3)*S2(3)*S3(3)*(+CI*(TMP18+TMP29+TMP30)))
      END


