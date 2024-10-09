C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     ProjP(2,1)*ProjP(4,3)
C     
      SUBROUTINE FFFF12_3(F1, F2, F4, COUP, M3, W3,F3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 F3(6)
      COMPLEX*16 F4(*)
      REAL*8 M3
      REAL*8 P3(0:3)
      COMPLEX*16 TMP4
      REAL*8 W3
      COMPLEX*16 DENOM
      F3(1) = +F1(1)+F2(1)+F4(1)
      F3(2) = +F1(2)+F2(2)+F4(2)
      P3(0) = -DBLE(F3(1))
      P3(1) = -DBLE(F3(2))
      P3(2) = -DIMAG(F3(2))
      P3(3) = -DIMAG(F3(1))
      TMP4 = (F1(5)*F2(5)+F1(6)*F2(6))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      F3(3)= DENOM*(-CI )* TMP4*(F4(5)*(P3(0)+P3(3))+F4(6)*(P3(1)+CI
     $ *(P3(2))))
      F3(4)= DENOM*CI * TMP4*(F4(5)*(-P3(1)+CI*(P3(2)))+F4(6)*(-P3(0)
     $ +P3(3)))
      F3(5)= DENOM*CI * F4(5)*TMP4*M3
      F3(6)= DENOM*CI * F4(6)*TMP4*M3
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     ProjP(2,1)*ProjP(4,3)
C     
      SUBROUTINE FFFF12_2_3(F1, F2, F4, COUP1, COUP2, M3, W3,F3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP1
      COMPLEX*16 COUP2
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 F3(6)
      COMPLEX*16 F4(*)
      COMPLEX*16 FTMP(6)
      REAL*8 M3
      REAL*8 P3(0:3)
      REAL*8 W3
      COMPLEX*16 DENOM
      INTEGER*4 I
      CALL FFFF12_3(F1,F2,F4,COUP1,M3,W3,F3)
      CALL FFFF2_3(F1,F2,F4,COUP2,M3,W3,FTMP)
      DO I = 3, 6
        F3(I) = F3(I) + FTMP(I)
      ENDDO
      END

