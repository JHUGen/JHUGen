C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     ProjP(2,1)*ProjP(4,3)
C     
      SUBROUTINE FFFF12_1(F2, F3, F4, COUP, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      COMPLEX*16 F1(6)
      COMPLEX*16 F2(*)
      COMPLEX*16 F3(*)
      COMPLEX*16 F4(*)
      REAL*8 M1
      REAL*8 P1(0:3)
      COMPLEX*16 TMP11
      REAL*8 W1
      COMPLEX*16 DENOM
      F1(1) = +F2(1)+F3(1)+F4(1)
      F1(2) = +F2(2)+F3(2)+F4(2)
      P1(0) = -DBLE(F1(1))
      P1(1) = -DBLE(F1(2))
      P1(2) = -DIMAG(F1(2))
      P1(3) = -DIMAG(F1(1))
      TMP11 = (F4(5)*F3(5)+F4(6)*F3(6))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      F1(3)= DENOM*(-CI )* TMP11*(F2(5)*(P1(0)+P1(3))+F2(6)*(P1(1)+CI
     $ *(P1(2))))
      F1(4)= DENOM*CI * TMP11*(F2(5)*(-P1(1)+CI*(P1(2)))+F2(6)*(-P1(0)
     $ +P1(3)))
      F1(5)= DENOM*CI * TMP11*F2(5)*M1
      F1(6)= DENOM*CI * TMP11*F2(6)*M1
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     ProjP(2,1)*ProjP(4,3)
C     
      SUBROUTINE FFFF12_2_1(F2, F3, F4, COUP1, COUP2, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP1
      COMPLEX*16 COUP2
      COMPLEX*16 F1(6)
      COMPLEX*16 F2(*)
      COMPLEX*16 F3(*)
      COMPLEX*16 F4(*)
      COMPLEX*16 FTMP(6)
      REAL*8 M1
      REAL*8 P1(0:3)
      REAL*8 W1
      COMPLEX*16 DENOM
      INTEGER*4 I
      CALL FFFF12_1(F2,F3,F4,COUP1,M1,W1,F1)
      CALL FFFF2_1(F2,F3,F4,COUP2,M1,W1,FTMP)
      DO I = 3, 6
        F1(I) = F1(I) + FTMP(I)
      ENDDO
      END


