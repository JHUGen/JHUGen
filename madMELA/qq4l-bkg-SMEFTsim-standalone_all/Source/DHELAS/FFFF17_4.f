C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-2,1)*ProjP(-3,3)
C     
      SUBROUTINE FFFF17_4(F1, F2, F3, COUP, M4, W4,F4)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 F3(*)
      COMPLEX*16 F4(6)
      REAL*8 M4
      REAL*8 P4(0:3)
      REAL*8 W4
      COMPLEX*16 DENOM
      F4(1) = +F1(1)+F2(1)+F3(1)
      F4(2) = +F1(2)+F2(2)+F3(2)
      P4(0) = -DBLE(F4(1))
      P4(1) = -DBLE(F4(2))
      P4(2) = -DIMAG(F4(2))
      P4(3) = -DIMAG(F4(1))
      DENOM = COUP/(P4(0)**2-P4(1)**2-P4(2)**2-P4(3)**2 - M4 * (M4 -CI
     $ * W4))
      F4(3)= DENOM*2D0 * CI * F1(3)*M4*(F2(5)*F3(5)+F2(6)*F3(6))
      F4(4)= DENOM*2D0 * CI * F1(4)*M4*(F2(5)*F3(5)+F2(6)*F3(6))
      F4(5)= DENOM*2D0 * CI*(F1(3)*(F2(5)*F3(5)*(P4(0)+P4(3))+F2(6)
     $ *F3(6)*(P4(0)+P4(3)))+F1(4)*(F2(5)*F3(5)*(P4(1)-CI*(P4(2)))
     $ +F2(6)*F3(6)*(P4(1)-CI*(P4(2)))))
      F4(6)= DENOM*2D0 * CI*(F1(3)*(F2(5)*F3(5)*(P4(1)+CI*(P4(2)))
     $ +F2(6)*F3(6)*(P4(1)+CI*(P4(2))))+F1(4)*(F2(5)*F3(5)*(P4(0)-P4(3)
     $ )+F2(6)*F3(6)*(P4(0)-P4(3))))
      END


