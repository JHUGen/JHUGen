
C     PY ((2, -2), (-11, -11, 11, 11)) : (2, -2, -11, 11, -11, 11) #
C      M0_ 1
C     PY ((2, -2), (-13, -13, 13, 13)) : (2, -2, -13, 13, -13, 13) #
C      M1_ 1
C     PY ((2, -2), (-15, -15, 15, 15)) : (2, -2, -15, 15, -15, 15) #
C      M2_ 1
C     PY ((1, -1), (-11, -11, 11, 11)) : (1, -1, -11, 11, -11, 11) #
C      M3_ 1
C     PY ((1, -1), (-13, -13, 13, 13)) : (1, -1, -13, 13, -13, 13) #
C      M4_ 1
C     PY ((1, -1), (-15, -15, 15, 15)) : (1, -1, -15, 15, -15, 15) #
C      M5_ 1
C     PY ((3, -3), (-11, -11, 11, 11)) : (3, -3, -11, 11, -11, 11) #
C      M6_ 1
C     PY ((3, -3), (-13, -13, 13, 13)) : (3, -3, -13, 13, -13, 13) #
C      M7_ 1
C     PY ((3, -3), (-15, -15, 15, 15)) : (3, -3, -15, 15, -15, 15) #
C      M8_ 1
C     PY ((4, -4), (-11, -11, 11, 11)) : (4, -4, -11, 11, -11, 11) #
C      M9_ 1
C     PY ((4, -4), (-13, -13, 13, 13)) : (4, -4, -13, 13, -13, 13) #
C      M10_ 1
C     PY ((4, -4), (-15, -15, 15, 15)) : (4, -4, -15, 15, -15, 15) #
C      M11_ 1
C     PY ((5, -5), (-11, -11, 11, 11)) : (5, -5, -11, 11, -11, 11) #
C      M12_ 1
C     PY ((5, -5), (-13, -13, 13, 13)) : (5, -5, -13, 13, -13, 13) #
C      M13_ 1
C     PY ((5, -5), (-15, -15, 15, 15)) : (5, -5, -15, 15, -15, 15) #
C      M14_ 1
C     PY ((2, -2), (-13, -11, 11, 13)) : (2, -2, -11, 11, -13, 13) #
C      M15_ 1
C     PY ((2, -2), (-15, -11, 11, 15)) : (2, -2, -11, 11, -15, 15) #
C      M16_ 1
C     PY ((2, -2), (-15, -13, 13, 15)) : (2, -2, -13, 13, -15, 15) #
C      M17_ 1
C     PY ((1, -1), (-13, -11, 11, 13)) : (1, -1, -11, 11, -13, 13) #
C      M18_ 1
C     PY ((1, -1), (-15, -11, 11, 15)) : (1, -1, -11, 11, -15, 15) #
C      M19_ 1
C     PY ((1, -1), (-15, -13, 13, 15)) : (1, -1, -13, 13, -15, 15) #
C      M20_ 1
C     PY ((3, -3), (-13, -11, 11, 13)) : (3, -3, -11, 11, -13, 13) #
C      M21_ 1
C     PY ((3, -3), (-15, -11, 11, 15)) : (3, -3, -11, 11, -15, 15) #
C      M22_ 1
C     PY ((3, -3), (-15, -13, 13, 15)) : (3, -3, -13, 13, -15, 15) #
C      M23_ 1
C     PY ((4, -4), (-13, -11, 11, 13)) : (4, -4, -11, 11, -13, 13) #
C      M24_ 1
C     PY ((4, -4), (-15, -11, 11, 15)) : (4, -4, -11, 11, -15, 15) #
C      M25_ 1
C     PY ((4, -4), (-15, -13, 13, 15)) : (4, -4, -13, 13, -15, 15) #
C      M26_ 1
C     PY ((5, -5), (-13, -11, 11, 13)) : (5, -5, -11, 11, -13, 13) #
C      M27_ 1
C     PY ((5, -5), (-15, -11, 11, 15)) : (5, -5, -11, 11, -15, 15) #
C      M28_ 1
C     PY ((5, -5), (-15, -13, 13, 15)) : (5, -5, -13, 13, -15, 15) #
C      M29_ 1
      SUBROUTINE SMATRIXHEL(PDGS, PROCID, NPDG, P, ALPHAS, SCALE2,
     $  NHEL, ANS)
      IMPLICIT NONE
C     ALPHAS is given at scale2 (SHOULD be different of 0 for loop
C      induced, ignore for LO)  

CF2PY double precision, intent(in), dimension(0:3,npdg) :: p
CF2PY integer, intent(in), dimension(npdg) :: pdgs
CF2PY integer, intent(in):: procid
CF2PY integer, intent(in) :: npdg
CF2PY double precision, intent(out) :: ANS
CF2PY double precision, intent(in) :: ALPHAS
CF2PY double precision, intent(in) :: SCALE2
      INTEGER PDGS(*)
      INTEGER NPDG, NHEL, PROCID
      DOUBLE PRECISION P(*)
      DOUBLE PRECISION ANS, ALPHAS, PI,SCALE2
      INCLUDE 'coupl.inc'

      PI = 3.141592653589793D0
      G = 2* DSQRT(ALPHAS*PI)
      CALL UPDATE_AS_PARAM()
C     if (scale2.ne.0d0) stop 1

      IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-11.EQ.PDGS(5).AND.11.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 0
        CALL M0_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 1
        CALL M1_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.-15.EQ.PDGS(3)
     $ .AND.15.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 2
        CALL M2_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-11.EQ.PDGS(5).AND.11.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 3
        CALL M3_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 4
        CALL M4_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.-15.EQ.PDGS(3)
     $ .AND.15.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 5
        CALL M5_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-11.EQ.PDGS(5).AND.11.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 6
        CALL M6_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 7
        CALL M7_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.-15.EQ.PDGS(3)
     $ .AND.15.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 8
        CALL M8_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-11.EQ.PDGS(5).AND.11.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 9
        CALL M9_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 10
        CALL M10_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.-15.EQ.PDGS(3)
     $ .AND.15.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 11
        CALL M11_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-11.EQ.PDGS(5).AND.11.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 12
        CALL M12_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 13
        CALL M13_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.-15.EQ.PDGS(3)
     $ .AND.15.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 14
        CALL M14_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 15
        CALL M15_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 16
        CALL M16_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 17
        CALL M17_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 18
        CALL M18_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 19
        CALL M19_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 20
        CALL M20_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 21
        CALL M21_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 22
        CALL M22_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 23
        CALL M23_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 24
        CALL M24_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 25
        CALL M25_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 26
        CALL M26_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-13.EQ.PDGS(5).AND.13.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 27
        CALL M27_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.-11.EQ.PDGS(3)
     $ .AND.11.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 28
        CALL M28_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.-13.EQ.PDGS(3)
     $ .AND.13.EQ.PDGS(4).AND.-15.EQ.PDGS(5).AND.15.EQ.PDGS(6)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 29
        CALL M29_SMATRIXHEL(P, NHEL, ANS)
      ENDIF

      RETURN
      END

      SUBROUTINE INITIALISE(PATH)
C     ROUTINE FOR F2PY to read the benchmark point.
      IMPLICIT NONE
      CHARACTER*512 PATH
CF2PY INTENT(IN) :: PATH
      CALL SETPARA(PATH)  !first call to setup the paramaters
      RETURN
      END


      SUBROUTINE CHANGE_PARA(NAME, VALUE)
      IMPLICIT NONE
CF2PY intent(in) :: name
CF2PY intent(in) :: value

      CHARACTER*512 NAME
      DOUBLE PRECISION VALUE

      LOGICAL M15_HELRESET
      COMMON /M15_HELRESET/ M15_HELRESET
      LOGICAL M13_HELRESET
      COMMON /M13_HELRESET/ M13_HELRESET
      LOGICAL M2_HELRESET
      COMMON /M2_HELRESET/ M2_HELRESET
      LOGICAL M12_HELRESET
      COMMON /M12_HELRESET/ M12_HELRESET
      LOGICAL M11_HELRESET
      COMMON /M11_HELRESET/ M11_HELRESET
      LOGICAL M25_HELRESET
      COMMON /M25_HELRESET/ M25_HELRESET
      LOGICAL M16_HELRESET
      COMMON /M16_HELRESET/ M16_HELRESET
      LOGICAL M5_HELRESET
      COMMON /M5_HELRESET/ M5_HELRESET
      LOGICAL M18_HELRESET
      COMMON /M18_HELRESET/ M18_HELRESET
      LOGICAL M14_HELRESET
      COMMON /M14_HELRESET/ M14_HELRESET
      LOGICAL M6_HELRESET
      COMMON /M6_HELRESET/ M6_HELRESET
      LOGICAL M7_HELRESET
      COMMON /M7_HELRESET/ M7_HELRESET
      LOGICAL M19_HELRESET
      COMMON /M19_HELRESET/ M19_HELRESET
      LOGICAL M4_HELRESET
      COMMON /M4_HELRESET/ M4_HELRESET
      LOGICAL M3_HELRESET
      COMMON /M3_HELRESET/ M3_HELRESET
      LOGICAL M21_HELRESET
      COMMON /M21_HELRESET/ M21_HELRESET
      LOGICAL M20_HELRESET
      COMMON /M20_HELRESET/ M20_HELRESET
      LOGICAL M0_HELRESET
      COMMON /M0_HELRESET/ M0_HELRESET
      LOGICAL M28_HELRESET
      COMMON /M28_HELRESET/ M28_HELRESET
      LOGICAL M26_HELRESET
      COMMON /M26_HELRESET/ M26_HELRESET
      LOGICAL M29_HELRESET
      COMMON /M29_HELRESET/ M29_HELRESET
      LOGICAL M27_HELRESET
      COMMON /M27_HELRESET/ M27_HELRESET
      LOGICAL M22_HELRESET
      COMMON /M22_HELRESET/ M22_HELRESET
      LOGICAL M10_HELRESET
      COMMON /M10_HELRESET/ M10_HELRESET
      LOGICAL M9_HELRESET
      COMMON /M9_HELRESET/ M9_HELRESET
      LOGICAL M8_HELRESET
      COMMON /M8_HELRESET/ M8_HELRESET
      LOGICAL M23_HELRESET
      COMMON /M23_HELRESET/ M23_HELRESET
      LOGICAL M1_HELRESET
      COMMON /M1_HELRESET/ M1_HELRESET
      LOGICAL M17_HELRESET
      COMMON /M17_HELRESET/ M17_HELRESET
      LOGICAL M24_HELRESET
      COMMON /M24_HELRESET/ M24_HELRESET

      INCLUDE '../Source/MODEL/input.inc'
      INCLUDE '../Source/MODEL/coupl.inc'

      M15_HELRESET = .TRUE.
      M13_HELRESET = .TRUE.
      M2_HELRESET = .TRUE.
      M12_HELRESET = .TRUE.
      M11_HELRESET = .TRUE.
      M25_HELRESET = .TRUE.
      M16_HELRESET = .TRUE.
      M5_HELRESET = .TRUE.
      M18_HELRESET = .TRUE.
      M14_HELRESET = .TRUE.
      M6_HELRESET = .TRUE.
      M7_HELRESET = .TRUE.
      M19_HELRESET = .TRUE.
      M4_HELRESET = .TRUE.
      M3_HELRESET = .TRUE.
      M21_HELRESET = .TRUE.
      M20_HELRESET = .TRUE.
      M0_HELRESET = .TRUE.
      M28_HELRESET = .TRUE.
      M26_HELRESET = .TRUE.
      M29_HELRESET = .TRUE.
      M27_HELRESET = .TRUE.
      M22_HELRESET = .TRUE.
      M10_HELRESET = .TRUE.
      M9_HELRESET = .TRUE.
      M8_HELRESET = .TRUE.
      M23_HELRESET = .TRUE.
      M1_HELRESET = .TRUE.
      M17_HELRESET = .TRUE.
      M24_HELRESET = .TRUE.

      SELECT CASE (NAME)
      CASE ('CKMlambda')
      MDL_CKMLAMBDA = VALUE
      CASE ('CKMBLOCK_2')
      MDL_CKMLAMBDA = VALUE
      CASE ('CKMA')
      MDL_CKMA = VALUE
      CASE ('CKMBLOCK_3')
      MDL_CKMA = VALUE
      CASE ('CKMrho')
      MDL_CKMRHO = VALUE
      CASE ('CKMBLOCK_4')
      MDL_CKMRHO = VALUE
      CASE ('CKMeta')
      MDL_CKMETA = VALUE
      CASE ('CKMBLOCK_5')
      MDL_CKMETA = VALUE
      CASE ('MD')
      MDL_MD = VALUE
      CASE ('MASS_1')
      MDL_MD = VALUE
      CASE ('MU')
      MDL_MU = VALUE
      CASE ('MASS_2')
      MDL_MU = VALUE
      CASE ('MS')
      MDL_MS = VALUE
      CASE ('MASS_3')
      MDL_MS = VALUE
      CASE ('MC')
      MDL_MC = VALUE
      CASE ('MASS_4')
      MDL_MC = VALUE
      CASE ('MB')
      MDL_MB = VALUE
      CASE ('MASS_5')
      MDL_MB = VALUE
      CASE ('MT')
      MDL_MT = VALUE
      CASE ('MASS_6')
      MDL_MT = VALUE
      CASE ('Me')
      MDL_ME = VALUE
      CASE ('MASS_11')
      MDL_ME = VALUE
      CASE ('MMU')
      MDL_MMU = VALUE
      CASE ('MASS_13')
      MDL_MMU = VALUE
      CASE ('MTA')
      MDL_MTA = VALUE
      CASE ('MASS_15')
      MDL_MTA = VALUE
      CASE ('MZ')
      MDL_MZ = VALUE
      CASE ('MASS_23')
      MDL_MZ = VALUE
      CASE ('MH')
      MDL_MH = VALUE
      CASE ('MASS_25')
      MDL_MH = VALUE
      CASE ('cG')
      MDL_CG = VALUE
      CASE ('SMEFT_1')
      MDL_CG = VALUE
      CASE ('cW')
      MDL_CW = VALUE
      CASE ('SMEFT_2')
      MDL_CW = VALUE
      CASE ('cH')
      MDL_CH = VALUE
      CASE ('SMEFT_3')
      MDL_CH = VALUE
      CASE ('cHbox')
      MDL_CHBOX = VALUE
      CASE ('SMEFT_4')
      MDL_CHBOX = VALUE
      CASE ('cHDD')
      MDL_CHDD = VALUE
      CASE ('SMEFT_5')
      MDL_CHDD = VALUE
      CASE ('cHG')
      MDL_CHG = VALUE
      CASE ('SMEFT_6')
      MDL_CHG = VALUE
      CASE ('cHW')
      MDL_CHW = VALUE
      CASE ('SMEFT_7')
      MDL_CHW = VALUE
      CASE ('cHB')
      MDL_CHB = VALUE
      CASE ('SMEFT_8')
      MDL_CHB = VALUE
      CASE ('cHWB')
      MDL_CHWB = VALUE
      CASE ('SMEFT_9')
      MDL_CHWB = VALUE
      CASE ('ceHRe')
      MDL_CEHRE = VALUE
      CASE ('SMEFT_10')
      MDL_CEHRE = VALUE
      CASE ('cuHRe')
      MDL_CUHRE = VALUE
      CASE ('SMEFT_11')
      MDL_CUHRE = VALUE
      CASE ('cdHRe')
      MDL_CDHRE = VALUE
      CASE ('SMEFT_12')
      MDL_CDHRE = VALUE
      CASE ('ceWRe')
      MDL_CEWRE = VALUE
      CASE ('SMEFT_13')
      MDL_CEWRE = VALUE
      CASE ('ceBRe')
      MDL_CEBRE = VALUE
      CASE ('SMEFT_14')
      MDL_CEBRE = VALUE
      CASE ('cuGRe')
      MDL_CUGRE = VALUE
      CASE ('SMEFT_15')
      MDL_CUGRE = VALUE
      CASE ('cuWRe')
      MDL_CUWRE = VALUE
      CASE ('SMEFT_16')
      MDL_CUWRE = VALUE
      CASE ('cuBRe')
      MDL_CUBRE = VALUE
      CASE ('SMEFT_17')
      MDL_CUBRE = VALUE
      CASE ('cdGRe')
      MDL_CDGRE = VALUE
      CASE ('SMEFT_18')
      MDL_CDGRE = VALUE
      CASE ('cdWRe')
      MDL_CDWRE = VALUE
      CASE ('SMEFT_19')
      MDL_CDWRE = VALUE
      CASE ('cdBRe')
      MDL_CDBRE = VALUE
      CASE ('SMEFT_20')
      MDL_CDBRE = VALUE
      CASE ('cHl1')
      MDL_CHL1 = VALUE
      CASE ('SMEFT_21')
      MDL_CHL1 = VALUE
      CASE ('cHl3')
      MDL_CHL3 = VALUE
      CASE ('SMEFT_22')
      MDL_CHL3 = VALUE
      CASE ('cHe')
      MDL_CHE = VALUE
      CASE ('SMEFT_23')
      MDL_CHE = VALUE
      CASE ('cHq1')
      MDL_CHQ1 = VALUE
      CASE ('SMEFT_24')
      MDL_CHQ1 = VALUE
      CASE ('cHq3')
      MDL_CHQ3 = VALUE
      CASE ('SMEFT_25')
      MDL_CHQ3 = VALUE
      CASE ('cHu')
      MDL_CHU = VALUE
      CASE ('SMEFT_26')
      MDL_CHU = VALUE
      CASE ('cHd')
      MDL_CHD = VALUE
      CASE ('SMEFT_27')
      MDL_CHD = VALUE
      CASE ('cHudRe')
      MDL_CHUDRE = VALUE
      CASE ('SMEFT_28')
      MDL_CHUDRE = VALUE
      CASE ('cll')
      MDL_CLL = VALUE
      CASE ('SMEFT_29')
      MDL_CLL = VALUE
      CASE ('cll1')
      MDL_CLL1 = VALUE
      CASE ('SMEFT_30')
      MDL_CLL1 = VALUE
      CASE ('cqq1')
      MDL_CQQ1 = VALUE
      CASE ('SMEFT_31')
      MDL_CQQ1 = VALUE
      CASE ('cqq11')
      MDL_CQQ11 = VALUE
      CASE ('SMEFT_32')
      MDL_CQQ11 = VALUE
      CASE ('cqq3')
      MDL_CQQ3 = VALUE
      CASE ('SMEFT_33')
      MDL_CQQ3 = VALUE
      CASE ('cqq31')
      MDL_CQQ31 = VALUE
      CASE ('SMEFT_34')
      MDL_CQQ31 = VALUE
      CASE ('clq1')
      MDL_CLQ1 = VALUE
      CASE ('SMEFT_35')
      MDL_CLQ1 = VALUE
      CASE ('clq3')
      MDL_CLQ3 = VALUE
      CASE ('SMEFT_36')
      MDL_CLQ3 = VALUE
      CASE ('cee')
      MDL_CEE = VALUE
      CASE ('SMEFT_37')
      MDL_CEE = VALUE
      CASE ('cuu')
      MDL_CUU = VALUE
      CASE ('SMEFT_38')
      MDL_CUU = VALUE
      CASE ('cuu1')
      MDL_CUU1 = VALUE
      CASE ('SMEFT_39')
      MDL_CUU1 = VALUE
      CASE ('cdd')
      MDL_CDD = VALUE
      CASE ('SMEFT_40')
      MDL_CDD = VALUE
      CASE ('cdd1')
      MDL_CDD1 = VALUE
      CASE ('SMEFT_41')
      MDL_CDD1 = VALUE
      CASE ('ceu')
      MDL_CEU = VALUE
      CASE ('SMEFT_42')
      MDL_CEU = VALUE
      CASE ('ced')
      MDL_CED = VALUE
      CASE ('SMEFT_43')
      MDL_CED = VALUE
      CASE ('cud1')
      MDL_CUD1 = VALUE
      CASE ('SMEFT_44')
      MDL_CUD1 = VALUE
      CASE ('cud8')
      MDL_CUD8 = VALUE
      CASE ('SMEFT_45')
      MDL_CUD8 = VALUE
      CASE ('cle')
      MDL_CLE = VALUE
      CASE ('SMEFT_46')
      MDL_CLE = VALUE
      CASE ('clu')
      MDL_CLU = VALUE
      CASE ('SMEFT_47')
      MDL_CLU = VALUE
      CASE ('cld')
      MDL_CLD = VALUE
      CASE ('SMEFT_48')
      MDL_CLD = VALUE
      CASE ('cqe')
      MDL_CQE = VALUE
      CASE ('SMEFT_49')
      MDL_CQE = VALUE
      CASE ('cqu1')
      MDL_CQU1 = VALUE
      CASE ('SMEFT_50')
      MDL_CQU1 = VALUE
      CASE ('cqu8')
      MDL_CQU8 = VALUE
      CASE ('SMEFT_51')
      MDL_CQU8 = VALUE
      CASE ('cqd1')
      MDL_CQD1 = VALUE
      CASE ('SMEFT_52')
      MDL_CQD1 = VALUE
      CASE ('cqd8')
      MDL_CQD8 = VALUE
      CASE ('SMEFT_53')
      MDL_CQD8 = VALUE
      CASE ('cledqRe')
      MDL_CLEDQRE = VALUE
      CASE ('SMEFT_54')
      MDL_CLEDQRE = VALUE
      CASE ('cquqd1Re')
      MDL_CQUQD1RE = VALUE
      CASE ('SMEFT_55')
      MDL_CQUQD1RE = VALUE
      CASE ('cquqd11Re')
      MDL_CQUQD11RE = VALUE
      CASE ('SMEFT_56')
      MDL_CQUQD11RE = VALUE
      CASE ('cquqd8Re')
      MDL_CQUQD8RE = VALUE
      CASE ('SMEFT_57')
      MDL_CQUQD8RE = VALUE
      CASE ('cquqd81Re')
      MDL_CQUQD81RE = VALUE
      CASE ('SMEFT_58')
      MDL_CQUQD81RE = VALUE
      CASE ('clequ1Re')
      MDL_CLEQU1RE = VALUE
      CASE ('SMEFT_59')
      MDL_CLEQU1RE = VALUE
      CASE ('clequ3Re')
      MDL_CLEQU3RE = VALUE
      CASE ('SMEFT_60')
      MDL_CLEQU3RE = VALUE
      CASE ('cGtil')
      MDL_CGTIL = VALUE
      CASE ('SMEFTCPV_1')
      MDL_CGTIL = VALUE
      CASE ('cWtil')
      MDL_CWTIL = VALUE
      CASE ('SMEFTCPV_2')
      MDL_CWTIL = VALUE
      CASE ('cHGtil')
      MDL_CHGTIL = VALUE
      CASE ('SMEFTCPV_3')
      MDL_CHGTIL = VALUE
      CASE ('cHWtil')
      MDL_CHWTIL = VALUE
      CASE ('SMEFTCPV_4')
      MDL_CHWTIL = VALUE
      CASE ('cHBtil')
      MDL_CHBTIL = VALUE
      CASE ('SMEFTCPV_5')
      MDL_CHBTIL = VALUE
      CASE ('cHWBtil')
      MDL_CHWBTIL = VALUE
      CASE ('SMEFTCPV_6')
      MDL_CHWBTIL = VALUE
      CASE ('ceWIm')
      MDL_CEWIM = VALUE
      CASE ('SMEFTCPV_7')
      MDL_CEWIM = VALUE
      CASE ('ceBIm')
      MDL_CEBIM = VALUE
      CASE ('SMEFTCPV_8')
      MDL_CEBIM = VALUE
      CASE ('cuGIm')
      MDL_CUGIM = VALUE
      CASE ('SMEFTCPV_9')
      MDL_CUGIM = VALUE
      CASE ('cuWIm')
      MDL_CUWIM = VALUE
      CASE ('SMEFTCPV_10')
      MDL_CUWIM = VALUE
      CASE ('cuBIm')
      MDL_CUBIM = VALUE
      CASE ('SMEFTCPV_11')
      MDL_CUBIM = VALUE
      CASE ('cdGIm')
      MDL_CDGIM = VALUE
      CASE ('SMEFTCPV_12')
      MDL_CDGIM = VALUE
      CASE ('cdWIm')
      MDL_CDWIM = VALUE
      CASE ('SMEFTCPV_13')
      MDL_CDWIM = VALUE
      CASE ('cdBIm')
      MDL_CDBIM = VALUE
      CASE ('SMEFTCPV_14')
      MDL_CDBIM = VALUE
      CASE ('cHudIm')
      MDL_CHUDIM = VALUE
      CASE ('SMEFTCPV_15')
      MDL_CHUDIM = VALUE
      CASE ('ceHIm')
      MDL_CEHIM = VALUE
      CASE ('SMEFTCPV_16')
      MDL_CEHIM = VALUE
      CASE ('cuHIm')
      MDL_CUHIM = VALUE
      CASE ('SMEFTCPV_17')
      MDL_CUHIM = VALUE
      CASE ('cdHIm')
      MDL_CDHIM = VALUE
      CASE ('SMEFTCPV_18')
      MDL_CDHIM = VALUE
      CASE ('cledqIm')
      MDL_CLEDQIM = VALUE
      CASE ('SMEFTCPV_19')
      MDL_CLEDQIM = VALUE
      CASE ('cquqd1Im')
      MDL_CQUQD1IM = VALUE
      CASE ('SMEFTCPV_20')
      MDL_CQUQD1IM = VALUE
      CASE ('cquqd8Im')
      MDL_CQUQD8IM = VALUE
      CASE ('SMEFTCPV_21')
      MDL_CQUQD8IM = VALUE
      CASE ('cquqd11Im')
      MDL_CQUQD11IM = VALUE
      CASE ('SMEFTCPV_22')
      MDL_CQUQD11IM = VALUE
      CASE ('cquqd81Im')
      MDL_CQUQD81IM = VALUE
      CASE ('SMEFTCPV_23')
      MDL_CQUQD81IM = VALUE
      CASE ('clequ1Im')
      MDL_CLEQU1IM = VALUE
      CASE ('SMEFTCPV_24')
      MDL_CLEQU1IM = VALUE
      CASE ('clequ3Im')
      MDL_CLEQU3IM = VALUE
      CASE ('SMEFTCPV_25')
      MDL_CLEQU3IM = VALUE
      CASE ('LambdaSMEFT')
      MDL_LAMBDASMEFT = VALUE
      CASE ('SMEFTCUTOFF_1')
      MDL_LAMBDASMEFT = VALUE
      CASE ('aEW')
      MDL_AEW = VALUE
      CASE ('SMINPUTS_1')
      MDL_AEW = VALUE
      CASE ('Gf')
      MDL_GF = VALUE
      CASE ('SMINPUTS_2')
      MDL_GF = VALUE
      CASE ('aS')
      AS = VALUE
      CASE ('SMINPUTS_3')
      AS = VALUE
      CASE ('linearPropCorrections')
      MDL_LINEARPROPCORRECTIONS = VALUE
      CASE ('SWITCHES_1')
      MDL_LINEARPROPCORRECTIONS = VALUE
      CASE ('ymdo')
      MDL_YMDO = VALUE
      CASE ('YUKAWA_1')
      MDL_YMDO = VALUE
      CASE ('ymup')
      MDL_YMUP = VALUE
      CASE ('YUKAWA_2')
      MDL_YMUP = VALUE
      CASE ('yms')
      MDL_YMS = VALUE
      CASE ('YUKAWA_3')
      MDL_YMS = VALUE
      CASE ('ymc')
      MDL_YMC = VALUE
      CASE ('YUKAWA_4')
      MDL_YMC = VALUE
      CASE ('ymb')
      MDL_YMB = VALUE
      CASE ('YUKAWA_5')
      MDL_YMB = VALUE
      CASE ('ymt')
      MDL_YMT = VALUE
      CASE ('YUKAWA_6')
      MDL_YMT = VALUE
      CASE ('yme')
      MDL_YME = VALUE
      CASE ('YUKAWA_11')
      MDL_YME = VALUE
      CASE ('ymm')
      MDL_YMM = VALUE
      CASE ('YUKAWA_13')
      MDL_YMM = VALUE
      CASE ('ymtau')
      MDL_YMTAU = VALUE
      CASE ('YUKAWA_15')
      MDL_YMTAU = VALUE
      CASE ('WT')
      MDL_WT = VALUE
      CASE ('DECAY_6')
      MDL_WT = VALUE
      CASE ('WZ')
      MDL_WZ = VALUE
      CASE ('DECAY_23')
      MDL_WZ = VALUE
      CASE ('WW')
      MDL_WW = VALUE
      CASE ('DECAY_24')
      MDL_WW = VALUE
      CASE ('WH')
      MDL_WH = VALUE
      CASE ('DECAY_25')
      MDL_WH = VALUE
      CASE DEFAULT
      WRITE(*,*) 'no parameter matching', NAME, VALUE
      END SELECT

      RETURN
      END

      SUBROUTINE UPDATE_ALL_COUP()
      IMPLICIT NONE
      CALL COUP()
      RETURN
      END


      SUBROUTINE GET_PDG_ORDER(PDG, ALLPROC)
      IMPLICIT NONE
CF2PY INTEGER, intent(out) :: PDG(30,6)
CF2PY INTEGER, intent(out) :: ALLPROC(30)
      INTEGER PDG(30,6), PDGS(30,6)
      INTEGER ALLPROC(30),PIDS(30)
      DATA PDGS/ 2,2,2,1,1,1,3,3,3,4,4,4,5,5,5,2,2,2,1,1,1,3,3,3,4,4,4
     $ ,5,5,5,-2,-2,-2,-1,-1,-1,-3,-3,-3,-4,-4,-4,-5,-5,-5,-2,-2,-2,-1
     $ ,-1,-1,-3,-3,-3,-4,-4,-4,-5,-5,-5,-11,-13,-15,-11,-13,-15,-11,
     $ -13,-15,-11,-13,-15,-11,-13,-15,-11,-11,-13,-11,-11,-13,-11,-11
     $ ,-13,-11,-11,-13,-11,-11,-13,11,13,15,11,13,15,11,13,15,11,13
     $ ,15,11,13,15,11,11,13,11,11,13,11,11,13,11,11,13,11,11,13,-11,
     $ -13,-15,-11,-13,-15,-11,-13,-15,-11,-13,-15,-11,-13,-15,-13,-15
     $ ,-15,-13,-15,-15,-13,-15,-15,-13,-15,-15,-13,-15,-15,11,13,15
     $ ,11,13,15,11,13,15,11,13,15,11,13,15,13,15,15,13,15,15,13,15,15
     $ ,13,15,15,13,15,15 /
      DATA PIDS/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
     $ ,1,1,1 /
      PDG = PDGS
      ALLPROC = PIDS
      RETURN
      END

      SUBROUTINE GET_PREFIX(PREFIX)
      IMPLICIT NONE
CF2PY CHARACTER*20, intent(out) :: PREFIX(30)
      CHARACTER*20 PREFIX(30),PREF(30)
      DATA PREF / 'M0_','M1_','M2_','M3_','M4_','M5_','M6_','M7_'
     $ ,'M8_','M9_','M10_','M11_','M12_','M13_','M14_','M15_','M16_'
     $ ,'M17_','M18_','M19_','M20_','M21_','M22_','M23_','M24_','M25_'
     $ ,'M26_','M27_','M28_','M29_'/
      PREFIX = PREF
      RETURN
      END



