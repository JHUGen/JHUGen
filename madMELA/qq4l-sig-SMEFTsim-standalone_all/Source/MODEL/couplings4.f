ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP4()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1171 = -((MDL_CEBRE*MDL_COMPLEXI*MDL_STH*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1261 = -((MDL_COMPLEXI*MDL_YUP)/MDL_SQRT__2)
      GC_1313 = -((MDL_CHBOX*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1314 = (MDL_CHDD*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(4.000000D+00*MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1315 = (MDL_CHL3*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1316 = -(MDL_CLL1*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(2.000000D+00*MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1317 = -((MDL_CUHIM*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1318 = (MDL_CUHRE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      END
