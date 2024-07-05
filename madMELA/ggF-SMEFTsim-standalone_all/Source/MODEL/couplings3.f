ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP3()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_225 = (4.000000D+00*MDL_COMPLEXI*MDL_GHGG1)/MDL_VEVHAT
      GC_229 = (-4.000000D+00*MDL_COMPLEXI*MDL_GHGG2)/(MDL_MT__EXP__2
     $ *MDL_VEVHAT)
      GC_234 = (-2.000000D+00*MDL_COMPLEXI*MDL_GHGG4)/(MDL_MT__EXP__2
     $ *MDL_VEVHAT)
      GC_238 = (MDL_COMPLEXI*MDL_GHGG5)/(MDL_MT__EXP__2*MDL_VEVHAT)
      END
