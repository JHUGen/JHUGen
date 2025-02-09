ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP2()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_319 = (-4.000000D+00*MDL_CHB*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_320 = (-4.000000D+00*MDL_CHBTIL*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_321 = (4.000000D+00*MDL_CHW*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_322 = (-4.000000D+00*MDL_CHWB*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_323 = (4.000000D+00*MDL_CHWB*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_324 = (-2.000000D+00*MDL_CHWBTIL*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_325 = (2.000000D+00*MDL_CHWBTIL*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_326 = (4.000000D+00*MDL_CHWTIL*MDL_CTH*MDL_COMPLEXI*MDL_STH
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_327 = (4.000000D+00*MDL_CHB*MDL_COMPLEXI*MDL_STH__EXP__2
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_328 = (-2.000000D+00*MDL_CHBTIL*MDL_COMPLEXI*MDL_STH__EXP__2
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_329 = (4.000000D+00*MDL_CHW*MDL_COMPLEXI*MDL_STH__EXP__2
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_330 = (-2.000000D+00*MDL_CHWTIL*MDL_COMPLEXI*MDL_STH__EXP__2
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_350 = (MDL_CHD*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_352 = (MDL_CHE*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_353 = (MDL_CHL1*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_354 = (MDL_CHQ1*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_355 = -(MDL_CHQ3*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_356 = (MDL_CHQ3*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_357 = (MDL_CHU*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_361 = (2.000000D+00*MDL_CHWB*MDL_EE*MDL_COMPLEXI
     $ *MDL_VEVHAT__EXP__2)/(3.000000D+00*MDL_LAMBDASMEFT__EXP__2
     $ *(1.000000D+00-2.000000D+00*MDL_STH__EXP__2))
      GC_363 = (MDL_CHWB*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(3.000000D+00*MDL_LAMBDASMEFT__EXP__2*(-1.000000D+00+2.000000D
     $ +00*MDL_STH__EXP__2))
      GC_364 = (MDL_CHWB*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(MDL_LAMBDASMEFT__EXP__2*(-1.000000D+00+2.000000D+00
     $ *MDL_STH__EXP__2))
      GC_383 = -(MDL_CHDD*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(8.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH*(
     $ -1.000000D+00+2.000000D+00*MDL_STH__EXP__2))
      GC_384 = (MDL_CHDD*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(8.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH*(
     $ -1.000000D+00+2.000000D+00*MDL_STH__EXP__2))
      GC_385 = -(MDL_CHL3*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2)
     $ /(2.000000D+00*MDL_CTH*MDL_LAMBDASMEFT__EXP__2*MDL_STH*(
     $ -1.000000D+00+2.000000D+00*MDL_STH__EXP__2))
      END
