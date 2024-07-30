ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP13()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1189 = -((MDL_CLEQU1IM*MDL_YC*MDL_YTAU)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_1190 = (MDL_CLEQU1IM*MDL_YC*MDL_YTAU)/MDL_LAMBDASMEFT__EXP__2
      GC_1194 = -((MDL_CLEQU1RE*MDL_COMPLEXI*MDL_YC*MDL_YTAU)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_1198 = -(MDL_CLEQU3IM*MDL_YC*MDL_YTAU)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1199 = -(MDL_CLEQU3IM*MDL_YC*MDL_YTAU)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1200 = (MDL_CLEQU3IM*MDL_YC*MDL_YTAU)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1201 = (MDL_CLEQU3IM*MDL_YC*MDL_YTAU)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1208 = -(MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YC*MDL_YTAU)
     $ /(4.000000D+00*MDL_LAMBDASMEFT__EXP__2)
      GC_1209 = (MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YC*MDL_YTAU)/(2.000000D
     $ +00*MDL_LAMBDASMEFT__EXP__2)
      GC_1216 = -((MDL_CLEDQIM*MDL_YDO*MDL_YTAU)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_1217 = (MDL_CLEDQIM*MDL_YDO*MDL_YTAU)/MDL_LAMBDASMEFT__EXP__2
      GC_1221 = (MDL_CLEDQRE*MDL_COMPLEXI*MDL_YDO*MDL_YTAU)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_1225 = -((MDL_CLEDQIM*MDL_YS*MDL_YTAU)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_1226 = (MDL_CLEDQIM*MDL_YS*MDL_YTAU)/MDL_LAMBDASMEFT__EXP__2
      GC_1230 = (MDL_CLEDQRE*MDL_COMPLEXI*MDL_YS*MDL_YTAU)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_1261 = -((MDL_COMPLEXI*MDL_YUP)/MDL_SQRT__2)
      GC_1262 = -((MDL_CTH*MDL_CUBIM*MDL_YUP)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_1263 = (MDL_CTH*MDL_CUBRE*MDL_COMPLEXI*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1270 = -((MDL_CTH*MDL_CUWIM*MDL_YUP)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_1272 = (MDL_CTH*MDL_CUWRE*MDL_COMPLEXI*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1284 = (MDL_CUBIM*MDL_STH*MDL_YUP)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2)
      GC_1285 = -((MDL_CUBRE*MDL_COMPLEXI*MDL_STH*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1286 = -((MDL_CUWIM*MDL_STH*MDL_YUP)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_1287 = (MDL_CUWRE*MDL_COMPLEXI*MDL_STH*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1288 = -((MDL_CTH*MDL_CUBIM*MDL_VEVHAT*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      END
