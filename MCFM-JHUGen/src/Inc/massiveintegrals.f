      complex(dp):: 
     & I41x2x3x4(2),I4m1x2x3x4(2),I41x2x4x3,
     & I3m12x3x4(2),I312x3x4(2),
     & I3m13x2x4(2),I313x2x4(2),
     & F41x2x3x4(2),F4m1x2x3x4(2),
     & I461x2x3x4(2),I46m1x2x3x4(2),F212(2),
     & I3m1x23x4,I3m2x3x41,I32x3x41,I31x23x4,I2m,I2h23,F2m23,Im2m,I23
      common/massiveintegrals/
     & I41x2x3x4,I4m1x2x3x4,I41x2x4x3,
     & I3m12x3x4,I312x3x4,
     & I3m13x2x4,I313x2x4,
     & F41x2x3x4,F4m1x2x3x4,
     & I461x2x3x4,I46m1x2x3x4,F212,
     & I3m1x23x4,I3m2x3x41,I32x3x41,I31x23x4,
     & I2m,I2h23,F2m23,Im2m,I23
!$omp threadprivate(/massiveintegrals/)
