      real(dp):: p25(4),p34(4),p134(4),p346(4),p1346(4)
      real(dp):: p16(4),p234(4),p345(4),p156(4),p235(4)
      real(dp):: p12(4),p26(4)
      real(dp):: s34,s25,s346,s134,s16,s234,s345,s156,s235
      real(dp):: s12,s26
      common/kininvc/p25,p34,p134,p1346,p346,p16,p234,p345,p235,
     &     s34,s25,s346,s134,s16,s234,s345,p156,s156,s235,p12,p26,s12,
     & s26
!$omp threadprivate(/kininvc/)
