      DATA (NHEL(IHEL,  1),IHEL=1,5) / +1, +1, +1, +1, +1/
      DATA (NHEL(IHEL,  2),IHEL=1,5) / -1, +1, +1, +1, +1/
      DATA (NHEL(IHEL,  3),IHEL=1,5) / +1, -1, +1, +1, +1/
      DATA (NHEL(IHEL,  4),IHEL=1,5) / +1, +1, -1, +1, +1/
      DATA (NHEL(IHEL,  5),IHEL=1,5) / +1, +1, +1, -1, +1/
      DATA (NHEL(IHEL,  6),IHEL=1,5) / +1, +1, +1, +1, -1/
      DATA (NHEL(IHEL,  7),IHEL=1,5) / -1, -1, +1, +1, +1/
      DATA (NHEL(IHEL,  8),IHEL=1,5) / -1, +1, -1, +1, +1/
      DATA (NHEL(IHEL,  9),IHEL=1,5) / -1, +1, +1, -1, +1/
      DATA (NHEL(IHEL, 10),IHEL=1,5) / -1, +1, +1, +1, -1/
      DATA (NHEL(IHEL, 11),IHEL=1,5) / +1, -1, -1, +1, +1/
      DATA (NHEL(IHEL, 12),IHEL=1,5) / +1, -1, +1, -1, +1/
      DATA (NHEL(IHEL, 13),IHEL=1,5) / +1, -1, +1, +1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,5) / +1, +1, -1, -1, +1/
      DATA (NHEL(IHEL, 15),IHEL=1,5) / +1, +1, -1, +1, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,5) / +1, +1, +1, -1, -1/
      save NHEL
!$omp threadprivate(NHEL)
