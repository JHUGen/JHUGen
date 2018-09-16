      complex(dp)::
     & J34x0,J34x1(4),J34x3(4,4,4),
     & J61x0,J61x1(4),J61x3(4,4,4),J61x5(4,4,4,4,4),
     & J52x0,J52x1(4),J52x2(4,4),J52x3(4,4,4),
     & J52x4(4,4,4,4),J52x5(4,4,4,4,4)

      common/current/J34x0,J34x1,J34x3,
     & J61x0,J61x1,J61x3,J61x5,
     & J52x0,J52x1,J52x2,J52x3,J52x4,J52x5
!$omp threadprivate(/current/)
