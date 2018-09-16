      logical FixBrHbb
      real(dp):: GamHbb,GamHbb0,GamHbb1,mb_eff
      common/FixBrHbbFlag/FixBrHbb
      common/hbbparams/GamHbb,GamHbb0,GamHbb1,mb_eff
!$omp threadprivate(/hbbparams/)
