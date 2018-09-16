      real(dp):: xmass,medmass,dm_lam
      real(dp):: medwidth
      real(dp):: gdm,g_dmx,g_dmq
      logical:: effective_th
      character*6 dm_mediator 
      real(dp):: dmL(5),dmR(5)
      logical:: yukawa_scal
      common/yuk_scal/yukawa_scal
      common/dm_params/xmass,medmass,dm_lam,medwidth
      common/dm_coup/dmL,dmR
      common/dm_med/dm_mediator
      common/effec_dm/effective_th
      common/dm_g/gdm,g_dmx,g_dmq
