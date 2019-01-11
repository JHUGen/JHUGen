c--- include file for controlling fragmentation/Frixione variables

c--- the following are set by the input file
      logical frag ! whether fragmentation process should be included
      character*8 fragset         ! label for fragmentation set
      double precision frag_scale ! fragmentation scale
      double precision cone_ang   ! cone size for Frixione isolation
      double precision epsilon_h  ! energy fraction for isolation
      double precision frag_scalestart ! frag. scale value in input file

c-- the following is used when computing fragmentation processes
      double precision z_frag ! energy fraction carried by photon
      logical rescale ! Indicates if p_part->1/z*p_gamma or not
c      double precision p_phys(mxpart,4) ! Physical momenta with jets
c                                        ! rescaled by factor of z to photons
      
      common/fraginputs/frag_scale,cone_ang,epsilon_h,frag_scalestart,
     & frag,fragset
      common/fragvars/z_frag,rescale
!===== logical variable to specific fragintmore 
      logical fragint_mode 
      common/fragint_mode/fragint_mode
!$omp threadprivate(/fragvars/)
