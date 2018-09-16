c--- include file for controlling fragmentation/Frixione variables

c--- the following are set by the input file
      logical:: frag ! whether fragmentation process should be included
      character*8 fragset         ! label for fragmentation set
      real(dp):: frag_scale ! fragmentation scale
      real(dp):: cone_ang   ! cone size for Frixione isolation
      real(dp):: epsilon_h  ! energy fraction for isolation
      real(dp):: frag_scalestart ! frag. scale value in input file
      real(dp):: n_pow   ! exponent for smooth-cone (Frixione) isolation

c-- the following is used when computing fragmentation processes
      real(dp):: z_frag ! energy fraction carried by photon
      logical:: rescale ! Indicates if p_part->1/z*p_gamma or not
c      real(dp):: p_phys(mxpart,4) ! Physical momenta with jets
c                                        ! rescaled by factor of z to photons
      
      common/fraginputs/frag_scale,cone_ang,epsilon_h,frag_scalestart,
     & n_pow,frag,fragset
      common/fragvars/z_frag,rescale
!===== logical:: variable to specific fragintmore 
      logical:: fragint_mode 
      common/fragint_mode/fragint_mode
!$omp threadprivate(/fragvars/)
