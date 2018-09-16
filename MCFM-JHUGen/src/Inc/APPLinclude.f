      logical:: creategrid
      integer:: nSubProcess
      common/grid/creategrid,nSubProcess

      real(dp):: weightb( -nf:nf, -nf:nf)
      real(dp):: weightv( -nf:nf, -nf:nf)
      real(dp):: weightv1( -nf:nf, -nf:nf)
      real(dp):: weightv2( -nf:nf, -nf:nf)
      real(dp):: weightr( 0:maxd , -nf:nf, -nf:nf)
      real(dp):: weightfactor
      common/gridweight/
     &     weightfactor,
     &     weightb,
     &     weightv,weightv1,weightv2,
     &     weightr

      integer:: contrib,dipole
      real(dp):: ag_xx1,ag_xx2,ag_x1z,ag_x2z,ag_scale,refwt,refwt2
      common/gridevent/
     &     ag_xx1,ag_xx2,ag_x1z,ag_x2z,
     &     ag_scale,refwt,refwt2,
     &     contrib,dipole
!$omp threadprivate(/gridweight/,/gridevent/)


