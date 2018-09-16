      logical creategrid
      integer nSubProcess
      common/grid/creategrid,nSubProcess

      double precision weightb( -nf:nf, -nf:nf)
      double precision weightv( -nf:nf, -nf:nf)
      double precision weightv1( -nf:nf, -nf:nf)
      double precision weightv2( -nf:nf, -nf:nf)
      double precision weightr( 0:maxd , -nf:nf, -nf:nf)
      double precision weightfactor
      common/gridweight/
     .     weightfactor,
     .     weightb,
     .     weightv,weightv1,weightv2,
     .     weightr
      
      integer contrib,dipole
      double precision ag_xx1,ag_xx2,ag_x1z,ag_x2z,ag_scale,refwt,refwt2
      common/gridevent/
     .     ag_xx1,ag_xx2,ag_x1z,ag_x2z,
     .     ag_scale,refwt,refwt2,
     .     contrib,dipole
!$omp threadprivate(/gridweight/,/gridevent/)
      

