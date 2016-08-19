      double complex function A41HAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex A41phiAQggmppp,A41phiAQggmpmm
      integer j1,j2,j3,j4
      
      A41HAQggmppp= 
     .   A41phiAQggmppp(j1,j2,j3,j4,za,zb)  
     .  -A41phiAQggmpmm(j2,j1,j4,j3,zb,za)
      
      return
      end
      
      double complex function A41HAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex A41phiAQggmpmm,A41phiAQggmppp
      integer j1,j2,j3,j4
      
      A41HAQggmpmm= 
     .    A41phiAQggmpmm(j1,j2,j3,j4,za,zb)
     .   -A41phiAQggmppp(j2,j1,j4,j3,zb,za)
      
      return
      end
      
      double complex function A41HAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex A41phiAQggmpmp
      integer j1,j2,j3,j4
      
      A41HAQggmpmp= 
     .    A41phiAQggmpmp(j1,j2,j3,j4,za,zb)
     .   -A41phiAQggmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
      
      double complex function A41HAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex A41phiAQggmppm
      integer j1,j2,j3,j4
      
      A41HAQggmppm= 
     .    A41phiAQggmppm(j1,j2,j3,j4,za,zb)
     .   -A41phiAQggmppm(j2,j1,j4,j3,zb,za)
      
      return
      end
      
