      function A41HAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41HAQggmppp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A41phiAQggmppp,A41phiAQggmpmm
      integer:: j1,j2,j3,j4

      A41HAQggmppp=
     &   A41phiAQggmppp(j1,j2,j3,j4,za,zb)
     &  -A41phiAQggmpmm(j2,j1,j4,j3,zb,za)

      return
      end

      function A41HAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41HAQggmpmm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A41phiAQggmpmm,A41phiAQggmppp
      integer:: j1,j2,j3,j4

      A41HAQggmpmm=
     &    A41phiAQggmpmm(j1,j2,j3,j4,za,zb)
     &   -A41phiAQggmppp(j2,j1,j4,j3,zb,za)

      return
      end

      function A41HAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41HAQggmpmp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A41phiAQggmpmp
      integer:: j1,j2,j3,j4

      A41HAQggmpmp=
     &    A41phiAQggmpmp(j1,j2,j3,j4,za,zb)
     &   -A41phiAQggmpmp(j2,j1,j4,j3,zb,za)

      return
      end

      function A41HAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41HAQggmppm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A41phiAQggmppm
      integer:: j1,j2,j3,j4

      A41HAQggmppm=
     &    A41phiAQggmppm(j1,j2,j3,j4,za,zb)
     &   -A41phiAQggmppm(j2,j1,j4,j3,zb,za)

      return
      end

