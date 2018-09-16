      function A43HAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmppp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmm
      integer:: j1,j2,j3,j4

      A43HAQggmppp=
c     &   A43phiAQggmppp(j1,j2,j3,j4,za,zb)  ! This term is zero
     &  -A43phiAQggmpmm(j2,j1,j4,j3,zb,za)

      return
      end

      function A43HAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmpmm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmm
      integer:: j1,j2,j3,j4

      A43HAQggmpmm=
     &    A43phiAQggmpmm(j1,j2,j3,j4,za,zb)
c     &   -A43phiAQggmppp(j2,j1,j4,j3,zb,za) ! This term is zero

      return
      end

      function A43HAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmpmp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmp
      integer:: j1,j2,j3,j4

      A43HAQggmpmp=
     &    A43phiAQggmpmp(j1,j2,j3,j4,za,zb)
     &   -A43phiAQggmpmp(j2,j1,j4,j3,zb,za)

      return
      end

      function A43HAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmppm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmppm
      integer:: j1,j2,j3,j4

      A43HAQggmppm=
     &    A43phiAQggmppm(j1,j2,j3,j4,za,zb)
     &   -A43phiAQggmppm(j2,j1,j4,j3,zb,za)

      return
      end

