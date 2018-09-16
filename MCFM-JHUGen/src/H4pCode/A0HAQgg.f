      function A0HAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0HAQggmppp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A0phidAQggmppp
      integer:: j1,j2,j3,j4

      A0HAQggmppp=
c     &  A0phiAQggmppp(j1,j2,j3,j4,za,zb)   ! This amplitude is zero
     &  +A0phidAQggmppp(j1,j2,j3,j4,za,zb)

      return
      end

      function A0HAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0HAQggmpmm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A0phiAQggmpmm
      integer:: j1,j2,j3,j4

      A0HAQggmpmm=
     &    A0phiAQggmpmm(j1,j2,j3,j4,za,zb)
c     &  +A0phidAQggmpmm(j1,j2,j3,j4,za,zb)   ! This amplitude is zero

      return
      end

      function A0HAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0HAQggmpmp

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A0phiAQggmpmp,A0phidAQggmpmp
      integer:: j1,j2,j3,j4

      A0HAQggmpmp=
     &    A0phiAQggmpmp(j1,j2,j3,j4,za,zb)
     &  +A0phidAQggmpmp(j1,j2,j3,j4,za,zb)

      return
      end

      function A0HAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0HAQggmppm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A0phiAQggmppm,A0phidAQggmppm
      integer:: j1,j2,j3,j4

      A0HAQggmppm=
     &    A0phiAQggmppm(j1,j2,j3,j4,za,zb)
     &  +A0phidAQggmppm(j1,j2,j3,j4,za,zb)

      return
      end

