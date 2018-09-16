      function A0phiggggmmpp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiggggmmpp

C----Expresssion of Eq. (2.14) of hep-ph/0804.4149v3
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phiggggmmpp=za(j1,j2)**4
     &             /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))
      return
      end

