      function A0phiggggpppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiggggpppp

C----Expresssion of Eq. (4) of hep-ph/0607139v2
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phiggggpppp=czip
      return
      end

