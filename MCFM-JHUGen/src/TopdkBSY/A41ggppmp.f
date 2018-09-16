      function A41ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A41ggppmp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: ALggppmp,ARggppmp,
     & BSYA1fggppmp,BSYA1Hggppmp
      integer:: e1,p2,p3,e4
      real(dp):: nlf,nhf
      parameter(nlf=5d0,nhf=1d0)

      A41ggppmp=ALggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     & -Ncinv**2*ARggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     & -BSYA1fggppmp(e1,p2,p3,e4,za,zb,zab,zba)*nlf/xn
     & -BSYA1Hggppmp(e1,p2,p3,e4,za,zb,zab,zba)*nhf/xn

      return
      end

