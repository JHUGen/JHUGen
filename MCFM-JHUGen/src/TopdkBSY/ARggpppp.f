      function ARggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: ARggpppp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], Eq. (92)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      real(dp):: mt2,xlog,beta
      complex(dp):: BSYA0ggpppp,BSYARggpppp,VR
      integer:: e1,p2,p3,e4

      mt2=mt**2
      beta=sqrt(1d0-4d0*mt2/s(p2,p3))
      xlog=log((1d0-beta)/(1d0+beta))/beta
      VR=cplx1(0.5d0*epinv)
     & -cplx1(epinv)*(1d0-2d0*mt2/s(p2,p3))*xlog
      ARggpppp=VR*BSYA0ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     &           +BSYARggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
