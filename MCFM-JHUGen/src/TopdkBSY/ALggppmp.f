      function ALggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: ALggppmp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (92)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: mt2
      complex(dp):: BSYA0ggppmp,BSYALggppmp,VL,lnrat
      integer:: e1,p2,p3,e4

      mt2=mt**2
      VL=-2d0*cplx1(epinv**2)+cplx1(0.5d0*epinv)
     & -epinv*(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2)))
     & -epinv*lnrat(musq,-s(p2,p3))
      ALggppmp=VL*BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &           +BSYALggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
