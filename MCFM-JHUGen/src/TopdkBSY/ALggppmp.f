      double complex function ALggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (92)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      double precision mt2
      double complex BSYA0ggppmp,BSYALggppmp,VL,lnrat
      integer e1,p2,p3,e4
 
      mt2=mt**2
      VL=-2d0*dcmplx(epinv**2)+dcmplx(0.5d0*epinv)
     & -epinv*(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2)))
     & -epinv*lnrat(musq,-s(p2,p3))
      ALggppmp=VL*BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &           +BSYALggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
