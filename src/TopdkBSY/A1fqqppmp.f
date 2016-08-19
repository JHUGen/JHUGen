      double complex function A1fqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (99)
      include 'constants.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'massiveintegrals.f'
      include 'epinv.f'
      double complex BSYA0qqppmp
      integer e1,p2,p3,e4
 
      A1fqqppmp=2d0/3d0*BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
     & *(epinv+I23-1d0/3d0)

      return
      end


