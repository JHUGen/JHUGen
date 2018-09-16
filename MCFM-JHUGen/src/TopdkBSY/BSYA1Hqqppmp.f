      function BSYA1Hqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA1Hqqppmp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], Eq. (102)
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'sprods_com.f'
      include 'massiveintegrals.f'
      complex(dp):: BSYA0qqppmp
      integer:: e1,p2,p3,e4

      BSYA1Hqqppmp=
     &  2d0/3d0*BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)*(
     & (2d0*mt**2/S(p2,p3)+1d0)*F2m23+Im2m-1d0/3d0)

      return
      end

