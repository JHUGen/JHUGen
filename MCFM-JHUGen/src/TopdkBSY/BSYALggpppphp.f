      double complex function BSYALggpppphp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (92), fully Badger-compliant
C---- higher-point (3- and 4-point) contributions only
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'massiveintegrals.f'
      double precision mt3
      integer e1,p2,p3,e4,j
 
C-----setup variable controlling integrals to be used,
C-----depending on whether p2=2 or 3
      j=p2-1

      mt3=mt**3
      BSYALggpppphp=
     & +I41x2x3x4(j)*mt3*za(e1,e4)*zb(p2,p3)**2

      return
      end
