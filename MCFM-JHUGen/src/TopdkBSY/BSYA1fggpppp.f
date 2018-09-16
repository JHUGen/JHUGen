      function BSYA1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA1fggpppp

C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (100),fully Badger-compliant
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      real(dp):: s23
      integer:: e1,p2,p3,e4

      s23=s(p2,p3)
      BSYA1fggpppp=
     & -2d0*mt*(za(e1,e4)*zab(p2,q1,p2)-za(p2,e1)*za(p3,e4)*zb(p2,p3))
     & /(za(p2,p3)**3*zb(p2,p3))*(s23/6d0)

      return
      end

