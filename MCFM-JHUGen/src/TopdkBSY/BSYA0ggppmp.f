      function BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA0ggppmp

C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (A6), fully BSY compliant
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'sprods_com.f'
      include 'qdef.f'
      integer:: e1,p2,p3,e4
      BSYA0ggppmp=mt*zab(p3,q1,p2)
     & *(za(e1,e4)*zab(p3,q1,p2)-zb(p2,p3)*za(p3,e1)*za(p3,e4))
     & /(s(p2,p3)*zab(p2,q1,p2))

      return
      end
