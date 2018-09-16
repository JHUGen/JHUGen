      function BSYA0qedpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA0qedpppp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], Eq. (A7)
C---- These are the twiddle functions arXiv:1101.5947 [hep-ph], Eq. (91)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      integer:: e1,p2,p3,e4
      BSYA0qedpppp=mt**3*za(e1,e4)*zb(p2,p3)**2
     & /(zab(p2,q1,p2)*zab(p3,q1,p3))

c--- flip sign to agree with BSY code
      BSYA0qedpppp=-BSYA0qedpppp

      return
      end
