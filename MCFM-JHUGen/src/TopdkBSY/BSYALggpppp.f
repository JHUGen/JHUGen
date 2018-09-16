      function BSYALggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYALggpppp

C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (92), fully Badger-compliant
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
      include 'massiveintegrals.f'
      real(dp):: s12,mt2
      complex(dp):: ze1xK12xK23xe4,BSYA0ggpppp,A0ggpppp
      integer:: e1,p2,p3,e4,j

C-----setup variable controlling integrals to be used,
C-----depending on whether p2=2 or 3
      j=p2-1

      mt2=mt**2
      s12=mt2+s(1,p2)

      A0ggpppp=BSYA0ggpppp(e1,p2,p3,e4,za,zb,zab,zba)

      ze1xK12xK23xe4=
     & zab(e1,q1,p2)*za(p2,e4)+zab(e1,q1,p3)*za(p3,e4)
     & +za(e1,p2)*zb(p2,p3)*za(p3,e4)

      BSYALggpppp =
     & +I41x2x3x4(j)*mt**3*za(e1,e4)*zb(p2,p3)**2
     & -F212(j)*(mt**3*zb(p2,p3)*(2d0*s12*za(e1,e4)-ze1xK12xK23xe4)
     &  /(za(p2,p3)*zab(p2,q1,p2)**2))
     & +0.5d0*A0ggpppp*(I2m-1d0)
     & -(mt*zb(p2,p3)*(za(e1,e4)*zab(p2,q1,p2)+ze1xK12xK23xe4)
     & /(2d0*za(p2,p3)*zab(p2,q1,p2))
     &  -mt*(za(e1,e4)*zab(p2,q1,p2)-za(p2,e1)*za(p3,e4)*zb(p2,p3))
     & /(3d0*za(p2,p3)**2))

      return
      end
