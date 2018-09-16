      function BSYARggpppphp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYARggpppphp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], Eq. (93), fully Badger compliant
C---- higher-point (3- and 4-point) contributions only
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
      real(dp):: s23,mt3
      complex(dp):: ze1xK12xK23xe4
      integer:: e1,p2,p3,e4,j

C-----setup variable controlling integrals to be used,
C-----depending on whether p2=2 or 3
      j=p2-1
      mt3=mt**3
      s23=s(p2,p3)


      ze1xK12xK23xe4=
     & zab(e1,q1,p2)*za(p2,e4)+zab(e1,q1,p3)*za(p3,e4)
     & +za(e1,p2)*zb(p2,p3)*za(p3,e4)

      BSYARggpppphp=(
     & - F4m1x2x3x4(j)*(
     &  +za(p3,e1)*za(p3,e4)*(2d0*mt**2+zab(p2,q1,p2))*zb(p2,p3)**2
     &  *mt3/2/za(p2,p3)/zab(p3,q1,p2)
     &  -za(p2,e1)*za(p2,e4)*(2d0*mt**2+zab(p2,q1,p2))*zb(p2,p3)**2
     & *mt3/2/za(p2,p3)/zab(p2,q1,p3)
     &   +((2d0*mt**2-s23)*za(e1,e4)+ze1xK12xK23xe4)*zb(p2,p3)
     &  *mt3/za(p2,p3))

     & -I3m1x23x4*(2d0*mt**2-s23)*za(e1,e4)*zb(p2,p3)*mt3
     & /za(p2,p3)/zab(p2,q1,p2)

     & +I3m12x3x4(j)
     & *((2d0*za(e1,e4)*za(p2,p3)+4d0*za(p2,e4)*za(p3,e1))*zb(p2,p3)
     & *mt3/za(p2,p3)**2
     &  +za(p3,e1)*za(p3,e4)*zab(p2,q1,p3)*zb(p2,p3)
     & *mt3/za(p2,p3)**2/zab(p2,q1,p2)
     &  -za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)*zb(p2,p3)
     & *mt3/za(p2,p3)**2/zab(p2,q1,p2)
     &  +za(p2,e1)*za(p2,e4)*zab(p2,q1,p2)*zb(p2,p3)
     & *mt3/za(p2,p3)**2/zab(p2,q1,p3)
     &  -za(p3,e1)*za(p3,e4)*zab(p2,q1,p2)*zb(p2,p3)
     & *mt3/za(p2,p3)**2/zab(p3,q1,p2) )

     & + (
     &    za(p2,e1)*za(p2,e4)*zb(p3,p2)*zb(p2,p3)
     & *mt3/2d0/za(p2,p3)/zab(p2,q1,p3)
     &   -za(p3,e1)*za(p3,e4)*zb(p3,p2)*zb(p2,p3)
     & *mt3/2d0/za(p2,p3)/zab(p3,q1,p2)
     &   -za(e1,e4)*zb(p2,p3)*mt3/za(p2,p3)
     &   )*I3m2x3x41)

c---- sign required in accordance with LHS of eq. (93)
      BSYARggpppphp=-BSYARggpppphp

      return
      end
