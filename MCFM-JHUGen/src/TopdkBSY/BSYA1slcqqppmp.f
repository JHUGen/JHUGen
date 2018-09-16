      function BSYA1slcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA1slcqqppmp

C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (99), fully Badger-compliant
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
      real(dp):: xbeta2,s23
      complex(dp):: BSYA0qqppmp
      integer:: e1,p2,p3,e4

      s23=s(p2,p3)
      xbeta2=1d0-4d0*mt**2/s23

      BSYA1slcqqppmp=-BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)*(
     & s23*I32x3x41+(s23-2d0*mt**2)*I3m1x23x4+1.5d0*I2h23
     & +(1d0+(2d0*zab(p2,q1,p2)+s23)/(2d0*s23*xbeta2))*F2m23
     &  +I2m-4d0
     &  -(2d0*zab(p2,q1,p2)+s23)/(s23*xbeta2))
     &  -(F2m23-2d0)*(
     & +mt**3*za(p3,e1)*za(p3,e4)/(za(p2,p3)*s23*xbeta2)
     & +mt*za(p2,e1)*za(p3,e4)*zab(p3,q1,p2)*(2d0*zab(p2,q1,p2)+s23)
     & /(za(p2,p3)*s23**2*xbeta2)
     & -mt*(za(p3,e1)*za(p3,e4)*zab(p2,q1,p2)**2
     &     +za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)**2)
     & /(za(p2,p3)*s23**2*xbeta2))

      return
      end
