      function BSYA1lcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA1lcqqppmp

C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (98), fully Badger-compliant
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
      real(dp):: xbeta2,xbeta4,s12,s23,mt2,mt3
      complex(dp):: BSYA0qqppmp,A0qqppmp,bracket
      integer:: e1,p2,p3,e4,j

      j=p2-1

      mt2=mt**2
      mt3=mt**3
      s12=mt2+s(1,p2)
      s23=s(p2,p3)
      xbeta2=1d0-4d0*mt**2/s23
      xbeta4=xbeta2**2

      A0qqppmp=BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      bracket=2d0*zab(p2,q1,p2)+s23

      BSYA1lcqqppmp=A0qqppmp*(
     & 0.5d0*s23*zab(p2,q1,p2)*I41x2x3x4(j)-0.5d0*s23*I32x3x41
     & +zab(p2,q1,p2)*I312x3x4(j)
     & -(s23+3d0*mt2*bracket/(s23*xbeta4))*I31x23x4
     & -8d0/3d0*I2m+23d0/9d0
     & +12d0*mt2*bracket/(s23**2*xbeta4)
     & -(7d0/6d0+bracket/(2d0*s23*xbeta2)
     & +6d0*mt2*bracket/(s23**2*xbeta4))*I2h23)

     & +mt*(za(p2,e1)*za(p2,e4)*zb(p2,p3)+za(e1,e4)*zab(p2,q1,p3))
     & /(2d0*zab(p2,q1,p3)**2)
     & *(zab(p2,q1,p2)**3*I41x2x3x4(j)-zab(p2,q1,p2)**2*I32x3x41
     & -2d0*zab(p2,q1,p2)**3/s23*I312x3x4(j)
     & +(s12*s23-2d0*mt2*zab(p2,q1,p2))*I31x23x4)

     & +(I31x23x4+2d0/s23*(I2h23-2d0))
     & *(6d0*mt3*(za(p3,e1)*za(p3,e4)*zab(p2,q1,p2)**2
     & +za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)**2)
     & /(za(p2,p3)*s23**2*xbeta4)
     & -6d0*mt3*za(p2,e1)*za(p3,e4)*zab(p3,q1,p2)
     & *bracket/(za(p2,p3)*s23**2*xbeta4)
     & -3d0*mt3*za(p3,e1)*za(p3,e4)/(2d0*za(p2,p3)*xbeta4))

     & -(I2h23-2d0)*(
     & mt*za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)*(zab(p2,q1,p2)+2d0*mt2)
     & /(za(p2,p3)*zab(p2,q1,p3)*s23*xbeta2)
     & +2d0*mt*za(p2,e1)*za(p3,e4)*zab(p3,q1,p2)/(za(p2,p3)*s23*xbeta2)
     & -mt*za(p3,e1)*za(p3,e4)*s12/(za(p2,p3)*s23*xbeta2))

     &  +I31x23x4*(
     &  mt*(za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)
     & +za(p3,e1)*za(p3,e4)*zab(p2,q1,p3))*zb(p2,p3)
     &  /(2d0*zab(p2,q1,p3) )
     & +mt3*za(p3,e1)*za(p3,e4)*(4d0*zab(p2,q1,p2)+s23)
     & /(2d0*za(p2,p3)*s23*xbeta2)
     & -mt3*za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)*bracket
     & /(za(p2,p3)*zab(p2,q1,p3)*s23*xbeta2)
     & -4d0*mt3*za(p2,e1)*za(p3,e4)*zab(p3,q1,p2)
     & /(za(p2,p3)*s23*xbeta2))

     &  + I2h23*(
     &  mt*(za(p3,e1)*za(p3,e4)*zab(p2,q1,p2)**2
     &  +za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)**2)
     & /(za(p2,p3)*s23**2*xbeta2)
     &  -mt*za(p2,e1)*za(p3,e4)*zab(p3,q1,p2)*bracket
     & /(za(p2,p3)*s23**2*xbeta2)
     & -mt3*za(p3,e1)*za(p3,e4)/(za(p2,p3)*s23*xbeta2))

     &  + F212(j)*(
     & (mt*za(p2,e1)*za(p2,e4)*zab(p3,q1,p2)**2)
     & /(s23*za(p2,p3)*zab(p2,q1,p2))
     & +(mt*za(p2,e1)*za(p2,e4)*zab(p2,q1,p2)*zab(p3,q1,p2))
     & /(s23*za(p2,p3)*zab(p2,q1,p3))
     & -(mt*za(p3,e4)*(za(p2,e1)*zab(p3,q1,p2)+mt2*za(p3,e1)))
     & /(za(p2,p3)*zab(p2,q1,p2)))
C---  overall sign change as required by the LHS of eq. (98)
      BSYA1lcqqppmp=-BSYA1lcqqppmp

      return
      end
