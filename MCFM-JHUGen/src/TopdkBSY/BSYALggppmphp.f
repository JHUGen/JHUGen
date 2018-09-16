      function BSYALggppmphp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYALggppmphp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], Eq. (94)
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
      real(dp):: xbeta2,s23,mt3
      complex(dp):: BSYA0ggppmp,AT0ggppmp,tmp
      integer:: e1,p2,p3,e4,j
      j=p2-1
      mt3=mt**3
      s23=s(p2,p3)
      xbeta2=1d0-4d0*mt**2/s23

c---- SB see changes in BSYALggppmp.f
      AT0ggppmp=-BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      tmp =  + I41x2x3x4(j) * (  - 1.D0/2.D0*za(e1,p3)*za(e4,p3
     &    )*zb(p3,p2)*zab(p3,q1,p2)*mt - zab(p2,q1,p2)*s23*AT0ggppmp )
      tmp = tmp + I31x23x4 * (  - 6.D0/(za(p3,p2))*za(
     &    e1,p2)*za(e4,p3)*zab(p3,q1,p2)**2*s23**(-2)*xbeta2**(-2)*mt3
     &     - 6.D0/(za(p3,p2))*za(e1,p3)*za(e4,p2)*zab(p3,q1,p2)**2*
     &    s23**(-2)*xbeta2**(-2)*mt3 + 1/(zab(p2,q1,p2))/(za(p3,p2))*
     &    za(p3,p2)*za(e4,e1)*zab(p3,q1,p2)**2*mt + 2.D0/(zab(p2,q1,p2)
     &    )/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)**2*mt - 1.D0/
     &    2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*
     &    za(e4,p3)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*mt + 1.D0/2.D0/(zab(
     &    p2,q1,p2))/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*
     &    zab(p2,q1,p2)*zab(p3,q1,p2)*s23*mt - 1/(zab(p2,q1,p2))/(zab(
     &    p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p3,q1,p2)*s23*
     &    mt3 - 1.D0/4.D0/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)
     &    *zab(p3,q1,p2)**2*s23**(-1)*mt*xbeta2**(-1) + 1.D0/4.D0/(zab(
     &    p2,q1,p3))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p2,q1,p3)*zab(p3
     &    ,q1,p2)*s23**(-1)*mt*xbeta2**(-1) - 1/(zab(p2,q1,p3))*za(e4,
     &    e1)*zab(p3,q1,p2)*mt3 )
      tmp = tmp + I31x23x4 * (  - 6.D0/(zab(p2,q1,p3))
     &    /(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)
     &    **2*s23**(-2)*xbeta2**(-2)*mt3 - 3.D0/4.D0/(zab(p2,q1,p3))/(
     &    za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,p2)**2*mt*
     &    xbeta2**(-2) + 1.D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2
     &    )*za(e4,p2)*zab(p3,q1,p2)**2*mt + 1/(zab(p2,q1,p3))/(za(p3,p2
     &    ))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt*
     &    xbeta2**(-1) + 1/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,
     &    p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt + 2.D0/(zab(p2,q1,p3))/(
     &    za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)*xbeta2**(-1)*mt3
     &     + 2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(
     &    p3,q1,p2)*s23*mt + 6.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)
     &    *za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*
     &    s23**(-2)*xbeta2**(-2)*mt3 + 3.D0/4.D0/(zab(p2,q1,p3))/(za(p3
     &    ,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt*
     &    xbeta2**(-2) )
      tmp = tmp + I31x23x4 * ( 2.D0/(zab(p2,q1,p3))/(
     &    zab(p2,q1,p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(p2,q1,p2)*
     &    zab(p3,q1,p2)*s23**(-1)*xbeta2**(-1)*mt3 + 1/(zab(p2,q1,p3))
     &    /(zab(p2,q1,p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(p3,q1,p2)*
     &    xbeta2**(-1)*mt3 )
      tmp = tmp + I32x3x41 * ( 1.D0/2.D0/(zab(p2,q1,p2)
     &    )*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p2)*mt )
      tmp = tmp + I312x3x4(j) * (  - 1/(za(p3,p2))*za(
     &    e1,p3)*za(e4,p3)*zab(p3,q1,p2)*mt )
      tmp = tmp + I461x2x3x4(j) * ( 1/(zab(p2,q1,p2))/(
     &    zab(p2,q1,p2))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p2)*
     &    mt3 + 2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3
     &    )*zb(p3,p2)*zab(p3,q1,p2)*mt3 - 1/(zab(p2,q1,p2))/(zab(p2,q1,
     &    p3))*za(e4,e1)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*mt + 1/(zab(p2,
     &    q1,p2))/(zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p3)*zab(p3,q1,p2)
     &    **2*mt - 2.D0/(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p2)*mt3 +
     &    1/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4
     &    ,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)*s23*mt + 1/(zab(p2,q1,p3))/(
     &    zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p3)*
     &    zab(p3,q1,p2)**2*mt )
      tmp = tmp + I2h23 * (  - 12.D0/(za(p3,p2))*za(e1,
     &    p2)*za(e4,p3)*zab(p3,q1,p2)**2*s23**(-3)*xbeta2**(-2)*mt3 -
     &    1/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)**2*s23**(-2)*
     &    mt*xbeta2**(-1) - 12.D0/(za(p3,p2))*za(e1,p3)*za(e4,p2)*zab(
     &    p3,q1,p2)**2*s23**(-3)*xbeta2**(-2)*mt3 - 1/(za(p3,p2))*za(e1
     &    ,p3)*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-2)*mt*xbeta2**(-1) -
     &    1/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p2)
     &    *s23**(-1)*mt + 1/(zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)*
     &    zab(p3,q1,p2)*s23**(-1)*mt - 12.D0/(zab(p2,q1,p3))/(za(p3,p2)
     &    )*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)**2*
     &    s23**(-3)*xbeta2**(-2)*mt3 - 1/(zab(p2,q1,p3))/(za(p3,p2))*
     &    za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)**2*s23**(-2)*
     &    mt*xbeta2**(-1) - 3.D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1
     &    ,p2)*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-1)*mt*xbeta2**(-2) + 2.
     &    D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,
     &    p2)*zab(p3,q1,p2)*s23**(-1)*mt*xbeta2**(-1) )
      tmp = tmp + I2h23 * ( 1/(zab(p2,q1,p3))/(za(p3,p2
     &    ))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)*mt*xbeta2**(-1) + 12.D0
     &    /(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,p2
     &    )*zab(p2,q1,p3)*zab(p3,q1,p2)*s23**(-3)*xbeta2**(-2)*mt3 +
     &    1/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,
     &    p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*s23**(-2)*mt*xbeta2**(-1) + 3.
     &    D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(
     &    p2,q1,p3)*zab(p3,q1,p2)*s23**(-1)*mt*xbeta2**(-2) - 1/(zab(p2
     &    ,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(
     &    p2,q1,p2)*zab(p3,q1,p2)*mt*xbeta2**(-1) - 2.D0/(zab(p2,q1,p3)
     &    )/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,
     &    p2)*xbeta2**(-1)*mt3 )

c---- SB see changes in BSYALggppmp.f
      BSYALggppmphp=-tmp

      return
      end
