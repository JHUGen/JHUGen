      subroutine vertices_bt2(mtsq,ep,
     &  vert25x8,vert25x9,vert25x10,
     &  vert25x11,vert25x12,vert25x13,vert25x14,vert25x15,
     &  vert25x16)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'poles.f'
      include 'scale.f'
      include 'masses.f'
      include 'decl_kininv.f'

      real(dp)::mtsq
      integer::ep
      complex(dp)::qlI2,qlI3,epqlI3,ep2qlI3
      complex(dp)::vert25x8,vert25x9
      complex(dp)::vert25x10,vert25x11,vert25x12,
     & vert25x13,vert25x14,vert25x15, vert25x16
      real(dp)::p5Dp234,denbt2

      mtsq=mt**2
      p5Dp234=(s16-s234-mtsq)/2d0
!      qlI2diffs234s34(ep)=qlI2(s234,0d0,0d0,musq,ep)
!     & -qlI2(s34,0d0,0d0,musq,ep)
     &
      denbt2=1d0/(mtsq*s234-p5Dp234**2)

      vert25x8= + qlI2(mtsq,zip,mtsq,musq,ep) * (  - 1.D0/2.D0*denbt2*
     &    s16 + 1.D0/8.D0*denbt2**2*s16**3 + 1.D0/2.D0*s234*denbt2 - 3.D
     &    0/8.D0*s234*denbt2**2*s16**2 + 3.D0/8.D0*s234**2*denbt2**2*
     &    s16 - 1.D0/8.D0*s234**3*denbt2**2 - mtsq*denbt2 - 3.D0/8.D0*
     &    mtsq*denbt2**2*s16**2 + 1.D0/2.D0*mtsq*s234*denbt2*s16**(-1)
     &     + mtsq*s234*denbt2**2*s16 - 5.D0/8.D0*mtsq*s234**2*denbt2**2
     &     - 1.D0/2.D0*mtsq**2*denbt2*s16**(-1) + 3.D0/8.D0*mtsq**2*
     &    denbt2**2*s16 - 5.D0/8.D0*mtsq**2*s234*denbt2**2 - 1.D0/8.D0*
     &    mtsq**3*denbt2**2 )
      vert25x8 = vert25x8 + qlI2(s234,zip,zip,musq,ep) * (  - denbt2*
     &    s16 + 3.D0/8.D0*s234*denbt2**2*s16**2 - 3.D0/4.D0*s234**2*
     &    denbt2**2*s16 + 3.D0/8.D0*s234**3*denbt2**2 + mtsq*denbt2 - 3.
     &    D0/4.D0*mtsq*s234*denbt2**2*s16 + 3.D0/4.D0*mtsq*s234**2*
     &    denbt2**2 + 3.D0/8.D0*mtsq**2*s234*denbt2**2 )
      vert25x8 = vert25x8 + qlI2(s16,zip,mtsq,musq,ep) * ( 3.D0/2.D0*
     &    denbt2*s16 - 1.D0/8.D0*denbt2**2*s16**3 - 1.D0/2.D0*s234*
     &    denbt2 + 3.D0/8.D0*s234**2*denbt2**2*s16 - 1.D0/4.D0*s234**3*
     &    denbt2**2 + 3.D0/8.D0*mtsq*denbt2**2*s16**2 - 1.D0/2.D0*mtsq*
     &    s234*denbt2*s16**(-1) - 1.D0/4.D0*mtsq*s234*denbt2**2*s16 - 1.
     &    D0/8.D0*mtsq*s234**2*denbt2**2 + 1.D0/2.D0*mtsq**2*denbt2*
     &    s16**(-1) - 3.D0/8.D0*mtsq**2*denbt2**2*s16 + 1.D0/4.D0*
     &    mtsq**2*s234*denbt2**2 + 1.D0/8.D0*mtsq**3*denbt2**2 )
      vert25x8 = vert25x8 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * ( denbt2*s16**2 - 2.D0*s234*denbt2*s16 - 1.D0/4.D0*s234*
     &    denbt2**2*s16**3 + s234**2*denbt2 + 3.D0/4.D0*s234**2*
     &    denbt2**2*s16**2 - 3.D0/4.D0*s234**3*denbt2**2*s16 + 1.D0/4.D0
     &    *s234**4*denbt2**2 - 2.D0*mtsq*denbt2*s16 + 3.D0/4.D0*mtsq*
     &    s234*denbt2**2*s16**2 - 5.D0/4.D0*mtsq*s234**2*denbt2**2*s16
     &     + 1.D0/2.D0*mtsq*s234**3*denbt2**2 + mtsq**2*denbt2 - 3.D0/4.
     &    D0*mtsq**2*s234*denbt2**2*s16 + 1.D0/2.D0*mtsq**2*s234**2*
     &    denbt2**2 + 1.D0/4.D0*mtsq**3*s234*denbt2**2 )
      vert25x8 = vert25x8 + 1.D0/2.D0*fp(ep)*denbt2*s16 - 1.D0/2.D0*fp(
     &    ep)*s234*denbt2 - fp(ep)*mtsq*denbt2 - 1.D0/2.D0*fp(ep)*mtsq*
     &    s234*denbt2*s16**(-1) + 1.D0/2.D0*fp(ep)*mtsq**2*denbt2*
     &    s16**(-1) - 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,
     &    ep)*mtsq*s234**2*denbt2**2*s16 + 1.D0/4.D0*epqlI3(s234,mtsq,
     &    s16,zip,zip,mtsq,musq,ep)*mtsq*s234**3*denbt2**2 + 1.D0/4.D0*
     &    epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*
     &    denbt2**2 - 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq
     &    ,ep)*mtsq*s234**2*denbt2**2*s16 + 1.D0/4.D0*ep2qlI3(s234,mtsq
     &    ,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**3*denbt2**2 + 1.D0/4.D0
     &    *ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*
     &    denbt2**2

      vert25x9= + qlI2(mtsq,zip,mtsq,musq,ep) * ( 1.D0/2.D0*denbt2*s16
     &     - 1.D0/8.D0*denbt2**2*s16**3 - 1.D0/2.D0*s234*denbt2 + 3.D0/
     &    8.D0*s234*denbt2**2*s16**2 - 3.D0/8.D0*s234**2*denbt2**2*s16
     &     + 1.D0/8.D0*s234**3*denbt2**2 - 1.D0/2.D0*mtsq*denbt2 - 3.D0/
     &    8.D0*mtsq*denbt2**2*s16**2 + 1.D0/2.D0*mtsq*s234*denbt2**2*
     &    s16 - 1.D0/8.D0*mtsq*s234**2*denbt2**2 + 9.D0/8.D0*mtsq**2*
     &    denbt2**2*s16 + 5.D0/8.D0*mtsq**2*s234*denbt2**2 - 5.D0/8.D0*
     &    mtsq**3*denbt2**2 )
      vert25x9 = vert25x9 + qlI2(s234,zip,zip,musq,ep) * ( 1.D0/2.D0*
     &    denbt2*s16 - 1.D0/4.D0*denbt2**2*s16**3 + 1.D0/2.D0*s234*
     &    denbt2 + 3.D0/8.D0*s234*denbt2**2*s16**2 - 1.D0/8.D0*s234**3*
     &    denbt2**2 - 1.D0/2.D0*mtsq*denbt2 + 3.D0/4.D0*mtsq*denbt2**2*
     &    s16**2 - 1.D0/2.D0*mtsq*s234*denbt2**2*s16 - 1.D0/4.D0*mtsq*
     &    s234**2*denbt2**2 - 3.D0/4.D0*mtsq**2*denbt2**2*s16 + 1.D0/8.D
     &    0*mtsq**2*s234*denbt2**2 + 1.D0/4.D0*mtsq**3*denbt2**2 )
      vert25x9 = vert25x9 + qlI2(s16,zip,mtsq,musq,ep) * (  - denbt2*
     &    s16 + 3.D0/8.D0*denbt2**2*s16**3 - 3.D0/4.D0*s234*denbt2**2*
     &    s16**2 + 3.D0/8.D0*s234**2*denbt2**2*s16 + mtsq*denbt2 - 3.D0/
     &    8.D0*mtsq*denbt2**2*s16**2 + 3.D0/8.D0*mtsq*s234**2*denbt2**2
     &     - 3.D0/8.D0*mtsq**2*denbt2**2*s16 - 3.D0/4.D0*mtsq**2*s234*
     &    denbt2**2 + 3.D0/8.D0*mtsq**3*denbt2**2 )
      vert25x9 = vert25x9 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * ( 1.D0/4.D0*denbt2**2*s16**4 - 3.D0/4.D0*s234*denbt2**2*
     &    s16**3 + 3.D0/4.D0*s234**2*denbt2**2*s16**2 - 1.D0/4.D0*
     &    s234**3*denbt2**2*s16 - mtsq*denbt2**2*s16**3 + 5.D0/4.D0*
     &    mtsq*s234*denbt2**2*s16**2 + 1.D0/4.D0*mtsq*s234**2*denbt2**2
     &    *s16 - 1.D0/2.D0*mtsq*s234**3*denbt2**2 + 3.D0/2.D0*mtsq**2*
     &    denbt2**2*s16**2 - 1.D0/4.D0*mtsq**2*s234*denbt2**2*s16 + 1.D0
     &    /2.D0*mtsq**2*s234**2*denbt2**2 - mtsq**3*denbt2**2*s16 - 1.D0
     &    /4.D0*mtsq**3*s234*denbt2**2 + 1.D0/4.D0*mtsq**4*denbt2**2 )
      vert25x9 = vert25x9 - 1.D0/2.D0*fp(ep)*denbt2*s16 + 1.D0/2.D0*fp(
     &    ep)*s234*denbt2 + 1.D0/2.D0*fp(ep)*mtsq*denbt2 + 1.D0/4.D0*
     &    epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2*
     &    denbt2**2*s16 - 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,
     &    musq,ep)*mtsq*s234**3*denbt2**2 + 1.D0/4.D0*epqlI3(s234,mtsq,
     &    s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*denbt2**2 + 1.D0/4.D
     &    0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2*
     &    denbt2**2*s16 - 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,
     &    musq,ep)*mtsq*s234**3*denbt2**2 + 1.D0/4.D0*ep2qlI3(s234,mtsq
     &    ,s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*denbt2**2

      vert25x10= + qlI2(mtsq,zip,mtsq,musq,ep) * (  - 1.D0 - 1.D0/16.D0
     &    *denbt2**2*s16**4 + 1.D0/4.D0*s234*denbt2**2*s16**3 - 3.D0/8.D
     &    0*s234**2*denbt2**2*s16**2 + 1.D0/4.D0*s234**3*denbt2**2*s16
     &     - 1.D0/16.D0*s234**4*denbt2**2 + mtsq*denbt2*s16 + 1.D0/4.D0
     &    *mtsq*denbt2**2*s16**3 + 1.D0/2.D0*mtsq*s234*denbt2 - 1.D0/4.D
     &    0*mtsq*s234*denbt2**2*s16**2 - 1.D0/4.D0*mtsq*s234**2*
     &    denbt2**2*s16 + 1.D0/4.D0*mtsq*s234**3*denbt2**2 - 3.D0/8.D0*
     &    mtsq**2*denbt2**2*s16**2 - 1.D0/4.D0*mtsq**2*s234*denbt2**2*
     &    s16 - 3.D0/8.D0*mtsq**2*s234**2*denbt2**2 + 1.D0/4.D0*mtsq**3
     &    *denbt2**2*s16 + 1.D0/4.D0*mtsq**3*s234*denbt2**2 - 1.D0/16.D0
     &    *mtsq**4*denbt2**2 )
      vert25x10 = vert25x10 + qlI2(s234,zip,zip,musq,ep) * ( 1.D0/2.D0*
     &    denbt2*s16**2 - 1.D0/2.D0*s234*denbt2*s16 - 1.D0/16.D0*s234*
     &    denbt2**2*s16**3 + 3.D0/16.D0*s234**2*denbt2**2*s16**2 - 3.D0/
     &    16.D0*s234**3*denbt2**2*s16 + 1.D0/16.D0*s234**4*denbt2**2 -
     &    1.D0/2.D0*mtsq*denbt2*s16 + 3.D0/16.D0*mtsq*s234*denbt2**2*
     &    s16**2 - 1.D0/8.D0*mtsq*s234**2*denbt2**2*s16 - 1.D0/16.D0*
     &    mtsq*s234**3*denbt2**2 - 3.D0/16.D0*mtsq**2*s234*denbt2**2*
     &    s16 - 1.D0/16.D0*mtsq**2*s234**2*denbt2**2 + 1.D0/16.D0*
     &    mtsq**3*s234*denbt2**2 )
      vert25x10 = vert25x10 + qlI2(s16,zip,mtsq,musq,ep) * (  - 1.D0/2.D
     &    0*denbt2*s16**2 + 1.D0/16.D0*denbt2**2*s16**4 + 1.D0/2.D0*
     &    s234*denbt2*s16 - 3.D0/16.D0*s234*denbt2**2*s16**3 + 3.D0/16.D
     &    0*s234**2*denbt2**2*s16**2 - 1.D0/16.D0*s234**3*denbt2**2*s16
     &     - 1.D0/2.D0*mtsq*denbt2*s16 - 1.D0/4.D0*mtsq*denbt2**2*
     &    s16**3 - 1.D0/2.D0*mtsq*s234*denbt2 + 1.D0/16.D0*mtsq*s234*
     &    denbt2**2*s16**2 + 3.D0/8.D0*mtsq*s234**2*denbt2**2*s16 - 3.D0
     &    /16.D0*mtsq*s234**3*denbt2**2 + 3.D0/8.D0*mtsq**2*denbt2**2*
     &    s16**2 + 7.D0/16.D0*mtsq**2*s234*denbt2**2*s16 + 7.D0/16.D0*
     &    mtsq**2*s234**2*denbt2**2 - 1.D0/4.D0*mtsq**3*denbt2**2*s16
     &     - 5.D0/16.D0*mtsq**3*s234*denbt2**2 + 1.D0/16.D0*mtsq**4*
     &    denbt2**2 )
      vert25x10 = vert25x10 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * (  - 1.D0/2.D0*denbt2*s16**3 + s234*denbt2*s16**2 - 1.D0/2.D0
     &    *s234**2*denbt2*s16 + mtsq*denbt2*s16**2 + 1.D0/2.D0*mtsq*
     &    s234**2*denbt2 + 1.D0/4.D0*mtsq*s234**2*denbt2**2*s16**2 - 1.D
     &    0/2.D0*mtsq*s234**3*denbt2**2*s16 + 1.D0/4.D0*mtsq*s234**4*
     &    denbt2**2 - 1.D0/2.D0*mtsq**2*denbt2*s16 - 1.D0/2.D0*mtsq**2*
     &    s234**2*denbt2**2*s16 - 1.D0/2.D0*mtsq**2*s234**3*denbt2**2
     &     + 1.D0/4.D0*mtsq**3*s234**2*denbt2**2 )
      vert25x10 = vert25x10 - 1.D0/4.D0*fp(ep)*denbt2*s16**2 + 1.D0/2.D0
     &    *fp(ep)*s234*denbt2*s16 - 1.D0/4.D0*fp(ep)*s234**2*denbt2 + 1.
     &    D0/2.D0*fp(ep)*mtsq*denbt2*s16 + 1.D0/2.D0*fp(ep)*mtsq*s234*
     &    denbt2 - 1.D0/4.D0*fp(ep)*mtsq**2*denbt2 + 1.D0/2.D0*epqlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2*denbt2 + 1.D0
     &    /8.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2
     &    *denbt2**2*s16**2 - 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,
     &    mtsq,musq,ep)*mtsq*s234**3*denbt2**2*s16 + 1.D0/8.D0*epqlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**4*denbt2**2 -
     &    1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2*
     &    s234**2*denbt2**2*s16 - 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,
     &    zip,mtsq,musq,ep)*mtsq**2*s234**3*denbt2**2 + 1.D0/8.D0*
     &    epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**3*s234**2*
     &    denbt2**2 + 1.D0/2.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq
     &    ,ep)*mtsq*s234**2*denbt2 + 1.D0/8.D0*ep2qlI3(s234,mtsq,s16,
     &    zip,zip,mtsq,musq,ep)*mtsq*s234**2*denbt2**2*s16**2
      vert25x10 = vert25x10 - 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,
     & mtsq,musq,ep)*mtsq*s234**3*denbt2**2*s16 + 1.D0/8.D0*ep2qlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**4*denbt2**2 -
     &    1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2
     &    *s234**2*denbt2**2*s16 - 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,
     &    zip,mtsq,musq,ep)*mtsq**2*s234**3*denbt2**2 + 1.D0/8.D0*
     &    ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**3*s234**2*
     &    denbt2**2

      vert25x11= + qlI2(mtsq,zip,mtsq,musq,ep) * ( 1.D0/8.D0*mt*
     &    denbt2**2*s16**3 + mt*s234*denbt2 + 1.D0/8.D0*mt*s234*
     &    denbt2**2*s16**2 - 5.D0/8.D0*mt*s234**2*denbt2**2*s16 + 3.D0/
     &    8.D0*mt*s234**3*denbt2**2 - 2.D0*mt*mtsq*denbt2 - 3.D0/8.D0*
     &    mt*mtsq*denbt2**2*s16**2 - 1.D0/8.D0*mt*mtsq*s234**2*
     &    denbt2**2 + 3.D0/8.D0*mt*mtsq**2*denbt2**2*s16 - 1.D0/8.D0*mt
     &    *mtsq**2*s234*denbt2**2 - 1.D0/8.D0*mt*mtsq**3*denbt2**2 )
      vert25x11 = vert25x11 + qlI2(s234,zip,zip,musq,ep) * (  - mt*
     &    denbt2*s16 + 3.D0/8.D0*mt*s234*denbt2**2*s16**2 - 3.D0/8.D0*
     &    mt*s234**3*denbt2**2 + mt*mtsq*denbt2 - 3.D0/4.D0*mt*mtsq*
     &    s234*denbt2**2*s16 + 3.D0/8.D0*mt*mtsq**2*s234*denbt2**2 )
      vert25x11 = vert25x11 + qlI2(s16,zip,mtsq,musq,ep) * ( mt*denbt2*
     &    s16 - 1.D0/8.D0*mt*denbt2**2*s16**3 - mt*s234*denbt2 - 1.D0/2.
     &    D0*mt*s234*denbt2**2*s16**2 + 5.D0/8.D0*mt*s234**2*denbt2**2*
     &    s16 + mt*mtsq*denbt2 + 3.D0/8.D0*mt*mtsq*denbt2**2*s16**2 + 3.
     &    D0/4.D0*mt*mtsq*s234*denbt2**2*s16 + 1.D0/8.D0*mt*mtsq*
     &    s234**2*denbt2**2 - 3.D0/8.D0*mt*mtsq**2*denbt2**2*s16 - 1.D0/
     &    4.D0*mt*mtsq**2*s234*denbt2**2 + 1.D0/8.D0*mt*mtsq**3*
     &    denbt2**2 )
      vert25x11 = vert25x11 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * ( mt*denbt2*s16**2 - 2.D0*mt*s234*denbt2*s16 - 1.D0/4.D0*mt*
     &    s234*denbt2**2*s16**3 + mt*s234**2*denbt2 + 1.D0/2.D0*mt*
     &    s234**2*denbt2**2*s16**2 - 1.D0/4.D0*mt*s234**3*denbt2**2*s16
     &     - 2.D0*mt*mtsq*denbt2*s16 + 3.D0/4.D0*mt*mtsq*s234*denbt2**2
     &    *s16**2 - 3.D0/4.D0*mt*mtsq*s234**2*denbt2**2*s16 - 1.D0/2.D0
     &    *mt*mtsq*s234**3*denbt2**2 + mt*mtsq**2*denbt2 - 3.D0/4.D0*mt
     &    *mtsq**2*s234*denbt2**2*s16 + 1.D0/4.D0*mt*mtsq**2*s234**2*
     &    denbt2**2 + 1.D0/4.D0*mt*mtsq**3*s234*denbt2**2 )
      vert25x11 = vert25x11 + fp(ep)*mt*s234*denbt2 - 1.D0/4.D0*epqlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mt*mtsq*s234**2*denbt2**2
     &    *s16 - 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*
     &    mt*mtsq*s234**3*denbt2**2 + 1.D0/4.D0*epqlI3(s234,mtsq,s16,
     &    zip,zip,mtsq,musq,ep)*mt*mtsq**2*s234**2*denbt2**2 - 1.D0/4.D0
     &    *ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mt*mtsq*s234**2*
     &    denbt2**2*s16 - 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,
     &    musq,ep)*mt*mtsq*s234**3*denbt2**2 + 1.D0/4.D0*ep2qlI3(s234,
     &    mtsq,s16,zip,zip,mtsq,musq,ep)*mt*mtsq**2*s234**2*denbt2**2

      vert25x12= + qlI2(mtsq,zip,mtsq,musq,ep) * (  - 1.D0/2.D0*mt*
     &    denbt2*s16 - 1.D0/4.D0*mt*denbt2**2*s16**3 - 1.D0/2.D0*mt*
     &    s234*denbt2 + 1.D0/4.D0*mt*s234*denbt2**2*s16**2 + 1.D0/4.D0*
     &    mt*s234**2*denbt2**2*s16 - 1.D0/4.D0*mt*s234**3*denbt2**2 - 3.
     &    D0/2.D0*mt*mtsq*denbt2 + 1.D0/2.D0*mt*mtsq*s234*denbt2**2*s16
     &     + 3.D0/4.D0*mt*mtsq**2*denbt2**2*s16 + 3.D0/4.D0*mt*mtsq**2*
     &    s234*denbt2**2 - 1.D0/2.D0*mt*mtsq**3*denbt2**2 )
      vert25x12 = vert25x12 + qlI2(s234,zip,zip,musq,ep) * (  - 1.D0/2.D
     &    0*mt*denbt2*s16 - 1.D0/4.D0*mt*denbt2**2*s16**3 + 1.D0/2.D0*
     &    mt*s234*denbt2 + 1.D0/4.D0*mt*s234**3*denbt2**2 + 1.D0/2.D0*
     &    mt*mtsq*denbt2 + 3.D0/4.D0*mt*mtsq*denbt2**2*s16**2 + 1.D0/4.D
     &    0*mt*mtsq*s234*denbt2**2*s16 - 1.D0/4.D0*mt*mtsq*s234**2*
     &    denbt2**2 - 3.D0/4.D0*mt*mtsq**2*denbt2**2*s16 - 1.D0/4.D0*mt
     &    *mtsq**2*s234*denbt2**2 + 1.D0/4.D0*mt*mtsq**3*denbt2**2 )
      vert25x12 = vert25x12 + qlI2(s16,zip,mtsq,musq,ep) * ( mt*denbt2*
     &    s16 + 1.D0/2.D0*mt*denbt2**2*s16**3 - 1.D0/4.D0*mt*s234*
     &    denbt2**2*s16**2 - 1.D0/4.D0*mt*s234**2*denbt2**2*s16 + mt*
     &    mtsq*denbt2 - 3.D0/4.D0*mt*mtsq*denbt2**2*s16**2 - 3.D0/4.D0*
     &    mt*mtsq*s234*denbt2**2*s16 + 1.D0/4.D0*mt*mtsq*s234**2*
     &    denbt2**2 - 1.D0/2.D0*mt*mtsq**2*s234*denbt2**2 + 1.D0/4.D0*
     &    mt*mtsq**3*denbt2**2 )
      vert25x12 = vert25x12 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * ( mt*denbt2*s16**2 + 1.D0/4.D0*mt*denbt2**2*s16**4 - mt*s234*
     &    denbt2*s16 - 1.D0/2.D0*mt*s234*denbt2**2*s16**3 + 1.D0/4.D0*
     &    mt*s234**2*denbt2**2*s16**2 - 2.D0*mt*mtsq*denbt2*s16 - mt*
     &    mtsq*denbt2**2*s16**3 - mt*mtsq*s234*denbt2 + 1.D0/2.D0*mt*
     &    mtsq*s234*denbt2**2*s16**2 + mt*mtsq*s234**2*denbt2**2*s16 +
     &    mt*mtsq**2*denbt2 + 3.D0/2.D0*mt*mtsq**2*denbt2**2*s16**2 + 1.
     &    D0/2.D0*mt*mtsq**2*s234*denbt2**2*s16 + 1.D0/4.D0*mt*mtsq**2*
     &    s234**2*denbt2**2 - mt*mtsq**3*denbt2**2*s16 - 1.D0/2.D0*mt*
     &    mtsq**3*s234*denbt2**2 + 1.D0/4.D0*mt*mtsq**4*denbt2**2 )
      vert25x12 = vert25x12 - 1.D0/2.D0*fp(ep)*mt*denbt2*s16 - 1.D0/2.D0
     &    *fp(ep)*mt*s234*denbt2 + 1.D0/2.D0*fp(ep)*mt*mtsq*denbt2 + 1.D
     &    0/2.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mt*mtsq*
     &    s234**2*denbt2**2*s16 + 1.D0/2.D0*ep2qlI3(s234,mtsq,s16,zip,
     &    zip,mtsq,musq,ep)*mt*mtsq*s234**2*denbt2**2*s16

      vert25x13= + qlI2(mtsq,zip,mtsq,musq,ep) * ( 1.D0/2.D0*mt*denbt2*
     &    s16 - 1.D0/2.D0*mt*s234*denbt2 + 1.D0/2.D0*mt*mtsq*denbt2 )
      vert25x13 = vert25x13 + qlI2(s234,zip,zip,musq,ep) * ( 1.D0/2.D0*
     &    mt*denbt2*s16 + 1.D0/2.D0*mt*s234*denbt2 - 1.D0/2.D0*mt*mtsq*
     &    denbt2 )
      vert25x13 = vert25x13 + qlI2(s16,zip,mtsq,musq,ep) * (  - mt*
     &    denbt2*s16 )
      vert25x13 = vert25x13 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * (  - 1.D0/2.D0*mt*denbt2*s16**2 + 1.D0/2.D0*mt*s234*denbt2*
     &    s16 + mt*mtsq*denbt2*s16 + 1.D0/2.D0*mt*mtsq*s234*denbt2 - 1.D
     &    0/2.D0*mt*mtsq**2*denbt2 )

      vert25x14= + qlI2(mtsq,zip,mtsq,musq,ep) * ( 1.D0/2.D0*denbt2*s16
     &     - 1.D0/8.D0*denbt2**2*s16**3 - 1.D0/2.D0*s234*denbt2 + 3.D0/
     &    8.D0*s234*denbt2**2*s16**2 - 3.D0/8.D0*s234**2*denbt2**2*s16
     &     + 1.D0/8.D0*s234**3*denbt2**2 + mtsq*denbt2 + 3.D0/8.D0*mtsq
     &    *denbt2**2*s16**2 - 1.D0/2.D0*mtsq*s234*denbt2*s16**(-1) -
     &    mtsq*s234*denbt2**2*s16 + 5.D0/8.D0*mtsq*s234**2*denbt2**2 +
     &    1.D0/2.D0*mtsq**2*denbt2*s16**(-1) - 3.D0/8.D0*mtsq**2*
     &    denbt2**2*s16 + 5.D0/8.D0*mtsq**2*s234*denbt2**2 + 1.D0/8.D0*
     &    mtsq**3*denbt2**2 )
      vert25x14 = vert25x14 + qlI2(s234,zip,zip,musq,ep) * ( denbt2*s16
     &     - 3.D0/8.D0*s234*denbt2**2*s16**2 + 3.D0/4.D0*s234**2*
     &    denbt2**2*s16 - 3.D0/8.D0*s234**3*denbt2**2 - mtsq*denbt2 + 3.
     &    D0/4.D0*mtsq*s234*denbt2**2*s16 - 3.D0/4.D0*mtsq*s234**2*
     &    denbt2**2 - 3.D0/8.D0*mtsq**2*s234*denbt2**2 )
      vert25x14 = vert25x14 + qlI2(s16,zip,mtsq,musq,ep) * (  - 3.D0/2.D
     &    0*denbt2*s16 + 1.D0/8.D0*denbt2**2*s16**3 + 1.D0/2.D0*s234*
     &    denbt2 - 3.D0/8.D0*s234**2*denbt2**2*s16 + 1.D0/4.D0*s234**3*
     &    denbt2**2 - 3.D0/8.D0*mtsq*denbt2**2*s16**2 + 1.D0/2.D0*mtsq*
     &    s234*denbt2*s16**(-1) + 1.D0/4.D0*mtsq*s234*denbt2**2*s16 + 1.
     &    D0/8.D0*mtsq*s234**2*denbt2**2 - 1.D0/2.D0*mtsq**2*denbt2*
     &    s16**(-1) + 3.D0/8.D0*mtsq**2*denbt2**2*s16 - 1.D0/4.D0*
     &    mtsq**2*s234*denbt2**2 - 1.D0/8.D0*mtsq**3*denbt2**2 )
      vert25x14 = vert25x14 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * (  - denbt2*s16**2 + 2.D0*s234*denbt2*s16 + 1.D0/4.D0*s234*
     &    denbt2**2*s16**3 - s234**2*denbt2 - 3.D0/4.D0*s234**2*
     &    denbt2**2*s16**2 + 3.D0/4.D0*s234**3*denbt2**2*s16 - 1.D0/4.D0
     &    *s234**4*denbt2**2 + 2.D0*mtsq*denbt2*s16 - 3.D0/4.D0*mtsq*
     &    s234*denbt2**2*s16**2 + 5.D0/4.D0*mtsq*s234**2*denbt2**2*s16
     &     - 1.D0/2.D0*mtsq*s234**3*denbt2**2 - mtsq**2*denbt2 + 3.D0/4.
     &    D0*mtsq**2*s234*denbt2**2*s16 - 1.D0/2.D0*mtsq**2*s234**2*
     &    denbt2**2 - 1.D0/4.D0*mtsq**3*s234*denbt2**2 )
      vert25x14 = vert25x14 - 1.D0/2.D0*fp(ep)*denbt2*s16 + 1.D0/2.D0*
     &    fp(ep)*s234*denbt2 + fp(ep)*mtsq*denbt2 + 1.D0/2.D0*fp(ep)*
     &    mtsq*s234*denbt2*s16**(-1) - 1.D0/2.D0*fp(ep)*mtsq**2*denbt2*
     &    s16**(-1) + 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,
     &    ep)*mtsq*s234**2*denbt2**2*s16 - 1.D0/4.D0*epqlI3(s234,mtsq,
     &    s16,zip,zip,mtsq,musq,ep)*mtsq*s234**3*denbt2**2 - 1.D0/4.D0*
     &    epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*
     &    denbt2**2 + 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq
     &    ,ep)*mtsq*s234**2*denbt2**2*s16 - 1.D0/4.D0*ep2qlI3(s234,mtsq
     &    ,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**3*denbt2**2 - 1.D0/4.D0
     &    *ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*
     &    denbt2**2

      vert25x15= + qlI2(mtsq,zip,mtsq,musq,ep) * (  - 1.D0/2.D0*denbt2*
     &    s16 + 1.D0/8.D0*denbt2**2*s16**3 + 1.D0/2.D0*s234*denbt2 - 3.D
     &    0/8.D0*s234*denbt2**2*s16**2 + 3.D0/8.D0*s234**2*denbt2**2*
     &    s16 - 1.D0/8.D0*s234**3*denbt2**2 + 1.D0/2.D0*mtsq*denbt2 + 3.
     &    D0/8.D0*mtsq*denbt2**2*s16**2 - 1.D0/2.D0*mtsq*s234*denbt2**2
     &    *s16 + 1.D0/8.D0*mtsq*s234**2*denbt2**2 - 9.D0/8.D0*mtsq**2*
     &    denbt2**2*s16 - 5.D0/8.D0*mtsq**2*s234*denbt2**2 + 5.D0/8.D0*
     &    mtsq**3*denbt2**2 )
      vert25x15 = vert25x15 + qlI2(s234,zip,zip,musq,ep) * (  - 1.D0/2.D
     &    0*denbt2*s16 + 1.D0/4.D0*denbt2**2*s16**3 - 1.D0/2.D0*s234*
     &    denbt2 - 3.D0/8.D0*s234*denbt2**2*s16**2 + 1.D0/8.D0*s234**3*
     &    denbt2**2 + 1.D0/2.D0*mtsq*denbt2 - 3.D0/4.D0*mtsq*denbt2**2*
     &    s16**2 + 1.D0/2.D0*mtsq*s234*denbt2**2*s16 + 1.D0/4.D0*mtsq*
     &    s234**2*denbt2**2 + 3.D0/4.D0*mtsq**2*denbt2**2*s16 - 1.D0/8.D
     &    0*mtsq**2*s234*denbt2**2 - 1.D0/4.D0*mtsq**3*denbt2**2 )
      vert25x15 = vert25x15 + qlI2(s16,zip,mtsq,musq,ep) * ( denbt2*s16
     &     - 3.D0/8.D0*denbt2**2*s16**3 + 3.D0/4.D0*s234*denbt2**2*
     &    s16**2 - 3.D0/8.D0*s234**2*denbt2**2*s16 - mtsq*denbt2 + 3.D0/
     &    8.D0*mtsq*denbt2**2*s16**2 - 3.D0/8.D0*mtsq*s234**2*denbt2**2
     &     + 3.D0/8.D0*mtsq**2*denbt2**2*s16 + 3.D0/4.D0*mtsq**2*s234*
     &    denbt2**2 - 3.D0/8.D0*mtsq**3*denbt2**2 )
      vert25x15 = vert25x15 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * (  - 1.D0/4.D0*denbt2**2*s16**4 + 3.D0/4.D0*s234*denbt2**2*
     &    s16**3 - 3.D0/4.D0*s234**2*denbt2**2*s16**2 + 1.D0/4.D0*
     &    s234**3*denbt2**2*s16 + mtsq*denbt2**2*s16**3 - 5.D0/4.D0*
     &    mtsq*s234*denbt2**2*s16**2 - 1.D0/4.D0*mtsq*s234**2*denbt2**2
     &    *s16 + 1.D0/2.D0*mtsq*s234**3*denbt2**2 - 3.D0/2.D0*mtsq**2*
     &    denbt2**2*s16**2 + 1.D0/4.D0*mtsq**2*s234*denbt2**2*s16 - 1.D0
     &    /2.D0*mtsq**2*s234**2*denbt2**2 + mtsq**3*denbt2**2*s16 + 1.D0
     &    /4.D0*mtsq**3*s234*denbt2**2 - 1.D0/4.D0*mtsq**4*denbt2**2 )
      vert25x15 = vert25x15 + 1.D0/2.D0*fp(ep)*denbt2*s16 - 1.D0/2.D0*
     &    fp(ep)*s234*denbt2 - 1.D0/2.D0*fp(ep)*mtsq*denbt2 - 1.D0/4.D0
     &    *epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2*
     &    denbt2**2*s16 + 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,
     &    musq,ep)*mtsq*s234**3*denbt2**2 - 1.D0/4.D0*epqlI3(s234,mtsq,
     &    s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*denbt2**2 - 1.D0/4.D
     &    0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2*
     &    denbt2**2*s16 + 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,
     &    musq,ep)*mtsq*s234**3*denbt2**2 - 1.D0/4.D0*ep2qlI3(s234,mtsq
     &    ,s16,zip,zip,mtsq,musq,ep)*mtsq**2*s234**2*denbt2**2

      vert25x16= + qlI2(mtsq,zip,mtsq,musq,ep) * ( 1.D0 + 1.D0/16.D0*
     &    denbt2**2*s16**4 - 1.D0/4.D0*s234*denbt2**2*s16**3 + 3.D0/8.D0
     &    *s234**2*denbt2**2*s16**2 - 1.D0/4.D0*s234**3*denbt2**2*s16
     &     + 1.D0/16.D0*s234**4*denbt2**2 - mtsq*denbt2*s16 - 1.D0/4.D0
     &    *mtsq*denbt2**2*s16**3 - 1.D0/2.D0*mtsq*s234*denbt2 + 1.D0/4.D
     &    0*mtsq*s234*denbt2**2*s16**2 + 1.D0/4.D0*mtsq*s234**2*
     &    denbt2**2*s16 - 1.D0/4.D0*mtsq*s234**3*denbt2**2 + 3.D0/8.D0*
     &    mtsq**2*denbt2**2*s16**2 + 1.D0/4.D0*mtsq**2*s234*denbt2**2*
     &    s16 + 3.D0/8.D0*mtsq**2*s234**2*denbt2**2 - 1.D0/4.D0*mtsq**3
     &    *denbt2**2*s16 - 1.D0/4.D0*mtsq**3*s234*denbt2**2 + 1.D0/16.D0
     &    *mtsq**4*denbt2**2 )
      vert25x16 = vert25x16 + qlI2(s234,zip,zip,musq,ep) * (  - 1.D0/2.D
     &    0*denbt2*s16**2 + 1.D0/2.D0*s234*denbt2*s16 + 1.D0/16.D0*s234
     &    *denbt2**2*s16**3 - 3.D0/16.D0*s234**2*denbt2**2*s16**2 + 3.D0
     &    /16.D0*s234**3*denbt2**2*s16 - 1.D0/16.D0*s234**4*denbt2**2
     &     + 1.D0/2.D0*mtsq*denbt2*s16 - 3.D0/16.D0*mtsq*s234*denbt2**2
     &    *s16**2 + 1.D0/8.D0*mtsq*s234**2*denbt2**2*s16 + 1.D0/16.D0*
     &    mtsq*s234**3*denbt2**2 + 3.D0/16.D0*mtsq**2*s234*denbt2**2*
     &    s16 + 1.D0/16.D0*mtsq**2*s234**2*denbt2**2 - 1.D0/16.D0*
     &    mtsq**3*s234*denbt2**2 )
      vert25x16 = vert25x16 + qlI2(s16,zip,mtsq,musq,ep) * ( 1.D0/2.D0*
     &    denbt2*s16**2 - 1.D0/16.D0*denbt2**2*s16**4 - 1.D0/2.D0*s234*
     &    denbt2*s16 + 3.D0/16.D0*s234*denbt2**2*s16**3 - 3.D0/16.D0*
     &    s234**2*denbt2**2*s16**2 + 1.D0/16.D0*s234**3*denbt2**2*s16
     &     + 1.D0/2.D0*mtsq*denbt2*s16 + 1.D0/4.D0*mtsq*denbt2**2*
     &    s16**3 + 1.D0/2.D0*mtsq*s234*denbt2 - 1.D0/16.D0*mtsq*s234*
     &    denbt2**2*s16**2 - 3.D0/8.D0*mtsq*s234**2*denbt2**2*s16 + 3.D0
     &    /16.D0*mtsq*s234**3*denbt2**2 - 3.D0/8.D0*mtsq**2*denbt2**2*
     &    s16**2 - 7.D0/16.D0*mtsq**2*s234*denbt2**2*s16 - 7.D0/16.D0*
     &    mtsq**2*s234**2*denbt2**2 + 1.D0/4.D0*mtsq**3*denbt2**2*s16
     &     + 5.D0/16.D0*mtsq**3*s234*denbt2**2 - 1.D0/16.D0*mtsq**4*
     &    denbt2**2 )
      vert25x16 = vert25x16 + qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)
     &  * ( 1.D0/2.D0*denbt2*s16**3 - s234*denbt2*s16**2 + 1.D0/2.D0*
     &    s234**2*denbt2*s16 - mtsq*denbt2*s16**2 - 1.D0/2.D0*mtsq*
     &    s234**2*denbt2 - 1.D0/4.D0*mtsq*s234**2*denbt2**2*s16**2 + 1.D
     &    0/2.D0*mtsq*s234**3*denbt2**2*s16 - 1.D0/4.D0*mtsq*s234**4*
     &    denbt2**2 + 1.D0/2.D0*mtsq**2*denbt2*s16 + 1.D0/2.D0*mtsq**2*
     &    s234**2*denbt2**2*s16 + 1.D0/2.D0*mtsq**2*s234**3*denbt2**2
     &     - 1.D0/4.D0*mtsq**3*s234**2*denbt2**2 )
      vert25x16 = vert25x16 + 1.D0/4.D0*fp(ep)*denbt2*s16**2 - 1.D0/2.D0
     &    *fp(ep)*s234*denbt2*s16 + 1.D0/4.D0*fp(ep)*s234**2*denbt2 - 1.
     &    D0/2.D0*fp(ep)*mtsq*denbt2*s16 - 1.D0/2.D0*fp(ep)*mtsq*s234*
     &    denbt2 + 1.D0/4.D0*fp(ep)*mtsq**2*denbt2 - 1.D0/2.D0*epqlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2*denbt2 - 1.D0
     &    /8.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**2
     &    *denbt2**2*s16**2 + 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,
     &    mtsq,musq,ep)*mtsq*s234**3*denbt2**2*s16 - 1.D0/8.D0*epqlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**4*denbt2**2 +
     &    1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2*
     &    s234**2*denbt2**2*s16 + 1.D0/4.D0*epqlI3(s234,mtsq,s16,zip,
     &    zip,mtsq,musq,ep)*mtsq**2*s234**3*denbt2**2 - 1.D0/8.D0*
     &    epqlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**3*s234**2*
     &    denbt2**2 - 1.D0/2.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq
     &    ,ep)*mtsq*s234**2*denbt2 - 1.D0/8.D0*ep2qlI3(s234,mtsq,s16,
     &    zip,zip,mtsq,musq,ep)*mtsq*s234**2*denbt2**2*s16**2
      vert25x16 = vert25x16 + 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,
     & mtsq,musq,ep)*mtsq*s234**3*denbt2**2*s16 - 1.D0/8.D0*ep2qlI3(
     &    s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq*s234**4*denbt2**2 +
     &    1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**2
     &    *s234**2*denbt2**2*s16 + 1.D0/4.D0*ep2qlI3(s234,mtsq,s16,zip,
     &    zip,mtsq,musq,ep)*mtsq**2*s234**3*denbt2**2 - 1.D0/8.D0*
     &    ep2qlI3(s234,mtsq,s16,zip,zip,mtsq,musq,ep)*mtsq**3*s234**2*
     &    denbt2**2


      return
      end

