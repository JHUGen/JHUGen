      subroutine TWZggNAB(p1,p2,p3,p4,p5,p6,p7,p8,NAB)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
C     Author: R.K. Ellis Feb, 2013
C     written by program WZggdiags.frm
C     These are the non-Abelian diagrams with both the W and Z coming
C     from 3 boson vertex off light line
C     Calculation is performed for LH quark-line (perforce because of W)
C     Calculation is performed for LH 56line
C     The two indices of NAB are the gluon helicities
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s34,s56,s78,s278,s178,
     & s3456
      complex(dp):: zab2,zba3,zba2,NAB(2,2),iza,izb
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zba3(p1,p2,p3,p4,p5)=
     & +zb(p1,p2)*za(p2,p5)
     & +zb(p1,p3)*za(p3,p5)
     & +zb(p1,p4)*za(p4,p5)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
C     end statement functions
      s34=s(p3,p4)
      s56=s(p5,p6)
      s78=s(p7,p8)
      s278=s3(p2,p7,p8)
      s178=s3(p1,p7,p8)
      s3456=s4(p3,p4,p5,p6)
      NAB(1,1)= + s78**(-1)*s34**(-1)*s56**(-1)*s3456**(-1)*s278**(-1)
     &  * (  - za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p3)*zb(p4,p6)*izb(p1
     &    ,p8)*zba3(p1,p2,p7,p8,p3) + za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(
     &    p1,p4)*zb(p3,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) - za(p2,p7)*
     &    za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p4) - za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p6)*zb(p4,p5
     &    )*izb(p1,p8)*zba3(p1,p2,p7,p8,p5) - za(p2,p7)*za(p3,p6)*za(p7
     &    ,p8)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5) +
     &    za(p2,p7)*za(p4,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p3) - za(p2,p8)*za(p3,p5)*za(p7,p8)*zb(p1,p3
     &    )*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,p3) + za(p2,p8)*za(p3
     &    ,p5)*za(p7,p8)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p3) - za(p2,p8)*za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*
     &    izb(p1,p7)*zba3(p1,p2,p7,p8,p4) - za(p2,p8)*za(p3,p5)*za(p7,
     &    p8)*zb(p1,p6)*zb(p4,p5)*izb(p1,p7)*zba3(p1,p2,p7,p8,p5) - za(
     &    p2,p8)*za(p3,p6)*za(p7,p8)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*
     &    zba3(p1,p2,p7,p8,p5) )
      NAB(1,1) = NAB(1,1) + s78**(-1)*s34**(-1)*s56**(-1)*s3456**(-1)*
     & s278**(-1) * ( za(p2,p8)*za(p4,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,p6)
     &    *izb(p1,p7)*zba3(p1,p2,p7,p8,p3) )
      NAB(1,1) = NAB(1,1) + s34**(-1)*s56**(-1)*s3456**(-1)*s278**(-1)
     &  * ( za(p2,p8)*za(p3,p5)*zb(p1,p2)*zb(p1,p3)*zb(p4,p6)*izb(p1,p7
     &    )*izb(p1,p8)*izb(p2,p7)*zba3(p1,p2,p7,p8,p3) - za(p2,p8)*za(
     &    p3,p5)*zb(p1,p2)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*izb(p1,p8)*
     &    izb(p2,p7)*zba3(p1,p2,p7,p8,p3) + za(p2,p8)*za(p3,p5)*zb(p1,
     &    p2)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*
     &    zba3(p1,p2,p7,p8,p4) + za(p2,p8)*za(p3,p5)*zb(p1,p2)*zb(p1,p6
     &    )*zb(p4,p5)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*zba3(p1,p2,p7,p8
     &    ,p5) + za(p2,p8)*za(p3,p6)*zb(p1,p2)*zb(p1,p6)*zb(p4,p6)*izb(
     &    p1,p7)*izb(p1,p8)*izb(p2,p7)*zba3(p1,p2,p7,p8,p5) - za(p2,p8)
     &    *za(p4,p5)*zb(p1,p2)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p1,p8
     &    )*izb(p2,p7)*zba3(p1,p2,p7,p8,p3) + za(p3,p5)*za(p7,p8)*zb(p1
     &    ,p3)*zb(p4,p6)*izb(p1,p8)*izb(p2,p7)*zba3(p1,p2,p7,p8,p3) -
     &    za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p3,p6)*izb(p1,p8)*izb(p2,p7)
     &    *zba3(p1,p2,p7,p8,p3) + za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,
     &    p6)*izb(p1,p8)*izb(p2,p7)*zba3(p1,p2,p7,p8,p4) )
      NAB(1,1) = NAB(1,1) + s34**(-1)*s56**(-1)*s3456**(-1)*s278**(-1)
     &  * ( za(p3,p5)*za(p7,p8)*zb(p1,p6)*zb(p4,p5)*izb(p1,p8)*izb(p2,
     &    p7)*zba3(p1,p2,p7,p8,p5) + za(p3,p6)*za(p7,p8)*zb(p1,p6)*zb(
     &    p4,p6)*izb(p1,p8)*izb(p2,p7)*zba3(p1,p2,p7,p8,p5) - za(p4,p5)
     &    *za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*izb(p2,p7)*zba3(p1,
     &    p2,p7,p8,p3) )

      NAB(2,2)= + s78**(-1)*s34**(-1)*s56**(-1)*s3456**(-1)*s178**(-1)
     &  * (  - za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*iza(p2
     &    ,p8)*zba3(p4,p1,p7,p8,p2) + za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(
     &    p4,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p3,p1,p7,p8,p2) - za(p2,p3)*
     &    za(p3,p5)*zb(p1,p8)*zb(p3,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p4,p1
     &    ,p7,p8,p2) + za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*zb(p7,p8
     &    )*iza(p2,p7)*zba3(p3,p1,p7,p8,p2) - za(p2,p3)*za(p4,p5)*zb(p1
     &    ,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2) -
     &    za(p2,p3)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*zb(p7,p8)*iza(p2,p7)*
     &    zba3(p4,p1,p7,p8,p2) + za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p4,p6
     &    )*zb(p7,p8)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2) + za(p2,p4)*za(p3
     &    ,p5)*zb(p1,p8)*zb(p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p4,p1,p7,
     &    p8,p2) + za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*zb(p7,p8)*
     &    iza(p2,p8)*zba3(p6,p1,p7,p8,p2) + za(p2,p5)*za(p3,p5)*zb(p1,
     &    p8)*zb(p4,p5)*zb(p7,p8)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2) + za(
     &    p2,p5)*za(p3,p6)*zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*
     &    zba3(p6,p1,p7,p8,p2) )
      NAB(2,2) = NAB(2,2) + s78**(-1)*s34**(-1)*s56**(-1)*s3456**(-1)*
     & s178**(-1) * ( za(p2,p5)*za(p3,p6)*zb(p1,p8)*zb(p4,p6)*zb(p7,p8)
     &    *iza(p2,p7)*zba3(p6,p1,p7,p8,p2) )
      NAB(2,2) = NAB(2,2) + s34**(-1)*s56**(-1)*s3456**(-1)*s178**(-1)
     &  * ( za(p1,p2)*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*iza(p1,p8
     &    )*iza(p2,p7)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2) - za(p1,p2)*za(
     &    p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*iza(p1,p8)*iza(p2,p7)*
     &    iza(p2,p8)*zba3(p3,p1,p7,p8,p2) + za(p1,p2)*za(p2,p3)*za(p4,
     &    p5)*zb(p1,p7)*zb(p4,p6)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*
     &    zba3(p4,p1,p7,p8,p2) - za(p1,p2)*za(p2,p4)*za(p3,p5)*zb(p1,p7
     &    )*zb(p4,p6)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*zba3(p4,p1,p7,p8
     &    ,p2) - za(p1,p2)*za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(
     &    p1,p8)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,p7,p8,p2) - za(p1,p2)
     &    *za(p2,p5)*za(p3,p6)*zb(p1,p7)*zb(p4,p6)*iza(p1,p8)*iza(p2,p7
     &    )*iza(p2,p8)*zba3(p6,p1,p7,p8,p2) + za(p2,p3)*za(p3,p5)*zb(p3
     &    ,p6)*zb(p7,p8)*iza(p1,p8)*iza(p2,p7)*zba3(p4,p1,p7,p8,p2) -
     &    za(p2,p3)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*iza(p1,p8)*iza(p2,p7)
     &    *zba3(p3,p1,p7,p8,p2) + za(p2,p3)*za(p4,p5)*zb(p4,p6)*zb(p7,
     &    p8)*iza(p1,p8)*iza(p2,p7)*zba3(p4,p1,p7,p8,p2) )
      NAB(2,2) = NAB(2,2) + s34**(-1)*s56**(-1)*s3456**(-1)*s178**(-1)
     &  * (  - za(p2,p4)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*iza(p1,p8)*iza(
     &    p2,p7)*zba3(p4,p1,p7,p8,p2) - za(p2,p5)*za(p3,p5)*zb(p4,p5)*
     &    zb(p7,p8)*iza(p1,p8)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2) - za(p2,
     &    p5)*za(p3,p6)*zb(p4,p6)*zb(p7,p8)*iza(p1,p8)*iza(p2,p7)*zba3(
     &    p6,p1,p7,p8,p2) )

      NAB(1,2)= + s34**(-1)*s56**(-1)*s3456**(-1)*s278**(-1) * (  - za(
     &    p2,p7)*za(p3,p5)*zb(p1,p3)*zb(p2,p8)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*zba2(p8,p2,p7,p3) + za(p2,p7)*za(p3,p5)
     &    *zb(p1,p4)*zb(p2,p8)*zb(p3,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*zba2(p8,p2,p7,p3) - za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,
     &    p8)*zb(p4,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*zba2(p8,p2,p7,
     &    p4) - za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p2,p8)*zb(p4,p5)*iza(
     &    p7,p8)*izb(p2,p7)*izb(p7,p8)*zba2(p8,p2,p7,p5) - za(p2,p7)*
     &    za(p3,p6)*zb(p1,p6)*zb(p2,p8)*zb(p4,p6)*iza(p7,p8)*izb(p2,p7)
     &    *izb(p7,p8)*zba2(p8,p2,p7,p5) + za(p2,p7)*za(p4,p5)*zb(p1,p4)
     &    *zb(p2,p8)*zb(p4,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*zba2(p8
     &    ,p2,p7,p3) )
      NAB(1,2) = NAB(1,2) + s34**(-1)*s56**(-1)*s3456**(-1)*s178**(-1)
     &  * ( za(p1,p7)*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p1,p8
     &    )*iza(p7,p8)*izb(p7,p8)*zba2(p4,p1,p8,p7) - za(p1,p7)*za(p2,
     &    p3)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8)*zba2(p3,p1,p8,p7) + za(p1,p7)*za(p2,p3)*za(p4,p5)*zb(
     &    p1,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba2(p4,p1,
     &    p8,p7) - za(p1,p7)*za(p2,p4)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba2(p4,p1,p8,p7) - za(p1,p7
     &    )*za(p2,p5)*za(p3,p5)*zb(p1,p8)*zb(p4,p5)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8)*zba2(p6,p1,p8,p7) - za(p1,p7)*za(p2,p5)*za(p3,
     &    p6)*zb(p1,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    zba2(p6,p1,p8,p7) )
      NAB(1,2) = NAB(1,2) + s34**(-1)*s56**(-1)*s3456**(-1) * (  - za(
     &    p1,p7)*za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p2,p8)*zb(p3,p6)*iza(
     &    p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - za(p1,p7)*za(p2,p3)
     &    *za(p4,p5)*zb(p1,p4)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8
     &    )*izb(p2,p7)*izb(p7,p8) + za(p1,p7)*za(p2,p5)*za(p3,p5)*zb(p1
     &    ,p6)*zb(p2,p8)*zb(p4,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) + za(p1,p7)*za(p2,p5)*za(p3,p6)*zb(p1,p6)*zb(p2,p8
     &    )*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + za(
     &    p1,p7)*za(p3,p5)*za(p3,p7)*zb(p1,p4)*zb(p3,p6)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7) - za(p1,p7)*za(p3,p5)*za(p5,p7)*zb(p1,
     &    p6)*zb(p4,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7) - za(p1,p7)*
     &    za(p3,p5)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7
     &    )*izb(p7,p8)*zab2(p2,p3,p4,p1) - za(p1,p7)*za(p3,p5)*zb(p4,p6
     &    )*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*zab2(p7,p3,p4,p1) - za(p1,
     &    p7)*za(p3,p6)*za(p5,p7)*zb(p1,p6)*zb(p4,p6)*iza(p1,p8)*iza(p7
     &    ,p8)*izb(p2,p7) )
      NAB(1,2) = NAB(1,2) + s34**(-1)*s56**(-1)*s3456**(-1) * ( za(p1,
     &    p7)*za(p3,p7)*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*iza(p1,p8)*iza(p7
     &    ,p8)*izb(p2,p7) - za(p2,p3)*za(p3,p5)*zb(p2,p8)*zb(p3,p6)*zb(
     &    p4,p8)*iza(p1,p8)*izb(p2,p7)*izb(p7,p8) - za(p2,p3)*za(p4,p5)
     &    *zb(p2,p8)*zb(p4,p6)*zb(p4,p8)*iza(p1,p8)*izb(p2,p7)*izb(p7,
     &    p8) + za(p2,p5)*za(p3,p5)*zb(p2,p8)*zb(p4,p5)*zb(p6,p8)*iza(
     &    p1,p8)*izb(p2,p7)*izb(p7,p8) + za(p2,p5)*za(p3,p6)*zb(p2,p8)*
     &    zb(p4,p6)*zb(p6,p8)*iza(p1,p8)*izb(p2,p7)*izb(p7,p8) + za(p3,
     &    p5)*za(p3,p7)*zb(p3,p6)*zb(p4,p8)*iza(p1,p8)*izb(p2,p7) - za(
     &    p3,p5)*za(p5,p7)*zb(p4,p5)*zb(p6,p8)*iza(p1,p8)*izb(p2,p7) +
     &    za(p3,p5)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*izb(p2,p7)*izb(p7,p8
     &    )*zab2(p2,p3,p4,p8) + za(p3,p5)*zb(p4,p6)*iza(p1,p8)*izb(p2,
     &    p7)*zab2(p7,p3,p4,p8) - za(p3,p6)*za(p5,p7)*zb(p4,p6)*zb(p6,
     &    p8)*iza(p1,p8)*izb(p2,p7) + za(p3,p7)*za(p4,p5)*zb(p4,p6)*zb(
     &    p4,p8)*iza(p1,p8)*izb(p2,p7) )

      NAB(2,1)= + s34**(-1)*s56**(-1)*s3456**(-1)*s278**(-1) * (  - za(
     &    p2,p8)**2*za(p3,p5)*zb(p1,p3)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p7,p8)*zba2(p7,p2,p8,p3) + za(p2,p8)**2*za(p3,p5)*zb(p1,
     &    p4)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba2(p7,p2,p8,
     &    p3) - za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p7,p8)*zba2(p7,p2,p8,p4) - za(p2,p8)**2*za(p3,
     &    p5)*zb(p1,p6)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    zba2(p7,p2,p8,p5) - za(p2,p8)**2*za(p3,p6)*zb(p1,p6)*zb(p4,p6
     &    )*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba2(p7,p2,p8,p5) + za(p2,
     &    p8)**2*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*zba2(p7,p2,p8,p3) )
      NAB(2,1) = NAB(2,1) + s34**(-1)*s56**(-1)*s3456**(-1)*s178**(-1)
     &  * ( za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*zb(p3,p6)*iza(p7,p8)*izb(
     &    p1,p8)*izb(p7,p8)*zba2(p4,p1,p7,p8) - za(p2,p3)*za(p3,p5)*zb(
     &    p1,p7)**2*zb(p4,p6)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba2(p3,
     &    p1,p7,p8) + za(p2,p3)*za(p4,p5)*zb(p1,p7)**2*zb(p4,p6)*iza(p7
     &    ,p8)*izb(p1,p8)*izb(p7,p8)*zba2(p4,p1,p7,p8) - za(p2,p4)*za(
     &    p3,p5)*zb(p1,p7)**2*zb(p4,p6)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8
     &    )*zba2(p4,p1,p7,p8) - za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p4,
     &    p5)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba2(p6,p1,p7,p8) - za(
     &    p2,p5)*za(p3,p6)*zb(p1,p7)**2*zb(p4,p6)*iza(p7,p8)*izb(p1,p8)
     &    *izb(p7,p8)*zba2(p6,p1,p7,p8) )
      NAB(2,1) = NAB(2,1) + s34**(-1)*s56**(-1)*s3456**(-1) * (  - za(
     &    p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(p3,p6)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - za(p2,p3)*za(p2,p8)
     &    *za(p4,p5)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8
     &    )*izb(p1,p8)*izb(p7,p8) + za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1
     &    ,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) + za(p2,p5)*za(p2,p8)*za(p3,p6)*zb(p1,p6)*zb(p1,p7
     &    )*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - za(
     &    p2,p8)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p1,p8)*izb(p7,p8)*zab2(p2,p3,p4,p1) )

      return
      end
