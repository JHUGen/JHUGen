      subroutine TWZggAB(p1,p2,p3,p4,p5,p6,p7,p8,AB)
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
C     These are the Abelian diagrams with both the W and Z coming
C     off light line, only three diagrams with W emitted before Z
C     as we proceed from p1 to p2
C     Calculation is performed for LH light-line (perforce because of W)
C     Calculation is performed for LH Z-dcay line
C     The two indices of AB are the gluon helicities
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s34,s56,s78,s134,s256,s278,s178,s2567
      complex(dp):: zba2,zba3,AB(2,2),iza,izb
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
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
      s134=s3(p1,p3,p4)
      s256=s3(p2,p5,p6)
      s278=s3(p2,p7,p8)
      s178=s3(p1,p7,p8)
      s2567=s4(p2,p5,p6,p7)
      AB(1,1)= + s78**(-1)*s34**(-1)*s56**(-1)*s134**(-1)*s278**(-1)
     &  * ( za(p2,p7)*za(p7,p8)*zb(p1,p4)*izb(p1,p8)*zba2(p6,p1,p4,p3)*
     &    zba3(p1,p2,p7,p8,p5) + za(p2,p8)*za(p7,p8)*zb(p1,p4)*izb(p1,
     &    p7)*zba2(p6,p1,p4,p3)*zba3(p1,p2,p7,p8,p5) )
      AB(1,1) = AB(1,1) + s78**(-1)*s34**(-1)*s56**(-1)*s134**(-1)*
     & s256**(-1) * ( za(p2,p5)*za(p7,p8)*zb(p1,p4)*izb(p1,p7)*zba2(p1,
     &    p1,p4,p3)*zba2(p6,p2,p5,p8) + za(p2,p5)*za(p7,p8)*zb(p1,p4)*
     &    izb(p1,p8)*zba2(p1,p1,p4,p3)*zba2(p6,p2,p5,p7) )
      AB(1,1) = AB(1,1) + s34**(-1)*s56**(-1)*s134**(-1)*s278**(-1)
     &  * (  - za(p2,p8)*zb(p1,p2)*zb(p1,p4)*izb(p1,p7)*izb(p1,p8)*izb(
     &    p2,p7)*zba2(p6,p1,p4,p3)*zba3(p1,p2,p7,p8,p5) - za(p7,p8)*zb(
     &    p1,p4)*izb(p1,p8)*izb(p2,p7)*zba2(p6,p1,p4,p3)*zba3(p1,p2,p7,
     &    p8,p5) )
      AB(1,1) = AB(1,1) + s34**(-1)*s56**(-1)*s134**(-1)*s256**(-1)*
     & s2567**(-1) * (  - za(p2,p5)*zb(p1,p4)*izb(p1,p7)*izb(p1,p8)*
     &    zba2(p1,p1,p4,p3)*zba2(p1,p3,p4,p8)*zba2(p6,p2,p5,p7) )
      AB(1,1) = AB(1,1) + s34**(-1)*s56**(-1)*s134**(-1)*s2567**(-1)
     &  * (  - za(p2,p5)*zb(p1,p2)*zb(p1,p4)*izb(p1,p7)*izb(p1,p8)*izb(
     &    p2,p7)*zba2(p1,p1,p4,p3)*zba3(p6,p2,p5,p7,p8) + za(p5,p7)*zb(
     &    p1,p4)*izb(p1,p8)*izb(p2,p7)*zba2(p1,p1,p4,p3)*zba3(p6,p2,p5,
     &    p7,p8) )

      AB(2,2)= + s78**(-1)*s34**(-1)*s56**(-1)*s134**(-1)*s256**(-1)
     &  * ( za(p2,p5)*zb(p1,p4)*zb(p7,p8)*iza(p2,p7)*zba2(p6,p5,p6,p2)*
     &    zba2(p8,p1,p4,p3) + za(p2,p5)*zb(p1,p4)*zb(p7,p8)*iza(p2,p8)*
     &    zba2(p6,p5,p6,p2)*zba2(p7,p1,p4,p3) )
      AB(2,2) = AB(2,2) + s78**(-1)*s34**(-1)*s56**(-1)*s178**(-1)*
     & s256**(-1) * ( za(p2,p5)*zb(p1,p7)*zb(p7,p8)*iza(p2,p8)*zba2(p6,
     &    p2,p5,p3)*zba3(p4,p1,p7,p8,p2) + za(p2,p5)*zb(p1,p8)*zb(p7,p8
     &    )*iza(p2,p7)*zba2(p6,p2,p5,p3)*zba3(p4,p1,p7,p8,p2) )
      AB(2,2) = AB(2,2) + s34**(-1)*s56**(-1)*s134**(-1)*s256**(-1)*
     & s2567**(-1) * ( za(p2,p5)*zb(p1,p4)*iza(p2,p7)*iza(p2,p8)*zba2(
     &    p6,p5,p6,p2)*zba2(p7,p5,p6,p2)*zba2(p8,p1,p4,p3) )
      AB(2,2) = AB(2,2) + s34**(-1)*s56**(-1)*s178**(-1)*s256**(-1)
     &  * (  - za(p1,p2)*za(p2,p5)*zb(p1,p7)*iza(p1,p8)*iza(p2,p7)*iza(
     &    p2,p8)*zba2(p6,p2,p5,p3)*zba3(p4,p1,p7,p8,p2) - za(p2,p5)*zb(
     &    p7,p8)*iza(p1,p8)*iza(p2,p7)*zba2(p6,p2,p5,p3)*zba3(p4,p1,p7,
     &    p8,p2) )
      AB(2,2) = AB(2,2) + s34**(-1)*s56**(-1)*s256**(-1)*s2567**(-1)
     &  * ( za(p1,p2)*za(p2,p5)*zb(p1,p4)*iza(p1,p8)*iza(p2,p7)*iza(p2,
     &    p8)*zba2(p6,p5,p6,p2)*zba3(p7,p2,p5,p6,p3) + za(p2,p5)*zb(p4,
     &    p8)*iza(p1,p8)*iza(p2,p7)*zba2(p6,p5,p6,p2)*zba3(p7,p2,p5,p6,
     &    p3) )

      AB(1,2)= + s34**(-1)*s56**(-1)*s134**(-1)*s278**(-1) * ( za(p2,p7
     &    )*zb(p1,p4)*zb(p2,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*zba2(
     &    p6,p1,p4,p3)*zba2(p8,p2,p7,p5) )
      AB(1,2) = AB(1,2) + s34**(-1)*s56**(-1)*s134**(-1)*s256**(-1)*
     & s2567**(-1) * ( za(p2,p5)*zb(p1,p4)*iza(p7,p8)*izb(p7,p8)*zba2(
     &    p6,p2,p5,p7)*zba2(p8,p1,p4,p3)*zba3(p8,p2,p5,p6,p7) )
      AB(1,2) = AB(1,2) + s34**(-1)*s56**(-1)*s134**(-1)*s2567**(-1)
     &  * ( za(p2,p5)*zb(p1,p4)*zb(p2,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*zba2(p6,p2,p5,p7)*zba2(p8,p1,p4,p3) - za(p5,p7)*zb(p1,p4)
     &    *iza(p7,p8)*izb(p2,p7)*zba2(p6,p2,p5,p7)*zba2(p8,p1,p4,p3) )
      AB(1,2) = AB(1,2) + s34**(-1)*s56**(-1)*s178**(-1)*s256**(-1)
     &  * (  - za(p1,p7)*za(p2,p5)*zb(p1,p8)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8)*zba2(p4,p1,p8,p7)*zba2(p6,p2,p5,p3) )
      AB(1,2) = AB(1,2) + s34**(-1)*s56**(-1)*s256**(-1)*s2567**(-1)
     &  * (  - za(p1,p7)*za(p2,p5)*zb(p1,p4)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8)*zba2(p6,p2,p5,p7)*zba2(p8,p1,p4,p3) - za(p2,p5)*zb(p4,
     &    p8)*iza(p1,p8)*izb(p7,p8)*zba2(p6,p2,p5,p7)*zba2(p8,p1,p4,p3)
     &     )
      AB(1,2) = AB(1,2) + s34**(-1)*s56**(-1)*s2567**(-1) * ( za(p1,p7)
     &    *za(p2,p5)*zb(p1,p4)*zb(p2,p8)*iza(p1,p8)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8)*zba3(p6,p2,p5,p7,p3) - za(p1,p7)*za(p5,p7)*zb(
     &    p1,p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*zba3(p6,p2,p5,p7,p3)
     &     + za(p2,p5)*zb(p2,p8)*zb(p4,p8)*iza(p1,p8)*izb(p2,p7)*izb(p7
     &    ,p8)*zba3(p6,p2,p5,p7,p3) - za(p5,p7)*zb(p4,p8)*iza(p1,p8)*
     &    izb(p2,p7)*zba3(p6,p2,p5,p7,p3) )

      AB(2,1)= + s34**(-1)*s56**(-1)*s134**(-1)*s278**(-1) * ( za(p2,p8
     &    )**2*zb(p1,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba2(p6,p1,p4
     &    ,p3)*zba2(p7,p2,p8,p5) )
      AB(2,1) = AB(2,1) + s34**(-1)*s56**(-1)*s134**(-1)*s256**(-1)*
     & s2567**(-1) * ( za(p2,p5)*zb(p1,p4)*iza(p7,p8)*izb(p7,p8)*zba2(
     &    p6,p2,p5,p8)*zba2(p7,p1,p4,p3)*zba3(p7,p2,p5,p6,p8) )
      AB(2,1) = AB(2,1) + s34**(-1)*s56**(-1)*s134**(-1)*s2567**(-1)
     &  * ( za(p2,p5)*za(p2,p8)*zb(p1,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*zba2(p7,p1,p4,p3)*zba3(p6,p2,p5,p7,p8) )
      AB(2,1) = AB(2,1) + s34**(-1)*s56**(-1)*s178**(-1)*s256**(-1)
     &  * (  - za(p2,p5)*zb(p1,p7)**2*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    zba2(p4,p1,p7,p8)*zba2(p6,p2,p5,p3) )
      AB(2,1) = AB(2,1) + s34**(-1)*s56**(-1)*s256**(-1)*s2567**(-1)
     &  * ( za(p2,p5)*zb(p1,p4)*zb(p1,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*zba2(p6,p2,p5,p8)*zba3(p7,p2,p5,p6,p3) )
      AB(2,1) = AB(2,1) + s34**(-1)*s56**(-1)*s2567**(-1) * ( za(p2,p5)
     &    *za(p2,p8)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*iza(p7,p8)*izb(p1,
     &    p8)*izb(p7,p8)*zba3(p6,p2,p5,p7,p3) )

      return
      end
