      function TWZbbnr1(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      complex(dp):: TWZbbnr1

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
C     Author: R.K. Ellis Feb, 2013
C     written by program WZbbdiags.frm
C     These are the non-resonant diagrams
C     Calculation is performed for LH light-line (perforce because of W)
C     Calculation is performed for LH bbbar-line
C     Calculation is performed for LH Z-dcay line
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s56,s78,s3456,s456,s178,s278
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
C     end statement functions
      s56=s(p5,p6)
      s78=s(p7,p8)
      s178=s3(p1,p7,p8)
      s278=s3(p2,p7,p8)
      s456=s3(p4,p5,p6)
      s3456=s4(p3,p4,p5,p6)
      TWZbbnr1= + s3456**(-1)*s56**(-1)*s78**(-1)*s456**(-1)*s178**(-1)
     &  * (  - za(p1,p7)*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p1,p8)*zb(p4,
     &    p6) + za(p1,p7)*za(p2,p3)*za(p5,p6)*zb(p1,p6)*zb(p1,p8)*zb(p4
     &    ,p6) - za(p2,p3)*za(p4,p5)*za(p7,p8)*zb(p1,p8)*zb(p4,p6)*zb(
     &    p4,p8) + za(p2,p3)*za(p5,p6)*za(p7,p8)*zb(p1,p8)*zb(p4,p6)*
     &    zb(p6,p8) )
      TWZbbnr1 = TWZbbnr1 + s3456**(-1)*s56**(-1)*s78**(-1)*s456**(-1)*
     & s278**(-1) * ( za(p2,p3)*za(p2,p7)*za(p4,p5)*zb(p1,p4)*zb(p2,p8)
     &    *zb(p4,p6) - za(p2,p3)*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p2,p8
     &    )*zb(p4,p6) - za(p2,p7)*za(p3,p7)*za(p4,p5)*zb(p1,p4)*zb(p4,
     &    p6)*zb(p7,p8) + za(p2,p7)*za(p3,p7)*za(p5,p6)*zb(p1,p6)*zb(p4
     &    ,p6)*zb(p7,p8) )

      return
      end
