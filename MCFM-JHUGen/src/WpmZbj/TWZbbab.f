      function TWZbbab(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      complex(dp):: TWZbbab
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
C     Author: R.K. Ellis Feb, 2013
C     written by program WZbbdiags.frm
C     These are the Abelian diagrams with both the W and Z coming
C     off light line, only three diagrams with Z emitted before W
C     Calculation is performed for LH light-line (perforce because of W)
C     Calculation is performed for LH bbbar-line
C     Calculation is performed for LH Z-dcay line
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s34,s56,s78,s134,s256,s278,s178
      complex(dp):: d1,d2,d3
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
C     end statement functions
      s34=s(p3,p4)
      s56=s(p5,p6)
      s78=s(p7,p8)
      s134=s3(p1,p3,p4)
      s256=s3(p2,p5,p6)
      s278=s3(p2,p7,p8)
      s178=s3(p1,p7,p8)
      d1= + s34**(-1)*s56**(-1)*s78**(-1)*s256**(-1)*s178**(-1) * ( za(
     &    p1,p7)*za(p2,p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p8)*zb(p2,p6) - 
     &    za(p1,p7)*za(p2,p5)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*zb(p5,p6)
     &     + za(p2,p3)*za(p2,p5)*za(p7,p8)*zb(p1,p8)*zb(p2,p6)*zb(p4,p8
     &    ) - za(p2,p5)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*zb(p4,p8)*zb(p5,
     &    p6) )

      d2= + s34**(-1)*s56**(-1)*s78**(-1)*s134**(-1)*s256**(-1) * ( za(
     &    p1,p3)*za(p2,p5)*za(p2,p7)*zb(p1,p4)*zb(p1,p8)*zb(p2,p6) + 
     &    za(p1,p3)*za(p2,p5)*za(p5,p7)*zb(p1,p4)*zb(p1,p8)*zb(p5,p6)
     &     - za(p2,p5)*za(p2,p7)*za(p3,p4)*zb(p1,p4)*zb(p2,p6)*zb(p4,p8
     &    ) - za(p2,p5)*za(p3,p4)*za(p5,p7)*zb(p1,p4)*zb(p4,p8)*zb(p5,
     &    p6) )

      d3= + s34**(-1)*s56**(-1)*s78**(-1)*s134**(-1)*s278**(-1) * ( za(
     &    p1,p3)*za(p2,p5)*za(p2,p7)*zb(p1,p4)*zb(p1,p6)*zb(p2,p8) - 
     &    za(p1,p3)*za(p2,p7)*za(p5,p7)*zb(p1,p4)*zb(p1,p6)*zb(p7,p8)
     &     - za(p2,p5)*za(p2,p7)*za(p3,p4)*zb(p1,p4)*zb(p2,p8)*zb(p4,p6
     &    ) + za(p2,p7)*za(p3,p4)*za(p5,p7)*zb(p1,p4)*zb(p4,p6)*zb(p7,
     &    p8) )

      TWZbbab=d1+d2+d3
      return
      end
