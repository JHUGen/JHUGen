c--- File written by FORM program dksrwz.frm on Thu Jul 26 20:12:12 CDT 2012
      subroutine dksrwz(p1,p2,p3,p4,p5,p6,p7,f)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'nwz.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: s12,s123,s124,s567
      complex(dp):: iza,izb,cc11,cc12,f(4,2,2)
C     first index of f controls the choice of diagrams
C     second index of f  is helicity of 5-6 line,
C     third is gluon helicity
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
      s123=s(p1,p2)+s(p2,p3)+s(p1,p3)
      s124=s(p1,p2)+s(p2,p4)+s(p1,p4)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)
      s12=s(p1,p2)
      f(1,1,1)= + s12**(-1)*s567**(-1) * (  - za(p1,p2)*za(p3,p5)*zb(p1
     &    ,p4)*zb(p1,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s124**(-1) +
     &    za(p1,p2)*za(p3,p7)*zb(p1,p4)*zb(p1,p6)*izb(p5,p7)*s124**(-1)
     &     + za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*zb(p5,p6)*izb(p5,
     &    p7)*izb(p6,p7)*s124**(-1) - za(p2,p4)*za(p3,p7)*zb(p1,p4)*zb(
     &    p4,p6)*izb(p5,p7)*s124**(-1) )

      f(1,1,2)= + s12**(-1)*s567**(-1) * ( za(p1,p2)*za(p3,p5)*za(p5,p6
     &    )*zb(p1,p4)*zb(p1,p6)*iza(p5,p7)*iza(p6,p7)*s124**(-1) + za(
     &    p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p6,p7)*s124**(-1) -
     &    za(p2,p4)*za(p3,p5)*za(p5,p6)*zb(p1,p4)*zb(p4,p6)*iza(p5,p7)*
     &    iza(p6,p7)*s124**(-1) - za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p4,
     &    p7)*iza(p6,p7)*s124**(-1) )

      f(1,2,1)= + s12**(-1)*s567**(-1) * (  - za(p1,p2)*za(p3,p6)*zb(p1
     &    ,p4)*zb(p1,p5)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s124**(-1) -
     &    za(p1,p2)*za(p3,p7)*zb(p1,p4)*zb(p1,p5)*izb(p6,p7)*s124**(-1)
     &     + za(p2,p4)*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*zb(p5,p6)*izb(p5,
     &    p7)*izb(p6,p7)*s124**(-1) + za(p2,p4)*za(p3,p7)*zb(p1,p4)*zb(
     &    p4,p5)*izb(p6,p7)*s124**(-1) )

      f(1,2,2)= + s12**(-1)*s567**(-1) * ( za(p1,p2)*za(p3,p6)*za(p5,p6
     &    )*zb(p1,p4)*zb(p1,p5)*iza(p5,p7)*iza(p6,p7)*s124**(-1) - za(
     &    p1,p2)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*iza(p5,p7)*s124**(-1) -
     &    za(p2,p4)*za(p3,p6)*za(p5,p6)*zb(p1,p4)*zb(p4,p5)*iza(p5,p7)*
     &    iza(p6,p7)*s124**(-1) + za(p2,p4)*za(p3,p6)*zb(p1,p4)*zb(p4,
     &    p7)*iza(p5,p7)*s124**(-1) )

      f(2,1,1)= + s12**(-1)*s567**(-1) * (  - za(p2,p3)*za(p2,p5)*zb(p1
     &    ,p2)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s123**(-1) +
     &    za(p2,p3)*za(p2,p7)*zb(p1,p2)*zb(p4,p6)*izb(p5,p7)*s123**(-1)
     &     - za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p4,p6)*zb(p5,p6)*izb(p5,
     &    p7)*izb(p6,p7)*s123**(-1) + za(p2,p3)*za(p3,p7)*zb(p1,p3)*zb(
     &    p4,p6)*izb(p5,p7)*s123**(-1) )

      f(2,1,2)= + s12**(-1)*s567**(-1) * ( za(p2,p3)*za(p2,p5)*za(p5,p6
     &    )*zb(p1,p2)*zb(p4,p6)*iza(p5,p7)*iza(p6,p7)*s123**(-1) + za(
     &    p2,p3)*za(p2,p5)*zb(p1,p2)*zb(p4,p7)*iza(p6,p7)*s123**(-1) +
     &    za(p2,p3)*za(p3,p5)*za(p5,p6)*zb(p1,p3)*zb(p4,p6)*iza(p5,p7)*
     &    iza(p6,p7)*s123**(-1) + za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p4,
     &    p7)*iza(p6,p7)*s123**(-1) )

      f(2,2,1)= + s12**(-1)*s567**(-1) * (  - za(p2,p3)*za(p2,p6)*zb(p1
     &    ,p2)*zb(p4,p5)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s123**(-1) -
     &    za(p2,p3)*za(p2,p7)*zb(p1,p2)*zb(p4,p5)*izb(p6,p7)*s123**(-1)
     &     - za(p2,p3)*za(p3,p6)*zb(p1,p3)*zb(p4,p5)*zb(p5,p6)*izb(p5,
     &    p7)*izb(p6,p7)*s123**(-1) - za(p2,p3)*za(p3,p7)*zb(p1,p3)*zb(
     &    p4,p5)*izb(p6,p7)*s123**(-1) )

      f(2,2,2)= + s12**(-1)*s567**(-1) * ( za(p2,p3)*za(p2,p6)*za(p5,p6
     &    )*zb(p1,p2)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)*s123**(-1) - za(
     &    p2,p3)*za(p2,p6)*zb(p1,p2)*zb(p4,p7)*iza(p5,p7)*s123**(-1) +
     &    za(p2,p3)*za(p3,p6)*za(p5,p6)*zb(p1,p3)*zb(p4,p5)*iza(p5,p7)*
     &    iza(p6,p7)*s123**(-1) - za(p2,p3)*za(p3,p6)*zb(p1,p3)*zb(p4,
     &    p7)*iza(p5,p7)*s123**(-1) )

      if (nwz == 1) then
      call dksrwzc(p1,p2,p3,p4,p5,p6,p7,cc11,cc12)
      f(3,1,1)=cc11
      f(3,1,2)=cc12
      f(3,2,1)=czip
      f(3,2,2)=czip
      call dksrwzc(p4,p3,p2,p1,p5,p6,p7,cc11,cc12)
      f(4,1,1)=cc11
      f(4,1,2)=cc12
      f(4,2,1)=czip
      f(4,2,2)=czip
      elseif (nwz == -1) then
      call dksrwzc(p4,p3,p2,p1,p5,p6,p7,cc11,cc12)
      f(3,1,1)=cc11
      f(3,1,2)=cc12
      f(3,2,1)=czip
      f(3,2,2)=czip
      call dksrwzc(p1,p2,p3,p4,p5,p6,p7,cc11,cc12)
      f(4,1,1)=cc11
      f(4,1,2)=cc12
      f(4,2,1)=czip
      f(4,2,2)=czip
      endif
      return
      end
