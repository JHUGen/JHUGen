      subroutine c7tree(p1,p2,p3,p4,p5,p6,p7,f)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      complex(dp):: iza,izb,f(4,2,2)
C     first index controls the choice of diagrams
C     second index of fa is helicity of 1-2 line,
C     third is gluon helicity
      real(dp):: s12,s34,s123,s124,s125,s126,s345,s346,s567
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
      s123=s(p1,p2)+s(p2,p3)+s(p1,p3)
      s124=s(p1,p2)+s(p2,p4)+s(p1,p4)
      s125=s(p1,p2)+s(p1,p5)+s(p2,p5)
      s126=s(p1,p2)+s(p1,p6)+s(p2,p6)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      s346=s(p3,p4)+s(p3,p6)+s(p4,p6)
      s34=s(p3,p4)
      s12=s(p1,p2)
      f(1,1,1)= + s567**(-1)*s12**(-1) * (  - za(p1,p2)*za(p3,p5)*zb(p1
     &    ,p4)*zb(p1,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s124**(-1) +
     &    za(p1,p2)*za(p3,p7)*zb(p1,p4)*zb(p1,p6)*izb(p5,p7)*s124**(-1)
     &     + za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*zb(p5,p6)*izb(p5,
     &    p7)*izb(p6,p7)*s124**(-1) - za(p2,p4)*za(p3,p7)*zb(p1,p4)*zb(
     &    p4,p6)*izb(p5,p7)*s124**(-1) )

      f(1,1,2)= + s567**(-1)*s12**(-1) * ( za(p1,p2)*za(p3,p5)*za(p5,p6
     &    )*zb(p1,p4)*zb(p1,p6)*iza(p5,p7)*iza(p6,p7)*s124**(-1) + za(
     &    p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p6,p7)*s124**(-1) -
     &    za(p2,p4)*za(p3,p5)*za(p5,p6)*zb(p1,p4)*zb(p4,p6)*iza(p5,p7)*
     &    iza(p6,p7)*s124**(-1) - za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p4,
     &    p7)*iza(p6,p7)*s124**(-1) )

      f(1,2,1)= + s567**(-1)*s12**(-1) * ( za(p1,p2)*za(p3,p5)*zb(p2,p4
     &    )*zb(p2,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s124**(-1) - za(
     &    p1,p2)*za(p3,p7)*zb(p2,p4)*zb(p2,p6)*izb(p5,p7)*s124**(-1) +
     &    za(p1,p4)*za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*
     &    izb(p6,p7)*s124**(-1) - za(p1,p4)*za(p3,p7)*zb(p2,p4)*zb(p4,
     &    p6)*izb(p5,p7)*s124**(-1) )

      f(1,2,2)= + s567**(-1)*s12**(-1) * (  - za(p1,p2)*za(p3,p5)*za(p5
     &    ,p6)*zb(p2,p4)*zb(p2,p6)*iza(p5,p7)*iza(p6,p7)*s124**(-1) -
     &    za(p1,p2)*za(p3,p5)*zb(p2,p4)*zb(p2,p7)*iza(p6,p7)*s124**(-1)
     &     - za(p1,p4)*za(p3,p5)*za(p5,p6)*zb(p2,p4)*zb(p4,p6)*iza(p5,
     &    p7)*iza(p6,p7)*s124**(-1) - za(p1,p4)*za(p3,p5)*zb(p2,p4)*zb(
     &    p4,p7)*iza(p6,p7)*s124**(-1) )

      f(2,1,1)= + s567**(-1)*s12**(-1) * (  - za(p2,p3)*za(p2,p5)*zb(p1
     &    ,p2)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s123**(-1) +
     &    za(p2,p3)*za(p2,p7)*zb(p1,p2)*zb(p4,p6)*izb(p5,p7)*s123**(-1)
     &     - za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p4,p6)*zb(p5,p6)*izb(p5,
     &    p7)*izb(p6,p7)*s123**(-1) + za(p2,p3)*za(p3,p7)*zb(p1,p3)*zb(
     &    p4,p6)*izb(p5,p7)*s123**(-1) )

      f(2,1,2)= + s567**(-1)*s12**(-1) * ( za(p2,p3)*za(p2,p5)*za(p5,p6
     &    )*zb(p1,p2)*zb(p4,p6)*iza(p5,p7)*iza(p6,p7)*s123**(-1) + za(
     &    p2,p3)*za(p2,p5)*zb(p1,p2)*zb(p4,p7)*iza(p6,p7)*s123**(-1) +
     &    za(p2,p3)*za(p3,p5)*za(p5,p6)*zb(p1,p3)*zb(p4,p6)*iza(p5,p7)*
     &    iza(p6,p7)*s123**(-1) + za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p4,
     &    p7)*iza(p6,p7)*s123**(-1) )

      f(2,2,1)= + s567**(-1)*s12**(-1) * ( za(p1,p3)*za(p1,p5)*zb(p1,p2
     &    )*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s123**(-1) - za(
     &    p1,p3)*za(p1,p7)*zb(p1,p2)*zb(p4,p6)*izb(p5,p7)*s123**(-1) -
     &    za(p1,p3)*za(p3,p5)*zb(p2,p3)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*
     &    izb(p6,p7)*s123**(-1) + za(p1,p3)*za(p3,p7)*zb(p2,p3)*zb(p4,
     &    p6)*izb(p5,p7)*s123**(-1) )

      f(2,2,2)= + s567**(-1)*s12**(-1) * (  - za(p1,p3)*za(p1,p5)*za(p5
     &    ,p6)*zb(p1,p2)*zb(p4,p6)*iza(p5,p7)*iza(p6,p7)*s123**(-1) -
     &    za(p1,p3)*za(p1,p5)*zb(p1,p2)*zb(p4,p7)*iza(p6,p7)*s123**(-1)
     &     + za(p1,p3)*za(p3,p5)*za(p5,p6)*zb(p2,p3)*zb(p4,p6)*iza(p5,
     &    p7)*iza(p6,p7)*s123**(-1) + za(p1,p3)*za(p3,p5)*zb(p2,p3)*zb(
     &    p4,p7)*iza(p6,p7)*s123**(-1) )

      f(3,1,1)= + s34**(-1)*s12**(-1) * (  - za(p2,p5)*za(p2,p7)*za(p3,
     &    p4)*zb(p1,p2)*zb(p4,p6)**2*izb(p6,p7)*s346**(-1)*s125**(-1)
     &     - za(p2,p5)*za(p3,p4)*za(p5,p7)*zb(p1,p5)*zb(p4,p6)**2*izb(
     &    p6,p7)*s346**(-1)*s125**(-1) + za(p2,p5)*za(p3,p4)*zb(p1,p4)*
     &    zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s346**(-1) + za(p2,
     &    p5)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6
     &    ,p7)*s346**(-1) - za(p2,p7)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)*
     &    izb(p5,p7)*s346**(-1) - za(p2,p7)*za(p3,p6)*zb(p1,p6)*zb(p4,
     &    p6)*izb(p5,p7)*s346**(-1) )

      f(3,1,2)= + s34**(-1)*s12**(-1) * (  - za(p2,p3)*za(p2,p5)*za(p5,
     &    p6)*zb(p1,p2)*zb(p4,p6)*iza(p5,p7)*iza(p6,p7)*s125**(-1) -
     &    za(p2,p3)*za(p2,p5)*zb(p1,p2)*zb(p4,p7)*iza(p6,p7)*s125**(-1)
     &     + za(p2,p5)**2*za(p3,p4)*zb(p1,p2)*zb(p4,p6)*zb(p4,p7)*iza(
     &    p5,p7)*s346**(-1)*s125**(-1) + za(p2,p5)**2*za(p3,p6)*zb(p1,
     &    p2)*zb(p4,p6)*zb(p6,p7)*iza(p5,p7)*s346**(-1)*s125**(-1) +
     &    za(p2,p5)*za(p3,p5)*za(p5,p6)*zb(p1,p5)*zb(p4,p6)*iza(p5,p7)*
     &    iza(p6,p7)*s125**(-1) + za(p2,p5)*za(p3,p5)*zb(p1,p5)*zb(p4,
     &    p7)*iza(p6,p7)*s125**(-1) )

      f(3,2,1)= + s34**(-1)*s12**(-1) * ( za(p1,p5)*za(p1,p7)*za(p3,p4)
     &    *zb(p1,p2)*zb(p4,p6)**2*izb(p6,p7)*s346**(-1)*s125**(-1) -
     &    za(p1,p5)*za(p3,p4)*za(p5,p7)*zb(p2,p5)*zb(p4,p6)**2*izb(p6,
     &    p7)*s346**(-1)*s125**(-1) + za(p1,p5)*za(p3,p4)*zb(p2,p4)*zb(
     &    p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s346**(-1) + za(p1,p5)
     &    *za(p3,p6)*zb(p2,p6)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7
     &    )*s346**(-1) - za(p1,p7)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*izb(p5
     &    ,p7)*s346**(-1) - za(p1,p7)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    izb(p5,p7)*s346**(-1) )

      f(3,2,2)= + s34**(-1)*s12**(-1) * ( za(p1,p3)*za(p1,p5)*za(p5,p6)
     &    *zb(p1,p2)*zb(p4,p6)*iza(p5,p7)*iza(p6,p7)*s125**(-1) + za(p1
     &    ,p3)*za(p1,p5)*zb(p1,p2)*zb(p4,p7)*iza(p6,p7)*s125**(-1) -
     &    za(p1,p5)**2*za(p3,p4)*zb(p1,p2)*zb(p4,p6)*zb(p4,p7)*iza(p5,
     &    p7)*s346**(-1)*s125**(-1) - za(p1,p5)**2*za(p3,p6)*zb(p1,p2)*
     &    zb(p4,p6)*zb(p6,p7)*iza(p5,p7)*s346**(-1)*s125**(-1) + za(p1,
     &    p5)*za(p3,p5)*za(p5,p6)*zb(p2,p5)*zb(p4,p6)*iza(p5,p7)*iza(p6
     &    ,p7)*s125**(-1) + za(p1,p5)*za(p3,p5)*zb(p2,p5)*zb(p4,p7)*
     &    iza(p6,p7)*s125**(-1) )

      f(4,1,1)= + s34**(-1)*s12**(-1) * (  - za(p1,p2)*za(p3,p5)*za(p3,
     &    p7)*zb(p1,p6)**2*zb(p3,p4)*izb(p6,p7)*s126**(-1)*s345**(-1)
     &     + za(p1,p2)*za(p3,p5)*za(p5,p7)*zb(p1,p6)**2*zb(p4,p5)*izb(
     &    p6,p7)*s126**(-1)*s345**(-1) + za(p1,p2)*za(p3,p5)*zb(p1,p4)*
     &    zb(p1,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s126**(-1) - za(p1,
     &    p2)*za(p3,p7)*zb(p1,p4)*zb(p1,p6)*izb(p5,p7)*s126**(-1) + za(
     &    p2,p6)*za(p3,p5)*zb(p1,p6)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*
     &    izb(p6,p7)*s126**(-1) - za(p2,p6)*za(p3,p7)*zb(p1,p6)*zb(p4,
     &    p6)*izb(p5,p7)*s126**(-1) )

      f(4,1,2)= + s34**(-1)*s12**(-1) * ( za(p1,p2)*za(p3,p5)**2*zb(p1,
     &    p6)*zb(p1,p7)*zb(p3,p4)*iza(p5,p7)*s126**(-1)*s345**(-1) -
     &    za(p2,p3)*za(p3,p5)*za(p5,p6)*zb(p1,p6)*zb(p3,p4)*iza(p5,p7)*
     &    iza(p6,p7)*s345**(-1) - za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,
     &    p4)*iza(p6,p7)*s345**(-1) + za(p2,p5)*za(p3,p5)*za(p5,p6)*zb(
     &    p1,p6)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)*s345**(-1) + za(p2,p5)
     &    *za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p6,p7)*s345**(-1) - za(p2,
     &    p6)*za(p3,p5)**2*zb(p1,p6)*zb(p3,p4)*zb(p6,p7)*iza(p5,p7)*
     &    s126**(-1)*s345**(-1) )

      f(4,2,1)= + s34**(-1)*s12**(-1) * ( za(p1,p2)*za(p3,p5)*za(p3,p7)
     &    *zb(p2,p6)**2*zb(p3,p4)*izb(p6,p7)*s126**(-1)*s345**(-1) -
     &    za(p1,p2)*za(p3,p5)*za(p5,p7)*zb(p2,p6)**2*zb(p4,p5)*izb(p6,
     &    p7)*s126**(-1)*s345**(-1) - za(p1,p2)*za(p3,p5)*zb(p2,p4)*zb(
     &    p2,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s126**(-1) + za(p1,p2)
     &    *za(p3,p7)*zb(p2,p4)*zb(p2,p6)*izb(p5,p7)*s126**(-1) + za(p1,
     &    p6)*za(p3,p5)*zb(p2,p6)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6
     &    ,p7)*s126**(-1) - za(p1,p6)*za(p3,p7)*zb(p2,p6)*zb(p4,p6)*
     &    izb(p5,p7)*s126**(-1) )

      f(4,2,2)= + s34**(-1)*s12**(-1) * (  - za(p1,p2)*za(p3,p5)**2*zb(
     &    p2,p6)*zb(p2,p7)*zb(p3,p4)*iza(p5,p7)*s126**(-1)*s345**(-1)
     &     - za(p1,p3)*za(p3,p5)*za(p5,p6)*zb(p2,p6)*zb(p3,p4)*iza(p5,
     &    p7)*iza(p6,p7)*s345**(-1) - za(p1,p3)*za(p3,p5)*zb(p2,p7)*zb(
     &    p3,p4)*iza(p6,p7)*s345**(-1) + za(p1,p5)*za(p3,p5)*za(p5,p6)*
     &    zb(p2,p6)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)*s345**(-1) + za(p1,
     &    p5)*za(p3,p5)*zb(p2,p7)*zb(p4,p5)*iza(p6,p7)*s345**(-1) - za(
     &    p1,p6)*za(p3,p5)**2*zb(p2,p6)*zb(p3,p4)*zb(p6,p7)*iza(p5,p7)*
     &    s126**(-1)*s345**(-1) )

      return
      end
