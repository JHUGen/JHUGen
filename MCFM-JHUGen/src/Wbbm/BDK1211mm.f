      function BDK1211mm(k1,k2,k3,k4,k5,k6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: BDK1211mm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'momwbbm.f'
      include 'Wbbmlabels.f'
      integer:: k1,k2,k3,k4,k5,k6,nu
      real(dp):: t,p2(4),p3(4),p23(4),p123(4),p234(4),p14(4),
     & p1234(4),s23,s123,s234,s1234,s14,s56
      complex(dp):: zab,bubrat
      real(dp):: IDELTA3,DELTA3,ddelta
      complex(dp):: L0,L1,Lnrat,i3m,Lsm1_2mh,zab2
C  statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
c      zba2(k1,k2,k3,k4)=za(k4,k2)*zb(k2,k1)+za(k4,k3)*zb(k3,k1)
      DELTA3(k1,k2,k3,k4,k5,k6)=
     & +s(k1,k2)**2+s(k3,k4)**2+s(k5,k6)**2
     & -2d0*(s(k1,k2)*s(k3,k4)+s(k3,k4)*s(k5,k6)+s(k5,k6)*s(k1,k2))
      ddelta(k1,k2,k3,k4,k5,k6)=s(k1,k2)-s(k3,k4)-s(k5,k6)
      IDELTA3(k1,k2,k3,k4,k5,k6)=1d0/(
     & s(k1,k2)**2+s(k3,k4)**2+s(k5,k6)**2
     & -2d0*(s(k1,k2)*s(k3,k4)+s(k1,k2)*s(k5,k6)+s(k5,k6)*s(k3,k4)))

      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
c      p12(nu)=mom(k1,nu)+p2(nu)
c      p34(nu)=p3(nu)+mom(k4,nu)
      p23(nu)=p2(nu)+p3(nu)
      p14(nu)=mom(k1,nu)+mom(k4,nu)
      p123(nu)=p23(nu)+mom(k1,nu)
      p234(nu)=p23(nu)+mom(k4,nu)
      p1234(nu)=p123(nu)+mom(k4,nu)
      enddo

      s14=p14(4)**2-p14(1)**2-p14(2)**2-p14(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2
      s123=p123(4)**2-p123(1)**2-p123(2)**2-p123(3)**2
      s234=p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2
c      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
c      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s1234=p1234(4)**2-p1234(1)**2-p1234(2)**2-p1234(3)**2
      s56=s1234
c      p1Dp3=mom(k1,4)*p3(4)-mom(k1,1)*p3(1)
c     &     -mom(k1,2)*p3(2)-mom(k1,3)*p3(3)
c      p2Dp4=mom(k4,4)*p2(4)-mom(k4,1)*p2(1)
c     &     -mom(k4,2)*p2(2)-mom(k4,3)*p2(3)
c      p2Dp3=p2(4)*p3(4)-p2(1)*p3(1)
c     &     -p2(2)*p3(2)-p2(3)*p3(3)
c      p1Dp4=mom(k1,4)*mom(k4,4)-mom(k1,1)*mom(k4,1)
c     &     -mom(k1,2)*mom(k4,2)-mom(k1,3)*mom(k4,3)
c      p1Dp2=mom(k1,4)*p2(4)-mom(k1,1)*p2(1)
c     &     -mom(k1,2)*p2(2)-mom(k1,3)*p2(3)
c      p3Dp4=mom(k4,4)*p3(4)-mom(k4,1)*p3(1)
c     &     -mom(k4,2)*p3(2)-mom(k4,3)*p3(3)
c      msq=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
c      mb=sqrt(msq)

c---------------------------- BOXES --------------------------------------


c--- this is the subleading (23,1,4) mm box with factor mb/zb(k2,k3) removed
      coeff(4,d23x1x4)=-cone/(s23*s1234*s123)
     *  *(za(k1,k2)*za(k5,k6)*za(k1,k3)
     &  *zb(k4,k5)**2*zb(k2,k3)/zab2(k1,k2,k3,k4)
     & +za(k2,k3)*zb(k5,k6)
     & *zab2(k6,k2,k3,k1)**2*zab2(k4,k1,k2,k3)*zab2(k4,k1,k3,k2)
     & /zab2(k4,k2,k3,k1)**3)

c--- this is the subleading (23,4,1) mm box with factor mb/zb(k2,k3) removed
      coeff(4,d23x4x1)=-cone/(s23*s1234*s234)
     & *(za(k1,k6)**2*za(k2,k3)
     & *zb(k2,k4)*zb(k3,k4)*zb(k5,k6)/zab2(k1,k2,k3,k4)
     & +za(k5,k6)*zb(k2,k3)
     & *zab2(k4,k2,k3,k5)**2*zab2(k3,k2,k4,k1)*zab2(k2,k3,k4,k1)
     & /zab2(k4,k2,k3,k1)**3)


C     This is a coeff(3,c14x23bis) with factor mb/zb(k2,k3) removed
      coeff(3,c14x23bis)=  + IDELTA3(k2,k3,k1,k4,k5,k6)**2 * ( 3.D0/2.D0
     &    /(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*s(k1,
     &    k4)*ddelta(k1,k4,k2,k3,k5,k6)*s123 - 3.D0/2.D0/(zab(k4,p23,k1
     &    ))*za(k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*s(k1,k4)*ddelta(k1,k4
     &    ,k2,k3,k5,k6)*s234 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(
     &    k2,k5)*zab(k2,p1234,k2)*DELTA3(k3,k2,k1,k4,k5,k6)*s123 + 1.D0/
     &    2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*
     &    DELTA3(k3,k2,k1,k4,k5,k6)*s234 - 3.D0/2.D0/(zab(k4,p23,k1))*
     &    za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s(k1,k4)*ddelta(k1,k4,k2
     &    ,k3,k5,k6)*s123 + 3.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,
     &    k5)*zab(k3,p1234,k3)*s(k1,k4)*ddelta(k1,k4,k2,k3,k5,k6)*s234
     &     + 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k3,
     &    p1234,k3)*DELTA3(k3,k2,k1,k4,k5,k6)*s123 - 1.D0/2.D0/(zab(k4,
     &    p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*DELTA3(k3,k2,k1
     &    ,k4,k5,k6)*s234 + 3.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,
     &    k5)*zab(k2,p1234,k2)*s(k1,k4)*ddelta(k1,k4,k2,k3,k5,k6)*s123
     &     )
      coeff(3,c14x23bis) = coeff(3,c14x23bis) + IDELTA3(k2,k3,k1,k4,k5,
     & k6)**2 * (  - 3.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*
     &    zab(k2,p1234,k2)*s(k1,k4)*ddelta(k1,k4,k2,k3,k5,k6)*s234 - 1.D
     &    0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*
     &    DELTA3(k3,k2,k1,k4,k5,k6)*s123 + 1.D0/2.D0/(zab(k4,p23,k1))*
     &    za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*DELTA3(k3,k2,k1,k4,k5,k6
     &    )*s234 - 3.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(
     &    k3,p1234,k3)*s(k1,k4)*ddelta(k1,k4,k2,k3,k5,k6)*s123 + 3.D0/2.
     &    D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k3,p1234,k3)*s(k1
     &    ,k4)*ddelta(k1,k4,k2,k3,k5,k6)*s234 + 1.D0/2.D0/(zab(k4,p23,
     &    k1))*za(k3,k6)*zb(k3,k5)*zab(k3,p1234,k3)*DELTA3(k3,k2,k1,k4,
     &    k5,k6)*s123 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*
     &    zab(k3,p1234,k3)*DELTA3(k3,k2,k1,k4,k5,k6)*s234 )
      coeff(3,c14x23bis) = coeff(3,c14x23bis) + IDELTA3(k2,k3,k1,k4,k5,
     & k6) * (  - 2.D0/(zab(k4,p23,k1))*za(k2,k3)*zb(k5,k6)*zab(k6,p14,
     &    k2)*zab(k6,p14,k3) + 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(
     &    k2,k5)*s(k1,k4)*s123 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*
     &    zb(k2,k5)*s(k1,k4)*s234 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6
     &    )*zb(k3,k5)*s(k1,k4)*s123 + 1.D0/2.D0/(zab(k4,p23,k1))*za(k3,
     &    k6)*zb(k3,k5)*s(k1,k4)*s234-1d0/(zab(k4,p23,k1))/(zab(k4,p23,
     &    k1))*za(k2,k3)*za(k4,k6)*zb(k1,k2)*zb(k3,k5)*ddelta(k1,k4,k2,
     &    k3,k5,k6)*ddelta(k2,k3,k1,k4,k5,k6)-1d0/(zab(k4,p23,k1))/(
     &    zab(k4,p23,k1))*za(k2,k3)*za(k4,k6)*zb(k1,k3)*zb(k2,k5)*
     &    ddelta(k1,k4,k2,k3,k5,k6)*ddelta(k2,k3,k1,k4,k5,k6)-1d0/(zab(
     &    k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k6)*za(k4,k6)*zb(k1,k2)*
     &    zb(k5,k6)*ddelta(k1,k4,k2,k3,k5,k6)*ddelta(k5,k6,k1,k4,k2,k3)
     &    +1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k6)*za(k4,k6)*
     &    zb(k1,k3)*zb(k5,k6)*ddelta(k1,k4,k2,k3,k5,k6)*ddelta(k5,k6,k1
     &    ,k4,k2,k3) )
      coeff(3,c14x23bis) = coeff(3,c14x23bis) + IDELTA3(k5,k6,k1,k4,k2,
     & k3)**2 * ( 3.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k2
     &    ,p1234,k2)*s(k1,k4)*ddelta(k1,k4,k5,k6,k2,k3)*s123 - 3.D0/2.D0
     &    /(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*s(k1,
     &    k4)*ddelta(k1,k4,k5,k6,k2,k3)*s234 - 1.D0/2.D0/(zab(k4,p23,k1
     &    ))*za(k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*DELTA3(k6,k5,k1,k4,k2
     &    ,k3)*s123 + 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*
     &    zab(k2,p1234,k2)*DELTA3(k6,k5,k1,k4,k2,k3)*s234 - 3.D0/2.D0/(
     &    zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s(k1,k4)
     &    *ddelta(k1,k4,k5,k6,k2,k3)*s123 + 3.D0/2.D0/(zab(k4,p23,k1))*
     &    za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s(k1,k4)*ddelta(k1,k4,k5
     &    ,k6,k2,k3)*s234 + 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,
     &    k5)*zab(k3,p1234,k3)*DELTA3(k6,k5,k1,k4,k2,k3)*s123 - 1.D0/2.D
     &    0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*
     &    DELTA3(k6,k5,k1,k4,k2,k3)*s234 + 3.D0/2.D0/(zab(k4,p23,k1))*
     &    za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*s(k1,k4)*ddelta(k1,k4,k5
     &    ,k6,k2,k3)*s123 )
      coeff(3,c14x23bis) = coeff(3,c14x23bis) + IDELTA3(k5,k6,k1,k4,k2,
     & k3)**2 * (  - 3.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*
     &    zab(k2,p1234,k2)*s(k1,k4)*ddelta(k1,k4,k5,k6,k2,k3)*s234 - 1.D
     &    0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*
     &    DELTA3(k6,k5,k1,k4,k2,k3)*s123 + 1.D0/2.D0/(zab(k4,p23,k1))*
     &    za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*DELTA3(k6,k5,k1,k4,k2,k3
     &    )*s234 - 3.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(
     &    k3,p1234,k3)*s(k1,k4)*ddelta(k1,k4,k5,k6,k2,k3)*s123 + 3.D0/2.
     &    D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k3,p1234,k3)*s(k1
     &    ,k4)*ddelta(k1,k4,k5,k6,k2,k3)*s234 + 1.D0/2.D0/(zab(k4,p23,
     &    k1))*za(k3,k6)*zb(k3,k5)*zab(k3,p1234,k3)*DELTA3(k6,k5,k1,k4,
     &    k2,k3)*s123 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*
     &    zab(k3,p1234,k3)*DELTA3(k6,k5,k1,k4,k2,k3)*s234 )
      coeff(3,c14x23bis) = coeff(3,c14x23bis) + IDELTA3(k5,k6,k1,k4,k2,
     & k3) * ( 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*s(k1,k4)*
     &    s123 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k2,k6)*zb(k2,k5)*s(k1,k4
     &    )*s234 - 1.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*s(k1,
     &    k4)*s123 + 1.D0/2.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*s(
     &    k1,k4)*s234 - 2.D0/(zab(k4,p23,k1))*za(k5,k6)*zb(k2,k3)*zab(
     &    k2,p14,k5)*zab(k3,p14,k5)-1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1
     &    ))*za(k2,k4)*za(k3,k6)*zb(k1,k5)*zb(k2,k3)*ddelta(k1,k4,k5,k6
     &    ,k2,k3)*ddelta(k2,k3,k1,k4,k5,k6)+1d0/(zab(k4,p23,k1))/(zab(
     &    k4,p23,k1))*za(k2,k4)*za(k5,k6)*zb(k1,k5)*zb(k2,k5)*ddelta(k1
     &    ,k4,k5,k6,k2,k3)*ddelta(k5,k6,k1,k4,k2,k3)-1d0/(zab(k4,p23,k1
     &    ))/(zab(k4,p23,k1))*za(k2,k6)*za(k3,k4)*zb(k1,k5)*zb(k2,k3)*
     &    ddelta(k1,k4,k5,k6,k2,k3)*ddelta(k2,k3,k1,k4,k5,k6)-1d0/(zab(
     &    k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k4)*za(k5,k6)*zb(k1,k5)*
     &    zb(k3,k5)*ddelta(k1,k4,k5,k6,k2,k3)*ddelta(k5,k6,k1,k4,k2,k3)
     &     )
      coeff(3,c14x23bis) = coeff(3,c14x23bis) + 2.D0/(zab(k4,p23,k1))*
     & za(k1,k2)*za(k2,k6)*zb(k1,k2)*zb(k4,k5)*zab(k4,p1234,k2)*
     & s123**(-2) + 2.D0/(zab(k4,p23,k1))*za(k1,k2)*za(k3,k6)*zb(k1,k3)
     &    *zb(k4,k5)*zab(k4,p1234,k2)*s123**(-2) - 2.D0/(zab(k4,p23,k1)
     &    )*za(k1,k3)*za(k2,k6)*zb(k1,k2)*zb(k4,k5)*zab(k4,p1234,k3)*
     &    s123**(-2) - 2.D0/(zab(k4,p23,k1))*za(k1,k3)*za(k3,k6)*zb(k1,
     &    k3)*zb(k4,k5)*zab(k4,p1234,k3)*s123**(-2) - 2.D0/(zab(k4,p23,
     &    k1))*za(k1,k6)*za(k2,k4)*zb(k2,k4)*zb(k2,k5)*zab(k2,p1234,k1)
     &    *s234**(-2) + 2.D0/(zab(k4,p23,k1))*za(k1,k6)*za(k2,k4)*zb(k2
     &    ,k5)*zb(k3,k4)*zab(k3,p1234,k1)*s234**(-2) - 2.D0/(zab(k4,p23
     &    ,k1))*za(k1,k6)*za(k3,k4)*zb(k2,k4)*zb(k3,k5)*zab(k2,p1234,k1
     &    )*s234**(-2) + 2.D0/(zab(k4,p23,k1))*za(k1,k6)*za(k3,k4)*zb(
     &    k3,k4)*zb(k3,k5)*zab(k3,p1234,k1)*s234**(-2)-1d0/(zab(k4,p23,
     &    k1))/(zab(k4,p23,k1))*za(k2,k4)*zb(k1,k5)*zab(k6,p234,k2) +
     &  1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k4)*zb(k1,k5)*zab(
     &    k6,p234,k3)
      coeff(3,c14x23bis) = coeff(3,c14x23bis)-1d0/(zab(k4,p23,k1))/(
     & zab(k4,p23,k1))*za(k4,k6)*zb(k1,k2)*zab(k2,p123,k5)+1d0/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k3)*zab(k3,p123,k5)


c---------------------------- BUBBLES ------------------------------------

      bubrat=  + Lnrat( - s123, - s56) * (  - 3.D0/2.D0/(za(k5,k6))/(
     &    zb(k2,k3))/(zab(k4,p23,k1))*zab(k6,p123,k2)*zab(k6,p123,k3)*
     &    s123**(-1) - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1
     &    ))/(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k2)*zab(k6,p123,k3) - 1.D0
     &    /2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))/(zab(k4,p23,k1
     &    ))*za(k4,k6)*zb(k1,k3)*zab(k6,p123,k2) )
      bubrat = bubrat + Lnrat( - s123, - s23) * (  - 3.D0/2.D0/(za(k5,
     &    k6))/(zb(k2,k3))/(zab(k4,p23,k1))*zab(k6,p123,k2)*zab(k6,p123
     &    ,k3)*s123**(-1) - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k2)*zab(k6,p123,k3)
     &     - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))/(zab(k4
     &    ,p23,k1))*za(k4,k6)*zb(k1,k3)*zab(k6,p123,k2) )
      bubrat = bubrat + Lnrat( - s56, - s14)*IDELTA3(k5,k6,k1,k4,k2,k3)
     & **2*ddelta(k2,k3,k5,k6,k1,k4) * ( 3.D0/(zab(k4,p23,k1))*za(k2,k6
     &    )*zb(k2,k5)*zab(k2,p1234,k2)*s123 - 3.D0/(zab(k4,p23,k1))*za(
     &    k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*s234 - 3.D0/(zab(k4,p23,k1)
     &    )*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s123 + 3.D0/(zab(k4,
     &    p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s234 + 3.D0/(
     &    zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*s123 - 3.
     &    D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*s234
     &     - 3.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k3,p1234,k3)
     &    *s123 + 3.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k3,
     &    p1234,k3)*s234 )
      bubrat = bubrat + Lnrat( - s56, - s14)*IDELTA3(k5,k6,k1,k4,k2,k3)
     & *ddelta(k1,k4,k5,k6,k2,k3) * ( -1d0/(za(k2,k3))/(zab(k4,p23,k1))
     &    /(zab(k4,p23,k1))*za(k2,k4)*za(k3,k6)*zb(k1,k5)*s234-1d0/(za(
     &    k2,k3))/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k6)*za(k3,k4)
     &    *zb(k1,k5)*s234-1d0/(zb(k2,k3))/(zab(k4,p23,k1))/(zab(k4,p23,
     &    k1))*za(k4,k6)*zb(k1,k2)*zb(k3,k5)*s123-1d0/(zb(k2,k3))/(zab(
     &    k4,p23,k1))/(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k3)*zb(k2,k5)*
     &    s123 )
      bubrat = bubrat + Lnrat( - s56, - s14)*IDELTA3(k5,k6,k1,k4,k2,k3)
     &  * ( -1d0/(za(k2,k3))/(zab(k4,p23,k1))*za(k1,k2)*za(k3,k6)*zb(k1
     &    ,k5)*s234-1d0/(za(k2,k3))/(zab(k4,p23,k1))*za(k1,k3)*za(k2,k6
     &    )*zb(k1,k5)*s234-1d0/(za(k2,k3))/(zab(k4,p23,k1))*za(k2,k4)*
     &    za(k3,k6)*zb(k4,k5)*s123-1d0/(za(k2,k3))/(zab(k4,p23,k1))*za(
     &    k2,k6)*za(k3,k4)*zb(k4,k5)*s123-1d0/(zb(k2,k3))/(zab(k4,p23,
     &    k1))*za(k1,k6)*zb(k1,k2)*zb(k3,k5)*s234-1d0/(zb(k2,k3))/(zab(
     &    k4,p23,k1))*za(k1,k6)*zb(k1,k3)*zb(k2,k5)*s234-1d0/(zb(k2,k3)
     &    )/(zab(k4,p23,k1))*za(k4,k6)*zb(k2,k4)*zb(k3,k5)*s123-1d0/(
     &    zb(k2,k3))/(zab(k4,p23,k1))*za(k4,k6)*zb(k2,k5)*zb(k3,k4)*
     &    s123+1d0/(zab(k4,p23,k1))*za(k1,k2)*za(k5,k6)*zb(k1,k5)*zb(k2
     &    ,k5)-1d0/(zab(k4,p23,k1))*za(k1,k3)*za(k5,k6)*zb(k1,k5)*zb(k3
     &    ,k5)-1d0/(zab(k4,p23,k1))*za(k1,k6)*za(k2,k6)*zb(k1,k2)*zb(k5
     &    ,k6)+1d0/(zab(k4,p23,k1))*za(k1,k6)*za(k3,k6)*zb(k1,k3)*zb(k5
     &    ,k6)+1d0/(zab(k4,p23,k1))*za(k2,k4)*za(k5,k6)*zb(k2,k5)*zb(k4
     &    ,k5) )
      bubrat = bubrat + Lnrat( - s56, - s14)*IDELTA3(k5,k6,k1,k4,k2,k3)
     &  * ( -1d0/(zab(k4,p23,k1))*za(k2,k6)*za(k4,k6)*zb(k2,k4)*zb(k5,
     &    k6)-1d0/(zab(k4,p23,k1))*za(k3,k4)*za(k5,k6)*zb(k3,k5)*zb(k4,
     &    k5)+1d0/(zab(k4,p23,k1))*za(k3,k6)*za(k4,k6)*zb(k3,k4)*zb(k5,
     &    k6) - 2.D0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k3)*za(k4,
     &    k6)*za(k5,k6)*zb(k1,k2)*zb(k3,k5)*zb(k5,k6) - 2.D0/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k3)*za(k4,k6)*za(k5,k6)*zb(k1
     &    ,k3)*zb(k2,k5)*zb(k5,k6) - 2.D0/(zab(k4,p23,k1))/(zab(k4,p23,
     &    k1))*za(k2,k4)*za(k3,k6)*za(k5,k6)*zb(k1,k5)*zb(k2,k3)*zb(k5,
     &    k6)+1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*za(k5,k6)
     &    *zb(k1,k5)*zb(k2,k5)*s123-1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1
     &    ))*za(k2,k4)*za(k5,k6)*zb(k1,k5)*zb(k2,k5)*s234 - 2.D0/(zab(
     &    k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k6)*za(k3,k4)*za(k5,k6)*
     &    zb(k1,k5)*zb(k2,k3)*zb(k5,k6)+1d0/(zab(k4,p23,k1))/(zab(k4,
     &    p23,k1))*za(k2,k6)*za(k4,k6)*zb(k1,k2)*zb(k5,k6)*s123-1d0/(
     &    zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k6)*za(k4,k6)*zb(k1,k2
     &    )*zb(k5,k6)*s234 )
      bubrat = bubrat + Lnrat( - s56, - s14)*IDELTA3(k5,k6,k1,k4,k2,k3)
     &  * ( -1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k4)*za(k5,k6)*
     &    zb(k1,k5)*zb(k3,k5)*s123+1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1)
     &    )*za(k3,k4)*za(k5,k6)*zb(k1,k5)*zb(k3,k5)*s234-1d0/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k3,k6)*za(k4,k6)*zb(k1,k3)*zb(k5
     &    ,k6)*s123+1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k6)*za(
     &    k4,k6)*zb(k1,k3)*zb(k5,k6)*s234 )
      bubrat = bubrat + Lnrat( - s234, - s56) * (  - 3.D0/2.D0/(za(k2,
     &    k3))/(zb(k5,k6))/(zab(k4,p23,k1))*zab(k2,p234,k5)*zab(k3,p234
     &    ,k5)*s234**(-1) + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*zb(k1,k5)*zab(k3,p234,k5)
     &     + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))/(zab(k4
     &    ,p23,k1))*za(k3,k4)*zb(k1,k5)*zab(k2,p234,k5) )
      bubrat = bubrat + Lnrat( - s234, - s23) * (  - 3.D0/2.D0/(za(k2,
     &    k3))/(zb(k5,k6))/(zab(k4,p23,k1))*zab(k2,p234,k5)*zab(k3,p234
     &    ,k5)*s234**(-1) + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*zb(k1,k5)*zab(k3,p234,k5)
     &     + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))/(zab(k4
     &    ,p23,k1))*za(k3,k4)*zb(k1,k5)*zab(k2,p234,k5) )
      bubrat = bubrat + Lnrat( - s23, - s14)*IDELTA3(k2,k3,k1,k4,k5,k6)
     & **2*ddelta(k5,k6,k2,k3,k1,k4) * ( 3.D0/(zab(k4,p23,k1))*za(k2,k6
     &    )*zb(k2,k5)*zab(k2,p1234,k2)*s123 - 3.D0/(zab(k4,p23,k1))*za(
     &    k2,k6)*zb(k2,k5)*zab(k2,p1234,k2)*s234 - 3.D0/(zab(k4,p23,k1)
     &    )*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s123 + 3.D0/(zab(k4,
     &    p23,k1))*za(k2,k6)*zb(k2,k5)*zab(k3,p1234,k3)*s234 + 3.D0/(
     &    zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*s123 - 3.
     &    D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k2,p1234,k2)*s234
     &     - 3.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k3,p1234,k3)
     &    *s123 + 3.D0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*zab(k3,
     &    p1234,k3)*s234 )
      bubrat = bubrat + Lnrat( - s23, - s14)*IDELTA3(k2,k3,k1,k4,k5,k6)
     & *ddelta(k1,k4,k2,k3,k5,k6) * ( -1d0/(za(k5,k6))/(zab(k4,p23,k1))
     &    /(zab(k4,p23,k1))*za(k2,k6)*za(k4,k6)*zb(k1,k2)*s123+1d0/(za(
     &    k5,k6))/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k6)*za(k4,k6)
     &    *zb(k1,k3)*s123+1d0/(zb(k5,k6))/(zab(k4,p23,k1))/(zab(k4,p23,
     &    k1))*za(k2,k4)*zb(k1,k5)*zb(k2,k5)*s234-1d0/(zb(k5,k6))/(zab(
     &    k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k4)*zb(k1,k5)*zb(k3,k5)*
     &    s234 )
      bubrat = bubrat + Lnrat( - s23, - s14)*IDELTA3(k2,k3,k1,k4,k5,k6)
     &  * ( -1d0/(za(k5,k6))/(zab(k4,p23,k1))*za(k1,k6)*za(k2,k6)*zb(k1
     &    ,k2)*s123+1d0/(za(k5,k6))/(zab(k4,p23,k1))*za(k1,k6)*za(k3,k6
     &    )*zb(k1,k3)*s123-1d0/(za(k5,k6))/(zab(k4,p23,k1))*za(k2,k6)*
     &    za(k4,k6)*zb(k2,k4)*s234+1d0/(za(k5,k6))/(zab(k4,p23,k1))*za(
     &    k3,k6)*za(k4,k6)*zb(k3,k4)*s234+1d0/(zb(k5,k6))/(zab(k4,p23,
     &    k1))*za(k1,k2)*zb(k1,k5)*zb(k2,k5)*s123-1d0/(zb(k5,k6))/(zab(
     &    k4,p23,k1))*za(k1,k3)*zb(k1,k5)*zb(k3,k5)*s123+1d0/(zb(k5,k6)
     &    )/(zab(k4,p23,k1))*za(k2,k4)*zb(k2,k5)*zb(k4,k5)*s234-1d0/(
     &    zb(k5,k6))/(zab(k4,p23,k1))*za(k3,k4)*zb(k3,k5)*zb(k4,k5)*
     &    s234-1d0/(zab(k4,p23,k1))*za(k1,k2)*za(k3,k6)*zb(k1,k5)*zb(k2
     &    ,k3)-1d0/(zab(k4,p23,k1))*za(k1,k3)*za(k2,k6)*zb(k1,k5)*zb(k2
     &    ,k3)-1d0/(zab(k4,p23,k1))*za(k1,k6)*za(k2,k3)*zb(k1,k2)*zb(k3
     &    ,k5)-1d0/(zab(k4,p23,k1))*za(k1,k6)*za(k2,k3)*zb(k1,k3)*zb(k2
     &    ,k5)-1d0/(zab(k4,p23,k1))*za(k2,k3)*za(k4,k6)*zb(k2,k4)*zb(k3
     &    ,k5) )
      bubrat = bubrat + Lnrat( - s23, - s14)*IDELTA3(k2,k3,k1,k4,k5,k6)
     &  * ( -1d0/(zab(k4,p23,k1))*za(k2,k3)*za(k4,k6)*zb(k2,k5)*zb(k3,
     &    k4)-1d0/(zab(k4,p23,k1))*za(k2,k4)*za(k3,k6)*zb(k2,k3)*zb(k4,
     &    k5)-1d0/(zab(k4,p23,k1))*za(k2,k6)*za(k3,k4)*zb(k2,k3)*zb(k4,
     &    k5) - 2.D0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k3)*za(k2,
     &    k6)*za(k4,k6)*zb(k1,k2)*zb(k2,k3)*zb(k5,k6) + 2.D0/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k3)*za(k3,k6)*za(k4,k6)*zb(k1
     &    ,k3)*zb(k2,k3)*zb(k5,k6)+1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1)
     &    )*za(k2,k3)*za(k4,k6)*zb(k1,k2)*zb(k3,k5)*s123-1d0/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k3)*za(k4,k6)*zb(k1,k2)*zb(k3
     &    ,k5)*s234+1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k3)*za(
     &    k4,k6)*zb(k1,k3)*zb(k2,k5)*s123-1d0/(zab(k4,p23,k1))/(zab(k4,
     &    p23,k1))*za(k2,k3)*za(k4,k6)*zb(k1,k3)*zb(k2,k5)*s234-1d0/(
     &    zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*za(k3,k6)*zb(k1,k5
     &    )*zb(k2,k3)*s123+1d0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,
     &    k4)*za(k3,k6)*zb(k1,k5)*zb(k2,k3)*s234 )
      bubrat = bubrat + Lnrat( - s23, - s14)*IDELTA3(k2,k3,k1,k4,k5,k6)
     &  * (  - 2.D0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*za(k5,
     &    k6)*zb(k1,k5)*zb(k2,k5)*s23-1d0/(zab(k4,p23,k1))/(zab(k4,p23,
     &    k1))*za(k2,k6)*za(k3,k4)*zb(k1,k5)*zb(k2,k3)*s123+1d0/(zab(k4
     &    ,p23,k1))/(zab(k4,p23,k1))*za(k2,k6)*za(k3,k4)*zb(k1,k5)*zb(
     &    k2,k3)*s234 + 2.D0/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k3,k4
     &    )*za(k5,k6)*zb(k1,k5)*zb(k3,k5)*s23 )
      bubrat = bubrat + Lnrat( - s14, - s56) * (  - 3.D0/2.D0/(za(k2,k3
     &    ))/(zb(k5,k6))/(zab(k4,p23,k1))*zab(k2,p234,k5)*zab(k3,p234,
     &    k5)*s234**(-1) + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*zb(k1,k5)*zab(k3,p234,k5)
     &     + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))/(zab(k4
     &    ,p23,k1))*za(k3,k4)*zb(k1,k5)*zab(k2,p234,k5) - 3.D0/2.D0/(
     &    za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))*zab(k6,p123,k2)*zab(
     &    k6,p123,k3)*s123**(-1) - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(
     &    zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k2)*zab(k6,
     &    p123,k3) - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))
     &    /(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k3)*zab(k6,p123,k2) )
      bubrat = bubrat + Lnrat( - s14, - s23) * (  - 3.D0/2.D0/(za(k2,k3
     &    ))/(zb(k5,k6))/(zab(k4,p23,k1))*zab(k2,p234,k5)*zab(k3,p234,
     &    k5)*s234**(-1) + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,
     &    p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*zb(k1,k5)*zab(k3,p234,k5)
     &     + 1.D0/2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))/(zab(k4
     &    ,p23,k1))*za(k3,k4)*zb(k1,k5)*zab(k2,p234,k5) - 3.D0/2.D0/(
     &    za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))*zab(k6,p123,k2)*zab(
     &    k6,p123,k3)*s123**(-1) - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(
     &    zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k2)*zab(k6,
     &    p123,k3) - 1.D0/2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))
     &    /(zab(k4,p23,k1))*za(k4,k6)*zb(k1,k3)*zab(k6,p123,k2) )
      bubrat = bubrat + L0( - s123, - s56) * ( 2.D0/(za(k5,k6))/(zb(k2,
     &    k3))/(zab(k4,p23,k1))*za(k4,k6)*zb(k2,k4)*zab(k6,p123,k3)*
     &    s56**(-1) + 2.D0/(za(k5,k6))/(zb(k2,k3))/(zab(k4,p23,k1))*za(
     &    k4,k6)*zb(k3,k4)*zab(k6,p123,k2)*s56**(-1)+1d0/(za(k5,k6))/(
     &    zb(k2,k3))/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k4,k6)**2*zb(
     &    k1,k2)*zb(k3,k4)*s56**(-1)*s123+1d0/(za(k5,k6))/(zb(k2,k3))/(
     &    zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k4,k6)**2*zb(k1,k3)*zb(k2
     &    ,k4)*s56**(-1)*s123 )
      bubrat = bubrat + L0( - s123, - s23) * (  - 2.D0/(za(k2,k3))/(za(
     &    k5,k6))/(zb(k2,k3))/(zb(k2,k3))/(zab(k4,p23,k1))*za(k1,k6)*
     &    zb(k1,k2)*zab(k6,p123,k3) - 2.D0/(za(k2,k3))/(za(k5,k6))/(zb(
     &    k2,k3))/(zb(k2,k3))/(zab(k4,p23,k1))*za(k1,k6)*zb(k1,k3)*zab(
     &    k6,p123,k2) - 2.D0/(za(k2,k3))/(za(k5,k6))/(zb(k2,k3))/(zb(k2
     &    ,k3))/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k1,k6)*za(k4,k6)*
     &    zb(k1,k2)*zb(k1,k3)*s123 )
      bubrat = bubrat + L0( - s234, - s56) * (  - 2.D0/(za(k2,k3))/(zb(
     &    k5,k6))/(zab(k4,p23,k1))*za(k1,k2)*zb(k1,k5)*zab(k3,p234,k5)*
     &    s56**(-1) - 2.D0/(za(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))*za(
     &    k1,k3)*zb(k1,k5)*zab(k2,p234,k5)*s56**(-1)+1d0/(za(k2,k3))/(
     &    zb(k5,k6))/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k1,k2)*za(k3,
     &    k4)*zb(k1,k5)**2*s234*s56**(-1)+1d0/(za(k2,k3))/(zb(k5,k6))/(
     &    zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k1,k3)*za(k2,k4)*zb(k1,k5
     &    )**2*s234*s56**(-1) )
      bubrat = bubrat + L0( - s234, - s23) * ( 2.D0/(za(k2,k3))/(za(k2,
     &    k3))/(zb(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))*za(k2,k4)*zb(k4
     &    ,k5)*zab(k3,p234,k5) + 2.D0/(za(k2,k3))/(za(k2,k3))/(zb(k2,k3
     &    ))/(zb(k5,k6))/(zab(k4,p23,k1))*za(k3,k4)*zb(k4,k5)*zab(k2,
     &    p234,k5) - 2.D0/(za(k2,k3))/(za(k2,k3))/(zb(k2,k3))/(zb(k5,k6
     &    ))/(zab(k4,p23,k1))/(zab(k4,p23,k1))*za(k2,k4)*za(k3,k4)*zb(
     &    k1,k5)*zb(k4,k5)*s234 )
      bubrat = bubrat + L1( - s56, - s123) *(1d0/(za(k5,k6))/(zb(k2,k3)
     &    )/(zab(k4,p23,k1))*za(k4,k6)**2*zb(k2,k4)*zb(k3,k4)*
     &    s123**(-1) )
      bubrat = bubrat + L1( - s56, - s234) *(1d0/(za(k2,k3))/(zb(k5,k6)
     &    )/(zab(k4,p23,k1))*za(k1,k2)*za(k1,k3)*zb(k1,k5)**2*
     &    s234**(-1) )
      bubrat = bubrat + L1( - s23, - s123) *(1d0/(za(k5,k6))/(zb(k2,k3)
     &    )/(zab(k4,p23,k1))*za(k1,k6)**2*zb(k1,k2)*zb(k1,k3)*
     &    s123**(-1) )
      bubrat = bubrat + L1( - s23, - s234) *(1d0/(za(k2,k3))/(zb(k5,k6)
     &    )/(zab(k4,p23,k1))*za(k2,k4)*za(k3,k4)*zb(k4,k5)**2*
     &    s234**(-1) )
      bubrat = bubrat + IDELTA3(k2,k3,k1,k4,k5,k6)*ddelta(k1,k4,k2,k3,
     & k5,k6) * ( 1.D0/2.D0/(za(k2,k3))/(za(k5,k6))/(zab(k4,p23,k1))*
     &    za(k2,k6)*za(k3,k6)*s123 - 1.D0/2.D0/(za(k2,k3))/(za(k5,k6))
     &    /(zab(k4,p23,k1))*za(k2,k6)*za(k3,k6)*s234 - 1.D0/2.D0/(zb(k2
     &    ,k3))/(zb(k5,k6))/(zab(k4,p23,k1))*zb(k2,k5)*zb(k3,k5)*s123
     &     + 1.D0/2.D0/(zb(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))*zb(k2,
     &    k5)*zb(k3,k5)*s234 )
      bubrat = bubrat + IDELTA3(k2,k3,k1,k4,k5,k6) *(1d0/(zab(k4,p23,k1
     &    ))*za(k2,k6)*zb(k2,k5)*s123-1d0/(zab(k4,p23,k1))*za(k2,k6)*
     &    zb(k2,k5)*s234-1d0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*s123
     &    +1d0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*s234 )
      bubrat = bubrat + IDELTA3(k5,k6,k1,k4,k2,k3)*ddelta(k1,k4,k5,k6,
     & k2,k3) * ( 1.D0/2.D0/(za(k2,k3))/(za(k5,k6))/(zab(k4,p23,k1))*
     &    za(k2,k6)*za(k3,k6)*s123 - 1.D0/2.D0/(za(k2,k3))/(za(k5,k6))
     &    /(zab(k4,p23,k1))*za(k2,k6)*za(k3,k6)*s234 - 1.D0/2.D0/(zb(k2
     &    ,k3))/(zb(k5,k6))/(zab(k4,p23,k1))*zb(k2,k5)*zb(k3,k5)*s123
     &     + 1.D0/2.D0/(zb(k2,k3))/(zb(k5,k6))/(zab(k4,p23,k1))*zb(k2,
     &    k5)*zb(k3,k5)*s234 )
      bubrat = bubrat + IDELTA3(k5,k6,k1,k4,k2,k3) *(1d0/(zab(k4,p23,k1
     &    ))*za(k2,k6)*zb(k2,k5)*s123-1d0/(zab(k4,p23,k1))*za(k2,k6)*
     &    zb(k2,k5)*s234-1d0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*s123
     &    +1d0/(zab(k4,p23,k1))*za(k3,k6)*zb(k3,k5)*s234 )


      BDK1211mm=
     & +Lsm1_2mh(s(k1,k4),t(k1,k2,k3),s(k2,k3),s(k5,k6))
     & *coeff(4,d23x1x4)*2d0
     & +Lsm1_2mh(s(k1,k4),t(k6,k5,k1),s(k6,k5),s(k2,k3))
     & *coeff(4,d23x4x1)*2d0
     & +i3m(s(k1,k4),s(k2,k3),s(k5,k6))*coeff(3,c14x23bis)

      BDK1211mm=BDK1211mm+bubrat


      return
      end
