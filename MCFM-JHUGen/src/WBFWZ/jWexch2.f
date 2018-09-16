      subroutine jW2exch2(p1,p2,p3,p4,p5,p6,p7,p8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer:: jdu1,p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,amp(2,2),rxw,
     & propW34,propz56,propW17,propw3456,
     & gameW(2),gamzeq(2,2,2),gamzee(2,2),gamzne(2,2)
      real(dp):: qn,t3,t4,s345,s356,s456,
     & s256,s3456,s34,s56,s17,s178,s234,
     & s568,s127,s348,twop34Dp1278,q3,q4,l3,l4
C     amp(h56)
      parameter(qn=0._dp)
C-----Begin statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      t4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s17=s(p1,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s127=t3(p1,p2,p7)
      s234=t3(p2,p3,p4)
      s256=t3(p2,p5,p6)
      s178=t3(p1,p7,p8)
      s345=t3(p3,p4,p5)
      s356=t3(p3,p5,p6)
      s456=t3(p4,p5,p6)
      s568=t3(p5,p6,p8)
      s348=t3(p3,p4,p8)
      s3456=t4(p3,p4,p5,p6)
      twop34Dp1278=s56-s3456-s34

      propW34=s34-cwmass2
      propz56=s56-czmass2
      propW17=s17-cwmass2
      propw3456=s3456-cwmass2

      q3=qn
      l3=ln
      q4=qe
      l4=le
      gamzeq(:,:,:)=czip
      do jdu1=1,2
      gamzeq(jdu1,1,1)=qe*Q(jdu1)/s56+le*L(jdu1)/propz56
      gamzeq(jdu1,2,1)=qe*Q(jdu1)/s56+re*L(jdu1)/propz56
      enddo
      gamzee(1,1)=q4*qe/s56+l4*le/propz56
      gamzee(1,2)=q4*qe/s56+l4*re/propz56
      gamzne(1,1)=q3*qe/s56+l3*le/propz56
      gamzne(1,2)=q3*qe/s56+l3*re/propz56
      gameW(1)=qe/s56+le*rxw/propz56
      gameW(2)=qe/s56+re*rxw/propz56

      amp(1,1)= + cxw**(-3)*propW17**(-1)*propW34**(-1)*propW3456**(-1)
     & *s127**(-1) * (  - za(p7,p1)*za(p3,p5)*zb(p1,p6)*zb(p1,p2)*zba2(
     &    p4,p3,p5,p8)*s345**(-1) + za(p7,p2)*za(p3,p5)*zb(p1,p2)*zb(p6
     &    ,p2)*zba2(p4,p3,p5,p8)*s345**(-1) )
      amp(1,1) = amp(1,1) + cxw**(-3)*propW17**(-1)*propW34**(-1)*
     & propW3456**(-1) * (  - half*za(p7,p8)*za(p3,p5)*zb(p1,p2)*
     &    zb(p4,p6)*cwmass2**(-1) )
      amp(1,1) = amp(1,1) + gamZee(1,1)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1)*s127**(-1) * ( two*za(p7,p1)*za(p3,p8)*za(p4,p5
     &    )*zb(p1,p4)*zb(p1,p2)*zb(p4,p6)*s456**(-1) - two*za(p7,p1)*
     &    za(p3,p8)*za(p5,p6)*zb(p1,p6)*zb(p1,p2)*zb(p4,p6)*s456**(-1)
     &     + two*za(p7,p2)*za(p3,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p2,p4,p6
     &    ,p5)*s456**(-1) )
      amp(1,1) = amp(1,1) + gamZee(1,1)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * ( za(p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*
     &    cwmass2**(-1) )
      amp(1,1) = amp(1,1) + gamZne(1,1)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1)*s127**(-1) * (  - two*za(p7,p1)*za(p3,p5)*za(p3
     &    ,p8)*zb(p1,p4)*zb(p1,p2)*zb(p3,p6)*s356**(-1) - two*za(p7,p1
     &    )*za(p3,p5)*za(p5,p8)*zb(p1,p4)*zb(p1,p2)*zb(p5,p6)*
     &    s356**(-1) - two*za(p7,p2)*za(p3,p5)*zb(p1,p2)*zb(p4,p2)*
     &    zba2(p6,p3,p5,p8)*s356**(-1) )
      amp(1,1) = amp(1,1) + gamZne(1,1)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * (  - za(p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*
     &    cwmass2**(-1) )
      amp(1,1) = amp(1,1) + gamZeq(1,1,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*s127**(-1) * ( two*za(p5,p8)*zb(p1,p2)*zba2(p4,p1
     &    ,p2,p7)*zba2(p6,p5,p8,p3)*s568**(-1) )
      amp(1,1) = amp(1,1) + gamZeq(1,1,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*s348**(-1) * ( two*za(p3,p8)*zb(p6,p2)*zba2(p1,p2
     &    ,p6,p5)*zba2(p4,p3,p8,p7)*s256**(-1) )
      amp(1,1) = amp(1,1) + gamZeq(2,1,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*s127**(-1)*s348**(-1) * ( two*za(p3,p8)*zb(p1,p2)
     &    *zba2(p4,p3,p8,p5)*zba2(p6,p1,p2,p7) )
      amp(1,1) = amp(1,1) + gameW(1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1)*s127**(-1) * ( two*za(p7,p1)*za(
     &    p3,p5)*zb(p1,p2)*zb(p4,p6)*zab2(p8,p3,p4,p1) + za(p7,p1)*za(
     &    p3,p5)*zb(p1,p2)*zb(p4,p6)*zab2(p8,p2,p7,p1) - two*za(p7,p1)
     &    *za(p5,p8)*zb(p1,p6)*zb(p1,p2)*zab2(p3,p5,p6,p4) + za(p7,p2)*
     &    za(p3,p5)*zb(p1,p2)*zb(p4,p6)*zab2(p8,p1,p7,p2) + two*za(p7,
     &    p2)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*zab2(p8,p3,p4,p2) + two*
     &    za(p7,p2)*za(p5,p8)*zb(p1,p2)*zb(p6,p2)*zab2(p3,p5,p6,p4) + 2.
     &    D0*za(p3,p8)*zb(p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p7) )
      amp(1,1) = amp(1,1) + gameW(1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1) * ( za(p7,p8)*za(p3,p5)*zb(p1,p2)*
     &    zb(p3,p6)*zab2(p3,p5,p6,p4)*cwmass2**(-1) - za(p7,p8)*za(p3,
     &    p5)*zb(p1,p2)*zb(p4,p5)*zab2(p5,p3,p4,p6)*cwmass2**(-1) - za(
     &    p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*cwmass2**(-1)*s3456 -
     &    za(p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*cwmass2**(-1)*
     &    twop34Dp1278 - za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p6)*zab2(
     &    p5,p3,p4,p6)*cwmass2**(-1) + za(p7,p8)*za(p4,p5)*zb(p1,p2)*
     &    zb(p4,p6)*zab2(p3,p5,p6,p4)*cwmass2**(-1) )
      amp(1,2)= + gamZee(1,2)*cxw**(-2)*propW17**(-1)*propW3456**(-1)*
     & s127**(-1) * ( two*za(p7,p1)*za(p3,p8)*za(p4,p6)*zb(p1,p4)*zb(
     &    p1,p2)*zb(p4,p5)*s456**(-1) + two*za(p7,p1)*za(p3,p8)*za(p5,
     &    p6)*zb(p1,p5)*zb(p1,p2)*zb(p4,p5)*s456**(-1) + two*za(p7,p2)
     &    *za(p3,p8)*zb(p1,p2)*zb(p4,p5)*zba2(p2,p4,p5,p6)*s456**(-1) )
      amp(1,2) = amp(1,2) + gamZee(1,2)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * ( za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*
     &    cwmass2**(-1) )
      amp(1,2) = amp(1,2) + gamZne(1,2)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1)*s127**(-1) * (  - two*za(p7,p1)*za(p3,p6)*za(p3
     &    ,p8)*zb(p1,p4)*zb(p1,p2)*zb(p3,p5)*s356**(-1) + two*za(p7,p1
     &    )*za(p3,p6)*za(p6,p8)*zb(p1,p4)*zb(p1,p2)*zb(p5,p6)*
     &    s356**(-1) - two*za(p7,p2)*za(p3,p6)*zb(p1,p2)*zb(p4,p2)*
     &    zba2(p5,p3,p6,p8)*s356**(-1) )
      amp(1,2) = amp(1,2) + gamZne(1,2)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * (  - za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*
     &    cwmass2**(-1) )
      amp(1,2) = amp(1,2) + gamZeq(1,2,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*s127**(-1) * ( two*za(p6,p8)*zb(p1,p2)*zba2(p4,p1
     &    ,p2,p7)*zba2(p5,p6,p8,p3)*s568**(-1) )
      amp(1,2) = amp(1,2) + gamZeq(1,2,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*s348**(-1) * ( two*za(p3,p8)*zb(p5,p2)*zba2(p1,p2
     &    ,p5,p6)*zba2(p4,p3,p8,p7)*s256**(-1) )
      amp(1,2) = amp(1,2) + gamZeq(2,2,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*s127**(-1)*s348**(-1) * ( two*za(p3,p8)*zb(p1,p2)
     &    *zba2(p4,p3,p8,p6)*zba2(p5,p1,p2,p7) )
      amp(1,2) = amp(1,2) + gameW(2)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1)*s127**(-1) * ( two*za(p7,p1)*za(
     &    p3,p6)*zb(p1,p2)*zb(p4,p5)*zab2(p8,p3,p4,p1) + za(p7,p1)*za(
     &    p3,p6)*zb(p1,p2)*zb(p4,p5)*zab2(p8,p2,p7,p1) - two*za(p7,p1)
     &    *za(p6,p8)*zb(p1,p5)*zb(p1,p2)*zab2(p3,p5,p6,p4) + za(p7,p2)*
     &    za(p3,p6)*zb(p1,p2)*zb(p4,p5)*zab2(p8,p1,p7,p2) + two*za(p7,
     &    p2)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*zab2(p8,p3,p4,p2) + two*
     &    za(p7,p2)*za(p6,p8)*zb(p1,p2)*zb(p5,p2)*zab2(p3,p5,p6,p4) + 2.
     &    D0*za(p3,p8)*zb(p1,p2)*zba2(p4,p1,p2,p7)*zba2(p5,p3,p4,p6) )
      amp(1,2) = amp(1,2) + gameW(2)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1) * (  - za(p7,p8)*za(p3,p5)*zb(p1,
     &    p2)*zb(p4,p5)*zba2(p5,p3,p4,p6)*cwmass2**(-1) + za(p7,p8)*za(
     &    p3,p6)*zb(p1,p2)*zb(p3,p5)*zab2(p3,p5,p6,p4)*cwmass2**(-1) -
     &    za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*cwmass2**(-1)*s3456
     &     - za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*cwmass2**(-1)*
     &    twop34Dp1278 - za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p6)*zba2(
     &    p5,p3,p4,p6)*cwmass2**(-1) + za(p7,p8)*za(p4,p6)*zb(p1,p2)*
     &    zb(p4,p5)*zab2(p3,p5,p6,p4)*cwmass2**(-1) )
      amp(2,1)= + cxw**(-3)*propW17**(-1)*propW34**(-1)*propW3456**(-1)
     &  * ( za(p7,p8)*za(p3,p5)*zb(p7,p1)*zb(p6,p2)*zba2(p4,p3,p5,p7)*
     &    s178**(-1)*s345**(-1) - za(p7,p8)*za(p3,p5)*zb(p1,p8)*zb(p6,
     &    p2)*zba2(p4,p3,p5,p8)*s178**(-1)*s345**(-1) + half*za(p7
     &    ,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*cwmass2**(-1) )
      amp(2,1) = amp(2,1) + gamZee(1,1)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * (  - two*za(p7,p3)*za(p7,p8)*zb(p7,p1)*zb(p4,
     &    p6)*zba2(p2,p4,p6,p5)*s178**(-1)*s456**(-1) - za(p7,p8)*za(p3
     &    ,p5)*zb(p1,p2)*zb(p4,p6)*cwmass2**(-1) - two*za(p7,p8)*za(p3
     &    ,p8)*zb(p1,p8)*zb(p4,p6)*zba2(p2,p4,p6,p5)*s178**(-1)*
     &    s456**(-1) )
      amp(2,1) = amp(2,1) + gamZne(1,1)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * (  - two*za(p7,p8)*za(p3,p5)*zb(p7,p1)*zb(p4,
     &    p2)*zba2(p6,p3,p5,p7)*s178**(-1)*s356**(-1) + two*za(p7,p8)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p2)*zba2(p6,p3,p5,p8)*s178**(-1)*
     &    s356**(-1) + za(p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*
     &    cwmass2**(-1) )
      amp(2,1) = amp(2,1) + gamZeq(1,1,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1) * ( two*za(p7,p8)*zb(p4,p2)*zba2(p1,p7,p8,p5)*
     &    zba2(p6,p2,p4,p3)*s178**(-1)*s234**(-1) )
      amp(2,1) = amp(2,1) + gamZeq(2,1,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1) * ( two*za(p7,p8)*zb(p6,p2)*zba2(p1,p7,p8,p3)*
     &    zba2(p4,p2,p6,p5)*s178**(-1)*s256**(-1) + two*za(p5,p8)*zb(
     &    p4,p2)*zba2(p1,p2,p4,p3)*zba2(p6,p5,p8,p7)*s234**(-1)*
     &    s568**(-1) )
      amp(2,1) = amp(2,1) + gameW(1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1) * (  - two*za(p7,p5)*za(p7,p8)*
     &    zb(p7,p1)*zb(p6,p2)*zab2(p3,p5,p6,p4)*s178**(-1) + za(p7,p8)*
     &    za(p3,p5)*zb(p7,p1)*zb(p4,p6)*zab2(p7,p1,p8,p2)*s178**(-1) +
     &    two*za(p7,p8)*za(p3,p5)*zb(p7,p1)*zb(p4,p6)*zab2(p7,p3,p4,p2
     &    )*s178**(-1) - za(p7,p8)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*zab2(
     &    p8,p1,p7,p2)*s178**(-1) - two*za(p7,p8)*za(p3,p5)*zb(p1,p8)*
     &    zb(p4,p6)*zab2(p8,p3,p4,p2)*s178**(-1) - za(p7,p8)*za(p3,p5)*
     &    zb(p1,p2)*zb(p3,p6)*zab2(p3,p5,p6,p4)*cwmass2**(-1) + za(p7,
     &    p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p5)*zab2(p5,p3,p4,p6)*
     &    cwmass2**(-1) + za(p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)*
     &    cwmass2**(-1)*s3456 + za(p7,p8)*za(p3,p5)*zb(p1,p2)*zb(p4,p6)
     &    *cwmass2**(-1)*twop34Dp1278 + za(p7,p8)*za(p3,p6)*zb(p1,p2)*
     &    zb(p4,p6)*zab2(p5,p3,p4,p6)*cwmass2**(-1) - za(p7,p8)*za(p4,
     &    p5)*zb(p1,p2)*zb(p4,p6)*zab2(p3,p5,p6,p4)*cwmass2**(-1) - two
     &    *za(p7,p8)*za(p5,p8)*zb(p1,p8)*zb(p6,p2)*zab2(p3,p5,p6,p4)*
     &    s178**(-1) )
      amp(2,1) = amp(2,1) + gameW(1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1) * (  - two*za(p7,p8)*zb(p4,p2)*
     &    zab2(p5,p3,p4,p6)*zba2(p1,p7,p8,p3)*s178**(-1) )
      amp(2,2)= + gamZee(1,2)*cxw**(-2)*propW17**(-1)*propW3456**(-1)
     &  * (  - two*za(p7,p3)*za(p7,p8)*zb(p7,p1)*zb(p4,p5)*zba2(p2,p4,
     &    p5,p6)*s178**(-1)*s456**(-1) - za(p7,p8)*za(p3,p6)*zb(p1,p2)*
     &    zb(p4,p5)*cwmass2**(-1) - two*za(p7,p8)*za(p3,p8)*zb(p1,p8)*
     &    zb(p4,p5)*zba2(p2,p4,p5,p6)*s178**(-1)*s456**(-1) )
      amp(2,2) = amp(2,2) + gamZne(1,2)*cxw**(-2)*propW17**(-1)*
     & propW3456**(-1) * (  - two*za(p7,p8)*za(p3,p6)*zb(p7,p1)*zb(p4,
     &    p2)*zba2(p5,p3,p6,p7)*s178**(-1)*s356**(-1) + two*za(p7,p8)*
     &    za(p3,p6)*zb(p1,p8)*zb(p4,p2)*zba2(p5,p3,p6,p8)*s178**(-1)*
     &    s356**(-1) + za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*
     &    cwmass2**(-1) )
      amp(2,2) = amp(2,2) + gamZeq(1,2,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1) * ( two*za(p7,p8)*zb(p4,p2)*zba2(p1,p7,p8,p6)*
     &    zba2(p5,p2,p4,p3)*s178**(-1)*s234**(-1) )
      amp(2,2) = amp(2,2) + gamZeq(2,2,1)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1) * ( two*za(p7,p8)*zb(p5,p2)*zba2(p1,p7,p8,p3)*
     &    zba2(p4,p2,p5,p6)*s178**(-1)*s256**(-1) + two*za(p6,p8)*zb(
     &    p4,p2)*zba2(p1,p2,p4,p3)*zba2(p5,p6,p8,p7)*s234**(-1)*
     &    s568**(-1) )
      amp(2,2) = amp(2,2) + gameW(2)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1) * (  - two*za(p7,p6)*za(p7,p8)*
     &    zb(p7,p1)*zb(p5,p2)*zab2(p3,p5,p6,p4)*s178**(-1) + za(p7,p8)*
     &    za(p3,p5)*zb(p1,p2)*zb(p4,p5)*zba2(p5,p3,p4,p6)*cwmass2**(-1)
     &     + za(p7,p8)*za(p3,p6)*zb(p7,p1)*zb(p4,p5)*zab2(p7,p1,p8,p2)*
     &    s178**(-1) + two*za(p7,p8)*za(p3,p6)*zb(p7,p1)*zb(p4,p5)*
     &    zab2(p7,p3,p4,p2)*s178**(-1) - za(p7,p8)*za(p3,p6)*zb(p1,p8)*
     &    zb(p4,p5)*zab2(p8,p1,p7,p2)*s178**(-1) - two*za(p7,p8)*za(p3
     &    ,p6)*zb(p1,p8)*zb(p4,p5)*zab2(p8,p3,p4,p2)*s178**(-1) - za(p7
     &    ,p8)*za(p3,p6)*zb(p1,p2)*zb(p3,p5)*zab2(p3,p5,p6,p4)*
     &    cwmass2**(-1) + za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)*
     &    cwmass2**(-1)*s3456 + za(p7,p8)*za(p3,p6)*zb(p1,p2)*zb(p4,p5)
     &    *cwmass2**(-1)*twop34Dp1278 + za(p7,p8)*za(p3,p6)*zb(p1,p2)*
     &    zb(p4,p6)*zba2(p5,p3,p4,p6)*cwmass2**(-1) - za(p7,p8)*za(p4,
     &    p6)*zb(p1,p2)*zb(p4,p5)*zab2(p3,p5,p6,p4)*cwmass2**(-1) - two
     &    *za(p7,p8)*za(p6,p8)*zb(p1,p8)*zb(p5,p2)*zab2(p3,p5,p6,p4)*
     &    s178**(-1) )
      amp(2,2) = amp(2,2) + gameW(2)*cxw**(-2)*propW17**(-1)*
     & propW34**(-1)*propW3456**(-1) * (  - two*za(p7,p8)*zb(p4,p2)*
     &    zba2(p1,p7,p8,p3)*zba2(p5,p3,p4,p6)*s178**(-1) )
      return
      end
