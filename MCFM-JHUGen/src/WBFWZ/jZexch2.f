      subroutine jZexch2(p1,p2,p3,p4,p5,p6,p7,p8,za,zb,amp)
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
      integer:: jdu1,jdu2,i1,i2,i3,i4,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,amp(2,2,2),rxw,gamZ56ne(2),
     & propw34,propz28,propZ56,propW3456,gameW(2),gamZ56ee(2),
     & gamZ56qe(2,2,2),gamZ28qq(2,2,2,2)
      real(dp):: qn,t3,t4,s134,s128,s345,s347,s356,s156,s456,
     & s278,s567,s3456,s28,s34,s56,
     & twop34Dp1278,twop56Dp1278,q3,q4,l3,l4
C     amp(jdu1,jdu2,h17,h28)
      parameter(qn=0d0)
C-----Begin statement functions
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s28=s(p2,p8)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s128=t3(p1,p2,p8)
      s134=t3(p1,p3,p4)
      s156=t3(p1,p5,p6)
      s278=t3(p2,p7,p8)
      s345=t3(p3,p4,p5)
      s347=t3(p3,p4,p7)
      s356=t3(p3,p5,p6)
      s456=t3(p4,p5,p6)
      s567=t3(p5,p6,p7)
      s3456=t4(p3,p4,p5,p6)
      twop34Dp1278=s56-s3456-s34
      twop56Dp1278=s34-s3456-s56

      propw34=s34-cwmass2
      propz56=s56-czmass2
      propz28=s28-czmass2
      propw3456=s3456-cwmass2

      q3=qn
      l3=ln
      q4=qe
      l4=le
      do jdu1=1,2
      gamZ56qe(jdu1,1,1)=Q(jdu1)*qe/s56+L(jdu1)*le/propZ56
      gamZ56qe(jdu1,1,2)=Q(jdu1)*qe/s56+L(jdu1)*re/propZ56
      gamZ56qe(jdu1,2,1)=Q(jdu1)*qe/s56+R(jdu1)*le/propZ56
      gamZ56qe(jdu1,2,2)=Q(jdu1)*qe/s56+R(jdu1)*re/propZ56

c      gamZ28qe(jdu1,1,1)=Q(jdu1)*q4/s28+L(jdu1)*l4/propZ28
c      gamZ28qe(jdu1,2,1)=Q(jdu1)*q4/s28+R(jdu1)*l4/propZ28

c      gamqW(jdu1,1)=Q(jdu1)/s56+L(jdu1)*rxw/propz56
c      gamqW(jdu1,2)=Q(jdu1)/s56+R(jdu1)*rxw/propz56

      do jdu2=1,2
      gamZ28qq(jdu1,jdu2,1,1)=
     & Q(jdu1)*Q(jdu2)/s28+L(jdu1)*L(jdu2)/propZ28
      gamZ28qq(jdu1,jdu2,1,2)=
     & Q(jdu1)*Q(jdu2)/s28+L(jdu1)*R(jdu2)/propZ28
      gamZ28qq(jdu1,jdu2,2,1)=
     & Q(jdu1)*Q(jdu2)/s28+R(jdu1)*L(jdu2)/propZ28
      gamZ28qq(jdu1,jdu2,2,2)=
     & Q(jdu1)*Q(jdu2)/s28+R(jdu1)*R(jdu2)/propZ28
      enddo

      enddo

      gameW(1)=qe/s56+le*rxw/propz56
      gameW(2)=qe/s56+re*rxw/propz56
c      gamnW(1)=      +ln*rxw/propz56
c      gamnW(2)=      +rn*rxw/propz56

      gamZ56ee(1)=q4*qe/s56+l4*le/propZ56
      gamZ56ee(2)=q4*qe/s56+l4*re/propZ56

      gamZ56ne(1)=q3*qe/s56+l3*le/propZ56
      gamZ56ne(2)=q3*qe/s56+l3*re/propZ56

      do jdu2=1,2
      amp(jdu2,1,1)= + gamZ56qe(1,1,1)*gamZ28qq(1,jdu2,1,1)*
     & propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p8)*zb(p1,p4)*zba2(p6,p1
     &    ,p4,p3)*zba2(p2,p7,p8,p5)*s134**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ56qe(1,1,1)*gamZ28qq(1,jdu2,1
     & ,1)*propW34**(-1) * ( 4.D0*za(p7,p5)*zb(p1,p4)*zba2(p6,p5,p7,p8)
     &    *zba2(p2,p1,p4,p3)*s134**(-1)*s567**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ56qe(1,1,1)*gamZ28qq(2,jdu2,1
     & ,1)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p5)*zb(p1,p2)*zba2(
     &    p6,p5,p7,p3)*zba2(p4,p1,p2,p8)*s567**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ56qe(2,1,1)*gamZ28qq(1,jdu2,1
     & ,1)*propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p8)*zb(p1,p6)*zba2(
     &    p4,p1,p6,p5)*zba2(p2,p7,p8,p3)*s156**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ56qe(2,1,1)*gamZ28qq(2,jdu2,1
     & ,1)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p2)*zba2(
     &    p6,p1,p2,p8)*zba2(p4,p3,p7,p5)*s347**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ56qe(2,1,1)*gamZ28qq(2,jdu2,1
     & ,1)*propW34**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p6)*zba2(p4,p3,p7,p8)
     &    *zba2(p2,p1,p6,p5)*s347**(-1)*s156**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*
     & cwmass2**(-1)*cxw**(-1)*propW34**(-1)*propW3456**(-1) * ( za(p7,
     &    p8)*za(p5,p3)**2*zb(p1,p2)*zb(p6,p3)*zb(p4,p5)*s345**(-1) -
     &    za(p7,p8)*za(p5,p3)*za(p3,p5)*zb(p1,p2)*zb(p5,p6)*zb(p4,p3)*
     &    s345**(-1) + za(p7,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p4,
     &    p3,p5,p4)*s345**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*cxw**(-1)*
     & propW34**(-1)*propW3456**(-1)*s278**(-1) * (  - 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p7,p2)*zb(p1,p6)*zba2(p4,p3,p5,p7)*s345**(-1) -
     &    2.D0*za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p8,p2)*zba2(p4,p3,p5,p8
     &    )*s345**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*gamZ56ee(1)*
     & cwmass2**(-1)*propW3456**(-1)*s456**(-1) * (  - 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p5,p4,p6,p5) + 2.D0*za(p7,
     &    p8)*za(p6,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)*zb(p4,p6) - 2.D0*
     &    za(p7,p8)*za(p6,p3)*za(p4,p5)*zb(p1,p2)*zb(p6,p4)**2 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*gamZ56ee(1)*
     & propW3456**(-1)*s278**(-1)*s456**(-1) * (  - 4.D0*za(p7,p3)*za(
     &    p7,p8)*zb(p7,p2)*zb(p6,p4)*zba2(p1,p4,p6,p5) + 4.D0*za(p7,p8)
     &    *za(p3,p8)*zb(p6,p4)*zb(p8,p2)*zba2(p1,p4,p6,p5) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*gamZ56ne(1)*
     & cwmass2**(-1)*propW3456**(-1)*s356**(-1) * ( 2.D0*za(p7,p8)*za(
     &    p5,p3)**2*zb(p1,p2)*zb(p6,p5)*zb(p3,p4) + 2.D0*za(p7,p8)*za(
     &    p5,p3)*za(p3,p5)*zb(p1,p2)*zb(p5,p4)*zb(p6,p3) + 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p3,p5,p6) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*gamZ56ne(1)*
     & propW3456**(-1)*s278**(-1)*s356**(-1) * ( 4.D0*za(p7,p8)*za(p5,
     &    p3)*zb(p7,p2)*zb(p1,p4)*zba2(p6,p3,p5,p7) + 4.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p1,p4)*zb(p8,p2)*zba2(p6,p3,p5,p8) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*gameW(1)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1) * (  - za(p7,p8)*za(
     &    p5,p3)*zb(p1,p2)*zb(p6,p4)*twop56Dp1278 + za(p7,p8)*za(p5,p3)
     &    *zb(p1,p2)*zb(p6,p4)*twop34Dp1278 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(1,jdu2,1,1)*gameW(1)*
     & propW34**(-1)*propW3456**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*za(
     &    p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,p8)
     &    *za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p7
     &    ,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1) - 2.D0*
     &    za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1) + 4.
     &    D0*za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,p5) -
     &    4.D0*za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,p8,p3)
     &     )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*
     & cwmass2**(-1)*cxw**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1)
     &  * ( za(p1,p8)*za(p5,p3)**2*zb(p1,p2)*zb(p6,p3)*zb(p4,p5)*zab2(
     &    p7,p5,p6,p1)*s345**(-1) + za(p1,p8)*za(p5,p3)**2*zb(p1,p2)*
     &    zb(p6,p3)*zb(p4,p5)*zab2(p7,p3,p4,p1)*s345**(-1) - za(p1,p8)*
     &    za(p5,p3)*za(p3,p5)*zb(p1,p2)*zb(p5,p6)*zb(p4,p3)*zab2(p7,p5,
     &    p6,p1)*s345**(-1) - za(p1,p8)*za(p5,p3)*za(p3,p5)*zb(p1,p2)*
     &    zb(p5,p6)*zb(p4,p3)*zab2(p7,p3,p4,p1)*s345**(-1) + za(p1,p8)*
     &    za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)*zba2(p4,p3,p5
     &    ,p4)*s345**(-1) + za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p1)*zba2(p4,p3,p5,p4)*s345**(-1) - za(p5,p3)**2
     &    *za(p8,p2)*zb(p1,p2)*zb(p6,p3)*zb(p4,p5)*zab2(p7,p5,p6,p2)*
     &    s345**(-1) - za(p5,p3)**2*za(p8,p2)*zb(p1,p2)*zb(p6,p3)*zb(p4
     &    ,p5)*zab2(p7,p3,p4,p2)*s345**(-1) + za(p5,p3)*za(p3,p5)*za(p8
     &    ,p2)*zb(p1,p2)*zb(p5,p6)*zb(p4,p3)*zab2(p7,p5,p6,p2)*
     &    s345**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*
     & cwmass2**(-1)*cxw**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1)
     &  * ( za(p5,p3)*za(p3,p5)*za(p8,p2)*zb(p1,p2)*zb(p5,p6)*zb(p4,p3)
     &    *zab2(p7,p3,p4,p2)*s345**(-1) - za(p5,p3)*za(p8,p2)*zb(p1,p2)
     &    *zb(p6,p4)*zab2(p7,p5,p6,p2)*zba2(p4,p3,p5,p4)*s345**(-1) -
     &    za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p2)*
     &    zba2(p4,p3,p5,p4)*s345**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*cxw**(-1)*
     & propW34**(-1)*propW3456**(-1)*s128**(-1) * ( 2.D0*za(p1,p8)*za(
     &    p5,p3)*zb(p1,p6)*zb(p1,p2)*zba2(p4,p3,p5,p7)*s345**(-1) + 2.D0
     &    *za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p2)*zba2(p4,p3,p5,p7)*
     &    s345**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gamZ56ee(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * (  - 2.D0*
     &    za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)*
     &    zba2(p5,p4,p6,p5) - 2.D0*za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,
     &    p4)*zab2(p7,p3,p4,p1)*zba2(p5,p4,p6,p5) + 2.D0*za(p1,p8)*za(
     &    p6,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)*zb(p4,p6)*zab2(p7,p5,p6,
     &    p1) + 2.D0*za(p1,p8)*za(p6,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)*
     &    zb(p4,p6)*zab2(p7,p3,p4,p1) - 2.D0*za(p1,p8)*za(p6,p3)*za(p4,
     &    p5)*zb(p1,p2)*zb(p6,p4)**2*zab2(p7,p5,p6,p1) - 2.D0*za(p1,p8)
     &    *za(p6,p3)*za(p4,p5)*zb(p1,p2)*zb(p6,p4)**2*zab2(p7,p3,p4,p1)
     &     + 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6
     &    ,p2)*zba2(p5,p4,p6,p5) + 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*
     &    zb(p6,p4)*zab2(p7,p3,p4,p2)*zba2(p5,p4,p6,p5) - 2.D0*za(p6,p5
     &    )*za(p3,p4)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zb(p4,p6)*zab2(p7,
     &    p5,p6,p2) - 2.D0*za(p6,p5)*za(p3,p4)*za(p8,p2)*zb(p1,p2)*zb(
     &    p6,p4)*zb(p4,p6)*zab2(p7,p3,p4,p2) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gamZ56ee(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * ( 2.D0*za(
     &    p6,p3)*za(p4,p5)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)**2*zab2(p7,p5,
     &    p6,p2) + 2.D0*za(p6,p3)*za(p4,p5)*za(p8,p2)*zb(p1,p2)*zb(p6,
     &    p4)**2*zab2(p7,p3,p4,p2) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gamZ56ee(1)*
     & propW3456**(-1)*s128**(-1)*s456**(-1) * ( 4.D0*za(p7,p3)*za(p1,
     &    p8)*zb(p1,p2)*zb(p6,p4)*zba2(p1,p4,p6,p5) - 4.D0*za(p7,p3)*
     &    za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zba2(p2,p4,p6,p5) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gamZ56ne(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * ( 2.D0*za(
     &    p1,p8)*za(p5,p3)**2*zb(p1,p2)*zb(p6,p5)*zb(p3,p4)*zab2(p7,p5,
     &    p6,p1) + 2.D0*za(p1,p8)*za(p5,p3)**2*zb(p1,p2)*zb(p6,p5)*zb(
     &    p3,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p1,p8)*za(p5,p3)*za(p3,p5)
     &    *zb(p1,p2)*zb(p5,p4)*zb(p6,p3)*zab2(p7,p5,p6,p1) + 2.D0*za(p1
     &    ,p8)*za(p5,p3)*za(p3,p5)*zb(p1,p2)*zb(p5,p4)*zb(p6,p3)*zab2(
     &    p7,p3,p4,p1) + 2.D0*za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*
     &    zab2(p7,p5,p6,p1)*zba2(p6,p3,p5,p6) + 2.D0*za(p1,p8)*za(p5,p3
     &    )*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1)*zba2(p6,p3,p5,p6) - 2.
     &    D0*za(p5,p3)**2*za(p8,p2)*zb(p1,p2)*zb(p6,p5)*zb(p3,p4)*zab2(
     &    p7,p5,p6,p2) - 2.D0*za(p5,p3)**2*za(p8,p2)*zb(p1,p2)*zb(p6,p5
     &    )*zb(p3,p4)*zab2(p7,p3,p4,p2) - 2.D0*za(p5,p3)*za(p3,p5)*za(
     &    p8,p2)*zb(p1,p2)*zb(p5,p4)*zb(p6,p3)*zab2(p7,p5,p6,p2) - 2.D0
     &    *za(p5,p3)*za(p3,p5)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zb(p6,p3)*
     &    zab2(p7,p3,p4,p2) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gamZ56ne(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * (  - 2.D0*
     &    za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2)*
     &    zba2(p6,p3,p5,p6) - 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,
     &    p4)*zab2(p7,p3,p4,p2)*zba2(p6,p3,p5,p6) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gamZ56ne(1)*
     & propW3456**(-1)*s128**(-1)*s356**(-1) * (  - 4.D0*za(p1,p8)*za(
     &    p5,p3)*zb(p1,p4)*zb(p1,p2)*zba2(p6,p3,p5,p7) - 4.D0*za(p5,p3)
     &    *za(p8,p2)*zb(p1,p2)*zb(p4,p2)*zba2(p6,p3,p5,p7) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gameW(1)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1) * (  -
     &    za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)*
     &    twop56Dp1278 + za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(
     &    p7,p5,p6,p1)*twop34Dp1278 - za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(
     &    p6,p4)*zab2(p7,p3,p4,p1)*twop56Dp1278 + za(p1,p8)*za(p5,p3)*
     &    zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1)*twop34Dp1278 + za(p5,p3
     &    )*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2)*
     &    twop56Dp1278 - za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(
     &    p7,p5,p6,p2)*twop34Dp1278 + za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(
     &    p6,p4)*zab2(p7,p3,p4,p2)*twop56Dp1278 - za(p5,p3)*za(p8,p2)*
     &    zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p2)*twop34Dp1278 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qq(2,jdu2,1,1)*gameW(1)*
     & propW34**(-1)*propW3456**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5)*
     &    zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8) + 4.D0*za(p7,p3
     &    )*zb(p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8) - 2.D0*za(p1,
     &    p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0*
     &    za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.
     &    D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2)
     &     - 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4
     &    ,p2) )
      amp(jdu2,1,2)= + gamZ56qe(1,1,2)*gamZ28qq(1,jdu2,1,1)*
     & propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p8)*zb(p1,p4)*zba2(p5,p1
     &    ,p4,p3)*zba2(p2,p7,p8,p6)*s134**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ56qe(1,1,2)*gamZ28qq(1,jdu2,1
     & ,1)*propW34**(-1) * ( 4.D0*za(p7,p6)*zb(p1,p4)*zba2(p5,p6,p7,p8)
     &    *zba2(p2,p1,p4,p3)*s134**(-1)*s567**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ56qe(1,1,2)*gamZ28qq(2,jdu2,1
     & ,1)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p6)*zb(p1,p2)*zba2(
     &    p5,p6,p7,p3)*zba2(p4,p1,p2,p8)*s567**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ56qe(2,1,2)*gamZ28qq(1,jdu2,1
     & ,1)*propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p8)*zb(p1,p5)*zba2(
     &    p4,p1,p5,p6)*zba2(p2,p7,p8,p3)*s156**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ56qe(2,1,2)*gamZ28qq(2,jdu2,1
     & ,1)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p2)*zba2(
     &    p5,p1,p2,p8)*zba2(p4,p3,p7,p6)*s347**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ56qe(2,1,2)*gamZ28qq(2,jdu2,1
     & ,1)*propW34**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p5)*zba2(p4,p3,p7,p8)
     &    *zba2(p2,p1,p5,p6)*s347**(-1)*s156**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(1,jdu2,1,1)*gamZ56ee(2)*
     & cwmass2**(-1)*propW3456**(-1)*s456**(-1) * ( 2.D0*za(p7,p8)*za(
     &    p5,p6)*za(p3,p4)*zb(p1,p2)*zb(p5,p4)*zb(p4,p5) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*za(p4,p6)*zb(p1,p2)*zb(p5,p4)**2 - 2.D0*za(p7,
     &    p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*zba2(p6,p4,p5,p6) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(1,jdu2,1,1)*gamZ56ee(2)*
     & propW3456**(-1)*s278**(-1)*s456**(-1) * (  - 4.D0*za(p7,p3)*za(
     &    p7,p8)*zb(p7,p2)*zb(p5,p4)*zba2(p1,p4,p5,p6) + 4.D0*za(p7,p8)
     &    *za(p3,p8)*zb(p5,p4)*zb(p8,p2)*zba2(p1,p4,p5,p6) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(1,jdu2,1,1)*gamZ56ne(2)*
     & cwmass2**(-1)*propW3456**(-1)*s356**(-1) * ( 2.D0*za(p7,p8)*za(
     &    p6,p3)**2*zb(p1,p2)*zb(p5,p6)*zb(p3,p4) + 2.D0*za(p7,p8)*za(
     &    p6,p3)*za(p3,p6)*zb(p1,p2)*zb(p5,p3)*zb(p6,p4) + 2.D0*za(p7,
     &    p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*zba2(p5,p3,p6,p5) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(1,jdu2,1,1)*gamZ56ne(2)*
     & propW3456**(-1)*s278**(-1)*s356**(-1) * ( 4.D0*za(p7,p8)*za(p6,
     &    p3)*zb(p7,p2)*zb(p1,p4)*zba2(p5,p3,p6,p7) + 4.D0*za(p7,p8)*
     &    za(p6,p3)*zb(p1,p4)*zb(p8,p2)*zba2(p5,p3,p6,p8) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(1,jdu2,1,1)*gameW(2)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1) * (  - za(p7,p8)*za(
     &    p6,p3)*zb(p1,p2)*zb(p5,p4)*twop56Dp1278 + za(p7,p8)*za(p6,p3)
     &    *zb(p1,p2)*zb(p5,p4)*twop34Dp1278 )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(1,jdu2,1,1)*gameW(2)*
     & propW34**(-1)*propW3456**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*za(
     &    p6,p3)*zb(p7,p2)*zb(p5,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,p8)
     &    *za(p6,p3)*zb(p7,p2)*zb(p5,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p7
     &    ,p8)*za(p6,p3)*zb(p5,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1) - 2.D0*
     &    za(p7,p8)*za(p6,p3)*zb(p5,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1) + 4.
     &    D0*za(p7,p8)*zb(p1,p5)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,p6) -
     &    4.D0*za(p7,p8)*zb(p1,p4)*zba2(p5,p3,p4,p6)*zba2(p2,p7,p8,p3)
     &     )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gamZ56ee(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * ( 2.D0*za(
     &    p1,p8)*za(p5,p6)*za(p3,p4)*zb(p1,p2)*zb(p5,p4)*zb(p4,p5)*
     &    zab2(p7,p5,p6,p1) + 2.D0*za(p1,p8)*za(p5,p6)*za(p3,p4)*zb(p1,
     &    p2)*zb(p5,p4)*zb(p4,p5)*zab2(p7,p3,p4,p1) - 2.D0*za(p1,p8)*
     &    za(p5,p3)*za(p4,p6)*zb(p1,p2)*zb(p5,p4)**2*zab2(p7,p5,p6,p1)
     &     - 2.D0*za(p1,p8)*za(p5,p3)*za(p4,p6)*zb(p1,p2)*zb(p5,p4)**2*
     &    zab2(p7,p3,p4,p1) - 2.D0*za(p1,p8)*za(p6,p3)*zb(p1,p2)*zb(p5,
     &    p4)*zab2(p7,p5,p6,p1)*zba2(p6,p4,p5,p6) - 2.D0*za(p1,p8)*za(
     &    p6,p3)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p3,p4,p1)*zba2(p6,p4,p5,p6
     &    ) - 2.D0*za(p5,p6)*za(p3,p4)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*
     &    zb(p4,p5)*zab2(p7,p5,p6,p2) - 2.D0*za(p5,p6)*za(p3,p4)*za(p8,
     &    p2)*zb(p1,p2)*zb(p5,p4)*zb(p4,p5)*zab2(p7,p3,p4,p2) + 2.D0*
     &    za(p5,p3)*za(p4,p6)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)**2*zab2(p7,
     &    p5,p6,p2) + 2.D0*za(p5,p3)*za(p4,p6)*za(p8,p2)*zb(p1,p2)*zb(
     &    p5,p4)**2*zab2(p7,p3,p4,p2) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gamZ56ee(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * ( 2.D0*za(
     &    p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p5,p6,p2)*zba2(
     &    p6,p4,p5,p6) + 2.D0*za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*
     &    zab2(p7,p3,p4,p2)*zba2(p6,p4,p5,p6) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gamZ56ee(2)*
     & propW3456**(-1)*s128**(-1)*s456**(-1) * ( 4.D0*za(p7,p3)*za(p1,
     &    p8)*zb(p1,p2)*zb(p5,p4)*zba2(p1,p4,p5,p6) - 4.D0*za(p7,p3)*
     &    za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zba2(p2,p4,p5,p6) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gamZ56ne(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * ( 2.D0*za(
     &    p1,p8)*za(p6,p3)**2*zb(p1,p2)*zb(p5,p6)*zb(p3,p4)*zab2(p7,p5,
     &    p6,p1) + 2.D0*za(p1,p8)*za(p6,p3)**2*zb(p1,p2)*zb(p5,p6)*zb(
     &    p3,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p1,p8)*za(p6,p3)*za(p3,p6)
     &    *zb(p1,p2)*zb(p5,p3)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0*za(p1
     &    ,p8)*za(p6,p3)*za(p3,p6)*zb(p1,p2)*zb(p5,p3)*zb(p6,p4)*zab2(
     &    p7,p3,p4,p1) + 2.D0*za(p1,p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*
     &    zab2(p7,p5,p6,p1)*zba2(p5,p3,p6,p5) + 2.D0*za(p1,p8)*za(p6,p3
     &    )*zb(p1,p2)*zb(p5,p4)*zab2(p7,p3,p4,p1)*zba2(p5,p3,p6,p5) - 2.
     &    D0*za(p6,p3)**2*za(p8,p2)*zb(p1,p2)*zb(p5,p6)*zb(p3,p4)*zab2(
     &    p7,p5,p6,p2) - 2.D0*za(p6,p3)**2*za(p8,p2)*zb(p1,p2)*zb(p5,p6
     &    )*zb(p3,p4)*zab2(p7,p3,p4,p2) - 2.D0*za(p6,p3)*za(p3,p6)*za(
     &    p8,p2)*zb(p1,p2)*zb(p5,p3)*zb(p6,p4)*zab2(p7,p5,p6,p2) - 2.D0
     &    *za(p6,p3)*za(p3,p6)*za(p8,p2)*zb(p1,p2)*zb(p5,p3)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p2) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gamZ56ne(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * (  - 2.D0*
     &    za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p5,p6,p2)*
     &    zba2(p5,p3,p6,p5) - 2.D0*za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,
     &    p4)*zab2(p7,p3,p4,p2)*zba2(p5,p3,p6,p5) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gamZ56ne(2)*
     & propW3456**(-1)*s128**(-1)*s356**(-1) * (  - 4.D0*za(p1,p8)*za(
     &    p6,p3)*zb(p1,p4)*zb(p1,p2)*zba2(p5,p3,p6,p7) - 4.D0*za(p6,p3)
     &    *za(p8,p2)*zb(p1,p2)*zb(p4,p2)*zba2(p5,p3,p6,p7) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gameW(2)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1) * (  -
     &    za(p1,p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p5,p6,p1)*
     &    twop56Dp1278 + za(p1,p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*zab2(
     &    p7,p5,p6,p1)*twop34Dp1278 - za(p1,p8)*za(p6,p3)*zb(p1,p2)*zb(
     &    p5,p4)*zab2(p7,p3,p4,p1)*twop56Dp1278 + za(p1,p8)*za(p6,p3)*
     &    zb(p1,p2)*zb(p5,p4)*zab2(p7,p3,p4,p1)*twop34Dp1278 + za(p6,p3
     &    )*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p5,p6,p2)*
     &    twop56Dp1278 - za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zab2(
     &    p7,p5,p6,p2)*twop34Dp1278 + za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(
     &    p5,p4)*zab2(p7,p3,p4,p2)*twop56Dp1278 - za(p6,p3)*za(p8,p2)*
     &    zb(p1,p2)*zb(p5,p4)*zab2(p7,p3,p4,p2)*twop34Dp1278 )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qq(2,jdu2,1,1)*gameW(2)*
     & propW34**(-1)*propW3456**(-1)*s128**(-1) * (  - 4.D0*za(p7,p6)*
     &    zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p5,p1,p2,p8) + 4.D0*za(p7,p3
     &    )*zb(p1,p2)*zba2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8) - 2.D0*za(p1,
     &    p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p5,p6,p1) + 2.D0*
     &    za(p1,p8)*za(p6,p3)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p3,p4,p1) + 2.
     &    D0*za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p5,p6,p2)
     &     - 2.D0*za(p6,p3)*za(p8,p2)*zb(p1,p2)*zb(p5,p4)*zab2(p7,p3,p4
     &    ,p2) )
      amp(jdu2,2,1)= + gamZ56qe(1,1,1)*gamZ28qq(1,jdu2,1,2)*
     & propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p2)*zb(p1,p4)*zba2(p6,p1
     &    ,p4,p3)*zba2(p8,p2,p7,p5)*s134**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ56qe(1,1,1)*gamZ28qq(1,jdu2,1
     & ,2)*propW34**(-1) * ( 4.D0*za(p7,p5)*zb(p1,p4)*zba2(p6,p5,p7,p2)
     &    *zba2(p8,p1,p4,p3)*s134**(-1)*s567**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ56qe(1,1,1)*gamZ28qq(2,jdu2,1
     & ,2)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p5)*zb(p1,p8)*zba2(
     &    p6,p5,p7,p3)*zba2(p4,p1,p8,p2)*s567**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ56qe(2,1,1)*gamZ28qq(1,jdu2,1
     & ,2)*propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p2)*zb(p1,p6)*zba2(
     &    p4,p1,p6,p5)*zba2(p8,p2,p7,p3)*s156**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ56qe(2,1,1)*gamZ28qq(2,jdu2,1
     & ,2)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p8)*zba2(
     &    p6,p1,p8,p2)*zba2(p4,p3,p7,p5)*s347**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ56qe(2,1,1)*gamZ28qq(2,jdu2,1
     & ,2)*propW34**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p6)*zba2(p4,p3,p7,p2)
     &    *zba2(p8,p1,p6,p5)*s347**(-1)*s156**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*
     & cwmass2**(-1)*cxw**(-1)*propW34**(-1)*propW3456**(-1) * ( za(p7,
     &    p2)*za(p5,p3)**2*zb(p1,p8)*zb(p6,p3)*zb(p4,p5)*s345**(-1) -
     &    za(p7,p2)*za(p5,p3)*za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zb(p4,p3)*
     &    s345**(-1) + za(p7,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zba2(p4,
     &    p3,p5,p4)*s345**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*cxw**(-1)*
     & propW34**(-1)*propW3456**(-1)*s278**(-1) * (  - 2.D0*za(p7,p2)*
     &    za(p5,p3)*zb(p7,p8)*zb(p1,p6)*zba2(p4,p3,p5,p7)*s345**(-1) +
     &    2.D0*za(p7,p2)*za(p5,p3)*zb(p1,p6)*zb(p8,p2)*zba2(p4,p3,p5,p2
     &    )*s345**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*gamZ56ee(1)*
     & cwmass2**(-1)*propW3456**(-1)*s456**(-1) * (  - 2.D0*za(p7,p2)*
     &    za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zba2(p5,p4,p6,p5) + 2.D0*za(p7,
     &    p2)*za(p6,p5)*za(p3,p4)*zb(p1,p8)*zb(p6,p4)*zb(p4,p6) - 2.D0*
     &    za(p7,p2)*za(p6,p3)*za(p4,p5)*zb(p1,p8)*zb(p6,p4)**2 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*gamZ56ee(1)*
     & propW3456**(-1)*s278**(-1)*s456**(-1) * (  - 4.D0*za(p7,p3)*za(
     &    p7,p2)*zb(p7,p8)*zb(p6,p4)*zba2(p1,p4,p6,p5) - 4.D0*za(p7,p2)
     &    *za(p3,p2)*zb(p6,p4)*zb(p8,p2)*zba2(p1,p4,p6,p5) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*gamZ56ne(1)*
     & cwmass2**(-1)*propW3456**(-1)*s356**(-1) * ( 2.D0*za(p7,p2)*za(
     &    p5,p3)**2*zb(p1,p8)*zb(p6,p5)*zb(p3,p4) + 2.D0*za(p7,p2)*za(
     &    p5,p3)*za(p3,p5)*zb(p1,p8)*zb(p5,p4)*zb(p6,p3) + 2.D0*za(p7,
     &    p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zba2(p6,p3,p5,p6) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*gamZ56ne(1)*
     & propW3456**(-1)*s278**(-1)*s356**(-1) * ( 4.D0*za(p7,p2)*za(p5,
     &    p3)*zb(p7,p8)*zb(p1,p4)*zba2(p6,p3,p5,p7) - 4.D0*za(p7,p2)*
     &    za(p5,p3)*zb(p1,p4)*zb(p8,p2)*zba2(p6,p3,p5,p2) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*gameW(1)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1) * (  - za(p7,p2)*za(
     &    p5,p3)*zb(p1,p8)*zb(p6,p4)*twop56Dp1278 + za(p7,p2)*za(p5,p3)
     &    *zb(p1,p8)*zb(p6,p4)*twop34Dp1278 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(1,jdu2,1,2)*gameW(1)*
     & propW34**(-1)*propW3456**(-1)*s278**(-1) * ( 2.D0*za(p7,p2)*za(
     &    p5,p3)*zb(p7,p8)*zb(p6,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,p2)
     &    *za(p5,p3)*zb(p7,p8)*zb(p6,p4)*zab2(p7,p3,p4,p1) - 2.D0*za(p7
     &    ,p2)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p2,p5,p6,p1) + 2.D0*
     &    za(p7,p2)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p2,p3,p4,p1) + 4.
     &    D0*za(p7,p2)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p5) -
     &    4.D0*za(p7,p2)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p8,p2,p7,p3)
     &     )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*
     & cwmass2**(-1)*cxw**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1)
     &  * ( za(p1,p2)*za(p5,p3)**2*zb(p1,p8)*zb(p6,p3)*zb(p4,p5)*zab2(
     &    p7,p5,p6,p1)*s345**(-1) + za(p1,p2)*za(p5,p3)**2*zb(p1,p8)*
     &    zb(p6,p3)*zb(p4,p5)*zab2(p7,p3,p4,p1)*s345**(-1) - za(p1,p2)*
     &    za(p5,p3)*za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zb(p4,p3)*zab2(p7,p5,
     &    p6,p1)*s345**(-1) - za(p1,p2)*za(p5,p3)*za(p3,p5)*zb(p1,p8)*
     &    zb(p5,p6)*zb(p4,p3)*zab2(p7,p3,p4,p1)*s345**(-1) + za(p1,p2)*
     &    za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p1)*zba2(p4,p3,p5
     &    ,p4)*s345**(-1) + za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p1)*zba2(p4,p3,p5,p4)*s345**(-1) + za(p5,p3)**2
     &    *za(p8,p2)*zb(p1,p8)*zb(p6,p3)*zb(p4,p5)*zab2(p7,p5,p6,p8)*
     &    s345**(-1) + za(p5,p3)**2*za(p8,p2)*zb(p1,p8)*zb(p6,p3)*zb(p4
     &    ,p5)*zab2(p7,p3,p4,p8)*s345**(-1) - za(p5,p3)*za(p3,p5)*za(p8
     &    ,p2)*zb(p1,p8)*zb(p5,p6)*zb(p4,p3)*zab2(p7,p5,p6,p8)*
     &    s345**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*
     & cwmass2**(-1)*cxw**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1)
     &  * (  - za(p5,p3)*za(p3,p5)*za(p8,p2)*zb(p1,p8)*zb(p5,p6)*zb(p4,
     &    p3)*zab2(p7,p3,p4,p8)*s345**(-1) + za(p5,p3)*za(p8,p2)*zb(p1,
     &    p8)*zb(p6,p4)*zab2(p7,p5,p6,p8)*zba2(p4,p3,p5,p4)*s345**(-1)
     &     + za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p3,p4,p8)*
     &    zba2(p4,p3,p5,p4)*s345**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*cxw**(-1)*
     & propW34**(-1)*propW3456**(-1)*s128**(-1) * ( 2.D0*za(p1,p2)*za(
     &    p5,p3)*zb(p1,p6)*zb(p1,p8)*zba2(p4,p3,p5,p7)*s345**(-1) - 2.D0
     &    *za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p8)*zba2(p4,p3,p5,p7)*
     &    s345**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gamZ56ee(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * (  - 2.D0*
     &    za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p1)*
     &    zba2(p5,p4,p6,p5) - 2.D0*za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,
     &    p4)*zab2(p7,p3,p4,p1)*zba2(p5,p4,p6,p5) + 2.D0*za(p1,p2)*za(
     &    p6,p5)*za(p3,p4)*zb(p1,p8)*zb(p6,p4)*zb(p4,p6)*zab2(p7,p5,p6,
     &    p1) + 2.D0*za(p1,p2)*za(p6,p5)*za(p3,p4)*zb(p1,p8)*zb(p6,p4)*
     &    zb(p4,p6)*zab2(p7,p3,p4,p1) - 2.D0*za(p1,p2)*za(p6,p3)*za(p4,
     &    p5)*zb(p1,p8)*zb(p6,p4)**2*zab2(p7,p5,p6,p1) - 2.D0*za(p1,p2)
     &    *za(p6,p3)*za(p4,p5)*zb(p1,p8)*zb(p6,p4)**2*zab2(p7,p3,p4,p1)
     &     - 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6
     &    ,p8)*zba2(p5,p4,p6,p5) - 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p8)*
     &    zb(p6,p4)*zab2(p7,p3,p4,p8)*zba2(p5,p4,p6,p5) + 2.D0*za(p6,p5
     &    )*za(p3,p4)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zb(p4,p6)*zab2(p7,
     &    p5,p6,p8) + 2.D0*za(p6,p5)*za(p3,p4)*za(p8,p2)*zb(p1,p8)*zb(
     &    p6,p4)*zb(p4,p6)*zab2(p7,p3,p4,p8) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gamZ56ee(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * (  - 2.D0*
     &    za(p6,p3)*za(p4,p5)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)**2*zab2(p7,
     &    p5,p6,p8) - 2.D0*za(p6,p3)*za(p4,p5)*za(p8,p2)*zb(p1,p8)*zb(
     &    p6,p4)**2*zab2(p7,p3,p4,p8) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gamZ56ee(1)*
     & propW3456**(-1)*s128**(-1)*s456**(-1) * ( 4.D0*za(p7,p3)*za(p1,
     &    p2)*zb(p1,p8)*zb(p6,p4)*zba2(p1,p4,p6,p5) + 4.D0*za(p7,p3)*
     &    za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zba2(p8,p4,p6,p5) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gamZ56ne(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * ( 2.D0*za(
     &    p1,p2)*za(p5,p3)**2*zb(p1,p8)*zb(p6,p5)*zb(p3,p4)*zab2(p7,p5,
     &    p6,p1) + 2.D0*za(p1,p2)*za(p5,p3)**2*zb(p1,p8)*zb(p6,p5)*zb(
     &    p3,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p1,p2)*za(p5,p3)*za(p3,p5)
     &    *zb(p1,p8)*zb(p5,p4)*zb(p6,p3)*zab2(p7,p5,p6,p1) + 2.D0*za(p1
     &    ,p2)*za(p5,p3)*za(p3,p5)*zb(p1,p8)*zb(p5,p4)*zb(p6,p3)*zab2(
     &    p7,p3,p4,p1) + 2.D0*za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*
     &    zab2(p7,p5,p6,p1)*zba2(p6,p3,p5,p6) + 2.D0*za(p1,p2)*za(p5,p3
     &    )*zb(p1,p8)*zb(p6,p4)*zab2(p7,p3,p4,p1)*zba2(p6,p3,p5,p6) + 2.
     &    D0*za(p5,p3)**2*za(p8,p2)*zb(p1,p8)*zb(p6,p5)*zb(p3,p4)*zab2(
     &    p7,p5,p6,p8) + 2.D0*za(p5,p3)**2*za(p8,p2)*zb(p1,p8)*zb(p6,p5
     &    )*zb(p3,p4)*zab2(p7,p3,p4,p8) + 2.D0*za(p5,p3)*za(p3,p5)*za(
     &    p8,p2)*zb(p1,p8)*zb(p5,p4)*zb(p6,p3)*zab2(p7,p5,p6,p8) + 2.D0
     &    *za(p5,p3)*za(p3,p5)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zb(p6,p3)*
     &    zab2(p7,p3,p4,p8) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gamZ56ne(1)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * ( 2.D0*za(
     &    p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p8)*zba2(
     &    p6,p3,p5,p6) + 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p8)*zba2(p6,p3,p5,p6) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gamZ56ne(1)*
     & propW3456**(-1)*s128**(-1)*s356**(-1) * (  - 4.D0*za(p1,p2)*za(
     &    p5,p3)*zb(p1,p4)*zb(p1,p8)*zba2(p6,p3,p5,p7) + 4.D0*za(p5,p3)
     &    *za(p8,p2)*zb(p1,p8)*zb(p4,p8)*zba2(p6,p3,p5,p7) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gameW(1)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1) * (  -
     &    za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p1)*
     &    twop56Dp1278 + za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zab2(
     &    p7,p5,p6,p1)*twop34Dp1278 - za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(
     &    p6,p4)*zab2(p7,p3,p4,p1)*twop56Dp1278 + za(p1,p2)*za(p5,p3)*
     &    zb(p1,p8)*zb(p6,p4)*zab2(p7,p3,p4,p1)*twop34Dp1278 - za(p5,p3
     &    )*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p8)*
     &    twop56Dp1278 + za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(
     &    p7,p5,p6,p8)*twop34Dp1278 - za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(
     &    p6,p4)*zab2(p7,p3,p4,p8)*twop56Dp1278 + za(p5,p3)*za(p8,p2)*
     &    zb(p1,p8)*zb(p6,p4)*zab2(p7,p3,p4,p8)*twop34Dp1278 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qq(2,jdu2,1,2)*gameW(1)*
     & propW34**(-1)*propW3456**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5)*
     &    zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p8,p2) + 4.D0*za(p7,p3
     &    )*zb(p1,p8)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p8,p2) - 2.D0*za(p1,
     &    p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0*
     &    za(p1,p2)*za(p5,p3)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p3,p4,p1) - 2.
     &    D0*za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p5,p6,p8)
     &     + 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p8)*zb(p6,p4)*zab2(p7,p3,p4
     &    ,p8) )
      amp(jdu2,2,2)= + gamZ56qe(1,1,2)*gamZ28qq(1,jdu2,1,2)*
     & propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p2)*zb(p1,p4)*zba2(p5,p1
     &    ,p4,p3)*zba2(p8,p2,p7,p6)*s134**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ56qe(1,1,2)*gamZ28qq(1,jdu2,1
     & ,2)*propW34**(-1) * ( 4.D0*za(p7,p6)*zb(p1,p4)*zba2(p5,p6,p7,p2)
     &    *zba2(p8,p1,p4,p3)*s134**(-1)*s567**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ56qe(1,1,2)*gamZ28qq(2,jdu2,1
     & ,2)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p6)*zb(p1,p8)*zba2(
     &    p5,p6,p7,p3)*zba2(p4,p1,p8,p2)*s567**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ56qe(2,1,2)*gamZ28qq(1,jdu2,1
     & ,2)*propW34**(-1)*s278**(-1) * ( 4.D0*za(p7,p2)*zb(p1,p5)*zba2(
     &    p4,p1,p5,p6)*zba2(p8,p2,p7,p3)*s156**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ56qe(2,1,2)*gamZ28qq(2,jdu2,1
     & ,2)*propW34**(-1)*s128**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p8)*zba2(
     &    p5,p1,p8,p2)*zba2(p4,p3,p7,p6)*s347**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ56qe(2,1,2)*gamZ28qq(2,jdu2,1
     & ,2)*propW34**(-1) * ( 4.D0*za(p7,p3)*zb(p1,p5)*zba2(p4,p3,p7,p2)
     &    *zba2(p8,p1,p5,p6)*s347**(-1)*s156**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(1,jdu2,1,2)*gamZ56ee(2)*
     & cwmass2**(-1)*propW3456**(-1)*s456**(-1) * ( 2.D0*za(p7,p2)*za(
     &    p5,p6)*za(p3,p4)*zb(p1,p8)*zb(p5,p4)*zb(p4,p5) - 2.D0*za(p7,
     &    p2)*za(p5,p3)*za(p4,p6)*zb(p1,p8)*zb(p5,p4)**2 - 2.D0*za(p7,
     &    p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*zba2(p6,p4,p5,p6) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(1,jdu2,1,2)*gamZ56ee(2)*
     & propW3456**(-1)*s278**(-1)*s456**(-1) * (  - 4.D0*za(p7,p3)*za(
     &    p7,p2)*zb(p7,p8)*zb(p5,p4)*zba2(p1,p4,p5,p6) - 4.D0*za(p7,p2)
     &    *za(p3,p2)*zb(p5,p4)*zb(p8,p2)*zba2(p1,p4,p5,p6) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(1,jdu2,1,2)*gamZ56ne(2)*
     & cwmass2**(-1)*propW3456**(-1)*s356**(-1) * ( 2.D0*za(p7,p2)*za(
     &    p6,p3)**2*zb(p1,p8)*zb(p5,p6)*zb(p3,p4) + 2.D0*za(p7,p2)*za(
     &    p6,p3)*za(p3,p6)*zb(p1,p8)*zb(p5,p3)*zb(p6,p4) + 2.D0*za(p7,
     &    p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*zba2(p5,p3,p6,p5) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(1,jdu2,1,2)*gamZ56ne(2)*
     & propW3456**(-1)*s278**(-1)*s356**(-1) * ( 4.D0*za(p7,p2)*za(p6,
     &    p3)*zb(p7,p8)*zb(p1,p4)*zba2(p5,p3,p6,p7) - 4.D0*za(p7,p2)*
     &    za(p6,p3)*zb(p1,p4)*zb(p8,p2)*zba2(p5,p3,p6,p2) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(1,jdu2,1,2)*gameW(2)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1) * (  - za(p7,p2)*za(
     &    p6,p3)*zb(p1,p8)*zb(p5,p4)*twop56Dp1278 + za(p7,p2)*za(p6,p3)
     &    *zb(p1,p8)*zb(p5,p4)*twop34Dp1278 )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(1,jdu2,1,2)*gameW(2)*
     & propW34**(-1)*propW3456**(-1)*s278**(-1) * ( 2.D0*za(p7,p2)*za(
     &    p6,p3)*zb(p7,p8)*zb(p5,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,p2)
     &    *za(p6,p3)*zb(p7,p8)*zb(p5,p4)*zab2(p7,p3,p4,p1) - 2.D0*za(p7
     &    ,p2)*za(p6,p3)*zb(p5,p4)*zb(p8,p2)*zab2(p2,p5,p6,p1) + 2.D0*
     &    za(p7,p2)*za(p6,p3)*zb(p5,p4)*zb(p8,p2)*zab2(p2,p3,p4,p1) + 4.
     &    D0*za(p7,p2)*zb(p1,p5)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p6) -
     &    4.D0*za(p7,p2)*zb(p1,p4)*zba2(p5,p3,p4,p6)*zba2(p8,p2,p7,p3)
     &     )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gamZ56ee(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * ( 2.D0*za(
     &    p1,p2)*za(p5,p6)*za(p3,p4)*zb(p1,p8)*zb(p5,p4)*zb(p4,p5)*
     &    zab2(p7,p5,p6,p1) + 2.D0*za(p1,p2)*za(p5,p6)*za(p3,p4)*zb(p1,
     &    p8)*zb(p5,p4)*zb(p4,p5)*zab2(p7,p3,p4,p1) - 2.D0*za(p1,p2)*
     &    za(p5,p3)*za(p4,p6)*zb(p1,p8)*zb(p5,p4)**2*zab2(p7,p5,p6,p1)
     &     - 2.D0*za(p1,p2)*za(p5,p3)*za(p4,p6)*zb(p1,p8)*zb(p5,p4)**2*
     &    zab2(p7,p3,p4,p1) - 2.D0*za(p1,p2)*za(p6,p3)*zb(p1,p8)*zb(p5,
     &    p4)*zab2(p7,p5,p6,p1)*zba2(p6,p4,p5,p6) - 2.D0*za(p1,p2)*za(
     &    p6,p3)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p3,p4,p1)*zba2(p6,p4,p5,p6
     &    ) + 2.D0*za(p5,p6)*za(p3,p4)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*
     &    zb(p4,p5)*zab2(p7,p5,p6,p8) + 2.D0*za(p5,p6)*za(p3,p4)*za(p8,
     &    p2)*zb(p1,p8)*zb(p5,p4)*zb(p4,p5)*zab2(p7,p3,p4,p8) - 2.D0*
     &    za(p5,p3)*za(p4,p6)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)**2*zab2(p7,
     &    p5,p6,p8) - 2.D0*za(p5,p3)*za(p4,p6)*za(p8,p2)*zb(p1,p8)*zb(
     &    p5,p4)**2*zab2(p7,p3,p4,p8) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gamZ56ee(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s456**(-1) * (  - 2.D0*
     &    za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p5,p6,p8)*
     &    zba2(p6,p4,p5,p6) - 2.D0*za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,
     &    p4)*zab2(p7,p3,p4,p8)*zba2(p6,p4,p5,p6) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gamZ56ee(2)*
     & propW3456**(-1)*s128**(-1)*s456**(-1) * ( 4.D0*za(p7,p3)*za(p1,
     &    p2)*zb(p1,p8)*zb(p5,p4)*zba2(p1,p4,p5,p6) + 4.D0*za(p7,p3)*
     &    za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zba2(p8,p4,p5,p6) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gamZ56ne(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * ( 2.D0*za(
     &    p1,p2)*za(p6,p3)**2*zb(p1,p8)*zb(p5,p6)*zb(p3,p4)*zab2(p7,p5,
     &    p6,p1) + 2.D0*za(p1,p2)*za(p6,p3)**2*zb(p1,p8)*zb(p5,p6)*zb(
     &    p3,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p1,p2)*za(p6,p3)*za(p3,p6)
     &    *zb(p1,p8)*zb(p5,p3)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0*za(p1
     &    ,p2)*za(p6,p3)*za(p3,p6)*zb(p1,p8)*zb(p5,p3)*zb(p6,p4)*zab2(
     &    p7,p3,p4,p1) + 2.D0*za(p1,p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*
     &    zab2(p7,p5,p6,p1)*zba2(p5,p3,p6,p5) + 2.D0*za(p1,p2)*za(p6,p3
     &    )*zb(p1,p8)*zb(p5,p4)*zab2(p7,p3,p4,p1)*zba2(p5,p3,p6,p5) + 2.
     &    D0*za(p6,p3)**2*za(p8,p2)*zb(p1,p8)*zb(p5,p6)*zb(p3,p4)*zab2(
     &    p7,p5,p6,p8) + 2.D0*za(p6,p3)**2*za(p8,p2)*zb(p1,p8)*zb(p5,p6
     &    )*zb(p3,p4)*zab2(p7,p3,p4,p8) + 2.D0*za(p6,p3)*za(p3,p6)*za(
     &    p8,p2)*zb(p1,p8)*zb(p5,p3)*zb(p6,p4)*zab2(p7,p5,p6,p8) + 2.D0
     &    *za(p6,p3)*za(p3,p6)*za(p8,p2)*zb(p1,p8)*zb(p5,p3)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p8) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gamZ56ne(2)*
     & cwmass2**(-1)*propW3456**(-1)*s128**(-1)*s356**(-1) * ( 2.D0*za(
     &    p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p5,p6,p8)*zba2(
     &    p5,p3,p6,p5) + 2.D0*za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*
     &    zab2(p7,p3,p4,p8)*zba2(p5,p3,p6,p5) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gamZ56ne(2)*
     & propW3456**(-1)*s128**(-1)*s356**(-1) * (  - 4.D0*za(p1,p2)*za(
     &    p6,p3)*zb(p1,p4)*zb(p1,p8)*zba2(p5,p3,p6,p7) + 4.D0*za(p6,p3)
     &    *za(p8,p2)*zb(p1,p8)*zb(p4,p8)*zba2(p5,p3,p6,p7) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gameW(2)*
     & cwmass2**(-1)*propW34**(-1)*propW3456**(-1)*s128**(-1) * (  -
     &    za(p1,p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p5,p6,p1)*
     &    twop56Dp1278 + za(p1,p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*zab2(
     &    p7,p5,p6,p1)*twop34Dp1278 - za(p1,p2)*za(p6,p3)*zb(p1,p8)*zb(
     &    p5,p4)*zab2(p7,p3,p4,p1)*twop56Dp1278 + za(p1,p2)*za(p6,p3)*
     &    zb(p1,p8)*zb(p5,p4)*zab2(p7,p3,p4,p1)*twop34Dp1278 - za(p6,p3
     &    )*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p5,p6,p8)*
     &    twop56Dp1278 + za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zab2(
     &    p7,p5,p6,p8)*twop34Dp1278 - za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(
     &    p5,p4)*zab2(p7,p3,p4,p8)*twop56Dp1278 + za(p6,p3)*za(p8,p2)*
     &    zb(p1,p8)*zb(p5,p4)*zab2(p7,p3,p4,p8)*twop34Dp1278 )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qq(2,jdu2,1,2)*gameW(2)*
     & propW34**(-1)*propW3456**(-1)*s128**(-1) * (  - 4.D0*za(p7,p6)*
     &    zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p5,p1,p8,p2) + 4.D0*za(p7,p3
     &    )*zb(p1,p8)*zba2(p5,p3,p4,p6)*zba2(p4,p1,p8,p2) - 2.D0*za(p1,
     &    p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p5,p6,p1) + 2.D0*
     &    za(p1,p2)*za(p6,p3)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p3,p4,p1) - 2.
     &    D0*za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p5,p6,p8)
     &     + 2.D0*za(p6,p3)*za(p8,p2)*zb(p1,p8)*zb(p5,p4)*zab2(p7,p3,p4
     &    ,p8) )
      enddo
      amp(:,:,:)=amp(:,:,:)/cxw

      return
      end
