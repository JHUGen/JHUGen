      subroutine amp2current(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
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
      integer:: jdu1,jdu2,h28,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,amp(2,2,2,2),rxw,
     & propw34,propw56,propz28,propz3456,gam28(2,2,2,2),
     & gamn3456(2,2,2),game3456(2,2,2),gamv(2,2),gamv17(2,2),
     & gamv28(2,2),game28(2,2),gamn28(2,2),propW2348,propW2568,
     & propZ17
      real(dp):: qn,t3,t4,s134,s128,s345,s347,s356,s156,s456,
     & s278,s346,s567,s3456,s28,s34,s56,s238,s248,s258,s268,s17,
     & p34Dp28,p34Dp56,p56Dp28,p17Dp28,p17Dp34,p17Dp56
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
      s28=s(i2,i8)
      s17=s(i1,i7)
      s34=s(i3,i4)
      s56=s(i5,i6)
      s128=t3(i1,i2,i8)
      s134=t3(i1,i3,i4)
      s156=t3(i1,i5,i6)
      s278=t3(i2,i7,i8)
      s345=t3(i3,i4,i5)
      s346=t3(i3,i4,i6)
      s347=t3(i3,i4,i7)
      s356=t3(i3,i5,i6)
      s456=t3(i4,i5,i6)
      s567=t3(i5,i6,i7)
      s238=t3(i2,i3,i8)
      s248=t3(i2,i4,i8)
      s258=t3(i2,i5,i8)
      s268=t3(i2,i6,i8)
      s3456=t4(i3,i4,i5,i6)
      p34Dp28=(s(i3,i2)+s(i3,i8)+s(i4,i2)+s(i4,i8))/2d0
      p34Dp56=(s(i3,i5)+s(i3,i6)+s(i4,i5)+s(i4,i6))/2d0
      p56Dp28=(s(i5,i2)+s(i5,i8)+s(i6,i2)+s(i6,i8))/2d0
      p17Dp28=(s(i1,i2)+s(i1,i8)+s(i7,i2)+s(i7,i8))/2d0
      p17Dp34=(s(i1,i3)+s(i1,i4)+s(i7,i3)+s(i7,i4))/2d0
      p17Dp56=(s(i1,i5)+s(i1,i6)+s(i7,i5)+s(i7,i6))/2d0

      propw34=s34-cwmass2
      propw56=s56-cwmass2
      propz28=s28-czmass2
      propz17=s17-czmass2
      propz3456=s3456-czmass2
      propw2348=t4(i2,i3,i4,i8)-cwmass2
      propw2568=t4(i2,i5,i6,i8)-cwmass2

      do jdu1=1,2
      game3456(jdu1,1,1)=Q(jdu1)*qe/s3456+L(jdu1)*le/propZ3456
      game3456(jdu1,1,2)=Q(jdu1)*qe/s3456+L(jdu1)*re/propZ3456
      game3456(jdu1,2,1)=Q(jdu1)*qe/s3456+R(jdu1)*le/propZ3456
      game3456(jdu1,2,2)=Q(jdu1)*qe/s3456+R(jdu1)*re/propZ3456

      gamV(jdu1,1)=Q(jdu1)/s3456+L(jdu1)*rxw/propZ3456
      gamV(jdu1,2)=Q(jdu1)/s3456+R(jdu1)*rxw/propZ3456

      gamV17(jdu1,1)=Q(jdu1)/s17+L(jdu1)*rxw/propZ17
      gamV17(jdu1,2)=Q(jdu1)/s17+R(jdu1)*rxw/propZ17

      gamV28(jdu1,1)=Q(jdu1)/s28+L(jdu1)*rxw/propZ28
      gamV28(jdu1,2)=Q(jdu1)/s28+R(jdu1)*rxw/propZ28

      gamn3456(jdu1,1,1)=Q(jdu1)*qn/s3456+L(jdu1)*ln/propZ3456
      gamn3456(jdu1,1,2)=Q(jdu1)*qn/s3456+L(jdu1)*rn/propZ3456
      gamn3456(jdu1,2,1)=Q(jdu1)*qn/s3456+R(jdu1)*ln/propZ3456
      gamn3456(jdu1,2,2)=Q(jdu1)*qn/s3456+R(jdu1)*rn/propZ3456

      game28(jdu1,1)=Q(jdu1)*qe/s28+L(jdu1)*le/propZ28
      game28(jdu1,2)=Q(jdu1)*qe/s28+R(jdu1)*le/propZ28
      gamn28(jdu1,1)=Q(jdu1)*qn/s28+L(jdu1)*ln/propZ28
      gamn28(jdu1,2)=Q(jdu1)*qn/s28+R(jdu1)*ln/propZ28
      do jdu2=1,2
      gam28(jdu1,jdu2,1,1)=Q(jdu1)*Q(jdu2)/s28+L(jdu1)*L(jdu2)/propZ28
      gam28(jdu1,jdu2,1,2)=Q(jdu1)*Q(jdu2)/s28+L(jdu1)*R(jdu2)/propZ28
      gam28(jdu1,jdu2,2,1)=Q(jdu1)*Q(jdu2)/s28+R(jdu1)*L(jdu2)/propZ28
      gam28(jdu1,jdu2,2,2)=Q(jdu1)*Q(jdu2)/s28+R(jdu1)*R(jdu2)/propZ28
      enddo
      enddo

      p1=i1
      p2=i2
      p3=i3
      p4=i4
      p5=i5
      p6=i6
      p7=i7
      p8=i8

      do h28=1,2
      if (h28 == 1) then
      p2=i2
      p8=i8
      elseif (h28 == 2) then
      p2=i8
      p8=i2
      endif
      do jdu2=1,2
      amp(1,jdu2,1,h28)= + gam28(1,jdu2,1,h28)*propw56**(-1)*
     & propw34**(-1)*cxw**(-1)*s128**(-1)*s567**(-1) * ( 2.D0*za(p7,p5)
     &    *zb(p1,p2)*zba2(p6,p5,p7,p3)*zba2(p4,p1,p2,p8) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gam28(1,jdu2,1,h28)*
     & propw56**(-1)*propw34**(-1)*cxw**(-1)*s278**(-1) * ( 2.D0*za(p7,
     &    p8)*zb(p1,p4)*zba2(p6,p1,p4,p3)*zba2(p2,p7,p8,p5)*s134**(-1)
     &     )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gam28(1,jdu2,1,h28)*
     & propw56**(-1)*s128**(-1) * ( 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,p2)*
     &    zb(p5,p6)*zba2(p4,p1,p2,p8)*gamn3456(1,1,1)*s356**(-1) + 4.D0
     &    *za(p7,p3)*za(p5,p6)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)*
     &    game3456(1,1,1)*s456**(-1) - 4.D0*za(p7,p3)*za(p5,p3)*zb(p1,
     &    p2)*zb(p6,p3)*zba2(p4,p1,p2,p8)*gamn3456(1,1,1)*s356**(-1) + 
     &    4.D0*za(p7,p3)*za(p5,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8
     &    )*game3456(1,1,1)*s456**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gam28(1,jdu2,1,h28)*
     & propw56**(-1)*s278**(-1) * (  - 4.D0*za(p7,p8)*za(p5,p6)*zb(p1,
     &    p6)*zb(p6,p4)*zba2(p2,p7,p8,p3)*game3456(1,1,1)*s456**(-1) - 
     &    4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p5,p6)*zba2(p2,p7,p8,p5
     &    )*gamn3456(1,1,1)*s356**(-1) + 4.D0*za(p7,p8)*za(p5,p3)*zb(p1
     &    ,p4)*zb(p6,p3)*zba2(p2,p7,p8,p3)*gamn3456(1,1,1)*s356**(-1)
     &     - 4.D0*za(p7,p8)*za(p5,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8
     &    ,p3)*game3456(1,1,1)*s456**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gam28(1,jdu2,1,h28)*
     & propw34**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,
     &    p2)*zb(p5,p4)*zba2(p6,p1,p2,p8)*game3456(1,1,1)*s345**(-1) + 
     &    4.D0*za(p7,p5)*za(p6,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8
     &    )*gamn3456(1,1,1)*s346**(-1) - 4.D0*za(p7,p5)*za(p3,p4)*zb(p1
     &    ,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)*gamn3456(1,1,1)*s346**(-1)
     &     - 4.D0*za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2
     &    ,p8)*game3456(1,1,1)*s345**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gam28(1,jdu2,1,h28)*
     & propw34**(-1)*s278**(-1) * ( 4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p6)*
     &    zb(p5,p4)*zba2(p2,p7,p8,p5)*game3456(1,1,1)*s345**(-1) + 4.D0
     &    *za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p3,p4)*zba2(p2,p7,p8,p3)*
     &    game3456(1,1,1)*s345**(-1) - 4.D0*za(p7,p8)*za(p6,p3)*zb(p1,
     &    p6)*zb(p6,p4)*zba2(p2,p7,p8,p5)*gamn3456(1,1,1)*s346**(-1) + 
     &    4.D0*za(p7,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5
     &    )*gamn3456(1,1,1)*s346**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw56**(-1)*propw34**(-1)*cxw**(-1)*s567**(-1) * ( 2.D0*za(p7,
     &    p5)*zb(p1,p4)*zba2(p6,p5,p7,p8)*zba2(p2,p1,p4,p3)*s134**(-1)
     &     )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv(1,1)*gam28(1,jdu2,1,
     & h28)*propw56**(-1)*propw34**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5
     &    )*zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8) + 4.D0*za(p7,
     &    p3)*zb(p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8) - 2.D0*za(
     &    p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0
     &    *za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 
     &    2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2
     &    ) - 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,
     &    p4,p2) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv(1,1)*gam28(1,jdu2,1,
     & h28)*propw56**(-1)*propw34**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*
     &    za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1) - 2.
     &    D0*za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1)
     &     + 4.D0*za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,
     &    p5) - 4.D0*za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,
     &    p8,p3) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2348**(-1)*cwmass2**(-1)*
     & cxw**(-1) * (  - za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*s28 + 
     &    za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*s34 )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2348**(-1)*cxw**(-1)*s567**(-1)
     &  * ( za(p7,p5)*za(p3,p8)*zb(p7,p6)*zb(p4,p2)*zab2(p7,p3,p4,p1)
     &     - za(p7,p5)*za(p3,p8)*zb(p7,p6)*zb(p4,p2)*zab2(p7,p2,p8,p1)
     &     + za(p7,p5)*za(p3,p8)*zb(p5,p6)*zb(p4,p2)*zab2(p5,p3,p4,p1)
     &     - za(p7,p5)*za(p3,p8)*zb(p5,p6)*zb(p4,p2)*zab2(p5,p2,p8,p1)
     &     + 2.D0*za(p7,p5)*zb(p1,p4)*zab2(p8,p3,p4,p2)*zba2(p6,p5,p7,
     &    p3) - 2.D0*za(p7,p5)*zb(p1,p2)*zab2(p3,p2,p8,p4)*zba2(p6,p5,
     &    p7,p8) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1)*
     & cxw**(-1) * (  - za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*s28 + 
     &    za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*s56 )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2568**(-1)*cxw**(-1) * ( 2.D0*
     &    za(p7,p5)*zb(p1,p4)*zab2(p8,p5,p6,p2)*zba2(p6,p1,p4,p3)*
     &    s134**(-1) - 2.D0*za(p7,p8)*zb(p1,p4)*zab2(p5,p2,p8,p6)*zba2(
     &    p2,p1,p4,p3)*s134**(-1) + za(p1,p3)*za(p5,p8)*zb(p1,p4)*zb(p6
     &    ,p2)*zab2(p7,p5,p6,p1)*s134**(-1) - za(p1,p3)*za(p5,p8)*zb(p1
     &    ,p4)*zb(p6,p2)*zab2(p7,p2,p8,p1)*s134**(-1) - za(p5,p8)*za(p3
     &    ,p4)*zb(p1,p4)*zb(p6,p2)*zab2(p7,p5,p6,p4)*s134**(-1) + za(p5
     &    ,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p2)*zab2(p7,p2,p8,p4)*
     &    s134**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cwmass2**(-1)*cxw**(-1) * ( za(p7,
     &    p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cxw**(-1)*s567**(-1)*s238**(-1)
     &  * ( 2.D0*za(p7,p5)*za(p3,p8)*zb(p1,p4)*zb(p3,p2)*zba2(p6,p5,p7,
     &    p3) + 2.D0*za(p7,p5)*za(p3,p8)*zb(p1,p4)*zb(p8,p2)*zba2(p6,p5
     &    ,p7,p8) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cwmass2**(-1)*cxw**(-1) * ( za(p7,
     &    p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cxw**(-1)*s268**(-1) * ( 2.D0*za(
     &    p7,p5)*za(p6,p8)*zb(p1,p4)*zb(p6,p2)*zba2(p6,p1,p4,p3)*
     &    s134**(-1) - 2.D0*za(p7,p5)*za(p8,p2)*zb(p1,p4)*zb(p6,p2)*
     &    zba2(p2,p1,p4,p3)*s134**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*gamv17(1
     & ,1)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * ( 2.D0*za(p7,
     &    p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p56Dp28 + 2.D0*za(p7,p5)*
     &    za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p34Dp56 - 2.D0*za(p7,p5)*za(p3,
     &    p8)*zb(p1,p6)*zb(p4,p2)*p17Dp28 - 2.D0*za(p7,p5)*za(p3,p8)*
     &    zb(p1,p6)*zb(p4,p2)*p17Dp34 )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*gamv17(1
     & ,1)*propw56**(-1)*propW2348**(-1)*s238**(-1) * ( 2.D0*za(p7,p5)*
     &    za(p3,p8)*zb(p1,p6)*zb(p3,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p7,
     &    p5)*za(p3,p8)*zb(p1,p6)*zb(p3,p2)*zab2(p3,p5,p6,p4) + 2.D0*
     &    za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p8,p2)*zab2(p8,p1,p7,p4) - 2.
     &    D0*za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p8,p2)*zab2(p8,p5,p6,p4)
     &     - 4.D0*za(p7,p3)*za(p3,p8)*zb(p1,p4)*zb(p3,p2)*zab2(p5,p1,p7
     &    ,p6) - 4.D0*za(p7,p8)*za(p3,p8)*zb(p1,p4)*zb(p8,p2)*zab2(p5,
     &    p1,p7,p6) + 4.D0*za(p5,p3)*za(p3,p8)*zb(p6,p4)*zb(p3,p2)*
     &    zab2(p7,p5,p6,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zab2(p7,p5,p6,p1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*gamv17(1
     & ,1)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * ( 2.D0*za(p7,
     &    p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp56 + 2.D0*za(p7,p3)*
     &    za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp28 - 2.D0*za(p7,p3)*za(p5,
     &    p8)*zb(p1,p4)*zb(p6,p2)*p17Dp28 - 2.D0*za(p7,p3)*za(p5,p8)*
     &    zb(p1,p4)*zb(p6,p2)*p17Dp56 )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + game28(jdu2,h28)*gamv17(1
     & ,1)*propw34**(-1)*propW2568**(-1)*s268**(-1) * (  - 4.D0*za(p7,
     &    p5)*za(p6,p8)*zb(p1,p6)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 4.D0*
     &    za(p7,p5)*za(p8,p2)*zb(p1,p2)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 2.
     &    D0*za(p7,p3)*za(p6,p8)*zb(p1,p4)*zb(p6,p2)*zab2(p5,p1,p7,p6)
     &     - 2.D0*za(p7,p3)*za(p6,p8)*zb(p1,p4)*zb(p6,p2)*zab2(p5,p3,p4
     &    ,p6) - 2.D0*za(p7,p3)*za(p8,p2)*zb(p1,p4)*zb(p6,p2)*zab2(p5,
     &    p1,p7,p2) + 2.D0*za(p7,p3)*za(p8,p2)*zb(p1,p4)*zb(p6,p2)*
     &    zab2(p5,p3,p4,p2) + 4.D0*za(p5,p3)*za(p6,p8)*zb(p6,p4)*zb(p6,
     &    p2)*zab2(p7,p3,p4,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zab2(p7,p3,p4,p1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cwmass2**(-1)*cxw**(-1) * (  - za(
     &    p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cxw**(-1)*s567**(-1)*s248**(-1)
     &  * (  - 2.D0*za(p7,p5)*za(p4,p8)*zb(p1,p4)*zb(p4,p2)*zba2(p6,p5,
     &    p7,p3) + 2.D0*za(p7,p5)*za(p8,p2)*zb(p1,p2)*zb(p4,p2)*zba2(p6
     &    ,p5,p7,p3) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cwmass2**(-1)*cxw**(-1) * (  - za(
     &    p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cxw**(-1)*s258**(-1) * (  - 2.D0*
     &    za(p7,p5)*za(p5,p8)*zb(p1,p4)*zb(p5,p2)*zba2(p6,p1,p4,p3)*
     &    s134**(-1) - 2.D0*za(p7,p8)*za(p5,p8)*zb(p1,p4)*zb(p8,p2)*
     &    zba2(p6,p1,p4,p3)*s134**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,1)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p56Dp28 - 2.D0*za(p7,p5)
     &    *za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p34Dp56 + 2.D0*za(p7,p5)*za(p3
     &    ,p8)*zb(p1,p6)*zb(p4,p2)*p17Dp28 + 2.D0*za(p7,p5)*za(p3,p8)*
     &    zb(p1,p6)*zb(p4,p2)*p17Dp34 )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,1)*propw56**(-1)*propW2348**(-1)*s248**(-1) * (  - 2.D0*za(p7,
     &    p5)*za(p4,p8)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p1,p7,p4) + 2.D0*
     &    za(p7,p5)*za(p4,p8)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p5,p6,p4) + 2.
     &    D0*za(p7,p5)*za(p8,p2)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p1,p7,p2)
     &     - 2.D0*za(p7,p5)*za(p8,p2)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p5,p6
     &    ,p2) + 4.D0*za(p7,p3)*za(p4,p8)*zb(p1,p4)*zb(p4,p2)*zab2(p5,
     &    p1,p7,p6) - 4.D0*za(p7,p3)*za(p8,p2)*zb(p1,p2)*zb(p4,p2)*
     &    zab2(p5,p1,p7,p6) - 4.D0*za(p5,p3)*za(p4,p8)*zb(p6,p4)*zb(p4,
     &    p2)*zab2(p7,p5,p6,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zab2(p7,p5,p6,p1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,1)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp56 - 2.D0*za(p7,p3)
     &    *za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp28 + 2.D0*za(p7,p3)*za(p5
     &    ,p8)*zb(p1,p4)*zb(p6,p2)*p17Dp28 + 2.D0*za(p7,p3)*za(p5,p8)*
     &    zb(p1,p4)*zb(p6,p2)*p17Dp56 )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,1)*propw34**(-1)*propW2568**(-1)*s258**(-1) * ( 4.D0*za(p7,p5)*
     &    za(p5,p8)*zb(p1,p6)*zb(p5,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p7,
     &    p3)*za(p5,p8)*zb(p1,p4)*zb(p5,p2)*zab2(p5,p1,p7,p6) + 2.D0*
     &    za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p5,p2)*zab2(p5,p3,p4,p6) - 2.
     &    D0*za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p8,p2)*zab2(p8,p1,p7,p6)
     &     + 2.D0*za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p8,p2)*zab2(p8,p3,p4
     &    ,p6) + 4.D0*za(p7,p8)*za(p5,p8)*zb(p1,p6)*zb(p8,p2)*zab2(p3,
     &    p1,p7,p4) - 4.D0*za(p5,p3)*za(p5,p8)*zb(p5,p2)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zab2(p7,p3,p4,p1) )
      amp(1,jdu2,2,h28)= + gam28(1,jdu2,2,h28)*propw56**(-1)*s128**(-1)
     &  * ( 4.D0*za(p1,p8)*za(p5,p6)*zb(p7,p6)*zb(p6,p4)*zab2(p3,p1,p8,
     &    p2)*game3456(1,2,1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p3)*zb(
     &    p7,p4)*zb(p5,p6)*zab2(p5,p1,p8,p2)*gamn3456(1,2,1)*s356**(-1)
     &     - 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p4)*zb(p6,p3)*zab2(p3,p1,p8
     &    ,p2)*gamn3456(1,2,1)*s356**(-1) + 4.D0*za(p1,p8)*za(p5,p4)*
     &    zb(p7,p4)*zb(p6,p4)*zab2(p3,p1,p8,p2)*game3456(1,2,1)*
     &    s456**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gam28(1,jdu2,2,h28)*
     & propw56**(-1)*s278**(-1) * (  - 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,
     &    p2)*zb(p5,p6)*zab2(p8,p2,p7,p4)*gamn3456(1,2,1)*s356**(-1) - 
     &    4.D0*za(p1,p3)*za(p5,p6)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6
     &    )*game3456(1,2,1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p3)*zb(p7
     &    ,p2)*zb(p6,p3)*zab2(p8,p2,p7,p4)*gamn3456(1,2,1)*s356**(-1)
     &     - 4.D0*za(p1,p3)*za(p5,p4)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7
     &    ,p4)*game3456(1,2,1)*s456**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gam28(1,jdu2,2,h28)*
     & propw34**(-1)*s128**(-1) * (  - 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,
     &    p6)*zb(p5,p4)*zab2(p5,p1,p8,p2)*game3456(1,2,1)*s345**(-1) - 
     &    4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p6)*zb(p3,p4)*zab2(p3,p1,p8,p2
     &    )*game3456(1,2,1)*s345**(-1) + 4.D0*za(p1,p8)*za(p6,p3)*zb(p7
     &    ,p6)*zb(p6,p4)*zab2(p5,p1,p8,p2)*gamn3456(1,2,1)*s346**(-1)
     &     - 4.D0*za(p1,p8)*za(p3,p4)*zb(p7,p4)*zb(p6,p4)*zab2(p5,p1,p8
     &    ,p2)*gamn3456(1,2,1)*s346**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gam28(1,jdu2,2,h28)*
     & propw34**(-1)*s278**(-1) * ( 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,p2)*
     &    zb(p5,p4)*zab2(p8,p2,p7,p6)*game3456(1,2,1)*s345**(-1) - 4.D0
     &    *za(p1,p5)*za(p6,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6)*
     &    gamn3456(1,2,1)*s346**(-1) + 4.D0*za(p1,p5)*za(p3,p4)*zb(p7,
     &    p2)*zb(p6,p4)*zab2(p8,p2,p7,p4)*gamn3456(1,2,1)*s346**(-1) + 
     &    4.D0*za(p1,p3)*za(p5,p3)*zb(p7,p2)*zb(p3,p4)*zab2(p8,p2,p7,p6
     &    )*game3456(1,2,1)*s345**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamv(1,2)*gam28(1,jdu2,2,
     & h28)*propw56**(-1)*propw34**(-1)*s128**(-1) * (  - 2.D0*za(p1,p8
     &    )*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) + 2.D0*za(
     &    p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) - 2.D0
     &    *za(p1,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p5,p6,p8) + 
     &    2.D0*za(p1,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p3,p4,p8
     &    ) - 4.D0*za(p1,p8)*zb(p7,p6)*zab2(p5,p1,p8,p2)*zab2(p3,p5,p6,
     &    p4) + 4.D0*za(p1,p8)*zb(p7,p4)*zab2(p5,p3,p4,p6)*zab2(p3,p1,
     &    p8,p2) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamv(1,2)*gam28(1,jdu2,2,
     & h28)*propw56**(-1)*propw34**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) + 4.D0*
     &    za(p1,p5)*zb(p7,p2)*zab2(p3,p5,p6,p4)*zab2(p8,p2,p7,p6) - 4.D0
     &    *za(p1,p3)*zb(p7,p2)*zab2(p5,p3,p4,p6)*zab2(p8,p2,p7,p4) - 2.D
     &    0*za(p5,p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p5,p6,p1)
     &     + 2.D0*za(p5,p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p3,p4
     &    ,p1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + game28(jdu2,h28)*gamv17(1
     & ,2)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * ( 2.D0*za(p1,
     &    p5)*za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p56Dp28 + 2.D0*za(p1,p5)*
     &    za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p34Dp56 - 2.D0*za(p1,p5)*za(p3,
     &    p8)*zb(p7,p6)*zb(p4,p2)*p17Dp28 - 2.D0*za(p1,p5)*za(p3,p8)*
     &    zb(p7,p6)*zb(p4,p2)*p17Dp34 )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + game28(jdu2,h28)*gamv17(1
     & ,2)*propw56**(-1)*propW2348**(-1)*s238**(-1) * ( 2.D0*za(p1,p5)*
     &    za(p3,p8)*zb(p7,p6)*zb(p3,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p1,
     &    p5)*za(p3,p8)*zb(p7,p6)*zb(p3,p2)*zab2(p3,p5,p6,p4) + 2.D0*
     &    za(p1,p5)*za(p3,p8)*zb(p7,p6)*zb(p8,p2)*zab2(p8,p1,p7,p4) - 2.
     &    D0*za(p1,p5)*za(p3,p8)*zb(p7,p6)*zb(p8,p2)*zab2(p8,p5,p6,p4)
     &     - 4.D0*za(p1,p3)*za(p3,p8)*zb(p7,p4)*zb(p3,p2)*zab2(p5,p1,p7
     &    ,p6) - 4.D0*za(p1,p8)*za(p3,p8)*zb(p7,p4)*zb(p8,p2)*zab2(p5,
     &    p1,p7,p6) + 4.D0*za(p5,p3)*za(p3,p8)*zb(p6,p4)*zb(p3,p2)*
     &    zba2(p7,p5,p6,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zba2(p7,p5,p6,p1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + game28(jdu2,h28)*gamv17(1
     & ,2)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * ( 2.D0*za(p1,
     &    p3)*za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp56 + 2.D0*za(p1,p3)*
     &    za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp28 - 2.D0*za(p1,p3)*za(p5,
     &    p8)*zb(p7,p4)*zb(p6,p2)*p17Dp28 - 2.D0*za(p1,p3)*za(p5,p8)*
     &    zb(p7,p4)*zb(p6,p2)*p17Dp56 )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + game28(jdu2,h28)*gamv17(1
     & ,2)*propw34**(-1)*propW2568**(-1)*s268**(-1) * (  - 4.D0*za(p1,
     &    p5)*za(p6,p8)*zb(p7,p6)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 4.D0*
     &    za(p1,p5)*za(p8,p2)*zb(p7,p2)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 2.
     &    D0*za(p1,p3)*za(p6,p8)*zb(p7,p4)*zb(p6,p2)*zab2(p5,p1,p7,p6)
     &     - 2.D0*za(p1,p3)*za(p6,p8)*zb(p7,p4)*zb(p6,p2)*zab2(p5,p3,p4
     &    ,p6) - 2.D0*za(p1,p3)*za(p8,p2)*zb(p7,p4)*zb(p6,p2)*zab2(p5,
     &    p1,p7,p2) + 2.D0*za(p1,p3)*za(p8,p2)*zb(p7,p4)*zb(p6,p2)*
     &    zab2(p5,p3,p4,p2) + 4.D0*za(p5,p3)*za(p6,p8)*zb(p6,p4)*zb(p6,
     &    p2)*zba2(p7,p3,p4,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zba2(p7,p3,p4,p1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,2)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p1,p5)*za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p56Dp28 - 2.D0*za(p1,p5)
     &    *za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p34Dp56 + 2.D0*za(p1,p5)*za(p3
     &    ,p8)*zb(p7,p6)*zb(p4,p2)*p17Dp28 + 2.D0*za(p1,p5)*za(p3,p8)*
     &    zb(p7,p6)*zb(p4,p2)*p17Dp34 )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,2)*propw56**(-1)*propW2348**(-1)*s248**(-1) * (  - 2.D0*za(p1,
     &    p5)*za(p4,p8)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p1,p7,p4) + 2.D0*
     &    za(p1,p5)*za(p4,p8)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p5,p6,p4) + 2.
     &    D0*za(p1,p5)*za(p8,p2)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p1,p7,p2)
     &     - 2.D0*za(p1,p5)*za(p8,p2)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p5,p6
     &    ,p2) + 4.D0*za(p1,p3)*za(p4,p8)*zb(p7,p4)*zb(p4,p2)*zab2(p5,
     &    p1,p7,p6) - 4.D0*za(p1,p3)*za(p8,p2)*zb(p7,p2)*zb(p4,p2)*
     &    zab2(p5,p1,p7,p6) - 4.D0*za(p5,p3)*za(p4,p8)*zb(p6,p4)*zb(p4,
     &    p2)*zba2(p7,p5,p6,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zba2(p7,p5,p6,p1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,2)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp56 - 2.D0*za(p1,p3)
     &    *za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp28 + 2.D0*za(p1,p3)*za(p5
     &    ,p8)*zb(p7,p4)*zb(p6,p2)*p17Dp28 + 2.D0*za(p1,p3)*za(p5,p8)*
     &    zb(p7,p4)*zb(p6,p2)*p17Dp56 )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(1
     & ,2)*propw34**(-1)*propW2568**(-1)*s258**(-1) * ( 4.D0*za(p1,p5)*
     &    za(p5,p8)*zb(p7,p6)*zb(p5,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p1,
     &    p3)*za(p5,p8)*zb(p7,p4)*zb(p5,p2)*zab2(p5,p1,p7,p6) + 2.D0*
     &    za(p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p5,p2)*zab2(p5,p3,p4,p6) - 2.
     &    D0*za(p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p8,p2)*zab2(p8,p1,p7,p6)
     &     + 2.D0*za(p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p8,p2)*zab2(p8,p3,p4
     &    ,p6) + 4.D0*za(p1,p8)*za(p5,p8)*zb(p7,p6)*zb(p8,p2)*zab2(p3,
     &    p1,p7,p4) - 4.D0*za(p5,p3)*za(p5,p8)*zb(p5,p2)*zb(p6,p4)*
     &    zba2(p7,p3,p4,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zba2(p7,p3,p4,p1) )
      amp(2,jdu2,1,h28)= + gam28(1,jdu2,1,h28)*propw56**(-1)*
     & propw34**(-1)*cxw**(-1) * ( 2.D0*za(p7,p3)*za(p7,p8)*za(p1,p5)*
     &    zb(p7,p4)*zb(p1,p6)*zb(p1,p2)*s347**(-1)*s156**(-1) - 2.D0*
     &    za(p7,p3)*za(p7,p8)*za(p5,p6)*zb(p7,p4)*zb(p1,p6)*zb(p6,p2)*
     &    s347**(-1)*s156**(-1) + 2.D0*za(p7,p3)*za(p1,p5)*za(p3,p8)*
     &    zb(p1,p6)*zb(p1,p2)*zb(p3,p4)*s347**(-1)*s156**(-1) - 2.D0*
     &    za(p7,p3)*za(p5,p6)*za(p3,p8)*zb(p1,p6)*zb(p6,p2)*zb(p3,p4)*
     &    s347**(-1)*s156**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw56**(-1)*propw34**(-1)*cxw**(-1)*s128**(-1) * (  - 2.D0*za(
     &    p7,p5)*za(p7,p3)*zb(p7,p4)*zb(p1,p2)*zba2(p6,p1,p2,p8)*
     &    s347**(-1) + 2.D0*za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*
     &    zba2(p6,p1,p2,p8)*s347**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw56**(-1)*propw34**(-1)*cxw**(-1)*s278**(-1) * (  - 2.D0*za(
     &    p7,p8)*za(p1,p5)*zb(p1,p6)*zb(p1,p4)*zba2(p2,p7,p8,p3)*
     &    s156**(-1) + 2.D0*za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*
     &    zba2(p2,p7,p8,p3)*s156**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw56**(-1)*s128**(-1) * ( 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,p2)*
     &    zb(p5,p6)*zba2(p4,p1,p2,p8)*gamn3456(2,1,1)*s356**(-1) + 4.D0
     &    *za(p7,p3)*za(p5,p6)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)*
     &    game3456(2,1,1)*s456**(-1) - 4.D0*za(p7,p3)*za(p5,p3)*zb(p1,
     &    p2)*zb(p6,p3)*zba2(p4,p1,p2,p8)*gamn3456(2,1,1)*s356**(-1) + 
     &    4.D0*za(p7,p3)*za(p5,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8
     &    )*game3456(2,1,1)*s456**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw56**(-1)*s278**(-1) * (  - 4.D0*za(p7,p8)*za(p5,p6)*zb(p1,
     &    p6)*zb(p6,p4)*zba2(p2,p7,p8,p3)*game3456(2,1,1)*s456**(-1) - 
     &    4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p5,p6)*zba2(p2,p7,p8,p5
     &    )*gamn3456(2,1,1)*s356**(-1) + 4.D0*za(p7,p8)*za(p5,p3)*zb(p1
     &    ,p4)*zb(p6,p3)*zba2(p2,p7,p8,p3)*gamn3456(2,1,1)*s356**(-1)
     &     - 4.D0*za(p7,p8)*za(p5,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8
     &    ,p3)*game3456(2,1,1)*s456**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw34**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,
     &    p2)*zb(p5,p4)*zba2(p6,p1,p2,p8)*game3456(2,1,1)*s345**(-1) + 
     &    4.D0*za(p7,p5)*za(p6,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8
     &    )*gamn3456(2,1,1)*s346**(-1) - 4.D0*za(p7,p5)*za(p3,p4)*zb(p1
     &    ,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)*gamn3456(2,1,1)*s346**(-1)
     &     - 4.D0*za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2
     &    ,p8)*game3456(2,1,1)*s345**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gam28(2,jdu2,1,h28)*
     & propw34**(-1)*s278**(-1) * ( 4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p6)*
     &    zb(p5,p4)*zba2(p2,p7,p8,p5)*game3456(2,1,1)*s345**(-1) + 4.D0
     &    *za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p3,p4)*zba2(p2,p7,p8,p3)*
     &    game3456(2,1,1)*s345**(-1) - 4.D0*za(p7,p8)*za(p6,p3)*zb(p1,
     &    p6)*zb(p6,p4)*zba2(p2,p7,p8,p5)*gamn3456(2,1,1)*s346**(-1) + 
     &    4.D0*za(p7,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5
     &    )*gamn3456(2,1,1)*s346**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv(2,1)*gam28(2,jdu2,1,
     & h28)*propw56**(-1)*propw34**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5
     &    )*zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8) + 4.D0*za(p7,
     &    p3)*zb(p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8) - 2.D0*za(
     &    p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0
     &    *za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 
     &    2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2
     &    ) - 2.D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,
     &    p4,p2) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv(2,1)*gam28(2,jdu2,1,
     & h28)*propw56**(-1)*propw34**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*
     &    za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1) - 2.
     &    D0*za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1)
     &     + 4.D0*za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,
     &    p5) - 4.D0*za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,
     &    p8,p3) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2348**(-1)*cwmass2**(-1)*
     & cxw**(-1) * ( za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*s28 - za(
     &    p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*s34 )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2348**(-1)*cxw**(-1) * ( 2.D0*
     &    za(p7,p3)*za(p1,p5)*zb(p1,p6)*zb(p1,p4)*zab2(p8,p3,p4,p2)*
     &    s156**(-1) - 2.D0*za(p7,p3)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*
     &    zab2(p8,p3,p4,p2)*s156**(-1) - 2.D0*za(p7,p8)*za(p1,p5)*zb(p1
     &    ,p6)*zb(p1,p2)*zab2(p3,p2,p8,p4)*s156**(-1) + 2.D0*za(p7,p8)*
     &    za(p5,p6)*zb(p1,p6)*zb(p6,p2)*zab2(p3,p2,p8,p4)*s156**(-1) - 
     &    za(p1,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*zab2(p7,p3,p4,p1)*
     &    s156**(-1) + za(p1,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*zab2(p7,
     &    p2,p8,p1)*s156**(-1) + za(p5,p6)*za(p3,p8)*zb(p1,p6)*zb(p4,p2
     &    )*zab2(p7,p3,p4,p6)*s156**(-1) - za(p5,p6)*za(p3,p8)*zb(p1,p6
     &    )*zb(p4,p2)*zab2(p7,p2,p8,p6)*s156**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1)*
     & cxw**(-1) * ( za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*s28 - za(
     &    p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*s56 )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv28(jdu2,h28)*
     & propw56**(-1)*propw34**(-1)*propW2568**(-1)*cxw**(-1) * ( 2.D0*
     &    za(p7,p5)*za(p7,p3)*zb(p7,p4)*zb(p1,p6)*zab2(p8,p5,p6,p2)*
     &    s347**(-1) - 2.D0*za(p7,p3)*za(p7,p8)*zb(p7,p4)*zb(p1,p2)*
     &    zab2(p5,p2,p8,p6)*s347**(-1) - 2.D0*za(p7,p3)*za(p5,p3)*zb(p1
     &    ,p6)*zb(p3,p4)*zab2(p8,p5,p6,p2)*s347**(-1) - za(p7,p3)*za(p5
     &    ,p8)*zb(p7,p4)*zb(p6,p2)*zab2(p7,p5,p6,p1)*s347**(-1) + za(p7
     &    ,p3)*za(p5,p8)*zb(p7,p4)*zb(p6,p2)*zab2(p7,p2,p8,p1)*
     &    s347**(-1) - za(p7,p3)*za(p5,p8)*zb(p6,p2)*zb(p3,p4)*zab2(p3,
     &    p5,p6,p1)*s347**(-1) + za(p7,p3)*za(p5,p8)*zb(p6,p2)*zb(p3,p4
     &    )*zab2(p3,p2,p8,p1)*s347**(-1) - 2.D0*za(p7,p3)*za(p3,p8)*zb(
     &    p1,p2)*zb(p3,p4)*zab2(p5,p2,p8,p6)*s347**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cwmass2**(-1)*cxw**(-1) * (  - za(
     &    p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cxw**(-1)*s238**(-1) * ( 2.D0*za(
     &    p7,p3)*za(p1,p5)*za(p3,p8)*zb(p1,p6)*zb(p1,p4)*zb(p3,p2)*
     &    s156**(-1) - 2.D0*za(p7,p3)*za(p5,p6)*za(p3,p8)*zb(p1,p6)*zb(
     &    p6,p4)*zb(p3,p2)*s156**(-1) + 2.D0*za(p7,p8)*za(p1,p5)*za(p3,
     &    p8)*zb(p1,p6)*zb(p1,p4)*zb(p8,p2)*s156**(-1) - 2.D0*za(p7,p8)
     &    *za(p5,p6)*za(p3,p8)*zb(p1,p6)*zb(p6,p4)*zb(p8,p2)*s156**(-1)
     &     )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cwmass2**(-1)*cxw**(-1) * (  - za(
     &    p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cxw**(-1)*s268**(-1) * ( 2.D0*za(
     &    p7,p5)*za(p7,p3)*za(p6,p8)*zb(p7,p4)*zb(p1,p6)*zb(p6,p2)*
     &    s347**(-1) - 2.D0*za(p7,p5)*za(p7,p3)*za(p8,p2)*zb(p7,p4)*zb(
     &    p1,p2)*zb(p6,p2)*s347**(-1) - 2.D0*za(p7,p3)*za(p5,p3)*za(p6,
     &    p8)*zb(p1,p6)*zb(p6,p2)*zb(p3,p4)*s347**(-1) + 2.D0*za(p7,p3)
     &    *za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p2)*zb(p3,p4)*s347**(-1)
     &     )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*gamv17(2
     & ,1)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * ( 2.D0*za(p7,
     &    p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p56Dp28 + 2.D0*za(p7,p5)*
     &    za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p34Dp56 - 2.D0*za(p7,p5)*za(p3,
     &    p8)*zb(p1,p6)*zb(p4,p2)*p17Dp28 - 2.D0*za(p7,p5)*za(p3,p8)*
     &    zb(p1,p6)*zb(p4,p2)*p17Dp34 )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*gamv17(2
     & ,1)*propw56**(-1)*propW2348**(-1)*s238**(-1) * ( 2.D0*za(p7,p5)*
     &    za(p3,p8)*zb(p1,p6)*zb(p3,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p7,
     &    p5)*za(p3,p8)*zb(p1,p6)*zb(p3,p2)*zab2(p3,p5,p6,p4) + 2.D0*
     &    za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p8,p2)*zab2(p8,p1,p7,p4) - 2.
     &    D0*za(p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p8,p2)*zab2(p8,p5,p6,p4)
     &     - 4.D0*za(p7,p3)*za(p3,p8)*zb(p1,p4)*zb(p3,p2)*zab2(p5,p1,p7
     &    ,p6) - 4.D0*za(p7,p8)*za(p3,p8)*zb(p1,p4)*zb(p8,p2)*zab2(p5,
     &    p1,p7,p6) + 4.D0*za(p5,p3)*za(p3,p8)*zb(p6,p4)*zb(p3,p2)*
     &    zab2(p7,p5,p6,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zab2(p7,p5,p6,p1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*gamv17(2
     & ,1)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * ( 2.D0*za(p7,
     &    p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp56 + 2.D0*za(p7,p3)*
     &    za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp28 - 2.D0*za(p7,p3)*za(p5,
     &    p8)*zb(p1,p4)*zb(p6,p2)*p17Dp28 - 2.D0*za(p7,p3)*za(p5,p8)*
     &    zb(p1,p4)*zb(p6,p2)*p17Dp56 )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + game28(jdu2,h28)*gamv17(2
     & ,1)*propw34**(-1)*propW2568**(-1)*s268**(-1) * (  - 4.D0*za(p7,
     &    p5)*za(p6,p8)*zb(p1,p6)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 4.D0*
     &    za(p7,p5)*za(p8,p2)*zb(p1,p2)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 2.
     &    D0*za(p7,p3)*za(p6,p8)*zb(p1,p4)*zb(p6,p2)*zab2(p5,p1,p7,p6)
     &     - 2.D0*za(p7,p3)*za(p6,p8)*zb(p1,p4)*zb(p6,p2)*zab2(p5,p3,p4
     &    ,p6) - 2.D0*za(p7,p3)*za(p8,p2)*zb(p1,p4)*zb(p6,p2)*zab2(p5,
     &    p1,p7,p2) + 2.D0*za(p7,p3)*za(p8,p2)*zb(p1,p4)*zb(p6,p2)*
     &    zab2(p5,p3,p4,p2) + 4.D0*za(p5,p3)*za(p6,p8)*zb(p6,p4)*zb(p6,
     &    p2)*zab2(p7,p3,p4,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zab2(p7,p3,p4,p1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cwmass2**(-1)*cxw**(-1) * ( za(p7,
     &    p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw56**(-1)*propW2348**(-1)*cxw**(-1)*s248**(-1) * (  - 2.D0*
     &    za(p7,p3)*za(p1,p5)*za(p4,p8)*zb(p1,p6)*zb(p1,p4)*zb(p4,p2)*
     &    s156**(-1) + 2.D0*za(p7,p3)*za(p1,p5)*za(p8,p2)*zb(p1,p6)*zb(
     &    p1,p2)*zb(p4,p2)*s156**(-1) + 2.D0*za(p7,p3)*za(p5,p6)*za(p4,
     &    p8)*zb(p1,p6)*zb(p6,p4)*zb(p4,p2)*s156**(-1) - 2.D0*za(p7,p3)
     &    *za(p5,p6)*za(p8,p2)*zb(p1,p6)*zb(p6,p2)*zb(p4,p2)*s156**(-1)
     &     )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cwmass2**(-1)*cxw**(-1) * ( za(p7,
     &    p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*
     & propw34**(-1)*propW2568**(-1)*cxw**(-1)*s258**(-1) * (  - 2.D0*
     &    za(p7,p5)*za(p7,p3)*za(p5,p8)*zb(p7,p4)*zb(p1,p6)*zb(p5,p2)*
     &    s347**(-1) - 2.D0*za(p7,p3)*za(p7,p8)*za(p5,p8)*zb(p7,p4)*zb(
     &    p1,p6)*zb(p8,p2)*s347**(-1) + 2.D0*za(p7,p3)*za(p5,p3)*za(p5,
     &    p8)*zb(p1,p6)*zb(p5,p2)*zb(p3,p4)*s347**(-1) - 2.D0*za(p7,p3)
     &    *za(p5,p8)*za(p3,p8)*zb(p1,p6)*zb(p3,p4)*zb(p8,p2)*s347**(-1)
     &     )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,1)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p7,p5)*za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p56Dp28 - 2.D0*za(p7,p5)
     &    *za(p3,p8)*zb(p1,p6)*zb(p4,p2)*p34Dp56 + 2.D0*za(p7,p5)*za(p3
     &    ,p8)*zb(p1,p6)*zb(p4,p2)*p17Dp28 + 2.D0*za(p7,p5)*za(p3,p8)*
     &    zb(p1,p6)*zb(p4,p2)*p17Dp34 )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,1)*propw56**(-1)*propW2348**(-1)*s248**(-1) * (  - 2.D0*za(p7,
     &    p5)*za(p4,p8)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p1,p7,p4) + 2.D0*
     &    za(p7,p5)*za(p4,p8)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p5,p6,p4) + 2.
     &    D0*za(p7,p5)*za(p8,p2)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p1,p7,p2)
     &     - 2.D0*za(p7,p5)*za(p8,p2)*zb(p1,p6)*zb(p4,p2)*zab2(p3,p5,p6
     &    ,p2) + 4.D0*za(p7,p3)*za(p4,p8)*zb(p1,p4)*zb(p4,p2)*zab2(p5,
     &    p1,p7,p6) - 4.D0*za(p7,p3)*za(p8,p2)*zb(p1,p2)*zb(p4,p2)*
     &    zab2(p5,p1,p7,p6) - 4.D0*za(p5,p3)*za(p4,p8)*zb(p6,p4)*zb(p4,
     &    p2)*zab2(p7,p5,p6,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zab2(p7,p5,p6,p1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,1)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp56 - 2.D0*za(p7,p3)
     &    *za(p5,p8)*zb(p1,p4)*zb(p6,p2)*p34Dp28 + 2.D0*za(p7,p3)*za(p5
     &    ,p8)*zb(p1,p4)*zb(p6,p2)*p17Dp28 + 2.D0*za(p7,p3)*za(p5,p8)*
     &    zb(p1,p4)*zb(p6,p2)*p17Dp56 )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,1)*propw34**(-1)*propW2568**(-1)*s258**(-1) * ( 4.D0*za(p7,p5)*
     &    za(p5,p8)*zb(p1,p6)*zb(p5,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p7,
     &    p3)*za(p5,p8)*zb(p1,p4)*zb(p5,p2)*zab2(p5,p1,p7,p6) + 2.D0*
     &    za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p5,p2)*zab2(p5,p3,p4,p6) - 2.
     &    D0*za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p8,p2)*zab2(p8,p1,p7,p6)
     &     + 2.D0*za(p7,p3)*za(p5,p8)*zb(p1,p4)*zb(p8,p2)*zab2(p8,p3,p4
     &    ,p6) + 4.D0*za(p7,p8)*za(p5,p8)*zb(p1,p6)*zb(p8,p2)*zab2(p3,
     &    p1,p7,p4) - 4.D0*za(p5,p3)*za(p5,p8)*zb(p5,p2)*zb(p6,p4)*
     &    zab2(p7,p3,p4,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zab2(p7,p3,p4,p1) )
      amp(2,jdu2,2,h28)= + gam28(2,jdu2,2,h28)*propw56**(-1)*s128**(-1)
     &  * ( 4.D0*za(p1,p8)*za(p5,p6)*zb(p7,p6)*zb(p6,p4)*zab2(p3,p1,p8,
     &    p2)*game3456(2,2,1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p3)*zb(
     &    p7,p4)*zb(p5,p6)*zab2(p5,p1,p8,p2)*gamn3456(2,2,1)*s356**(-1)
     &     - 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p4)*zb(p6,p3)*zab2(p3,p1,p8
     &    ,p2)*gamn3456(2,2,1)*s356**(-1) + 4.D0*za(p1,p8)*za(p5,p4)*
     &    zb(p7,p4)*zb(p6,p4)*zab2(p3,p1,p8,p2)*game3456(2,2,1)*
     &    s456**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gam28(2,jdu2,2,h28)*
     & propw56**(-1)*s278**(-1) * (  - 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,
     &    p2)*zb(p5,p6)*zab2(p8,p2,p7,p4)*gamn3456(2,2,1)*s356**(-1) - 
     &    4.D0*za(p1,p3)*za(p5,p6)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6
     &    )*game3456(2,2,1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p3)*zb(p7
     &    ,p2)*zb(p6,p3)*zab2(p8,p2,p7,p4)*gamn3456(2,2,1)*s356**(-1)
     &     - 4.D0*za(p1,p3)*za(p5,p4)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7
     &    ,p4)*game3456(2,2,1)*s456**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gam28(2,jdu2,2,h28)*
     & propw34**(-1)*s128**(-1) * (  - 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,
     &    p6)*zb(p5,p4)*zab2(p5,p1,p8,p2)*game3456(2,2,1)*s345**(-1) - 
     &    4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p6)*zb(p3,p4)*zab2(p3,p1,p8,p2
     &    )*game3456(2,2,1)*s345**(-1) + 4.D0*za(p1,p8)*za(p6,p3)*zb(p7
     &    ,p6)*zb(p6,p4)*zab2(p5,p1,p8,p2)*gamn3456(2,2,1)*s346**(-1)
     &     - 4.D0*za(p1,p8)*za(p3,p4)*zb(p7,p4)*zb(p6,p4)*zab2(p5,p1,p8
     &    ,p2)*gamn3456(2,2,1)*s346**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gam28(2,jdu2,2,h28)*
     & propw34**(-1)*s278**(-1) * ( 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,p2)*
     &    zb(p5,p4)*zab2(p8,p2,p7,p6)*game3456(2,2,1)*s345**(-1) - 4.D0
     &    *za(p1,p5)*za(p6,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6)*
     &    gamn3456(2,2,1)*s346**(-1) + 4.D0*za(p1,p5)*za(p3,p4)*zb(p7,
     &    p2)*zb(p6,p4)*zab2(p8,p2,p7,p4)*gamn3456(2,2,1)*s346**(-1) + 
     &    4.D0*za(p1,p3)*za(p5,p3)*zb(p7,p2)*zb(p3,p4)*zab2(p8,p2,p7,p6
     &    )*game3456(2,2,1)*s345**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamv(2,2)*gam28(2,jdu2,2,
     & h28)*propw56**(-1)*propw34**(-1)*s128**(-1) * (  - 2.D0*za(p1,p8
     &    )*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) + 2.D0*za(
     &    p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) - 2.D0
     &    *za(p1,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p5,p6,p8) + 
     &    2.D0*za(p1,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p3,p4,p8
     &    ) - 4.D0*za(p1,p8)*zb(p7,p6)*zab2(p5,p1,p8,p2)*zab2(p3,p5,p6,
     &    p4) + 4.D0*za(p1,p8)*zb(p7,p4)*zab2(p5,p3,p4,p6)*zab2(p3,p1,
     &    p8,p2) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamv(2,2)*gam28(2,jdu2,2,
     & h28)*propw56**(-1)*propw34**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) + 4.D0*
     &    za(p1,p5)*zb(p7,p2)*zab2(p3,p5,p6,p4)*zab2(p8,p2,p7,p6) - 4.D0
     &    *za(p1,p3)*zb(p7,p2)*zab2(p5,p3,p4,p6)*zab2(p8,p2,p7,p4) - 2.D
     &    0*za(p5,p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p5,p6,p1)
     &     + 2.D0*za(p5,p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p3,p4
     &    ,p1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + game28(jdu2,h28)*gamv17(2
     & ,2)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * ( 2.D0*za(p1,
     &    p5)*za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p56Dp28 + 2.D0*za(p1,p5)*
     &    za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p34Dp56 - 2.D0*za(p1,p5)*za(p3,
     &    p8)*zb(p7,p6)*zb(p4,p2)*p17Dp28 - 2.D0*za(p1,p5)*za(p3,p8)*
     &    zb(p7,p6)*zb(p4,p2)*p17Dp34 )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + game28(jdu2,h28)*gamv17(2
     & ,2)*propw56**(-1)*propW2348**(-1)*s238**(-1) * ( 2.D0*za(p1,p5)*
     &    za(p3,p8)*zb(p7,p6)*zb(p3,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p1,
     &    p5)*za(p3,p8)*zb(p7,p6)*zb(p3,p2)*zab2(p3,p5,p6,p4) + 2.D0*
     &    za(p1,p5)*za(p3,p8)*zb(p7,p6)*zb(p8,p2)*zab2(p8,p1,p7,p4) - 2.
     &    D0*za(p1,p5)*za(p3,p8)*zb(p7,p6)*zb(p8,p2)*zab2(p8,p5,p6,p4)
     &     - 4.D0*za(p1,p3)*za(p3,p8)*zb(p7,p4)*zb(p3,p2)*zab2(p5,p1,p7
     &    ,p6) - 4.D0*za(p1,p8)*za(p3,p8)*zb(p7,p4)*zb(p8,p2)*zab2(p5,
     &    p1,p7,p6) + 4.D0*za(p5,p3)*za(p3,p8)*zb(p6,p4)*zb(p3,p2)*
     &    zba2(p7,p5,p6,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zba2(p7,p5,p6,p1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + game28(jdu2,h28)*gamv17(2
     & ,2)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * ( 2.D0*za(p1,
     &    p3)*za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp56 + 2.D0*za(p1,p3)*
     &    za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp28 - 2.D0*za(p1,p3)*za(p5,
     &    p8)*zb(p7,p4)*zb(p6,p2)*p17Dp28 - 2.D0*za(p1,p3)*za(p5,p8)*
     &    zb(p7,p4)*zb(p6,p2)*p17Dp56 )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + game28(jdu2,h28)*gamv17(2
     & ,2)*propw34**(-1)*propW2568**(-1)*s268**(-1) * (  - 4.D0*za(p1,
     &    p5)*za(p6,p8)*zb(p7,p6)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 4.D0*
     &    za(p1,p5)*za(p8,p2)*zb(p7,p2)*zb(p6,p2)*zab2(p3,p1,p7,p4) + 2.
     &    D0*za(p1,p3)*za(p6,p8)*zb(p7,p4)*zb(p6,p2)*zab2(p5,p1,p7,p6)
     &     - 2.D0*za(p1,p3)*za(p6,p8)*zb(p7,p4)*zb(p6,p2)*zab2(p5,p3,p4
     &    ,p6) - 2.D0*za(p1,p3)*za(p8,p2)*zb(p7,p4)*zb(p6,p2)*zab2(p5,
     &    p1,p7,p2) + 2.D0*za(p1,p3)*za(p8,p2)*zb(p7,p4)*zb(p6,p2)*
     &    zab2(p5,p3,p4,p2) + 4.D0*za(p5,p3)*za(p6,p8)*zb(p6,p4)*zb(p6,
     &    p2)*zba2(p7,p3,p4,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zba2(p7,p3,p4,p1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,2)*propw56**(-1)*propW2348**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p1,p5)*za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p56Dp28 - 2.D0*za(p1,p5)
     &    *za(p3,p8)*zb(p7,p6)*zb(p4,p2)*p34Dp56 + 2.D0*za(p1,p5)*za(p3
     &    ,p8)*zb(p7,p6)*zb(p4,p2)*p17Dp28 + 2.D0*za(p1,p5)*za(p3,p8)*
     &    zb(p7,p6)*zb(p4,p2)*p17Dp34 )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,2)*propw56**(-1)*propW2348**(-1)*s248**(-1) * (  - 2.D0*za(p1,
     &    p5)*za(p4,p8)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p1,p7,p4) + 2.D0*
     &    za(p1,p5)*za(p4,p8)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p5,p6,p4) + 2.
     &    D0*za(p1,p5)*za(p8,p2)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p1,p7,p2)
     &     - 2.D0*za(p1,p5)*za(p8,p2)*zb(p7,p6)*zb(p4,p2)*zab2(p3,p5,p6
     &    ,p2) + 4.D0*za(p1,p3)*za(p4,p8)*zb(p7,p4)*zb(p4,p2)*zab2(p5,
     &    p1,p7,p6) - 4.D0*za(p1,p3)*za(p8,p2)*zb(p7,p2)*zb(p4,p2)*
     &    zab2(p5,p1,p7,p6) - 4.D0*za(p5,p3)*za(p4,p8)*zb(p6,p4)*zb(p4,
     &    p2)*zba2(p7,p5,p6,p1) + 4.D0*za(p5,p3)*za(p8,p2)*zb(p6,p2)*
     &    zb(p4,p2)*zba2(p7,p5,p6,p1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,2)*propw34**(-1)*propW2568**(-1)*cwmass2**(-1) * (  - 2.D0*za(
     &    p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp56 - 2.D0*za(p1,p3)
     &    *za(p5,p8)*zb(p7,p4)*zb(p6,p2)*p34Dp28 + 2.D0*za(p1,p3)*za(p5
     &    ,p8)*zb(p7,p4)*zb(p6,p2)*p17Dp28 + 2.D0*za(p1,p3)*za(p5,p8)*
     &    zb(p7,p4)*zb(p6,p2)*p17Dp56 )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamn28(jdu2,h28)*gamv17(2
     & ,2)*propw34**(-1)*propW2568**(-1)*s258**(-1) * ( 4.D0*za(p1,p5)*
     &    za(p5,p8)*zb(p7,p6)*zb(p5,p2)*zab2(p3,p1,p7,p4) - 2.D0*za(p1,
     &    p3)*za(p5,p8)*zb(p7,p4)*zb(p5,p2)*zab2(p5,p1,p7,p6) + 2.D0*
     &    za(p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p5,p2)*zab2(p5,p3,p4,p6) - 2.
     &    D0*za(p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p8,p2)*zab2(p8,p1,p7,p6)
     &     + 2.D0*za(p1,p3)*za(p5,p8)*zb(p7,p4)*zb(p8,p2)*zab2(p8,p3,p4
     &    ,p6) + 4.D0*za(p1,p8)*za(p5,p8)*zb(p7,p6)*zb(p8,p2)*zab2(p3,
     &    p1,p7,p4) - 4.D0*za(p5,p3)*za(p5,p8)*zb(p5,p2)*zb(p6,p4)*
     &    zba2(p7,p3,p4,p1) + 4.D0*za(p5,p8)*za(p3,p8)*zb(p6,p4)*zb(p8,
     &    p2)*zba2(p7,p3,p4,p1) )
      enddo
      enddo
      amp(:,:,:,:)=amp(:,:,:,:)/cxw
      return
      end
