      subroutine amp2currentstrong(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer:: jdu1,jdu2,h28,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,amp(2,2,2,2),rxw,propZ3456,
     & propw34,propw56,gamn3456(2,2,2),game3456(2,2,2),gamv(2,2)
      real(dp):: qn,t3,t4,s134,s128,s345,s347,s356,s156,s456,
     & s278,s346,s567,s3456,s28,s34,s56
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
      s3456=t4(i3,i4,i5,i6)

      propw34=s34-cwmass2
      propw56=s56-cwmass2
      propZ3456=s3456-czmass2

      do jdu1=1,2
      game3456(jdu1,1,1)=Q(jdu1)*qe/s3456+L(jdu1)*le/propZ3456
      game3456(jdu1,1,2)=Q(jdu1)*qe/s3456+L(jdu1)*re/propZ3456
      game3456(jdu1,2,1)=Q(jdu1)*qe/s3456+R(jdu1)*le/propZ3456
      game3456(jdu1,2,2)=Q(jdu1)*qe/s3456+R(jdu1)*re/propZ3456

      gamV(jdu1,1)=Q(jdu1)/s3456+L(jdu1)*rxw/propZ3456
      gamV(jdu1,2)=Q(jdu1)/s3456+R(jdu1)*rxw/propZ3456

      gamn3456(jdu1,1,1)=Q(jdu1)*qn/s3456+L(jdu1)*ln/propZ3456
      gamn3456(jdu1,1,2)=Q(jdu1)*qn/s3456+L(jdu1)*rn/propZ3456
      gamn3456(jdu1,2,1)=Q(jdu1)*qn/s3456+R(jdu1)*ln/propZ3456
      gamn3456(jdu1,2,2)=Q(jdu1)*qn/s3456+R(jdu1)*rn/propZ3456
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
      amp(1,jdu2,1,h28)= + propw56**(-1)*propw34**(-1)*cxw**(-1)*
     & s28**(-1)*s128**(-1) * (  - 2.D0*za(p7,p5)*za(p7,p3)*zb(p7,p6)*
     &    zb(p1,p2)*zba2(p4,p1,p2,p8)*s567**(-1) - 2.D0*za(p7,p5)*za(p5
     &    ,p3)*zb(p1,p2)*zb(p5,p6)*zba2(p4,p1,p2,p8)*s567**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + propw56**(-1)*
     & propw34**(-1)*cxw**(-1)*s28**(-1)*s278**(-1) * (  - 2.D0*za(p7,
     &    p8)*za(p1,p3)*zb(p1,p6)*zb(p1,p4)*zba2(p2,p7,p8,p5)*
     &    s134**(-1) - 2.D0*za(p7,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*
     &    zba2(p2,p7,p8,p5)*s134**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + propw56**(-1)*
     & propw34**(-1)*cxw**(-1)*s28**(-1) * ( 2.D0*za(p7,p5)*za(p7,p8)*
     &    za(p1,p3)*zb(p7,p6)*zb(p1,p4)*zb(p1,p2)*s134**(-1)*s567**(-1)
     &     - 2.D0*za(p7,p5)*za(p7,p8)*za(p3,p4)*zb(p7,p6)*zb(p1,p4)*zb(
     &    p4,p2)*s134**(-1)*s567**(-1) + 2.D0*za(p7,p5)*za(p1,p3)*za(p5
     &    ,p8)*zb(p1,p4)*zb(p1,p2)*zb(p5,p6)*s134**(-1)*s567**(-1) - 2.D
     &    0*za(p7,p5)*za(p5,p8)*za(p3,p4)*zb(p1,p4)*zb(p5,p6)*zb(p4,p2)
     &    *s134**(-1)*s567**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + propw56**(-1)*s28**(-1)*
     & s128**(-1) * ( 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p6)*
     &    zba2(p4,p1,p2,p8)*gamn3456(1,1,1)*s356**(-1) + 4.D0*za(p7,p3)
     &    *za(p5,p6)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)*game3456(1,1
     &    ,1)*s456**(-1) - 4.D0*za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p6,p3)
     &    *zba2(p4,p1,p2,p8)*gamn3456(1,1,1)*s356**(-1) + 4.D0*za(p7,p3
     &    )*za(p5,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)*game3456(1,
     &    1,1)*s456**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + propw56**(-1)*s28**(-1)*
     & s278**(-1) * (  - 4.D0*za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*
     &    zba2(p2,p7,p8,p3)*game3456(1,1,1)*s456**(-1) - 4.D0*za(p7,p8)
     &    *za(p5,p3)*zb(p1,p4)*zb(p5,p6)*zba2(p2,p7,p8,p5)*gamn3456(1,1
     &    ,1)*s356**(-1) + 4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p6,p3)
     &    *zba2(p2,p7,p8,p3)*gamn3456(1,1,1)*s356**(-1) - 4.D0*za(p7,p8
     &    )*za(p5,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p3)*game3456(1,
     &    1,1)*s456**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + propw34**(-1)*s28**(-1)*
     & s128**(-1) * (  - 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p4)*
     &    zba2(p6,p1,p2,p8)*game3456(1,1,1)*s345**(-1) + 4.D0*za(p7,p5)
     &    *za(p6,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)*gamn3456(1,1
     &    ,1)*s346**(-1) - 4.D0*za(p7,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)
     &    *zba2(p4,p1,p2,p8)*gamn3456(1,1,1)*s346**(-1) - 4.D0*za(p7,p3
     &    )*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2,p8)*game3456(1,
     &    1,1)*s345**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + propw34**(-1)*s28**(-1)*
     & s278**(-1) * ( 4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p5,p4)*
     &    zba2(p2,p7,p8,p5)*game3456(1,1,1)*s345**(-1) + 4.D0*za(p7,p8)
     &    *za(p5,p3)*zb(p1,p6)*zb(p3,p4)*zba2(p2,p7,p8,p3)*game3456(1,1
     &    ,1)*s345**(-1) - 4.D0*za(p7,p8)*za(p6,p3)*zb(p1,p6)*zb(p6,p4)
     &    *zba2(p2,p7,p8,p5)*gamn3456(1,1,1)*s346**(-1) + 4.D0*za(p7,p8
     &    )*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5)*gamn3456(1,
     &    1,1)*s346**(-1) )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv(1,1)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5)*zb(p1,
     &    p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8) + 4.D0*za(p7,p3)*zb(
     &    p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8) - 2.D0*za(p1,p8)*
     &    za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0*za(p1,
     &    p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*
     &    za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2) - 2.
     &    D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p2)
     &     )
      amp(1,jdu2,1,h28) = amp(1,jdu2,1,h28) + gamv(1,1)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*za(p5,p3)*
     &    zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,p8)*za(p5,
     &    p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1) + 4.D0*
     &    za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,p5) - 4.D0
     &    *za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,p8,p3) )
      amp(1,jdu2,2,h28)= + propw56**(-1)*s28**(-1)*s128**(-1) * ( 4.D0*
     &    za(p1,p8)*za(p5,p6)*zb(p7,p6)*zb(p6,p4)*zab2(p3,p1,p8,p2)*
     &    game3456(1,2,1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,
     &    p4)*zb(p5,p6)*zab2(p5,p1,p8,p2)*gamn3456(1,2,1)*s356**(-1) - 
     &    4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p4)*zb(p6,p3)*zab2(p3,p1,p8,p2
     &    )*gamn3456(1,2,1)*s356**(-1) + 4.D0*za(p1,p8)*za(p5,p4)*zb(p7
     &    ,p4)*zb(p6,p4)*zab2(p3,p1,p8,p2)*game3456(1,2,1)*s456**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + propw56**(-1)*s28**(-1)*
     & s278**(-1) * (  - 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,p2)*zb(p5,p6)*
     &    zab2(p8,p2,p7,p4)*gamn3456(1,2,1)*s356**(-1) - 4.D0*za(p1,p3)
     &    *za(p5,p6)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6)*game3456(1,2
     &    ,1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p3)*zb(p7,p2)*zb(p6,p3)
     &    *zab2(p8,p2,p7,p4)*gamn3456(1,2,1)*s356**(-1) - 4.D0*za(p1,p3
     &    )*za(p5,p4)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p4)*game3456(1,
     &    2,1)*s456**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + propw34**(-1)*s28**(-1)*
     & s128**(-1) * (  - 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p6)*zb(p5,p4)*
     &    zab2(p5,p1,p8,p2)*game3456(1,2,1)*s345**(-1) - 4.D0*za(p1,p8)
     &    *za(p5,p3)*zb(p7,p6)*zb(p3,p4)*zab2(p3,p1,p8,p2)*game3456(1,2
     &    ,1)*s345**(-1) + 4.D0*za(p1,p8)*za(p6,p3)*zb(p7,p6)*zb(p6,p4)
     &    *zab2(p5,p1,p8,p2)*gamn3456(1,2,1)*s346**(-1) - 4.D0*za(p1,p8
     &    )*za(p3,p4)*zb(p7,p4)*zb(p6,p4)*zab2(p5,p1,p8,p2)*gamn3456(1,
     &    2,1)*s346**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + propw34**(-1)*s28**(-1)*
     & s278**(-1) * ( 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,p2)*zb(p5,p4)*
     &    zab2(p8,p2,p7,p6)*game3456(1,2,1)*s345**(-1) - 4.D0*za(p1,p5)
     &    *za(p6,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6)*gamn3456(1,2
     &    ,1)*s346**(-1) + 4.D0*za(p1,p5)*za(p3,p4)*zb(p7,p2)*zb(p6,p4)
     &    *zab2(p8,p2,p7,p4)*gamn3456(1,2,1)*s346**(-1) + 4.D0*za(p1,p3
     &    )*za(p5,p3)*zb(p7,p2)*zb(p3,p4)*zab2(p8,p2,p7,p6)*game3456(1,
     &    2,1)*s345**(-1) )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamv(1,2)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s128**(-1) * (  - 2.D0*za(p1,p8)*za(p5,
     &    p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) + 2.D0*za(p1,p8)*
     &    za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) - 2.D0*za(p1,
     &    p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p5,p6,p8) + 2.D0*
     &    za(p1,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p3,p4,p8) - 4.
     &    D0*za(p1,p8)*zb(p7,p6)*zab2(p5,p1,p8,p2)*zab2(p3,p5,p6,p4) + 
     &    4.D0*za(p1,p8)*zb(p7,p4)*zab2(p5,p3,p4,p6)*zab2(p3,p1,p8,p2)
     &     )
      amp(1,jdu2,2,h28) = amp(1,jdu2,2,h28) + gamv(1,2)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*za(p5,p3)*
     &    zb(p7,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) - 2.D0*za(p7,p8)*za(p5,
     &    p3)*zb(p7,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) + 4.D0*za(p1,p5)*
     &    zb(p7,p2)*zab2(p3,p5,p6,p4)*zab2(p8,p2,p7,p6) - 4.D0*za(p1,p3
     &    )*zb(p7,p2)*zab2(p5,p3,p4,p6)*zab2(p8,p2,p7,p4) - 2.D0*za(p5,
     &    p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p5,p6,p1) + 2.D0*
     &    za(p5,p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p3,p4,p1) )
      amp(2,jdu2,1,h28)= + propw56**(-1)*propw34**(-1)*cxw**(-1)*
     & s28**(-1)*s128**(-1) * (  - 2.D0*za(p7,p5)*za(p7,p3)*zb(p7,p4)*
     &    zb(p1,p2)*zba2(p6,p1,p2,p8)*s347**(-1) + 2.D0*za(p7,p3)*za(p5
     &    ,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2,p8)*s347**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + propw56**(-1)*
     & propw34**(-1)*cxw**(-1)*s28**(-1)*s278**(-1) * (  - 2.D0*za(p7,
     &    p8)*za(p1,p5)*zb(p1,p6)*zb(p1,p4)*zba2(p2,p7,p8,p3)*
     &    s156**(-1) + 2.D0*za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*
     &    zba2(p2,p7,p8,p3)*s156**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + propw56**(-1)*
     & propw34**(-1)*cxw**(-1)*s28**(-1) * ( 2.D0*za(p7,p3)*za(p7,p8)*
     &    za(p1,p5)*zb(p7,p4)*zb(p1,p6)*zb(p1,p2)*s347**(-1)*s156**(-1)
     &     - 2.D0*za(p7,p3)*za(p7,p8)*za(p5,p6)*zb(p7,p4)*zb(p1,p6)*zb(
     &    p6,p2)*s347**(-1)*s156**(-1) + 2.D0*za(p7,p3)*za(p1,p5)*za(p3
     &    ,p8)*zb(p1,p6)*zb(p1,p2)*zb(p3,p4)*s347**(-1)*s156**(-1) - 2.D
     &    0*za(p7,p3)*za(p5,p6)*za(p3,p8)*zb(p1,p6)*zb(p6,p2)*zb(p3,p4)
     &    *s347**(-1)*s156**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + propw56**(-1)*s28**(-1)*
     & s128**(-1) * ( 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p6)*
     &    zba2(p4,p1,p2,p8)*gamn3456(2,1,1)*s356**(-1) + 4.D0*za(p7,p3)
     &    *za(p5,p6)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)*game3456(2,1
     &    ,1)*s456**(-1) - 4.D0*za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p6,p3)
     &    *zba2(p4,p1,p2,p8)*gamn3456(2,1,1)*s356**(-1) + 4.D0*za(p7,p3
     &    )*za(p5,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)*game3456(2,
     &    1,1)*s456**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + propw56**(-1)*s28**(-1)*
     & s278**(-1) * (  - 4.D0*za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*
     &    zba2(p2,p7,p8,p3)*game3456(2,1,1)*s456**(-1) - 4.D0*za(p7,p8)
     &    *za(p5,p3)*zb(p1,p4)*zb(p5,p6)*zba2(p2,p7,p8,p5)*gamn3456(2,1
     &    ,1)*s356**(-1) + 4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p6,p3)
     &    *zba2(p2,p7,p8,p3)*gamn3456(2,1,1)*s356**(-1) - 4.D0*za(p7,p8
     &    )*za(p5,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p3)*game3456(2,
     &    1,1)*s456**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + propw34**(-1)*s28**(-1)*
     & s128**(-1) * (  - 4.D0*za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p4)*
     &    zba2(p6,p1,p2,p8)*game3456(2,1,1)*s345**(-1) + 4.D0*za(p7,p5)
     &    *za(p6,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)*gamn3456(2,1
     &    ,1)*s346**(-1) - 4.D0*za(p7,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)
     &    *zba2(p4,p1,p2,p8)*gamn3456(2,1,1)*s346**(-1) - 4.D0*za(p7,p3
     &    )*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2,p8)*game3456(2,
     &    1,1)*s345**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + propw34**(-1)*s28**(-1)*
     & s278**(-1) * ( 4.D0*za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p5,p4)*
     &    zba2(p2,p7,p8,p5)*game3456(2,1,1)*s345**(-1) + 4.D0*za(p7,p8)
     &    *za(p5,p3)*zb(p1,p6)*zb(p3,p4)*zba2(p2,p7,p8,p3)*game3456(2,1
     &    ,1)*s345**(-1) - 4.D0*za(p7,p8)*za(p6,p3)*zb(p1,p6)*zb(p6,p4)
     &    *zba2(p2,p7,p8,p5)*gamn3456(2,1,1)*s346**(-1) + 4.D0*za(p7,p8
     &    )*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5)*gamn3456(2,
     &    1,1)*s346**(-1) )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv(2,1)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s128**(-1) * (  - 4.D0*za(p7,p5)*zb(p1,
     &    p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8) + 4.D0*za(p7,p3)*zb(
     &    p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8) - 2.D0*za(p1,p8)*
     &    za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) + 2.D0*za(p1,
     &    p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*
     &    za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2) - 2.
     &    D0*za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p2)
     &     )
      amp(2,jdu2,1,h28) = amp(2,jdu2,1,h28) + gamv(2,1)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*za(p5,p3)*
     &    zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p7,p8)*za(p5,
     &    p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1) + 2.D0*za(p7,p8)*
     &    za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1) - 2.D0*za(p7,
     &    p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1) + 4.D0*
     &    za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,p5) - 4.D0
     &    *za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,p8,p3) )
      amp(2,jdu2,2,h28)= + propw56**(-1)*s28**(-1)*s128**(-1) * ( 4.D0*
     &    za(p1,p8)*za(p5,p6)*zb(p7,p6)*zb(p6,p4)*zab2(p3,p1,p8,p2)*
     &    game3456(2,2,1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,
     &    p4)*zb(p5,p6)*zab2(p5,p1,p8,p2)*gamn3456(2,2,1)*s356**(-1) - 
     &    4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p4)*zb(p6,p3)*zab2(p3,p1,p8,p2
     &    )*gamn3456(2,2,1)*s356**(-1) + 4.D0*za(p1,p8)*za(p5,p4)*zb(p7
     &    ,p4)*zb(p6,p4)*zab2(p3,p1,p8,p2)*game3456(2,2,1)*s456**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + propw56**(-1)*s28**(-1)*
     & s278**(-1) * (  - 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,p2)*zb(p5,p6)*
     &    zab2(p8,p2,p7,p4)*gamn3456(2,2,1)*s356**(-1) - 4.D0*za(p1,p3)
     &    *za(p5,p6)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6)*game3456(2,2
     &    ,1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p3)*zb(p7,p2)*zb(p6,p3)
     &    *zab2(p8,p2,p7,p4)*gamn3456(2,2,1)*s356**(-1) - 4.D0*za(p1,p3
     &    )*za(p5,p4)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p4)*game3456(2,
     &    2,1)*s456**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + propw34**(-1)*s28**(-1)*
     & s128**(-1) * (  - 4.D0*za(p1,p8)*za(p5,p3)*zb(p7,p6)*zb(p5,p4)*
     &    zab2(p5,p1,p8,p2)*game3456(2,2,1)*s345**(-1) - 4.D0*za(p1,p8)
     &    *za(p5,p3)*zb(p7,p6)*zb(p3,p4)*zab2(p3,p1,p8,p2)*game3456(2,2
     &    ,1)*s345**(-1) + 4.D0*za(p1,p8)*za(p6,p3)*zb(p7,p6)*zb(p6,p4)
     &    *zab2(p5,p1,p8,p2)*gamn3456(2,2,1)*s346**(-1) - 4.D0*za(p1,p8
     &    )*za(p3,p4)*zb(p7,p4)*zb(p6,p4)*zab2(p5,p1,p8,p2)*gamn3456(2,
     &    2,1)*s346**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + propw34**(-1)*s28**(-1)*
     & s278**(-1) * ( 4.D0*za(p1,p5)*za(p5,p3)*zb(p7,p2)*zb(p5,p4)*
     &    zab2(p8,p2,p7,p6)*game3456(2,2,1)*s345**(-1) - 4.D0*za(p1,p5)
     &    *za(p6,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p8,p2,p7,p6)*gamn3456(2,2
     &    ,1)*s346**(-1) + 4.D0*za(p1,p5)*za(p3,p4)*zb(p7,p2)*zb(p6,p4)
     &    *zab2(p8,p2,p7,p4)*gamn3456(2,2,1)*s346**(-1) + 4.D0*za(p1,p3
     &    )*za(p5,p3)*zb(p7,p2)*zb(p3,p4)*zab2(p8,p2,p7,p6)*game3456(2,
     &    2,1)*s345**(-1) )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamv(2,2)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s128**(-1) * (  - 2.D0*za(p1,p8)*za(p5,
     &    p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) + 2.D0*za(p1,p8)*
     &    za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) - 2.D0*za(p1,
     &    p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p5,p6,p8) + 2.D0*
     &    za(p1,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zba2(p7,p3,p4,p8) - 4.
     &    D0*za(p1,p8)*zb(p7,p6)*zab2(p5,p1,p8,p2)*zab2(p3,p5,p6,p4) + 
     &    4.D0*za(p1,p8)*zb(p7,p4)*zab2(p5,p3,p4,p6)*zab2(p3,p1,p8,p2)
     &     )
      amp(2,jdu2,2,h28) = amp(2,jdu2,2,h28) + gamv(2,2)*propw56**(-1)*
     & propw34**(-1)*s28**(-1)*s278**(-1) * ( 2.D0*za(p7,p8)*za(p5,p3)*
     &    zb(p7,p2)*zb(p6,p4)*zba2(p7,p5,p6,p1) - 2.D0*za(p7,p8)*za(p5,
     &    p3)*zb(p7,p2)*zb(p6,p4)*zba2(p7,p3,p4,p1) + 4.D0*za(p1,p5)*
     &    zb(p7,p2)*zab2(p3,p5,p6,p4)*zab2(p8,p2,p7,p6) - 4.D0*za(p1,p3
     &    )*zb(p7,p2)*zab2(p5,p3,p4,p6)*zab2(p8,p2,p7,p4) - 2.D0*za(p5,
     &    p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p5,p6,p1) + 2.D0*
     &    za(p5,p3)*za(p8,p2)*zb(p7,p2)*zb(p6,p4)*zba2(p2,p3,p4,p1) )
      enddo
      enddo
      amp(:,:,:,:)=amp(:,:,:,:)/cxw
      return
      end
