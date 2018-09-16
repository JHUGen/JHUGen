      subroutine qq2l2nuggamp(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      integer:: jdu,j78,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8,h1,h7,h8
      complex(dp)::ratxw,zab2,zba2,coupel(2,2),coupnu(2,2),
     & coupw(2,2),
     & adr78(2,2,2,2),adr87(2,2,2,2),a78(2,2,2,2),
     & asr78(2,2,2,2),asr87(2,2,2,2),a87(2,2,2,2),aq(2,2,2,2),
     & propw34,propw56,propz3456,propa3456,
     & zba3,zba4,iza,izb,sign
      real(dp):: t3,t4,s356,s456,s34,s56,
     & s278,s178,s156,s234,s345,s346,s2347,s3456,msq(2)
C-----Begin statement functions
      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      zba3(i1,i2,i3,i4,i5)=
     & zb(i1,i2)*za(i2,i5)+zb(i1,i3)*za(i3,i5)+zb(i1,i4)*za(i4,i5)
      zba4(i1,i2,i3,i4,i5,i6)=
     & +zb(i1,i2)*za(i2,i6)+zb(i1,i3)*za(i3,i6)
     & +zb(i1,i4)*za(i4,i6)+zb(i1,i5)*za(i5,i6)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      p1=i1
      p2=i2
      s3456=t4(i3,i4,i5,i6)
      propz3456=s3456-czmass2
      propa3456=s3456
      ratxw=sqrt((cone-cxw)/cxw)
      do jdu=1,2
      coupel(jdu,1)=Q(jdu)*qe/propa3456+L(jdu)*le/propz3456
      coupel(jdu,2)=Q(jdu)*qe/propa3456+R(jdu)*le/propz3456
      coupnu(jdu,1)=L(jdu)*ln/propz3456
      coupnu(jdu,2)=R(jdu)*ln/propz3456
      coupw(jdu,1)=ratxw*L(jdu)/propz3456+Q(jdu)/propa3456
      coupw(jdu,2)=ratxw*R(jdu)/propz3456+Q(jdu)/propa3456
      enddo

      adr78=czip
      adr87=czip
      asr78=czip
      asr87=czip

      do jdu=1,2
      p3=i3
      p4=i4
      p5=i5
      p6=i6
      p7=i7
      p8=i8
      s34=s(p3,p4)
      s56=s(p5,p6)
      propw34=s34-cwmass2
      propw56=s56-cwmass2
      s278=t3(p2,p7,p8)
      s178=t3(p1,p7,p8)
      s156=t3(p1,p5,p6)
      s234=t3(p2,p3,p4)
      s356=t3(p3,p5,p6)
      s345=t3(p3,p4,p5)
      s346=t3(p3,p4,p6)
      s456=t3(p4,p5,p6)

      do j78=1,2
      if (j78 == 1) then
      p7=i7
      p8=i8
      s2347=t4(p2,p3,p4,p7)

      asr78(jdu,1,1,1)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)*izb(
     &    p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*za(
     &    p2,p7)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*za(p2,p8)*za(p3,p4)*
     &    zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*s346**(-1) + 2.D0*za(p2,p8)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*
     &    za(p3,p4)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*
     &    za(p3,p6)*zb(p1,p6)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) )
      asr78(jdu,1,1,1) = asr78(jdu,1,1,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*
     &    s356**(-1) - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*
     &    izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s356**(-1) - 2.D0*
     &    za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p3)*s356**(-1) - 2.D0*za(p2,p8)*za(p3,p5)*
     &    zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*s356**(-1) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*zba2(p1,p2,
     &    p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*
     &    s356**(-1) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*
     &    s356**(-1) )
      asr78(jdu,1,1,1) = asr78(jdu,1,1,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(
     &    p3,p4)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*s345**(-1)
     &     - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s345**(-1) + 2.D0*za(p2,p8)*
     &    za(p3,p5)*zb(p1,p6)*zb(p3,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1,
     &    p2,p7,p8,p3)*s345**(-1) - 2.D0*za(p2,p8)*za(p3,p5)*zb(p1,p6)*
     &    zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*
     &    s345**(-1) + 2.D0*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*
     &    s345**(-1) - 2.D0*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*
     &    s345**(-1) )
      asr78(jdu,1,1,1) = asr78(jdu,1,1,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * ( 2.D0*za(p2,p7)*za(p4,p5)*
     &    zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3
     &    ) - 2.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p8)*za(p4,p5)*
     &    zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3
     &    ) - 2.D0*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p4,p5)*zb(p1,p4)*
     &    zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p3) - 2.D0*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) )
      asr78(jdu,1,1,2)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p4)*zb(p1,p8)*
     &    zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p6)*zb(p1,p8)*zb(
     &    p4,p6)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p4)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s346**(-1) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s346**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p4)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(
     &    p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p3,p6)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(
     &    p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(
     &    p1,p8)*zb(p3,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,
     &    p8)*zb(p5,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(
     &    p1,p4)*zb(p3,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,
     &    p4)*zb(p5,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s356**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p3,p5)*zb(p3,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *s356**(-1) - 2.D0*za(p3,p5)*zb(p5,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *s356**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p3,p4)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p4,p5)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p3,p4)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p5)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s345**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p5)*zb(p3,p4)*zba2(p6,p1,p8,p7)*zba2(
     &    p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p3,p5)*zb(p4,p5)*zba2(p6,p1,p8,p7)*zba2(
     &    p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)
     &    *za(p4,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p5)
     &    *zb(p2,p8)*zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8) + 2.D0*za(p2,p7)*za(p5,p6)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) )
      asr78(jdu,1,1,2) = asr78(jdu,1,1,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p4,p5)*zb(p4,p6)*zba2(p4,p1,p8
     &    ,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8) - 2.D0*za(p5,p6)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(p8,
     &    p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      asr78(jdu,1,2,1)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p7)**2*zb(p4,p6
     &    )*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p2,p5)*za(p3,p6)*zb(p1,p7)**2*zb(p4,p6)*
     &    zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s346**(-1)
     &     )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p4)*zb(p1,
     &    p4)*zb(p4,p6)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*s346**(-1) - 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p6)*zb(p4,
     &    p6)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p4)*zb(p1,p4)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p6)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p5,p6)*
     &    zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s356**(-1)
     &     )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*
     &    zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s356**(-1)
     &     )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*
     &    zb(p1,p7)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s356**(-1) - 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p4)
     &    *zb(p1,p7)*zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)
     &    **2*zb(p3,p4)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p4,
     &    p5)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,
     &    p6)*zb(p3,p4)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p6)*zb(p4,
     &    p5)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p3,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)
     &    *zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,
     &    p8)*izb(p7,p8) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)**2*zb(p4,
     &    p6)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)**2*za(p4,
     &    p5)*zb(p1,p4)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,
     &    p8)*izb(p7,p8) + 2.D0*za(p2,p8)**2*za(p5,p6)*zb(p1,p6)*zb(p4,
     &    p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      asr78(jdu,1,2,1) = asr78(jdu,1,2,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p4,p5)*zb(
     &    p1,p4)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p3)*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     )
      asr78(jdu,1,2,2)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s346**(-1) - 2.D0*
     &    za(p2,p5)*za(p3,p4)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *zba3(p4,p1,p7,p8,p2)*s346**(-1) + 2.D0*za(p2,p5)*za(p3,p4)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s346**(-1) - 2.D0*za(p2,p5)*za(p3,p6)*
     &    zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2
     &    )*s346**(-1) - 2.D0*za(p2,p5)*za(p3,p6)*zb(p1,p8)*zb(p4,p6)*
     &    iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s346**(-1) + 2.D0*
     &    za(p2,p5)*za(p3,p6)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*s346**(-1) )
      asr78(jdu,1,2,2) = asr78(jdu,1,2,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(
     &    p3,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s356**(-1)
     &     + 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p2,p7)*
     &    iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s356**(-1) - 2.D0*za(p2,p3)*
     &    za(p3,p5)*zb(p3,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*
     &    iza(p1,p8)*zba3(p4,p1,p7,p8,p2)*s356**(-1) + 2.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p7)*zb(p5,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2)*s356**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p8)*
     &    zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*
     &    s356**(-1) - 2.D0*za(p2,p5)*za(p3,p5)*zb(p5,p6)*zba2(p7,p1,p8
     &    ,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)*
     &    s356**(-1) )
      asr78(jdu,1,2,2) = asr78(jdu,1,2,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*
     &    zb(p3,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*
     &    s345**(-1) - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p4)*
     &    iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) + 2.D0*
     &    za(p2,p3)*za(p3,p5)*zb(p3,p4)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) + 2.D0*
     &    za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,p8)
     &    *zba3(p6,p1,p7,p8,p2)*s345**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*
     &    zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2
     &    )*s345**(-1) - 2.D0*za(p2,p5)*za(p3,p5)*zb(p4,p5)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*
     &    s345**(-1) )
      asr78(jdu,1,2,2) = asr78(jdu,1,2,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)
     &    *zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,
     &    p2) - 2.D0*za(p2,p3)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)
     &    *iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(
     &    p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) + 2.D0*za(
     &    p2,p3)*za(p5,p6)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p6,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p4,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1
     &    ,p7,p8,p2) )
      asr78(jdu,2,1,1)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p1,p7)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*s346**(-1) - 2.D0*za(
     &    p1,p7)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1
     &    ,p8)*izb(p7,p8)*s346**(-1) - 2.D0*za(p1,p8)*za(p3,p4)*zb(p2,
     &    p4)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p1,p8)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*s346**(-1) - 2.D0*za(
     &    p3,p4)*za(p7,p8)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8)*s346**(-1) - 2.D0*za(p3,p6)*za(p7,p8)*zb(p2,
     &    p6)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*
     &    s346**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p1)
     &    **2*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,
     &    p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p1)*zb(p4,p6)
     &    *zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*s346**(-1) + 2.D0*
     &    za(p1,p5)*za(p3,p4)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*
     &    izb(p1,p7)*izb(p7,p8)*s346**(-1) - 2.D0*za(p1,p5)*za(p3,p6)*
     &    zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p1)*
     &    zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*s346**(-1)
     &     + 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*s346**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p4)*za(p5,p8)*zb(p2,p1)*zb(p2,p4)*zb(
     &    p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*s346**(-1) + 2.D0*za(
     &    p3,p6)*za(p5,p8)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*s346**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(
     &    p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*s356**(-1) + 2.
     &    D0*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*
     &    izb(p1,p8)*izb(p7,p8)*s356**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*
     &    zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*s356**(-1) + 2.D0*za(
     &    p3,p5)*za(p7,p8)*zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8)*s356**(-1) + 2.D0*za(p3,p5)*za(p7,p8)*zb(p2,
     &    p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*
     &    s356**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)**2*
     &    zb(p3,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    s356**(-1) - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p6)*
     &    zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*s356**(-1) - 2.D0*za(
     &    p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p4)*izb(p1
     &    ,p7)*izb(p7,p8)*s356**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,
     &    p1)**2*zb(p5,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(
     &    p1,p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p5,
     &    p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*s356**(-1) - 2.D0
     &    *za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p5,p6)*zab2(p8,p2,p7,p4)*
     &    izb(p1,p7)*izb(p7,p8)*s356**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,p4)*
     &    zb(p3,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*s356**(-1) - 2.D0*
     &    za(p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p4)*zb(p5,p6)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*s356**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p7)*za(p3,p5)*zb(p2,p6)*
     &    zb(p3,p4)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*s345**(-1)
     &     + 2.D0*za(p1,p7)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)*zab2(p5,p7,p8
     &    ,p1)*izb(p1,p8)*izb(p7,p8)*s345**(-1) - 2.D0*za(p1,p8)*za(p3,
     &    p5)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)
     &    *zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*s345**(-1) - 2.D0*
     &    za(p3,p5)*za(p7,p8)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p7,p8,p1)*
     &    izb(p1,p7)*izb(p1,p8)*s345**(-1) + 2.D0*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p5)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*
     &    s345**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)
     &    **2*zb(p3,p4)*zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,
     &    p8)*s345**(-1) + 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p4)
     &    *zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*s345**(-1) + 2.D0*
     &    za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p4)*zab2(p8,p2,p7,p6)*
     &    izb(p1,p7)*izb(p7,p8)*s345**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*
     &    zb(p2,p1)**2*zb(p4,p5)*zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*s345**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p1)*
     &    zb(p4,p5)*zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*s345**(-1)
     &     - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p4,p5)*zab2(p8,p2,p7
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*s345**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,p6)*zb(
     &    p3,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*s345**(-1) - 2.D0*za(
     &    p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p6)*zb(p4,p5)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*s345**(-1) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p1,p7)*za(p4,p5)
     &    *zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)
     &     + 2.D0*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8
     &    ,p1)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(p1,p8)*za(p4,p5)*zb(p2,
     &    p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,p8) + 2.D0*
     &    za(p1,p8)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*
     &    izb(p1,p7)*izb(p7,p8) - 2.D0*za(p4,p5)*za(p7,p8)*zb(p2,p4)*
     &    zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p1,p8) + 2.D0*za(
     &    p5,p6)*za(p7,p8)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p1,p3)*za(p4,p5)
     &    *zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,
     &    p7)*izb(p1,p8) + 2.D0*za(p1,p3)*za(p4,p5)*zb(p2,p1)*zb(p4,p6)
     &    *zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8) + 2.D0*za(p1,p3)*za(
     &    p4,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(
     &    p7,p8) + 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p1)**2*zb(p4,p6)*
     &    zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8) - 2.D0*za(
     &    p1,p3)*za(p5,p6)*zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1
     &    ,p8)*izb(p7,p8) - 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p1)*zb(p4,p6
     &    )*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p7,p8) )
      asr78(jdu,2,1,1) = asr78(jdu,2,1,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p3,p8)*za(p4,p5)*zb(p2,p1)*zb(
     &    p2,p4)*zb(p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8) - 2.D0*za(
     &    p3,p8)*za(p5,p6)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8) )
      asr78(jdu,2,2,1)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p1,p8)*za(p3,p4)*zb(p2,p4)*zb(p1,p7)*
     &    zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p1,p8)*za(p3,p6)*zb(p2,p6)*zb(p1,p7)*zb(
     &    p4,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p4)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s346**(-1) - 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p6)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s346**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p4)*zb(p4,p6)*zab2(p5,p1,p8,p7)*zab2(
     &    p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p3,p6)*zb(p4,p6)*zab2(p5,p1,p8,p7)*zab2(
     &    p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(
     &    p1,p7)*zb(p3,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p1,
     &    p7)*zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p8)*za(p1,p3)*za(p3,p5)*zb(
     &    p2,p7)*zb(p3,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p5)*zb(p2,
     &    p7)*zb(p5,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p3,p5)*zb(p3,p6)*zab2(p3,p1,p8,p7)*
     &    zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *s356**(-1) - 2.D0*za(p3,p5)*zb(p5,p6)*zab2(p5,p1,p8,p7)*
     &    zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *s356**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p3,p4)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p4,p5)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8)*s345**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)*za(p1,p3)*za(p3,p5)*
     &    zb(p2,p7)*zb(p3,p4)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p5)*
     &    zb(p2,p7)*zb(p4,p5)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p5)*zb(p3,p4)*zab2(p3,p1,p8,p7)*zab2(
     &    p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p3,p5)*zb(p4,p5)*zab2(p5,p1,p8,p7)*zab2(
     &    p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p1,p8)*za(p4,p5)
     &    *zb(p2,p4)*zb(p1,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*
     &    izb(p1,p8)*izb(p7,p8) + 2.D0*za(p1,p8)*za(p5,p6)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)*za(p1,p3)
     &    *za(p4,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p8)*za(p1,p3)*za(p5,p6)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      asr78(jdu,2,2,1) = asr78(jdu,2,2,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p4,p5)*zb(p4,p6)*zab2(p3,p1,p8
     &    ,p7)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(
     &    p7,p8) - 2.D0*za(p5,p6)*zb(p4,p6)*zab2(p3,p1,p8,p7)*zab2(p8,
     &    p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      asr78(jdu,2,1,2)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p1,p7)**2*za(p3,p4)*zb(p2,p4)*zb(p4,p6
     &    )*zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p1,p7)**2*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s346**(-1)
     &     )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p8)
     &    **2*zb(p4,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s346**(-1) - 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p8)**2*zb(p4,
     &    p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p4)*zb(p2,p4)*zb(
     &    p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p6)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*
     &    zb(p3,p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s356**(-1)
     &     )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)**2*
     &    zb(p3,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p8)**2*zb(p5,p6)*
     &    zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*s356**(-1)
     &     )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*
     &    zb(p2,p8)*zb(p3,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p4)
     &    *zb(p2,p8)*zb(p5,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s356**(-1) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,
     &    p6)*zb(p3,p4)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,p6)*zb(p4,
     &    p5)*zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)
     &    **2*zb(p3,p4)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p8)**2*zb(p4,
     &    p5)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p3,p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p4,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p1,p7)**2*za(p4,
     &    p5)*zb(p2,p4)*zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8) + 2.D0*za(p1,p7)**2*za(p5,p6)*zb(p2,p6)*zb(p4,
     &    p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p1,p3)*za(p4,p5)
     &    *zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8) + 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p8)**2*zb(p4,
     &    p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      asr78(jdu,2,1,2) = asr78(jdu,2,1,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p1,p3)*za(p1,p7)*za(p4,p5)*zb(
     &    p2,p4)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) - 2.D0*za(p1,p3)*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     )
      asr78(jdu,2,2,2)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p1)**2*za(p3,p4)*zb(p2,p4)*zb(p4,p6
     &    )*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    s346**(-1) - 2.D0*za(p2,p1)**2*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s346**(-1)
     &     + 2.D0*za(p2,p1)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p1,p7
     &    ,p8)*iza(p2,p7)*iza(p7,p8)*s346**(-1) + 2.D0*za(p2,p1)*za(p3,
     &    p4)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,
     &    p8)*s346**(-1) + 2.D0*za(p2,p1)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)
     &    *zab2(p5,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*s346**(-1) + 2.D0*
     &    za(p2,p1)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p5,p1,p8,p7)*
     &    iza(p2,p8)*iza(p7,p8)*s346**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p7)*zb(
     &    p4,p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*s346**(-1) + 2.
     &    D0*za(p1,p5)*za(p3,p4)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p4)*
     &    iza(p2,p7)*iza(p7,p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p4)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8)*
     &    s346**(-1) + 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p7)*zb(p4,p6)*
     &    zab2(p2,p7,p8,p6)*iza(p2,p8)*iza(p7,p8)*s346**(-1) + 2.D0*za(
     &    p1,p5)*za(p3,p6)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p6)*iza(p2
     &    ,p7)*iza(p7,p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p6)*zb(p4,
     &    p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p2,p8)*
     &    s346**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p2,p1)*za(p1,p5)*za(p3,p4)*zb(p4,p6)*
     &    zb(p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s346**(-1) - 2.D0*
     &    za(p2,p1)*za(p1,p5)*za(p3,p6)*zb(p4,p6)*zb(p6,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*s346**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*
     &    zb(p3,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    s356**(-1) + 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s356**(-1)
     &     - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(p3,p6)*zab2(p3,p1,p7
     &    ,p8)*iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*za(p2,p1)*za(p3,
     &    p5)*zb(p2,p4)*zb(p3,p6)*zab2(p3,p1,p8,p7)*iza(p2,p8)*iza(p7,
     &    p8)*s356**(-1) - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)
     &    *zab2(p5,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*
     &    za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*zab2(p5,p1,p8,p7)*
     &    iza(p2,p8)*iza(p7,p8)*s356**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*s356**(-1)
     &     - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(p3,p6)*zab2(p2,p7,p8
     &    ,p4)*iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*za(p1,p3)*za(p3,
     &    p5)*zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,
     &    p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*zb(p5,p6)
     &    *zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*s356**(-1) - 2.D0*
     &    za(p1,p5)*za(p3,p5)*zb(p2,p8)*zb(p5,p6)*zab2(p2,p7,p8,p4)*
     &    iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8)*
     &    s356**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * ( 2.D0*za(p2,p1)*za(p1,p3)*za(p3,p5)*zb(p3,p6)*zb(
     &    p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s356**(-1) + 2.D0*za(
     &    p2,p1)*za(p1,p5)*za(p3,p5)*zb(p4,p7)*zb(p5,p6)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*s356**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,
     &    p6)*zb(p3,p4)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,
     &    p8)*s345**(-1) + 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p6)*zb(p4,
     &    p5)*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    s345**(-1) + 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p3,p4)*
     &    zab2(p3,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*s345**(-1) + 2.D0*za(
     &    p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p1,p8,p7)*iza(p2
     &    ,p8)*iza(p7,p8)*s345**(-1) - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,
     &    p6)*zb(p4,p5)*zab2(p5,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*
     &    s345**(-1) - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*s345**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*zb(
     &    p3,p4)*zab2(p2,p7,p8,p6)*iza(p2,p8)*iza(p7,p8)*s345**(-1) + 2.
     &    D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(p3,p4)*zab2(p2,p7,p8,p6)*
     &    iza(p2,p7)*iza(p7,p8)*s345**(-1) + 2.D0*za(p1,p3)*za(p3,p5)*
     &    zb(p3,p4)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p2,p8)*
     &    s345**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*zb(p4,p5)*
     &    zab2(p2,p7,p8,p6)*iza(p2,p8)*iza(p7,p8)*s345**(-1) - 2.D0*za(
     &    p1,p5)*za(p3,p5)*zb(p2,p8)*zb(p4,p5)*zab2(p2,p7,p8,p6)*iza(p2
     &    ,p7)*iza(p7,p8)*s345**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p4,
     &    p5)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p2,p8)*
     &    s345**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p2,p1)*za(p1,p3)*za(p3,p5)*zb(p3,p4)*
     &    zb(p6,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s345**(-1) + 2.D0*
     &    za(p2,p1)*za(p1,p5)*za(p3,p5)*zb(p4,p5)*zb(p6,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*s345**(-1) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p1)**2*za(p4,
     &    p5)*zb(p2,p4)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,
     &    p8)*iza(p1,p8) + 2.D0*za(p2,p1)**2*za(p5,p6)*zb(p2,p6)*zb(p4,
     &    p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8) + 2.D0
     &    *za(p2,p1)*za(p4,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p3,p1,p7,p8)*
     &    iza(p2,p7)*iza(p7,p8) + 2.D0*za(p2,p1)*za(p4,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p8)*iza(p7,p8) - 2.D0*za(
     &    p2,p1)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p2
     &    ,p7)*iza(p7,p8) - 2.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*zb(p4,p6
     &    )*zab2(p3,p1,p8,p7)*iza(p2,p8)*iza(p7,p8) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * ( 2.D0*za(p1,p3)*za(p4,p5)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)
     &     + 2.D0*za(p1,p3)*za(p4,p5)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8
     &    ,p4)*iza(p2,p7)*iza(p7,p8) + 2.D0*za(p1,p3)*za(p4,p5)*zb(p4,
     &    p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8) - 2.D0*
     &    za(p1,p3)*za(p5,p6)*zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,p6)*
     &    iza(p2,p8)*iza(p7,p8) - 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p8)*
     &    zb(p4,p6)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p7,p8) - 2.D0*za(
     &    p1,p3)*za(p5,p6)*zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2
     &    ,p7)*iza(p2,p8) )
      asr78(jdu,2,2,2) = asr78(jdu,2,2,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * (  - 2.D0*za(p2,p1)*za(p1,p3)*za(p4,p5)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8) + 2.D0*
     &    za(p2,p1)*za(p1,p3)*za(p5,p6)*zb(p4,p6)*zb(p6,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8) )
      elseif (j78 == 2) then
      p7=i8
      p8=i7
      s2347=t4(p2,p3,p4,p7)

      asr87(jdu,1,1,1)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)*izb(
     &    p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*za(
     &    p2,p7)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*za(p2,p8)*za(p3,p4)*
     &    zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*s346**(-1) + 2.D0*za(p2,p8)*za(p3,p6)*zb(p1,p6)*zb(p4,p6)*
     &    izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*
     &    za(p3,p4)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) + 2.D0*
     &    za(p3,p6)*zb(p1,p6)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s346**(-1) )
      asr87(jdu,1,1,1) = asr87(jdu,1,1,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*
     &    s356**(-1) - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*
     &    izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s356**(-1) - 2.D0*
     &    za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p3)*s356**(-1) - 2.D0*za(p2,p8)*za(p3,p5)*
     &    zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*s356**(-1) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*zba2(p1,p2,
     &    p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*
     &    s356**(-1) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*
     &    s356**(-1) )
      asr87(jdu,1,1,1) = asr87(jdu,1,1,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(
     &    p3,p4)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*s345**(-1)
     &     - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s345**(-1) + 2.D0*za(p2,p8)*
     &    za(p3,p5)*zb(p1,p6)*zb(p3,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1,
     &    p2,p7,p8,p3)*s345**(-1) - 2.D0*za(p2,p8)*za(p3,p5)*zb(p1,p6)*
     &    zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*
     &    s345**(-1) + 2.D0*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*
     &    s345**(-1) - 2.D0*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*
     &    s345**(-1) )
      asr87(jdu,1,1,1) = asr87(jdu,1,1,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * ( 2.D0*za(p2,p7)*za(p4,p5)*
     &    zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3
     &    ) - 2.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p8)*za(p4,p5)*
     &    zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3
     &    ) - 2.D0*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p4,p5)*zb(p1,p4)*
     &    zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p3) - 2.D0*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) )
      asr87(jdu,1,2,1)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p4)*zb(p1,p8)*
     &    zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p6)*zb(p1,p8)*zb(
     &    p4,p6)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p4)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s346**(-1) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s346**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p4)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(
     &    p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p3,p6)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(
     &    p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(
     &    p1,p8)*zb(p3,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,
     &    p8)*zb(p5,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(
     &    p1,p4)*zb(p3,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,
     &    p4)*zb(p5,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s356**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p3,p5)*zb(p3,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *s356**(-1) - 2.D0*za(p3,p5)*zb(p5,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *s356**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p3,p4)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p4,p5)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p3,p4)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p5)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*s345**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p5)*zb(p3,p4)*zba2(p6,p1,p8,p7)*zba2(
     &    p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p3,p5)*zb(p4,p5)*zba2(p6,p1,p8,p7)*zba2(
     &    p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)
     &    *za(p4,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p5)
     &    *zb(p2,p8)*zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8) + 2.D0*za(p2,p7)*za(p5,p6)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) )
      asr87(jdu,1,2,1) = asr87(jdu,1,2,1) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p4,p5)*zb(p4,p6)*zba2(p4,p1,p8
     &    ,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8) - 2.D0*za(p5,p6)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(p8,
     &    p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      asr87(jdu,1,1,2)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p7)**2*zb(p4,p6
     &    )*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p2,p5)*za(p3,p6)*zb(p1,p7)**2*zb(p4,p6)*
     &    zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s346**(-1)
     &     )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p4)*zb(p1,
     &    p4)*zb(p4,p6)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*s346**(-1) - 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p6)*zb(p4,
     &    p6)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupnu(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p4)*zb(p1,p4)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p6)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p5,p6)*
     &    zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s356**(-1)
     &     )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*
     &    zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s356**(-1)
     &     )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*
     &    zb(p1,p7)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s356**(-1) - 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p4)
     &    *zb(p1,p7)*zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)
     &    **2*zb(p3,p4)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p4,
     &    p5)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,
     &    p6)*zb(p3,p4)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p6)*zb(p4,
     &    p5)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p3,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)
     &    *zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,
     &    p8)*izb(p7,p8) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)**2*zb(p4,
     &    p6)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)**2*za(p4,
     &    p5)*zb(p1,p4)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,
     &    p8)*izb(p7,p8) + 2.D0*za(p2,p8)**2*za(p5,p6)*zb(p1,p6)*zb(p4,
     &    p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      asr87(jdu,1,1,2) = asr87(jdu,1,1,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p4,p5)*zb(
     &    p1,p4)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p3)*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     )
      asr87(jdu,1,2,2)= + coupnu(jdu,1)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s346**(-1) - 2.D0*
     &    za(p2,p5)*za(p3,p4)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *zba3(p4,p1,p7,p8,p2)*s346**(-1) + 2.D0*za(p2,p5)*za(p3,p4)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s346**(-1) - 2.D0*za(p2,p5)*za(p3,p6)*
     &    zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2
     &    )*s346**(-1) - 2.D0*za(p2,p5)*za(p3,p6)*zb(p1,p8)*zb(p4,p6)*
     &    iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s346**(-1) + 2.D0*
     &    za(p2,p5)*za(p3,p6)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*s346**(-1) )
      asr87(jdu,1,2,2) = asr87(jdu,1,2,2) + coupnu(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(
     &    p3,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s356**(-1)
     &     + 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p2,p7)*
     &    iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s356**(-1) - 2.D0*za(p2,p3)*
     &    za(p3,p5)*zb(p3,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*
     &    iza(p1,p8)*zba3(p4,p1,p7,p8,p2)*s356**(-1) + 2.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p7)*zb(p5,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2)*s356**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p8)*
     &    zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*
     &    s356**(-1) - 2.D0*za(p2,p5)*za(p3,p5)*zb(p5,p6)*zba2(p7,p1,p8
     &    ,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)*
     &    s356**(-1) )
      asr87(jdu,1,2,2) = asr87(jdu,1,2,2) + coupel(jdu,1)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*
     &    zb(p3,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*
     &    s345**(-1) - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p4)*
     &    iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) + 2.D0*
     &    za(p2,p3)*za(p3,p5)*zb(p3,p4)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) + 2.D0*
     &    za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,p8)
     &    *zba3(p6,p1,p7,p8,p2)*s345**(-1) + 2.D0*za(p2,p5)*za(p3,p5)*
     &    zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2
     &    )*s345**(-1) - 2.D0*za(p2,p5)*za(p3,p5)*zb(p4,p5)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*
     &    s345**(-1) )
      asr87(jdu,1,2,2) = asr87(jdu,1,2,2) + coupel(jdu,1)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)
     &    *zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,
     &    p2) - 2.D0*za(p2,p3)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)
     &    *iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(
     &    p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) + 2.D0*za(
     &    p2,p3)*za(p5,p6)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p6,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p4,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1
     &    ,p7,p8,p2) )
      asr87(jdu,2,1,1)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p1,p7)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*s346**(-1) - 2.D0*za(
     &    p1,p7)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1
     &    ,p8)*izb(p7,p8)*s346**(-1) - 2.D0*za(p1,p8)*za(p3,p4)*zb(p2,
     &    p4)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p1,p8)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*s346**(-1) - 2.D0*za(
     &    p3,p4)*za(p7,p8)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8)*s346**(-1) - 2.D0*za(p3,p6)*za(p7,p8)*zb(p2,
     &    p6)*zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*
     &    s346**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p1)
     &    **2*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,
     &    p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p1)*zb(p4,p6)
     &    *zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*s346**(-1) + 2.D0*
     &    za(p1,p5)*za(p3,p4)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*
     &    izb(p1,p7)*izb(p7,p8)*s346**(-1) - 2.D0*za(p1,p5)*za(p3,p6)*
     &    zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p1)*
     &    zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*s346**(-1)
     &     + 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*s346**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p4)*za(p5,p8)*zb(p2,p1)*zb(p2,p4)*zb(
     &    p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*s346**(-1) + 2.D0*za(
     &    p3,p6)*za(p5,p8)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*s346**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(
     &    p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*s356**(-1) + 2.
     &    D0*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*
     &    izb(p1,p8)*izb(p7,p8)*s356**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*
     &    zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*s356**(-1) + 2.D0*za(
     &    p3,p5)*za(p7,p8)*zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8)*s356**(-1) + 2.D0*za(p3,p5)*za(p7,p8)*zb(p2,
     &    p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*
     &    s356**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)**2*
     &    zb(p3,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    s356**(-1) - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p6)*
     &    zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*s356**(-1) - 2.D0*za(
     &    p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p4)*izb(p1
     &    ,p7)*izb(p7,p8)*s356**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,
     &    p1)**2*zb(p5,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(
     &    p1,p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p5,
     &    p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*s356**(-1) - 2.D0
     &    *za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p5,p6)*zab2(p8,p2,p7,p4)*
     &    izb(p1,p7)*izb(p7,p8)*s356**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,p4)*
     &    zb(p3,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*s356**(-1) - 2.D0*
     &    za(p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p4)*zb(p5,p6)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*s356**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p7)*za(p3,p5)*zb(p2,p6)*
     &    zb(p3,p4)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*s345**(-1)
     &     + 2.D0*za(p1,p7)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)*zab2(p5,p7,p8
     &    ,p1)*izb(p1,p8)*izb(p7,p8)*s345**(-1) - 2.D0*za(p1,p8)*za(p3,
     &    p5)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)
     &    *zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*s345**(-1) - 2.D0*
     &    za(p3,p5)*za(p7,p8)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p7,p8,p1)*
     &    izb(p1,p7)*izb(p1,p8)*s345**(-1) + 2.D0*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p5)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*
     &    s345**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)
     &    **2*zb(p3,p4)*zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,
     &    p8)*s345**(-1) + 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p4)
     &    *zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*s345**(-1) + 2.D0*
     &    za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p3,p4)*zab2(p8,p2,p7,p6)*
     &    izb(p1,p7)*izb(p7,p8)*s345**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*
     &    zb(p2,p1)**2*zb(p4,p5)*zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*s345**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p1)*
     &    zb(p4,p5)*zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*s345**(-1)
     &     - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p4,p5)*zab2(p8,p2,p7
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*s345**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,p6)*zb(
     &    p3,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*s345**(-1) - 2.D0*za(
     &    p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p6)*zb(p4,p5)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*s345**(-1) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p1,p7)*za(p4,p5)
     &    *zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)
     &     + 2.D0*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8
     &    ,p1)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(p1,p8)*za(p4,p5)*zb(p2,
     &    p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,p8) + 2.D0*
     &    za(p1,p8)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*
     &    izb(p1,p7)*izb(p7,p8) - 2.D0*za(p4,p5)*za(p7,p8)*zb(p2,p4)*
     &    zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p1,p8) + 2.D0*za(
     &    p5,p6)*za(p7,p8)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p1,p3)*za(p4,p5)
     &    *zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,
     &    p7)*izb(p1,p8) + 2.D0*za(p1,p3)*za(p4,p5)*zb(p2,p1)*zb(p4,p6)
     &    *zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8) + 2.D0*za(p1,p3)*za(
     &    p4,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(
     &    p7,p8) + 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p1)**2*zb(p4,p6)*
     &    zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8) - 2.D0*za(
     &    p1,p3)*za(p5,p6)*zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1
     &    ,p8)*izb(p7,p8) - 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p1)*zb(p4,p6
     &    )*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p7,p8) )
      asr87(jdu,2,1,1) = asr87(jdu,2,1,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p3,p8)*za(p4,p5)*zb(p2,p1)*zb(
     &    p2,p4)*zb(p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8) - 2.D0*za(
     &    p3,p8)*za(p5,p6)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8) )
      asr87(jdu,2,1,2)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p1,p8)*za(p3,p4)*zb(p2,p4)*zb(p1,p7)*
     &    zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p1,p8)*za(p3,p6)*zb(p2,p6)*zb(p1,p7)*zb(
     &    p4,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p4)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s346**(-1) - 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p6)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s346**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p4)*zb(p4,p6)*zab2(p5,p1,p8,p7)*zab2(
     &    p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p3,p6)*zb(p4,p6)*zab2(p5,p1,p8,p7)*zab2(
     &    p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(
     &    p1,p7)*zb(p3,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p1,
     &    p7)*zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p2,p8)*za(p1,p3)*za(p3,p5)*zb(
     &    p2,p7)*zb(p3,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(
     &    p7,p8)*s356**(-1) + 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p5)*zb(p2,
     &    p7)*zb(p5,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8)*s356**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p3,p5)*zb(p3,p6)*zab2(p3,p1,p8,p7)*
     &    zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *s356**(-1) - 2.D0*za(p3,p5)*zb(p5,p6)*zab2(p5,p1,p8,p7)*
     &    zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *s356**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p3,p4)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p1,p8)*za(p3,p5)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p4,p5)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8)*s345**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)*za(p1,p3)*za(p3,p5)*
     &    zb(p2,p7)*zb(p3,p4)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) + 2.D0*za(p2,p8)*za(p1,p5)*za(p3,p5)*
     &    zb(p2,p7)*zb(p4,p5)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*s345**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p3,p5)*zb(p3,p4)*zab2(p3,p1,p8,p7)*zab2(
     &    p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p3,p5)*zb(p4,p5)*zab2(p5,p1,p8,p7)*zab2(
     &    p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p1,p8)*za(p4,p5)
     &    *zb(p2,p4)*zb(p1,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*
     &    izb(p1,p8)*izb(p7,p8) + 2.D0*za(p1,p8)*za(p5,p6)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)*za(p1,p3)
     &    *za(p4,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p8)*za(p1,p3)*za(p5,p6)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      asr87(jdu,2,1,2) = asr87(jdu,2,1,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p4,p5)*zb(p4,p6)*zab2(p3,p1,p8
     &    ,p7)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(
     &    p7,p8) - 2.D0*za(p5,p6)*zb(p4,p6)*zab2(p3,p1,p8,p7)*zab2(p8,
     &    p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      asr87(jdu,2,2,1)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p1,p7)**2*za(p3,p4)*zb(p2,p4)*zb(p4,p6
     &    )*zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s346**(-1) - 2.D0*za(p1,p7)**2*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s346**(-1)
     &     )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p8)
     &    **2*zb(p4,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s346**(-1) - 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p8)**2*zb(p4,
     &    p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p4)*zb(p2,p4)*zb(
     &    p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) + 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p6)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s346**(-1) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*
     &    zb(p3,p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s356**(-1)
     &     )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)**2*
     &    zb(p3,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s356**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p8)**2*zb(p5,p6)*
     &    zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*s356**(-1)
     &     )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*
     &    zb(p2,p8)*zb(p3,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p4)
     &    *zb(p2,p8)*zb(p5,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s356**(-1) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,
     &    p6)*zb(p3,p4)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p1,p7)**2*za(p3,p5)*zb(p2,p6)*zb(p4,
     &    p5)*zab2(p5,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)
     &    **2*zb(p3,p4)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s345**(-1) + 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p8)**2*zb(p4,
     &    p5)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * ( 2.D0*za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p3,p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) - 2.D0*za(p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p4,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p1,p7)**2*za(p4,
     &    p5)*zb(p2,p4)*zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8) + 2.D0*za(p1,p7)**2*za(p5,p6)*zb(p2,p6)*zb(p4,
     &    p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * (  - 2.D0*za(p1,p3)*za(p4,p5)
     &    *zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8) + 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p8)**2*zb(p4,
     &    p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      asr87(jdu,2,2,1) = asr87(jdu,2,2,1) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * ( 2.D0*za(p1,p3)*za(p1,p7)*za(p4,p5)*zb(
     &    p2,p4)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) - 2.D0*za(p1,p3)*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(
     &    p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     )
      asr87(jdu,2,2,2)= + coupnu(jdu,2)*propw34**(-1)*cxw**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p1)**2*za(p3,p4)*zb(p2,p4)*zb(p4,p6
     &    )*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    s346**(-1) - 2.D0*za(p2,p1)**2*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s346**(-1)
     &     + 2.D0*za(p2,p1)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p1,p7
     &    ,p8)*iza(p2,p7)*iza(p7,p8)*s346**(-1) + 2.D0*za(p2,p1)*za(p3,
     &    p4)*zb(p2,p4)*zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,
     &    p8)*s346**(-1) + 2.D0*za(p2,p1)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)
     &    *zab2(p5,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*s346**(-1) + 2.D0*
     &    za(p2,p1)*za(p3,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p5,p1,p8,p7)*
     &    iza(p2,p8)*iza(p7,p8)*s346**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p5)*za(p3,p4)*zb(p2,p7)*zb(
     &    p4,p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*s346**(-1) + 2.
     &    D0*za(p1,p5)*za(p3,p4)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p4)*
     &    iza(p2,p7)*iza(p7,p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p4)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8)*
     &    s346**(-1) + 2.D0*za(p1,p5)*za(p3,p6)*zb(p2,p7)*zb(p4,p6)*
     &    zab2(p2,p7,p8,p6)*iza(p2,p8)*iza(p7,p8)*s346**(-1) + 2.D0*za(
     &    p1,p5)*za(p3,p6)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p6)*iza(p2
     &    ,p7)*iza(p7,p8)*s346**(-1) + 2.D0*za(p1,p5)*za(p3,p6)*zb(p4,
     &    p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p2,p8)*
     &    s346**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupnu(jdu,2)*propw34**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p2,p1)*za(p1,p5)*za(p3,p4)*zb(p4,p6)*
     &    zb(p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s346**(-1) - 2.D0*
     &    za(p2,p1)*za(p1,p5)*za(p3,p6)*zb(p4,p6)*zb(p6,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*s346**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1) * ( 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*
     &    zb(p3,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    s356**(-1) + 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s356**(-1)
     &     - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(p3,p6)*zab2(p3,p1,p7
     &    ,p8)*iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*za(p2,p1)*za(p3,
     &    p5)*zb(p2,p4)*zb(p3,p6)*zab2(p3,p1,p8,p7)*iza(p2,p8)*iza(p7,
     &    p8)*s356**(-1) - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)
     &    *zab2(p5,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*
     &    za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(p5,p6)*zab2(p5,p1,p8,p7)*
     &    iza(p2,p8)*iza(p7,p8)*s356**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1) * (  - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*s356**(-1)
     &     - 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(p3,p6)*zab2(p2,p7,p8
     &    ,p4)*iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*za(p1,p3)*za(p3,
     &    p5)*zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,
     &    p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*zb(p5,p6)
     &    *zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*s356**(-1) - 2.D0*
     &    za(p1,p5)*za(p3,p5)*zb(p2,p8)*zb(p5,p6)*zab2(p2,p7,p8,p4)*
     &    iza(p2,p7)*iza(p7,p8)*s356**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8)*
     &    s356**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupnu(jdu,2)*propw56**(-1)
     & *cxw**(-1) * ( 2.D0*za(p2,p1)*za(p1,p3)*za(p3,p5)*zb(p3,p6)*zb(
     &    p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s356**(-1) + 2.D0*za(
     &    p2,p1)*za(p1,p5)*za(p3,p5)*zb(p4,p7)*zb(p5,p6)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*s356**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s178**(-1) * (  - 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,
     &    p6)*zb(p3,p4)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,
     &    p8)*s345**(-1) + 2.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p6)*zb(p4,
     &    p5)*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    s345**(-1) + 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p3,p4)*
     &    zab2(p3,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*s345**(-1) + 2.D0*za(
     &    p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p3,p4)*zab2(p3,p1,p8,p7)*iza(p2
     &    ,p8)*iza(p7,p8)*s345**(-1) - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,
     &    p6)*zb(p4,p5)*zab2(p5,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*
     &    s345**(-1) - 2.D0*za(p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*s345**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1)*s278**(-1) * ( 2.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*zb(
     &    p3,p4)*zab2(p2,p7,p8,p6)*iza(p2,p8)*iza(p7,p8)*s345**(-1) + 2.
     &    D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(p3,p4)*zab2(p2,p7,p8,p6)*
     &    iza(p2,p7)*iza(p7,p8)*s345**(-1) + 2.D0*za(p1,p3)*za(p3,p5)*
     &    zb(p3,p4)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p2,p8)*
     &    s345**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*zb(p4,p5)*
     &    zab2(p2,p7,p8,p6)*iza(p2,p8)*iza(p7,p8)*s345**(-1) - 2.D0*za(
     &    p1,p5)*za(p3,p5)*zb(p2,p8)*zb(p4,p5)*zab2(p2,p7,p8,p6)*iza(p2
     &    ,p7)*iza(p7,p8)*s345**(-1) - 2.D0*za(p1,p5)*za(p3,p5)*zb(p4,
     &    p5)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p2,p8)*
     &    s345**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupel(jdu,2)*propw34**(-1)
     & *cxw**(-1) * (  - 2.D0*za(p2,p1)*za(p1,p3)*za(p3,p5)*zb(p3,p4)*
     &    zb(p6,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*s345**(-1) + 2.D0*
     &    za(p2,p1)*za(p1,p5)*za(p3,p5)*zb(p4,p5)*zb(p6,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*s345**(-1) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p1)**2*za(p4,
     &    p5)*zb(p2,p4)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,
     &    p8)*iza(p1,p8) + 2.D0*za(p2,p1)**2*za(p5,p6)*zb(p2,p6)*zb(p4,
     &    p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8) + 2.D0
     &    *za(p2,p1)*za(p4,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p3,p1,p7,p8)*
     &    iza(p2,p7)*iza(p7,p8) + 2.D0*za(p2,p1)*za(p4,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p8)*iza(p7,p8) - 2.D0*za(
     &    p2,p1)*za(p5,p6)*zb(p2,p6)*zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p2
     &    ,p7)*iza(p7,p8) - 2.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*zb(p4,p6
     &    )*zab2(p3,p1,p8,p7)*iza(p2,p8)*iza(p7,p8) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s278**(-1)*s456**(-1) * ( 2.D0*za(p1,p3)*za(p4,p5)*
     &    zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)
     &     + 2.D0*za(p1,p3)*za(p4,p5)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8
     &    ,p4)*iza(p2,p7)*iza(p7,p8) + 2.D0*za(p1,p3)*za(p4,p5)*zb(p4,
     &    p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8) - 2.D0*
     &    za(p1,p3)*za(p5,p6)*zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,p6)*
     &    iza(p2,p8)*iza(p7,p8) - 2.D0*za(p1,p3)*za(p5,p6)*zb(p2,p8)*
     &    zb(p4,p6)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p7,p8) - 2.D0*za(
     &    p1,p3)*za(p5,p6)*zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2
     &    ,p7)*iza(p2,p8) )
      asr87(jdu,2,2,2) = asr87(jdu,2,2,2) + coupel(jdu,2)*propw56**(-1)
     & *cxw**(-1)*s456**(-1) * (  - 2.D0*za(p2,p1)*za(p1,p3)*za(p4,p5)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8) + 2.D0*
     &    za(p2,p1)*za(p1,p3)*za(p5,p6)*zb(p4,p6)*zb(p6,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8) )
      endif
      enddo
      enddo

      do jdu=1,2
      if (jdu==2) then
      p3=i3
      p4=i4
      p5=i5
      p6=i6
      sign=cone
      elseif (jdu==1) then
      p3=i5
      p4=i6
      p5=i3
      p6=i4
      sign=-cone
      endif

      p7=i7
      p8=i8
      s34=s(p3,p4)
      s56=s(p5,p6)
      propw34=s34-cwmass2
      propw56=s56-cwmass2
      s278=t3(p2,p7,p8)
      s178=t3(p1,p7,p8)
      s156=t3(p1,p5,p6)
      s234=t3(p2,p3,p4)
      s356=t3(p3,p5,p6)
      s345=t3(p3,p4,p5)
      s346=t3(p3,p4,p6)
      s456=t3(p4,p5,p6)

      do j78=1,2
      if (j78 == 1) then
      p7=i7
      p8=i8
      s2347=t4(p2,p3,p4,p7)
      adr78(jdu,1,1,1)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s278**(-1)*s156**(-1) * (  - za(p2,p7)*zb(p1,p6)*zba2(p4,p1,p6,
     &    p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - za(p2,p8)*
     &    zb(p1,p6)*zba2(p4,p1,p6,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,
     &    p7,p8,p3) - zb(p1,p6)*zba2(p1,p2,p7,p8)*zba2(p4,p1,p6,p5)*
     &    izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) )
      adr78(jdu,1,1,1) = adr78(jdu,1,1,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * (  - za(p2,p3)*za(p5,p6)*zb(
     &    p1,p6)**2*zba2(p4,p2,p3,p7)*izb(p1,p7)*izb(p1,p8)*zba4(p1,p2,
     &    p3,p4,p7,p8)*s2347**(-1) + za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*
     &    zba2(p4,p2,p3,p7)*izb(p1,p8)*izb(p7,p8) + za(p2,p3)*za(p5,p6)
     &    *zb(p1,p6)**2*zba2(p4,p2,p3,p8)*izb(p1,p7)*izb(p7,p8) )
      adr78(jdu,1,1,1) = adr78(jdu,1,1,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1) * ( za(p5,p6)*zb(p1,p6)**2*zba2(p1,p2,p7,
     &    p3)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*
     &    s2347**(-1) )
      adr78(jdu,1,1,1) = adr78(jdu,1,1,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * ( za(p2,p7)*za(p3,p5)*zb(p1
     &    ,p3)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*
     &    sign + za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p4)*sign - za(p2,p7)*za(p3,p5)*
     &    zb(p1,p5)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*sign - za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p6)*sign + 2.D0*za(p2,p7)*zb(p1,
     &    p4)*zab2(p5,p3,p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,
     &    p3)*sign - 2.D0*za(p2,p7)*zb(p1,p6)*zab2(p3,p5,p6,p4)*izb(p1,
     &    p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*sign + za(p2,p8)*za(p3,p5
     &    )*zb(p1,p3)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,
     &    p3)*sign + za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)
     &    *izb(p7,p8)*zba3(p1,p2,p7,p8,p4)*sign - za(p2,p8)*za(p3,p5)*
     &    zb(p1,p5)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*sign )
      adr78(jdu,1,1,1) = adr78(jdu,1,1,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p8)*za(p3,p5)*
     &    zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p6
     &    )*sign + 2.D0*za(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*izb(p1,p7
     &    )*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*sign - 2.D0*za(p2,p8)*zb(p1
     &    ,p6)*zab2(p3,p5,p6,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8
     &    ,p5)*sign + za(p3,p5)*zb(p1,p3)*zb(p4,p6)*zba2(p1,p2,p7,p8)*
     &    izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*sign + 
     &    za(p3,p5)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*sign - za(p3,p5)*
     &    zb(p1,p5)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*
     &    izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*sign - za(p3,p5)*zb(p1,p6)*
     &    zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p6)*sign + 2.D0*zb(p1,p4)*zab2(p5,p3,p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3)*sign )
      adr78(jdu,1,1,1) = adr78(jdu,1,1,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - 2.D0*zb(p1,p6)*zab2(p3
     &    ,p5,p6,p4)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)
     &    *zba3(p1,p2,p7,p8,p5)*sign )
      adr78(jdu,1,1,2)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s178**(-1)*s234**(-1) * (  - za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(
     &    p4,p2,p3,p5)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8) )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s278**(-1)*s156**(-1) * ( za(p2,p7)*zb(p2,p8)*zb(p1,
     &    p6)*zba2(p4,p1,p6,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)
     &    *izb(p7,p8) )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * ( za(p2,p3)*zb(p1,p6)*zba2(p4
     &    ,p2,p3,p7)*zba2(p8,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2
     &    ,p3,p4,p7)*s2347**(-1) )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1) * (  - zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(p8
     &    ,p2,p7,p3)*zba2(p8,p1,p6,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *s2347**(-1) )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s234**(-1) * (  - za(p2,p3)*zba2(p4,p2,p3,p7)*zba2(p6
     &    ,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba4(p8,p2,p3,p4,
     &    p7,p5)*s2347**(-1) )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2) * ( zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*
     &    s2347**(-1) )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - za(p2,p3)*za(p1,p7)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p3,p1,p8,p7)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8)*sign - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,
     &    p8)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)
     &    *izb(p7,p8)*sign - za(p2,p4)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*
     &    zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    sign + za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zba2(p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + 2.D0
     &    *za(p2,p5)*za(p1,p7)*zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p6,p1,
     &    p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + za(p2,p6)*za(
     &    p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p6,p1,p8,p7)*iza(p1
     &    ,p8)*iza(p7,p8)*izb(p7,p8)*sign )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p7)*za(p3,p5)*
     &    zb(p2,p8)*zb(p1,p3)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*sign - za(p2,p7)*za(p3,p5)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*sign + za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p5)*zb(
     &    p4,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    sign + za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p6)*zb(p4,p6)*
     &    zba2(p8,p2,p7,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - 2.D0
     &    *za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p8,p2,
     &    p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign + 2.D0*za(p2,p7)
     &    *zb(p2,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p5)*iza(
     &    p7,p8)*izb(p2,p7)*izb(p7,p8)*sign )
      adr78(jdu,1,1,2) = adr78(jdu,1,1,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p3,p5)*zb(p4,p6)*zba2(p3,p1,p8,p7
     &    )*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*sign + za(p3,p5)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,
     &    p7,p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(
     &    p3,p5)*zb(p4,p6)*zba2(p5,p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p1,
     &    p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(p3,p5)*zb(p4,
     &    p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p6)*iza(p1,p8)*iza(p7,p8)
     &    *izb(p2,p7)*izb(p7,p8)*sign - 2.D0*zab2(p3,p5,p6,p4)*zba2(p6,
     &    p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*sign + 2.D0*zab2(p5,p3,p4,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *sign )
      adr78(jdu,1,2,1)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s178**(-1)*s234**(-1) * (  - za(p2,p3)*zb(p1,p7)**2*zba2(p4,p2,
     &    p3,p5)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s278**(-1)*s156**(-1) * ( za(p2,p8)**2*zb(p1,p6)*
     &    zba2(p4,p1,p6,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * ( za(p2,p3)*zb(p1,p6)*zba2(p4
     &    ,p2,p3,p8)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p7,p2
     &    ,p3,p4,p8)*s2347**(-1) )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1) * ( za(p2,p3)*za(p2,p8)*zb(p1,p6)*zba2(p7,
     &    p1,p6,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,
     &    p8)*s2347**(-1) )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s234**(-1) * ( za(p2,p3)*zb(p1,p6)*zb(p1,p7)*zba2(p4,
     &    p2,p3,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p2,p3,p4,
     &    p5)*s2347**(-1) )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2) * ( za(p2,p3)*za(p2,p8)*zb(p1,p6)*zb(p1,p7)*iza(p2,p7
     &    )*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*
     &    s2347**(-1) )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - za(p2,p3)*za(p3,p5)*
     &    zb(p1,p7)**2*zb(p4,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p1,p8
     &    )*izb(p7,p8)*sign - 2.D0*za(p2,p3)*zb(p1,p7)**2*zab2(p5,p3,p4
     &    ,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign
     &     - za(p2,p4)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,
     &    p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign + za(p2,p5)*za(p3,
     &    p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p5,p1,p7,p8)*iza(p7,p8)*izb(
     &    p1,p8)*izb(p7,p8)*sign + 2.D0*za(p2,p5)*zb(p1,p7)**2*zab2(p3,
     &    p5,p6,p4)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    sign + za(p2,p6)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p6,p1,
     &    p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p8)**2*za(p3,p5)
     &    *zb(p1,p3)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*sign - za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*
     &    zba2(p7,p2,p8,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + za(
     &    p2,p8)**2*za(p3,p5)*zb(p1,p5)*zb(p4,p6)*zba2(p7,p2,p8,p5)*
     &    iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + za(p2,p8)**2*za(p3,p5
     &    )*zb(p1,p6)*zb(p4,p6)*zba2(p7,p2,p8,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p7,p8)*sign - 2.D0*za(p2,p8)**2*zb(p1,p4)*zab2(p5,p3,p4,
     &    p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign
     &     + 2.D0*za(p2,p8)**2*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p7,p2,
     &    p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign )
      adr78(jdu,1,2,1) = adr78(jdu,1,2,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,
     &    p3)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(
     &    p7,p8)*sign + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*
     &    zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *sign + za(p2,p4)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(
     &    p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(
     &    p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p5)*zb(p1,p7)*zb(p4,p6)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - 2.D0*za(p2,p5)
     &    *za(p2,p8)*zb(p1,p6)*zb(p1,p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(p2,p6)*za(p2,p8)*
     &    za(p3,p5)*zb(p1,p6)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p1,p8)*izb(p7,p8)*sign )
      adr78(jdu,1,2,2)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s178**(-1)*s234**(-1) * (  - za(p2,p3)*zb(p1,p7)*zba2(p4,p2,p3,
     &    p5)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) - za(p2,p3)*
     &    zb(p1,p8)*zba2(p4,p2,p3,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,
     &    p7,p8,p2) + za(p2,p3)*zba2(p4,p2,p3,p5)*zba2(p7,p1,p8,p2)*
     &    iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2) )
      adr78(jdu,1,2,2) = adr78(jdu,1,2,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * (  - za(p2,p3)**2*zb(p1,p6)*
     &    zb(p3,p4)*zba2(p7,p1,p6,p5)*iza(p2,p8)*iza(p7,p8) + za(p2,p3)
     &    **2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*iza(p2,p8)*s2347**(-1) - za(p2,p3)**2*zb(p1,p6)*
     &    zb(p3,p4)*zba2(p8,p1,p6,p5)*iza(p2,p7)*iza(p7,p8) )
      adr78(jdu,1,2,2) = adr78(jdu,1,2,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s234**(-1) * (  - za(p2,p3)**2*zb(p3,p4)*zba2(p6,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p7,p2,p3,p4,p5)*
     &    s2347**(-1) )
      adr78(jdu,1,2,2) = adr78(jdu,1,2,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - za(p2,p3)*za(p3,p5)*
     &    zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p3,p1,p7,p8,p2
     &    )*sign - za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*zba3(p3,p1,p7,p8,p2)*sign + za(p2,p3)*za(p3,p5)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p3,p1,p7,p8,p2)*sign - 2.D0*za(p2,p3)*zb(p1,p7)*zab2(p5,
     &    p3,p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*sign - 2.
     &    D0*za(p2,p3)*zb(p1,p8)*zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8
     &    )*zba3(p4,p1,p7,p8,p2)*sign + 2.D0*za(p2,p3)*zab2(p5,p3,p4,p6
     &    )*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,
     &    p1,p7,p8,p2)*sign - za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*sign - za(p2,p4)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2)*sign + za(p2,p4)*za(p3,p5)*zb(p4,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)*
     &    sign )
      adr78(jdu,1,2,2) = adr78(jdu,1,2,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * ( za(p2,p5)*za(p3,p5)*zb(p1
     &    ,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p5,p1,p7,p8,p2)*
     &    sign + za(p2,p5)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*zba3(p5,p1,p7,p8,p2)*sign - za(p2,p5)*za(p3,p5)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p5,p1,p7,p8,p2)*sign + 2.D0*za(p2,p5)*zb(p1,p7)*zab2(p3,
     &    p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*sign + 2.
     &    D0*za(p2,p5)*zb(p1,p8)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p7,p8
     &    )*zba3(p6,p1,p7,p8,p2)*sign - 2.D0*za(p2,p5)*zab2(p3,p5,p6,p4
     &    )*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,
     &    p1,p7,p8,p2)*sign + za(p2,p6)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*sign + za(p2,p6)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p6,
     &    p1,p7,p8,p2)*sign - za(p2,p6)*za(p3,p5)*zb(p4,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*
     &    sign )
      adr78(jdu,2,1,1)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p1,p7)*za(p3,p5)*zb(p2,p3)*zb(p4,
     &    p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*sign - za(p1,p7)*
     &    za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*
     &    izb(p7,p8)*sign + za(p1,p7)*za(p3,p5)*zb(p2,p5)*zb(p4,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*sign + za(p1,p7)*za(
     &    p3,p5)*zb(p2,p6)*zb(p4,p6)*zab2(p6,p7,p8,p1)*izb(p1,p8)*izb(
     &    p7,p8)*sign - 2.D0*za(p1,p7)*zb(p2,p4)*zab2(p3,p7,p8,p1)*
     &    zab2(p5,p3,p4,p6)*izb(p1,p8)*izb(p7,p8)*sign + 2.D0*za(p1,p7)
     &    *zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p7,p8,p1)*izb(p1,p8)*
     &    izb(p7,p8)*sign - za(p1,p8)*za(p3,p5)*zb(p2,p3)*zb(p4,p6)*
     &    zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*sign - za(p1,p8)*za(
     &    p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*izb(
     &    p7,p8)*sign + za(p1,p8)*za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5
     &    ,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*sign + za(p1,p8)*za(p3,p5)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p6,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*
     &    sign )
      adr78(jdu,2,1,1) = adr78(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p8)*zb(p2,
     &    p4)*zab2(p3,p7,p8,p1)*zab2(p5,p3,p4,p6)*izb(p1,p7)*izb(p7,p8)
     &    *sign + 2.D0*za(p1,p8)*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p7
     &    ,p8,p1)*izb(p1,p7)*izb(p7,p8)*sign - za(p3,p5)*za(p7,p8)*zb(
     &    p2,p3)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*sign
     &     - za(p3,p5)*za(p7,p8)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p7,p8,p1)*
     &    izb(p1,p7)*izb(p1,p8)*sign + za(p3,p5)*za(p7,p8)*zb(p2,p5)*
     &    zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*sign + za(
     &    p3,p5)*za(p7,p8)*zb(p2,p6)*zb(p4,p6)*zab2(p6,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8)*sign - 2.D0*za(p7,p8)*zb(p2,p4)*zab2(p3,p7,p8
     &    ,p1)*zab2(p5,p3,p4,p6)*izb(p1,p7)*izb(p1,p8)*sign + 2.D0*za(
     &    p7,p8)*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p7,p8,p1)*izb(p1,
     &    p7)*izb(p1,p8)*sign )
      adr78(jdu,2,1,1) = adr78(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p1,p3)*za(p3,p5)*
     &    zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*sign + za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*
     &    zab2(p7,p2,p8,p3)*izb(p1,p8)*izb(p7,p8)*sign + za(p1,p3)*za(
     &    p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*izb(
     &    p7,p8)*sign - 2.D0*za(p1,p3)*zb(p2,p1)**2*zab2(p5,p3,p4,p6)*
     &    zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign + 2.D0
     &    *za(p1,p3)*zb(p2,p1)*zab2(p5,p3,p4,p6)*zab2(p7,p2,p8,p4)*izb(
     &    p1,p8)*izb(p7,p8)*sign + 2.D0*za(p1,p3)*zb(p2,p1)*zab2(p5,p3,
     &    p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p7,p8)*sign - za(p1,
     &    p4)*za(p3,p5)*zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p2
     &    ,p7)*izb(p1,p7)*izb(p1,p8)*sign + za(p1,p4)*za(p3,p5)*zb(p2,
     &    p1)*zb(p4,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*sign + 
     &    za(p1,p4)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*
     &    izb(p1,p7)*izb(p7,p8)*sign )
      adr78(jdu,2,1,1) = adr78(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * ( za(p1,p5)*za(p3,p5)*zb(p2
     &    ,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p5)*izb(p2,p7)*izb(p1,p7)*
     &    izb(p1,p8)*sign - za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*
     &    zab2(p7,p2,p8,p5)*izb(p1,p8)*izb(p7,p8)*sign - za(p1,p5)*za(
     &    p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p5)*izb(p1,p7)*izb(
     &    p7,p8)*sign + 2.D0*za(p1,p5)*zb(p2,p1)**2*zab2(p3,p5,p6,p4)*
     &    zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign - 2.D0
     &    *za(p1,p5)*zb(p2,p1)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p6)*izb(
     &    p1,p8)*izb(p7,p8)*sign - 2.D0*za(p1,p5)*zb(p2,p1)*zab2(p3,p5,
     &    p6,p4)*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p7,p8)*sign + za(p1,
     &    p6)*za(p3,p5)*zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p2
     &    ,p7)*izb(p1,p7)*izb(p1,p8)*sign - za(p1,p6)*za(p3,p5)*zb(p2,
     &    p1)*zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*sign - 
     &    za(p1,p6)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p6)*
     &    izb(p1,p7)*izb(p7,p8)*sign )
      adr78(jdu,2,1,1) = adr78(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,
     &    p3)*zb(p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign + za(p3,
     &    p5)*za(p4,p8)*zb(p2,p1)*zb(p2,p4)*zb(p4,p6)*izb(p2,p7)*izb(p1
     &    ,p7)*izb(p1,p8)*sign - za(p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p5
     &    )*zb(p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign - za(p3,p5)
     &    *za(p6,p8)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*sign + 2.D0*za(p3,p8)*zb(p2,p1)*zb(p2,p4)*zab2(
     &    p5,p3,p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign - 2.D0*za(
     &    p5,p8)*zb(p2,p1)*zb(p2,p6)*zab2(p3,p5,p6,p4)*izb(p2,p7)*izb(
     &    p1,p7)*izb(p1,p8)*sign )
      adr78(jdu,2,2,1)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p1,p8)*za(p3,p5)*zb(p2,p3)*zb(p1,
     &    p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*sign - za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p1,p7)*zb(p4,p6)*
     &    zab2(p4,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign + za(
     &    p1,p8)*za(p3,p5)*zb(p2,p5)*zb(p1,p7)*zb(p4,p6)*zab2(p5,p1,p8,
     &    p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign + za(p1,p8)*za(p3,
     &    p5)*zb(p2,p6)*zb(p1,p7)*zb(p4,p6)*zab2(p6,p1,p8,p7)*iza(p7,p8
     &    )*izb(p1,p8)*izb(p7,p8)*sign - 2.D0*za(p1,p8)*zb(p2,p4)*zb(p1
     &    ,p7)*zab2(p3,p1,p8,p7)*zab2(p5,p3,p4,p6)*iza(p7,p8)*izb(p1,p8
     &    )*izb(p7,p8)*sign + 2.D0*za(p1,p8)*zb(p2,p6)*zb(p1,p7)*zab2(
     &    p3,p5,p6,p4)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*sign )
      adr78(jdu,2,2,1) = adr78(jdu,2,2,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p8)*za(p1,p3)*
     &    za(p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p3)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p7,p8)*sign - 2.D0*za(p2,p8)*za(p1,p3)*zb(p2,
     &    p7)*zab2(p5,p3,p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p7,p8)*sign - za(p2,p8)*za(p1,p4)*za(p3,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    sign + za(p2,p8)*za(p1,p5)*za(p3,p5)*zb(p2,p7)*zb(p4,p6)*
     &    zab2(p8,p2,p7,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + 2.D0
     &    *za(p2,p8)*za(p1,p5)*zb(p2,p7)*zab2(p3,p5,p6,p4)*zab2(p8,p2,
     &    p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + za(p2,p8)*za(
     &    p1,p6)*za(p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p2
     &    ,p7)*iza(p7,p8)*izb(p7,p8)*sign )
      adr78(jdu,2,2,1) = adr78(jdu,2,2,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p3,p5)*zb(p4,p6)*zab2(p3,p1,p8,p7
     &    )*zab2(p8,p2,p7,p3)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*sign + za(p3,p5)*zb(p4,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,
     &    p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(
     &    p3,p5)*zb(p4,p6)*zab2(p5,p1,p8,p7)*zab2(p8,p2,p7,p5)*iza(p2,
     &    p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(p3,p5)*zb(p4,
     &    p6)*zab2(p6,p1,p8,p7)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p1,p8)*izb(p7,p8)*sign + 2.D0*zab2(p3,p1,p8,p7)*zab2(p5,
     &    p3,p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8)*sign - 2.D0*zab2(p3,p5,p6,p4)*zab2(p5,p1,p8,p7)*
     &    zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *sign )
      adr78(jdu,2,1,2)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p1,p7)**2*za(p3,p5)*zb(p2,p3)*zb(
     &    p4,p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    sign - za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p1,
     &    p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + za(p1,p7)**2*
     &    za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5,p1,p7,p8)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8)*sign + za(p1,p7)**2*za(p3,p5)*zb(p2,p6)
     &    *zb(p4,p6)*zab2(p6,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &    *sign - 2.D0*za(p1,p7)**2*zb(p2,p4)*zab2(p3,p1,p7,p8)*zab2(p5
     &    ,p3,p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + 2.D0*za(p1
     &    ,p7)**2*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p7,p8)*iza(p1,
     &    p8)*iza(p7,p8)*izb(p7,p8)*sign )
      adr78(jdu,2,1,2) = adr78(jdu,2,1,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p1,p3)*za(p3,p5)*
     &    zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p2,p7
     &    )*izb(p7,p8)*sign - 2.D0*za(p1,p3)*zb(p2,p8)**2*zab2(p5,p3,p4
     &    ,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign
     &     - za(p1,p4)*za(p3,p5)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,
     &    p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign + za(p1,p5)*za(p3,
     &    p5)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p5)*iza(p7,p8)*izb(
     &    p2,p7)*izb(p7,p8)*sign + 2.D0*za(p1,p5)*zb(p2,p8)**2*zab2(p3,
     &    p5,p6,p4)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    sign + za(p1,p6)*za(p3,p5)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,
     &    p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign )
      adr78(jdu,2,1,2) = adr78(jdu,2,1,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,
     &    p3)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8)*sign + 2.D0*za(p1,p3)*za(p1,p7)*zb(p2,p4)*zb(p2,p8)*
     &    zab2(p5,p3,p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *sign + za(p1,p4)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p2,p8)*zb(
     &    p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(
     &    p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p5)*zb(p2,p8)*zb(p4,p6)*iza(
     &    p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - 2.D0*za(p1,p5)
     &    *za(p1,p7)*zb(p2,p6)*zb(p2,p8)*zab2(p3,p5,p6,p4)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(p1,p6)*za(p1,p7)*
     &    za(p3,p5)*zb(p2,p6)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)
     &    *izb(p2,p7)*izb(p7,p8)*sign )
      adr78(jdu,2,2,2)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p2,p1)**2*za(p3,p5)*zb(p2,p3)*zb(
     &    p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    sign - za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p1,
     &    p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign + za(p2,p1)**2*
     &    za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*sign + za(p2,p1)**2*za(p3,p5)*zb(p2,p6)
     &    *zb(p4,p6)*zab2(p6,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)
     &    *sign - 2.D0*za(p2,p1)**2*zb(p2,p4)*zab2(p3,p1,p8,p7)*zab2(p5
     &    ,p3,p4,p6)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign + 2.D0*za(p2
     &    ,p1)**2*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p8,p7)*iza(p2,
     &    p7)*iza(p2,p8)*iza(p1,p8)*sign + za(p2,p1)*za(p3,p5)*zb(p2,p3
     &    )*zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*sign + 
     &    za(p2,p1)*za(p3,p5)*zb(p2,p3)*zb(p4,p6)*zab2(p3,p1,p8,p7)*
     &    iza(p2,p8)*iza(p7,p8)*sign + za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zab2(p4,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*sign )
      adr78(jdu,2,2,2) = adr78(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * ( za(p2,p1)*za(p3,p5)*zb(p2
     &    ,p4)*zb(p4,p6)*zab2(p4,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign
     &     - za(p2,p1)*za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5,p1,p7,p8)*
     &    iza(p2,p7)*iza(p7,p8)*sign - za(p2,p1)*za(p3,p5)*zb(p2,p5)*
     &    zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign - za(
     &    p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p4,p6)*zab2(p6,p1,p7,p8)*iza(p2
     &    ,p7)*iza(p7,p8)*sign - za(p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p4,p6
     &    )*zab2(p6,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign + 2.D0*za(p2,
     &    p1)*zb(p2,p4)*zab2(p3,p1,p7,p8)*zab2(p5,p3,p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*sign + 2.D0*za(p2,p1)*zb(p2,p4)*zab2(p3,p1,p8,p7)*
     &    zab2(p5,p3,p4,p6)*iza(p2,p8)*iza(p7,p8)*sign - 2.D0*za(p2,p1)
     &    *zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p7,p8)*iza(p2,p7)*
     &    iza(p7,p8)*sign - 2.D0*za(p2,p1)*zb(p2,p6)*zab2(p3,p5,p6,p4)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign )
      adr78(jdu,2,2,2) = adr78(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * ( za(p1,p3)*za(p3,p5)*zb(p2
     &    ,p7)*zb(p4,p6)*zab2(p2,p7,p8,p3)*iza(p2,p8)*iza(p7,p8)*sign
     &     + za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p3)*
     &    iza(p2,p7)*iza(p7,p8)*sign + za(p1,p3)*za(p3,p5)*zb(p4,p6)*
     &    zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*iza(p2,p8)*sign + 2.D0
     &    *za(p1,p3)*zb(p2,p7)*zab2(p2,p7,p8,p4)*zab2(p5,p3,p4,p6)*iza(
     &    p2,p8)*iza(p7,p8)*sign + 2.D0*za(p1,p3)*zb(p2,p8)*zab2(p2,p7,
     &    p8,p4)*zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*sign + 2.D0*
     &    za(p1,p3)*zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p5,p3,p4,p6)*iza(
     &    p2,p7)*iza(p2,p8)*sign + za(p1,p4)*za(p3,p5)*zb(p2,p7)*zb(p4,
     &    p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*sign + za(p1,p4)*
     &    za(p3,p5)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p4)*iza(p2,p7)*
     &    iza(p7,p8)*sign + za(p1,p4)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*
     &    zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8)*sign - za(p1,p5)*za(
     &    p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,p5)*iza(p2,p8)*iza(
     &    p7,p8)*sign )
      adr78(jdu,2,2,2) = adr78(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p1,p5)*za(p3,p5)*
     &    zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p5)*iza(p2,p7)*iza(p7,p8)*
     &    sign - za(p1,p5)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,
     &    p5)*iza(p2,p7)*iza(p2,p8)*sign - 2.D0*za(p1,p5)*zb(p2,p7)*
     &    zab2(p2,p7,p8,p6)*zab2(p3,p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*
     &    sign - 2.D0*za(p1,p5)*zb(p2,p8)*zab2(p2,p7,p8,p6)*zab2(p3,p5,
     &    p6,p4)*iza(p2,p7)*iza(p7,p8)*sign - 2.D0*za(p1,p5)*zb(p7,p8)*
     &    zab2(p2,p7,p8,p6)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p2,p8)*
     &    sign - za(p1,p6)*za(p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,
     &    p6)*iza(p2,p8)*iza(p7,p8)*sign - za(p1,p6)*za(p3,p5)*zb(p2,p8
     &    )*zb(p4,p6)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p7,p8)*sign - 
     &    za(p1,p6)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*
     &    iza(p2,p7)*iza(p2,p8)*sign )
      adr78(jdu,2,2,2) = adr78(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * (  - za(p2,p1)*za(p1,p3)*za(p3,p5)*zb(
     &    p3,p7)*zb(p4,p6)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign - 2.D0
     &    *za(p2,p1)*za(p1,p3)*zb(p4,p7)*zab2(p5,p3,p4,p6)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*sign - za(p2,p1)*za(p1,p4)*za(p3,p5)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign + 
     &    za(p2,p1)*za(p1,p5)*za(p3,p5)*zb(p4,p6)*zb(p5,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*sign + 2.D0*za(p2,p1)*za(p1,p5)*zb(p6,
     &    p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign
     &     + za(p2,p1)*za(p1,p6)*za(p3,p5)*zb(p4,p6)*zb(p6,p7)*iza(p2,
     &    p7)*iza(p2,p8)*iza(p1,p8)*sign )
      elseif (j78 == 2) then
      p7=i8
      p8=i7
      s2347=t4(p2,p3,p4,p7)
      adr87(jdu,1,1,1)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s278**(-1)*s156**(-1) * (  - za(p2,p7)*zb(p1,p6)*zba2(p4,p1,p6,
     &    p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - za(p2,p8)*
     &    zb(p1,p6)*zba2(p4,p1,p6,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,
     &    p7,p8,p3) - zb(p1,p6)*zba2(p1,p2,p7,p8)*zba2(p4,p1,p6,p5)*
     &    izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) )
      adr87(jdu,1,1,1) = adr87(jdu,1,1,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * (  - za(p2,p3)*za(p5,p6)*zb(
     &    p1,p6)**2*zba2(p4,p2,p3,p7)*izb(p1,p7)*izb(p1,p8)*zba4(p1,p2,
     &    p3,p4,p7,p8)*s2347**(-1) + za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*
     &    zba2(p4,p2,p3,p7)*izb(p1,p8)*izb(p7,p8) + za(p2,p3)*za(p5,p6)
     &    *zb(p1,p6)**2*zba2(p4,p2,p3,p8)*izb(p1,p7)*izb(p7,p8) )
      adr87(jdu,1,1,1) = adr87(jdu,1,1,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1) * ( za(p5,p6)*zb(p1,p6)**2*zba2(p1,p2,p7,
     &    p3)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*
     &    s2347**(-1) )
      adr87(jdu,1,1,1) = adr87(jdu,1,1,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * ( za(p2,p7)*za(p3,p5)*zb(p1
     &    ,p3)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*
     &    sign + za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p4)*sign - za(p2,p7)*za(p3,p5)*
     &    zb(p1,p5)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*sign - za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*
     &    izb(p7,p8)*zba3(p1,p2,p7,p8,p6)*sign + 2.D0*za(p2,p7)*zb(p1,
     &    p4)*zab2(p5,p3,p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,
     &    p3)*sign - 2.D0*za(p2,p7)*zb(p1,p6)*zab2(p3,p5,p6,p4)*izb(p1,
     &    p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*sign + za(p2,p8)*za(p3,p5
     &    )*zb(p1,p3)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,
     &    p3)*sign + za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)
     &    *izb(p7,p8)*zba3(p1,p2,p7,p8,p4)*sign - za(p2,p8)*za(p3,p5)*
     &    zb(p1,p5)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5
     &    )*sign )
      adr87(jdu,1,1,1) = adr87(jdu,1,1,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p8)*za(p3,p5)*
     &    zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p6
     &    )*sign + 2.D0*za(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*izb(p1,p7
     &    )*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*sign - 2.D0*za(p2,p8)*zb(p1
     &    ,p6)*zab2(p3,p5,p6,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8
     &    ,p5)*sign + za(p3,p5)*zb(p1,p3)*zb(p4,p6)*zba2(p1,p2,p7,p8)*
     &    izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*sign + 
     &    za(p3,p5)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*sign - za(p3,p5)*
     &    zb(p1,p5)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*
     &    izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*sign - za(p3,p5)*zb(p1,p6)*
     &    zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p6)*sign + 2.D0*zb(p1,p4)*zab2(p5,p3,p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3)*sign )
      adr87(jdu,1,1,1) = adr87(jdu,1,1,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - 2.D0*zb(p1,p6)*zab2(p3
     &    ,p5,p6,p4)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)
     &    *zba3(p1,p2,p7,p8,p5)*sign )
      adr87(jdu,1,2,1)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s178**(-1)*s234**(-1) * (  - za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(
     &    p4,p2,p3,p5)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8) )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s278**(-1)*s156**(-1) * ( za(p2,p7)*zb(p2,p8)*zb(p1,
     &    p6)*zba2(p4,p1,p6,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)
     &    *izb(p7,p8) )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * ( za(p2,p3)*zb(p1,p6)*zba2(p4
     &    ,p2,p3,p7)*zba2(p8,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2
     &    ,p3,p4,p7)*s2347**(-1) )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1) * (  - zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(p8
     &    ,p2,p7,p3)*zba2(p8,p1,p6,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *s2347**(-1) )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s234**(-1) * (  - za(p2,p3)*zba2(p4,p2,p3,p7)*zba2(p6
     &    ,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba4(p8,p2,p3,p4,
     &    p7,p5)*s2347**(-1) )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2) * ( zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*
     &    s2347**(-1) )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - za(p2,p3)*za(p1,p7)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p3,p1,p8,p7)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8)*sign - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,
     &    p8)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)
     &    *izb(p7,p8)*sign - za(p2,p4)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*
     &    zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    sign + za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zba2(p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + 2.D0
     &    *za(p2,p5)*za(p1,p7)*zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p6,p1,
     &    p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + za(p2,p6)*za(
     &    p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p6,p1,p8,p7)*iza(p1
     &    ,p8)*iza(p7,p8)*izb(p7,p8)*sign )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p7)*za(p3,p5)*
     &    zb(p2,p8)*zb(p1,p3)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*sign - za(p2,p7)*za(p3,p5)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*sign + za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p5)*zb(
     &    p4,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    sign + za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p6)*zb(p4,p6)*
     &    zba2(p8,p2,p7,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - 2.D0
     &    *za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p8,p2,
     &    p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign + 2.D0*za(p2,p7)
     &    *zb(p2,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p5)*iza(
     &    p7,p8)*izb(p2,p7)*izb(p7,p8)*sign )
      adr87(jdu,1,2,1) = adr87(jdu,1,2,1) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p3,p5)*zb(p4,p6)*zba2(p3,p1,p8,p7
     &    )*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*sign + za(p3,p5)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,
     &    p7,p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(
     &    p3,p5)*zb(p4,p6)*zba2(p5,p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p1,
     &    p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(p3,p5)*zb(p4,
     &    p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p6)*iza(p1,p8)*iza(p7,p8)
     &    *izb(p2,p7)*izb(p7,p8)*sign - 2.D0*zab2(p3,p5,p6,p4)*zba2(p6,
     &    p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*sign + 2.D0*zab2(p5,p3,p4,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *sign )
      adr87(jdu,1,1,2)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s178**(-1)*s234**(-1) * (  - za(p2,p3)*zb(p1,p7)**2*zba2(p4,p2,
     &    p3,p5)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s278**(-1)*s156**(-1) * ( za(p2,p8)**2*zb(p1,p6)*
     &    zba2(p4,p1,p6,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * ( za(p2,p3)*zb(p1,p6)*zba2(p4
     &    ,p2,p3,p8)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p7,p2
     &    ,p3,p4,p8)*s2347**(-1) )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1) * ( za(p2,p3)*za(p2,p8)*zb(p1,p6)*zba2(p7,
     &    p1,p6,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,
     &    p8)*s2347**(-1) )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s234**(-1) * ( za(p2,p3)*zb(p1,p6)*zb(p1,p7)*zba2(p4,
     &    p2,p3,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p2,p3,p4,
     &    p5)*s2347**(-1) )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2) * ( za(p2,p3)*za(p2,p8)*zb(p1,p6)*zb(p1,p7)*iza(p2,p7
     &    )*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*
     &    s2347**(-1) )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - za(p2,p3)*za(p3,p5)*
     &    zb(p1,p7)**2*zb(p4,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p1,p8
     &    )*izb(p7,p8)*sign - 2.D0*za(p2,p3)*zb(p1,p7)**2*zab2(p5,p3,p4
     &    ,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign
     &     - za(p2,p4)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,
     &    p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign + za(p2,p5)*za(p3,
     &    p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p5,p1,p7,p8)*iza(p7,p8)*izb(
     &    p1,p8)*izb(p7,p8)*sign + 2.D0*za(p2,p5)*zb(p1,p7)**2*zab2(p3,
     &    p5,p6,p4)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    sign + za(p2,p6)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p6,p1,
     &    p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p8)**2*za(p3,p5)
     &    *zb(p1,p3)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*sign - za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*
     &    zba2(p7,p2,p8,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + za(
     &    p2,p8)**2*za(p3,p5)*zb(p1,p5)*zb(p4,p6)*zba2(p7,p2,p8,p5)*
     &    iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + za(p2,p8)**2*za(p3,p5
     &    )*zb(p1,p6)*zb(p4,p6)*zba2(p7,p2,p8,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p7,p8)*sign - 2.D0*za(p2,p8)**2*zb(p1,p4)*zab2(p5,p3,p4,
     &    p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign
     &     + 2.D0*za(p2,p8)**2*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p7,p2,
     &    p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign )
      adr87(jdu,1,1,2) = adr87(jdu,1,1,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,
     &    p3)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(
     &    p7,p8)*sign + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*
     &    zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *sign + za(p2,p4)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(
     &    p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(
     &    p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p5)*zb(p1,p7)*zb(p4,p6)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - 2.D0*za(p2,p5)
     &    *za(p2,p8)*zb(p1,p6)*zb(p1,p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(p2,p6)*za(p2,p8)*
     &    za(p3,p5)*zb(p1,p6)*zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p1,p8)*izb(p7,p8)*sign )
      adr87(jdu,1,2,2)= + propw34**(-1)*propw56**(-1)*cxw**(-2)*
     & s178**(-1)*s234**(-1) * (  - za(p2,p3)*zb(p1,p7)*zba2(p4,p2,p3,
     &    p5)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) - za(p2,p3)*
     &    zb(p1,p8)*zba2(p4,p2,p3,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,
     &    p7,p8,p2) + za(p2,p3)*zba2(p4,p2,p3,p5)*zba2(p7,p1,p8,p2)*
     &    iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2) )
      adr87(jdu,1,2,2) = adr87(jdu,1,2,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s156**(-1)*s234**(-1) * (  - za(p2,p3)**2*zb(p1,p6)*
     &    zb(p3,p4)*zba2(p7,p1,p6,p5)*iza(p2,p8)*iza(p7,p8) + za(p2,p3)
     &    **2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*iza(p2,p8)*s2347**(-1) - za(p2,p3)**2*zb(p1,p6)*
     &    zb(p3,p4)*zba2(p8,p1,p6,p5)*iza(p2,p7)*iza(p7,p8) )
      adr87(jdu,1,2,2) = adr87(jdu,1,2,2) + propw34**(-1)*propw56**(-1)
     & *cxw**(-2)*s234**(-1) * (  - za(p2,p3)**2*zb(p3,p4)*zba2(p6,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p7,p2,p3,p4,p5)*
     &    s2347**(-1) )
      adr87(jdu,1,2,2) = adr87(jdu,1,2,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - za(p2,p3)*za(p3,p5)*
     &    zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p3,p1,p7,p8,p2
     &    )*sign - za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*zba3(p3,p1,p7,p8,p2)*sign + za(p2,p3)*za(p3,p5)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p3,p1,p7,p8,p2)*sign - 2.D0*za(p2,p3)*zb(p1,p7)*zab2(p5,
     &    p3,p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*sign - 2.
     &    D0*za(p2,p3)*zb(p1,p8)*zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8
     &    )*zba3(p4,p1,p7,p8,p2)*sign + 2.D0*za(p2,p3)*zab2(p5,p3,p4,p6
     &    )*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,
     &    p1,p7,p8,p2)*sign - za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*sign - za(p2,p4)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2)*sign + za(p2,p4)*za(p3,p5)*zb(p4,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)*
     &    sign )
      adr87(jdu,1,2,2) = adr87(jdu,1,2,2) + coupw(jdu,1)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * ( za(p2,p5)*za(p3,p5)*zb(p1
     &    ,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p5,p1,p7,p8,p2)*
     &    sign + za(p2,p5)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*zba3(p5,p1,p7,p8,p2)*sign - za(p2,p5)*za(p3,p5)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p5,p1,p7,p8,p2)*sign + 2.D0*za(p2,p5)*zb(p1,p7)*zab2(p3,
     &    p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*sign + 2.
     &    D0*za(p2,p5)*zb(p1,p8)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p7,p8
     &    )*zba3(p6,p1,p7,p8,p2)*sign - 2.D0*za(p2,p5)*zab2(p3,p5,p6,p4
     &    )*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,
     &    p1,p7,p8,p2)*sign + za(p2,p6)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*sign + za(p2,p6)*
     &    za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p6,
     &    p1,p7,p8,p2)*sign - za(p2,p6)*za(p3,p5)*zb(p4,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)*
     &    sign )
      adr87(jdu,2,1,1)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p1,p7)*za(p3,p5)*zb(p2,p3)*zb(p4,
     &    p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*sign - za(p1,p7)*
     &    za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*
     &    izb(p7,p8)*sign + za(p1,p7)*za(p3,p5)*zb(p2,p5)*zb(p4,p6)*
     &    zab2(p5,p7,p8,p1)*izb(p1,p8)*izb(p7,p8)*sign + za(p1,p7)*za(
     &    p3,p5)*zb(p2,p6)*zb(p4,p6)*zab2(p6,p7,p8,p1)*izb(p1,p8)*izb(
     &    p7,p8)*sign - 2.D0*za(p1,p7)*zb(p2,p4)*zab2(p3,p7,p8,p1)*
     &    zab2(p5,p3,p4,p6)*izb(p1,p8)*izb(p7,p8)*sign + 2.D0*za(p1,p7)
     &    *zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p7,p8,p1)*izb(p1,p8)*
     &    izb(p7,p8)*sign - za(p1,p8)*za(p3,p5)*zb(p2,p3)*zb(p4,p6)*
     &    zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*sign - za(p1,p8)*za(
     &    p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*izb(
     &    p7,p8)*sign + za(p1,p8)*za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5
     &    ,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*sign + za(p1,p8)*za(p3,p5)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p6,p7,p8,p1)*izb(p1,p7)*izb(p7,p8)*
     &    sign )
      adr87(jdu,2,1,1) = adr87(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * (  - 2.D0*za(p1,p8)*zb(p2,
     &    p4)*zab2(p3,p7,p8,p1)*zab2(p5,p3,p4,p6)*izb(p1,p7)*izb(p7,p8)
     &    *sign + 2.D0*za(p1,p8)*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p7
     &    ,p8,p1)*izb(p1,p7)*izb(p7,p8)*sign - za(p3,p5)*za(p7,p8)*zb(
     &    p2,p3)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*sign
     &     - za(p3,p5)*za(p7,p8)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p7,p8,p1)*
     &    izb(p1,p7)*izb(p1,p8)*sign + za(p3,p5)*za(p7,p8)*zb(p2,p5)*
     &    zb(p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*izb(p1,p8)*sign + za(
     &    p3,p5)*za(p7,p8)*zb(p2,p6)*zb(p4,p6)*zab2(p6,p7,p8,p1)*izb(p1
     &    ,p7)*izb(p1,p8)*sign - 2.D0*za(p7,p8)*zb(p2,p4)*zab2(p3,p7,p8
     &    ,p1)*zab2(p5,p3,p4,p6)*izb(p1,p7)*izb(p1,p8)*sign + 2.D0*za(
     &    p7,p8)*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p7,p8,p1)*izb(p1,
     &    p7)*izb(p1,p8)*sign )
      adr87(jdu,2,1,1) = adr87(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p1,p3)*za(p3,p5)*
     &    zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*sign + za(p1,p3)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*
     &    zab2(p7,p2,p8,p3)*izb(p1,p8)*izb(p7,p8)*sign + za(p1,p3)*za(
     &    p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*izb(
     &    p7,p8)*sign - 2.D0*za(p1,p3)*zb(p2,p1)**2*zab2(p5,p3,p4,p6)*
     &    zab2(p8,p2,p7,p4)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign + 2.D0
     &    *za(p1,p3)*zb(p2,p1)*zab2(p5,p3,p4,p6)*zab2(p7,p2,p8,p4)*izb(
     &    p1,p8)*izb(p7,p8)*sign + 2.D0*za(p1,p3)*zb(p2,p1)*zab2(p5,p3,
     &    p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p7,p8)*sign - za(p1,
     &    p4)*za(p3,p5)*zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p2
     &    ,p7)*izb(p1,p7)*izb(p1,p8)*sign + za(p1,p4)*za(p3,p5)*zb(p2,
     &    p1)*zb(p4,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*izb(p7,p8)*sign + 
     &    za(p1,p4)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*
     &    izb(p1,p7)*izb(p7,p8)*sign )
      adr87(jdu,2,1,1) = adr87(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * ( za(p1,p5)*za(p3,p5)*zb(p2
     &    ,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p5)*izb(p2,p7)*izb(p1,p7)*
     &    izb(p1,p8)*sign - za(p1,p5)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*
     &    zab2(p7,p2,p8,p5)*izb(p1,p8)*izb(p7,p8)*sign - za(p1,p5)*za(
     &    p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p5)*izb(p1,p7)*izb(
     &    p7,p8)*sign + 2.D0*za(p1,p5)*zb(p2,p1)**2*zab2(p3,p5,p6,p4)*
     &    zab2(p8,p2,p7,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign - 2.D0
     &    *za(p1,p5)*zb(p2,p1)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p6)*izb(
     &    p1,p8)*izb(p7,p8)*sign - 2.D0*za(p1,p5)*zb(p2,p1)*zab2(p3,p5,
     &    p6,p4)*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p7,p8)*sign + za(p1,
     &    p6)*za(p3,p5)*zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p2
     &    ,p7)*izb(p1,p7)*izb(p1,p8)*sign - za(p1,p6)*za(p3,p5)*zb(p2,
     &    p1)*zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*izb(p7,p8)*sign - 
     &    za(p1,p6)*za(p3,p5)*zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p6)*
     &    izb(p1,p7)*izb(p7,p8)*sign )
      adr87(jdu,2,1,1) = adr87(jdu,2,1,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,
     &    p3)*zb(p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign + za(p3,
     &    p5)*za(p4,p8)*zb(p2,p1)*zb(p2,p4)*zb(p4,p6)*izb(p2,p7)*izb(p1
     &    ,p7)*izb(p1,p8)*sign - za(p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p5
     &    )*zb(p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign - za(p3,p5)
     &    *za(p6,p8)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*sign + 2.D0*za(p3,p8)*zb(p2,p1)*zb(p2,p4)*zab2(
     &    p5,p3,p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*sign - 2.D0*za(
     &    p5,p8)*zb(p2,p1)*zb(p2,p6)*zab2(p3,p5,p6,p4)*izb(p2,p7)*izb(
     &    p1,p7)*izb(p1,p8)*sign )
      adr87(jdu,2,1,2)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p1,p8)*za(p3,p5)*zb(p2,p3)*zb(p1,
     &    p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*sign - za(p1,p8)*za(p3,p5)*zb(p2,p4)*zb(p1,p7)*zb(p4,p6)*
     &    zab2(p4,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign + za(
     &    p1,p8)*za(p3,p5)*zb(p2,p5)*zb(p1,p7)*zb(p4,p6)*zab2(p5,p1,p8,
     &    p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign + za(p1,p8)*za(p3,
     &    p5)*zb(p2,p6)*zb(p1,p7)*zb(p4,p6)*zab2(p6,p1,p8,p7)*iza(p7,p8
     &    )*izb(p1,p8)*izb(p7,p8)*sign - 2.D0*za(p1,p8)*zb(p2,p4)*zb(p1
     &    ,p7)*zab2(p3,p1,p8,p7)*zab2(p5,p3,p4,p6)*iza(p7,p8)*izb(p1,p8
     &    )*izb(p7,p8)*sign + 2.D0*za(p1,p8)*zb(p2,p6)*zb(p1,p7)*zab2(
     &    p3,p5,p6,p4)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*sign )
      adr87(jdu,2,1,2) = adr87(jdu,2,1,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p2,p8)*za(p1,p3)*
     &    za(p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p3)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p7,p8)*sign - 2.D0*za(p2,p8)*za(p1,p3)*zb(p2,
     &    p7)*zab2(p5,p3,p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p7,p8)*sign - za(p2,p8)*za(p1,p4)*za(p3,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    sign + za(p2,p8)*za(p1,p5)*za(p3,p5)*zb(p2,p7)*zb(p4,p6)*
     &    zab2(p8,p2,p7,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + 2.D0
     &    *za(p2,p8)*za(p1,p5)*zb(p2,p7)*zab2(p3,p5,p6,p4)*zab2(p8,p2,
     &    p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*sign + za(p2,p8)*za(
     &    p1,p6)*za(p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p2
     &    ,p7)*iza(p7,p8)*izb(p7,p8)*sign )
      adr87(jdu,2,1,2) = adr87(jdu,2,1,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p3,p5)*zb(p4,p6)*zab2(p3,p1,p8,p7
     &    )*zab2(p8,p2,p7,p3)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8)*sign + za(p3,p5)*zb(p4,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,
     &    p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(
     &    p3,p5)*zb(p4,p6)*zab2(p5,p1,p8,p7)*zab2(p8,p2,p7,p5)*iza(p2,
     &    p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*sign - za(p3,p5)*zb(p4,
     &    p6)*zab2(p6,p1,p8,p7)*zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)
     &    *izb(p1,p8)*izb(p7,p8)*sign + 2.D0*zab2(p3,p1,p8,p7)*zab2(p5,
     &    p3,p4,p6)*zab2(p8,p2,p7,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8)*sign - 2.D0*zab2(p3,p5,p6,p4)*zab2(p5,p1,p8,p7)*
     &    zab2(p8,p2,p7,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &    *sign )
      adr87(jdu,2,2,1)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p1,p7)**2*za(p3,p5)*zb(p2,p3)*zb(
     &    p4,p6)*zab2(p3,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    sign - za(p1,p7)**2*za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p1,
     &    p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + za(p1,p7)**2*
     &    za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5,p1,p7,p8)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8)*sign + za(p1,p7)**2*za(p3,p5)*zb(p2,p6)
     &    *zb(p4,p6)*zab2(p6,p1,p7,p8)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &    *sign - 2.D0*za(p1,p7)**2*zb(p2,p4)*zab2(p3,p1,p7,p8)*zab2(p5
     &    ,p3,p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*sign + 2.D0*za(p1
     &    ,p7)**2*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p7,p8)*iza(p1,
     &    p8)*iza(p7,p8)*izb(p7,p8)*sign )
      adr87(jdu,2,2,1) = adr87(jdu,2,2,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p1,p3)*za(p3,p5)*
     &    zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p2,p7
     &    )*izb(p7,p8)*sign - 2.D0*za(p1,p3)*zb(p2,p8)**2*zab2(p5,p3,p4
     &    ,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign
     &     - za(p1,p4)*za(p3,p5)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,
     &    p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign + za(p1,p5)*za(p3,
     &    p5)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p5)*iza(p7,p8)*izb(
     &    p2,p7)*izb(p7,p8)*sign + 2.D0*za(p1,p5)*zb(p2,p8)**2*zab2(p3,
     &    p5,p6,p4)*zab2(p7,p2,p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    sign + za(p1,p6)*za(p3,p5)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,
     &    p8,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign )
      adr87(jdu,2,2,1) = adr87(jdu,2,2,1) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * ( za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,
     &    p3)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8)*sign + 2.D0*za(p1,p3)*za(p1,p7)*zb(p2,p4)*zb(p2,p8)*
     &    zab2(p5,p3,p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &    *sign + za(p1,p4)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p2,p8)*zb(
     &    p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(
     &    p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p5)*zb(p2,p8)*zb(p4,p6)*iza(
     &    p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - 2.D0*za(p1,p5)
     &    *za(p1,p7)*zb(p2,p6)*zb(p2,p8)*zab2(p3,p5,p6,p4)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*sign - za(p1,p6)*za(p1,p7)*
     &    za(p3,p5)*zb(p2,p6)*zb(p2,p8)*zb(p4,p6)*iza(p1,p8)*iza(p7,p8)
     &    *izb(p2,p7)*izb(p7,p8)*sign )
      adr87(jdu,2,2,2)= + coupw(jdu,2)*propw34**(-1)*propw56**(-1)*
     & cxw**(-1)*s178**(-1) * (  - za(p2,p1)**2*za(p3,p5)*zb(p2,p3)*zb(
     &    p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    sign - za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p4,p6)*zab2(p4,p1,
     &    p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign + za(p2,p1)**2*
     &    za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*sign + za(p2,p1)**2*za(p3,p5)*zb(p2,p6)
     &    *zb(p4,p6)*zab2(p6,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)
     &    *sign - 2.D0*za(p2,p1)**2*zb(p2,p4)*zab2(p3,p1,p8,p7)*zab2(p5
     &    ,p3,p4,p6)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign + 2.D0*za(p2
     &    ,p1)**2*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p8,p7)*iza(p2,
     &    p7)*iza(p2,p8)*iza(p1,p8)*sign + za(p2,p1)*za(p3,p5)*zb(p2,p3
     &    )*zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*sign + 
     &    za(p2,p1)*za(p3,p5)*zb(p2,p3)*zb(p4,p6)*zab2(p3,p1,p8,p7)*
     &    iza(p2,p8)*iza(p7,p8)*sign + za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zab2(p4,p1,p7,p8)*iza(p2,p7)*iza(p7,p8)*sign )
      adr87(jdu,2,2,2) = adr87(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s178**(-1) * ( za(p2,p1)*za(p3,p5)*zb(p2
     &    ,p4)*zb(p4,p6)*zab2(p4,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign
     &     - za(p2,p1)*za(p3,p5)*zb(p2,p5)*zb(p4,p6)*zab2(p5,p1,p7,p8)*
     &    iza(p2,p7)*iza(p7,p8)*sign - za(p2,p1)*za(p3,p5)*zb(p2,p5)*
     &    zb(p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign - za(
     &    p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p4,p6)*zab2(p6,p1,p7,p8)*iza(p2
     &    ,p7)*iza(p7,p8)*sign - za(p2,p1)*za(p3,p5)*zb(p2,p6)*zb(p4,p6
     &    )*zab2(p6,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign + 2.D0*za(p2,
     &    p1)*zb(p2,p4)*zab2(p3,p1,p7,p8)*zab2(p5,p3,p4,p6)*iza(p2,p7)*
     &    iza(p7,p8)*sign + 2.D0*za(p2,p1)*zb(p2,p4)*zab2(p3,p1,p8,p7)*
     &    zab2(p5,p3,p4,p6)*iza(p2,p8)*iza(p7,p8)*sign - 2.D0*za(p2,p1)
     &    *zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p7,p8)*iza(p2,p7)*
     &    iza(p7,p8)*sign - 2.D0*za(p2,p1)*zb(p2,p6)*zab2(p3,p5,p6,p4)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p8)*iza(p7,p8)*sign )
      adr87(jdu,2,2,2) = adr87(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * ( za(p1,p3)*za(p3,p5)*zb(p2
     &    ,p7)*zb(p4,p6)*zab2(p2,p7,p8,p3)*iza(p2,p8)*iza(p7,p8)*sign
     &     + za(p1,p3)*za(p3,p5)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p3)*
     &    iza(p2,p7)*iza(p7,p8)*sign + za(p1,p3)*za(p3,p5)*zb(p4,p6)*
     &    zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*iza(p2,p8)*sign + 2.D0
     &    *za(p1,p3)*zb(p2,p7)*zab2(p2,p7,p8,p4)*zab2(p5,p3,p4,p6)*iza(
     &    p2,p8)*iza(p7,p8)*sign + 2.D0*za(p1,p3)*zb(p2,p8)*zab2(p2,p7,
     &    p8,p4)*zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*sign + 2.D0*
     &    za(p1,p3)*zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p5,p3,p4,p6)*iza(
     &    p2,p7)*iza(p2,p8)*sign + za(p1,p4)*za(p3,p5)*zb(p2,p7)*zb(p4,
     &    p6)*zab2(p2,p7,p8,p4)*iza(p2,p8)*iza(p7,p8)*sign + za(p1,p4)*
     &    za(p3,p5)*zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p4)*iza(p2,p7)*
     &    iza(p7,p8)*sign + za(p1,p4)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*
     &    zab2(p2,p7,p8,p4)*iza(p2,p7)*iza(p2,p8)*sign - za(p1,p5)*za(
     &    p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,p5)*iza(p2,p8)*iza(
     &    p7,p8)*sign )
      adr87(jdu,2,2,2) = adr87(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1)*s278**(-1) * (  - za(p1,p5)*za(p3,p5)*
     &    zb(p2,p8)*zb(p4,p6)*zab2(p2,p7,p8,p5)*iza(p2,p7)*iza(p7,p8)*
     &    sign - za(p1,p5)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,
     &    p5)*iza(p2,p7)*iza(p2,p8)*sign - 2.D0*za(p1,p5)*zb(p2,p7)*
     &    zab2(p2,p7,p8,p6)*zab2(p3,p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*
     &    sign - 2.D0*za(p1,p5)*zb(p2,p8)*zab2(p2,p7,p8,p6)*zab2(p3,p5,
     &    p6,p4)*iza(p2,p7)*iza(p7,p8)*sign - 2.D0*za(p1,p5)*zb(p7,p8)*
     &    zab2(p2,p7,p8,p6)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p2,p8)*
     &    sign - za(p1,p6)*za(p3,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p2,p7,p8,
     &    p6)*iza(p2,p8)*iza(p7,p8)*sign - za(p1,p6)*za(p3,p5)*zb(p2,p8
     &    )*zb(p4,p6)*zab2(p2,p7,p8,p6)*iza(p2,p7)*iza(p7,p8)*sign - 
     &    za(p1,p6)*za(p3,p5)*zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*
     &    iza(p2,p7)*iza(p2,p8)*sign )
      adr87(jdu,2,2,2) = adr87(jdu,2,2,2) + coupw(jdu,2)*propw34**(-1)*
     & propw56**(-1)*cxw**(-1) * (  - za(p2,p1)*za(p1,p3)*za(p3,p5)*zb(
     &    p3,p7)*zb(p4,p6)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign - 2.D0
     &    *za(p2,p1)*za(p1,p3)*zb(p4,p7)*zab2(p5,p3,p4,p6)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*sign - za(p2,p1)*za(p1,p4)*za(p3,p5)*
     &    zb(p4,p6)*zb(p4,p7)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign + 
     &    za(p2,p1)*za(p1,p5)*za(p3,p5)*zb(p4,p6)*zb(p5,p7)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*sign + 2.D0*za(p2,p1)*za(p1,p5)*zb(p6,
     &    p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*sign
     &     + za(p2,p1)*za(p1,p6)*za(p3,p5)*zb(p4,p6)*zb(p6,p7)*iza(p2,
     &    p7)*iza(p2,p8)*iza(p1,p8)*sign )
      endif
      enddo
      enddo
      a78=adr78+asr78
      a87=adr87+asr87
      aq=a78+a87

      msq(:)=zip
      do jdu=1,2
      do h1=1,2
      do h7=1,2
      do h8=1,2
C23456789012345678901234567890123456789012345678901234567890123456789012
      msq(jdu)=msq(jdu)+gsq**2*esq**4*V*xn*(
     & +real(a78(jdu,h1,h7,h8)*conjg(a78(jdu,h1,h7,h8)))
     & +real(a87(jdu,h1,h7,h8)*conjg(a87(jdu,h1,h7,h8)))
     & -real(aq(jdu,h1,h7,h8)*conjg(aq(jdu,h1,h7,h8)))
     &  /xn**2)
      enddo
      enddo
      enddo
      enddo





      return
      end
