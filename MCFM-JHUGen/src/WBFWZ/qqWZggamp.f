      subroutine qqWZggamp(i1,i2,i3,i4,i5,i6,i7,i8,
     & b7,b8,za,zb,amp78,amp87)
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
      integer:: jdu,jdd,j78,i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,rxw,zba3,zba4,iza,izb,
     & amp78(2,2,2),amp87(2,2,2),bmp78(2,2,2),bmp87(2,2,2),
     & propw34,propz56,propw3456,gamzew56(2),
     & gamzee56(2,2),gamzne56(2,2),gamzqe56(2,2,2)
      real(dp):: t3,t4,s356,s456,
     & s34,s56,s278,s178,s156,s234,s2347,s3456,
     & s1348,s2567,s134,s256,s345
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
     &             +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      p1=i1
      p2=i2
      p5=i5
      p6=i6
      p3=i3
      p4=i4

      rxw=sqrt((cone-cxw)/cxw)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s278=t3(p2,i7,i8)
      s178=t3(p1,i7,i8)
      s134=t3(p1,p3,p4)
      s156=t3(p1,p5,p6)
      s256=t3(p2,p5,p6)
      s234=t3(p2,p3,p4)
      s345=t3(p3,p4,p5)
      s356=t3(p3,p5,p6)
      s456=t3(p4,p5,p6)
      s3456=t4(p3,p4,p5,p6)

      propw3456=s3456-cwmass2
      propw34=s34-cwmass2
      propz56=s56-czmass2


      gamzew56(1)=qe/s56+rxw*le/propz56
      gamzew56(2)=qe/s56+rxw*re/propz56
      gamzee56(1,1)=qe**2/s56+le**2/propz56
      gamzee56(1,2)=qe**2/s56+le*re/propz56
      gamzee56(2,1)=qe**2/s56+re*le/propz56
      gamzee56(2,2)=qe**2/s56+re**2/propz56
      gamzne56(1,1)=ln*le/propz56
      gamzne56(1,2)=ln*re/propz56
      gamzne56(2,1)=rn*le/propz56
      gamzne56(2,2)=rn*re/propz56
      do jdu=1,2
      gamzqe56(jdu,1,1)=Q(jdu)*qe/s56+L(jdu)*le/propz56
      gamzqe56(jdu,1,2)=Q(jdu)*qe/s56+L(jdu)*re/propz56
      gamzqe56(jdu,2,1)=Q(jdu)*qe/s56+R(jdu)*le/propz56
      gamzqe56(jdu,2,2)=Q(jdu)*qe/s56+R(jdu)*re/propz56
      enddo

      jdu=2
      jdd=1

      do j78=1,2
      if (j78 == 1) then
      p7=i7
      p8=i8
      elseif (j78 == 2) then
      p7=i8
      p8=i7
      endif

      s2347=t4(p2,p3,p4,p7)
      s2567=t4(p2,p5,p6,p7)
      s1348=t4(p1,p3,p4,p8)

      if (j78==1) then
C---  amp78(h56,h7,h8)
      amp78(1,1,1)= + gamzqe56(jdu,1,1)*propw34**(-1)*s278**(-1) * ( 
     &     - 2.D0*za(p2,p7)*zb(p1,p6)*zba2(p4,p1,p6,p5)*izb(p1,p8)*izb(
     &    p7,p8)*zba3(p1,p2,p7,p8,p3)*s156**(-1) - 2.D0*za(p2,p8)*zb(p1
     &    ,p6)*zba2(p4,p1,p6,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8
     &    ,p3)*s156**(-1) - 2.D0*zb(p1,p6)*zba2(p1,p2,p7,p8)*zba2(p4,p1
     &    ,p6,p5)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &    *s156**(-1) )
      amp78(1,1,1) = amp78(1,1,1) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*zba2(p4,p2,p3,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba4(p1,p2,p3,p4,p7,p8)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p6)**2
     &    *zba2(p4,p2,p3,p7)*izb(p1,p8)*izb(p7,p8)*s156**(-1)*
     &    s234**(-1) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*zba2(p4,p2
     &    ,p3,p8)*izb(p1,p7)*izb(p7,p8)*s156**(-1)*s234**(-1) + 2.D0*
     &    za(p5,p6)*zb(p1,p6)**2*zba2(p1,p2,p7,p3)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*s2347**(-1)*s156**(-1) )
      amp78(1,1,1) = amp78(1,1,1) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s278**(-1) * (  - 2.D0*za(p2,p7)*zb(p1,p4)*zba2(p6,p1,p4,p3)*
     &    izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s134**(-1) - 2.D0*
     &    za(p2,p8)*zb(p1,p4)*zba2(p6,p1,p4,p3)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s134**(-1) - 2.D0*zb(p1,p4)*zba2(p1,p2,
     &    p7,p8)*zba2(p6,p1,p4,p3)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s134**(-1) )
      amp78(1,1,1) = amp78(1,1,1) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p4)**2*zba2(p1,p3,p4,p8)*
     &    zba2(p6,p2,p5,p7)*izb(p1,p7)*izb(p1,p8)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p4)**2
     &    *zba2(p6,p2,p5,p7)*izb(p1,p8)*izb(p7,p8)*s134**(-1)*
     &    s256**(-1) + 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p4)**2*zba2(p6,p2
     &    ,p5,p8)*izb(p1,p7)*izb(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*
     &    za(p3,p4)*zb(p1,p4)**2*zba2(p1,p2,p7,p5)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p6,p2,p5,p7,p8)*s2567**(-1)*s134**(-1) )
      amp78(1,1,1) = amp78(1,1,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p3
     &    )*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8
     &    )*zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p7)*zb(p1,p4)*zab2(p5,p3,
     &    p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p7)*zb(p1,p6)*zab2(p3,p5,p6,p4)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5) + 2.D0*za(p2,p8)*za(p3,p5)*zb(p1,p3)*zb(
     &    p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(
     &    p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(p2,
     &    p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1
     &    ,p2,p7,p8,p5) + 2.D0*za(p3,p5)*zb(p1,p3)*zb(p4,p6)*zba2(p1,p2
     &    ,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &     + 2.D0*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(
     &    p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4) )
      amp78(1,1,1) = amp78(1,1,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*zb(p1,p4)*zab2(p5,p3,p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p1,p2,p7,
     &    p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5) )
      amp78(1,2,1)= + gamzqe56(jdu,1,1)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)**2*zba2(p4,p2,p3,p5)*zba2(p6,p1,
     &    p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s234**(-1) )
      amp78(1,2,1) = amp78(1,2,1) + gamzqe56(jdu,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p6)*zba2(p4,p1,p6,p5)*
     &    zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s156**(-1)
     &     )
      amp78(1,2,1) = amp78(1,2,1) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p6)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*
     &    s2347**(-1) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p6)*zba2(p7,p1,
     &    p6,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p8)*
     &    s2347**(-1)*s156**(-1) + 2.D0*za(p2,p3)*zb(p1,p6)*zb(p1,p7)*
     &    zba2(p4,p2,p3,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p2
     &    ,p3,p4,p5)*s2347**(-1)*s234**(-1) + 2.D0*za(p2,p3)*zb(p1,p6)*
     &    zba2(p4,p2,p3,p8)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p2,p3,p4,p8)*s2347**(-1)*s156**(-1)*s234**(-1) )
      amp78(1,2,1) = amp78(1,2,1) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*zb(p1,p7)**2*zba2(p4,p1,p7,p8)*
     &    zba2(p6,p2,p5,p3)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s256**(-1)
     &     )
      amp78(1,2,1) = amp78(1,2,1) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p4)*zba2(p6,p1,p4,p3)*
     &    zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s134**(-1)
     &     )
      amp78(1,2,1) = amp78(1,2,1) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p5)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p6,p2,p5,p7,p3)*
     &    s2567**(-1) + 2.D0*za(p2,p5)*za(p2,p8)*zb(p1,p4)*zba2(p7,p1,
     &    p4,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p6,p2,p5,p7,p8)*
     &    s2567**(-1)*s134**(-1) - 2.D0*za(p2,p5)*zb(p1,p4)*zb(p1,p7)*
     &    zba2(p6,p2,p5,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p1
     &    ,p4,p8,p3)*s1348**(-1)*s256**(-1) - 2.D0*za(p2,p5)*zb(p1,p4)*
     &    zba2(p6,p2,p5,p8)*zba2(p7,p1,p4,p3)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p1,p3,p4,p8)*s1348**(-1)*s134**(-1)*s256**(-1) )
      amp78(1,2,1) = amp78(1,2,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1
     &    ,p7)**2*zb(p4,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p3)*zb(p1,p7)**2*zab2(p5,p3,p4,p6)*
     &    zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p4)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,p8)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*za(p2,p5)*zb(p1,p7)**
     &    2*zab2(p3,p5,p6,p4)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) )
      amp78(1,2,1) = amp78(1,2,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p5)*
     &    zb(p1,p3)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*
     &    zba2(p7,p2,p8,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p8)**2*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p7,p2,p8,p3)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p8)**2*zb(p1,p6)*
     &    zab2(p3,p5,p6,p4)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      amp78(1,2,1) = amp78(1,2,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p3)
     &    *zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*zab2(p5,p3
     &    ,p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*
     &    za(p2,p4)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(p2,p5)*
     &    za(p2,p8)*zb(p1,p6)*zb(p1,p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      amp78(1,1,2)= + gamzqe56(jdu,1,1)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(p4,p2,p3,p5)*zba2(
     &    p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s234**(-1) )
      amp78(1,1,2) = amp78(1,1,2) + gamzqe56(jdu,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p6)*zba2(p4,p1,p6,
     &    p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s156**(-1) )
      amp78(1,1,2) = amp78(1,1,2) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(p8,p1,p6,p5
     &    )*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2,p3,p4,p7)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)*zba2(p4,p2,p3,p7)*
     &    zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba4(p8,p2
     &    ,p3,p4,p7,p5)*s2347**(-1)*s234**(-1) - 2.D0*zb(p1,p6)*zba2(p4
     &    ,p2,p3,p7)*zba2(p8,p2,p7,p3)*zba2(p8,p1,p6,p5)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*s2347**(-1)*s156**(-1) + 2.D0*zba2(p6,
     &    p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*s2347**(-1) )
      amp78(1,1,2) = amp78(1,1,2) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p1,p7)*zb(p1,p8)*zba2(p4,p1,
     &    p8,p7)*zba2(p6,p2,p5,p3)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s256**(-1) )
      amp78(1,1,2) = amp78(1,1,2) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zba2(p6,p1,p4,
     &    p3)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s134**(-1) )
      amp78(1,1,2) = amp78(1,1,2) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p5)*zb(p1,p4)*zba2(p6,p2,p5,p7)*zba2(p8,p1,p4
     &    ,p3)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p1,p3,p4,p7)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p5)*zba2(p4,p1,p8,p7)*
     &    zba2(p6,p2,p5,p7)*zba2(p8,p1,p4,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s1348**(-1)*s256**(-1) - 2.D0*zb(p1,p4)*zba2(p6,p2
     &    ,p5,p7)*zba2(p8,p2,p7,p5)*zba2(p8,p1,p4,p3)*iza(p7,p8)*izb(p2
     &    ,p7)*izb(p7,p8)*s2567**(-1)*s134**(-1) + 2.D0*zba2(p4,p1,p8,
     &    p7)*zba2(p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7
     &    ,p8)*zba3(p6,p2,p5,p7,p3)*s2567**(-1) )
      amp78(1,1,2) = amp78(1,1,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p3
     &    ,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p3,p1,p8,p7)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8) - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zab2(p5,
     &    p3,p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p4)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(
     &    p2,p5)*za(p1,p7)*zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p8,p7
     &    )*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      amp78(1,1,2) = amp78(1,1,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p2
     &    ,p8)*zb(p1,p3)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8) - 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p4)
     &    *zb(p4,p6)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     - 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(
     &    p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*za(p2,p7
     &    )*zb(p2,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p5)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      amp78(1,1,2) = amp78(1,1,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p3,p5)*zb(p4,p6)*zba2(p3,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     + 2.D0*za(p3,p5)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - 2.D0*zab2(
     &    p3,p5,p6,p4)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*zab2(p5,p3,p4,p6)*
     &    zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8) )
      amp78(1,2,2)= + gamzqe56(jdu,1,1)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)*zba2(p4,p2,p3,p5)*iza(p2,p8)*iza(
     &    p7,p8)*zba3(p6,p1,p7,p8,p2)*s234**(-1) - 2.D0*za(p2,p3)*zb(p1
     &    ,p8)*zba2(p4,p2,p3,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8
     &    ,p2)*s234**(-1) + 2.D0*za(p2,p3)*zba2(p4,p2,p3,p5)*zba2(p7,p1
     &    ,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)
     &    *s234**(-1) )
      amp78(1,2,2) = amp78(1,2,2) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p1,p6,p5)*
     &    iza(p2,p8)*iza(p7,p8)*s156**(-1)*s234**(-1) + 2.D0*za(p2,p3)
     &    **2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*iza(p2,p8)*s2347**(-1)*s156**(-1)*s234**(-1) - 2.D0
     &    *za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zba2(p8,p1,p6,p5)*iza(p2,p7
     &    )*iza(p7,p8)*s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)**2*zb(p3,
     &    p4)*zba2(p6,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(
     &    p7,p2,p3,p4,p5)*s2347**(-1)*s234**(-1) )
      amp78(1,2,2) = amp78(1,2,2) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*zb(p1,p7)*zba2(p6,p2,p5,p3)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s256**(-1) - 2.D0*
     &    za(p2,p5)*zb(p1,p8)*zba2(p6,p2,p5,p3)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) + 2.D0*za(p2,p5)*zba2(p6,p2,
     &    p5,p3)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) )
      amp78(1,2,2) = amp78(1,2,2) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p5)**2*zb(p1,p4)*zb(p5,p6)*zba2(p7,p1,p4,p3)*
     &    iza(p2,p8)*iza(p7,p8)*s134**(-1)*s256**(-1) - 2.D0*za(p2,p5)
     &    **2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*iza(p2,p7)*iza(p2,
     &    p8)*zba4(p7,p1,p3,p4,p8,p2)*s1348**(-1)*s134**(-1)*s256**(-1)
     &     - 2.D0*za(p2,p5)**2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*
     &    iza(p2,p7)*iza(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*za(p2,p5)
     &    **2*zb(p5,p6)*zba2(p4,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,
     &    p8)*zba3(p7,p1,p4,p8,p3)*s1348**(-1)*s256**(-1) )
      amp78(1,2,2) = amp78(1,2,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1
     &    ,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p3,p1,p7,p8,p2) - 2.
     &    D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,
     &    p8)*zba3(p3,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p3,p5)*zb(p4,p6)
     &    *zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p3,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p3)*zb(p1,p7)*zab2(p5,p3,p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*
     &    zb(p1,p8)*zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,
     &    p7,p8,p2) + 2.D0*za(p2,p3)*zab2(p5,p3,p4,p6)*zba2(p7,p1,p8,p2
     &    )*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0
     &    *za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8
     &    )*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p4)*za(p3,p5)*zb(p1,p8)*
     &    zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p4)*za(p3,p5)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) )
      amp78(1,2,2) = amp78(1,2,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * ( 2.D0*za(p2,p5)*zb(p1,p7)*zab2(p3,
     &    p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p5)*zb(p1,p8)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p6,p1,p7,p8,p2) - 2.D0*za(p2,p5)*zab2(p3,p5,p6,p4)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,
     &    p8,p2) )
      amp78(2,1,1)= + gamzqe56(jdu,1,2)*propw34**(-1)*s278**(-1) * ( 
     &     - 2.D0*za(p2,p7)*zb(p1,p5)*zba2(p4,p1,p5,p6)*izb(p1,p8)*izb(
     &    p7,p8)*zba3(p1,p2,p7,p8,p3)*s156**(-1) - 2.D0*za(p2,p8)*zb(p1
     &    ,p5)*zba2(p4,p1,p5,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8
     &    ,p3)*s156**(-1) - 2.D0*zb(p1,p5)*zba2(p1,p2,p7,p8)*zba2(p4,p1
     &    ,p5,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &    *s156**(-1) )
      amp78(2,1,1) = amp78(2,1,1) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p5)**2*zba2(p4,p2,p3,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba4(p1,p2,p3,p4,p7,p8)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p5)**2
     &    *zba2(p4,p2,p3,p7)*izb(p1,p8)*izb(p7,p8)*s156**(-1)*
     &    s234**(-1) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p5)**2*zba2(p4,p2
     &    ,p3,p8)*izb(p1,p7)*izb(p7,p8)*s156**(-1)*s234**(-1) - 2.D0*
     &    za(p5,p6)*zb(p1,p5)**2*zba2(p1,p2,p7,p3)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*s2347**(-1)*s156**(-1) )
      amp78(2,1,1) = amp78(2,1,1) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s278**(-1) * (  - 2.D0*za(p2,p7)*zb(p1,p4)*zba2(p5,p1,p4,p3)*
     &    izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p6)*s134**(-1) - 2.D0*
     &    za(p2,p8)*zb(p1,p4)*zba2(p5,p1,p4,p3)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p6)*s134**(-1) - 2.D0*zb(p1,p4)*zba2(p1,p2,
     &    p7,p8)*zba2(p5,p1,p4,p3)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p6)*s134**(-1) )
      amp78(2,1,1) = amp78(2,1,1) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*zba2(p1,p3,p4,p8)*
     &    zba2(p5,p2,p6,p7)*izb(p1,p7)*izb(p1,p8)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p4)**2
     &    *zba2(p5,p2,p6,p7)*izb(p1,p8)*izb(p7,p8)*s134**(-1)*
     &    s256**(-1) + 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*zba2(p5,p2
     &    ,p6,p8)*izb(p1,p7)*izb(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*
     &    za(p3,p4)*zb(p1,p4)**2*zba2(p1,p2,p7,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p5,p2,p6,p7,p8)*s2567**(-1)*s134**(-1) )
      amp78(2,1,1) = amp78(2,1,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p3
     &    )*zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0
     &    *za(p2,p7)*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,p8)*izb(p7,p8
     &    )*zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p7)*zb(p1,p4)*zba2(p5,p3,
     &    p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p7)*zb(p1,p5)*zab2(p3,p5,p6,p4)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p6) + 2.D0*za(p2,p8)*za(p3,p6)*zb(p1,p3)*zb(
     &    p4,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(
     &    p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p8)*zb(p1,p4)*zba2(p5,p3,p4
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(p2,
     &    p8)*zb(p1,p5)*zab2(p3,p5,p6,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1
     &    ,p2,p7,p8,p6) + 2.D0*za(p3,p6)*zb(p1,p3)*zb(p4,p5)*zba2(p1,p2
     &    ,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &     + 2.D0*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*zba2(p1,p2,p7,p8)*izb(
     &    p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4) )
      amp78(2,1,1) = amp78(2,1,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*zb(p1,p4)*zba2(p1,p2,p7,p8)*
     &    zba2(p5,p3,p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*zb(p1,p5)*zab2(p3,p5,p6,p4)*zba2(p1,p2,p7,
     &    p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p6) )
      amp78(2,2,1)= + gamzqe56(jdu,1,2)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)**2*zba2(p4,p2,p3,p6)*zba2(p5,p1,
     &    p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s234**(-1) )
      amp78(2,2,1) = amp78(2,2,1) + gamzqe56(jdu,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p5)*zba2(p4,p1,p5,p6)*
     &    zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s156**(-1)
     &     )
      amp78(2,2,1) = amp78(2,2,1) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p5)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p6)*
     &    s2347**(-1) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p5)*zba2(p7,p1,
     &    p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p8)*
     &    s2347**(-1)*s156**(-1) + 2.D0*za(p2,p3)*zb(p1,p5)*zb(p1,p7)*
     &    zba2(p4,p2,p3,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p2
     &    ,p3,p4,p6)*s2347**(-1)*s234**(-1) + 2.D0*za(p2,p3)*zb(p1,p5)*
     &    zba2(p4,p2,p3,p8)*zba2(p7,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p2,p3,p4,p8)*s2347**(-1)*s156**(-1)*s234**(-1) )
      amp78(2,2,1) = amp78(2,2,1) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p6)*zb(p1,p7)**2*zba2(p4,p1,p7,p8)*
     &    zba2(p5,p2,p6,p3)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s256**(-1)
     &     )
      amp78(2,2,1) = amp78(2,2,1) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p4)*zba2(p5,p1,p4,p3)*
     &    zba2(p7,p2,p8,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s134**(-1)
     &     )
      amp78(2,2,1) = amp78(2,2,1) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p6)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p5,p2,p6,p7,p3)*
     &    s2567**(-1) + 2.D0*za(p2,p6)*za(p2,p8)*zb(p1,p4)*zba2(p7,p1,
     &    p4,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p5,p2,p6,p7,p8)*
     &    s2567**(-1)*s134**(-1) - 2.D0*za(p2,p6)*zb(p1,p4)*zb(p1,p7)*
     &    zba2(p5,p2,p6,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p1
     &    ,p4,p8,p3)*s1348**(-1)*s256**(-1) - 2.D0*za(p2,p6)*zb(p1,p4)*
     &    zba2(p5,p2,p6,p8)*zba2(p7,p1,p4,p3)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p1,p3,p4,p8)*s1348**(-1)*s134**(-1)*s256**(-1) )
      amp78(2,2,1) = amp78(2,2,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p6)*zb(p1
     &    ,p7)**2*zb(p4,p5)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p3)*zb(p1,p7)**2*zba2(p4,p1,p7,p8)*
     &    zba2(p5,p3,p4,p6)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p4)*za(p3,p6)*zb(p1,p7)**2*zb(p4,p5)*zba2(p4,p1,p7,p8)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*za(p2,p6)*zb(p1,p7)**
     &    2*zab2(p3,p5,p6,p4)*zba2(p5,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) )
      amp78(2,2,1) = amp78(2,2,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p6)*
     &    zb(p1,p3)*zb(p4,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*
     &    zba2(p7,p2,p8,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p8)**2*zb(p1,p4)*zba2(p5,p3,p4,p6)*zba2(p7,p2,p8,p3)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p8)**2*zb(p1,p5)*
     &    zab2(p3,p5,p6,p4)*zba2(p7,p2,p8,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      amp78(2,2,1) = amp78(2,2,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p6)*zb(p1,p3)
     &    *zb(p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*zba2(p5,p3
     &    ,p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*
     &    za(p2,p4)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(p2,p6)*
     &    za(p2,p8)*zb(p1,p5)*zb(p1,p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      amp78(2,1,2)= + gamzqe56(jdu,1,2)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(p4,p2,p3,p6)*zba2(
     &    p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s234**(-1) )
      amp78(2,1,2) = amp78(2,1,2) + gamzqe56(jdu,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p5)*zba2(p4,p1,p5,
     &    p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s156**(-1) )
      amp78(2,1,2) = amp78(2,1,2) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*zb(p1,p5)*zba2(p4,p2,p3,p7)*zba2(p8,p1,p5,p6
     &    )*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2,p3,p4,p7)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)*zba2(p4,p2,p3,p7)*
     &    zba2(p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba4(p8,p2
     &    ,p3,p4,p7,p6)*s2347**(-1)*s234**(-1) - 2.D0*zb(p1,p5)*zba2(p4
     &    ,p2,p3,p7)*zba2(p8,p2,p7,p3)*zba2(p8,p1,p5,p6)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*s2347**(-1)*s156**(-1) + 2.D0*zba2(p5,
     &    p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*zba3(p4,p2,p3,p7,p6)*s2347**(-1) )
      amp78(2,1,2) = amp78(2,1,2) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p6)*za(p1,p7)*zb(p1,p8)*zba2(p4,p1,
     &    p8,p7)*zba2(p5,p2,p6,p3)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s256**(-1) )
      amp78(2,1,2) = amp78(2,1,2) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zba2(p5,p1,p4,
     &    p3)*zba2(p8,p2,p7,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s134**(-1) )
      amp78(2,1,2) = amp78(2,1,2) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p6)*zb(p1,p4)*zba2(p5,p2,p6,p7)*zba2(p8,p1,p4
     &    ,p3)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p1,p3,p4,p7)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p5,p2,p6,p7)*zba2(p8,p1,p4,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s1348**(-1)*s256**(-1) - 2.D0*zb(p1,p4)*zba2(p5,p2
     &    ,p6,p7)*zba2(p8,p2,p7,p6)*zba2(p8,p1,p4,p3)*iza(p7,p8)*izb(p2
     &    ,p7)*izb(p7,p8)*s2567**(-1)*s134**(-1) + 2.D0*zba2(p4,p1,p8,
     &    p7)*zba2(p8,p2,p7,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7
     &    ,p8)*zba3(p5,p2,p6,p7,p3)*s2567**(-1) )
      amp78(2,1,2) = amp78(2,1,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p3
     &    ,p6)*zb(p1,p8)*zb(p4,p5)*zba2(p3,p1,p8,p7)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8) - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(p4,
     &    p1,p8,p7)*zba2(p5,p3,p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p4)*za(p1,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*
     &    zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(
     &    p2,p6)*za(p1,p7)*zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p5,p1,p8,p7
     &    )*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      amp78(2,1,2) = amp78(2,1,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2
     &    ,p8)*zb(p1,p3)*zb(p4,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*zb(p1,p4)
     &    *zb(p4,p5)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     - 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zba2(p5,p3,p4,p6)*zba2(
     &    p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*za(p2,p7
     &    )*zb(p2,p8)*zb(p1,p5)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p6)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      amp78(2,1,2) = amp78(2,1,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p3,p6)*zb(p4,p5)*zba2(p3,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     + 2.D0*za(p3,p6)*zb(p4,p5)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - 2.D0*zab2(
     &    p3,p5,p6,p4)*zba2(p5,p1,p8,p7)*zba2(p8,p2,p7,p6)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*zba2(p4,p1,p8,p7)*
     &    zba2(p5,p3,p4,p6)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8) )
      amp78(2,2,2)= + gamzqe56(jdu,1,2)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)*zba2(p4,p2,p3,p6)*iza(p2,p8)*iza(
     &    p7,p8)*zba3(p5,p1,p7,p8,p2)*s234**(-1) - 2.D0*za(p2,p3)*zb(p1
     &    ,p8)*zba2(p4,p2,p3,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p5,p1,p7,p8
     &    ,p2)*s234**(-1) + 2.D0*za(p2,p3)*zba2(p4,p2,p3,p6)*zba2(p7,p1
     &    ,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p5,p1,p7,p8,p2)
     &    *s234**(-1) )
      amp78(2,2,2) = amp78(2,2,2) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p3)**2*zb(p1,p5)*zb(p3,p4)*zba2(p7,p1,p5,p6)*
     &    iza(p2,p8)*iza(p7,p8)*s156**(-1)*s234**(-1) + 2.D0*za(p2,p3)
     &    **2*zb(p1,p5)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p5,p6)*
     &    iza(p2,p7)*iza(p2,p8)*s2347**(-1)*s156**(-1)*s234**(-1) - 2.D0
     &    *za(p2,p3)**2*zb(p1,p5)*zb(p3,p4)*zba2(p8,p1,p5,p6)*iza(p2,p7
     &    )*iza(p7,p8)*s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)**2*zb(p3,
     &    p4)*zba2(p5,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(
     &    p7,p2,p3,p4,p6)*s2347**(-1)*s234**(-1) )
      amp78(2,2,2) = amp78(2,2,2) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p6)*zb(p1,p7)*zba2(p5,p2,p6,p3)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s256**(-1) - 2.D0*
     &    za(p2,p6)*zb(p1,p8)*zba2(p5,p2,p6,p3)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) + 2.D0*za(p2,p6)*zba2(p5,p2,
     &    p6,p3)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) )
      amp78(2,2,2) = amp78(2,2,2) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p6)**2*zb(p1,p4)*zb(p5,p6)*zba2(p7,p1,p4,p3)*
     &    iza(p2,p8)*iza(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*za(p2,p6)
     &    **2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*iza(p2,p7)*iza(p2,
     &    p8)*zba4(p7,p1,p3,p4,p8,p2)*s1348**(-1)*s134**(-1)*s256**(-1)
     &     + 2.D0*za(p2,p6)**2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*
     &    iza(p2,p7)*iza(p7,p8)*s134**(-1)*s256**(-1) - 2.D0*za(p2,p6)
     &    **2*zb(p5,p6)*zba2(p4,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,
     &    p8)*zba3(p7,p1,p4,p8,p3)*s1348**(-1)*s256**(-1) )
      amp78(2,2,2) = amp78(2,2,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p6)*zb(p1
     &    ,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,p8)*zba3(p3,p1,p7,p8,p2) - 2.
     &    D0*za(p2,p3)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,
     &    p8)*zba3(p3,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p3,p6)*zb(p4,p5)
     &    *zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p3,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p3)*zb(p1,p7)*zba2(p5,p3,p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*
     &    zb(p1,p8)*zba2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,
     &    p7,p8,p2) + 2.D0*za(p2,p3)*zba2(p5,p3,p4,p6)*zba2(p7,p1,p8,p2
     &    )*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0
     &    *za(p2,p4)*za(p3,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,p8
     &    )*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p4)*za(p3,p6)*zb(p1,p8)*
     &    zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p4)*za(p3,p6)*zb(p4,p5)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) )
      amp78(2,2,2) = amp78(2,2,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * ( 2.D0*za(p2,p6)*zb(p1,p7)*zab2(p3,
     &    p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p5,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p6)*zb(p1,p8)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p5,p1,p7,p8,p2) - 2.D0*za(p2,p6)*zab2(p3,p5,p6,p4)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p5,p1,p7,
     &    p8,p2) )
C---  bmp78(h56,h7,h8)
      bmp78(1,1,1)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s278**(-1) * ( za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*izb(p1,p8
     &    )*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*s345**(-1) - za(p2,p7)*za(
     &    p3,p5)*zb(p1,p6)*zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,
     &    p7,p8,p5)*s345**(-1) + za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p3,p4
     &    )*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*s345**(-1) - za(
     &    p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s345**(-1) + za(p3,p5)*zb(p1,p6)*zb(p3,
     &    p4)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(
     &    p1,p2,p7,p8,p3)*s345**(-1) - za(p3,p5)*zb(p1,p6)*zb(p4,p5)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p5)*s345**(-1) )
      bmp78(1,1,1) = bmp78(1,1,1) + gamzee56(1,1)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * ( 2.D0*za(p2,p7)*za(p4,p5)*zb(p1,p4)*zb(
     &    p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p8)*za(p4,p5)*zb(p1,p4)*zb(
     &    p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p3) + 2.D0*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) )
      bmp78(1,1,1) = bmp78(1,1,1) + gamzne56(1,1)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*
     &    za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*izb(p1,p8)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p5) - 2.D0*za(p2,p8)*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*
     &    za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p5) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5) )
      bmp78(1,2,1)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s178**(-1) * (  - za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*zb(p3,p4)*
     &    zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s345**(-1)
     &     + za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p5)*zba2(p6,p1,p7,
     &    p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s345**(-1) )
      bmp78(1,2,1) = bmp78(1,2,1) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1)*s278**(-1) * (  - za(p2,p8)**2*za(p3,p5)*zb(p1,p6)*zb(
     &    p3,p4)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) + za(p2,p8)**2*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*zba2(
     &    p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s345**(-1) )
      bmp78(1,2,1) = bmp78(1,2,1) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1) * ( za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p1,p7)*
     &    zb(p3,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) - za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p1,p7
     &    )*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      bmp78(1,2,1) = bmp78(1,2,1) + gamzee56(1,1)*propw3456**(-1)*
     & s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)*zb(p1,p7)
     &    **2*zb(p4,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)**2*zb(p4,p6)*zba2(p6
     &    ,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(1,2,1) = bmp78(1,2,1) + gamzee56(1,1)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)**2*za(p4,p5)*zb(p1,
     &    p4)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p8)**2*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*zba2(p7
     &    ,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp78(1,2,1) = bmp78(1,2,1) + gamzee56(1,1)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p4,p5)*zb(p1,p4)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p3)*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p1,p7)*zb(
     &    p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(1,2,1) = bmp78(1,2,1) + gamzne56(1,1)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p5,p6)*zba2(p4,p1
     &    ,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(1,2,1) = bmp78(1,2,1) + gamzne56(1,1)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)
     &     + 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*zba2(p7,p2
     &    ,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp78(1,2,1) = bmp78(1,2,1) + gamzne56(1,1)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*
     &    zb(p1,p7)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) - 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(1,1,2)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s178**(-1) * (  - za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,
     &    p4)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) + za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p5
     &    )*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) )
      bmp78(1,1,2) = bmp78(1,1,2) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1)*s278**(-1) * (  - za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,
     &    p6)*zb(p3,p4)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s345**(-1) + za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p6)*zb(
     &    p4,p5)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      bmp78(1,1,2) = bmp78(1,1,2) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1) * ( za(p3,p5)*zb(p3,p4)*zba2(p6,p1,p8,p7)*zba2(p8,p2,
     &    p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*s345**(-1)
     &     - za(p3,p5)*zb(p4,p5)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p5)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*s345**(-1) )
      bmp78(1,1,2) = bmp78(1,1,2) + gamzee56(1,1)*propw3456**(-1)*
     & s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8) + 2.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*zb(p1,p8)*zb(
     &    p4,p6)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp78(1,1,2) = bmp78(1,1,2) + gamzee56(1,1)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p5)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) + 2.D0*za(p2,p7)*za(p5,p6)*zb(p2,p8)*zb(p1,p6)*zb(
     &    p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(1,1,2) = bmp78(1,1,2) + gamzee56(1,1)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p4,p5)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(
     &    p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - 2.D
     &    0*za(p5,p6)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(1,1,2) = bmp78(1,1,2) + gamzne56(1,1)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(
     &    p1,p8)*zb(p3,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8) + 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p5,
     &    p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp78(1,1,2) = bmp78(1,1,2) + gamzne56(1,1)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(
     &    p1,p4)*zb(p3,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8) + 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p4)*zb(p5,
     &    p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(1,1,2) = bmp78(1,1,2) + gamzne56(1,1)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p3,p5)*zb(p3,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     - 2.D0*za(p3,p5)*zb(p5,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(1,2,2)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s178**(-1) * (  - za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p4)*iza(p2
     &    ,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) - za(p2,p3)*
     &    za(p3,p5)*zb(p1,p8)*zb(p3,p4)*iza(p2,p7)*iza(p7,p8)*zba3(p6,
     &    p1,p7,p8,p2)*s345**(-1) + za(p2,p3)*za(p3,p5)*zb(p3,p4)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,
     &    p8,p2)*s345**(-1) + za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) + za(p2
     &    ,p5)*za(p3,p5)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p6,p1,p7,p8,p2)*s345**(-1) - za(p2,p5)*za(p3,p5)*zb(p4,
     &    p5)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(
     &    p6,p1,p7,p8,p2)*s345**(-1) )
      bmp78(1,2,2) = bmp78(1,2,2) + gamzee56(1,1)*propw3456**(-1)*
     & s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)*zb(p1,p7)*
     &    zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*
     &    za(p2,p3)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p4,p5)*zb(p4,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1
     &    ,p7,p8,p2) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) + 2.D0*za(p2,p3)*
     &    za(p5,p6)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p6,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p4,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)
     &     )
      bmp78(1,2,2) = bmp78(1,2,2) + gamzne56(1,1)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(
     &    p3,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(
     &    p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p3,p5)*zb(p3,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1
     &    ,p7,p8,p2) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p5,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p8)*zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p5)*za(p3,p5)*zb(p5,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)
     &     )
      bmp78(2,1,1)= + gamzee56(1,2)*propw3456**(-1)*s278**(-1)*
     & s456**(-1) * ( 2.D0*za(p2,p7)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(
     &    p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p7)*za(p5
     &    ,p6)*zb(p1,p5)*zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,
     &    p8,p3) + 2.D0*za(p2,p8)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,
     &    p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p8)*za(p5,p6
     &    )*zb(p1,p5)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,
     &    p3) + 2.D0*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*zba2(p1,p2,p7,p8)*
     &    izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*
     &    za(p5,p6)*zb(p1,p5)*zb(p4,p5)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) )
      bmp78(2,1,1) = bmp78(2,1,1) + gamzne56(1,2)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p4)*
     &    zb(p3,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*
     &    za(p2,p7)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*izb(p1,p8)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p6) - 2.D0*za(p2,p8)*za(p3,p6)*zb(p1,p4)*
     &    zb(p3,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*
     &    za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p6) - 2.D0*za(p3,p6)*zb(p1,p4)*zb(p3,p5)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) + 2.D0*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p6) )
      bmp78(2,2,1)= + gamzee56(1,2)*propw3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p6)*zb(p1,p7)**2*zb(p4,p5
     &    )*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*
     &    za(p2,p3)*za(p5,p6)*zb(p1,p7)**2*zb(p4,p5)*zba2(p5,p1,p7,p8)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(2,2,1) = bmp78(2,2,1) + gamzee56(1,2)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)**2*za(p4,p6)*zb(p1,
     &    p4)*zb(p4,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8) - 2.D0*za(p2,p8)**2*za(p5,p6)*zb(p1,p5)*zb(p4,p5)*zba2(p7
     &    ,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp78(2,2,1) = bmp78(2,2,1) + gamzee56(1,2)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p4,p6)*zb(p1,p4)*zb(
     &    p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     + 2.D0*za(p2,p3)*za(p2,p8)*za(p5,p6)*zb(p1,p5)*zb(p1,p7)*zb(
     &    p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(2,2,1) = bmp78(2,2,1) + gamzne56(1,2)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p6)*zb(p1,p7)**2*
     &    zb(p3,p5)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p6)*za(p3,p6)*zb(p1,p7)**2*zb(p5,p6)*zba2(p4,p1
     &    ,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(2,2,1) = bmp78(2,2,1) + gamzne56(1,2)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*
     &    zb(p3,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*zba2(p7,p2
     &    ,p8,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp78(2,2,1) = bmp78(2,2,1) + gamzne56(1,2)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*
     &    zb(p1,p7)*zb(p3,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p6)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp78(2,1,2)= + gamzee56(1,2)*propw3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p4,p6)*zb(p1,p8)*
     &    zb(p4,p5)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*zb(p1,p8)*zb(p4,p5)*
     &    zba2(p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp78(2,1,2) = bmp78(2,1,2) + gamzee56(1,2)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p6)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) - 2.D0*za(p2,p7)*za(p5,p6)*zb(p2,p8)*zb(p1,p5)*zb(
     &    p4,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(2,1,2) = bmp78(2,1,2) + gamzee56(1,2)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p4,p6)*zb(p4,p5)*zba2(p4,p1,p8,p7)*zba2(
     &    p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D
     &    0*za(p5,p6)*zb(p4,p5)*zba2(p5,p1,p8,p7)*zba2(p8,p2,p7,p3)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(2,1,2) = bmp78(2,1,2) + gamzne56(1,2)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p6)*zb(
     &    p1,p8)*zb(p3,p5)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8) - 2.D0*za(p2,p6)*za(p1,p7)*za(p3,p6)*zb(p1,p8)*zb(p5,
     &    p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp78(2,1,2) = bmp78(2,1,2) + gamzne56(1,2)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*zb(
     &    p1,p4)*zb(p3,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*zb(p1,p4)*zb(p5,
     &    p6)*zba2(p8,p2,p7,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(2,1,2) = bmp78(2,1,2) + gamzne56(1,2)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p3,p6)*zb(p3,p5)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     + 2.D0*za(p3,p6)*zb(p5,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp78(2,2,2)= + gamzee56(1,2)*propw3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p6)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*
     &    za(p4,p6)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p4,p6)*zb(p4,p5)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)
     &     - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*
     &    iza(p7,p8)*zba3(p5,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p5,p6)*
     &    zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p5,p1,p7,p8,p2
     &    ) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p4,p5)*zba2(p7,p1,p8,p2)*iza(
     &    p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p5,p1,p7,p8,p2) )
      bmp78(2,2,2) = bmp78(2,2,2) + gamzne56(1,2)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p6)*zb(p1,p7)*zb(
     &    p3,p5)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(
     &    p2,p3)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p3,p6)*zb(p3,p5)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1
     &    ,p7,p8,p2) - 2.D0*za(p2,p6)*za(p3,p6)*zb(p1,p7)*zb(p5,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p6)*
     &    za(p3,p6)*zb(p1,p8)*zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2) + 2.D0*za(p2,p6)*za(p3,p6)*zb(p5,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)
     &     )
      elseif (j78==2) then
      amp87(1,1,1)= + gamzqe56(jdu,1,1)*propw34**(-1)*s278**(-1) * ( 
     &     - 2.D0*za(p2,p7)*zb(p1,p6)*zba2(p4,p1,p6,p5)*izb(p1,p8)*izb(
     &    p7,p8)*zba3(p1,p2,p7,p8,p3)*s156**(-1) - 2.D0*za(p2,p8)*zb(p1
     &    ,p6)*zba2(p4,p1,p6,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8
     &    ,p3)*s156**(-1) - 2.D0*zb(p1,p6)*zba2(p1,p2,p7,p8)*zba2(p4,p1
     &    ,p6,p5)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &    *s156**(-1) )
      amp87(1,1,1) = amp87(1,1,1) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*zba2(p4,p2,p3,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba4(p1,p2,p3,p4,p7,p8)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p6)**2
     &    *zba2(p4,p2,p3,p7)*izb(p1,p8)*izb(p7,p8)*s156**(-1)*
     &    s234**(-1) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p6)**2*zba2(p4,p2
     &    ,p3,p8)*izb(p1,p7)*izb(p7,p8)*s156**(-1)*s234**(-1) + 2.D0*
     &    za(p5,p6)*zb(p1,p6)**2*zba2(p1,p2,p7,p3)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*s2347**(-1)*s156**(-1) )
      amp87(1,1,1) = amp87(1,1,1) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s278**(-1) * (  - 2.D0*za(p2,p7)*zb(p1,p4)*zba2(p6,p1,p4,p3)*
     &    izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p5)*s134**(-1) - 2.D0*
     &    za(p2,p8)*zb(p1,p4)*zba2(p6,p1,p4,p3)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s134**(-1) - 2.D0*zb(p1,p4)*zba2(p1,p2,
     &    p7,p8)*zba2(p6,p1,p4,p3)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s134**(-1) )
      amp87(1,1,1) = amp87(1,1,1) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p4)**2*zba2(p1,p3,p4,p8)*
     &    zba2(p6,p2,p5,p7)*izb(p1,p7)*izb(p1,p8)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p4)**2
     &    *zba2(p6,p2,p5,p7)*izb(p1,p8)*izb(p7,p8)*s134**(-1)*
     &    s256**(-1) + 2.D0*za(p2,p5)*za(p3,p4)*zb(p1,p4)**2*zba2(p6,p2
     &    ,p5,p8)*izb(p1,p7)*izb(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*
     &    za(p3,p4)*zb(p1,p4)**2*zba2(p1,p2,p7,p5)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p6,p2,p5,p7,p8)*s2567**(-1)*s134**(-1) )
      amp87(1,1,1) = amp87(1,1,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p3
     &    )*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8
     &    )*zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p7)*zb(p1,p4)*zab2(p5,p3,
     &    p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p7)*zb(p1,p6)*zab2(p3,p5,p6,p4)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5) + 2.D0*za(p2,p8)*za(p3,p5)*zb(p1,p3)*zb(
     &    p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(
     &    p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(p2,
     &    p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1
     &    ,p2,p7,p8,p5) + 2.D0*za(p3,p5)*zb(p1,p3)*zb(p4,p6)*zba2(p1,p2
     &    ,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &     + 2.D0*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(
     &    p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4) )
      amp87(1,1,1) = amp87(1,1,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*zb(p1,p4)*zab2(p5,p3,p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p1,p2,p7,
     &    p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5) )
      amp87(1,1,2)= + gamzqe56(jdu,1,1)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)**2*zba2(p4,p2,p3,p5)*zba2(p6,p1,
     &    p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s234**(-1) )
      amp87(1,1,2) = amp87(1,1,2) + gamzqe56(jdu,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p6)*zba2(p4,p1,p6,p5)*
     &    zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s156**(-1)
     &     )
      amp87(1,1,2) = amp87(1,1,2) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p6)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*
     &    s2347**(-1) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p6)*zba2(p7,p1,
     &    p6,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p8)*
     &    s2347**(-1)*s156**(-1) + 2.D0*za(p2,p3)*zb(p1,p6)*zb(p1,p7)*
     &    zba2(p4,p2,p3,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p2
     &    ,p3,p4,p5)*s2347**(-1)*s234**(-1) + 2.D0*za(p2,p3)*zb(p1,p6)*
     &    zba2(p4,p2,p3,p8)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p2,p3,p4,p8)*s2347**(-1)*s156**(-1)*s234**(-1) )
      amp87(1,1,2) = amp87(1,1,2) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*zb(p1,p7)**2*zba2(p4,p1,p7,p8)*
     &    zba2(p6,p2,p5,p3)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s256**(-1)
     &     )
      amp87(1,1,2) = amp87(1,1,2) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p4)*zba2(p6,p1,p4,p3)*
     &    zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s134**(-1)
     &     )
      amp87(1,1,2) = amp87(1,1,2) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p5)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p6,p2,p5,p7,p3)*
     &    s2567**(-1) + 2.D0*za(p2,p5)*za(p2,p8)*zb(p1,p4)*zba2(p7,p1,
     &    p4,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p6,p2,p5,p7,p8)*
     &    s2567**(-1)*s134**(-1) - 2.D0*za(p2,p5)*zb(p1,p4)*zb(p1,p7)*
     &    zba2(p6,p2,p5,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p1
     &    ,p4,p8,p3)*s1348**(-1)*s256**(-1) - 2.D0*za(p2,p5)*zb(p1,p4)*
     &    zba2(p6,p2,p5,p8)*zba2(p7,p1,p4,p3)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p1,p3,p4,p8)*s1348**(-1)*s134**(-1)*s256**(-1) )
      amp87(1,1,2) = amp87(1,1,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1
     &    ,p7)**2*zb(p4,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p3)*zb(p1,p7)**2*zab2(p5,p3,p4,p6)*
     &    zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p4)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,p8)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*za(p2,p5)*zb(p1,p7)**
     &    2*zab2(p3,p5,p6,p4)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) )
      amp87(1,1,2) = amp87(1,1,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p5)*
     &    zb(p1,p3)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p4,p6)*
     &    zba2(p7,p2,p8,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p8)**2*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p7,p2,p8,p3)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p8)**2*zb(p1,p6)*
     &    zab2(p3,p5,p6,p4)*zba2(p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      amp87(1,1,2) = amp87(1,1,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p3)
     &    *zb(p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*zab2(p5,p3
     &    ,p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*
     &    za(p2,p4)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(p2,p5)*
     &    za(p2,p8)*zb(p1,p6)*zb(p1,p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      amp87(1,2,1)= + gamzqe56(jdu,1,1)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(p4,p2,p3,p5)*zba2(
     &    p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s234**(-1) )
      amp87(1,2,1) = amp87(1,2,1) + gamzqe56(jdu,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p6)*zba2(p4,p1,p6,
     &    p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s156**(-1) )
      amp87(1,2,1) = amp87(1,2,1) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(p8,p1,p6,p5
     &    )*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2,p3,p4,p7)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)*zba2(p4,p2,p3,p7)*
     &    zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba4(p8,p2
     &    ,p3,p4,p7,p5)*s2347**(-1)*s234**(-1) - 2.D0*zb(p1,p6)*zba2(p4
     &    ,p2,p3,p7)*zba2(p8,p2,p7,p3)*zba2(p8,p1,p6,p5)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*s2347**(-1)*s156**(-1) + 2.D0*zba2(p6,
     &    p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*s2347**(-1) )
      amp87(1,2,1) = amp87(1,2,1) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*za(p1,p7)*zb(p1,p8)*zba2(p4,p1,
     &    p8,p7)*zba2(p6,p2,p5,p3)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s256**(-1) )
      amp87(1,2,1) = amp87(1,2,1) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zba2(p6,p1,p4,
     &    p3)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s134**(-1) )
      amp87(1,2,1) = amp87(1,2,1) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p5)*zb(p1,p4)*zba2(p6,p2,p5,p7)*zba2(p8,p1,p4
     &    ,p3)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p1,p3,p4,p7)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p5)*zba2(p4,p1,p8,p7)*
     &    zba2(p6,p2,p5,p7)*zba2(p8,p1,p4,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s1348**(-1)*s256**(-1) - 2.D0*zb(p1,p4)*zba2(p6,p2
     &    ,p5,p7)*zba2(p8,p2,p7,p5)*zba2(p8,p1,p4,p3)*iza(p7,p8)*izb(p2
     &    ,p7)*izb(p7,p8)*s2567**(-1)*s134**(-1) + 2.D0*zba2(p4,p1,p8,
     &    p7)*zba2(p8,p2,p7,p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7
     &    ,p8)*zba3(p6,p2,p5,p7,p3)*s2567**(-1) )
      amp87(1,2,1) = amp87(1,2,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p3
     &    ,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p3,p1,p8,p7)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8) - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zab2(p5,
     &    p3,p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p4)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(
     &    p2,p5)*za(p1,p7)*zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p8,p7
     &    )*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      amp87(1,2,1) = amp87(1,2,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p2
     &    ,p8)*zb(p1,p3)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8) - 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p4)
     &    *zb(p4,p6)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     - 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(
     &    p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*za(p2,p7
     &    )*zb(p2,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p5)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      amp87(1,2,1) = amp87(1,2,1) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p3,p5)*zb(p4,p6)*zba2(p3,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     + 2.D0*za(p3,p5)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - 2.D0*zab2(
     &    p3,p5,p6,p4)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*zab2(p5,p3,p4,p6)*
     &    zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8) )
      amp87(1,2,2)= + gamzqe56(jdu,1,1)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)*zba2(p4,p2,p3,p5)*iza(p2,p8)*iza(
     &    p7,p8)*zba3(p6,p1,p7,p8,p2)*s234**(-1) - 2.D0*za(p2,p3)*zb(p1
     &    ,p8)*zba2(p4,p2,p3,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p6,p1,p7,p8
     &    ,p2)*s234**(-1) + 2.D0*za(p2,p3)*zba2(p4,p2,p3,p5)*zba2(p7,p1
     &    ,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)
     &    *s234**(-1) )
      amp87(1,2,2) = amp87(1,2,2) + gamzqe56(jdu,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p1,p6,p5)*
     &    iza(p2,p8)*iza(p7,p8)*s156**(-1)*s234**(-1) + 2.D0*za(p2,p3)
     &    **2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*iza(p2,p8)*s2347**(-1)*s156**(-1)*s234**(-1) - 2.D0
     &    *za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zba2(p8,p1,p6,p5)*iza(p2,p7
     &    )*iza(p7,p8)*s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)**2*zb(p3,
     &    p4)*zba2(p6,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(
     &    p7,p2,p3,p4,p5)*s2347**(-1)*s234**(-1) )
      amp87(1,2,2) = amp87(1,2,2) + gamzqe56(jdd,1,1)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p5)*zb(p1,p7)*zba2(p6,p2,p5,p3)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s256**(-1) - 2.D0*
     &    za(p2,p5)*zb(p1,p8)*zba2(p6,p2,p5,p3)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) + 2.D0*za(p2,p5)*zba2(p6,p2,
     &    p5,p3)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) )
      amp87(1,2,2) = amp87(1,2,2) + gamzqe56(jdd,1,1)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p5)**2*zb(p1,p4)*zb(p5,p6)*zba2(p7,p1,p4,p3)*
     &    iza(p2,p8)*iza(p7,p8)*s134**(-1)*s256**(-1) - 2.D0*za(p2,p5)
     &    **2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*iza(p2,p7)*iza(p2,
     &    p8)*zba4(p7,p1,p3,p4,p8,p2)*s1348**(-1)*s134**(-1)*s256**(-1)
     &     - 2.D0*za(p2,p5)**2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*
     &    iza(p2,p7)*iza(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*za(p2,p5)
     &    **2*zb(p5,p6)*zba2(p4,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,
     &    p8)*zba3(p7,p1,p4,p8,p3)*s1348**(-1)*s256**(-1) )
      amp87(1,2,2) = amp87(1,2,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p5)*zb(p1
     &    ,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p3,p1,p7,p8,p2) - 2.
     &    D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,
     &    p8)*zba3(p3,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p3,p5)*zb(p4,p6)
     &    *zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p3,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p3)*zb(p1,p7)*zab2(p5,p3,p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*
     &    zb(p1,p8)*zab2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,
     &    p7,p8,p2) + 2.D0*za(p2,p3)*zab2(p5,p3,p4,p6)*zba2(p7,p1,p8,p2
     &    )*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0
     &    *za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p4,p6)*iza(p2,p8)*iza(p7,p8
     &    )*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p4)*za(p3,p5)*zb(p1,p8)*
     &    zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p4)*za(p3,p5)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) )
      amp87(1,2,2) = amp87(1,2,2) + gamzew56(1)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * ( 2.D0*za(p2,p5)*zb(p1,p7)*zab2(p3,
     &    p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p5)*zb(p1,p8)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p6,p1,p7,p8,p2) - 2.D0*za(p2,p5)*zab2(p3,p5,p6,p4)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,
     &    p8,p2) )
      amp87(2,1,1)= + gamzqe56(jdu,1,2)*propw34**(-1)*s278**(-1) * ( 
     &     - 2.D0*za(p2,p7)*zb(p1,p5)*zba2(p4,p1,p5,p6)*izb(p1,p8)*izb(
     &    p7,p8)*zba3(p1,p2,p7,p8,p3)*s156**(-1) - 2.D0*za(p2,p8)*zb(p1
     &    ,p5)*zba2(p4,p1,p5,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8
     &    ,p3)*s156**(-1) - 2.D0*zb(p1,p5)*zba2(p1,p2,p7,p8)*zba2(p4,p1
     &    ,p5,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &    *s156**(-1) )
      amp87(2,1,1) = amp87(2,1,1) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p5)**2*zba2(p4,p2,p3,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba4(p1,p2,p3,p4,p7,p8)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p5)**2
     &    *zba2(p4,p2,p3,p7)*izb(p1,p8)*izb(p7,p8)*s156**(-1)*
     &    s234**(-1) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p5)**2*zba2(p4,p2
     &    ,p3,p8)*izb(p1,p7)*izb(p7,p8)*s156**(-1)*s234**(-1) - 2.D0*
     &    za(p5,p6)*zb(p1,p5)**2*zba2(p1,p2,p7,p3)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*s2347**(-1)*s156**(-1) )
      amp87(2,1,1) = amp87(2,1,1) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s278**(-1) * (  - 2.D0*za(p2,p7)*zb(p1,p4)*zba2(p5,p1,p4,p3)*
     &    izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p6)*s134**(-1) - 2.D0*
     &    za(p2,p8)*zb(p1,p4)*zba2(p5,p1,p4,p3)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p6)*s134**(-1) - 2.D0*zb(p1,p4)*zba2(p1,p2,
     &    p7,p8)*zba2(p5,p1,p4,p3)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*
     &    zba3(p1,p2,p7,p8,p6)*s134**(-1) )
      amp87(2,1,1) = amp87(2,1,1) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*zba2(p1,p3,p4,p8)*
     &    zba2(p5,p2,p6,p7)*izb(p1,p7)*izb(p1,p8)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p4)**2
     &    *zba2(p5,p2,p6,p7)*izb(p1,p8)*izb(p7,p8)*s134**(-1)*
     &    s256**(-1) + 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*zba2(p5,p2
     &    ,p6,p8)*izb(p1,p7)*izb(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*
     &    za(p3,p4)*zb(p1,p4)**2*zba2(p1,p2,p7,p6)*izb(p2,p7)*izb(p1,p7
     &    )*izb(p1,p8)*zba3(p5,p2,p6,p7,p8)*s2567**(-1)*s134**(-1) )
      amp87(2,1,1) = amp87(2,1,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p3
     &    )*zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0
     &    *za(p2,p7)*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,p8)*izb(p7,p8
     &    )*zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p7)*zb(p1,p4)*zba2(p5,p3,
     &    p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p7)*zb(p1,p5)*zab2(p3,p5,p6,p4)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p6) + 2.D0*za(p2,p8)*za(p3,p6)*zb(p1,p3)*zb(
     &    p4,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(
     &    p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p4) + 2.D0*za(p2,p8)*zb(p1,p4)*zba2(p5,p3,p4
     &    ,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(p2,
     &    p8)*zb(p1,p5)*zab2(p3,p5,p6,p4)*izb(p1,p7)*izb(p7,p8)*zba3(p1
     &    ,p2,p7,p8,p6) + 2.D0*za(p3,p6)*zb(p1,p3)*zb(p4,p5)*zba2(p1,p2
     &    ,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)
     &     + 2.D0*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*zba2(p1,p2,p7,p8)*izb(
     &    p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4) )
      amp87(2,1,1) = amp87(2,1,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * ( 2.D0*zb(p1,p4)*zba2(p1,p2,p7,p8)*
     &    zba2(p5,p3,p4,p6)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*zb(p1,p5)*zab2(p3,p5,p6,p4)*zba2(p1,p2,p7,
     &    p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p6) )
      amp87(2,1,2)= + gamzqe56(jdu,1,2)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)**2*zba2(p4,p2,p3,p6)*zba2(p5,p1,
     &    p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s234**(-1) )
      amp87(2,1,2) = amp87(2,1,2) + gamzqe56(jdu,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p5)*zba2(p4,p1,p5,p6)*
     &    zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s156**(-1)
     &     )
      amp87(2,1,2) = amp87(2,1,2) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p5)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p6)*
     &    s2347**(-1) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p5)*zba2(p7,p1,
     &    p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p8)*
     &    s2347**(-1)*s156**(-1) + 2.D0*za(p2,p3)*zb(p1,p5)*zb(p1,p7)*
     &    zba2(p4,p2,p3,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p2
     &    ,p3,p4,p6)*s2347**(-1)*s234**(-1) + 2.D0*za(p2,p3)*zb(p1,p5)*
     &    zba2(p4,p2,p3,p8)*zba2(p7,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p2,p3,p4,p8)*s2347**(-1)*s156**(-1)*s234**(-1) )
      amp87(2,1,2) = amp87(2,1,2) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p6)*zb(p1,p7)**2*zba2(p4,p1,p7,p8)*
     &    zba2(p5,p2,p6,p3)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s256**(-1)
     &     )
      amp87(2,1,2) = amp87(2,1,2) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p8)**2*zb(p1,p4)*zba2(p5,p1,p4,p3)*
     &    zba2(p7,p2,p8,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s134**(-1)
     &     )
      amp87(2,1,2) = amp87(2,1,2) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p6)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p5,p2,p6,p7,p3)*
     &    s2567**(-1) + 2.D0*za(p2,p6)*za(p2,p8)*zb(p1,p4)*zba2(p7,p1,
     &    p4,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p5,p2,p6,p7,p8)*
     &    s2567**(-1)*s134**(-1) - 2.D0*za(p2,p6)*zb(p1,p4)*zb(p1,p7)*
     &    zba2(p5,p2,p6,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*zba3(p7,p1
     &    ,p4,p8,p3)*s1348**(-1)*s256**(-1) - 2.D0*za(p2,p6)*zb(p1,p4)*
     &    zba2(p5,p2,p6,p8)*zba2(p7,p1,p4,p3)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p7,p1,p3,p4,p8)*s1348**(-1)*s134**(-1)*s256**(-1) )
      amp87(2,1,2) = amp87(2,1,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p6)*zb(p1
     &    ,p7)**2*zb(p4,p5)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p3)*zb(p1,p7)**2*zba2(p4,p1,p7,p8)*
     &    zba2(p5,p3,p4,p6)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p4)*za(p3,p6)*zb(p1,p7)**2*zb(p4,p5)*zba2(p4,p1,p7,p8)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*za(p2,p6)*zb(p1,p7)**
     &    2*zab2(p3,p5,p6,p4)*zba2(p5,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*
     &    izb(p7,p8) )
      amp87(2,1,2) = amp87(2,1,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p8)**2*za(p3,p6)*
     &    zb(p1,p3)*zb(p4,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) - 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*zb(p4,p5)*
     &    zba2(p7,p2,p8,p4)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) - 2.D0*za(
     &    p2,p8)**2*zb(p1,p4)*zba2(p5,p3,p4,p6)*zba2(p7,p2,p8,p3)*iza(
     &    p2,p7)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(p2,p8)**2*zb(p1,p5)*
     &    zab2(p3,p5,p6,p4)*zba2(p7,p2,p8,p6)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8) )
      amp87(2,1,2) = amp87(2,1,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p6)*zb(p1,p3)
     &    *zb(p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)*zba2(p5,p3
     &    ,p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) + 2.D0*
     &    za(p2,p4)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*za(p2,p6)*
     &    za(p2,p8)*zb(p1,p5)*zb(p1,p7)*zab2(p3,p5,p6,p4)*iza(p2,p7)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      amp87(2,2,1)= + gamzqe56(jdu,1,2)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(p4,p2,p3,p6)*zba2(
     &    p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*s234**(-1) )
      amp87(2,2,1) = amp87(2,2,1) + gamzqe56(jdu,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p5)*zba2(p4,p1,p5,
     &    p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s156**(-1) )
      amp87(2,2,1) = amp87(2,2,1) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p3)*zb(p1,p5)*zba2(p4,p2,p3,p7)*zba2(p8,p1,p5,p6
     &    )*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2,p3,p4,p7)*s2347**(-1)*
     &    s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)*zba2(p4,p2,p3,p7)*
     &    zba2(p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*zba4(p8,p2
     &    ,p3,p4,p7,p6)*s2347**(-1)*s234**(-1) - 2.D0*zb(p1,p5)*zba2(p4
     &    ,p2,p3,p7)*zba2(p8,p2,p7,p3)*zba2(p8,p1,p5,p6)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*s2347**(-1)*s156**(-1) + 2.D0*zba2(p5,
     &    p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8)*zba3(p4,p2,p3,p7,p6)*s2347**(-1) )
      amp87(2,2,1) = amp87(2,2,1) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p6)*za(p1,p7)*zb(p1,p8)*zba2(p4,p1,
     &    p8,p7)*zba2(p5,p2,p6,p3)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s256**(-1) )
      amp87(2,2,1) = amp87(2,2,1) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s278**(-1) * ( 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zba2(p5,p1,p4,
     &    p3)*zba2(p8,p2,p7,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s134**(-1) )
      amp87(2,2,1) = amp87(2,2,1) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p6)*zb(p1,p4)*zba2(p5,p2,p6,p7)*zba2(p8,p1,p4
     &    ,p3)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p1,p3,p4,p7)*s1348**(-1)*
     &    s134**(-1)*s256**(-1) + 2.D0*za(p2,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p5,p2,p6,p7)*zba2(p8,p1,p4,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s1348**(-1)*s256**(-1) - 2.D0*zb(p1,p4)*zba2(p5,p2
     &    ,p6,p7)*zba2(p8,p2,p7,p6)*zba2(p8,p1,p4,p3)*iza(p7,p8)*izb(p2
     &    ,p7)*izb(p7,p8)*s2567**(-1)*s134**(-1) + 2.D0*zba2(p4,p1,p8,
     &    p7)*zba2(p8,p2,p7,p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7
     &    ,p8)*zba3(p5,p2,p6,p7,p3)*s2567**(-1) )
      amp87(2,2,1) = amp87(2,2,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p3
     &    ,p6)*zb(p1,p8)*zb(p4,p5)*zba2(p3,p1,p8,p7)*iza(p1,p8)*iza(p7,
     &    p8)*izb(p7,p8) - 2.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)*zba2(p4,
     &    p1,p8,p7)*zba2(p5,p3,p4,p6)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p4)*za(p1,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*
     &    zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) + 2.D0*za(
     &    p2,p6)*za(p1,p7)*zb(p1,p8)*zab2(p3,p5,p6,p4)*zba2(p5,p1,p8,p7
     &    )*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      amp87(2,2,1) = amp87(2,2,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s278**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2
     &    ,p8)*zb(p1,p3)*zb(p4,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,
     &    p7)*izb(p7,p8) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*zb(p1,p4)
     &    *zb(p4,p5)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     - 2.D0*za(p2,p7)*zb(p2,p8)*zb(p1,p4)*zba2(p5,p3,p4,p6)*zba2(
     &    p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*za(p2,p7
     &    )*zb(p2,p8)*zb(p1,p5)*zab2(p3,p5,p6,p4)*zba2(p8,p2,p7,p6)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      amp87(2,2,1) = amp87(2,2,1) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p3,p6)*zb(p4,p5)*zba2(p3,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     + 2.D0*za(p3,p6)*zb(p4,p5)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - 2.D0*zab2(
     &    p3,p5,p6,p4)*zba2(p5,p1,p8,p7)*zba2(p8,p2,p7,p6)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D0*zba2(p4,p1,p8,p7)*
     &    zba2(p5,p3,p4,p6)*zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8) )
      amp87(2,2,2)= + gamzqe56(jdu,1,2)*propw34**(-1)*s178**(-1) * ( 
     &     - 2.D0*za(p2,p3)*zb(p1,p7)*zba2(p4,p2,p3,p6)*iza(p2,p8)*iza(
     &    p7,p8)*zba3(p5,p1,p7,p8,p2)*s234**(-1) - 2.D0*za(p2,p3)*zb(p1
     &    ,p8)*zba2(p4,p2,p3,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p5,p1,p7,p8
     &    ,p2)*s234**(-1) + 2.D0*za(p2,p3)*zba2(p4,p2,p3,p6)*zba2(p7,p1
     &    ,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p5,p1,p7,p8,p2)
     &    *s234**(-1) )
      amp87(2,2,2) = amp87(2,2,2) + gamzqe56(jdu,1,2)*propw34**(-1)
     &  * (  - 2.D0*za(p2,p3)**2*zb(p1,p5)*zb(p3,p4)*zba2(p7,p1,p5,p6)*
     &    iza(p2,p8)*iza(p7,p8)*s156**(-1)*s234**(-1) + 2.D0*za(p2,p3)
     &    **2*zb(p1,p5)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p5,p6)*
     &    iza(p2,p7)*iza(p2,p8)*s2347**(-1)*s156**(-1)*s234**(-1) - 2.D0
     &    *za(p2,p3)**2*zb(p1,p5)*zb(p3,p4)*zba2(p8,p1,p5,p6)*iza(p2,p7
     &    )*iza(p7,p8)*s156**(-1)*s234**(-1) - 2.D0*za(p2,p3)**2*zb(p3,
     &    p4)*zba2(p5,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(
     &    p7,p2,p3,p4,p6)*s2347**(-1)*s234**(-1) )
      amp87(2,2,2) = amp87(2,2,2) + gamzqe56(jdd,1,2)*propw34**(-1)*
     & s178**(-1) * (  - 2.D0*za(p2,p6)*zb(p1,p7)*zba2(p5,p2,p6,p3)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2)*s256**(-1) - 2.D0*
     &    za(p2,p6)*zb(p1,p8)*zba2(p5,p2,p6,p3)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) + 2.D0*za(p2,p6)*zba2(p5,p2,
     &    p6,p3)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s256**(-1) )
      amp87(2,2,2) = amp87(2,2,2) + gamzqe56(jdd,1,2)*propw34**(-1)
     &  * ( 2.D0*za(p2,p6)**2*zb(p1,p4)*zb(p5,p6)*zba2(p7,p1,p4,p3)*
     &    iza(p2,p8)*iza(p7,p8)*s134**(-1)*s256**(-1) + 2.D0*za(p2,p6)
     &    **2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*iza(p2,p7)*iza(p2,
     &    p8)*zba4(p7,p1,p3,p4,p8,p2)*s1348**(-1)*s134**(-1)*s256**(-1)
     &     + 2.D0*za(p2,p6)**2*zb(p1,p4)*zb(p5,p6)*zba2(p8,p1,p4,p3)*
     &    iza(p2,p7)*iza(p7,p8)*s134**(-1)*s256**(-1) - 2.D0*za(p2,p6)
     &    **2*zb(p5,p6)*zba2(p4,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,
     &    p8)*zba3(p7,p1,p4,p8,p3)*s1348**(-1)*s256**(-1) )
      amp87(2,2,2) = amp87(2,2,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * (  - 2.D0*za(p2,p3)*za(p3,p6)*zb(p1
     &    ,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,p8)*zba3(p3,p1,p7,p8,p2) - 2.
     &    D0*za(p2,p3)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,
     &    p8)*zba3(p3,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p3,p6)*zb(p4,p5)
     &    *zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p3,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p3)*zb(p1,p7)*zba2(p5,p3,p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*
     &    zb(p1,p8)*zba2(p5,p3,p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,
     &    p7,p8,p2) + 2.D0*za(p2,p3)*zba2(p5,p3,p4,p6)*zba2(p7,p1,p8,p2
     &    )*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0
     &    *za(p2,p4)*za(p3,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*iza(p7,p8
     &    )*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p4)*za(p3,p6)*zb(p1,p8)*
     &    zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p4)*za(p3,p6)*zb(p4,p5)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2) )
      amp87(2,2,2) = amp87(2,2,2) + gamzew56(2)*propw34**(-1)*
     & propw3456**(-1)*s178**(-1) * ( 2.D0*za(p2,p6)*zb(p1,p7)*zab2(p3,
     &    p5,p6,p4)*iza(p2,p8)*iza(p7,p8)*zba3(p5,p1,p7,p8,p2) + 2.D0*
     &    za(p2,p6)*zb(p1,p8)*zab2(p3,p5,p6,p4)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p5,p1,p7,p8,p2) - 2.D0*za(p2,p6)*zab2(p3,p5,p6,p4)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p5,p1,p7,
     &    p8,p2) )
      bmp87(1,1,1)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s278**(-1) * ( za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*izb(p1,p8
     &    )*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*s345**(-1) - za(p2,p7)*za(
     &    p3,p5)*zb(p1,p6)*zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,
     &    p7,p8,p5)*s345**(-1) + za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p3,p4
     &    )*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3)*s345**(-1) - za(
     &    p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p5)*s345**(-1) + za(p3,p5)*zb(p1,p6)*zb(p3,
     &    p4)*zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(
     &    p1,p2,p7,p8,p3)*s345**(-1) - za(p3,p5)*zb(p1,p6)*zb(p4,p5)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p5)*s345**(-1) )
      bmp87(1,1,1) = bmp87(1,1,1) + gamzee56(1,1)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * ( 2.D0*za(p2,p7)*za(p4,p5)*zb(p1,p4)*zb(
     &    p4,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p8)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p8)*za(p4,p5)*zb(p1,p4)*zb(
     &    p4,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*za(
     &    p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*izb(p7,p8)*
     &    zba3(p1,p2,p7,p8,p3) + 2.D0*za(p4,p5)*zb(p1,p4)*zb(p4,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) )
      bmp87(1,1,1) = bmp87(1,1,1) + gamzne56(1,1)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*
     &    za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*izb(p1,p8)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p5) - 2.D0*za(p2,p8)*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) - 2.D0*
     &    za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p5) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) - 2.D0*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5) )
      bmp87(1,1,2)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s178**(-1) * (  - za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*zb(p3,p4)*
     &    zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s345**(-1)
     &     + za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p4,p5)*zba2(p6,p1,p7,
     &    p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*s345**(-1) )
      bmp87(1,1,2) = bmp87(1,1,2) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1)*s278**(-1) * (  - za(p2,p8)**2*za(p3,p5)*zb(p1,p6)*zb(
     &    p3,p4)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) + za(p2,p8)**2*za(p3,p5)*zb(p1,p6)*zb(p4,p5)*zba2(
     &    p7,p2,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*s345**(-1) )
      bmp87(1,1,2) = bmp87(1,1,2) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1) * ( za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p1,p7)*
     &    zb(p3,p4)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) - za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p6)*zb(p1,p7
     &    )*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    s345**(-1) )
      bmp87(1,1,2) = bmp87(1,1,2) + gamzee56(1,1)*propw3456**(-1)*
     & s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)*zb(p1,p7)
     &    **2*zb(p4,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)**2*zb(p4,p6)*zba2(p6
     &    ,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(1,1,2) = bmp87(1,1,2) + gamzee56(1,1)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)**2*za(p4,p5)*zb(p1,
     &    p4)*zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p8)**2*za(p5,p6)*zb(p1,p6)*zb(p4,p6)*zba2(p7
     &    ,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp87(1,1,2) = bmp87(1,1,2) + gamzee56(1,1)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p4,p5)*zb(p1,p4)*zb(
     &    p1,p7)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p3)*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p1,p7)*zb(
     &    p4,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(1,1,2) = bmp87(1,1,2) + gamzne56(1,1)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)**2*zb(p5,p6)*zba2(p4,p1
     &    ,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(1,1,2) = bmp87(1,1,2) + gamzne56(1,1)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*
     &    zb(p3,p6)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)
     &     + 2.D0*za(p2,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p5,p6)*zba2(p7,p2
     &    ,p8,p5)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp87(1,1,2) = bmp87(1,1,2) + gamzne56(1,1)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*
     &    zb(p1,p7)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) - 2.D0*za(p2,p5)*za(p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(1,2,1)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s178**(-1) * (  - za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,
     &    p4)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) + za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p5
     &    )*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s345**(-1) )
      bmp87(1,2,1) = bmp87(1,2,1) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1)*s278**(-1) * (  - za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,
     &    p6)*zb(p3,p4)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,
     &    p8)*s345**(-1) + za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p6)*zb(
     &    p4,p5)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    s345**(-1) )
      bmp87(1,2,1) = bmp87(1,2,1) + propw34**(-1)*propw3456**(-1)*
     & cxw**(-1) * ( za(p3,p5)*zb(p3,p4)*zba2(p6,p1,p8,p7)*zba2(p8,p2,
     &    p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*s345**(-1)
     &     - za(p3,p5)*zb(p4,p5)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p5)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*s345**(-1) )
      bmp87(1,2,1) = bmp87(1,2,1) + gamzee56(1,1)*propw3456**(-1)*
     & s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p7,p8) + 2.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*zb(p1,p8)*zb(
     &    p4,p6)*zba2(p6,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp87(1,2,1) = bmp87(1,2,1) + gamzee56(1,1)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p5)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) + 2.D0*za(p2,p7)*za(p5,p6)*zb(p2,p8)*zb(p1,p6)*zb(
     &    p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(1,2,1) = bmp87(1,2,1) + gamzee56(1,1)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p4,p5)*zb(p4,p6)*zba2(p4,p1,p8,p7)*zba2(
     &    p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) - 2.D
     &    0*za(p5,p6)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(1,2,1) = bmp87(1,2,1) + gamzne56(1,1)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(
     &    p1,p8)*zb(p3,p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8) + 2.D0*za(p2,p5)*za(p1,p7)*za(p3,p5)*zb(p1,p8)*zb(p5,
     &    p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp87(1,2,1) = bmp87(1,2,1) + gamzne56(1,1)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(
     &    p1,p4)*zb(p3,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8) + 2.D0*za(p2,p7)*za(p3,p5)*zb(p2,p8)*zb(p1,p4)*zb(p5,
     &    p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(1,2,1) = bmp87(1,2,1) + gamzne56(1,1)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p3,p5)*zb(p3,p6)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     - 2.D0*za(p3,p5)*zb(p5,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p5)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(1,2,2)= + propw34**(-1)*propw3456**(-1)*cxw**(-1)*
     & s178**(-1) * (  - za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p4)*iza(p2
     &    ,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) - za(p2,p3)*
     &    za(p3,p5)*zb(p1,p8)*zb(p3,p4)*iza(p2,p7)*iza(p7,p8)*zba3(p6,
     &    p1,p7,p8,p2)*s345**(-1) + za(p2,p3)*za(p3,p5)*zb(p3,p4)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,
     &    p8,p2)*s345**(-1) + za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2)*s345**(-1) + za(p2
     &    ,p5)*za(p3,p5)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p6,p1,p7,p8,p2)*s345**(-1) - za(p2,p5)*za(p3,p5)*zb(p4,
     &    p5)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(
     &    p6,p1,p7,p8,p2)*s345**(-1) )
      bmp87(1,2,2) = bmp87(1,2,2) + gamzee56(1,1)*propw3456**(-1)*
     & s178**(-1)*s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p5)*zb(p1,p7)*
     &    zb(p4,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*
     &    za(p2,p3)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)
     &    *zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p4,p5)*zb(p4,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1
     &    ,p7,p8,p2) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(p4,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p6,p1,p7,p8,p2) + 2.D0*za(p2,p3)*
     &    za(p5,p6)*zb(p1,p8)*zb(p4,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p6,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p5,p6)*zb(p4,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p6,p1,p7,p8,p2)
     &     )
      bmp87(1,2,2) = bmp87(1,2,2) + gamzne56(1,1)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(
     &    p3,p6)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(
     &    p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p3,p5)*zb(p3,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1
     &    ,p7,p8,p2) + 2.D0*za(p2,p5)*za(p3,p5)*zb(p1,p7)*zb(p5,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p8)*zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2) - 2.D0*za(p2,p5)*za(p3,p5)*zb(p5,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)
     &     )
      bmp87(2,1,1)= + gamzee56(1,2)*propw3456**(-1)*s278**(-1)*
     & s456**(-1) * ( 2.D0*za(p2,p7)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(
     &    p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p7)*za(p5
     &    ,p6)*zb(p1,p5)*zb(p4,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,
     &    p8,p3) + 2.D0*za(p2,p8)*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*izb(p1,
     &    p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*za(p2,p8)*za(p5,p6
     &    )*zb(p1,p5)*zb(p4,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,
     &    p3) + 2.D0*za(p4,p6)*zb(p1,p4)*zb(p4,p5)*zba2(p1,p2,p7,p8)*
     &    izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*
     &    za(p5,p6)*zb(p1,p5)*zb(p4,p5)*zba2(p1,p2,p7,p8)*izb(p2,p7)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3) )
      bmp87(2,1,1) = bmp87(2,1,1) + gamzne56(1,2)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * (  - 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p4)*
     &    zb(p3,p5)*izb(p1,p8)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*
     &    za(p2,p7)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*izb(p1,p8)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p6) - 2.D0*za(p2,p8)*za(p3,p6)*zb(p1,p4)*
     &    zb(p3,p5)*izb(p1,p7)*izb(p7,p8)*zba3(p1,p2,p7,p8,p3) + 2.D0*
     &    za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*izb(p7,p8)
     &    *zba3(p1,p2,p7,p8,p6) - 2.D0*za(p3,p6)*zb(p1,p4)*zb(p3,p5)*
     &    zba2(p1,p2,p7,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3) + 2.D0*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*zba2(p1,p2,p7
     &    ,p8)*izb(p2,p7)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p6) )
      bmp87(2,1,2)= + gamzee56(1,2)*propw3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p6)*zb(p1,p7)**2*zb(p4,p5
     &    )*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) - 2.D0*
     &    za(p2,p3)*za(p5,p6)*zb(p1,p7)**2*zb(p4,p5)*zba2(p5,p1,p7,p8)*
     &    iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(2,1,2) = bmp87(2,1,2) + gamzee56(1,2)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p8)**2*za(p4,p6)*zb(p1,
     &    p4)*zb(p4,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,
     &    p8) - 2.D0*za(p2,p8)**2*za(p5,p6)*zb(p1,p5)*zb(p4,p5)*zba2(p7
     &    ,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp87(2,1,2) = bmp87(2,1,2) + gamzee56(1,2)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p2,p3)*za(p2,p8)*za(p4,p6)*zb(p1,p4)*zb(
     &    p1,p7)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     + 2.D0*za(p2,p3)*za(p2,p8)*za(p5,p6)*zb(p1,p5)*zb(p1,p7)*zb(
     &    p4,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(2,1,2) = bmp87(2,1,2) + gamzne56(1,2)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p6)*zb(p1,p7)**2*
     &    zb(p3,p5)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p6)*za(p3,p6)*zb(p1,p7)**2*zb(p5,p6)*zba2(p4,p1
     &    ,p7,p8)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(2,1,2) = bmp87(2,1,2) + gamzne56(1,2)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*
     &    zb(p3,p5)*zba2(p7,p2,p8,p3)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p8)**2*za(p3,p6)*zb(p1,p4)*zb(p5,p6)*zba2(p7,p2
     &    ,p8,p6)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8) )
      bmp87(2,1,2) = bmp87(2,1,2) + gamzne56(1,2)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p2,p3)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*
     &    zb(p1,p7)*zb(p3,p5)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,
     &    p8) + 2.D0*za(p2,p6)*za(p2,p8)*za(p3,p6)*zb(p1,p4)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8) )
      bmp87(2,2,1)= + gamzee56(1,2)*propw3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p3)*za(p1,p7)*za(p4,p6)*zb(p1,p8)*
     &    zb(p4,p5)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8)
     &     - 2.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*zb(p1,p8)*zb(p4,p5)*
     &    zba2(p5,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp87(2,2,1) = bmp87(2,2,1) + gamzee56(1,2)*propw3456**(-1)*
     & s278**(-1)*s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p6)*zb(p2,p8)*
     &    zb(p1,p4)*zb(p4,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*
     &    izb(p7,p8) - 2.D0*za(p2,p7)*za(p5,p6)*zb(p2,p8)*zb(p1,p5)*zb(
     &    p4,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(2,2,1) = bmp87(2,2,1) + gamzee56(1,2)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p4,p6)*zb(p4,p5)*zba2(p4,p1,p8,p7)*zba2(
     &    p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) + 2.D
     &    0*za(p5,p6)*zb(p4,p5)*zba2(p5,p1,p8,p7)*zba2(p8,p2,p7,p3)*
     &    iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(2,2,1) = bmp87(2,2,1) + gamzne56(1,2)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p1,p7)*za(p3,p6)*zb(
     &    p1,p8)*zb(p3,p5)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(
     &    p7,p8) - 2.D0*za(p2,p6)*za(p1,p7)*za(p3,p6)*zb(p1,p8)*zb(p5,
     &    p6)*zba2(p4,p1,p8,p7)*iza(p1,p8)*iza(p7,p8)*izb(p7,p8) )
      bmp87(2,2,1) = bmp87(2,2,1) + gamzne56(1,2)*propw3456**(-1)*
     & s278**(-1)*s356**(-1) * ( 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*zb(
     &    p1,p4)*zb(p3,p5)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p2,p7)*izb(
     &    p7,p8) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p2,p8)*zb(p1,p4)*zb(p5,
     &    p6)*zba2(p8,p2,p7,p6)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(2,2,1) = bmp87(2,2,1) + gamzne56(1,2)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p3,p6)*zb(p3,p5)*zba2(p4,p1,p8,p7)*
     &    zba2(p8,p2,p7,p3)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)
     &     + 2.D0*za(p3,p6)*zb(p5,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,
     &    p6)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8) )
      bmp87(2,2,2)= + gamzee56(1,2)*propw3456**(-1)*s178**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p3)*za(p4,p6)*zb(p1,p7)*zb(p4,p5)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*
     &    za(p4,p6)*zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2) + 2.D0*za(p2,p3)*za(p4,p6)*zb(p4,p5)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)
     &     - 2.D0*za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(p4,p5)*iza(p2,p8)*
     &    iza(p7,p8)*zba3(p5,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p5,p6)*
     &    zb(p1,p8)*zb(p4,p5)*iza(p2,p7)*iza(p7,p8)*zba3(p5,p1,p7,p8,p2
     &    ) + 2.D0*za(p2,p3)*za(p5,p6)*zb(p4,p5)*zba2(p7,p1,p8,p2)*iza(
     &    p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p5,p1,p7,p8,p2) )
      bmp87(2,2,2) = bmp87(2,2,2) + gamzne56(1,2)*propw3456**(-1)*
     & s178**(-1)*s356**(-1) * ( 2.D0*za(p2,p3)*za(p3,p6)*zb(p1,p7)*zb(
     &    p3,p5)*iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) + 2.D0*za(
     &    p2,p3)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*iza(p2,p7)*iza(p7,p8)*
     &    zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p3)*za(p3,p6)*zb(p3,p5)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1
     &    ,p7,p8,p2) - 2.D0*za(p2,p6)*za(p3,p6)*zb(p1,p7)*zb(p5,p6)*
     &    iza(p2,p8)*iza(p7,p8)*zba3(p4,p1,p7,p8,p2) - 2.D0*za(p2,p6)*
     &    za(p3,p6)*zb(p1,p8)*zb(p5,p6)*iza(p2,p7)*iza(p7,p8)*zba3(p4,
     &    p1,p7,p8,p2) + 2.D0*za(p2,p6)*za(p3,p6)*zb(p5,p6)*zba2(p7,p1,
     &    p8,p2)*iza(p2,p7)*iza(p2,p8)*iza(p1,p8)*zba3(p4,p1,p7,p8,p2)
     &     )
      endif
      enddo  !j78
c--- Overall factor
      amp78(:,:,:)=(amp78(:,:,:)+bmp78(:,:,:))/cxw
      amp87(:,:,:)=(amp87(:,:,:)+bmp87(:,:,:))/cxw
      return
      end
