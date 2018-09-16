      subroutine ampmidWZ(p1,p2,p3,p4,p5,p6,p7,p8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cmplxmass.f'
      include 'masses.f'
      include 'runstring.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'WWbits.f'
      integer:: jdu1,jdu2,i1,i2,i3,i4,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,amp(2,2,2),sqwmass,rxw,cosW,
     & propw34,propz56,propz28,propw17,
     & propw3456,proph1347,gamZ56ne(2,2),
     & gamZ56ee(2,2),gamZ28nq(2,2,2),gamZ28eq(2,2,2),gamZ28qW(2,2),
     & gamZ56eW(2),gamZ28eqZ(2,2,2),gamZ28nqZ(2,2,2)
      real(dp):: qn,t3,t4,s34,s56,s17,s28,s137,s147,s167,
     & s238,s248,s258,s268,s456,s345,s356,s3456,s1347,
     & twop17Dp3456,twop34Dp3456,q3,q4,l3,l4
      parameter(qn=zip)
C-----Begin statement functions
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
      include 'cplx.h'
C-----end statement functions

c--- special fix for Madgraph check
      if (index(runstring,'mad') > 0) then
        sqwmass=cplx2(wmass**2,zip)
      else
        sqwmass=cwmass2
      endif

      rxw=sqrt((cone-cxw)/cxw)
      cosW=sqrt((cone-cxw))
      s3456=t4(p3,p4,p5,p6)
      s1347=t4(p1,p3,p4,p7)
      s137=t3(p1,p3,p7)
      s147=t3(p1,p4,p7)
      s167=t3(p1,p6,p7)
      s238=t3(p2,p3,p8)
      s248=t3(p2,p4,p8)
      s258=t3(p2,p5,p8)
      s268=t3(p2,p6,p8)
      s345=t3(p3,p4,p5)
      s356=t3(p3,p5,p6)
      s456=t3(p4,p5,p6)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s17=s(p1,p7)
      s28=s(p2,p8)

      twop17Dp3456=s28-s17-s3456
      twop34Dp3456=s34-s56+s3456

      proph1347=cplx2(s1347-hmass**2,hmass*hwidth)
      propw34=s34-cwmass2
      propz56=s56-czmass2
      propw3456=s3456-cwmass2
      propw17=s17-cwmass2
      propz28=s28-czmass2

      q3=qn
      l3=ln
      q4=qe
      l4=le
      do jdu1=1,2
      gamZ28nq(jdu1,1,1)=q3*Q(jdu1)/s28+l3*L(jdu1)/propZ28
      gamZ28nq(jdu1,1,2)=q3*Q(jdu1)/s28+l3*R(jdu1)/propZ28

      gamZ28eq(jdu1,1,1)=q4*Q(jdu1)/s28+l4*L(jdu1)/propZ28
      gamZ28eq(jdu1,1,2)=q4*Q(jdu1)/s28+l4*R(jdu1)/propZ28

      gamZ28nqZ(jdu1,1,1)=qn*Q(jdu1)/s28+ln*L(jdu1)/propZ28
      gamZ28nqZ(jdu1,1,2)=qn*Q(jdu1)/s28+ln*R(jdu1)/propZ28

      gamZ28eqZ(jdu1,1,1)=qe*Q(jdu1)/s28+le*L(jdu1)/propZ28
      gamZ28eqZ(jdu1,1,2)=qe*Q(jdu1)/s28+le*R(jdu1)/propZ28

      gamZ28qW(jdu1,1)=Q(jdu1)/s28+rxw*L(jdu1)/propZ28
      gamZ28qW(jdu1,2)=Q(jdu1)/s28+rxw*R(jdu1)/propZ28

      gamZ56eW(1)=qe/s56+rxw*le/propZ56
      gamZ56eW(2)=qe/s56+rxw*re/propZ56
      enddo


      gamZ56ne(1,1)=q3*qe/s56+l3*le/propZ56
      gamZ56ne(1,2)=q3*qe/s56+l3*re/propZ56

      gamZ56ee(1,1)=q4*qe/s56+l4*le/propZ56
      gamZ56ee(1,2)=q4*qe/s56+l4*re/propZ56

      do jdu2=1,2
      amp(jdu2,1,1)= + L(jdu2)*cxw**(-2)*Hbit*cosW**(-2)*propw17**(-1)*
     & propw34**(-1)*propZ28**(-1) * (  - 2.D0*za(p3,p7)*za(p5,p8)*zb(
     &    p1,p4)*zb(p2,p6)*le*sqwmass*proph1347**(-1)*propZ56**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*cwmass2**(-1)*
     & cxw**(-2)*Bbit*propw17**(-1)*propw34**(-1)*propw3456**(-1)*
     & s345**(-1) * ( za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p3,p4)*zab2(p3,
     &    p4,p5,p6)*s3456 + za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p3,p4)*
     &    zab2(p3,p4,p5,p6)*twop17Dp3456 - za(p3,p5)*za(p7,p8)*zb(p1,p2
     &    )*zb(p4,p5)*zab2(p5,p3,p4,p6)*s3456 - za(p3,p5)*za(p7,p8)*zb(
     &    p1,p2)*zb(p4,p5)*zab2(p5,p3,p4,p6)*twop17Dp3456 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*cxw**(-2)*Bbit*
     & propw17**(-1)*propw34**(-1)*propw3456**(-1)*s345**(-1) * ( 2.D0*
     &    za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p6)*zba2(p4,p3,p5,p1) -
     &    za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p3,p4)*zab2(p3,p4,p5,p6) +
     &    za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*zab2(p5,p3,p4,p6) - 2.
     &    D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p6,p7)*zba2(p4,p3,p5,p7)
     &     + 2.D0*za(p3,p5)*zb(p1,p6)*zab2(p8,p1,p7,p2)*zba2(p4,p3,p5,
     &    p7) - 2.D0*za(p3,p5)*zb(p2,p6)*zab2(p7,p2,p8,p1)*zba2(p4,p3,
     &    p5,p8) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56eW(1)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*
     & propw3456**(-1) * (  - 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p6)*s3456**2 + 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    twop34Dp3456*s3456 - 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4
     &    ,p6)*twop17Dp3456*s3456 + 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*twop17Dp3456*twop34Dp3456 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56eW(1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*propw3456**(-1) * (
     &     - 4.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*s1347 + 2.D0*
     &    za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*s3456 + 4.D0*za(p3,p5
     &    )*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*s34 + 4.D0*za(p3,p5)*za(p7,p8
     &    )*zb(p1,p2)*zb(p4,p6)*s17 - 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2
     &    )*zb(p4,p6)*twop34Dp3456 + 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)
     &    *zb(p4,p6)*twop17Dp3456 + 4.D0*za(p3,p5)*zb(p4,p6)*zab2(p7,p2
     &    ,p8,p1)*zab2(p8,p3,p4,p2) - 4.D0*za(p3,p5)*zb(p4,p6)*zab2(p7,
     &    p3,p4,p1)*zab2(p8,p1,p7,p2) - 4.D0*za(p3,p7)*zb(p1,p4)*zab2(
     &    p5,p3,p4,p6)*zab2(p8,p1,p7,p2) + 4.D0*za(p3,p8)*zb(p2,p4)*
     &    zab2(p5,p3,p4,p6)*zab2(p7,p2,p8,p1) + 4.D0*za(p5,p7)*zb(p1,p6
     &    )*zab2(p3,p5,p6,p4)*zab2(p8,p1,p7,p2) - 4.D0*za(p5,p8)*zb(p2,
     &    p6)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1) - 4.D0*za(p7,p8)*zb(
     &    p1,p2)*zab2(p3,p1,p7,p4)*zab2(p5,p3,p4,p6) + 4.D0*za(p7,p8)*
     &    zb(p1,p2)*zab2(p3,p5,p6,p4)*zab2(p5,p1,p7,p6) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56eW(1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1) * ( 2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6) - 4.D0*za(p3,p7)*za(p5,p8)*zb(
     &    p1,p4)*zb(p2,p6) + 2.D0*za(p3,p8)*za(p5,p7)*zb(p1,p6)*zb(p2,
     &    p4) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56ne(1,1)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p3,p6)*
     &    zab2(p3,p5,p6,p4)*s3456 - 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p3,p6)*zab2(p3,p5,p6,p4)*twop17Dp3456 - 2.D0*za(p3,p5)*za(
     &    p7,p8)*zb(p1,p2)*zb(p5,p6)*zab2(p5,p3,p6,p4)*s3456 - 2.D0*za(
     &    p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p5,p6)*zab2(p5,p3,p6,p4)*
     &    twop17Dp3456 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s356**(-1) * (  - 4.
     &    D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*zba2(p6,p3,p5,p1)
     &     + 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p3,p6)*zab2(p3,p5,p6
     &    ,p4) + 4.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p7)*zba2(p6,
     &    p3,p5,p7) + 2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p5,p6)*
     &    zab2(p5,p3,p6,p4) - 4.D0*za(p3,p5)*zb(p1,p4)*zab2(p8,p1,p7,p2
     &    )*zba2(p6,p3,p5,p7) + 4.D0*za(p3,p5)*zb(p2,p4)*zab2(p7,p2,p8,
     &    p1)*zba2(p6,p3,p5,p8) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56ee(1,1)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p4,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    zab2(p3,p5,p6,p4)*s3456 + 2.D0*za(p4,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*zab2(p3,p5,p6,p4)*twop17Dp3456 - 2.D0*za(p5,p6)*za(
     &    p7,p8)*zb(p1,p2)*zb(p4,p6)*zab2(p3,p4,p5,p6)*s3456 - 2.D0*za(
     &    p5,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zab2(p3,p4,p5,p6)*
     &    twop17Dp3456 )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28qW(jdu2,1)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s456**(-1) * ( 4.D0
     &    *za(p1,p3)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p1,p4,p6,p5) -
     &    4.D0*za(p3,p7)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p7,p4,p6,p5
     &    ) - 4.D0*za(p3,p7)*zb(p4,p6)*zab2(p8,p1,p7,p2)*zba2(p1,p4,p6,
     &    p5) + 4.D0*za(p3,p8)*zb(p4,p6)*zab2(p7,p2,p8,p1)*zba2(p2,p4,
     &    p6,p5) - 2.D0*za(p4,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zab2(p3
     &    ,p5,p6,p4) + 2.D0*za(p5,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    zab2(p3,p4,p5,p6) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28eq(jdu2,1,1)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s356**(-1) * (  - 4.D0*za(p3,p5)*
     &    zb(p2,p4)*zba2(p1,p2,p4,p8)*zba2(p6,p3,p5,p7)*s248**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28eq(jdu2,1,1)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * ( 4.D0*za(p3,p7)*zb(p4
     &    ,p6)*zba2(p1,p3,p7,p8)*zba2(p2,p4,p6,p5)*s137**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28eq(jdu2,1,1)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1) * (  - 4.D0*za(p3,p7)*zb(p2,p4)*
     &    zba2(p1,p3,p7,p5)*zba2(p6,p2,p4,p8)*s137**(-1)*s248**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28nq(jdu2,1,1)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1)*s356**(-1) * (  - 4.D0*
     &    za(p3,p5)*zb(p1,p4)*zba2(p2,p1,p4,p7)*zba2(p6,p3,p5,p8) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28nq(jdu2,1,1)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1) * (  - 4.D0*za(p3,p8)*
     &    zb(p1,p4)*zba2(p2,p3,p8,p5)*zba2(p6,p1,p4,p7)*s238**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28nq(jdu2,1,1)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * ( 4.D0*za(p3,p8)*zb(p4
     &    ,p6)*zba2(p1,p4,p6,p5)*zba2(p2,p3,p8,p7)*s238**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28eqZ(jdu2,1,1)*cxw**(-2)*
     & Bbit*propw17**(-1)*propw34**(-1)*s258**(-1) * (  - 2.D0*za(p5,p8
     &    )*zb(p1,p6)*zba2(p2,p5,p8,p3)*zba2(p4,p1,p6,p7)*s167**(-1) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28eqZ(jdu2,1,1)*cxw**(-2)*
     & Bbit*propw17**(-1)*propw34**(-1)*s268**(-1)*s345**(-1) * ( 2.D0*
     &    za(p3,p5)*zb(p2,p6)*zba2(p1,p2,p6,p8)*zba2(p4,p3,p5,p7) )
      amp(jdu2,1,1) = amp(jdu2,1,1) + gamZ28nqZ(jdu2,1,1)*cxw**(-2)*
     & Bbit*propw17**(-1)*propw34**(-1)*s345**(-1) * ( 2.D0*za(p3,p5)*
     &    zb(p1,p6)*zba2(p2,p1,p6,p7)*zba2(p4,p3,p5,p8)*s167**(-1) )
      amp(jdu2,1,2)= + L(jdu2)*cxw**(-2)*Hbit*cosW**(-2)*propw17**(-1)*
     & propw34**(-1)*propZ28**(-1) * (  - 2.D0*za(p3,p7)*za(p6,p8)*zb(
     &    p1,p4)*zb(p2,p5)*re*sqwmass*proph1347**(-1)*propZ56**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56eW(2)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*
     & propw3456**(-1) * (  - 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p5)*s3456**2 + 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*
     &    twop34Dp3456*s3456 - 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p4
     &    ,p5)*twop17Dp3456*s3456 + 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p5)*twop17Dp3456*twop34Dp3456 )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56eW(2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*propw3456**(-1) * (
     &     - 4.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*s1347 + 2.D0*
     &    za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*s3456 + 4.D0*za(p3,p6
     &    )*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*s34 + 4.D0*za(p3,p6)*za(p7,p8
     &    )*zb(p1,p2)*zb(p4,p5)*s17 - 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2
     &    )*zb(p4,p5)*twop34Dp3456 + 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)
     &    *zb(p4,p5)*twop17Dp3456 + 4.D0*za(p3,p6)*zb(p4,p5)*zab2(p7,p2
     &    ,p8,p1)*zab2(p8,p3,p4,p2) - 4.D0*za(p3,p6)*zb(p4,p5)*zab2(p7,
     &    p3,p4,p1)*zab2(p8,p1,p7,p2) - 4.D0*za(p3,p7)*zb(p1,p4)*zab2(
     &    p8,p1,p7,p2)*zba2(p5,p3,p4,p6) + 4.D0*za(p3,p8)*zb(p2,p4)*
     &    zab2(p7,p2,p8,p1)*zba2(p5,p3,p4,p6) + 4.D0*za(p6,p7)*zb(p1,p5
     &    )*zab2(p3,p5,p6,p4)*zab2(p8,p1,p7,p2) - 4.D0*za(p6,p8)*zb(p2,
     &    p5)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1) - 4.D0*za(p7,p8)*zb(
     &    p1,p2)*zab2(p3,p1,p7,p4)*zba2(p5,p3,p4,p6) + 4.D0*za(p7,p8)*
     &    zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p5,p1,p7,p6) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56eW(2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1) * ( 2.D0*za(p3,p6)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p5) - 4.D0*za(p3,p7)*za(p6,p8)*zb(
     &    p1,p4)*zb(p2,p5) + 2.D0*za(p3,p8)*za(p6,p7)*zb(p1,p5)*zb(p2,
     &    p4) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56ne(1,2)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s356**(-1) * (  - 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p3,p5)*
     &    zab2(p3,p5,p6,p4)*s3456 - 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*
     &    zb(p3,p5)*zab2(p3,p5,p6,p4)*twop17Dp3456 + 2.D0*za(p3,p6)*za(
     &    p7,p8)*zb(p1,p2)*zb(p5,p6)*zab2(p6,p3,p5,p4)*s3456 + 2.D0*za(
     &    p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p5,p6)*zab2(p6,p3,p5,p4)*
     &    twop17Dp3456 )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s356**(-1) * (  - 4.
     &    D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*zba2(p5,p3,p6,p1)
     &     + 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p3,p5)*zab2(p3,p5,p6
     &    ,p4) + 4.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p7)*zba2(p5,
     &    p3,p6,p7) - 2.D0*za(p3,p6)*za(p7,p8)*zb(p1,p2)*zb(p5,p6)*
     &    zab2(p6,p3,p5,p4) - 4.D0*za(p3,p6)*zb(p1,p4)*zab2(p8,p1,p7,p2
     &    )*zba2(p5,p3,p6,p7) + 4.D0*za(p3,p6)*zb(p2,p4)*zab2(p7,p2,p8,
     &    p1)*zba2(p5,p3,p6,p8) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56ee(1,2)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s456**(-1) * ( 2.D0*za(p4,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*
     &    zab2(p3,p5,p6,p4)*s3456 + 2.D0*za(p4,p6)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p5)*zab2(p3,p5,p6,p4)*twop17Dp3456 + 2.D0*za(p5,p6)*za(
     &    p7,p8)*zb(p1,p2)*zb(p4,p5)*zab2(p3,p4,p6,p5)*s3456 + 2.D0*za(
     &    p5,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*zab2(p3,p4,p6,p5)*
     &    twop17Dp3456 )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28qW(jdu2,1)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s456**(-1) * ( 4.D0
     &    *za(p1,p3)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*zba2(p1,p4,p5,p6) -
     &    4.D0*za(p3,p7)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*zba2(p7,p4,p5,p6
     &    ) - 4.D0*za(p3,p7)*zb(p4,p5)*zab2(p8,p1,p7,p2)*zba2(p1,p4,p5,
     &    p6) + 4.D0*za(p3,p8)*zb(p4,p5)*zab2(p7,p2,p8,p1)*zba2(p2,p4,
     &    p5,p6) - 2.D0*za(p4,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*zab2(p3
     &    ,p5,p6,p4) - 2.D0*za(p5,p6)*za(p7,p8)*zb(p1,p2)*zb(p4,p5)*
     &    zab2(p3,p4,p6,p5) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28eq(jdu2,1,1)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s356**(-1) * (  - 4.D0*za(p3,p6)*
     &    zb(p2,p4)*zba2(p1,p2,p4,p8)*zba2(p5,p3,p6,p7)*s248**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28eq(jdu2,1,1)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * ( 4.D0*za(p3,p7)*zb(p4
     &    ,p5)*zba2(p1,p3,p7,p8)*zba2(p2,p4,p5,p6)*s137**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28eq(jdu2,1,1)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1) * (  - 4.D0*za(p3,p7)*zb(p2,p4)*
     &    zba2(p1,p3,p7,p6)*zba2(p5,p2,p4,p8)*s137**(-1)*s248**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28nq(jdu2,1,1)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1)*s356**(-1) * (  - 4.D0*
     &    za(p3,p6)*zb(p1,p4)*zba2(p2,p1,p4,p7)*zba2(p5,p3,p6,p8) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28nq(jdu2,1,1)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1) * (  - 4.D0*za(p3,p8)*
     &    zb(p1,p4)*zba2(p2,p3,p8,p6)*zba2(p5,p1,p4,p7)*s238**(-1) )
      amp(jdu2,1,2) = amp(jdu2,1,2) + gamZ28nq(jdu2,1,1)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * ( 4.D0*za(p3,p8)*zb(p4
     &    ,p5)*zba2(p1,p4,p5,p6)*zba2(p2,p3,p8,p7)*s238**(-1) )
      amp(jdu2,2,1)= + cxw**(-2)*Hbit*cosW**(-2)*propw17**(-1)*
     & propw34**(-1)*propZ28**(-1) * (  - 2.D0*R(jdu2)*za(p2,p5)*za(p3,
     &    p7)*zb(p1,p4)*zb(p6,p8)*le*sqwmass*proph1347**(-1)*
     &    propZ56**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*cwmass2**(-1)*
     & cxw**(-2)*Bbit*propw17**(-1)*propw34**(-1)*propw3456**(-1)*
     & s345**(-1) * (  - za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,p4)*zab2(
     &    p3,p4,p5,p6)*s3456 - za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,p4)*
     &    zab2(p3,p4,p5,p6)*twop17Dp3456 + za(p2,p7)*za(p3,p5)*zb(p1,p8
     &    )*zb(p4,p5)*zab2(p5,p3,p4,p6)*s3456 + za(p2,p7)*za(p3,p5)*zb(
     &    p1,p8)*zb(p4,p5)*zab2(p5,p3,p4,p6)*twop17Dp3456 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*cxw**(-2)*Bbit*
     & propw17**(-1)*propw34**(-1)*propw3456**(-1)*s345**(-1) * (  - 2.D
     &    0*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p1,p8)*zba2(p4,p3,p5,p1)
     &     + za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,p4)*zab2(p3,p4,p5,p6)
     &     - za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p5)*zab2(p5,p3,p4,p6)
     &     + 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p6,p7)*zba2(p4,p3,p5
     &    ,p7) + 2.D0*za(p3,p5)*zb(p1,p6)*zba2(p4,p3,p5,p7)*zba2(p8,p1,
     &    p7,p2) + 2.D0*za(p3,p5)*zb(p6,p8)*zab2(p7,p2,p8,p1)*zba2(p4,
     &    p3,p5,p2) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56eW(1)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)
     &    *s3456**2 - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*
     &    twop34Dp3456*s3456 + 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4
     &    ,p6)*twop17Dp3456*s3456 - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*
     &    zb(p4,p6)*twop17Dp3456*twop34Dp3456 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56eW(1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*propw3456**(-1) * ( 4.
     &    D0*za(p2,p3)*zb(p4,p8)*zab2(p5,p3,p4,p6)*zab2(p7,p2,p8,p1) -
     &    4.D0*za(p2,p5)*zb(p6,p8)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1)
     &     + 4.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*s1347 - 2.D0*
     &    za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*s3456 - 4.D0*za(p2,p7
     &    )*za(p3,p5)*zb(p1,p8)*zb(p4,p6)*s34 - 4.D0*za(p2,p7)*za(p3,p5
     &    )*zb(p1,p8)*zb(p4,p6)*s17 + 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8
     &    )*zb(p4,p6)*twop34Dp3456 - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)
     &    *zb(p4,p6)*twop17Dp3456 + 4.D0*za(p2,p7)*zb(p1,p8)*zab2(p3,p1
     &    ,p7,p4)*zab2(p5,p3,p4,p6) - 4.D0*za(p2,p7)*zb(p1,p8)*zab2(p3,
     &    p5,p6,p4)*zab2(p5,p1,p7,p6) + 4.D0*za(p3,p5)*zb(p4,p6)*zab2(
     &    p7,p2,p8,p1)*zba2(p8,p3,p4,p2) - 4.D0*za(p3,p5)*zb(p4,p6)*
     &    zab2(p7,p3,p4,p1)*zba2(p8,p1,p7,p2) - 4.D0*za(p3,p7)*zb(p1,p4
     &    )*zab2(p5,p3,p4,p6)*zba2(p8,p1,p7,p2) + 4.D0*za(p5,p7)*zb(p1,
     &    p6)*zab2(p3,p5,p6,p4)*zba2(p8,p1,p7,p2) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56eW(1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1) * ( 2.D0*za(p2,p3)*
     &    za(p5,p7)*zb(p1,p6)*zb(p4,p8) - 4.D0*za(p2,p5)*za(p3,p7)*zb(
     &    p1,p4)*zb(p6,p8) - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,
     &    p6) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56ne(1,1)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s356**(-1) * ( 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*
     &    zab2(p3,p5,p6,p4)*s3456 + 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*
     &    zb(p3,p6)*zab2(p3,p5,p6,p4)*twop17Dp3456 + 2.D0*za(p2,p7)*za(
     &    p3,p5)*zb(p1,p8)*zb(p5,p6)*zab2(p5,p3,p6,p4)*s3456 + 2.D0*za(
     &    p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zab2(p5,p3,p6,p4)*
     &    twop17Dp3456 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s356**(-1) * ( 4.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*zba2(p6,p3,p5,p1) -
     &    2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*zab2(p3,p5,p6,p4
     &    ) - 4.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p4,p7)*zba2(p6,p3,
     &    p5,p7) - 2.D0*za(p2,p7)*za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zab2(p5
     &    ,p3,p6,p4) - 4.D0*za(p3,p5)*zb(p1,p4)*zba2(p6,p3,p5,p7)*zba2(
     &    p8,p1,p7,p2) - 4.D0*za(p3,p5)*zb(p4,p8)*zab2(p7,p2,p8,p1)*
     &    zba2(p6,p3,p5,p2) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56ee(1,1)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zab2(p3,p5,p6,p4)*s3456 - 2.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*
     &    zb(p4,p6)*zab2(p3,p5,p6,p4)*twop17Dp3456 + 2.D0*za(p2,p7)*za(
     &    p5,p6)*zb(p1,p8)*zb(p4,p6)*zab2(p3,p4,p5,p6)*s3456 + 2.D0*za(
     &    p2,p7)*za(p5,p6)*zb(p1,p8)*zb(p4,p6)*zab2(p3,p4,p5,p6)*
     &    twop17Dp3456 )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28qW(jdu2,2)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s456**(-1) * (  - 4.
     &    D0*za(p1,p3)*za(p2,p7)*zb(p1,p8)*zb(p4,p6)*zba2(p1,p4,p6,p5)
     &     - 4.D0*za(p2,p3)*zb(p4,p6)*zab2(p7,p2,p8,p1)*zba2(p8,p4,p6,
     &    p5) + 4.D0*za(p2,p7)*za(p3,p7)*zb(p1,p8)*zb(p4,p6)*zba2(p7,p4
     &    ,p6,p5) + 2.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*zab2(
     &    p3,p5,p6,p4) - 2.D0*za(p2,p7)*za(p5,p6)*zb(p1,p8)*zb(p4,p6)*
     &    zab2(p3,p4,p5,p6) - 4.D0*za(p3,p7)*zb(p4,p6)*zba2(p1,p4,p6,p5
     &    )*zba2(p8,p1,p7,p2) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28eq(jdu2,1,2)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s356**(-1) * ( 4.D0*za(p3,p5)*zb(p4
     &    ,p8)*zba2(p1,p4,p8,p2)*zba2(p6,p3,p5,p7)*s248**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28eq(jdu2,1,2)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * ( 4.D0*za(p3,p7)*zb(p4
     &    ,p6)*zba2(p1,p3,p7,p2)*zba2(p8,p4,p6,p5)*s137**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28eq(jdu2,1,2)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1) * ( 4.D0*za(p3,p7)*zb(p4,p8)*zba2(
     &    p1,p3,p7,p5)*zba2(p6,p4,p8,p2)*s137**(-1)*s248**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28nq(jdu2,1,2)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1)*s356**(-1) * (  - 4.D0*
     &    za(p3,p5)*zb(p1,p4)*zba2(p6,p3,p5,p2)*zba2(p8,p1,p4,p7) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28nq(jdu2,1,2)*gamZ56ne(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1) * ( 4.D0*za(p2,p3)*zb(p1
     &    ,p4)*zba2(p6,p1,p4,p7)*zba2(p8,p2,p3,p5)*s238**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28nq(jdu2,1,2)*gamZ56ee(1,1)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * (  - 4.D0*za(p2,p3)*
     &    zb(p4,p6)*zba2(p1,p4,p6,p5)*zba2(p8,p2,p3,p7)*s238**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28eqZ(jdu2,1,2)*cxw**(-2)*
     & Bbit*propw17**(-1)*propw34**(-1)*s258**(-1) * ( 2.D0*za(p2,p5)*
     &    zb(p1,p6)*zba2(p4,p1,p6,p7)*zba2(p8,p2,p5,p3)*s167**(-1) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28eqZ(jdu2,1,2)*cxw**(-2)*
     & Bbit*propw17**(-1)*propw34**(-1)*s268**(-1)*s345**(-1) * (  - 2.D
     &    0*za(p3,p5)*zb(p6,p8)*zba2(p1,p6,p8,p2)*zba2(p4,p3,p5,p7) )
      amp(jdu2,2,1) = amp(jdu2,2,1) + gamZ28nqZ(jdu2,1,2)*cxw**(-2)*
     & Bbit*propw17**(-1)*propw34**(-1)*s345**(-1) * ( 2.D0*za(p3,p5)*
     &    zb(p1,p6)*zba2(p4,p3,p5,p2)*zba2(p8,p1,p6,p7)*s167**(-1) )
      amp(jdu2,2,2)= + cxw**(-2)*Hbit*cosW**(-2)*propw17**(-1)*
     & propw34**(-1)*propZ28**(-1) * (  - 2.D0*R(jdu2)*za(p2,p6)*za(p3,
     &    p7)*zb(p1,p4)*zb(p5,p8)*re*sqwmass*proph1347**(-1)*
     &    propZ56**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56eW(2)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*
     & propw3456**(-1) * ( 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)
     &    *s3456**2 - 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*
     &    twop34Dp3456*s3456 + 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4
     &    ,p5)*twop17Dp3456*s3456 - 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*
     &    zb(p4,p5)*twop17Dp3456*twop34Dp3456 )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56eW(2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1)*propw3456**(-1) * ( 4.
     &    D0*za(p2,p3)*zb(p4,p8)*zab2(p7,p2,p8,p1)*zba2(p5,p3,p4,p6) -
     &    4.D0*za(p2,p6)*zb(p5,p8)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1)
     &     + 4.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*s1347 - 2.D0*
     &    za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*s3456 - 4.D0*za(p2,p7
     &    )*za(p3,p6)*zb(p1,p8)*zb(p4,p5)*s34 - 4.D0*za(p2,p7)*za(p3,p6
     &    )*zb(p1,p8)*zb(p4,p5)*s17 + 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8
     &    )*zb(p4,p5)*twop34Dp3456 - 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)
     &    *zb(p4,p5)*twop17Dp3456 + 4.D0*za(p2,p7)*zb(p1,p8)*zab2(p3,p1
     &    ,p7,p4)*zba2(p5,p3,p4,p6) - 4.D0*za(p2,p7)*zb(p1,p8)*zab2(p3,
     &    p5,p6,p4)*zba2(p5,p1,p7,p6) + 4.D0*za(p3,p6)*zb(p4,p5)*zab2(
     &    p7,p2,p8,p1)*zba2(p8,p3,p4,p2) - 4.D0*za(p3,p6)*zb(p4,p5)*
     &    zab2(p7,p3,p4,p1)*zba2(p8,p1,p7,p2) - 4.D0*za(p3,p7)*zb(p1,p4
     &    )*zba2(p5,p3,p4,p6)*zba2(p8,p1,p7,p2) + 4.D0*za(p6,p7)*zb(p1,
     &    p5)*zab2(p3,p5,p6,p4)*zba2(p8,p1,p7,p2) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56eW(2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw34**(-1) * ( 2.D0*za(p2,p3)*
     &    za(p6,p7)*zb(p1,p5)*zb(p4,p8) - 4.D0*za(p2,p6)*za(p3,p7)*zb(
     &    p1,p4)*zb(p5,p8) - 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,
     &    p5) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56ne(1,2)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s356**(-1) * ( 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*
     &    zab2(p3,p5,p6,p4)*s3456 + 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*
     &    zb(p3,p5)*zab2(p3,p5,p6,p4)*twop17Dp3456 - 2.D0*za(p2,p7)*za(
     &    p3,p6)*zb(p1,p8)*zb(p5,p6)*zab2(p6,p3,p5,p4)*s3456 - 2.D0*za(
     &    p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p5,p6)*zab2(p6,p3,p5,p4)*
     &    twop17Dp3456 )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s356**(-1) * ( 4.D0
     &    *za(p2,p7)*za(p3,p6)*zb(p1,p4)*zb(p1,p8)*zba2(p5,p3,p6,p1) -
     &    2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p3,p5)*zab2(p3,p5,p6,p4
     &    ) - 4.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p4,p7)*zba2(p5,p3,
     &    p6,p7) + 2.D0*za(p2,p7)*za(p3,p6)*zb(p1,p8)*zb(p5,p6)*zab2(p6
     &    ,p3,p5,p4) - 4.D0*za(p3,p6)*zb(p1,p4)*zba2(p5,p3,p6,p7)*zba2(
     &    p8,p1,p7,p2) - 4.D0*za(p3,p6)*zb(p4,p8)*zab2(p7,p2,p8,p1)*
     &    zba2(p5,p3,p6,p2) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56ee(1,2)*
     & cwmass2**(-1)*cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*
     & s456**(-1) * (  - 2.D0*za(p2,p7)*za(p4,p6)*zb(p1,p8)*zb(p4,p5)*
     &    zab2(p3,p5,p6,p4)*s3456 - 2.D0*za(p2,p7)*za(p4,p6)*zb(p1,p8)*
     &    zb(p4,p5)*zab2(p3,p5,p6,p4)*twop17Dp3456 - 2.D0*za(p2,p7)*za(
     &    p5,p6)*zb(p1,p8)*zb(p4,p5)*zab2(p3,p4,p6,p5)*s3456 - 2.D0*za(
     &    p2,p7)*za(p5,p6)*zb(p1,p8)*zb(p4,p5)*zab2(p3,p4,p6,p5)*
     &    twop17Dp3456 )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28qW(jdu2,2)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*propw3456**(-1)*s456**(-1) * (  - 4.
     &    D0*za(p1,p3)*za(p2,p7)*zb(p1,p8)*zb(p4,p5)*zba2(p1,p4,p5,p6)
     &     - 4.D0*za(p2,p3)*zb(p4,p5)*zab2(p7,p2,p8,p1)*zba2(p8,p4,p5,
     &    p6) + 4.D0*za(p2,p7)*za(p3,p7)*zb(p1,p8)*zb(p4,p5)*zba2(p7,p4
     &    ,p5,p6) + 2.D0*za(p2,p7)*za(p4,p6)*zb(p1,p8)*zb(p4,p5)*zab2(
     &    p3,p5,p6,p4) + 2.D0*za(p2,p7)*za(p5,p6)*zb(p1,p8)*zb(p4,p5)*
     &    zab2(p3,p4,p6,p5) - 4.D0*za(p3,p7)*zb(p4,p5)*zba2(p1,p4,p5,p6
     &    )*zba2(p8,p1,p7,p2) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28eq(jdu2,1,2)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s356**(-1) * ( 4.D0*za(p3,p6)*zb(p4
     &    ,p8)*zba2(p1,p4,p8,p2)*zba2(p5,p3,p6,p7)*s248**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28eq(jdu2,1,2)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * ( 4.D0*za(p3,p7)*zb(p4
     &    ,p5)*zba2(p1,p3,p7,p2)*zba2(p8,p4,p5,p6)*s137**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28eq(jdu2,1,2)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1) * ( 4.D0*za(p3,p7)*zb(p4,p8)*zba2(
     &    p1,p3,p7,p6)*zba2(p5,p4,p8,p2)*s137**(-1)*s248**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28nq(jdu2,1,2)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1)*s356**(-1) * (  - 4.D0*
     &    za(p3,p6)*zb(p1,p4)*zba2(p5,p3,p6,p2)*zba2(p8,p1,p4,p7) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28nq(jdu2,1,2)*gamZ56ne(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s147**(-1) * ( 4.D0*za(p2,p3)*zb(p1
     &    ,p4)*zba2(p5,p1,p4,p7)*zba2(p8,p2,p3,p6)*s238**(-1) )
      amp(jdu2,2,2) = amp(jdu2,2,2) + gamZ28nq(jdu2,1,2)*gamZ56ee(1,2)*
     & cxw**(-1)*Bbit*propw17**(-1)*s456**(-1) * (  - 4.D0*za(p2,p3)*
     &    zb(p4,p5)*zba2(p1,p4,p5,p6)*zba2(p8,p2,p3,p7)*s238**(-1) )
      enddo

      return
      end
