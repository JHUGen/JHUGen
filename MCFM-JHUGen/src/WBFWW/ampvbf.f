      subroutine ampvbf(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'runstring.f'
      include 'zcouple.f'
      include 'WWbits.f'
      integer:: jdu1,jdu2,h28,h17,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,amp(2,2,2,2),sqzmass,
     & propw34,propw56,propz28,propz17,propH,propw1347,propw1567,
     & gamZ17(2,2),gamz28(2,2),ZZ17(2,2),ZZ28(2,2),rxw
      real(dp):: t4,s17,s28,s34,s56,s3456,s1347,s1567
C     amp(jdu1,jdu2,h17,h28)
C-----Begin statement functions
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

c--- special fix for Madgraph check
      if (index(runstring,'mad') > 0) then
        sqzmass=cplx2(zmass**2,zip)
      else
        sqzmass=czmass2
      endif

      rxw=sqrt((cone-cxw)/cxw)
      s17=s(i1,i7)
      s28=s(i2,i8)
      s34=s(i3,i4)
      s56=s(i5,i6)
      s1347=t4(i1,i3,i4,i7)
      s1567=t4(i1,i5,i6,i7)
      s3456=t4(i3,i4,i5,i6)

      propw1347=s1347-cwmass2
      propw1567=s1567-cwmass2
      propw34=s34-cwmass2
      propw56=s56-cwmass2
      propz17=s17-czmass2
      propz28=s28-czmass2
      propH=cplx2(s3456-hmass**2,hmass*hwidth)

      do jdu1=1,2
      gamz17(jdu1,1)=Q(jdu1)/s(i1,i7)+rxW*L(jdu1)/propz17
      gamz17(jdu1,2)=Q(jdu1)/s(i1,i7)+rxW*R(jdu1)/propz17
      gamz28(jdu1,1)=Q(jdu1)/s(i2,i8)+rxW*L(jdu1)/propz28
      gamz28(jdu1,2)=Q(jdu1)/s(i2,i8)+rxW*R(jdu1)/propz28
      ZZ17(jdu1,1)=L(jdu1)/propz17
      ZZ17(jdu1,2)=R(jdu1)/propz17
      ZZ28(jdu1,1)=L(jdu1)/propz28
      ZZ28(jdu1,2)=R(jdu1)/propz28
      enddo

      p3=i3
      p4=i4
      p5=i5
      p6=i6

      do h17=1,2
      if (h17 == 1) then
      p1=i1
      p7=i7
      elseif (h17 == 2) then
      p1=i7
      p7=i1
      endif
      do h28=1,2
      if (h28 == 1) then
      p2=i2
      p8=i8
      elseif (h28 == 2) then
      p2=i8
      p8=i2
      endif
      do jdu1=1,2
      do jdu2=1,2
      amp(jdu1,jdu2,h17,h28)= + propw56**(-1)*propw34**(-1)*cxw**(-2)*
     & Hbit * (  - 2.D0*za(p3,p5)*za(p7,p8)*zb(p4,p6)*zb(p1,p2)*ZZ17(
     &    jdu1,h17)*ZZ28(jdu2,h28)*propH**(-1)*sqzmass )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamz17(jdu1,h17
     & )*gamz28(jdu2,h28)*propw56**(-1)*propw34**(-1)*propw1347**(-1)*
     & cxw**(-1)*Bbit * ( 4.D0*za(p3,p5)*zb(p4,p6)*zab2(p7,p3,p4,p1)*
     &    zab2(p8,p5,p6,p2) - 2.D0*za(p3,p7)*za(p5,p8)*zb(p4,p1)*zb(p6,
     &    p2)*cwmass2**(-1)*s17*s28 + 2.D0*za(p3,p7)*za(p5,p8)*zb(p4,p1
     &    )*zb(p6,p2)*cwmass2**(-1)*s56*s17 + 2.D0*za(p3,p7)*za(p5,p8)*
     &    zb(p4,p1)*zb(p6,p2)*cwmass2**(-1)*s34*s28 - 2.D0*za(p3,p7)*
     &    za(p5,p8)*zb(p4,p1)*zb(p6,p2)*cwmass2**(-1)*s34*s56 - 2.D0*
     &    za(p3,p7)*za(p5,p8)*zb(p4,p1)*zb(p6,p2)*s3456 + 2.D0*za(p3,p7
     &    )*za(p5,p8)*zb(p4,p1)*zb(p6,p2)*s1567 - 2.D0*za(p3,p7)*zb(p4,
     &    p1)*zab2(p5,p3,p4,p6)*zab2(p8,p5,p6,p2) + 2.D0*za(p3,p7)*zb(
     &    p4,p1)*zab2(p5,p1,p7,p6)*zab2(p8,p5,p6,p2) + 2.D0*za(p3,p7)*
     &    zb(p4,p1)*zab2(p5,p2,p8,p6)*zab2(p8,p3,p4,p2) - 2.D0*za(p3,p7
     &    )*zb(p4,p1)*zab2(p5,p2,p8,p6)*zab2(p8,p1,p7,p2) - 4.D0*za(p3,
     &    p8)*zb(p4,p2)*zab2(p5,p2,p8,p6)*zab2(p7,p3,p4,p1) - 4.D0*za(
     &    p5,p7)*zb(p6,p1)*zab2(p3,p1,p7,p4)*zab2(p8,p5,p6,p2) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamz17(jdu1,h17
     & )*gamz28(jdu2,h28)*propw56**(-1)*propw34**(-1)*propw1347**(-1)*
     & cxw**(-1)*Bbit * (  - 2.D0*za(p5,p8)*zb(p6,p2)*zab2(p3,p5,p6,p4)
     &    *zab2(p7,p3,p4,p1) + 2.D0*za(p5,p8)*zb(p6,p2)*zab2(p3,p1,p7,
     &    p4)*zab2(p7,p5,p6,p1) - 2.D0*za(p5,p8)*zb(p6,p2)*zab2(p3,p1,
     &    p7,p4)*zab2(p7,p2,p8,p1) + 2.D0*za(p5,p8)*zb(p6,p2)*zab2(p3,
     &    p2,p8,p4)*zab2(p7,p3,p4,p1) + 4.D0*za(p7,p8)*zb(p1,p2)*zab2(
     &    p3,p1,p7,p4)*zab2(p5,p2,p8,p6) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamz17(jdu1,h17
     & )*gamz28(jdu2,h28)*propw56**(-1)*propw34**(-1)*propw1567**(-1)*
     & cxw**(-1)*Bbit * ( 4.D0*za(p3,p5)*zb(p4,p6)*zab2(p7,p5,p6,p1)*
     &    zab2(p8,p3,p4,p2) - 4.D0*za(p3,p7)*zb(p4,p1)*zab2(p5,p1,p7,p6
     &    )*zab2(p8,p3,p4,p2) - 2.D0*za(p3,p8)*za(p5,p7)*zb(p4,p2)*zb(
     &    p6,p1)*cwmass2**(-1)*s17*s28 + 2.D0*za(p3,p8)*za(p5,p7)*zb(p4
     &    ,p2)*zb(p6,p1)*cwmass2**(-1)*s56*s28 + 2.D0*za(p3,p8)*za(p5,
     &    p7)*zb(p4,p2)*zb(p6,p1)*cwmass2**(-1)*s34*s17 - 2.D0*za(p3,p8
     &    )*za(p5,p7)*zb(p4,p2)*zb(p6,p1)*cwmass2**(-1)*s34*s56 - 2.D0*
     &    za(p3,p8)*za(p5,p7)*zb(p4,p2)*zb(p6,p1)*s3456 + 2.D0*za(p3,p8
     &    )*za(p5,p7)*zb(p4,p2)*zb(p6,p1)*s1347 - 2.D0*za(p3,p8)*zb(p4,
     &    p2)*zab2(p5,p3,p4,p6)*zab2(p7,p5,p6,p1) + 2.D0*za(p3,p8)*zb(
     &    p4,p2)*zab2(p5,p1,p7,p6)*zab2(p7,p3,p4,p1) - 2.D0*za(p3,p8)*
     &    zb(p4,p2)*zab2(p5,p1,p7,p6)*zab2(p7,p2,p8,p1) + 2.D0*za(p3,p8
     &    )*zb(p4,p2)*zab2(p5,p2,p8,p6)*zab2(p7,p5,p6,p1) - 2.D0*za(p5,
     &    p7)*zb(p6,p1)*zab2(p3,p5,p6,p4)*zab2(p8,p3,p4,p2) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamz17(jdu1,h17
     & )*gamz28(jdu2,h28)*propw56**(-1)*propw34**(-1)*propw1567**(-1)*
     & cxw**(-1)*Bbit * ( 2.D0*za(p5,p7)*zb(p6,p1)*zab2(p3,p1,p7,p4)*
     &    zab2(p8,p3,p4,p2) + 2.D0*za(p5,p7)*zb(p6,p1)*zab2(p3,p2,p8,p4
     &    )*zab2(p8,p5,p6,p2) - 2.D0*za(p5,p7)*zb(p6,p1)*zab2(p3,p2,p8,
     &    p4)*zab2(p8,p1,p7,p2) - 4.D0*za(p5,p8)*zb(p6,p2)*zab2(p3,p2,
     &    p8,p4)*zab2(p7,p5,p6,p1) + 4.D0*za(p7,p8)*zb(p1,p2)*zab2(p3,
     &    p2,p8,p4)*zab2(p5,p1,p7,p6) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamz17(jdu1,h17
     & )*gamz28(jdu2,h28)*propw56**(-1)*propw34**(-1)*cxw**(-1)*Bbit
     &  * (  - 4.D0*za(p3,p5)*za(p7,p8)*zb(p4,p6)*zb(p1,p2) + 2.D0*za(
     &    p3,p7)*za(p5,p8)*zb(p4,p1)*zb(p6,p2) + 2.D0*za(p3,p8)*za(p5,
     &    p7)*zb(p4,p2)*zb(p6,p1) )
      enddo
      enddo
      enddo
      enddo
      return
      end
