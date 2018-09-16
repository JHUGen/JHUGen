      subroutine amptwodiags(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8,iswap
      complex(dp):: amp(2,2,2,2),propW34,propW56,
     & propW1347,propW1567,tmp
      real(dp):: t3,t4,s34,s56,s1347,s1567,
     & s134,s156,s256,s347,s348,s567
C-----Begin statement functions
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      p1=i1
      p2=i2
      p7=i7
      p8=i8

      amp(:,:,:,:)=czip
c--- need to check corresponding result for (2,1,1,1)
      do iswap=1,2
      if (iswap == 1) then
        p3=i3
        p4=i4
        p5=i5
        p6=i6
      else
        p3=i5
        p4=i6
        p5=i3
        p6=i4
      endif
      s34=s(p3,p4)
      s56=s(p5,p6)
      s1347=t4(p1,p3,p4,p7)
      s1567=t4(p1,p5,p6,p7)
      s134=t3(p1,p3,p4)
      s156=t3(p1,p5,p6)
      s256=t3(p2,p5,p6)
      s347=t3(p3,p4,p7)
      s348=t3(p3,p4,p8)
      s567=t3(p5,p6,p7)
      propW34=s34-cwmass2
      propW56=s56-cwmass2
      propW1347=s1347-cwmass2
      propW1567=s1567-cwmass2
      tmp= + s134**(-1)*s256**(-1)*propW34**(-1)*propW56**(-1)*
     & propW1347**(-1) * (  - za(p7,p8)*za(p1,p3)*za(p2,p5)*zb(p1,p2)*
     &    zb(p1,p4)*zb(p2,p6) + za(p7,p8)*za(p1,p3)*za(p5,p6)*zb(p1,p6)
     &    *zb(p1,p4)*zb(p2,p6) - za(p7,p8)*za(p2,p5)*za(p3,p4)*zb(p1,p4
     &    )*zb(p2,p6)*zb(p2,p4) + za(p7,p8)*za(p5,p6)*za(p3,p4)*zb(p1,
     &    p4)*zb(p2,p6)*zb(p6,p4) )
      tmp = tmp + s348**(-1)*s567**(-1)*propW34**(-1)*propW56**(-1)*
     & propW1567**(-1) * (  - za(p7,p8)*za(p7,p5)*za(p8,p3)*zb(p7,p6)*
     &    zb(p1,p2)*zb(p8,p4) - za(p7,p5)*za(p7,p3)*za(p8,p3)*zb(p7,p6)
     &    *zb(p1,p2)*zb(p3,p4) + za(p7,p5)*za(p8,p5)*za(p8,p3)*zb(p1,p2
     &    )*zb(p8,p4)*zb(p5,p6) - za(p7,p5)*za(p8,p3)*za(p5,p3)*zb(p1,
     &    p2)*zb(p5,p6)*zb(p3,p4) )
      tmp = tmp + propW34**(-1)*propW56**(-1)*propW1347**(-1) * ( 1.D0/
     &    2.D0*za(p7,p3)*za(p8,p5)*zb(p1,p4)*zb(p2,p6)*cwmass2**(-1) )
      tmp = tmp + propW34**(-1)*propW56**(-1)*propW1567**(-1) * ( 1.D0/
     &    2.D0*za(p7,p5)*za(p8,p3)*zb(p1,p6)*zb(p2,p4)*cwmass2**(-1) )
      if (iswap == 1) amp(1,2,1,1)=tmp
      if (iswap == 2) amp(2,1,1,1)=tmp
      tmp= + s256**(-1)*propW34**(-1)*propW56**(-1)*propW1347**(-1)
     &  * ( za(p7,p8)*za(p7,p3)*za(p2,p5)*zb(p7,p4)*zb(p1,p2)*zb(p2,p6)
     &    *s347**(-1) - za(p7,p8)*za(p7,p3)*za(p5,p6)*zb(p7,p4)*zb(p1,
     &    p6)*zb(p2,p6)*s347**(-1) - za(p7,p3)*za(p8,p3)*za(p2,p5)*zb(
     &    p1,p2)*zb(p2,p6)*zb(p3,p4)*s347**(-1) + za(p7,p3)*za(p8,p3)*
     &    za(p5,p6)*zb(p1,p6)*zb(p2,p6)*zb(p3,p4)*s347**(-1) )
      tmp = tmp + s348**(-1)*propW34**(-1)*propW56**(-1)*
     & propW1567**(-1) * ( za(p7,p8)*za(p1,p5)*za(p8,p3)*zb(p1,p2)*zb(
     &    p1,p6)*zb(p8,p4)*s156**(-1) + za(p7,p8)*za(p8,p3)*za(p5,p6)*
     &    zb(p1,p6)*zb(p8,p4)*zb(p2,p6)*s156**(-1) + za(p7,p3)*za(p1,p5
     &    )*za(p8,p3)*zb(p1,p2)*zb(p1,p6)*zb(p3,p4)*s156**(-1) + za(p7,
     &    p3)*za(p8,p3)*za(p5,p6)*zb(p1,p6)*zb(p2,p6)*zb(p3,p4)*
     &    s156**(-1) )
      tmp = tmp + propW34**(-1)*propW56**(-1)*propW1347**(-1) * (  - 1.D
     &    0/2.D0*za(p7,p3)*za(p8,p5)*zb(p1,p4)*zb(p2,p6)*cwmass2**(-1)
     &     )
      tmp = tmp + propW34**(-1)*propW56**(-1)*propW1567**(-1) * (  - 1.D
     &    0/2.D0*za(p7,p5)*za(p8,p3)*zb(p1,p6)*zb(p2,p4)*cwmass2**(-1)
     &     )
      if (iswap == 1) amp(2,2,1,1)=tmp
      if (iswap == 2) amp(1,1,1,1)=tmp
      enddo

      amp(:,:,:,:)=amp(:,:,:,:)/cxw**3

      return
      end
