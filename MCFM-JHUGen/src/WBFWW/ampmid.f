      subroutine ampmid(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cmplxmass.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer:: jdu1,jdu2,h28,h17,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zba2,amp(2,2,2,2),
     & propw34,propw56,propz28,propz17,
     & gam17e(2,2),gam17n(2,2),gam28e(2,2),gam28n(2,2)
      real(dp):: t3,
     & s456,s345,s356,s346,s137,s147,s157,s167,s238,s248,s258,s268
C     amp(jdu1,jdu2,h17,h28)
C-----Begin statement functions
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
C-----end statement functions
      s456=t3(i4,i5,i6)
      s345=t3(i3,i4,i5)
      s356=t3(i3,i5,i6)
      s346=t3(i3,i4,i6)
      s137=t3(i1,i3,i7)
      s147=t3(i1,i4,i7)
      s157=t3(i1,i5,i7)
      s167=t3(i1,i6,i7)
      s238=t3(i2,i3,i8)
      s248=t3(i2,i4,i8)
      s258=t3(i2,i5,i8)
      s268=t3(i2,i6,i8)

      propw34=s(i3,i4)-cwmass2
      propw56=s(i5,i6)-cwmass2
      propz17=s(i1,i7)-czmass2
      propz28=s(i2,i8)-czmass2

      do jdu1=1,2
      gam17e(jdu1,1)=Q(jdu1)*qe/s(i1,i7)+L(jdu1)*le/propz17
      gam17e(jdu1,2)=Q(jdu1)*qe/s(i1,i7)+R(jdu1)*le/propz17
      gam17n(jdu1,1)=L(jdu1)*ln/propz17
      gam17n(jdu1,2)=R(jdu1)*ln/propz17
      gam28e(jdu1,1)=Q(jdu1)*qe/s(i2,i8)+L(jdu1)*le/propz28
      gam28e(jdu1,2)=Q(jdu1)*qe/s(i2,i8)+R(jdu1)*le/propz28
      gam28n(jdu1,1)=L(jdu1)*ln/propz28
      gam28n(jdu1,2)=R(jdu1)*ln/propz28
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
      amp(jdu1,jdu2,h17,h28)= + gam17e(jdu1,h17)*gam28e(jdu2,h28) * (
     &     - 4.D0*za(p3,p5)*zb(p6,p1)*zba2(p4,p3,p5,p8)*zba2(p2,p1,p6,
     &    p7)*propw34**(-1)*s345**(-1)*s167**(-1) - 4.D0*za(p3,p5)*zb(
     &    p6,p2)*zba2(p4,p3,p5,p7)*zba2(p1,p2,p6,p8)*propw34**(-1)*
     &    s345**(-1)*s268**(-1) + 4.D0*za(p3,p7)*zb(p4,p6)*zba2(p1,p3,
     &    p7,p8)*zba2(p2,p4,p6,p5)*propw56**(-1)*s456**(-1)*s137**(-1)
     &     + 4.D0*za(p3,p8)*zb(p4,p6)*zba2(p1,p4,p6,p5)*zba2(p2,p3,p8,
     &    p7)*propw56**(-1)*s456**(-1)*s238**(-1) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gam17e(jdu1,h17
     & )*gam28n(jdu2,h28) * ( 4.D0*za(p3,p7)*zb(p4,p2)*zba2(p6,p2,p4,p8
     &    )*zba2(p1,p3,p7,p5)*propw56**(-1)*s137**(-1)*s248**(-1) + 4.D0
     &    *za(p5,p8)*zb(p6,p1)*zba2(p4,p1,p6,p7)*zba2(p2,p5,p8,p3)*
     &    propw34**(-1)*s167**(-1)*s258**(-1) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gam17n(jdu1,h17
     & )*gam28e(jdu2,h28) * ( 4.D0*za(p3,p8)*zb(p4,p1)*zba2(p6,p1,p4,p7
     &    )*zba2(p2,p3,p8,p5)*propw56**(-1)*s147**(-1)*s238**(-1) + 4.D0
     &    *za(p5,p7)*zb(p6,p2)*zba2(p4,p2,p6,p8)*zba2(p1,p5,p7,p3)*
     &    propw34**(-1)*s157**(-1)*s268**(-1) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gam17n(jdu1,h17
     & )*gam28n(jdu2,h28) * ( 4.D0*za(p3,p5)*zb(p4,p1)*zba2(p6,p3,p5,p8
     &    )*zba2(p2,p1,p4,p7)*propw56**(-1)*s356**(-1)*s147**(-1) + 4.D0
     &    *za(p3,p5)*zb(p4,p2)*zba2(p6,p3,p5,p7)*zba2(p1,p2,p4,p8)*
     &    propw56**(-1)*s356**(-1)*s248**(-1) - 4.D0*za(p5,p7)*zb(p4,p6
     &    )*zba2(p1,p5,p7,p8)*zba2(p2,p4,p6,p3)*propw34**(-1)*
     &    s346**(-1)*s157**(-1) - 4.D0*za(p5,p8)*zb(p4,p6)*zba2(p1,p4,
     &    p6,p3)*zba2(p2,p5,p8,p7)*propw34**(-1)*s346**(-1)*s258**(-1)
     &     )
      amp(jdu1,jdu2,h17,h28)=amp(jdu1,jdu2,h17,h28)/cxw
      enddo
      enddo
      enddo
      enddo
      return
      end
