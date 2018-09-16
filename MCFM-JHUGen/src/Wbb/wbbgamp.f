      subroutine wbbgamp(j1,j2,j3,j4,j5,j6,j7,rmsq) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,h1,h2
      complex(dp):: xa(2,2),xb(2,2),xi(2,2),xf(2,2),t2
      complex(dp):: xd1(2,2),xd2(2,2),xd3(2,2),
     &               xd4(2,2),xd5(2,2),xd6(2,2),
     &               xd7(2,2),xd8(2,2),xd9(2,2),
     &               xd10(2,2),xd11(2,2),xd12(2,2)
      real(dp):: s45,s345,s2345,s267,s2367,s245,rmsq,prop

      s45=s(j4,j5)
      s245=s(j2,j4)+s(j2,j5)+s(j4,j5)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      s2367=s(j2,j3)+s(j3,j6)+s(j3,j7)+s267

C****************************************************************
C      t2=zb(j1,j2)*za(j4,j2)+zb(j1,j3)*za(j4,j3)
C  choose gauge vector jb=j4

      xd1(1,2)=t2(j4,j2,j3,j5)*t2(j4,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,j4)*zb(j3,j2)*s45*s2345)

      xd2(1,2)=t2(j4,j2,j3,j6)*t2(j7,j1,j4,j5)
     & *zb(j4,j1)/(zb(j3,j4)*zb(j3,j2)*s45*s2367)

      xd3(1,2)=-za(j6,j2)*t2(j7,j2,j6,j3)*zb(j4,j1)**2*za(j5,j1)
     & /(zb(j3,j4)*s45*s267*s2367)

      xd4(1,2)=-za(j5,j2)*t2(j4,j2,j5,j3)*t2(j4,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,j4)*s45*s245*s2345)

      xd5(1,2)=za(j6,j2)*t2(j7,j2,j6,j5)*zb(j4,j1)
     & *zb(j4,j1)/(zb(j3,j4)*zb(j3,j1)*s45*s267)

      xd6(1,2)=za(j5,j2)*t2(j4,j2,j5,j6)*zb(j7,j1)
     & *zb(j4,j1)/(zb(j3,j4)*zb(j3,j1)*s45*s245)

      xd7(1,2)=zb(j7,j1)/(zb(j3,j4)*s45*s345*s2345)
     & *t2(j4,j1,j7,j6)*za(j3,j5)
     & *(za(j3,j2)*zb(j3,j4)-za(j5,j2)*zb(j4,j5))

      xd8(1,2)=za(j6,j2)*za(j3,j5)*zb(j4,j1)/(zb(j3,j4)*s45*s345*s267)
     & *(-t2(j7,j2,j6,j3)*zb(j3,j4)+t2(j7,j2,j6,j5)*zb(j4,j5))

      xd9(1,2)=czip

      xd10(1,2)=-zb(j7,j1)*t2(j4,j1,j7,j6)
     & *(za(j5,j2)*zb(j5,j4)+za(j3,j2)*zb(j3,j4))
     & /(zb(j3,j4)*zb(j5,j3)*s345*s2345)

      xd11(1,2)=czip

      xd12(1,2)=za(j6,j2)*zb(j4,j1)
     & *(t2(j7,j2,j6,j3)*zb(j3,j4)+t2(j7,j2,j6,j5)*zb(j5,j4))
     & /(zb(j3,j4)*zb(j5,j3)*s345*s267)

C***********************************************************************

      xd1(2,1)=
     & -za(j4,j2)*za(j4,j2)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,j4)*za(j3,j2)*s45*s2345)

      xd2(2,1)=
     & -za(j4,j2)*za(j6,j2)*t2(j7,j1,j5,j4)
     & *zb(j5,j1)/(za(j3,j4)*za(j3,j2)*s45*s2367)

      xd3(2,1)=
     & +za(j6,j2)*t2(j7,j2,j6,j4)*t2(j3,j1,j5,j4)
     & *zb(j5,j1)/(za(j3,j4)*s45*s267*s2367)

      xd4(2,1)=za(j4,j2)**2*zb(j5,j2)*t2(j3,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,j4)*s45*s245*s2345)

      xd5(2,1)=za(j6,j2)*t2(j7,j2,j6,j4)*t2(j5,j1,j3,j4)
     & *zb(j3,j1)/(za(j3,j4)*s(j1,j3)*s45*s267)

      xd6(2,1)=
     & -za(j4,j2)*t2(j5,j2,j4,j6)*t2(j7,j1,j3,j4)
     & /(za(j3,j4)*za(j3,j1)*s45*s245)

      xd7(2,1)=-zb(j7,j1)*za(j4,j2)*zb(j3,j5)/(za(j3,j4)*s45*s345*s2345)
     & *(+t2(j3,j1,j7,j6)*za(j3,j4)-t2(j5,j1,j7,j6)*za(j4,j5))

      xd8(2,1)=-za(j6,j2)*zb(j3,j5)*t2(j7,j2,j6,j4)
     & /(za(j3,j4)*s45*s345*s267)
     & *(-zb(j3,j1)*za(j3,j4)+zb(j5,j1)*za(j4,j5))

      xd9(2,1)=czip

      xd10(2,1)=za(j2,j4)*zb(j7,j1)
     & *(za(j6,j1)*t2(j1,j3,j5,j4)+za(j6,j7)*t2(j7,j3,j5,j4))
     & /(za(j3,j4)*za(j3,j5)*s345*s2345)

      xd11(2,1)=czip

      xd12(2,1)=-za(j6,j2)*t2(j7,j2,j6,j4)
     & *(za(j4,j3)*zb(j3,j1)+za(j4,j5)*zb(j5,j1))
     & /(za(j3,j4)*za(j3,j5)*s345*s267)


*******************************************************************************

C   Choose gauge vector jb=j5

      xd1(1,1)=t2(j5,j2,j3,j4)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,j5)*zb(j3,j2)*s45*s2345)

      xd2(1,1)=t2(j5,j2,j3,j6)*t2(j7,j1,j5,j4)
     & *zb(j5,j1)/(zb(j3,j5)*zb(j3,j2)*s45*s2367)

      xd3(1,1)=-za(j6,j2)*t2(j7,j2,j6,j3)*zb(j5,j1)**2*za(j4,j1)
     & /(zb(j3,j5)*s45*s267*s2367)

      xd4(1,1)=
     & -za(j4,j2)*t2(j5,j2,j4,j3)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,j5)*s45*s245*s2345)

      xd5(1,1)=za(j6,j2)*t2(j7,j2,j6,j4)*zb(j5,j1)
     & *zb(j5,j1)/(zb(j3,j5)*zb(j3,j1)*s45*s267)

      xd6(1,1)=za(j4,j2)*t2(j5,j2,j4,j6)*zb(j7,j1)
     & *zb(j5,j1)/(zb(j3,j5)*zb(j3,j1)*s45*s245)

      xd7(1,1)=
     & zb(j7,j1)/(zb(j3,j5)*s45*s345*s2345)
     & *t2(j5,j1,j7,j6)*za(j3,j4)
     * *(za(j3,j2)*zb(j3,j5)-za(j4,j2)*zb(j5,j4))


      xd8(1,1)=
     & za(j6,j2)*zb(j5,j1)*za(j3,j4)/(zb(j3,j5)*s45*s345*s267)
     & *(-t2(j7,j2,j6,j3)*zb(j3,j5)+t2(j7,j2,j6,j4)*zb(j5,j4))

      xd9(1,1)=-t2(j5,j3,j4,j2)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,j5)*zb(j3,j4)*s345*s2345)

      xd10(1,1)=czip

      xd11(1,1)=-za(j6,j2)*zb(j5,j1)
     & *(zb(j5,j3)*t2(j7,j2,j6,j3)+zb(j5,j4)*t2(j7,j2,j6,j4))
     & /(zb(j3,j5)*zb(j3,j4)*s345*s267)

      xd12(1,1)=czip

C**********************************************************************

      xd1(2,2)=-za(j5,j2)*za(j5,j2)*t2(j4,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,j5)*za(j3,j2)*s45*s2345)

      xd2(2,2)=-za(j5,j2)*za(j6,j2)*t2(j7,j1,j4,j5)
     & *zb(j4,j1)/(za(j3,j5)*za(j3,j2)*s45*s2367)

      xd3(2,2)=za(j6,j2)*t2(j7,j2,j6,j5)*t2(j3,j1,j4,j5)
     & *zb(j4,j1)/(za(j3,j5)*s45*s267*s2367)

      xd4(2,2)=za(j5,j2)**2*zb(j4,j2)*t2(j3,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,j5)*s45*s245*s2345)

      xd5(2,2)=za(j6,j2)*t2(j7,j2,j6,j5)*t2(j4,j1,j3,j5)
     & *zb(j3,j1)/(za(j3,j5)*s(j1,j3)*s45*s267)

      xd6(2,2)=-za(j5,j2)*t2(j4,j2,j5,j6)*t2(j7,j1,j3,j5)
     & /(za(j3,j5)*za(j3,j1)*s45*s245)

      xd7(2,2)=-zb(j7,j1)*za(j5,j2)*zb(j3,j4)/(za(j3,j5)*s45*s345*s2345)
     &  *(+t2(j3,j1,j7,j6)*za(j3,j5)-t2(j4,j1,j7,j6)*za(j5,j4))

      xd8(2,2)=-za(j6,j2)/(za(j3,j5)*s45*s345*s267)
     & *t2(j7,j2,j6,j5)*zb(j3,j4)
     & *(zb(j3,j1)*za(j5,j3)+zb(j4,j1)*za(j5,j4))

      xd9(2,2)=-za(j5,j2)*zb(j7,j1)
     & *(za(j3,j5)*t2(j3,j1,j7,j6)+za(j4,j5)*t2(j4,j1,j7,j6))
     & /(za(j3,j5)*za(j4,j3)*s345*s2345)

      xd10(2,2)=czip

      xd11(2,2)=za(j6,j2)*t2(j7,j2,j6,j5)
     & *(za(j3,j5)*zb(j3,j1)+za(j4,j5)*zb(j4,j1))
     & /(za(j3,j5)*za(j4,j3)*s345*s267)

      xd12(2,2)=czip


C****************************************************************

      rmsq=0._dp

      do h1=1,2
      do h2=1,2
      xi(h1,h2)=+xd1(h1,h2)+xd2(h1,h2)+xd3(h1,h2)
     &             +xd4(h1,h2)+xd5(h1,h2)+xd6(h1,h2)
      xf(h1,h2)=+xd9(h1,h2)+xd10(h1,h2)+xd11(h1,h2)+xd12(h1,h2)
      xa(h1,h2)=+xd1(h1,h2)+xd2(h1,h2)+xd3(h1,h2)
     &             +xd10(h1,h2)+xd12(h1,h2)-xd7(h1,h2)-xd8(h1,h2)
      xb(h1,h2)=+xd4(h1,h2)+xd5(h1,h2)+xd6(h1,h2)
     &             +xd9(h1,h2)+xd11(h1,h2)+xd7(h1,h2)+xd8(h1,h2)

      rmsq=rmsq
     & +V*xn/eight*(abs(xa(h1,h2))**2+abs(xb(h1,h2))**2)
     & +V/(eight*xn)*(abs(xi(h1,h2))**2+abs(xf(h1,h2))**2
     & -two*(abs(xi(h1,h2)+xf(h1,h2)))**2)

      enddo
      enddo

      prop=((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2)
      rmsq=32._dp*rmsq/prop
      return
      end

