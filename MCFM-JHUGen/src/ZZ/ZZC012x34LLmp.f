      subroutine ZZC012x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, October 2013
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: k1,k2,k3,k4,k5,k6
      integer:: h3,h5,j1,j2,j3,j4,j5,j6,itot
      real(dp):: mt,IDelta,delta,t,s134,s234
      complex(dp):: izab2,zab2,amp0,Xmp(2,2,2),Xpm(2,2,2)
      parameter(itot=1)

C---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      izab2(k1,k2,k3,k4)=cone/zab2(k1,k2,k3,k4)
      delta(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
C---end statement functions

      IDelta=one/(s(j3,j4)**2+s(j1,j2)**2+s(j5,j6)**2
     & -two*(+s(j3,j4)*s(j1,j2)+s(j3,j4)*s(j5,j6)+s(j5,j6)*s(j1,j2)))

      s134=t(j1,j3,j4)
      s234=t(j2,j3,j4)
      k1=j1
      k2=j2

      do h3=1,2
      do h5=1,2
      if (h3 == 1) then
        k3=j3
        k4=j4
      elseif (h3 == 2) then
        k3=j4
        k4=j3
      endif
      if (h5 == 1) then
        k5=j5
        k6=j6
      elseif (h5 == 2) then
        k5=j6
        k6=j5
      endif

      amp0=
     & +six*IDelta**2*s(k1,k2)*zab2(k1,k3,k4,k2)
     &     *zab2(k3,k1,k2,k4)*zab2(k5,k1,k2,k6)/zab2(k2,k3,k4,k1)
     &     *delta(k1,k2,k3,k4,k5,k6)

     & +two*IDelta/zab2(k2,k3,k4,k1)*(
     &  +zab2(k1,k3,k4,k2)*izab2(k2,k3,k4,k1)**2
     &   *(s134-s234)*(za(k2,k5)**2*za(k3,k4)*zb(k1,k4)**2*zb(k5,k6)
     &                -za(k2,k3)**2*za(k5,k6)*zb(k1,k6)**2*zb(k3,k4))
     &  +(s134-s234)*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)
     &   *(za(k2,k3)*zb(k1,k6)*zab2(k5,k1,k3,k4)
     &    +za(k2,k5)*zb(k1,k4)*zab2(k3,k1,k5,k6))
     & -three*zab2(k1,k3,k4,k2)*(za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)
     &                         +za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4))
     &  +za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*zab2(k1,k3,k4,k2)**2
     &   /zab2(k2,k3,k4,k1)
     &  +za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)*zab2(k2,k3,k4,k1))

     & -za(k1,k3)*zb(k2,k6)*zab2(k5,k1,k3,k4)/s134*izab2(k2,k3,k4,k1)
     & +za(k1,k5)*zb(k2,k4)*zab2(k3,k1,k5,k6)/s234*izab2(k2,k3,k4,k1)


c--- This is raw Form output: use numerical relation instead
c      amp2= + ( 2.D0/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k3)*
c     &    za(k1,k5)*zb(k2,k4)*zb(k2,k6) )
c      amp2=amp2 + izab2(k1,k3,k4,k1) * (  - 1/(s(k1,k2))/(s(k5,k6)
c     &    )*za(k1,k3)**2*za(k1,k5)**2*zb(k1,k2)*zb(k5,k6)*iza(k1,k2)*
c     &    iza(k3,k4) )
c      amp2=amp2 + izab2(k1,k3,k4,k1)*izab2(k2,k3,k4,k1) * (  - 1/(
c     &    s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(k1,k5)*za(k3,k4)
c     &    *zb(k1,k4)**2*zb(k2,k6)*zab2(k1,k3,k4,k2) + 2.D0/(s(k1,k2))/(
c     &    s(k5,k6))*za(k1,k3)*za(k1,k5)**2*zb(k1,k2)*zb(k1,k4)*zb(k5,k6
c     &    ) )
c      amp2=amp2 + izab2(k1,k3,k4,k1)*izab2(k2,k3,k4,k1)**2 * ( 2.D0
c     &    /(s(k3,k4))/(s(k5,k6))*za(k1,k5)*za(k3,k4)*zb(k1,k4)**2*zb(k1
c     &    ,k6)*izb(k1,k2)*zab2(k1,k3,k4,k2)*s134 - 2.D0/(s(k5,k6
c     &    ))*za(k1,k5)**2*za(k3,k4)*zb(k1,k4)**2*zb(k5,k6) - 1/(s(k5,k6
c     &    ))*za(k1,k5)*za(k3,k4)*zb(k1,k4)**2*zb(k1,k6)*izb(k1,k2)*
c     &    zab2(k1,k3,k4,k2) )
c      amp2=amp2 + izab2(k1,k5,k6,k1) * (  - 1/(s(k1,k2))/(s(k5,k6)
c     &    )*za(k1,k3)**2*za(k1,k5)**2*zb(k1,k2)*zb(k5,k6)*iza(k1,k2)*
c     &    iza(k3,k4) )
c      amp2=amp2 + izab2(k1,k5,k6,k1)*izab2(k2,k3,k4,k1) * ( 2.D0/(
c     &    s(k1,k2))*za(k1,k3)**2*za(k1,k5)*zb(k1,k2)*zb(k1,k6)*iza(k3,
c     &    k4) + 1/(s(k1,k2))/(s(k3,k4))*za(k1,k2)*za(k1,k3)*zb(k1,k6)**
c     &    2*zb(k2,k4)*izb(k5,k6)*zab2(k1,k3,k4,k2) )
c      amp2=amp2 + izab2(k1,k5,k6,k1)*izab2(k2,k3,k4,k1)**2 * ( 2.D0
c     &    *za(k1,k3)**2*za(k5,k6)*zb(k1,k6)**2*iza(k3,k4) + 1/(s(k3,k4)
c     &    )*za(k1,k3)*za(k5,k6)*zb(k1,k4)*zb(k1,k6)**2*izb(k1,k2)*zab2(
c     &    k1,k3,k4,k2) + 2.D0/(s(k3,k4))*za(k1,k3)*zb(k1,k4)*zb(k1,k6)
c     &    **2*izb(k1,k2)*izb(k5,k6)*zab2(k1,k3,k4,k2)*s234 )
c      amp2=amp2 + izab2(k2,k3,k4,k1) * ( 1/(s(k1,k2))/(s(k3,k4))/(
c     &    s(k5,k6))*za(k1,k2)*za(k2,k3)*zb(k2,k4)*zb(k2,k6)*zab2(k5,k3,
c     &    k4,k2) - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(k2,
c     &    k5)*zb(k2,k4)*zb(k2,k6)*zab2(k3,k5,k6,k2) + 1/(s(k1,k2))/(s(
c     &    k3,k4))/(s(k5,k6))*za(k1,k2)*za(k3,k5)*zb(k2,k4)*zb(k2,k6)*t(
c     &    k1,k3,k4) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(
c     &    k3,k5)*zb(k2,k4)*zb(k2,k6)*s234 - 1/(s(k1,k2))/(s(k3,
c     &    k4))/(s(k5,k6))*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zab2(
c     &    k1,k3,k4,k6) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k3)*
c     &    za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zab2(k1,k5,k6,k4) - 1/(s(k1,k2)
c     &    )/(s(k3,k4))/(s(k5,k6))*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k4,
c     &    k6)*s134 - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k3
c     &    )*za(k1,k5)*zb(k1,k2)*zb(k4,k6)*s234 - 2.D0/(s(k1,k2))
c     &    /(s(k3,k4))/(s(k5,k6))*za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6
c     &    )*zab2(k1,k3,k4,k2) - 2.D0/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*
c     &    za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*zab2(k1,k3,k4,k2) )
c      amp2=amp2 + izab2(k2,k3,k4,k1) * ( 1/(s(k3,k4))*za(k1,k3)*
c     &    zb(k2,k6)*zb(k4,k6)*izb(k5,k6) - 1/(s(k3,k4))/(s(k5,k6))*za(
c     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k4,k6) - 1/(s(k3,k4))/(s(k5,k6)
c     &    )*za(k1,k5)*za(k3,k4)*zb(k2,k4)*zb(k4,k6) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2 * (  - 4.D0*za(k3,k5)*zb(
c     &    k4,k6) - 2.D0*iza(k3,k4)*izb(k5,k6)*zab2(k3,k1,k4,k6)**2 +
c     &    1/(s(k1,k2))/(s(k3,k4))*za(k1,k2)*za(k2,k3)*zb(k1,k6)**2*zb(
c     &    k2,k4)*izb(k5,k6)*zab2(k1,k3,k4,k2) + 1/(s(k1,k2))/(s(k3,k4))
c     &    *za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6)*zab2(k2,k5,k6,k2) +
c     &    1/(s(k1,k2))/(s(k3,k4))*za(k1,k3)*zb(k1,k4)*zb(k1,k6)**2*izb(
c     &    k1,k2)*izb(k5,k6)*zab2(k1,k3,k4,k2)*zab2(k2,k5,k6,k2) + 1/(s(
c     &    k1,k2))/(s(k3,k4))*za(k2,k3)*za(k2,k5)*zb(k2,k4)*zb(k2,k6)*
c     &    zab2(k1,k5,k6,k1) + 1/(s(k1,k2))/(s(k3,k4))*za(k2,k3)*zb(k1,
c     &    k4)*zb(k1,k6)**2*izb(k1,k2)*izb(k5,k6)*zab2(k1,k3,k4,k2)**2
c     &     - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(k1,k3)*za(
c     &    k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*zab2(k1,k3,k4,k2) + 1/(
c     &    s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(k1,k5)*za(k2,k3)
c     &    *zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*zab2(k1,k3,k4,k2) + 1/(s(k1,k2
c     &    ))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(k1,
c     &    k4)**2*zb(k2,k6)*zab2(k1,k3,k4,k2) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2 * (  - 1/(s(k1,k2))/(s(k3,
c     &    k4))/(s(k5,k6))*za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zab2(
c     &    k1,k3,k4,k2)*zab2(k2,k1,k5,k6) + 1/(s(k1,k2))/(s(k3,k4))/(s(
c     &    k5,k6))*za(k1,k3)*zb(k1,k4)*zab2(k1,k3,k4,k6)*zab2(k2,k3,k4,
c     &    k2)*zab2(k5,k3,k4,k1) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*
c     &    za(k1,k5)*za(k2,k3)*zb(k1,k2)*zb(k1,k6)*zab2(k1,k3,k4,k2)*
c     &    zab2(k2,k1,k3,k4) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,
c     &    k5)*za(k3,k4)*zb(k1,k4)**2*zb(k1,k6)*izb(k1,k2)*zab2(k1,k3,k4
c     &    ,k2)*zab2(k2,k3,k4,k2) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*
c     &    za(k1,k5)*zb(k1,k6)*zab2(k1,k5,k6,k4)*zab2(k2,k5,k6,k2)*zab2(
c     &    k3,k5,k6,k1) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k2,k3)*
c     &    za(k2,k5)*zb(k1,k4)*zb(k2,k6)*zab2(k1,k3,k4,k1)*zab2(k1,k3,k4
c     &    ,k2) - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)
c     &    *zb(k1,k4)*iza(k1,k2)*zab2(k1,k3,k4,k2)**2*zab2(k2,k1,k5,k6)
c     &     - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(
c     &    k1,k6)*zb(k2,k4)*zab2(k1,k3,k4,k2)*zab2(k1,k5,k6,k1) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2 * (  - 1/(s(k1,k2))/(s(k3,
c     &    k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(k1,k6)*iza(k1,k2)*
c     &    zab2(k1,k3,k4,k2)**2*zab2(k2,k1,k3,k4) - 1/(s(k1,k2))/(s(k3,
c     &    k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(k2,k4)*iza(k1,k2)*
c     &    zab2(k1,k3,k4,k2)*zab2(k1,k5,k6,k1)*zab2(k2,k1,k5,k6) + 1/(s(
c     &    k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(k2,k6)*
c     &    iza(k1,k2)*zab2(k1,k3,k4,k1)*zab2(k1,k3,k4,k2)*zab2(k2,k1,k3,
c     &    k4) + 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k2,k3)*zb(k2,k4)*
c     &    zab2(k1,k3,k4,k1)*zab2(k2,k3,k4,k6)*zab2(k5,k3,k4,k2) - 1/(s(
c     &    k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k2,k5)*za(k3,k4)*zb(k1,k4)**
c     &    2*zb(k1,k6)*izb(k1,k2)*zab2(k1,k3,k4,k2)**2 + 1/(s(k1,k2))/(
c     &    s(k3,k4))/(s(k5,k6))*za(k2,k5)*zb(k2,k6)*zab2(k1,k5,k6,k1)*
c     &    zab2(k2,k5,k6,k4)*zab2(k3,k5,k6,k2) + 1/(s(k1,k2))/(s(k5,k6))
c     &    *za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6)*zab2(k2,k3,k4,k2) +
c     &    1/(s(k1,k2))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(k2,k4)*zb(k2,
c     &    k6)*zab2(k1,k3,k4,k1) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2 * ( 1/(s(k3,k4))*za(k1,k3)
c     &    *za(k1,k5)*zb(k1,k4)*zb(k1,k6) + 1/(s(k3,k4))*za(k2,k3)*za(k2
c     &    ,k5)*zb(k2,k4)*zb(k2,k6) - 1.D0/2.D0/(s(k3,k4))*za(k2,k3)*zb(
c     &    k1,k6)*zb(k4,k6)*izb(k5,k6)*zab2(k1,k3,k4,k2) - 3.D0/2.D0/(s(
c     &    k3,k4))/(s(k5,k6))*za(k1,k2)*za(k3,k5)*zb(k1,k4)*zb(k1,k6)*
c     &    zab2(k1,k3,k4,k2) - 1/(s(k3,k4))/(s(k5,k6))*za(k1,k3)*za(k2,
c     &    k5)*zb(k1,k4)*zb(k2,k6)*s134 - 1/(s(k3,k4))/(s(k5,k6))
c     &    *za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)*s234 + 1/(s(
c     &    k3,k4))/(s(k5,k6))*za(k1,k3)*zb(k1,k4)*zab2(k1,k3,k4,k6)*
c     &    zab2(k5,k3,k4,k1) - 1/(s(k3,k4))/(s(k5,k6))*za(k1,k5)*za(k2,
c     &    k3)*zb(k1,k6)*zb(k2,k4)*s134 - 1/(s(k3,k4))/(s(k5,k6))
c     &    *za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*s234 + 1/(s(
c     &    k3,k4))/(s(k5,k6))*za(k1,k5)*zb(k1,k6)*zab2(k1,k5,k6,k4)*
c     &    zab2(k3,k5,k6,k1) + 1.D0/2.D0/(s(k3,k4))/(s(k5,k6))*za(k2,k3)
c     &    *za(k2,k5)*zb(k1,k2)*zb(k4,k6)*zab2(k1,k3,k4,k2) + 1/(s(k3,k4
c     &    ))/(s(k5,k6))*za(k2,k3)*zb(k2,k4)*zab2(k2,k3,k4,k6)*zab2(k5,
c     &    k3,k4,k2) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2 * ( 1.D0/2.D0/(s(k3,k4))/(
c     &    s(k5,k6))*za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6)*zab2(k1,k3,
c     &    k4,k2) + 1/(s(k3,k4))/(s(k5,k6))*za(k2,k5)*zb(k2,k6)*zab2(k2,
c     &    k5,k6,k4)*zab2(k3,k5,k6,k2) - 2.D0/(s(k3,k4))/(s(k5,k6))*za(
c     &    k3,k4)*zb(k5,k6)*zab2(k5,k2,k3,k4)**2 - 1.D0/2.D0/(s(k3,k4))
c     &    /(s(k5,k6))*za(k3,k5)*zb(k1,k4)*zab2(k1,k3,k4,k2)*zab2(k2,k1,
c     &    k5,k6) - 1.D0/2.D0/(s(k3,k4))/(s(k5,k6))*za(k3,k5)*zb(k1,k6)*
c     &    zab2(k1,k3,k4,k2)*zab2(k2,k1,k3,k4) + 1/(s(k5,k6))*za(k1,k3)*
c     &    za(k1,k5)*zb(k1,k4)*zb(k1,k6) + 1/(s(k5,k6))*za(k2,k3)*za(k2,
c     &    k5)*zb(k2,k4)*zb(k2,k6) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**3 * (  - 2.D0/(s(k3,k4))/(s(
c     &    k5,k6))*za(k2,k3)*za(k2,k5)*zb(k1,k4)*iza(k1,k2)*zab2(k1,k3,
c     &    k4,k2)*zab2(k2,k1,k5,k6)*s134 - 2.D0/(s(k3,k4))/(s(k5,
c     &    k6))*za(k2,k3)*za(k2,k5)*zb(k1,k6)*iza(k1,k2)*zab2(k1,k3,k4,
c     &    k2)*zab2(k2,k1,k3,k4)*s234 - 2.D0/(s(k3,k4))/(s(k5,k6)
c     &    )*za(k2,k3)*za(k5,k6)*zb(k1,k4)*zb(k1,k6)**2*izb(k1,k2)*zab2(
c     &    k1,k3,k4,k2)*s234 - 2.D0/(s(k3,k4))/(s(k5,k6))*za(k2,
c     &    k5)*za(k3,k4)*zb(k1,k4)**2*zb(k1,k6)*izb(k1,k2)*zab2(k1,k3,k4
c     &    ,k2)*s134 )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2*izab2(k2,k3,k4,k2) * ( 2.D0
c     &    *za(k1,k2)*za(k2,k3)*zb(k1,k4)*zb(k2,k6)**2*izb(k5,k6) + 2.D0
c     &    *za(k2,k3)*zb(k2,k6)**2*izb(k5,k6)*zab2(k2,k1,k3,k4) + 2.D0/(
c     &    s(k3,k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)*
c     &    zab2(k1,k3,k4,k2)*s234 + 2.D0/(s(k3,k4))/(s(k5,k6))*
c     &    za(k2,k3)*za(k2,k5)*zb(k2,k6)*iza(k1,k2)*zab2(k1,k3,k4,k2)*
c     &    zab2(k2,k1,k3,k4)*s234 - 1/(s(k5,k6))*za(k2,k3)*za(k2,
c     &    k5)*zb(k1,k4)*zb(k2,k6)*zab2(k1,k3,k4,k2) - 1/(s(k5,k6))*za(
c     &    k2,k3)*za(k2,k5)*zb(k2,k6)*iza(k1,k2)*zab2(k1,k3,k4,k2)*zab2(
c     &    k2,k1,k3,k4) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)**2*izab2(k2,k5,k6,k2) * (  -
c     &    2.D0/(s(k3,k4))*za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(k1,k6)*zb(k2
c     &    ,k4)**2 + 1/(s(k3,k4))*za(k2,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,k4
c     &    )*zab2(k1,k3,k4,k2) + 1/(s(k3,k4))*za(k2,k3)*za(k2,k5)*zb(k2,
c     &    k4)*iza(k1,k2)*zab2(k1,k3,k4,k2)*zab2(k2,k1,k5,k6) - 2.D0/(s(
c     &    k3,k4))*za(k2,k5)*za(k3,k4)*zb(k2,k4)**2*zab2(k2,k1,k5,k6) -
c     &    2.D0/(s(k3,k4))/(s(k5,k6))*za(k2,k3)*za(k2,k5)*zb(k1,k6)*zb(
c     &    k2,k4)*zab2(k1,k3,k4,k2)*s134 - 2.D0/(s(k3,k4))/(s(k5,
c     &    k6))*za(k2,k3)*za(k2,k5)*zb(k2,k4)*iza(k1,k2)*zab2(k1,k3,k4,
c     &    k2)*zab2(k2,k1,k5,k6)*s134 )
c      amp2=amp2 + izab2(k2,k3,k4,k1)*izab2(k2,k3,k4,k2) * (  - 2.D0
c     &    /(s(k1,k2))*za(k1,k2)*za(k2,k3)*zb(k2,k4)*zb(k2,k6)**2*izb(k5
c     &    ,k6) - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)*za(k1,k5)
c     &    *za(k2,k3)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*zab2(k1,k3,k4,k2) -
c     &    1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k5)*za(k2,k3)*zb(k1,
c     &    k2)*zb(k2,k6)*zab2(k1,k3,k4,k2)*zab2(k2,k1,k3,k4) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)*izab2(k2,k5,k6,k2) * (  - 2.D0
c     &    /(s(k1,k2))/(s(k3,k4))*za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(k2,k4
c     &    )**2*zb(k2,k6) - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k2)
c     &    *za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4)*zab2(k1,k3
c     &    ,k4,k2) - 1/(s(k1,k2))/(s(k3,k4))/(s(k5,k6))*za(k1,k3)*za(k2,
c     &    k5)*zb(k1,k2)*zb(k2,k4)*zab2(k1,k3,k4,k2)*zab2(k2,k1,k5,k6) )
c      amp2=amp2 + izab2(k2,k3,k4,k1)*IDelta * (  - 8.D0*za(k3,k5)*
c     &    zb(k4,k6)*zab2(k1,k3,k4,k2) + 2.D0/(s(k3,k4))*s(k1,k2)*za(k3,
c     &    k4)*zb(k4,k6)**2*izb(k5,k6)*zab2(k1,k3,k4,k2) - 2.D0/(s(k3,k4
c     &    ))*s(k3,k4)*za(k3,k4)*zb(k4,k6)**2*izb(k5,k6)*zab2(k1,k3,k4,
c     &    k2) - 2.D0/(s(k3,k4))*s(k5,k6)*za(k3,k4)*zb(k4,k6)**2*izb(k5,
c     &    k6)*zab2(k1,k3,k4,k2) + 2.D0/(s(k5,k6))*s(k1,k2)*za(k3,k5)**2
c     &    *zb(k5,k6)*iza(k3,k4)*zab2(k1,k3,k4,k2) - 2.D0/(s(k5,k6))*s(
c     &    k3,k4)*za(k3,k5)**2*zb(k5,k6)*iza(k3,k4)*zab2(k1,k3,k4,k2) -
c     &    2.D0/(s(k5,k6))*s(k5,k6)*za(k3,k5)**2*zb(k5,k6)*iza(k3,k4)*
c     &    zab2(k1,k3,k4,k2) )
c      amp2=amp2 + izab2(k2,k3,k4,k2) * (  - 1/(s(k1,k2))/(s(k3,k4)
c     &    )*za(k1,k2)*za(k3,k4)*zb(k2,k4)**2*zb(k2,k6)**2*izb(k1,k2)*
c     &    izb(k5,k6) )
c      amp2=amp2 + izab2(k2,k5,k6,k2) * (  - 1/(s(k1,k2))/(s(k3,k4)
c     &    )*za(k1,k2)*za(k3,k4)*zb(k2,k4)**2*zb(k2,k6)**2*izb(k1,k2)*
c     &    izb(k5,k6) )

c      Xmp(h3,h5,itot)=amp0+amp2*mtsq
c      Xmp(h3,h5,irat)=amp2*(+half)

c--- Only return the massless coefficient
      Xmp(h3,h5,itot)=amp0

      enddo
      enddo

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5,:)=conjg(Xmp(3-h3,3-h5,:))
      enddo
      enddo

      return
      end
