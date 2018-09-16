      subroutine upper_parttri(p,upper_tri,first)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'alpha1.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'nwz.f'
      include 'decl_kininv.f'
      integer:: k1,k2,ep,eta
      complex(dp):: prW,prq,upper_tri(2,2,-2:0),iprZ,iza,izb
      complex(dp):: vert25x1,vert25x2,vert25x3,
     & vert16x2,vert16x3,vert16x4,vert16x5,vert16x6,vert16x7,
     & vert16x8,vert16x9,vert16x10,vert16x11,vert16x12,vert16x13
      complex(dp):: cprop,facuLl,facdLl
      real(dp):: p(mxpart,4),q(mxpart,4),mtsq,mwsq
      real(dp):: p5Deta,omal
      integer:: j3,p1,p2,p3,p4,k5,e5,p6
      integer:: epmin
      logical:: first
      parameter(p1=1,p2=2,k5=5,p6=7,e5=6)

c----statement function
      include 'cplx.h'
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      prq(s16)=cone/cplx1(s16)
      iza(k1,k2)=cone/za(k1,k2)
      izb(k1,k2)=cone/zb(k1,k2)
c----end statement function

      omal=1d0-alpha1

C---choose auxiliary vector
      eta=2
      mtsq=mt**2
      mwsq=wmass**2
      p5Deta=p(5,4)*p(eta,4)
     &  -p(5,1)*p(eta,1)-p(5,2)*p(eta,2)-p(5,3)*p(eta,3)
      q(:,:)=p(:,:)
      q(7,:)=p(6,:)
      q(5,:)=p(5,:)-0.5d0*mtsq*p(eta,:)/p5Deta
      q(e5,:)=0.5d0*mtsq*p(eta,:)/p5Deta

      if (nwz == +1) then
      call spinoru(7,q,za,zb)
      elseif (nwz == -1) then
      call spinoru(7,q,zb,za)
      endif


c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=cplx1(1d0/sqrt((s34-zmass**2)**2+(zmass*zwidth)**2))
      iprZ=cplx1(s34-zmass**2)
      if (first) then
         epmin=-2
      else
         epmin=-1
      endif

c      write(*,*) 'epmin in upper_tri', epmin

      do ep=epmin,0
      do j3=1,2
      if (j3 == 1) then
        p3=3
        p4=4
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*le)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*le)
      elseif (j3 == 2) then
        p3=4
        p4=3
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*re)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*re)
      endif
 
      call upper_vertices(mtsq,ep,facuLl,facdLl,
     & vert25x1,vert25x2,vert25x3,vert16x2,vert16x3,vert16x4,
     & vert16x5,vert16x6,vert16x7,vert16x8,vert16x9,vert16x10,
     & vert16x11,vert16x12,vert16x13)
      upper_tri(j3,1,ep)= + prW(s25)*vert25x1*facdLl*s346**(-1) * (  - 
     &    2.D0*za(p3,p6)*za(p3,k5)*zb(p1,p2)*zb(p3,p4) + 2.D0*za(p3,p6)
     &    *za(p6,k5)*zb(p1,p2)*zb(p4,p6) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x1*
     & facuLl*s134**(-1) * (  - 2.D0*za(p1,p3)*za(p6,k5)*zb(p1,p2)*zb(
     &    p1,p4) - 2.D0*za(p3,p4)*za(p6,k5)*zb(p1,p4)*zb(p2,p4) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x1*omal*
     & mwsq**(-1)*facdLl * ( za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    mtsq )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x1*omal*
     & mwsq**(-1)*facuLl * (  - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5
     &    )*mtsq )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x2*
     & facdLl*mt*s346**(-1) * ( za(p2,p3)*za(p3,p6)*zb(p1,p2)*zb(p2,e5)
     &    *zb(p3,p4)*izb(k5,e5) - za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,
     &    e5)*zb(p4,p6)*izb(k5,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x2*
     & facuLl*mt*s134**(-1) * ( za(p1,p3)*za(p2,p6)*zb(p1,p2)*zb(p1,p4)
     &    *zb(p2,e5)*izb(k5,e5) + za(p2,p6)*za(p3,p4)*zb(p1,p4)*zb(p2,
     &    p4)*zb(p2,e5)*izb(k5,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x2*omal*
     & mwsq**(-1)*facdLl*mt * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)
     &    *izb(k5,e5)*s25 - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*
     &    izb(k5,e5)*mtsq )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x2*omal*
     & mwsq**(-1)*facuLl*mt * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,
     &    e5)*izb(k5,e5)*s25 + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*
     &    izb(k5,e5)*mtsq )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x3*
     & facdLl*mt * (  - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x3*
     & facuLl*mt * ( za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x3*omal*
     & mwsq**(-1)*facdLl*mt * ( za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5
     &    )*s25 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert25x3*omal*
     & mwsq**(-1)*facuLl*mt * (  - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5
     &    ,e5)*s25 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x2*
     & facdLl * (  - za(p1,k5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x2*omal*
     & mwsq**(-1)*facdLl * (  - 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(p1,p4)
     &    *zb(p2,e5)*s25 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x2*omal*
     & mwsq**(-1)*facdLl*s346 * ( 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(p1,
     &    p4)*zb(p2,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x3*
     & facdLl * ( za(p1,k5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4) + za(p3,p6)*
     &    za(k5,e5)*zb(p1,p4)*zb(p2,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x3*omal*
     & mwsq**(-1)*facdLl * (  - 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(p1,p4)
     &    *zb(p2,e5)*s25 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x3*omal*
     & mwsq**(-1)*facdLl*s346 * (  - 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(
     &    p1,p4)*zb(p2,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x4*
     & facdLl*s346**(-1) * (  - 2.D0*za(p3,p6)*za(p3,k5)*zb(p1,p2)*zb(
     &    p3,p4) + 2.D0*za(p3,p6)*za(p6,k5)*zb(p1,p2)*zb(p4,p6) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x4*omal*
     & mwsq**(-1)*facdLl*s346**(-1) * ( za(p3,p2)*za(p3,p6)*za(k5,e5)*
     &    zb(p2,p1)*zb(p2,e5)*zb(p3,p4) + za(p3,p6)*za(p3,k5)*za(k5,e5)
     &    *zb(p2,e5)*zb(p3,p4)*zb(k5,p1) + za(p3,p6)*za(p3,e5)*za(k5,e5
     &    )*zb(p2,e5)*zb(p3,p4)*zb(e5,p1) - za(p3,p6)*za(p6,p2)*za(k5,
     &    e5)*zb(p2,p1)*zb(p2,e5)*zb(p4,p6) - za(p3,p6)*za(p6,k5)*za(k5
     &    ,e5)*zb(p2,e5)*zb(p4,p6)*zb(k5,p1) - za(p3,p6)*za(p6,e5)*za(
     &    k5,e5)*zb(p2,e5)*zb(p4,p6)*zb(e5,p1) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x5*
     & s346**(-1) * (  - za(p3,p6)**2*za(p4,k5)*zb(p1,p2)*zb(p3,p4)*zb(
     &    p4,p6) - za(p3,p6)**2*za(p6,k5)*zb(p1,p2)*zb(p3,p6)*zb(p4,p6)
     &     + za(p3,p6)*za(p3,k5)*za(p4,p6)*zb(p1,p2)*zb(p3,p4)*zb(p4,p6
     &    ) - za(p3,p6)*za(p4,p6)*za(p6,k5)*zb(p1,p2)*zb(p4,p6)**2 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x5*omal*
     & mwsq**(-1)*s346**(-1) * (  - 1.D0/2.D0*za(p3,p2)*za(p3,p6)*za(p4
     &    ,p6)*za(k5,e5)*zb(p2,p1)*zb(p2,e5)*zb(p3,p4)*zb(p4,p6) + 1.D0/
     &    2.D0*za(p3,p6)**2*za(p4,p2)*za(k5,e5)*zb(p2,p1)*zb(p2,e5)*zb(
     &    p3,p4)*zb(p4,p6) + 1.D0/2.D0*za(p3,p6)**2*za(p4,k5)*za(k5,e5)
     &    *zb(p2,e5)*zb(p3,p4)*zb(p4,p6)*zb(k5,p1) + 1.D0/2.D0*za(p3,p6
     &    )**2*za(p4,e5)*za(k5,e5)*zb(p2,e5)*zb(p3,p4)*zb(p4,p6)*zb(e5,
     &    p1) + 1.D0/2.D0*za(p3,p6)**2*za(p6,p2)*za(k5,e5)*zb(p2,p1)*
     &    zb(p2,e5)*zb(p3,p6)*zb(p4,p6) + 1.D0/2.D0*za(p3,p6)**2*za(p6,
     &    k5)*za(k5,e5)*zb(p2,e5)*zb(p3,p6)*zb(p4,p6)*zb(k5,p1) + 1.D0/
     &    2.D0*za(p3,p6)**2*za(p6,e5)*za(k5,e5)*zb(p2,e5)*zb(p3,p6)*zb(
     &    p4,p6)*zb(e5,p1) - 1.D0/2.D0*za(p3,p6)*za(p3,k5)*za(p4,p6)*
     &    za(k5,e5)*zb(p2,e5)*zb(p3,p4)*zb(p4,p6)*zb(k5,p1) - 1.D0/2.D0
     &    *za(p3,p6)*za(p3,e5)*za(p4,p6)*za(k5,e5)*zb(p2,e5)*zb(p3,p4)*
     &    zb(p4,p6)*zb(e5,p1) + 1.D0/2.D0*za(p3,p6)*za(p4,p6)*za(p6,p2)
     &    *za(k5,e5)*zb(p2,p1)*zb(p2,e5)*zb(p4,p6)**2 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x5*omal*
     & mwsq**(-1)*s346**(-1) * ( 1.D0/2.D0*za(p3,p6)*za(p4,p6)*za(p6,k5
     &    )*za(k5,e5)*zb(p2,e5)*zb(p4,p6)**2*zb(k5,p1) + 1.D0/2.D0*za(
     &    p3,p6)*za(p4,p6)*za(p6,e5)*za(k5,e5)*zb(p2,e5)*zb(p4,p6)**2*
     &    zb(e5,p1) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x7*
     & s346**(-1) * (  - 2.D0*za(p3,p6)*za(p3,k5)*zb(p1,p2)*zb(p3,p4)
     &     + 2.D0*za(p3,p6)*za(p6,k5)*zb(p1,p2)*zb(p4,p6) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x7*omal*
     & mwsq**(-1)*s346**(-1) * ( za(p3,p2)*za(p3,p6)*za(k5,e5)*zb(p2,p1
     &    )*zb(p2,e5)*zb(p3,p4) + za(p3,p6)*za(p3,k5)*za(k5,e5)*zb(p2,
     &    e5)*zb(p3,p4)*zb(k5,p1) + za(p3,p6)*za(p3,e5)*za(k5,e5)*zb(p2
     &    ,e5)*zb(p3,p4)*zb(e5,p1) - za(p3,p6)*za(p6,p2)*za(k5,e5)*zb(
     &    p2,p1)*zb(p2,e5)*zb(p4,p6) - za(p3,p6)*za(p6,k5)*za(k5,e5)*
     &    zb(p2,e5)*zb(p4,p6)*zb(k5,p1) - za(p3,p6)*za(p6,e5)*za(k5,e5)
     &    *zb(p2,e5)*zb(p4,p6)*zb(e5,p1) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x8*
     & s134**(-1) * ( za(p1,p3)**2*za(p6,k5)*zb(p1,p2)*zb(p1,p3)*zb(p1,
     &    p4) + za(p1,p3)*za(p1,p4)*za(p6,k5)*zb(p1,p2)*zb(p1,p4)**2 + 
     &    za(p1,p3)*za(p3,p4)*za(p6,k5)*zb(p1,p3)*zb(p1,p4)*zb(p2,p4)
     &     - za(p1,p3)*za(p3,p4)*za(p6,k5)*zb(p1,p4)**2*zb(p2,p3) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x8*omal*
     & mwsq**(-1)*s134**(-1) * (  - 1.D0/2.D0*za(p1,p3)**2*za(p6,p2)*
     &    za(k5,e5)*zb(p1,p3)*zb(p1,p4)*zb(p2,p1)*zb(p2,e5) - 1.D0/2.D0
     &    *za(p1,p3)**2*za(p6,k5)*za(k5,e5)*zb(p1,p3)*zb(p1,p4)*zb(p2,
     &    e5)*zb(k5,p1) - 1.D0/2.D0*za(p1,p3)**2*za(p6,e5)*za(k5,e5)*
     &    zb(p1,p3)*zb(p1,p4)*zb(p2,e5)*zb(e5,p1) - 1.D0/2.D0*za(p1,p3)
     &    *za(p1,p4)*za(p6,p2)*za(k5,e5)*zb(p1,p4)**2*zb(p2,p1)*zb(p2,
     &    e5) - 1.D0/2.D0*za(p1,p3)*za(p1,p4)*za(p6,k5)*za(k5,e5)*zb(p1
     &    ,p4)**2*zb(p2,e5)*zb(k5,p1) - 1.D0/2.D0*za(p1,p3)*za(p1,p4)*
     &    za(p6,e5)*za(k5,e5)*zb(p1,p4)**2*zb(p2,e5)*zb(e5,p1) + 1.D0/2.
     &    D0*za(p1,p3)*za(p3,p4)*za(p6,p2)*za(k5,e5)*zb(p1,p3)*zb(p1,p4
     &    )*zb(p2,p4)*zb(p2,e5) - 1.D0/2.D0*za(p1,p3)*za(p3,p4)*za(p6,
     &    p2)*za(k5,e5)*zb(p1,p4)**2*zb(p2,p3)*zb(p2,e5) + 1.D0/2.D0*
     &    za(p1,p3)*za(p3,p4)*za(p6,k5)*za(k5,e5)*zb(p1,p3)*zb(p1,p4)*
     &    zb(p2,e5)*zb(k5,p4) - 1.D0/2.D0*za(p1,p3)*za(p3,p4)*za(p6,k5)
     &    *za(k5,e5)*zb(p1,p4)**2*zb(p2,e5)*zb(k5,p3) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x8*omal*
     & mwsq**(-1)*s134**(-1) * ( 1.D0/2.D0*za(p1,p3)*za(p3,p4)*za(p6,e5
     &    )*za(k5,e5)*zb(p1,p3)*zb(p1,p4)*zb(p2,e5)*zb(e5,p4) - 1.D0/2.D
     &    0*za(p1,p3)*za(p3,p4)*za(p6,e5)*za(k5,e5)*zb(p1,p4)**2*zb(p2,
     &    e5)*zb(e5,p3) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x10*
     & s134**(-1) * (  - 2.D0*za(p1,p3)*za(p6,k5)*zb(p1,p2)*zb(p1,p4)
     &     - 2.D0*za(p3,p4)*za(p6,k5)*zb(p1,p4)*zb(p2,p4) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x10*omal
     & *mwsq**(-1)*s134**(-1) * ( za(p1,p3)*za(p6,p2)*za(k5,e5)*zb(p1,
     &    p4)*zb(p2,p1)*zb(p2,e5) + za(p1,p3)*za(p6,k5)*za(k5,e5)*zb(p1
     &    ,p4)*zb(p2,e5)*zb(k5,p1) + za(p1,p3)*za(p6,e5)*za(k5,e5)*zb(
     &    p1,p4)*zb(p2,e5)*zb(e5,p1) - za(p3,p4)*za(p6,p2)*za(k5,e5)*
     &    zb(p1,p4)*zb(p2,p4)*zb(p2,e5) - za(p3,p4)*za(p6,k5)*za(k5,e5)
     &    *zb(p1,p4)*zb(p2,e5)*zb(k5,p4) - za(p3,p4)*za(p6,e5)*za(k5,e5
     &    )*zb(p1,p4)*zb(p2,e5)*zb(e5,p4) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x11*
     & facuLl * (  - za(p3,p6)*za(p6,k5)*zb(p1,p4)*zb(p2,p6) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x11*omal
     & *mwsq**(-1)*facuLl * ( 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(p1,p4)*
     &    zb(p2,e5)*s25 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x11*omal
     & *mwsq**(-1)*facuLl*s134 * (  - 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(
     &    p1,p4)*zb(p2,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x12*
     & facuLl * ( za(p1,k5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4) - za(p3,p6)*
     &    za(p3,k5)*zb(p1,p4)*zb(p2,p3) - za(p3,p6)*za(p4,k5)*zb(p1,p4)
     &    *zb(p2,p4) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x12*omal
     & *mwsq**(-1)*facuLl * ( 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(p1,p4)*
     &    zb(p2,e5)*s25 )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x12*omal
     & *mwsq**(-1)*facuLl*s134 * ( 1.D0/2.D0*za(p3,p6)*za(k5,e5)*zb(p1,
     &    p4)*zb(p2,e5) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x13*
     & facuLl*s134**(-1) * (  - 2.D0*za(p1,p3)*za(p6,k5)*zb(p1,p2)*zb(
     &    p1,p4) - 2.D0*za(p3,p4)*za(p6,k5)*zb(p1,p4)*zb(p2,p4) )
      upper_tri(j3,1,ep) = upper_tri(j3,1,ep) + prW(s25)*vert16x13*omal
     & *mwsq**(-1)*facuLl*s134**(-1) * ( za(p1,p3)*za(p6,p2)*za(k5,e5)*
     &    zb(p1,p4)*zb(p2,p1)*zb(p2,e5) + za(p1,p3)*za(p6,k5)*za(k5,e5)
     &    *zb(p1,p4)*zb(p2,e5)*zb(k5,p1) + za(p1,p3)*za(p6,e5)*za(k5,e5
     &    )*zb(p1,p4)*zb(p2,e5)*zb(e5,p1) - za(p3,p4)*za(p6,p2)*za(k5,
     &    e5)*zb(p1,p4)*zb(p2,p4)*zb(p2,e5) - za(p3,p4)*za(p6,k5)*za(k5
     &    ,e5)*zb(p1,p4)*zb(p2,e5)*zb(k5,p4) - za(p3,p4)*za(p6,e5)*za(
     &    k5,e5)*zb(p1,p4)*zb(p2,e5)*zb(e5,p4) )

      upper_tri(j3,2,ep)= + prW(s25)*vert25x1*facdLl*mt*s346**(-1) * ( 
     &    2.D0*za(p3,p6)*za(p3,e5)*zb(p1,p2)*zb(p3,p4)*iza(k5,e5) - 2.D0
     &    *za(p3,p6)*za(p6,e5)*zb(p1,p2)*zb(p4,p6)*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x1*
     & facuLl*mt*s134**(-1) * ( 2.D0*za(p1,p3)*za(p6,e5)*zb(p1,p2)*zb(
     &    p1,p4)*iza(k5,e5) + 2.D0*za(p3,p4)*za(p6,e5)*zb(p1,p4)*zb(p2,
     &    p4)*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x1*omal*
     & mwsq**(-1)*facdLl*mt * (  - za(p3,p6)*zb(p1,p4)*zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x1*omal*
     & mwsq**(-1)*facuLl*mt * ( za(p3,p6)*zb(p1,p4)*zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x2*
     & facdLl*s346**(-1) * (  - za(p2,p3)*za(p3,p6)*zb(p1,p2)*zb(p2,k5)
     &    *zb(p3,p4) + za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,k5)*zb(p4,p6
     &    ) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x2*
     & facuLl*s134**(-1) * (  - za(p1,p3)*za(p2,p6)*zb(p1,p2)*zb(p1,p4)
     &    *zb(p2,k5) - za(p2,p6)*za(p3,p4)*zb(p1,p4)*zb(p2,p4)*zb(p2,k5
     &    ) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x2*omal*
     & mwsq**(-1)*facdLl * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)
     &    *s25 + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*mtsq )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x2*omal*
     & mwsq**(-1)*facuLl * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*
     &    s25 - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*mtsq )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x3*
     & facdLl * ( za(p3,p6)*zb(p1,p4)*zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x3*
     & facuLl * (  - za(p3,p6)*zb(p1,p4)*zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x3*omal*
     & mwsq**(-1)*facdLl * (  - za(p3,p6)*zb(p1,p4)*zb(p2,k5)*s25 )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert25x3*omal*
     & mwsq**(-1)*facuLl * ( za(p3,p6)*zb(p1,p4)*zb(p2,k5)*s25 )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x2*
     & facdLl*mt * ( za(p1,e5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*iza(k5,e5)
     &     )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x2*omal*
     & mwsq**(-1)*facdLl*mt * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,
     &    k5)*s25 )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x2*omal*
     & mwsq**(-1)*facdLl*mt*s346 * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(
     &    p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x3*
     & facdLl*mt * (  - za(p1,e5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*iza(k5,
     &    e5) + za(p3,p6)*zb(p1,p4)*zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x3*omal*
     & mwsq**(-1)*facdLl*mt * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,
     &    k5)*s25 )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x3*omal*
     & mwsq**(-1)*facdLl*mt*s346 * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*
     &    zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x4*
     & facdLl*mt*s346**(-1) * ( 2.D0*za(p3,p6)*za(p3,e5)*zb(p1,p2)*zb(
     &    p3,p4)*iza(k5,e5) - 2.D0*za(p3,p6)*za(p6,e5)*zb(p1,p2)*zb(p4,
     &    p6)*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x4*omal*
     & mwsq**(-1)*facdLl*mt*s346**(-1) * ( za(p3,p2)*za(p3,p6)*zb(p2,p1
     &    )*zb(p2,k5)*zb(p3,p4) + za(p3,p6)*za(p3,k5)*zb(p2,k5)*zb(p3,
     &    p4)*zb(k5,p1) + za(p3,p6)*za(p3,e5)*zb(p2,k5)*zb(p3,p4)*zb(e5
     &    ,p1) - za(p3,p6)*za(p6,p2)*zb(p2,p1)*zb(p2,k5)*zb(p4,p6) - 
     &    za(p3,p6)*za(p6,k5)*zb(p2,k5)*zb(p4,p6)*zb(k5,p1) - za(p3,p6)
     &    *za(p6,e5)*zb(p2,k5)*zb(p4,p6)*zb(e5,p1) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x5*mt*
     & s346**(-1) * ( za(p3,p6)**2*za(p4,e5)*zb(p1,p2)*zb(p3,p4)*zb(p4,
     &    p6)*iza(k5,e5) + za(p3,p6)**2*za(p6,e5)*zb(p1,p2)*zb(p3,p6)*
     &    zb(p4,p6)*iza(k5,e5) - za(p3,p6)*za(p3,e5)*za(p4,p6)*zb(p1,p2
     &    )*zb(p3,p4)*zb(p4,p6)*iza(k5,e5) + za(p3,p6)*za(p4,p6)*za(p6,
     &    e5)*zb(p1,p2)*zb(p4,p6)**2*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x5*omal*
     & mwsq**(-1)*mt*s346**(-1) * (  - 1.D0/2.D0*za(p3,p2)*za(p3,p6)*
     &    za(p4,p6)*zb(p2,p1)*zb(p2,k5)*zb(p3,p4)*zb(p4,p6) + 1.D0/2.D0
     &    *za(p3,p6)**2*za(p4,p2)*zb(p2,p1)*zb(p2,k5)*zb(p3,p4)*zb(p4,
     &    p6) + 1.D0/2.D0*za(p3,p6)**2*za(p4,k5)*zb(p2,k5)*zb(p3,p4)*
     &    zb(p4,p6)*zb(k5,p1) + 1.D0/2.D0*za(p3,p6)**2*za(p4,e5)*zb(p2,
     &    k5)*zb(p3,p4)*zb(p4,p6)*zb(e5,p1) + 1.D0/2.D0*za(p3,p6)**2*
     &    za(p6,p2)*zb(p2,p1)*zb(p2,k5)*zb(p3,p6)*zb(p4,p6) + 1.D0/2.D0
     &    *za(p3,p6)**2*za(p6,k5)*zb(p2,k5)*zb(p3,p6)*zb(p4,p6)*zb(k5,
     &    p1) + 1.D0/2.D0*za(p3,p6)**2*za(p6,e5)*zb(p2,k5)*zb(p3,p6)*
     &    zb(p4,p6)*zb(e5,p1) - 1.D0/2.D0*za(p3,p6)*za(p3,k5)*za(p4,p6)
     &    *zb(p2,k5)*zb(p3,p4)*zb(p4,p6)*zb(k5,p1) - 1.D0/2.D0*za(p3,p6
     &    )*za(p3,e5)*za(p4,p6)*zb(p2,k5)*zb(p3,p4)*zb(p4,p6)*zb(e5,p1)
     &     + 1.D0/2.D0*za(p3,p6)*za(p4,p6)*za(p6,p2)*zb(p2,p1)*zb(p2,k5
     &    )*zb(p4,p6)**2 + 1.D0/2.D0*za(p3,p6)*za(p4,p6)*za(p6,k5)*zb(
     &    p2,k5)*zb(p4,p6)**2*zb(k5,p1) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x5*omal*
     & mwsq**(-1)*mt*s346**(-1) * ( 1.D0/2.D0*za(p3,p6)*za(p4,p6)*za(p6
     &    ,e5)*zb(p2,k5)*zb(p4,p6)**2*zb(e5,p1) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x7*mt*
     & s346**(-1) * ( 2.D0*za(p3,p6)*za(p3,e5)*zb(p1,p2)*zb(p3,p4)*iza(
     &    k5,e5) - 2.D0*za(p3,p6)*za(p6,e5)*zb(p1,p2)*zb(p4,p6)*iza(k5,
     &    e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x7*omal*
     & mwsq**(-1)*mt*s346**(-1) * ( za(p3,p2)*za(p3,p6)*zb(p2,p1)*zb(p2
     &    ,k5)*zb(p3,p4) + za(p3,p6)*za(p3,k5)*zb(p2,k5)*zb(p3,p4)*zb(
     &    k5,p1) + za(p3,p6)*za(p3,e5)*zb(p2,k5)*zb(p3,p4)*zb(e5,p1) - 
     &    za(p3,p6)*za(p6,p2)*zb(p2,p1)*zb(p2,k5)*zb(p4,p6) - za(p3,p6)
     &    *za(p6,k5)*zb(p2,k5)*zb(p4,p6)*zb(k5,p1) - za(p3,p6)*za(p6,e5
     &    )*zb(p2,k5)*zb(p4,p6)*zb(e5,p1) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x8*mt*
     & s134**(-1) * (  - za(p1,p3)**2*za(p6,e5)*zb(p1,p2)*zb(p1,p3)*zb(
     &    p1,p4)*iza(k5,e5) - za(p1,p3)*za(p1,p4)*za(p6,e5)*zb(p1,p2)*
     &    zb(p1,p4)**2*iza(k5,e5) - za(p1,p3)*za(p3,p4)*za(p6,e5)*zb(p1
     &    ,p3)*zb(p1,p4)*zb(p2,p4)*iza(k5,e5) + za(p1,p3)*za(p3,p4)*za(
     &    p6,e5)*zb(p1,p4)**2*zb(p2,p3)*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x8*omal*
     & mwsq**(-1)*mt*s134**(-1) * (  - 1.D0/2.D0*za(p1,p3)**2*za(p6,p2)
     &    *zb(p1,p3)*zb(p1,p4)*zb(p2,p1)*zb(p2,k5) - 1.D0/2.D0*za(p1,p3
     &    )**2*za(p6,k5)*zb(p1,p3)*zb(p1,p4)*zb(p2,k5)*zb(k5,p1) - 1.D0/
     &    2.D0*za(p1,p3)**2*za(p6,e5)*zb(p1,p3)*zb(p1,p4)*zb(p2,k5)*zb(
     &    e5,p1) - 1.D0/2.D0*za(p1,p3)*za(p1,p4)*za(p6,p2)*zb(p1,p4)**2
     &    *zb(p2,p1)*zb(p2,k5) - 1.D0/2.D0*za(p1,p3)*za(p1,p4)*za(p6,k5
     &    )*zb(p1,p4)**2*zb(p2,k5)*zb(k5,p1) - 1.D0/2.D0*za(p1,p3)*za(
     &    p1,p4)*za(p6,e5)*zb(p1,p4)**2*zb(p2,k5)*zb(e5,p1) + 1.D0/2.D0
     &    *za(p1,p3)*za(p3,p4)*za(p6,p2)*zb(p1,p3)*zb(p1,p4)*zb(p2,p4)*
     &    zb(p2,k5) - 1.D0/2.D0*za(p1,p3)*za(p3,p4)*za(p6,p2)*zb(p1,p4)
     &    **2*zb(p2,p3)*zb(p2,k5) + 1.D0/2.D0*za(p1,p3)*za(p3,p4)*za(p6
     &    ,k5)*zb(p1,p3)*zb(p1,p4)*zb(p2,k5)*zb(k5,p4) - 1.D0/2.D0*za(
     &    p1,p3)*za(p3,p4)*za(p6,k5)*zb(p1,p4)**2*zb(p2,k5)*zb(k5,p3)
     &     + 1.D0/2.D0*za(p1,p3)*za(p3,p4)*za(p6,e5)*zb(p1,p3)*zb(p1,p4
     &    )*zb(p2,k5)*zb(e5,p4) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x8*omal*
     & mwsq**(-1)*mt*s134**(-1) * (  - 1.D0/2.D0*za(p1,p3)*za(p3,p4)*
     &    za(p6,e5)*zb(p1,p4)**2*zb(p2,k5)*zb(e5,p3) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x10*mt*
     & s134**(-1) * ( 2.D0*za(p1,p3)*za(p6,e5)*zb(p1,p2)*zb(p1,p4)*iza(
     &    k5,e5) + 2.D0*za(p3,p4)*za(p6,e5)*zb(p1,p4)*zb(p2,p4)*iza(k5,
     &    e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x10*omal
     & *mwsq**(-1)*mt*s134**(-1) * ( za(p1,p3)*za(p6,p2)*zb(p1,p4)*zb(
     &    p2,p1)*zb(p2,k5) + za(p1,p3)*za(p6,k5)*zb(p1,p4)*zb(p2,k5)*
     &    zb(k5,p1) + za(p1,p3)*za(p6,e5)*zb(p1,p4)*zb(p2,k5)*zb(e5,p1)
     &     - za(p3,p4)*za(p6,p2)*zb(p1,p4)*zb(p2,p4)*zb(p2,k5) - za(p3,
     &    p4)*za(p6,k5)*zb(p1,p4)*zb(p2,k5)*zb(k5,p4) - za(p3,p4)*za(p6
     &    ,e5)*zb(p1,p4)*zb(p2,k5)*zb(e5,p4) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x11*
     & facuLl*mt * ( za(p3,p6)*za(p6,e5)*zb(p1,p4)*zb(p2,p6)*iza(k5,e5)
     &     )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x11*omal
     & *mwsq**(-1)*facuLl*mt * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5
     &    )*s25 )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x11*omal
     & *mwsq**(-1)*facuLl*mt*s134 * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*
     &    zb(p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x12*
     & facuLl*mt * (  - za(p1,e5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*iza(k5,
     &    e5) + za(p3,p6)*za(p3,e5)*zb(p1,p4)*zb(p2,p3)*iza(k5,e5) + 
     &    za(p3,p6)*za(p4,e5)*zb(p1,p4)*zb(p2,p4)*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x12*omal
     & *mwsq**(-1)*facuLl*mt * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5
     &    )*s25 )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x12*omal
     & *mwsq**(-1)*facuLl*mt*s134 * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(
     &    p2,k5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x13*
     & facuLl*mt*s134**(-1) * ( 2.D0*za(p1,p3)*za(p6,e5)*zb(p1,p2)*zb(
     &    p1,p4)*iza(k5,e5) + 2.D0*za(p3,p4)*za(p6,e5)*zb(p1,p4)*zb(p2,
     &    p4)*iza(k5,e5) )
      upper_tri(j3,2,ep) = upper_tri(j3,2,ep) + prW(s25)*vert16x13*omal
     & *mwsq**(-1)*facuLl*mt*s134**(-1) * ( za(p1,p3)*za(p6,p2)*zb(p1,
     &    p4)*zb(p2,p1)*zb(p2,k5) + za(p1,p3)*za(p6,k5)*zb(p1,p4)*zb(p2
     &    ,k5)*zb(k5,p1) + za(p1,p3)*za(p6,e5)*zb(p1,p4)*zb(p2,k5)*zb(
     &    e5,p1) - za(p3,p4)*za(p6,p2)*zb(p1,p4)*zb(p2,p4)*zb(p2,k5) - 
     &    za(p3,p4)*za(p6,k5)*zb(p1,p4)*zb(p2,k5)*zb(k5,p4) - za(p3,p4)
     &    *za(p6,e5)*zb(p1,p4)*zb(p2,k5)*zb(e5,p4) )


      enddo
      enddo

      upper_tri=upper_tri*cprop

      return
      end




      subroutine upper_vertices(mtsq,ep,
     & facuLl,facdLl,vert25x1,vert25x2,vert25x3,
     & vert16x2,vert16x3,vert16x4,vert16x5,vert16x6,vert16x7,
     & vert16x8,vert16x9,vert16x10,vert16x11,vert16x12,vert16x13)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'poles.f'
      include 'scale.f'
      include 'masses.f'
      include 'decl_kininv.f'
      
      real(dp):: mtsq
      integer:: ep
      complex(dp):: facuLl,facdLl
      complex(dp):: vert25x1,vert25x2,vert25x3,vert16x2,vert16x3,
     & vert16x4,vert16x5,vert16x6,vert16x7,vert16x8,vert16x9,
     & vert16x10,vert16x11,vert16x12,vert16x13
      complex(dp):: qlI2diffs346s25(-2:0),qlI2,qlI3,qlI2diff(-2:0)
      complex(dp):: qlI2diffs134s25(-2:0)
      complex(dp):: qlI2x25,qlI2x34,qlI2x134,qlI2x346
      real(dp):: p2Dp5
      
      qlI2x25=qlI2(s25,0d0,0d0,musq,ep)
      qlI2x34=qlI2(s34,0d0,0d0,musq,ep)
      qlI2x134=qlI2(s134,0d0,0d0,musq,ep)
      qlI2x346=qlI2(s346,0d0,0d0,musq,ep)
      
      qlI2diffs346s25(ep)=qlI2x346-qlI2x25
      qlI2diffs134s25(ep)=qlI2x134-qlI2x25
      qlI2diff(ep)=qlI2(s25,0d0,mtsq,musq,ep)-
     &qlI2(mtsq,0d0,mtsq,musq,ep)
      p2Dp5=0.5d0*(s25-mtsq)
      
      vert25x1=
     &  +4d0*qlI3(0d0,mtsq,s25,0d0,0d0,mtsq,musq,ep)*p2Dp5
     &  +qlI2(mtsq,0d0,mtsq,musq,ep)
     &  +qlI2diff(ep)*(4d0-0.5d0/p2Dp5*s25+1.5d0/p2Dp5*mtsq)+fp(ep)
      vert25x2=-mt/p2Dp5*qlI2diff(ep)
      vert25x3=-2d0*mt/s25*(fp(ep)+qlI2diff(ep))
      vert16x2=+ (                                                    
     &     - 2.D0/(s25 - s346)*fp(ep)				    
     &     - 4.D0/(s25 - s346)*qlI2x25		    
     &  - 4.D0/(s25 - s346)*qlI3(0d0,s346,s25,0d0,0d0,0d0,musq,ep)*s25  
   
     &     - 6.D0/(s25 - s346)*qlI2diffs346s25(ep)		    
     &     + 10.D0/(s25 - s346)/(s25 - s346)*qlI2diffs346s25(ep)*s25  
     &     )							    
      vert16x3= + (						    
     &     - 2.D0/(s25 - s346)*qlI2diffs346s25(ep)		    
     &     )							    
      vert16x4= + (  						    
     &     + fp(ep)						    
     &     + qlI2x25				    
     &     + 2.D0*qlI3(0d0,s346,s25,0d0,0d0,0d0,musq,ep)*s25		    
     &     + qlI2diffs346s25(ep)					    
     &     - 3.D0/(s25 - s346)*qlI2diffs346s25(ep)*s25                
     &     )

      vert16x5 =  + (
     &     + 4.D0*qlI3(0d0,s34,s346,0d0,0d0,0d0,musq,ep)*facdLl    
     &     - 2.D0/(s346 - s34)*fp(ep)*facdLl		     
     &     - 10.D0/(s346 - s34)*qlI2x34*facdLl    
     &     + 6.D0/(s346 - s34)*qlI2x346*facdLl    
     &     - 4.D0/(s346 - s34)*qlI3(0d0,s34,s346,0d0,0d0,0d0,musq,ep)*
     &    facdLl*s346						     
     &     + 10.D0/(s346 - s34)/(s346 - s34)*qlI2x34* 
     &    facdLl*s346					     
     &     - 10.D0/(s346 - s34)/(s346 - s34)*qlI2x346*
     &    facdLl*s346					     
     &     )							     
      vert16x6 = + (						     
     &     - 2.D0/(s346 - s34)*qlI2x34*facdLl
     &     + 2.D0/(s346 - s34)*qlI2x346*facdLl    
     &     )							     
      vert16x7= +  (						     
     & + fp(ep)*facdLl				     
     & + 3.D0*qlI2x34*facdLl	     
     & - 2.D0*qlI2x346*facdLl	     
     & + 2.D0*qlI3(0d0,s34,s346,0d0,0d0,0d0,musq,ep)*facdLl*s34
     & - 3.D0/(s346 - s34)*qlI2x34*facdLl
     &*s346							     
     & + 3.D0/(s346 - s34)*qlI2x346*facdLl    
     &*s346
     & )
      vert16x8=- (
     & + 2.D0/(s134 - s34)*fp(ep)*facuLl                 
     & + 4.D0/(s134 - s34)*qlI2x134*facuLl 
     & + 4.D0/(s134 - s34)*qlI3(s34,0d0,s134,0d0,0d0,0d0,musq,ep)
     &*facuLl*s34						      
     & - 10.D0/(s134 - s34)/(s134 - s34)*qlI2x34*facuLl*s34					      
     & + 10.D0/(s134 - s34)/(s134 - s34)*qlI2x134*facuLl*s34					      
     & )							      
      vert16x9 =-  (						      
     & + 2.D0/(s134 - s34)*qlI2x34*facuLl  
     & - 2.D0/(s134 - s34)*qlI2x134*facuLl 
     &     )							      
      vert16x10=-(						      
     & - fp(ep)*facuLl				      
     & - qlI2x134*facuLl		      
     & - 2.D0*qlI3(s34,0d0,s134,0d0,0d0,0d0,musq,ep)*facuLl*s34  
     & + 3.D0/(s134 - s34)*qlI2x34*facuLl* 
     &s34							      
     & - 3.D0/(s134 - s34)*qlI2x134*facuLl 
     &*s34)


      vert16x11= +(                                                   
     &     - 2.D0/(s25 - s134)*fp(ep)
     &     - 4.D0/(s25 - s134)*qlI2x25
     &     - 4.D0/(s25 - s134)*qlI3(s134,0d0,s25,0d0,0d0,0d0,musq,ep)
     & *s25
     &     - 6.D0/(s25 - s134)*qlI2diffs134s25(ep)
     &     + 10.D0/(s25 - s134)/(s25 - s134)*qlI2diffs134s25(ep)
     &     *s25)							  
									  
      vert16x12=+(							  
     &     - 2.D0/(s25 - s134)*qlI2diffs134s25(ep)		  
     &     )							  
      vert16x13=+(						  
     &     + fp(ep)						  
     &     + qlI2x25				  
     &     + 2.D0*qlI3(s134,0d0,s25,0d0,0d0,0d0,musq,ep)*s25	  
     &     + qlI2diffs134s25(ep)					  
     &     - 3.D0/(s25 - s134)*qlI2diffs134s25(ep)*s25              
     &     )

      return
      end



