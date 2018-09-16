      subroutine lowerdk_parttri(q,ylower,first)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'poles.f'
      include 'masses.f'
      include 'scale.f'
      include 'zcouple.f'
      include 'nwz.f'
      include 'decl_kininv.f'
      integer:: k1,k2,ep,epmin
      complex(dp):: prW,prt,ylower(2,-2:0),
     & qlI2,qlI3,izb,
     & vert16x1,vert25x5,vert25x6,
     & vert25x7,vert25x8,vert25x9,vert25x10,vert25x11,vert25x12,
     & vert25x13,vert25x14,vert25x15,vert25x16,vert25x17,vert25x18,
     & vert25x19,vert25x20,vert25x21,vert25x22,vert25x23,
     & vert25x24,vert25x25,vert25x26,vert25x27,vert25x28,
     & vert25x29,vert25x30,vert25x31,vert25x32,vert25x16a
      complex(dp):: facuLl,facuRl,facdLl,cprop,
     & iprZ
      real(dp):: q(mxpart,4),mtsq
      logical:: first
      integer:: j3,p1,p2,p3,p4,k5,e5,p6
      parameter(p1=1,p2=2,k5=5,p6=7,e5=6)

c----statement function
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      prt(s345)=cone/cplx2(s345-mt**2,zip)
      izb(k1,k2)=cone/zb(k1,k2)
c----end statement function

      mtsq=mt**2
      
      if (nwz == +1) then
      call spinoru(7,q,za,zb)
      elseif (nwz == -1) then
      call spinoru(7,q,zb,za)
      endif
      
c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=cplx1(1d0/sqrt((s34-zmass**2)**2+(zmass*zwidth)**2))
      cprop=cprop/cplx2(zip,mt*twidth)
      iprZ=cplx1(s34-zmass**2)
c--- only compute poles for checking on first call
      if (first) then
         epmin=-2
      else
         epmin=-1
      endif

c      write(*,*) 'epmin in lowerdk_tri', epmin
      do ep=epmin,0
 
      vert16x1=
     &  + 2d0*qlI3(s16,0d0,0d0,0d0,0d0,0d0,musq,ep)*s16
     &  + 3d0*qlI2(s16,0d0,0d0,musq,ep)+fp(ep)

      do j3=1,2
      if (j3 == 1) then
        p3=3
        p4=4
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*le)
        facuRl=cplx1(Qu*q1)*iprZ/s34+cplx1(R(2)*le)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*le)
      elseif (j3 == 2) then
        p3=4
        p4=3
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*re)
        facuRl=cplx1(Qu*q1)*iprZ/s34+cplx1(R(2)*re)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*re)
      endif
      call vertices_bt1(ep,facdLl,vert25x5,vert25x6
     & ,vert25x7)
      call vertices_bt2(mtsq,ep,
     & vert25x8,vert25x9,vert25x10,vert25x11,vert25x12,
     & vert25x13,vert25x14,vert25x15,vert25x16)
      call vertices_tt1(mtsq,ep,
     & vert25x16a,vert25x17,vert25x18,vert25x19,vert25x20)
      call vertices_tt2(mtsq,ep,facuLl,facuRl,
     & vert25x21,vert25x22,vert25x23,vert25x24,vert25x25,
     & vert25x26,vert25x27,vert25x28,vert25x29,vert25x30,
     & vert25x31,vert25x32)
      ylower(j3,ep)= + prW(s16)*prt(s345)*facuRl*mt * ( za(p1,p3)*za(p2
     &    ,p6)*zb(p1,p2)**2*zb(p4,e5)*izb(k5,e5)*vert25x16a + za(p2,p6)
     &    *za(p3,p6)*zb(p1,p2)*zb(p2,p6)*zb(p4,e5)*izb(k5,e5)*
     &    vert25x16a )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*facuLl*mt * ( 
     &    za(p2,p6)*za(p3,k5)*zb(p1,p2)*zb(p2,p4)*vert25x16a )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert16x1*
     & facuRl*mt**2 * ( 2.D0*za(p3,p6)*zb(p1,p2)*zb(p4,e5)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert16x1*
     & facuLl * (  - 2.D0*za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,p4) - 2.D0
     &    *za(p2,p6)*za(p3,k5)*zb(p1,p2)*zb(p2,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x17*
     & facuRl*mt * (  - za(p1,p3)*za(p6,p2)*zb(p1,p2)*zb(p2,p1)*zb(p4,
     &    e5)*izb(k5,e5) - za(p3,p6)*za(p6,p2)*zb(p2,p1)*zb(p2,p6)*zb(
     &    p4,e5)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x17*
     & facuLl*mt * (  - za(p3,k5)*za(p6,p2)*zb(p2,p1)*zb(p2,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x18*
     & facuRl*mt**2 * ( za(p1,p3)*za(p2,p6)*zb(p1,p2)**2*zb(p4,e5)*izb(
     &    k5,e5) + za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,p6)*zb(p4,e5)*
     &    izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x18*
     & facuLl * (  - za(p1,p2)*za(p2,p6)*za(p3,k5)*zb(p1,p2)**2*zb(p2,
     &    p4) + za(p1,p6)*za(p2,p6)*za(p3,k5)*zb(p1,p2)**2*zb(p4,p6) - 
     &    za(p1,p6)*za(p2,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,p4)*zb(p2,p6)
     &     - za(p2,p6)**2*za(p3,k5)*zb(p1,p2)*zb(p2,p4)*zb(p2,p6) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x19*
     & facuRl*mt**2 * (  - za(p1,p3)*za(p6,p2)*zb(p1,p2)*zb(p2,p1)*zb(
     &    p4,e5)*izb(k5,e5) - za(p3,p6)*za(p6,p2)*zb(p2,p1)*zb(p2,p6)*
     &    zb(p4,e5)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x19*
     & facuLl * ( za(p1,p2)*za(p3,k5)*za(p6,p2)*zb(p1,p2)*zb(p2,p1)*zb(
     &    p2,p4) - za(p1,p6)*za(p3,k5)*za(p6,p2)*zb(p1,p2)*zb(p2,p1)*
     &    zb(p4,p6) + za(p1,p6)*za(p3,k5)*za(p6,p2)*zb(p1,p4)*zb(p2,p1)
     &    *zb(p2,p6) + za(p2,p6)*za(p3,k5)*za(p6,p2)*zb(p2,p1)*zb(p2,p4
     &    )*zb(p2,p6) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x20*
     & facuRl*mt**2 * ( 2.D0*za(p3,p6)*zb(p1,p2)*zb(p4,e5)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x20*
     & facuLl * (  - 2.D0*za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,p4) - 2.D0
     &    *za(p2,p6)*za(p3,k5)*zb(p1,p2)*zb(p2,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x21*mt
     &  * ( za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,e5)*zb(k5,p4)*izb(k5,e5
     &    ) + za(p1,p6)*za(p3,e5)*zb(p1,p2)*zb(p1,e5)*zb(e5,p4)*izb(k5,
     &    e5) + za(p2,p6)*za(p3,k5)*zb(p1,p2)*zb(p2,e5)*zb(k5,p4)*izb(
     &    k5,e5) + za(p2,p6)*za(p3,e5)*zb(p1,p2)*zb(p2,e5)*zb(e5,p4)*
     &    izb(k5,e5) + za(p3,k5)*za(p6,k5)*zb(p1,p2)*zb(k5,p4) + za(p3,
     &    e5)*za(p6,k5)*zb(p1,p2)*zb(e5,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x24*mt
     &  * ( 2.D0*za(p1,p6)*za(p3,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,e5)*izb(
     &    k5,e5) + 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p2)*zb(p2,p4)*zb(p4,
     &    e5)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x23*mt
     &  * ( 2.D0*za(p1,p6)*za(p3,p4)*zb(p1,p2)*zb(p1,p4)*zb(p4,e5)*izb(
     &    k5,e5) + 2.D0*za(p2,p6)*za(p3,p4)*zb(p1,p2)*zb(p2,p4)*zb(p4,
     &    e5)*izb(k5,e5) - 2.D0*za(p3,p6)*za(k5,p3)*zb(p1,p2)*zb(p3,p4)
     &     )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x25*
     & mt**2 * ( za(p3,k5)*za(p6,p3)*zb(p1,p2)*zb(p3,e5)*zb(k5,p4)*izb(
     &    k5,e5) + za(p3,k5)*za(p6,p4)*zb(p1,p2)*zb(p4,e5)*zb(k5,p4)*
     &    izb(k5,e5) + za(p3,e5)*za(p6,p3)*zb(p1,p2)*zb(p3,e5)*zb(e5,p4
     &    )*izb(k5,e5) + za(p3,e5)*za(p6,p4)*zb(p1,p2)*zb(p4,e5)*zb(e5,
     &    p4)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x25 * ( 
     &    za(p1,p6)*za(p3,k5)*za(k5,p3)*zb(p1,p2)*zb(p3,p1)*zb(k5,p4)
     &     + za(p1,p6)*za(p3,k5)*za(k5,p4)*zb(p1,p2)*zb(p4,p1)*zb(k5,p4
     &    ) + za(p1,p6)*za(p3,e5)*za(k5,p3)*zb(p1,p2)*zb(p3,p1)*zb(e5,
     &    p4) + za(p1,p6)*za(p3,e5)*za(k5,p4)*zb(p1,p2)*zb(p4,p1)*zb(e5
     &    ,p4) + za(p2,p6)*za(p3,k5)*za(k5,p3)*zb(p1,p2)*zb(p3,p2)*zb(
     &    k5,p4) + za(p2,p6)*za(p3,k5)*za(k5,p4)*zb(p1,p2)*zb(p4,p2)*
     &    zb(k5,p4) + za(p2,p6)*za(p3,e5)*za(k5,p3)*zb(p1,p2)*zb(p3,p2)
     &    *zb(e5,p4) + za(p2,p6)*za(p3,e5)*za(k5,p4)*zb(p1,p2)*zb(p4,p2
     &    )*zb(e5,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x27 * ( 
     &    za(p1,p6)*za(p3,k5)*za(k5,p3)*zb(p1,p2)*zb(p3,p1)*zb(k5,p4)
     &     + za(p1,p6)*za(p3,k5)*za(k5,p4)*zb(p1,p2)*zb(p4,p1)*zb(k5,p4
     &    ) + za(p1,p6)*za(p3,e5)*za(k5,p3)*zb(p1,p2)*zb(p3,p1)*zb(e5,
     &    p4) + za(p1,p6)*za(p3,e5)*za(k5,p4)*zb(p1,p2)*zb(p4,p1)*zb(e5
     &    ,p4) + za(p2,p6)*za(p3,k5)*za(k5,p3)*zb(p1,p2)*zb(p3,p2)*zb(
     &    k5,p4) + za(p2,p6)*za(p3,k5)*za(k5,p4)*zb(p1,p2)*zb(p4,p2)*
     &    zb(k5,p4) + za(p2,p6)*za(p3,e5)*za(k5,p3)*zb(p1,p2)*zb(p3,p2)
     &    *zb(e5,p4) + za(p2,p6)*za(p3,e5)*za(k5,p4)*zb(p1,p2)*zb(p4,p2
     &    )*zb(e5,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x29*
     & mt**2 * (  - 2.D0*za(p3,p6)*zb(p1,p2)*zb(p4,e5)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x29 * ( 
     &    2.D0*za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,p4) + 2.D0*za(p2,p6)
     &    *za(p3,k5)*zb(p1,p2)*zb(p2,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x30 * ( 
     &    2.D0*za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,p4) + 2.D0*za(p2,p6)
     &    *za(p3,k5)*zb(p1,p2)*zb(p2,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*prt(s345)*vert25x31*mt
     &  * ( za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,e5)*zb(k5,p4)*izb(k5,e5
     &    ) + za(p1,p6)*za(p3,e5)*zb(p1,p2)*zb(p1,e5)*zb(e5,p4)*izb(k5,
     &    e5) + za(p2,p6)*za(p3,k5)*zb(p1,p2)*zb(p2,e5)*zb(k5,p4)*izb(
     &    k5,e5) + za(p2,p6)*za(p3,e5)*zb(p1,p2)*zb(p2,e5)*zb(e5,p4)*
     &    izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert16x1*facdLl*
     & s234**(-1) * (  - 2.D0*za(p3,p2)*za(p6,k5)*zb(p2,p1)*zb(p2,p4)
     &     - 2.D0*za(p3,p4)*za(p6,k5)*zb(p2,p4)*zb(p4,p1) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert25x5*s234**(-1) * ( 
     &     - za(p2,p3)**2*za(p6,k5)*zb(p1,p2)*zb(p2,p3)*zb(p2,p4) - za(
     &    p2,p3)*za(p2,p4)*za(p6,k5)*zb(p1,p2)*zb(p2,p4)**2 - za(p2,p3)
     &    *za(p3,p4)*za(p6,k5)*zb(p1,p3)*zb(p2,p4)**2 + za(p2,p3)*za(p3
     &    ,p4)*za(p6,k5)*zb(p1,p4)*zb(p2,p3)*zb(p2,p4) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert25x7*s234**(-1) * ( 
     &    2.D0*za(p3,p2)*za(p6,k5)*zb(p2,p1)*zb(p2,p4) + 2.D0*za(p3,p4)
     &    *za(p6,k5)*zb(p2,p4)*zb(p4,p1) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert25x9*facdLl*
     & s234**(-1) * (  - za(p1,k5)*za(p2,p3)*za(p6,k5)*zb(p1,p2)*zb(p2,
     &    p4)*zb(k5,p1) - za(p1,k5)*za(p2,p3)*za(p6,e5)*zb(p1,p2)*zb(p2
     &    ,p4)*zb(e5,p1) + za(p1,k5)*za(p3,p4)*za(p6,k5)*zb(p1,p4)*zb(
     &    p2,p4)*zb(k5,p1) + za(p1,k5)*za(p3,p4)*za(p6,e5)*zb(p1,p4)*
     &    zb(p2,p4)*zb(e5,p1) + za(p2,p3)*za(p6,k5)**2*zb(p2,p4)*zb(p2,
     &    p6)*zb(k5,p1) + za(p2,p3)*za(p6,k5)*za(p6,e5)*zb(p2,p4)*zb(p2
     &    ,p6)*zb(e5,p1) - za(p3,p4)*za(p6,k5)**2*zb(p2,p4)*zb(p4,p6)*
     &    zb(k5,p1) - za(p3,p4)*za(p6,k5)*za(p6,e5)*zb(p2,p4)*zb(p4,p6)
     &    *zb(e5,p1) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert25x10*facdLl*
     & s234**(-1) * (  - 2.D0*za(p3,p2)*za(p6,k5)*zb(p2,p1)*zb(p2,p4)
     &     - 2.D0*za(p3,p4)*za(p6,k5)*zb(p2,p4)*zb(p4,p1) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert25x12*facdLl*mt*
     & s234**(-1) * (  - za(p3,p2)*za(p6,k5)*zb(p2,p4)*zb(p2,e5)*zb(k5,
     &    p1)*izb(k5,e5) - za(p3,p2)*za(p6,e5)*zb(p2,p4)*zb(p2,e5)*zb(
     &    e5,p1)*izb(k5,e5) - za(p3,p4)*za(p6,k5)*zb(p2,p4)*zb(p4,e5)*
     &    zb(k5,p1)*izb(k5,e5) - za(p3,p4)*za(p6,e5)*zb(p2,p4)*zb(p4,e5
     &    )*zb(e5,p1)*izb(k5,e5) )
      ylower(j3,ep) = ylower(j3,ep) + prW(s16)*vert25x13*facdLl*mt*
     & s234**(-1) * (  - 2.D0*za(p3,p2)*za(p6,p1)*zb(p1,e5)*zb(p2,p1)*
     &    zb(p2,p4)*izb(k5,e5) - 2.D0*za(p3,p4)*za(p6,p1)*zb(p1,e5)*zb(
     &    p2,p4)*zb(p4,p1)*izb(k5,e5) )

      enddo
      enddo
      ylower=ylower*cprop
      return
      end
      
