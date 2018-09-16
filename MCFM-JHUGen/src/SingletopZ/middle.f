      subroutine middle(p,ymiddle)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'poles.f'
      include 'alpha1.f'
      include 'masses.f'
      include 'scale.f'
      include 'zcouple.f'
      include 'nwz.f'
      integer:: j3,k1,k2,ep,eta
      complex(dp):: prW,ymiddle(2,2,-2:0),
     & qlI2,qlI3,qlI2diff(-2:0),iza,izb,
     & vert25x1,vert25x2,vert25x3,vert16x1,cprop,iprZ
      complex(dp):: facuLl,facdLl,facLdiff
      real(dp):: p(mxpart,4),q(mxpart,4),mtsq,mwsq
      real(dp):: p12(4),p26(4),p16(4),p34(4),p25(4),
     & s16,s25,s34,s12,s26,p5Deta,p2Dp5,he,omal
      integer:: p1,p2,p3,p4,k5,e5,p6
      parameter(p1=1,p2=2,k5=5,p6=7,e5=6)

c----statement functions
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      iza(k1,k2)=cone/za(k1,k2)
      izb(k1,k2)=cone/zb(k1,k2)
c----end statement functions

      omal=1d0-alpha1

C---choose auxiliary vector
      eta=2
      mtsq=mt**2
      mwsq=wmass**2
      p12(:)=p(1,:)+p(2,:)
      p16(:)=p(1,:)+p(6,:)
      p26(:)=p(2,:)+p(6,:)
      p25(:)=p(2,:)+p(5,:)
      p34(:)=p(3,:)+p(4,:)
      s26=p26(4)**2-p26(1)**2-p26(2)**2-p26(3)**2
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
      s25=p25(4)**2-p25(1)**2-p25(2)**2-p25(3)**2
      p2Dp5=0.5d0*(s25-mtsq)
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

      do ep=-2,0
      qlI2diff(ep)=qlI2(s25,0d0,mtsq,musq,ep)
      enddo
      do ep=-2,0
      qlI2diff(ep)=qlI2diff(ep)-qlI2(mtsq,0d0,mtsq,musq,ep)
      vert25x1=
     &  -4d0*qlI3(0d0,mtsq,s25,0d0,0d0,mtsq,musq,ep)*p2Dp5
     &  -qlI2(mtsq,0d0,mtsq,musq,ep)
     &  -qlI2diff(ep)*(4d0-0.5d0/p2Dp5*s25+1.5d0/p2Dp5*mtsq)-fp(ep)
      vert25x2=mt/p2Dp5*qlI2diff(ep)
      vert25x3=2d0*mt/s25*(fp(ep)+qlI2diff(ep))
      vert16x1=
     &  - 2d0*qlI3(s16,0d0,0d0,0d0,0d0,0d0,musq,ep)*s16
     &  - 3d0*qlI2(s16,0d0,0d0,musq,ep)-fp(ep)
      do j3=1,2
      if (j3 == 1) then
        p3=3
        p4=4
        facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*le)*s34
        facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*le)*s34
        he=le
      elseif (j3 == 2) then
        p3=4
        p4=3
        facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*re)*s34
        facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*re)*s34
        he=re
      endif
      facLdiff=facuLl-facdLl

      ymiddle(j3,1,ep)= + vert25x1*s34**(-1)*facLdiff * ( 2.D0*za(p1,p3
     &    )*za(p6,k5)*zb(p1,p2)*zb(p1,p4) + 2.D0*za(p1,k5)*za(p3,p6)*
     &    zb(p1,p2)*zb(p1,p4) - 2.D0*za(p3,p6)*za(p3,k5)*zb(p1,p3)*zb(
     &    p2,p4) + 2.D0*za(p3,p6)*za(p6,k5)*zb(p1,p2)*zb(p4,p6) - 2.D0*
     &    za(p3,p6)*za(p6,k5)*zb(p1,p4)*zb(p2,p6) - za(p3,p6)*zb(p1,p4)
     &    *zb(p2,e5)*izb(k5,e5)*mtsq - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*
     &    izb(k5,e5)*omal*mwsq**(-1)*mtsq*s16 - 2.D0*za(p3,k5)*za(p4,p6
     &    )*zb(p1,p4)*zb(p2,p4) )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert25x1*facLdiff * ( za(p3
     &    ,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal*mwsq**(-1)*mtsq )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert25x2*mt*s34**(-1)*
     & facLdiff * (  - za(p1,p3)*za(p2,p6)*zb(p1,p2)*zb(p1,p4)*zb(p2,e5
     &    )*izb(k5,e5) + za(p2,p3)*za(p3,p6)*zb(p1,p3)*zb(p2,p4)*zb(p2,
     &    e5)*izb(k5,e5) + za(p2,p3)*za(p4,p6)*zb(p1,p4)*zb(p2,p4)*zb(
     &    p2,e5)*izb(k5,e5) - za(p2,p6)*za(p3,p6)*zb(p1,p2)*zb(p2,e5)*
     &    zb(p4,p6)*izb(k5,e5) - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,
     &    e5)*s26 - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    s25 + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*mtsq
     &     - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*s12 - 1.D0/2.D0*
     &    za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal*mwsq**(-1)*s16*
     &    s25 + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal
     &    *mwsq**(-1)*mtsq*s16 )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert25x2*mt*facLdiff * ( 1.D
     &    0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal*
     &    mwsq**(-1)*s25 - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(
     &    k5,e5)*omal*mwsq**(-1)*mtsq )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert25x3*mt*s34**(-1)*
     & facLdiff * ( za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*s16 - za(
     &    p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal*mwsq**(-1)*s16*s25
     &     )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert25x3*mt*facLdiff * (
     &     - za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5) + za(p3,p6)*zb(p1
     &    ,p4)*zb(p2,e5)*izb(k5,e5)*omal*mwsq**(-1)*s25 )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert16x1*s34**(-1)*facLdiff
     &  * ( 2.D0*za(p1,p3)*za(p6,k5)*zb(p1,p2)*zb(p1,p4) + 2.D0*za(p1,
     &    k5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4) - 2.D0*za(p3,p6)*za(p3,k5)*
     &    zb(p1,p3)*zb(p2,p4) + 2.D0*za(p3,p6)*za(p6,k5)*zb(p1,p2)*zb(
     &    p4,p6) - 2.D0*za(p3,p6)*za(p6,k5)*zb(p1,p4)*zb(p2,p6) - za(p3
     &    ,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*mtsq - za(p3,p6)*zb(p1,p4
     &    )*zb(p2,e5)*izb(k5,e5)*omal*mwsq**(-1)*mtsq*s16 - 2.D0*za(p3,
     &    k5)*za(p4,p6)*zb(p1,p4)*zb(p2,p4) )
      ymiddle(j3,1,ep) = ymiddle(j3,1,ep) + vert16x1*facLdiff * ( za(p3
     &    ,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal*mwsq**(-1)*mtsq )

      ymiddle(j3,2,ep)= + vert25x1*mt*s34**(-1)*facLdiff * (  - 2.D0*
     &    za(p1,p3)*za(p6,e5)*zb(p1,p2)*zb(p1,p4)*iza(k5,e5) - 2.D0*za(
     &    p1,e5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*iza(k5,e5) + 2.D0*za(p3,
     &    p6)*za(p3,e5)*zb(p1,p3)*zb(p2,p4)*iza(k5,e5) - 2.D0*za(p3,p6)
     &    *za(p6,e5)*zb(p1,p2)*zb(p4,p6)*iza(k5,e5) + 2.D0*za(p3,p6)*
     &    za(p6,e5)*zb(p1,p4)*zb(p2,p6)*iza(k5,e5) + za(p3,p6)*zb(p1,p4
     &    )*zb(p2,k5) + za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*mwsq**(-1)*
     &    s16 + 2.D0*za(p3,e5)*za(p4,p6)*zb(p1,p4)*zb(p2,p4)*iza(k5,e5)
     &     )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert25x1*mt*facLdiff * (
     &     - za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*mwsq**(-1) )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert25x2*s34**(-1)*facLdiff
     &  * ( za(p1,p3)*za(p2,p6)*zb(p1,p2)*zb(p1,p4)*zb(p2,k5) - za(p2,
     &    p3)*za(p3,p6)*zb(p1,p3)*zb(p2,p4)*zb(p2,k5) - za(p2,p3)*za(p4
     &    ,p6)*zb(p1,p4)*zb(p2,p4)*zb(p2,k5) + za(p2,p6)*za(p3,p6)*zb(
     &    p1,p2)*zb(p2,k5)*zb(p4,p6) + za(p3,p6)*zb(p1,p4)*zb(p2,k5)*
     &    s26 + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*s25 - 1.D0/2.D0
     &    *za(p3,p6)*zb(p1,p4)*zb(p2,k5)*mtsq + za(p3,p6)*zb(p1,p4)*zb(
     &    p2,k5)*s12 + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*
     &    mwsq**(-1)*s16*s25 - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*
     &    omal*mwsq**(-1)*mtsq*s16 )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert25x2*facLdiff * (  - 1.D
     &    0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*mwsq**(-1)*s25 + 1.D
     &    0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*mwsq**(-1)*mtsq )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert25x3*s34**(-1)*facLdiff
     &  * (  - za(p3,p6)*zb(p1,p4)*zb(p2,k5)*s16 + za(p3,p6)*zb(p1,p4)*
     &    zb(p2,k5)*omal*mwsq**(-1)*s16*s25 )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert25x3*facLdiff * ( za(p3
     &    ,p6)*zb(p1,p4)*zb(p2,k5) - za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal
     &    *mwsq**(-1)*s25 )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert16x1*mt*s34**(-1)*
     & facLdiff * (  - 2.D0*za(p1,p3)*za(p6,e5)*zb(p1,p2)*zb(p1,p4)*
     &    iza(k5,e5) - 2.D0*za(p1,e5)*za(p3,p6)*zb(p1,p2)*zb(p1,p4)*
     &    iza(k5,e5) + 2.D0*za(p3,p6)*za(p3,e5)*zb(p1,p3)*zb(p2,p4)*
     &    iza(k5,e5) - 2.D0*za(p3,p6)*za(p6,e5)*zb(p1,p2)*zb(p4,p6)*
     &    iza(k5,e5) + 2.D0*za(p3,p6)*za(p6,e5)*zb(p1,p4)*zb(p2,p6)*
     &    iza(k5,e5) + za(p3,p6)*zb(p1,p4)*zb(p2,k5) + za(p3,p6)*zb(p1,
     &    p4)*zb(p2,k5)*omal*mwsq**(-1)*s16 + 2.D0*za(p3,e5)*za(p4,p6)*
     &    zb(p1,p4)*zb(p2,p4)*iza(k5,e5) )
      ymiddle(j3,2,ep) = ymiddle(j3,2,ep) + vert16x1*mt*facLdiff * (
     &     - za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*mwsq**(-1) )


      ymiddle(j3,1,ep)=ymiddle(j3,1,ep)*cprop*prW(s16)*prW(s25)
      ymiddle(j3,2,ep)=ymiddle(j3,2,ep)*cprop*prW(s16)*prW(s25)

      enddo
      enddo
      return
      end
