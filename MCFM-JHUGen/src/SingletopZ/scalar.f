      subroutine scalar(p,yscalar)
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
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'nwz.f'
      integer:: j3,k1,k2,ep,eta
      complex(dp):: prW,yscalar(2,2,-2:0),
     & qlI2,qlI3,qlI2diff(-2:0),izb,
     & vert25x4,vert16x1,cprop,iprZ
      real(dp):: p(mxpart,4),q(mxpart,4),mtsq
      real(dp):: p16(4),p34(4),p25(4),
     & s16,s25,s34,p5Deta,p2Dp5,sinW,cosW,he
      real(dp):: Zm(-2:0)
      integer:: p1,p2,p3,p4,k5,e5,p6
      parameter(p1=1,p2=2,k5=5,p6=7,e5=6)

c----statement functions
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      izb(k1,k2)=cone/zb(k1,k2)
c----end statement functions

      sinW=sqrt(xw)
      cosW=sqrt(1d0-xw)
C---choose auxiliary vector
      eta=2
      mtsq=mt**2
      p16(:)=p(1,:)+p(6,:)
      p25(:)=p(2,:)+p(5,:)
      p34(:)=p(3,:)+p(4,:)
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
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

      Zm(-2)=0d0
      Zm(-1)=-3d0/2d0
      Zm(0)=(-3d0*log(musq/mtsq)-5d0)/2d0
      do ep=-2,0
      qlI2diff(ep)=qlI2(s25,0d0,mtsq,musq,ep)
      enddo
      do ep=-2,0
      qlI2diff(ep)=qlI2diff(ep)-qlI2(mtsq,0d0,mtsq,musq,ep)
      vert16x1=
     &  - 2d0*qlI3(s16,0d0,0d0,0d0,0d0,0d0,musq,ep)*s16
     &  - 3d0*qlI2(s16,0d0,0d0,musq,ep)-fp(ep)
      vert25x4=
     &  -4d0*qlI3(0d0,mtsq,s25,0d0,0d0,mtsq,musq,ep)*p2Dp5
     &  +2d0*qlI2(mtsq,0d0,mtsq,musq,ep)
     &  -qlI2diff(ep)*mt**2/p2Dp5
      do j3=1,2
      if (j3 == 1) then
        p3=3
        p4=4
        he=le
      elseif (j3 == 2) then
        p3=4
        p4=3
        he=re
      endif

      yscalar(j3,1,ep)= + prW(s16)*prW(s25)*alpha1*mtsq * (  - 2.D0*
     &    iprZ*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*Zm(ep)*
     &    s34**(-1)*qe - iprZ*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    vert25x4*s34**(-1)*qe + 2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*
     &    izb(k5,e5)*Zm(ep)*cosW**(-1)*he*sinW + za(p3,p6)*zb(p1,p4)*
     &    zb(p2,e5)*izb(k5,e5)*vert25x4*cosW**(-1)*he*sinW )
      yscalar(j3,1,ep) = yscalar(j3,1,ep) + prW(s16)*prW(s25)*vert16x1*
     & alpha1*mtsq * (  - iprZ*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)
     &    *s34**(-1)*qe + za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    cosW**(-1)*he*sinW )

      yscalar(j3,2,ep)= + prW(s16)*prW(s25)*alpha1*mt * ( 2.D0*iprZ*za(
     &    p3,p6)*zb(p1,p4)*zb(p2,k5)*Zm(ep)*s34**(-1)*qe + iprZ*za(p3,
     &    p6)*zb(p1,p4)*zb(p2,k5)*vert25x4*s34**(-1)*qe - 2.D0*za(p3,p6
     &    )*zb(p1,p4)*zb(p2,k5)*Zm(ep)*cosW**(-1)*he*sinW - za(p3,p6)*
     &    zb(p1,p4)*zb(p2,k5)*vert25x4*cosW**(-1)*he*sinW )
      yscalar(j3,2,ep) = yscalar(j3,2,ep) + prW(s16)*prW(s25)*vert16x1*
     & alpha1*mt * ( iprZ*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*s34**(-1)*qe - 
     &    za(p3,p6)*zb(p1,p4)*zb(p2,k5)*cosW**(-1)*he*sinW )


      yscalar(j3,1,ep)=yscalar(j3,1,ep)*cprop
      yscalar(j3,2,ep)=yscalar(j3,2,ep)*cprop

      enddo
      enddo
      return
      end
