      subroutine extra(p,yextra)
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
      include 'nwz.f'
      integer:: k1,k2,ep,eta
      complex(dp):: prW,yextra(2,2,-2:0),
     & qlI2,qlI3,qlI2diff(-2:0),iza,izb,
     & vert25x1,vert25x2,vert25x3,vert16x1,iprZ,cprop
      real(dp):: p(mxpart,4),q(mxpart,4),mtsq,mwsq
      real(dp):: p16(4),p34(4),p235(4),p25(4),
     & s16,s25,s34,s235,p5Deta,p2Dp5,omal
      integer:: p1,p2,p3,p4,k5,e5,p6
      parameter(p1=1,p2=2,p3=3,p4=4,k5=5,p6=7,e5=6)

c----statement functions
      prW(s16)=cone/cplx1(s16-wmass**2)
      iza(k1,k2)=cone/za(k1,k2)
      izb(k1,k2)=cone/zb(k1,k2)
c----end statement functions

      omal=1d0-alpha1

c--- amplitude only exists for left-handed lepton couplings
      yextra(2,:,:)=czip
C---choose auxiliary vector
      eta=2
      mtsq=mt**2
      mwsq=wmass**2
      p235(:)=p(2,:)+p(3,:)+p(5,:)
      p16(:)=p(1,:)+p(6,:)
      p25(:)=p(2,:)+p(5,:)
      p34(:)=p(3,:)+p(4,:)
      s235=p235(4)**2-p235(1)**2-p235(2)**2-p235(3)**2
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

      yextra(1,1,ep)= + prW(s16)*prW(s25)*vert25x1*xw**(-1)*s235**(-1)
     &  * (  - za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,p4) + za(p3,k5)*za(
     &    p4,p6)*zb(p1,p4)*zb(p2,p4) )
      yextra(1,1,ep) = yextra(1,1,ep) + prW(s16)*prW(s25)*vert25x1*
     & xw**(-1) * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    omal*mwsq**(-1)*mtsq )
      yextra(1,1,ep) = yextra(1,1,ep) + prW(s16)*prW(s25)*vert25x2*mt*
     & xw**(-1)*s235**(-1) * ( 1.D0/2.D0*za(p1,p6)*za(p2,p3)*zb(p1,p2)*
     &    zb(p1,p4)*zb(p2,e5)*izb(k5,e5) - 1.D0/2.D0*za(p2,p3)*za(p4,p6
     &    )*zb(p1,p4)*zb(p2,p4)*zb(p2,e5)*izb(k5,e5) )
      yextra(1,1,ep) = yextra(1,1,ep) + prW(s16)*prW(s25)*vert25x2*mt*
     & xw**(-1) * ( 1.D0/4.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    omal*mwsq**(-1)*s25 - 1.D0/4.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)
     &    *izb(k5,e5)*omal*mwsq**(-1)*mtsq )
      yextra(1,1,ep) = yextra(1,1,ep) + prW(s16)*prW(s25)*vert25x3*mt*
     & xw**(-1) * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,
     &    e5) + 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*omal
     &    *mwsq**(-1)*s25 )
      yextra(1,1,ep) = yextra(1,1,ep) + prW(s16)*prW(s25)*vert16x1*
     & xw**(-1)*s235**(-1) * (  - za(p1,p6)*za(p3,k5)*zb(p1,p2)*zb(p1,
     &    p4) + za(p3,k5)*za(p4,p6)*zb(p1,p4)*zb(p2,p4) )
      yextra(1,1,ep) = yextra(1,1,ep) + prW(s16)*prW(s25)*vert16x1*
     & xw**(-1) * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    omal*mwsq**(-1)*mtsq )

      yextra(1,2,ep)= + prW(s16)*prW(s25)*vert25x1*mt*xw**(-1)*
     & s235**(-1) * ( za(p1,p6)*za(p3,e5)*zb(p1,p2)*zb(p1,p4)*iza(k5,e5
     &    ) - za(p3,e5)*za(p4,p6)*zb(p1,p4)*zb(p2,p4)*iza(k5,e5) )
      yextra(1,2,ep) = yextra(1,2,ep) + prW(s16)*prW(s25)*vert25x1*mt*
     & xw**(-1) * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*
     &    mwsq**(-1) )
      yextra(1,2,ep) = yextra(1,2,ep) + prW(s16)*prW(s25)*vert25x2*
     & xw**(-1)*s235**(-1) * (  - 1.D0/2.D0*za(p1,p6)*za(p2,p3)*zb(p1,
     &    p2)*zb(p1,p4)*zb(p2,k5) + 1.D0/2.D0*za(p2,p3)*za(p4,p6)*zb(p1
     &    ,p4)*zb(p2,p4)*zb(p2,k5) )
      yextra(1,2,ep) = yextra(1,2,ep) + prW(s16)*prW(s25)*vert25x2*
     & xw**(-1) * (  - 1.D0/4.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*
     &    mwsq**(-1)*s25 + 1.D0/4.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal
     &    *mwsq**(-1)*mtsq )
      yextra(1,2,ep) = yextra(1,2,ep) + prW(s16)*prW(s25)*vert25x3*
     & xw**(-1) * ( 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5) - 1.D0/2.D0
     &    *za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*mwsq**(-1)*s25 )
      yextra(1,2,ep) = yextra(1,2,ep) + prW(s16)*prW(s25)*vert16x1*mt*
     & xw**(-1)*s235**(-1) * ( za(p1,p6)*za(p3,e5)*zb(p1,p2)*zb(p1,p4)*
     &    iza(k5,e5) - za(p3,e5)*za(p4,p6)*zb(p1,p4)*zb(p2,p4)*iza(k5,
     &    e5) )
      yextra(1,2,ep) = yextra(1,2,ep) + prW(s16)*prW(s25)*vert16x1*mt*
     & xw**(-1) * (  - 1.D0/2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,k5)*omal*
     &    mwsq**(-1) )


      yextra(1,1,ep)=yextra(1,1,ep)*iprZ*cprop
      yextra(1,2,ep)=yextra(1,2,ep)*iprZ*cprop
      enddo

      return
      end
