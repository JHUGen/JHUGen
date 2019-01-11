      subroutine scalardk(q,yscalar)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'poles.f'
      include 'alpha1.f'
      include 'masses.f'
      include 'scale.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'nwz.f'
      include 'decl_kininv.f'
      integer j3,k1,k2,ep
      double complex prW,yscalar(2,-2:0),
     & qlI2,qlI3,qlI2diff(-2:0),izb,
     & vert25x4,vert16x1,cprop,iprZ
      double precision q(mxpart,4),mtsq
      double precision p2Dp5,sinW,cosW,he
      double precision Zm(-2:0)
      integer p1,p2,p3,p4,k5,e5,p6
      parameter(p1=1,p2=2,k5=5,p6=7,e5=6)

c----statement functions
      prW(s16)=cone/dcmplx(s16-wmass**2,zip)
      izb(k1,k2)=cone/zb(k1,k2)
c----end statement functions

      sinW=dsqrt(xw)
      cosW=dsqrt(1d0-xw)
      mtsq=mt**2
      p2Dp5=1d0/2d0*(s25-mtsq)
      
      if (nwz .eq. +1) then
      call spinoru(7,q,za,zb)
      elseif (nwz .eq. -1) then
      call spinoru(7,q,zb,za)
      endif
      

c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=dcmplx(1d0/dsqrt((s34-zmass**2)**2+(zmass*zwidth)**2))
      iprZ=dcmplx(s34-zmass**2)
      cprop=cprop/dcmplx(zip,mt*twidth)

      Zm(-2)=0d0
      Zm(-1)=-3d0/2d0
      Zm(0)=(-3d0*dlog(musq/mtsq)-5d0)/2d0

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
      if (j3 .eq. 1) then
        p3=3
        p4=4
        he=le
      elseif (j3 .eq. 2) then
        p3=4
        p4=3
        he=re
      endif

      yscalar(j3,ep)= + prW(s16)*prW(s25)*alpha1*mtsq * (  - 2.D0*iprZ*
     &    za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*Zm(ep)*s34**(-1)*qe
     &     - iprZ*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*vert25x4*
     &    s34**(-1)*qe + 2.D0*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    Zm(ep)*cosW**(-1)*he*sinW + za(p3,p6)*zb(p1,p4)*zb(p2,e5)*
     &    izb(k5,e5)*vert25x4*cosW**(-1)*he*sinW )
      yscalar(j3,ep) = yscalar(j3,ep) + prW(s16)*prW(s25)*vert16x1*
     & alpha1*mtsq * (  - iprZ*za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)
     &    *s34**(-1)*qe + za(p3,p6)*zb(p1,p4)*zb(p2,e5)*izb(k5,e5)*
     &    cosW**(-1)*he*sinW )


      yscalar(j3,ep)=yscalar(j3,ep)*cprop

      enddo
      enddo
      return
      end
