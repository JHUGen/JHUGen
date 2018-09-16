      function ubhtdsq(p1,p2,p3,p4,k5,p6,e5)
      implicit none
      include 'types.f'
      real(dp):: ubhtdsq

C     Matrix element squared for u(p1)+b(p2)->h(p3,p4)+t(p5)+d(p6)
C     Matrix element is constructed so that dependence on p5
C     only enters through k5 and e5 in wave function
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'anomHiggs.f'
      real(dp):: s25,s16,s345,mw
      integer:: p1,p2,p3,p4,k5,p6,e5
      complex(dp):: prW,prt,iza,izb,Ampw(2),Ampf(2),Ampt(2),zab2
      zab2(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      prt(s16)=cone/cplx2(s16-mt**2,zip*mt*twidth)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
C
      mw=wmass
C-----since 1,2,3,4,6 are unchanged in demassification we may calculate
C-----if we only use these vectors.
      s16=s(p1,p6)
      s25=s(p1,p3)+s(p1,p4)+s(p1,p6)+s(p3,p4)+s(p3,p6)+s(p4,p6)
      s345=s(p1,p2)+s(p1,p6)+s(p2,p6)
c      do j=0,1
c      al=real(j,dp)
c      omal=1d0-al
      Ampw(2)= + prW(s16)*prW(s25)*iza(k5,e5)*mt*mw * ( za(p6,e5)*zb(p1
     &    ,p2) )

      Ampf(2)= + prW(s16)*prW(s25)*mt*mw**(-1) * (  - 1.D0/2.D0*zb(p2,
     &    k5)*zab2(p6,p3,p4,p1) )

      Ampt(2)= + prW(s16)*prt(s345)*mt*mw**(-1) * (  - 1.D0/2.D0*zb(p1,
     &    p2)*zab2(p6,p1,p2,k5) )
      Ampt(2) = Ampt(2) + prW(s16)*prt(s345)*iza(k5,e5)*mt**3*mw**(-1)
     &  * ( 1.D0/2.D0*za(p6,e5)*zb(p1,p2) )

      Ampw(1)= + prW(s16)*prW(s25)*mw * (  - za(p6,k5)*zb(p1,p2) )

      Ampf(1)= + prW(s16)*prW(s25)*izb(k5,e5)*mt**2*mw**(-1) * ( 1.D0/2.
     &    D0*zb(p2,e5)*zab2(p6,p3,p4,p1) )

      Ampt(1)= + prW(s16)*prt(s345)*mt**2*mw**(-1) * (  - 1.D0/2.D0*za(
     &    p6,k5)*zb(p1,p2) )
      Ampt(1) = Ampt(1) + prW(s16)*prt(s345)*izb(k5,e5)*mt**2*mw**(-1)
     &  * ( 1.D0/2.D0*zb(p1,p2)*zab2(p6,p1,p2,e5) )

c --   include anomalous Yukawa and EW couplings for the Higgs
      Ampt=cttH*Ampt
      Ampf=cWWH*Ampf
      Ampw=cWWH*Ampw

      ubhtdsq=
     & +abs(Ampw(1)+Ampf(1)+Ampt(1))**2
     & +abs(Ampw(2)+Ampf(2)+Ampt(2))**2
c      write(6,*) 'omal,ubht',omal,ubhtdsq
c      enddo
c      pause
      return
      end

c --   include anomalous Yukawa and EW couplings for the Higgs

