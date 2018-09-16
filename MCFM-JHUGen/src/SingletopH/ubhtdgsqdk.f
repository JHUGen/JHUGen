      subroutine ubhtdgsqdk(p1,p2,p3,p4,k5,p6,p7,e5,
     &  mdecay,ampsq16,ampsq25)
      implicit none
      include 'types.f'

C     Matrix element squared for u(p1)+b(p2)->h(p3,p4)+t(p5)+d(p6)+g(p7)
C     split into two contributions:
C     radiation from 16 and 25 lines (ampsq16 and ampsq25)
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
      integer:: p1,p2,p3,p4,k5,p6,p7,e5,jg,icol,t16,t25
      parameter(t16=1,t25=2)
      real(dp):: sz,s16,s167,s1346,s13467,s12346,s126,s1267
      real(dp):: mw,ampsq,theta,ampsq16,ampsq25
      complex(dp):: prW,prt,iza,izb,Amp(2,2)
      complex(dp):: zab2,zba2,mdecay
c--- definitions of propagators
      theta(sz)=half+half*sign(one,sz)
      prW(sz)=cone/cplx2(sz-wmass**2,theta(sz)*wmass*wwidth)
      prt(sz)=cone/cplx2(sz-mt**2,zip)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
c--- definitions of compound spinors
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
C   End statement functions
      mw=wmass

      s16=s(p1,p6)
      s167=s(p1,p6)+s(p1,p7)+s(p6,p7)
      s12346=s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p1,p6)
     &               +s(p2,p3)+s(p2,p4)+s(p2,p6)
     &                        +s(p3,p4)+s(p3,p6)
     &                                 +s(p4,p6)
      s1346=s(p1,p3)+s(p1,p4)+s(p1,p6)
     &              +s(p3,p4)+s(p3,p6)
     &                       +s(p4,p6)
      s13467=s(p1,p3)+s(p1,p4)+s(p1,p6)+s(p1,p7)
     &               +s(p3,p4)+s(p3,p6)+s(p3,p7)
     &                        +s(p4,p6)+s(p4,p7)
     &                                 +s(p6,p7)
      s126=s(p1,p2)+s(p1,p6)+s(p2,p6)
      s1267=s(p1,p2)+s(p1,p6)+s(p1,p7)
     &              +s(p2,p6)+s(p2,p7)
     &                       +s(p6,p7)

      Amp(1,t16)= + prW(s167)*prW(s13467) * (  - 1.D0/2.D0*izb(p1,p7)*
     &    izb(p6,p7)*izb(k5,e5)*zb(p1,p6)*zb(p2,e5)*zab2(p6,p3,p4,p1)*
     &    mt**2*mw**(-1)*cWWH + izb(p1,p7)*izb(p6,p7)*za(p6,k5)*zb(p1,
     &    p2)*zb(p1,p6)*mw*cWWH - 1.D0/2.D0*izb(p6,p7)*izb(k5,e5)*zb(p2
     &    ,e5)*zab2(p7,p3,p4,p1)*mt**2*mw**(-1)*cWWH + izb(p6,p7)*za(p7
     &    ,k5)*zb(p1,p2)*mw*cWWH )
      Amp(1,t16) = Amp(1,t16) + prW(s167)*prt(s1267) * ( 1.D0/2.D0*izb(
     &    p1,p7)*izb(p6,p7)*izb(k5,e5)*zb(p1,p2)*zb(p1,p6)*zba2(e5,p3,
     &    p4,p6)*mt**2*mw**(-1)*cttH + izb(p1,p7)*izb(p6,p7)*za(p6,k5)*
     &    zb(p1,p2)*zb(p1,p6)*mt**2*mw**(-1)*cttH + 1.D0/2.D0*izb(p6,p7
     &    )*izb(k5,e5)*zb(p1,p2)*zba2(e5,p3,p4,p7)*mt**2*mw**(-1)*cttH
     &     + izb(p6,p7)*za(p7,k5)*zb(p1,p2)*mt**2*mw**(-1)*cttH )

      Amp(2,t16)= + prW(s167)*prW(s13467) * ( 1.D0/2.D0*iza(p1,p7)*iza(
     &    p6,p7)*izb(k5,e5)*za(p1,p6)*zb(p2,e5)*zab2(p6,p3,p4,p1)*mt**2
     &    *mw**(-1)*cWWH - iza(p1,p7)*iza(p6,p7)*za(p1,p6)*za(p6,k5)*
     &    zb(p1,p2)*mw*cWWH - 1.D0/2.D0*iza(p1,p7)*izb(k5,e5)*zb(p2,e5)
     &    *zab2(p6,p3,p4,p7)*mt**2*mw**(-1)*cWWH - iza(p1,p7)*za(p6,k5)
     &    *zb(p2,p7)*mw*cWWH )
      Amp(2,t16) = Amp(2,t16) + prW(s167)*prt(s1267) * (  - 1.D0/2.D0*
     &    iza(p1,p7)*iza(p6,p7)*izb(k5,e5)*za(p1,p6)*zb(p1,p2)*zba2(e5,
     &    p3,p4,p6)*mt**2*mw**(-1)*cttH - iza(p1,p7)*iza(p6,p7)*za(p1,
     &    p6)*za(p6,k5)*zb(p1,p2)*mt**2*mw**(-1)*cttH - 1.D0/2.D0*iza(
     &    p1,p7)*izb(k5,e5)*zb(p2,p7)*zba2(e5,p3,p4,p6)*mt**2*mw**(-1)*
     &    cttH - iza(p1,p7)*za(p6,k5)*zb(p2,p7)*mt**2*mw**(-1)*cttH )


      Amp(1,t25)= + prW(s16)*prW(s1346)*prt(s12346) * ( 1.D0/2.D0*izb(
     &    p2,p7)*izb(k5,e5)*za(p1,p7)*zb(p1,p2)*zb(p2,e5)*zab2(p6,p3,p4
     &    ,p1)*mt**2*mw**(-1)*cWWH - izb(p2,p7)*izb(k5,e5)*za(p6,p7)*
     &    zb(p1,p2)*zb(p2,e5)*mt**2*mw*cWWH - 1.D0/2.D0*izb(p2,p7)*izb(
     &    k5,e5)*za(p6,p7)*zb(p2,p6)*zb(p2,e5)*zab2(p6,p3,p4,p1)*mt**2*
     &    mw**(-1)*cWWH - 1.D0/2.D0*izb(p2,p7)*izb(k5,e5)*zb(p2,e5)*
     &    zab2(p6,p3,p4,p1)*zab2(p7,p3,p4,p2)*mt**2*mw**(-1)*cWWH -
     &    izb(p2,p7)*za(p1,p6)*za(p7,k5)*zb(p1,p2)**2*mw*cWWH + izb(p2,
     &    p7)*za(p7,k5)*zb(p1,p2)*zba2(p2,p3,p4,p6)*mw*cWWH )
      Amp(1,t25) = Amp(1,t25) + prW(s16)*prt(s1267)*prt(s126) * (  - 1.D
     &    0/2.D0*izb(p2,p7)*izb(k5,e5)*za(p1,p6)*zb(p1,p2)**2*zba2(e5,
     &    p3,p4,p7)*mt**2*mw**(-1)*cttH - izb(p2,p7)*izb(k5,e5)*za(p6,
     &    p7)*zb(p1,p2)*zb(p2,e5)*mt**4*mw**(-1)*cttH - izb(p2,p7)*za(
     &    p1,p6)*za(p7,k5)*zb(p1,p2)**2*mt**2*mw**(-1)*cttH - 1.D0/2.D0
     &    *izb(p2,p7)*za(p6,p7)*zb(p1,p2)*zab2(k5,p3,p4,p2)*mt**2*
     &    mw**(-1)*cttH )
      Amp(1,t25) = Amp(1,t25) + prW(s16)*prt(s12346)*prt(s126) * ( 1.D0/
     &    2.D0*izb(p2,p7)*izb(k5,e5)*za(p1,p6)*zb(p1,p2)*zb(p2,e5)*
     &    zab2(p7,p3,p4,p1)*mt**2*mw**(-1)*cttH + 1.D0/2.D0*izb(p2,p7)*
     &    izb(k5,e5)*za(p2,p6)*zb(p1,p2)*zb(p2,e5)*zab2(p7,p3,p4,p2)*
     &    mt**2*mw**(-1)*cttH - izb(p2,p7)*izb(k5,e5)*za(p6,p7)*zb(p1,
     &    p2)*zb(p2,e5)*mt**4*mw**(-1)*cttH - izb(p2,p7)*za(p1,p6)*za(
     &    p7,k5)*zb(p1,p2)**2*mt**2*mw**(-1)*cttH + 1.D0/2.D0*izb(p2,p7
     &    )*za(p7,k5)*zb(p1,p2)*zba2(p2,p3,p4,p6)*mt**2*mw**(-1)*cttH )
      Amp(1,t25) = Amp(1,t25) + prW(s16)*prt(s12346) * (  - 1.D0/2.D0*
     &    izb(p2,p7)*izb(k5,e5)*za(p6,p7)*zb(p1,p2)*zb(p2,e5)*mt**2*
     &    mw**(-1)*cttH )

      Amp(2,t25)= + prW(s16)*prW(s1346)*prt(s12346) * (  - 1.D0/2.D0*
     &    iza(p2,p7)*izb(k5,e5)*za(p1,p2)*zb(p1,p2)*zb(p7,e5)*zab2(p6,
     &    p3,p4,p1)*mt**2*mw**(-1)*cWWH - iza(p2,p7)*izb(k5,e5)*za(p2,
     &    p6)*zb(p1,p2)*zb(p7,e5)*mt**2*mw*cWWH - 1.D0/2.D0*iza(p2,p7)*
     &    izb(k5,e5)*za(p2,p6)*zb(p2,p6)*zb(p7,e5)*zab2(p6,p3,p4,p1)*
     &    mt**2*mw**(-1)*cWWH + 1.D0/2.D0*iza(p2,p7)*izb(k5,e5)*zb(p7,
     &    e5)*zab2(p2,p3,p4,p2)*zab2(p6,p3,p4,p1)*mt**2*mw**(-1)*cWWH
     &     + iza(p2,p7)*za(p1,p6)*za(p2,k5)*zb(p1,p2)*zb(p1,p7)*mw*cWWH
     &     + iza(p2,p7)*za(p2,p6)*za(p2,k5)*zb(p1,p2)*zb(p2,p7)*mw*cWWH
     &     - iza(p2,p7)*za(p2,k5)*zb(p1,p2)*zba2(p7,p3,p4,p6)*mw*cWWH
     &     - 1.D0/2.D0*iza(p2,p7)*za(p2,k5)*zb(p2,p7)*zab2(p6,p3,p4,p1)
     &    *mt**2*mw**(-1)*cWWH )
      Amp(2,t25) = Amp(2,t25) + prW(s16)*prW(s1346) * (  - 1.D0/2.D0*
     &    iza(p2,p7)*izb(k5,e5)*zb(p7,e5)*zab2(p6,p3,p4,p1)*mt**2*
     &    mw**(-1)*cWWH + iza(p2,p7)*za(p6,k5)*zb(p1,p7)*mw*cWWH )
      Amp(2,t25) = Amp(2,t25) + prW(s16)*prt(s1267)*prt(s126) * ( 1.D0/
     &    2.D0*iza(p2,p7)*izb(k5,e5)*za(p1,p6)*zb(p1,p2)*zb(p1,p7)*
     &    zba2(e5,p3,p4,p2)*mt**2*mw**(-1)*cttH + 1.D0/2.D0*iza(p2,p7)*
     &    izb(k5,e5)*za(p2,p6)*zb(p1,p2)*zb(p2,p7)*zba2(e5,p3,p4,p2)*
     &    mt**2*mw**(-1)*cttH - iza(p2,p7)*izb(k5,e5)*za(p2,p6)*zb(p1,
     &    p2)*zb(p7,e5)*mt**4*mw**(-1)*cttH + iza(p2,p7)*za(p1,p6)*za(
     &    p2,k5)*zb(p1,p2)*zb(p1,p7)*mt**2*mw**(-1)*cttH + iza(p2,p7)*
     &    za(p2,p6)*za(p2,k5)*zb(p1,p2)*zb(p2,p7)*mt**2*mw**(-1)*cttH
     &     - 1.D0/2.D0*iza(p2,p7)*za(p2,p6)*zb(p1,p2)*zab2(k5,p3,p4,p7)
     &    *mt**2*mw**(-1)*cttH )
      Amp(2,t25) = Amp(2,t25) + prW(s16)*prt(s1267) * ( 1.D0/2.D0*iza(
     &    p2,p7)*izb(k5,e5)*zb(p1,p7)*zba2(e5,p3,p4,p6)*mt**2*mw**(-1)*
     &    cttH + iza(p2,p7)*za(p6,k5)*zb(p1,p7)*mt**2*mw**(-1)*cttH )
      Amp(2,t25) = Amp(2,t25) + prW(s16)*prt(s12346)*prt(s126) * ( 1.D0/
     &    2.D0*iza(p2,p7)*izb(k5,e5)*za(p1,p2)*za(p2,p6)*zb(p1,p2)**2*
     &    zb(p7,e5)*mt**2*mw**(-1)*cttH + 1.D0/2.D0*iza(p2,p7)*izb(k5,
     &    e5)*za(p1,p6)*za(p2,p6)*zb(p1,p2)*zb(p1,p6)*zb(p7,e5)*mt**2*
     &    mw**(-1)*cttH - 1.D0/2.D0*iza(p2,p7)*izb(k5,e5)*za(p1,p6)*zb(
     &    p1,p2)*zb(p7,e5)*zab2(p2,p3,p4,p1)*mt**2*mw**(-1)*cttH + 1.D0/
     &    2.D0*iza(p2,p7)*izb(k5,e5)*za(p2,p6)**2*zb(p1,p2)*zb(p2,p6)*
     &    zb(p7,e5)*mt**2*mw**(-1)*cttH - 1.D0/2.D0*iza(p2,p7)*izb(k5,
     &    e5)*za(p2,p6)*zb(p1,p2)*zb(p7,e5)*mt**4*mw**(-1)*cttH - 1.D0/
     &    2.D0*iza(p2,p7)*izb(k5,e5)*za(p2,p6)*zb(p1,p2)*zb(p7,e5)*
     &    zab2(p2,p3,p4,p2)*mt**2*mw**(-1)*cttH + iza(p2,p7)*za(p1,p6)*
     &    za(p2,k5)*zb(p1,p2)*zb(p1,p7)*mt**2*mw**(-1)*cttH + iza(p2,p7
     &    )*za(p2,p6)*za(p2,k5)*zb(p1,p2)*zb(p2,p7)*mt**2*mw**(-1)*cttH
     &     - 1.D0/2.D0*iza(p2,p7)*za(p2,k5)*zb(p1,p2)*zba2(p7,p3,p4,p6)
     &    *mt**2*mw**(-1)*cttH )


      ampsq16=0d0
      ampsq25=0d0
      do jg=1,2
      ampsq16=ampsq16+abs(Amp(jg,t16)*mdecay)**2
      ampsq25=ampsq25+abs(Amp(jg,t25)*mdecay)**2
      enddo

      return
      end
