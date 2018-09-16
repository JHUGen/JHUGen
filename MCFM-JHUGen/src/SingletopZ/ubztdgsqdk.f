      subroutine ubztdgsqdk(p1,p2,pl,pa,k5,p6,p7,e5,
     &  mdecay,ampsq16,ampsq25)
      implicit none
C     Matrix element squared for
C     split into two contributions:
C     radiation from 16 and 25 lines (ampsq16 and ampsq25)
C      u(-p1)+b(p2) -> e^-(p3)+e^+(p4)+t(nu(p5)+e(p6)+b(p7))+d(p8)+g(p9)
C     Matrix element is constructed so that dependence on p5
C     only enters through k5 and e5 in wave function
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer p1,p2,p3,p4,p5,k5,p6,p7,e5,jl,jg,icol,t16,t25,pl,pa,j
      parameter(t16=1,t25=2)
      double precision s34,s1467,s16,s167,s1346,s13467,s12346,s126,s1267
      double precision s146,s134,s234,s1347,s2347,mton2mwsq,ampsq
      double precision s346,s3467,theta,cprop,ampsq16,ampsq25
      double complex prW,prt,iza,izb,Amp(2,2,2),nr(2,2)
      double complex zab2,zab3,zab4,mdecay
      double complex facuLl,facuRl,facdLl,iprZ
      double complex prWs16,prWs167,prWs1346,prWs13467,
     & prts126,prts1267,prts12346
c--- definitions of propagators
      theta(s34)=half+half*sign(one,s34)
      prW(s34)=cone/dcmplx(s34-wmass**2,theta(s34)*wmass*wwidth)
      prt(s34)=cone/dcmplx(s34-mt**2,zip)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
c--- definitions of compound spinors
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)
     &                    +za(p1,p4)*zb(p4,p5)
      zab4(p1,p2,p3,p4,p5,p6)=za(p1,p2)*zb(p2,p6)+za(p1,p3)*zb(p3,p6)
     &                       +za(p1,p4)*zb(p4,p6)+za(p1,p5)*zb(p5,p6)
C   End statement functions
      mton2mwsq=mt/(2d0*wmass**2)

      p3=pl
      p4=pa
      s34=s(p3,p4)
      s16=s(p1,p6)
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s167=s(p1,p6)+s(p1,p7)+s(p6,p7)
      s346=s(p3,p4)+s(p3,p6)+s(p4,p6)
      s126=s(p1,p2)+s(p1,p6)+s(p2,p6)
      s1267=s(p1,p2)+s(p1,p6)+s(p1,p7)
     &              +s(p2,p6)+s(p2,p7)
     &                       +s(p6,p7)

      s3467=s(p3,p4)+s(p3,p6)+s(p3,p7)
     &              +s(p4,p6)+s(p4,p7)
     &                       +s(p6,p7)

      s1346=s(p1,p3)+s(p1,p4)+s(p1,p6)
     &              +s(p3,p4)+s(p3,p6)
     &                       +s(p4,p6)
      s1347=s(p1,p3)+s(p1,p4)+s(p1,p7)
     &              +s(p3,p4)+s(p3,p7)
     &                       +s(p4,p7)
      s2347=s(p2,p3)+s(p2,p4)+s(p2,p7)
     &              +s(p3,p4)+s(p3,p7)
     &                       +s(p4,p7)
      s13467=s(p1,p3)+s(p1,p4)+s(p1,p6)+s(p1,p7)
     &               +s(p3,p4)+s(p3,p6)+s(p3,p7)
     &                        +s(p4,p6)+s(p4,p7)
     &                                 +s(p6,p7)
      s12346=s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p1,p6)
     &               +s(p2,p3)+s(p2,p4)+s(p2,p6)
     &                        +s(p3,p4)+s(p3,p6)
     &                                 +s(p4,p6)

      prWs16=prW(s16)
      prWs167=prW(s167)
      prWs1346=prW(s1346)
      prWs13467=prW(s13467)
      prts126=prt(s126)
      prts1267=prt(s1267)
      prts12346=prt(s12346)

c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=1d0/((s34-zmass**2)**2+(zmass*zwidth)**2)
      iprZ=dcmplx(s34-zmass**2)

      cprop=cprop/(mt*twidth)**2

      do j=1,2
      if (j .eq. 1) then
        p3=pl
        p4=pa
        facuLl=dcmplx(Qu*q1)*iprZ+dcmplx(L(2)*le)*s34
        facuRl=dcmplx(Qu*q1)*iprZ+dcmplx(R(2)*le)*s34
        facdLl=dcmplx(Qd*q1)*iprZ+dcmplx(L(1)*le)*s34
      else
        p3=pa
        p4=pl
        facuLl=dcmplx(Qu*q1)*iprZ+dcmplx(L(2)*re)*s34
        facuRl=dcmplx(Qu*q1)*iprZ+dcmplx(R(2)*re)*s34
        facdLl=dcmplx(Qd*q1)*iprZ+dcmplx(L(1)*re)*s34
      endif

      s146=s(p1,p4)+s(p1,p6)+s(p4,p6)
      s1467=s(p1,p4)+s(p1,p6)+s(p1,p7)
     &              +s(p4,p6)+s(p4,p7)
     &                       +s(p6,p7)

      Amp(j,1,t16)= + izb(p1,p7)*prWs13467*s34**(-1)*facdLl*s346**(-1)
     &  * ( za(p6,p3)*zb(p1,p2)*zab2(p7,p3,p6,p4)*zab4(k5,p3,p4,p6,p7,
     &    p1)*s3467**(-1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*prWs13467*s34**(-1)*
     & facuLl*s134**(-1) * (  - za(p6,k5)*za(p3,p4)*zb(p1,p4)**2*zab3(
     &    p7,p1,p3,p4,p2)*s1347**(-1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs13467*
     & s34**(-1)*facdLl * (  - zb(p1,p2)*zab2(p3,p6,p7,p1)*zab3(k5,p3,
     &    p6,p7,p4)*s3467**(-1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs13467*
     & s34**(-1)*facuLl*s134**(-1) * ( zb(p1,p4)*zab2(k5,p6,p7,p1)*
     &    zab2(p3,p1,p4,p2) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & s34**(-1)*facdLl*s234**(-1) * (  - zb(p2,p4)*zab2(k5,p6,p7,p1)*
     &    zab2(p3,p2,p4,p1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zb(p1,p6)*
     &    zab3(p6,p1,p2,p7,p4) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facdLl * (  - za(k5,p3)*zb(p1,p3)*zb(p2,p4)*
     &    zab2(p3,p6,p7,p1) - za(k5,p3)*zb(p1,p4)*zb(p2,p4)*zab2(p4,p6,
     &    p7,p1) - zb(p1,p2)*zab2(k5,p6,p7,p1)*zab3(p3,p1,p6,p7,p4) + 
     &    zb(p1,p4)*zab2(p3,p6,p7,p1)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facuLl * ( za(k5,p3)*zb(p1,p3)*zb(p2,p4)*
     &    zab2(p3,p6,p7,p1) + za(k5,p3)*zb(p1,p4)*zb(p2,p4)*zab2(p4,p6,
     &    p7,p1) + zb(p1,p2)*zab2(k5,p6,p7,p1)*zab3(p3,p1,p6,p7,p4) - 
     &    zb(p1,p4)*zab2(p3,p6,p7,p1)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(e5,k5)*mt
     & *mton2mwsq*prWs167*prWs13467*facdLl * ( zb(p1,p4)*zb(p2,e5)*
     &    zab2(p3,p6,p7,p1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(e5,k5)*mt
     & *mton2mwsq*prWs167*prWs13467*facuLl * (  - zb(p1,p4)*zb(p2,e5)*
     &    zab2(p3,p6,p7,p1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(e5,k5)*
     & mt**2*prWs167*prts1267*s34**(-1)*facuRl * (  - zb(p1,p2)*zb(e5,
     &    p4)*zab2(p3,p6,p7,p1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(e5,k5)*
     & mt**2*prWs167*prWs13467*s34**(-1)*facdLl * (  - zb(p1,p4)*zb(p2,
     &    e5)*zab2(p3,p6,p7,p1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(e5,k5)*
     & mt**2*prWs167*prWs13467*s34**(-1)*facuLl * ( zb(p1,p4)*zb(p2,e5)
     &    *zab2(p3,p6,p7,p1) )
      Amp(j,1,t16) = Amp(j,1,t16) + izb(p6,p7)*prWs167*prts1267*
     & s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zab3(p7,p1,p2,p6,p4)
     &     )

      Amp(j,2,t16)= + iza(p1,p7)*prWs167*s34**(-1)*facdLl*s234**(-1)
     &  * (  - za(p6,k5)*zb(p2,p4)*zab2(p3,p2,p4,p7) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*prWs167*prWs13467*
     & s34**(-1)*facdLl * ( za(k5,p3)*zb(p2,p4)*zab2(p6,p3,p4,p7) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*prWs167*prWs13467*
     & s34**(-1)*facuLl * (  - za(k5,p3)*zb(p2,p4)*zab2(p6,p3,p4,p7) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs13467*
     & s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*zab2(p6,p1,p7,p2)*
     &    zab2(k5,p3,p6,p4) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs13467*
     & s34**(-1)*facuLl * ( za(p6,k5)*zab2(p6,p1,p7,p4)*zab3(p3,p1,p4,
     &    p7,p2)*s1347**(-1) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & s34**(-1)*facdLl*s234**(-1) * ( za(p1,p6)*za(p6,k5)*zb(p2,p4)*
     &    zab2(p3,p2,p4,p1) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zab2(p6,p1,p7,p2)*
     &    zab3(p6,p1,p2,p7,p4) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facdLl * (  - za(p1,p6)*za(k5,p3)*zb(p2,p4)*
     &    zab2(p6,p3,p4,p1) - za(p6,k5)*zab2(p6,p1,p7,p2)*zab3(p3,p1,p6
     &    ,p7,p4) + za(p6,p3)*zab2(p6,p1,p7,p4)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facuLl * ( za(p1,p6)*za(k5,p3)*zb(p2,p4)*
     &    zab2(p6,p3,p4,p1) + za(p6,k5)*zab2(p6,p1,p7,p2)*zab3(p3,p1,p6
     &    ,p7,p4) - za(p6,p3)*zab2(p6,p1,p7,p4)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(e5,k5)*mt
     & *mton2mwsq*prWs167*prWs13467*facdLl * ( za(p6,p3)*zb(p2,e5)*
     &    zab2(p6,p1,p7,p4) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(e5,k5)*mt
     & *mton2mwsq*prWs167*prWs13467*facuLl * (  - za(p6,p3)*zb(p2,e5)*
     &    zab2(p6,p1,p7,p4) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(e5,k5)*
     & mt**2*prWs167*prts1267*s34**(-1)*facuRl * (  - za(p6,p3)*zb(e5,
     &    p4)*zab2(p6,p1,p7,p2) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(e5,k5)*
     & mt**2*prWs167*prWs13467*s34**(-1)*facdLl * (  - za(p6,p3)*zb(p2,
     &    e5)*zab2(p6,p1,p7,p4) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(e5,k5)*
     & mt**2*prWs167*prWs13467*s34**(-1)*facuLl * ( za(p6,p3)*zb(p2,e5)
     &    *zab2(p6,p1,p7,p4) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p6,p7)*prWs13467*s34**(-1)*
     & facdLl*s346**(-1) * (  - za(p6,p3)**2*zb(p1,p2)*zb(p3,p4)*zab3(
     &    k5,p3,p4,p6,p7)*s3467**(-1) )
      Amp(j,2,t16) = Amp(j,2,t16) + iza(p6,p7)*prWs13467*s34**(-1)*
     & facuLl*s134**(-1) * (  - za(p6,k5)*zb(p1,p4)*zab2(p3,p1,p4,p7)*
     &    zab4(p6,p1,p3,p4,p7,p2)*s1347**(-1) )

      Amp(j,1,t25)= + izb(p2,p7)*prWs1346*prts12346*s34**(-1)*facdLl*
     & s346**(-1) * ( za(p6,p3)*za(p7,k5)*zb(p1,p2)*zb(p6,p4)*zab3(p6,
     &    p1,p3,p4,p2) + za(p6,p3)*za(p7,k5)*zb(p1,p2)*zb(p3,p4)*zab3(
     &    p3,p1,p4,p6,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs1346*prts12346*
     & s34**(-1)*facuLl*s134**(-1) * ( za(p7,k5)*zb(p1,p4)*zab2(p3,p1,
     &    p4,p2)*zab3(p6,p1,p3,p4,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs16*s34**(-1)*facdLl*
     & s234**(-1)*s2347**(-1) * ( za(p6,k5)*za(p3,p4)*zb(p2,p4)**2*
     &    zab3(p7,p2,p3,p4,p1) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs16*prts12346*
     & s34**(-1)*facdLl*s234**(-1) * (  - za(p7,k5)*zb(p2,p4)*zab2(p3,
     &    p2,p4,p1)*zab3(p6,p1,p3,p4,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs16*prts126*prts12346
     & *s34**(-1)*facuLl * ( za(p7,k5)*zb(p1,p2)*zab2(p6,p1,p2,p4)*
     &    zab3(p3,p1,p4,p6,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs16*prts126*prts1267*
     & s34**(-1)*facuLl * ( za(p1,p6)*za(k5,p3)*zb(p1,p2)**2*zab3(p7,p1
     &    ,p2,p6,p4) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs16*prWs1346*
     & prts12346*s34**(-1)*facdLl * (  - za(p6,p3)*za(p7,k5)*zb(p1,p2)*
     &    zb(p1,p4)*zab3(p1,p3,p4,p6,p2) + za(p6,p3)*za(p7,k5)*zb(p1,p4
     &    )*zb(p2,p6)*zab3(p6,p1,p3,p4,p2) - za(p7,k5)*zb(p1,p2)*zab2(
     &    p3,p1,p6,p4)*zab3(p6,p1,p3,p4,p2) - za(p7,k5)*zb(p2,p4)*zab2(
     &    p6,p3,p4,p1)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*prWs16*prWs1346*
     & prts12346*s34**(-1)*facuLl * ( za(p6,p3)*za(p7,k5)*zb(p1,p2)*zb(
     &    p1,p4)*zab3(p1,p3,p4,p6,p2) - za(p6,p3)*za(p7,k5)*zb(p1,p4)*
     &    zb(p2,p6)*zab3(p6,p1,p3,p4,p2) + za(p7,k5)*zb(p1,p2)*zab2(p3,
     &    p1,p6,p4)*zab3(p6,p1,p3,p4,p2) + za(p7,k5)*zb(p2,p4)*zab2(p6,
     &    p3,p4,p1)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*mt**2*prWs16*prts126*
     & prts12346*s34**(-1)*facuRl * (  - za(p6,p3)*za(p7,k5)*zb(p1,p2)*
     &    zb(p2,p4) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*mt**2*prWs16*prts126*
     & prts1267*s34**(-1)*facuLl * (  - za(p6,p7)*za(k5,p3)*zb(p1,p2)*
     &    zb(p2,p4) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt*mton2mwsq*
     & prWs16*prWs1346*prts12346*facdLl * ( za(p6,p3)*zb(p1,p4)*zb(p2,
     &    e5)*zab4(p7,p1,p3,p4,p6,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt*mton2mwsq*
     & prWs16*prWs1346*prts12346*facuLl * (  - za(p6,p3)*zb(p1,p4)*zb(
     &    p2,e5)*zab4(p7,p1,p3,p4,p6,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*
     &    zb(p1,p2)*zb(p2,e5)*zab2(p7,p3,p6,p4) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facuLl*s134**(-1) * ( za(p6,p7)*zb(
     &    p1,p4)*zb(p2,e5)*zab2(p3,p1,p4,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*prWs16*
     & prts12346*s34**(-1)*facdLl*s234**(-1) * (  - za(p6,p7)*zb(p2,e5)
     &    *zb(p2,p4)*zab2(p3,p2,p4,p1) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*prWs16*
     & prts126*prts12346*s34**(-1)*facuRl * ( za(p6,p3)*zb(p1,p2)*zb(p2
     &    ,e5)*zab4(p7,p1,p2,p3,p6,p4) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*prWs16*
     & prts126*prts12346*s34**(-1)*facuLl * (  - za(p7,p3)*zb(p1,p2)*
     &    zb(p2,e5)*zab2(p6,p1,p2,p4) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*prWs16*
     & prts126*prts1267*s34**(-1)*facuRl * ( za(p1,p6)*za(p7,p3)*zb(p1,
     &    p2)**2*zb(e5,p4) - za(p6,p7)*zb(p1,p2)*zb(e5,p4)*zab3(p3,p1,
     &    p6,p7,p2) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*prWs16*
     & prWs1346*prts12346*s34**(-1)*facdLl * (  - za(p6,p7)*zb(p1,p2)*
     &    zb(p2,e5)*zab2(p3,p1,p6,p4) + za(p6,p3)*zb(p1,p4)*zb(p2,e5)*
     &    zab2(p7,p1,p6,p2) - za(p6,p3)*zb(p1,p4)*zb(p2,e5)*zab4(p7,p1,
     &    p3,p4,p6,p2) + za(p7,p3)*zb(p2,e5)*zb(p2,p4)*zab2(p6,p3,p4,p1
     &    ) )
      Amp(j,1,t25) = Amp(j,1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*prWs16*
     & prWs1346*prts12346*s34**(-1)*facuLl * ( za(p6,p7)*zb(p1,p2)*zb(
     &    p2,e5)*zab2(p3,p1,p6,p4) - za(p6,p3)*zb(p1,p4)*zb(p2,e5)*
     &    zab2(p7,p1,p6,p2) + za(p6,p3)*zb(p1,p4)*zb(p2,e5)*zab4(p7,p1,
     &    p3,p4,p6,p2) - za(p7,p3)*zb(p2,e5)*zb(p2,p4)*zab2(p6,p3,p4,p1
     &    ) )

      Amp(j,2,t25)= + iza(p2,p7)*prWs1346*s34**(-1)*facuLl*s134**(-1)
     &  * ( za(p6,k5)*zb(p1,p4)*zab2(p3,p1,p4,p7) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*prWs16*prWs1346*
     & s34**(-1)*facdLl * ( za(p6,p3)*zb(p1,p4)*zab2(k5,p1,p6,p7) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*prWs16*prWs1346*
     & s34**(-1)*facuLl * (  - za(p6,p3)*zb(p1,p4)*zab2(k5,p1,p6,p7) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs1346*
     & s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*zab2(k5,p2,p7,p1)*
     &    zab2(k5,p3,p6,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs1346*
     & s34**(-1)*facuLl*s134**(-1) * ( za(p2,k5)*za(p6,k5)*zb(p1,p4)*
     &    zab2(p3,p1,p4,p2) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & s34**(-1)*facdLl*s2347**(-1) * ( za(p6,k5)*zab2(k5,p2,p7,p4)*
     &    zab3(p3,p2,p4,p7,p1) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zab2(k5,p2,p7,p1)*
     &    zab3(p6,p1,p2,p7,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & prWs1346*s34**(-1)*facdLl * ( za(p2,k5)*za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p1,p6,p2) - za(p6,k5)*zab2(k5,p2,p7,p1)*zab2(p3,p1,p6
     &    ,p4) - za(k5,p3)*zab2(p6,p3,p4,p1)*zab2(k5,p2,p7,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & prWs1346*s34**(-1)*facuLl * (  - za(p2,k5)*za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p1,p6,p2) + za(p6,k5)*zab2(k5,p2,p7,p1)*zab2(p3,p1,p6
     &    ,p4) + za(k5,p3)*zab2(p6,p3,p4,p1)*zab2(k5,p2,p7,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(e5,k5)*mt
     & *mton2mwsq*prWs16*prWs1346*facdLl * (  - za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p2,p7,e5) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(e5,k5)*mt
     & *mton2mwsq*prWs16*prWs1346*facuLl * ( za(p6,p3)*zb(p1,p4)*zab2(
     &    k5,p2,p7,e5) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(e5,k5)*
     & mt**2*prWs16*prts1267*s34**(-1)*facuRl * (  - za(p6,p3)*zb(e5,p4
     &    )*zab2(k5,p2,p7,p1) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(e5,k5)*
     & mt**2*prWs16*prWs1346*s34**(-1)*facdLl * ( za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p2,p7,e5) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(e5,k5)*
     & mt**2*prWs16*prWs1346*s34**(-1)*facuLl * (  - za(p6,p3)*zb(p1,p4
     &    )*zab2(k5,p2,p7,e5) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*prWs16*s34**(-1)*facdLl*
     & s234**(-1)*s2347**(-1) * (  - za(p6,k5)*zb(p2,p4)*zab2(p3,p2,p4,
     &    p7)*zab4(k5,p2,p3,p4,p7,p1) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*prWs16*prts126*prts1267*
     & s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zab2(p6,p1,p2,p7)*
     &    zab4(k5,p1,p2,p6,p7,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*mt**2*prWs16*prts126*
     & prts1267*s34**(-1)*facuLl * (  - za(p6,k5)*za(k5,p3)*zb(p1,p2)*
     &    zb(p7,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt*mton2mwsq*
     & prWs16*prWs1346*prts12346*facdLl * ( za(p6,p3)*zb(p1,p4)*zb(p7,
     &    e5)*zab4(k5,p1,p3,p4,p6,p2) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt*mton2mwsq*
     & prWs16*prWs1346*prts12346*facuLl * (  - za(p6,p3)*zb(p1,p4)*zb(
     &    p7,e5)*zab4(k5,p1,p3,p4,p6,p2) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*
     &    zb(p1,p2)*zb(p7,e5)*zab2(k5,p3,p6,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facuLl*s134**(-1) * ( za(p6,k5)*zb(
     &    p1,p4)*zb(p7,e5)*zab2(p3,p1,p4,p2) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*prWs16*
     & prts12346*s34**(-1)*facdLl*s234**(-1) * (  - za(p6,k5)*zb(p2,p4)
     &    *zb(p7,e5)*zab2(p3,p2,p4,p1) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*prWs16*
     & prts126*prts12346*s34**(-1)*facuRl * ( za(p6,p3)*zb(p1,p2)*zb(p7
     &    ,e5)*zab4(k5,p1,p2,p3,p6,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*prWs16*
     & prts126*prts12346*s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*
     &    zb(p7,e5)*zab2(p6,p1,p2,p4) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*prWs16*
     & prts126*prts1267*s34**(-1)*facuRl * (  - za(p6,k5)*zb(p1,p2)*zb(
     &    e5,p4)*zab3(p3,p1,p2,p6,p7) - za(k5,p3)*zb(p1,p2)*zb(e5,p4)*
     &    zab2(p6,p1,p2,p7) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*prWs16*
     & prWs1346*prts12346*s34**(-1)*facdLl * (  - za(p6,k5)*zb(p1,p2)*
     &    zb(p7,e5)*zab2(p3,p1,p6,p4) + za(p6,p3)*zb(p1,p4)*zb(p7,e5)*
     &    zab2(k5,p1,p6,p2) - za(p6,p3)*zb(p1,p4)*zb(p7,e5)*zab4(k5,p1,
     &    p3,p4,p6,p2) + za(k5,p3)*zb(p2,p4)*zb(p7,e5)*zab2(p6,p3,p4,p1
     &    ) )
      Amp(j,2,t25) = Amp(j,2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*prWs16*
     & prWs1346*prts12346*s34**(-1)*facuLl * ( za(p6,k5)*zb(p1,p2)*zb(
     &    p7,e5)*zab2(p3,p1,p6,p4) - za(p6,p3)*zb(p1,p4)*zb(p7,e5)*
     &    zab2(k5,p1,p6,p2) + za(p6,p3)*zb(p1,p4)*zb(p7,e5)*zab4(k5,p1,
     &    p3,p4,p6,p2) - za(k5,p3)*zb(p2,p4)*zb(p7,e5)*zab2(p6,p3,p4,p1
     &    ) )


c--- Non-resonant diagrams for left-handed lepton helicity only
      if (j. eq. 1) then
      nr(1,t16)= + izb(p1,p7)*izb(p6,p7)*xw**(-1)*iprZ*prWs167*
     & prWs13467 * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p6)*zb(p1,p4)*zab3(p6,p1
     &    ,p4,p7,p2)*s1467**(-1) )
      nr(1,t16) = nr(1,t16) + izb(p1,p7)*izb(p6,p7)*izb(e5,k5)*mt*
     & mton2mwsq*xw**(-1)*iprZ*prWs167*prWs13467 * (  - 1.D0/2.D0*zb(p1
     &    ,p4)*zb(p2,e5)*zab2(p3,p6,p7,p1) )
      nr(1,t16) = nr(1,t16) + izb(p6,p7)*xw**(-1)*iprZ*prWs167*
     & prWs13467 * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p4)*zab3(p7,p1,p4,p6,p2)
     &    *s1467**(-1) )

      nr(2,t16)= + iza(p1,p7)*iza(p6,p7)*xw**(-1)*iprZ*prWs167*
     & prWs13467 * ( 1.D0/2.D0*za(k5,p3)*zab2(p6,p1,p7,p4)*zab3(p6,p1,
     &    p4,p7,p2)*s1467**(-1) )
      nr(2,t16) = nr(2,t16) + iza(p1,p7)*iza(p6,p7)*izb(e5,k5)*mt*
     & mton2mwsq*xw**(-1)*iprZ*prWs167*prWs13467 * (  - 1.D0/2.D0*za(p6
     &    ,p3)*zb(p2,e5)*zab2(p6,p1,p7,p4) )

      nr(1,t25)= + izb(p2,p7)*xw**(-1)*iprZ*prWs16*prWs1346*prts12346*
     & s146**(-1) * (  - 1.D0/2.D0*za(p7,k5)*zb(p1,p4)*zab2(p6,p1,p4,p2
     &    )*zab3(p3,p1,p4,p6,p2) )
      nr(1,t25) = nr(1,t25) + izb(p2,p7)*izb(e5,k5)*mt*mton2mwsq*
     & xw**(-1)*iprZ*prWs16*prWs1346*prts12346 * (  - 1.D0/2.D0*za(p6,
     &    p3)*zb(p1,p4)*zb(p2,e5)*zab4(p7,p1,p3,p4,p6,p2) )
      nr(1,t25) = nr(1,t25) + izb(p2,p7)*izb(e5,k5)*mt**2*xw**(-1)*iprZ
     & *prWs16*prWs1346*prts12346*s146**(-1) * ( 1.D0/2.D0*za(p7,p3)*
     &    zb(p1,p4)*zb(p2,e5)*zab2(p6,p1,p4,p2) )

      nr(2,t25)= + iza(p2,p7)*xw**(-1)*iprZ*prWs16*prWs1346*s146**(-1)
     &  * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p4)*zab2(p6,p1,p4,p7) )
      nr(2,t25) = nr(2,t25) + iza(p2,p7)*iza(p7,k5)*xw**(-1)*iprZ*
     & prWs16*prWs1346*s146**(-1) * ( 1.D0/2.D0*za(p2,k5)*za(k5,p3)*zb(
     &    p1,p4)*zab2(p6,p1,p4,p2) )
      nr(2,t25) = nr(2,t25) + iza(p2,p7)*iza(p7,k5)*izb(e5,k5)*mt*
     & mton2mwsq*xw**(-1)*iprZ*prWs16*prWs1346 * ( 1.D0/2.D0*za(p6,p3)*
     &    zb(p1,p4)*zab2(k5,p2,p7,e5) )
      nr(2,t25) = nr(2,t25) + iza(p7,k5)*izb(e5,k5)*mt*mton2mwsq*
     & xw**(-1)*iprZ*prWs16*prWs1346*prts12346 * (  - 1.D0/2.D0*za(p6,
     &    p3)*zb(p1,p4)*zb(p7,e5)*zab4(k5,p1,p3,p4,p6,p2) )
      nr(2,t25) = nr(2,t25) + iza(p7,k5)*izb(e5,k5)*mt**2*xw**(-1)*iprZ
     & *prWs16*prWs1346*prts12346*s146**(-1) * ( 1.D0/2.D0*za(k5,p3)*
     &    zb(p1,p4)*zb(p7,e5)*zab2(p6,p1,p4,p2) )

      do jg=1,2
      do icol=1,2
      Amp(1,jg,icol)=Amp(1,jg,icol)+nr(jg,icol)
      enddo
      enddo
      endif

      enddo

      ampsq16=0d0
      ampsq25=0d0
      do jl=1,2
      do jg=1,2
      ampsq16=ampsq16+cdabs(Amp(jl,jg,t16)*mdecay)**2
      ampsq25=ampsq25+cdabs(Amp(jl,jg,t25)*mdecay)**2
      enddo
      enddo
      ampsq16=ampsq16*cprop
      ampsq25=ampsq25*cprop

      return
      end
