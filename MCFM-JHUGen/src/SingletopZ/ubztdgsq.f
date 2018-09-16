      subroutine ubztdgsq(p1,p2,pl,pa,k5a,k5b,p6,p7,ampsq16,ampsq25)
      implicit none
      include 'types.f'

C     Matrix element squared for u(p1)+b(p2)->Z(p3,p4)+t(p5)+d(p6)+g(p7)
C     split into two contributions:
C     radiation from 16 and 25 lines (ampsq16 and ampsq25)
C     Matrix element is constructed so that dependence on p5 only enters
C     through k5; p5 is made massless wrt both pl (k5a) and pa (k5b)
C     in order to facilitate easy calculation of both lepton helicities
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,p5,k5,k5a,k5b,p6,p7,jl,jg,
     & icol,t16,t25,pl,pa,j,jt
      parameter(t16=1,t25=2)
      real(dp):: s34,s1467,s16,s167,s1346,s13467,s12346,s126,s1267
      real(dp):: s146,s134,s234,s1347,s2347,mton2mwsq,ampsq
      real(dp):: s346,s3467,theta,cprop,ampsq16,ampsq25
      complex(dp):: prW,prt,iza,izb,Amp(2,2,2,2),nr(2,2,2)
      complex(dp):: zab2,zab3,zab4
      complex(dp):: facuLl,facuRl,facdLl,iprZ
      complex(dp):: prWs16,prWs167,prWs1346,prWs13467,
     & prts126,prts1267,prts12346
c--- definitions of propagators
      theta(s34)=half+half*sign(one,s34)
      prW(s34)=cone/cplx2(s34-wmass**2,theta(s34)*wmass*wwidth)
      prt(s34)=cone/cplx2(s34-mt**2,zip)
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
      iprZ=cplx1(s34-zmass**2)
      do j=1,2
      if (j == 1) then
        p3=pl
        p4=pa
        k5=k5a    ! massless wrt p3=pl
        facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*le)*s34
        facuRl=cplx1(Qu*q1)*iprZ+cplx1(R(2)*le)*s34
        facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*le)*s34
      else
        p3=pa
        p4=pl
        k5=k5b    ! massless wrt p3=pa
        facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*re)*s34
        facuRl=cplx1(Qu*q1)*iprZ+cplx1(R(2)*re)*s34
        facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*re)*s34
      endif

      s146=s(p1,p4)+s(p1,p6)+s(p4,p6)
      s1467=s(p1,p4)+s(p1,p6)+s(p1,p7)
     &              +s(p4,p6)+s(p4,p7)
     &                       +s(p6,p7)

c--- ordering of labels is as follows:
c--- lepton helicity, top spin, gluon helicity, color
      Amp(j,1,1,t16)= + izb(p1,p7)*prWs13467*s34**(-1)*facdLl*
     & s346**(-1) * ( za(p6,p3)*zb(p1,p2)*zab2(p7,p3,p6,p4)*zab4(k5,p3,
     &    p4,p6,p7,p1)*s3467**(-1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*prWs13467*s34**(-1)*
     & facuLl*s134**(-1) * (  - za(p6,k5)*za(p3,p4)*zb(p1,p4)**2*zab3(
     &    p7,p1,p3,p4,p2)*s1347**(-1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs13467
     & *s34**(-1)*facdLl * (  - zb(p1,p2)*zab2(p3,p6,p7,p1)*zab3(k5,p3,
     &    p6,p7,p4)*s3467**(-1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs13467
     & *s34**(-1)*facuLl*s134**(-1) * ( zb(p1,p4)*zab2(k5,p6,p7,p1)*
     &    zab2(p3,p1,p4,p2) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & s34**(-1)*facdLl*s234**(-1) * (  - zb(p2,p4)*zab2(k5,p6,p7,p1)*
     &    zab2(p3,p2,p4,p1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zb(p1,p6)*
     &    zab3(p6,p1,p2,p7,p4) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facdLl * (  - za(k5,p3)*zb(p1,p3)*zb(p2,p4)*
     &    zab2(p3,p6,p7,p1) - za(k5,p3)*zb(p1,p4)*zb(p2,p4)*zab2(p4,p6,
     &    p7,p1) - zb(p1,p2)*zab2(k5,p6,p7,p1)*zab3(p3,p1,p6,p7,p4) +
     &    zb(p1,p4)*zab2(p3,p6,p7,p1)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facuLl * ( za(k5,p3)*zb(p1,p3)*zb(p2,p4)*
     &    zab2(p3,p6,p7,p1) + za(k5,p3)*zb(p1,p4)*zb(p2,p4)*zab2(p4,p6,
     &    p7,p1) + zb(p1,p2)*zab2(k5,p6,p7,p1)*zab3(p3,p1,p6,p7,p4) -
     &    zb(p1,p4)*zab2(p3,p6,p7,p1)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(k5,p3
     & )*mt*mton2mwsq*prWs167*prWs13467*facdLl * (  - zb(p1,p4)*zb(p2,
     &    p3)*zab2(p3,p6,p7,p1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(k5,p3
     & )*mt*mton2mwsq*prWs167*prWs13467*facuLl * ( zb(p1,p4)*zb(p2,p3)*
     &    zab2(p3,p6,p7,p1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(k5,p3
     & )*mt**2*prWs167*prts1267*s34**(-1)*facuRl * ( zb(p1,p2)*zb(p3,p4
     &    )*zab2(p3,p6,p7,p1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(k5,p3
     & )*mt**2*prWs167*prWs13467*s34**(-1)*facdLl * ( zb(p1,p4)*zb(p2,
     &    p3)*zab2(p3,p6,p7,p1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(k5,p3
     & )*mt**2*prWs167*prWs13467*s34**(-1)*facuLl * (  - zb(p1,p4)*zb(
     &    p2,p3)*zab2(p3,p6,p7,p1) )
      Amp(j,1,1,t16) = Amp(j,1,1,t16) + izb(p6,p7)*prWs167*prts1267*
     & s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zab3(p7,p1,p2,p6,p4)
     &     )

      Amp(j,1,2,t16)= + iza(p1,p7)*prWs167*s34**(-1)*facdLl*s234**(-1)
     &  * (  - za(p6,k5)*zb(p2,p4)*zab2(p3,p2,p4,p7) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*prWs167*prWs13467*
     & s34**(-1)*facdLl * ( za(k5,p3)*zb(p2,p4)*zab2(p6,p3,p4,p7) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*prWs167*prWs13467*
     & s34**(-1)*facuLl * (  - za(k5,p3)*zb(p2,p4)*zab2(p6,p3,p4,p7) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs13467
     & *s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*zab2(p6,p1,p7,p2)*
     &    zab2(k5,p3,p6,p4) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs13467
     & *s34**(-1)*facuLl * ( za(p6,k5)*zab2(p6,p1,p7,p4)*zab3(p3,p1,p4,
     &    p7,p2)*s1347**(-1) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & s34**(-1)*facdLl*s234**(-1) * ( za(p1,p6)*za(p6,k5)*zb(p2,p4)*
     &    zab2(p3,p2,p4,p1) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zab2(p6,p1,p7,p2)*
     &    zab3(p6,p1,p2,p7,p4) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facdLl * (  - za(p1,p6)*za(k5,p3)*zb(p2,p4)*
     &    zab2(p6,p3,p4,p1) - za(p6,k5)*zab2(p6,p1,p7,p2)*zab3(p3,p1,p6
     &    ,p7,p4) + za(p6,p3)*zab2(p6,p1,p7,p4)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*prWs167*
     & prWs13467*s34**(-1)*facuLl * ( za(p1,p6)*za(k5,p3)*zb(p2,p4)*
     &    zab2(p6,p3,p4,p1) + za(p6,k5)*zab2(p6,p1,p7,p2)*zab3(p3,p1,p6
     &    ,p7,p4) - za(p6,p3)*zab2(p6,p1,p7,p4)*zab3(k5,p1,p6,p7,p2) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(k5,p3
     & )*mt*mton2mwsq*prWs167*prWs13467*facdLl * (  - za(p6,p3)*zb(p2,
     &    p3)*zab2(p6,p1,p7,p4) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(k5,p3
     & )*mt*mton2mwsq*prWs167*prWs13467*facuLl * ( za(p6,p3)*zb(p2,p3)*
     &    zab2(p6,p1,p7,p4) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(k5,p3
     & )*mt**2*prWs167*prts1267*s34**(-1)*facuRl * ( za(p6,p3)*zb(p3,p4
     &    )*zab2(p6,p1,p7,p2) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(k5,p3
     & )*mt**2*prWs167*prWs13467*s34**(-1)*facdLl * ( za(p6,p3)*zb(p2,
     &    p3)*zab2(p6,p1,p7,p4) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(k5,p3
     & )*mt**2*prWs167*prWs13467*s34**(-1)*facuLl * (  - za(p6,p3)*zb(
     &    p2,p3)*zab2(p6,p1,p7,p4) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p6,p7)*prWs13467*s34**(-1)*
     & facdLl*s346**(-1) * (  - za(p6,p3)**2*zb(p1,p2)*zb(p3,p4)*zab3(
     &    k5,p3,p4,p6,p7)*s3467**(-1) )
      Amp(j,1,2,t16) = Amp(j,1,2,t16) + iza(p6,p7)*prWs13467*s34**(-1)*
     & facuLl*s134**(-1) * (  - za(p6,k5)*zb(p1,p4)*zab2(p3,p1,p4,p7)*
     &    zab4(p6,p1,p3,p4,p7,p2)*s1347**(-1) )

      Amp(j,2,1,t16)= + iza(k5,p3)*izb(p1,p7)*mt*prWs13467*s34**(-1)*
     & facdLl*s346**(-1) * (  - za(p6,p3)*zb(p1,p2)*zab2(p7,p3,p6,p4)*
     &    zab4(p3,p3,p4,p6,p7,p1)*s3467**(-1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + iza(k5,p3)*izb(p1,p7)*mt*
     & prWs13467*s34**(-1)*facuLl*s134**(-1) * ( za(p6,p3)*za(p3,p4)*
     &    zb(p1,p4)**2*zab3(p7,p1,p3,p4,p2)*s1347**(-1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + iza(k5,p3)*izb(p1,p7)*izb(p6,p7
     & )*mt*prWs13467*s34**(-1)*facdLl * ( zb(p1,p2)*zab2(p3,p6,p7,p1)*
     &    zab3(p3,p3,p6,p7,p4)*s3467**(-1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + iza(k5,p3)*izb(p1,p7)*izb(p6,p7
     & )*mt*prWs13467*s34**(-1)*facuLl*s134**(-1) * (  - zb(p1,p4)*
     &    zab2(p3,p1,p4,p2)*zab2(p3,p6,p7,p1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + iza(k5,p3)*izb(p1,p7)*izb(p6,p7
     & )*mt*prWs167*s34**(-1)*facdLl*s234**(-1) * ( zb(p2,p4)*zab2(p3,
     &    p2,p4,p1)*zab2(p3,p6,p7,p1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + iza(k5,p3)*izb(p1,p7)*izb(p6,p7
     & )*mt*prWs167*prWs13467*s34**(-1)*facdLl * ( zb(p1,p2)*zab2(p3,p6
     &    ,p7,p1)*zab3(p3,p1,p6,p7,p4) - zb(p1,p4)*zab2(p3,p6,p7,p1)*
     &    zab3(p3,p1,p6,p7,p2) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + iza(k5,p3)*izb(p1,p7)*izb(p6,p7
     & )*mt*prWs167*prWs13467*s34**(-1)*facuLl * (  - zb(p1,p2)*zab2(p3
     &    ,p6,p7,p1)*zab3(p3,p1,p6,p7,p4) + zb(p1,p4)*zab2(p3,p6,p7,p1)
     &    *zab3(p3,p1,p6,p7,p2) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + izb(p1,p7)*izb(p6,p7)*mton2mwsq
     & *prWs167*prWs13467*facdLl * ( zb(p1,p4)*zb(p2,k5)*zab2(p3,p6,p7,
     &    p1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + izb(p1,p7)*izb(p6,p7)*mton2mwsq
     & *prWs167*prWs13467*facuLl * (  - zb(p1,p4)*zb(p2,k5)*zab2(p3,p6,
     &    p7,p1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + izb(p1,p7)*izb(p6,p7)*mt*
     & prWs167*prts1267*s34**(-1)*facuRl * (  - zb(p1,p2)*zb(k5,p4)*
     &    zab2(p3,p6,p7,p1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + izb(p1,p7)*izb(p6,p7)*mt*
     & prWs167*prWs13467*s34**(-1)*facdLl * (  - zb(p1,p4)*zb(p2,k5)*
     &    zab2(p3,p6,p7,p1) )
      Amp(j,2,1,t16) = Amp(j,2,1,t16) + izb(p1,p7)*izb(p6,p7)*mt*
     & prWs167*prWs13467*s34**(-1)*facuLl * ( zb(p1,p4)*zb(p2,k5)*zab2(
     &    p3,p6,p7,p1) )

      Amp(j,2,2,t16)= + iza(p1,p7)*iza(p6,p7)*mton2mwsq*prWs167*
     & prWs13467*facdLl * ( za(p6,p3)*zb(p2,k5)*zab2(p6,p1,p7,p4) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*mton2mwsq
     & *prWs167*prWs13467*facuLl * (  - za(p6,p3)*zb(p2,k5)*zab2(p6,p1,
     &    p7,p4) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*mt*
     & prWs167*prts1267*s34**(-1)*facuRl * (  - za(p6,p3)*zb(k5,p4)*
     &    zab2(p6,p1,p7,p2) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*mt*
     & prWs167*prWs13467*s34**(-1)*facdLl * (  - za(p6,p3)*zb(p2,k5)*
     &    zab2(p6,p1,p7,p4) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*mt*
     & prWs167*prWs13467*s34**(-1)*facuLl * ( za(p6,p3)*zb(p2,k5)*zab2(
     &    p6,p1,p7,p4) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*iza(k5,p3
     & )*mt*prWs13467*s34**(-1)*facdLl*s346**(-1) * ( za(p6,p3)*zab2(p6
     &    ,p1,p7,p2)*zab2(p3,p3,p6,p4) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*iza(k5,p3
     & )*mt*prWs13467*s34**(-1)*facuLl * (  - za(p6,p3)*zab2(p6,p1,p7,
     &    p4)*zab3(p3,p1,p4,p7,p2)*s1347**(-1) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*iza(k5,p3
     & )*mt*prWs167*s34**(-1)*facdLl*s234**(-1) * (  - za(p1,p6)*za(p6,
     &    p3)*zb(p2,p4)*zab2(p3,p2,p4,p1) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*iza(k5,p3
     & )*mt*prWs167*prWs13467*s34**(-1)*facdLl * ( za(p6,p3)*zab2(p6,p1
     &    ,p7,p2)*zab3(p3,p1,p6,p7,p4) - za(p6,p3)*zab2(p6,p1,p7,p4)*
     &    zab3(p3,p1,p6,p7,p2) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(p6,p7)*iza(k5,p3
     & )*mt*prWs167*prWs13467*s34**(-1)*facuLl * (  - za(p6,p3)*zab2(p6
     &    ,p1,p7,p2)*zab3(p3,p1,p6,p7,p4) + za(p6,p3)*zab2(p6,p1,p7,p4)
     &    *zab3(p3,p1,p6,p7,p2) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p1,p7)*iza(k5,p3)*mt*
     & prWs167*s34**(-1)*facdLl*s234**(-1) * ( za(p6,p3)*zb(p2,p4)*
     &    zab2(p3,p2,p4,p7) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p6,p7)*iza(k5,p3)*mt*
     & prWs13467*s34**(-1)*facdLl*s346**(-1) * ( za(p6,p3)**2*zb(p1,p2)
     &    *zb(p3,p4)*zab3(p3,p3,p4,p6,p7)*s3467**(-1) )
      Amp(j,2,2,t16) = Amp(j,2,2,t16) + iza(p6,p7)*iza(k5,p3)*mt*
     & prWs13467*s34**(-1)*facuLl*s134**(-1) * ( za(p6,p3)*zb(p1,p4)*
     &    zab2(p3,p1,p4,p7)*zab4(p6,p1,p3,p4,p7,p2)*s1347**(-1) )


      Amp(j,1,1,t25)= + izb(p2,p7)*prWs1346*prts12346*s34**(-1)*facdLl*
     & s346**(-1) * ( za(p6,p3)*za(p7,k5)*zb(p1,p2)*zb(p6,p4)*zab3(p6,
     &    p1,p3,p4,p2) + za(p6,p3)*za(p7,k5)*zb(p1,p2)*zb(p3,p4)*zab3(
     &    p3,p1,p4,p6,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs1346*prts12346*
     & s34**(-1)*facuLl*s134**(-1) * ( za(p7,k5)*zb(p1,p4)*zab2(p3,p1,
     &    p4,p2)*zab3(p6,p1,p3,p4,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs16*s34**(-1)*
     & facdLl*s234**(-1)*s2347**(-1) * ( za(p6,k5)*za(p3,p4)*zb(p2,p4)
     &    **2*zab3(p7,p2,p3,p4,p1) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs16*prts12346*
     & s34**(-1)*facdLl*s234**(-1) * (  - za(p7,k5)*zb(p2,p4)*zab2(p3,
     &    p2,p4,p1)*zab3(p6,p1,p3,p4,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs16*prts126*
     & prts12346*s34**(-1)*facuLl * ( za(p7,k5)*zb(p1,p2)*zab2(p6,p1,p2
     &    ,p4)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs16*prts126*
     & prts1267*s34**(-1)*facuLl * ( za(p1,p6)*za(k5,p3)*zb(p1,p2)**2*
     &    zab3(p7,p1,p2,p6,p4) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs16*prWs1346*
     & prts12346*s34**(-1)*facdLl * (  - za(p6,p3)*za(p7,k5)*zb(p1,p2)*
     &    zb(p1,p4)*zab3(p1,p3,p4,p6,p2) + za(p6,p3)*za(p7,k5)*zb(p1,p4
     &    )*zb(p2,p6)*zab3(p6,p1,p3,p4,p2) - za(p7,k5)*zb(p1,p2)*zab2(
     &    p3,p1,p6,p4)*zab3(p6,p1,p3,p4,p2) - za(p7,k5)*zb(p2,p4)*zab2(
     &    p6,p3,p4,p1)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*prWs16*prWs1346*
     & prts12346*s34**(-1)*facuLl * ( za(p6,p3)*za(p7,k5)*zb(p1,p2)*zb(
     &    p1,p4)*zab3(p1,p3,p4,p6,p2) - za(p6,p3)*za(p7,k5)*zb(p1,p4)*
     &    zb(p2,p6)*zab3(p6,p1,p3,p4,p2) + za(p7,k5)*zb(p1,p2)*zab2(p3,
     &    p1,p6,p4)*zab3(p6,p1,p3,p4,p2) + za(p7,k5)*zb(p2,p4)*zab2(p6,
     &    p3,p4,p1)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*mt**2*prWs16*prts126
     & *prts12346*s34**(-1)*facuRl * (  - za(p6,p3)*za(p7,k5)*zb(p1,p2)
     &    *zb(p2,p4) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*mt**2*prWs16*prts126
     & *prts1267*s34**(-1)*facuLl * (  - za(p6,p7)*za(k5,p3)*zb(p1,p2)*
     &    zb(p2,p4) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt*
     & mton2mwsq*prWs16*prWs1346*prts12346*facdLl * (  - za(p6,p3)*zb(
     &    p1,p4)*zb(p2,p3)*zab4(p7,p1,p3,p4,p6,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt*
     & mton2mwsq*prWs16*prWs1346*prts12346*facuLl * ( za(p6,p3)*zb(p1,
     &    p4)*zb(p2,p3)*zab4(p7,p1,p3,p4,p6,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facdLl*s346**(-1) * ( za(p6,p3)*zb(
     &    p1,p2)*zb(p2,p3)*zab2(p7,p3,p6,p4) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facuLl*s134**(-1) * (  - za(p6,p7)*
     &    zb(p1,p4)*zb(p2,p3)*zab2(p3,p1,p4,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs16*prts12346*s34**(-1)*facdLl*s234**(-1) * ( za(p6,p7)*zb(p2
     &    ,p3)*zb(p2,p4)*zab2(p3,p2,p4,p1) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs16*prts126*prts12346*s34**(-1)*facuRl * (  - za(p6,p3)*zb(p1
     &    ,p2)*zb(p2,p3)*zab4(p7,p1,p2,p3,p6,p4) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs16*prts126*prts12346*s34**(-1)*facuLl * ( za(p7,p3)*zb(p1,p2
     &    )*zb(p2,p3)*zab2(p6,p1,p2,p4) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs16*prts126*prts1267*s34**(-1)*facuRl * (  - za(p1,p6)*za(p7,
     &    p3)*zb(p1,p2)**2*zb(p3,p4) + za(p6,p7)*zb(p1,p2)*zb(p3,p4)*
     &    zab3(p3,p1,p6,p7,p2) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs16*prWs1346*prts12346*s34**(-1)*facdLl * ( za(p6,p7)*zb(p1,
     &    p2)*zb(p2,p3)*zab2(p3,p1,p6,p4) - za(p6,p3)*zb(p1,p4)*zb(p2,
     &    p3)*zab2(p7,p1,p6,p2) + za(p6,p3)*zb(p1,p4)*zb(p2,p3)*zab4(p7
     &    ,p1,p3,p4,p6,p2) - za(p7,p3)*zb(p2,p3)*zb(p2,p4)*zab2(p6,p3,
     &    p4,p1) )
      Amp(j,1,1,t25) = Amp(j,1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*
     & prWs16*prWs1346*prts12346*s34**(-1)*facuLl * (  - za(p6,p7)*zb(
     &    p1,p2)*zb(p2,p3)*zab2(p3,p1,p6,p4) + za(p6,p3)*zb(p1,p4)*zb(
     &    p2,p3)*zab2(p7,p1,p6,p2) - za(p6,p3)*zb(p1,p4)*zb(p2,p3)*
     &    zab4(p7,p1,p3,p4,p6,p2) + za(p7,p3)*zb(p2,p3)*zb(p2,p4)*zab2(
     &    p6,p3,p4,p1) )

      Amp(j,1,2,t25)= + iza(p2,p7)*prWs1346*s34**(-1)*facuLl*s134**(-1)
     &  * ( za(p6,k5)*zb(p1,p4)*zab2(p3,p1,p4,p7) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*prWs16*prWs1346*
     & s34**(-1)*facdLl * ( za(p6,p3)*zb(p1,p4)*zab2(k5,p1,p6,p7) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*prWs16*prWs1346*
     & s34**(-1)*facuLl * (  - za(p6,p3)*zb(p1,p4)*zab2(k5,p1,p6,p7) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs1346*
     & s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*zab2(k5,p2,p7,p1)*
     &    zab2(k5,p3,p6,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs1346*
     & s34**(-1)*facuLl*s134**(-1) * ( za(p2,k5)*za(p6,k5)*zb(p1,p4)*
     &    zab2(p3,p1,p4,p2) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & s34**(-1)*facdLl*s2347**(-1) * ( za(p6,k5)*zab2(k5,p2,p7,p4)*
     &    zab3(p3,p2,p4,p7,p1) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zab2(k5,p2,p7,p1)*
     &    zab3(p6,p1,p2,p7,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & prWs1346*s34**(-1)*facdLl * ( za(p2,k5)*za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p1,p6,p2) - za(p6,k5)*zab2(k5,p2,p7,p1)*zab2(p3,p1,p6
     &    ,p4) - za(k5,p3)*zab2(p6,p3,p4,p1)*zab2(k5,p2,p7,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*prWs16*
     & prWs1346*s34**(-1)*facuLl * (  - za(p2,k5)*za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p1,p6,p2) + za(p6,k5)*zab2(k5,p2,p7,p1)*zab2(p3,p1,p6
     &    ,p4) + za(k5,p3)*zab2(p6,p3,p4,p1)*zab2(k5,p2,p7,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(k5,p3
     & )*mt*mton2mwsq*prWs16*prWs1346*facdLl * ( za(p6,p3)*zb(p1,p4)*
     &    zab2(k5,p2,p7,p3) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(k5,p3
     & )*mt*mton2mwsq*prWs16*prWs1346*facuLl * (  - za(p6,p3)*zb(p1,p4)
     &    *zab2(k5,p2,p7,p3) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(k5,p3
     & )*mt**2*prWs16*prts1267*s34**(-1)*facuRl * ( za(p6,p3)*zb(p3,p4)
     &    *zab2(k5,p2,p7,p1) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(k5,p3
     & )*mt**2*prWs16*prWs1346*s34**(-1)*facdLl * (  - za(p6,p3)*zb(p1,
     &    p4)*zab2(k5,p2,p7,p3) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(k5,p3
     & )*mt**2*prWs16*prWs1346*s34**(-1)*facuLl * ( za(p6,p3)*zb(p1,p4)
     &    *zab2(k5,p2,p7,p3) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*prWs16*s34**(-1)*
     & facdLl*s234**(-1)*s2347**(-1) * (  - za(p6,k5)*zb(p2,p4)*zab2(p3
     &    ,p2,p4,p7)*zab4(k5,p2,p3,p4,p7,p1) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*prWs16*prts126*
     & prts1267*s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zab2(p6,p1,
     &    p2,p7)*zab4(k5,p1,p2,p6,p7,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*mt**2*prWs16*prts126
     & *prts1267*s34**(-1)*facuLl * (  - za(p6,k5)*za(k5,p3)*zb(p1,p2)*
     &    zb(p7,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt*
     & mton2mwsq*prWs16*prWs1346*prts12346*facdLl * (  - za(p6,p3)*zb(
     &    p1,p4)*zb(p7,p3)*zab4(k5,p1,p3,p4,p6,p2) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt*
     & mton2mwsq*prWs16*prWs1346*prts12346*facuLl * ( za(p6,p3)*zb(p1,
     &    p4)*zb(p7,p3)*zab4(k5,p1,p3,p4,p6,p2) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facdLl*s346**(-1) * ( za(p6,p3)*zb(
     &    p1,p2)*zb(p7,p3)*zab2(k5,p3,p6,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs1346*prts12346*s34**(-1)*facuLl*s134**(-1) * (  - za(p6,k5)*
     &    zb(p1,p4)*zb(p7,p3)*zab2(p3,p1,p4,p2) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs16*prts12346*s34**(-1)*facdLl*s234**(-1) * ( za(p6,k5)*zb(p2
     &    ,p4)*zb(p7,p3)*zab2(p3,p2,p4,p1) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs16*prts126*prts12346*s34**(-1)*facuRl * (  - za(p6,p3)*zb(p1
     &    ,p2)*zb(p7,p3)*zab4(k5,p1,p2,p3,p6,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs16*prts126*prts12346*s34**(-1)*facuLl * ( za(k5,p3)*zb(p1,p2
     &    )*zb(p7,p3)*zab2(p6,p1,p2,p4) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs16*prts126*prts1267*s34**(-1)*facuRl * ( za(p6,k5)*zb(p1,p2)
     &    *zb(p3,p4)*zab3(p3,p1,p2,p6,p7) + za(k5,p3)*zb(p1,p2)*zb(p3,
     &    p4)*zab2(p6,p1,p2,p7) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs16*prWs1346*prts12346*s34**(-1)*facdLl * ( za(p6,k5)*zb(p1,
     &    p2)*zb(p7,p3)*zab2(p3,p1,p6,p4) - za(p6,p3)*zb(p1,p4)*zb(p7,
     &    p3)*zab2(k5,p1,p6,p2) + za(p6,p3)*zb(p1,p4)*zb(p7,p3)*zab4(k5
     &    ,p1,p3,p4,p6,p2) - za(k5,p3)*zb(p2,p4)*zb(p7,p3)*zab2(p6,p3,
     &    p4,p1) )
      Amp(j,1,2,t25) = Amp(j,1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*
     & prWs16*prWs1346*prts12346*s34**(-1)*facuLl * (  - za(p6,k5)*zb(
     &    p1,p2)*zb(p7,p3)*zab2(p3,p1,p6,p4) + za(p6,p3)*zb(p1,p4)*zb(
     &    p7,p3)*zab2(k5,p1,p6,p2) - za(p6,p3)*zb(p1,p4)*zb(p7,p3)*
     &    zab4(k5,p1,p3,p4,p6,p2) + za(k5,p3)*zb(p2,p4)*zb(p7,p3)*zab2(
     &    p6,p3,p4,p1) )

      Amp(j,2,1,t25)= + iza(k5,p3)*izb(p2,p7)*mt*prWs1346*prts12346*
     & s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*za(p7,p3)*zb(p1,p2)
     &    *zb(p6,p4)*zab3(p6,p1,p3,p4,p2) - za(p6,p3)*za(p7,p3)*zb(p1,
     &    p2)*zb(p3,p4)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt*
     & prWs1346*prts12346*s34**(-1)*facuLl*s134**(-1) * (  - za(p7,p3)*
     &    zb(p1,p4)*zab2(p3,p1,p4,p2)*zab3(p6,p1,p3,p4,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt*prWs16
     & *s34**(-1)*facdLl*s234**(-1)*s2347**(-1) * (  - za(p6,p3)*za(p3,
     &    p4)*zb(p2,p4)**2*zab3(p7,p2,p3,p4,p1) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt*prWs16
     & *prts12346*s34**(-1)*facdLl*s234**(-1) * ( za(p7,p3)*zb(p2,p4)*
     &    zab2(p3,p2,p4,p1)*zab3(p6,p1,p3,p4,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt*prWs16
     & *prts126*prts12346*s34**(-1)*facuLl * (  - za(p7,p3)*zb(p1,p2)*
     &    zab2(p6,p1,p2,p4)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt*prWs16
     & *prWs1346*prts12346*s34**(-1)*facdLl * ( za(p6,p3)*za(p7,p3)*zb(
     &    p1,p2)*zb(p1,p4)*zab3(p1,p3,p4,p6,p2) - za(p6,p3)*za(p7,p3)*
     &    zb(p1,p4)*zb(p2,p6)*zab3(p6,p1,p3,p4,p2) + za(p7,p3)*zb(p1,p2
     &    )*zab2(p3,p1,p6,p4)*zab3(p6,p1,p3,p4,p2) + za(p7,p3)*zb(p2,p4
     &    )*zab2(p6,p3,p4,p1)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt*prWs16
     & *prWs1346*prts12346*s34**(-1)*facuLl * (  - za(p6,p3)*za(p7,p3)*
     &    zb(p1,p2)*zb(p1,p4)*zab3(p1,p3,p4,p6,p2) + za(p6,p3)*za(p7,p3
     &    )*zb(p1,p4)*zb(p2,p6)*zab3(p6,p1,p3,p4,p2) - za(p7,p3)*zb(p1,
     &    p2)*zab2(p3,p1,p6,p4)*zab3(p6,p1,p3,p4,p2) - za(p7,p3)*zb(p2,
     &    p4)*zab2(p6,p3,p4,p1)*zab3(p3,p1,p4,p6,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + iza(k5,p3)*izb(p2,p7)*mt**3*
     & prWs16*prts126*prts12346*s34**(-1)*facuRl * ( za(p6,p3)*za(p7,p3
     &    )*zb(p1,p2)*zb(p2,p4) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mton2mwsq*prWs16*
     & prWs1346*prts12346*facdLl * ( za(p6,p3)*zb(p1,p4)*zb(p2,k5)*
     &    zab4(p7,p1,p3,p4,p6,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mton2mwsq*prWs16*
     & prWs1346*prts12346*facuLl * (  - za(p6,p3)*zb(p1,p4)*zb(p2,k5)*
     &    zab4(p7,p1,p3,p4,p6,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs1346*
     & prts12346*s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*zb(p1,p2)
     &    *zb(p2,k5)*zab2(p7,p3,p6,p4) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs1346*
     & prts12346*s34**(-1)*facuLl*s134**(-1) * ( za(p6,p7)*zb(p1,p4)*
     &    zb(p2,k5)*zab2(p3,p1,p4,p2) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs16*prts12346*
     & s34**(-1)*facdLl*s234**(-1) * (  - za(p6,p7)*zb(p2,k5)*zb(p2,p4)
     &    *zab2(p3,p2,p4,p1) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs16*prts126*
     & prts12346*s34**(-1)*facuRl * ( za(p6,p3)*zb(p1,p2)*zb(p2,k5)*
     &    zab4(p7,p1,p2,p3,p6,p4) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs16*prts126*
     & prts12346*s34**(-1)*facuLl * (  - za(p7,p3)*zb(p1,p2)*zb(p2,k5)*
     &    zab2(p6,p1,p2,p4) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs16*prts126*
     & prts1267*s34**(-1)*facuRl * ( za(p1,p6)*za(p7,p3)*zb(p1,p2)**2*
     &    zb(k5,p4) - za(p6,p7)*zb(p1,p2)*zb(k5,p4)*zab3(p3,p1,p6,p7,p2
     &    ) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs16*prWs1346*
     & prts12346*s34**(-1)*facdLl * (  - za(p6,p7)*zb(p1,p2)*zb(p2,k5)*
     &    zab2(p3,p1,p6,p4) + za(p6,p3)*zb(p1,p4)*zb(p2,k5)*zab2(p7,p1,
     &    p6,p2) - za(p6,p3)*zb(p1,p4)*zb(p2,k5)*zab4(p7,p1,p3,p4,p6,p2
     &    ) + za(p7,p3)*zb(p2,k5)*zb(p2,p4)*zab2(p6,p3,p4,p1) )
      Amp(j,2,1,t25) = Amp(j,2,1,t25) + izb(p2,p7)*mt*prWs16*prWs1346*
     & prts12346*s34**(-1)*facuLl * ( za(p6,p7)*zb(p1,p2)*zb(p2,k5)*
     &    zab2(p3,p1,p6,p4) - za(p6,p3)*zb(p1,p4)*zb(p2,k5)*zab2(p7,p1,
     &    p6,p2) + za(p6,p3)*zb(p1,p4)*zb(p2,k5)*zab4(p7,p1,p3,p4,p6,p2
     &    ) - za(p7,p3)*zb(p2,k5)*zb(p2,p4)*zab2(p6,p3,p4,p1) )

      Amp(j,2,2,t25)= + iza(p2,p7)*iza(p7,k5)*mton2mwsq*prWs16*prWs1346
     & *facdLl * (  - za(p6,p3)*zb(p1,p4)*zab2(k5,p2,p7,k5) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*mton2mwsq
     & *prWs16*prWs1346*facuLl * ( za(p6,p3)*zb(p1,p4)*zab2(k5,p2,p7,k5
     &    ) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*mt*prWs16
     & *prts1267*s34**(-1)*facuRl * (  - za(p6,p3)*zb(k5,p4)*zab2(k5,p2
     &    ,p7,p1) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*mt*prWs16
     & *prWs1346*s34**(-1)*facdLl * ( za(p6,p3)*zb(p1,p4)*zab2(k5,p2,p7
     &    ,k5) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*mt*prWs16
     & *prWs1346*s34**(-1)*facuLl * (  - za(p6,p3)*zb(p1,p4)*zab2(k5,p2
     &    ,p7,k5) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*iza(k5,p3
     & )*mt*prWs1346*s34**(-1)*facdLl*s346**(-1) * ( za(p6,p3)*zab2(k5,
     &    p2,p7,p1)*zab2(p3,p3,p6,p4) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*iza(k5,p3
     & )*mt*prWs1346*s34**(-1)*facuLl*s134**(-1) * (  - za(p2,k5)*za(p6
     &    ,p3)*zb(p1,p4)*zab2(p3,p1,p4,p2) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*iza(k5,p3
     & )*mt*prWs16*s34**(-1)*facdLl*s2347**(-1) * (  - za(p6,p3)*zab2(
     &    k5,p2,p7,p4)*zab3(p3,p2,p4,p7,p1) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*iza(k5,p3
     & )*mt*prWs16*prWs1346*s34**(-1)*facdLl * (  - za(p2,k5)*za(p6,p3)
     &    *zb(p1,p4)*zab2(p3,p1,p6,p2) + za(p6,p3)*zab2(k5,p2,p7,p1)*
     &    zab2(p3,p1,p6,p4) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(p7,k5)*iza(k5,p3
     & )*mt*prWs16*prWs1346*s34**(-1)*facuLl * ( za(p2,k5)*za(p6,p3)*
     &    zb(p1,p4)*zab2(p3,p1,p6,p2) - za(p6,p3)*zab2(k5,p2,p7,p1)*
     &    zab2(p3,p1,p6,p4) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(k5,p3)*mt*
     & prWs1346*s34**(-1)*facuLl*s134**(-1) * (  - za(p6,p3)*zb(p1,p4)*
     &    zab2(p3,p1,p4,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(k5,p3)*mt*prWs16
     & *prWs1346*s34**(-1)*facdLl * (  - za(p6,p3)*zb(p1,p4)*zab2(p3,p1
     &    ,p6,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p2,p7)*iza(k5,p3)*mt*prWs16
     & *prWs1346*s34**(-1)*facuLl * ( za(p6,p3)*zb(p1,p4)*zab2(p3,p1,p6
     &    ,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mton2mwsq*prWs16*
     & prWs1346*prts12346*facdLl * ( za(p6,p3)*zb(p1,p4)*zb(p7,k5)*
     &    zab4(k5,p1,p3,p4,p6,p2) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mton2mwsq*prWs16*
     & prWs1346*prts12346*facuLl * (  - za(p6,p3)*zb(p1,p4)*zb(p7,k5)*
     &    zab4(k5,p1,p3,p4,p6,p2) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs1346*
     & prts12346*s34**(-1)*facdLl*s346**(-1) * (  - za(p6,p3)*zb(p1,p2)
     &    *zb(p6,p4)*zab4(p6,p1,p2,p3,p4,p7) - za(p6,p3)*zb(p1,p2)*zb(
     &    p7,k5)*zab2(k5,p3,p6,p4) - za(p6,p3)*zb(p1,p2)*zb(p3,p4)*
     &    zab4(p3,p1,p2,p4,p6,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs1346*
     & prts12346*s34**(-1)*facuLl*s134**(-1) * ( za(p6,k5)*zb(p1,p4)*
     &    zb(p7,k5)*zab2(p3,p1,p4,p2) - zb(p1,p4)*zab2(p3,p1,p4,p2)*
     &    zab4(p6,p1,p2,p3,p4,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs16*prts12346*
     & s34**(-1)*facdLl*s234**(-1) * (  - za(p6,k5)*zb(p2,p4)*zb(p7,k5)
     &    *zab2(p3,p2,p4,p1) + zb(p2,p4)*zab2(p3,p2,p4,p1)*zab4(p6,p1,
     &    p2,p3,p4,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs16*prts126*
     & prts12346*s34**(-1)*facuRl * ( za(p6,p3)*zb(p1,p2)*zb(p7,k5)*
     &    zab4(k5,p1,p2,p3,p6,p4) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs16*prts126*
     & prts12346*s34**(-1)*facuLl * (  - za(k5,p3)*zb(p1,p2)*zb(p7,k5)*
     &    zab2(p6,p1,p2,p4) - zb(p1,p2)*zab2(p6,p1,p2,p4)*zab4(p3,p1,p2
     &    ,p4,p6,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs16*prts126*
     & prts1267*s34**(-1)*facuRl * (  - za(p6,k5)*zb(p1,p2)*zb(k5,p4)*
     &    zab3(p3,p1,p2,p6,p7) - za(k5,p3)*zb(p1,p2)*zb(k5,p4)*zab2(p6,
     &    p1,p2,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs16*prWs1346*
     & prts12346*s34**(-1)*facdLl * (  - za(p6,k5)*zb(p1,p2)*zb(p7,k5)*
     &    zab2(p3,p1,p6,p4) + za(p6,p3)*zb(p1,p2)*zb(p1,p4)*zab4(p1,p2,
     &    p3,p4,p6,p7) - za(p6,p3)*zb(p1,p4)*zb(p2,p6)*zab4(p6,p1,p2,p3
     &    ,p4,p7) + za(p6,p3)*zb(p1,p4)*zb(p7,k5)*zab2(k5,p1,p6,p2) -
     &    za(p6,p3)*zb(p1,p4)*zb(p7,k5)*zab4(k5,p1,p3,p4,p6,p2) + za(k5
     &    ,p3)*zb(p2,p4)*zb(p7,k5)*zab2(p6,p3,p4,p1) + zb(p1,p2)*zab2(
     &    p3,p1,p6,p4)*zab4(p6,p1,p2,p3,p4,p7) + zb(p2,p4)*zab2(p6,p3,
     &    p4,p1)*zab4(p3,p1,p2,p4,p6,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt*prWs16*prWs1346*
     & prts12346*s34**(-1)*facuLl * ( za(p6,k5)*zb(p1,p2)*zb(p7,k5)*
     &    zab2(p3,p1,p6,p4) - za(p6,p3)*zb(p1,p2)*zb(p1,p4)*zab4(p1,p2,
     &    p3,p4,p6,p7) + za(p6,p3)*zb(p1,p4)*zb(p2,p6)*zab4(p6,p1,p2,p3
     &    ,p4,p7) - za(p6,p3)*zb(p1,p4)*zb(p7,k5)*zab2(k5,p1,p6,p2) +
     &    za(p6,p3)*zb(p1,p4)*zb(p7,k5)*zab4(k5,p1,p3,p4,p6,p2) - za(k5
     &    ,p3)*zb(p2,p4)*zb(p7,k5)*zab2(p6,p3,p4,p1) - zb(p1,p2)*zab2(
     &    p3,p1,p6,p4)*zab4(p6,p1,p2,p3,p4,p7) - zb(p2,p4)*zab2(p6,p3,
     &    p4,p1)*zab4(p3,p1,p2,p4,p6,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt**2*mton2mwsq*
     & prWs16*prWs1346*prts12346*facdLl * (  - za(p6,p3)*zb(p1,p4)*zb(
     &    p2,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt**2*mton2mwsq*
     & prWs16*prWs1346*prts12346*facuLl * ( za(p6,p3)*zb(p1,p4)*zb(p2,
     &    p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt**3*prWs16*prts126
     & *prts12346*s34**(-1)*facuRl * ( za(p6,p3)*zb(p1,p2)*zb(p7,p4) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt**3*prWs16*
     & prWs1346*prts12346*s34**(-1)*facdLl * ( za(p6,p3)*zb(p1,p4)*zb(
     &    p2,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*mt**3*prWs16*
     & prWs1346*prts12346*s34**(-1)*facuLl * (  - za(p6,p3)*zb(p1,p4)*
     &    zb(p2,p7) )
      Amp(j,2,2,t25) = Amp(j,2,2,t25) + iza(p7,k5)*iza(k5,p3)*mt*prWs16
     & *s34**(-1)*facdLl*s234**(-1)*s2347**(-1) * ( za(p6,p3)*zb(p2,p4)
     &    *zab2(p3,p2,p4,p7)*zab4(k5,p2,p3,p4,p7,p1) )


c--- Non-resonant diagrams for left-handed lepton helicity only
      if (j. eq. 1) then
      nr(1,1,t16)= + izb(p1,p7)*izb(p6,p7)*xw**(-1)*iprZ*prWs167*
     & prWs13467 * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p6)*zb(p1,p4)*zab3(p6,p1
     &    ,p4,p7,p2)*s1467**(-1) )
      nr(1,1,t16) = nr(1,1,t16) + izb(p1,p7)*izb(p6,p7)*izb(k5,p3)*mt*
     & mton2mwsq*xw**(-1)*iprZ*prWs167*prWs13467 * ( 1.D0/2.D0*zb(p1,p4
     &    )*zb(p2,p3)*zab2(p3,p6,p7,p1) )
      nr(1,1,t16) = nr(1,1,t16) + izb(p6,p7)*xw**(-1)*iprZ*prWs167*
     & prWs13467 * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p4)*zab3(p7,p1,p4,p6,p2)
     &    *s1467**(-1) )

      nr(1,2,t16)= + iza(p1,p7)*iza(p6,p7)*xw**(-1)*iprZ*prWs167*
     & prWs13467 * ( 1.D0/2.D0*za(k5,p3)*zab2(p6,p1,p7,p4)*zab3(p6,p1,
     &    p4,p7,p2)*s1467**(-1) )
      nr(1,2,t16) = nr(1,2,t16) + iza(p1,p7)*iza(p6,p7)*izb(k5,p3)*mt*
     & mton2mwsq*xw**(-1)*iprZ*prWs167*prWs13467 * ( 1.D0/2.D0*za(p6,p3
     &    )*zb(p2,p3)*zab2(p6,p1,p7,p4) )

      nr(2,1,t16)= + izb(p1,p7)*izb(p6,p7)*mton2mwsq*xw**(-1)*iprZ*
     & prWs167*prWs13467 * (  - 1.D0/2.D0*zb(p1,p4)*zb(p2,k5)*zab2(p3,
     &    p6,p7,p1) )

      nr(2,2,t16)= + iza(p1,p7)*iza(p6,p7)*mton2mwsq*xw**(-1)*iprZ*
     & prWs167*prWs13467 * (  - 1.D0/2.D0*za(p6,p3)*zb(p2,k5)*zab2(p6,
     &    p1,p7,p4) )

      nr(1,1,t25)= + izb(p2,p7)*xw**(-1)*iprZ*prWs16*prWs1346*prts12346
     & *s146**(-1) * (  - 1.D0/2.D0*za(p7,k5)*zb(p1,p4)*zab2(p6,p1,p4,
     &    p2)*zab3(p3,p1,p4,p6,p2) )
      nr(1,1,t25) = nr(1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt*mton2mwsq*
     & xw**(-1)*iprZ*prWs16*prWs1346*prts12346 * ( 1.D0/2.D0*za(p6,p3)*
     &    zb(p1,p4)*zb(p2,p3)*zab4(p7,p1,p3,p4,p6,p2) )
      nr(1,1,t25) = nr(1,1,t25) + izb(p2,p7)*izb(k5,p3)*mt**2*xw**(-1)*
     & iprZ*prWs16*prWs1346*prts12346*s146**(-1) * (  - 1.D0/2.D0*za(p7
     &    ,p3)*zb(p1,p4)*zb(p2,p3)*zab2(p6,p1,p4,p2) )

      nr(1,2,t25)= + iza(p2,p7)*xw**(-1)*iprZ*prWs16*prWs1346*
     & s146**(-1) * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p4)*zab2(p6,p1,p4,p7) )
      nr(1,2,t25) = nr(1,2,t25) + iza(p2,p7)*iza(p7,k5)*xw**(-1)*iprZ*
     & prWs16*prWs1346*s146**(-1) * ( 1.D0/2.D0*za(p2,k5)*za(k5,p3)*zb(
     &    p1,p4)*zab2(p6,p1,p4,p2) )
      nr(1,2,t25) = nr(1,2,t25) + iza(p2,p7)*iza(p7,k5)*izb(k5,p3)*mt*
     & mton2mwsq*xw**(-1)*iprZ*prWs16*prWs1346 * (  - 1.D0/2.D0*za(p6,
     &    p3)*zb(p1,p4)*zab2(k5,p2,p7,p3) )
      nr(1,2,t25) = nr(1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt*mton2mwsq*
     & xw**(-1)*iprZ*prWs16*prWs1346*prts12346 * ( 1.D0/2.D0*za(p6,p3)*
     &    zb(p1,p4)*zb(p7,p3)*zab4(k5,p1,p3,p4,p6,p2) )
      nr(1,2,t25) = nr(1,2,t25) + iza(p7,k5)*izb(k5,p3)*mt**2*xw**(-1)*
     & iprZ*prWs16*prWs1346*prts12346*s146**(-1) * (  - 1.D0/2.D0*za(k5
     &    ,p3)*zb(p1,p4)*zb(p7,p3)*zab2(p6,p1,p4,p2) )

      nr(2,1,t25)= + iza(k5,p3)*izb(p2,p7)*mt*xw**(-1)*iprZ*prWs16*
     & prWs1346*prts12346*s146**(-1) * ( 1.D0/2.D0*za(p7,p3)*zb(p1,p4)*
     &    zab2(p6,p1,p4,p2)*zab3(p3,p1,p4,p6,p2) )
      nr(2,1,t25) = nr(2,1,t25) + izb(p2,p7)*mton2mwsq*xw**(-1)*iprZ*
     & prWs16*prWs1346*prts12346 * (  - 1.D0/2.D0*za(p6,p3)*zb(p1,p4)*
     &    zb(p2,k5)*zab4(p7,p1,p3,p4,p6,p2) )
      nr(2,1,t25) = nr(2,1,t25) + izb(p2,p7)*mt*xw**(-1)*iprZ*prWs16*
     & prWs1346*prts12346*s146**(-1) * ( 1.D0/2.D0*za(p7,p3)*zb(p1,p4)*
     &    zb(p2,k5)*zab2(p6,p1,p4,p2) )

      nr(2,2,t25)= + iza(p2,p7)*iza(p7,k5)*mton2mwsq*xw**(-1)*iprZ*
     & prWs16*prWs1346 * ( 1.D0/2.D0*za(p6,p3)*zb(p1,p4)*zab2(k5,p2,p7,
     &    k5) )
      nr(2,2,t25) = nr(2,2,t25) + iza(p7,k5)*mton2mwsq*xw**(-1)*iprZ*
     & prWs16*prWs1346*prts12346 * (  - 1.D0/2.D0*za(p6,p3)*zb(p1,p4)*
     &    zb(p7,k5)*zab4(k5,p1,p3,p4,p6,p2) )
      nr(2,2,t25) = nr(2,2,t25) + iza(p7,k5)*mt*xw**(-1)*iprZ*prWs16*
     & prWs1346*prts12346*s146**(-1) * ( 1.D0/2.D0*za(k5,p3)*zb(p1,p4)*
     &    zb(p7,k5)*zab2(p6,p1,p4,p2) + 1.D0/2.D0*zb(p1,p4)*zab2(p6,p1,
     &    p4,p2)*zab4(p3,p1,p2,p4,p6,p7) )
      nr(2,2,t25) = nr(2,2,t25) + iza(p7,k5)*mt**2*mton2mwsq*xw**(-1)*
     & iprZ*prWs16*prWs1346*prts12346 * ( 1.D0/2.D0*za(p6,p3)*zb(p1,p4)
     &    *zb(p2,p7) )

      do jt=1,2
      do jg=1,2
      do icol=1,2
      Amp(1,jt,jg,icol)=Amp(1,jt,jg,icol)+nr(jt,jg,icol)
      enddo
      enddo
      enddo
      endif

      enddo

      ampsq16=0d0
      ampsq25=0d0
      do jl=1,2
      do jt=1,2
      do jg=1,2
      ampsq16=ampsq16+abs(Amp(jl,jt,jg,t16))**2
      ampsq25=ampsq25+abs(Amp(jl,jt,jg,t25))**2
      enddo
      enddo
      enddo
      ampsq16=ampsq16*cprop
      ampsq25=ampsq25*cprop

      return
      end
