      subroutine box5(k1,k2,k3,k4,k5,k6,za,zb,app,apm,amp,amm)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      real(dp):: mtsq,s134,rat,mtsqps134
      integer:: k1,k2,k3,k4,k5,k6,i1,i2
      complex(dp):: app,apm,amp,amm,zab2,izab2,iza,izb
c--- statement functions
      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      izab2(k1,k2,k3,k4)=cone/(za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4))
c--- end statement functions

      mtsq=mt**2
      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      rat=mtsq/s(k1,k2)
      mtsqps134=mtsq+s134
      app =  + iza(k1,k2)*mtsq * (  - 2._dp*za(k1,k3)*za(k1,k5)*zb(k1,k2
     &    )*zb(k1,k4)*zb(k1,k6) + 2._dp*za(k1,k3)*za(k2,k5)*zb(k1,k2)*
     &    zb(k1,k4)*zb(k2,k6)*rat + 2._dp*za(k1,k5)*za(k2,k3)*zb(k1,k2)*
     &    zb(k1,k6)*zb(k2,k4)*rat + 2._dp*za(k1,k5)*za(k3,k4)*zb(k1,k2)*
     &    zb(k1,k4)*zb(k4,k6) )
      app = app + iza(k1,k2)*izab2(k1,k3,k4,k2)*mtsq*mtsqps134 * ( za(
     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) + za(k1,k3)*
     &    za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4) - za(k1,k5)*za(k3,k4)
     &    *zb(k1,k2)*zb(k2,k4)*zb(k4,k6) )
      app = app + iza(k1,k2)*izab2(k1,k3,k4,k2)**2*mtsq*mtsqps134**2
     &  * (  - za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6) )
      app = app + iza(k1,k2)*izab2(k1,k3,k4,k2)*zab2(k2,k3,k4,k1)*mtsq
     &  * (  - 2._dp*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*
     &    rat )
      app = app + iza(k1,k2)*izab2(k2,k3,k4,k1)*mtsq*mtsqps134 * ( za(
     &    k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k1,k5)*
     &    za(k2,k3)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k2,k5)*za(k3,k4)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      app = app + iza(k1,k2)*izab2(k2,k3,k4,k1)**2*mtsq*mtsqps134**2
     &  * (  - za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) )
      app = app + iza(k1,k2)*izab2(k2,k3,k4,k1)*zab2(k1,k3,k4,k2)*mtsq
     &  * (  - 2._dp*za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*
     &    rat )
      amm =  + izb(k1,k2)*mtsq * ( 2._dp*za(k1,k2)*za(k1,k3)*za(k2,k5)*
     &    zb(k1,k4)*zb(k2,k6)*rat + 2._dp*za(k1,k2)*za(k1,k3)*za(k2,k5)*
     &    zb(k1,k6)*zb(k2,k4) + 2._dp*za(k1,k2)*za(k1,k5)*za(k2,k3)*zb(
     &    k1,k6)*zb(k2,k4)*rat - 2._dp*za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(
     &    k2,k4)*zb(k4,k6) )
      amm = amm + izb(k1,k2)*izab2(k1,k3,k4,k2)*mtsq*mtsqps134 * ( za(
     &    k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4) - za(k1,k2)*
     &    za(k1,k3)*za(k2,k5)*zb(k2,k4)*zb(k2,k6) - za(k1,k2)*za(k1,k5)
     &    *za(k3,k4)*zb(k2,k4)*zb(k4,k6) )
      amm = amm + izb(k1,k2)*izab2(k1,k3,k4,k2)**2*mtsq*mtsqps134**2
     &  * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6) )
      amm = amm + izb(k1,k2)*izab2(k1,k3,k4,k2)*zab2(k2,k3,k4,k1)*mtsq
     &  * (  - 2._dp*za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)*
     &    rat )
      amm = amm + izb(k1,k2)*izab2(k2,k3,k4,k1)*mtsq*mtsqps134 * ( za(
     &    k1,k2)*za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*
     &    za(k2,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,k4) - za(k1,k2)*za(k2,k5)
     &    *za(k3,k4)*zb(k1,k4)*zb(k4,k6) )
      amm = amm + izb(k1,k2)*izab2(k2,k3,k4,k1)**2*mtsq*mtsqps134**2
     &  * (  - za(k1,k2)*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6) )
      amm = amm + izb(k1,k2)*izab2(k2,k3,k4,k1)*zab2(k1,k3,k4,k2)*mtsq
     &  * (  - 2._dp*za(k1,k2)*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*
     &    rat )
      apm =  + mtsq * (  - 2._dp*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)
     & - 2._dp*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*rat )
      apm = apm + izab2(k1,k3,k4,k2)*mtsqps134 * ( za(k1,k2)*za(k1,k3)*
     &    za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*za(k2,k5)
     &    *za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)*mtsq*mtsqps134 * (  - za(k1,k3)*
     &    za(k2,k5)*zb(k1,k4)*zb(k1,k6) - za(k1,k5)*za(k2,k3)*zb(k1,k4)
     &    *zb(k1,k6) + za(k2,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,k4) + za(k2,
     &    k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**2*mtsqps134**2 * ( za(k1,k2)*za(
     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*
     &    za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) - za(k1,k2)
     &    *za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4) - za(k1,k2
     &    )*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6) + za(k1,
     &    k2)*za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**2*mtsq*mtsqps134**2 * ( za(k1,k3)
     &    *za(k2,k5)*zb(k1,k4)*zb(k2,k6) + za(k1,k5)*za(k2,k3)*zb(k1,k6
     &    )*zb(k2,k4) )
      apm = apm + izab2(k1,k3,k4,k2)**3*mtsqps134**3 * (  - za(k1,k2)*
     &    za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) - za(k1,k2)
     &    *za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4) + za(k1,k2
     &    )*za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6) + za(k1,
     &    k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**4*mtsqps134**4 * ( za(k1,k2)*za(
     &    k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**3*zab2(k2,k3,k4,k1)*mtsq*
     & mtsqps134**2 * (  - 4._dp*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)
     & )
      apm = apm + izab2(k1,k3,k4,k2)**2*zab2(k2,k3,k4,k1)*mtsq*
     & mtsqps134 * ( 3._dp*za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k2,k6) + 3._dp
     &    *za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4) - 3._dp*za(k1,k3)*za(
     &    k2,k5)*zb(k2,k4)*zb(k2,k6) - 3._dp*za(k1,k5)*za(k3,k4)*zb(k2,
     &    k4)*zb(k4,k6) )
      apm = apm + izab2(k1,k3,k4,k2)**2*zab2(k2,k3,k4,k1)**2*mtsq * ( 
     & - 2._dp*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)*rat )
      apm = apm + izab2(k1,k3,k4,k2)*zab2(k2,k3,k4,k1)*mtsq * (  - 2._dp
     &    *za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6) + 2._dp*za(k1,k3)*za(
     &    k2,k5)*zb(k1,k4)*zb(k2,k6) + 2._dp*za(k1,k3)*za(k2,k5)*zb(k1,
     &    k4)*zb(k2,k6)*rat + 2._dp*za(k1,k3)*za(k2,k5)*zb(k1,k6)*zb(k2,
     &    k4) + 2._dp*za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*rat + 2._dp
     &    *za(k1,k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6) - 2._dp*za(k2,k5)*za(
     &    k3,k4)*zb(k2,k4)*zb(k4,k6) )
      amp =  + mtsq * (  - 2._dp*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6)
     &    *rat )
      amp = amp + izab2(k2,k3,k4,k1)*mtsq*mtsqps134 * (  - za(k1,k3)*
     &    za(k1,k5)*zb(k1,k6)*zb(k2,k4) + za(k1,k5)*za(k3,k4)*zb(k2,k4)
     &    *zb(k4,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**2*mtsq*mtsqps134**2 * ( za(k1,k3)
     &    *za(k2,k5)*zb(k1,k4)*zb(k2,k6) + za(k1,k5)*za(k2,k3)*zb(k1,k6
     &    )*zb(k2,k4) )
      amp = amp + izab2(k2,k3,k4,k1)**3*mtsqps134**3 * (  - za(k1,k2)*
     &    za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k1,k2)
     &    *za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**4*mtsqps134**4 * ( za(k1,k2)*za(
     &    k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**3*zab2(k1,k3,k4,k2)*mtsq*
     & mtsqps134**2 * (  - 4._dp*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6))
      amp = amp + izab2(k2,k3,k4,k1)**2*zab2(k1,k3,k4,k2)*mtsq*
     & mtsqps134 * ( 3._dp*za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6) - 3._dp
     &    *za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6) )
      amp = amp + izab2(k2,k3,k4,k1)**2*zab2(k1,k3,k4,k2)**2*mtsq * ( 
     &     - 2._dp*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*rat )
      amp = amp + izab2(k2,k3,k4,k1)*zab2(k1,k3,k4,k2)*mtsq * ( 2._dp*
     &    za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)*rat + 2._dp*za(k1,k5)*
     &    za(k2,k3)*zb(k1,k6)*zb(k2,k4)*rat )
      app=app*(-half*im)/(s(k3,k4)*s(k5,k6))
      amm=amm*(-half*im)/(s(k3,k4)*s(k5,k6))
      apm=apm*(-half*im)/(s(k3,k4)*s(k5,k6))
      amp=amp*(-half*im)/(s(k3,k4)*s(k5,k6))
      return
      end
