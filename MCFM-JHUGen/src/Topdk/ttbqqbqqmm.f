      double complex function ttbqqbqqmm(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer k1,k2,k3,k4,k5,k6,k7
      double precision s129,mtsq
      s129=s(k1,k2)+s(k1,k3)+s(k2,k3)
      mtsq=mt**2
      ttbqqbqqmm =  + xn**(-1)*s129**(-1) * ( 1/(zb(k1,k3))/(zb(k5,k3))
     &    *za(k2,k5)*za(k6,k7)*zb(k1,k5)*zb(k1,k6)*zb(k4,k5) - 1/(zb(k2
     &    ,k3))*za(k5,k3)*za(k6,k7)*zb(k1,k6)*zb(k4,k5) - 1/(zb(k2,k3))
     &    /(zb(k5,k3))*za(k2,k5)*za(k6,k7)*zb(k1,k6)*zb(k2,k5)*zb(k4,k5
     &    ) )
      ttbqqbqqmm = ttbqqbqqmm + xn**(-1)*mtsq*s129**(-1) * (  - 1/(zb(
     &    k1,k3))/(zb(k5,k3))*za(k2,k7)*zb(k1,k4)*zb(k1,k5) + 1/(zb(k2,
     &    k3))*za(k7,k3)*zb(k1,k4) + 1/(zb(k2,k3))/(zb(k5,k3))*za(k2,k7
     &    )*zb(k1,k4)*zb(k2,k5) )

      return
      end
