      double complex function ttbqqbsqpp(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer k1,k2,k3,k4,k5,k6,k7
      double precision s129,s3459,mtsq
      s129=s(k1,k2)+s(k1,k3)+s(k2,k3)
      s3459=s(k4,k3)+s(k5,k3)
      mtsq=mt**2
      ttbqqbsqpp =  + s129**(-1) * ( 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,
     &    k2))*za(k1,k5)**2*za(k6,k7)*zb(k1,k3)*zb(k2,k6)*zb(k4,k5) - 
     &    1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)**2*za(k6,k7)*
     &    zb(k1,k6)*zb(k2,k3)*zb(k4,k5) + 1/(za(k1,k3))/(za(k5,k3))*za(
     &    k1,k5)**2*za(k6,k7)*zb(k2,k6)*zb(k4,k5) )
      ttbqqbsqpp = ttbqqbsqpp + mtsq*s129**(-1) * (  - 1/(za(k1,k2))/(
     &    za(k5,k3))/(zb(k1,k2))*za(k1,k3)*za(k5,k7)*zb(k2,k3)*zb(k4,k3
     &    ) - 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k1,k7)
     &    *zb(k1,k3)*zb(k2,k4) + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*
     &    za(k1,k5)*za(k1,k7)*zb(k1,k4)*zb(k2,k3) + 1/(za(k1,k2))/(za(
     &    k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k2,k7)*zb(k2,k3)*zb(k2,k4)
     &     - 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k7)*za(k2,k5)*
     &    zb(k2,k3)*zb(k2,k4) - 1/(za(k1,k3))/(za(k5,k3))*za(k1,k5)*za(
     &    k1,k7)*zb(k2,k4) )
      ttbqqbsqpp = ttbqqbsqpp + mtsq*s3459**(-1) * ( 1/(za(k1,k2))/(za(
     &    k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k6,k7)*zb(k2,k6)*zb(k4,k3)
     &     + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k7)*za(k4,k5)*
     &    zb(k2,k4)*zb(k4,k3) - 1/(za(k1,k2))/(zb(k1,k2))*za(k1,k7)*zb(
     &    k2,k3)*zb(k4,k3) )

      return
      end
