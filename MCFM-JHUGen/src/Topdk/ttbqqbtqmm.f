      function ttbqqbtqmm(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'types.f'
      complex(dp):: ttbqqbtqmm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      real(dp):: s129,s6789,mtsq
      s129=s(k1,k2)+s(k1,k3)+s(k2,k3)
      s6789=s(k3,k6)+s(k3,k7)
      mtsq=mt**2
      ttbqqbtqmm =  + s129**(-1) * ( 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,
     &    k3))*za(k1,k3)*za(k2,k5)*za(k6,k7)*zb(k1,k5)*zb(k1,k6)*zb(k4,
     &    k5) - 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k1,k5)*za(k2,
     &    k3)*za(k6,k7)*zb(k1,k5)*zb(k1,k6)*zb(k4,k5) - 1/(za(k1,k2))/(
     &    zb(k1,k2))/(zb(k5,k3))*za(k2,k3)*za(k2,k5)*za(k6,k7)*zb(k1,k5
     &    )*zb(k2,k6)*zb(k4,k5) + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))
     &    *za(k2,k3)*za(k2,k5)*za(k6,k7)*zb(k1,k6)*zb(k2,k5)*zb(k4,k5)
     &     + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k2,k3)*za(k5,k3)*
     &    za(k6,k7)*zb(k1,k3)*zb(k4,k5)*zb(k5,k6) + 1/(zb(k2,k3))*za(k5
     &    ,k3)*za(k6,k7)*zb(k1,k6)*zb(k4,k5) + 1/(zb(k2,k3))/(zb(k5,k3)
     &    )*za(k2,k5)*za(k6,k7)*zb(k1,k6)*zb(k2,k5)*zb(k4,k5) )
      ttbqqbtqmm = ttbqqbtqmm + s6789**(-1) * ( 1/(za(k1,k2))/(zb(k1,k2
     &    ))/(zb(k5,k3))*za(k2,k5)*za(k6,k3)*za(k6,k7)*zb(k1,k6)*zb(k4,
     &    k5)*zb(k5,k6) + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k2,
     &    k5)*za(k6,k7)*za(k7,k3)*zb(k1,k7)*zb(k4,k5)*zb(k5,k6) )
      ttbqqbtqmm = ttbqqbtqmm + mtsq*s129**(-1) * (  - 1/(za(k1,k2))/(
     &    zb(k1,k2))/(zb(k5,k3))*za(k1,k3)*za(k2,k7)*zb(k1,k4)*zb(k1,k5
     &    ) + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k1,k7)*za(k2,k3)
     &    *zb(k1,k4)*zb(k1,k5) - 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*
     &    za(k2,k3)*za(k2,k7)*zb(k1,k4)*zb(k2,k5) + 1/(za(k1,k2))/(zb(
     &    k1,k2))/(zb(k5,k3))*za(k2,k3)*za(k2,k7)*zb(k1,k5)*zb(k2,k4)
     &     + 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k2,k3)*za(k7,k3)*
     &    zb(k1,k3)*zb(k4,k5) - 1/(zb(k2,k3))*za(k7,k3)*zb(k1,k4) - 1/(
     &    zb(k2,k3))/(zb(k5,k3))*za(k2,k7)*zb(k1,k4)*zb(k2,k5) )
      ttbqqbtqmm = ttbqqbtqmm + mtsq*s6789**(-1) * (  - 1/(za(k1,k2))/(
     &    zb(k1,k2))*za(k2,k3)*za(k7,k3)*zb(k1,k4) - 1/(za(k1,k2))/(zb(
     &    k1,k2))/(zb(k5,k3))*za(k2,k3)*za(k6,k7)*zb(k1,k4)*zb(k5,k6)
     &     - 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*za(k2,k5)*za(k7,k3)*
     &    zb(k1,k5)*zb(k4,k5) - 1/(za(k1,k2))/(zb(k1,k2))/(zb(k5,k3))*
     &    za(k2,k6)*za(k7,k3)*zb(k1,k4)*zb(k5,k6) - 1/(za(k1,k2))/(zb(
     &    k1,k2))/(zb(k5,k3))*za(k2,k7)*za(k7,k3)*zb(k1,k4)*zb(k5,k7) )

      return
      end
