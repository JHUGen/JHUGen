      function ttbgggmpm(i1,i2,i3,i4,i5,i6,i7)
      implicit none
      include 'types.f'
      complex(dp):: ttbgggmpm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      real(dp):: s129,s1345,s6789,mtsq
      s129=s(i1,i2)+s(i1,i3)+s(i2,i3)
      s1345=s(i1,i4)+s(i1,i5)
      s6789=s(i3,i6)+s(i3,i7)
      mtsq=mt**2
      ttbgggmpm =  + s129**(-1) * (  - 1/(za(i1,i2))/(za(i2,i3))/(zb(i1
     &    ,i2))*za(i1,i3)**2*za(i5,i3)*za(i6,i7)*zb(i2,i6)*zb(i4,i5) -
     &    1/(za(i1,i2))/(za(i2,i3))/(zb(i2,i3))*za(i1,i3)**2*za(i1,i5)*
     &    za(i6,i7)*zb(i2,i6)*zb(i4,i5) )
      ttbgggmpm = ttbgggmpm + s1345**(-1)*s6789**(-1) * (  - 1/(za(i1,
     &    i2))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i4)*za(i1,i5)*za(i6,i3)*
     &    za(i6,i7)*zb(i2,i4)*zb(i2,i6)**2*zb(i4,i5) - 1/(za(i1,i2))/(
     &    zb(i1,i2))/(zb(i2,i3))*za(i1,i4)*za(i1,i5)*za(i6,i7)*za(i7,i3
     &    )*zb(i2,i4)*zb(i2,i6)*zb(i2,i7)*zb(i4,i5) - 1/(za(i1,i2))/(
     &    zb(i1,i2))/(zb(i2,i3))*za(i1,i5)**2*za(i6,i3)*za(i6,i7)*zb(i2
     &    ,i5)*zb(i2,i6)**2*zb(i4,i5) - 1/(za(i1,i2))/(zb(i1,i2))/(zb(
     &    i2,i3))*za(i1,i5)**2*za(i6,i7)*za(i7,i3)*zb(i2,i5)*zb(i2,i6)*
     &    zb(i2,i7)*zb(i4,i5) )
      ttbgggmpm = ttbgggmpm + s1345**(-1) * ( 1/(za(i1,i2))/(za(i2,i3))
     &    /(zb(i1,i2))/(zb(i2,i3))*za(i1,i3)*za(i1,i5)*za(i4,i3)*za(i6,
     &    i7)*zb(i2,i4)*zb(i2,i6)*zb(i4,i5) + 1/(za(i1,i2))/(za(i2,i3))
     &    /(zb(i1,i2))/(zb(i2,i3))*za(i1,i3)*za(i1,i5)*za(i5,i3)*za(i6,
     &    i7)*zb(i2,i5)*zb(i2,i6)*zb(i4,i5) - 1/(za(i1,i2))/(za(i2,i3))
     &    /(zb(i2,i3))*za(i1,i3)**2*za(i1,i5)*za(i6,i7)*zb(i2,i6)*zb(i4
     &    ,i5) )
      ttbgggmpm = ttbgggmpm + mtsq*s129**(-1) * ( 1/(za(i1,i2))/(za(i2,
     &    i3))/(zb(i1,i2))*za(i1,i3)**2*za(i7,i3)*zb(i2,i4) + 1/(za(i1,
     &    i2))/(za(i2,i3))/(zb(i2,i3))*za(i1,i3)**2*za(i1,i7)*zb(i2,i4)
     &     )
      ttbgggmpm = ttbgggmpm + mtsq*s1345**(-1)*s6789**(-1) * ( 1/(za(i1
     &    ,i2))/(zb(i1,i2))*za(i1,i3)*za(i1,i4)*za(i7,i3)*zb(i2,i4)**2
     &     + 1/(za(i1,i2))/(zb(i1,i2))*za(i1,i3)*za(i1,i5)*za(i7,i3)*
     &    zb(i2,i4)*zb(i2,i5) + 1/(za(i1,i2))/(zb(i1,i2))/(zb(i2,i3))*
     &    za(i1,i3)*za(i1,i4)*za(i6,i7)*zb(i2,i4)**2*zb(i2,i6) + 1/(za(
     &    i1,i2))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i3)*za(i1,i5)*za(i6,i7)
     &    *zb(i2,i4)*zb(i2,i5)*zb(i2,i6) + 1/(za(i1,i2))/(zb(i1,i2))/(
     &    zb(i2,i3))*za(i1,i4)*za(i1,i6)*za(i7,i3)*zb(i2,i4)**2*zb(i2,
     &    i6) + 1/(za(i1,i2))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i4)*za(i1,
     &    i7)*za(i7,i3)*zb(i2,i4)**2*zb(i2,i7) + 1/(za(i1,i2))/(zb(i1,
     &    i2))/(zb(i2,i3))*za(i1,i5)*za(i1,i6)*za(i7,i3)*zb(i2,i4)*zb(
     &    i2,i5)*zb(i2,i6) + 1/(za(i1,i2))/(zb(i1,i2))/(zb(i2,i3))*za(
     &    i1,i5)*za(i1,i7)*za(i7,i3)*zb(i2,i4)*zb(i2,i5)*zb(i2,i7) )
      ttbgggmpm = ttbgggmpm + mtsq*s1345**(-1) * ( 1/(za(i1,i2))/(za(i2
     &    ,i3))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i3)**2*za(i6,i7)*zb(i2,i4
     &    )*zb(i2,i6) - 1/(za(i1,i2))/(za(i2,i3))/(zb(i1,i2))/(zb(i2,i3
     &    ))*za(i1,i3)*za(i1,i4)*za(i7,i3)*zb(i2,i4)**2 - 1/(za(i1,i2))
     &    /(za(i2,i3))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i3)*za(i1,i5)*za(
     &    i7,i3)*zb(i2,i4)*zb(i2,i5) )

      return
      end
