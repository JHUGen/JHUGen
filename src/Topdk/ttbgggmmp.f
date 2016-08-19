      double complex function ttbgggmmp(i1,i2,i3,i4,i5,i6,i7)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer i1,i2,i3,i4,i5,i6,i7
      double precision s129,s1345,s6789,mtsq
      s129=s(i1,i2)+s(i1,i3)+s(i2,i3)
      s1345=s(i1,i4)+s(i1,i5)
      s6789=s(i3,i6)+s(i3,i7)
      mtsq=mt**2
      ttbgggmmp =  + s129**(-1) * ( 1/(za(i2,i3))/(za(i6,i3))/(zb(i1,i3
     &    ))*za(i1,i2)*za(i2,i5)*za(i2,i6)*za(i6,i7)*zb(i4,i5)*zb(i6,i3
     &    ) + 1/(za(i2,i3))/(za(i6,i3))/(zb(i2,i3))*za(i1,i2)*za(i1,i5)
     &    *za(i2,i6)*za(i6,i7)*zb(i4,i5)*zb(i6,i3) + 1/(za(i6,i3))/(zb(
     &    i1,i2))*za(i1,i5)*za(i2,i6)*za(i6,i7)*zb(i4,i5)*zb(i6,i3) + 
     &    1/(za(i6,i3))/(zb(i1,i2))*za(i1,i6)*za(i2,i5)*za(i6,i7)*zb(i4
     &    ,i5)*zb(i6,i3) + 1/(za(i6,i3))/(zb(i1,i2))/(zb(i1,i3))*za(i2,
     &    i5)*za(i2,i6)*za(i6,i7)*zb(i2,i3)*zb(i4,i5)*zb(i6,i3) + 1/(
     &    za(i6,i3))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i5)*za(i1,i6)*za(i6,
     &    i7)*zb(i1,i3)*zb(i4,i5)*zb(i6,i3) )
      ttbgggmmp = ttbgggmmp + s1345**(-1)*s6789**(-1) * ( 1/(za(i6,i3))
     &    /(zb(i1,i3))/(zb(i2,i3))*za(i1,i5)*za(i2,i4)*za(i6,i7)**2*zb(
     &    i4,i3)*zb(i4,i5)*zb(i6,i3)*zb(i7,i3) + 1/(za(i6,i3))/(zb(i1,
     &    i3))/(zb(i2,i3))*za(i1,i5)*za(i2,i5)*za(i6,i7)**2*zb(i4,i5)*
     &    zb(i5,i3)*zb(i6,i3)*zb(i7,i3) - 1/(za(i6,i3))/(zb(i2,i3))*za(
     &    i1,i2)*za(i1,i5)*za(i6,i7)**2*zb(i4,i5)*zb(i6,i3)*zb(i7,i3) )
      ttbgggmmp = ttbgggmmp + s1345**(-1) * (  - 1/(za(i2,i3))/(za(i6,
     &    i3))/(zb(i1,i3))/(zb(i2,i3))*za(i1,i5)*za(i2,i4)*za(i2,i6)*
     &    za(i6,i7)*zb(i4,i3)*zb(i4,i5)*zb(i6,i3) - 1/(za(i2,i3))/(za(
     &    i6,i3))/(zb(i1,i3))/(zb(i2,i3))*za(i1,i5)*za(i2,i5)*za(i2,i6)
     &    *za(i6,i7)*zb(i4,i5)*zb(i5,i3)*zb(i6,i3) + 1/(za(i2,i3))/(za(
     &    i6,i3))/(zb(i2,i3))*za(i1,i2)*za(i1,i5)*za(i2,i6)*za(i6,i7)*
     &    zb(i4,i5)*zb(i6,i3) )
      ttbgggmmp = ttbgggmmp + s6789**(-1) * ( 1/(za(i6,i3))/(zb(i1,i2))
     &    /(zb(i1,i3))*za(i2,i5)*za(i6,i7)**2*zb(i4,i5)*zb(i6,i3)*zb(i7
     &    ,i3) + 1/(za(i6,i3))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i5)*za(i6,
     &    i7)**2*zb(i4,i5)*zb(i6,i3)*zb(i7,i3) )
      ttbgggmmp = ttbgggmmp + mtsq*s129**(-1) * (  - 1/(za(i2,i3))/(za(
     &    i6,i3))/(zb(i1,i3))*za(i1,i2)*za(i2,i6)*za(i2,i7)*zb(i4,i3)
     &     - 1/(za(i2,i3))/(za(i6,i3))/(zb(i2,i3))*za(i1,i2)*za(i1,i7)*
     &    za(i2,i6)*zb(i4,i3) - 1/(za(i6,i3))/(zb(i1,i2))*za(i1,i6)*za(
     &    i2,i7)*zb(i4,i3) - 1/(za(i6,i3))/(zb(i1,i2))*za(i1,i7)*za(i2,
     &    i6)*zb(i4,i3) - 1/(za(i6,i3))/(zb(i1,i2))/(zb(i1,i3))*za(i2,
     &    i6)*za(i2,i7)*zb(i2,i3)*zb(i4,i3) - 1/(za(i6,i3))/(zb(i1,i2))
     &    /(zb(i2,i3))*za(i1,i6)*za(i1,i7)*zb(i1,i3)*zb(i4,i3) )
      ttbgggmmp = ttbgggmmp + mtsq*s1345**(-1)*s6789**(-1) * (  - 1/(
     &    za(i6,i3))/(zb(i1,i3))/(zb(i2,i3))*za(i1,i2)*za(i6,i7)**2*zb(
     &    i4,i3)*zb(i6,i3)*zb(i7,i3) - 1/(za(i6,i3))/(zb(i1,i3))/(zb(i2
     &    ,i3))*za(i1,i4)*za(i2,i7)*za(i6,i7)*zb(i4,i3)**2*zb(i7,i3) - 
     &    1/(za(i6,i3))/(zb(i1,i3))/(zb(i2,i3))*za(i1,i5)*za(i2,i7)*za(
     &    i6,i7)*zb(i4,i3)*zb(i5,i3)*zb(i7,i3) )
      ttbgggmmp = ttbgggmmp + mtsq*s1345**(-1) * ( 1/(za(i2,i3))/(za(i6
     &    ,i3))/(zb(i1,i3))/(zb(i2,i3))*za(i1,i2)*za(i2,i6)*za(i6,i7)*
     &    zb(i4,i3)*zb(i6,i3) + 1/(za(i2,i3))/(za(i6,i3))/(zb(i1,i3))/(
     &    zb(i2,i3))*za(i1,i4)*za(i2,i6)*za(i2,i7)*zb(i4,i3)**2 + 1/(
     &    za(i2,i3))/(za(i6,i3))/(zb(i1,i3))/(zb(i2,i3))*za(i1,i5)*za(
     &    i2,i6)*za(i2,i7)*zb(i4,i3)*zb(i5,i3) )
      ttbgggmmp = ttbgggmmp + mtsq*s6789**(-1) * (  - 1/(za(i6,i3))/(
     &    zb(i1,i2))/(zb(i1,i3))*za(i2,i7)*za(i6,i7)*zb(i4,i3)*zb(i7,i3
     &    ) - 1/(za(i6,i3))/(zb(i1,i2))/(zb(i2,i3))*za(i1,i7)*za(i6,i7)
     &    *zb(i4,i3)*zb(i7,i3) )

      return
      end
