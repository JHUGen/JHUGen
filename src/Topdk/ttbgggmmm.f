      double complex function ttbgggmmm(i1,i2,i3,i4,i5,i6,i7)
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
      ttbgggmmm =  + mtsq*s129**(-1) * ( 1/(zb(i1,i2))/(zb(i1,i6))*za(
     &    i2,i3)*za(i7,i3)*zb(i4,i6) + 1/(zb(i1,i2))/(zb(i1,i6))/(zb(i6
     &    ,i3))*za(i2,i3)*za(i2,i7)*zb(i2,i6)*zb(i4,i6) + 1/(zb(i1,i2))
     &    /(zb(i2,i6))*za(i1,i3)*za(i7,i3)*zb(i4,i6) + 1/(zb(i1,i2))/(
     &    zb(i2,i6))/(zb(i6,i3))*za(i1,i3)*za(i1,i7)*zb(i1,i6)*zb(i4,i6
     &    ) + 1/(zb(i1,i2))/(zb(i6,i3))*za(i1,i3)*za(i2,i7)*zb(i4,i6)
     &     + 1/(zb(i1,i2))/(zb(i6,i3))*za(i1,i7)*za(i2,i3)*zb(i4,i6) + 
     &    1/(zb(i1,i6))/(zb(i2,i3))*za(i1,i2)*za(i7,i3)*zb(i4,i6) - 1/(
     &    zb(i1,i6))/(zb(i2,i3))*za(i1,i3)*za(i2,i7)*zb(i4,i6) - 1/(zb(
     &    i1,i6))/(zb(i2,i3))/(zb(i2,i6))*za(i1,i3)*za(i7,i3)*zb(i4,i6)
     &    *zb(i6,i3) + 1/(zb(i1,i6))/(zb(i2,i3))/(zb(i6,i3))*za(i1,i2)*
     &    za(i2,i7)*zb(i2,i6)*zb(i4,i6) - 1/(zb(i2,i3))/(zb(i2,i6))*za(
     &    i1,i3)*za(i1,i7)*zb(i4,i6) + 1/(zb(i2,i3))/(zb(i6,i3))*za(i1,
     &    i2)*za(i1,i7)*zb(i4,i6) )
      ttbgggmmm = ttbgggmmm + mtsq*s1345**(-1)*s6789**(-1) * ( 1/(zb(i1
     &    ,i6))/(zb(i2,i6))*za(i1,i4)*za(i2,i3)*za(i7,i3)*zb(i4,i6)**2
     &     + 1/(zb(i1,i6))/(zb(i2,i6))*za(i1,i5)*za(i2,i3)*za(i7,i3)*
     &    zb(i4,i6)*zb(i5,i6) + 1/(zb(i1,i6))/(zb(i2,i6))/(zb(i6,i3))*
     &    za(i1,i4)*za(i2,i7)*za(i7,i3)*zb(i4,i6)**2*zb(i6,i7) + 1/(zb(
     &    i1,i6))/(zb(i2,i6))/(zb(i6,i3))*za(i1,i5)*za(i2,i7)*za(i7,i3)
     &    *zb(i4,i6)*zb(i5,i6)*zb(i6,i7) )
      ttbgggmmm = ttbgggmmm + mtsq*s1345**(-1) * (  - 1/(zb(i1,i6))/(
     &    zb(i2,i3))/(zb(i2,i6))*za(i1,i4)*za(i7,i3)*zb(i4,i6)**2 - 1/(
     &    zb(i1,i6))/(zb(i2,i3))/(zb(i2,i6))*za(i1,i5)*za(i7,i3)*zb(i4,
     &    i6)*zb(i5,i6) - 1/(zb(i1,i6))/(zb(i2,i3))/(zb(i6,i3))*za(i1,
     &    i4)*za(i2,i7)*zb(i4,i6)**2 - 1/(zb(i1,i6))/(zb(i2,i3))/(zb(i6
     &    ,i3))*za(i1,i5)*za(i2,i7)*zb(i4,i6)*zb(i5,i6) )
      ttbgggmmm = ttbgggmmm + mtsq*s6789**(-1) * ( 1/(zb(i1,i2))/(zb(i1
     &    ,i6))*za(i2,i3)*za(i7,i3)*zb(i4,i6) + 1/(zb(i1,i2))/(zb(i1,i6
     &    ))/(zb(i6,i3))*za(i2,i7)*za(i7,i3)*zb(i4,i6)*zb(i6,i7) + 1/(
     &    zb(i1,i2))/(zb(i2,i6))*za(i1,i3)*za(i7,i3)*zb(i4,i6) + 1/(zb(
     &    i1,i2))/(zb(i2,i6))/(zb(i6,i3))*za(i1,i7)*za(i7,i3)*zb(i4,i6)
     &    *zb(i6,i7) )

      return
      end
