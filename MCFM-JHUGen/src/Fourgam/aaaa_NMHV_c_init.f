!===== T. Dennen, May 2014
!===== Three-mass triangle coefficients for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,-)ga(i5,+)ga(i6,+)
      subroutine aaaa_NMHV_c_init(i1,i2,i3,i4,i5,i6,za,zb,aaaa_NMHV_c)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::i1,i2,i3,i4,i5,i6,i7,i8
      complex(dp)::aaaa_NMHV_c(90), zab2, zaa2, zbb2, zab3, t
      complex(dp)::sqrtd
      complex(dp)::sqrtd1245, sqrtd1246, sqrtd2345, sqrtd2346 
      complex(dp)::sqrtd2435, sqrtd2436

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zaa2(i1,i2,i3,i4,i5,i6) = zab2(i1,i2,i3,i4)*za(i4,i6) + 
     &   zab2(i1,i2,i3,i5)*za(i5,i6)
      zbb2(i1,i2,i3,i4,i5,i6) = zb(i1,i2)*zab2(i2,i4,i5,i6) +
     &   zb(i1,i3)*zab2(i3,i4,i5,i6)
      zab3(i1,i2,i3,i4,i5,i6,i7,i8)=zaa2(i1,i2,i3,i4,i5,i6)*zb(i6,i8) +
     &   zaa2(i1,i2,i3,i4,i5,i7)*zb(i7,i8)
      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      sqrtd(i1,i2,i3,i4)=sqrt(cplx2(
     & (s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4))**2-4*s(i1,i2)*s(i3,i4),
     & zip))

      sqrtd1245 = sqrtd(i1,i2,i4,i5)
      sqrtd1246 = sqrtd(i1,i2,i4,i6)
      sqrtd2345 = sqrtd(i2,i3,i4,i5)
      sqrtd2346 = sqrtd(i2,i3,i4,i6)
      sqrtd2435 = sqrtd(i2,i4,i3,i5)
      sqrtd2436 = sqrtd(i2,i4,i3,i6)


      aaaa_NMHV_c(:) = czip

      aaaa_NMHV_c(2) = (
     &  (s(i2,i5)*za(i1,i2)*za(i4,i6)*zb(i2,i1)*zb(i6,i1))/
     &  (za(i2,i5)*za(i2,i6)*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)) + 
     & (s(i2,i5)*za(i1,i2)*za(i3,i6)*zb(i2,i1)*zb(i6,i1))/
     &  (za(i2,i5)*za(i2,i6)*za(i6,i5)*zb(i3,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (za(i1,i2)*za(i1,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)**2*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)) - 
     & (s(i2,i3)*za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i3,i6)*zb(i3,i1))**2*zb(i3,i2)*
     &    zb(i4,i1)) - (za(i1,i2)*za(i1,i6)*za(i3,i6)*za(i4,i6)*
     &    zb(i2,i1)**2*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i3,i6)*zb(i3,i1))**2*zb(i3,i2)*
     &    zb(i4,i1)) - (za(i1,i3)*za(i1,i6)*za(i3,i6)*za(i4,i6)*
     &    zb(i2,i1)*zb(i3,i1)*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i3,i6)*zb(i3,i1))**2*zb(i3,i2)*
     &    zb(i4,i1)) + (za(i1,i2)*za(i1,i6)*za(i3,i6)*za(i4,i6)*
     &    zb(i2,i1)*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)**2*za(i6,i5)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i3,i6)*zb(i3,i1))*zb(i3,i2)*
     &    zb(i4,i1)) + (za(i1,i2)*za(i1,i6)*za(i3,i6)*zb(i2,i1)*
     &    zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)**2*za(i6,i5)*zb(i3,i1)*zb(i4,i1)*zb(i4,i2)) - 
     & (s(i2,i4)*za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i4,i6)*zb(i4,i1))**2*zb(i4,i2))
     &  - (za(i1,i2)*za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)**2*
     &    zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i4,i6)*zb(i4,i1))**2*zb(i4,i2))
     &  - (za(i1,i4)*za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)*zb(i4,i1)*
     &    zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i4,i6)*zb(i4,i1))**2*zb(i4,i2))
     &  + (za(i1,i2)*za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i6)**2*za(i6,i5)*zb(i3,i1)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i4,i6)*zb(i4,i1))*zb(i4,i2)) + 
     & (za(i1,i6)*za(i3,i4)*za(i3,i6)*zb(i2,i1)*zb(i5,i3)*zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i4,i6)*zb(i4,i1))*zb(i4,i2)) - 
     & (za(i1,i6)*za(i3,i4)*za(i4,i6)*zb(i2,i1)*zb(i5,i4)*zb(i6,i1))/
     &  (za(i2,i6)*za(i6,i5)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i3,i6)*zb(i3,i1))*zb(i3,i2)*
     &    zb(i4,i1)) + (s(i2,i5)*za(i4,i6)*zb(i6,i1)**2)/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)**2*zb(i4,i1)) + 
     & (s(i2,i5)*za(i3,i6)*zb(i3,i2)*zb(i6,i1)**2)/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)**2*zb(i4,i1)*zb(i4,i2)) + 
     & (za(i1,i6)*za(i4,i6)*zb(i5,i1)*zb(i6,i1)**2)/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)**2*zb(i4,i1)) + 
     & (za(i1,i6)*za(i3,i6)*zb(i3,i2)*zb(i5,i1)*zb(i6,i1)**2)/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)**2*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i3)*zb(i6,i1)**2)/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)**2*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i4,i6)*zb(i4,i1))*zb(i4,i2)) - 
     & (za(i1,i6)*za(i3,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i4)*zb(i6,i1)**2)/
     &  (za(i2,i6)*za(i6,i5)*
     &    (-(za(i2,i6)*zb(i2,i1)) - za(i3,i6)*zb(i3,i1))*zb(i3,i2)*
     &    zb(i4,i1)**2) + (za(i2,i6)*za(i4,i6)*zb(i2,i1)*zb(i6,i1)**2*
     &    zb(i6,i3))/(za(i2,i5)*za(i6,i5)*zb(i3,i1)**3*zb(i4,i1)) + 
     & (za(i2,i6)*za(i3,i6)*zb(i2,i1)*zb(i3,i2)*zb(i6,i1)**2*zb(i6,i3))/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)**3*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i4,i6)*zb(i5,i1)*zb(i6,i1)**2*zb(i6,i3))/
     &  (za(i2,i5)*zb(i3,i1)**3*zb(i4,i1)) - 
     & (za(i3,i6)*zb(i3,i2)*zb(i5,i1)*zb(i6,i1)**2*zb(i6,i3))/
     &  (za(i2,i5)*zb(i3,i1)**3*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i3,i4)*za(i3,i6)*zb(i5,i2)*zb(i6,i1)*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*za(i2,i5)*zb(i3,i1)**2*zb(i4,i2)) + 
     & (za(i3,i6)*za(i4,i6)*zb(i5,i2)*zb(i6,i1)**2*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*za(i2,i5)*zb(i3,i1)**3*zb(i4,i2)) - 
     & (s(i2,i5)*za(i4,i6)*zb(i2,i1)*zb(i6,i1)*zb(i6,i4))/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**2) - 
     & (s(i2,i5)*za(i3,i6)*zb(i2,i1)*zb(i6,i1)*zb(i6,i4))/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)) - 
     & (za(i1,i6)*za(i4,i6)*zb(i2,i1)*zb(i5,i1)*zb(i6,i1)*zb(i6,i4))/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**2) - 
     & (za(i1,i6)*za(i3,i6)*zb(i2,i1)*zb(i5,i1)*zb(i6,i1)*zb(i6,i4))/
     &  (za(i2,i6)*za(i6,i5)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)) + 
     & (za(i2,i6)*za(i4,i6)*zb(i2,i1)*zb(i6,i1)**2*zb(i6,i4))/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)**2*zb(i4,i1)**2) + 
     & (za(i2,i6)*za(i3,i6)*zb(i2,i1)*zb(i3,i2)*zb(i6,i1)**2*zb(i6,i4))/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)**2*zb(i4,i1)**2*zb(i4,i2)) - 
     & (za(i4,i6)*zb(i5,i1)*zb(i6,i1)**2*zb(i6,i4))/
     &  (za(i2,i5)*zb(i3,i1)**2*zb(i4,i1)**2) - 
     & (za(i3,i6)*zb(i3,i2)*zb(i5,i1)*zb(i6,i1)**2*zb(i6,i4))/
     &  (za(i2,i5)*zb(i3,i1)**2*zb(i4,i1)**2*zb(i4,i2)) - 
     & (za(i2,i6)*za(i4,i6)*zb(i2,i1)**2*zb(i6,i1)*zb(i6,i4)**2)/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**3) - 
     & (za(i2,i6)*za(i3,i6)*zb(i2,i1)**2*zb(i6,i1)*zb(i6,i4)**2)/
     &  (za(i2,i5)*za(i6,i5)*zb(i3,i1)*zb(i4,i1)**3*zb(i4,i2)) + 
     & (za(i4,i6)*zb(i2,i1)*zb(i5,i1)*zb(i6,i1)*zb(i6,i4)**2)/
     &  (za(i2,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**3) + 
     & (za(i3,i6)*zb(i2,i1)*zb(i5,i1)*zb(i6,i1)*zb(i6,i4)**2)/
     &  (za(i2,i5)*zb(i3,i1)*zb(i4,i1)**3*zb(i4,i2)) + 
     & (za(i3,i4)*za(i4,i6)*zb(i5,i2)*zb(i6,i1)*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*za(i2,i5)*zb(i3,i2)*zb(i4,i1)**2) + 
     & (za(i3,i6)*za(i4,i6)*zb(i5,i2)*zb(i6,i1)**2*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*za(i2,i5)*zb(i3,i2)*zb(i4,i1)**3)
     &)


      aaaa_NMHV_c(4) = (
     &  (s(i2,i6)*za(i1,i2)*za(i4,i5)*zb(i2,i1)*zb(i5,i1))/
     &  (za(i2,i5)*za(i2,i6)*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)) + 
     & (s(i2,i6)*za(i1,i2)*za(i3,i5)*zb(i2,i1)*zb(i5,i1))/
     &  (za(i2,i5)*za(i2,i6)*za(i5,i6)*zb(i3,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (s(i2,i6)*za(i4,i5)*zb(i5,i1)**2)/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)**2*zb(i4,i1)) + 
     & (s(i2,i6)*za(i3,i5)*zb(i3,i2)*zb(i5,i1)**2)/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)**2*zb(i4,i1)*zb(i4,i2)) + 
     & (za(i2,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)**2*zb(i5,i3))/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)**3*zb(i4,i1)) + 
     & (za(i2,i5)*za(i3,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i1)**2*zb(i5,i3))/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)**3*zb(i4,i1)*zb(i4,i2)) - 
     & (s(i2,i6)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*zb(i5,i4))/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**2) - 
     & (s(i2,i6)*za(i3,i5)*zb(i2,i1)*zb(i5,i1)*zb(i5,i4))/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)) + 
     & (za(i2,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)**2*zb(i5,i4))/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)**2*zb(i4,i1)**2) + 
     & (za(i2,i5)*za(i3,i5)*zb(i2,i1)*zb(i3,i2)*zb(i5,i1)**2*zb(i5,i4))/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)**2*zb(i4,i1)**2*zb(i4,i2)) - 
     & (za(i2,i5)*za(i4,i5)*zb(i2,i1)**2*zb(i5,i1)*zb(i5,i4)**2)/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**3) - 
     & (za(i2,i5)*za(i3,i5)*zb(i2,i1)**2*zb(i5,i1)*zb(i5,i4)**2)/
     &  (za(i2,i6)*za(i5,i6)*zb(i3,i1)*zb(i4,i1)**3*zb(i4,i2)) + 
     & (za(i1,i2)*za(i1,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)**2*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)) - 
     & (s(i2,i3)*za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i3,i5)*zb(i3,i1))**2*zb(i3,i2)*
     &    zb(i4,i1)) - (za(i1,i2)*za(i1,i5)*za(i3,i5)*za(i4,i5)*
     &    zb(i2,i1)**2*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i3,i5)*zb(i3,i1))**2*zb(i3,i2)*
     &    zb(i4,i1)) - (za(i1,i3)*za(i1,i5)*za(i3,i5)*za(i4,i5)*
     &    zb(i2,i1)*zb(i3,i1)*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i3,i5)*zb(i3,i1))**2*zb(i3,i2)*
     &    zb(i4,i1)) + (za(i1,i2)*za(i1,i5)*za(i3,i5)*za(i4,i5)*
     &    zb(i2,i1)*zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)**2*za(i5,i6)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i3,i5)*zb(i3,i1))*zb(i3,i2)*
     &    zb(i4,i1)) + (za(i1,i2)*za(i1,i5)*za(i3,i5)*zb(i2,i1)*
     &    zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)**2*za(i5,i6)*zb(i3,i1)*zb(i4,i1)*zb(i4,i2)) - 
     & (s(i2,i4)*za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i4,i5)*zb(i4,i1))**2*zb(i4,i2))
     &  - (za(i1,i2)*za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)**2*
     &    zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i4,i5)*zb(i4,i1))**2*zb(i4,i2))
     &  - (za(i1,i4)*za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i4,i1)*
     &    zb(i5,i1)*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i4,i5)*zb(i4,i1))**2*zb(i4,i2))
     &  + (za(i1,i2)*za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i5)**2*za(i5,i6)*zb(i3,i1)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i4,i5)*zb(i4,i1))*zb(i4,i2)) + 
     & (za(i1,i5)*za(i4,i5)*zb(i5,i1)**2*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)**2*zb(i4,i1)) + 
     & (za(i1,i5)*za(i3,i5)*zb(i3,i2)*zb(i5,i1)**2*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)**2*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i4,i5)*zb(i5,i1)**2*zb(i5,i3)*zb(i6,i1))/
     &  (za(i2,i6)*zb(i3,i1)**3*zb(i4,i1)) - 
     & (za(i3,i5)*zb(i3,i2)*zb(i5,i1)**2*zb(i5,i3)*zb(i6,i1))/
     &  (za(i2,i6)*zb(i3,i1)**3*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i1,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*zb(i5,i4)*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**2) - 
     & (za(i1,i5)*za(i3,i5)*zb(i2,i1)*zb(i5,i1)*zb(i5,i4)*zb(i6,i1))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)) - 
     & (za(i4,i5)*zb(i5,i1)**2*zb(i5,i4)*zb(i6,i1))/
     &  (za(i2,i6)*zb(i3,i1)**2*zb(i4,i1)**2) - 
     & (za(i3,i5)*zb(i3,i2)*zb(i5,i1)**2*zb(i5,i4)*zb(i6,i1))/
     &  (za(i2,i6)*zb(i3,i1)**2*zb(i4,i1)**2*zb(i4,i2)) + 
     & (za(i4,i5)*zb(i2,i1)*zb(i5,i1)*zb(i5,i4)**2*zb(i6,i1))/
     &  (za(i2,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**3) + 
     & (za(i3,i5)*zb(i2,i1)*zb(i5,i1)*zb(i5,i4)**2*zb(i6,i1))/
     &  (za(i2,i6)*zb(i3,i1)*zb(i4,i1)**3*zb(i4,i2)) - 
     & (za(i3,i4)*za(i3,i5)*zb(i5,i1)*zb(i5,i3)**2*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zb(i3,i1)**2*zb(i4,i2)) + 
     & (za(i3,i5)*za(i4,i5)*zb(i5,i1)**2*zb(i5,i3)**2*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zb(i3,i1)**3*zb(i4,i2)) + 
     & (za(i3,i4)*za(i4,i5)*zb(i5,i1)*zb(i5,i4)**2*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zb(i3,i2)*zb(i4,i1)**2) + 
     & (za(i3,i5)*za(i4,i5)*zb(i5,i1)**2*zb(i5,i4)**2*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zb(i3,i2)*zb(i4,i1)**3) + 
     & (za(i1,i5)*za(i3,i4)*za(i3,i5)*zb(i2,i1)*zb(i5,i1)*zb(i6,i3))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i4,i5)*zb(i4,i1))*zb(i4,i2)) - 
     & (za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)**2*zb(i6,i3))/
     &  (za(i2,i5)*za(i5,i6)*zb(i3,i1)**2*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i4,i5)*zb(i4,i1))*zb(i4,i2)) - 
     & (za(i1,i5)*za(i3,i4)*za(i4,i5)*zb(i2,i1)*zb(i5,i1)*zb(i6,i4))/
     &  (za(i2,i5)*za(i5,i6)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i3,i5)*zb(i3,i1))*zb(i3,i2)*
     &    zb(i4,i1)) - (za(i1,i5)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*
     &    zb(i5,i1)**2*zb(i6,i4))/
     &  (za(i2,i5)*za(i5,i6)*
     &    (-(za(i2,i5)*zb(i2,i1)) - za(i3,i5)*zb(i3,i1))*zb(i3,i2)*
     &    zb(i4,i1)**2)
     &)


      aaaa_NMHV_c(5) = (
     &  -((za(i1,i6)**2*za(i4,i6)*
     &      (2*s(i1,i5)**2 + zab2(i4,i1,i5,i4)**2 + 
     &        2*s(i1,i5)*(zab2(i1,i4,i6,i1) + zab2(i5,i4,i6,i5)) + 
     &        zab2(i6,i1,i5,i6)**2)*zb(i4,i2)**2*zb(i6,i4))/
     &    (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i4)**3*zb(i3,i2)*zb(i4,i3))
     &    ) - (za(i4,i6)*(za(i1,i4)**2*zb(i4,i2)**2 + 
     &      2*za(i1,i4)*za(i1,i5)*zb(i4,i2)*zb(i5,i2) + 
     &      2*za(i1,i5)**2*zb(i5,i2)**2)*zb(i6,i4))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i4)*zb(i3,i2)*zb(i4,i3)) - 
     & (za(i1,i6)*za(i3,i6)*za(i4,i6)*
     &    (2*s(i3,i5) + zab2(i3,i4,i6,i3) + zab2(i5,i4,i6,i5))*
     &    zb(i5,i4)*zb(i6,i4))/
     &  (za(i2,i6)*za(i5,i6)*zab2(i6,i3,i5,i4)**2*zb(i4,i3)) + 
     & (za(i3,i6)*za(i4,i6)*(2*s(i1,i2) + zab2(i1,i4,i6,i1) + 
     &      zab2(i2,i4,i6,i2))*zb(i4,i2)*zb(i5,i4)*zb(i6,i4))/
     &  (za(i5,i6)*zab2(i6,i1,i2,i4)**2*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i1,i6)*za(i4,i6)*(za(i3,i4)*zb(i4,i2) + 
     &      2*za(i3,i5)*zb(i5,i2))*zb(i5,i4)*zb(i6,i4))/
     &  (za(i2,i6)*za(i5,i6)*zab2(i6,i3,i5,i4)*zb(i4,i2)*zb(i4,i3)) - 
     & (za(i4,i6)*za(i5,i6)*(s(i2,i5) + zab2(i4,i2,i5,i4))**2*
     &    zb(i5,i4)**2*zb(i6,i4))/
     &  (za(i2,i5)*zab2(i6,i2,i5,i4)**3*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i2,i6)*za(i4,i6)*(s(i2,i5) + zab2(i4,i2,i5,i4))**2*zb(i4,i2)*
     &    (za(i2,i6)*zb(i4,i2) - 2*za(i5,i6)*zb(i5,i4))*zb(i6,i4))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i2,i5,i4)**3*zb(i3,i1)*zb(i4,i3))
     &  - (za(i1,i6)*za(i4,i6)*
     &    (za(i5,i6)*(za(i1,i4)*zb(i4,i2)**2*zb(i4,i3)*
     &          (za(i4,i5)*zb(i4,i3) + za(i5,i6)*zb(i6,i3)) + 
     &         za(i1,i6)*za(i5,i6)*
     &          (zb(i4,i3)**2*zb(i6,i2)**2 + zb(i4,i2)**2*zb(i6,i3)**2))
     &        + za(i1,i5)*zb(i4,i3)*
     &       (za(i4,i5)*zb(i4,i2)*zb(i4,i3)*
     &          (za(i4,i6)*zb(i4,i2) + 2*za(i5,i6)*zb(i5,i2)) + 
     &         za(i5,i6)*(za(i4,i6)*zb(i4,i2)**2*zb(i6,i3) + 
     &            2*za(i5,i6)*zb(i5,i2)*
     &             (zb(i4,i3)*zb(i6,i2) + zb(i4,i2)*zb(i6,i3)))))*
     &    zb(i6,i4))/
     &  (za(i1,i5)*za(i5,i6)**3*zab2(i6,i1,i5,i4)*zb(i3,i2)*
     &    zb(i4,i3)**3) - (za(i4,i6)*zb(i6,i4)*
     &    (-(zb(i4,i1)*zb(i4,i2)*zb(i4,i3)*
     &         (-2*za(i3,i6)*za(i5,i6)*zab2(i1,i3,i4,i5) + 
     &           za(i1,i6)*za(i3,i5)*za(i4,i6)*zb(i5,i4))) + 
     &      za(i1,i6)*za(i3,i6)*za(i5,i6)*
     &       (zb(i4,i1)*zb(i4,i2)*zb(i5,i3) + 
     &         zb(i2,i1)*zb(i4,i3)*zb(i5,i4))*zb(i6,i4)))/
     &  (za(i1,i6)*za(i5,i6)**2*zab2(i6,i1,i2,i4)*zb(i4,i1)**2*
     &    zb(i4,i3)**2) + (za(i1,i6)*za(i4,i6)*zb(i4,i2)*zb(i6,i4)*
     &    (za(i1,i4)**2*za(i5,i6)*zb(i4,i1)*zb(i4,i2)*zb(i4,i3) + 
     &      za(i1,i5)*za(i4,i5)*za(i4,i6)*zb(i4,i2)*zb(i4,i3)*
     &       zb(i5,i4) + 2*za(i1,i5)*za(i4,i5)*za(i5,i6)*zb(i4,i3)*
     &       zb(i5,i2)*zb(i5,i4) + 
     &      2*za(i1,i5)*za(i1,i6)*za(i5,i6)*zb(i4,i3)*zb(i5,i2)*
     &       zb(i6,i1) + 2*za(i1,i6)**2*za(i5,i6)*zb(i4,i3)*zb(i6,i1)*
     &       zb(i6,i2) + za(i1,i6)*za(i4,i5)*za(i5,i6)*zb(i4,i2)*
     &       zb(i5,i4)*zb(i6,i3) + 
     &      za(i1,i4)*(za(i1,i5)*zb(i4,i1)*zb(i4,i3)*
     &          (za(i4,i6)*zb(i4,i2) + 2*za(i5,i6)*zb(i5,i2)) + 
     &         za(i5,i6)*zb(i4,i2)*
     &          (za(i4,i5)*zb(i4,i3)*zb(i5,i4) + 
     &            za(i1,i6)*zb(i4,i1)*zb(i6,i3))) + 
     &      s(i1,i5)*(za(i1,i5)*zb(i4,i3)*
     &          (za(i4,i6)*zb(i4,i2) + 4*za(i5,i6)*zb(i5,i2)) + 
     &         za(i5,i6)*(zb(i4,i3)*
     &             (za(i1,i4)*zb(i4,i2) + 2*za(i1,i6)*zb(i6,i2)) + 
     &            za(i1,i6)*zb(i4,i2)*zb(i6,i3))) + 
     &      2*za(i5,i6)**2*zab2(i1,i5,i6,i2)*zb(i4,i3)*zb(i6,i5)))/
     &  (za(i1,i5)*za(i5,i6)**2*zab2(i6,i1,i5,i4)**2*zb(i3,i2)*
     &    zb(i4,i3)**2) - (za(i4,i6)*zb(i6,i4)*
     &    (za(i2,i6)*(za(i2,i5)*za(i4,i6)*zb(i4,i2)*zb(i4,i3)*
     &          (za(i4,i5)*zb(i4,i2) - 2*za(i5,i6)*zb(i6,i2)) + 
     &         za(i5,i6)*(zb(i4,i3)*
     &             (2*za(i5,i6)*zb(i6,i2)*
     &                (s(i2,i5) + za(i2,i6)*zb(i6,i2)) + 
     &               za(i2,i4)*zb(i4,i2)*
     &                (za(i4,i5)*zb(i4,i2) - 2*za(i5,i6)*zb(i6,i2))) + 
     &            za(i4,i5)*zb(i4,i2)*
     &             (za(i2,i6)*zb(i4,i2) - 2*za(i5,i6)*zb(i5,i4))*
     &             zb(i6,i3))) + 
     &      2*za(i5,i6)**3*zb(i4,i3)*
     &       (s(i2,i5) - za(i2,i4)*zb(i4,i2) + 2*za(i2,i6)*zb(i6,i2))*
     &       zb(i6,i5)))/
     &  (za(i2,i5)*za(i5,i6)**3*zab2(i6,i2,i5,i4)*zb(i3,i1)*
     &    zb(i4,i3)**2) - (za(i4,i6)*zb(i6,i4)*
     &    (s(i2,i5)**2*zb(i4,i3) + 
     &      za(i2,i4)**2*zb(i4,i2)**2*zb(i4,i3) + 
     &      2*za(i2,i4)*za(i5,i6)*zb(i4,i3)*zb(i5,i4)*zb(i6,i2) + 
     &      za(i5,i6)*(2*za(i5,i6)*zb(i4,i3)*zb(i6,i5)**2 + 
     &         za(i4,i5)*zb(i5,i4)*
     &          (zb(i5,i4)*zb(i6,i3) - 2*zb(i4,i3)*zb(i6,i5)))))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i2,i5,i4)*zb(i3,i1)*zb(i4,i3)**2)
     &  - (za(i4,i6)*zb(i5,i4)*zb(i6,i4)*
     &    (za(i4,i5)**2*zb(i4,i3)*zb(i5,i4)**2 - 
     &      za(i2,i4)*za(i5,i6)*zb(i4,i2)*
     &       (zb(i5,i3)*zb(i6,i4) + zb(i4,i3)*zb(i6,i5)) - 
     &      za(i4,i5)*za(i5,i6)*zb(i5,i4)*
     &       (zb(i5,i3)*zb(i6,i4) + zb(i4,i3)*zb(i6,i5)) + 
     &      s(i2,i5)*(za(i4,i5)*zb(i4,i3)*zb(i5,i4) - 
     &         za(i5,i6)*(zb(i5,i3)*zb(i6,i4) + zb(i4,i3)*zb(i6,i5)))))/
     &  (za(i2,i5)*zab2(i6,i2,i5,i4)**2*zb(i3,i1)*zb(i4,i3)**2) - 
     & (za(i3,i6)*za(i4,i6)*zb(i6,i4)*
     &    (za(i1,i4)*za(i2,i6)*za(i5,i6)*zb(i4,i2)*zb(i4,i3)*
     &       zb(i5,i4) + za(i1,i6)*
     &       (-(za(i2,i5)*za(i4,i6)*zb(i4,i2)*zb(i4,i3)*zb(i5,i4)) + 
     &         za(i2,i6)*za(i5,i6)*
     &          (zb(i5,i4)*(-(zb(i4,i3)*zb(i6,i2)) + 
     &               zb(i3,i2)*zb(i6,i4)) + 
     &            zb(i4,i2)*zb(i4,i3)*zb(i6,i5)))))/
     &  (za(i2,i6)**2*za(i5,i6)**2*zab2(i6,i3,i5,i4)*zb(i4,i2)*
     &    zb(i4,i3)**2) + (za(i4,i6)*zb(i6,i4)*
     &    (za(i2,i4)**2*za(i5,i6)*zb(i4,i2)**2*zb(i4,i3)*
     &       (za(i2,i6)*zb(i4,i2) - 2*za(i5,i6)*zb(i5,i4)) + 
     &      za(i2,i6)*za(i4,i5)*zb(i5,i4)*
     &       (za(i2,i5)*za(i4,i6)*zb(i4,i2)**2*zb(i4,i3) - 
     &         za(i5,i6)*(za(i2,i6)*zb(i4,i2)*
     &             (zb(i4,i3)*zb(i6,i2) - zb(i3,i2)*zb(i6,i4)) + 
     &            2*za(i5,i6)*
     &             (zb(i3,i2)*zb(i5,i4)*zb(i6,i4) + 
     &               zb(i4,i2)*zb(i4,i3)*zb(i6,i5)))) + 
     &      za(i2,i4)*zb(i4,i2)*
     &       (za(i2,i5)*za(i2,i6)*za(i4,i6)*zb(i4,i2)**2*zb(i4,i3) + 
     &         za(i5,i6)*(-3*za(i4,i5)*za(i5,i6)*zb(i4,i3)*
     &             zb(i5,i4)**2 + 
     &            za(i2,i6)**2*zb(i4,i2)*
     &             (-(zb(i4,i3)*zb(i6,i2)) + zb(i3,i2)*zb(i6,i4)) + 
     &            za(i2,i6)*(za(i4,i5)*zb(i4,i2)*zb(i4,i3)*zb(i5,i4) - 
     &               2*za(i5,i6)*
     &                (zb(i3,i2)*zb(i5,i4)*zb(i6,i4) + 
     &                  zb(i4,i2)*zb(i4,i3)*zb(i6,i5))))) + 
     &      s(i2,i5)*(za(i2,i5)*za(i2,i6)*za(i4,i6)*zb(i4,i2)**2*
     &          zb(i4,i3) + za(i5,i6)*
     &          (za(i2,i4)*zb(i4,i2)*zb(i4,i3)*
     &             (za(i2,i6)*zb(i4,i2) - 2*za(i5,i6)*zb(i5,i4)) - 
     &            za(i2,i6)*(za(i2,i6)*zb(i4,i2)*
     &                (zb(i4,i3)*zb(i6,i2) - zb(i3,i2)*zb(i6,i4)) + 
     &               2*za(i5,i6)*
     &                (zb(i3,i2)*zb(i5,i4)*zb(i6,i4) + 
     &                  zb(i4,i2)*zb(i4,i3)*zb(i6,i5)))))))/
     &  (za(i2,i5)*za(i5,i6)**2*zab2(i6,i2,i5,i4)**2*zb(i3,i1)*
     &    zb(i4,i3)**2) - (za(i4,i6)*zb(i6,i4)*
     &    (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)**3*za(i5,i6)**3*zb(i3,i1)*
     &       zb(i4,i2)**2*zb(i4,i3)**3*zb(i5,i2)*
     &       (-(za(i1,i4)**2*zab2(i3,i2,i5,i4)*zb(i4,i1)**2) + 
     &         za(i1,i6)*zab2(i1,i4,i6,i1)*
     &          (za(i2,i3)*zb(i2,i1) - za(i3,i5)*zb(i5,i1))*zb(i6,i4))
     &       + t(i2,i3,i5)*(za(i1,i3)*za(i1,i6)**2*za(i2,i5)*
     &          za(i5,i6)**3*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)**3*
     &          (-(za(i2,i4)*za(i4,i6)*zab2(i2,i1,i3,i5)*
     &               zb(i4,i2)**2) - 
     &            za(i2,i6)*za(i4,i6)*zab2(i2,i1,i3,i5)*zb(i4,i2)*
     &             zb(i6,i2) + 
     &            za(i2,i6)**2*zab2(i6,i1,i3,i5)*zb(i6,i2)**2) + 
     &         t(i2,i4,i6)*(-(za(i1,i5)*za(i1,i6)**2*za(i2,i6)*
     &               (-(za(i2,i5)**2*za(i4,i6)*za(i5,i6)**2*zb(i3,i2)*
     &                    zb(i4,i1)**3*zb(i4,i2)*zb(i4,i3)**2*zb(i5,i2)*
     &                    zb(i5,i4)) + 
     &                 za(i2,i6)**2*za(i5,i6)**2*zb(i4,i2)**2*
     &                  (-2*za(i2,i4)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i2)*
     &                     zb(i4,i3)*zb(i6,i3) - 
     &                    2*za(i4,i5)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)*
     &                     zb(i5,i4)*zb(i6,i3) + 
     &                    2*za(i2,i6)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)*
     &                     zb(i6,i2)*zb(i6,i3) + 
     &                    za(i2,i6)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)**2*
     &                     zb(i6,i3)**2 - 
     &                    za(i5,i6)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)*
     &                     zb(i5,i4)*zb(i6,i3)**2 + 
     &                    za(i2,i6)*zb(i2,i1)**2*zb(i3,i1)*zb(i4,i3)**2*
     &                     zb(i6,i4)**2 + 
     &                    za(i5,i6)*zb(i2,i1)*zb(i3,i1)*zb(i4,i3)**2*
     &                     zb(i5,i1)*zb(i6,i4)**2 + 
     &                    2*za(i5,i6)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)*
     &                     zb(i6,i3)*zb(i6,i5)) + 
     &                 za(i2,i5)*za(i2,i6)*zb(i4,i1)*zb(i4,i3)*
     &                  (za(i2,i6)*zb(i4,i2)*
     &                     (za(i4,i5)*zb(i4,i1)*zb(i4,i2)*zb(i4,i3)*
     &                        (za(i4,i6)*zb(i4,i2)*
     &                        (zb(i3,i2)*zb(i4,i1) + 
     &                       zb(i3,i1)*zb(i4,i2)) + 
     &                        za(i5,i6)*zb(i3,i2)*zb(i4,i1)*zb(i5,i2))
     &                        + za(i5,i6)*
     &                        (za(i5,i6)*zb(i3,i2)*zb(i4,i1)**2*
     &                        zb(i5,i2)*
     &                        (zb(i4,i3)*zb(i6,i2) - 
     &                       zb(i3,i2)*zb(i6,i4)) - 
     &                        za(i4,i6)*zb(i2,i1)*zb(i4,i2)*zb(i4,i3)*
     &                        (zb(i4,i1)*zb(i4,i3)*zb(i6,i2) + 
     &                        zb(i3,i1)*zb(i4,i2)*zb(i6,i4)))) + 
     &                    za(i5,i6)*zb(i3,i2)*zb(i4,i1)**2*
     &                     (zb(i5,i4)*
     &                        (-(zb(i4,i2)*zb(i4,i3)*
     &                        (-(za(i4,i5)*zab2(i6,i4,i5,i2)) + 
     &                        za(i4,i6)*za(i5,i6)*zb(i6,i2))) + 
     &                        za(i5,i6)**2*zb(i3,i2)*zb(i5,i2)*zb(i6,i4)
     &                        ) + 
     &                       za(i5,i6)**2*zb(i4,i2)*zb(i4,i3)*zb(i5,i2)*
     &                        zb(i6,i5))))) + 
     &            za(i1,i5)**2*za(i2,i6)**3*za(i5,i6)*zb(i3,i1)*
     &             zb(i4,i1)*zb(i4,i2)**2*zb(i4,i3)*zb(i5,i2)*
     &             (za(i1,i4)*za(i5,i6)*zab2(i6,i2,i5,i4)*zb(i4,i1)*
     &                zb(i4,i3) + 
     &               za(i1,i6)*
     &                (za(i2,i5)*za(i4,i6)*zb(i4,i1)*zb(i4,i2)*
     &                   zb(i4,i3) - 
     &                  za(i5,i6)*
     &                   (za(i2,i6)*
     &                      (zb(i4,i1)*zb(i4,i2)*zb(i6,i3) + 
     &                        zb(i2,i1)*zb(i4,i3)*zb(i6,i4)) + 
     &                     za(i5,i6)*
     &                      (-(zb(i3,i1)*zb(i5,i4)*zb(i6,i4)) + 
     &                        zb(i4,i1)*zb(i4,i3)*zb(i6,i5))))) + 
     &            za(i1,i6)**2*za(i5,i6)**2*zb(i3,i2)*zb(i4,i1)**3*
     &             (za(i1,i2)*za(i2,i6)**2*za(i5,i6)*zab2(i6,i2,i5,i4)*
     &                zb(i4,i2)**2*zb(i6,i3)**2 + 
     &               za(i1,i6)*za(i2,i5)*zb(i4,i3)*
     &                (za(i2,i4)*za(i2,i5)*za(i4,i6)*zb(i4,i2)**2*
     &                   zb(i4,i3)*zb(i5,i4) + 
     &                  za(i2,i6)*
     &                   (za(i2,i6)*zb(i3,i2)*zb(i6,i4)*
     &                      ((2*za(i2,i6)*zb(i4,i2) - 
     &                        za(i5,i6)*zb(i5,i4))*zb(i6,i2) + 
     &                        za(i5,i6)*zb(i4,i2)*zb(i6,i5)) - 
     &                     za(i2,i5)*za(i4,i6)*zb(i4,i2)*
     &                      (zb(i3,i2)*zb(i5,i4)*zb(i6,i4) + 
     &                        zb(i4,i2)*zb(i4,i3)*zb(i6,i5)))))))))/
     &  (t(i2,i3,i5)*t(i2,i4,i6)*za(i1,i5)*za(i1,i6)**2*za(i2,i5)*
     &    za(i2,i6)**3*za(i5,i6)**3*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**3*
     &    zb(i4,i2)**2*zb(i4,i3)**3)
     &)


      aaaa_NMHV_c(6) = (
     &  -((za(i1,i4)*(-(t(i2,i3,i6)*za(i1,i4)*za(i1,i5)*za(i2,i6)*
     &          zab2(i3,i1,i4,i6)*zb(i4,i3)*zb(i5,i2)) + 
     &       t(i2,i3,i5)*(t(i2,i3,i6)*t(i2,i5,i6)*za(i1,i2)*
     &           zab2(i1,i3,i4,i2) - 
     &          za(i1,i4)*za(i1,i6)*za(i2,i5)*zab2(i3,i1,i4,i5)*
     &           zb(i4,i3)*zb(i6,i2))))/
     &   (t(i2,i3,i5)*t(i2,i3,i6)*za(i1,i5)*za(i1,i6)*za(i2,i5)*
     &     za(i2,i6)*zb(i3,i2)*zb(i4,i3)))
     &)


      aaaa_NMHV_c(7) = (
     &  -((za(i1,i5)**2*za(i4,i5)*
     &      (2*s(i1,i6)**2 + zab2(i4,i1,i6,i4)**2 + 
     &        zab2(i5,i1,i6,i5)**2 + 
     &        2*s(i1,i6)*(zab2(i1,i4,i5,i1) + zab2(i6,i4,i5,i6)))*
     &      zb(i4,i2)**2*zb(i5,i4))/
     &    (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i4)**3*zb(i3,i2)*zb(i4,i3))
     &    ) - (za(i4,i5)*zb(i5,i4)*
     &    (za(i1,i4)**2*zb(i4,i2)**2 + 
     &      2*za(i1,i4)*za(i1,i6)*zb(i4,i2)*zb(i6,i2) + 
     &      2*za(i1,i6)**2*zb(i6,i2)**2))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i4)*zb(i3,i2)*zb(i4,i3)) - 
     & (za(i1,i5)*za(i4,i5)*zb(i5,i4)*
     &    (za(i6,i5)*(za(i1,i4)*zb(i4,i2)**2*zb(i4,i3)*
     &          (za(i4,i6)*zb(i4,i3) + za(i6,i5)*zb(i5,i3)) + 
     &         za(i1,i5)*za(i6,i5)*
     &          (zb(i4,i3)**2*zb(i5,i2)**2 + zb(i4,i2)**2*zb(i5,i3)**2))
     &        + za(i1,i6)*zb(i4,i3)*
     &       (za(i4,i6)*zb(i4,i2)*zb(i4,i3)*
     &          (za(i4,i5)*zb(i4,i2) + 2*za(i6,i5)*zb(i6,i2)) + 
     &         za(i6,i5)*(za(i4,i5)*zb(i4,i2)**2*zb(i5,i3) + 
     &            2*za(i6,i5)*
     &             (zb(i4,i3)*zb(i5,i2) + zb(i4,i2)*zb(i5,i3))*zb(i6,i2)
     &            ))))/
     &  (za(i1,i6)*za(i6,i5)**3*zab2(i5,i1,i6,i4)*zb(i3,i2)*
     &    zb(i4,i3)**3) - (za(i1,i5)*za(i3,i5)*za(i4,i5)*
     &    (2*s(i3,i6) + zab2(i3,i4,i5,i3) + zab2(i6,i4,i5,i6))*
     &    zb(i5,i4)*zb(i6,i4))/
     &  (za(i2,i5)*za(i6,i5)*zab2(i5,i3,i6,i4)**2*zb(i4,i3)) + 
     & (za(i3,i5)*za(i4,i5)*(2*s(i1,i2) + zab2(i1,i4,i5,i1) + 
     &      zab2(i2,i4,i5,i2))*zb(i4,i2)*zb(i5,i4)*zb(i6,i4))/
     &  (za(i6,i5)*zab2(i5,i1,i2,i4)**2*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i1,i5)*za(i4,i5)*zb(i5,i4)*
     &    (za(i3,i4)*zb(i4,i2) + 2*za(i3,i6)*zb(i6,i2))*zb(i6,i4))/
     &  (za(i2,i5)*za(i6,i5)*zab2(i5,i3,i6,i4)*zb(i4,i2)*zb(i4,i3)) - 
     & (za(i4,i5)*za(i6,i5)*(s(i2,i6) + zab2(i4,i2,i6,i4))**2*zb(i5,i4)*
     &    zb(i6,i4)**2)/
     &  (za(i2,i6)*zab2(i5,i2,i6,i4)**3*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i2,i5)*za(i4,i5)*(s(i2,i6) + zab2(i4,i2,i6,i4))**2*zb(i4,i2)*
     &    zb(i5,i4)*(za(i2,i5)*zb(i4,i2) - 2*za(i6,i5)*zb(i6,i4)))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i2,i6,i4)**3*zb(i3,i1)*zb(i4,i3))
     &  - (za(i4,i5)*zb(i5,i4)*zb(i6,i4)*
     &    (-(za(i2,i4)*za(i6,i5)*zb(i4,i2)*
     &         (zb(i4,i3)*zb(i5,i6) + zb(i5,i4)*zb(i6,i3))) - 
     &      za(i4,i6)*za(i6,i5)*
     &       (zb(i4,i3)*zb(i5,i6) + zb(i5,i4)*zb(i6,i3))*zb(i6,i4) + 
     &      za(i4,i6)**2*zb(i4,i3)*zb(i6,i4)**2 + 
     &      s(i2,i6)*(-(za(i6,i5)*
     &            (zb(i4,i3)*zb(i5,i6) + zb(i5,i4)*zb(i6,i3))) + 
     &         za(i4,i6)*zb(i4,i3)*zb(i6,i4))))/
     &  (za(i2,i6)*zab2(i5,i2,i6,i4)**2*zb(i3,i1)*zb(i4,i3)**2) - 
     & (za(i4,i5)*zb(i5,i4)*(-(zb(i4,i1)*zb(i4,i2)*zb(i4,i3)*
     &         (-2*za(i3,i5)*za(i6,i5)*zab2(i1,i3,i4,i6) + 
     &           za(i1,i5)*za(i3,i6)*za(i4,i5)*zb(i6,i4))) + 
     &      za(i1,i5)*za(i3,i5)*za(i6,i5)*zb(i5,i4)*
     &       (zb(i4,i1)*zb(i4,i2)*zb(i6,i3) + 
     &         zb(i2,i1)*zb(i4,i3)*zb(i6,i4))))/
     &  (za(i1,i5)*za(i6,i5)**2*zab2(i5,i1,i2,i4)*zb(i4,i1)**2*
     &    zb(i4,i3)**2) + (za(i1,i5)*za(i4,i5)*zb(i4,i2)*zb(i5,i4)*
     &    (za(i1,i4)**2*za(i6,i5)*zb(i4,i1)*zb(i4,i2)*zb(i4,i3) + 
     &      2*za(i1,i5)**2*za(i6,i5)*zb(i4,i3)*zb(i5,i1)*zb(i5,i2) + 
     &      2*za(i6,i5)**2*zab2(i1,i5,i6,i2)*zb(i4,i3)*zb(i5,i6) + 
     &      2*za(i1,i5)*za(i1,i6)*za(i6,i5)*zb(i4,i3)*zb(i5,i1)*
     &       zb(i6,i2) + s(i1,i6)*
     &       (za(i6,i5)*(zb(i4,i3)*
     &             (za(i1,i4)*zb(i4,i2) + 2*za(i1,i5)*zb(i5,i2)) + 
     &            za(i1,i5)*zb(i4,i2)*zb(i5,i3)) + 
     &         za(i1,i6)*zb(i4,i3)*
     &          (za(i4,i5)*zb(i4,i2) + 4*za(i6,i5)*zb(i6,i2))) + 
     &      za(i1,i6)*za(i4,i5)*za(i4,i6)*zb(i4,i2)*zb(i4,i3)*
     &       zb(i6,i4) + za(i1,i5)*za(i4,i6)*za(i6,i5)*zb(i4,i2)*
     &       zb(i5,i3)*zb(i6,i4) + 
     &      2*za(i1,i6)*za(i4,i6)*za(i6,i5)*zb(i4,i3)*zb(i6,i2)*
     &       zb(i6,i4) + za(i1,i4)*
     &       (za(i1,i6)*zb(i4,i1)*zb(i4,i3)*
     &          (za(i4,i5)*zb(i4,i2) + 2*za(i6,i5)*zb(i6,i2)) + 
     &         za(i6,i5)*zb(i4,i2)*
     &          (za(i1,i5)*zb(i4,i1)*zb(i5,i3) + 
     &            za(i4,i6)*zb(i4,i3)*zb(i6,i4)))))/
     &  (za(i1,i6)*za(i6,i5)**2*zab2(i5,i1,i6,i4)**2*zb(i3,i2)*
     &    zb(i4,i3)**2) - (za(i4,i5)*zb(i5,i4)*
     &    (s(i2,i6)**2*zb(i4,i3) + 
     &      za(i2,i4)**2*zb(i4,i2)**2*zb(i4,i3) + 
     &      2*za(i2,i4)*za(i6,i5)*zb(i4,i3)*zb(i5,i2)*zb(i6,i4) + 
     &      za(i6,i5)*(2*za(i6,i5)*zb(i4,i3)*zb(i5,i6)**2 + 
     &         za(i4,i6)*zb(i6,i4)*
     &          (-2*zb(i4,i3)*zb(i5,i6) + zb(i5,i3)*zb(i6,i4)))))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i2,i6,i4)*zb(i3,i1)*zb(i4,i3)**2)
     &  - (za(i3,i5)*za(i4,i5)*zb(i5,i4)*
     &    (za(i1,i4)*za(i2,i5)*za(i6,i5)*zb(i4,i2)*zb(i4,i3)*
     &       zb(i6,i4) + za(i1,i5)*
     &       (-(za(i2,i6)*za(i4,i5)*zb(i4,i2)*zb(i4,i3)*zb(i6,i4)) + 
     &         za(i2,i5)*za(i6,i5)*
     &          (zb(i4,i2)*zb(i4,i3)*zb(i5,i6) + 
     &            (-(zb(i4,i3)*zb(i5,i2)) + zb(i3,i2)*zb(i5,i4))*
     &             zb(i6,i4)))))/
     &  (za(i2,i5)**2*za(i6,i5)**2*zab2(i5,i3,i6,i4)*zb(i4,i2)*
     &    zb(i4,i3)**2) - (za(i4,i5)*zb(i5,i4)*
     &    (2*za(i6,i5)**3*zb(i4,i3)*
     &       (s(i2,i6) - za(i2,i4)*zb(i4,i2) + 2*za(i2,i5)*zb(i5,i2))*
     &       zb(i5,i6) + za(i2,i5)*
     &       (za(i2,i6)*za(i4,i5)*zb(i4,i2)*zb(i4,i3)*
     &          (za(i4,i6)*zb(i4,i2) - 2*za(i6,i5)*zb(i5,i2)) + 
     &         za(i6,i5)*(zb(i4,i3)*
     &             (2*za(i6,i5)*zb(i5,i2)*
     &                (s(i2,i6) + za(i2,i5)*zb(i5,i2)) + 
     &               za(i2,i4)*zb(i4,i2)*
     &                (za(i4,i6)*zb(i4,i2) - 2*za(i6,i5)*zb(i5,i2))) + 
     &            za(i4,i6)*zb(i4,i2)*zb(i5,i3)*
     &             (za(i2,i5)*zb(i4,i2) - 2*za(i6,i5)*zb(i6,i4))))))/
     &  (za(i2,i6)*za(i6,i5)**3*zab2(i5,i2,i6,i4)*zb(i3,i1)*
     &    zb(i4,i3)**2) + (za(i4,i5)*zb(i5,i4)*
     &    (za(i2,i4)**2*za(i6,i5)*zb(i4,i2)**2*zb(i4,i3)*
     &       (za(i2,i5)*zb(i4,i2) - 2*za(i6,i5)*zb(i6,i4)) + 
     &      za(i2,i5)*za(i4,i6)*zb(i6,i4)*
     &       (za(i2,i6)*za(i4,i5)*zb(i4,i2)**2*zb(i4,i3) - 
     &         za(i6,i5)*(za(i2,i5)*zb(i4,i2)*
     &             (zb(i4,i3)*zb(i5,i2) - zb(i3,i2)*zb(i5,i4)) + 
     &            2*za(i6,i5)*
     &             (zb(i4,i2)*zb(i4,i3)*zb(i5,i6) + 
     &               zb(i3,i2)*zb(i5,i4)*zb(i6,i4)))) + 
     &      za(i2,i4)*zb(i4,i2)*
     &       (za(i2,i5)*za(i2,i6)*za(i4,i5)*zb(i4,i2)**2*zb(i4,i3) + 
     &         za(i6,i5)*(za(i2,i5)**2*zb(i4,i2)*
     &             (-(zb(i4,i3)*zb(i5,i2)) + zb(i3,i2)*zb(i5,i4)) - 
     &            3*za(i4,i6)*za(i6,i5)*zb(i4,i3)*zb(i6,i4)**2 + 
     &            za(i2,i5)*(za(i4,i6)*zb(i4,i2)*zb(i4,i3)*zb(i6,i4) - 
     &               2*za(i6,i5)*
     &                (zb(i4,i2)*zb(i4,i3)*zb(i5,i6) + 
     &                  zb(i3,i2)*zb(i5,i4)*zb(i6,i4))))) + 
     &      s(i2,i6)*(za(i2,i5)*za(i2,i6)*za(i4,i5)*zb(i4,i2)**2*
     &          zb(i4,i3) + za(i6,i5)*
     &          (za(i2,i4)*zb(i4,i2)*zb(i4,i3)*
     &             (za(i2,i5)*zb(i4,i2) - 2*za(i6,i5)*zb(i6,i4)) - 
     &            za(i2,i5)*(za(i2,i5)*zb(i4,i2)*
     &                (zb(i4,i3)*zb(i5,i2) - zb(i3,i2)*zb(i5,i4)) + 
     &               2*za(i6,i5)*
     &                (zb(i4,i2)*zb(i4,i3)*zb(i5,i6) + 
     &                  zb(i3,i2)*zb(i5,i4)*zb(i6,i4)))))))/
     &  (za(i2,i6)*za(i6,i5)**2*zab2(i5,i2,i6,i4)**2*zb(i3,i1)*
     &    zb(i4,i3)**2) - (za(i4,i5)*zb(i5,i4)*
     &    (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)**3*za(i6,i5)**3*zb(i3,i1)*
     &       zb(i4,i2)**2*zb(i4,i3)**3*
     &       (-(za(i1,i4)**2*zab2(i3,i2,i6,i4)*zb(i4,i1)**2) + 
     &         za(i1,i5)*zab2(i1,i4,i5,i1)*zb(i5,i4)*
     &          (za(i2,i3)*zb(i2,i1) - za(i3,i6)*zb(i6,i1)))*zb(i6,i2)
     &       + t(i2,i3,i6)*(za(i1,i3)*za(i1,i5)**2*za(i2,i6)*
     &          za(i6,i5)**3*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)**3*
     &          (-(za(i2,i4)*za(i4,i5)*zab2(i2,i1,i3,i6)*
     &               zb(i4,i2)**2) - 
     &            za(i2,i5)*za(i4,i5)*zab2(i2,i1,i3,i6)*zb(i4,i2)*
     &             zb(i5,i2) + 
     &            za(i2,i5)**2*zab2(i5,i1,i3,i6)*zb(i5,i2)**2) + 
     &         t(i2,i4,i5)*(-(za(i1,i5)**2*za(i1,i6)*za(i2,i5)*
     &               (-(za(i2,i6)**2*za(i4,i5)*za(i6,i5)**2*zb(i3,i2)*
     &                    zb(i4,i1)**3*zb(i4,i2)*zb(i4,i3)**2*zb(i6,i2)*
     &                    zb(i6,i4)) + 
     &                 za(i2,i5)**2*za(i6,i5)**2*zb(i4,i2)**2*
     &                  (-2*za(i2,i4)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i2)*
     &                     zb(i4,i3)*zb(i5,i3) + 
     &                    2*za(i2,i5)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)*
     &                     zb(i5,i2)*zb(i5,i3) + 
     &                    za(i2,i5)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)**2*
     &                     zb(i5,i3)**2 + 
     &                    za(i2,i5)*zb(i2,i1)**2*zb(i3,i1)*zb(i4,i3)**2*
     &                     zb(i5,i4)**2 + 
     &                    2*za(i6,i5)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)*
     &                     zb(i5,i3)*zb(i5,i6) + 
     &                    za(i6,i5)*zb(i2,i1)*zb(i3,i1)*zb(i4,i3)**2*
     &                     zb(i5,i4)**2*zb(i6,i1) - 
     &                    2*za(i4,i6)*zb(i3,i2)*zb(i4,i1)**3*zb(i4,i3)*
     &                     zb(i5,i3)*zb(i6,i4) - 
     &                    za(i6,i5)*zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)*
     &                     zb(i5,i3)**2*zb(i6,i4)) + 
     &                 za(i2,i5)*za(i2,i6)*zb(i4,i1)*zb(i4,i3)*
     &                  (za(i2,i5)*zb(i4,i2)*
     &                     (za(i4,i6)*zb(i4,i1)*zb(i4,i2)*zb(i4,i3)*
     &                        (za(i4,i5)*zb(i4,i2)*
     &                        (zb(i3,i2)*zb(i4,i1) + 
     &                       zb(i3,i1)*zb(i4,i2)) + 
     &                        za(i6,i5)*zb(i3,i2)*zb(i4,i1)*zb(i6,i2))
     &                        + za(i6,i5)*
     &                        (-(za(i4,i5)*zb(i2,i1)*zb(i4,i2)*
     &                        zb(i4,i3)*
     &                        (zb(i4,i1)*zb(i4,i3)*zb(i5,i2) + 
     &                        zb(i3,i1)*zb(i4,i2)*zb(i5,i4))) + 
     &                        za(i6,i5)*zb(i3,i2)*zb(i4,i1)**2*
     &                        (zb(i4,i3)*zb(i5,i2) - 
     &                       zb(i3,i2)*zb(i5,i4))*zb(i6,i2))) + 
     &                    za(i6,i5)*zb(i3,i2)*zb(i4,i1)**2*
     &                     (za(i6,i5)**2*zb(i4,i2)*zb(i4,i3)*zb(i5,i6)*
     &                        zb(i6,i2) + 
     &                       (-(zb(i4,i2)*zb(i4,i3)*
     &                        (-(za(i4,i6)*zab2(i5,i4,i6,i2)) + 
     &                        za(i4,i5)*za(i6,i5)*zb(i5,i2))) + 
     &                        za(i6,i5)**2*zb(i3,i2)*zb(i5,i4)*zb(i6,i2)
     &                        )*zb(i6,i4))))) + 
     &            za(i1,i6)**2*za(i2,i5)**3*za(i6,i5)*zb(i3,i1)*
     &             zb(i4,i1)*zb(i4,i2)**2*zb(i4,i3)*zb(i6,i2)*
     &             (za(i1,i4)*za(i6,i5)*zab2(i5,i2,i6,i4)*zb(i4,i1)*
     &                zb(i4,i3) + 
     &               za(i1,i5)*
     &                (za(i2,i6)*za(i4,i5)*zb(i4,i1)*zb(i4,i2)*
     &                   zb(i4,i3) - 
     &                  za(i6,i5)*
     &                   (za(i2,i5)*
     &                      (zb(i4,i1)*zb(i4,i2)*zb(i5,i3) + 
     &                        zb(i2,i1)*zb(i4,i3)*zb(i5,i4)) + 
     &                     za(i6,i5)*
     &                      (zb(i4,i1)*zb(i4,i3)*zb(i5,i6) - 
     &                        zb(i3,i1)*zb(i5,i4)*zb(i6,i4))))) + 
     &            za(i1,i5)**2*za(i6,i5)**2*zb(i3,i2)*zb(i4,i1)**3*
     &             (za(i1,i2)*za(i2,i5)**2*za(i6,i5)*zab2(i5,i2,i6,i4)*
     &                zb(i4,i2)**2*zb(i5,i3)**2 + 
     &               za(i1,i5)*za(i2,i6)*zb(i4,i3)*
     &                (za(i2,i4)*za(i2,i6)*za(i4,i5)*zb(i4,i2)**2*
     &                   zb(i4,i3)*zb(i6,i4) + 
     &                  za(i2,i5)*
     &                   (-(za(i2,i6)*za(i4,i5)*zb(i4,i2)*
     &                        (zb(i4,i2)*zb(i4,i3)*zb(i5,i6) + 
     &                        zb(i3,i2)*zb(i5,i4)*zb(i6,i4))) + 
     &                     za(i2,i5)*zb(i3,i2)*zb(i5,i4)*
     &                      (za(i6,i5)*zb(i4,i2)*zb(i5,i6) + 
     &                        zb(i5,i2)*
     &                        (2*za(i2,i5)*zb(i4,i2) - 
     &                        za(i6,i5)*zb(i6,i4))))))))))/
     &  (t(i2,i3,i6)*t(i2,i4,i5)*za(i1,i5)**2*za(i1,i6)*za(i2,i5)**3*
     &    za(i2,i6)*za(i6,i5)**3*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)**3*
     &    zb(i4,i2)**2*zb(i4,i3)**3)
     &)


      aaaa_NMHV_c(8) = (
     &  -((za(i1,i6)**2*za(i3,i6)*
     &      (2*s(i1,i5)**2 + zab2(i3,i1,i5,i3)**2 + 
     &        2*s(i1,i5)*(zab2(i1,i3,i6,i1) + zab2(i5,i3,i6,i5)) + 
     &        zab2(i6,i1,i5,i6)**2)*zb(i3,i2)**2*zb(i6,i3))/
     &    (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)**3*zb(i3,i4)*zb(i4,i2))
     &    ) - (za(i3,i6)*(za(i1,i3)**2*zb(i3,i2)**2 + 
     &      2*za(i1,i3)*za(i1,i5)*zb(i3,i2)*zb(i5,i2) + 
     &      2*za(i1,i5)**2*zb(i5,i2)**2)*zb(i6,i3))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)*zb(i3,i4)*zb(i4,i2)) - 
     & (za(i1,i6)*za(i3,i6)*za(i4,i6)*
     &    (2*s(i4,i5) + zab2(i4,i3,i6,i4) + zab2(i5,i3,i6,i5))*
     &    zb(i5,i3)*zb(i6,i3))/
     &  (za(i2,i6)*za(i5,i6)*zab2(i6,i4,i5,i3)**2*zb(i3,i4)) + 
     & (za(i3,i6)*za(i4,i6)*(2*s(i1,i2) + zab2(i1,i3,i6,i1) + 
     &      zab2(i2,i3,i6,i2))*zb(i3,i2)*zb(i5,i3)*zb(i6,i3))/
     &  (za(i5,i6)*zab2(i6,i1,i2,i3)**2*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i1,i6)*za(i3,i6)*(za(i4,i3)*zb(i3,i2) + 
     &      2*za(i4,i5)*zb(i5,i2))*zb(i5,i3)*zb(i6,i3))/
     &  (za(i2,i6)*za(i5,i6)*zab2(i6,i4,i5,i3)*zb(i3,i2)*zb(i3,i4)) - 
     & (za(i3,i6)*za(i5,i6)*(s(i2,i5) + zab2(i3,i2,i5,i3))**2*
     &    zb(i5,i3)**2*zb(i6,i3))/
     &  (za(i2,i5)*zab2(i6,i2,i5,i3)**3*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i2,i6)*za(i3,i6)*(s(i2,i5) + zab2(i3,i2,i5,i3))**2*zb(i3,i2)*
     &    (za(i2,i6)*zb(i3,i2) - 2*za(i5,i6)*zb(i5,i3))*zb(i6,i3))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i2,i5,i3)**3*zb(i3,i4)*zb(i4,i1))
     &  - (za(i3,i6)*zb(i6,i3)*
     &    (-(zb(i3,i1)*zb(i3,i2)*zb(i3,i4)*
     &         (-2*za(i4,i6)*za(i5,i6)*zab2(i1,i3,i4,i5) + 
     &           za(i1,i6)*za(i3,i6)*za(i4,i5)*zb(i5,i3))) + 
     &      za(i1,i6)*za(i4,i6)*za(i5,i6)*
     &       (zb(i2,i1)*zb(i3,i4)*zb(i5,i3) + 
     &         zb(i3,i1)*zb(i3,i2)*zb(i5,i4))*zb(i6,i3)))/
     &  (za(i1,i6)*za(i5,i6)**2*zab2(i6,i1,i2,i3)*zb(i3,i1)**2*
     &    zb(i3,i4)**2) - (za(i1,i6)*za(i3,i6)*zb(i6,i3)*
     &    (za(i5,i6)*(za(i1,i3)*zb(i3,i2)**2*zb(i3,i4)*
     &          (za(i3,i5)*zb(i3,i4) + za(i5,i6)*zb(i6,i4)) + 
     &         za(i1,i6)*za(i5,i6)*
     &          (zb(i3,i4)**2*zb(i6,i2)**2 + zb(i3,i2)**2*zb(i6,i4)**2))
     &        + za(i1,i5)*zb(i3,i4)*
     &       (za(i3,i5)*zb(i3,i2)*zb(i3,i4)*
     &          (za(i3,i6)*zb(i3,i2) + 2*za(i5,i6)*zb(i5,i2)) + 
     &         za(i5,i6)*(za(i3,i6)*zb(i3,i2)**2*zb(i6,i4) + 
     &            2*za(i5,i6)*zb(i5,i2)*
     &             (zb(i3,i4)*zb(i6,i2) + zb(i3,i2)*zb(i6,i4))))))/
     &  (za(i1,i5)*za(i5,i6)**3*zab2(i6,i1,i5,i3)*zb(i3,i4)**3*
     &    zb(i4,i2)) + (za(i1,i6)*za(i3,i6)*zb(i3,i2)*zb(i6,i3)*
     &    (za(i1,i3)**2*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i3,i4) + 
     &      za(i1,i5)*za(i3,i5)*za(i3,i6)*zb(i3,i2)*zb(i3,i4)*
     &       zb(i5,i3) + 2*za(i1,i5)*za(i3,i5)*za(i5,i6)*zb(i3,i4)*
     &       zb(i5,i2)*zb(i5,i3) + 
     &      2*za(i1,i5)*za(i1,i6)*za(i5,i6)*zb(i3,i4)*zb(i5,i2)*
     &       zb(i6,i1) + 2*za(i1,i6)**2*za(i5,i6)*zb(i3,i4)*zb(i6,i1)*
     &       zb(i6,i2) + za(i1,i6)*za(i3,i5)*za(i5,i6)*zb(i3,i2)*
     &       zb(i5,i3)*zb(i6,i4) + 
     &      za(i1,i3)*(za(i1,i5)*zb(i3,i1)*zb(i3,i4)*
     &          (za(i3,i6)*zb(i3,i2) + 2*za(i5,i6)*zb(i5,i2)) + 
     &         za(i5,i6)*zb(i3,i2)*
     &          (za(i3,i5)*zb(i3,i4)*zb(i5,i3) + 
     &            za(i1,i6)*zb(i3,i1)*zb(i6,i4))) + 
     &      s(i1,i5)*(za(i1,i5)*zb(i3,i4)*
     &          (za(i3,i6)*zb(i3,i2) + 4*za(i5,i6)*zb(i5,i2)) + 
     &         za(i5,i6)*(zb(i3,i4)*
     &             (za(i1,i3)*zb(i3,i2) + 2*za(i1,i6)*zb(i6,i2)) + 
     &            za(i1,i6)*zb(i3,i2)*zb(i6,i4))) + 
     &      2*za(i5,i6)**2*zab2(i1,i5,i6,i2)*zb(i3,i4)*zb(i6,i5)))/
     &  (za(i1,i5)*za(i5,i6)**2*zab2(i6,i1,i5,i3)**2*zb(i3,i4)**2*
     &    zb(i4,i2)) - (za(i3,i6)*zb(i6,i3)*
     &    (za(i2,i6)*(za(i2,i5)*za(i3,i6)*zb(i3,i2)*zb(i3,i4)*
     &          (za(i3,i5)*zb(i3,i2) - 2*za(i5,i6)*zb(i6,i2)) + 
     &         za(i5,i6)*(zb(i3,i4)*
     &             (2*za(i5,i6)*zb(i6,i2)*
     &                (s(i2,i5) + za(i2,i6)*zb(i6,i2)) + 
     &               za(i2,i3)*zb(i3,i2)*
     &                (za(i3,i5)*zb(i3,i2) - 2*za(i5,i6)*zb(i6,i2))) + 
     &            za(i3,i5)*zb(i3,i2)*
     &             (za(i2,i6)*zb(i3,i2) - 2*za(i5,i6)*zb(i5,i3))*
     &             zb(i6,i4))) + 
     &      2*za(i5,i6)**3*zb(i3,i4)*
     &       (s(i2,i5) - za(i2,i3)*zb(i3,i2) + 2*za(i2,i6)*zb(i6,i2))*
     &       zb(i6,i5)))/
     &  (za(i2,i5)*za(i5,i6)**3*zab2(i6,i2,i5,i3)*zb(i3,i4)**2*
     &    zb(i4,i1)) - (za(i3,i6)*zb(i6,i3)*
     &    (s(i2,i5)**2*zb(i3,i4) + 
     &      za(i2,i3)**2*zb(i3,i2)**2*zb(i3,i4) + 
     &      2*za(i2,i3)*za(i5,i6)*zb(i3,i4)*zb(i5,i3)*zb(i6,i2) + 
     &      za(i5,i6)*(2*za(i5,i6)*zb(i3,i4)*zb(i6,i5)**2 + 
     &         za(i3,i5)*zb(i5,i3)*
     &          (zb(i5,i3)*zb(i6,i4) - 2*zb(i3,i4)*zb(i6,i5)))))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i2,i5,i3)*zb(i3,i4)**2*zb(i4,i1))
     &  - (za(i3,i6)*zb(i5,i3)*zb(i6,i3)*
     &    (za(i3,i5)**2*zb(i3,i4)*zb(i5,i3)**2 - 
     &      za(i2,i3)*za(i5,i6)*zb(i3,i2)*
     &       (zb(i5,i4)*zb(i6,i3) + zb(i3,i4)*zb(i6,i5)) - 
     &      za(i3,i5)*za(i5,i6)*zb(i5,i3)*
     &       (zb(i5,i4)*zb(i6,i3) + zb(i3,i4)*zb(i6,i5)) + 
     &      s(i2,i5)*(za(i3,i5)*zb(i3,i4)*zb(i5,i3) - 
     &         za(i5,i6)*(zb(i5,i4)*zb(i6,i3) + zb(i3,i4)*zb(i6,i5)))))/
     &  (za(i2,i5)*zab2(i6,i2,i5,i3)**2*zb(i3,i4)**2*zb(i4,i1)) - 
     & (za(i3,i6)*za(i4,i6)*zb(i6,i3)*
     &    (za(i1,i3)*za(i2,i6)*za(i5,i6)*zb(i3,i2)*zb(i3,i4)*
     &       zb(i5,i3) + za(i1,i6)*
     &       (-(za(i2,i5)*za(i3,i6)*zb(i3,i2)*zb(i3,i4)*zb(i5,i3)) + 
     &         za(i2,i6)*za(i5,i6)*
     &          (zb(i5,i3)*(-(zb(i3,i4)*zb(i6,i2)) + 
     &               zb(i4,i2)*zb(i6,i3)) + 
     &            zb(i3,i2)*zb(i3,i4)*zb(i6,i5)))))/
     &  (za(i2,i6)**2*za(i5,i6)**2*zab2(i6,i4,i5,i3)*zb(i3,i2)*
     &    zb(i3,i4)**2) + (za(i3,i6)*zb(i6,i3)*
     &    (za(i2,i3)**2*za(i5,i6)*zb(i3,i2)**2*zb(i3,i4)*
     &       (za(i2,i6)*zb(i3,i2) - 2*za(i5,i6)*zb(i5,i3)) + 
     &      za(i2,i6)*za(i3,i5)*zb(i5,i3)*
     &       (za(i2,i5)*za(i3,i6)*zb(i3,i2)**2*zb(i3,i4) - 
     &         za(i5,i6)*(za(i2,i6)*zb(i3,i2)*
     &             (zb(i3,i4)*zb(i6,i2) - zb(i4,i2)*zb(i6,i3)) + 
     &            2*za(i5,i6)*
     &             (zb(i4,i2)*zb(i5,i3)*zb(i6,i3) + 
     &               zb(i3,i2)*zb(i3,i4)*zb(i6,i5)))) + 
     &      za(i2,i3)*zb(i3,i2)*
     &       (za(i2,i5)*za(i2,i6)*za(i3,i6)*zb(i3,i2)**2*zb(i3,i4) + 
     &         za(i5,i6)*(-3*za(i3,i5)*za(i5,i6)*zb(i3,i4)*
     &             zb(i5,i3)**2 + 
     &            za(i2,i6)**2*zb(i3,i2)*
     &             (-(zb(i3,i4)*zb(i6,i2)) + zb(i4,i2)*zb(i6,i3)) + 
     &            za(i2,i6)*(za(i3,i5)*zb(i3,i2)*zb(i3,i4)*zb(i5,i3) - 
     &               2*za(i5,i6)*
     &                (zb(i4,i2)*zb(i5,i3)*zb(i6,i3) + 
     &                  zb(i3,i2)*zb(i3,i4)*zb(i6,i5))))) + 
     &      s(i2,i5)*(za(i2,i5)*za(i2,i6)*za(i3,i6)*zb(i3,i2)**2*
     &          zb(i3,i4) + za(i5,i6)*
     &          (za(i2,i3)*zb(i3,i2)*zb(i3,i4)*
     &             (za(i2,i6)*zb(i3,i2) - 2*za(i5,i6)*zb(i5,i3)) - 
     &            za(i2,i6)*(za(i2,i6)*zb(i3,i2)*
     &                (zb(i3,i4)*zb(i6,i2) - zb(i4,i2)*zb(i6,i3)) + 
     &               2*za(i5,i6)*
     &                (zb(i4,i2)*zb(i5,i3)*zb(i6,i3) + 
     &                  zb(i3,i2)*zb(i3,i4)*zb(i6,i5)))))))/
     &  (za(i2,i5)*za(i5,i6)**2*zab2(i6,i2,i5,i3)**2*zb(i3,i4)**2*
     &    zb(i4,i1)) - (za(i3,i6)*zb(i6,i3)*
     &    (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)**3*za(i5,i6)**3*zb(i3,i2)**2*
     &       zb(i3,i4)**3*zb(i4,i1)*zb(i5,i2)*
     &       (-(za(i1,i3)**2*zab2(i4,i2,i5,i3)*zb(i3,i1)**2) + 
     &         za(i1,i6)*zab2(i1,i3,i6,i1)*
     &          (za(i2,i4)*zb(i2,i1) - za(i4,i5)*zb(i5,i1))*zb(i6,i3))
     &       + t(i2,i4,i5)*(za(i1,i4)*za(i1,i6)**2*za(i2,i5)*
     &          za(i5,i6)**3*zb(i3,i1)**3*zb(i3,i4)**3*zb(i4,i2)*
     &          (-(za(i2,i3)*za(i3,i6)*zab2(i2,i1,i4,i5)*
     &               zb(i3,i2)**2) - 
     &            za(i2,i6)*za(i3,i6)*zab2(i2,i1,i4,i5)*zb(i3,i2)*
     &             zb(i6,i2) + 
     &            za(i2,i6)**2*zab2(i6,i1,i4,i5)*zb(i6,i2)**2) + 
     &         t(i2,i3,i6)*(-(za(i1,i5)*za(i1,i6)**2*za(i2,i6)*
     &               (-(za(i2,i5)**2*za(i3,i6)*za(i5,i6)**2*
     &                    zb(i3,i1)**3*zb(i3,i2)*zb(i3,i4)**2*zb(i4,i2)*
     &                    zb(i5,i2)*zb(i5,i3)) + 
     &                 za(i2,i6)**2*za(i5,i6)**2*zb(i3,i2)**2*
     &                  (za(i2,i6)*zb(i2,i1)**2*zb(i3,i4)**2*zb(i4,i1)*
     &                     zb(i6,i3)**2 + 
     &                    za(i5,i6)*zb(i2,i1)*zb(i3,i4)**2*zb(i4,i1)*
     &                     zb(i5,i1)*zb(i6,i3)**2 - 
     &                    2*za(i2,i3)*zb(i3,i1)**3*zb(i3,i2)*zb(i3,i4)*
     &                     zb(i4,i2)*zb(i6,i4) - 
     &                    2*za(i3,i5)*zb(i3,i1)**3*zb(i3,i4)*zb(i4,i2)*
     &                     zb(i5,i3)*zb(i6,i4) + 
     &                    2*za(i2,i6)*zb(i3,i1)**3*zb(i3,i4)*zb(i4,i2)*
     &                     zb(i6,i2)*zb(i6,i4) + 
     &                    za(i2,i6)*zb(i3,i1)**2*zb(i3,i2)**2*zb(i4,i1)*
     &                     zb(i6,i4)**2 - 
     &                    za(i5,i6)*zb(i3,i1)**2*zb(i3,i2)*zb(i4,i1)*
     &                     zb(i5,i3)*zb(i6,i4)**2 + 
     &                    2*za(i5,i6)*zb(i3,i1)**3*zb(i3,i4)*zb(i4,i2)*
     &                     zb(i6,i4)*zb(i6,i5)) + 
     &                 za(i2,i5)*za(i2,i6)*zb(i3,i1)*zb(i3,i4)*
     &                  (za(i2,i6)*zb(i3,i2)*
     &                     (za(i3,i5)*zb(i3,i1)*zb(i3,i2)*zb(i3,i4)*
     &                        (za(i3,i6)*zb(i3,i2)*
     &                        (zb(i3,i2)*zb(i4,i1) + 
     &                       zb(i3,i1)*zb(i4,i2)) + 
     &                        za(i5,i6)*zb(i3,i1)*zb(i4,i2)*zb(i5,i2))
     &                        + za(i5,i6)*
     &                        (-(za(i3,i6)*zb(i2,i1)*zb(i3,i2)*
     &                        zb(i3,i4)*
     &                        (zb(i3,i1)*zb(i3,i4)*zb(i6,i2) + 
     &                        zb(i3,i2)*zb(i4,i1)*zb(i6,i3))) + 
     &                        za(i5,i6)*zb(i3,i1)**2*zb(i4,i2)*
     &                        zb(i5,i2)*
     &                        (zb(i3,i4)*zb(i6,i2) - 
     &                       zb(i4,i2)*zb(i6,i3)))) + 
     &                    za(i5,i6)*zb(i3,i1)**2*zb(i4,i2)*
     &                     (zb(i5,i3)*
     &                        (-(zb(i3,i2)*zb(i3,i4)*
     &                        (-(za(i3,i5)*zab2(i6,i3,i5,i2)) + 
     &                        za(i3,i6)*za(i5,i6)*zb(i6,i2))) + 
     &                        za(i5,i6)**2*zb(i4,i2)*zb(i5,i2)*zb(i6,i3)
     &                        ) + 
     &                       za(i5,i6)**2*zb(i3,i2)*zb(i3,i4)*zb(i5,i2)*
     &                        zb(i6,i5))))) + 
     &            za(i1,i5)**2*za(i2,i6)**3*za(i5,i6)*zb(i3,i1)*
     &             zb(i3,i2)**2*zb(i3,i4)*zb(i4,i1)*zb(i5,i2)*
     &             (za(i1,i3)*za(i5,i6)*zab2(i6,i2,i5,i3)*zb(i3,i1)*
     &                zb(i3,i4) + 
     &               za(i1,i6)*
     &                (za(i2,i5)*za(i3,i6)*zb(i3,i1)*zb(i3,i2)*
     &                   zb(i3,i4) - 
     &                  za(i5,i6)*
     &                   (za(i2,i6)*
     &                      (zb(i2,i1)*zb(i3,i4)*zb(i6,i3) + 
     &                        zb(i3,i1)*zb(i3,i2)*zb(i6,i4)) + 
     &                     za(i5,i6)*
     &                      (-(zb(i4,i1)*zb(i5,i3)*zb(i6,i3)) + 
     &                        zb(i3,i1)*zb(i3,i4)*zb(i6,i5))))) + 
     &            za(i1,i6)**2*za(i5,i6)**2*zb(i3,i1)**3*zb(i4,i2)*
     &             (za(i1,i2)*za(i2,i6)**2*za(i5,i6)*zab2(i6,i2,i5,i3)*
     &                zb(i3,i2)**2*zb(i6,i4)**2 + 
     &               za(i1,i6)*za(i2,i5)*zb(i3,i4)*
     &                (za(i2,i3)*za(i2,i5)*za(i3,i6)*zb(i3,i2)**2*
     &                   zb(i3,i4)*zb(i5,i3) + 
     &                  za(i2,i6)*
     &                   (za(i2,i6)*zb(i4,i2)*zb(i6,i3)*
     &                      ((2*za(i2,i6)*zb(i3,i2) - 
     &                        za(i5,i6)*zb(i5,i3))*zb(i6,i2) + 
     &                        za(i5,i6)*zb(i3,i2)*zb(i6,i5)) - 
     &                     za(i2,i5)*za(i3,i6)*zb(i3,i2)*
     &                      (zb(i4,i2)*zb(i5,i3)*zb(i6,i3) + 
     &                        zb(i3,i2)*zb(i3,i4)*zb(i6,i5)))))))))/
     &  (t(i2,i3,i6)*t(i2,i4,i5)*za(i1,i5)*za(i1,i6)**2*za(i2,i5)*
     &    za(i2,i6)**3*za(i5,i6)**3*zb(i3,i1)**3*zb(i3,i2)**2*
     &    zb(i3,i4)**3*zb(i4,i1)*zb(i4,i2))
     &)


      aaaa_NMHV_c(9) = (
     &  -((za(i1,i3)*(-(t(i2,i4,i6)*za(i1,i3)*za(i1,i5)*za(i2,i6)*
     &          zab2(i4,i1,i3,i6)*zb(i3,i4)*zb(i5,i2)) + 
     &       t(i2,i4,i5)*(t(i2,i4,i6)*t(i2,i5,i6)*za(i1,i2)*
     &           zab2(i1,i3,i4,i2) - 
     &          za(i1,i3)*za(i1,i6)*za(i2,i5)*zab2(i4,i1,i3,i5)*
     &           zb(i3,i4)*zb(i6,i2))))/
     &   (t(i2,i4,i5)*t(i2,i4,i6)*za(i1,i5)*za(i1,i6)*za(i2,i5)*
     &     za(i2,i6)*zb(i3,i4)*zb(i4,i2)))
     &)


      aaaa_NMHV_c(10) = (
     &  -((za(i1,i5)**2*za(i3,i5)*
     &      (2*s(i1,i6)**2 + zab2(i3,i1,i6,i3)**2 + 
     &        zab2(i5,i1,i6,i5)**2 + 
     &        2*s(i1,i6)*(zab2(i1,i3,i5,i1) + zab2(i6,i3,i5,i6)))*
     &      zb(i3,i2)**2*zb(i5,i3))/
     &    (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i3)**3*zb(i3,i4)*zb(i4,i2))
     &    ) - (za(i3,i5)*zb(i5,i3)*
     &    (za(i1,i3)**2*zb(i3,i2)**2 + 
     &      2*za(i1,i3)*za(i1,i6)*zb(i3,i2)*zb(i6,i2) + 
     &      2*za(i1,i6)**2*zb(i6,i2)**2))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i3)*zb(i3,i4)*zb(i4,i2)) - 
     & (za(i1,i5)*za(i3,i5)*zb(i5,i3)*
     &    (za(i6,i5)*(za(i1,i3)*zb(i3,i2)**2*zb(i3,i4)*
     &          (za(i3,i6)*zb(i3,i4) + za(i6,i5)*zb(i5,i4)) + 
     &         za(i1,i5)*za(i6,i5)*
     &          (zb(i3,i4)**2*zb(i5,i2)**2 + zb(i3,i2)**2*zb(i5,i4)**2))
     &        + za(i1,i6)*zb(i3,i4)*
     &       (za(i3,i6)*zb(i3,i2)*zb(i3,i4)*
     &          (za(i3,i5)*zb(i3,i2) + 2*za(i6,i5)*zb(i6,i2)) + 
     &         za(i6,i5)*(za(i3,i5)*zb(i3,i2)**2*zb(i5,i4) + 
     &            2*za(i6,i5)*
     &             (zb(i3,i4)*zb(i5,i2) + zb(i3,i2)*zb(i5,i4))*zb(i6,i2)
     &            ))))/
     &  (za(i1,i6)*za(i6,i5)**3*zab2(i5,i1,i6,i3)*zb(i3,i4)**3*
     &    zb(i4,i2)) - (za(i1,i5)*za(i3,i5)*za(i4,i5)*
     &    (2*s(i4,i6) + zab2(i4,i3,i5,i4) + zab2(i6,i3,i5,i6))*
     &    zb(i5,i3)*zb(i6,i3))/
     &  (za(i2,i5)*za(i6,i5)*zab2(i5,i4,i6,i3)**2*zb(i3,i4)) + 
     & (za(i3,i5)*za(i4,i5)*(2*s(i1,i2) + zab2(i1,i3,i5,i1) + 
     &      zab2(i2,i3,i5,i2))*zb(i3,i2)*zb(i5,i3)*zb(i6,i3))/
     &  (za(i6,i5)*zab2(i5,i1,i2,i3)**2*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i1,i5)*za(i3,i5)*zb(i5,i3)*
     &    (za(i4,i3)*zb(i3,i2) + 2*za(i4,i6)*zb(i6,i2))*zb(i6,i3))/
     &  (za(i2,i5)*za(i6,i5)*zab2(i5,i4,i6,i3)*zb(i3,i2)*zb(i3,i4)) - 
     & (za(i3,i5)*za(i6,i5)*(s(i2,i6) + zab2(i3,i2,i6,i3))**2*zb(i5,i3)*
     &    zb(i6,i3)**2)/
     &  (za(i2,i6)*zab2(i5,i2,i6,i3)**3*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i2,i5)*za(i3,i5)*(s(i2,i6) + zab2(i3,i2,i6,i3))**2*zb(i3,i2)*
     &    zb(i5,i3)*(za(i2,i5)*zb(i3,i2) - 2*za(i6,i5)*zb(i6,i3)))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i2,i6,i3)**3*zb(i3,i4)*zb(i4,i1))
     &  + (za(i1,i5)*za(i3,i5)*zb(i3,i2)*zb(i5,i3)*
     &    (za(i1,i3)**2*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i3,i4) + 
     &      2*za(i1,i5)**2*za(i6,i5)*zb(i3,i4)*zb(i5,i1)*zb(i5,i2) + 
     &      2*za(i6,i5)**2*zab2(i1,i5,i6,i2)*zb(i3,i4)*zb(i5,i6) + 
     &      2*za(i1,i5)*za(i1,i6)*za(i6,i5)*zb(i3,i4)*zb(i5,i1)*
     &       zb(i6,i2) + s(i1,i6)*
     &       (za(i6,i5)*(zb(i3,i4)*
     &             (za(i1,i3)*zb(i3,i2) + 2*za(i1,i5)*zb(i5,i2)) + 
     &            za(i1,i5)*zb(i3,i2)*zb(i5,i4)) + 
     &         za(i1,i6)*zb(i3,i4)*
     &          (za(i3,i5)*zb(i3,i2) + 4*za(i6,i5)*zb(i6,i2))) + 
     &      za(i1,i6)*za(i3,i5)*za(i3,i6)*zb(i3,i2)*zb(i3,i4)*
     &       zb(i6,i3) + za(i1,i5)*za(i3,i6)*za(i6,i5)*zb(i3,i2)*
     &       zb(i5,i4)*zb(i6,i3) + 
     &      2*za(i1,i6)*za(i3,i6)*za(i6,i5)*zb(i3,i4)*zb(i6,i2)*
     &       zb(i6,i3) + za(i1,i3)*
     &       (za(i1,i6)*zb(i3,i1)*zb(i3,i4)*
     &          (za(i3,i5)*zb(i3,i2) + 2*za(i6,i5)*zb(i6,i2)) + 
     &         za(i6,i5)*zb(i3,i2)*
     &          (za(i1,i5)*zb(i3,i1)*zb(i5,i4) + 
     &            za(i3,i6)*zb(i3,i4)*zb(i6,i3)))))/
     &  (za(i1,i6)*za(i6,i5)**2*zab2(i5,i1,i6,i3)**2*zb(i3,i4)**2*
     &    zb(i4,i2)) - (za(i3,i5)*za(i4,i5)*zb(i5,i3)*
     &    (za(i1,i3)*za(i2,i5)*za(i6,i5)*zb(i3,i2)*zb(i3,i4)*
     &       zb(i6,i3) + za(i1,i5)*
     &       (-(za(i2,i6)*za(i3,i5)*zb(i3,i2)*zb(i3,i4)*zb(i6,i3)) + 
     &         za(i2,i5)*za(i6,i5)*
     &          (zb(i3,i2)*zb(i3,i4)*zb(i5,i6) + 
     &            (-(zb(i3,i4)*zb(i5,i2)) + zb(i4,i2)*zb(i5,i3))*
     &             zb(i6,i3)))))/
     &  (za(i2,i5)**2*za(i6,i5)**2*zab2(i5,i4,i6,i3)*zb(i3,i2)*
     &    zb(i3,i4)**2) - (za(i3,i5)*zb(i5,i3)*
     &    (s(i2,i6)**2*zb(i3,i4) + 
     &      za(i2,i3)**2*zb(i3,i2)**2*zb(i3,i4) + 
     &      2*za(i2,i3)*za(i6,i5)*zb(i3,i4)*zb(i5,i2)*zb(i6,i3) + 
     &      za(i6,i5)*(2*za(i6,i5)*zb(i3,i4)*zb(i5,i6)**2 + 
     &         za(i3,i6)*zb(i6,i3)*
     &          (-2*zb(i3,i4)*zb(i5,i6) + zb(i5,i4)*zb(i6,i3)))))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i2,i6,i3)*zb(i3,i4)**2*zb(i4,i1))
     &  - (za(i3,i5)*zb(i5,i3)*
     &    (2*za(i6,i5)**3*zb(i3,i4)*
     &       (s(i2,i6) - za(i2,i3)*zb(i3,i2) + 2*za(i2,i5)*zb(i5,i2))*
     &       zb(i5,i6) + za(i2,i5)*
     &       (za(i2,i6)*za(i3,i5)*zb(i3,i2)*zb(i3,i4)*
     &          (za(i3,i6)*zb(i3,i2) - 2*za(i6,i5)*zb(i5,i2)) + 
     &         za(i6,i5)*(zb(i3,i4)*
     &             (2*za(i6,i5)*zb(i5,i2)*
     &                (s(i2,i6) + za(i2,i5)*zb(i5,i2)) + 
     &               za(i2,i3)*zb(i3,i2)*
     &                (za(i3,i6)*zb(i3,i2) - 2*za(i6,i5)*zb(i5,i2))) + 
     &            za(i3,i6)*zb(i3,i2)*zb(i5,i4)*
     &             (za(i2,i5)*zb(i3,i2) - 2*za(i6,i5)*zb(i6,i3))))))/
     &  (za(i2,i6)*za(i6,i5)**3*zab2(i5,i2,i6,i3)*zb(i3,i4)**2*
     &    zb(i4,i1)) + (za(i3,i5)*zb(i5,i3)*
     &    (za(i2,i3)**2*za(i6,i5)*zb(i3,i2)**2*zb(i3,i4)*
     &       (za(i2,i5)*zb(i3,i2) - 2*za(i6,i5)*zb(i6,i3)) + 
     &      za(i2,i5)*za(i3,i6)*zb(i6,i3)*
     &       (za(i2,i6)*za(i3,i5)*zb(i3,i2)**2*zb(i3,i4) - 
     &         za(i6,i5)*(za(i2,i5)*zb(i3,i2)*
     &             (zb(i3,i4)*zb(i5,i2) - zb(i4,i2)*zb(i5,i3)) + 
     &            2*za(i6,i5)*
     &             (zb(i3,i2)*zb(i3,i4)*zb(i5,i6) + 
     &               zb(i4,i2)*zb(i5,i3)*zb(i6,i3)))) + 
     &      za(i2,i3)*zb(i3,i2)*
     &       (za(i2,i5)*za(i2,i6)*za(i3,i5)*zb(i3,i2)**2*zb(i3,i4) + 
     &         za(i6,i5)*(za(i2,i5)**2*zb(i3,i2)*
     &             (-(zb(i3,i4)*zb(i5,i2)) + zb(i4,i2)*zb(i5,i3)) - 
     &            3*za(i3,i6)*za(i6,i5)*zb(i3,i4)*zb(i6,i3)**2 + 
     &            za(i2,i5)*(za(i3,i6)*zb(i3,i2)*zb(i3,i4)*zb(i6,i3) - 
     &               2*za(i6,i5)*
     &                (zb(i3,i2)*zb(i3,i4)*zb(i5,i6) + 
     &                  zb(i4,i2)*zb(i5,i3)*zb(i6,i3))))) + 
     &      s(i2,i6)*(za(i2,i5)*za(i2,i6)*za(i3,i5)*zb(i3,i2)**2*
     &          zb(i3,i4) + za(i6,i5)*
     &          (za(i2,i3)*zb(i3,i2)*zb(i3,i4)*
     &             (za(i2,i5)*zb(i3,i2) - 2*za(i6,i5)*zb(i6,i3)) - 
     &            za(i2,i5)*(za(i2,i5)*zb(i3,i2)*
     &                (zb(i3,i4)*zb(i5,i2) - zb(i4,i2)*zb(i5,i3)) + 
     &               2*za(i6,i5)*
     &                (zb(i3,i2)*zb(i3,i4)*zb(i5,i6) + 
     &                  zb(i4,i2)*zb(i5,i3)*zb(i6,i3)))))))/
     &  (za(i2,i6)*za(i6,i5)**2*zab2(i5,i2,i6,i3)**2*zb(i3,i4)**2*
     &    zb(i4,i1)) - (za(i3,i5)*zb(i5,i3)*
     &    (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)**3*za(i6,i5)**3*zb(i3,i2)**2*
     &       zb(i3,i4)**3*zb(i4,i1)*
     &       (-(za(i1,i3)**2*zab2(i4,i2,i6,i3)*zb(i3,i1)**2) + 
     &         za(i1,i5)*zab2(i1,i3,i5,i1)*zb(i5,i3)*
     &          (za(i2,i4)*zb(i2,i1) - za(i4,i6)*zb(i6,i1)))*zb(i6,i2)
     &       + t(i2,i4,i6)*(za(i1,i4)*za(i1,i5)**2*za(i2,i6)*
     &          za(i6,i5)**3*zb(i3,i1)**3*zb(i3,i4)**3*zb(i4,i2)*
     &          (-(za(i2,i3)*za(i3,i5)*zab2(i2,i1,i4,i6)*
     &               zb(i3,i2)**2) - 
     &            za(i2,i5)*za(i3,i5)*zab2(i2,i1,i4,i6)*zb(i3,i2)*
     &             zb(i5,i2) + 
     &            za(i2,i5)**2*zab2(i5,i1,i4,i6)*zb(i5,i2)**2) + 
     &         t(i2,i3,i5)*(-(za(i1,i5)**2*za(i1,i6)*za(i2,i5)*
     &               (-(za(i2,i6)**2*za(i3,i5)*za(i6,i5)**2*
     &                    zb(i3,i1)**3*zb(i3,i2)*zb(i3,i4)**2*zb(i4,i2)*
     &                    zb(i6,i2)*zb(i6,i3)) + 
     &                 za(i2,i5)**2*za(i6,i5)**2*zb(i3,i2)**2*
     &                  (za(i2,i5)*zb(i2,i1)**2*zb(i3,i4)**2*zb(i4,i1)*
     &                     zb(i5,i3)**2 - 
     &                    2*za(i2,i3)*zb(i3,i1)**3*zb(i3,i2)*zb(i3,i4)*
     &                     zb(i4,i2)*zb(i5,i4) + 
     &                    2*za(i2,i5)*zb(i3,i1)**3*zb(i3,i4)*zb(i4,i2)*
     &                     zb(i5,i2)*zb(i5,i4) + 
     &                    za(i2,i5)*zb(i3,i1)**2*zb(i3,i2)**2*zb(i4,i1)*
     &                     zb(i5,i4)**2 + 
     &                    2*za(i6,i5)*zb(i3,i1)**3*zb(i3,i4)*zb(i4,i2)*
     &                     zb(i5,i4)*zb(i5,i6) + 
     &                    za(i6,i5)*zb(i2,i1)*zb(i3,i4)**2*zb(i4,i1)*
     &                     zb(i5,i3)**2*zb(i6,i1) - 
     &                    2*za(i3,i6)*zb(i3,i1)**3*zb(i3,i4)*zb(i4,i2)*
     &                     zb(i5,i4)*zb(i6,i3) - 
     &                    za(i6,i5)*zb(i3,i1)**2*zb(i3,i2)*zb(i4,i1)*
     &                     zb(i5,i4)**2*zb(i6,i3)) + 
     &                 za(i2,i5)*za(i2,i6)*zb(i3,i1)*zb(i3,i4)*
     &                  (za(i2,i5)*zb(i3,i2)*
     &                     (za(i3,i6)*zb(i3,i1)*zb(i3,i2)*zb(i3,i4)*
     &                        (za(i3,i5)*zb(i3,i2)*
     &                        (zb(i3,i2)*zb(i4,i1) + 
     &                       zb(i3,i1)*zb(i4,i2)) + 
     &                        za(i6,i5)*zb(i3,i1)*zb(i4,i2)*zb(i6,i2))
     &                        + za(i6,i5)*
     &                        (-(za(i3,i5)*zb(i2,i1)*zb(i3,i2)*
     &                        zb(i3,i4)*
     &                        (zb(i3,i1)*zb(i3,i4)*zb(i5,i2) + 
     &                        zb(i3,i2)*zb(i4,i1)*zb(i5,i3))) + 
     &                        za(i6,i5)*zb(i3,i1)**2*zb(i4,i2)*
     &                        (zb(i3,i4)*zb(i5,i2) - 
     &                       zb(i4,i2)*zb(i5,i3))*zb(i6,i2))) + 
     &                    za(i6,i5)*zb(i3,i1)**2*zb(i4,i2)*
     &                     (za(i6,i5)**2*zb(i3,i2)*zb(i3,i4)*zb(i5,i6)*
     &                        zb(i6,i2) + 
     &                       (-(zb(i3,i2)*zb(i3,i4)*
     &                        (-(za(i3,i6)*zab2(i5,i3,i6,i2)) + 
     &                        za(i3,i5)*za(i6,i5)*zb(i5,i2))) + 
     &                        za(i6,i5)**2*zb(i4,i2)*zb(i5,i3)*zb(i6,i2)
     &                        )*zb(i6,i3))))) + 
     &            za(i1,i6)**2*za(i2,i5)**3*za(i6,i5)*zb(i3,i1)*
     &             zb(i3,i2)**2*zb(i3,i4)*zb(i4,i1)*zb(i6,i2)*
     &             (za(i1,i3)*za(i6,i5)*zab2(i5,i2,i6,i3)*zb(i3,i1)*
     &                zb(i3,i4) + 
     &               za(i1,i5)*
     &                (za(i2,i6)*za(i3,i5)*zb(i3,i1)*zb(i3,i2)*
     &                   zb(i3,i4) - 
     &                  za(i6,i5)*
     &                   (za(i2,i5)*
     &                      (zb(i2,i1)*zb(i3,i4)*zb(i5,i3) + 
     &                        zb(i3,i1)*zb(i3,i2)*zb(i5,i4)) + 
     &                     za(i6,i5)*
     &                      (zb(i3,i1)*zb(i3,i4)*zb(i5,i6) - 
     &                        zb(i4,i1)*zb(i5,i3)*zb(i6,i3))))) + 
     &            za(i1,i5)**2*za(i6,i5)**2*zb(i3,i1)**3*zb(i4,i2)*
     &             (za(i1,i2)*za(i2,i5)**2*za(i6,i5)*zab2(i5,i2,i6,i3)*
     &                zb(i3,i2)**2*zb(i5,i4)**2 + 
     &               za(i1,i5)*za(i2,i6)*zb(i3,i4)*
     &                (za(i2,i3)*za(i2,i6)*za(i3,i5)*zb(i3,i2)**2*
     &                   zb(i3,i4)*zb(i6,i3) + 
     &                  za(i2,i5)*
     &                   (-(za(i2,i6)*za(i3,i5)*zb(i3,i2)*
     &                        (zb(i3,i2)*zb(i3,i4)*zb(i5,i6) + 
     &                        zb(i4,i2)*zb(i5,i3)*zb(i6,i3))) + 
     &                     za(i2,i5)*zb(i4,i2)*zb(i5,i3)*
     &                      (za(i6,i5)*zb(i3,i2)*zb(i5,i6) + 
     &                        zb(i5,i2)*
     &                        (2*za(i2,i5)*zb(i3,i2) - 
     &                        za(i6,i5)*zb(i6,i3))))))))))/
     &  (t(i2,i3,i5)*t(i2,i4,i6)*za(i1,i5)**2*za(i1,i6)*za(i2,i5)**3*
     &    za(i2,i6)*za(i6,i5)**3*zb(i3,i1)**3*zb(i3,i2)**2*zb(i3,i4)**3*
     &    zb(i4,i1)*zb(i4,i2)) - 
     & (za(i3,i5)*zb(i5,i3)*(-(zb(i3,i1)*zb(i3,i2)*zb(i3,i4)*
     &         (-2*za(i4,i5)*za(i6,i5)*zab2(i1,i3,i4,i6) + 
     &           za(i1,i5)*za(i3,i5)*za(i4,i6)*zb(i6,i3))) + 
     &      za(i1,i5)*za(i4,i5)*za(i6,i5)*zb(i5,i3)*
     &       (zb(i2,i1)*zb(i3,i4)*zb(i6,i3) + 
     &         zb(i3,i1)*zb(i3,i2)*zb(i6,i4))))/
     &  (za(i1,i5)*za(i6,i5)**2*zab2(i5,i1,i2,i3)*zb(i3,i1)**2*
     &    zb(i3,i4)**2) - (za(i3,i5)*zb(i5,i3)*zb(i6,i3)*
     &    (za(i3,i6)**2*zb(i3,i4)*zb(i6,i3)**2 - 
     &      za(i2,i3)*za(i6,i5)*zb(i3,i2)*
     &       (zb(i3,i4)*zb(i5,i6) + zb(i5,i3)*zb(i6,i4)) - 
     &      za(i3,i6)*za(i6,i5)*zb(i6,i3)*
     &       (zb(i3,i4)*zb(i5,i6) + zb(i5,i3)*zb(i6,i4)) + 
     &      s(i2,i6)*(za(i3,i6)*zb(i3,i4)*zb(i6,i3) - 
     &         za(i6,i5)*(zb(i3,i4)*zb(i5,i6) + zb(i5,i3)*zb(i6,i4)))))/
     &  (za(i2,i6)*zab2(i5,i2,i6,i3)**2*zb(i3,i4)**2*zb(i4,i1))
     &)


      aaaa_NMHV_c(12) = (
     &  (t(i1,i3,i4)*zb(i2,i1)*zb(i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)) + 
     & (t(i1,i3,i4)*za(i1,i6)*zb(i2,i1)*zb(i6,i2)**2)/
     &  (za(i1,i5)*za(i6,i5)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i1,i4)*za(i2,i3)*zb(i5,i2)*zb(i6,i2)**2)/
     &  (t(i1,i4,i5)*za(i1,i5)*zb(i3,i2)*zb(i4,i1)) - 
     & (za(i1,i3)*za(i2,i4)*zb(i5,i2)*zb(i6,i2)**2)/
     &  (t(i1,i3,i5)*za(i1,i5)*zb(i3,i1)*zb(i4,i2)) + 
     & (za(i1,i4)*za(i3,i6)*zb(i5,i6)*zb(i6,i2)**2)/
     &  (t(i1,i4,i5)*za(i1,i5)*zb(i3,i2)*zb(i4,i1)) + 
     & (za(i1,i3)*za(i4,i6)*zb(i5,i6)*zb(i6,i2)**2)/
     &  (t(i1,i3,i5)*za(i1,i5)*zb(i3,i1)*zb(i4,i2))
     &)


      aaaa_NMHV_c(13) = (
     &  (t(i1,i3,i4)*za(i1,i5)*zb(i2,i1)*zb(i5,i2)**2)/
     &  (za(i1,i6)*za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)) + 
     & (t(i1,i3,i4)*zb(i2,i1)*zb(i5,i2)*zb(i6,i2))/
     &  (za(i5,i6)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)) - 
     & (za(i1,i4)*za(i2,i3)*zb(i5,i2)**2*zb(i6,i2))/
     &  (t(i1,i4,i6)*za(i1,i6)*zb(i3,i2)*zb(i4,i1)) - 
     & (za(i1,i3)*za(i2,i4)*zb(i5,i2)**2*zb(i6,i2))/
     &  (t(i1,i3,i6)*za(i1,i6)*zb(i3,i1)*zb(i4,i2)) + 
     & (za(i1,i4)*za(i3,i5)*zb(i5,i2)**2*zb(i6,i5))/
     &  (t(i1,i4,i6)*za(i1,i6)*zb(i3,i2)*zb(i4,i1)) + 
     & (za(i1,i3)*za(i4,i5)*zb(i5,i2)**2*zb(i6,i5))/
     &  (t(i1,i3,i6)*za(i1,i6)*zb(i3,i1)*zb(i4,i2))
     &)


      aaaa_NMHV_c(14) = (
     &  (za(i1,i5)*za(i2,i3)*za(i2,i4)**2*za(i4,i5)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)**3*za(i2,i6)*zb(i3,i1)) + 
     & (za(i1,i5)*za(i2,i3)*za(i2,i4)**2*za(i4,i6)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)**2*zb(i3,i1)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i4,i6)**2*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**3*zb(i3,i1)) - 
     & (s(i1,i3)*za(i1,i5)*za(i2,i4)**2*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)*zb(i3,i1)*zb(i4,i3)) + 
     & (s(i1,i3)*za(i1,i2)*za(i2,i4)*za(i4,i6)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**2*zb(i3,i1)*zb(i4,i3)) - 
     & (s(i1,i3)*za(i1,i2)*za(i2,i4)*zb(i2,i1)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)*zb(i3,i1)*zb(i4,i1)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i1,i5)*za(i2,i4)**2*za(i4,i5)*zb(i4,i1)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)**3*za(i2,i6)*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i1,i5)*za(i2,i4)**2*za(i4,i6)*zb(i4,i1)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)**2*zb(i3,i1)*zb(i4,i3)) + 
     & (za(i1,i2)**2*za(i2,i4)*za(i4,i6)**2*zb(i4,i1)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**3*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)*zb(i4,i1)**2*zb(i4,i3)) - 
     & (za(i1,i5)*za(i2,i3)*za(i2,i4)**2*zb(i4,i2)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i4,i6)*zb(i4,i2)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**2*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i2,i3)*za(i2,i4)**2*za(i4,i5)*zb(i6,i4))/
     &  (za(i2,i5)**3*za(i2,i6)*zb(i3,i1)) + 
     & (za(i2,i3)*za(i2,i4)**2*za(i4,i6)*zb(i6,i4))/
     &  (za(i2,i5)**2*za(i2,i6)**2*zb(i3,i1)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i4,i6)**2*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**3*zb(i3,i1)) - 
     & (s(i1,i3)*za(i2,i4)**2*zb(i6,i4))/
     &  (za(i2,i5)**2*za(i2,i6)*zb(i3,i1)*zb(i4,i3)) + 
     & (s(i1,i3)*za(i1,i2)*za(i2,i4)*za(i4,i6)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**2*zb(i3,i1)*zb(i4,i3)) - 
     & (s(i1,i3)*za(i1,i2)*za(i2,i4)*zb(i2,i1)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*zb(i3,i1)*zb(i4,i1)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i2,i4)**2*za(i4,i5)*zb(i4,i1)*zb(i6,i4))/
     &  (za(i2,i5)**3*za(i2,i6)*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i2,i4)**2*za(i4,i6)*zb(i4,i1)*zb(i6,i4))/
     &  (za(i2,i5)**2*za(i2,i6)**2*zb(i3,i1)*zb(i4,i3)) + 
     & (za(i1,i2)**2*za(i2,i4)*za(i4,i6)**2*zb(i4,i1)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**3*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*zb(i4,i1)**2*zb(i4,i3)) - 
     & (za(i2,i3)*za(i2,i4)**2*zb(i4,i2)*zb(i6,i4))/
     &  (za(i2,i5)**2*za(i2,i6)*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i4,i6)*zb(i4,i2)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**2*zb(i4,i1)*zb(i4,i3)) - 
     & (za(i1,i3)*za(i2,i4)**2*za(i4,i5)**2*zb(i5,i4)*zb(i6,i4))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)**3*zb(i3,i1)) - 
     & (za(i1,i3)*za(i2,i4)**2*za(i4,i6)**2*zb(i5,i4)*zb(i6,i4))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)**3*zb(i3,i1)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i4)*zb(i4,i1)**2*zb(i4,i3))
     &  + (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i4)*zb(i4,i1)**2*zb(i4,i3))
     &  + (s(i1,i5)*za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i4)**2*zb(i4,i1)*zb(i4,i3))
     &  - (za(i1,i2)*za(i2,i4)**2*za(i3,i6)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i6)**2*zab2(i2,i1,i5,i4)*zb(i4,i1)*zb(i4,i3))
     &  + (s(i1,i6)*za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i4)**2*zb(i4,i1)*zb(i4,i3))
     &  - (za(i1,i2)*za(i2,i4)**2*za(i3,i5)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i1,i6)*za(i2,i5)**2*zab2(i2,i1,i6,i4)*zb(i4,i1)*zb(i4,i3))
     &  + (za(i1,i2)**2*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2)*
     &    zb(i5,i4)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i4)**2*zb(i4,i1)*zb(i4,i3))
     &  + (za(i1,i2)**2*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2)*
     &    zb(i5,i4)*zb(i6,i4))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i4)**2*zb(i4,i1)*zb(i4,i3))
     &  + (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i2,i5)*zb(i4,i2)*zb(i5,i2)*
     &    zb(i5,i4)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i4)**2*zb(i4,i1)*zb(i4,i3))
     &  + (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i2,i6)*zb(i4,i2)*zb(i5,i4)*
     &    zb(i6,i2)*zb(i6,i4))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i4)**2*zb(i4,i1)*zb(i4,i3))
     &  - (za(i1,i3)*za(i2,i4)*za(i4,i5)**2*zb(i5,i4)*zb(i6,i5))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)**2*zb(i3,i1)) - 
     & (za(i1,i2)*za(i2,i4)*za(i3,i5)*zb(i4,i2)*zb(i5,i4)*zb(i6,i5))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i4)*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i1,i3)*za(i2,i4)*za(i4,i6)**2*zb(i6,i4)*zb(i6,i5))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)**2*zb(i3,i1)) + 
     & (za(i1,i2)*za(i2,i4)*za(i3,i6)*zb(i4,i2)*zb(i6,i4)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i4)*zb(i4,i1)*zb(i4,i3))
     &)


      aaaa_NMHV_c(15) = (
     &  -((za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*zb(i5,i3))/
     &    (za(i1,i6)*za(i2,i5)*za(i2,i6)*zb(i3,i1)**2*zb(i3,i4))) - 
     & (za(i1,i5)*za(i2,i3)**2*za(i2,i4)*zb(i3,i2)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i3,i6)*zb(i3,i2)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**2*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i1,i5)*za(i2,i3)**2*za(i2,i4)*za(i3,i5)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)**3*za(i2,i6)*zb(i4,i1)) + 
     & (za(i1,i5)*za(i2,i3)**2*za(i2,i4)*za(i3,i6)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)**2*zb(i4,i1)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i3,i6)**2*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**3*zb(i4,i1)) - 
     & (s(i1,i4)*za(i1,i5)*za(i2,i3)**2*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)*zb(i3,i4)*zb(i4,i1)) + 
     & (s(i1,i4)*za(i1,i2)*za(i2,i3)*za(i3,i6)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**2*zb(i3,i4)*zb(i4,i1)) - 
     & (s(i1,i4)*za(i1,i2)*za(i2,i3)*zb(i2,i1)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)*zb(i3,i1)*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i1,i2)*za(i1,i5)*za(i2,i3)**2*za(i3,i5)*zb(i3,i1)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)**3*za(i2,i6)*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i1,i2)*za(i1,i5)*za(i2,i3)**2*za(i3,i6)*zb(i3,i1)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i2,i6)**2*zb(i3,i4)*zb(i4,i1)) + 
     & (za(i1,i2)**2*za(i2,i3)*za(i3,i6)**2*zb(i3,i1)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)**3*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*zb(i3,i1)**2*zb(i3,i4)) - 
     & (za(i2,i3)**2*za(i2,i4)*zb(i3,i2)*zb(i6,i3))/
     &  (za(i2,i5)**2*za(i2,i6)*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i3,i6)*zb(i3,i2)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**2*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i2,i3)**2*za(i2,i4)*za(i3,i5)*zb(i6,i3))/
     &  (za(i2,i5)**3*za(i2,i6)*zb(i4,i1)) + 
     & (za(i2,i3)**2*za(i2,i4)*za(i3,i6)*zb(i6,i3))/
     &  (za(i2,i5)**2*za(i2,i6)**2*zb(i4,i1)) - 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i3,i6)**2*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**3*zb(i4,i1)) - 
     & (s(i1,i4)*za(i2,i3)**2*zb(i6,i3))/
     &  (za(i2,i5)**2*za(i2,i6)*zb(i3,i4)*zb(i4,i1)) + 
     & (s(i1,i4)*za(i1,i2)*za(i2,i3)*za(i3,i6)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**2*zb(i3,i4)*zb(i4,i1)) - 
     & (s(i1,i4)*za(i1,i2)*za(i2,i3)*zb(i2,i1)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*zb(i3,i1)*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i1,i2)*za(i2,i3)**2*za(i3,i5)*zb(i3,i1)*zb(i6,i3))/
     &  (za(i2,i5)**3*za(i2,i6)*zb(i3,i4)*zb(i4,i1)) - 
     & (za(i1,i2)*za(i2,i3)**2*za(i3,i6)*zb(i3,i1)*zb(i6,i3))/
     &  (za(i2,i5)**2*za(i2,i6)**2*zb(i3,i4)*zb(i4,i1)) + 
     & (za(i1,i2)**2*za(i2,i3)*za(i3,i6)**2*zb(i3,i1)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)**3*zb(i3,i4)*zb(i4,i1)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i3)*zb(i3,i1)**2*zb(i3,i4))
     &  + (za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i3)*zb(i3,i1)**2*zb(i3,i4))
     &  + (s(i1,i5)*za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i3)**2*zb(i3,i1)*zb(i3,i4))
     &  - (za(i1,i2)*za(i2,i3)**2*za(i4,i6)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i6)**2*zab2(i2,i1,i5,i3)*zb(i3,i1)*zb(i3,i4))
     &  + (s(i1,i6)*za(i1,i2)*za(i2,i3)*za(i2,i4)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i3)**2*zb(i3,i1)*zb(i3,i4))
     &  - (za(i1,i2)*za(i2,i3)**2*za(i4,i5)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i1,i6)*za(i2,i5)**2*zab2(i2,i1,i6,i3)*zb(i3,i1)*zb(i3,i4))
     &  + (za(i1,i2)**2*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &    zb(i5,i3)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i3)**2*zb(i3,i1)*zb(i3,i4))
     &  + (za(i1,i2)**2*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &    zb(i5,i3)*zb(i6,i3))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i3)**2*zb(i3,i1)*zb(i3,i4))
     &  - (za(i1,i4)*za(i2,i3)**2*za(i3,i5)**2*zb(i5,i3)*zb(i6,i3))/
     &  (t(i1,i4,i6)*za(i1,i6)*za(i2,i5)**3*zb(i4,i1)) - 
     & (za(i1,i4)*za(i2,i3)**2*za(i3,i6)**2*zb(i5,i3)*zb(i6,i3))/
     &  (t(i1,i4,i5)*za(i1,i5)*za(i2,i6)**3*zb(i4,i1)) + 
     & (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i5,i2)*
     &    zb(i5,i3)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i3)**2*zb(i3,i1)*zb(i3,i4))
     &  + (za(i1,i2)*za(i2,i3)*za(i2,i4)*za(i2,i6)*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i2)*zb(i6,i3))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i3)**2*zb(i3,i1)*zb(i3,i4))
     &  - (za(i1,i2)*za(i2,i3)*za(i4,i5)*zb(i3,i2)*zb(i5,i3)*zb(i6,i5))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i2,i1,i6,i3)*zb(i3,i1)*zb(i3,i4)) - 
     & (za(i1,i4)*za(i2,i3)*za(i3,i5)**2*zb(i5,i3)*zb(i6,i5))/
     &  (t(i1,i4,i6)*za(i1,i6)*za(i2,i5)**2*zb(i4,i1)) + 
     & (za(i1,i2)*za(i2,i3)*za(i4,i6)*zb(i3,i2)*zb(i6,i3)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i2,i1,i5,i3)*zb(i3,i1)*zb(i3,i4)) + 
     & (za(i1,i4)*za(i2,i3)*za(i3,i6)**2*zb(i6,i3)*zb(i6,i5))/
     &  (t(i1,i4,i5)*za(i1,i5)*za(i2,i6)**2*zb(i4,i1))
     &)


      aaaa_NMHV_c(16) = (
     &  (t(i4,i5,i6)*zab2(i1,i3,i2,i1)*zab2(i4,i3,i2,i1)**2*
     &   zb(i2,i1)**2)/
     & (zab2(i5,i3,i2,i1)**2*zab2(i6,i3,i2,i1)**2*zb(i3,i1)*zb(i3,i2))
     &)


      aaaa_NMHV_c(18) = (
     &  (2*t(i4,i5,i6)*zab2(i3,i1,i2,i3)*zab2(i4,i1,i2,i3)**2*
     &   zb(i3,i2))/
     & (zab2(i5,i1,i2,i3)**2*zab2(i6,i1,i2,i3)**2*zb(i3,i1))
     &)


      aaaa_NMHV_c(19) = (
     &  (-2*t(i3,i4,i5)*za(i1,i6)*zab2(i6,i1,i2,i5)**2*
     &   zab2(i6,i1,i2,i6))/
     & (za(i2,i6)*zab2(i6,i1,i2,i3)**2*zab2(i6,i1,i2,i4)**2)
     &)


      aaaa_NMHV_c(21) = (
     &  (t(i1,i2,i3)*za(i4,i6)**2*zab2(i6,i4,i5,i2)*zab2(i6,i4,i5,i6)*
     &    zb(i2,i1))/
     &  (za(i5,i6)**2*zab2(i6,i4,i5,i1)**2*zab2(i6,i4,i5,i3)*zb(i3,i2))
     &  + (za(i4,i6)**2*zab2(i5,i1,i3,i2)*zab2(i6,i4,i5,i6)*zb(i2,i1))/
     &  (za(i5,i6)**3*zab2(i6,i4,i5,i1)*zb(i3,i1)*zb(i3,i2)) + 
     & (za(i4,i5)*za(i4,i6)*zab2(i6,i4,i5,i2)*zab2(i6,i4,i5,i6)*
     &    zb(i2,i1))/
     &  (za(i5,i6)**3*zab2(i6,i4,i5,i1)*zb(i3,i1)*zb(i3,i2)) + 
     & (t(i1,i2,i3)*za(i4,i6)**2*zab2(i6,i4,i5,i6)*zb(i3,i2))/
     &  (za(i5,i6)**2*zab2(i6,i4,i5,i3)**2*zb(i3,i1)) + 
     & (za(i4,i5)*za(i4,i6)*zab2(i6,i4,i5,i6)*zb(i5,i2))/
     &  (za(i5,i6)**2*zab2(i6,i4,i5,i3)*zb(i3,i1)) - 
     & (za(i4,i5)*za(i4,i6)*zab2(i6,i4,i5,i2)**2*zab2(i6,i4,i5,i6)*
     &    zb(i5,i3))/
     &  (za(i5,i6)**2*zab2(i6,i4,i5,i1)*zab2(i6,i4,i5,i3)**2*zb(i3,i2))
     &  + (za(i4,i6)**2*zab2(i6,i4,i5,i6)*zb(i6,i2))/
     &  (za(i5,i6)**2*zab2(i6,i4,i5,i3)*zb(i3,i1)) - 
     & (za(i4,i6)**2*zab2(i6,i4,i5,i2)**2*zab2(i6,i4,i5,i6)*zb(i6,i3))/
     &  (za(i5,i6)**2*zab2(i6,i4,i5,i1)*zab2(i6,i4,i5,i3)**2*zb(i3,i2))
     &)


      aaaa_NMHV_c(22) = (
     &  (-2*t(i3,i4,i6)*za(i1,i5)*zab2(i5,i1,i2,i5)*
     &   zab2(i5,i1,i2,i6)**2)/
     & (za(i2,i5)*zab2(i5,i1,i2,i3)**2*zab2(i5,i1,i2,i4)**2)
     &)


      aaaa_NMHV_c(23) = (
     &  (t(i1,i2,i3)*za(i4,i5)**2*zab2(i5,i4,i6,i2)*zab2(i5,i4,i6,i5)*
     &    zb(i2,i1))/
     &  (za(i6,i5)**2*zab2(i5,i4,i6,i1)**2*zab2(i5,i4,i6,i3)*zb(i3,i2))
     &  + (za(i4,i5)*za(i4,i6)*zab2(i5,i4,i6,i2)*zab2(i5,i4,i6,i5)*
     &    zb(i2,i1))/
     &  (za(i6,i5)**3*zab2(i5,i4,i6,i1)*zb(i3,i1)*zb(i3,i2)) + 
     & (za(i4,i5)**2*zab2(i5,i4,i6,i5)*zab2(i6,i1,i3,i2)*zb(i2,i1))/
     &  (za(i6,i5)**3*zab2(i5,i4,i6,i1)*zb(i3,i1)*zb(i3,i2)) + 
     & (t(i1,i2,i3)*za(i4,i5)**2*zab2(i5,i4,i6,i5)*zb(i3,i2))/
     &  (za(i6,i5)**2*zab2(i5,i4,i6,i3)**2*zb(i3,i1)) + 
     & (za(i4,i5)**2*zab2(i5,i4,i6,i5)*zb(i5,i2))/
     &  (za(i6,i5)**2*zab2(i5,i4,i6,i3)*zb(i3,i1)) - 
     & (za(i4,i5)**2*zab2(i5,i4,i6,i2)**2*zab2(i5,i4,i6,i5)*zb(i5,i3))/
     &  (za(i6,i5)**2*zab2(i5,i4,i6,i1)*zab2(i5,i4,i6,i3)**2*zb(i3,i2))
     &  + (za(i4,i5)*za(i4,i6)*zab2(i5,i4,i6,i5)*zb(i6,i2))/
     &  (za(i6,i5)**2*zab2(i5,i4,i6,i3)*zb(i3,i1)) - 
     & (za(i4,i5)*za(i4,i6)*zab2(i5,i4,i6,i2)**2*zab2(i5,i4,i6,i5)*
     &    zb(i6,i3))/
     &  (za(i6,i5)**2*zab2(i5,i4,i6,i1)*zab2(i5,i4,i6,i3)**2*zb(i3,i2))
     &)


      aaaa_NMHV_c(24) = (
     &  -((zab2(i1,i4,i6,i1)*zab2(i3,i4,i6,i1)**2*zb(i6,i1)*
     &     (-(za(i1,i2)*zb(i4,i1)*zb(i6,i1)) + 
     &       (zab2(i2,i4,i6,i1) + za(i2,i4)*zb(i4,i1))*zb(i6,i4)))/
     &   (t(i2,i3,i5)*za(i2,i5)*zab2(i2,i4,i6,i1)*zab2(i5,i4,i6,i1)*
     &     zb(i4,i1)**3))
     &)


      aaaa_NMHV_c(25) = (
     &  (2*t(i3,i5,i6)*zab2(i3,i1,i2,i4)**2*zab2(i4,i1,i2,i4)*
     &   zb(i4,i2))/
     & (zab2(i5,i1,i2,i4)**2*zab2(i6,i1,i2,i4)**2*zb(i4,i1))
     &)


      aaaa_NMHV_c(26) = (
     &  -((zab2(i1,i4,i5,i1)*zab2(i3,i4,i5,i1)**2*zb(i5,i1)*
     &     (-(za(i1,i2)*zb(i4,i1)*zb(i5,i1)) + 
     &       (zab2(i2,i4,i5,i1) + za(i2,i4)*zb(i4,i1))*zb(i5,i4)))/
     &   (t(i2,i3,i6)*za(i2,i6)*zab2(i2,i4,i5,i1)*zab2(i6,i4,i5,i1)*
     &     zb(i4,i1)**3))
     &)


      aaaa_NMHV_c(27) = (
     &  (t(i3,i5,i6)*zab2(i1,i4,i2,i1)*zab2(i3,i4,i2,i1)**2*
     &   zb(i2,i1)**2)/
     & (zab2(i5,i4,i2,i1)**2*zab2(i6,i4,i2,i1)**2*zb(i4,i1)*zb(i4,i2))
     &)


      aaaa_NMHV_c(29) = (
     &  (t(i1,i2,i4)*za(i3,i6)**2*zab2(i6,i3,i5,i2)*zab2(i6,i3,i5,i6)*
     &    zb(i2,i1))/
     &  (za(i5,i6)**2*zab2(i6,i3,i5,i1)**2*zab2(i6,i3,i5,i4)*zb(i4,i2))
     &  + (za(i3,i6)**2*zab2(i5,i1,i4,i2)*zab2(i6,i3,i5,i6)*zb(i2,i1))/
     &  (za(i5,i6)**3*zab2(i6,i3,i5,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (za(i3,i5)*za(i3,i6)*zab2(i6,i3,i5,i2)*zab2(i6,i3,i5,i6)*
     &    zb(i2,i1))/
     &  (za(i5,i6)**3*zab2(i6,i3,i5,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (t(i1,i2,i4)*za(i3,i6)**2*zab2(i6,i3,i5,i6)*zb(i4,i2))/
     &  (za(i5,i6)**2*zab2(i6,i3,i5,i4)**2*zb(i4,i1)) + 
     & (za(i3,i5)*za(i3,i6)*zab2(i6,i3,i5,i6)*zb(i5,i2))/
     &  (za(i5,i6)**2*zab2(i6,i3,i5,i4)*zb(i4,i1)) - 
     & (za(i3,i5)*za(i3,i6)*zab2(i6,i3,i5,i2)**2*zab2(i6,i3,i5,i6)*
     &    zb(i5,i4))/
     &  (za(i5,i6)**2*zab2(i6,i3,i5,i1)*zab2(i6,i3,i5,i4)**2*zb(i4,i2))
     &  + (za(i3,i6)**2*zab2(i6,i3,i5,i6)*zb(i6,i2))/
     &  (za(i5,i6)**2*zab2(i6,i3,i5,i4)*zb(i4,i1)) - 
     & (za(i3,i6)**2*zab2(i6,i3,i5,i2)**2*zab2(i6,i3,i5,i6)*zb(i6,i4))/
     &  (za(i5,i6)**2*zab2(i6,i3,i5,i1)*zab2(i6,i3,i5,i4)**2*zb(i4,i2))
     &)


      aaaa_NMHV_c(30) = (
     &  (t(i1,i2,i4)*za(i3,i5)**2*zab2(i5,i3,i6,i2)*zab2(i5,i3,i6,i5)*
     &    zb(i2,i1))/
     &  (za(i6,i5)**2*zab2(i5,i3,i6,i1)**2*zab2(i5,i3,i6,i4)*zb(i4,i2))
     &  + (za(i3,i5)*za(i3,i6)*zab2(i5,i3,i6,i2)*zab2(i5,i3,i6,i5)*
     &    zb(i2,i1))/
     &  (za(i6,i5)**3*zab2(i5,i3,i6,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (za(i3,i5)**2*zab2(i5,i3,i6,i5)*zab2(i6,i1,i4,i2)*zb(i2,i1))/
     &  (za(i6,i5)**3*zab2(i5,i3,i6,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (t(i1,i2,i4)*za(i3,i5)**2*zab2(i5,i3,i6,i5)*zb(i4,i2))/
     &  (za(i6,i5)**2*zab2(i5,i3,i6,i4)**2*zb(i4,i1)) + 
     & (za(i3,i5)**2*zab2(i5,i3,i6,i5)*zb(i5,i2))/
     &  (za(i6,i5)**2*zab2(i5,i3,i6,i4)*zb(i4,i1)) - 
     & (za(i3,i5)**2*zab2(i5,i3,i6,i2)**2*zab2(i5,i3,i6,i5)*zb(i5,i4))/
     &  (za(i6,i5)**2*zab2(i5,i3,i6,i1)*zab2(i5,i3,i6,i4)**2*zb(i4,i2))
     &  + (za(i3,i5)*za(i3,i6)*zab2(i5,i3,i6,i5)*zb(i6,i2))/
     &  (za(i6,i5)**2*zab2(i5,i3,i6,i4)*zb(i4,i1)) - 
     & (za(i3,i5)*za(i3,i6)*zab2(i5,i3,i6,i2)**2*zab2(i5,i3,i6,i5)*
     &    zb(i6,i4))/
     &  (za(i6,i5)**2*zab2(i5,i3,i6,i1)*zab2(i5,i3,i6,i4)**2*zb(i4,i2))
     &)


      aaaa_NMHV_c(31) = (
     &  -((zab2(i1,i3,i6,i1)*zab2(i4,i3,i6,i1)**2*zb(i6,i1)*
     &     (-(za(i1,i2)*zb(i3,i1)*zb(i6,i1)) + 
     &       (zab2(i2,i3,i6,i1) + za(i2,i3)*zb(i3,i1))*zb(i6,i3)))/
     &   (t(i2,i4,i5)*za(i2,i5)*zab2(i2,i3,i6,i1)*zab2(i5,i3,i6,i1)*
     &     zb(i3,i1)**3))
     &)


      aaaa_NMHV_c(32) = (
     &  -((zab2(i1,i3,i5,i1)*zab2(i4,i3,i5,i1)**2*zb(i5,i1)*
     &     (-(za(i1,i2)*zb(i3,i1)*zb(i5,i1)) + 
     &       (zab2(i2,i3,i5,i1) + za(i2,i3)*zb(i3,i1))*zb(i5,i3)))/
     &   (t(i2,i4,i6)*za(i2,i6)*zab2(i2,i3,i5,i1)*zab2(i6,i3,i5,i1)*
     &     zb(i3,i1)**3))
     &)


      aaaa_NMHV_c(34) = (
     &  -((za(i1,i2)*zab2(i1,i5,i2,i4)*zab2(i3,i4,i6,i3)*zb(i6,i3)**2)/
     &    (za(i1,i5)*za(i2,i5)*zab2(i2,i4,i6,i3)*zb(i4,i3)**3)) - 
     & (t(i1,i2,i5)*za(i1,i5)*zab2(i3,i4,i6,i3)*zb(i6,i3)**2)/
     &  (za(i2,i5)*zab2(i5,i4,i6,i3)**2*zb(i4,i3)**2) - 
     & (za(i3,i5)*zab2(i1,i4,i6,i3)**2*zab2(i3,i4,i6,i3)*zb(i6,i3)**2)/
     &  (za(i1,i5)*zab2(i2,i4,i6,i3)*zab2(i5,i4,i6,i3)**2*zb(i4,i3)**2)
     &  - (za(i1,i3)*zab2(i3,i4,i6,i3)*zb(i6,i3)**2)/
     &  (za(i2,i5)*zab2(i5,i4,i6,i3)*zb(i4,i3)**2) + 
     & (t(i1,i2,i5)*za(i1,i2)*zab2(i1,i4,i6,i3)*zab2(i3,i4,i6,i3)*
     &    zb(i6,i3)**2)/
     &  (za(i1,i5)*zab2(i2,i4,i6,i3)**2*zab2(i5,i4,i6,i3)*zb(i4,i3)**2)
     &  - (za(i1,i2)*zab2(i1,i4,i6,i3)*zab2(i3,i4,i6,i3)*zb(i6,i3)*
     &    zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i5)*zab2(i2,i4,i6,i3)*zb(i4,i3)**3) - 
     & (za(i4,i5)*zab2(i1,i4,i6,i3)**2*zab2(i3,i4,i6,i3)*zb(i6,i3)*
     &    zb(i6,i4))/
     &  (za(i1,i5)*zab2(i2,i4,i6,i3)*zab2(i5,i4,i6,i3)**2*zb(i4,i3)**2)
     &  - (za(i1,i4)*zab2(i3,i4,i6,i3)*zb(i6,i3)*zb(i6,i4))/
     &  (za(i2,i5)*zab2(i5,i4,i6,i3)*zb(i4,i3)**2)
     &)


      aaaa_NMHV_c(36) = (
     &  -((za(i1,i2)*zab2(i1,i3,i6,i4)*zab2(i4,i3,i6,i4)*zb(i6,i3)*
     &      zb(i6,i4))/
     &    (za(i1,i5)*za(i2,i5)*zab2(i2,i3,i6,i4)*zb(i3,i4)**3)) - 
     & (za(i3,i5)*zab2(i1,i3,i6,i4)**2*zab2(i4,i3,i6,i4)*zb(i6,i3)*
     &    zb(i6,i4))/
     &  (za(i1,i5)*zab2(i2,i3,i6,i4)*zab2(i5,i3,i6,i4)**2*zb(i3,i4)**2)
     &  - (za(i1,i3)*zab2(i4,i3,i6,i4)*zb(i6,i3)*zb(i6,i4))/
     &  (za(i2,i5)*zab2(i5,i3,i6,i4)*zb(i3,i4)**2) - 
     & (za(i1,i2)*zab2(i1,i5,i2,i3)*zab2(i4,i3,i6,i4)*zb(i6,i4)**2)/
     &  (za(i1,i5)*za(i2,i5)*zab2(i2,i3,i6,i4)*zb(i3,i4)**3) - 
     & (t(i1,i2,i5)*za(i1,i5)*zab2(i4,i3,i6,i4)*zb(i6,i4)**2)/
     &  (za(i2,i5)*zab2(i5,i3,i6,i4)**2*zb(i3,i4)**2) - 
     & (za(i4,i5)*zab2(i1,i3,i6,i4)**2*zab2(i4,i3,i6,i4)*zb(i6,i4)**2)/
     &  (za(i1,i5)*zab2(i2,i3,i6,i4)*zab2(i5,i3,i6,i4)**2*zb(i3,i4)**2)
     &  - (za(i1,i4)*zab2(i4,i3,i6,i4)*zb(i6,i4)**2)/
     &  (za(i2,i5)*zab2(i5,i3,i6,i4)*zb(i3,i4)**2) + 
     & (t(i1,i2,i5)*za(i1,i2)*zab2(i1,i3,i6,i4)*zab2(i4,i3,i6,i4)*
     &    zb(i6,i4)**2)/
     &  (za(i1,i5)*zab2(i2,i3,i6,i4)**2*zab2(i5,i3,i6,i4)*zb(i3,i4)**2)
     &)


      aaaa_NMHV_c(39) = (
     &  -((za(i1,i2)*zab2(i1,i6,i2,i4)*zab2(i3,i4,i5,i3)*zb(i5,i3)**2)/
     &    (za(i1,i6)*za(i2,i6)*zab2(i2,i4,i5,i3)*zb(i4,i3)**3)) - 
     & (t(i1,i2,i6)*za(i1,i6)*zab2(i3,i4,i5,i3)*zb(i5,i3)**2)/
     &  (za(i2,i6)*zab2(i6,i4,i5,i3)**2*zb(i4,i3)**2) - 
     & (za(i3,i6)*zab2(i1,i4,i5,i3)**2*zab2(i3,i4,i5,i3)*zb(i5,i3)**2)/
     &  (za(i1,i6)*zab2(i2,i4,i5,i3)*zab2(i6,i4,i5,i3)**2*zb(i4,i3)**2)
     &  - (za(i1,i3)*zab2(i3,i4,i5,i3)*zb(i5,i3)**2)/
     &  (za(i2,i6)*zab2(i6,i4,i5,i3)*zb(i4,i3)**2) + 
     & (t(i1,i2,i6)*za(i1,i2)*zab2(i1,i4,i5,i3)*zab2(i3,i4,i5,i3)*
     &    zb(i5,i3)**2)/
     &  (za(i1,i6)*zab2(i2,i4,i5,i3)**2*zab2(i6,i4,i5,i3)*zb(i4,i3)**2)
     &  - (za(i1,i2)*zab2(i1,i4,i5,i3)*zab2(i3,i4,i5,i3)*zb(i5,i3)*
     &    zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i6)*zab2(i2,i4,i5,i3)*zb(i4,i3)**3) - 
     & (za(i4,i6)*zab2(i1,i4,i5,i3)**2*zab2(i3,i4,i5,i3)*zb(i5,i3)*
     &    zb(i5,i4))/
     &  (za(i1,i6)*zab2(i2,i4,i5,i3)*zab2(i6,i4,i5,i3)**2*zb(i4,i3)**2)
     &  - (za(i1,i4)*zab2(i3,i4,i5,i3)*zb(i5,i3)*zb(i5,i4))/
     &  (za(i2,i6)*zab2(i6,i4,i5,i3)*zb(i4,i3)**2)
     &)


      aaaa_NMHV_c(41) = (
     &  -((za(i1,i2)*zab2(i1,i3,i5,i4)*zab2(i4,i3,i5,i4)*zb(i5,i3)*
     &      zb(i5,i4))/
     &    (za(i1,i6)*za(i2,i6)*zab2(i2,i3,i5,i4)*zb(i3,i4)**3)) - 
     & (za(i3,i6)*zab2(i1,i3,i5,i4)**2*zab2(i4,i3,i5,i4)*zb(i5,i3)*
     &    zb(i5,i4))/
     &  (za(i1,i6)*zab2(i2,i3,i5,i4)*zab2(i6,i3,i5,i4)**2*zb(i3,i4)**2)
     &  - (za(i1,i3)*zab2(i4,i3,i5,i4)*zb(i5,i3)*zb(i5,i4))/
     &  (za(i2,i6)*zab2(i6,i3,i5,i4)*zb(i3,i4)**2) - 
     & (za(i1,i2)*zab2(i1,i6,i2,i3)*zab2(i4,i3,i5,i4)*zb(i5,i4)**2)/
     &  (za(i1,i6)*za(i2,i6)*zab2(i2,i3,i5,i4)*zb(i3,i4)**3) - 
     & (t(i1,i2,i6)*za(i1,i6)*zab2(i4,i3,i5,i4)*zb(i5,i4)**2)/
     &  (za(i2,i6)*zab2(i6,i3,i5,i4)**2*zb(i3,i4)**2) - 
     & (za(i4,i6)*zab2(i1,i3,i5,i4)**2*zab2(i4,i3,i5,i4)*zb(i5,i4)**2)/
     &  (za(i1,i6)*zab2(i2,i3,i5,i4)*zab2(i6,i3,i5,i4)**2*zb(i3,i4)**2)
     &  - (za(i1,i4)*zab2(i4,i3,i5,i4)*zb(i5,i4)**2)/
     &  (za(i2,i6)*zab2(i6,i3,i5,i4)*zb(i3,i4)**2) + 
     & (t(i1,i2,i6)*za(i1,i2)*zab2(i1,i3,i5,i4)*zab2(i4,i3,i5,i4)*
     &    zb(i5,i4)**2)/
     &  (za(i1,i6)*zab2(i2,i3,i5,i4)**2*zab2(i6,i3,i5,i4)*zb(i3,i4)**2)
     &)


      aaaa_NMHV_c(43) = (
     &  (za(i1,i3)**2*zab2(i6,i1,i3,i6)*zb(i5,i2)**2)/
     & (t(i2,i4,i5)*za(i1,i6)*zab2(i6,i1,i3,i4)*zb(i4,i2))
     &)


      aaaa_NMHV_c(44) = (
     &  (za(i1,i3)**2*zab2(i5,i1,i3,i5)*zb(i6,i2)**2)/
     & (t(i2,i4,i6)*za(i1,i5)*zab2(i5,i1,i3,i4)*zb(i4,i2))
     &)


      aaaa_NMHV_c(45) = (
     &  (t(i2,i5,i6)**2*zab2(i2,i1,i3,i4)*zab2(i4,i1,i3,i4))/
     & (za(i2,i5)*za(i2,i6)*zab2(i5,i1,i3,i4)*zab2(i6,i1,i3,i4)*
     &   zb(i3,i1)*zb(i4,i3))
     &)


      aaaa_NMHV_c(47) = (
     &  -((t(i1,i3,i4)**2*zab2(i6,i5,i2,i1)*zab2(i6,i5,i2,i6))/
     &   (za(i2,i5)*za(i5,i6)*zab2(i6,i5,i2,i3)*zab2(i6,i5,i2,i4)*
     &     zb(i3,i1)*zb(i4,i1)))
     &)


      aaaa_NMHV_c(48) = (
     &  -((t(i1,i3,i4)**2*zab2(i5,i6,i2,i1)*zab2(i5,i6,i2,i5))/
     &   (za(i2,i6)*za(i6,i5)*zab2(i5,i6,i2,i3)*zab2(i5,i6,i2,i4)*
     &     zb(i3,i1)*zb(i4,i1)))
     &)


      aaaa_NMHV_c(49) = (
     &  -((za(i2,i4)*zab2(i2,i4,i6,i2)*zab2(i2,i4,i6,i5)**2*
     &     (-(za(i2,i4)*za(i2,i6)*zb(i2,i1)) + 
     &       za(i4,i6)*(zab2(i2,i4,i6,i1) + za(i2,i6)*zb(i6,i1))))/
     &   (t(i1,i3,i5)*za(i2,i6)**3*zab2(i2,i4,i6,i1)*zab2(i2,i4,i6,i3)*
     &     zb(i3,i1)))
     &)


      aaaa_NMHV_c(50) = (
     &  (za(i4,i6)**2*zab2(i2,i1,i3,i5)*zab2(i2,i3,i5,i1)*
     &    zab2(i6,i2,i4,i5)*zab2(i6,i2,i4,i6))/
     &  (t(i1,i3,i5)*za(i2,i6)**3*zab2(i6,i2,i4,i1)*zab2(i6,i2,i4,i3)*
     &    zb(i3,i1)) + (za(i4,i6)**2*zab2(i2,i3,i5,i1)*
     &    zab2(i6,i2,i4,i5)*zab2(i6,i2,i4,i6)*zb(i5,i3))/
     &  (za(i2,i6)**2*zab2(i6,i2,i4,i1)*zab2(i6,i2,i4,i3)**2*zb(i3,i1))
     &  + (t(i1,i3,i5)*za(i4,i6)**2*zab2(i6,i2,i4,i6)*zb(i5,i3)**2)/
     &  (za(i2,i6)*zab2(i6,i2,i4,i3)**3*zb(i3,i1))
     &)


      aaaa_NMHV_c(51) = (
     &  -((za(i1,i3)**2*zab2(i4,i6,i2,i4)*zb(i6,i2)**2)/
     &   (t(i1,i3,i5)*za(i1,i5)*zab2(i5,i6,i2,i4)*zb(i4,i2)))
     &)


      aaaa_NMHV_c(52) = (
     &  -((za(i2,i4)*zab2(i2,i4,i5,i2)*zab2(i2,i4,i5,i6)**2*
     &     (-(za(i2,i4)*za(i2,i5)*zb(i2,i1)) + 
     &       za(i4,i5)*(zab2(i2,i4,i5,i1) + za(i2,i5)*zb(i5,i1))))/
     &   (t(i1,i3,i6)*za(i2,i5)**3*zab2(i2,i4,i5,i1)*zab2(i2,i4,i5,i3)*
     &     zb(i3,i1)))
     &)


      aaaa_NMHV_c(53) = (
     &  (za(i4,i5)**2*zab2(i2,i1,i3,i6)*zab2(i2,i3,i6,i1)*
     &    zab2(i5,i2,i4,i5)*zab2(i5,i2,i4,i6))/
     &  (t(i1,i3,i6)*za(i2,i5)**3*zab2(i5,i2,i4,i1)*zab2(i5,i2,i4,i3)*
     &    zb(i3,i1)) + (za(i4,i5)**2*zab2(i2,i3,i6,i1)*
     &    zab2(i5,i2,i4,i5)*zab2(i5,i2,i4,i6)*zb(i6,i3))/
     &  (za(i2,i5)**2*zab2(i5,i2,i4,i1)*zab2(i5,i2,i4,i3)**2*zb(i3,i1))
     &  + (t(i1,i3,i6)*za(i4,i5)**2*zab2(i5,i2,i4,i5)*zb(i6,i3)**2)/
     &  (za(i2,i5)*zab2(i5,i2,i4,i3)**3*zb(i3,i1))
     &)


      aaaa_NMHV_c(54) = (
     &  -((za(i1,i3)**2*zab2(i4,i5,i2,i4)*zb(i5,i2)**2)/
     &   (t(i1,i3,i6)*za(i1,i6)*zab2(i6,i5,i2,i4)*zb(i4,i2)))
     &)


      aaaa_NMHV_c(56) = (
     &  (za(i1,i4)**2*zab2(i6,i1,i4,i6)*zb(i5,i2)**2)/
     & (t(i2,i3,i5)*za(i1,i6)*zab2(i6,i1,i4,i3)*zb(i3,i2))
     &)


      aaaa_NMHV_c(57) = (
     &  (za(i1,i4)**2*zab2(i5,i1,i4,i5)*zb(i6,i2)**2)/
     & (t(i2,i3,i6)*za(i1,i5)*zab2(i5,i1,i4,i3)*zb(i3,i2))
     &)


      aaaa_NMHV_c(58) = (
     &  (t(i2,i5,i6)**2*zab2(i2,i1,i4,i3)*zab2(i3,i1,i4,i3))/
     & (za(i2,i5)*za(i2,i6)*zab2(i5,i1,i4,i3)*zab2(i6,i1,i4,i3)*
     &   zb(i3,i4)*zb(i4,i1))
     &)


      aaaa_NMHV_c(59) = (
     &  -((za(i2,i3)*zab2(i2,i3,i6,i2)*zab2(i2,i3,i6,i5)**2*
     &     (-(za(i2,i3)*za(i2,i6)*zb(i2,i1)) + 
     &       za(i3,i6)*(zab2(i2,i3,i6,i1) + za(i2,i6)*zb(i6,i1))))/
     &   (t(i1,i4,i5)*za(i2,i6)**3*zab2(i2,i3,i6,i1)*zab2(i2,i3,i6,i4)*
     &     zb(i4,i1)))
     &)


      aaaa_NMHV_c(60) = (
     &  (za(i3,i6)**2*zab2(i2,i1,i4,i5)*zab2(i2,i4,i5,i1)*
     &    zab2(i6,i2,i3,i5)*zab2(i6,i2,i3,i6))/
     &  (t(i1,i4,i5)*za(i2,i6)**3*zab2(i6,i2,i3,i1)*zab2(i6,i2,i3,i4)*
     &    zb(i4,i1)) + (za(i3,i6)**2*zab2(i2,i4,i5,i1)*
     &    zab2(i6,i2,i3,i5)*zab2(i6,i2,i3,i6)*zb(i5,i4))/
     &  (za(i2,i6)**2*zab2(i6,i2,i3,i1)*zab2(i6,i2,i3,i4)**2*zb(i4,i1))
     &  + (t(i1,i4,i5)*za(i3,i6)**2*zab2(i6,i2,i3,i6)*zb(i5,i4)**2)/
     &  (za(i2,i6)*zab2(i6,i2,i3,i4)**3*zb(i4,i1))
     &)


      aaaa_NMHV_c(61) = (
     &  -((za(i1,i4)**2*zab2(i3,i6,i2,i3)*zb(i6,i2)**2)/
     &   (t(i1,i4,i5)*za(i1,i5)*zab2(i5,i6,i2,i3)*zb(i3,i2)))
     &)


      aaaa_NMHV_c(62) = (
     &  -((za(i2,i3)*zab2(i2,i3,i5,i2)*zab2(i2,i3,i5,i6)**2*
     &     (-(za(i2,i3)*za(i2,i5)*zb(i2,i1)) + 
     &       za(i3,i5)*(zab2(i2,i3,i5,i1) + za(i2,i5)*zb(i5,i1))))/
     &   (t(i1,i4,i6)*za(i2,i5)**3*zab2(i2,i3,i5,i1)*zab2(i2,i3,i5,i4)*
     &     zb(i4,i1)))
     &)


      aaaa_NMHV_c(63) = (
     &  (za(i3,i5)**2*zab2(i2,i1,i4,i6)*zab2(i2,i4,i6,i1)*
     &    zab2(i5,i2,i3,i5)*zab2(i5,i2,i3,i6))/
     &  (t(i1,i4,i6)*za(i2,i5)**3*zab2(i5,i2,i3,i1)*zab2(i5,i2,i3,i4)*
     &    zb(i4,i1)) + (za(i3,i5)**2*zab2(i2,i4,i6,i1)*
     &    zab2(i5,i2,i3,i5)*zab2(i5,i2,i3,i6)*zb(i6,i4))/
     &  (za(i2,i5)**2*zab2(i5,i2,i3,i1)*zab2(i5,i2,i3,i4)**2*zb(i4,i1))
     &  + (t(i1,i4,i6)*za(i3,i5)**2*zab2(i5,i2,i3,i5)*zb(i6,i4)**2)/
     &  (za(i2,i5)*zab2(i5,i2,i3,i4)**3*zb(i4,i1))
     &)


      aaaa_NMHV_c(64) = (
     &  -((za(i1,i4)**2*zab2(i3,i5,i2,i3)*zb(i5,i2)**2)/
     &   (t(i1,i4,i6)*za(i1,i6)*zab2(i6,i5,i2,i3)*zb(i3,i2)))
     &)


      aaaa_NMHV_c(65) = (
     &  -((t(i3,i4,i6)*za(i1,i2)**2*zab2(i2,i1,i5,i2)*
     &     zab2(i2,i1,i5,i6)**2)/
     &   (za(i1,i5)*za(i2,i5)*zab2(i2,i1,i5,i3)**2*zab2(i2,i1,i5,i4)**2)
     &   )
     &)


      aaaa_NMHV_c(66) = (
     &  -((zab2(i1,i2,i3,i4)**2*zab2(i6,i1,i5,i2)**3*zab2(i6,i1,i5,i6))/
     &    (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)*zab2(i6,i1,i5,i4)**3*
     &      zb(i3,i2)*zb(i4,i2))) + 
     & (t(i1,i5,i6)**2*za(i1,i6)**2*zab2(i6,i1,i5,i2)**2*zb(i6,i3))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)**3*zab2(i6,i1,i5,i4)*
     &    zb(i4,i2)) + (t(i1,i5,i6)**2*za(i1,i6)**2*
     &    zab2(i6,i1,i5,i2)**2*zb(i6,i4))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)**2*zab2(i6,i1,i5,i4)**2*
     &    zb(i4,i2)) - (2*t(i1,i5,i6)*za(i1,i6)*zab2(i6,i1,i5,i2)**2*
     &    zb(i6,i5))/
     &  (za(i5,i6)*zab2(i6,i1,i5,i3)**2*zab2(i6,i1,i5,i4)*zb(i4,i2))
     &)


      aaaa_NMHV_c(67) = (
     &  -((za(i3,i6)*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)**3*
     &      zb(i5,i4)**2)/
     &    (t(i2,i3,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**3*zb(i4,i1))
     &    ) - (za(i3,i4)*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)**2*
     &    zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**2*zb(i4,i1))
     &  - (za(i3,i4)**2*zab2(i4,i1,i5,i4)*zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i6,i1,i5,i4)*zb(i4,i1)) - 
     & (za(i3,i5)*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)**2*zb(i5,i1)*
     &    zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**2*
     &    zb(i4,i1)**2) - (2*za(i3,i4)*za(i3,i5)*zab2(i4,i1,i5,i4)*
     &    zb(i5,i1)*zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i6,i1,i5,i4)*zb(i4,i1)**2) - 
     & (2*za(i1,i5)*za(i3,i6)*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)**2*
     &    zb(i5,i1)*zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**3*zb(i4,i1))
     &  - (za(i1,i5)*za(i3,i4)*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)*
     &    zb(i5,i1)*zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**2*zb(i4,i1))
     &  - (za(i2,i5)*za(i3,i5)*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)*
     &    zb(i5,i1)**2*zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)*
     &    zb(i4,i1)**3) - (za(i1,i5)*za(i2,i5)*za(i3,i6)*
     &    zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i4)*zb(i5,i1)**2*zb(i5,i4)**2)
     &   /(t(i2,i3,i6)*za(i2,i6)*zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**2*
     &    zb(i4,i1)**2) - (za(i1,i5)**2*za(i3,i6)**2*zab2(i4,i1,i5,i4)*
     &    zb(i5,i1)**2*zb(i5,i4)**2)/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i6,i1,i5,i4)**3*zb(i4,i1))
     &)


      aaaa_NMHV_c(68) = (
     &  -((za(i4,i6)*zab2(i3,i1,i5,i3)**3*zab2(i4,i1,i5,i3)*
     &      zb(i5,i3)**2)/
     &    (t(i2,i4,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**3*zb(i3,i1))
     &    ) - (za(i4,i3)*zab2(i3,i1,i5,i3)**2*zab2(i4,i1,i5,i3)*
     &    zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**2*zb(i3,i1))
     &  - (za(i4,i3)**2*zab2(i3,i1,i5,i3)*zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i6,i1,i5,i3)*zb(i3,i1)) - 
     & (za(i4,i5)*zab2(i3,i1,i5,i3)**2*zab2(i4,i1,i5,i3)*zb(i5,i1)*
     &    zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**2*
     &    zb(i3,i1)**2) - (2*za(i4,i3)*za(i4,i5)*zab2(i3,i1,i5,i3)*
     &    zb(i5,i1)*zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i6,i1,i5,i3)*zb(i3,i1)**2) - 
     & (2*za(i1,i5)*za(i4,i6)*zab2(i3,i1,i5,i3)**2*zab2(i4,i1,i5,i3)*
     &    zb(i5,i1)*zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**3*zb(i3,i1))
     &  - (za(i1,i5)*za(i4,i3)*zab2(i3,i1,i5,i3)*zab2(i4,i1,i5,i3)*
     &    zb(i5,i1)*zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**2*zb(i3,i1))
     &  - (za(i2,i5)*za(i4,i5)*zab2(i3,i1,i5,i3)*zab2(i4,i1,i5,i3)*
     &    zb(i5,i1)**2*zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)*
     &    zb(i3,i1)**3) - (za(i1,i5)*za(i2,i5)*za(i4,i6)*
     &    zab2(i3,i1,i5,i3)*zab2(i4,i1,i5,i3)*zb(i5,i1)**2*zb(i5,i3)**2)
     &   /(t(i2,i4,i6)*za(i2,i6)*zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**2*
     &    zb(i3,i1)**2) - (za(i1,i5)**2*za(i4,i6)**2*zab2(i3,i1,i5,i3)*
     &    zb(i5,i1)**2*zb(i5,i3)**2)/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i6,i1,i5,i3)**3*zb(i3,i1))
     &)


      aaaa_NMHV_c(70) = (
     &  (zab2(i1,i3,i2,i4)**3*zab2(i4,i3,i2,i4)*zab2(i6,i5,i1,i2)**2)/
     &  (za(i1,i5)*za(i1,i6)*zab2(i5,i3,i2,i4)*zab2(i6,i3,i2,i4)**3*
     &    zb(i3,i2)*zb(i4,i3)) + 
     & (2*t(i2,i3,i4)*za(i3,i4)*zab2(i1,i3,i2,i4)**2*zb(i4,i2))/
     &  (za(i1,i6)*zab2(i5,i3,i2,i4)**2*zab2(i6,i3,i2,i4)*zb(i4,i3)) + 
     & (t(i2,i3,i4)**2*za(i4,i6)*zab2(i1,i3,i2,i4)**2*zb(i4,i2)**2)/
     &  (za(i1,i6)*zab2(i5,i3,i2,i4)**2*zab2(i6,i3,i2,i4)**2*zb(i3,i2)*
     &    zb(i4,i3)) + (t(i2,i3,i4)**2*za(i4,i5)*zab2(i1,i3,i2,i4)**2*
     &    zb(i4,i2)**2)/
     &  (za(i1,i6)*zab2(i5,i3,i2,i4)**3*zab2(i6,i3,i2,i4)*zb(i3,i2)*
     &    zb(i4,i3))
     &)


      aaaa_NMHV_c(71) = (
     &  (2*t(i2,i3,i4)*za(i4,i3)*zab2(i1,i4,i2,i3)**2*zb(i3,i2))/
     &  (za(i1,i6)*zab2(i5,i4,i2,i3)**2*zab2(i6,i4,i2,i3)*zb(i3,i4)) + 
     & (zab2(i1,i4,i2,i3)**3*zab2(i3,i4,i2,i3)*zab2(i6,i5,i1,i2)**2)/
     &  (za(i1,i5)*za(i1,i6)*zab2(i5,i4,i2,i3)*zab2(i6,i4,i2,i3)**3*
     &    zb(i3,i4)*zb(i4,i2)) + 
     & (t(i2,i3,i4)**2*za(i3,i6)*zab2(i1,i4,i2,i3)**2*zb(i3,i2)**2)/
     &  (za(i1,i6)*zab2(i5,i4,i2,i3)**2*zab2(i6,i4,i2,i3)**2*zb(i3,i4)*
     &    zb(i4,i2)) + (t(i2,i3,i4)**2*za(i3,i5)*zab2(i1,i4,i2,i3)**2*
     &    zb(i3,i2)**2)/
     &  (za(i1,i6)*zab2(i5,i4,i2,i3)**3*zab2(i6,i4,i2,i3)*zb(i3,i4)*
     &    zb(i4,i2))
     &)


      aaaa_NMHV_c(72) = (
     &  -((t(i3,i4,i5)*za(i1,i2)**2*zab2(i2,i1,i6,i2)*
     &     zab2(i2,i1,i6,i5)**2)/
     &   (za(i1,i6)*za(i2,i6)*zab2(i2,i1,i6,i3)**2*zab2(i2,i1,i6,i4)**2)
     &   )
     &)


      aaaa_NMHV_c(73) = (
     &  -((zab2(i1,i2,i3,i4)**2*zab2(i5,i1,i6,i2)**3*zab2(i5,i1,i6,i5))/
     &    (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i3)*zab2(i5,i1,i6,i4)**3*
     &      zb(i3,i2)*zb(i4,i2))) + 
     & (t(i1,i5,i6)**2*za(i1,i5)**2*zab2(i5,i1,i6,i2)**2*zb(i5,i3))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i3)**3*zab2(i5,i1,i6,i4)*
     &    zb(i4,i2)) + (t(i1,i5,i6)**2*za(i1,i5)**2*
     &    zab2(i5,i1,i6,i2)**2*zb(i5,i4))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i1,i6,i3)**2*zab2(i5,i1,i6,i4)**2*
     &    zb(i4,i2)) - (2*t(i1,i5,i6)*za(i1,i5)*zab2(i5,i1,i6,i2)**2*
     &    zb(i5,i6))/
     &  (za(i6,i5)*zab2(i5,i1,i6,i3)**2*zab2(i5,i1,i6,i4)*zb(i4,i2))
     &)


      aaaa_NMHV_c(74) = (
     &  -((za(i3,i5)*zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)**3*
     &      zb(i6,i4)**2)/
     &    (t(i2,i3,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**3*zb(i4,i1))
     &    ) - (za(i3,i4)*zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)**2*
     &    zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**2*zb(i4,i1))
     &  - (za(i3,i4)**2*zab2(i4,i1,i6,i4)*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i5,i1,i6,i4)*zb(i4,i1)) - 
     & (za(i3,i6)*zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)**2*zb(i6,i1)*
     &    zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**2*
     &    zb(i4,i1)**2) - (2*za(i3,i4)*za(i3,i6)*zab2(i4,i1,i6,i4)*
     &    zb(i6,i1)*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i5,i1,i6,i4)*zb(i4,i1)**2) - 
     & (2*za(i1,i6)*za(i3,i5)*zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)**2*
     &    zb(i6,i1)*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**3*zb(i4,i1))
     &  - (za(i1,i6)*za(i3,i4)*zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)*
     &    zb(i6,i1)*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**2*zb(i4,i1))
     &  - (za(i2,i6)*za(i3,i6)*zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)*
     &    zb(i6,i1)**2*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)*
     &    zb(i4,i1)**3) - (za(i1,i6)*za(i2,i6)*za(i3,i5)*
     &    zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i4)*zb(i6,i1)**2*zb(i6,i4)**2)
     &   /(t(i2,i3,i5)*za(i2,i5)*zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**2*
     &    zb(i4,i1)**2) - (za(i1,i6)**2*za(i3,i5)**2*zab2(i4,i1,i6,i4)*
     &    zb(i6,i1)**2*zb(i6,i4)**2)/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i5,i1,i6,i4)**3*zb(i4,i1))
     &)


      aaaa_NMHV_c(75) = (
     &  -((za(i4,i5)*zab2(i3,i1,i6,i3)**3*zab2(i4,i1,i6,i3)*
     &      zb(i6,i3)**2)/
     &    (t(i2,i4,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**3*zb(i3,i1))
     &    ) - (za(i4,i3)*zab2(i3,i1,i6,i3)**2*zab2(i4,i1,i6,i3)*
     &    zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**2*zb(i3,i1))
     &  - (za(i4,i3)**2*zab2(i3,i1,i6,i3)*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i5,i1,i6,i3)*zb(i3,i1)) - 
     & (za(i4,i6)*zab2(i3,i1,i6,i3)**2*zab2(i4,i1,i6,i3)*zb(i6,i1)*
     &    zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**2*
     &    zb(i3,i1)**2) - (2*za(i4,i3)*za(i4,i6)*zab2(i3,i1,i6,i3)*
     &    zb(i6,i1)*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i5,i1,i6,i3)*zb(i3,i1)**2) - 
     & (2*za(i1,i6)*za(i4,i5)*zab2(i3,i1,i6,i3)**2*zab2(i4,i1,i6,i3)*
     &    zb(i6,i1)*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**3*zb(i3,i1))
     &  - (za(i1,i6)*za(i4,i3)*zab2(i3,i1,i6,i3)*zab2(i4,i1,i6,i3)*
     &    zb(i6,i1)*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**2*zb(i3,i1))
     &  - (za(i2,i6)*za(i4,i6)*zab2(i3,i1,i6,i3)*zab2(i4,i1,i6,i3)*
     &    zb(i6,i1)**2*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)*
     &    zb(i3,i1)**3) - (za(i1,i6)*za(i2,i6)*za(i4,i5)*
     &    zab2(i3,i1,i6,i3)*zab2(i4,i1,i6,i3)*zb(i6,i1)**2*zb(i6,i3)**2)
     &   /(t(i2,i4,i5)*za(i2,i5)*zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**2*
     &    zb(i3,i1)**2) - (za(i1,i6)**2*za(i4,i5)**2*zab2(i3,i1,i6,i3)*
     &    zb(i6,i1)**2*zb(i6,i3)**2)/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i5,i1,i6,i3)**3*zb(i3,i1))
     &)


      aaaa_NMHV_c(77) = (
     &  ((-(sqrtd1246*za(i1,i2)) + zaa2(i1,i1,i2,i4,i6,i2) + 
     &      zaa2(i2,i1,i2,i4,i6,i1))*
     &    (-(sqrtd1246*za(i3,i5)) + zaa2(i3,i1,i2,i4,i6,i5) + 
     &       zaa2(i5,i1,i2,i4,i6,i3))**2*
     &    (-(sqrtd1246*za(i4,i6)) + zaa2(i4,i1,i2,i4,i6,i6) + 
     &       zaa2(i6,i1,i2,i4,i6,i4))**2)/
     &  (16.*zaa2(i2,i1,i2,i4,i6,i2)*zaa2(i5,i1,i2,i4,i6,i5)**2*
     &    zaa2(i6,i1,i2,i4,i6,i6)**2) + 
     & ((sqrtd1246*za(i1,i2) + zaa2(i1,i1,i2,i4,i6,i2) + 
     &      zaa2(i2,i1,i2,i4,i6,i1))*
     &    (sqrtd1246*za(i3,i5) + zaa2(i3,i1,i2,i4,i6,i5) + 
     &       zaa2(i5,i1,i2,i4,i6,i3))**2*
     &    (sqrtd1246*za(i4,i6) + zaa2(i4,i1,i2,i4,i6,i6) + 
     &       zaa2(i6,i1,i2,i4,i6,i4))**2)/
     &  (16.*zaa2(i2,i1,i2,i4,i6,i2)*zaa2(i5,i1,i2,i4,i6,i5)**2*
     &    zaa2(i6,i1,i2,i4,i6,i6)**2)
     &)


      aaaa_NMHV_c(78) = (
     &  ((-(sqrtd1245*za(i1,i2)) + zaa2(i1,i1,i2,i4,i5,i2) + 
     &      zaa2(i2,i1,i2,i4,i5,i1))*
     &    (-(sqrtd1245*za(i3,i5)) + zaa2(i3,i1,i2,i4,i5,i5) + 
     &       zaa2(i5,i1,i2,i4,i5,i3))**2*
     &    (-(sqrtd1245*za(i4,i6)) + zaa2(i4,i1,i2,i4,i5,i6) + 
     &       zaa2(i6,i1,i2,i4,i5,i4))**2)/
     &  (16.*zaa2(i2,i1,i2,i4,i5,i2)*zaa2(i5,i1,i2,i4,i5,i5)**2*
     &    zaa2(i6,i1,i2,i4,i5,i6)**2) + 
     & ((sqrtd1245*za(i1,i2) + zaa2(i1,i1,i2,i4,i5,i2) + 
     &      zaa2(i2,i1,i2,i4,i5,i1))*
     &    (sqrtd1245*za(i3,i5) + zaa2(i3,i1,i2,i4,i5,i5) + 
     &       zaa2(i5,i1,i2,i4,i5,i3))**2*
     &    (sqrtd1245*za(i4,i6) + zaa2(i4,i1,i2,i4,i5,i6) + 
     &       zaa2(i6,i1,i2,i4,i5,i4))**2)/
     &  (16.*zaa2(i2,i1,i2,i4,i5,i2)*zaa2(i5,i1,i2,i4,i5,i5)**2*
     &    zaa2(i6,i1,i2,i4,i5,i6)**2)
     &)


      aaaa_NMHV_c(85) = (
     &  (-sqrtd2346 + s(i1,i5) + s(i2,i3) - s(i4,i6))*
     &  ((-(s(i4,i6)*(-(sqrtd2346*za(i1,i2)) + 
     &             zaa2(i1,i2,i3,i4,i6,i2) + zaa2(i2,i2,i3,i4,i6,i1))*
     &           zaa2(i4,i2,i3,i4,i6,i4)*
     &           (-(sqrtd2346*za(i1,i5)) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &             zaa2(i5,i2,i3,i4,i6,i1))*
     &           (-(sqrtd2346*za(i3,i6)) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &              zaa2(i6,i2,i3,i4,i6,i3))**2)/
     &        (32.*sqrtd2346**2*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &          zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**2) - 
     &       (s(i4,i6)*za(i4,i5)*
     &          (-(sqrtd2346*za(i1,i2)) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &            zaa2(i2,i2,i3,i4,i6,i1))*
     &          (-(sqrtd2346*za(i1,i5)) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i1))*
     &          (-(sqrtd2346*za(i3,i5)) + zaa2(i3,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i3))*
     &          (-(sqrtd2346*za(i3,i6)) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &            zaa2(i6,i2,i3,i4,i6,i3))*
     &          (-(sqrtd2346*za(i4,i6)) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &            zaa2(i6,i2,i3,i4,i6,i4)))/
     &        (64.*sqrtd2346*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &          zaa2(i5,i2,i3,i4,i6,i5)**2*zaa2(i6,i2,i3,i4,i6,i6)**2))/
     &     (sqrtd2346/2. + (s(i1,i5) - s(i2,i3) - s(i4,i6))/2.) + 
     &    (sqrtd2346/2. + (s(i1,i5) - s(i2,i3) - s(i4,i6))/2.)*
     &     ((za(i3,i6)*(sqrtd2346*za(i1,i2) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &            zaa2(i2,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i1,i5) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i3,i6) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &            zaa2(i6,i2,i3,i4,i6,i3))*
     &          (sqrtd2346*za(i4,i6) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &             zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &        (64.*sqrtd2346*s(i2,i3)*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &          zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**3) - 
     &       (za(i1,i2)*(-(sqrtd2346*za(i1,i2)) + 
     &            zaa2(i1,i2,i3,i4,i6,i2) + zaa2(i2,i2,i3,i4,i6,i1))*
     &          (-(sqrtd2346*za(i3,i2)) + zaa2(i2,i2,i3,i4,i6,i3) + 
     &            zaa2(i3,i2,i3,i4,i6,i2))*
     &          (-(sqrtd2346*za(i3,i5)) + zaa2(i3,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i3))*
     &          (-(sqrtd2346*za(i4,i6)) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &             zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &        (64.*sqrtd2346*s(i2,i3)*za(i1,i5)*
     &          zaa2(i2,i2,i3,i4,i6,i2)**2*zaa2(i5,i2,i3,i4,i6,i5)*
     &          zaa2(i6,i2,i3,i4,i6,i6)**2) - 
     &       ((sqrtd2346*za(i1,i2) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &            zaa2(i2,i2,i3,i4,i6,i1))*zaa2(i3,i2,i3,i4,i6,i3)*
     &          (sqrtd2346*za(i1,i5) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i4,i6) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &             zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &        (32.*sqrtd2346**2*s(i2,i3)*za(i1,i5)*
     &          zaa2(i2,i2,i3,i4,i6,i2)*zaa2(i5,i2,i3,i4,i6,i5)*
     &          zaa2(i6,i2,i3,i4,i6,i6)**2)) - 
     &    (zaa2(i1,i2,i3,i4,i6,i6)*
     &       (-(sqrtd2346*za(i1,i2)) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &         zaa2(i2,i2,i3,i4,i6,i1))*
     &       (-(sqrtd2346*za(i3,i5)) + zaa2(i3,i2,i3,i4,i6,i5) + 
     &         zaa2(i5,i2,i3,i4,i6,i3))*
     &       (-(sqrtd2346*za(i3,i6)) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &         zaa2(i6,i2,i3,i4,i6,i3))*
     &       (-(sqrtd2346*za(i4,i6)) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &          zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &     (64.*sqrtd2346*s(i2,i3)*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &       zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**3)) + 
     & (-sqrtd2346 + s(i1,i5) - s(i2,i3) + s(i4,i6))*
     &  ((sqrtd2346/2. + (s(i1,i5) - s(i2,i3) - s(i4,i6))/2.)**2*
     &     ((za(i4,i6)*(sqrtd2346*za(i1,i2) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &            zaa2(i2,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i1,i5) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i3,i6) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &             zaa2(i6,i2,i3,i4,i6,i3))**2*
     &          (sqrtd2346*za(i4,i6) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &            zaa2(i6,i2,i3,i4,i6,i4)))/
     &        (64.*sqrtd2346*s(i2,i3)*s(i4,i6)*za(i1,i5)*
     &          zaa2(i2,i2,i3,i4,i6,i2)*zaa2(i5,i2,i3,i4,i6,i5)*
     &          zaa2(i6,i2,i3,i4,i6,i6)**3) + 
     &       (za(i4,i5)*(sqrtd2346*za(i1,i2) + 
     &            zaa2(i1,i2,i3,i4,i6,i2) + zaa2(i2,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i1,i5) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i1))*
     &          (sqrtd2346*za(i3,i5) + zaa2(i3,i2,i3,i4,i6,i5) + 
     &            zaa2(i5,i2,i3,i4,i6,i3))*
     &          (sqrtd2346*za(i3,i6) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &            zaa2(i6,i2,i3,i4,i6,i3))*
     &          (sqrtd2346*za(i4,i6) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &            zaa2(i6,i2,i3,i4,i6,i4)))/
     &        (64.*sqrtd2346*s(i2,i3)*s(i4,i6)*za(i1,i5)*
     &          zaa2(i2,i2,i3,i4,i6,i2)*zaa2(i5,i2,i3,i4,i6,i5)**2*
     &          zaa2(i6,i2,i3,i4,i6,i6)**2)) - 
     &    (za(i3,i6)*(-(sqrtd2346*za(i1,i2)) + 
     &         zaa2(i1,i2,i3,i4,i6,i2) + zaa2(i2,i2,i3,i4,i6,i1))*
     &       (-(sqrtd2346*za(i1,i5)) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &         zaa2(i5,i2,i3,i4,i6,i1))*
     &       (-(sqrtd2346*za(i3,i6)) + zaa2(i3,i2,i3,i4,i6,i6) + 
     &         zaa2(i6,i2,i3,i4,i6,i3))*
     &       (-(sqrtd2346*za(i4,i6)) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &          zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &     (64.*sqrtd2346*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &       zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**3) - 
     &    ((-(sqrtd2346*za(i1,i2)) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &         zaa2(i2,i2,i3,i4,i6,i1))*zaa2(i3,i2,i3,i4,i6,i3)*
     &       (-(sqrtd2346*za(i1,i5)) + zaa2(i1,i2,i3,i4,i6,i5) + 
     &         zaa2(i5,i2,i3,i4,i6,i1))*
     &       (-(sqrtd2346*za(i4,i6)) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &          zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &     (32.*sqrtd2346**2*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &       zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**2) + 
     &    (za(i1,i2)*(sqrtd2346*za(i1,i2) + zaa2(i1,i2,i3,i4,i6,i2) + 
     &         zaa2(i2,i2,i3,i4,i6,i1))*
     &       (sqrtd2346*za(i3,i2) + zaa2(i2,i2,i3,i4,i6,i3) + 
     &         zaa2(i3,i2,i3,i4,i6,i2))*
     &       (sqrtd2346*za(i3,i5) + zaa2(i3,i2,i3,i4,i6,i5) + 
     &         zaa2(i5,i2,i3,i4,i6,i3))*
     &       (sqrtd2346*za(i4,i6) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &          zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &     (64.*sqrtd2346*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)**2*
     &       zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**2) - 
     &    (zaa2(i1,i2,i3,i4,i6,i1)*
     &       (sqrtd2346*za(i3,i2) + zaa2(i2,i2,i3,i4,i6,i3) + 
     &         zaa2(i3,i2,i3,i4,i6,i2))*
     &       (sqrtd2346*za(i3,i5) + zaa2(i3,i2,i3,i4,i6,i5) + 
     &         zaa2(i5,i2,i3,i4,i6,i3))*
     &       (sqrtd2346*za(i4,i6) + zaa2(i4,i2,i3,i4,i6,i6) + 
     &          zaa2(i6,i2,i3,i4,i6,i4))**2)/
     &     (32.*sqrtd2346**2*za(i1,i5)*zaa2(i2,i2,i3,i4,i6,i2)*
     &       zaa2(i5,i2,i3,i4,i6,i5)*zaa2(i6,i2,i3,i4,i6,i6)**2))
     &)


      aaaa_NMHV_c(86) = (
     &  (-sqrtd2436 + s(i1,i5) + s(i2,i4) - s(i3,i6))*
     &  ((-(s(i3,i6)*za(i3,i5)*
     &           (-(sqrtd2436*za(i1,i2)) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &             zaa2(i2,i2,i4,i3,i6,i1))*
     &           (-(sqrtd2436*za(i1,i5)) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &             zaa2(i5,i2,i4,i3,i6,i1))*
     &           (-(sqrtd2436*za(i3,i5)) + zaa2(i3,i2,i4,i3,i6,i5) + 
     &             zaa2(i5,i2,i4,i3,i6,i3))*
     &           (-(sqrtd2436*za(i4,i6)) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &              zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &        (64.*sqrtd2436*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &          zaa2(i5,i2,i4,i3,i6,i5)**2*zaa2(i6,i2,i4,i3,i6,i6)**2)
     &        - (s(i3,i6)*(-(sqrtd2436*za(i1,i2)) + 
     &            zaa2(i1,i2,i4,i3,i6,i2) + zaa2(i2,i2,i4,i3,i6,i1))*
     &          zaa2(i3,i2,i4,i3,i6,i3)*
     &          (-(sqrtd2436*za(i1,i5)) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i1))*
     &          (-(sqrtd2436*za(i4,i6)) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &             zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &        (32.*sqrtd2436**2*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &          zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**2))/
     &     (sqrtd2436/2. + (s(i1,i5) - s(i2,i4) - s(i3,i6))/2.) + 
     &    (sqrtd2436/2. + (s(i1,i5) - s(i2,i4) - s(i3,i6))/2.)*
     &     ((za(i4,i6)*(sqrtd2436*za(i1,i2) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &            zaa2(i2,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i1,i5) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i3,i6) + zaa2(i3,i2,i4,i3,i6,i6) + 
     &             zaa2(i6,i2,i4,i3,i6,i3))**2*
     &          (sqrtd2436*za(i4,i6) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &            zaa2(i6,i2,i4,i3,i6,i4)))/
     &        (64.*sqrtd2436*s(i2,i4)*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &          zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**3) - 
     &       ((sqrtd2436*za(i1,i2) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &            zaa2(i2,i2,i4,i3,i6,i1))*zaa2(i4,i2,i4,i3,i6,i4)*
     &          (sqrtd2436*za(i1,i5) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i3,i6) + zaa2(i3,i2,i4,i3,i6,i6) + 
     &             zaa2(i6,i2,i4,i3,i6,i3))**2)/
     &        (32.*sqrtd2436**2*s(i2,i4)*za(i1,i5)*
     &          zaa2(i2,i2,i4,i3,i6,i2)*zaa2(i5,i2,i4,i3,i6,i5)*
     &          zaa2(i6,i2,i4,i3,i6,i6)**2) - 
     &       (za(i1,i2)*(-(sqrtd2436*za(i1,i2)) + 
     &            zaa2(i1,i2,i4,i3,i6,i2) + zaa2(i2,i2,i4,i3,i6,i1))*
     &          (-(sqrtd2436*za(i3,i2)) + zaa2(i2,i2,i4,i3,i6,i3) + 
     &            zaa2(i3,i2,i4,i3,i6,i2))*
     &          (-(sqrtd2436*za(i3,i5)) + zaa2(i3,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i3))*
     &          (-(sqrtd2436*za(i4,i6)) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &             zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &        (64.*sqrtd2436*s(i2,i4)*za(i1,i5)*
     &          zaa2(i2,i2,i4,i3,i6,i2)**2*zaa2(i5,i2,i4,i3,i6,i5)*
     &          zaa2(i6,i2,i4,i3,i6,i6)**2)) - 
     &    (zaa2(i1,i2,i4,i3,i6,i6)*
     &       (-(sqrtd2436*za(i1,i2)) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &         zaa2(i2,i2,i4,i3,i6,i1))*
     &       (-(sqrtd2436*za(i3,i5)) + zaa2(i3,i2,i4,i3,i6,i5) + 
     &         zaa2(i5,i2,i4,i3,i6,i3))*
     &       (-(sqrtd2436*za(i3,i6)) + zaa2(i3,i2,i4,i3,i6,i6) + 
     &         zaa2(i6,i2,i4,i3,i6,i3))*
     &       (-(sqrtd2436*za(i4,i6)) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &          zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &     (64.*sqrtd2436*s(i2,i4)*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &       zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**3)) + 
     & (-sqrtd2436 + s(i1,i5) - s(i2,i4) + s(i3,i6))*
     &  ((sqrtd2436/2. + (s(i1,i5) - s(i2,i4) - s(i3,i6))/2.)**2*
     &     ((za(i3,i6)*(sqrtd2436*za(i1,i2) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &            zaa2(i2,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i1,i5) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i3,i6) + zaa2(i3,i2,i4,i3,i6,i6) + 
     &            zaa2(i6,i2,i4,i3,i6,i3))*
     &          (sqrtd2436*za(i4,i6) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &             zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &        (64.*sqrtd2436*s(i2,i4)*s(i3,i6)*za(i1,i5)*
     &          zaa2(i2,i2,i4,i3,i6,i2)*zaa2(i5,i2,i4,i3,i6,i5)*
     &          zaa2(i6,i2,i4,i3,i6,i6)**3) + 
     &       (za(i3,i5)*(sqrtd2436*za(i1,i2) + 
     &            zaa2(i1,i2,i4,i3,i6,i2) + zaa2(i2,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i1,i5) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i1))*
     &          (sqrtd2436*za(i3,i5) + zaa2(i3,i2,i4,i3,i6,i5) + 
     &            zaa2(i5,i2,i4,i3,i6,i3))*
     &          (sqrtd2436*za(i4,i6) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &             zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &        (64.*sqrtd2436*s(i2,i4)*s(i3,i6)*za(i1,i5)*
     &          zaa2(i2,i2,i4,i3,i6,i2)*zaa2(i5,i2,i4,i3,i6,i5)**2*
     &          zaa2(i6,i2,i4,i3,i6,i6)**2)) - 
     &    (za(i4,i6)*(-(sqrtd2436*za(i1,i2)) + 
     &         zaa2(i1,i2,i4,i3,i6,i2) + zaa2(i2,i2,i4,i3,i6,i1))*
     &       (-(sqrtd2436*za(i1,i5)) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &         zaa2(i5,i2,i4,i3,i6,i1))*
     &       (-(sqrtd2436*za(i3,i6)) + zaa2(i3,i2,i4,i3,i6,i6) + 
     &          zaa2(i6,i2,i4,i3,i6,i3))**2*
     &       (-(sqrtd2436*za(i4,i6)) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &         zaa2(i6,i2,i4,i3,i6,i4)))/
     &     (64.*sqrtd2436*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &       zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**3) - 
     &    ((-(sqrtd2436*za(i1,i2)) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &         zaa2(i2,i2,i4,i3,i6,i1))*zaa2(i4,i2,i4,i3,i6,i4)*
     &       (-(sqrtd2436*za(i1,i5)) + zaa2(i1,i2,i4,i3,i6,i5) + 
     &         zaa2(i5,i2,i4,i3,i6,i1))*
     &       (-(sqrtd2436*za(i3,i6)) + zaa2(i3,i2,i4,i3,i6,i6) + 
     &          zaa2(i6,i2,i4,i3,i6,i3))**2)/
     &     (32.*sqrtd2436**2*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &       zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**2) + 
     &    (za(i1,i2)*(sqrtd2436*za(i1,i2) + zaa2(i1,i2,i4,i3,i6,i2) + 
     &         zaa2(i2,i2,i4,i3,i6,i1))*
     &       (sqrtd2436*za(i3,i2) + zaa2(i2,i2,i4,i3,i6,i3) + 
     &         zaa2(i3,i2,i4,i3,i6,i2))*
     &       (sqrtd2436*za(i3,i5) + zaa2(i3,i2,i4,i3,i6,i5) + 
     &         zaa2(i5,i2,i4,i3,i6,i3))*
     &       (sqrtd2436*za(i4,i6) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &          zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &     (64.*sqrtd2436*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)**2*
     &       zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**2) - 
     &    (zaa2(i1,i2,i4,i3,i6,i1)*
     &       (sqrtd2436*za(i3,i2) + zaa2(i2,i2,i4,i3,i6,i3) + 
     &         zaa2(i3,i2,i4,i3,i6,i2))*
     &       (sqrtd2436*za(i3,i5) + zaa2(i3,i2,i4,i3,i6,i5) + 
     &         zaa2(i5,i2,i4,i3,i6,i3))*
     &       (sqrtd2436*za(i4,i6) + zaa2(i4,i2,i4,i3,i6,i6) + 
     &          zaa2(i6,i2,i4,i3,i6,i4))**2)/
     &     (32.*sqrtd2436**2*za(i1,i5)*zaa2(i2,i2,i4,i3,i6,i2)*
     &       zaa2(i5,i2,i4,i3,i6,i5)*zaa2(i6,i2,i4,i3,i6,i6)**2))
     &)


      aaaa_NMHV_c(88) = (
     &  (-sqrtd2345 + s(i1,i6) + s(i2,i3) - s(i4,i5))*
     &  ((-(s(i4,i5)*za(i4,i6)*
     &           (-(sqrtd2345*za(i1,i2)) + zaa2(i1,i2,i3,i4,i5,i2) + 
     &             zaa2(i2,i2,i3,i4,i5,i1))*
     &           (-(sqrtd2345*za(i1,i5)) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &             zaa2(i5,i2,i3,i4,i5,i1))*
     &           (-(sqrtd2345*za(i3,i5)) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &             zaa2(i5,i2,i3,i4,i5,i3))*
     &           (-(sqrtd2345*za(i3,i6)) + zaa2(i3,i2,i3,i4,i5,i6) + 
     &             zaa2(i6,i2,i3,i4,i5,i3))*
     &           (-(sqrtd2345*za(i4,i6)) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &             zaa2(i6,i2,i3,i4,i5,i4)))/
     &        (64.*sqrtd2345*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &          zaa2(i5,i2,i3,i4,i5,i5)**2*zaa2(i6,i2,i3,i4,i5,i6)**2)
     &        - (s(i4,i5)*(-(sqrtd2345*za(i1,i2)) + 
     &            zaa2(i1,i2,i3,i4,i5,i2) + zaa2(i2,i2,i3,i4,i5,i1))*
     &          zaa2(i4,i2,i3,i4,i5,i4)*
     &          (-(sqrtd2345*za(i1,i5)) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i1))*
     &          (-(sqrtd2345*za(i3,i5)) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i3))*
     &          (-(sqrtd2345*za(i3,i6)) + zaa2(i3,i2,i3,i4,i5,i6) + 
     &            zaa2(i6,i2,i3,i4,i5,i3)))/
     &        (32.*sqrtd2345**2*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &          zaa2(i5,i2,i3,i4,i5,i5)**2*zaa2(i6,i2,i3,i4,i5,i6)))/
     &     (sqrtd2345/2. + (s(i1,i6) - s(i2,i3) - s(i4,i5))/2.) + 
     &    (sqrtd2345/2. + (s(i1,i6) - s(i2,i3) - s(i4,i5))/2.)*
     &     (-(za(i1,i2)*(-(sqrtd2345*za(i1,i2)) + 
     &             zaa2(i1,i2,i3,i4,i5,i2) + zaa2(i2,i2,i3,i4,i5,i1))*
     &           (-(sqrtd2345*za(i3,i2)) + zaa2(i2,i2,i3,i4,i5,i3) + 
     &             zaa2(i3,i2,i3,i4,i5,i2))*
     &           (-(sqrtd2345*za(i3,i5)) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &             zaa2(i5,i2,i3,i4,i5,i3))*
     &           (-(sqrtd2345*za(i4,i5)) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &             zaa2(i5,i2,i3,i4,i5,i4))*
     &           (-(sqrtd2345*za(i4,i6)) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &             zaa2(i6,i2,i3,i4,i5,i4)))/
     &        (64.*sqrtd2345*s(i2,i3)*za(i1,i6)*
     &          zaa2(i2,i2,i3,i4,i5,i2)**2*zaa2(i5,i2,i3,i4,i5,i5)**2*
     &          zaa2(i6,i2,i3,i4,i5,i6)) + 
     &       (za(i3,i5)*(sqrtd2345*za(i1,i2) + 
     &            zaa2(i1,i2,i3,i4,i5,i2) + zaa2(i2,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i1,i5) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i3,i5) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i3))*
     &          (sqrtd2345*za(i4,i5) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i4))*
     &          (sqrtd2345*za(i4,i6) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &            zaa2(i6,i2,i3,i4,i5,i4)))/
     &        (64.*sqrtd2345*s(i2,i3)*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &          zaa2(i5,i2,i3,i4,i5,i5)**3*zaa2(i6,i2,i3,i4,i5,i6)) - 
     &       ((sqrtd2345*za(i1,i2) + zaa2(i1,i2,i3,i4,i5,i2) + 
     &            zaa2(i2,i2,i3,i4,i5,i1))*zaa2(i3,i2,i3,i4,i5,i3)*
     &          (sqrtd2345*za(i1,i5) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i4,i5) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i4))*
     &          (sqrtd2345*za(i4,i6) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &            zaa2(i6,i2,i3,i4,i5,i4)))/
     &        (32.*sqrtd2345**2*s(i2,i3)*za(i1,i6)*
     &          zaa2(i2,i2,i3,i4,i5,i2)*zaa2(i5,i2,i3,i4,i5,i5)**2*
     &          zaa2(i6,i2,i3,i4,i5,i6))) - 
     &    (zaa2(i1,i2,i3,i4,i5,i5)*
     &       (-(sqrtd2345*za(i1,i2)) + zaa2(i1,i2,i3,i4,i5,i2) + 
     &         zaa2(i2,i2,i3,i4,i5,i1))*
     &       (-(sqrtd2345*za(i3,i5)) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &          zaa2(i5,i2,i3,i4,i5,i3))**2*
     &       (-(sqrtd2345*za(i4,i5)) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i4))*
     &       (-(sqrtd2345*za(i4,i6)) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &         zaa2(i6,i2,i3,i4,i5,i4)))/
     &     (64.*sqrtd2345*s(i2,i3)*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &       zaa2(i5,i2,i3,i4,i5,i5)**3*zaa2(i6,i2,i3,i4,i5,i6))) + 
     & (-sqrtd2345 + s(i1,i6) - s(i2,i3) + s(i4,i5))*
     &  ((sqrtd2345/2. + (s(i1,i6) - s(i2,i3) - s(i4,i5))/2.)**2*
     &     ((za(i4,i6)*(sqrtd2345*za(i1,i2) + zaa2(i1,i2,i3,i4,i5,i2) + 
     &            zaa2(i2,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i1,i5) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i3,i5) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i3))*
     &          (sqrtd2345*za(i3,i6) + zaa2(i3,i2,i3,i4,i5,i6) + 
     &            zaa2(i6,i2,i3,i4,i5,i3))*
     &          (sqrtd2345*za(i4,i6) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &            zaa2(i6,i2,i3,i4,i5,i4)))/
     &        (64.*sqrtd2345*s(i2,i3)*s(i4,i5)*za(i1,i6)*
     &          zaa2(i2,i2,i3,i4,i5,i2)*zaa2(i5,i2,i3,i4,i5,i5)**2*
     &          zaa2(i6,i2,i3,i4,i5,i6)**2) + 
     &       (za(i4,i5)*(sqrtd2345*za(i1,i2) + 
     &            zaa2(i1,i2,i3,i4,i5,i2) + zaa2(i2,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i1,i5) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &            zaa2(i5,i2,i3,i4,i5,i1))*
     &          (sqrtd2345*za(i3,i5) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &             zaa2(i5,i2,i3,i4,i5,i3))**2*
     &          (sqrtd2345*za(i4,i6) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &            zaa2(i6,i2,i3,i4,i5,i4)))/
     &        (64.*sqrtd2345*s(i2,i3)*s(i4,i5)*za(i1,i6)*
     &          zaa2(i2,i2,i3,i4,i5,i2)*zaa2(i5,i2,i3,i4,i5,i5)**3*
     &          zaa2(i6,i2,i3,i4,i5,i6))) - 
     &    (za(i3,i5)*(-(sqrtd2345*za(i1,i2)) + 
     &         zaa2(i1,i2,i3,i4,i5,i2) + zaa2(i2,i2,i3,i4,i5,i1))*
     &       (-(sqrtd2345*za(i1,i5)) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i1))*
     &       (-(sqrtd2345*za(i3,i5)) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i3))*
     &       (-(sqrtd2345*za(i4,i5)) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i4))*
     &       (-(sqrtd2345*za(i4,i6)) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &         zaa2(i6,i2,i3,i4,i5,i4)))/
     &     (64.*sqrtd2345*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &       zaa2(i5,i2,i3,i4,i5,i5)**3*zaa2(i6,i2,i3,i4,i5,i6)) - 
     &    ((-(sqrtd2345*za(i1,i2)) + zaa2(i1,i2,i3,i4,i5,i2) + 
     &         zaa2(i2,i2,i3,i4,i5,i1))*zaa2(i3,i2,i3,i4,i5,i3)*
     &       (-(sqrtd2345*za(i1,i5)) + zaa2(i1,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i1))*
     &       (-(sqrtd2345*za(i4,i5)) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i4))*
     &       (-(sqrtd2345*za(i4,i6)) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &         zaa2(i6,i2,i3,i4,i5,i4)))/
     &     (32.*sqrtd2345**2*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &       zaa2(i5,i2,i3,i4,i5,i5)**2*zaa2(i6,i2,i3,i4,i5,i6)) + 
     &    (za(i1,i2)*(sqrtd2345*za(i1,i2) + zaa2(i1,i2,i3,i4,i5,i2) + 
     &         zaa2(i2,i2,i3,i4,i5,i1))*
     &       (sqrtd2345*za(i3,i2) + zaa2(i2,i2,i3,i4,i5,i3) + 
     &         zaa2(i3,i2,i3,i4,i5,i2))*
     &       (sqrtd2345*za(i3,i5) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i3))*
     &       (sqrtd2345*za(i4,i5) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i4))*
     &       (sqrtd2345*za(i4,i6) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &         zaa2(i6,i2,i3,i4,i5,i4)))/
     &     (64.*sqrtd2345*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)**2*
     &       zaa2(i5,i2,i3,i4,i5,i5)**2*zaa2(i6,i2,i3,i4,i5,i6)) - 
     &    (zaa2(i1,i2,i3,i4,i5,i1)*
     &       (sqrtd2345*za(i3,i2) + zaa2(i2,i2,i3,i4,i5,i3) + 
     &         zaa2(i3,i2,i3,i4,i5,i2))*
     &       (sqrtd2345*za(i3,i5) + zaa2(i3,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i3))*
     &       (sqrtd2345*za(i4,i5) + zaa2(i4,i2,i3,i4,i5,i5) + 
     &         zaa2(i5,i2,i3,i4,i5,i4))*
     &       (sqrtd2345*za(i4,i6) + zaa2(i4,i2,i3,i4,i5,i6) + 
     &         zaa2(i6,i2,i3,i4,i5,i4)))/
     &     (32.*sqrtd2345**2*za(i1,i6)*zaa2(i2,i2,i3,i4,i5,i2)*
     &       zaa2(i5,i2,i3,i4,i5,i5)**2*zaa2(i6,i2,i3,i4,i5,i6)))
     &)


      aaaa_NMHV_c(89) = (
     &  (-sqrtd2435 + s(i1,i6) + s(i2,i4) - s(i3,i5))*
     &  ((-(s(i3,i5)*za(i3,i6)*
     &           (-(sqrtd2435*za(i1,i2)) + zaa2(i1,i2,i4,i3,i5,i2) + 
     &             zaa2(i2,i2,i4,i3,i5,i1))*
     &           (-(sqrtd2435*za(i1,i5)) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &             zaa2(i5,i2,i4,i3,i5,i1))*
     &           (-(sqrtd2435*za(i3,i5)) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &             zaa2(i5,i2,i4,i3,i5,i3))*
     &           (-(sqrtd2435*za(i4,i6)) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &              zaa2(i6,i2,i4,i3,i5,i4))**2)/
     &        (64.*sqrtd2435*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &          zaa2(i5,i2,i4,i3,i5,i5)**2*zaa2(i6,i2,i4,i3,i5,i6)**2)
     &        - (s(i3,i5)*(-(sqrtd2435*za(i1,i2)) + 
     &            zaa2(i1,i2,i4,i3,i5,i2) + zaa2(i2,i2,i4,i3,i5,i1))*
     &          zaa2(i3,i2,i4,i3,i5,i3)*
     &          (-(sqrtd2435*za(i1,i5)) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i1))*
     &          (-(sqrtd2435*za(i4,i5)) + zaa2(i4,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i4))*
     &          (-(sqrtd2435*za(i4,i6)) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &            zaa2(i6,i2,i4,i3,i5,i4)))/
     &        (32.*sqrtd2435**2*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &          zaa2(i5,i2,i4,i3,i5,i5)**2*zaa2(i6,i2,i4,i3,i5,i6)))/
     &     (sqrtd2435/2. + (s(i1,i6) - s(i2,i4) - s(i3,i5))/2.) + 
     &    (sqrtd2435/2. + (s(i1,i6) - s(i2,i4) - s(i3,i5))/2.)*
     &     (-((sqrtd2435*za(i1,i2) + zaa2(i1,i2,i4,i3,i5,i2) + 
     &             zaa2(i2,i2,i4,i3,i5,i1))*zaa2(i4,i2,i4,i3,i5,i4)*
     &           (sqrtd2435*za(i1,i5) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &             zaa2(i5,i2,i4,i3,i5,i1))*
     &           (sqrtd2435*za(i3,i5) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &             zaa2(i5,i2,i4,i3,i5,i3))*
     &           (sqrtd2435*za(i3,i6) + zaa2(i3,i2,i4,i3,i5,i6) + 
     &             zaa2(i6,i2,i4,i3,i5,i3)))/
     &        (32.*sqrtd2435**2*s(i2,i4)*za(i1,i6)*
     &          zaa2(i2,i2,i4,i3,i5,i2)*zaa2(i5,i2,i4,i3,i5,i5)**2*
     &          zaa2(i6,i2,i4,i3,i5,i6)) - 
     &       (za(i1,i2)*(-(sqrtd2435*za(i1,i2)) + 
     &            zaa2(i1,i2,i4,i3,i5,i2) + zaa2(i2,i2,i4,i3,i5,i1))*
     &          (-(sqrtd2435*za(i3,i2)) + zaa2(i2,i2,i4,i3,i5,i3) + 
     &            zaa2(i3,i2,i4,i3,i5,i2))*
     &          (-(sqrtd2435*za(i3,i5)) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i3))*
     &          (-(sqrtd2435*za(i4,i5)) + zaa2(i4,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i4))*
     &          (-(sqrtd2435*za(i4,i6)) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &            zaa2(i6,i2,i4,i3,i5,i4)))/
     &        (64.*sqrtd2435*s(i2,i4)*za(i1,i6)*
     &          zaa2(i2,i2,i4,i3,i5,i2)**2*zaa2(i5,i2,i4,i3,i5,i5)**2*
     &          zaa2(i6,i2,i4,i3,i5,i6)) + 
     &       (za(i4,i5)*(sqrtd2435*za(i1,i2) + 
     &            zaa2(i1,i2,i4,i3,i5,i2) + zaa2(i2,i2,i4,i3,i5,i1))*
     &          (sqrtd2435*za(i1,i5) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i1))*
     &          (sqrtd2435*za(i3,i5) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &             zaa2(i5,i2,i4,i3,i5,i3))**2*
     &          (sqrtd2435*za(i4,i6) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &            zaa2(i6,i2,i4,i3,i5,i4)))/
     &        (64.*sqrtd2435*s(i2,i4)*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &          zaa2(i5,i2,i4,i3,i5,i5)**3*zaa2(i6,i2,i4,i3,i5,i6))) - 
     &    (zaa2(i1,i2,i4,i3,i5,i5)*
     &       (-(sqrtd2435*za(i1,i2)) + zaa2(i1,i2,i4,i3,i5,i2) + 
     &         zaa2(i2,i2,i4,i3,i5,i1))*
     &       (-(sqrtd2435*za(i3,i5)) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &          zaa2(i5,i2,i4,i3,i5,i3))**2*
     &       (-(sqrtd2435*za(i4,i5)) + zaa2(i4,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i4))*
     &       (-(sqrtd2435*za(i4,i6)) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &         zaa2(i6,i2,i4,i3,i5,i4)))/
     &     (64.*sqrtd2435*s(i2,i4)*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &       zaa2(i5,i2,i4,i3,i5,i5)**3*zaa2(i6,i2,i4,i3,i5,i6))) + 
     & (-sqrtd2435 + s(i1,i6) - s(i2,i4) + s(i3,i5))*
     &  ((sqrtd2435/2. + (s(i1,i6) - s(i2,i4) - s(i3,i5))/2.)**2*
     &     ((za(i3,i6)*(sqrtd2435*za(i1,i2) + zaa2(i1,i2,i4,i3,i5,i2) + 
     &            zaa2(i2,i2,i4,i3,i5,i1))*
     &          (sqrtd2435*za(i1,i5) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i1))*
     &          (sqrtd2435*za(i3,i5) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i3))*
     &          (sqrtd2435*za(i4,i6) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &             zaa2(i6,i2,i4,i3,i5,i4))**2)/
     &        (64.*sqrtd2435*s(i2,i4)*s(i3,i5)*za(i1,i6)*
     &          zaa2(i2,i2,i4,i3,i5,i2)*zaa2(i5,i2,i4,i3,i5,i5)**2*
     &          zaa2(i6,i2,i4,i3,i5,i6)**2) + 
     &       (za(i3,i5)*(sqrtd2435*za(i1,i2) + 
     &            zaa2(i1,i2,i4,i3,i5,i2) + zaa2(i2,i2,i4,i3,i5,i1))*
     &          (sqrtd2435*za(i1,i5) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i1))*
     &          (sqrtd2435*za(i3,i5) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i3))*
     &          (sqrtd2435*za(i4,i5) + zaa2(i4,i2,i4,i3,i5,i5) + 
     &            zaa2(i5,i2,i4,i3,i5,i4))*
     &          (sqrtd2435*za(i4,i6) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &            zaa2(i6,i2,i4,i3,i5,i4)))/
     &        (64.*sqrtd2435*s(i2,i4)*s(i3,i5)*za(i1,i6)*
     &          zaa2(i2,i2,i4,i3,i5,i2)*zaa2(i5,i2,i4,i3,i5,i5)**3*
     &          zaa2(i6,i2,i4,i3,i5,i6))) - 
     &    ((-(sqrtd2435*za(i1,i2)) + zaa2(i1,i2,i4,i3,i5,i2) + 
     &         zaa2(i2,i2,i4,i3,i5,i1))*zaa2(i4,i2,i4,i3,i5,i4)*
     &       (-(sqrtd2435*za(i1,i5)) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i1))*
     &       (-(sqrtd2435*za(i3,i5)) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i3))*
     &       (-(sqrtd2435*za(i3,i6)) + zaa2(i3,i2,i4,i3,i5,i6) + 
     &         zaa2(i6,i2,i4,i3,i5,i3)))/
     &     (32.*sqrtd2435**2*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &       zaa2(i5,i2,i4,i3,i5,i5)**2*zaa2(i6,i2,i4,i3,i5,i6)) - 
     &    (za(i4,i5)*(-(sqrtd2435*za(i1,i2)) + 
     &         zaa2(i1,i2,i4,i3,i5,i2) + zaa2(i2,i2,i4,i3,i5,i1))*
     &       (-(sqrtd2435*za(i1,i5)) + zaa2(i1,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i1))*
     &       (-(sqrtd2435*za(i3,i5)) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &          zaa2(i5,i2,i4,i3,i5,i3))**2*
     &       (-(sqrtd2435*za(i4,i6)) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &         zaa2(i6,i2,i4,i3,i5,i4)))/
     &     (64.*sqrtd2435*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &       zaa2(i5,i2,i4,i3,i5,i5)**3*zaa2(i6,i2,i4,i3,i5,i6)) + 
     &    (za(i1,i2)*(sqrtd2435*za(i1,i2) + zaa2(i1,i2,i4,i3,i5,i2) + 
     &         zaa2(i2,i2,i4,i3,i5,i1))*
     &       (sqrtd2435*za(i3,i2) + zaa2(i2,i2,i4,i3,i5,i3) + 
     &         zaa2(i3,i2,i4,i3,i5,i2))*
     &       (sqrtd2435*za(i3,i5) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i3))*
     &       (sqrtd2435*za(i4,i5) + zaa2(i4,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i4))*
     &       (sqrtd2435*za(i4,i6) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &         zaa2(i6,i2,i4,i3,i5,i4)))/
     &     (64.*sqrtd2435*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)**2*
     &       zaa2(i5,i2,i4,i3,i5,i5)**2*zaa2(i6,i2,i4,i3,i5,i6)) - 
     &    (zaa2(i1,i2,i4,i3,i5,i1)*
     &       (sqrtd2435*za(i3,i2) + zaa2(i2,i2,i4,i3,i5,i3) + 
     &         zaa2(i3,i2,i4,i3,i5,i2))*
     &       (sqrtd2435*za(i3,i5) + zaa2(i3,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i3))*
     &       (sqrtd2435*za(i4,i5) + zaa2(i4,i2,i4,i3,i5,i5) + 
     &         zaa2(i5,i2,i4,i3,i5,i4))*
     &       (sqrtd2435*za(i4,i6) + zaa2(i4,i2,i4,i3,i5,i6) + 
     &         zaa2(i6,i2,i4,i3,i5,i4)))/
     &     (32.*sqrtd2435**2*za(i1,i6)*zaa2(i2,i2,i4,i3,i5,i2)*
     &       zaa2(i5,i2,i4,i3,i5,i5)**2*zaa2(i6,i2,i4,i3,i5,i6)))
     &)


      return
      end
