!===== T. Dennen, May 2014
!===== Bubble coefficients for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,-)ga(i5,+)ga(i6,+)
      subroutine aaaa_NMHV_b_init(i1,i2,i3,i4,i5,i6,za,zb,aaaa_NMHV_b)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7,i8
      complex(dp):: aaaa_NMHV_b(25), zab2, zaa2, zbb2, zab3, t
      complex(dp):: sqrtd
      complex(dp):: sqrtd1325, sqrtd1326, sqrtd1425, sqrtd1426
      complex(dp):: sqrtd1523, sqrtd1524, sqrtd1536, sqrtd1546
      complex(dp):: sqrtd1623, sqrtd1624, sqrtd1635, sqrtd1645
      complex(dp):: sqrtd2354, sqrtd2364, sqrtd2453, sqrtd2463

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

      sqrtd1325 = sqrtd(i1,i3,i2,i5)
      sqrtd1326 = sqrtd(i1,i3,i2,i6)
      sqrtd1425 = sqrtd(i1,i4,i2,i5)
      sqrtd1426 = sqrtd(i1,i4,i2,i6)
      sqrtd1523 = sqrtd(i1,i5,i2,i3)
      sqrtd1524 = sqrtd(i1,i5,i2,i4)
      sqrtd1536 = sqrtd(i1,i5,i3,i6)
      sqrtd1546 = sqrtd(i1,i5,i4,i6)
      sqrtd1623 = sqrtd(i1,i6,i2,i3)
      sqrtd1624 = sqrtd(i1,i6,i2,i4)
      sqrtd1635 = sqrtd(i1,i6,i3,i5)
      sqrtd1645 = sqrtd(i1,i6,i4,i5)
      sqrtd2354 = sqrtd(i2,i3,i5,i4)
      sqrtd2364 = sqrtd(i2,i3,i6,i4)
      sqrtd2453 = sqrtd(i2,i4,i5,i3)
      sqrtd2463 = sqrtd(i2,i4,i6,i3)


      aaaa_NMHV_b(:) = czip

      aaaa_NMHV_b(6) = (
     &  (za(i1,i3)**2*zab2(i4,i1,i5,i3)**2*zb(i5,i1)**2)/
     &  (t(i1,i3,i5)*za(i2,i6)*zab2(i3,i1,i5,i3)**2*zab2(i6,i1,i5,i3)*
     &    zb(i3,i1)) + (2*za(i2,i4)**2*zab2(i2,i1,i3,i5)*zb(i5,i2))/
     &  (za(i2,i6)**2*zab2(i2,i1,i5,i3)*zab2(i2,i4,i6,i2)*zb(i3,i1)) + 
     & (2*t(i1,i3,i5)*za(i1,i3)*za(i4,i6)*zab2(i4,i1,i5,i3)*zb(i5,i1)*
     &    zb(i5,i3))/
     &  (zab2(i2,i1,i5,i3)*zab2(i3,i1,i5,i3)*zab2(i6,i1,i5,i3)**2*
     &    zab2(i6,i3,i5,i1)) - 
     & (2*za(i1,i3)*zab2(i2,i3,i5,i1)*zab2(i4,i1,i5,i3)*
     &    zab2(i4,i3,i5,i1)*zb(i5,i1)*zb(i5,i3))/
     &  (t(i1,i3,i5)*za(i2,i6)*zab2(i2,i1,i5,i3)*zab2(i3,i1,i5,i3)*
     &    zab2(i6,i3,i5,i1)*zb(i3,i1)**2) - 
     & (2*za(i1,i3)*zab2(i4,i3,i5,i1)**2*zb(i5,i1)*zb(i5,i3))/
     &  (t(i1,i3,i5)*za(i2,i6)*zab2(i1,i3,i5,i1)*zab2(i6,i3,i5,i1)*
     &    zb(i3,i1)**2) + (2*za(i1,i3)*za(i3,i4)*zab2(i4,i1,i5,i3)*
     &    zb(i5,i1)*zb(i5,i3))/
     &  (za(i2,i6)*zab2(i3,i1,i5,i3)**2*zab2(i6,i1,i5,i3)*zb(i3,i1)) + 
     & (2*za(i4,i6)**2*zab2(i2,i1,i3,i5)*zab2(i2,i3,i5,i1)*zb(i6,i5))/
     &  (za(i2,i6)**2*zab2(i2,i1,i5,i3)*zab2(i6,i2,i4,i6)*
     &    zab2(i6,i3,i5,i1)*zb(i3,i1)) + 
     & (2*t(i1,i3,i5)**2*za(i4,i6)**2*zb(i5,i3)*zb(i6,i5))/
     &  (zab2(i2,i1,i5,i3)*zab2(i6,i1,i5,i3)**2*zab2(i6,i2,i4,i6)*
     &    zab2(i6,i3,i5,i1)) - 
     & (t(i1,i3,i5)*za(i4,i6)**2*zb(i6,i5)**2)/
     &  (za(i2,i6)*zab2(i6,i1,i5,i3)*zab2(i6,i2,i4,i6)**2*zb(i3,i1))
     &)


      aaaa_NMHV_b(7) = (
     &  (2*za(i4,i5)**2*zab2(i2,i1,i3,i6)*zab2(i2,i3,i6,i1)*zb(i5,i6))/
     &  (za(i2,i5)**2*zab2(i2,i1,i6,i3)*zab2(i5,i2,i4,i5)*
     &    zab2(i5,i3,i6,i1)*zb(i3,i1)) - 
     & (t(i1,i3,i6)*za(i4,i5)**2*zb(i5,i6)**2)/
     &  (za(i2,i5)*zab2(i5,i1,i6,i3)*zab2(i5,i2,i4,i5)**2*zb(i3,i1)) + 
     & (za(i1,i3)**2*zab2(i4,i1,i6,i3)**2*zb(i6,i1)**2)/
     &  (t(i1,i3,i6)*za(i2,i5)*zab2(i3,i1,i6,i3)**2*zab2(i5,i1,i6,i3)*
     &    zb(i3,i1)) + (2*za(i2,i4)**2*zab2(i2,i1,i3,i6)*zb(i6,i2))/
     &  (za(i2,i5)**2*zab2(i2,i1,i6,i3)*zab2(i2,i4,i5,i2)*zb(i3,i1)) + 
     & (2*t(i1,i3,i6)**2*za(i4,i5)**2*zb(i5,i6)*zb(i6,i3))/
     &  (zab2(i2,i1,i6,i3)*zab2(i5,i1,i6,i3)**2*zab2(i5,i2,i4,i5)*
     &    zab2(i5,i3,i6,i1)) + 
     & (2*t(i1,i3,i6)*za(i1,i3)*za(i4,i5)*zab2(i4,i1,i6,i3)*zb(i6,i1)*
     &    zb(i6,i3))/
     &  (zab2(i2,i1,i6,i3)*zab2(i3,i1,i6,i3)*zab2(i5,i1,i6,i3)**2*
     &    zab2(i5,i3,i6,i1)) - 
     & (2*za(i1,i3)*zab2(i2,i3,i6,i1)*zab2(i4,i1,i6,i3)*
     &    zab2(i4,i3,i6,i1)*zb(i6,i1)*zb(i6,i3))/
     &  (t(i1,i3,i6)*za(i2,i5)*zab2(i2,i1,i6,i3)*zab2(i3,i1,i6,i3)*
     &    zab2(i5,i3,i6,i1)*zb(i3,i1)**2) - 
     & (2*za(i1,i3)*zab2(i4,i3,i6,i1)**2*zb(i6,i1)*zb(i6,i3))/
     &  (t(i1,i3,i6)*za(i2,i5)*zab2(i1,i3,i6,i1)*zab2(i5,i3,i6,i1)*
     &    zb(i3,i1)**2) + (2*za(i1,i3)*za(i3,i4)*zab2(i4,i1,i6,i3)*
     &    zb(i6,i1)*zb(i6,i3))/
     &  (za(i2,i5)*zab2(i3,i1,i6,i3)**2*zab2(i5,i1,i6,i3)*zb(i3,i1))
     &)


      aaaa_NMHV_b(8) = (
     &  (za(i1,i4)**2*zab2(i3,i1,i5,i4)**2*zb(i5,i1)**2)/
     &  (t(i1,i4,i5)*za(i2,i6)*zab2(i4,i1,i5,i4)**2*zab2(i6,i1,i5,i4)*
     &    zb(i4,i1)) + (2*za(i2,i3)**2*zab2(i2,i1,i4,i5)*zb(i5,i2))/
     &  (za(i2,i6)**2*zab2(i2,i1,i5,i4)*zab2(i2,i3,i6,i2)*zb(i4,i1)) + 
     & (2*t(i1,i4,i5)*za(i1,i4)*za(i3,i6)*zab2(i3,i1,i5,i4)*zb(i5,i1)*
     &    zb(i5,i4))/
     &  (zab2(i2,i1,i5,i4)*zab2(i4,i1,i5,i4)*zab2(i6,i1,i5,i4)**2*
     &    zab2(i6,i4,i5,i1)) - 
     & (2*za(i1,i4)*zab2(i3,i4,i5,i1)**2*zb(i5,i1)*zb(i5,i4))/
     &  (t(i1,i4,i5)*za(i2,i6)*zab2(i1,i4,i5,i1)*zab2(i6,i4,i5,i1)*
     &    zb(i4,i1)**2) - (2*za(i1,i4)*zab2(i2,i4,i5,i1)*
     &    zab2(i3,i1,i5,i4)*zab2(i3,i4,i5,i1)*zb(i5,i1)*zb(i5,i4))/
     &  (t(i1,i4,i5)*za(i2,i6)*zab2(i2,i1,i5,i4)*zab2(i4,i1,i5,i4)*
     &    zab2(i6,i4,i5,i1)*zb(i4,i1)**2) + 
     & (2*za(i1,i4)*za(i4,i3)*zab2(i3,i1,i5,i4)*zb(i5,i1)*zb(i5,i4))/
     &  (za(i2,i6)*zab2(i4,i1,i5,i4)**2*zab2(i6,i1,i5,i4)*zb(i4,i1)) + 
     & (2*za(i3,i6)**2*zab2(i2,i1,i4,i5)*zab2(i2,i4,i5,i1)*zb(i6,i5))/
     &  (za(i2,i6)**2*zab2(i2,i1,i5,i4)*zab2(i6,i2,i3,i6)*
     &    zab2(i6,i4,i5,i1)*zb(i4,i1)) + 
     & (2*t(i1,i4,i5)**2*za(i3,i6)**2*zb(i5,i4)*zb(i6,i5))/
     &  (zab2(i2,i1,i5,i4)*zab2(i6,i1,i5,i4)**2*zab2(i6,i2,i3,i6)*
     &    zab2(i6,i4,i5,i1)) - 
     & (t(i1,i4,i5)*za(i3,i6)**2*zb(i6,i5)**2)/
     &  (za(i2,i6)*zab2(i6,i1,i5,i4)*zab2(i6,i2,i3,i6)**2*zb(i4,i1))
     &)


      aaaa_NMHV_b(9) = (
     &  (2*za(i3,i5)**2*zab2(i2,i1,i4,i6)*zab2(i2,i4,i6,i1)*zb(i5,i6))/
     &  (za(i2,i5)**2*zab2(i2,i1,i6,i4)*zab2(i5,i2,i3,i5)*
     &    zab2(i5,i4,i6,i1)*zb(i4,i1)) - 
     & (t(i1,i4,i6)*za(i3,i5)**2*zb(i5,i6)**2)/
     &  (za(i2,i5)*zab2(i5,i1,i6,i4)*zab2(i5,i2,i3,i5)**2*zb(i4,i1)) + 
     & (za(i1,i4)**2*zab2(i3,i1,i6,i4)**2*zb(i6,i1)**2)/
     &  (t(i1,i4,i6)*za(i2,i5)*zab2(i4,i1,i6,i4)**2*zab2(i5,i1,i6,i4)*
     &    zb(i4,i1)) + (2*za(i2,i3)**2*zab2(i2,i1,i4,i6)*zb(i6,i2))/
     &  (za(i2,i5)**2*zab2(i2,i1,i6,i4)*zab2(i2,i3,i5,i2)*zb(i4,i1)) + 
     & (2*t(i1,i4,i6)**2*za(i3,i5)**2*zb(i5,i6)*zb(i6,i4))/
     &  (zab2(i2,i1,i6,i4)*zab2(i5,i1,i6,i4)**2*zab2(i5,i2,i3,i5)*
     &    zab2(i5,i4,i6,i1)) + 
     & (2*t(i1,i4,i6)*za(i1,i4)*za(i3,i5)*zab2(i3,i1,i6,i4)*zb(i6,i1)*
     &    zb(i6,i4))/
     &  (zab2(i2,i1,i6,i4)*zab2(i4,i1,i6,i4)*zab2(i5,i1,i6,i4)**2*
     &    zab2(i5,i4,i6,i1)) - 
     & (2*za(i1,i4)*zab2(i3,i4,i6,i1)**2*zb(i6,i1)*zb(i6,i4))/
     &  (t(i1,i4,i6)*za(i2,i5)*zab2(i1,i4,i6,i1)*zab2(i5,i4,i6,i1)*
     &    zb(i4,i1)**2) - (2*za(i1,i4)*zab2(i2,i4,i6,i1)*
     &    zab2(i3,i1,i6,i4)*zab2(i3,i4,i6,i1)*zb(i6,i1)*zb(i6,i4))/
     &  (t(i1,i4,i6)*za(i2,i5)*zab2(i2,i1,i6,i4)*zab2(i4,i1,i6,i4)*
     &    zab2(i5,i4,i6,i1)*zb(i4,i1)**2) + 
     & (2*za(i1,i4)*za(i4,i3)*zab2(i3,i1,i6,i4)*zb(i6,i1)*zb(i6,i4))/
     &  (za(i2,i5)*zab2(i4,i1,i6,i4)**2*zab2(i5,i1,i6,i4)*zb(i4,i1))
     &)


      aaaa_NMHV_b(10) = (
     &  (2*zab2(i1,i5,i6,i3)**3*zab2(i3,i2,i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i3,i2,i4,i3)*zab2(i5,i1,i6,i3)*
     &    zab2(i6,i1,i5,i3)*zb(i4,i3)**2) + 
     & (2*zab2(i1,i5,i6,i4)**3*zab2(i4,i2,i3,i2))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i4,i2,i3,i4)*zab2(i5,i1,i6,i4)*
     &    zab2(i6,i1,i5,i4)*zb(i4,i3)**2) - 
     & (2*zab2(i1,i5,i6,i3)**2*zab2(i1,i5,i6,i4)*zab2(i3,i2,i4,i2)*
     &    zb(i3,i2))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i3,i2,i4,i3)*zab2(i5,i1,i6,i3)*
     &    zab2(i6,i1,i5,i3)*zb(i4,i2)*zb(i4,i3)**2) - 
     & (2*zab2(i1,i5,i6,i3)*zab2(i1,i5,i6,i4)**2*zab2(i4,i2,i3,i2)*
     &    zb(i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i4,i2,i3,i4)*zab2(i5,i1,i6,i4)*
     &    zab2(i6,i1,i5,i4)*zb(i3,i2)*zb(i4,i3)**2) - 
     & (zab2(i1,i5,i6,i4)**3*zab2(i4,i2,i3,i2)**2)/
     &  (za(i1,i5)*za(i1,i6)*zab2(i4,i2,i3,i4)**2*zab2(i5,i1,i6,i4)*
     &    zab2(i6,i1,i5,i4)*zb(i3,i2)*zb(i4,i3)) + 
     & (zab2(i1,i5,i6,i3)**3*zab2(i3,i2,i4,i2)**2)/
     &  (za(i1,i5)*za(i1,i6)*zab2(i3,i2,i4,i3)**2*zab2(i5,i1,i6,i3)*
     &    zab2(i6,i1,i5,i3)*zb(i4,i2)*zb(i4,i3)) - 
     & (2*t(i1,i5,i6)*zab2(i1,i5,i6,i3)**2*zab2(i3,i2,i4,i2)*zb(i3,i2))/
     &  (za(i1,i5)*zab2(i3,i2,i4,i3)*zab2(i5,i1,i6,i3)*
     &    zab2(i6,i1,i5,i3)**2*zb(i4,i2)*zb(i4,i3)) - 
     & (2*t(i1,i5,i6)*zab2(i1,i5,i6,i3)**2*zab2(i3,i2,i4,i2)*zb(i3,i2))/
     &  (za(i1,i6)*zab2(i3,i2,i4,i3)*zab2(i5,i1,i6,i3)**2*
     &    zab2(i6,i1,i5,i3)*zb(i4,i2)*zb(i4,i3)) + 
     & (2*t(i1,i5,i6)*zab2(i1,i5,i6,i4)**2*zab2(i4,i2,i3,i2)*zb(i4,i2))/
     &  (za(i1,i5)*zab2(i4,i2,i3,i4)*zab2(i5,i1,i6,i4)*
     &    zab2(i6,i1,i5,i4)**2*zb(i3,i2)*zb(i4,i3)) + 
     & (2*t(i1,i5,i6)*zab2(i1,i5,i6,i4)**2*zab2(i4,i2,i3,i2)*zb(i4,i2))/
     &  (za(i1,i6)*zab2(i4,i2,i3,i4)*zab2(i5,i1,i6,i4)**2*
     &    zab2(i6,i1,i5,i4)*zb(i3,i2)*zb(i4,i3)) - 
     & (2*t(i1,i5,i6)*za(i1,i5)*zab2(i1,i5,i6,i4)*zab2(i5,i1,i6,i2)**2*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i5,i6)*zab2(i5,i1,i6,i3)*zab2(i5,i1,i6,i4)**2*
     &    zab2(i5,i1,i6,i5)*zb(i3,i2)*zb(i4,i2)) - 
     & (2*t(i1,i5,i6)*za(i1,i5)*zab2(i1,i5,i6,i3)*zab2(i5,i1,i6,i2)**2*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i5,i6)*zab2(i5,i1,i6,i3)**2*zab2(i5,i1,i6,i4)*
     &    zab2(i5,i1,i6,i5)*zb(i3,i2)*zb(i4,i2)) + 
     & (t(i1,i5,i6)**2*za(i1,i5)**2*zab2(i5,i1,i6,i2)*zb(i5,i2)**2)/
     &  (za(i1,i6)*za(i5,i6)*zab2(i5,i1,i6,i3)*zab2(i5,i1,i6,i4)*
     &    zab2(i5,i1,i6,i5)**2*zb(i3,i2)*zb(i4,i2)) + 
     & (2*t(i1,i5,i6)*za(i1,i6)*zab2(i1,i5,i6,i4)*zab2(i6,i1,i5,i2)**2*
     &    zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)*zab2(i6,i1,i5,i4)**2*
     &    zab2(i6,i1,i5,i6)*zb(i3,i2)*zb(i4,i2)) + 
     & (2*t(i1,i5,i6)*za(i1,i6)*zab2(i1,i5,i6,i3)*zab2(i6,i1,i5,i2)**2*
     &    zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)**2*zab2(i6,i1,i5,i4)*
     &    zab2(i6,i1,i5,i6)*zb(i3,i2)*zb(i4,i2)) - 
     & (t(i1,i5,i6)**2*za(i1,i6)**2*zab2(i6,i1,i5,i2)*zb(i6,i2)**2)/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i1,i5,i3)*zab2(i6,i1,i5,i4)*
     &    zab2(i6,i1,i5,i6)**2*zb(i3,i2)*zb(i4,i2))
     &)


      aaaa_NMHV_b(13) = (
     &  (2*zab2(i1,i4,i6,i3)*zab2(i3,i4,i6,i5)*zab2(i4,i4,i6,i3)**2)/
     &  (za(i1,i5)*zab2(i2,i4,i6,i3)*zab2(i3,i4,i6,i3)*
     &    zab2(i6,i4,i6,i3)**2*zb(i3,i1)) + 
     & (2*zab2(i1,i4,i6,i3)*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2)/
     &  (za(i1,i5)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i3)**2*zb(i3,i1)) + 
     & (2*zab2(i3,i4,i6,i5)*zab2(i4,i4,i6,i3)**2*zab2(i6,i4,i6,i2))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i6,i4,i6,i1)*
     &    zab2(i6,i4,i6,i3)**2*zb(i3,i2)) + 
     & (2*zab2(i2,i4,i6,i3)*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*
     &    zab2(i6,i4,i6,i2))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)**2*zb(i3,i2)) - 
     & (2*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*zab2(i5,i4,i6,i2))/
     &  (za(i5,i6)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i1)*
     &    zab2(i5,i4,i6,i3)*zab2(i6,i4,i6,i3)*zb(i3,i2)) + 
     & (2*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*zab2(i6,i4,i6,i2))/
     &  (za(i5,i6)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)*zb(i3,i2)) + 
     & (2*zab2(i1,i4,i6,i5)*zab2(i4,i4,i6,i1)**2*zb(i2,i1))/
     &  (za(i2,i5)*zab2(i1,i4,i6,i1)*zab2(i6,i4,i6,i1)**2*zb(i3,i1)*
     &    zb(i3,i2)) + (2*zab2(i1,i4,i6,i2)*zab2(i2,i4,i6,i1)*
     &    zab2(i4,i4,i6,i1)**2*zb(i2,i1))/
     &  (za(i2,i5)*zab2(i1,i4,i6,i1)*zab2(i5,i4,i6,i1)*
     &    zab2(i6,i4,i6,i1)**2*zb(i3,i1)*zb(i3,i2)) - 
     & (2*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*zb(i2,i1))/
     &  (za(i5,i6)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i1)*zb(i3,i1)*zb(i3,i2)) + 
     & (2*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*zb(i2,i1))/
     &  (za(i5,i6)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i1)*
     &    zab2(i6,i4,i6,i3)*zb(i3,i1)*zb(i3,i2)) - 
     & (2*zab2(i3,i4,i6,i5)*zab2(i4,i4,i6,i3)**2*zb(i2,i1))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i6,i4,i6,i1)*
     &    zab2(i6,i4,i6,i3)*zb(i3,i1)*zb(i3,i2)) - 
     & (2*zab2(i2,i4,i6,i3)*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*
     &    zb(i2,i1))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)*zb(i3,i1)*zb(i3,i2)) - 
     & (2*za(i1,i3)*za(i2,i4)**2*zab2(i1,i4,i6,i2)*zb(i5,i1))/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i2,i4,i6,i2)*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i4,i6)*zab2(i1,i4,i6,i6)*zab2(i4,i4,i6,i6)*
     &    zb(i5,i1))/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)*zab2(i6,i4,i6,i6)**2*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i2,i4)*za(i4,i6)*zab2(i1,i4,i6,i6)*
     &    zb(i5,i1))/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i6,i4,i6,i6)*zb(i3,i1))
     &   + (2*s(i4,i6)*za(i1,i2)*za(i2,i4)**2*zb(i5,i2))/
     &  (za(i1,i5)*za(i2,i6)**2*zab2(i2,i4,i6,i2)*zab2(i2,i4,i6,i3)*
     &    zb(i3,i1)) - (2*s(i4,i6)*za(i4,i5)**2*zb(i5,i2))/
     &  (za(i5,i6)**2*zab2(i5,i4,i6,i3)*zab2(i5,i4,i6,i5)*zb(i3,i1)) - 
     & (2*s(i4,i6)*za(i2,i5)*za(i4,i5)**2*zab2(i5,i4,i6,i2)*zb(i5,i2))/
     &  (za(i5,i6)**2*zaa2(i5,i1,i3,i2,i5,i5)*zab2(i5,i4,i6,i3)*
     &    zab2(i5,i4,i6,i5)*zb(i3,i1)) + 
     & (2*s(i4,i6)*za(i1,i5)*za(i4,i5)**2*zab2(i5,i4,i6,i2)*zb(i5,i2))/
     &  (za(i5,i6)**2*zaa2(i5,i1,i5,i2,i3,i5)*zab2(i5,i4,i6,i3)*
     &    zab2(i5,i4,i6,i5)*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i4,i5)**2*zab2(i5,i4,i6,i2)*zb(i5,i2))/
     &  (za(i5,i6)**2*zab2(i5,i4,i6,i1)*zab2(i5,i4,i6,i3)*
     &    zab2(i5,i4,i6,i5)*zb(i3,i2)) + 
     & (2*za(i2,i3)*zab2(i1,i4,i6,i2)*zab2(i4,i4,i6,i1)**2*zb(i5,i2))/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i1,i4,i6,i1)*zab2(i6,i4,i6,i1)**2*
     &    zb(i3,i2)) - (2*za(i3,i5)*zab2(i1,i4,i6,i5)*
     &    zab2(i4,i4,i6,i1)**2*zb(i5,i2))/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i1,i4,i6,i1)*zab2(i6,i4,i6,i1)**2*
     &    zb(i3,i2)) - (2*za(i1,i3)*za(i2,i4)**2*zab2(i3,i4,i6,i2)*
     &    zb(i5,i3))/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i2,i4,i6,i2)*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i4,i6)*zab2(i3,i4,i6,i6)*zab2(i4,i4,i6,i6)*
     &    zb(i5,i3))/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)*zab2(i6,i4,i6,i6)**2*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i2,i4)*za(i4,i6)*zab2(i3,i4,i6,i6)*
     &    zb(i5,i3))/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i6,i4,i6,i6)*zb(i3,i1))
     &   - (2*s(i4,i6)*za(i1,i6)*za(i4,i6)*zab2(i4,i4,i6,i6)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)**2*
     &    zb(i3,i1)) - (2*s(i4,i6)*za(i1,i6)*za(i4,i6)*
     &    zab2(i4,i4,i6,i3)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i4,i6,i3)**2*zab2(i6,i4,i6,i6)*
     &    zb(i3,i1)) - (2*s(i4,i6)*za(i2,i6)**2*za(i4,i6)*
     &    zab2(i4,i4,i6,i3)*zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i3,i2,i5,i6)*
     &    zab2(i6,i4,i6,i3)**2*zab2(i6,i4,i6,i6)*zb(i3,i1)) + 
     & (2*s(i4,i6)*za(i4,i6)**2*zb(i6,i2))/
     &  (za(i5,i6)**2*zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*zb(i3,i1)) + 
     & (2*s(i4,i6)*za(i2,i6)*za(i4,i6)**2*zaa2(i2,i2,i5,i1,i3,i6)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i3,i2,i5,i6)**2*
     &    zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*zb(i3,i1)) - 
     & (2*s(i4,i6)*za(i2,i6)**2*za(i4,i6)*zaa2(i4,i1,i3,i2,i5,i6)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i3,i2,i5,i6)**2*
     &    zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*zb(i3,i1)) + 
     & (2*s(i4,i6)*za(i2,i6)*za(i4,i6)**2*zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i5,i6)**2*zaa2(i6,i1,i3,i2,i5,i6)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) - 
     & (4*s(i4,i6)*za(i2,i6)*za(i4,i6)*zab2(i4,i4,i6,i3)*
     &    zab2(i6,i4,i6,i5)*zb(i6,i2))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)*zab2(i6,i4,i6,i3)**2*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) + 
     & (4*s(i4,i6)*za(i4,i6)**2*zaa2(i2,i2,i5,i1,i3,i6)*
     &    zab2(i6,i4,i6,i5)*zb(i6,i2))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) - 
     & (4*s(i4,i6)*za(i2,i6)*za(i4,i6)*zaa2(i4,i1,i3,i2,i5,i6)*
     &    zab2(i6,i4,i6,i5)*zb(i6,i2))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) - 
     & (2*s(i4,i6)*za(i2,i6)*za(i4,i6)*zab2(i4,i4,i6,i6)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)**2*zb(i3,i2)) + 
     & (2*s(i4,i6)*za(i1,i6)**2*za(i4,i6)*zab2(i4,i4,i6,i3)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i3,i6)*
     &    zab2(i6,i4,i6,i3)**2*zab2(i6,i4,i6,i6)*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i2,i6)*za(i4,i6)*zab2(i4,i4,i6,i3)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)**2*
     &    zab2(i6,i4,i6,i6)*zb(i3,i2)) + 
     & (2*s(i4,i6)*za(i4,i6)**2*zab2(i5,i4,i6,i2)*zb(i6,i2))/
     &  (za(i5,i6)**2*zab2(i5,i4,i6,i1)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i1,i6)*za(i4,i6)**2*zaa2(i1,i2,i3,i1,i5,i6)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i3,i6)**2*
     &    zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*zb(i3,i2)) + 
     & (2*s(i4,i6)*za(i1,i6)**2*za(i4,i6)*zaa2(i4,i1,i5,i2,i3,i6)*
     &    zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i3,i6)**2*
     &    zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i1,i6)*za(i4,i6)**2*zab2(i6,i4,i6,i2)*zb(i6,i2))/
     &  (za(i5,i6)**2*zaa2(i6,i1,i5,i2,i3,i6)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i2)) - 
     & (2*s(i4,i6)**2*za(i2,i6)*za(i4,i6)**2*zb(i2,i1)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i4,i6,i1)**2*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i2)) + 
     & (2*s(i4,i6)**2*za(i4,i6)**2*zb(i2,i1)*zb(i6,i2))/
     &  (za(i5,i6)*zab2(i5,i4,i6,i1)*zab2(i6,i4,i6,i1)*
     &    zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i2,i3)*za(i4,i6)*zab2(i4,i4,i6,i6)*zb(i5,i2)*
     &    zb(i6,i2))/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i6)**2*
     &    zb(i3,i2)) - (2*s(i4,i6)*za(i2,i3)*za(i4,i6)*
     &    zab2(i4,i4,i6,i1)*zb(i5,i2)*zb(i6,i2))/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i6,i4,i6,i1)**2*zab2(i6,i4,i6,i6)*
     &    zb(i3,i2)) + (2*s(i4,i6)**2*za(i2,i6)**2*za(i4,i6)**2*
     &    zb(i6,i2)**2)/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i3,i2,i5,i6)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)**2*zb(i3,i1)) - 
     & (2*s(i4,i6)**2*za(i1,i6)**2*za(i4,i6)**2*zb(i6,i2)**2)/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i3,i6)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)**2*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i1,i6)*za(i4,i6)*zab2(i4,i4,i6,i6)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)**2*
     &    zb(i3,i1)) - (2*s(i4,i6)*za(i1,i6)*za(i4,i6)*
     &    zab2(i4,i4,i6,i3)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i6,i4,i6,i3)**2*zab2(i6,i4,i6,i6)*
     &    zb(i3,i1)) + (2*s(i4,i6)*za(i1,i2)*za(i4,i6)**2*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)**2*zab2(i6,i4,i6,i3)*zab2(i6,i4,i6,i6)*
     &    zb(i3,i1)) - (2*s(i4,i6)*za(i4,i6)**2*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i5)*zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)*zab2(i6,i4,i6,i3)**2*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) - 
     & (2*s(i4,i6)*za(i4,i6)*za(i5,i6)*zaa2(i4,i1,i3,i2,i5,i6)*
     &    zab2(i6,i4,i6,i5)*zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) + 
     & (2*s(i4,i6)*za(i4,i6)*za(i5,i6)*zaa2(i4,i2,i5,i1,i3,i6)*
     &    zab2(i6,i4,i6,i5)*zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i1)) - 
     & (2*s(i4,i6)*za(i4,i6)*zab2(i4,i4,i6,i6)*zab2(i6,i4,i6,i2)*
     &    zb(i6,i5))/
     &  (za(i2,i5)*zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)**2*zb(i3,i2)) - 
     & (2*s(i4,i6)*za(i4,i6)*zab2(i4,i4,i6,i3)*zab2(i6,i4,i6,i2)*
     &    zb(i6,i5))/
     &  (za(i2,i5)*zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i3)**2*
     &    zab2(i6,i4,i6,i6)*zb(i3,i2)) - 
     & (2*s(i4,i6)**2*za(i4,i6)**2*zb(i2,i1)*zb(i6,i5))/
     &  (za(i2,i5)*zab2(i6,i4,i6,i1)**2*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)*zb(i3,i2)) + 
     & (2*s(i4,i6)*za(i3,i5)*za(i4,i6)*zab2(i4,i4,i6,i6)*zb(i5,i2)*
     &    zb(i6,i5))/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i6,i4,i6,i1)*zab2(i6,i4,i6,i6)**2*
     &    zb(i3,i2)) + (2*s(i4,i6)*za(i3,i5)*za(i4,i6)*
     &    zab2(i4,i4,i6,i1)*zb(i5,i2)*zb(i6,i5))/
     &  (t(i2,i3,i5)*za(i2,i5)*zab2(i6,i4,i6,i1)**2*zab2(i6,i4,i6,i6)*
     &    zb(i3,i2)) + (4*s(i4,i6)**2*za(i2,i6)*za(i4,i6)**2*zb(i6,i2)*
     &    zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)**2*zb(i3,i1)) + 
     & (2*s(i4,i6)**2*za(i4,i6)**2*za(i5,i6)*zb(i6,i5)**2)/
     &  (za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)*zab2(i6,i4,i6,i3)*
     &    zab2(i6,i4,i6,i6)**2*zb(i3,i1)) - 
     & (2*zab2(i1,i4,i6,i3)**2*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2)/
     &  (za(i1,i5)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i3)**2*zbb2(i3,i2,i3,i1,i5,i3)) + 
     & ((sqrtd1523*za(i1,i5) + zaa2(i1,i1,i5,i2,i3,i5) - 
     &      zaa2(i1,i2,i3,i1,i5,i5))*
     &    (sqrtd1523*za(i1,i6) + zaa2(i1,i1,i5,i2,i3,i6) - 
     &      zaa2(i1,i2,i3,i1,i5,i6))*
     &    (sqrtd1523*za(i4,i6) + zaa2(i4,i1,i5,i2,i3,i6) - 
     &      zaa2(i4,i2,i3,i1,i5,i6))*
     &    (sqrtd1523*zab2(i4,i4,i6,i2) + 
     &      zab3(i4,i1,i5,i2,i3,i4,i6,i2) - 
     &      zab3(i4,i2,i3,i1,i5,i4,i6,i2))*
     &    (sqrtd1523*zb(i3,i2) + zbb2(i3,i1,i5,i2,i3,i2) - 
     &      zbb2(i3,i2,i3,i1,i5,i2)))/
     &  (16.*sqrtd1523**2*za(i1,i5)*zaa2(i5,i1,i5,i2,i3,i5)*
     &    zaa2(i6,i1,i5,i2,i3,i6)**2*zb(i3,i2)*zbb2(i3,i2,i3,i1,i5,i3))
     &  - ((sqrtd1523*za(i1,i5) - zaa2(i1,i1,i5,i2,i3,i5) + 
     &      zaa2(i1,i2,i3,i1,i5,i5))*
     &    (sqrtd1523*za(i1,i6) - zaa2(i1,i1,i5,i2,i3,i6) + 
     &      zaa2(i1,i2,i3,i1,i5,i6))*
     &    (sqrtd1523*za(i4,i6) - zaa2(i4,i1,i5,i2,i3,i6) + 
     &      zaa2(i4,i2,i3,i1,i5,i6))*
     &    (sqrtd1523*zab2(i4,i4,i6,i2) - 
     &      zab3(i4,i1,i5,i2,i3,i4,i6,i2) + 
     &      zab3(i4,i2,i3,i1,i5,i4,i6,i2))*
     &    (sqrtd1523*zb(i3,i2) - zbb2(i3,i1,i5,i2,i3,i2) + 
     &      zbb2(i3,i2,i3,i1,i5,i2)))/
     &  (16.*sqrtd1523**2*za(i1,i5)*zaa2(i5,i1,i5,i2,i3,i5)*
     &    zaa2(i6,i1,i5,i2,i3,i6)**2*zb(i3,i2)*zbb2(i3,i2,i3,i1,i5,i3))
     &  + (2*zab2(i2,i4,i6,i3)**2*zab2(i3,i4,i6,i2)*
     &    zab2(i4,i4,i6,i3)**2*zb(i3,i2))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i5,i4,i6,i3)*
     &    zab2(i6,i4,i6,i3)**2*zb(i3,i1)*zbb2(i3,i2,i5,i1,i3,i3)) - 
     & (4*zab2(i2,i4,i6,i3)*zab2(i3,i4,i6,i2)*zab2(i4,i4,i6,i3)**2*
     &    zb(i5,i3))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i6,i4,i6,i3)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i5,i1,i3,i3)) - 
     & (2*zab2(i3,i4,i6,i5)*zab2(i4,i4,i6,i3)**2*zab2(i5,i4,i6,i3)*
     &    zb(i5,i3))/
     &  (za(i2,i5)*zab2(i3,i4,i6,i3)*zab2(i6,i4,i6,i3)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i5,i1,i3,i3)) - 
     & ((sqrtd1325*za(i2,i5) + zaa2(i2,i1,i3,i2,i5,i5) - 
     &      zaa2(i2,i2,i5,i1,i3,i5))*
     &    (sqrtd1325*za(i2,i6) + zaa2(i2,i1,i3,i2,i5,i6) - 
     &      zaa2(i2,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*za(i4,i6) + zaa2(i4,i1,i3,i2,i5,i6) - 
     &      zaa2(i4,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*zab2(i4,i4,i6,i2) + 
     &      zab3(i4,i1,i3,i2,i5,i4,i6,i2) - 
     &      zab3(i4,i2,i5,i1,i3,i4,i6,i2))*
     &    (sqrtd1325*zb(i3,i2) + zbb2(i3,i1,i3,i2,i5,i2) - 
     &      zbb2(i3,i2,i5,i1,i3,i2)))/
     &  (16.*sqrtd1325**2*za(i2,i5)*zaa2(i5,i1,i3,i2,i5,i5)*
     &    zaa2(i6,i1,i3,i2,i5,i6)**2*zb(i3,i1)*zbb2(i3,i2,i5,i1,i3,i3))
     &  + ((sqrtd1325*za(i2,i5) - zaa2(i2,i1,i3,i2,i5,i5) + 
     &      zaa2(i2,i2,i5,i1,i3,i5))*
     &    (sqrtd1325*za(i2,i6) - zaa2(i2,i1,i3,i2,i5,i6) + 
     &      zaa2(i2,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*za(i4,i6) - zaa2(i4,i1,i3,i2,i5,i6) + 
     &      zaa2(i4,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*zab2(i4,i4,i6,i2) - 
     &      zab3(i4,i1,i3,i2,i5,i4,i6,i2) + 
     &      zab3(i4,i2,i5,i1,i3,i4,i6,i2))*
     &    (sqrtd1325*zb(i3,i2) - zbb2(i3,i1,i3,i2,i5,i2) + 
     &      zbb2(i3,i2,i5,i1,i3,i2)))/
     &  (16.*sqrtd1325**2*za(i2,i5)*zaa2(i5,i1,i3,i2,i5,i5)*
     &    zaa2(i6,i1,i3,i2,i5,i6)**2*zb(i3,i1)*zbb2(i3,i2,i5,i1,i3,i3))
     &  + ((sqrtd1325*za(i2,i6) - zaa2(i2,i1,i3,i2,i5,i6) + 
     &      zaa2(i2,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*za(i4,i6) - zaa2(i4,i1,i3,i2,i5,i6) + 
     &      zaa2(i4,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*zab2(i4,i4,i6,i2) - 
     &      zab3(i4,i1,i3,i2,i5,i4,i6,i2) + 
     &      zab3(i4,i2,i5,i1,i3,i4,i6,i2))*
     &    (sqrtd1325*zb(i5,i3) + zbb2(i5,i1,i3,i2,i5,i3) - 
     &      zbb2(i5,i2,i5,i1,i3,i3)))/
     &  (4.*sqrtd1325**2*za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i5,i1,i3,i3)) + 
     & ((sqrtd1325*za(i4,i6) - zaa2(i4,i1,i3,i2,i5,i6) + 
     &       zaa2(i4,i2,i5,i1,i3,i6))**2*
     &    (sqrtd1325*zab2(i5,i4,i6,i5) - 
     &      zab3(i5,i1,i3,i2,i5,i4,i6,i5) + 
     &      zab3(i5,i2,i5,i1,i3,i4,i6,i5))*
     &    (sqrtd1325*zb(i5,i3) + zbb2(i5,i1,i3,i2,i5,i3) - 
     &      zbb2(i5,i2,i5,i1,i3,i3)))/
     &  (8.*sqrtd1325**2*za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i5,i1,i3,i3)) + 
     & ((sqrtd1325*za(i2,i6) + zaa2(i2,i1,i3,i2,i5,i6) - 
     &      zaa2(i2,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*za(i4,i6) + zaa2(i4,i1,i3,i2,i5,i6) - 
     &      zaa2(i4,i2,i5,i1,i3,i6))*
     &    (sqrtd1325*zab2(i4,i4,i6,i2) + 
     &      zab3(i4,i1,i3,i2,i5,i4,i6,i2) - 
     &      zab3(i4,i2,i5,i1,i3,i4,i6,i2))*
     &    (sqrtd1325*zb(i5,i3) - zbb2(i5,i1,i3,i2,i5,i3) + 
     &      zbb2(i5,i2,i5,i1,i3,i3)))/
     &  (4.*sqrtd1325**2*za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i5,i1,i3,i3)) + 
     & ((sqrtd1325*za(i4,i6) + zaa2(i4,i1,i3,i2,i5,i6) - 
     &       zaa2(i4,i2,i5,i1,i3,i6))**2*
     &    (sqrtd1325*zab2(i5,i4,i6,i5) + 
     &      zab3(i5,i1,i3,i2,i5,i4,i6,i5) - 
     &      zab3(i5,i2,i5,i1,i3,i4,i6,i5))*
     &    (sqrtd1325*zb(i5,i3) - zbb2(i5,i1,i3,i2,i5,i3) + 
     &      zbb2(i5,i2,i5,i1,i3,i3)))/
     &  (8.*sqrtd1325**2*za(i2,i5)*zaa2(i6,i1,i3,i2,i5,i6)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i5,i1,i3,i3))
     &)


      aaaa_NMHV_b(14) = (
     &  (2*zab2(i1,i4,i5,i3)*zab2(i3,i4,i5,i6)*zab2(i4,i4,i5,i3)**2)/
     &  (za(i1,i6)*zab2(i2,i4,i5,i3)*zab2(i3,i4,i5,i3)*
     &    zab2(i5,i4,i5,i3)**2*zb(i3,i1)) + 
     & (2*zab2(i1,i4,i5,i3)*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2)/
     &  (za(i1,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i6,i4,i5,i3)*zb(i3,i1)) + 
     & (2*zab2(i3,i4,i5,i6)*zab2(i4,i4,i5,i3)**2*zab2(i5,i4,i5,i2))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i1)*
     &    zab2(i5,i4,i5,i3)**2*zb(i3,i2)) + 
     & (2*zab2(i2,i4,i5,i3)*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*
     &    zab2(i5,i4,i5,i2))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i1)*
     &    zab2(i5,i4,i5,i3)**2*zab2(i6,i4,i5,i3)*zb(i3,i2)) + 
     & (2*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*zab2(i5,i4,i5,i2))/
     &  (za(i6,i5)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i1)*
     &    zab2(i5,i4,i5,i3)*zab2(i6,i4,i5,i3)*zb(i3,i2)) - 
     & (2*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*zab2(i6,i4,i5,i2))/
     &  (za(i6,i5)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)*
     &    zab2(i6,i4,i5,i1)*zab2(i6,i4,i5,i3)*zb(i3,i2)) + 
     & (2*zab2(i1,i4,i5,i6)*zab2(i4,i4,i5,i1)**2*zb(i2,i1))/
     &  (za(i2,i6)*zab2(i1,i4,i5,i1)*zab2(i5,i4,i5,i1)**2*zb(i3,i1)*
     &    zb(i3,i2)) - (2*zab2(i3,i4,i5,i6)*zab2(i4,i4,i5,i3)**2*
     &    zb(i2,i1))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i1)*
     &    zab2(i5,i4,i5,i3)*zb(i3,i1)*zb(i3,i2)) + 
     & (2*zab2(i1,i4,i5,i2)*zab2(i2,i4,i5,i1)*zab2(i4,i4,i5,i1)**2*
     &    zb(i2,i1))/
     &  (za(i2,i6)*zab2(i1,i4,i5,i1)*zab2(i5,i4,i5,i1)**2*
     &    zab2(i6,i4,i5,i1)*zb(i3,i1)*zb(i3,i2)) + 
     & (2*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*zb(i2,i1))/
     &  (za(i6,i5)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)*
     &    zab2(i6,i4,i5,i1)*zb(i3,i1)*zb(i3,i2)) - 
     & (2*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*zb(i2,i1))/
     &  (za(i6,i5)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i1)*
     &    zab2(i6,i4,i5,i3)*zb(i3,i1)*zb(i3,i2)) - 
     & (2*zab2(i2,i4,i5,i3)*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*
     &    zb(i2,i1))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i1)*
     &    zab2(i5,i4,i5,i3)*zab2(i6,i4,i5,i3)*zb(i3,i1)*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i1,i5)*za(i4,i5)*zab2(i4,i4,i5,i5)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)**2*
     &    zb(i3,i1)) - (2*s(i4,i5)*za(i1,i5)*za(i4,i5)*
     &    zab2(i4,i4,i5,i3)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i4,i5,i3)**2*zab2(i5,i4,i5,i5)*
     &    zb(i3,i1)) - (2*s(i4,i5)*za(i2,i5)**2*za(i4,i5)*
     &    zab2(i4,i4,i5,i3)*zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i3,i2,i6,i5)*
     &    zab2(i5,i4,i5,i3)**2*zab2(i5,i4,i5,i5)*zb(i3,i1)) + 
     & (2*s(i4,i5)*za(i4,i5)**2*zb(i5,i2))/
     &  (za(i6,i5)**2*zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*zb(i3,i1)) + 
     & (2*s(i4,i5)*za(i2,i5)*za(i4,i5)**2*zaa2(i2,i2,i6,i1,i3,i5)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i3,i2,i6,i5)**2*
     &    zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*zb(i3,i1)) - 
     & (2*s(i4,i5)*za(i2,i5)**2*za(i4,i5)*zaa2(i4,i1,i3,i2,i6,i5)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i3,i2,i6,i5)**2*
     &    zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*zb(i3,i1)) + 
     & (2*s(i4,i5)*za(i2,i5)*za(i4,i5)**2*zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i6,i5)**2*zaa2(i5,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) - 
     & (4*s(i4,i5)*za(i2,i5)*za(i4,i5)*zab2(i4,i4,i5,i3)*
     &    zab2(i5,i4,i5,i6)*zb(i5,i2))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) + 
     & (4*s(i4,i5)*za(i4,i5)**2*zaa2(i2,i2,i6,i1,i3,i5)*
     &    zab2(i5,i4,i5,i6)*zb(i5,i2))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) - 
     & (4*s(i4,i5)*za(i2,i5)*za(i4,i5)*zaa2(i4,i1,i3,i2,i6,i5)*
     &    zab2(i5,i4,i5,i6)*zb(i5,i2))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) - 
     & (2*s(i4,i5)*za(i2,i5)*za(i4,i5)*zab2(i4,i4,i5,i5)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)**2*zb(i3,i2)) + 
     & (2*s(i4,i5)*za(i1,i5)**2*za(i4,i5)*zab2(i4,i4,i5,i3)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i3,i5)*
     &    zab2(i5,i4,i5,i3)**2*zab2(i5,i4,i5,i5)*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i2,i5)*za(i4,i5)*zab2(i4,i4,i5,i3)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i5,i4,i5,i5)*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i1,i5)*za(i4,i5)**2*zaa2(i1,i2,i3,i1,i6,i5)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i3,i5)**2*
     &    zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*zb(i3,i2)) + 
     & (2*s(i4,i5)*za(i1,i5)**2*za(i4,i5)*zaa2(i4,i1,i6,i2,i3,i5)*
     &    zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i3,i5)**2*
     &    zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i1,i5)*za(i4,i5)**2*zab2(i5,i4,i5,i2)*zb(i5,i2))/
     &  (za(i6,i5)**2*zaa2(i5,i1,i6,i2,i3,i5)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i2)) + 
     & (2*s(i4,i5)*za(i4,i5)**2*zab2(i6,i4,i5,i2)*zb(i5,i2))/
     &  (za(i6,i5)**2*zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*
     &    zab2(i6,i4,i5,i1)*zb(i3,i2)) - 
     & (2*s(i4,i5)**2*za(i2,i5)*za(i4,i5)**2*zb(i2,i1)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i4,i5,i1)**2*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i2)) + 
     & (2*s(i4,i5)**2*za(i4,i5)**2*zb(i2,i1)*zb(i5,i2))/
     &  (za(i6,i5)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zab2(i6,i4,i5,i1)*zb(i3,i2)) + 
     & (2*s(i4,i5)**2*za(i2,i5)**2*za(i4,i5)**2*zb(i5,i2)**2)/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)**2*zb(i3,i1)) - 
     & (2*s(i4,i5)**2*za(i1,i5)**2*za(i4,i5)**2*zb(i5,i2)**2)/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i3,i5)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)**2*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i1,i5)*za(i4,i5)*zab2(i4,i4,i5,i5)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)**2*
     &    zb(i3,i1)) - (2*s(i4,i5)*za(i1,i5)*za(i4,i5)*
     &    zab2(i4,i4,i5,i3)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i5,i4,i5,i3)**2*zab2(i5,i4,i5,i5)*
     &    zb(i3,i1)) + (2*s(i4,i5)*za(i1,i2)*za(i4,i5)**2*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i5)**2*zab2(i5,i4,i5,i3)*zab2(i5,i4,i5,i5)*
     &    zb(i3,i1)) - (2*s(i4,i5)*za(i4,i5)*za(i6,i5)*
     &    zaa2(i4,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i6)*zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) + 
     & (2*s(i4,i5)*za(i4,i5)*za(i6,i5)*zaa2(i4,i2,i6,i1,i3,i5)*
     &    zab2(i5,i4,i5,i6)*zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) - 
     & (2*s(i4,i5)*za(i4,i5)**2*zab2(i5,i4,i5,i6)*zab2(i6,i4,i5,i3)*
     &    zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i5,i4,i5,i5)*zb(i3,i1)) - 
     & (2*s(i4,i5)*za(i4,i5)*zab2(i4,i4,i5,i5)*zab2(i5,i4,i5,i2)*
     &    zb(i5,i6))/
     &  (za(i2,i6)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)**2*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i4,i5)*zab2(i4,i4,i5,i3)*zab2(i5,i4,i5,i2)*
     &    zb(i5,i6))/
     &  (za(i2,i6)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i5,i4,i5,i5)*zb(i3,i2)) - 
     & (2*s(i4,i5)**2*za(i4,i5)**2*zb(i2,i1)*zb(i5,i6))/
     &  (za(i2,i6)*zab2(i5,i4,i5,i1)**2*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)*zb(i3,i2)) + 
     & (4*s(i4,i5)**2*za(i2,i5)*za(i4,i5)**2*zb(i5,i2)*zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)**2*zb(i3,i1)) + 
     & (2*s(i4,i5)**2*za(i4,i5)**2*za(i6,i5)*zb(i5,i6)**2)/
     &  (za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)*zab2(i5,i4,i5,i3)*
     &    zab2(i5,i4,i5,i5)**2*zb(i3,i1)) - 
     & (2*za(i1,i3)*za(i2,i4)**2*zab2(i1,i4,i5,i2)*zb(i6,i1))/
     &  (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i2,i4,i5,i2)*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i4,i5)*zab2(i1,i4,i5,i5)*zab2(i4,i4,i5,i5)*
     &    zb(i6,i1))/
     &  (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)*zab2(i5,i4,i5,i5)**2*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i2,i4)*za(i4,i5)*zab2(i1,i4,i5,i5)*
     &    zb(i6,i1))/
     &  (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i5,i4,i5,i5)*zb(i3,i1))
     &   + (2*s(i4,i5)*za(i1,i2)*za(i2,i4)**2*zb(i6,i2))/
     &  (za(i1,i6)*za(i2,i5)**2*zab2(i2,i4,i5,i2)*zab2(i2,i4,i5,i3)*
     &    zb(i3,i1)) - (2*s(i4,i5)*za(i4,i6)**2*zb(i6,i2))/
     &  (za(i6,i5)**2*zab2(i6,i4,i5,i3)*zab2(i6,i4,i5,i6)*zb(i3,i1)) - 
     & (2*s(i4,i5)*za(i2,i6)*za(i4,i6)**2*zab2(i6,i4,i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)**2*zaa2(i6,i1,i3,i2,i6,i6)*zab2(i6,i4,i5,i3)*
     &    zab2(i6,i4,i5,i6)*zb(i3,i1)) + 
     & (2*za(i2,i3)*zab2(i1,i4,i5,i2)*zab2(i4,i4,i5,i1)**2*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i1,i4,i5,i1)*zab2(i5,i4,i5,i1)**2*
     &    zb(i3,i2)) - (2*za(i3,i6)*zab2(i1,i4,i5,i6)*
     &    zab2(i4,i4,i5,i1)**2*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i1,i4,i5,i1)*zab2(i5,i4,i5,i1)**2*
     &    zb(i3,i2)) + (2*s(i4,i5)*za(i1,i6)*za(i4,i6)**2*
     &    zab2(i6,i4,i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)**2*zaa2(i6,i1,i6,i2,i3,i6)*zab2(i6,i4,i5,i3)*
     &    zab2(i6,i4,i5,i6)*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i4,i6)**2*zab2(i6,i4,i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)**2*zab2(i6,i4,i5,i1)*zab2(i6,i4,i5,i3)*
     &    zab2(i6,i4,i5,i6)*zb(i3,i2)) - 
     & (2*s(i4,i5)*za(i2,i3)*za(i4,i5)*zab2(i4,i4,i5,i5)*zb(i5,i2)*
     &    zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i5)**2*
     &    zb(i3,i2)) - (2*s(i4,i5)*za(i2,i3)*za(i4,i5)*
     &    zab2(i4,i4,i5,i1)*zb(i5,i2)*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i5,i4,i5,i1)**2*zab2(i5,i4,i5,i5)*
     &    zb(i3,i2)) + (2*s(i4,i5)*za(i3,i6)*za(i4,i5)*
     &    zab2(i4,i4,i5,i5)*zb(i5,i6)*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i5,i4,i5,i1)*zab2(i5,i4,i5,i5)**2*
     &    zb(i3,i2)) + (2*s(i4,i5)*za(i3,i6)*za(i4,i5)*
     &    zab2(i4,i4,i5,i1)*zb(i5,i6)*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*zab2(i5,i4,i5,i1)**2*zab2(i5,i4,i5,i5)*
     &    zb(i3,i2)) - (2*za(i1,i3)*za(i2,i4)**2*zab2(i3,i4,i5,i2)*
     &    zb(i6,i3))/
     &  (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i2,i4,i5,i2)*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i4,i5)*zab2(i3,i4,i5,i5)*zab2(i4,i4,i5,i5)*
     &    zb(i6,i3))/
     &  (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)*zab2(i5,i4,i5,i5)**2*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i2,i4)*za(i4,i5)*zab2(i3,i4,i5,i5)*
     &    zb(i6,i3))/
     &  (t(i2,i4,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i5,i4,i5,i5)*zb(i3,i1))
     &   - (2*zab2(i1,i4,i5,i3)**2*zab2(i3,i4,i5,i2)*
     &    zab2(i4,i4,i5,i3)**2)/
     &  (za(i1,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i6,i4,i5,i3)*zbb2(i3,i2,i3,i1,i6,i3)) + 
     & ((sqrtd1623*za(i1,i5) + zaa2(i1,i1,i6,i2,i3,i5) - 
     &      zaa2(i1,i2,i3,i1,i6,i5))*
     &    (sqrtd1623*za(i1,i6) + zaa2(i1,i1,i6,i2,i3,i6) - 
     &      zaa2(i1,i2,i3,i1,i6,i6))*
     &    (sqrtd1623*za(i4,i5) + zaa2(i4,i1,i6,i2,i3,i5) - 
     &      zaa2(i4,i2,i3,i1,i6,i5))*
     &    (sqrtd1623*zab2(i4,i4,i5,i2) + 
     &      zab3(i4,i1,i6,i2,i3,i4,i5,i2) - 
     &      zab3(i4,i2,i3,i1,i6,i4,i5,i2))*
     &    (sqrtd1623*zb(i3,i2) + zbb2(i3,i1,i6,i2,i3,i2) - 
     &      zbb2(i3,i2,i3,i1,i6,i2)))/
     &  (16.*sqrtd1623**2*za(i1,i6)*zaa2(i5,i1,i6,i2,i3,i5)**2*
     &    zaa2(i6,i1,i6,i2,i3,i6)*zb(i3,i2)*zbb2(i3,i2,i3,i1,i6,i3)) - 
     & ((sqrtd1623*za(i1,i5) - zaa2(i1,i1,i6,i2,i3,i5) + 
     &      zaa2(i1,i2,i3,i1,i6,i5))*
     &    (sqrtd1623*za(i1,i6) - zaa2(i1,i1,i6,i2,i3,i6) + 
     &      zaa2(i1,i2,i3,i1,i6,i6))*
     &    (sqrtd1623*za(i4,i5) - zaa2(i4,i1,i6,i2,i3,i5) + 
     &      zaa2(i4,i2,i3,i1,i6,i5))*
     &    (sqrtd1623*zab2(i4,i4,i5,i2) - 
     &      zab3(i4,i1,i6,i2,i3,i4,i5,i2) + 
     &      zab3(i4,i2,i3,i1,i6,i4,i5,i2))*
     &    (sqrtd1623*zb(i3,i2) - zbb2(i3,i1,i6,i2,i3,i2) + 
     &      zbb2(i3,i2,i3,i1,i6,i2)))/
     &  (16.*sqrtd1623**2*za(i1,i6)*zaa2(i5,i1,i6,i2,i3,i5)**2*
     &    zaa2(i6,i1,i6,i2,i3,i6)*zb(i3,i2)*zbb2(i3,i2,i3,i1,i6,i3)) + 
     & (2*zab2(i2,i4,i5,i3)**2*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*
     &    zb(i3,i2))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)**2*
     &    zab2(i6,i4,i5,i3)*zb(i3,i1)*zbb2(i3,i2,i6,i1,i3,i3)) - 
     & (4*zab2(i2,i4,i5,i3)*zab2(i3,i4,i5,i2)*zab2(i4,i4,i5,i3)**2*
     &    zb(i6,i3))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i6,i1,i3,i3)) - 
     & (2*zab2(i3,i4,i5,i6)*zab2(i4,i4,i5,i3)**2*zab2(i6,i4,i5,i3)*
     &    zb(i6,i3))/
     &  (za(i2,i6)*zab2(i3,i4,i5,i3)*zab2(i5,i4,i5,i3)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i6,i1,i3,i3)) - 
     & ((sqrtd1326*za(i2,i5) + zaa2(i2,i1,i3,i2,i6,i5) - 
     &      zaa2(i2,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*za(i2,i6) + zaa2(i2,i1,i3,i2,i6,i6) - 
     &      zaa2(i2,i2,i6,i1,i3,i6))*
     &    (sqrtd1326*za(i4,i5) + zaa2(i4,i1,i3,i2,i6,i5) - 
     &      zaa2(i4,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*zab2(i4,i4,i5,i2) + 
     &      zab3(i4,i1,i3,i2,i6,i4,i5,i2) - 
     &      zab3(i4,i2,i6,i1,i3,i4,i5,i2))*
     &    (sqrtd1326*zb(i3,i2) + zbb2(i3,i1,i3,i2,i6,i2) - 
     &      zbb2(i3,i2,i6,i1,i3,i2)))/
     &  (16.*sqrtd1326**2*za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*
     &    zaa2(i6,i1,i3,i2,i6,i6)*zb(i3,i1)*zbb2(i3,i2,i6,i1,i3,i3)) + 
     & ((sqrtd1326*za(i2,i5) - zaa2(i2,i1,i3,i2,i6,i5) + 
     &      zaa2(i2,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*za(i2,i6) - zaa2(i2,i1,i3,i2,i6,i6) + 
     &      zaa2(i2,i2,i6,i1,i3,i6))*
     &    (sqrtd1326*za(i4,i5) - zaa2(i4,i1,i3,i2,i6,i5) + 
     &      zaa2(i4,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*zab2(i4,i4,i5,i2) - 
     &      zab3(i4,i1,i3,i2,i6,i4,i5,i2) + 
     &      zab3(i4,i2,i6,i1,i3,i4,i5,i2))*
     &    (sqrtd1326*zb(i3,i2) - zbb2(i3,i1,i3,i2,i6,i2) + 
     &      zbb2(i3,i2,i6,i1,i3,i2)))/
     &  (16.*sqrtd1326**2*za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*
     &    zaa2(i6,i1,i3,i2,i6,i6)*zb(i3,i1)*zbb2(i3,i2,i6,i1,i3,i3)) + 
     & ((sqrtd1326*za(i2,i5) - zaa2(i2,i1,i3,i2,i6,i5) + 
     &      zaa2(i2,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*za(i4,i5) - zaa2(i4,i1,i3,i2,i6,i5) + 
     &      zaa2(i4,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*zab2(i4,i4,i5,i2) - 
     &      zab3(i4,i1,i3,i2,i6,i4,i5,i2) + 
     &      zab3(i4,i2,i6,i1,i3,i4,i5,i2))*
     &    (sqrtd1326*zb(i6,i3) + zbb2(i6,i1,i3,i2,i6,i3) - 
     &      zbb2(i6,i2,i6,i1,i3,i3)))/
     &  (4.*sqrtd1326**2*za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i6,i1,i3,i3)) + 
     & ((sqrtd1326*za(i4,i5) - zaa2(i4,i1,i3,i2,i6,i5) + 
     &       zaa2(i4,i2,i6,i1,i3,i5))**2*
     &    (sqrtd1326*zab2(i6,i4,i5,i6) - 
     &      zab3(i6,i1,i3,i2,i6,i4,i5,i6) + 
     &      zab3(i6,i2,i6,i1,i3,i4,i5,i6))*
     &    (sqrtd1326*zb(i6,i3) + zbb2(i6,i1,i3,i2,i6,i3) - 
     &      zbb2(i6,i2,i6,i1,i3,i3)))/
     &  (8.*sqrtd1326**2*za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i6,i1,i3,i3)) + 
     & ((sqrtd1326*za(i2,i5) + zaa2(i2,i1,i3,i2,i6,i5) - 
     &      zaa2(i2,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*za(i4,i5) + zaa2(i4,i1,i3,i2,i6,i5) - 
     &      zaa2(i4,i2,i6,i1,i3,i5))*
     &    (sqrtd1326*zab2(i4,i4,i5,i2) + 
     &      zab3(i4,i1,i3,i2,i6,i4,i5,i2) - 
     &      zab3(i4,i2,i6,i1,i3,i4,i5,i2))*
     &    (sqrtd1326*zb(i6,i3) - zbb2(i6,i1,i3,i2,i6,i3) + 
     &      zbb2(i6,i2,i6,i1,i3,i3)))/
     &  (4.*sqrtd1326**2*za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i6,i1,i3,i3)) + 
     & ((sqrtd1326*za(i4,i5) + zaa2(i4,i1,i3,i2,i6,i5) - 
     &       zaa2(i4,i2,i6,i1,i3,i5))**2*
     &    (sqrtd1326*zab2(i6,i4,i5,i6) + 
     &      zab3(i6,i1,i3,i2,i6,i4,i5,i6) - 
     &      zab3(i6,i2,i6,i1,i3,i4,i5,i6))*
     &    (sqrtd1326*zb(i6,i3) - zbb2(i6,i1,i3,i2,i6,i3) + 
     &      zbb2(i6,i2,i6,i1,i3,i3)))/
     &  (8.*sqrtd1326**2*za(i2,i6)*zaa2(i5,i1,i3,i2,i6,i5)**2*zb(i3,i1)*
     &    zbb2(i3,i2,i6,i1,i3,i3))
     &)


      aaaa_NMHV_b(15) = (
     &  (2*zab2(i1,i3,i6,i4)*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i5))/
     &  (za(i1,i5)*zab2(i2,i3,i6,i4)*zab2(i4,i3,i6,i4)*
     &    zab2(i6,i3,i6,i4)**2*zb(i4,i1)) + 
     & (2*zab2(i1,i3,i6,i4)*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2))/
     &  (za(i1,i5)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i4)**2*zb(i4,i1)) + 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i5)*zab2(i6,i3,i6,i2))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i6,i3,i6,i1)*
     &    zab2(i6,i3,i6,i4)**2*zb(i4,i2)) + 
     & (2*zab2(i2,i3,i6,i4)*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*
     &    zab2(i6,i3,i6,i2))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)**2*zb(i4,i2)) - 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*zab2(i5,i3,i6,i2))/
     &  (za(i5,i6)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i1)*
     &    zab2(i5,i3,i6,i4)*zab2(i6,i3,i6,i4)*zb(i4,i2)) + 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*zab2(i6,i3,i6,i2))/
     &  (za(i5,i6)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)*zb(i4,i2)) + 
     & (2*zab2(i1,i3,i6,i5)*zab2(i3,i3,i6,i1)**2*zb(i2,i1))/
     &  (za(i2,i5)*zab2(i1,i3,i6,i1)*zab2(i6,i3,i6,i1)**2*zb(i4,i1)*
     &    zb(i4,i2)) + (2*zab2(i1,i3,i6,i2)*zab2(i2,i3,i6,i1)*
     &    zab2(i3,i3,i6,i1)**2*zb(i2,i1))/
     &  (za(i2,i5)*zab2(i1,i3,i6,i1)*zab2(i5,i3,i6,i1)*
     &    zab2(i6,i3,i6,i1)**2*zb(i4,i1)*zb(i4,i2)) - 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*zb(i2,i1))/
     &  (za(i5,i6)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*zb(i2,i1))/
     &  (za(i5,i6)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i1)*
     &    zab2(i6,i3,i6,i4)*zb(i4,i1)*zb(i4,i2)) - 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i5)*zb(i2,i1))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i6,i3,i6,i1)*
     &    zab2(i6,i3,i6,i4)*zb(i4,i1)*zb(i4,i2)) - 
     & (2*zab2(i2,i3,i6,i4)*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*
     &    zb(i2,i1))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)*zb(i4,i1)*zb(i4,i2)) - 
     & (2*za(i1,i4)*za(i2,i3)**2*zab2(i1,i3,i6,i2)*zb(i5,i1))/
     &  (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i2,i3,i6,i2)*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i3,i6)*zab2(i1,i3,i6,i6)*zab2(i3,i3,i6,i6)*
     &    zb(i5,i1))/
     &  (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)*zab2(i6,i3,i6,i6)**2*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i2,i3)*za(i3,i6)*zab2(i1,i3,i6,i6)*
     &    zb(i5,i1))/
     &  (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i6,i3,i6,i6)*zb(i4,i1))
     &   + (2*s(i3,i6)*za(i1,i2)*za(i2,i3)**2*zb(i5,i2))/
     &  (za(i1,i5)*za(i2,i6)**2*zab2(i2,i3,i6,i2)*zab2(i2,i3,i6,i4)*
     &    zb(i4,i1)) - (2*s(i3,i6)*za(i3,i5)**2*zb(i5,i2))/
     &  (za(i5,i6)**2*zab2(i5,i3,i6,i4)*zab2(i5,i3,i6,i5)*zb(i4,i1)) - 
     & (2*s(i3,i6)*za(i2,i5)*za(i3,i5)**2*zab2(i5,i3,i6,i2)*zb(i5,i2))/
     &  (za(i5,i6)**2*zaa2(i5,i1,i4,i2,i5,i5)*zab2(i5,i3,i6,i4)*
     &    zab2(i5,i3,i6,i5)*zb(i4,i1)) + 
     & (2*s(i3,i6)*za(i1,i5)*za(i3,i5)**2*zab2(i5,i3,i6,i2)*zb(i5,i2))/
     &  (za(i5,i6)**2*zaa2(i5,i1,i5,i2,i4,i5)*zab2(i5,i3,i6,i4)*
     &    zab2(i5,i3,i6,i5)*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i3,i5)**2*zab2(i5,i3,i6,i2)*zb(i5,i2))/
     &  (za(i5,i6)**2*zab2(i5,i3,i6,i1)*zab2(i5,i3,i6,i4)*
     &    zab2(i5,i3,i6,i5)*zb(i4,i2)) + 
     & (2*za(i2,i4)*zab2(i1,i3,i6,i2)*zab2(i3,i3,i6,i1)**2*zb(i5,i2))/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i1,i3,i6,i1)*zab2(i6,i3,i6,i1)**2*
     &    zb(i4,i2)) - (2*za(i4,i5)*zab2(i1,i3,i6,i5)*
     &    zab2(i3,i3,i6,i1)**2*zb(i5,i2))/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i1,i3,i6,i1)*zab2(i6,i3,i6,i1)**2*
     &    zb(i4,i2)) - (2*za(i1,i4)*za(i2,i3)**2*zab2(i4,i3,i6,i2)*
     &    zb(i5,i4))/
     &  (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i2,i3,i6,i2)*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i3,i6)*zab2(i3,i3,i6,i6)*zab2(i4,i3,i6,i6)*
     &    zb(i5,i4))/
     &  (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)*zab2(i6,i3,i6,i6)**2*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i2,i3)*za(i3,i6)*zab2(i4,i3,i6,i6)*
     &    zb(i5,i4))/
     &  (t(i2,i3,i6)*za(i1,i5)*za(i2,i6)**2*zab2(i6,i3,i6,i6)*zb(i4,i1))
     &   - (2*s(i3,i6)*za(i1,i6)*za(i3,i6)*zab2(i3,i3,i6,i6)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)**2*
     &    zb(i4,i1)) - (2*s(i3,i6)*za(i1,i6)*za(i3,i6)*
     &    zab2(i3,i3,i6,i4)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zab2(i6,i3,i6,i4)**2*zab2(i6,i3,i6,i6)*
     &    zb(i4,i1)) - (2*s(i3,i6)*za(i2,i6)**2*za(i3,i6)*
     &    zab2(i3,i3,i6,i4)*zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i4,i2,i5,i6)*
     &    zab2(i6,i3,i6,i4)**2*zab2(i6,i3,i6,i6)*zb(i4,i1)) + 
     & (2*s(i3,i6)*za(i3,i6)**2*zb(i6,i2))/
     &  (za(i5,i6)**2*zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*zb(i4,i1)) + 
     & (2*s(i3,i6)*za(i2,i6)*za(i3,i6)**2*zaa2(i2,i2,i5,i1,i4,i6)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i4,i2,i5,i6)**2*
     &    zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*zb(i4,i1)) - 
     & (2*s(i3,i6)*za(i2,i6)**2*za(i3,i6)*zaa2(i3,i1,i4,i2,i5,i6)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i4,i2,i5,i6)**2*
     &    zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*zb(i4,i1)) + 
     & (2*s(i3,i6)*za(i2,i6)*za(i3,i6)**2*zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i5,i6)**2*zaa2(i6,i1,i4,i2,i5,i6)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) - 
     & (4*s(i3,i6)*za(i2,i6)*za(i3,i6)*zab2(i3,i3,i6,i4)*
     &    zab2(i6,i3,i6,i5)*zb(i6,i2))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)*zab2(i6,i3,i6,i4)**2*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) + 
     & (4*s(i3,i6)*za(i3,i6)**2*zaa2(i2,i2,i5,i1,i4,i6)*
     &    zab2(i6,i3,i6,i5)*zb(i6,i2))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) - 
     & (4*s(i3,i6)*za(i2,i6)*za(i3,i6)*zaa2(i3,i1,i4,i2,i5,i6)*
     &    zab2(i6,i3,i6,i5)*zb(i6,i2))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) - 
     & (2*s(i3,i6)*za(i2,i6)*za(i3,i6)*zab2(i3,i3,i6,i6)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)**2*zb(i4,i2)) + 
     & (2*s(i3,i6)*za(i1,i6)**2*za(i3,i6)*zab2(i3,i3,i6,i4)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i4,i6)*
     &    zab2(i6,i3,i6,i4)**2*zab2(i6,i3,i6,i6)*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i2,i6)*za(i3,i6)*zab2(i3,i3,i6,i4)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)**2*
     &    zab2(i6,i3,i6,i6)*zb(i4,i2)) + 
     & (2*s(i3,i6)*za(i3,i6)**2*zab2(i5,i3,i6,i2)*zb(i6,i2))/
     &  (za(i5,i6)**2*zab2(i5,i3,i6,i1)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i1,i6)*za(i3,i6)**2*zaa2(i1,i2,i4,i1,i5,i6)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i4,i6)**2*
     &    zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*zb(i4,i2)) + 
     & (2*s(i3,i6)*za(i1,i6)**2*za(i3,i6)*zaa2(i3,i1,i5,i2,i4,i6)*
     &    zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i4,i6)**2*
     &    zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i1,i6)*za(i3,i6)**2*zab2(i6,i3,i6,i2)*zb(i6,i2))/
     &  (za(i5,i6)**2*zaa2(i6,i1,i5,i2,i4,i6)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i2)) - 
     & (2*s(i3,i6)**2*za(i2,i6)*za(i3,i6)**2*zb(i2,i1)*zb(i6,i2))/
     &  (za(i2,i5)*za(i5,i6)*zab2(i6,i3,i6,i1)**2*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i2)) + 
     & (2*s(i3,i6)**2*za(i3,i6)**2*zb(i2,i1)*zb(i6,i2))/
     &  (za(i5,i6)*zab2(i5,i3,i6,i1)*zab2(i6,i3,i6,i1)*
     &    zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i2,i4)*za(i3,i6)*zab2(i3,i3,i6,i6)*zb(i5,i2)*
     &    zb(i6,i2))/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i6)**2*
     &    zb(i4,i2)) - (2*s(i3,i6)*za(i2,i4)*za(i3,i6)*
     &    zab2(i3,i3,i6,i1)*zb(i5,i2)*zb(i6,i2))/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i6,i3,i6,i1)**2*zab2(i6,i3,i6,i6)*
     &    zb(i4,i2)) + (2*s(i3,i6)**2*za(i2,i6)**2*za(i3,i6)**2*
     &    zb(i6,i2)**2)/
     &  (za(i2,i5)*za(i5,i6)*zaa2(i6,i1,i4,i2,i5,i6)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)**2*zb(i4,i1)) - 
     & (2*s(i3,i6)**2*za(i1,i6)**2*za(i3,i6)**2*zb(i6,i2)**2)/
     &  (za(i1,i5)*za(i5,i6)*zaa2(i6,i1,i5,i2,i4,i6)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)**2*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i1,i6)*za(i3,i6)*zab2(i3,i3,i6,i6)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)**2*
     &    zb(i4,i1)) - (2*s(i3,i6)*za(i1,i6)*za(i3,i6)*
     &    zab2(i3,i3,i6,i4)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zab2(i6,i3,i6,i4)**2*zab2(i6,i3,i6,i6)*
     &    zb(i4,i1)) + (2*s(i3,i6)*za(i1,i2)*za(i3,i6)**2*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)**2*zab2(i6,i3,i6,i4)*zab2(i6,i3,i6,i6)*
     &    zb(i4,i1)) - (2*s(i3,i6)*za(i3,i6)**2*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i5)*zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)*zab2(i6,i3,i6,i4)**2*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) - 
     & (2*s(i3,i6)*za(i3,i6)*za(i5,i6)*zaa2(i3,i1,i4,i2,i5,i6)*
     &    zab2(i6,i3,i6,i5)*zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) + 
     & (2*s(i3,i6)*za(i3,i6)*za(i5,i6)*zaa2(i3,i2,i5,i1,i4,i6)*
     &    zab2(i6,i3,i6,i5)*zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i1)) - 
     & (2*s(i3,i6)*za(i3,i6)*zab2(i3,i3,i6,i6)*zab2(i6,i3,i6,i2)*
     &    zb(i6,i5))/
     &  (za(i2,i5)*zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)**2*zb(i4,i2)) - 
     & (2*s(i3,i6)*za(i3,i6)*zab2(i3,i3,i6,i4)*zab2(i6,i3,i6,i2)*
     &    zb(i6,i5))/
     &  (za(i2,i5)*zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i4)**2*
     &    zab2(i6,i3,i6,i6)*zb(i4,i2)) - 
     & (2*s(i3,i6)**2*za(i3,i6)**2*zb(i2,i1)*zb(i6,i5))/
     &  (za(i2,i5)*zab2(i6,i3,i6,i1)**2*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)*zb(i4,i2)) + 
     & (2*s(i3,i6)*za(i3,i6)*za(i4,i5)*zab2(i3,i3,i6,i6)*zb(i5,i2)*
     &    zb(i6,i5))/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i6,i3,i6,i1)*zab2(i6,i3,i6,i6)**2*
     &    zb(i4,i2)) + (2*s(i3,i6)*za(i3,i6)*za(i4,i5)*
     &    zab2(i3,i3,i6,i1)*zb(i5,i2)*zb(i6,i5))/
     &  (t(i2,i4,i5)*za(i2,i5)*zab2(i6,i3,i6,i1)**2*zab2(i6,i3,i6,i6)*
     &    zb(i4,i2)) + (4*s(i3,i6)**2*za(i2,i6)*za(i3,i6)**2*zb(i6,i2)*
     &    zb(i6,i5))/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)**2*zb(i4,i1)) + 
     & (2*s(i3,i6)**2*za(i3,i6)**2*za(i5,i6)*zb(i6,i5)**2)/
     &  (za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)*zab2(i6,i3,i6,i4)*
     &    zab2(i6,i3,i6,i6)**2*zb(i4,i1)) - 
     & (2*zab2(i1,i3,i6,i4)**2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2))/
     &  (za(i1,i5)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i4)**2*zbb2(i4,i2,i4,i1,i5,i4)) + 
     & ((sqrtd1524*za(i1,i5) + zaa2(i1,i1,i5,i2,i4,i5) - 
     &      zaa2(i1,i2,i4,i1,i5,i5))*
     &    (sqrtd1524*za(i1,i6) + zaa2(i1,i1,i5,i2,i4,i6) - 
     &      zaa2(i1,i2,i4,i1,i5,i6))*
     &    (sqrtd1524*za(i3,i6) + zaa2(i3,i1,i5,i2,i4,i6) - 
     &      zaa2(i3,i2,i4,i1,i5,i6))*
     &    (sqrtd1524*zab2(i3,i3,i6,i2) + 
     &      zab3(i3,i1,i5,i2,i4,i3,i6,i2) - 
     &      zab3(i3,i2,i4,i1,i5,i3,i6,i2))*
     &    (sqrtd1524*zb(i4,i2) + zbb2(i4,i1,i5,i2,i4,i2) - 
     &      zbb2(i4,i2,i4,i1,i5,i2)))/
     &  (16.*sqrtd1524**2*za(i1,i5)*zaa2(i5,i1,i5,i2,i4,i5)*
     &    zaa2(i6,i1,i5,i2,i4,i6)**2*zb(i4,i2)*zbb2(i4,i2,i4,i1,i5,i4))
     &  - ((sqrtd1524*za(i1,i5) - zaa2(i1,i1,i5,i2,i4,i5) + 
     &      zaa2(i1,i2,i4,i1,i5,i5))*
     &    (sqrtd1524*za(i1,i6) - zaa2(i1,i1,i5,i2,i4,i6) + 
     &      zaa2(i1,i2,i4,i1,i5,i6))*
     &    (sqrtd1524*za(i3,i6) - zaa2(i3,i1,i5,i2,i4,i6) + 
     &      zaa2(i3,i2,i4,i1,i5,i6))*
     &    (sqrtd1524*zab2(i3,i3,i6,i2) - 
     &      zab3(i3,i1,i5,i2,i4,i3,i6,i2) + 
     &      zab3(i3,i2,i4,i1,i5,i3,i6,i2))*
     &    (sqrtd1524*zb(i4,i2) - zbb2(i4,i1,i5,i2,i4,i2) + 
     &      zbb2(i4,i2,i4,i1,i5,i2)))/
     &  (16.*sqrtd1524**2*za(i1,i5)*zaa2(i5,i1,i5,i2,i4,i5)*
     &    zaa2(i6,i1,i5,i2,i4,i6)**2*zb(i4,i2)*zbb2(i4,i2,i4,i1,i5,i4))
     &  + (2*zab2(i2,i3,i6,i4)**2*zab2(i3,i3,i6,i4)**2*
     &    zab2(i4,i3,i6,i2)*zb(i4,i2))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i5,i3,i6,i4)*
     &    zab2(i6,i3,i6,i4)**2*zb(i4,i1)*zbb2(i4,i2,i5,i1,i4,i4)) - 
     & (4*zab2(i2,i3,i6,i4)*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i2)*
     &    zb(i5,i4))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i6,i3,i6,i4)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i5,i1,i4,i4)) - 
     & (2*zab2(i3,i3,i6,i4)**2*zab2(i4,i3,i6,i5)*zab2(i5,i3,i6,i4)*
     &    zb(i5,i4))/
     &  (za(i2,i5)*zab2(i4,i3,i6,i4)*zab2(i6,i3,i6,i4)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i5,i1,i4,i4)) - 
     & ((sqrtd1425*za(i2,i5) + zaa2(i2,i1,i4,i2,i5,i5) - 
     &      zaa2(i2,i2,i5,i1,i4,i5))*
     &    (sqrtd1425*za(i2,i6) + zaa2(i2,i1,i4,i2,i5,i6) - 
     &      zaa2(i2,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*za(i3,i6) + zaa2(i3,i1,i4,i2,i5,i6) - 
     &      zaa2(i3,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*zab2(i3,i3,i6,i2) + 
     &      zab3(i3,i1,i4,i2,i5,i3,i6,i2) - 
     &      zab3(i3,i2,i5,i1,i4,i3,i6,i2))*
     &    (sqrtd1425*zb(i4,i2) + zbb2(i4,i1,i4,i2,i5,i2) - 
     &      zbb2(i4,i2,i5,i1,i4,i2)))/
     &  (16.*sqrtd1425**2*za(i2,i5)*zaa2(i5,i1,i4,i2,i5,i5)*
     &    zaa2(i6,i1,i4,i2,i5,i6)**2*zb(i4,i1)*zbb2(i4,i2,i5,i1,i4,i4))
     &  + ((sqrtd1425*za(i2,i5) - zaa2(i2,i1,i4,i2,i5,i5) + 
     &      zaa2(i2,i2,i5,i1,i4,i5))*
     &    (sqrtd1425*za(i2,i6) - zaa2(i2,i1,i4,i2,i5,i6) + 
     &      zaa2(i2,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*za(i3,i6) - zaa2(i3,i1,i4,i2,i5,i6) + 
     &      zaa2(i3,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*zab2(i3,i3,i6,i2) - 
     &      zab3(i3,i1,i4,i2,i5,i3,i6,i2) + 
     &      zab3(i3,i2,i5,i1,i4,i3,i6,i2))*
     &    (sqrtd1425*zb(i4,i2) - zbb2(i4,i1,i4,i2,i5,i2) + 
     &      zbb2(i4,i2,i5,i1,i4,i2)))/
     &  (16.*sqrtd1425**2*za(i2,i5)*zaa2(i5,i1,i4,i2,i5,i5)*
     &    zaa2(i6,i1,i4,i2,i5,i6)**2*zb(i4,i1)*zbb2(i4,i2,i5,i1,i4,i4))
     &  + ((sqrtd1425*za(i2,i6) - zaa2(i2,i1,i4,i2,i5,i6) + 
     &      zaa2(i2,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*za(i3,i6) - zaa2(i3,i1,i4,i2,i5,i6) + 
     &      zaa2(i3,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*zab2(i3,i3,i6,i2) - 
     &      zab3(i3,i1,i4,i2,i5,i3,i6,i2) + 
     &      zab3(i3,i2,i5,i1,i4,i3,i6,i2))*
     &    (sqrtd1425*zb(i5,i4) + zbb2(i5,i1,i4,i2,i5,i4) - 
     &      zbb2(i5,i2,i5,i1,i4,i4)))/
     &  (4.*sqrtd1425**2*za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i5,i1,i4,i4)) + 
     & ((sqrtd1425*za(i3,i6) - zaa2(i3,i1,i4,i2,i5,i6) + 
     &       zaa2(i3,i2,i5,i1,i4,i6))**2*
     &    (sqrtd1425*zab2(i5,i3,i6,i5) - 
     &      zab3(i5,i1,i4,i2,i5,i3,i6,i5) + 
     &      zab3(i5,i2,i5,i1,i4,i3,i6,i5))*
     &    (sqrtd1425*zb(i5,i4) + zbb2(i5,i1,i4,i2,i5,i4) - 
     &      zbb2(i5,i2,i5,i1,i4,i4)))/
     &  (8.*sqrtd1425**2*za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i5,i1,i4,i4)) + 
     & ((sqrtd1425*za(i2,i6) + zaa2(i2,i1,i4,i2,i5,i6) - 
     &      zaa2(i2,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*za(i3,i6) + zaa2(i3,i1,i4,i2,i5,i6) - 
     &      zaa2(i3,i2,i5,i1,i4,i6))*
     &    (sqrtd1425*zab2(i3,i3,i6,i2) + 
     &      zab3(i3,i1,i4,i2,i5,i3,i6,i2) - 
     &      zab3(i3,i2,i5,i1,i4,i3,i6,i2))*
     &    (sqrtd1425*zb(i5,i4) - zbb2(i5,i1,i4,i2,i5,i4) + 
     &      zbb2(i5,i2,i5,i1,i4,i4)))/
     &  (4.*sqrtd1425**2*za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i5,i1,i4,i4)) + 
     & ((sqrtd1425*za(i3,i6) + zaa2(i3,i1,i4,i2,i5,i6) - 
     &       zaa2(i3,i2,i5,i1,i4,i6))**2*
     &    (sqrtd1425*zab2(i5,i3,i6,i5) + 
     &      zab3(i5,i1,i4,i2,i5,i3,i6,i5) - 
     &      zab3(i5,i2,i5,i1,i4,i3,i6,i5))*
     &    (sqrtd1425*zb(i5,i4) - zbb2(i5,i1,i4,i2,i5,i4) + 
     &      zbb2(i5,i2,i5,i1,i4,i4)))/
     &  (8.*sqrtd1425**2*za(i2,i5)*zaa2(i6,i1,i4,i2,i5,i6)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i5,i1,i4,i4))
     &)


      aaaa_NMHV_b(16) = (
     &  (2*zab2(i1,i3,i5,i4)*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i6))/
     &  (za(i1,i6)*zab2(i2,i3,i5,i4)*zab2(i4,i3,i5,i4)*
     &    zab2(i5,i3,i5,i4)**2*zb(i4,i1)) + 
     & (2*zab2(i1,i3,i5,i4)*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2))/
     &  (za(i1,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i6,i3,i5,i4)*zb(i4,i1)) + 
     & (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i6)*zab2(i5,i3,i5,i2))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i1)*
     &    zab2(i5,i3,i5,i4)**2*zb(i4,i2)) + 
     & (2*zab2(i2,i3,i5,i4)*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*
     &    zab2(i5,i3,i5,i2))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i1)*
     &    zab2(i5,i3,i5,i4)**2*zab2(i6,i3,i5,i4)*zb(i4,i2)) + 
     & (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*zab2(i5,i3,i5,i2))/
     &  (za(i6,i5)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i1)*
     &    zab2(i5,i3,i5,i4)*zab2(i6,i3,i5,i4)*zb(i4,i2)) - 
     & (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*zab2(i6,i3,i5,i2))/
     &  (za(i6,i5)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)*
     &    zab2(i6,i3,i5,i1)*zab2(i6,i3,i5,i4)*zb(i4,i2)) + 
     & (2*zab2(i1,i3,i5,i6)*zab2(i3,i3,i5,i1)**2*zb(i2,i1))/
     &  (za(i2,i6)*zab2(i1,i3,i5,i1)*zab2(i5,i3,i5,i1)**2*zb(i4,i1)*
     &    zb(i4,i2)) - (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i6)*
     &    zb(i2,i1))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i1)*
     &    zab2(i5,i3,i5,i4)*zb(i4,i1)*zb(i4,i2)) + 
     & (2*zab2(i1,i3,i5,i2)*zab2(i2,i3,i5,i1)*zab2(i3,i3,i5,i1)**2*
     &    zb(i2,i1))/
     &  (za(i2,i6)*zab2(i1,i3,i5,i1)*zab2(i5,i3,i5,i1)**2*
     &    zab2(i6,i3,i5,i1)*zb(i4,i1)*zb(i4,i2)) + 
     & (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*zb(i2,i1))/
     &  (za(i6,i5)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)*
     &    zab2(i6,i3,i5,i1)*zb(i4,i1)*zb(i4,i2)) - 
     & (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*zb(i2,i1))/
     &  (za(i6,i5)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i1)*
     &    zab2(i6,i3,i5,i4)*zb(i4,i1)*zb(i4,i2)) - 
     & (2*zab2(i2,i3,i5,i4)*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*
     &    zb(i2,i1))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i1)*
     &    zab2(i5,i3,i5,i4)*zab2(i6,i3,i5,i4)*zb(i4,i1)*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i1,i5)*za(i3,i5)*zab2(i3,i3,i5,i5)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)**2*
     &    zb(i4,i1)) - (2*s(i3,i5)*za(i1,i5)*za(i3,i5)*
     &    zab2(i3,i3,i5,i4)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zab2(i5,i3,i5,i4)**2*zab2(i5,i3,i5,i5)*
     &    zb(i4,i1)) - (2*s(i3,i5)*za(i2,i5)**2*za(i3,i5)*
     &    zab2(i3,i3,i5,i4)*zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i4,i2,i6,i5)*
     &    zab2(i5,i3,i5,i4)**2*zab2(i5,i3,i5,i5)*zb(i4,i1)) + 
     & (2*s(i3,i5)*za(i3,i5)**2*zb(i5,i2))/
     &  (za(i6,i5)**2*zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*zb(i4,i1)) + 
     & (2*s(i3,i5)*za(i2,i5)*za(i3,i5)**2*zaa2(i2,i2,i6,i1,i4,i5)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i4,i2,i6,i5)**2*
     &    zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*zb(i4,i1)) - 
     & (2*s(i3,i5)*za(i2,i5)**2*za(i3,i5)*zaa2(i3,i1,i4,i2,i6,i5)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i4,i2,i6,i5)**2*
     &    zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*zb(i4,i1)) + 
     & (2*s(i3,i5)*za(i2,i5)*za(i3,i5)**2*zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i6,i5)**2*zaa2(i5,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) - 
     & (4*s(i3,i5)*za(i2,i5)*za(i3,i5)*zab2(i3,i3,i5,i4)*
     &    zab2(i5,i3,i5,i6)*zb(i5,i2))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) + 
     & (4*s(i3,i5)*za(i3,i5)**2*zaa2(i2,i2,i6,i1,i4,i5)*
     &    zab2(i5,i3,i5,i6)*zb(i5,i2))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) - 
     & (4*s(i3,i5)*za(i2,i5)*za(i3,i5)*zaa2(i3,i1,i4,i2,i6,i5)*
     &    zab2(i5,i3,i5,i6)*zb(i5,i2))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) - 
     & (2*s(i3,i5)*za(i2,i5)*za(i3,i5)*zab2(i3,i3,i5,i5)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)**2*zb(i4,i2)) + 
     & (2*s(i3,i5)*za(i1,i5)**2*za(i3,i5)*zab2(i3,i3,i5,i4)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i4,i5)*
     &    zab2(i5,i3,i5,i4)**2*zab2(i5,i3,i5,i5)*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i2,i5)*za(i3,i5)*zab2(i3,i3,i5,i4)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i5,i3,i5,i5)*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i1,i5)*za(i3,i5)**2*zaa2(i1,i2,i4,i1,i6,i5)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i4,i5)**2*
     &    zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*zb(i4,i2)) + 
     & (2*s(i3,i5)*za(i1,i5)**2*za(i3,i5)*zaa2(i3,i1,i6,i2,i4,i5)*
     &    zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i4,i5)**2*
     &    zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i1,i5)*za(i3,i5)**2*zab2(i5,i3,i5,i2)*zb(i5,i2))/
     &  (za(i6,i5)**2*zaa2(i5,i1,i6,i2,i4,i5)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i2)) + 
     & (2*s(i3,i5)*za(i3,i5)**2*zab2(i6,i3,i5,i2)*zb(i5,i2))/
     &  (za(i6,i5)**2*zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*
     &    zab2(i6,i3,i5,i1)*zb(i4,i2)) - 
     & (2*s(i3,i5)**2*za(i2,i5)*za(i3,i5)**2*zb(i2,i1)*zb(i5,i2))/
     &  (za(i2,i6)*za(i6,i5)*zab2(i5,i3,i5,i1)**2*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i2)) + 
     & (2*s(i3,i5)**2*za(i3,i5)**2*zb(i2,i1)*zb(i5,i2))/
     &  (za(i6,i5)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zab2(i6,i3,i5,i1)*zb(i4,i2)) + 
     & (2*s(i3,i5)**2*za(i2,i5)**2*za(i3,i5)**2*zb(i5,i2)**2)/
     &  (za(i2,i6)*za(i6,i5)*zaa2(i5,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)**2*zb(i4,i1)) - 
     & (2*s(i3,i5)**2*za(i1,i5)**2*za(i3,i5)**2*zb(i5,i2)**2)/
     &  (za(i1,i6)*za(i6,i5)*zaa2(i5,i1,i6,i2,i4,i5)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)**2*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i1,i5)*za(i3,i5)*zab2(i3,i3,i5,i5)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)**2*
     &    zb(i4,i1)) - (2*s(i3,i5)*za(i1,i5)*za(i3,i5)*
     &    zab2(i3,i3,i5,i4)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i5)*zab2(i5,i3,i5,i4)**2*zab2(i5,i3,i5,i5)*
     &    zb(i4,i1)) + (2*s(i3,i5)*za(i1,i2)*za(i3,i5)**2*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i5)**2*zab2(i5,i3,i5,i4)*zab2(i5,i3,i5,i5)*
     &    zb(i4,i1)) - (2*s(i3,i5)*za(i3,i5)*za(i6,i5)*
     &    zaa2(i3,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i6)*zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) + 
     & (2*s(i3,i5)*za(i3,i5)*za(i6,i5)*zaa2(i3,i2,i6,i1,i4,i5)*
     &    zab2(i5,i3,i5,i6)*zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) - 
     & (2*s(i3,i5)*za(i3,i5)**2*zab2(i5,i3,i5,i6)*zab2(i6,i3,i5,i4)*
     &    zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i5,i3,i5,i5)*zb(i4,i1)) - 
     & (2*s(i3,i5)*za(i3,i5)*zab2(i3,i3,i5,i5)*zab2(i5,i3,i5,i2)*
     &    zb(i5,i6))/
     &  (za(i2,i6)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)**2*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i3,i5)*zab2(i3,i3,i5,i4)*zab2(i5,i3,i5,i2)*
     &    zb(i5,i6))/
     &  (za(i2,i6)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i5,i3,i5,i5)*zb(i4,i2)) - 
     & (2*s(i3,i5)**2*za(i3,i5)**2*zb(i2,i1)*zb(i5,i6))/
     &  (za(i2,i6)*zab2(i5,i3,i5,i1)**2*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)*zb(i4,i2)) + 
     & (4*s(i3,i5)**2*za(i2,i5)*za(i3,i5)**2*zb(i5,i2)*zb(i5,i6))/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)**2*zb(i4,i1)) + 
     & (2*s(i3,i5)**2*za(i3,i5)**2*za(i6,i5)*zb(i5,i6)**2)/
     &  (za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)*zab2(i5,i3,i5,i4)*
     &    zab2(i5,i3,i5,i5)**2*zb(i4,i1)) - 
     & (2*za(i1,i4)*za(i2,i3)**2*zab2(i1,i3,i5,i2)*zb(i6,i1))/
     &  (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i2,i3,i5,i2)*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i3,i5)*zab2(i1,i3,i5,i5)*zab2(i3,i3,i5,i5)*
     &    zb(i6,i1))/
     &  (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)*zab2(i5,i3,i5,i5)**2*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i2,i3)*za(i3,i5)*zab2(i1,i3,i5,i5)*
     &    zb(i6,i1))/
     &  (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i5,i3,i5,i5)*zb(i4,i1))
     &   + (2*s(i3,i5)*za(i1,i2)*za(i2,i3)**2*zb(i6,i2))/
     &  (za(i1,i6)*za(i2,i5)**2*zab2(i2,i3,i5,i2)*zab2(i2,i3,i5,i4)*
     &    zb(i4,i1)) - (2*s(i3,i5)*za(i3,i6)**2*zb(i6,i2))/
     &  (za(i6,i5)**2*zab2(i6,i3,i5,i4)*zab2(i6,i3,i5,i6)*zb(i4,i1)) - 
     & (2*s(i3,i5)*za(i2,i6)*za(i3,i6)**2*zab2(i6,i3,i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)**2*zaa2(i6,i1,i4,i2,i6,i6)*zab2(i6,i3,i5,i4)*
     &    zab2(i6,i3,i5,i6)*zb(i4,i1)) + 
     & (2*za(i2,i4)*zab2(i1,i3,i5,i2)*zab2(i3,i3,i5,i1)**2*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i1,i3,i5,i1)*zab2(i5,i3,i5,i1)**2*
     &    zb(i4,i2)) - (2*za(i4,i6)*zab2(i1,i3,i5,i6)*
     &    zab2(i3,i3,i5,i1)**2*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i1,i3,i5,i1)*zab2(i5,i3,i5,i1)**2*
     &    zb(i4,i2)) + (2*s(i3,i5)*za(i1,i6)*za(i3,i6)**2*
     &    zab2(i6,i3,i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)**2*zaa2(i6,i1,i6,i2,i4,i6)*zab2(i6,i3,i5,i4)*
     &    zab2(i6,i3,i5,i6)*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i3,i6)**2*zab2(i6,i3,i5,i2)*zb(i6,i2))/
     &  (za(i6,i5)**2*zab2(i6,i3,i5,i1)*zab2(i6,i3,i5,i4)*
     &    zab2(i6,i3,i5,i6)*zb(i4,i2)) - 
     & (2*s(i3,i5)*za(i2,i4)*za(i3,i5)*zab2(i3,i3,i5,i5)*zb(i5,i2)*
     &    zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i5)**2*
     &    zb(i4,i2)) - (2*s(i3,i5)*za(i2,i4)*za(i3,i5)*
     &    zab2(i3,i3,i5,i1)*zb(i5,i2)*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i5,i3,i5,i1)**2*zab2(i5,i3,i5,i5)*
     &    zb(i4,i2)) + (2*s(i3,i5)*za(i3,i5)*za(i4,i6)*
     &    zab2(i3,i3,i5,i5)*zb(i5,i6)*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i5,i3,i5,i1)*zab2(i5,i3,i5,i5)**2*
     &    zb(i4,i2)) + (2*s(i3,i5)*za(i3,i5)*za(i4,i6)*
     &    zab2(i3,i3,i5,i1)*zb(i5,i6)*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i6)*zab2(i5,i3,i5,i1)**2*zab2(i5,i3,i5,i5)*
     &    zb(i4,i2)) - (2*za(i1,i4)*za(i2,i3)**2*zab2(i4,i3,i5,i2)*
     &    zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i2,i3,i5,i2)*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i3,i5)*zab2(i3,i3,i5,i5)*zab2(i4,i3,i5,i5)*
     &    zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)*zab2(i5,i3,i5,i5)**2*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i2,i3)*za(i3,i5)*zab2(i4,i3,i5,i5)*
     &    zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i1,i6)*za(i2,i5)**2*zab2(i5,i3,i5,i5)*zb(i4,i1))
     &   - (2*zab2(i1,i3,i5,i4)**2*zab2(i3,i3,i5,i4)**2*
     &    zab2(i4,i3,i5,i2))/
     &  (za(i1,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i6,i3,i5,i4)*zbb2(i4,i2,i4,i1,i6,i4)) + 
     & ((sqrtd1624*za(i1,i5) + zaa2(i1,i1,i6,i2,i4,i5) - 
     &      zaa2(i1,i2,i4,i1,i6,i5))*
     &    (sqrtd1624*za(i1,i6) + zaa2(i1,i1,i6,i2,i4,i6) - 
     &      zaa2(i1,i2,i4,i1,i6,i6))*
     &    (sqrtd1624*za(i3,i5) + zaa2(i3,i1,i6,i2,i4,i5) - 
     &      zaa2(i3,i2,i4,i1,i6,i5))*
     &    (sqrtd1624*zab2(i3,i3,i5,i2) + 
     &      zab3(i3,i1,i6,i2,i4,i3,i5,i2) - 
     &      zab3(i3,i2,i4,i1,i6,i3,i5,i2))*
     &    (sqrtd1624*zb(i4,i2) + zbb2(i4,i1,i6,i2,i4,i2) - 
     &      zbb2(i4,i2,i4,i1,i6,i2)))/
     &  (16.*sqrtd1624**2*za(i1,i6)*zaa2(i5,i1,i6,i2,i4,i5)**2*
     &    zaa2(i6,i1,i6,i2,i4,i6)*zb(i4,i2)*zbb2(i4,i2,i4,i1,i6,i4)) - 
     & ((sqrtd1624*za(i1,i5) - zaa2(i1,i1,i6,i2,i4,i5) + 
     &      zaa2(i1,i2,i4,i1,i6,i5))*
     &    (sqrtd1624*za(i1,i6) - zaa2(i1,i1,i6,i2,i4,i6) + 
     &      zaa2(i1,i2,i4,i1,i6,i6))*
     &    (sqrtd1624*za(i3,i5) - zaa2(i3,i1,i6,i2,i4,i5) + 
     &      zaa2(i3,i2,i4,i1,i6,i5))*
     &    (sqrtd1624*zab2(i3,i3,i5,i2) - 
     &      zab3(i3,i1,i6,i2,i4,i3,i5,i2) + 
     &      zab3(i3,i2,i4,i1,i6,i3,i5,i2))*
     &    (sqrtd1624*zb(i4,i2) - zbb2(i4,i1,i6,i2,i4,i2) + 
     &      zbb2(i4,i2,i4,i1,i6,i2)))/
     &  (16.*sqrtd1624**2*za(i1,i6)*zaa2(i5,i1,i6,i2,i4,i5)**2*
     &    zaa2(i6,i1,i6,i2,i4,i6)*zb(i4,i2)*zbb2(i4,i2,i4,i1,i6,i4)) + 
     & (2*zab2(i2,i3,i5,i4)**2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*
     &    zb(i4,i2))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)**2*
     &    zab2(i6,i3,i5,i4)*zb(i4,i1)*zbb2(i4,i2,i6,i1,i4,i4)) - 
     & (4*zab2(i2,i3,i5,i4)*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i2)*
     &    zb(i6,i4))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i6,i1,i4,i4)) - 
     & (2*zab2(i3,i3,i5,i4)**2*zab2(i4,i3,i5,i6)*zab2(i6,i3,i5,i4)*
     &    zb(i6,i4))/
     &  (za(i2,i6)*zab2(i4,i3,i5,i4)*zab2(i5,i3,i5,i4)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i6,i1,i4,i4)) - 
     & ((sqrtd1426*za(i2,i5) + zaa2(i2,i1,i4,i2,i6,i5) - 
     &      zaa2(i2,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*za(i2,i6) + zaa2(i2,i1,i4,i2,i6,i6) - 
     &      zaa2(i2,i2,i6,i1,i4,i6))*
     &    (sqrtd1426*za(i3,i5) + zaa2(i3,i1,i4,i2,i6,i5) - 
     &      zaa2(i3,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*zab2(i3,i3,i5,i2) + 
     &      zab3(i3,i1,i4,i2,i6,i3,i5,i2) - 
     &      zab3(i3,i2,i6,i1,i4,i3,i5,i2))*
     &    (sqrtd1426*zb(i4,i2) + zbb2(i4,i1,i4,i2,i6,i2) - 
     &      zbb2(i4,i2,i6,i1,i4,i2)))/
     &  (16.*sqrtd1426**2*za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*
     &    zaa2(i6,i1,i4,i2,i6,i6)*zb(i4,i1)*zbb2(i4,i2,i6,i1,i4,i4)) + 
     & ((sqrtd1426*za(i2,i5) - zaa2(i2,i1,i4,i2,i6,i5) + 
     &      zaa2(i2,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*za(i2,i6) - zaa2(i2,i1,i4,i2,i6,i6) + 
     &      zaa2(i2,i2,i6,i1,i4,i6))*
     &    (sqrtd1426*za(i3,i5) - zaa2(i3,i1,i4,i2,i6,i5) + 
     &      zaa2(i3,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*zab2(i3,i3,i5,i2) - 
     &      zab3(i3,i1,i4,i2,i6,i3,i5,i2) + 
     &      zab3(i3,i2,i6,i1,i4,i3,i5,i2))*
     &    (sqrtd1426*zb(i4,i2) - zbb2(i4,i1,i4,i2,i6,i2) + 
     &      zbb2(i4,i2,i6,i1,i4,i2)))/
     &  (16.*sqrtd1426**2*za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*
     &    zaa2(i6,i1,i4,i2,i6,i6)*zb(i4,i1)*zbb2(i4,i2,i6,i1,i4,i4)) + 
     & ((sqrtd1426*za(i2,i5) - zaa2(i2,i1,i4,i2,i6,i5) + 
     &      zaa2(i2,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*za(i3,i5) - zaa2(i3,i1,i4,i2,i6,i5) + 
     &      zaa2(i3,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*zab2(i3,i3,i5,i2) - 
     &      zab3(i3,i1,i4,i2,i6,i3,i5,i2) + 
     &      zab3(i3,i2,i6,i1,i4,i3,i5,i2))*
     &    (sqrtd1426*zb(i6,i4) + zbb2(i6,i1,i4,i2,i6,i4) - 
     &      zbb2(i6,i2,i6,i1,i4,i4)))/
     &  (4.*sqrtd1426**2*za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i6,i1,i4,i4)) + 
     & ((sqrtd1426*za(i3,i5) - zaa2(i3,i1,i4,i2,i6,i5) + 
     &       zaa2(i3,i2,i6,i1,i4,i5))**2*
     &    (sqrtd1426*zab2(i6,i3,i5,i6) - 
     &      zab3(i6,i1,i4,i2,i6,i3,i5,i6) + 
     &      zab3(i6,i2,i6,i1,i4,i3,i5,i6))*
     &    (sqrtd1426*zb(i6,i4) + zbb2(i6,i1,i4,i2,i6,i4) - 
     &      zbb2(i6,i2,i6,i1,i4,i4)))/
     &  (8.*sqrtd1426**2*za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i6,i1,i4,i4)) + 
     & ((sqrtd1426*za(i2,i5) + zaa2(i2,i1,i4,i2,i6,i5) - 
     &      zaa2(i2,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*za(i3,i5) + zaa2(i3,i1,i4,i2,i6,i5) - 
     &      zaa2(i3,i2,i6,i1,i4,i5))*
     &    (sqrtd1426*zab2(i3,i3,i5,i2) + 
     &      zab3(i3,i1,i4,i2,i6,i3,i5,i2) - 
     &      zab3(i3,i2,i6,i1,i4,i3,i5,i2))*
     &    (sqrtd1426*zb(i6,i4) - zbb2(i6,i1,i4,i2,i6,i4) + 
     &      zbb2(i6,i2,i6,i1,i4,i4)))/
     &  (4.*sqrtd1426**2*za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i6,i1,i4,i4)) + 
     & ((sqrtd1426*za(i3,i5) + zaa2(i3,i1,i4,i2,i6,i5) - 
     &       zaa2(i3,i2,i6,i1,i4,i5))**2*
     &    (sqrtd1426*zab2(i6,i3,i5,i6) + 
     &      zab3(i6,i1,i4,i2,i6,i3,i5,i6) - 
     &      zab3(i6,i2,i6,i1,i4,i3,i5,i6))*
     &    (sqrtd1426*zb(i6,i4) - zbb2(i6,i1,i4,i2,i6,i4) + 
     &      zbb2(i6,i2,i6,i1,i4,i4)))/
     &  (8.*sqrtd1426**2*za(i2,i6)*zaa2(i5,i1,i4,i2,i6,i5)**2*zb(i4,i1)*
     &    zbb2(i4,i2,i6,i1,i4,i4))
     &)


      aaaa_NMHV_b(21) = (
     &  (2*s(i1,i3)*za(i2,i4)**2*zaa2(i1,i5,i6,i2,i4,i2)*
     &    zab2(i1,i2,i4,i2))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i2,i4,i2)*zab2(i2,i2,i4,i3)*zb(i3,i1)) + 
     & (2*za(i2,i4)**2*zaa2(i1,i5,i6,i2,i4,i2)*zaa2(i2,i1,i3,i2,i4,i5)*
     &    zab2(i1,i2,i4,i2))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)**2*za(i2,i6)*
     &    zab2(i2,i2,i4,i2)*zab2(i2,i2,i4,i3)*zb(i3,i1)) - 
     & (2*za(i2,i4)*zaa2(i1,i5,i6,i2,i4,i2)*zaa2(i2,i1,i3,i2,i4,i2)*
     &    zab2(i1,i2,i4,i2)*zab2(i4,i2,i4,i1))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i2,i4,i1)*zab2(i2,i2,i4,i2)*zab2(i2,i2,i4,i3)*
     &    zb(i3,i1)) - (2*za(i2,i4)*zaa2(i1,i5,i6,i2,i4,i2)*
     &    zaa2(i2,i1,i3,i2,i4,i2)*zab2(i4,i2,i4,i2))/
     &  (s(i2,i4)*za(i1,i5)*za(i2,i5)*za(i2,i6)**2*zab2(i2,i2,i4,i2)*
     &    zab2(i2,i2,i4,i3)*zb(i3,i1)) + 
     & (za(i1,i2)*zaa2(i1,i5,i6,i2,i4,i2)*zaa2(i2,i1,i3,i2,i4,i2)*
     &    zab2(i4,i2,i4,i2)**2)/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i2,i4,i2)**2*zab2(i2,i2,i4,i3)*zb(i3,i1)) + 
     & (za(i1,i3)*zab2(i2,i1,i3,i6)*zab2(i2,i2,i4,i5)*
     &    zab2(i4,i2,i4,i2)**2)/
     &  (s(i2,i4)*t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*zab2(i2,i2,i4,i2)**2*
     &    zb(i3,i1)) + (za(i1,i3)*zab2(i2,i1,i3,i5)*zab2(i2,i2,i4,i6)*
     &    zab2(i4,i2,i4,i2)**2)/
     &  (s(i2,i4)*t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*zab2(i2,i2,i4,i2)**2*
     &    zb(i3,i1)) - (zaa2(i1,i5,i6,i2,i4,i5)*zaa2(i5,i1,i3,i2,i4,i5)*
     &    zab2(i4,i2,i4,i5)**2)/
     &  (s(i2,i4)*za(i1,i6)*za(i2,i5)*za(i5,i6)*zab2(i5,i2,i4,i3)*
     &    zab2(i5,i2,i4,i5)**2*zb(i3,i1)) - 
     & (za(i1,i3)*zab2(i4,i2,i4,i5)**2*zab2(i5,i1,i3,i6))/
     &  (s(i2,i4)*t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*zab2(i5,i2,i4,i5)*
     &    zb(i3,i1)) + (2*za(i4,i5)**2*zaa2(i5,i1,i3,i2,i4,i5)*
     &    zab2(i1,i2,i4,i5)*zab2(i1,i5,i6,i3))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*zab2(i5,i2,i4,i3)**2*
     &    zab2(i5,i2,i4,i5)*zb(i3,i1)) - 
     & (2*s(i1,i3)*za(i4,i5)**2*zaa2(i1,i5,i6,i2,i4,i5)*
     &    zab2(i1,i2,i4,i5))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*
     &    zab2(i5,i2,i4,i3)*zab2(i5,i2,i4,i5)*zb(i3,i1)) + 
     & (2*za(i4,i5)*za(i4,i6)*zaa2(i1,i5,i6,i2,i4,i5)*
     &    zaa2(i5,i1,i3,i2,i4,i5)*zab2(i1,i2,i4,i5))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)**2*
     &    zab2(i5,i2,i4,i3)*zab2(i5,i2,i4,i5)*zb(i3,i1)) + 
     & (2*za(i2,i4)*za(i4,i5)*zaa2(i1,i5,i6,i2,i4,i5)*
     &    zaa2(i5,i1,i3,i2,i4,i5)*zab2(i1,i2,i4,i5))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)**2*za(i5,i6)*
     &    zab2(i5,i2,i4,i3)*zab2(i5,i2,i4,i5)*zb(i3,i1)) - 
     & (2*za(i4,i5)*zaa2(i1,i5,i6,i2,i4,i5)*zaa2(i5,i1,i3,i2,i4,i5)*
     &    zab2(i4,i2,i4,i5))/
     &  (s(i2,i4)*za(i1,i5)*za(i2,i5)*za(i5,i6)**2*zab2(i5,i2,i4,i3)*
     &    zab2(i5,i2,i4,i5)*zb(i3,i1)) - 
     & (2*za(i4,i5)*zaa2(i1,i5,i6,i2,i4,i5)*zaa2(i5,i1,i3,i2,i4,i5)*
     &    zab2(i1,i2,i4,i5)*zab2(i4,i2,i4,i1))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*
     &    zab2(i5,i2,i4,i1)*zab2(i5,i2,i4,i3)*zab2(i5,i2,i4,i5)*
     &    zb(i3,i1)) + (2*s(i1,i3)*za(i1,i3)*zab2(i1,i5,i6,i3)*
     &    zab2(i4,i2,i4,i3)**2)/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i4,i3)*zab2(i3,i2,i4,i3)*
     &    zab2(i5,i2,i4,i3)*zab2(i6,i2,i4,i3)*zb(i3,i1)) + 
     & (zaa2(i1,i5,i6,i2,i4,i6)*zaa2(i6,i1,i3,i2,i4,i6)*
     &    zab2(i4,i2,i4,i6)**2)/
     &  (s(i2,i4)*za(i1,i5)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i4,i3)*
     &    zab2(i6,i2,i4,i6)**2*zb(i3,i1)) - 
     & (za(i1,i3)*zab2(i4,i2,i4,i6)**2*zab2(i6,i1,i3,i5))/
     &  (s(i2,i4)*t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*zab2(i6,i2,i4,i6)*
     &    zb(i3,i1)) - (2*za(i4,i6)*zaa2(i1,i5,i6,i2,i4,i6)*
     &    zab2(i4,i2,i4,i6)*zab2(i6,i1,i3,i3))/
     &  (za(i1,i5)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i4,i3)**2*
     &    zab2(i6,i2,i4,i6)*zb(i3,i1)) + 
     & (2*s(i1,i3)*za(i4,i6)**2*zaa2(i1,i5,i6,i2,i4,i6)*
     &    zab2(i1,i2,i4,i6))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i5,i6)*
     &    zab2(i6,i2,i4,i3)*zab2(i6,i2,i4,i6)*zb(i3,i1)) + 
     & (2*za(i4,i6)**2*zaa2(i1,i5,i6,i2,i4,i6)*zaa2(i5,i2,i4,i1,i3,i6)*
     &    zab2(i1,i2,i4,i6))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i5,i6)**2*
     &    zab2(i6,i2,i4,i3)*zab2(i6,i2,i4,i6)*zb(i3,i1)) + 
     & (2*za(i4,i6)*zaa2(i1,i5,i6,i2,i4,i5)*zaa2(i6,i1,i3,i2,i4,i6)*
     &    zab2(i4,i2,i4,i6))/
     &  (s(i2,i4)*za(i1,i5)*za(i2,i6)*za(i5,i6)**2*zab2(i6,i2,i4,i3)*
     &    zab2(i6,i2,i4,i6)*zb(i3,i1)) - 
     & (2*za(i2,i4)*zaa2(i1,i5,i6,i2,i4,i6)*zaa2(i6,i1,i3,i2,i4,i6)*
     &    zab2(i4,i2,i4,i6))/
     &  (s(i2,i4)*za(i1,i5)*za(i2,i6)**2*za(i5,i6)*zab2(i6,i2,i4,i3)*
     &    zab2(i6,i2,i4,i6)*zb(i3,i1)) + 
     & (2*za(i4,i6)*zaa2(i1,i5,i6,i2,i4,i6)*zaa2(i6,i1,i3,i2,i4,i6)*
     &    zab2(i1,i2,i4,i6)*zab2(i4,i2,i4,i1))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i5,i6)*
     &    zab2(i6,i2,i4,i1)*zab2(i6,i2,i4,i3)*zab2(i6,i2,i4,i6)*
     &    zb(i3,i1)) + (2*za(i1,i3)*zab2(i1,i5,i6,i3)*
     &    zab2(i4,i2,i4,i3)**2*zab3(i5,i2,i4,i1,i3,i2,i4,i3))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i4,i3)*zab2(i3,i2,i4,i3)*
     &    zab2(i5,i2,i4,i3)**2*zab2(i6,i2,i4,i3)*zb(i3,i1)) - 
     & (2*za(i4,i5)**2*zaa2(i1,i5,i6,i2,i4,i5)*zab2(i1,i2,i4,i5)*
     &    zab3(i5,i2,i4,i1,i3,i2,i4,i5))/
     &  (s(i2,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*
     &    zab2(i5,i2,i4,i3)*zab2(i5,i2,i4,i5)**2*zb(i3,i1)) + 
     & (2*za(i4,i6)*zaa2(i1,i5,i6,i2,i4,i6)*zab2(i4,i2,i4,i6)*
     &    zab3(i6,i2,i4,i1,i3,i2,i4,i6))/
     &  (s(i2,i4)*za(i1,i5)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i4,i3)*
     &    zab2(i6,i2,i4,i6)**2*zb(i3,i1)) + 
     & (2*s(i2,i4)*za(i1,i2)*za(i2,i3)*za(i2,i4)**2*zab2(i2,i2,i4,i6)*
     &    zb(i5,i2))/
     &  (za(i1,i5)*za(i2,i6)*zaa2(i2,i1,i5,i3,i6,i2)*zab2(i2,i2,i4,i1)*
     &    zab2(i2,i2,i4,i2)*zab2(i2,i2,i4,i3)) + 
     & (2*s(i2,i4)*za(i1,i2)*za(i2,i3)*za(i2,i4)**2*zab2(i2,i2,i4,i6)*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i5)*zaa2(i2,i1,i6,i3,i5,i2)*zab2(i2,i2,i4,i1)*
     &    zab2(i2,i2,i4,i2)*zab2(i2,i2,i4,i3)) - 
     & (2*za(i1,i3)*za(i2,i4)*za(i4,i5)*zab2(i2,i1,i3,i6)*zb(i5,i2))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)**2*zab2(i2,i2,i4,i2)*zb(i3,i1))
     &   - (2*za(i1,i3)*za(i2,i4)*za(i4,i6)*zab2(i2,i1,i3,i5)*
     &    zb(i6,i2))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)**2*zab2(i2,i2,i4,i2)*zb(i3,i1))
     &   - (2*s(i2,i4)*za(i1,i6)*za(i3,i6)*za(i4,i6)**2*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zaa2(i6,i1,i5,i3,i6,i6)*zab2(i6,i2,i4,i1)*
     &    zab2(i6,i2,i4,i3)) - 
     & (2*zab2(i1,i2,i4,i5)*zab2(i3,i2,i4,i1)*zab2(i4,i2,i4,i1)**2*
     &    zb(i6,i1))/
     &  (s(i2,i4)*za(i1,i6)*zab2(i2,i2,i4,i1)*zab2(i5,i2,i4,i1)*
     &    zb(i3,i1)*zbb2(i1,i3,i5,i1,i6,i1)) - 
     & (2*zab2(i1,i2,i4,i5)*zab2(i3,i2,i4,i1)*zab2(i4,i2,i4,i1)**2*
     &    zb(i6,i1))/
     &  (s(i2,i4)*za(i1,i5)*zab2(i2,i2,i4,i1)*zab2(i6,i2,i4,i1)*
     &    zb(i3,i1)*zbb2(i1,i3,i6,i1,i5,i1)) + 
     & (2*za(i1,i3)*zab2(i1,i5,i6,i3)*zab2(i4,i2,i4,i1)*
     &    zab2(i4,i2,i4,i3)*zbb2(i3,i2,i4,i1,i3,i3))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i4,i3)*zab2(i3,i2,i4,i3)*
     &    zab2(i5,i2,i4,i3)*zab2(i6,i2,i4,i3)*zb(i3,i1)**2) - 
     & (2*s(i2,i4)*za(i3,i4)*zab2(i1,i5,i6,i3)*zab2(i4,i2,i4,i3)*
     &    zbb2(i3,i2,i4,i1,i3,i3))/
     &  (za(i1,i5)*zab2(i2,i2,i4,i3)*zab2(i3,i2,i4,i3)*
     &    zab2(i5,i2,i4,i3)*zab2(i6,i2,i4,i3)**2*zb(i3,i1)) + 
     & (s(i2,i4)*za(i3,i4)**2*zab2(i1,i2,i4,i3)*zab2(i1,i5,i6,i3)*
     &    zbb2(i3,i2,i4,i1,i3,i3))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i4,i3)*zab2(i3,i2,i4,i3)**2*
     &    zab2(i5,i2,i4,i3)*zab2(i6,i2,i4,i3)*zb(i3,i1)) + 
     & (2*zab2(i1,i2,i4,i3)*zab2(i3,i2,i4,i5)*zab2(i4,i2,i4,i3)**2*
     &    zb(i6,i3))/
     &  (s(i2,i4)*za(i1,i6)*zab2(i2,i2,i4,i3)*zab2(i5,i2,i4,i3)*
     &    zb(i3,i1)*zbb2(i3,i3,i5,i1,i6,i3)) + 
     & ((sqrtd1635*za(i4,i5) + zaa2(i4,i1,i6,i3,i5,i5) - 
     &      zaa2(i4,i3,i5,i1,i6,i5))*
     &    (sqrtd1635*zab2(i1,i2,i4,i1) + 
     &      zab3(i1,i1,i6,i3,i5,i2,i4,i1) - 
     &      zab3(i1,i3,i5,i1,i6,i2,i4,i1))*
     &    (sqrtd1635*zab2(i2,i2,i4,i6) - 
     &      zab3(i2,i1,i6,i3,i5,i2,i4,i6) + 
     &      zab3(i2,i3,i5,i1,i6,i2,i4,i6))*
     &    (sqrtd1635*zab2(i3,i2,i4,i3) + 
     &      zab3(i3,i1,i6,i3,i5,i2,i4,i3) - 
     &      zab3(i3,i3,i5,i1,i6,i2,i4,i3))*
     &    (sqrtd1635*zab2(i4,i2,i4,i5) + 
     &      zab3(i4,i1,i6,i3,i5,i2,i4,i5) - 
     &      zab3(i4,i3,i5,i1,i6,i2,i4,i5)))/
     &  (16.*sqrtd1635**2*s(i2,i4)**2*za(i1,i6)*zaa2(i2,i1,i6,i3,i5,i2)*
     &    zaa2(i5,i1,i6,i3,i5,i5)*zbb2(i1,i3,i5,i1,i6,i1)*
     &    zbb2(i3,i3,i5,i1,i6,i3)) - 
     & ((sqrtd1635*za(i4,i5) - zaa2(i4,i1,i6,i3,i5,i5) + 
     &      zaa2(i4,i3,i5,i1,i6,i5))*
     &    (sqrtd1635*zab2(i1,i2,i4,i1) - 
     &      zab3(i1,i1,i6,i3,i5,i2,i4,i1) + 
     &      zab3(i1,i3,i5,i1,i6,i2,i4,i1))*
     &    (sqrtd1635*zab2(i2,i2,i4,i6) + 
     &      zab3(i2,i1,i6,i3,i5,i2,i4,i6) - 
     &      zab3(i2,i3,i5,i1,i6,i2,i4,i6))*
     &    (sqrtd1635*zab2(i3,i2,i4,i3) - 
     &      zab3(i3,i1,i6,i3,i5,i2,i4,i3) + 
     &      zab3(i3,i3,i5,i1,i6,i2,i4,i3))*
     &    (sqrtd1635*zab2(i4,i2,i4,i5) - 
     &      zab3(i4,i1,i6,i3,i5,i2,i4,i5) + 
     &      zab3(i4,i3,i5,i1,i6,i2,i4,i5)))/
     &  (16.*sqrtd1635**2*s(i2,i4)**2*za(i1,i6)*zaa2(i2,i1,i6,i3,i5,i2)*
     &    zaa2(i5,i1,i6,i3,i5,i5)*zbb2(i1,i3,i5,i1,i6,i1)*
     &    zbb2(i3,i3,i5,i1,i6,i3)) + 
     & (2*zab2(i1,i2,i4,i3)*zab2(i3,i2,i4,i5)*zab2(i4,i2,i4,i3)**2*
     &    zb(i6,i3))/
     &  (s(i2,i4)*za(i1,i5)*zab2(i2,i2,i4,i3)*zab2(i6,i2,i4,i3)*
     &    zb(i3,i1)*zbb2(i3,i3,i6,i1,i5,i3)) + 
     & ((sqrtd1536*za(i2,i4) + zaa2(i2,i1,i5,i3,i6,i4) - 
     &      zaa2(i2,i3,i6,i1,i5,i4))*
     &    (sqrtd1536*zab2(i1,i2,i4,i1) - 
     &      zab3(i1,i1,i5,i3,i6,i2,i4,i1) + 
     &      zab3(i1,i3,i6,i1,i5,i2,i4,i1))*
     &    (sqrtd1536*zab2(i3,i2,i4,i3) - 
     &      zab3(i3,i1,i5,i3,i6,i2,i4,i3) + 
     &      zab3(i3,i3,i6,i1,i5,i2,i4,i3))*
     &    (sqrtd1536*zab2(i4,i2,i4,i5) - 
     &      zab3(i4,i1,i5,i3,i6,i2,i4,i5) + 
     &      zab3(i4,i3,i6,i1,i5,i2,i4,i5))*
     &    (sqrtd1536*zab2(i6,i2,i4,i6) + 
     &      zab3(i6,i1,i5,i3,i6,i2,i4,i6) - 
     &      zab3(i6,i3,i6,i1,i5,i2,i4,i6)))/
     &  (16.*sqrtd1536**2*s(i2,i4)**2*za(i1,i5)*zaa2(i2,i1,i5,i3,i6,i2)*
     &    zaa2(i6,i1,i5,i3,i6,i6)*zbb2(i1,i3,i6,i1,i5,i1)*
     &    zbb2(i3,i3,i6,i1,i5,i3)) - 
     & ((sqrtd1536*za(i2,i4) - zaa2(i2,i1,i5,i3,i6,i4) + 
     &      zaa2(i2,i3,i6,i1,i5,i4))*
     &    (sqrtd1536*zab2(i1,i2,i4,i1) + 
     &      zab3(i1,i1,i5,i3,i6,i2,i4,i1) - 
     &      zab3(i1,i3,i6,i1,i5,i2,i4,i1))*
     &    (sqrtd1536*zab2(i3,i2,i4,i3) + 
     &      zab3(i3,i1,i5,i3,i6,i2,i4,i3) - 
     &      zab3(i3,i3,i6,i1,i5,i2,i4,i3))*
     &    (sqrtd1536*zab2(i4,i2,i4,i5) + 
     &      zab3(i4,i1,i5,i3,i6,i2,i4,i5) - 
     &      zab3(i4,i3,i6,i1,i5,i2,i4,i5))*
     &    (sqrtd1536*zab2(i6,i2,i4,i6) - 
     &      zab3(i6,i1,i5,i3,i6,i2,i4,i6) + 
     &      zab3(i6,i3,i6,i1,i5,i2,i4,i6)))/
     &  (16.*sqrtd1536**2*s(i2,i4)**2*za(i1,i5)*zaa2(i2,i1,i5,i3,i6,i2)*
     &    zaa2(i6,i1,i5,i3,i6,i6)*zbb2(i1,i3,i6,i1,i5,i1)*
     &    zbb2(i3,i3,i6,i1,i5,i3))
     &)


      aaaa_NMHV_b(23) = (
     &  (2*s(i1,i4)*za(i2,i3)**2*zaa2(i1,i5,i6,i2,i3,i2)*
     &    zab2(i1,i2,i3,i2))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i2,i3,i2)*zab2(i2,i2,i3,i4)*zb(i4,i1)) + 
     & (2*za(i2,i3)**2*zaa2(i1,i5,i6,i2,i3,i2)*zaa2(i2,i1,i4,i2,i3,i5)*
     &    zab2(i1,i2,i3,i2))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)**2*za(i2,i6)*
     &    zab2(i2,i2,i3,i2)*zab2(i2,i2,i3,i4)*zb(i4,i1)) - 
     & (2*za(i2,i3)*zaa2(i1,i5,i6,i2,i3,i2)*zaa2(i2,i1,i4,i2,i3,i2)*
     &    zab2(i1,i2,i3,i2)*zab2(i3,i2,i3,i1))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i2,i3,i1)*zab2(i2,i2,i3,i2)*zab2(i2,i2,i3,i4)*
     &    zb(i4,i1)) - (2*za(i2,i3)*zaa2(i1,i5,i6,i2,i3,i2)*
     &    zaa2(i2,i1,i4,i2,i3,i2)*zab2(i3,i2,i3,i2))/
     &  (s(i2,i3)*za(i1,i5)*za(i2,i5)*za(i2,i6)**2*zab2(i2,i2,i3,i2)*
     &    zab2(i2,i2,i3,i4)*zb(i4,i1)) + 
     & (za(i1,i2)*zaa2(i1,i5,i6,i2,i3,i2)*zaa2(i2,i1,i4,i2,i3,i2)*
     &    zab2(i3,i2,i3,i2)**2)/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i2,i3,i2)**2*zab2(i2,i2,i3,i4)*zb(i4,i1)) + 
     & (za(i1,i4)*zab2(i2,i1,i4,i6)*zab2(i2,i2,i3,i5)*
     &    zab2(i3,i2,i3,i2)**2)/
     &  (s(i2,i3)*t(i1,i4,i6)*za(i1,i6)*za(i2,i5)*zab2(i2,i2,i3,i2)**2*
     &    zb(i4,i1)) + (za(i1,i4)*zab2(i2,i1,i4,i5)*zab2(i2,i2,i3,i6)*
     &    zab2(i3,i2,i3,i2)**2)/
     &  (s(i2,i3)*t(i1,i4,i5)*za(i1,i5)*za(i2,i6)*zab2(i2,i2,i3,i2)**2*
     &    zb(i4,i1)) - (zaa2(i1,i5,i6,i2,i3,i5)*zaa2(i5,i1,i4,i2,i3,i5)*
     &    zab2(i3,i2,i3,i5)**2)/
     &  (s(i2,i3)*za(i1,i6)*za(i2,i5)*za(i5,i6)*zab2(i5,i2,i3,i4)*
     &    zab2(i5,i2,i3,i5)**2*zb(i4,i1)) - 
     & (za(i1,i4)*zab2(i3,i2,i3,i5)**2*zab2(i5,i1,i4,i6))/
     &  (s(i2,i3)*t(i1,i4,i6)*za(i1,i6)*za(i2,i5)*zab2(i5,i2,i3,i5)*
     &    zb(i4,i1)) + (2*za(i3,i5)**2*zaa2(i5,i1,i4,i2,i3,i5)*
     &    zab2(i1,i2,i3,i5)*zab2(i1,i5,i6,i4))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*zab2(i5,i2,i3,i4)**2*
     &    zab2(i5,i2,i3,i5)*zb(i4,i1)) - 
     & (2*s(i1,i4)*za(i3,i5)**2*zaa2(i1,i5,i6,i2,i3,i5)*
     &    zab2(i1,i2,i3,i5))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*
     &    zab2(i5,i2,i3,i4)*zab2(i5,i2,i3,i5)*zb(i4,i1)) + 
     & (2*za(i3,i5)*za(i3,i6)*zaa2(i1,i5,i6,i2,i3,i5)*
     &    zaa2(i5,i1,i4,i2,i3,i5)*zab2(i1,i2,i3,i5))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)**2*
     &    zab2(i5,i2,i3,i4)*zab2(i5,i2,i3,i5)*zb(i4,i1)) + 
     & (2*za(i2,i3)*za(i3,i5)*zaa2(i1,i5,i6,i2,i3,i5)*
     &    zaa2(i5,i1,i4,i2,i3,i5)*zab2(i1,i2,i3,i5))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)**2*za(i5,i6)*
     &    zab2(i5,i2,i3,i4)*zab2(i5,i2,i3,i5)*zb(i4,i1)) - 
     & (2*za(i3,i5)*zaa2(i1,i5,i6,i2,i3,i5)*zaa2(i5,i1,i4,i2,i3,i5)*
     &    zab2(i3,i2,i3,i5))/
     &  (s(i2,i3)*za(i1,i5)*za(i2,i5)*za(i5,i6)**2*zab2(i5,i2,i3,i4)*
     &    zab2(i5,i2,i3,i5)*zb(i4,i1)) - 
     & (2*za(i3,i5)*zaa2(i1,i5,i6,i2,i3,i5)*zaa2(i5,i1,i4,i2,i3,i5)*
     &    zab2(i1,i2,i3,i5)*zab2(i3,i2,i3,i1))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*
     &    zab2(i5,i2,i3,i1)*zab2(i5,i2,i3,i4)*zab2(i5,i2,i3,i5)*
     &    zb(i4,i1)) + (2*s(i1,i4)*za(i1,i4)*zab2(i1,i5,i6,i4)*
     &    zab2(i3,i2,i3,i4)**2)/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i3,i4)*zab2(i4,i2,i3,i4)*
     &    zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)*zb(i4,i1)) + 
     & (zaa2(i1,i5,i6,i2,i3,i6)*zaa2(i6,i1,i4,i2,i3,i6)*
     &    zab2(i3,i2,i3,i6)**2)/
     &  (s(i2,i3)*za(i1,i5)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i3,i4)*
     &    zab2(i6,i2,i3,i6)**2*zb(i4,i1)) - 
     & (za(i1,i4)*zab2(i3,i2,i3,i6)**2*zab2(i6,i1,i4,i5))/
     &  (s(i2,i3)*t(i1,i4,i5)*za(i1,i5)*za(i2,i6)*zab2(i6,i2,i3,i6)*
     &    zb(i4,i1)) - (2*za(i3,i6)*zaa2(i1,i5,i6,i2,i3,i6)*
     &    zab2(i3,i2,i3,i6)*zab2(i6,i1,i4,i4))/
     &  (za(i1,i5)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i3,i4)**2*
     &    zab2(i6,i2,i3,i6)*zb(i4,i1)) + 
     & (2*s(i1,i4)*za(i3,i6)**2*zaa2(i1,i5,i6,i2,i3,i6)*
     &    zab2(i1,i2,i3,i6))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i5,i6)*
     &    zab2(i6,i2,i3,i4)*zab2(i6,i2,i3,i6)*zb(i4,i1)) + 
     & (2*za(i3,i6)**2*zaa2(i1,i5,i6,i2,i3,i6)*zaa2(i5,i2,i3,i1,i4,i6)*
     &    zab2(i1,i2,i3,i6))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i5,i6)**2*
     &    zab2(i6,i2,i3,i4)*zab2(i6,i2,i3,i6)*zb(i4,i1)) + 
     & (2*za(i3,i6)*zaa2(i1,i5,i6,i2,i3,i5)*zaa2(i6,i1,i4,i2,i3,i6)*
     &    zab2(i3,i2,i3,i6))/
     &  (s(i2,i3)*za(i1,i5)*za(i2,i6)*za(i5,i6)**2*zab2(i6,i2,i3,i4)*
     &    zab2(i6,i2,i3,i6)*zb(i4,i1)) - 
     & (2*za(i2,i3)*zaa2(i1,i5,i6,i2,i3,i6)*zaa2(i6,i1,i4,i2,i3,i6)*
     &    zab2(i3,i2,i3,i6))/
     &  (s(i2,i3)*za(i1,i5)*za(i2,i6)**2*za(i5,i6)*zab2(i6,i2,i3,i4)*
     &    zab2(i6,i2,i3,i6)*zb(i4,i1)) + 
     & (2*za(i3,i6)*zaa2(i1,i5,i6,i2,i3,i6)*zaa2(i6,i1,i4,i2,i3,i6)*
     &    zab2(i1,i2,i3,i6)*zab2(i3,i2,i3,i1))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i5,i6)*
     &    zab2(i6,i2,i3,i1)*zab2(i6,i2,i3,i4)*zab2(i6,i2,i3,i6)*
     &    zb(i4,i1)) + (2*za(i1,i4)*zab2(i1,i5,i6,i4)*
     &    zab2(i3,i2,i3,i4)**2*zab3(i5,i2,i3,i1,i4,i2,i3,i4))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i3,i4)*zab2(i4,i2,i3,i4)*
     &    zab2(i5,i2,i3,i4)**2*zab2(i6,i2,i3,i4)*zb(i4,i1)) - 
     & (2*za(i3,i5)**2*zaa2(i1,i5,i6,i2,i3,i5)*zab2(i1,i2,i3,i5)*
     &    zab3(i5,i2,i3,i1,i4,i2,i3,i5))/
     &  (s(i2,i3)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i5,i6)*
     &    zab2(i5,i2,i3,i4)*zab2(i5,i2,i3,i5)**2*zb(i4,i1)) + 
     & (2*za(i3,i6)*zaa2(i1,i5,i6,i2,i3,i6)*zab2(i3,i2,i3,i6)*
     &    zab3(i6,i2,i3,i1,i4,i2,i3,i6))/
     &  (s(i2,i3)*za(i1,i5)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i3,i4)*
     &    zab2(i6,i2,i3,i6)**2*zb(i4,i1)) + 
     & (2*s(i2,i3)*za(i1,i2)*za(i2,i3)**2*za(i2,i4)*zab2(i2,i2,i3,i6)*
     &    zb(i5,i2))/
     &  (za(i1,i5)*za(i2,i6)*zaa2(i2,i1,i5,i4,i6,i2)*zab2(i2,i2,i3,i1)*
     &    zab2(i2,i2,i3,i2)*zab2(i2,i2,i3,i4)) + 
     & (2*s(i2,i3)*za(i1,i2)*za(i2,i3)**2*za(i2,i4)*zab2(i2,i2,i3,i6)*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i5)*zaa2(i2,i1,i6,i4,i5,i2)*zab2(i2,i2,i3,i1)*
     &    zab2(i2,i2,i3,i2)*zab2(i2,i2,i3,i4)) - 
     & (2*za(i1,i4)*za(i2,i3)*za(i3,i5)*zab2(i2,i1,i4,i6)*zb(i5,i2))/
     &  (t(i1,i4,i6)*za(i1,i6)*za(i2,i5)**2*zab2(i2,i2,i3,i2)*zb(i4,i1))
     &   - (2*za(i1,i4)*za(i2,i3)*za(i3,i6)*zab2(i2,i1,i4,i5)*
     &    zb(i6,i2))/
     &  (t(i1,i4,i5)*za(i1,i5)*za(i2,i6)**2*zab2(i2,i2,i3,i2)*zb(i4,i1))
     &   - (2*s(i2,i3)*za(i1,i6)*za(i3,i6)**2*za(i4,i6)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*zaa2(i6,i1,i5,i4,i6,i6)*zab2(i6,i2,i3,i1)*
     &    zab2(i6,i2,i3,i4)) - 
     & (2*zab2(i1,i2,i3,i5)*zab2(i3,i2,i3,i1)**2*zab2(i4,i2,i3,i1)*
     &    zb(i6,i1))/
     &  (s(i2,i3)*za(i1,i6)*zab2(i2,i2,i3,i1)*zab2(i5,i2,i3,i1)*
     &    zb(i4,i1)*zbb2(i1,i4,i5,i1,i6,i1)) - 
     & (2*zab2(i1,i2,i3,i5)*zab2(i3,i2,i3,i1)**2*zab2(i4,i2,i3,i1)*
     &    zb(i6,i1))/
     &  (s(i2,i3)*za(i1,i5)*zab2(i2,i2,i3,i1)*zab2(i6,i2,i3,i1)*
     &    zb(i4,i1)*zbb2(i1,i4,i6,i1,i5,i1)) + 
     & (2*za(i1,i4)*zab2(i1,i5,i6,i4)*zab2(i3,i2,i3,i1)*
     &    zab2(i3,i2,i3,i4)*zbb2(i4,i2,i3,i1,i4,i4))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i3,i4)*zab2(i4,i2,i3,i4)*
     &    zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)*zb(i4,i1)**2) - 
     & (2*s(i2,i3)*za(i4,i3)*zab2(i1,i5,i6,i4)*zab2(i3,i2,i3,i4)*
     &    zbb2(i4,i2,i3,i1,i4,i4))/
     &  (za(i1,i5)*zab2(i2,i2,i3,i4)*zab2(i4,i2,i3,i4)*
     &    zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)**2*zb(i4,i1)) + 
     & (s(i2,i3)*za(i4,i3)**2*zab2(i1,i2,i3,i4)*zab2(i1,i5,i6,i4)*
     &    zbb2(i4,i2,i3,i1,i4,i4))/
     &  (za(i1,i5)*za(i1,i6)*zab2(i2,i2,i3,i4)*zab2(i4,i2,i3,i4)**2*
     &    zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)*zb(i4,i1)) + 
     & (2*zab2(i1,i2,i3,i4)*zab2(i3,i2,i3,i4)**2*zab2(i4,i2,i3,i5)*
     &    zb(i6,i4))/
     &  (s(i2,i3)*za(i1,i6)*zab2(i2,i2,i3,i4)*zab2(i5,i2,i3,i4)*
     &    zb(i4,i1)*zbb2(i4,i4,i5,i1,i6,i4)) + 
     & ((sqrtd1645*za(i3,i5) + zaa2(i3,i1,i6,i4,i5,i5) - 
     &      zaa2(i3,i4,i5,i1,i6,i5))*
     &    (sqrtd1645*zab2(i1,i2,i3,i1) + 
     &      zab3(i1,i1,i6,i4,i5,i2,i3,i1) - 
     &      zab3(i1,i4,i5,i1,i6,i2,i3,i1))*
     &    (sqrtd1645*zab2(i2,i2,i3,i6) - 
     &      zab3(i2,i1,i6,i4,i5,i2,i3,i6) + 
     &      zab3(i2,i4,i5,i1,i6,i2,i3,i6))*
     &    (sqrtd1645*zab2(i3,i2,i3,i5) + 
     &      zab3(i3,i1,i6,i4,i5,i2,i3,i5) - 
     &      zab3(i3,i4,i5,i1,i6,i2,i3,i5))*
     &    (sqrtd1645*zab2(i4,i2,i3,i4) + 
     &      zab3(i4,i1,i6,i4,i5,i2,i3,i4) - 
     &      zab3(i4,i4,i5,i1,i6,i2,i3,i4)))/
     &  (16.*sqrtd1645**2*s(i2,i3)**2*za(i1,i6)*zaa2(i2,i1,i6,i4,i5,i2)*
     &    zaa2(i5,i1,i6,i4,i5,i5)*zbb2(i1,i4,i5,i1,i6,i1)*
     &    zbb2(i4,i4,i5,i1,i6,i4)) - 
     & ((sqrtd1645*za(i3,i5) - zaa2(i3,i1,i6,i4,i5,i5) + 
     &      zaa2(i3,i4,i5,i1,i6,i5))*
     &    (sqrtd1645*zab2(i1,i2,i3,i1) - 
     &      zab3(i1,i1,i6,i4,i5,i2,i3,i1) + 
     &      zab3(i1,i4,i5,i1,i6,i2,i3,i1))*
     &    (sqrtd1645*zab2(i2,i2,i3,i6) + 
     &      zab3(i2,i1,i6,i4,i5,i2,i3,i6) - 
     &      zab3(i2,i4,i5,i1,i6,i2,i3,i6))*
     &    (sqrtd1645*zab2(i3,i2,i3,i5) - 
     &      zab3(i3,i1,i6,i4,i5,i2,i3,i5) + 
     &      zab3(i3,i4,i5,i1,i6,i2,i3,i5))*
     &    (sqrtd1645*zab2(i4,i2,i3,i4) - 
     &      zab3(i4,i1,i6,i4,i5,i2,i3,i4) + 
     &      zab3(i4,i4,i5,i1,i6,i2,i3,i4)))/
     &  (16.*sqrtd1645**2*s(i2,i3)**2*za(i1,i6)*zaa2(i2,i1,i6,i4,i5,i2)*
     &    zaa2(i5,i1,i6,i4,i5,i5)*zbb2(i1,i4,i5,i1,i6,i1)*
     &    zbb2(i4,i4,i5,i1,i6,i4)) + 
     & (2*zab2(i1,i2,i3,i4)*zab2(i3,i2,i3,i4)**2*zab2(i4,i2,i3,i5)*
     &    zb(i6,i4))/
     &  (s(i2,i3)*za(i1,i5)*zab2(i2,i2,i3,i4)*zab2(i6,i2,i3,i4)*
     &    zb(i4,i1)*zbb2(i4,i4,i6,i1,i5,i4)) + 
     & ((sqrtd1546*za(i2,i3) + zaa2(i2,i1,i5,i4,i6,i3) - 
     &      zaa2(i2,i4,i6,i1,i5,i3))*
     &    (sqrtd1546*zab2(i1,i2,i3,i1) - 
     &      zab3(i1,i1,i5,i4,i6,i2,i3,i1) + 
     &      zab3(i1,i4,i6,i1,i5,i2,i3,i1))*
     &    (sqrtd1546*zab2(i3,i2,i3,i5) - 
     &      zab3(i3,i1,i5,i4,i6,i2,i3,i5) + 
     &      zab3(i3,i4,i6,i1,i5,i2,i3,i5))*
     &    (sqrtd1546*zab2(i4,i2,i3,i4) - 
     &      zab3(i4,i1,i5,i4,i6,i2,i3,i4) + 
     &      zab3(i4,i4,i6,i1,i5,i2,i3,i4))*
     &    (sqrtd1546*zab2(i6,i2,i3,i6) + 
     &      zab3(i6,i1,i5,i4,i6,i2,i3,i6) - 
     &      zab3(i6,i4,i6,i1,i5,i2,i3,i6)))/
     &  (16.*sqrtd1546**2*s(i2,i3)**2*za(i1,i5)*zaa2(i2,i1,i5,i4,i6,i2)*
     &    zaa2(i6,i1,i5,i4,i6,i6)*zbb2(i1,i4,i6,i1,i5,i1)*
     &    zbb2(i4,i4,i6,i1,i5,i4)) - 
     & ((sqrtd1546*za(i2,i3) - zaa2(i2,i1,i5,i4,i6,i3) + 
     &      zaa2(i2,i4,i6,i1,i5,i3))*
     &    (sqrtd1546*zab2(i1,i2,i3,i1) + 
     &      zab3(i1,i1,i5,i4,i6,i2,i3,i1) - 
     &      zab3(i1,i4,i6,i1,i5,i2,i3,i1))*
     &    (sqrtd1546*zab2(i3,i2,i3,i5) + 
     &      zab3(i3,i1,i5,i4,i6,i2,i3,i5) - 
     &      zab3(i3,i4,i6,i1,i5,i2,i3,i5))*
     &    (sqrtd1546*zab2(i4,i2,i3,i4) + 
     &      zab3(i4,i1,i5,i4,i6,i2,i3,i4) - 
     &      zab3(i4,i4,i6,i1,i5,i2,i3,i4))*
     &    (sqrtd1546*zab2(i6,i2,i3,i6) - 
     &      zab3(i6,i1,i5,i4,i6,i2,i3,i6) + 
     &      zab3(i6,i4,i6,i1,i5,i2,i3,i6)))/
     &  (16.*sqrtd1546**2*s(i2,i3)**2*za(i1,i5)*zaa2(i2,i1,i5,i4,i6,i2)*
     &    zaa2(i6,i1,i5,i4,i6,i6)*zbb2(i1,i4,i6,i1,i5,i1)*
     &    zbb2(i4,i4,i6,i1,i5,i4))
     &)


      aaaa_NMHV_b(24) = (
     &  (2*za(i3,i2)*zab2(i2,i5,i1,i5)**2*zab2(i2,i5,i1,i6)*
     &    zab2(i4,i5,i1,i2))/
     &  (s(i1,i5)*za(i6,i2)*zaa2(i2,i6,i4,i2,i3,i2)*zab2(i2,i5,i1,i1)*
     &    zab2(i2,i5,i1,i4)*zb(i2,i3)) - 
     & (2*za(i3,i6)*zab2(i4,i5,i1,i6)*zab2(i6,i5,i1,i2)*
     &    zab2(i6,i5,i1,i5)**2)/
     &  (s(i1,i5)*za(i6,i2)*zaa2(i6,i6,i4,i2,i3,i6)*zab2(i6,i5,i1,i1)*
     &    zab2(i6,i5,i1,i4)*zb(i2,i3)) + 
     & (2*za(i3,i2)*zab2(i2,i5,i1,i5)**2*zab2(i2,i5,i1,i6)*
     &    zab2(i4,i5,i1,i2))/
     &  (s(i1,i5)*za(i6,i2)*zaa2(i2,i6,i3,i2,i4,i2)*zab2(i2,i5,i1,i1)*
     &    zab2(i2,i5,i1,i3)*zb(i2,i4)) - 
     & (2*za(i3,i6)*zab2(i4,i5,i1,i6)*zab2(i6,i5,i1,i2)*
     &    zab2(i6,i5,i1,i5)**2)/
     &  (s(i1,i5)*za(i6,i2)*zaa2(i6,i6,i3,i2,i4,i6)*zab2(i6,i5,i1,i1)*
     &    zab2(i6,i5,i1,i3)*zb(i2,i4)) - 
     & (zab2(i1,i5,i1,i5)**2*zab2(i3,i6,i2,i1)*zab2(i4,i5,i1,i1)*
     &    zb(i2,i6))/
     &  (s(i1,i5)*t(i2,i6,i3)*za(i6,i2)*zab2(i1,i5,i1,i1)**2*zb(i1,i4)*
     &    zb(i2,i3)) + (zab2(i3,i6,i2,i4)*zab2(i4,i5,i1,i5)**2*
     &    zb(i2,i6))/
     &  (s(i1,i5)*t(i2,i6,i3)*za(i6,i2)*zab2(i4,i5,i1,i4)*zb(i1,i4)*
     &    zb(i2,i3)) - (zab2(i1,i5,i1,i5)**2*zab2(i3,i5,i1,i1)*
     &    zab2(i4,i6,i2,i1)*zb(i2,i6))/
     &  (s(i1,i5)*t(i2,i6,i4)*za(i6,i2)*zab2(i1,i5,i1,i1)**2*zb(i1,i3)*
     &    zb(i2,i4)) + (zab2(i3,i5,i1,i5)**2*zab2(i4,i6,i2,i3)*
     &    zb(i2,i6))/
     &  (s(i1,i5)*t(i2,i6,i4)*za(i6,i2)*zab2(i3,i5,i1,i3)*zb(i1,i3)*
     &    zb(i2,i4)) - (2*zaa2(i6,i1,i5,i2,i6,i6)*zab2(i2,i5,i1,i5)*
     &    zab2(i6,i3,i4,i2)*zab2(i6,i5,i1,i5)*zb(i2,i6))/
     &  (za(i6,i2)**2*zab2(i6,i5,i1,i1)*zab2(i6,i5,i1,i3)*
     &    zab2(i6,i5,i1,i4)*zab2(i6,i5,i1,i6)*zb(i2,i3)*zb(i2,i4)) - 
     & (2*s(i2,i6)*zab2(i6,i3,i4,i2)*zab2(i6,i5,i1,i5)**2*zb(i2,i6))/
     &  (za(i6,i2)*zab2(i6,i5,i1,i1)*zab2(i6,i5,i1,i3)*
     &    zab2(i6,i5,i1,i4)*zab2(i6,i5,i1,i6)*zb(i2,i3)*zb(i2,i4)) - 
     & (2*zab2(i6,i3,i4,i2)*zab2(i6,i5,i1,i5)**2*
     &    zab3(i6,i5,i1,i6,i2,i5,i1,i4)*zb(i2,i6))/
     &  (za(i6,i2)*zab2(i6,i5,i1,i1)*zab2(i6,i5,i1,i3)*
     &    zab2(i6,i5,i1,i4)**2*zab2(i6,i5,i1,i6)*zb(i2,i3)*zb(i2,i4)) + 
     & (2*za(i3,i1)*zab2(i4,i6,i2,i1)*zb(i1,i5)*zb(i2,i6)*zb(i5,i3))/
     &  (t(i2,i6,i4)*za(i6,i2)*zab2(i1,i5,i1,i1)*zb(i1,i3)**2*zb(i2,i4))
     &   + (2*za(i4,i1)*zab2(i3,i6,i2,i1)*zb(i1,i5)*zb(i2,i6)*
     &    zb(i5,i4))/
     &  (t(i2,i6,i3)*za(i6,i2)*zab2(i1,i5,i1,i1)*zb(i1,i4)**2*zb(i2,i3))
     &   + (2*s(i1,i5)*zaa2(i6,i1,i5,i2,i6,i6)*zab2(i6,i3,i4,i2)*
     &    zab2(i6,i5,i1,i5)*zb(i6,i5))/
     &  (za(i6,i2)*zab2(i6,i5,i1,i1)*zab2(i6,i5,i1,i3)**2*
     &    zab2(i6,i5,i1,i4)*zab2(i6,i5,i1,i6)*zb(i2,i4)) - 
     & (s(i1,i5)*zaa2(i6,i1,i5,i2,i6,i6)*zab2(i6,i3,i4,i2)*
     &    zab2(i6,i5,i1,i2)*zb(i6,i5)**2)/
     &  (za(i6,i2)*zab2(i6,i5,i1,i1)*zab2(i6,i5,i1,i3)*
     &    zab2(i6,i5,i1,i4)*zab2(i6,i5,i1,i6)**2*zb(i2,i3)*zb(i2,i4)) - 
     & (2*s(i1,i5)*za(i4,i1)*zab2(i3,i5,i1,i1)*zb(i1,i5)**2*zb(i1,i6)*
     &    zb(i2,i1))/
     &  (zab2(i1,i5,i1,i1)*zab2(i2,i5,i1,i1)*zab2(i6,i5,i1,i1)*
     &    zb(i1,i4)*zb(i2,i3)*zbb2(i1,i2,i3,i6,i4,i1)) - 
     & (2*s(i1,i5)*za(i4,i1)*zab2(i3,i5,i1,i1)*zb(i1,i5)**2*zb(i1,i6)*
     &    zb(i2,i1))/
     &  (zab2(i1,i5,i1,i1)*zab2(i2,i5,i1,i1)*zab2(i6,i5,i1,i1)*
     &    zb(i1,i3)*zb(i2,i4)*zbb2(i1,i2,i4,i6,i3,i1)) - 
     & (2*s(i2,i6)*zab2(i1,i5,i1,i2)*zb(i1,i5)**2*
     &    zbb2(i2,i4,i3,i1,i5,i1))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i1,i5,i1,i1)*zab2(i6,i5,i1,i1)*
     &    zb(i1,i3)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)) + 
     & (2*zab2(i1,i5,i1,i5)*zb(i1,i5)*zbb2(i1,i2,i6,i1,i5,i1)*
     &    zbb2(i2,i4,i3,i1,i5,i1))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i1,i5,i1,i1)*zab2(i6,i5,i1,i1)*
     &    zb(i1,i3)**2*zb(i1,i4)*zb(i2,i4)) + 
     & (2*zab2(i1,i5,i1,i2)*zab2(i2,i5,i1,i5)*zb(i1,i5)*
     &    zbb2(i1,i2,i6,i1,i5,i1)*zbb2(i2,i4,i3,i1,i5,i1))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i1,i5,i1,i1)*zab2(i2,i5,i1,i1)*
     &    zab2(i6,i5,i1,i1)*zb(i1,i3)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)) - 
     & (zab2(i1,i5,i1,i5)**2*zb(i2,i1)*zbb2(i1,i2,i6,i1,i5,i1)*
     &    zbb2(i2,i4,i3,i1,i5,i1))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i1,i5,i1,i1)**2*zab2(i6,i5,i1,i1)*
     &    zb(i1,i3)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)) - 
     & (2*zab2(i1,i5,i1,i2)*zb(i1,i5)**2*zbb2(i1,i2,i6,i1,i5,i4)*
     &    zbb2(i2,i4,i3,i1,i5,i1))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i1,i5,i1,i1)*zab2(i6,i5,i1,i1)*
     &    zb(i1,i3)*zb(i1,i4)**2*zb(i2,i3)*zb(i2,i4)) + 
     & (2*zab2(i3,i5,i1,i5)*zab2(i6,i6,i2,i3)*zb(i5,i3)*
     &    zbb2(i2,i4,i3,i1,i5,i3))/
     &  (za(i6,i2)*zab2(i3,i5,i1,i3)*zab2(i6,i5,i1,i3)**2*zb(i1,i3)*
     &    zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i3,i5,i1,i5)*zab3(i3,i5,i1,i6,i2,i5,i1,i3)*zb(i5,i3)*
     &    zbb2(i2,i4,i3,i1,i5,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i3,i5,i1,i3)**2*zab2(i6,i5,i1,i3)*
     &    zb(i1,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*s(i2,i6)*zab2(i3,i5,i1,i2)*zb(i5,i3)**2*
     &    zbb2(i2,i4,i3,i1,i5,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i3,i5,i1,i3)*zab2(i6,i5,i1,i3)*
     &    zb(i1,i3)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*s(i2,i6)*zab2(i4,i5,i1,i2)*zb(i5,i4)**2*
     &    zbb2(i2,i4,i3,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i4,i5,i1,i4)*zab2(i6,i5,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*zab2(i4,i5,i1,i2)*zab3(i4,i5,i1,i6,i2,i5,i1,i4)*zb(i5,i4)**2*
     &    zbb2(i2,i4,i3,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i4,i5,i1,i4)**2*zab2(i6,i5,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*s(i1,i5)*za(i3,i4)*zb(i2,i3)*zb(i5,i3)**2*zb(i6,i3))/
     &  (zab2(i2,i5,i1,i3)*zab2(i6,i5,i1,i3)*zb(i1,i3)*zb(i2,i4)*
     &    zbb2(i3,i2,i4,i6,i3,i3)) - 
     & ((sqrtd2463*zab2(i2,i5,i1,i2) - zab3(i2,i5,i1,i3,i6,i4,i2,i2) + 
     &      zab3(i2,i5,i1,i4,i2,i3,i6,i2))*
     &    (sqrtd2463*zab2(i3,i5,i1,i3) + 
     &      zab3(i3,i5,i1,i3,i6,i4,i2,i3) - 
     &      zab3(i3,i5,i1,i4,i2,i3,i6,i3))*
     &    (sqrtd2463*zab2(i4,i5,i1,i5) - 
     &      zab3(i4,i5,i1,i3,i6,i4,i2,i5) + 
     &      zab3(i4,i5,i1,i4,i2,i3,i6,i5))*
     &    (sqrtd2463*zab2(i6,i5,i1,i6) - 
     &      zab3(i6,i5,i1,i3,i6,i4,i2,i6) + 
     &      zab3(i6,i5,i1,i4,i2,i3,i6,i6))*
     &    (sqrtd2463*zb(i1,i5) + zbb2(i1,i2,i4,i6,i3,i5) - 
     &      zbb2(i1,i6,i3,i2,i4,i5)))/
     &  (16.*sqrtd2463**2*s(i1,i5)**2*zaa2(i2,i6,i3,i2,i4,i2)*
     &    zaa2(i6,i6,i3,i2,i4,i6)*zb(i2,i4)*zbb2(i1,i2,i4,i6,i3,i1)*
     &    zbb2(i3,i2,i4,i6,i3,i3)) + 
     & ((sqrtd2463*zab2(i2,i5,i1,i2) + zab3(i2,i5,i1,i3,i6,i4,i2,i2) - 
     &      zab3(i2,i5,i1,i4,i2,i3,i6,i2))*
     &    (sqrtd2463*zab2(i3,i5,i1,i3) - 
     &      zab3(i3,i5,i1,i3,i6,i4,i2,i3) + 
     &      zab3(i3,i5,i1,i4,i2,i3,i6,i3))*
     &    (sqrtd2463*zab2(i4,i5,i1,i5) + 
     &      zab3(i4,i5,i1,i3,i6,i4,i2,i5) - 
     &      zab3(i4,i5,i1,i4,i2,i3,i6,i5))*
     &    (sqrtd2463*zab2(i6,i5,i1,i6) + 
     &      zab3(i6,i5,i1,i3,i6,i4,i2,i6) - 
     &      zab3(i6,i5,i1,i4,i2,i3,i6,i6))*
     &    (sqrtd2463*zb(i1,i5) - zbb2(i1,i2,i4,i6,i3,i5) + 
     &      zbb2(i1,i6,i3,i2,i4,i5)))/
     &  (16.*sqrtd2463**2*s(i1,i5)**2*zaa2(i2,i6,i3,i2,i4,i2)*
     &    zaa2(i6,i6,i3,i2,i4,i6)*zb(i2,i4)*zbb2(i1,i2,i4,i6,i3,i1)*
     &    zbb2(i3,i2,i4,i6,i3,i3)) - 
     & (zab2(i3,i5,i1,i5)**2*zbb2(i2,i4,i3,i1,i5,i3)*
     &    zbb2(i3,i2,i6,i1,i5,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i3,i5,i1,i3)**2*zab2(i6,i5,i1,i3)*
     &    zb(i1,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*zab2(i3,i5,i1,i5)*zb(i1,i5)*zbb2(i2,i4,i3,i1,i5,i3)*
     &    zbb2(i3,i2,i6,i1,i5,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i3,i5,i1,i3)*zab2(i6,i5,i1,i3)*
     &    zb(i1,i3)**2*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i2,i5,i1,i5)*zab2(i3,i5,i1,i2)*zb(i5,i3)*
     &    zbb2(i2,i4,i3,i1,i5,i3)*zbb2(i3,i2,i6,i1,i5,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i2,i5,i1,i3)*zab2(i3,i5,i1,i3)*
     &    zab2(i6,i5,i1,i3)*zb(i1,i3)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i3,i5,i1,i5)*zb(i5,i3)*zbb2(i2,i4,i3,i1,i5,i4)*
     &    zbb2(i3,i2,i6,i1,i5,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i3,i5,i1,i3)*zab2(i6,i5,i1,i3)*
     &    zb(i1,i3)*zb(i2,i4)*zb(i4,i3)**2) - 
     & (2*zab2(i3,i5,i1,i2)*zb(i5,i3)**2*zbb2(i2,i4,i3,i1,i5,i3)*
     &    zbb2(i4,i1,i5,i2,i6,i3))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i3,i5,i1,i3)*zab2(i6,i5,i1,i3)*
     &    zb(i1,i3)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)**2) - 
     & (2*zab2(i4,i5,i1,i2)*zab2(i6,i3,i4,i2)*zb(i5,i4)**2*
     &    zbb2(i4,i2,i6,i1,i5,i4))/
     &  (za(i6,i2)*zab2(i4,i5,i1,i4)*zab2(i6,i5,i1,i4)**2*zb(i1,i4)*
     &    zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (zab2(i4,i5,i1,i5)**2*zbb2(i2,i4,i3,i1,i5,i4)*
     &    zbb2(i4,i2,i6,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i4,i5,i1,i4)**2*zab2(i6,i5,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i4,i3)) + 
     & (2*zab2(i4,i5,i1,i5)*zb(i5,i4)*zbb2(i2,i4,i3,i1,i5,i4)*
     &    zbb2(i4,i2,i6,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i4,i5,i1,i4)*zab2(i6,i5,i1,i4)*
     &    zb(i1,i4)*zb(i2,i4)*zb(i4,i3)**2) + 
     & (2*zab2(i2,i5,i1,i5)*zab2(i4,i5,i1,i2)*zb(i5,i4)*
     &    zbb2(i2,i4,i3,i1,i5,i4)*zbb2(i4,i2,i6,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i2,i5,i1,i4)*zab2(i4,i5,i1,i4)*
     &    zab2(i6,i5,i1,i4)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i4,i5,i1,i2)*zb(i1,i5)*zb(i5,i4)*zbb2(i2,i4,i3,i1,i5,i4)*
     &    zbb2(i4,i2,i6,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i4,i5,i1,i4)*zab2(i6,i5,i1,i4)*
     &    zb(i1,i4)**2*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i4,i5,i1,i2)*zb(i5,i3)*zb(i5,i4)*zbb2(i2,i4,i3,i1,i5,i4)*
     &    zbb2(i4,i2,i6,i1,i5,i4))/
     &  (s(i1,i5)*za(i6,i2)*zab2(i4,i5,i1,i4)*zab2(i6,i5,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)**2) - 
     & ((sqrtd2364*zab2(i2,i5,i1,i2) - zab3(i2,i5,i1,i3,i2,i4,i6,i2) + 
     &      zab3(i2,i5,i1,i4,i6,i3,i2,i2))*
     &    (sqrtd2364*zab2(i3,i5,i1,i1) + 
     &      zab3(i3,i5,i1,i3,i2,i4,i6,i1) - 
     &      zab3(i3,i5,i1,i4,i6,i3,i2,i1))*
     &    (sqrtd2364*zab2(i4,i5,i1,i5) - 
     &      zab3(i4,i5,i1,i3,i2,i4,i6,i5) + 
     &      zab3(i4,i5,i1,i4,i6,i3,i2,i5))*
     &    (sqrtd2364*zab2(i6,i5,i1,i6) - 
     &      zab3(i6,i5,i1,i3,i2,i4,i6,i6) + 
     &      zab3(i6,i5,i1,i4,i6,i3,i2,i6))*
     &    (sqrtd2364*zb(i5,i4) + zbb2(i5,i2,i3,i6,i4,i4) - 
     &      zbb2(i5,i6,i4,i2,i3,i4)))/
     &  (16.*sqrtd2364**2*s(i1,i5)**2*zaa2(i2,i6,i4,i2,i3,i2)*
     &    zaa2(i6,i6,i4,i2,i3,i6)*zb(i2,i3)*zbb2(i1,i2,i3,i6,i4,i1)*
     &    zbb2(i4,i2,i3,i6,i4,i4)) + 
     & ((sqrtd2364*zab2(i2,i5,i1,i2) + zab3(i2,i5,i1,i3,i2,i4,i6,i2) - 
     &      zab3(i2,i5,i1,i4,i6,i3,i2,i2))*
     &    (sqrtd2364*zab2(i3,i5,i1,i1) - 
     &      zab3(i3,i5,i1,i3,i2,i4,i6,i1) + 
     &      zab3(i3,i5,i1,i4,i6,i3,i2,i1))*
     &    (sqrtd2364*zab2(i4,i5,i1,i5) + 
     &      zab3(i4,i5,i1,i3,i2,i4,i6,i5) - 
     &      zab3(i4,i5,i1,i4,i6,i3,i2,i5))*
     &    (sqrtd2364*zab2(i6,i5,i1,i6) + 
     &      zab3(i6,i5,i1,i3,i2,i4,i6,i6) - 
     &      zab3(i6,i5,i1,i4,i6,i3,i2,i6))*
     &    (sqrtd2364*zb(i5,i4) - zbb2(i5,i2,i3,i6,i4,i4) + 
     &      zbb2(i5,i6,i4,i2,i3,i4)))/
     &  (16.*sqrtd2364**2*s(i1,i5)**2*zaa2(i2,i6,i4,i2,i3,i2)*
     &    zaa2(i6,i6,i4,i2,i3,i6)*zb(i2,i3)*zbb2(i1,i2,i3,i6,i4,i1)*
     &    zbb2(i4,i2,i3,i6,i4,i4))
     &)


      aaaa_NMHV_b(25) = (
     &  (2*za(i3,i2)*zab2(i2,i6,i1,i5)*zab2(i2,i6,i1,i6)**2*
     &    zab2(i4,i6,i1,i2))/
     &  (s(i1,i6)*za(i5,i2)*zaa2(i2,i5,i4,i2,i3,i2)*zab2(i2,i6,i1,i1)*
     &    zab2(i2,i6,i1,i4)*zb(i2,i3)) - 
     & (2*za(i3,i5)*zab2(i4,i6,i1,i5)*zab2(i5,i6,i1,i2)*
     &    zab2(i5,i6,i1,i6)**2)/
     &  (s(i1,i6)*za(i5,i2)*zaa2(i5,i5,i4,i2,i3,i5)*zab2(i5,i6,i1,i1)*
     &    zab2(i5,i6,i1,i4)*zb(i2,i3)) + 
     & (2*za(i3,i2)*zab2(i2,i6,i1,i5)*zab2(i2,i6,i1,i6)**2*
     &    zab2(i4,i6,i1,i2))/
     &  (s(i1,i6)*za(i5,i2)*zaa2(i2,i5,i3,i2,i4,i2)*zab2(i2,i6,i1,i1)*
     &    zab2(i2,i6,i1,i3)*zb(i2,i4)) - 
     & (2*za(i3,i5)*zab2(i4,i6,i1,i5)*zab2(i5,i6,i1,i2)*
     &    zab2(i5,i6,i1,i6)**2)/
     &  (s(i1,i6)*za(i5,i2)*zaa2(i5,i5,i3,i2,i4,i5)*zab2(i5,i6,i1,i1)*
     &    zab2(i5,i6,i1,i3)*zb(i2,i4)) - 
     & (zab2(i1,i6,i1,i6)**2*zab2(i3,i5,i2,i1)*zab2(i4,i6,i1,i1)*
     &    zb(i2,i5))/
     &  (s(i1,i6)*t(i2,i5,i3)*za(i5,i2)*zab2(i1,i6,i1,i1)**2*zb(i1,i4)*
     &    zb(i2,i3)) + (zab2(i3,i5,i2,i4)*zab2(i4,i6,i1,i6)**2*
     &    zb(i2,i5))/
     &  (s(i1,i6)*t(i2,i5,i3)*za(i5,i2)*zab2(i4,i6,i1,i4)*zb(i1,i4)*
     &    zb(i2,i3)) - (zab2(i1,i6,i1,i6)**2*zab2(i3,i6,i1,i1)*
     &    zab2(i4,i5,i2,i1)*zb(i2,i5))/
     &  (s(i1,i6)*t(i2,i5,i4)*za(i5,i2)*zab2(i1,i6,i1,i1)**2*zb(i1,i3)*
     &    zb(i2,i4)) + (zab2(i3,i6,i1,i6)**2*zab2(i4,i5,i2,i3)*
     &    zb(i2,i5))/
     &  (s(i1,i6)*t(i2,i5,i4)*za(i5,i2)*zab2(i3,i6,i1,i3)*zb(i1,i3)*
     &    zb(i2,i4)) - (2*zaa2(i5,i1,i6,i2,i5,i5)*zab2(i2,i6,i1,i6)*
     &    zab2(i5,i3,i4,i2)*zab2(i5,i6,i1,i6)*zb(i2,i5))/
     &  (za(i5,i2)**2*zab2(i5,i6,i1,i1)*zab2(i5,i6,i1,i3)*
     &    zab2(i5,i6,i1,i4)*zab2(i5,i6,i1,i5)*zb(i2,i3)*zb(i2,i4)) - 
     & (2*s(i2,i5)*zab2(i5,i3,i4,i2)*zab2(i5,i6,i1,i6)**2*zb(i2,i5))/
     &  (za(i5,i2)*zab2(i5,i6,i1,i1)*zab2(i5,i6,i1,i3)*
     &    zab2(i5,i6,i1,i4)*zab2(i5,i6,i1,i5)*zb(i2,i3)*zb(i2,i4)) - 
     & (2*zab2(i5,i3,i4,i2)*zab2(i5,i6,i1,i6)**2*
     &    zab3(i5,i6,i1,i5,i2,i6,i1,i4)*zb(i2,i5))/
     &  (za(i5,i2)*zab2(i5,i6,i1,i1)*zab2(i5,i6,i1,i3)*
     &    zab2(i5,i6,i1,i4)**2*zab2(i5,i6,i1,i5)*zb(i2,i3)*zb(i2,i4)) + 
     & (2*s(i1,i6)*zaa2(i5,i1,i6,i2,i5,i5)*zab2(i5,i3,i4,i2)*
     &    zab2(i5,i6,i1,i6)*zb(i5,i6))/
     &  (za(i5,i2)*zab2(i5,i6,i1,i1)*zab2(i5,i6,i1,i3)**2*
     &    zab2(i5,i6,i1,i4)*zab2(i5,i6,i1,i5)*zb(i2,i4)) - 
     & (s(i1,i6)*zaa2(i5,i1,i6,i2,i5,i5)*zab2(i5,i3,i4,i2)*
     &    zab2(i5,i6,i1,i2)*zb(i5,i6)**2)/
     &  (za(i5,i2)*zab2(i5,i6,i1,i1)*zab2(i5,i6,i1,i3)*
     &    zab2(i5,i6,i1,i4)*zab2(i5,i6,i1,i5)**2*zb(i2,i3)*zb(i2,i4)) + 
     & (2*za(i3,i1)*zab2(i4,i5,i2,i1)*zb(i1,i6)*zb(i2,i5)*zb(i6,i3))/
     &  (t(i2,i5,i4)*za(i5,i2)*zab2(i1,i6,i1,i1)*zb(i1,i3)**2*zb(i2,i4))
     &   + (2*za(i4,i1)*zab2(i3,i5,i2,i1)*zb(i1,i6)*zb(i2,i5)*
     &    zb(i6,i4))/
     &  (t(i2,i5,i3)*za(i5,i2)*zab2(i1,i6,i1,i1)*zb(i1,i4)**2*zb(i2,i3))
     &   - (2*s(i1,i6)*za(i4,i1)*zab2(i3,i6,i1,i1)*zb(i1,i5)*
     &    zb(i1,i6)**2*zb(i2,i1))/
     &  (zab2(i1,i6,i1,i1)*zab2(i2,i6,i1,i1)*zab2(i5,i6,i1,i1)*
     &    zb(i1,i4)*zb(i2,i3)*zbb2(i1,i2,i3,i5,i4,i1)) - 
     & (2*s(i1,i6)*za(i4,i1)*zab2(i3,i6,i1,i1)*zb(i1,i5)*zb(i1,i6)**2*
     &    zb(i2,i1))/
     &  (zab2(i1,i6,i1,i1)*zab2(i2,i6,i1,i1)*zab2(i5,i6,i1,i1)*
     &    zb(i1,i3)*zb(i2,i4)*zbb2(i1,i2,i4,i5,i3,i1)) - 
     & (2*s(i2,i5)*zab2(i1,i6,i1,i2)*zb(i1,i6)**2*
     &    zbb2(i2,i4,i3,i1,i6,i1))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i1,i6,i1,i1)*zab2(i5,i6,i1,i1)*
     &    zb(i1,i3)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)) + 
     & (2*zab2(i1,i6,i1,i6)*zb(i1,i6)*zbb2(i1,i2,i5,i1,i6,i1)*
     &    zbb2(i2,i4,i3,i1,i6,i1))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i1,i6,i1,i1)*zab2(i5,i6,i1,i1)*
     &    zb(i1,i3)**2*zb(i1,i4)*zb(i2,i4)) + 
     & (2*zab2(i1,i6,i1,i2)*zab2(i2,i6,i1,i6)*zb(i1,i6)*
     &    zbb2(i1,i2,i5,i1,i6,i1)*zbb2(i2,i4,i3,i1,i6,i1))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i1,i6,i1,i1)*zab2(i2,i6,i1,i1)*
     &    zab2(i5,i6,i1,i1)*zb(i1,i3)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)) - 
     & (zab2(i1,i6,i1,i6)**2*zb(i2,i1)*zbb2(i1,i2,i5,i1,i6,i1)*
     &    zbb2(i2,i4,i3,i1,i6,i1))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i1,i6,i1,i1)**2*zab2(i5,i6,i1,i1)*
     &    zb(i1,i3)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)) - 
     & (2*zab2(i1,i6,i1,i2)*zb(i1,i6)**2*zbb2(i1,i2,i5,i1,i6,i4)*
     &    zbb2(i2,i4,i3,i1,i6,i1))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i1,i6,i1,i1)*zab2(i5,i6,i1,i1)*
     &    zb(i1,i3)*zb(i1,i4)**2*zb(i2,i3)*zb(i2,i4)) + 
     & (2*zab2(i3,i6,i1,i6)*zab2(i5,i5,i2,i3)*zb(i6,i3)*
     &    zbb2(i2,i4,i3,i1,i6,i3))/
     &  (za(i5,i2)*zab2(i3,i6,i1,i3)*zab2(i5,i6,i1,i3)**2*zb(i1,i3)*
     &    zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i3,i6,i1,i6)*zab3(i3,i6,i1,i5,i2,i6,i1,i3)*zb(i6,i3)*
     &    zbb2(i2,i4,i3,i1,i6,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i3,i6,i1,i3)**2*zab2(i5,i6,i1,i3)*
     &    zb(i1,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*s(i2,i5)*zab2(i3,i6,i1,i2)*zb(i6,i3)**2*
     &    zbb2(i2,i4,i3,i1,i6,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i3,i6,i1,i3)*zab2(i5,i6,i1,i3)*
     &    zb(i1,i3)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*s(i2,i5)*zab2(i4,i6,i1,i2)*zb(i6,i4)**2*
     &    zbb2(i2,i4,i3,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i4,i6,i1,i4)*zab2(i5,i6,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*zab2(i4,i6,i1,i2)*zab3(i4,i6,i1,i5,i2,i6,i1,i4)*zb(i6,i4)**2*
     &    zbb2(i2,i4,i3,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i4,i6,i1,i4)**2*zab2(i5,i6,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*s(i1,i6)*za(i3,i4)*zb(i2,i3)*zb(i5,i3)*zb(i6,i3)**2)/
     &  (zab2(i2,i6,i1,i3)*zab2(i5,i6,i1,i3)*zb(i1,i3)*zb(i2,i4)*
     &    zbb2(i3,i2,i4,i5,i3,i3)) - 
     & ((sqrtd2453*zab2(i2,i6,i1,i2) - zab3(i2,i6,i1,i3,i5,i4,i2,i2) + 
     &      zab3(i2,i6,i1,i4,i2,i3,i5,i2))*
     &    (sqrtd2453*zab2(i3,i6,i1,i3) + 
     &      zab3(i3,i6,i1,i3,i5,i4,i2,i3) - 
     &      zab3(i3,i6,i1,i4,i2,i3,i5,i3))*
     &    (sqrtd2453*zab2(i4,i6,i1,i6) - 
     &      zab3(i4,i6,i1,i3,i5,i4,i2,i6) + 
     &      zab3(i4,i6,i1,i4,i2,i3,i5,i6))*
     &    (sqrtd2453*zab2(i5,i6,i1,i5) - 
     &      zab3(i5,i6,i1,i3,i5,i4,i2,i5) + 
     &      zab3(i5,i6,i1,i4,i2,i3,i5,i5))*
     &    (sqrtd2453*zb(i1,i6) + zbb2(i1,i2,i4,i5,i3,i6) - 
     &      zbb2(i1,i5,i3,i2,i4,i6)))/
     &  (16.*sqrtd2453**2*s(i1,i6)**2*zaa2(i2,i5,i3,i2,i4,i2)*
     &    zaa2(i5,i5,i3,i2,i4,i5)*zb(i2,i4)*zbb2(i1,i2,i4,i5,i3,i1)*
     &    zbb2(i3,i2,i4,i5,i3,i3)) + 
     & ((sqrtd2453*zab2(i2,i6,i1,i2) + zab3(i2,i6,i1,i3,i5,i4,i2,i2) - 
     &      zab3(i2,i6,i1,i4,i2,i3,i5,i2))*
     &    (sqrtd2453*zab2(i3,i6,i1,i3) - 
     &      zab3(i3,i6,i1,i3,i5,i4,i2,i3) + 
     &      zab3(i3,i6,i1,i4,i2,i3,i5,i3))*
     &    (sqrtd2453*zab2(i4,i6,i1,i6) + 
     &      zab3(i4,i6,i1,i3,i5,i4,i2,i6) - 
     &      zab3(i4,i6,i1,i4,i2,i3,i5,i6))*
     &    (sqrtd2453*zab2(i5,i6,i1,i5) + 
     &      zab3(i5,i6,i1,i3,i5,i4,i2,i5) - 
     &      zab3(i5,i6,i1,i4,i2,i3,i5,i5))*
     &    (sqrtd2453*zb(i1,i6) - zbb2(i1,i2,i4,i5,i3,i6) + 
     &      zbb2(i1,i5,i3,i2,i4,i6)))/
     &  (16.*sqrtd2453**2*s(i1,i6)**2*zaa2(i2,i5,i3,i2,i4,i2)*
     &    zaa2(i5,i5,i3,i2,i4,i5)*zb(i2,i4)*zbb2(i1,i2,i4,i5,i3,i1)*
     &    zbb2(i3,i2,i4,i5,i3,i3)) - 
     & (zab2(i3,i6,i1,i6)**2*zbb2(i2,i4,i3,i1,i6,i3)*
     &    zbb2(i3,i2,i5,i1,i6,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i3,i6,i1,i3)**2*zab2(i5,i6,i1,i3)*
     &    zb(i1,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (2*zab2(i3,i6,i1,i6)*zb(i1,i6)*zbb2(i2,i4,i3,i1,i6,i3)*
     &    zbb2(i3,i2,i5,i1,i6,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i3,i6,i1,i3)*zab2(i5,i6,i1,i3)*
     &    zb(i1,i3)**2*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i2,i6,i1,i6)*zab2(i3,i6,i1,i2)*zb(i6,i3)*
     &    zbb2(i2,i4,i3,i1,i6,i3)*zbb2(i3,i2,i5,i1,i6,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i2,i6,i1,i3)*zab2(i3,i6,i1,i3)*
     &    zab2(i5,i6,i1,i3)*zb(i1,i3)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i3,i6,i1,i6)*zb(i6,i3)*zbb2(i2,i4,i3,i1,i6,i4)*
     &    zbb2(i3,i2,i5,i1,i6,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i3,i6,i1,i3)*zab2(i5,i6,i1,i3)*
     &    zb(i1,i3)*zb(i2,i4)*zb(i4,i3)**2) - 
     & (2*zab2(i3,i6,i1,i2)*zb(i6,i3)**2*zbb2(i2,i4,i3,i1,i6,i3)*
     &    zbb2(i4,i1,i6,i2,i5,i3))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i3,i6,i1,i3)*zab2(i5,i6,i1,i3)*
     &    zb(i1,i3)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)**2) - 
     & (2*zab2(i4,i6,i1,i2)*zab2(i5,i3,i4,i2)*zb(i6,i4)**2*
     &    zbb2(i4,i2,i5,i1,i6,i4))/
     &  (za(i5,i2)*zab2(i4,i6,i1,i4)*zab2(i5,i6,i1,i4)**2*zb(i1,i4)*
     &    zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) + 
     & (zab2(i4,i6,i1,i6)**2*zbb2(i2,i4,i3,i1,i6,i4)*
     &    zbb2(i4,i2,i5,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i4,i6,i1,i4)**2*zab2(i5,i6,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i4,i3)) + 
     & (2*zab2(i4,i6,i1,i6)*zb(i6,i4)*zbb2(i2,i4,i3,i1,i6,i4)*
     &    zbb2(i4,i2,i5,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i4,i6,i1,i4)*zab2(i5,i6,i1,i4)*
     &    zb(i1,i4)*zb(i2,i4)*zb(i4,i3)**2) + 
     & (2*zab2(i2,i6,i1,i6)*zab2(i4,i6,i1,i2)*zb(i6,i4)*
     &    zbb2(i2,i4,i3,i1,i6,i4)*zbb2(i4,i2,i5,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i2,i6,i1,i4)*zab2(i4,i6,i1,i4)*
     &    zab2(i5,i6,i1,i4)*zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i4,i6,i1,i2)*zb(i1,i6)*zb(i6,i4)*zbb2(i2,i4,i3,i1,i6,i4)*
     &    zbb2(i4,i2,i5,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i4,i6,i1,i4)*zab2(i5,i6,i1,i4)*
     &    zb(i1,i4)**2*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)) - 
     & (2*zab2(i4,i6,i1,i2)*zb(i6,i3)*zb(i6,i4)*zbb2(i2,i4,i3,i1,i6,i4)*
     &    zbb2(i4,i2,i5,i1,i6,i4))/
     &  (s(i1,i6)*za(i5,i2)*zab2(i4,i6,i1,i4)*zab2(i5,i6,i1,i4)*
     &    zb(i1,i4)*zb(i2,i3)*zb(i2,i4)*zb(i4,i3)**2) - 
     & ((sqrtd2354*zab2(i2,i6,i1,i2) - zab3(i2,i6,i1,i3,i2,i4,i5,i2) + 
     &      zab3(i2,i6,i1,i4,i5,i3,i2,i2))*
     &    (sqrtd2354*zab2(i3,i6,i1,i1) + 
     &      zab3(i3,i6,i1,i3,i2,i4,i5,i1) - 
     &      zab3(i3,i6,i1,i4,i5,i3,i2,i1))*
     &    (sqrtd2354*zab2(i4,i6,i1,i6) - 
     &      zab3(i4,i6,i1,i3,i2,i4,i5,i6) + 
     &      zab3(i4,i6,i1,i4,i5,i3,i2,i6))*
     &    (sqrtd2354*zab2(i5,i6,i1,i5) - 
     &      zab3(i5,i6,i1,i3,i2,i4,i5,i5) + 
     &      zab3(i5,i6,i1,i4,i5,i3,i2,i5))*
     &    (sqrtd2354*zb(i6,i4) + zbb2(i6,i2,i3,i5,i4,i4) - 
     &      zbb2(i6,i5,i4,i2,i3,i4)))/
     &  (16.*sqrtd2354**2*s(i1,i6)**2*zaa2(i2,i5,i4,i2,i3,i2)*
     &    zaa2(i5,i5,i4,i2,i3,i5)*zb(i2,i3)*zbb2(i1,i2,i3,i5,i4,i1)*
     &    zbb2(i4,i2,i3,i5,i4,i4)) + 
     & ((sqrtd2354*zab2(i2,i6,i1,i2) + zab3(i2,i6,i1,i3,i2,i4,i5,i2) - 
     &      zab3(i2,i6,i1,i4,i5,i3,i2,i2))*
     &    (sqrtd2354*zab2(i3,i6,i1,i1) - 
     &      zab3(i3,i6,i1,i3,i2,i4,i5,i1) + 
     &      zab3(i3,i6,i1,i4,i5,i3,i2,i1))*
     &    (sqrtd2354*zab2(i4,i6,i1,i6) + 
     &      zab3(i4,i6,i1,i3,i2,i4,i5,i6) - 
     &      zab3(i4,i6,i1,i4,i5,i3,i2,i6))*
     &    (sqrtd2354*zab2(i5,i6,i1,i5) + 
     &      zab3(i5,i6,i1,i3,i2,i4,i5,i5) - 
     &      zab3(i5,i6,i1,i4,i5,i3,i2,i5))*
     &    (sqrtd2354*zb(i6,i4) - zbb2(i6,i2,i3,i5,i4,i4) + 
     &      zbb2(i6,i5,i4,i2,i3,i4)))/
     &  (16.*sqrtd2354**2*s(i1,i6)**2*zaa2(i2,i5,i4,i2,i3,i2)*
     &    zaa2(i5,i5,i4,i2,i3,i5)*zb(i2,i3)*zbb2(i1,i2,i3,i5,i4,i1)*
     &    zbb2(i4,i2,i3,i5,i4,i4))
     &)


      return
      end
