!===== T. Dennen, May 2014
!===== Box coefficients for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,-)ga(i5,+)ga(i6,+)
      subroutine aaaa_NMHV_d_init(i1,i2,i3,i4,i5,i6,za,zb,aaaa_NMHV_d)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: aaaa_NMHV_d(195), zab2, t

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      aaaa_NMHV_d(:) = czip

      aaaa_NMHV_d(1) = (
     &  (s(i1,i2)*t(i4,i5,i6)*za(i2,i3)*zab2(i4,i5,i6,i3)**2*
     &   zb(i2,i1)**2)/
     & (zab2(i5,i4,i6,i1)*zab2(i5,i4,i6,i3)*zab2(i6,i4,i5,i1)*
     &   zab2(i6,i4,i5,i3)*zb(i3,i1))
     & )

      aaaa_NMHV_d(2) = (
     &  -((s(i1,i2)*t(i3,i4,i5)*za(i1,i2)**2*zab2(i6,i3,i4,i5)**2*
     &     zb(i6,i1))/
     &   (za(i2,i6)*zab2(i2,i3,i5,i4)*zab2(i2,i4,i5,i3)*
     &     zab2(i6,i3,i5,i4)*zab2(i6,i4,i5,i3)))
     & )

      aaaa_NMHV_d(5) = (
     &  -((s(i1,i2)*t(i3,i4,i6)*za(i1,i2)**2*zab2(i5,i3,i4,i6)**2*
     &     zb(i5,i1))/
     &   (za(i2,i5)*zab2(i2,i3,i6,i4)*zab2(i2,i4,i6,i3)*
     &     zab2(i5,i3,i6,i4)*zab2(i5,i4,i6,i3)))
     & )

      aaaa_NMHV_d(8) = (
     &  (za(i4,i6)*za(i6,i1)*zab2(i2,i3,i5,i1)*zab2(i3,i5,i2,i1)**2*
     &   zb(i1,i6)*zb(i6,i4)**3)/
     & (t(i2,i3,i5)*za(i5,i2)*zab2(i2,i3,i5,i4)*zab2(i5,i3,i2,i1)*
     &   zb(i1,i4)**3)
     & )

      aaaa_NMHV_d(9) = (
     &  (za(i4,i5)*za(i4,i6)**3*zab2(i5,i1,i3,i2)**2*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i5,i6)**3*zab2(i5,i1,i2,i3)*zab2(i6,i1,i3,i2)*zb(i3,i1)) - 
     & (za(i4,i5)**3*za(i4,i6)*zab2(i6,i1,i3,i2)**2*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i5,i6)**3*zab2(i5,i1,i3,i2)*zab2(i6,i1,i2,i3)*zb(i3,i1)) - 
     & (za(i4,i5)*za(i4,i6)**3*zab2(i5,i1,i3,i2)**2*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i5,i6)**3*zab2(i5,i2,i3,i1)*zab2(i6,i1,i2,i3)*zb(i3,i2)) + 
     & (za(i4,i5)**3*za(i4,i6)*zab2(i6,i1,i3,i2)**2*zb(i5,i4)*
     &    zb(i6,i4))/
     &  (za(i5,i6)**3*zab2(i5,i1,i2,i3)*zab2(i6,i2,i3,i1)*zb(i3,i2))
     & )

      aaaa_NMHV_d(10) = (
     &  (s(i1,i2)*t(i3,i5,i6)*za(i1,i4)*zab2(i3,i6,i5,i2)**2*zb(i4,i2))/
     & (zab2(i5,i3,i6,i2)*zab2(i5,i3,i6,i4)*zab2(i6,i3,i5,i2)*
     &   zab2(i6,i3,i5,i4))
     & )

      aaaa_NMHV_d(11) = (
     &  -((za(i1,i4)**3*za(i4,i6)*zab2(i1,i3,i5,i2)*zb(i4,i1)*
     &     zb(i5,i2)**2*zb(i6,i4))/
     &   (t(i2,i3,i5)*za(i1,i6)*zab2(i1,i5,i2,i3)*zab2(i6,i3,i5,i2)*
     &     zb(i3,i2)))
     & )

      aaaa_NMHV_d(12) = (
     &  (za(i4,i5)*za(i5,i1)*zab2(i2,i3,i6,i1)*zab2(i3,i6,i2,i1)**2*
     &   zb(i1,i5)*zb(i5,i4)**3)/
     & (t(i2,i3,i6)*za(i6,i2)*zab2(i2,i3,i6,i4)*zab2(i6,i3,i2,i1)*
     &   zb(i1,i4)**3)
     & )

      aaaa_NMHV_d(13) = (
     &  -((za(i1,i4)**3*za(i4,i5)*zab2(i1,i3,i6,i2)*zb(i4,i1)*zb(i5,i4)*
     &     zb(i6,i2)**2)/
     &   (t(i2,i3,i6)*za(i1,i5)*zab2(i1,i6,i2,i3)*zab2(i5,i3,i6,i2)*
     &     zb(i3,i2)))
     & )

      aaaa_NMHV_d(14) = (
     &  (s(i1,i2)*t(i3,i5,i6)*za(i2,i4)*zab2(i3,i5,i6,i4)**2*
     &   zb(i2,i1)**2)/
     & (zab2(i5,i3,i6,i1)*zab2(i5,i3,i6,i4)*zab2(i6,i3,i5,i1)*
     &   zab2(i6,i3,i5,i4)*zb(i4,i1))
     & )

      aaaa_NMHV_d(17) = (
     &  (za(i3,i6)*za(i6,i1)*zab2(i2,i4,i5,i1)*zab2(i4,i5,i2,i1)**2*
     &   zb(i1,i6)*zb(i6,i3)**3)/
     & (t(i2,i4,i5)*za(i5,i2)*zab2(i2,i4,i5,i3)*zab2(i5,i4,i2,i1)*
     &   zb(i1,i3)**3)
     & )

      aaaa_NMHV_d(18) = (
     &  (za(i3,i5)*za(i3,i6)**3*zab2(i5,i1,i4,i2)**2*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i5,i6)**3*zab2(i5,i1,i2,i4)*zab2(i6,i1,i4,i2)*zb(i4,i1)) - 
     & (za(i3,i5)**3*za(i3,i6)*zab2(i6,i1,i4,i2)**2*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i5,i6)**3*zab2(i5,i1,i4,i2)*zab2(i6,i1,i2,i4)*zb(i4,i1)) - 
     & (za(i3,i5)*za(i3,i6)**3*zab2(i5,i1,i4,i2)**2*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i5,i6)**3*zab2(i5,i2,i4,i1)*zab2(i6,i1,i2,i4)*zb(i4,i2)) + 
     & (za(i3,i5)**3*za(i3,i6)*zab2(i6,i1,i4,i2)**2*zb(i5,i3)*
     &    zb(i6,i3))/
     &  (za(i5,i6)**3*zab2(i5,i1,i2,i4)*zab2(i6,i2,i4,i1)*zb(i4,i2))
     & )

      aaaa_NMHV_d(19) = (
     &  (s(i1,i2)*t(i4,i5,i6)*za(i1,i3)*zab2(i4,i6,i5,i2)**2*zb(i3,i2))/
     & (zab2(i5,i4,i6,i2)*zab2(i5,i4,i6,i3)*zab2(i6,i4,i5,i2)*
     &   zab2(i6,i4,i5,i3))
     & )

      aaaa_NMHV_d(20) = (
     &  -((za(i1,i3)**3*za(i3,i6)*zab2(i1,i4,i5,i2)*zb(i3,i1)*
     &     zb(i5,i2)**2*zb(i6,i3))/
     &   (t(i2,i4,i5)*za(i1,i6)*zab2(i1,i5,i2,i4)*zab2(i6,i4,i5,i2)*
     &     zb(i4,i2)))
     & )

      aaaa_NMHV_d(21) = (
     &  (za(i3,i5)*za(i5,i1)*zab2(i2,i4,i6,i1)*zab2(i4,i6,i2,i1)**2*
     &   zb(i1,i5)*zb(i5,i3)**3)/
     & (t(i2,i4,i6)*za(i6,i2)*zab2(i2,i4,i6,i3)*zab2(i6,i4,i2,i1)*
     &   zb(i1,i3)**3)
     & )

      aaaa_NMHV_d(22) = (
     &  -((za(i1,i3)**3*za(i3,i5)*zab2(i1,i4,i6,i2)*zb(i3,i1)*zb(i5,i3)*
     &     zb(i6,i2)**2)/
     &   (t(i2,i4,i6)*za(i1,i5)*zab2(i1,i6,i2,i4)*zab2(i5,i4,i6,i2)*
     &     zb(i4,i2)))
     & )

      aaaa_NMHV_d(23) = (
     &  -((s(i1,i2)*t(i3,i4,i6)*za(i1,i5)*zab2(i1,i3,i4,i6)**2*
     &     zb(i5,i2))/
     &   (zab2(i1,i3,i6,i4)*zab2(i1,i4,i6,i3)*zab2(i5,i3,i6,i4)*
     &     zab2(i5,i4,i6,i3)))
     & )

      aaaa_NMHV_d(25) = (
     &  -((za(i3,i6)*za(i4,i6)*zab2(i1,i5,i2,i4)**2*zb(i6,i3)**3*
     &      zb(i6,i4))/
     &    (za(i1,i5)*zab2(i2,i5,i1,i4)*zab2(i5,i1,i2,i3)*zb(i4,i3)**3))
     &  + (za(i3,i6)*za(i4,i6)*zab2(i1,i5,i2,i4)**2*zb(i6,i3)**3*
     &    zb(i6,i4))/
     &  (za(i2,i5)*zab2(i1,i5,i2,i3)*zab2(i5,i1,i2,i4)*zb(i4,i3)**3) - 
     & (za(i3,i6)*za(i4,i6)*zab2(i1,i5,i2,i3)**2*zb(i6,i3)*
     &    zb(i6,i4)**3)/
     &  (za(i2,i5)*zab2(i1,i5,i2,i4)*zab2(i5,i1,i2,i3)*zb(i4,i3)**3) + 
     & (za(i3,i6)*za(i4,i6)*zab2(i1,i5,i2,i3)**2*zb(i6,i3)*
     &    zb(i6,i4)**3)/
     &  (za(i1,i5)*zab2(i2,i5,i1,i3)*zab2(i5,i1,i2,i4)*zb(i4,i3)**3)
     & )

      aaaa_NMHV_d(29) = (
     &  -((s(i1,i2)*t(i3,i4,i5)*za(i1,i6)*zab2(i1,i3,i4,i5)**2*
     &     zb(i6,i2))/
     &   (zab2(i1,i3,i5,i4)*zab2(i1,i4,i5,i3)*zab2(i6,i3,i5,i4)*
     &     zab2(i6,i4,i5,i3)))
     & )

      aaaa_NMHV_d(31) = (
     &  -((za(i3,i5)*za(i4,i5)*zab2(i1,i6,i2,i4)**2*zb(i5,i3)**3*
     &      zb(i5,i4))/
     &    (za(i1,i6)*zab2(i2,i6,i1,i4)*zab2(i6,i1,i2,i3)*zb(i4,i3)**3))
     &  + (za(i3,i5)*za(i4,i5)*zab2(i1,i6,i2,i4)**2*zb(i5,i3)**3*
     &    zb(i5,i4))/
     &  (za(i2,i6)*zab2(i1,i6,i2,i3)*zab2(i6,i1,i2,i4)*zb(i4,i3)**3) - 
     & (za(i3,i5)*za(i4,i5)*zab2(i1,i6,i2,i3)**2*zb(i5,i3)*
     &    zb(i5,i4)**3)/
     &  (za(i2,i6)*zab2(i1,i6,i2,i4)*zab2(i6,i1,i2,i3)*zb(i4,i3)**3) + 
     & (za(i3,i5)*za(i4,i5)*zab2(i1,i6,i2,i3)**2*zb(i5,i3)*
     &    zb(i5,i4)**3)/
     &  (za(i1,i6)*zab2(i2,i6,i1,i3)*zab2(i6,i1,i2,i4)*zb(i4,i3)**3)
     & )

      aaaa_NMHV_d(42) = (
     &  -((za(i2,i4)*za(i4,i6)**3*zab2(i2,i1,i3,i5)**2*
     &     zab2(i2,i3,i5,i1)*zb(i4,i2)*zb(i6,i4))/
     &   (t(i1,i3,i5)*za(i2,i6)**3*zab2(i2,i1,i5,i3)*zab2(i6,i3,i5,i1)*
     &     zb(i3,i1)))
     & )

      aaaa_NMHV_d(43) = (
     &  (za(i1,i3)**2*za(i2,i6)*za(i4,i6)*zab2(i1,i3,i5,i2)*
     &   zb(i6,i2)**3*zb(i6,i4))/
     & (t(i1,i3,i5)*za(i1,i5)*zab2(i1,i3,i5,i4)*zab2(i5,i1,i3,i2)*
     &   zb(i4,i2))
     & )

      aaaa_NMHV_d(46) = (
     &  -((za(i2,i4)*za(i4,i5)**3*zab2(i2,i1,i3,i6)**2*
     &     zab2(i2,i3,i6,i1)*zb(i4,i2)*zb(i5,i4))/
     &   (t(i1,i3,i6)*za(i2,i5)**3*zab2(i2,i1,i6,i3)*zab2(i5,i3,i6,i1)*
     &     zb(i3,i1)))
     & )

      aaaa_NMHV_d(47) = (
     &  (za(i1,i3)**2*za(i2,i5)*za(i4,i5)*zab2(i1,i3,i6,i2)*
     &   zb(i5,i2)**3*zb(i5,i4))/
     & (t(i1,i3,i6)*za(i1,i6)*zab2(i1,i3,i6,i4)*zab2(i6,i1,i3,i2)*
     &   zb(i4,i2))
     & )

      aaaa_NMHV_d(51) = (
     &  -((za(i2,i3)*za(i3,i6)**3*zab2(i2,i1,i4,i5)**2*
     &     zab2(i2,i4,i5,i1)*zb(i3,i2)*zb(i6,i3))/
     &   (t(i1,i4,i5)*za(i2,i6)**3*zab2(i2,i1,i5,i4)*zab2(i6,i4,i5,i1)*
     &     zb(i4,i1)))
     & )

      aaaa_NMHV_d(52) = (
     &  (za(i1,i4)**2*za(i2,i6)*za(i3,i6)*zab2(i1,i4,i5,i2)*
     &   zb(i6,i2)**3*zb(i6,i3))/
     & (t(i1,i4,i5)*za(i1,i5)*zab2(i1,i4,i5,i3)*zab2(i5,i1,i4,i2)*
     &   zb(i3,i2))
     & )

      aaaa_NMHV_d(54) = (
     &  -((za(i2,i3)*za(i3,i5)**3*zab2(i2,i1,i4,i6)**2*
     &     zab2(i2,i4,i6,i1)*zb(i3,i2)*zb(i5,i3))/
     &   (t(i1,i4,i6)*za(i2,i5)**3*zab2(i2,i1,i6,i4)*zab2(i5,i4,i6,i1)*
     &     zb(i4,i1)))
     & )

      aaaa_NMHV_d(55) = (
     &  (za(i1,i4)**2*za(i2,i5)*za(i3,i5)*zab2(i1,i4,i6,i2)*
     &   zb(i5,i2)**3*zb(i5,i3))/
     & (t(i1,i4,i6)*za(i1,i6)*zab2(i1,i4,i6,i3)*zab2(i6,i1,i4,i2)*
     &   zb(i3,i2))
     & )

      aaaa_NMHV_d(107) = (
     &  -((s(i1,i6)*t(i1,i2,i3)**3*za(i4,i6)**2*zb(i2,i1)**2)/
     &   (za(i5,i6)*zab2(i5,i3,i2,i1)*zab2(i6,i1,i2,i3)*
     &     zab2(i6,i3,i2,i1)**2*zb(i3,i2)))
     & )

      aaaa_NMHV_d(110) = (
     &  -((s(i1,i5)*t(i1,i2,i3)**3*za(i4,i5)**2*zb(i2,i1)**2)/
     &   (za(i6,i5)*zab2(i5,i1,i2,i3)*zab2(i5,i3,i2,i1)**2*
     &     zab2(i6,i3,i2,i1)*zb(i3,i2)))
     & )

      aaaa_NMHV_d(112) = (
     &  (s(i1,i2)*zab2(i2,i1,i4,i6)**2*zab2(i3,i4,i6,i1)**2)/
     &  (za(i2,i5)*zab2(i2,i3,i5,i4)*zab2(i2,i4,i6,i1)*
     &    zab2(i5,i4,i6,i1)*zb(i4,i1)) + 
     & (s(i1,i2)*za(i1,i4)**2*zab2(i1,i4,i6,i2)*zb(i5,i2)**2)/
     &  (za(i1,i6)*zab2(i1,i4,i6,i3)*zab2(i6,i1,i4,i2)*zb(i3,i2))
     & )

      aaaa_NMHV_d(113) = (
     &  (s(i3,i5)*t(i1,i2,i3)*
     &   (zab2(i4,i1,i2,i3)**2*zab2(i5,i1,i3,i2)**3*zb(i3,i1) + 
     &     t(i1,i2,i3)**2*za(i4,i5)**2*zab2(i5,i2,i3,i1)*zb(i3,i2)**3))/
     & (za(i6,i5)*zab2(i5,i1,i2,i3)**2*zab2(i5,i1,i3,i2)*
     &   zab2(i5,i2,i3,i1)*zab2(i6,i1,i2,i3)*zb(i3,i1)*zb(i3,i2))
     & )

      aaaa_NMHV_d(114) = (
     &  (s(i4,i6)*t(i1,i2,i6)*
     &   (-(za(i2,i6)*zab2(i1,i6,i2,i4)**3*zab2(i6,i1,i2,i5)**2) - 
     &     t(i1,i2,i6)**2*za(i1,i6)**3*zab2(i2,i6,i1,i4)*zb(i5,i4)**2))/
     & (za(i1,i6)*za(i2,i6)*zab2(i1,i6,i2,i4)*zab2(i2,i6,i1,i4)*
     &   zab2(i6,i1,i2,i3)*zab2(i6,i1,i2,i4)**2*zb(i4,i3))
     & )

      aaaa_NMHV_d(116) = (
     &  (s(i4,i6)*t(i1,i2,i4)*
     &   (zab2(i3,i1,i2,i4)**2*zab2(i6,i1,i4,i2)**3*zb(i4,i1) + 
     &     t(i1,i2,i4)**2*za(i3,i6)**2*zab2(i6,i2,i4,i1)*zb(i4,i2)**3))/
     & (za(i5,i6)*zab2(i5,i1,i2,i4)*zab2(i6,i1,i2,i4)**2*
     &   zab2(i6,i1,i4,i2)*zab2(i6,i2,i4,i1)*zb(i4,i1)*zb(i4,i2))
     & )

      aaaa_NMHV_d(117) = (
     &  (s(i1,i2)*zab2(i2,i1,i4,i5)**2*zab2(i3,i4,i5,i1)**2)/
     &  (za(i2,i6)*zab2(i2,i3,i6,i4)*zab2(i2,i4,i5,i1)*
     &    zab2(i6,i4,i5,i1)*zb(i4,i1)) + 
     & (s(i1,i2)*za(i1,i4)**2*zab2(i1,i4,i5,i2)*zb(i6,i2)**2)/
     &  (za(i1,i5)*zab2(i1,i4,i5,i3)*zab2(i5,i1,i4,i2)*zb(i3,i2))
     & )

      aaaa_NMHV_d(118) = (
     &  (s(i3,i6)*t(i1,i2,i3)*
     &   (zab2(i4,i1,i2,i3)**2*zab2(i6,i1,i3,i2)**3*zb(i3,i1) + 
     &     t(i1,i2,i3)**2*za(i4,i6)**2*zab2(i6,i2,i3,i1)*zb(i3,i2)**3))/
     & (za(i5,i6)*zab2(i5,i1,i2,i3)*zab2(i6,i1,i2,i3)**2*
     &   zab2(i6,i1,i3,i2)*zab2(i6,i2,i3,i1)*zb(i3,i1)*zb(i3,i2))
     & )

      aaaa_NMHV_d(119) = (
     &  (s(i4,i5)*t(i1,i2,i5)*
     &   (-(za(i2,i5)*zab2(i1,i5,i2,i4)**3*zab2(i5,i1,i2,i6)**2) - 
     &     t(i1,i2,i5)**2*za(i1,i5)**3*zab2(i2,i5,i1,i4)*zb(i6,i4)**2))/
     & (za(i1,i5)*za(i2,i5)*zab2(i1,i5,i2,i4)*zab2(i2,i5,i1,i4)*
     &   zab2(i5,i1,i2,i3)*zab2(i5,i1,i2,i4)**2*zb(i4,i3))
     & )

      aaaa_NMHV_d(120) = (
     &  (s(i4,i5)*t(i1,i2,i4)*
     &   (zab2(i3,i1,i2,i4)**2*zab2(i5,i1,i4,i2)**3*zb(i4,i1) + 
     &     t(i1,i2,i4)**2*za(i3,i5)**2*zab2(i5,i2,i4,i1)*zb(i4,i2)**3))/
     & (za(i6,i5)*zab2(i5,i1,i2,i4)**2*zab2(i5,i1,i4,i2)*
     &   zab2(i5,i2,i4,i1)*zab2(i6,i1,i2,i4)*zb(i4,i1)*zb(i4,i2))
     & )

      aaaa_NMHV_d(121) = (
     &  -((s(i1,i6)*t(i1,i2,i4)**3*za(i3,i6)**2*zb(i2,i1)**2)/
     &   (za(i5,i6)*zab2(i5,i4,i2,i1)*zab2(i6,i1,i2,i4)*
     &     zab2(i6,i4,i2,i1)**2*zb(i4,i2)))
     & )

      aaaa_NMHV_d(123) = (
     &  -((s(i1,i5)*t(i1,i2,i4)**3*za(i3,i5)**2*zb(i2,i1)**2)/
     &   (za(i6,i5)*zab2(i5,i1,i2,i4)*zab2(i5,i4,i2,i1)**2*
     &     zab2(i6,i4,i2,i1)*zb(i4,i2)))
     & )

      aaaa_NMHV_d(124) = (
     &  (s(i1,i2)*zab2(i2,i1,i3,i6)**2*zab2(i4,i3,i6,i1)**2)/
     &  (za(i2,i5)*zab2(i2,i3,i6,i1)*zab2(i2,i4,i5,i3)*
     &    zab2(i5,i3,i6,i1)*zb(i3,i1)) + 
     & (s(i1,i2)*za(i1,i3)**2*zab2(i1,i3,i6,i2)*zb(i5,i2)**2)/
     &  (za(i1,i6)*zab2(i1,i3,i6,i4)*zab2(i6,i1,i3,i2)*zb(i4,i2))
     & )

      aaaa_NMHV_d(125) = (
     &  (s(i3,i6)*t(i1,i2,i6)*
     &   (-(za(i2,i6)*zab2(i1,i6,i2,i3)**3*zab2(i6,i1,i2,i5)**2) - 
     &     t(i1,i2,i6)**2*za(i1,i6)**3*zab2(i2,i6,i1,i3)*zb(i5,i3)**2))/
     & (za(i1,i6)*za(i2,i6)*zab2(i1,i6,i2,i3)*zab2(i2,i6,i1,i3)*
     &   zab2(i6,i1,i2,i3)**2*zab2(i6,i1,i2,i4)*zb(i3,i4))
     & )

      aaaa_NMHV_d(127) = (
     &  (s(i1,i2)*zab2(i2,i1,i3,i5)**2*zab2(i4,i3,i5,i1)**2)/
     &  (za(i2,i6)*zab2(i2,i3,i5,i1)*zab2(i2,i4,i6,i3)*
     &    zab2(i6,i3,i5,i1)*zb(i3,i1)) + 
     & (s(i1,i2)*za(i1,i3)**2*zab2(i1,i3,i5,i2)*zb(i6,i2)**2)/
     &  (za(i1,i5)*zab2(i1,i3,i5,i4)*zab2(i5,i1,i3,i2)*zb(i4,i2))
     & )

      aaaa_NMHV_d(128) = (
     &  (s(i3,i5)*t(i1,i2,i5)*
     &   (-(za(i2,i5)*zab2(i1,i5,i2,i3)**3*zab2(i5,i1,i2,i6)**2) - 
     &     t(i1,i2,i5)**2*za(i1,i5)**3*zab2(i2,i5,i1,i3)*zb(i6,i3)**2))/
     & (za(i1,i5)*za(i2,i5)*zab2(i1,i5,i2,i3)*zab2(i2,i5,i1,i3)*
     &   zab2(i5,i1,i2,i3)**2*zab2(i5,i1,i2,i4)*zb(i3,i4))
     & )

      aaaa_NMHV_d(130) = (
     &  -((t(i3,i4,i6)*za(i1,i4)*zab2(i1,i4,i3,i6)**2*zb(i4,i1))/
     &   (za(i2,i5)*zab2(i1,i4,i6,i3)*zab2(i5,i3,i6,i4)*zb(i3,i4)))
     & )

      aaaa_NMHV_d(131) = (
     &  -((t(i3,i4,i6)*za(i1,i3)*zab2(i1,i3,i4,i6)**2*zb(i3,i1))/
     &   (za(i2,i5)*zab2(i1,i3,i6,i4)*zab2(i5,i4,i6,i3)*zb(i4,i3)))
     & )

      aaaa_NMHV_d(134) = (
     &  -((t(i3,i4,i5)*za(i1,i4)*zab2(i1,i4,i3,i5)**2*zb(i4,i1))/
     &   (za(i2,i6)*zab2(i1,i4,i5,i3)*zab2(i6,i3,i5,i4)*zb(i3,i4)))
     & )

      aaaa_NMHV_d(135) = (
     &  -((t(i3,i4,i5)*za(i1,i3)*zab2(i1,i3,i4,i5)**2*zb(i3,i1))/
     &   (za(i2,i6)*zab2(i1,i3,i5,i4)*zab2(i6,i4,i5,i3)*zb(i4,i3)))
     & )

      aaaa_NMHV_d(141) = (
     &  (s(i2,i5)*t(i1,i2,i3)*zab2(i4,i1,i3,i2)**2)/
     & (za(i5,i6)*zab2(i5,i1,i2,i3)*zab2(i6,i1,i3,i2)*zb(i3,i1))
     & )

      aaaa_NMHV_d(142) = (
     &  -((s(i2,i5)**2*s(i4,i6)*t(i2,i4,i5)*za(i1,i3)**2)/
     &   (za(i2,i5)**2*zab2(i1,i2,i5,i4)*zab2(i6,i2,i5,i4)*
     &     zab2(i6,i4,i5,i2)))
     & )

      aaaa_NMHV_d(143) = (
     &  (s(i4,i6)*t(i2,i5,i6)**3)/
     & (za(i2,i5)*za(i5,i6)*zab2(i6,i2,i5,i4)*zb(i3,i1)*zb(i4,i3))
     & )

      aaaa_NMHV_d(145) = (
     &  (s(i2,i6)*t(i1,i2,i3)*zab2(i4,i1,i3,i2)**2)/
     & (za(i6,i5)*zab2(i5,i1,i3,i2)*zab2(i6,i1,i2,i3)*zb(i3,i1))
     & )

      aaaa_NMHV_d(146) = (
     &  -((s(i2,i6)**2*s(i4,i5)*t(i2,i4,i6)*za(i1,i3)**2)/
     &   (za(i2,i6)**2*zab2(i1,i2,i6,i4)*zab2(i5,i2,i6,i4)*
     &     zab2(i5,i4,i6,i2)))
     & )

      aaaa_NMHV_d(147) = (
     &  (s(i4,i5)*t(i2,i5,i6)**3)/
     & (za(i2,i6)*za(i6,i5)*zab2(i5,i2,i6,i4)*zb(i3,i1)*zb(i4,i3))
     & )

      aaaa_NMHV_d(161) = (
     &  (s(i2,i5)*t(i1,i2,i4)*zab2(i3,i1,i4,i2)**2)/
     & (za(i5,i6)*zab2(i5,i1,i2,i4)*zab2(i6,i1,i4,i2)*zb(i4,i1))
     & )

      aaaa_NMHV_d(162) = (
     &  -((s(i2,i5)**2*s(i3,i6)*t(i2,i3,i5)*za(i1,i4)**2)/
     &   (za(i2,i5)**2*zab2(i1,i2,i5,i3)*zab2(i6,i2,i5,i3)*
     &     zab2(i6,i3,i5,i2)))
     & )

      aaaa_NMHV_d(163) = (
     &  (s(i3,i6)*t(i2,i5,i6)**3)/
     & (za(i2,i5)*za(i5,i6)*zab2(i6,i2,i5,i3)*zb(i3,i4)*zb(i4,i1))
     & )

      aaaa_NMHV_d(164) = (
     &  (s(i2,i6)*t(i1,i2,i4)*zab2(i3,i1,i4,i2)**2)/
     & (za(i6,i5)*zab2(i5,i1,i4,i2)*zab2(i6,i1,i2,i4)*zb(i4,i1))
     & )

      aaaa_NMHV_d(165) = (
     &  -((s(i2,i6)**2*s(i3,i5)*t(i2,i3,i6)*za(i1,i4)**2)/
     &   (za(i2,i6)**2*zab2(i1,i2,i6,i3)*zab2(i5,i2,i6,i3)*
     &     zab2(i5,i3,i6,i2)))
     & )

      aaaa_NMHV_d(166) = (
     &  (s(i3,i5)*t(i2,i5,i6)**3)/
     & (za(i2,i6)*za(i6,i5)*zab2(i5,i2,i6,i3)*zb(i3,i4)*zb(i4,i1))
     & )

      aaaa_NMHV_d(172) = (
     &  -((s(i2,i3)*t(i1,i2,i5)**3*za(i1,i2)**2*zb(i6,i3)**2)/
     &   (za(i1,i5)*zab2(i2,i1,i5,i3)**2*zab2(i2,i1,i5,i4)*
     &     zab2(i5,i1,i2,i3)*zb(i4,i3)))
     & )

      aaaa_NMHV_d(173) = (
     &  -((s(i4,i6)*t(i2,i3,i4)*za(i2,i3)*zab2(i1,i2,i3,i4)**2*
     &     zab2(i6,i3,i4,i2)**2)/
     &   (s(i2,i3)*za(i1,i5)*za(i5,i6)*zab2(i6,i2,i3,i4)**3*zb(i4,i3)))
     & )

      aaaa_NMHV_d(174) = (
     &  (s(i4,i6)*t(i2,i3,i6)**3*za(i3,i6)**2*zb(i5,i4)**2)/
     & (zab2(i2,i3,i6,i4)*zab2(i6,i2,i3,i1)*zab2(i6,i2,i3,i4)**3)
     & )

      aaaa_NMHV_d(175) = (
     &  -((s(i2,i4)*t(i1,i2,i5)**3*za(i1,i2)**2*zb(i6,i4)**2)/
     &   (za(i1,i5)*zab2(i2,i1,i5,i3)*zab2(i2,i1,i5,i4)**2*
     &     zab2(i5,i1,i2,i4)*zb(i3,i4)))
     & )

      aaaa_NMHV_d(176) = (
     &  -((s(i3,i6)*t(i2,i3,i4)*za(i2,i4)*zab2(i1,i2,i4,i3)**2*
     &     zab2(i6,i4,i3,i2)**2)/
     &   (s(i2,i4)*za(i1,i5)*za(i5,i6)*zab2(i6,i2,i4,i3)**3*zb(i3,i4)))
     & )

      aaaa_NMHV_d(177) = (
     &  (s(i3,i6)*t(i2,i4,i6)**3*za(i4,i6)**2*zb(i5,i3)**2)/
     & (zab2(i2,i4,i6,i3)*zab2(i6,i2,i4,i1)*zab2(i6,i2,i4,i3)**3)
     & )

      aaaa_NMHV_d(184) = (
     &  -((s(i2,i3)*t(i1,i2,i6)**3*za(i1,i2)**2*zb(i5,i3)**2)/
     &   (za(i1,i6)*zab2(i2,i1,i6,i3)**2*zab2(i2,i1,i6,i4)*
     &     zab2(i6,i1,i2,i3)*zb(i4,i3)))
     & )

      aaaa_NMHV_d(185) = (
     &  -((s(i4,i5)*t(i2,i3,i4)*za(i2,i3)*zab2(i1,i2,i3,i4)**2*
     &     zab2(i5,i3,i4,i2)**2)/
     &   (s(i2,i3)*za(i1,i6)*za(i6,i5)*zab2(i5,i2,i3,i4)**3*zb(i4,i3)))
     & )

      aaaa_NMHV_d(186) = (
     &  (s(i4,i5)*t(i2,i3,i5)**3*za(i3,i5)**2*zb(i6,i4)**2)/
     & (zab2(i2,i3,i5,i4)*zab2(i5,i2,i3,i1)*zab2(i5,i2,i3,i4)**3)
     & )

      aaaa_NMHV_d(187) = (
     &  -((s(i2,i4)*t(i1,i2,i6)**3*za(i1,i2)**2*zb(i5,i4)**2)/
     &   (za(i1,i6)*zab2(i2,i1,i6,i3)*zab2(i2,i1,i6,i4)**2*
     &     zab2(i6,i1,i2,i4)*zb(i3,i4)))
     & )

      aaaa_NMHV_d(188) = (
     &  -((s(i3,i5)*t(i2,i3,i4)*za(i2,i4)*zab2(i1,i2,i4,i3)**2*
     &     zab2(i5,i4,i3,i2)**2)/
     &   (s(i2,i4)*za(i1,i6)*za(i6,i5)*zab2(i5,i2,i4,i3)**3*zb(i3,i4)))
     & )

      aaaa_NMHV_d(189) = (
     &  (s(i3,i5)*t(i2,i4,i5)**3*za(i4,i5)**2*zb(i6,i3)**2)/
     & (zab2(i2,i4,i5,i3)*zab2(i5,i2,i4,i1)*zab2(i5,i2,i4,i3)**3)
     & )

      return
      end
