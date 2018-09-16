!===== T. Dennen, May 2014
!===== Box coefficients for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,+)ga(i5,+)ga(i6,+)
      subroutine aaaa_MHV_d_init(i1,i2,i3,i4,i5,i6,za,zb,
     & aaaa_MHV_d)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: aaaa_MHV_d(195), zab2

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      aaaa_MHV_d(:) = czip

      aaaa_MHV_d(2) = (
     &  (s(i1,i2)*za(i1,i2)**2*za(i3,i6)**2*zb(i6,i1))/
     & (za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i6,i4)*za(i6,i5))
     & )

      aaaa_MHV_d(5) = (
     &  (s(i1,i2)*za(i1,i2)**2*za(i3,i5)**2*zb(i5,i1))/
     & (za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i4)*za(i5,i6))
     & )

      aaaa_MHV_d(10) = (
     &  (s(i1,i2)*za(i1,i2)**2*za(i3,i4)**2*zb(i4,i1))/
     & (za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i4,i6))
     & )

      aaaa_MHV_d(14) = (
     &  (s(i1,i2)*za(i1,i3)**2*za(i1,i4)*zb(i4,i2))/
     & (za(i1,i5)*za(i1,i6)*za(i4,i5)*za(i4,i6))
     & )

      aaaa_MHV_d(18) = (
     &  -((s(i3,i5)*s(i3,i6)*za(i1,i6)**2*za(i3,i5)**2)/
     &    (za(i1,i5)*za(i2,i4)*za(i4,i6)*za(i5,i6)**3)) - 
     & (s(i3,i5)*s(i3,i6)*za(i1,i5)**2*za(i3,i6)**2)/
     &  (za(i1,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6)**3) - 
     & (s(i3,i5)*s(i3,i6)*za(i1,i6)**2*za(i3,i5)**2)/
     &  (za(i1,i4)*za(i2,i6)*za(i5,i4)*za(i5,i6)**3) - 
     & (s(i3,i5)*s(i3,i6)*za(i1,i5)**2*za(i3,i6)**2)/
     &  (za(i1,i6)*za(i2,i4)*za(i5,i4)*za(i5,i6)**3)
     & )

      aaaa_MHV_d(20) = (
     &  (s(i1,i3)*s(i3,i6)*za(i1,i3)**2*za(i2,i6))/
     & (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i6,i4)*za(i6,i5))
     & )

      aaaa_MHV_d(22) = (
     &  (s(i1,i3)*s(i3,i5)*za(i1,i3)**2*za(i2,i5))/
     & (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i5,i4)*za(i5,i6))
     & )

      aaaa_MHV_d(23) = (
     &  (s(i1,i2)*za(i1,i3)**2*za(i1,i5)*zb(i5,i2))/
     & (za(i1,i4)*za(i1,i6)*za(i5,i4)*za(i5,i6))
     & )

      aaaa_MHV_d(26) = (
     &  -((s(i3,i4)*s(i3,i6)*za(i1,i6)**2*za(i3,i4)**2)/
     &    (za(i1,i5)*za(i2,i6)*za(i4,i5)*za(i4,i6)**3)) - 
     & (s(i3,i4)*s(i3,i6)*za(i1,i4)**2*za(i3,i6)**2)/
     &  (za(i1,i6)*za(i2,i5)*za(i4,i5)*za(i4,i6)**3) - 
     & (s(i3,i4)*s(i3,i6)*za(i1,i6)**2*za(i3,i4)**2)/
     &  (za(i1,i4)*za(i2,i5)*za(i4,i6)**3*za(i5,i6)) - 
     & (s(i3,i4)*s(i3,i6)*za(i1,i4)**2*za(i3,i6)**2)/
     &  (za(i1,i5)*za(i2,i4)*za(i4,i6)**3*za(i5,i6))
     & )

      aaaa_MHV_d(28) = (
     &  (s(i1,i3)*s(i3,i4)*za(i1,i3)**2*za(i2,i4))/
     & (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i4,i6))
     & )

      aaaa_MHV_d(29) = (
     &  (s(i1,i2)*za(i1,i3)**2*za(i1,i6)*zb(i6,i2))/
     & (za(i1,i4)*za(i1,i5)*za(i6,i4)*za(i6,i5))
     & )

      aaaa_MHV_d(32) = (
     &  -((s(i3,i4)*s(i3,i5)*za(i1,i5)**2*za(i3,i4)**2)/
     &    (za(i1,i6)*za(i2,i5)*za(i4,i5)**3*za(i4,i6))) - 
     & (s(i3,i4)*s(i3,i5)*za(i1,i4)**2*za(i3,i5)**2)/
     &  (za(i1,i5)*za(i2,i6)*za(i4,i5)**3*za(i4,i6)) - 
     & (s(i3,i4)*s(i3,i5)*za(i1,i5)**2*za(i3,i4)**2)/
     &  (za(i1,i4)*za(i2,i6)*za(i4,i5)**3*za(i6,i5)) - 
     & (s(i3,i4)*s(i3,i5)*za(i1,i4)**2*za(i3,i5)**2)/
     &  (za(i1,i6)*za(i2,i4)*za(i4,i5)**3*za(i6,i5))
     & )

      aaaa_MHV_d(51) = (
     &  (s(i2,i3)*s(i3,i6)*za(i1,i2)**2*za(i1,i6)*za(i3,i6)**2)/
     & (za(i1,i4)*za(i1,i5)*za(i2,i6)**3*za(i6,i4)*za(i6,i5))
     & )

      aaaa_MHV_d(54) = (
     &  (s(i2,i3)*s(i3,i5)*za(i1,i2)**2*za(i1,i5)*za(i3,i5)**2)/
     & (za(i1,i4)*za(i1,i6)*za(i2,i5)**3*za(i5,i4)*za(i5,i6))
     & )

      aaaa_MHV_d(58) = (
     &  (s(i2,i3)*s(i3,i4)*za(i1,i2)**2*za(i1,i4)*za(i3,i4)**2)/
     & (za(i1,i5)*za(i1,i6)*za(i2,i4)**3*za(i4,i5)*za(i4,i6))
     & )

      aaaa_MHV_d(68) = (
     &  (za(i1,i6)**2*za(i3,i4)**2*zab2(i4,i1,i2,i6)*zab2(i6,i3,i5,i4))/
     &  (za(i1,i4)*za(i2,i6)*za(i4,i5)*za(i6,i4)**2*za(i6,i5)) + 
     & (za(i1,i4)**2*za(i3,i6)**2*zab2(i4,i1,i2,i6)*zab2(i6,i3,i5,i4))/
     &  (za(i1,i6)*za(i2,i4)*za(i4,i5)*za(i6,i4)**2*za(i6,i5))
     & )

      aaaa_MHV_d(69) = (
     &  (za(i1,i5)**2*za(i3,i4)**2*zab2(i4,i3,i6,i5)*zab2(i5,i1,i2,i4))/
     &  (za(i1,i4)*za(i2,i5)*za(i4,i5)**2*za(i4,i6)*za(i5,i6)) + 
     & (za(i1,i4)**2*za(i3,i5)**2*zab2(i4,i3,i6,i5)*zab2(i5,i1,i2,i4))/
     &  (za(i1,i5)*za(i2,i4)*za(i4,i5)**2*za(i4,i6)*za(i5,i6))
     & )

      aaaa_MHV_d(70) = (
     &  (za(i1,i3)**2*zab2(i1,i2,i4,i5)*zab2(i5,i3,i6,i1))/
     & (za(i1,i6)*za(i2,i4)*za(i4,i5)*za(i6,i5))
     & )

      aaaa_MHV_d(71) = (
     &  (za(i1,i3)**2*zab2(i1,i2,i4,i6)*zab2(i6,i3,i5,i1))/
     & (za(i1,i5)*za(i2,i4)*za(i4,i6)*za(i5,i6))
     & )

      aaaa_MHV_d(73) = (
     &  (za(i1,i6)**2*za(i3,i5)**2*zab2(i5,i3,i4,i6)*zab2(i6,i1,i2,i5))/
     &  (za(i1,i5)*za(i2,i6)*za(i5,i4)*za(i5,i6)**2*za(i6,i4)) + 
     & (za(i1,i5)**2*za(i3,i6)**2*zab2(i5,i3,i4,i6)*zab2(i6,i1,i2,i5))/
     &  (za(i1,i6)*za(i2,i5)*za(i5,i4)*za(i5,i6)**2*za(i6,i4))
     & )

      aaaa_MHV_d(74) = (
     &  (za(i1,i3)**2*zab2(i1,i2,i5,i4)*zab2(i4,i3,i6,i1))/
     & (za(i1,i6)*za(i2,i5)*za(i5,i4)*za(i6,i4))
     & )

      aaaa_MHV_d(75) = (
     &  (za(i1,i3)**2*zab2(i1,i2,i5,i6)*zab2(i6,i3,i4,i1))/
     & (za(i1,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6))
     & )

      aaaa_MHV_d(77) = (
     &  (za(i1,i3)**2*zab2(i1,i2,i6,i4)*zab2(i4,i3,i5,i1))/
     & (za(i1,i5)*za(i2,i6)*za(i5,i4)*za(i6,i4))
     & )

      aaaa_MHV_d(78) = (
     &  (za(i1,i3)**2*zab2(i1,i2,i6,i5)*zab2(i5,i3,i4,i1))/
     & (za(i1,i4)*za(i2,i6)*za(i4,i5)*za(i6,i5))
     & )

      aaaa_MHV_d(88) = (
     &  (za(i1,i2)**2*za(i3,i6)**2*zab2(i2,i3,i5,i6)*zab2(i6,i1,i4,i2))/
     & (za(i1,i4)*za(i2,i5)*za(i2,i6)**2*za(i4,i6)*za(i5,i6))
     & )

      aaaa_MHV_d(89) = (
     &  (za(i1,i2)**2*za(i3,i5)**2*zab2(i2,i3,i6,i5)*zab2(i5,i1,i4,i2))/
     & (za(i1,i4)*za(i2,i5)**2*za(i2,i6)*za(i4,i5)*za(i6,i5))
     & )

      aaaa_MHV_d(94) = (
     &  (za(i1,i2)**2*za(i3,i6)**2*zab2(i2,i3,i4,i6)*zab2(i6,i1,i5,i2))/
     & (za(i1,i5)*za(i2,i4)*za(i2,i6)**2*za(i4,i6)*za(i5,i6))
     & )

      aaaa_MHV_d(95) = (
     &  (za(i1,i2)**2*za(i3,i4)**2*zab2(i2,i3,i6,i4)*zab2(i4,i1,i5,i2))/
     & (za(i1,i5)*za(i2,i4)**2*za(i2,i6)*za(i5,i4)*za(i6,i4))
     & )

      aaaa_MHV_d(100) = (
     &  (za(i1,i2)**2*za(i3,i5)**2*zab2(i2,i3,i4,i5)*zab2(i5,i1,i6,i2))/
     & (za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i4,i5)*za(i6,i5))
     & )

      aaaa_MHV_d(101) = (
     &  (za(i1,i2)**2*za(i3,i4)**2*zab2(i2,i3,i5,i4)*zab2(i4,i1,i6,i2))/
     & (za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i5,i4)*za(i6,i4))
     & )

      return
      end
