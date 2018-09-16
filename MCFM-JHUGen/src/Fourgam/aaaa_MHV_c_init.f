!===== T. Dennen, May 2014
!===== Triangle coefficients for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,+)ga(i5,+)ga(i6,+)
      subroutine aaaa_MHV_c_init(i1,i2,i3,i4,i5,i6,za,zb,
     & aaaa_MHV_c)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
c      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: aaaa_MHV_c(90), zab2


      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)


      aaaa_MHV_c(:) = czip

      aaaa_MHV_c(2) = (
     &  -((za(i1,i2)**2*za(i3,i6)**2*zb(i6,i1))/
     &   (za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i6,i4)*za(i6,i5)))
     &)


      aaaa_MHV_c(4) = (
     &  -((za(i1,i2)**2*za(i3,i5)**2*zb(i5,i1))/
     &   (za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i4)*za(i5,i6)))
     &)


      aaaa_MHV_c(6) = (
     &  -((za(i1,i2)**2*za(i3,i4)**2*zb(i4,i1))/
     &   (za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i4,i6)))
     &)


      aaaa_MHV_c(8) = (
     &  (-2*za(i1,i2)*za(i1,i3)*za(i3,i4)*za(i3,i6)*zb(i6,i3))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i4,i6)**2) - 
     & (za(i1,i2)**2*za(i1,i6)*za(i3,i4)**2*za(i3,i6)*zb(i6,i3))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i6)**2)
     &   - (za(i1,i2)**2*za(i1,i6)*za(i3,i5)**2*za(i3,i6)*zb(i6,i3))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i6)**2)
     &   - (za(i1,i2)*za(i1,i3)*za(i1,i5)*za(i3,i6)**2*zb(i6,i3))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i5,i6)**2) - 
     & (za(i1,i6)*za(i3,i5)*za(i3,i6)**2*zb(i6,i3))/
     &  (za(i2,i5)*za(i4,i6)**2*za(i5,i6)**2) - 
     & (za(i1,i6)*za(i3,i6)**3*zb(i6,i3))/
     &  (za(i2,i6)*za(i4,i6)**2*za(i5,i6)**2) - 
     & (za(i1,i2)*za(i1,i3)*za(i1,i6)*za(i3,i5)*za(i3,i6)*za(i4,i5)*
     &    zb(i6,i3))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6)**2)
     &   - (za(i1,i4)*za(i2,i3)*za(i3,i6)**2*zb(i6,i3))/
     &  (za(i2,i4)*za(i2,i5)*za(i4,i6)**2*za(i5,i6)) + 
     & (za(i1,i3)*za(i1,i6)**2*za(i2,i3)*za(i3,i6)*zb(i6,i3))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i4,i6)*za(i5,i6)) - 
     & (za(i1,i3)**2*za(i2,i6)*za(i3,i6)*zb(i6,i3))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6)) + 
     & (za(i1,i2)*za(i1,i6)**2*za(i2,i3)*za(i3,i6)**2*zb(i6,i3))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)**3*za(i4,i6)*za(i5,i6))
     &)


      aaaa_MHV_c(9) = (
     &  -((za(i1,i2)**2*za(i1,i3)**3*zb(i3,i1))/
     &   (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)))
     &)


      aaaa_MHV_c(10) = (
     &  (-2*za(i1,i2)*za(i1,i3)*za(i3,i4)*za(i3,i5)*zb(i5,i3))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i4,i5)**2) - 
     & (za(i1,i2)**2*za(i1,i5)*za(i3,i4)**2*za(i3,i5)*zb(i5,i3))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)**2)
     &   - (za(i1,i2)*za(i1,i3)*za(i1,i6)*za(i3,i5)**2*zb(i5,i3))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i6,i5)**2) - 
     & (za(i1,i2)**2*za(i1,i5)*za(i3,i5)*za(i3,i6)**2*zb(i5,i3))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i6,i5)**2)
     &   - (za(i1,i5)*za(i3,i5)**3*zb(i5,i3))/
     &  (za(i2,i5)*za(i4,i5)**2*za(i6,i5)**2) - 
     & (za(i1,i5)*za(i3,i5)**2*za(i3,i6)*zb(i5,i3))/
     &  (za(i2,i6)*za(i4,i5)**2*za(i6,i5)**2) - 
     & (za(i1,i2)*za(i1,i3)*za(i1,i5)*za(i3,i5)*za(i3,i6)*za(i4,i6)*
     &    zb(i5,i3))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i4,i5)*za(i6,i5)**2)
     &   - (za(i1,i4)*za(i2,i3)*za(i3,i5)**2*zb(i5,i3))/
     &  (za(i2,i4)*za(i2,i6)*za(i4,i5)**2*za(i6,i5)) + 
     & (za(i1,i3)*za(i1,i5)**2*za(i2,i3)*za(i3,i5)*zb(i5,i3))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i4,i5)*za(i6,i5)) - 
     & (za(i1,i3)**2*za(i2,i5)*za(i3,i5)*zb(i5,i3))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i5)*za(i6,i5)) + 
     & (za(i1,i2)*za(i1,i5)**2*za(i2,i3)*za(i3,i5)**2*zb(i5,i3))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)**3*za(i4,i5)*za(i6,i5))
     &)


      aaaa_MHV_c(11) = (
     &  -((za(i1,i2)*za(i1,i3)*za(i1,i5)*za(i3,i4)**2*zb(i4,i3))/
     &    (za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i5,i4)**2)) - 
     & (za(i1,i2)**2*za(i1,i4)*za(i3,i4)*za(i3,i5)**2*zb(i4,i3))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i4)**2)
     &   - (2*za(i1,i2)*za(i1,i3)*za(i3,i4)*za(i3,i6)*zb(i4,i3))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i6,i4)**2) - 
     & (za(i1,i2)**2*za(i1,i4)*za(i3,i4)*za(i3,i6)**2*zb(i4,i3))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i6,i4)**2)
     &   - (za(i1,i4)*za(i3,i4)**3*zb(i4,i3))/
     &  (za(i2,i4)*za(i5,i4)**2*za(i6,i4)**2) - 
     & (za(i1,i4)*za(i3,i4)**2*za(i3,i5)*zb(i4,i3))/
     &  (za(i2,i5)*za(i5,i4)**2*za(i6,i4)**2) - 
     & (za(i1,i6)*za(i2,i3)*za(i3,i4)**2*zb(i4,i3))/
     &  (za(i2,i5)*za(i2,i6)*za(i5,i4)*za(i6,i4)**2) + 
     & (za(i1,i3)*za(i1,i4)**2*za(i2,i3)*za(i3,i4)*zb(i4,i3))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i5,i4)*za(i6,i4)) - 
     & (za(i1,i3)**2*za(i2,i4)*za(i3,i4)*zb(i4,i3))/
     &  (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i5,i4)*za(i6,i4)) + 
     & (za(i1,i2)*za(i1,i4)**2*za(i2,i3)*za(i3,i4)**2*zb(i4,i3))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**3*za(i5,i4)*za(i6,i4)) - 
     & (za(i1,i2)*za(i1,i3)*za(i1,i4)*za(i3,i4)*za(i3,i5)*za(i6,i5)*
     &    zb(i4,i3))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i5,i4)**2*za(i6,i4))
     &)


      aaaa_MHV_c(12) = (
     &  -((za(i1,i3)**2*za(i1,i6)*zb(i6,i2))/
     &   (za(i1,i4)*za(i1,i5)*za(i6,i4)*za(i6,i5)))
     &)


      aaaa_MHV_c(13) = (
     &  -((za(i1,i3)**2*za(i1,i5)*zb(i5,i2))/
     &   (za(i1,i4)*za(i1,i6)*za(i5,i4)*za(i5,i6)))
     &)


      aaaa_MHV_c(14) = (
     &  -((za(i1,i3)**2*za(i1,i4)*zb(i4,i2))/
     &   (za(i1,i5)*za(i1,i6)*za(i4,i5)*za(i4,i6)))
     &)


      aaaa_MHV_c(15) = (
     &  -((za(i1,i2)**2*za(i1,i4)*za(i2,i3)**3*zb(i3,i2))/
     &    (za(i1,i5)*za(i1,i6)*za(i2,i4)**3*za(i2,i5)*za(i2,i6))) + 
     & (za(i1,i2)**3*za(i2,i3)**2*za(i3,i5)*zb(i3,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**3*za(i2,i6)) + 
     & (za(i1,i2)**3*za(i2,i3)**2*za(i3,i5)*zb(i3,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i2,i5)**2*za(i2,i6)) + 
     & (za(i1,i2)**3*za(i2,i3)**2*za(i3,i6)*zb(i3,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)**2) + 
     & (za(i1,i2)**3*za(i2,i3)**2*za(i3,i6)*zb(i3,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)**2) - 
     & (za(i1,i2)**4*za(i2,i3)*za(i3,i6)**2*zb(i3,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)**3)
     &)


      aaaa_MHV_c(19) = (
     &  (2*za(i1,i6)*za(i3,i6)**2*zab2(i6,i1,i2,i6))/
     & (za(i2,i6)*za(i6,i4)**2*za(i6,i5)**2)
     &)


      aaaa_MHV_c(22) = (
     &  (2*za(i1,i5)*za(i3,i5)**2*zab2(i5,i1,i2,i5))/
     & (za(i2,i5)*za(i5,i4)**2*za(i5,i6)**2)
     &)


      aaaa_MHV_c(25) = (
     &  (2*za(i1,i4)*za(i3,i4)**2*zab2(i4,i1,i2,i4))/
     & (za(i2,i4)*za(i4,i5)**2*za(i4,i6)**2)
     &)


      aaaa_MHV_c(29) = (
     &  -((za(i1,i2)*za(i3,i6)**2*zab2(i1,i3,i5,i6))/
     &    (za(i1,i4)*za(i2,i4)*za(i2,i6)*za(i5,i6)**2)) - 
     & (za(i1,i6)**2*za(i3,i6)*zab2(i3,i3,i5,i6))/
     &  (za(i1,i4)*za(i2,i6)*za(i4,i6)*za(i5,i6)**2) + 
     & (za(i1,i6)*za(i3,i6)**2*zab2(i4,i3,i5,i6))/
     &  (za(i2,i4)*za(i4,i6)**2*za(i5,i6)**2) + 
     & (za(i1,i2)*za(i1,i3)*za(i3,i6)*zab2(i6,i3,i5,i6))/
     &  (za(i1,i4)*za(i2,i4)*za(i2,i6)*za(i5,i6)**2) + 
     & (za(i1,i6)**2*za(i2,i4)*za(i3,i6)**2*zab2(i6,i3,i5,i6))/
     &  (za(i1,i4)*za(i2,i6)**2*za(i4,i6)**2*za(i5,i6)**2) + 
     & (za(i1,i3)*za(i3,i6)*zab2(i6,i3,i5,i6))/
     &  (za(i2,i4)*za(i4,i6)*za(i5,i6)**2)
     &)


      aaaa_MHV_c(30) = (
     &  -((za(i1,i2)*za(i3,i5)**2*zab2(i1,i3,i6,i5))/
     &    (za(i1,i4)*za(i2,i4)*za(i2,i5)*za(i6,i5)**2)) - 
     & (za(i1,i5)**2*za(i3,i5)*zab2(i3,i3,i6,i5))/
     &  (za(i1,i4)*za(i2,i5)*za(i4,i5)*za(i6,i5)**2) + 
     & (za(i1,i5)*za(i3,i5)**2*zab2(i4,i3,i6,i5))/
     &  (za(i2,i4)*za(i4,i5)**2*za(i6,i5)**2) + 
     & (za(i1,i2)*za(i1,i3)*za(i3,i5)*zab2(i5,i3,i6,i5))/
     &  (za(i1,i4)*za(i2,i4)*za(i2,i5)*za(i6,i5)**2) + 
     & (za(i1,i5)**2*za(i2,i4)*za(i3,i5)**2*zab2(i5,i3,i6,i5))/
     &  (za(i1,i4)*za(i2,i5)**2*za(i4,i5)**2*za(i6,i5)**2) + 
     & (za(i1,i3)*za(i3,i5)*zab2(i5,i3,i6,i5))/
     &  (za(i2,i4)*za(i4,i5)*za(i6,i5)**2)
     &)


      aaaa_MHV_c(35) = (
     &  -((za(i1,i2)*za(i3,i6)**2*zab2(i1,i3,i4,i6))/
     &    (za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i4,i6)**2)) - 
     & (za(i1,i6)**2*za(i3,i6)*zab2(i3,i3,i4,i6))/
     &  (za(i1,i5)*za(i2,i6)*za(i4,i6)**2*za(i5,i6)) + 
     & (za(i1,i6)*za(i3,i6)**2*zab2(i5,i3,i4,i6))/
     &  (za(i2,i5)*za(i4,i6)**2*za(i5,i6)**2) + 
     & (za(i1,i2)*za(i1,i3)*za(i3,i6)*zab2(i6,i3,i4,i6))/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i4,i6)**2) + 
     & (za(i1,i6)**2*za(i2,i5)*za(i3,i6)**2*zab2(i6,i3,i4,i6))/
     &  (za(i1,i5)*za(i2,i6)**2*za(i4,i6)**2*za(i5,i6)**2) + 
     & (za(i1,i3)*za(i3,i6)*zab2(i6,i3,i4,i6))/
     &  (za(i2,i5)*za(i4,i6)**2*za(i5,i6))
     &)


      aaaa_MHV_c(36) = (
     &  -((za(i1,i2)*za(i3,i4)**2*zab2(i1,i3,i6,i4))/
     &    (za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i6,i4)**2)) - 
     & (za(i1,i4)**2*za(i3,i4)*zab2(i3,i3,i6,i4))/
     &  (za(i1,i5)*za(i2,i4)*za(i5,i4)*za(i6,i4)**2) + 
     & (za(i1,i2)*za(i1,i3)*za(i3,i4)*zab2(i4,i3,i6,i4))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i6,i4)**2) + 
     & (za(i1,i4)**2*za(i2,i5)*za(i3,i4)**2*zab2(i4,i3,i6,i4))/
     &  (za(i1,i5)*za(i2,i4)**2*za(i5,i4)**2*za(i6,i4)**2) + 
     & (za(i1,i3)*za(i3,i4)*zab2(i4,i3,i6,i4))/
     &  (za(i2,i5)*za(i5,i4)*za(i6,i4)**2) + 
     & (za(i1,i4)*za(i3,i4)**2*zab2(i5,i3,i6,i4))/
     &  (za(i2,i5)*za(i5,i4)**2*za(i6,i4)**2)
     &)


      aaaa_MHV_c(40) = (
     &  -((za(i1,i2)*za(i3,i5)**2*zab2(i1,i3,i4,i5))/
     &    (za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i4,i5)**2)) - 
     & (za(i1,i5)**2*za(i3,i5)*zab2(i3,i3,i4,i5))/
     &  (za(i1,i6)*za(i2,i5)*za(i4,i5)**2*za(i6,i5)) + 
     & (za(i1,i2)*za(i1,i3)*za(i3,i5)*zab2(i5,i3,i4,i5))/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i4,i5)**2) + 
     & (za(i1,i5)**2*za(i2,i6)*za(i3,i5)**2*zab2(i5,i3,i4,i5))/
     &  (za(i1,i6)*za(i2,i5)**2*za(i4,i5)**2*za(i6,i5)**2) + 
     & (za(i1,i3)*za(i3,i5)*zab2(i5,i3,i4,i5))/
     &  (za(i2,i6)*za(i4,i5)**2*za(i6,i5)) + 
     & (za(i1,i5)*za(i3,i5)**2*zab2(i6,i3,i4,i5))/
     &  (za(i2,i6)*za(i4,i5)**2*za(i6,i5)**2)
     &)


      aaaa_MHV_c(41) = (
     &  -((za(i1,i2)*za(i3,i4)**2*zab2(i1,i3,i5,i4))/
     &    (za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i5,i4)**2)) - 
     & (za(i1,i4)**2*za(i3,i4)*zab2(i3,i3,i5,i4))/
     &  (za(i1,i6)*za(i2,i4)*za(i5,i4)**2*za(i6,i4)) + 
     & (za(i1,i2)*za(i1,i3)*za(i3,i4)*zab2(i4,i3,i5,i4))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i5,i4)**2) + 
     & (za(i1,i4)**2*za(i2,i6)*za(i3,i4)**2*zab2(i4,i3,i5,i4))/
     &  (za(i1,i6)*za(i2,i4)**2*za(i5,i4)**2*za(i6,i4)**2) + 
     & (za(i1,i3)*za(i3,i4)*zab2(i4,i3,i5,i4))/
     &  (za(i2,i6)*za(i5,i4)**2*za(i6,i4)) + 
     & (za(i1,i4)*za(i3,i4)**2*zab2(i6,i3,i5,i4))/
     &  (za(i2,i6)*za(i5,i4)**2*za(i6,i4)**2)
     &)


      aaaa_MHV_c(43) = (
     &  -((za(i1,i3)**2*za(i2,i6)*zab2(i6,i1,i3,i6))/
     &   (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i6)*za(i6,i5)))
     &)


      aaaa_MHV_c(44) = (
     &  -((za(i1,i3)**2*za(i2,i5)*zab2(i5,i1,i3,i5))/
     &   (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i5)*za(i5,i6)))
     &)


      aaaa_MHV_c(45) = (
     &  -((za(i1,i3)**2*za(i2,i4)*zab2(i4,i1,i3,i4))/
     &   (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i6)*za(i5,i4)))
     &)


      aaaa_MHV_c(47) = (
     &  (za(i1,i3)**2*zab2(i6,i2,i5,i6))/
     & (za(i1,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6))
     &)


      aaaa_MHV_c(48) = (
     &  (za(i1,i3)**2*zab2(i5,i2,i6,i5))/
     & (za(i1,i4)*za(i2,i6)*za(i4,i5)*za(i6,i5))
     &)


      aaaa_MHV_c(50) = (
     &  (za(i1,i3)**2*zab2(i6,i2,i4,i6))/
     & (za(i1,i5)*za(i2,i4)*za(i4,i6)*za(i5,i6))
     &)


      aaaa_MHV_c(51) = (
     &  (za(i1,i3)**2*zab2(i4,i2,i6,i4))/
     & (za(i1,i5)*za(i2,i6)*za(i5,i4)*za(i6,i4))
     &)


      aaaa_MHV_c(53) = (
     &  (za(i1,i3)**2*zab2(i5,i2,i4,i5))/
     & (za(i1,i6)*za(i2,i4)*za(i4,i5)*za(i6,i5))
     &)


      aaaa_MHV_c(54) = (
     &  (za(i1,i3)**2*zab2(i4,i2,i5,i4))/
     & (za(i1,i6)*za(i2,i5)*za(i5,i4)*za(i6,i4))
     &)


      aaaa_MHV_c(55) = (
     &  (za(i1,i2)**2*za(i2,i3)**2*zab2(i2,i1,i4,i2))/
     & (za(i1,i4)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)**2)
     &)


      aaaa_MHV_c(56) = (
     &  -((za(i1,i2)**2*za(i3,i6)**2*zab2(i6,i1,i4,i6))/
     &   (za(i1,i4)*za(i2,i5)*za(i2,i6)**2*za(i4,i6)*za(i6,i5)))
     &)


      aaaa_MHV_c(57) = (
     &  -((za(i1,i2)**2*za(i3,i5)**2*zab2(i5,i1,i4,i5))/
     &   (za(i1,i4)*za(i2,i5)**2*za(i2,i6)*za(i4,i5)*za(i5,i6)))
     &)


      aaaa_MHV_c(59) = (
     &  (za(i1,i2)**2*za(i2,i3)**2*zab2(i2,i3,i6,i2))/
     &  (za(i1,i5)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)**2) - 
     & (za(i1,i2)**3*za(i2,i3)*za(i3,i5)*zab2(i2,i3,i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)**2) - 
     & (za(i1,i2)**3*za(i2,i3)*za(i3,i6)*zab2(i2,i3,i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)**3)
     &)


      aaaa_MHV_c(60) = (
     &  (za(i1,i2)**2*za(i1,i6)*za(i3,i6)**2*zab2(i6,i2,i3,i6))/
     & (za(i1,i4)*za(i1,i5)*za(i2,i6)**3*za(i4,i6)*za(i5,i6))
     &)


      aaaa_MHV_c(62) = (
     &  (za(i1,i2)**2*za(i2,i3)**2*zab2(i2,i3,i5,i2))/
     &  (za(i1,i6)*za(i2,i4)**2*za(i2,i5)**2*za(i2,i6)) - 
     & (za(i1,i2)**3*za(i2,i3)*za(i3,i5)*zab2(i2,i3,i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**3*za(i2,i6)) - 
     & (za(i1,i2)**3*za(i2,i3)*za(i3,i6)*zab2(i2,i3,i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)**2)
     &)


      aaaa_MHV_c(63) = (
     &  (za(i1,i2)**2*za(i1,i5)*za(i3,i5)**2*zab2(i5,i2,i3,i5))/
     & (za(i1,i4)*za(i1,i6)*za(i2,i5)**3*za(i4,i5)*za(i6,i5))
     &)


      aaaa_MHV_c(65) = (
     &  (za(i1,i2)**2*za(i2,i3)**2*zab2(i2,i1,i5,i2))/
     & (za(i1,i5)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)**2)
     &)


      aaaa_MHV_c(66) = (
     &  -((za(i1,i2)**2*za(i3,i6)**2*zab2(i6,i1,i5,i6))/
     &   (za(i1,i5)*za(i2,i4)*za(i2,i6)**2*za(i5,i6)*za(i6,i4)))
     &)


      aaaa_MHV_c(67) = (
     &  -((za(i1,i2)**2*za(i3,i4)**2*zab2(i4,i1,i5,i4))/
     &   (za(i1,i5)*za(i2,i4)**2*za(i2,i6)*za(i4,i6)*za(i5,i4)))
     &)


      aaaa_MHV_c(69) = (
     &  (za(i1,i2)**2*za(i2,i3)**2*zab2(i2,i3,i4,i2))/
     &  (za(i1,i6)*za(i2,i4)**2*za(i2,i5)**2*za(i2,i6)) - 
     & (za(i1,i2)**3*za(i2,i3)*za(i3,i4)*zab2(i2,i3,i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**3*za(i2,i5)*za(i2,i6)) - 
     & (za(i1,i2)**3*za(i2,i3)*za(i3,i6)*zab2(i2,i3,i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)**2)
     &)


      aaaa_MHV_c(70) = (
     &  (za(i1,i2)**2*za(i1,i4)*za(i3,i4)**2*zab2(i4,i2,i3,i4))/
     & (za(i1,i5)*za(i1,i6)*za(i2,i4)**3*za(i5,i4)*za(i6,i4))
     &)


      aaaa_MHV_c(72) = (
     &  (za(i1,i2)**2*za(i2,i3)**2*zab2(i2,i1,i6,i2))/
     & (za(i1,i6)*za(i2,i4)**2*za(i2,i5)**2*za(i2,i6))
     &)


      aaaa_MHV_c(73) = (
     &  -((za(i1,i2)**2*za(i3,i5)**2*zab2(i5,i1,i6,i5))/
     &   (za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i5,i4)*za(i6,i5)))
     &)


      aaaa_MHV_c(74) = (
     &  -((za(i1,i2)**2*za(i3,i4)**2*zab2(i4,i1,i6,i4))/
     &   (za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i4,i5)*za(i6,i4)))
     &)


      return
      end
