!===== T. Dennen, May 2014
!===== Rational terms for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,+)ga(i5,+)ga(i6,+)
      subroutine aaaa_MHV_r_init(i1,i2,i3,i4,i5,i6,za,zb,aaaa_MHV_r)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: aaaa_MHV_r, zab2, t

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      aaaa_MHV_r = czip

      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i3)*za(i2,i5)*zab2(i5,i1,i3,i5))/
     & (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i5,i4)*za(i5,i6)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i3)*za(i2,i6)*zab2(i6,i1,i3,i6))/
     & (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i6,i4)*za(i6,i5)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i3)*za(i2,i4)*zab2(i4,i1,i3,i4))/
     & (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i4,i6)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i2,i5)*za(i6,i5)**2*zab2(i2,i3,i6,i2)*zab2(i5,i1,i4,i3)*
     &      zab2(i5,i1,i4,i6)*zab2(i6,i1,i4,i5)*
     &      (za(i2,i3)*(za(i4,i5)*zab2(i1,i1,i4,i5)*
     &            (-(za(i4,i5)*zab2(i1,i1,i4,i5)) - 
     &              za(i1,i5)*zab2(i4,i1,i4,i5))*zab2(i5,i1,i4,i4)*
     &            zb(i5,i2) + 
     &           za(i1,i4)*za(i1,i5)*zab2(i4,i1,i4,i5)*
     &            zab2(i5,i1,i4,i2)*zab2(i5,i1,i4,i5)*zb(i5,i4)) + 
     &        za(i3,i6)*(-(za(i1,i4)*za(i1,i5)*zab2(i4,i1,i4,i5)*
     &              zab2(i5,i1,i4,i5)*zab2(i5,i1,i4,i6)*zb(i5,i4)) - 
     &           za(i4,i5)*zab2(i1,i1,i4,i5)*
     &            (-(za(i4,i5)*zab2(i1,i1,i4,i5)) - 
     &              za(i1,i5)*zab2(i4,i1,i4,i5))*zab2(i5,i1,i4,i4)*
     &            zb(i5,i6)))*zb(i6,i2)) + 
     &   t(i2,i3,i6)*(-(za(i2,i5)*za(i4,i5)**2*za(i6,i5)**2*
     &         zab2(i1,i1,i4,i5)**2*zab2(i3,i1,i4,i5)*zab2(i5,i1,i4,i3)*
     &         zab2(i5,i1,i4,i4)*zab2(i5,i1,i4,i6)*
     &         (za(i3,i6)*zb(i3,i2)*zb(i5,i6) + 
     &           za(i2,i6)*zb(i5,i2)*zb(i6,i2))) - 
     &      za(i1,i5)*za(i6,i5)*zab2(i4,i1,i4,i5)*zab2(i5,i1,i4,i6)*
     &       zab2(i6,i1,i4,i5)*
     &       (za(i2,i3)*za(i2,i5)*za(i4,i5)*za(i6,i5)*zab2(i1,i1,i4,i5)*
     &          zab2(i5,i1,i4,i4)*zab2(i5,i1,i4,i6)*zb(i3,i2)*zb(i5,i2)
     &          - za(i2,i5)*za(i3,i6)*za(i4,i5)*za(i6,i5)*
     &          zab2(i1,i1,i4,i5)*zab2(i5,i1,i4,i4)*zab2(i5,i1,i4,i6)*
     &          zb(i3,i2)*zb(i5,i6) + 
     &         za(i3,i5)*(za(i1,i4)*za(i2,i5)*zab2(i5,i1,i4,i3)*
     &             zab2(i5,i1,i4,i5)*zb(i5,i4)*
     &             (-(za(i3,i6)*zab2(i5,i1,i4,i6)*zb(i3,i2)) - 
     &               za(i2,i6)*zab2(i5,i1,i4,i2)*zb(i6,i2)) + 
     &            za(i2,i6)*za(i4,i5)*zab2(i5,i1,i4,i4)*zb(i5,i2)*
     &             (za(i1,i2)*zab2(i5,i1,i4,i5)*zab2(i5,i1,i4,i6)*
     &                zb(i3,i2) + 
     &               za(i2,i5)*zab2(i1,i1,i4,i5)*zab2(i5,i1,i4,i3)*
     &                zb(i6,i2)))) - 
     &      za(i1,i5)**2*za(i4,i5)*zab2(i4,i1,i4,i5)*zab2(i5,i1,i4,i4)*
     &       zab2(i6,i1,i4,i5)*
     &       (za(i2,i6)*zb(i5,i6)*
     &          (-((za(i2,i5)*za(i6,i5)*zab2(i3,i1,i4,i5)*
     &                  zab2(i5,i1,i4,i2) + 
     &                 za(i3,i5)*za(i3,i6)*zab2(i5,i1,i4,i3)*
     &                  zab2(i5,i1,i4,i5))*zab2(i5,i1,i4,i6)*zb(i3,i2))
     &             - za(i2,i5)*za(i3,i6)*zab2(i5,i1,i4,i2)*
     &             zab2(i5,i1,i4,i3)*zab2(i5,i1,i4,i5)*zb(i6,i2)) + 
     &         za(i2,i3)*za(i6,i5)**2*zab2(i3,i1,i4,i5)*
     &          zab2(i5,i1,i4,i5)*zab2(i5,i1,i4,i6)*zb(i3,i2)*zb(i6,i3))
     &      ))/
     & (t(i2,i3,i6)*za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i3,i6)*
     &   za(i4,i5)**2*za(i6,i5)**2*zab2(i4,i1,i4,i5)*zab2(i5,i1,i4,i3)*
     &   zab2(i5,i1,i4,i4)*zab2(i5,i1,i4,i5)*zab2(i5,i1,i4,i6)*
     &   zab2(i6,i1,i4,i5)*zb(i3,i2)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i2,i6)*za(i5,i6)**2*zab2(i2,i3,i5,i2)*zab2(i5,i1,i4,i6)*
     &      zab2(i6,i1,i4,i3)*zab2(i6,i1,i4,i5)*zb(i5,i2)*
     &      (za(i2,i3)*(za(i4,i6)*zab2(i1,i1,i4,i6)*
     &            (-(za(i4,i6)*zab2(i1,i1,i4,i6)) - 
     &              za(i1,i6)*zab2(i4,i1,i4,i6))*zab2(i6,i1,i4,i4)*
     &            zb(i6,i2) + 
     &           za(i1,i4)*za(i1,i6)*zab2(i4,i1,i4,i6)*
     &            zab2(i6,i1,i4,i2)*zab2(i6,i1,i4,i6)*zb(i6,i4)) + 
     &        za(i3,i5)*(-(za(i1,i4)*za(i1,i6)*zab2(i4,i1,i4,i6)*
     &              zab2(i6,i1,i4,i5)*zab2(i6,i1,i4,i6)*zb(i6,i4)) - 
     &           za(i4,i6)*zab2(i1,i1,i4,i6)*
     &            (-(za(i4,i6)*zab2(i1,i1,i4,i6)) - 
     &              za(i1,i6)*zab2(i4,i1,i4,i6))*zab2(i6,i1,i4,i4)*
     &            zb(i6,i5)))) + 
     &   t(i2,i3,i5)*(-(za(i2,i6)*za(i4,i6)**2*za(i5,i6)**2*
     &         zab2(i1,i1,i4,i6)**2*zab2(i3,i1,i4,i6)*zab2(i6,i1,i4,i3)*
     &         zab2(i6,i1,i4,i4)*zab2(i6,i1,i4,i5)*
     &         (za(i2,i5)*zb(i5,i2)*zb(i6,i2) + 
     &           za(i3,i5)*zb(i3,i2)*zb(i6,i5))) - 
     &      za(i1,i6)*za(i5,i6)*zab2(i4,i1,i4,i6)*zab2(i5,i1,i4,i6)*
     &       zab2(i6,i1,i4,i5)*
     &       (za(i2,i3)*za(i2,i6)*za(i4,i6)*za(i5,i6)*zab2(i1,i1,i4,i6)*
     &          zab2(i6,i1,i4,i4)*zab2(i6,i1,i4,i5)*zb(i3,i2)*zb(i6,i2)
     &          + za(i3,i6)*(za(i2,i5)*za(i4,i6)*zab2(i6,i1,i4,i4)*
     &             (za(i1,i2)*zab2(i6,i1,i4,i5)*zab2(i6,i1,i4,i6)*
     &                zb(i3,i2) + 
     &               za(i2,i6)*zab2(i1,i1,i4,i6)*zab2(i6,i1,i4,i3)*
     &                zb(i5,i2))*zb(i6,i2) + 
     &            za(i1,i4)*za(i2,i6)*zab2(i6,i1,i4,i3)*
     &             zab2(i6,i1,i4,i6)*
     &             (-(za(i3,i5)*zab2(i6,i1,i4,i5)*zb(i3,i2)) - 
     &               za(i2,i5)*zab2(i6,i1,i4,i2)*zb(i5,i2))*zb(i6,i4))
     &          - za(i2,i6)*za(i3,i5)*za(i4,i6)*za(i5,i6)*
     &          zab2(i1,i1,i4,i6)*zab2(i6,i1,i4,i4)*zab2(i6,i1,i4,i5)*
     &          zb(i3,i2)*zb(i6,i5)) - 
     &      za(i1,i6)**2*za(i4,i6)*zab2(i4,i1,i4,i6)*zab2(i5,i1,i4,i6)*
     &       zab2(i6,i1,i4,i4)*
     &       (za(i2,i3)*za(i5,i6)**2*zab2(i3,i1,i4,i6)*
     &          zab2(i6,i1,i4,i5)*zab2(i6,i1,i4,i6)*zb(i3,i2)*zb(i5,i3)
     &          + za(i2,i5)*(-(zab2(i6,i1,i4,i5)*
     &               (za(i2,i6)*za(i5,i6)*zab2(i3,i1,i4,i6)*
     &                  zab2(i6,i1,i4,i2) + 
     &                 za(i3,i5)*za(i3,i6)*zab2(i6,i1,i4,i3)*
     &                  zab2(i6,i1,i4,i6))*zb(i3,i2)) - 
     &            za(i2,i6)*za(i3,i5)*zab2(i6,i1,i4,i2)*
     &             zab2(i6,i1,i4,i3)*zab2(i6,i1,i4,i6)*zb(i5,i2))*
     &          zb(i6,i5))))/
     & (t(i2,i3,i5)*za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i3,i5)*
     &   za(i4,i6)**2*za(i5,i6)**2*zab2(i4,i1,i4,i6)*zab2(i5,i1,i4,i6)*
     &   zab2(i6,i1,i4,i3)*zab2(i6,i1,i4,i4)*zab2(i6,i1,i4,i5)*
     &   zab2(i6,i1,i4,i6)*zb(i3,i2)*zb(i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i2,i6)*za(i4,i6)**2*zab2(i2,i3,i4,i2)*zab2(i4,i1,i5,i6)*
     &      zab2(i6,i1,i5,i3)*zab2(i6,i1,i5,i4)*zb(i4,i2)*
     &      (za(i2,i3)*(za(i5,i6)*zab2(i1,i1,i5,i6)*
     &            (-(za(i5,i6)*zab2(i1,i1,i5,i6)) - 
     &              za(i1,i6)*zab2(i5,i1,i5,i6))*zab2(i6,i1,i5,i5)*
     &            zb(i6,i2) + 
     &           za(i1,i5)*za(i1,i6)*zab2(i5,i1,i5,i6)*
     &            zab2(i6,i1,i5,i2)*zab2(i6,i1,i5,i6)*zb(i6,i5)) + 
     &        za(i3,i4)*(-(za(i5,i6)*zab2(i1,i1,i5,i6)*
     &              (-(za(i5,i6)*zab2(i1,i1,i5,i6)) - 
     &                za(i1,i6)*zab2(i5,i1,i5,i6))*zab2(i6,i1,i5,i5)*
     &              zb(i6,i4)) - 
     &           za(i1,i5)*za(i1,i6)*zab2(i5,i1,i5,i6)*
     &            zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i6)*zb(i6,i5)))) + 
     &   t(i2,i3,i4)*(-(za(i2,i6)*za(i4,i6)**2*za(i5,i6)**2*
     &         zab2(i1,i1,i5,i6)**2*zab2(i3,i1,i5,i6)*zab2(i6,i1,i5,i3)*
     &         zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i5)*
     &         (za(i2,i4)*zb(i4,i2)*zb(i6,i2) + 
     &           za(i3,i4)*zb(i3,i2)*zb(i6,i4))) - 
     &      za(i1,i6)**2*za(i5,i6)*zab2(i4,i1,i5,i6)*zab2(i5,i1,i5,i6)*
     &       zab2(i6,i1,i5,i5)*
     &       (za(i2,i3)*za(i4,i6)**2*zab2(i3,i1,i5,i6)*
     &          zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i6)*zb(i3,i2)*zb(i4,i3)
     &          + za(i2,i4)*(-(zab2(i6,i1,i5,i4)*
     &               (za(i2,i6)*za(i4,i6)*zab2(i3,i1,i5,i6)*
     &                  zab2(i6,i1,i5,i2) + 
     &                 za(i3,i4)*za(i3,i6)*zab2(i6,i1,i5,i3)*
     &                  zab2(i6,i1,i5,i6))*zb(i3,i2)) - 
     &            za(i2,i6)*za(i3,i4)*zab2(i6,i1,i5,i2)*
     &             zab2(i6,i1,i5,i3)*zab2(i6,i1,i5,i6)*zb(i4,i2))*
     &          zb(i6,i4)) - 
     &      za(i1,i6)*za(i4,i6)*zab2(i4,i1,i5,i6)*zab2(i5,i1,i5,i6)*
     &       zab2(i6,i1,i5,i4)*
     &       (za(i2,i3)*za(i2,i6)*za(i4,i6)*za(i5,i6)*zab2(i1,i1,i5,i6)*
     &          zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i5)*zb(i3,i2)*zb(i6,i2)
     &          - za(i2,i6)*za(i3,i4)*za(i4,i6)*za(i5,i6)*
     &          zab2(i1,i1,i5,i6)*zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i5)*
     &          zb(i3,i2)*zb(i6,i4) + 
     &         za(i3,i6)*(za(i2,i4)*za(i5,i6)*zab2(i6,i1,i5,i5)*
     &             (za(i1,i2)*zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i6)*
     &                zb(i3,i2) + 
     &               za(i2,i6)*zab2(i1,i1,i5,i6)*zab2(i6,i1,i5,i3)*
     &                zb(i4,i2))*zb(i6,i2) + 
     &            za(i1,i5)*za(i2,i6)*zab2(i6,i1,i5,i3)*
     &             zab2(i6,i1,i5,i6)*
     &             (-(za(i3,i4)*zab2(i6,i1,i5,i4)*zb(i3,i2)) - 
     &               za(i2,i4)*zab2(i6,i1,i5,i2)*zb(i4,i2))*zb(i6,i5))))
     &   )/(t(i2,i3,i4)*za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i3,i4)*
     &   za(i4,i6)**2*za(i5,i6)**2*zab2(i4,i1,i5,i6)*zab2(i5,i1,i5,i6)*
     &   zab2(i6,i1,i5,i3)*zab2(i6,i1,i5,i4)*zab2(i6,i1,i5,i5)*
     &   zab2(i6,i1,i5,i6)*zb(i3,i2)*zb(i4,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i2,i4)*za(i6,i4)**2*zab2(i2,i3,i6,i2)*zab2(i4,i1,i5,i3)*
     &      zab2(i4,i1,i5,i6)*zab2(i6,i1,i5,i4)*
     &      (za(i2,i3)*(za(i5,i4)*zab2(i1,i1,i5,i4)*zab2(i4,i1,i5,i5)*
     &            (-(za(i5,i4)*zab2(i1,i1,i5,i4)) - 
     &              za(i1,i4)*zab2(i5,i1,i5,i4))*zb(i4,i2) + 
     &           za(i1,i4)*za(i1,i5)*zab2(i4,i1,i5,i2)*
     &            zab2(i4,i1,i5,i4)*zab2(i5,i1,i5,i4)*zb(i4,i5)) + 
     &        za(i3,i6)*(-(za(i1,i4)*za(i1,i5)*zab2(i4,i1,i5,i4)*
     &              zab2(i4,i1,i5,i6)*zab2(i5,i1,i5,i4)*zb(i4,i5)) - 
     &           za(i5,i4)*zab2(i1,i1,i5,i4)*zab2(i4,i1,i5,i5)*
     &            (-(za(i5,i4)*zab2(i1,i1,i5,i4)) - 
     &              za(i1,i4)*zab2(i5,i1,i5,i4))*zb(i4,i6)))*zb(i6,i2))
     &    + t(i2,i3,i6)*(-(za(i2,i4)*za(i5,i4)**2*za(i6,i4)**2*
     &         zab2(i1,i1,i5,i4)**2*zab2(i3,i1,i5,i4)*zab2(i4,i1,i5,i3)*
     &         zab2(i4,i1,i5,i5)*zab2(i4,i1,i5,i6)*
     &         (za(i3,i6)*zb(i3,i2)*zb(i4,i6) + 
     &           za(i2,i6)*zb(i4,i2)*zb(i6,i2))) - 
     &      za(i1,i4)*za(i6,i4)*zab2(i4,i1,i5,i6)*zab2(i5,i1,i5,i4)*
     &       zab2(i6,i1,i5,i4)*
     &       (za(i2,i3)*za(i2,i4)*za(i5,i4)*za(i6,i4)*zab2(i1,i1,i5,i4)*
     &          zab2(i4,i1,i5,i5)*zab2(i4,i1,i5,i6)*zb(i3,i2)*zb(i4,i2)
     &          - za(i2,i4)*za(i3,i6)*za(i5,i4)*za(i6,i4)*
     &          zab2(i1,i1,i5,i4)*zab2(i4,i1,i5,i5)*zab2(i4,i1,i5,i6)*
     &          zb(i3,i2)*zb(i4,i6) + 
     &         za(i3,i4)*(za(i1,i5)*za(i2,i4)*zab2(i4,i1,i5,i3)*
     &             zab2(i4,i1,i5,i4)*zb(i4,i5)*
     &             (-(za(i3,i6)*zab2(i4,i1,i5,i6)*zb(i3,i2)) - 
     &               za(i2,i6)*zab2(i4,i1,i5,i2)*zb(i6,i2)) + 
     &            za(i2,i6)*za(i5,i4)*zab2(i4,i1,i5,i5)*zb(i4,i2)*
     &             (za(i1,i2)*zab2(i4,i1,i5,i4)*zab2(i4,i1,i5,i6)*
     &                zb(i3,i2) + 
     &               za(i2,i4)*zab2(i1,i1,i5,i4)*zab2(i4,i1,i5,i3)*
     &                zb(i6,i2)))) - 
     &      za(i1,i4)**2*za(i5,i4)*zab2(i4,i1,i5,i5)*zab2(i5,i1,i5,i4)*
     &       zab2(i6,i1,i5,i4)*
     &       (za(i2,i6)*zb(i4,i6)*
     &          (-((za(i2,i4)*za(i6,i4)*zab2(i3,i1,i5,i4)*
     &                  zab2(i4,i1,i5,i2) + 
     &                 za(i3,i4)*za(i3,i6)*zab2(i4,i1,i5,i3)*
     &                  zab2(i4,i1,i5,i4))*zab2(i4,i1,i5,i6)*zb(i3,i2))
     &             - za(i2,i4)*za(i3,i6)*zab2(i4,i1,i5,i2)*
     &             zab2(i4,i1,i5,i3)*zab2(i4,i1,i5,i4)*zb(i6,i2)) + 
     &         za(i2,i3)*za(i6,i4)**2*zab2(i3,i1,i5,i4)*
     &          zab2(i4,i1,i5,i4)*zab2(i4,i1,i5,i6)*zb(i3,i2)*zb(i6,i3))
     &      ))/
     & (t(i2,i3,i6)*za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i3,i6)*
     &   za(i5,i4)**2*za(i6,i4)**2*zab2(i4,i1,i5,i3)*zab2(i4,i1,i5,i4)*
     &   zab2(i4,i1,i5,i5)*zab2(i4,i1,i5,i6)*zab2(i5,i1,i5,i4)*
     &   zab2(i6,i1,i5,i4)*zb(i3,i2)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i2,i4)*za(i5,i4)**2*zab2(i2,i3,i5,i2)*zab2(i4,i1,i6,i3)*
     &      zab2(i4,i1,i6,i5)*zab2(i5,i1,i6,i4)*
     &      (za(i2,i3)*(za(i6,i4)*zab2(i1,i1,i6,i4)*zab2(i4,i1,i6,i6)*
     &            (-(za(i6,i4)*zab2(i1,i1,i6,i4)) - 
     &              za(i1,i4)*zab2(i6,i1,i6,i4))*zb(i4,i2) + 
     &           za(i1,i4)*za(i1,i6)*zab2(i4,i1,i6,i2)*
     &            zab2(i4,i1,i6,i4)*zab2(i6,i1,i6,i4)*zb(i4,i6)) + 
     &        za(i3,i5)*(-(za(i6,i4)*zab2(i1,i1,i6,i4)*
     &              zab2(i4,i1,i6,i6)*
     &              (-(za(i6,i4)*zab2(i1,i1,i6,i4)) - 
     &                za(i1,i4)*zab2(i6,i1,i6,i4))*zb(i4,i5)) - 
     &           za(i1,i4)*za(i1,i6)*zab2(i4,i1,i6,i4)*
     &            zab2(i4,i1,i6,i5)*zab2(i6,i1,i6,i4)*zb(i4,i6)))*
     &      zb(i5,i2)) + t(i2,i3,i5)*
     &    (-(za(i2,i4)*za(i5,i4)**2*za(i6,i4)**2*zab2(i1,i1,i6,i4)**2*
     &         zab2(i3,i1,i6,i4)*zab2(i4,i1,i6,i3)*zab2(i4,i1,i6,i5)*
     &         zab2(i4,i1,i6,i6)*
     &         (za(i3,i5)*zb(i3,i2)*zb(i4,i5) + 
     &           za(i2,i5)*zb(i4,i2)*zb(i5,i2))) - 
     &      za(i1,i4)*za(i5,i4)*zab2(i4,i1,i6,i5)*zab2(i5,i1,i6,i4)*
     &       zab2(i6,i1,i6,i4)*
     &       (za(i2,i3)*za(i2,i4)*za(i5,i4)*za(i6,i4)*zab2(i1,i1,i6,i4)*
     &          zab2(i4,i1,i6,i5)*zab2(i4,i1,i6,i6)*zb(i3,i2)*zb(i4,i2)
     &          - za(i2,i4)*za(i3,i5)*za(i5,i4)*za(i6,i4)*
     &          zab2(i1,i1,i6,i4)*zab2(i4,i1,i6,i5)*zab2(i4,i1,i6,i6)*
     &          zb(i3,i2)*zb(i4,i5) + 
     &         za(i3,i4)*(za(i1,i6)*za(i2,i4)*zab2(i4,i1,i6,i3)*
     &             zab2(i4,i1,i6,i4)*zb(i4,i6)*
     &             (-(za(i3,i5)*zab2(i4,i1,i6,i5)*zb(i3,i2)) - 
     &               za(i2,i5)*zab2(i4,i1,i6,i2)*zb(i5,i2)) + 
     &            za(i2,i5)*za(i6,i4)*zab2(i4,i1,i6,i6)*zb(i4,i2)*
     &             (za(i1,i2)*zab2(i4,i1,i6,i4)*zab2(i4,i1,i6,i5)*
     &                zb(i3,i2) + 
     &               za(i2,i4)*zab2(i1,i1,i6,i4)*zab2(i4,i1,i6,i3)*
     &                zb(i5,i2)))) - 
     &      za(i1,i4)**2*za(i6,i4)*zab2(i4,i1,i6,i6)*zab2(i5,i1,i6,i4)*
     &       zab2(i6,i1,i6,i4)*
     &       (za(i2,i5)*zb(i4,i5)*
     &          (-((za(i2,i4)*za(i5,i4)*zab2(i3,i1,i6,i4)*
     &                  zab2(i4,i1,i6,i2) + 
     &                 za(i3,i4)*za(i3,i5)*zab2(i4,i1,i6,i3)*
     &                  zab2(i4,i1,i6,i4))*zab2(i4,i1,i6,i5)*zb(i3,i2))
     &             - za(i2,i4)*za(i3,i5)*zab2(i4,i1,i6,i2)*
     &             zab2(i4,i1,i6,i3)*zab2(i4,i1,i6,i4)*zb(i5,i2)) + 
     &         za(i2,i3)*za(i5,i4)**2*zab2(i3,i1,i6,i4)*
     &          zab2(i4,i1,i6,i4)*zab2(i4,i1,i6,i5)*zb(i3,i2)*zb(i5,i3))
     &      ))/
     & (t(i2,i3,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i3,i5)*
     &   za(i5,i4)**2*za(i6,i4)**2*zab2(i4,i1,i6,i3)*zab2(i4,i1,i6,i4)*
     &   zab2(i4,i1,i6,i5)*zab2(i4,i1,i6,i6)*zab2(i5,i1,i6,i4)*
     &   zab2(i6,i1,i6,i4)*zb(i3,i2)*zb(i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i2,i5)*za(i4,i5)**2*zab2(i2,i3,i4,i2)*zab2(i4,i1,i6,i5)*
     &      zab2(i5,i1,i6,i3)*zab2(i5,i1,i6,i4)*zb(i4,i2)*
     &      (za(i2,i3)*(za(i6,i5)*zab2(i1,i1,i6,i5)*zab2(i5,i1,i6,i6)*
     &            (-(za(i6,i5)*zab2(i1,i1,i6,i5)) - 
     &              za(i1,i5)*zab2(i6,i1,i6,i5))*zb(i5,i2) + 
     &           za(i1,i5)*za(i1,i6)*zab2(i5,i1,i6,i2)*
     &            zab2(i5,i1,i6,i5)*zab2(i6,i1,i6,i5)*zb(i5,i6)) + 
     &        za(i3,i4)*(-(za(i6,i5)*zab2(i1,i1,i6,i5)*
     &              zab2(i5,i1,i6,i6)*
     &              (-(za(i6,i5)*zab2(i1,i1,i6,i5)) - 
     &                za(i1,i5)*zab2(i6,i1,i6,i5))*zb(i5,i4)) - 
     &           za(i1,i5)*za(i1,i6)*zab2(i5,i1,i6,i4)*
     &            zab2(i5,i1,i6,i5)*zab2(i6,i1,i6,i5)*zb(i5,i6)))) + 
     &   t(i2,i3,i4)*(-(za(i2,i5)*za(i4,i5)**2*za(i6,i5)**2*
     &         zab2(i1,i1,i6,i5)**2*zab2(i3,i1,i6,i5)*zab2(i5,i1,i6,i3)*
     &         zab2(i5,i1,i6,i4)*zab2(i5,i1,i6,i6)*
     &         (za(i2,i4)*zb(i4,i2)*zb(i5,i2) + 
     &           za(i3,i4)*zb(i3,i2)*zb(i5,i4))) - 
     &      za(i1,i5)**2*za(i6,i5)*zab2(i4,i1,i6,i5)*zab2(i5,i1,i6,i6)*
     &       zab2(i6,i1,i6,i5)*
     &       (za(i2,i3)*za(i4,i5)**2*zab2(i3,i1,i6,i5)*
     &          zab2(i5,i1,i6,i4)*zab2(i5,i1,i6,i5)*zb(i3,i2)*zb(i4,i3)
     &          + za(i2,i4)*(-(zab2(i5,i1,i6,i4)*
     &               (za(i2,i5)*za(i4,i5)*zab2(i3,i1,i6,i5)*
     &                  zab2(i5,i1,i6,i2) + 
     &                 za(i3,i4)*za(i3,i5)*zab2(i5,i1,i6,i3)*
     &                  zab2(i5,i1,i6,i5))*zb(i3,i2)) - 
     &            za(i2,i5)*za(i3,i4)*zab2(i5,i1,i6,i2)*
     &             zab2(i5,i1,i6,i3)*zab2(i5,i1,i6,i5)*zb(i4,i2))*
     &          zb(i5,i4)) - 
     &      za(i1,i5)*za(i4,i5)*zab2(i4,i1,i6,i5)*zab2(i5,i1,i6,i4)*
     &       zab2(i6,i1,i6,i5)*
     &       (za(i2,i3)*za(i2,i5)*za(i4,i5)*za(i6,i5)*zab2(i1,i1,i6,i5)*
     &          zab2(i5,i1,i6,i4)*zab2(i5,i1,i6,i6)*zb(i3,i2)*zb(i5,i2)
     &          - za(i2,i5)*za(i3,i4)*za(i4,i5)*za(i6,i5)*
     &          zab2(i1,i1,i6,i5)*zab2(i5,i1,i6,i4)*zab2(i5,i1,i6,i6)*
     &          zb(i3,i2)*zb(i5,i4) + 
     &         za(i3,i5)*(za(i2,i4)*za(i6,i5)*zab2(i5,i1,i6,i6)*
     &             (za(i1,i2)*zab2(i5,i1,i6,i4)*zab2(i5,i1,i6,i5)*
     &                zb(i3,i2) + 
     &               za(i2,i5)*zab2(i1,i1,i6,i5)*zab2(i5,i1,i6,i3)*
     &                zb(i4,i2))*zb(i5,i2) + 
     &            za(i1,i6)*za(i2,i5)*zab2(i5,i1,i6,i3)*
     &             zab2(i5,i1,i6,i5)*
     &             (-(za(i3,i4)*zab2(i5,i1,i6,i4)*zb(i3,i2)) - 
     &               za(i2,i4)*zab2(i5,i1,i6,i2)*zb(i4,i2))*zb(i5,i6))))
     &   )/(t(i2,i3,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i3,i4)*
     &   za(i4,i5)**2*za(i6,i5)**2*zab2(i4,i1,i6,i5)*zab2(i5,i1,i6,i3)*
     &   zab2(i5,i1,i6,i4)*zab2(i5,i1,i6,i5)*zab2(i5,i1,i6,i6)*
     &   zab2(i6,i1,i6,i5)*zb(i3,i2)*zb(i4,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (t(i2,i5,i6)*za(i1,i2)*zab2(i1,i1,i4,i3)*zab2(i3,i1,i4,i3))/
     & (za(i1,i4)*za(i2,i5)*za(i2,i6)*zab2(i4,i1,i4,i3)*
     &   zab2(i5,i1,i4,i3)*zab2(i6,i1,i4,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (t(i2,i4,i6)*za(i1,i2)*zab2(i1,i1,i5,i3)*zab2(i3,i1,i5,i3))/
     & (za(i1,i5)*za(i2,i4)*za(i2,i6)*zab2(i4,i1,i5,i3)*
     &   zab2(i5,i1,i5,i3)*zab2(i6,i1,i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (t(i2,i4,i5)*za(i1,i2)*zab2(i1,i1,i6,i3)*zab2(i3,i1,i6,i3))/
     & (za(i1,i6)*za(i2,i4)*za(i2,i5)*zab2(i4,i1,i6,i3)*
     &   zab2(i5,i1,i6,i3)*zab2(i6,i1,i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i1,i6)*za(i2,i6)*za(i3,i6)**2*za(i6,i4)*zab2(i4,i2,i3,i6)*
     &      zab2(i5,i2,i3,i6)*zab2(i6,i2,i3,i4)*
     &      (za(i6,i5)**2*zab2(i1,i2,i3,i6)*zab2(i6,i2,i3,i5) + 
     &        za(i1,i5)*za(i6,i4)*zab2(i6,i2,i3,i4)*zab2(i6,i2,i3,i6))*
     &      zb(i3,i2)*zb(i5,i6)*zb(i6,i3)) - 
     &   za(i6,i5)*zab2(i6,i2,i3,i5)*
     &    (za(i1,i6)*za(i2,i6)*za(i3,i6)**2*zab2(i4,i2,i3,i6)*
     &       zab2(i5,i2,i3,i6)*
     &       (za(i6,i4)**2*zab2(i1,i2,i3,i6)*zab2(i6,i2,i3,i4) + 
     &         za(i1,i4)*za(i6,i5)*zab2(i6,i2,i3,i5)*zab2(i6,i2,i3,i6))*
     &       zb(i3,i2)*zb(i4,i6)*zb(i6,i3) - 
     &      za(i6,i4)*zab2(i6,i2,i3,i4)*
     &       (za(i2,i6)**2*za(i6,i4)*za(i6,i5)*zab2(i1,i2,i3,i6)**3*
     &          zab2(i6,i2,i3,i1)*zb(i6,i2)**2 - 
     &         za(i1,i6)**2*zab2(i4,i2,i3,i6)*zab2(i5,i2,i3,i6)*
     &          zab2(i6,i2,i3,i6)*
     &          (-(za(i2,i6)*za(i3,i6)*zab2(i1,i2,i3,i6)*zb(i2,i1)) + 
     &            za(i1,i6)*za(i2,i3)*zab2(i6,i2,i3,i1)*zb(i6,i2))*
     &          zb(i6,i3))))/
     & (za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i6,i4)**2*za(i6,i5)**2*
     &   zab2(i4,i2,i3,i6)*zab2(i5,i2,i3,i6)*zab2(i6,i2,i3,i1)*
     &   zab2(i6,i2,i3,i4)*zab2(i6,i2,i3,i5)*zab2(i6,i2,i3,i6)*
     &   zb(i3,i2)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i5,i6)*zab2(i5,i2,i3,i6)*
     &      (za(i1,i5)*za(i2,i5)*za(i3,i5)**2*zab2(i4,i2,i3,i5)*
     &         (za(i5,i4)**2*zab2(i1,i2,i3,i5)*zab2(i5,i2,i3,i4) + 
     &           za(i1,i4)*za(i5,i6)*zab2(i5,i2,i3,i5)*zab2(i5,i2,i3,i6)
     &           )*zab2(i6,i2,i3,i5)*zb(i3,i2)*zb(i4,i5)*zb(i5,i3) - 
     &        za(i5,i4)*zab2(i5,i2,i3,i4)*
     &         (za(i2,i5)**2*za(i5,i4)*za(i5,i6)*zab2(i1,i2,i3,i5)**3*
     &            zab2(i5,i2,i3,i1)*zb(i5,i2)**2 - 
     &           za(i1,i5)**2*zab2(i4,i2,i3,i5)*zab2(i5,i2,i3,i5)*
     &            zab2(i6,i2,i3,i5)*
     &            (-(za(i2,i5)*za(i3,i5)*zab2(i1,i2,i3,i5)*zb(i2,i1)) + 
     &              za(i1,i5)*za(i2,i3)*zab2(i5,i2,i3,i1)*zb(i5,i2))*
     &            zb(i5,i3)))) - 
     &   za(i1,i5)*za(i2,i5)*za(i3,i5)**2*za(i5,i4)*zab2(i4,i2,i3,i5)*
     &    zab2(i5,i2,i3,i4)*(za(i1,i6)*za(i5,i4)*zab2(i5,i2,i3,i4)*
     &       zab2(i5,i2,i3,i5) + 
     &      za(i5,i6)**2*zab2(i1,i2,i3,i5)*zab2(i5,i2,i3,i6))*
     &    zab2(i6,i2,i3,i5)*zb(i3,i2)*zb(i5,i3)*zb(i6,i5))/
     & (za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i5,i4)**2*za(i5,i6)**2*
     &   zab2(i4,i2,i3,i5)*zab2(i5,i2,i3,i1)*zab2(i5,i2,i3,i4)*
     &   zab2(i5,i2,i3,i5)*zab2(i5,i2,i3,i6)*zab2(i6,i2,i3,i5)*
     &   zb(i3,i2)*zb(i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-(za(i4,i6)*zab2(i4,i2,i3,i6)*
     &      (-(za(i4,i5)*zab2(i4,i2,i3,i5)*
     &           (za(i2,i4)**2*za(i4,i5)*za(i4,i6)*zab2(i1,i2,i3,i4)**3*
     &              zab2(i4,i2,i3,i1)*zb(i4,i2)**2 - 
     &             za(i1,i4)**2*zab2(i4,i2,i3,i4)*zab2(i5,i2,i3,i4)*
     &              zab2(i6,i2,i3,i4)*
     &              (-(za(i2,i4)*za(i3,i4)*zab2(i1,i2,i3,i4)*
     &                   zb(i2,i1)) + 
     &                za(i1,i4)*za(i2,i3)*zab2(i4,i2,i3,i1)*zb(i4,i2))*
     &              zb(i4,i3))) + 
     &        za(i1,i4)*za(i2,i4)*za(i3,i4)**2*
     &         (za(i4,i5)**2*zab2(i1,i2,i3,i4)*zab2(i4,i2,i3,i5) + 
     &           za(i1,i5)*za(i4,i6)*zab2(i4,i2,i3,i4)*zab2(i4,i2,i3,i6)
     &           )*zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)*zb(i3,i2)*
     &         zb(i4,i3)*zb(i5,i4))) - 
     &   za(i1,i4)*za(i2,i4)*za(i3,i4)**2*za(i4,i5)*zab2(i4,i2,i3,i5)*
     &    (za(i1,i6)*za(i4,i5)*zab2(i4,i2,i3,i4)*zab2(i4,i2,i3,i5) + 
     &      za(i4,i6)**2*zab2(i1,i2,i3,i4)*zab2(i4,i2,i3,i6))*
     &    zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)*zb(i3,i2)*zb(i4,i3)*
     &    zb(i6,i4))/
     & (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i4,i5)**2*za(i4,i6)**2*
     &   zab2(i4,i2,i3,i1)*zab2(i4,i2,i3,i4)*zab2(i4,i2,i3,i5)*
     &   zab2(i4,i2,i3,i6)*zab2(i5,i2,i3,i4)*zab2(i6,i2,i3,i4)*
     &   zb(i3,i2)*zb(i4,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i6)*zab2(i1,i2,i4,i6)*zab2(i6,i2,i4,i5))/
     &    (za(i1,i5)*za(i2,i4)*za(i4,i6)*za(i5,i6)*zab2(i6,i2,i4,i3)*
     &      zb(i5,i3))) + (zab2(i6,i2,i4,i5)**2*zb(i6,i4))/
     &  (t(i1,i3,i5)*za(i2,i6)*za(i4,i6)*zab2(i6,i2,i4,i3)*zb(i3,i1)) + 
     & (za(i1,i3)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i4)*za(i4,i6)*zb(i3,i1)) + 
     & (t(i1,i3,i5)*za(i1,i6)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i4)*za(i4,i6)*za(i5,i6)*zb(i3,i1)*zb(i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (zab2(i5,i2,i4,i6)**2*zb(i5,i4))/
     &  (t(i1,i3,i6)*za(i2,i5)*za(i4,i5)*zab2(i5,i2,i4,i3)*zb(i3,i1)) + 
     & (za(i1,i3)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i4)*za(i4,i5)*zb(i3,i1)) - 
     & (za(i1,i5)*zab2(i1,i2,i4,i5)*zab2(i5,i2,i4,i6))/
     &  (za(i1,i6)*za(i2,i4)*za(i4,i5)*za(i6,i5)*zab2(i5,i2,i4,i3)*
     &    zb(i6,i3)) + (t(i1,i3,i6)*za(i1,i5)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i4)*za(i4,i5)*za(i6,i5)*zb(i3,i1)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i6)*zab2(i1,i2,i5,i6)*zab2(i6,i2,i5,i4))/
     &    (za(i1,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6)*zab2(i6,i2,i5,i3)*
     &      zb(i4,i3))) + (za(i1,i3)*zb(i6,i4))/
     &  (za(i1,i4)*za(i2,i5)*za(i5,i6)*zb(i3,i1)) + 
     & (t(i1,i3,i4)*za(i1,i6)*zb(i6,i4))/
     &  (za(i1,i4)*za(i2,i5)*za(i4,i6)*za(i5,i6)*zb(i3,i1)*zb(i4,i3)) + 
     & (zab2(i6,i2,i5,i4)**2*zb(i6,i5))/
     &  (t(i1,i3,i4)*za(i2,i6)*za(i5,i6)*zab2(i6,i2,i5,i3)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (zab2(i4,i2,i5,i6)**2*zb(i4,i5))/
     &  (t(i1,i3,i6)*za(i2,i4)*za(i5,i4)*zab2(i4,i2,i5,i3)*zb(i3,i1)) + 
     & (za(i1,i3)*zb(i4,i6))/
     &  (za(i1,i6)*za(i2,i5)*za(i5,i4)*zb(i3,i1)) - 
     & (za(i1,i4)*zab2(i1,i2,i5,i4)*zab2(i4,i2,i5,i6))/
     &  (za(i1,i6)*za(i2,i5)*za(i5,i4)*za(i6,i4)*zab2(i4,i2,i5,i3)*
     &    zb(i6,i3)) + (t(i1,i3,i6)*za(i1,i4)*zb(i4,i6))/
     &  (za(i1,i6)*za(i2,i5)*za(i5,i4)*za(i6,i4)*zb(i3,i1)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i5)*zab2(i1,i2,i6,i5)*zab2(i5,i2,i6,i4))/
     &    (za(i1,i4)*za(i2,i6)*za(i4,i5)*za(i6,i5)*zab2(i5,i2,i6,i3)*
     &      zb(i4,i3))) + (za(i1,i3)*zb(i5,i4))/
     &  (za(i1,i4)*za(i2,i6)*za(i6,i5)*zb(i3,i1)) + 
     & (t(i1,i3,i4)*za(i1,i5)*zb(i5,i4))/
     &  (za(i1,i4)*za(i2,i6)*za(i4,i5)*za(i6,i5)*zb(i3,i1)*zb(i4,i3)) + 
     & (zab2(i5,i2,i6,i4)**2*zb(i5,i6))/
     &  (t(i1,i3,i4)*za(i2,i5)*za(i6,i5)*zab2(i5,i2,i6,i3)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i3)*zb(i4,i5))/
     &  (za(i1,i5)*za(i2,i6)*za(i6,i4)*zb(i3,i1)) + 
     & (zab2(i4,i2,i6,i5)**2*zb(i4,i6))/
     &  (t(i1,i3,i5)*za(i2,i4)*za(i6,i4)*zab2(i4,i2,i6,i3)*zb(i3,i1)) - 
     & (za(i1,i4)*zab2(i1,i2,i6,i4)*zab2(i4,i2,i6,i5))/
     &  (za(i1,i5)*za(i2,i6)*za(i5,i4)*za(i6,i4)*zab2(i4,i2,i6,i3)*
     &    zb(i5,i3)) + (t(i1,i3,i5)*za(i1,i4)*zb(i4,i5))/
     &  (za(i1,i5)*za(i2,i6)*za(i5,i4)*za(i6,i4)*zb(i3,i1)*zb(i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i2)*za(i1,i6)*za(i2,i3)*zab2(i3,i2,i4,i3)*
     &      zab2(i5,i2,i4,i3)*zb(i4,i2)*zb(i5,i1))/
     &    (s(i1,i5)*s(i2,i3)*za(i2,i4)**2*za(i5,i6)**2*
     &      zab2(i6,i2,i4,i3)*zb(i3,i4))) + 
     & (za(i1,i2)*za(i1,i6)*za(i2,i3)*zab2(i3,i2,i4,i3)*zb(i4,i2)*
     &    zb(i6,i1))/
     &  (s(i1,i6)*s(i2,i3)*za(i2,i4)**2*za(i5,i6)**2*zb(i3,i4)) + 
     & (za(i1,i2)*za(i2,i3)*zab2(i1,i2,i4,i3)*zab2(i3,i2,i4,i3)*
     &    zb(i4,i2)*zb(i6,i1))/
     &  (s(i1,i6)*s(i2,i3)*za(i2,i4)**2*za(i5,i6)*zab2(i5,i2,i4,i3)*
     &    zb(i3,i4))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i2)*za(i1,i6)*za(i2,i3)*zab2(i3,i2,i5,i3)*
     &      zab2(i4,i2,i5,i3)*zb(i4,i1)*zb(i5,i2))/
     &    (s(i1,i4)*s(i2,i3)*za(i2,i5)**2*za(i4,i6)**2*
     &      zab2(i6,i2,i5,i3)*zb(i3,i5))) + 
     & (za(i1,i2)*za(i1,i6)*za(i2,i3)*zab2(i3,i2,i5,i3)*zb(i5,i2)*
     &    zb(i6,i1))/
     &  (s(i1,i6)*s(i2,i3)*za(i2,i5)**2*za(i4,i6)**2*zb(i3,i5)) + 
     & (za(i1,i2)*za(i2,i3)*zab2(i1,i2,i5,i3)*zab2(i3,i2,i5,i3)*
     &    zb(i5,i2)*zb(i6,i1))/
     &  (s(i1,i6)*s(i2,i3)*za(i2,i5)**2*za(i4,i6)*zab2(i4,i2,i5,i3)*
     &    zb(i3,i5))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i2)*za(i1,i4)*za(i2,i3)*zab2(i3,i2,i6,i3)*zb(i4,i1)*
     &    zb(i6,i2))/
     &  (s(i1,i4)*s(i2,i3)*za(i2,i6)**2*za(i5,i4)**2*zb(i3,i6)) + 
     & (za(i1,i2)*za(i2,i3)*zab2(i1,i2,i6,i3)*zab2(i3,i2,i6,i3)*
     &    zb(i4,i1)*zb(i6,i2))/
     &  (s(i1,i4)*s(i2,i3)*za(i2,i6)**2*za(i5,i4)*zab2(i5,i2,i6,i3)*
     &    zb(i3,i6)) - (za(i1,i2)*za(i1,i4)*za(i2,i3)*zab2(i3,i2,i6,i3)*
     &    zab2(i5,i2,i6,i3)*zb(i5,i1)*zb(i6,i2))/
     &  (s(i1,i5)*s(i2,i3)*za(i2,i6)**2*za(i5,i4)**2*zab2(i4,i2,i6,i3)*
     &    zb(i3,i6))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (2*za(i1,i3)**2*za(i3,i5)*zb(i6,i5))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i2,i5)*za(i3,i4)*za(i5,i6)) + 
     & (2*za(i1,i3)**2*za(i3,i6)*zb(i6,i5))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i2,i6)*za(i3,i4)*za(i5,i6)) + 
     & (2*za(i1,i3)**2*za(i1,i5)*zb(i4,i1)*zb(i6,i5))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i2,i5)*za(i3,i4)*za(i5,i6)*zb(i4,i3))
     &  + (2*za(i1,i3)**2*za(i1,i6)*zb(i4,i1)*zb(i6,i5))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i2,i6)*za(i3,i4)*za(i5,i6)*zb(i4,i3))
     &  + (2*za(i3,i6)*zb(i4,i2)*zb(i5,i2)*zb(i6,i5))/
     &  (za(i3,i4)*za(i4,i6)*za(i5,i6)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1)) - 
     & (2*za(i2,i3)*zb(i4,i2)**2*zb(i5,i2)*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i3,i4)*za(i5,i6)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1))
     &  + (2*za(i1,i6)**2*zb(i5,i2)**2*zb(i6,i5))/
     &  (za(i1,i4)*za(i4,i6)*za(i5,i6)*zab2(i6,i2,i3,i5)*zb(i3,i2)*
     &    zb(i5,i3)) + (2*za(i3,i6)**2*zb(i5,i4)*zb(i6,i5))/
     &  (za(i2,i6)*za(i3,i4)*za(i4,i6)*za(i5,i6)*zb(i4,i3)*zb(i5,i1)) + 
     & (2*zb(i4,i2)**2*zb(i5,i4)*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i5,i6)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1)) + 
     & (2*za(i1,i6)**2*zb(i5,i4)*zb(i6,i5))/
     &  (za(i1,i4)*za(i2,i6)*za(i4,i6)*za(i5,i6)*zb(i4,i3)*zb(i5,i3)) + 
     & (2*za(i3,i5)*zb(i4,i2)*zb(i6,i2)*zb(i6,i5))/
     &  (za(i3,i4)*za(i4,i5)*za(i5,i6)*zb(i3,i2)*zb(i4,i3)*zb(i6,i1)) - 
     & (2*za(i2,i3)*zb(i4,i2)**2*zb(i6,i2)*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i3,i4)*za(i5,i6)*zb(i3,i2)*zb(i4,i3)*zb(i6,i1))
     &  + (2*za(i1,i5)**2*zb(i6,i2)**2*zb(i6,i5))/
     &  (za(i1,i4)*za(i4,i5)*za(i5,i6)*zab2(i5,i2,i3,i6)*zb(i3,i2)*
     &    zb(i6,i3)) + (2*za(i3,i5)**2*zb(i6,i4)*zb(i6,i5))/
     &  (za(i2,i5)*za(i3,i4)*za(i4,i5)*za(i5,i6)*zb(i4,i3)*zb(i6,i1)) + 
     & (2*zb(i4,i2)**2*zb(i6,i4)*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i5,i6)*zb(i3,i2)*zb(i4,i3)*zb(i6,i1)) + 
     & (2*za(i1,i5)**2*zb(i6,i4)*zb(i6,i5))/
     &  (za(i1,i4)*za(i2,i5)*za(i4,i5)*za(i5,i6)*zb(i4,i3)*zb(i6,i3)) + 
     & (2*(t(i2,i3,i4)*za(i1,i3)*za(i2,i4)*za(i4,i5)*za(i4,i6)*
     &       zb(i4,i1)*(-((za(i2,i6)*za(i3,i5) + za(i2,i5)*za(i3,i6))*
     &            zb(i3,i4)) + 
     &         (za(i1,i6)*za(i2,i5) + za(i1,i5)*za(i2,i6))*zb(i4,i1))*
     &       zb(i5,i1)*zb(i5,i3)*zb(i6,i1)*zb(i6,i3) + 
     &      t(i1,i3,i4)*(-(za(i2,i3)*za(i2,i5)*za(i2,i6)*za(i4,i5)*
     &            za(i4,i6)*zb(i3,i1)*zb(i4,i2)*zb(i5,i3)*zb(i6,i3)*
     &            (za(i2,i3)*
     &               (zb(i5,i2)*zb(i6,i1) + zb(i5,i1)*zb(i6,i2)) + 
     &              za(i4,i3)*
     &               (zb(i5,i4)*zb(i6,i1) + zb(i5,i1)*zb(i6,i4)))) + 
     &         t(i2,i3,i4)*(za(i2,i5)*za(i2,i6)*za(i4,i3)*za(i4,i6)*
     &             zb(i5,i1)*zb(i5,i3)*
     &             (za(i1,i5)*zb(i3,i4)*zb(i6,i1)**2 + 
     &               zab2(i5,i2,i4,i6)*zb(i3,i1)*zb(i6,i4)) + 
     &            zb(i6,i3)*((-(zb(i5,i3)*
     &                     ((-(za(i2,i3)*
     &                        (za(i2,i6)*za(i3,i5) + 
     &                       za(i2,i5)*za(i3,i6))*za(i4,i5)*za(i4,i6)*
     &                        zb(i3,i4)) + 
     &                        za(i2,i4)*
     &                        (za(i1,i6)*za(i2,i5)*za(i3,i6)*
     &                       za(i4,i5) + 
     &                        za(i1,i5)*za(i2,i6)*za(i3,i5)*za(i4,i6))*
     &                        zb(i4,i1))*zb(i5,i1) + 
     &                       za(i2,i4)*za(i2,i5)*za(i3,i6)**2*za(i4,i5)*
     &                        zb(i3,i1)*zb(i5,i4))) + 
     &                  za(i2,i5)*za(i2,i6)*za(i4,i3)*za(i4,i5)*
     &                   (za(i1,i6)*zb(i3,i4)*zb(i5,i1)**2 + 
     &                     zab2(i6,i2,i4,i5)*zb(i3,i1)*zb(i5,i4)))*
     &                zb(i6,i1) - 
     &               za(i2,i4)*za(i2,i6)*za(i3,i5)**2*za(i4,i6)*
     &                zb(i3,i1)*zb(i5,i1)*zb(i5,i3)*zb(i6,i4)))))*
     &    zb(i6,i5))/
     &  (t(i1,i3,i4)*t(i2,i3,i4)*za(i2,i4)*za(i2,i5)*za(i2,i6)*
     &    za(i4,i3)*za(i4,i5)*za(i4,i6)*za(i5,i6)*zb(i3,i1)*zb(i3,i4)*
     &    zb(i5,i1)*zb(i5,i3)*zb(i6,i1)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (2*za(i1,i3)**2*za(i3,i4)*zb(i4,i5))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i3,i6)*za(i5,i4)) + 
     & (2*za(i1,i3)**2*za(i3,i5)*zb(i4,i5))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*za(i3,i6)*za(i5,i4)) + 
     & (2*za(i1,i5)**2*zb(i4,i2)**2*zb(i4,i5))/
     &  (za(i1,i6)*za(i5,i4)*za(i6,i5)*zab2(i5,i2,i3,i4)*zb(i3,i2)*
     &    zb(i4,i3)) + (2*za(i1,i4)**2*zb(i4,i5)*zb(i5,i2)**2)/
     &  (za(i1,i6)*za(i5,i4)*za(i6,i4)*zab2(i4,i2,i3,i5)*zb(i3,i2)*
     &    zb(i5,i3)) + (2*zb(i4,i5)*
     &    (t(i2,i3,i6)*za(i1,i3)*za(i2,i6)*za(i6,i4)*za(i6,i5)*
     &       zb(i4,i1)*zb(i4,i3)*zb(i5,i1)*zb(i5,i3)*zb(i6,i1)*
     &       (-((za(i2,i5)*za(i3,i4) + za(i2,i4)*za(i3,i5))*
     &            zb(i3,i6)) + 
     &         (za(i1,i5)*za(i2,i4) + za(i1,i4)*za(i2,i5))*zb(i6,i1)) + 
     &      t(i1,i3,i6)*(t(i2,i3,i6)*
     &          (za(i2,i4)*za(i2,i5)*za(i6,i3)*za(i6,i4)*
     &             (za(i1,i5)*zb(i3,i6)*zb(i4,i1)**2 + 
     &               zab2(i5,i2,i6,i4)*zb(i3,i1)*zb(i4,i6))*zb(i5,i1)*
     &             zb(i5,i3) + 
     &            zb(i4,i3)*(-(za(i2,i4)*za(i2,i6)*za(i3,i5)**2*
     &                  za(i6,i4)*zb(i3,i1)*zb(i4,i6)*zb(i5,i1)*
     &                  zb(i5,i3)) + 
     &               zb(i4,i1)*
     &                (za(i2,i4)*za(i2,i5)*za(i6,i3)*za(i6,i5)*
     &                   (za(i1,i4)*zb(i3,i6)*zb(i5,i1)**2 + 
     &                     zab2(i4,i2,i6,i5)*zb(i3,i1)*zb(i5,i6)) - 
     &                  zb(i5,i3)*
     &                   (za(i2,i5)*za(i2,i6)*za(i3,i4)**2*za(i6,i5)*
     &                      zb(i3,i1)*zb(i5,i6) + 
     &                     zb(i5,i1)*
     &                      (-(za(i2,i3)*
     &                        (za(i2,i5)*za(i3,i4) + 
     &                       za(i2,i4)*za(i3,i5))*za(i6,i4)*za(i6,i5)*
     &                        zb(i3,i6)) + 
     &                        za(i2,i6)*
     &                        (za(i1,i5)*za(i2,i4)*za(i3,i5)*
     &                       za(i6,i4) + 
     &                        za(i1,i4)*za(i2,i5)*za(i3,i4)*za(i6,i5))*
     &                        zb(i6,i1)))))) - 
     &         za(i2,i3)*za(i2,i4)*za(i2,i5)*za(i6,i4)*za(i6,i5)*
     &          zb(i3,i1)*zb(i4,i3)*zb(i5,i3)*
     &          (za(i2,i3)*(zb(i4,i2)*zb(i5,i1) + 
     &               zb(i4,i1)*zb(i5,i2)) + 
     &            za(i6,i3)*(zb(i4,i6)*zb(i5,i1) + zb(i4,i1)*zb(i5,i6)))
     &           *zb(i6,i2))))/
     &  (t(i1,i3,i6)*t(i2,i3,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*
     &    za(i5,i4)*za(i6,i3)*za(i6,i4)*za(i6,i5)*zb(i3,i1)*zb(i3,i6)*
     &    zb(i4,i1)*zb(i4,i3)*zb(i5,i1)*zb(i5,i3)) + 
     & (2*za(i3,i5)**2*zb(i4,i5)*zb(i4,i6))/
     &  (za(i2,i5)*za(i3,i6)*za(i5,i4)*za(i6,i5)*zb(i4,i1)*zb(i6,i3)) + 
     & (2*za(i1,i5)**2*zb(i4,i5)*zb(i4,i6))/
     &  (za(i1,i6)*za(i2,i5)*za(i5,i4)*za(i6,i5)*zb(i4,i3)*zb(i6,i3)) + 
     & (2*za(i3,i4)**2*zb(i4,i5)*zb(i5,i6))/
     &  (za(i2,i4)*za(i3,i6)*za(i5,i4)*za(i6,i4)*zb(i5,i1)*zb(i6,i3)) + 
     & (2*za(i1,i4)**2*zb(i4,i5)*zb(i5,i6))/
     &  (za(i1,i6)*za(i2,i4)*za(i5,i4)*za(i6,i4)*zb(i5,i3)*zb(i6,i3)) + 
     & (2*za(i1,i3)**2*za(i1,i4)*zb(i4,i5)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i3,i6)*za(i5,i4)*zb(i6,i3))
     &  + (2*za(i1,i3)**2*za(i1,i5)*zb(i4,i5)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*za(i3,i6)*za(i5,i4)*zb(i6,i3))
     &  + (2*za(i3,i5)*zb(i4,i2)*zb(i4,i5)*zb(i6,i2))/
     &  (za(i3,i6)*za(i5,i4)*za(i6,i5)*zb(i3,i2)*zb(i4,i1)*zb(i6,i3)) + 
     & (2*za(i3,i4)*zb(i4,i5)*zb(i5,i2)*zb(i6,i2))/
     &  (za(i3,i6)*za(i5,i4)*za(i6,i4)*zb(i3,i2)*zb(i5,i1)*zb(i6,i3)) - 
     & (2*za(i2,i3)*zb(i4,i2)*zb(i4,i5)*zb(i6,i2)**2)/
     &  (t(i2,i3,i6)*za(i3,i6)*za(i5,i4)*zb(i3,i2)*zb(i4,i1)*zb(i6,i3))
     &  + (2*zb(i4,i5)*zb(i4,i6)*zb(i6,i2)**2)/
     &  (t(i2,i3,i6)*za(i5,i4)*zb(i3,i2)*zb(i4,i1)*zb(i6,i3)) - 
     & (2*za(i2,i3)*zb(i4,i5)*zb(i5,i2)*zb(i6,i2)**2)/
     &  (t(i2,i3,i6)*za(i3,i6)*za(i5,i4)*zb(i3,i2)*zb(i5,i1)*zb(i6,i3))
     &  + (2*zb(i4,i5)*zb(i5,i6)*zb(i6,i2)**2)/
     &  (t(i2,i3,i6)*za(i5,i4)*zb(i3,i2)*zb(i5,i1)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (2*za(i1,i3)**2*za(i3,i4)*zb(i6,i4))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i4)*za(i3,i5)*za(i4,i6)) + 
     & (2*za(i1,i3)**2*za(i3,i6)*zb(i6,i4))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*za(i3,i5)*za(i4,i6)) + 
     & (2*za(i1,i6)**2*zb(i4,i2)**2*zb(i6,i4))/
     &  (za(i1,i5)*za(i4,i6)*za(i5,i6)*zab2(i6,i2,i3,i4)*zb(i3,i2)*
     &    zb(i4,i3)) + (2*za(i3,i6)**2*zb(i4,i5)*zb(i6,i4))/
     &  (za(i2,i6)*za(i3,i5)*za(i4,i6)*za(i5,i6)*zb(i4,i1)*zb(i5,i3)) + 
     & (2*za(i1,i6)**2*zb(i4,i5)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i6)*za(i4,i6)*za(i5,i6)*zb(i4,i3)*zb(i5,i3)) + 
     & (2*za(i1,i3)**2*za(i1,i4)*zb(i5,i1)*zb(i6,i4))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i4)*za(i3,i5)*za(i4,i6)*zb(i5,i3))
     &  + (2*za(i1,i3)**2*za(i1,i6)*zb(i5,i1)*zb(i6,i4))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*za(i3,i5)*za(i4,i6)*zb(i5,i3))
     &  + (2*za(i3,i6)*zb(i4,i2)*zb(i5,i2)*zb(i6,i4))/
     &  (za(i3,i5)*za(i4,i6)*za(i5,i6)*zb(i3,i2)*zb(i4,i1)*zb(i5,i3)) - 
     & (2*za(i2,i3)*zb(i4,i2)*zb(i5,i2)**2*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i3,i5)*za(i4,i6)*zb(i3,i2)*zb(i4,i1)*zb(i5,i3))
     &  + (2*zb(i4,i5)*zb(i5,i2)**2*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i4,i6)*zb(i3,i2)*zb(i4,i1)*zb(i5,i3)) + 
     & (2*za(i3,i4)*zb(i5,i2)*zb(i6,i2)*zb(i6,i4))/
     &  (za(i3,i5)*za(i4,i6)*za(i5,i4)*zb(i3,i2)*zb(i5,i3)*zb(i6,i1)) - 
     & (2*za(i2,i3)*zb(i5,i2)**2*zb(i6,i2)*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i3,i5)*za(i4,i6)*zb(i3,i2)*zb(i5,i3)*zb(i6,i1))
     &  + (2*za(i1,i4)**2*zb(i6,i2)**2*zb(i6,i4))/
     &  (za(i1,i5)*za(i4,i6)*za(i5,i4)*zab2(i4,i2,i3,i6)*zb(i3,i2)*
     &    zb(i6,i3)) + (2*za(i3,i4)**2*zb(i6,i4)*zb(i6,i5))/
     &  (za(i2,i4)*za(i3,i5)*za(i4,i6)*za(i5,i4)*zb(i5,i3)*zb(i6,i1)) + 
     & (2*zb(i5,i2)**2*zb(i6,i4)*zb(i6,i5))/
     &  (t(i2,i3,i5)*za(i4,i6)*zb(i3,i2)*zb(i5,i3)*zb(i6,i1)) + 
     & (2*za(i1,i4)**2*zb(i6,i4)*zb(i6,i5))/
     &  (za(i1,i5)*za(i2,i4)*za(i4,i6)*za(i5,i4)*zb(i5,i3)*zb(i6,i3)) + 
     & (2*zb(i6,i4)*(t(i2,i3,i5)*za(i1,i3)*za(i2,i5)*za(i5,i4)*
     &       za(i5,i6)*zb(i4,i1)*zb(i4,i3)*zb(i5,i1)*
     &       (-((za(i2,i6)*za(i3,i4) + za(i2,i4)*za(i3,i6))*
     &            zb(i3,i5)) + 
     &         (za(i1,i6)*za(i2,i4) + za(i1,i4)*za(i2,i6))*zb(i5,i1))*
     &       zb(i6,i1)*zb(i6,i3) + 
     &      t(i1,i3,i5)*(-(za(i2,i3)*za(i2,i4)*za(i2,i6)*za(i5,i4)*
     &            za(i5,i6)*zb(i3,i1)*zb(i4,i3)*zb(i5,i2)*zb(i6,i3)*
     &            (za(i2,i3)*
     &               (zb(i4,i2)*zb(i6,i1) + zb(i4,i1)*zb(i6,i2)) + 
     &              za(i5,i3)*
     &               (zb(i4,i5)*zb(i6,i1) + zb(i4,i1)*zb(i6,i5)))) + 
     &         t(i2,i3,i5)*(za(i2,i4)*za(i2,i6)*za(i5,i3)*za(i5,i6)*
     &             zb(i4,i1)*zb(i4,i3)*
     &             (za(i1,i4)*zb(i3,i5)*zb(i6,i1)**2 + 
     &               zab2(i4,i2,i5,i6)*zb(i3,i1)*zb(i6,i5)) + 
     &            zb(i6,i3)*((za(i2,i4)*za(i2,i6)*za(i5,i3)*za(i5,i4)*
     &                   (za(i1,i6)*zb(i3,i5)*zb(i4,i1)**2 + 
     &                     zab2(i6,i2,i5,i4)*zb(i3,i1)*zb(i4,i5)) - 
     &                  zb(i4,i3)*
     &                   (za(i2,i4)*za(i2,i5)*za(i3,i6)**2*za(i5,i4)*
     &                      zb(i3,i1)*zb(i4,i5) + 
     &                     zb(i4,i1)*
     &                      (-(za(i2,i3)*
     &                        (za(i2,i6)*za(i3,i4) + 
     &                       za(i2,i4)*za(i3,i6))*za(i5,i4)*za(i5,i6)*
     &                        zb(i3,i5)) + 
     &                        za(i2,i5)*
     &                        (za(i1,i6)*za(i2,i4)*za(i3,i6)*
     &                       za(i5,i4) + 
     &                        za(i1,i4)*za(i2,i6)*za(i3,i4)*za(i5,i6))*
     &                        zb(i5,i1))))*zb(i6,i1) - 
     &               za(i2,i5)*za(i2,i6)*za(i3,i4)**2*za(i5,i6)*
     &                zb(i3,i1)*zb(i4,i1)*zb(i4,i3)*zb(i6,i5))))))/
     &  (t(i1,i3,i5)*t(i2,i3,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)*
     &    za(i4,i6)*za(i5,i3)*za(i5,i4)*za(i5,i6)*zb(i3,i1)*zb(i3,i5)*
     &    zb(i4,i1)*zb(i4,i3)*zb(i6,i1)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (2*(za(i1,i2)*((za(i1,i5)*za(i5,i3)*
     &           (-(((s(i2,i5) + s(i4,i5))*za(i1,i4))/
     &                (za(i2,i4)*zab2(i5,i2,i4,i3))) - 
     &             ((s(i2,i5) + s(i5,i6))*za(i1,i6))/
     &              (za(i2,i6)*zab2(i5,i2,i6,i3))))/
     &         (za(i2,i5)*za(i5,i4)*za(i5,i6)) - 
     &        (za(i1,i3)*zab2(i6,i2,i5,i1)*zb(i3,i5))/
     &         (za(i2,i4)*za(i2,i6)*za(i5,i6)*zb(i3,i1)**2) + 
     &        (za(i1,i3)**2*zb(i3,i5))/
     &         (za(i1,i5)*za(i2,i4)*za(i2,i6)*zb(i3,i1)) + 
     &        (za(i1,i5)**2*za(i5,i3)*zb(i5,i2))/
     &         (za(i2,i5)**2*za(i5,i4)*za(i5,i6)*zb(i3,i2))) + 
     &     (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i5,i3)*
     &         zb(i3,i1)**2*zb(i4,i5) + 
     &        za(i1,i3)*za(i1,i6)*za(i2,i4)*za(i5,i4)*zb(i3,i5)*
     &         zb(i4,i3)*(-(za(i2,i5)*zb(i2,i1)) + za(i5,i6)*zb(i6,i1))
     &         + za(i1,i4)*(za(i1,i3)*za(i2,i6)*za(i5,i4)**2*zb(i3,i5)*
     &            zb(i4,i1)*zb(i4,i3) + 
     &           za(i1,i6)*za(i5,i3)*zb(i3,i1)*
     &            (za(i2,i4)*zab2(i5,i2,i6,i3)*zb(i4,i5) - 
     &              za(i2,i6)*za(i5,i4)*zb(i4,i3)*zb(i6,i5))))/
     &      (za(i2,i4)*za(i2,i6)*za(i5,i4)**2*za(i5,i6)*zb(i3,i1)**2*
     &        zb(i4,i3))))/(za(i1,i4)*za(i1,i6))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (2*((za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i6,i3)*
     &         zb(i3,i1)**2*zb(i4,i6) + 
     &        za(i1,i3)*za(i1,i5)*za(i2,i4)*za(i6,i4)*zb(i3,i6)*
     &         zb(i4,i3)*(-(za(i2,i6)*zb(i2,i1)) + za(i6,i5)*zb(i5,i1))
     &         + za(i1,i4)*(za(i1,i3)*za(i2,i5)*za(i6,i4)**2*zb(i3,i6)*
     &            zb(i4,i1)*zb(i4,i3) + 
     &           za(i1,i5)*za(i6,i3)*zb(i3,i1)*
     &            (za(i2,i4)*zab2(i6,i2,i5,i3)*zb(i4,i6) - 
     &              za(i2,i5)*za(i6,i4)*zb(i4,i3)*zb(i5,i6))))/
     &      (za(i2,i4)*za(i2,i5)*za(i6,i4)**2*za(i6,i5)*zb(i3,i1)**2*
     &        zb(i4,i3)) + za(i1,i2)*
     &      ((za(i1,i6)*za(i6,i3)*
     &           (-(((s(i2,i6) + s(i4,i6))*za(i1,i4))/
     &                (za(i2,i4)*zab2(i6,i2,i4,i3))) - 
     &             ((s(i2,i6) + s(i5,i6))*za(i1,i5))/
     &              (za(i2,i5)*zab2(i6,i2,i5,i3))))/
     &         (za(i2,i6)*za(i6,i4)*za(i6,i5)) - 
     &        (za(i1,i3)*zab2(i5,i2,i6,i1)*zb(i3,i6))/
     &         (za(i2,i4)*za(i2,i5)*za(i6,i5)*zb(i3,i1)**2) + 
     &        (za(i1,i3)**2*zb(i3,i6))/
     &         (za(i1,i6)*za(i2,i4)*za(i2,i5)*zb(i3,i1)) + 
     &        (za(i1,i6)**2*za(i6,i3)*zb(i6,i2))/
     &         (za(i2,i6)**2*za(i6,i4)*za(i6,i5)*zb(i3,i2)))))/
     & (za(i1,i4)*za(i1,i5))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (2*(za(i1,i2)*((za(i1,i4)*za(i4,i3)*
     &           (-(((s(i2,i4) + s(i4,i5))*za(i1,i5))/
     &                (za(i2,i5)*zab2(i4,i2,i5,i3))) - 
     &             ((s(i2,i4) + s(i4,i6))*za(i1,i6))/
     &              (za(i2,i6)*zab2(i4,i2,i6,i3))))/
     &         (za(i2,i4)*za(i4,i5)*za(i4,i6)) - 
     &        (za(i1,i3)*zab2(i6,i2,i4,i1)*zb(i3,i4))/
     &         (za(i2,i5)*za(i2,i6)*za(i4,i6)*zb(i3,i1)**2) + 
     &        (za(i1,i3)**2*zb(i3,i4))/
     &         (za(i1,i4)*za(i2,i5)*za(i2,i6)*zb(i3,i1)) + 
     &        (za(i1,i4)**2*za(i4,i3)*zb(i4,i2))/
     &         (za(i2,i4)**2*za(i4,i5)*za(i4,i6)*zb(i3,i2))) + 
     &     (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i4,i3)*
     &         zb(i3,i1)**2*zb(i5,i4) + 
     &        za(i1,i3)*za(i1,i6)*za(i2,i5)*za(i4,i5)*zb(i3,i4)*
     &         zb(i5,i3)*(-(za(i2,i4)*zb(i2,i1)) + za(i4,i6)*zb(i6,i1))
     &         + za(i1,i5)*(za(i1,i3)*za(i2,i6)*za(i4,i5)**2*zb(i3,i4)*
     &            zb(i5,i1)*zb(i5,i3) + 
     &           za(i1,i6)*za(i4,i3)*zb(i3,i1)*
     &            (za(i2,i5)*zab2(i4,i2,i6,i3)*zb(i5,i4) - 
     &              za(i2,i6)*za(i4,i5)*zb(i5,i3)*zb(i6,i4))))/
     &      (za(i2,i5)*za(i2,i6)*za(i4,i5)**2*za(i4,i6)*zb(i3,i1)**2*
     &        zb(i5,i3))))/(za(i1,i5)*za(i1,i6))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (-2*za(i1,i2)*za(i1,i3)**3*zb(i5,i1))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i3,i4)) + 
     & (2*za(i1,i3)**2*za(i2,i3)*zb(i3,i4)*zb(i5,i1))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i2)*za(i1,i3)**2*zb(i4,i1)*zb(i5,i1))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i6)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i2)*za(i1,i3)**2*zb(i4,i1)*zb(i5,i1))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**3*za(i2,i3)*zb(i5,i3))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i3,i4)) - 
     & (2*za(i1,i3)**2*za(i2,i3)*zb(i4,i1)*zb(i5,i3))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i6)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**3*za(i2,i4)*zb(i5,i4))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i3,i4)) - 
     & (2*za(i1,i3)**2*za(i2,i4)*zb(i4,i1)*zb(i5,i4))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i6)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**2*za(i2,i5)*zb(i5,i1)*zb(i5,i4))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)**2*za(i4,i3)*zb(i3,i1)) - 
     & (za(i1,i2)*za(i1,i3)**2*za(i2,i3)*zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i3,i4))
     &  - (za(i1,i3)**2*za(i1,i6)*za(i2,i3)*zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i6)*za(i3,i4)*za(i6,i5))
     &  - (za(i1,i2)*za(i2,i3)**2*zb(i6,i2))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1)) + 
     & (za(i1,i6)*za(i2,i3)*za(i3,i6)*zb(i6,i2))/
     &  (za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i4,i6)*za(i5,i6)*zb(i3,i1)) + 
     & (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i4,i1)*zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1))
     &  - (za(i1,i3)*za(i1,i6)*za(i2,i3)*zb(i4,i1)*zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i5,i6)*zb(i3,i1))
     &  - (za(i1,i2)**2*za(i2,i3)*zb(i4,i1)*zb(i6,i2))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1)*
     &    zb(i3,i4)) + (za(i1,i2)*za(i1,i6)*za(i3,i6)*zb(i4,i1)*
     &    zb(i6,i2))/
     &  (za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i4,i6)*za(i5,i6)*zb(i3,i1)*
     &    zb(i3,i4)) + (za(i1,i2)**2*za(i1,i3)*zb(i4,i1)**2*zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1)*
     &    zb(i3,i4)) - (za(i1,i2)*za(i1,i3)*za(i1,i6)*zb(i4,i1)**2*
     &    zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i5,i6)*zb(i3,i1)*
     &    zb(i3,i4)) + (za(i1,i2)*za(i1,i6)*za(i3,i6)**2*zb(i4,i2)*
     &    zb(i6,i2))/
     &  (za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i3,i4,i2)*zb(i3,i4)) + 
     & (za(i1,i2)**2*za(i1,i3)**2*zb(i4,i1)*zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i3,i4)*
     &    zb(i4,i3)) + (za(i1,i2)*za(i1,i3)**2*za(i1,i6)*zb(i4,i1)*
     &    zb(i6,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i5)*za(i2,i6)*za(i3,i4)*za(i6,i5)*
     &    zb(i4,i3)) - (za(i1,i2)*za(i1,i3)*za(i1,i6)*zb(i4,i2)*
     &    zb(i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)*za(i3,i4)*za(i6,i5)*zb(i3,i2)*
     &    zb(i4,i3)) - (za(i1,i2)*za(i1,i6)*za(i3,i6)*zb(i4,i2)**2*
     &    zb(i6,i2))/
     &  (za(i1,i5)*za(i2,i6)*za(i3,i4)*za(i6,i5)*zab2(i6,i3,i4,i2)*
     &    zb(i3,i2)*zb(i4,i3)) - 
     & (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i5,i1)*zb(i6,i2))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i3)*zb(i3,i1))
     &  - (za(i1,i2)*za(i1,i3)*za(i3,i6)*zb(i5,i1)*zb(i6,i2))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i4,i6)*zb(i3,i1))
     &  + (za(i1,i3)*za(i2,i3)**2*zb(i5,i3)*zb(i6,i2))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i3)*zb(i3,i1))
     &  + (za(i1,i3)*za(i2,i3)*za(i3,i6)*zb(i5,i3)*zb(i6,i2))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*za(i4,i3)*za(i4,i6)*zb(i3,i1))
     &  - (za(i1,i2)**2*za(i2,i3)**2*zb(i6,i2)*zb(i6,i4))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i3)*
     &    zab2(i2,i3,i4,i6)*zb(i3,i4)) - 
     & (za(i1,i2)**2*za(i1,i3)*zb(i6,i2)*zb(i6,i4))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i3,i4)*zb(i4,i3)*
     &    zb(i6,i3)) - (za(i1,i2)**2*za(i2,i3)*zb(i6,i2)*zb(i6,i4)**2)/
     &  (za(i1,i5)*za(i2,i5)*za(i2,i6)*za(i3,i4)*zab2(i2,i3,i4,i6)*
     &    zb(i4,i3)*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i2)*za(i1,i3)**2*za(i2,i3)*zb(i4,i2))/
     &    (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i6)*
     &      za(i3,i5))) - (za(i1,i3)**2*za(i1,i4)*za(i2,i3)*zb(i4,i2))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i3,i5)*za(i4,i6))
     &  - (za(i1,i2)*za(i2,i3)**2*zb(i4,i2))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i3)*zb(i3,i1)) + 
     & (za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i2))/
     &  (za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i5,i4)*za(i6,i4)*zb(i3,i1)) - 
     & (za(i1,i2)**2*za(i2,i3)**2*zb(i4,i2)*zb(i4,i5))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i3)*
     &    zab2(i2,i3,i5,i4)*zb(i3,i5)) + 
     & (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i4,i2)*zb(i5,i1))/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i5,i3)*zb(i3,i1))
     &  - (za(i1,i3)*za(i1,i4)*za(i2,i3)*zb(i4,i2)*zb(i5,i1))/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i6,i4)*zb(i3,i1))
     &  - (za(i1,i2)**2*za(i2,i3)*zb(i4,i2)*zb(i5,i1))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i5,i3)*zb(i3,i1)*
     &    zb(i3,i5)) + (za(i1,i2)*za(i1,i4)*za(i3,i4)*zb(i4,i2)*
     &    zb(i5,i1))/
     &  (za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i5,i4)*za(i6,i4)*zb(i3,i1)*
     &    zb(i3,i5)) + (za(i1,i2)**2*za(i1,i3)*zb(i4,i2)*zb(i5,i1)**2)/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i5,i3)*zb(i3,i1)*
     &    zb(i3,i5)) - (za(i1,i2)*za(i1,i3)*za(i1,i4)*zb(i4,i2)*
     &    zb(i5,i1)**2)/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i6,i4)*zb(i3,i1)*
     &    zb(i3,i5)) + (za(i1,i2)*za(i1,i4)*za(i3,i4)**2*zb(i4,i2)*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i5,i4)*za(i6,i4)*
     &    zab2(i4,i3,i5,i2)*zb(i3,i5)) - 
     & (za(i1,i2)**2*za(i1,i3)*zb(i4,i2)*zb(i4,i5))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i3,i5)*zb(i4,i3)*
     &    zb(i5,i3)) - (za(i1,i2)**2*za(i2,i3)*zb(i4,i2)*zb(i4,i5)**2)/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i3,i5)*zab2(i2,i3,i5,i4)*
     &    zb(i4,i3)*zb(i5,i3)) + 
     & (za(i1,i2)**2*za(i1,i3)**2*zb(i4,i2)*zb(i5,i1))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i6)*za(i3,i5)*
     &    zb(i5,i3)) + (za(i1,i2)*za(i1,i3)**2*za(i1,i4)*zb(i4,i2)*
     &    zb(i5,i1))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i3,i5)*za(i4,i6)*
     &    zb(i5,i3)) - (za(i1,i2)*za(i1,i3)*za(i1,i4)*zb(i4,i2)*
     &    zb(i5,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i3,i5)*za(i4,i6)*zb(i3,i2)*
     &    zb(i5,i3)) - (za(i1,i2)*za(i1,i4)*za(i3,i4)*zb(i4,i2)*
     &    zb(i5,i2)**2)/
     &  (za(i1,i6)*za(i2,i4)*za(i3,i5)*za(i4,i6)*zab2(i4,i3,i5,i2)*
     &    zb(i3,i2)*zb(i5,i3)) - 
     & (2*za(i1,i2)*za(i1,i3)**3*zb(i6,i1))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i3,i5)) + 
     & (2*za(i1,i3)**2*za(i2,i3)*zb(i3,i5)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)**2*za(i5,i3)*zb(i3,i1)) - 
     & (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i4,i2)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i5,i3)*zb(i3,i1))
     &  - (za(i1,i2)*za(i1,i3)*za(i3,i4)*zb(i4,i2)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i5,i4)*zb(i3,i1))
     &  + (2*za(i1,i2)*za(i1,i3)**2*zb(i5,i1)*zb(i6,i1))/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)**2*za(i5,i3)*zb(i3,i1)) + 
     & (2*za(i1,i2)*za(i1,i3)**2*zb(i5,i1)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)**2*za(i5,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**3*za(i2,i3)*zb(i6,i3))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i3,i5)) + 
     & (za(i1,i3)*za(i2,i3)**2*zb(i4,i2)*zb(i6,i3))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i5,i3)*zb(i3,i1))
     &  + (za(i1,i3)*za(i2,i3)*za(i3,i4)*zb(i4,i2)*zb(i6,i3))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i5,i3)*za(i5,i4)*zb(i3,i1))
     &  - (2*za(i1,i3)**2*za(i2,i3)*zb(i5,i1)*zb(i6,i3))/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)**2*za(i5,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**3*za(i2,i5)*zb(i6,i5))/
     &  (t(i1,i3,i5)*za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i3,i5)) - 
     & (2*za(i1,i3)**2*za(i2,i5)*zb(i5,i1)*zb(i6,i5))/
     &  (t(i1,i3,i5)*za(i1,i6)*za(i2,i4)**2*za(i5,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**2*za(i2,i6)*zb(i6,i1)*zb(i6,i5))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)**2*za(i5,i3)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i2)*za(i1,i3)**2*za(i2,i3)*zb(i5,i2))/
     &    (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     &      za(i3,i4))) - (za(i1,i3)**2*za(i1,i5)*za(i2,i3)*zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6))
     &  - (za(i1,i2)*za(i2,i3)**2*zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1)) + 
     & (za(i1,i5)*za(i2,i3)*za(i3,i5)*zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i4,i5)*za(i6,i5)*zb(i3,i1)) + 
     & (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i4,i1)*zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1))
     &  - (za(i1,i3)*za(i1,i5)*za(i2,i3)*zb(i4,i1)*zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i6,i5)*zb(i3,i1))
     &  - (za(i1,i2)**2*za(i2,i3)*zb(i4,i1)*zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1)*
     &    zb(i3,i4)) + (za(i1,i2)*za(i1,i5)*za(i3,i5)*zb(i4,i1)*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i4,i5)*za(i6,i5)*zb(i3,i1)*
     &    zb(i3,i4)) + (za(i1,i2)**2*za(i1,i3)*zb(i4,i1)**2*zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i4,i3)*zb(i3,i1)*
     &    zb(i3,i4)) - (za(i1,i2)*za(i1,i3)*za(i1,i5)*zb(i4,i1)**2*
     &    zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i6,i5)*zb(i3,i1)*
     &    zb(i3,i4)) + (za(i1,i2)*za(i1,i5)*za(i3,i5)**2*zb(i4,i2)*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i4,i5)*za(i6,i5)*
     &    zab2(i5,i3,i4,i2)*zb(i3,i4)) + 
     & (za(i1,i2)**2*za(i1,i3)**2*zb(i4,i1)*zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i3,i4)*
     &    zb(i4,i3)) + (za(i1,i2)*za(i1,i3)**2*za(i1,i5)*zb(i4,i1)*
     &    zb(i5,i2))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6)*
     &    zb(i4,i3)) - (za(i1,i2)*za(i1,i3)*za(i1,i5)*zb(i4,i2)*
     &    zb(i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6)*zb(i3,i2)*
     &    zb(i4,i3)) - (za(i1,i2)*za(i1,i5)*za(i3,i5)*zb(i4,i2)**2*
     &    zb(i5,i2))/
     &  (za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6)*zab2(i5,i3,i4,i2)*
     &    zb(i3,i2)*zb(i4,i3)) - 
     & (za(i1,i2)**2*za(i2,i3)**2*zb(i5,i2)*zb(i5,i4))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)*za(i4,i3)*
     &    zab2(i2,i3,i4,i5)*zb(i3,i4)) - 
     & (za(i1,i2)**2*za(i1,i3)*zb(i5,i2)*zb(i5,i4))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i3,i4)*zb(i4,i3)*
     &    zb(i5,i3)) - (za(i1,i2)**2*za(i2,i3)*zb(i5,i2)*zb(i5,i4)**2)/
     &  (za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i3,i4)*zab2(i2,i3,i4,i5)*
     &    zb(i4,i3)*zb(i5,i3)) - 
     & (2*za(i1,i2)*za(i1,i3)**3*zb(i6,i1))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i3,i4)) + 
     & (2*za(i1,i3)**2*za(i2,i3)*zb(i3,i4)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i2)*za(i1,i3)**2*zb(i4,i1)*zb(i6,i1))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i2)*za(i1,i3)**2*zb(i4,i1)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)**2*za(i4,i3)*zb(i3,i1)) - 
     & (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i5,i2)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i3)*zb(i3,i1))
     &  - (za(i1,i2)*za(i1,i3)*za(i3,i5)*zb(i5,i2)*zb(i6,i1))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i4,i5)*zb(i3,i1))
     &  + (2*za(i1,i3)**3*za(i2,i3)*zb(i6,i3))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i3,i4)) - 
     & (2*za(i1,i3)**2*za(i2,i3)*zb(i4,i1)*zb(i6,i3))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)**2*za(i4,i3)*zb(i3,i1)) + 
     & (za(i1,i3)*za(i2,i3)**2*zb(i5,i2)*zb(i6,i3))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i3)*zb(i3,i1))
     &  + (za(i1,i3)*za(i2,i3)*za(i3,i5)*zb(i5,i2)*zb(i6,i3))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*za(i4,i3)*za(i4,i5)*zb(i3,i1))
     &  + (2*za(i1,i3)**3*za(i2,i4)*zb(i6,i4))/
     &  (t(i1,i3,i4)*za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i3,i4)) - 
     & (2*za(i1,i3)**2*za(i2,i4)*zb(i4,i1)*zb(i6,i4))/
     &  (t(i1,i3,i4)*za(i1,i6)*za(i2,i5)**2*za(i4,i3)*zb(i3,i1)) + 
     & (2*za(i1,i3)**2*za(i2,i6)*zb(i6,i1)*zb(i6,i4))/
     &  (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)**2*za(i4,i3)*zb(i3,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i2)**3*za(i2,i3)*zb(i4,i2))/
     &    (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)*
     &      zb(i4,i3))) - (za(i1,i2)**3*za(i2,i3)*zb(i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)*zb(i5,i3))
     &   - (za(i1,i2)**3*za(i2,i3)*zb(i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)**2*zb(i6,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i3)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i1))/
     &  (za(i2,i5)*za(i2,i6)*za(i3,i4)*za(i3,i5)*za(i4,i6)*zb(i3,i1)**2)
     &   - (za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     &  (za(i2,i6)*za(i3,i5)*za(i4,i5)*za(i4,i6)*zb(i3,i1)) - 
     & (za(i1,i3)*za(i2,i3)*za(i4,i5)*zb(i4,i1)*zb(i5,i1))/
     &  (za(i2,i5)*za(i2,i6)*za(i3,i4)*za(i3,i5)*za(i4,i6)*zb(i3,i1)**2)
     &   - (za(i1,i4)*za(i3,i4)*zb(i5,i4))/
     &  (za(i2,i6)*za(i4,i5)**2*za(i4,i6)*zb(i5,i3)) + 
     & (za(i1,i4)*zb(i5,i1)*zb(i5,i4))/
     &  (za(i2,i6)*za(i4,i5)*za(i4,i6)*zb(i3,i1)*zb(i5,i3)) - 
     & (za(i1,i4)*za(i2,i4)*zb(i2,i1)*zb(i5,i2)*zb(i5,i4))/
     &  (za(i2,i6)*za(i4,i5)**2*za(i4,i6)*zb(i3,i2)*zb(i5,i1)*zb(i5,i3))
     &   - (za(i1,i3)*za(i2,i3)*zb(i4,i1)*zb(i5,i2)*zb(i6,i1))/
     &  (t(i2,i5,i6)*za(i2,i6)*za(i3,i4)*za(i3,i5)*zb(i3,i1)**2) + 
     & (za(i1,i4)*za(i2,i3)*zb(i3,i1)*zb(i5,i2)*zb(i5,i4)*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*za(i4,i5)**2*zb(i3,i2)*zb(i5,i1)*
     &    zb(i5,i3)) - (za(i1,i3)*za(i2,i3)*zb(i4,i1)*zb(i6,i1)*
     &    zb(i6,i2))/
     &  (t(i2,i5,i6)*za(i2,i5)*za(i3,i4)*za(i3,i5)*zb(i3,i1)**2) + 
     & (za(i1,i4)*zb(i5,i2)*zb(i5,i4)*zb(i6,i1)*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i4,i5)**2*zb(i3,i2)*zb(i5,i1)*zb(i5,i3)) + 
     & (za(i1,i4)*za(i2,i4)*zb(i2,i1)*zb(i5,i1)*zb(i6,i4))/
     &  (za(i2,i5)*za(i4,i5)*za(i4,i6)**2*zb(i3,i1)*zb(i5,i3)*zb(i6,i1))
     &   - (za(i1,i4)*zb(i5,i1)**2*zb(i6,i4))/
     &  (za(i2,i5)*za(i4,i6)**2*zb(i3,i1)*zb(i5,i3)*zb(i6,i1)) - 
     & (za(i1,i4)*za(i2,i3)**2*zb(i2,i1)*zb(i5,i2)*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i2,i5)*za(i3,i5)*za(i4,i6)**2*zb(i5,i3)*
     &    zb(i6,i1)) - (za(i1,i4)*za(i3,i4)*zb(i2,i1)*zb(i5,i2)*
     &    zb(i6,i4))/
     &  (za(i3,i5)*za(i4,i5)*za(i4,i6)**2*zb(i3,i2)*zb(i5,i3)*zb(i6,i1))
     &   + (za(i1,i4)*za(i2,i3)*zb(i5,i1)*zb(i5,i2)*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i2,i5)*za(i4,i6)**2*zb(i5,i3)*zb(i6,i1)) - 
     & (za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i5,i2)**2*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i3,i5)*za(i4,i6)**2*zb(i3,i2)*zb(i5,i3)*
     &    zb(i6,i1)) + (za(i1,i4)*zb(i5,i1)*zb(i5,i2)**2*zb(i6,i4))/
     &  (t(i2,i3,i5)*za(i4,i6)**2*zb(i3,i2)*zb(i5,i3)*zb(i6,i1)) - 
     & (za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i5,i4)*zb(i6,i5))/
     &  (t(i2,i3,i6)*za(i2,i6)*za(i4,i5)**2*zb(i5,i1)*zb(i5,i3)) - 
     & (za(i1,i3)*za(i2,i3)*za(i5,i6)*zb(i4,i1)*zb(i6,i1)*zb(i6,i5))/
     &  (t(i2,i5,i6)*za(i2,i5)*za(i2,i6)*za(i3,i4)*za(i3,i5)*
     &    zb(i3,i1)**2) + (za(i1,i4)*za(i3,i6)*zb(i5,i4)*zb(i6,i1)*
     &    zb(i6,i5))/
     &  (t(i2,i3,i6)*za(i2,i6)*za(i4,i5)**2*zb(i5,i1)*zb(i5,i3))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i5)*za(i3,i5)*zb(i4,i5))/
     &    (za(i2,i6)*za(i5,i4)**2*za(i5,i6)*zb(i4,i3))) + 
     & (za(i1,i5)*zb(i4,i1)*zb(i4,i5))/
     &  (za(i2,i6)*za(i5,i4)*za(i5,i6)*zb(i3,i1)*zb(i4,i3)) - 
     & (za(i1,i5)*za(i2,i5)*zb(i2,i1)*zb(i4,i2)*zb(i4,i5))/
     &  (za(i2,i6)*za(i5,i4)**2*za(i5,i6)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3))
     &   + (za(i1,i3)*za(i2,i3)*za(i2,i5)*zb(i2,i1)*zb(i5,i1))/
     &  (za(i2,i4)*za(i2,i6)*za(i3,i4)*za(i3,i5)*za(i5,i6)*zb(i3,i1)**2)
     &   - (za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     &  (za(i2,i6)*za(i3,i4)*za(i5,i4)*za(i5,i6)*zb(i3,i1)) - 
     & (za(i1,i3)*za(i2,i3)*za(i5,i4)*zb(i4,i1)*zb(i5,i1))/
     &  (za(i2,i4)*za(i2,i6)*za(i3,i4)*za(i3,i5)*za(i5,i6)*zb(i3,i1)**2)
     &   - (za(i1,i3)*za(i2,i3)*zb(i4,i2)*zb(i5,i1)*zb(i6,i1))/
     &  (t(i2,i4,i6)*za(i2,i6)*za(i3,i4)*za(i3,i5)*zb(i3,i1)**2) + 
     & (za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i4,i2)*zb(i4,i5)*zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i2,i6)*za(i5,i4)**2*zb(i3,i2)*zb(i4,i1)*
     &    zb(i4,i3)) + (za(i1,i5)*zb(i4,i2)*zb(i4,i5)*zb(i6,i1)*
     &    zb(i6,i2))/
     &  (t(i2,i3,i6)*za(i5,i4)**2*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) - 
     & (za(i1,i3)*za(i2,i3)*zb(i5,i1)*zb(i6,i1)*zb(i6,i2))/
     &  (t(i2,i4,i6)*za(i2,i4)*za(i3,i4)*za(i3,i5)*zb(i3,i1)**2) - 
     & (za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i4,i5)*zb(i6,i4))/
     &  (t(i2,i3,i6)*za(i2,i6)*za(i5,i4)**2*zb(i4,i1)*zb(i4,i3)) + 
     & (za(i1,i5)*za(i3,i6)*zb(i4,i5)*zb(i6,i1)*zb(i6,i4))/
     &  (t(i2,i3,i6)*za(i2,i6)*za(i5,i4)**2*zb(i4,i1)*zb(i4,i3)) - 
     & (za(i1,i3)*za(i2,i3)*za(i4,i6)*zb(i5,i1)*zb(i6,i1)*zb(i6,i4))/
     &  (t(i2,i4,i6)*za(i2,i4)*za(i2,i6)*za(i3,i4)*za(i3,i5)*
     &    zb(i3,i1)**2) + (za(i1,i5)*za(i2,i5)*zb(i2,i1)*zb(i4,i1)*
     &    zb(i6,i5))/
     &  (za(i2,i4)*za(i5,i4)*za(i5,i6)**2*zb(i3,i1)*zb(i4,i3)*zb(i6,i1))
     &   - (za(i1,i5)*zb(i4,i1)**2*zb(i6,i5))/
     &  (za(i2,i4)*za(i5,i6)**2*zb(i3,i1)*zb(i4,i3)*zb(i6,i1)) - 
     & (za(i1,i5)*za(i2,i3)**2*zb(i2,i1)*zb(i4,i2)*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i2,i4)*za(i3,i4)*za(i5,i6)**2*zb(i4,i3)*
     &    zb(i6,i1)) - (za(i1,i5)*za(i3,i5)*zb(i2,i1)*zb(i4,i2)*
     &    zb(i6,i5))/
     &  (za(i3,i4)*za(i5,i4)*za(i5,i6)**2*zb(i3,i2)*zb(i4,i3)*zb(i6,i1))
     &   + (za(i1,i5)*za(i2,i3)*zb(i4,i1)*zb(i4,i2)*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i2,i4)*za(i5,i6)**2*zb(i4,i3)*zb(i6,i1)) - 
     & (za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i4,i2)**2*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i3,i4)*za(i5,i6)**2*zb(i3,i2)*zb(i4,i3)*
     &    zb(i6,i1)) + (za(i1,i5)*zb(i4,i1)*zb(i4,i2)**2*zb(i6,i5))/
     &  (t(i2,i3,i4)*za(i5,i6)**2*zb(i3,i2)*zb(i4,i3)*zb(i6,i1))
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  (za(i1,i6)*za(i2,i6)*zb(i2,i1)*zb(i4,i6)*zb(i5,i1))/
     &  (za(i2,i5)*za(i6,i4)**2*za(i6,i5)*zb(i3,i1)*zb(i4,i1)*zb(i5,i3))
     &   - (za(i1,i6)*zb(i4,i6)*zb(i5,i1)**2)/
     &  (za(i2,i5)*za(i6,i4)**2*zb(i3,i1)*zb(i4,i1)*zb(i5,i3)) - 
     & (za(i1,i6)*za(i2,i3)**2*zb(i2,i1)*zb(i4,i6)*zb(i5,i2))/
     &  (t(i2,i3,i5)*za(i2,i5)*za(i3,i5)*za(i6,i4)**2*zb(i4,i1)*
     &    zb(i5,i3)) - (za(i1,i6)*za(i3,i6)*zb(i2,i1)*zb(i4,i6)*
     &    zb(i5,i2))/
     &  (za(i3,i5)*za(i6,i4)**2*za(i6,i5)*zb(i3,i2)*zb(i4,i1)*zb(i5,i3))
     &   + (za(i1,i6)*za(i2,i3)*zb(i4,i6)*zb(i5,i1)*zb(i5,i2))/
     &  (t(i2,i3,i5)*za(i2,i5)*za(i6,i4)**2*zb(i4,i1)*zb(i5,i3)) - 
     & (za(i1,i6)*za(i2,i3)*zb(i2,i1)*zb(i4,i6)*zb(i5,i2)**2)/
     &  (t(i2,i3,i5)*za(i3,i5)*za(i6,i4)**2*zb(i3,i2)*zb(i4,i1)*
     &    zb(i5,i3)) + (za(i1,i6)*zb(i4,i6)*zb(i5,i1)*zb(i5,i2)**2)/
     &  (t(i2,i3,i5)*za(i6,i4)**2*zb(i3,i2)*zb(i4,i1)*zb(i5,i3)) - 
     & (za(i1,i6)*za(i3,i6)*zb(i5,i6))/
     &  (za(i2,i4)*za(i6,i4)*za(i6,i5)**2*zb(i5,i3)) - 
     & (za(i1,i6)*za(i2,i3)*zb(i2,i1)*zb(i4,i5)*zb(i5,i6))/
     &  (t(i2,i3,i4)*za(i2,i4)*za(i6,i5)**2*zb(i5,i1)*zb(i5,i3)) + 
     & (za(i1,i6)*za(i3,i4)*zb(i4,i1)*zb(i4,i5)*zb(i5,i6))/
     &  (t(i2,i3,i4)*za(i2,i4)*za(i6,i5)**2*zb(i5,i1)*zb(i5,i3)) + 
     & (za(i1,i6)*zb(i5,i1)*zb(i5,i6))/
     &  (za(i2,i4)*za(i6,i4)*za(i6,i5)*zb(i3,i1)*zb(i5,i3)) - 
     & (za(i1,i6)*za(i2,i6)*zb(i2,i1)*zb(i5,i2)*zb(i5,i6))/
     &  (za(i2,i4)*za(i6,i4)*za(i6,i5)**2*zb(i3,i2)*zb(i5,i1)*zb(i5,i3))
     &   + (za(i1,i6)*za(i2,i3)*zb(i3,i1)*zb(i4,i2)*zb(i5,i2)*
     &    zb(i5,i6))/
     &  (t(i2,i3,i4)*za(i2,i4)*za(i6,i5)**2*zb(i3,i2)*zb(i5,i1)*
     &    zb(i5,i3)) + (za(i1,i6)*zb(i4,i1)*zb(i4,i2)*zb(i5,i2)*
     &    zb(i5,i6))/
     &  (t(i2,i3,i4)*za(i6,i5)**2*zb(i3,i2)*zb(i5,i1)*zb(i5,i3)) + 
     & (za(i1,i3)*za(i2,i3)*za(i2,i6)*zb(i2,i1)*zb(i6,i1))/
     &  (za(i2,i4)*za(i2,i5)*za(i3,i5)*za(i3,i6)*za(i6,i4)*zb(i3,i1)**2)
     &   - (za(i1,i3)*za(i3,i6)*zb(i6,i1))/
     &  (za(i2,i4)*za(i3,i5)*za(i6,i4)*za(i6,i5)*zb(i3,i1)) - 
     & (za(i1,i3)*za(i2,i3)*zb(i4,i1)*zb(i4,i2)*zb(i6,i1))/
     &  (t(i2,i4,i5)*za(i2,i5)*za(i3,i5)*za(i3,i6)*zb(i3,i1)**2) - 
     & (za(i1,i3)*za(i2,i3)*za(i5,i4)*zb(i4,i1)*zb(i4,i5)*zb(i6,i1))/
     &  (t(i2,i4,i5)*za(i2,i4)*za(i2,i5)*za(i3,i5)*za(i3,i6)*
     &    zb(i3,i1)**2) - (za(i1,i3)*za(i2,i3)*za(i6,i5)*zb(i5,i1)*
     &    zb(i6,i1))/
     &  (za(i2,i4)*za(i2,i5)*za(i3,i5)*za(i3,i6)*za(i6,i4)*zb(i3,i1)**2)
     &   - (za(i1,i3)*za(i2,i3)*zb(i4,i1)*zb(i5,i2)*zb(i6,i1))/
     &  (t(i2,i4,i5)*za(i2,i4)*za(i3,i5)*za(i3,i6)*zb(i3,i1)**2)
     & )


      aaaa_MHV_r = aaaa_MHV_r + (
     &  -((za(i1,i2)**2*za(i1,i3)**2)/
     &   (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)))
     & )


      return
      end
