      subroutine reals2(u,g1,d,t1,b1,g2,mq,ma,za,zb,amp)
************************************************************************
*                                                                      *
* AUTHORS: R. FREDERIX AND F. TRAMONTANO                               *
* DATE  : 7/16/2008                                                    *
*                                                                      *
************************************************************************
c--- 0  ->  u~ + g1 + d + t1 + b1~ + g2  (t-channel single-top)
c--- where the g1 (g2) is the gluon attached to the massive (massless) quark
c--- line. Arguments of amp are helicities of gluon1, top, bottom and gluon2,
c--- respectively. Spin of top and bottom quarks should be projected on gluon1.
c--- mq is the mass of the top quark and ma the mass of the bottom quark.
 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer u,g1,g2,d,b1,t1
      double precision g1g2b,g1g2t,g2b,g2t,mq,ma,mq2,ma2
      double complex amp(2,2,2,2)
 
 
      mq2 = mq**2
      ma2 = ma**2
      g2b=s(g2,b1)/2d0+s(g1,g2)*ma2/2d0/s(g1,b1)
      g2t=s(g2,t1)/2d0+s(g1,g2)*mq2/2d0/s(g1,t1)
      g1g2b=s(g1,g2)/2d0+s(g1,b1)/2d0+g2b
      g1g2t=s(g1,g2)/2d0+s(g1,t1)/2d0+g2t
 
      amp(1,1,1,1)=
     &        (mq*(2*g1g2b*za(g1,b1)*za(g1,t1)*zb(g1,b1)*zb(g1,t1)**2*
     &      (-(g2b*za(d,g1)*za(g1,t1)*
     &           (za(g1,g2)*zb(g1,g2) + za(g2,t1)*zb(g2,t1))*
     &           zb(u,b1)) - 
     &        g1g2t*(za(d,t1)*za(g1,g2) - za(d,g2)*za(g1,t1))*
     &         zb(g2,b1)*(za(g1,b1)*zb(u,b1) + za(g1,g2)*zb(u,g2))
     &        ) + s(g1,t1)*
     &      (g1g2b*g2b*za(g1,b1)*za(g1,t1)*
     &         (za(d,t1)*za(g1,g2) - za(d,g2)*za(g1,t1))*
     &         zb(g1,b1)*zb(g1,t1)*zb(g2,t1)*zb(u,b1) + 
     &        za(d,g1)*(-(g1g2t*za(g1,b1)**2*za(g1,t1)*
     &              zb(g1,b1)**2*zb(g1,t1)*
     &              ((dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,b1) - 
     &                za(g1,g2)*zb(g2,b1)*zb(u,g1))) - 
     &           g1g2t*ma2*za(g1,g2)*za(g1,t1)*zb(g1,t1)*
     &            ((-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*zb(g2,b1)*
     &               zb(u,g1) + g2b*zb(g1,b1)*zb(u,g2)) + 
     &           za(g1,b1)*zb(g1,b1)*
     &            (g2b*za(g1,t1)*zb(g1,t1)*
     &               (g1g2b*za(g1,t1)*zb(g1,t1) + 
     &                 g1g2t*za(g2,b1)*zb(g2,b1))*zb(u,b1) + 
     &              g1g2t*za(g1,g2)**2*za(g1,t1)*zb(g1,g2)*
     &               zb(g1,t1)*zb(g2,b1)*zb(u,g1) + 
     &              za(g1,g2)*
     &               (zb(g1,g2)*
     &                  (dcmplx(g1g2b*g2b*mq2) - 
     &                    za(g1,t1)*zb(g1,t1)*
     &                   (g1g2b*dcmplx(g2b)+ g1g2t*za(g2,b1)*zb(g2,b1))
     &                    )*zb(u,b1) - 
     &                 g1g2t*g2b*za(g1,t1)*zb(g1,t1)*
     &                  (zb(g2,b1)*zb(u,g1) + zb(g1,b1)*zb(u,g2)))
     &              )))))/
     & (2d0*za(g1,b1)*za(g1,g2)**2*za(g1,t1)**2*zb(g1,b1)*
     &   zb(g1,t1))
      amp(1,1,1,2)=
     &        (mq*(-2*g1g2b*zb(g1,t1)*
     &      (za(d,t1)*za(g1,g2)*
     &         (-(g1g2t*za(g2,b1)*zb(g1,b1)) + 
     &           g2b*za(g2,t1)*zb(g1,t1)) + 
     &        za(d,g2)*(g2b*za(g1,g2)*za(g2,t1)*zb(g1,g2) + 
     &           za(g1,t1)*
     &            (g1g2t*za(g2,b1)*zb(g1,b1) - 
     &              g2b*za(g2,t1)*zb(g1,t1))))*zb(u,b1) + 
     &     g1g2t*s(g1,t1)*za(d,g1)*za(g2,b1)*zb(g1,b1)**2*
     &      (za(g2,b1)*zb(u,b1) - za(g1,g2)*zb(u,g1))))/
     & (2d0*za(g1,g2)*za(g1,t1)*zb(g1,g2))
      amp(1,1,2,1)=
     &        -(ma*mq*(2*g1g2b*za(g1,b1)*za(g1,t1)*zb(g1,b1)*
     &       zb(g1,t1)**2*(-(g1g2t*za(d,t1)*za(g1,g2)*
     &            (za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &              za(g1,b1)*
     &               (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2))))
     &          + za(g1,t1)*
     &          (g2b*za(d,g1)*
     &             (za(g1,g2)*zb(g1,g2) + za(g2,t1)*zb(g2,t1))*
     &             zb(u,g1) + 
     &            g1g2t*za(d,g2)*
     &             (za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &               za(g1,b1)*
     &                (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2)))))
     &        + s(g1,t1)*(g1g2b*g2b*za(g1,b1)*za(g1,t1)*
     &          (-(za(d,t1)*za(g1,g2)) + za(d,g2)*za(g1,t1))*
     &          zb(g1,b1)*zb(g1,t1)*zb(g2,t1)*zb(u,g1) + 
     &         za(d,g1)*(g1g2t*ma2*za(g1,g2)*za(g1,t1)*zb(g1,g2)*
     &            (dcmplx(g2b) - za(g1,g2)*zb(g1,g2))*zb(g1,t1)*zb(u,g1)
     &             + g1g2t*za(g1,b1)**2*za(g1,t1)*zb(g1,b1)**2*
     &            zb(g1,t1)*(dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,g1)
     &             - za(g1,b1)*zb(g1,b1)*
     &             (g1g2b*g2b*mq2*za(g1,g2)*zb(g1,g2)*zb(u,g1) + 
     &               g1g2b*g2b*za(g1,t1)**2*zb(g1,t1)**2*
     &                zb(u,g1) + 
     &               za(g1,t1)*zb(g1,t1)*
     &                (-(za(g1,g2)*zb(g1,g2)*
     &                     ((g1g2b - g1g2t)*dcmplx(g2b)+ 
     &                      g1g2t*za(g1,g2)*zb(g1,g2))*zb(u,g1))
     &                   + g1g2t*za(g2,b1)*
     &                  (zb(g1,g2)*(-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*
     &                      zb(u,b1) + g2b*zb(g1,b1)*zb(u,g2))))))
     &      ))/
     & (2d0*za(g1,b1)*za(g1,g2)**2*za(g1,t1)**2*zb(g1,b1)**2*
     &   zb(g1,t1))
      amp(1,1,2,2)=
     &        (ma*mq*(g1g2t*s(g1,t1)*za(d,g1)*za(g2,b1)**2*
     &      zb(g1,b1)**2 - 
     &     2*g1g2b*zb(g1,t1)*
     &      (za(d,t1)*za(g1,g2)*
     &         (-(g1g2t*za(g2,b1)*zb(g1,b1)) + 
     &           g2b*za(g2,t1)*zb(g1,t1)) + 
     &        za(d,g2)*(g2b*za(g1,g2)*za(g2,t1)*zb(g1,g2) + 
     &           za(g1,t1)*
     &            (g1g2t*za(g2,b1)*zb(g1,b1) - 
     &              g2b*za(g2,t1)*zb(g1,t1)))))*zb(u,g1))/
     & (2d0*za(g1,g2)*za(g1,t1)*zb(g1,b1)*zb(g1,g2))
      amp(1,2,1,1)=
     &        (2*g1g2b*za(g1,b1)*za(g1,t1)*za(g2,t1)*zb(g1,b1)*
     &    zb(g1,t1)**2*(g2b*za(d,g1)*za(g1,t1)*zb(g1,g2)*
     &       zb(u,b1) + za(d,t1)*
     &       ((g1g2t*za(g1,b1)*zb(g2,b1) - 
     &            g2b*za(g1,t1)*zb(g2,t1))*zb(u,b1) + 
     &         g1g2t*za(g1,g2)*zb(g2,b1)*zb(u,g2))) + 
     &   s(g1,t1)*(g1g2b*g2b*za(d,g1)*za(g1,b1)*za(g2,t1)*
     &       zb(g1,b1)*zb(g1,g2)*(-dcmplx(mq2) + za(g1,t1)*zb(g1,t1))*
     &       zb(u,b1) + g1g2b*g2b*za(d,g2)*za(g1,b1)*za(g1,t1)*
     &       zb(g1,b1)*zb(g1,g2)*(dcmplx(mq2) + za(g1,t1)*zb(g1,t1))*
     &       zb(u,b1) - za(d,t1)*za(g1,t1)*zb(g1,t1)*
     &       (g1g2t*za(g1,b1)**2*zb(g1,b1)**2*
     &          ((dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,b1) - 
     &            za(g1,g2)*zb(g2,b1)*zb(u,g1)) + 
     &         g1g2t*ma2*za(g1,g2)*
     &          ((-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*zb(g2,b1)*
     &             zb(u,g1) + g2b*zb(g1,b1)*zb(u,g2)) + 
     &         za(g1,b1)*zb(g1,b1)*
     &          (-(g1g2b*g2b*za(g1,t1)*zb(g1,t1)*zb(u,b1)) + 
     &            g1g2t*za(g2,b1)*(-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*
     &             zb(g2,b1)*zb(u,b1) + 
     &            g1g2b*g2b*za(g2,t1)*zb(g2,t1)*zb(u,b1) + 
     &            g1g2t*za(g1,g2)*
     &             ((dcmplx(g2b) - za(g1,g2)*zb(g1,g2))*zb(g2,b1)*
     &                zb(u,g1) + g2b*zb(g1,b1)*zb(u,g2))))))/
     & (2d0*za(g1,b1)*za(g1,g2)**2*za(g1,t1)*zb(g1,b1)*zb(g1,t1))
      amp(1,2,1,2)=
     &        (2*g1g2b*za(g2,t1)*zb(g1,t1)*
     &    (g2b*za(d,g2)*za(g2,t1)*zb(g1,g2) + 
     &      za(d,t1)*(-(g1g2t*za(g2,b1)*zb(g1,b1)) + 
     &         g2b*za(g2,t1)*zb(g1,t1)))*zb(u,b1) + 
     &   g1g2t*s(g1,t1)*za(d,t1)*za(g2,b1)*zb(g1,b1)**2*
     &    (za(g2,b1)*zb(u,b1) - za(g1,g2)*zb(u,g1)))/
     & (2d0*za(g1,g2)*zb(g1,g2))
      amp(1,2,2,1)=
     &        (ma*(2*g1g2b*za(g1,b1)*za(g1,t1)*za(g2,t1)*zb(g1,b1)*
     &      zb(g1,t1)**2*(g2b*za(d,g1)*za(g1,t1)*zb(g1,g2)*
     &         zb(u,g1) - za(d,t1)*
     &         (g2b*za(g1,t1)*zb(g2,t1)*zb(u,g1) + 
     &           g1g2t*za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &           g1g2t*za(g1,b1)*
     &            (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2)))) + 
     &     s(g1,t1)*(g1g2b*g2b*za(g1,b1)*zb(g1,b1)*zb(g1,g2)*
     &       (za(d,g1)*za(g2,t1)*(-dcmplx(mq2) + za(g1,t1)*zb(g1,t1)) + 
     &         za(d,g2)*za(g1,t1)*(dcmplx(mq2) + za(g1,t1)*zb(g1,t1)))*
     &         zb(u,g1) + za(d,t1)*za(g1,t1)*zb(g1,t1)*
     &         (g1g2t*ma2*za(g1,g2)*zb(g1,g2)*
     &            (-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*zb(u,g1) - 
     &           g1g2t*za(g1,b1)**2*zb(g1,b1)**2*
     &            (dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,g1) + 
     &           za(g1,b1)*zb(g1,b1)*
     &            ((g1g2t*g2b*za(g1,g2)*zb(g1,g2) - 
     &                 g1g2t*za(g1,g2)**2*zb(g1,g2)**2 + 
     &                 g1g2b*g2b*
     &                  (za(g1,t1)*zb(g1,t1) - 
     &                    za(g2,t1)*zb(g2,t1)))*zb(u,g1) + 
     &              g1g2t*za(g2,b1)*
     &               (zb(g1,g2)*(-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*
     &                  zb(u,b1) + g2b*zb(g1,b1)*zb(u,g2)))))))/
     & (2d0*za(g1,b1)*za(g1,g2)**2*za(g1,t1)*zb(g1,b1)**2*
     &   zb(g1,t1))
      amp(1,2,2,2)=
     &        (ma*(g1g2t*s(g1,t1)*za(d,t1)*za(g2,b1)**2*zb(g1,b1)**2 + 
     &     2*g1g2b*za(g2,t1)*zb(g1,t1)*
     &      (g2b*za(d,g2)*za(g2,t1)*zb(g1,g2) + 
     &        za(d,t1)*(-(g1g2t*za(g2,b1)*zb(g1,b1)) + 
     &           g2b*za(g2,t1)*zb(g1,t1))))*zb(u,g1))/
     & (2d0*za(g1,g2)*zb(g1,b1)*zb(g1,g2))
      amp(2,1,1,1)=
     &        (mq*za(d,g1)*(g1g2t*s(g1,t1)*za(g1,b1)*zb(g2,b1)**2*
     &      (za(g1,b1)*zb(u,b1) + za(g1,g2)*zb(u,g2)) + 
     &     2*g1g2b*za(g1,t1)*zb(g2,t1)*
     &      (-(g1g2t*za(g1,b1)*zb(g2,b1)*zb(u,b1)) + 
     &        g2b*za(g1,t1)*zb(g2,t1)*zb(u,b1) - 
     &        g1g2t*za(g1,g2)*zb(g2,b1)*zb(u,g2))))/
     & (2d0*za(g1,g2)*za(g1,t1)*zb(g1,g2))
      amp(2,1,1,2)=
     &        (mq*(2*g1g2b*za(g1,b1)*za(g1,t1)**2*zb(g1,b1)*zb(g1,t1)*
     &      (g1g2t*za(d,g1)*za(g2,b1)*zb(g1,b1) + 
     &        g2b*za(d,g2)*za(g1,g2)*zb(g1,g2) + 
     &        g2b*(za(d,t1)*za(g1,g2) - za(d,g2)*za(g1,t1))*
     &         zb(g1,t1))*zb(g2,t1)*zb(u,b1) + 
     &     s(g1,t1)*(g1g2b*g2b*za(g1,b1)*za(g1,t1)*
     &         (za(d,t1)*za(g1,g2) - za(d,g2)*za(g1,t1))*
     &         zb(g1,b1)*zb(g1,t1)*zb(g2,t1)*zb(u,b1) + 
     &        za(d,g1)*(-(g1g2t*za(g1,b1)**2*za(g1,t1)*
     &              zb(g1,b1)**2*zb(g1,t1)*
     &              (dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,b1)) + 
     &           g1g2t*ma2*za(g1,g2)*za(g1,t1)*zb(g1,t1)*
     &            (g2b*zb(g2,b1)*zb(u,g1) + 
     &              zb(g1,b1)*(-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*
     &               zb(u,g2)) - 
     &           za(g1,b1)*zb(g1,b1)*
     &            (-(g2b*za(g1,t1)*zb(g1,t1)*
     &                 (g1g2b*za(g1,t1)*zb(g1,t1) + 
     &                   g1g2t*za(g2,b1)*zb(g2,b1))*zb(u,b1)) + 
     &              za(g1,g2)*
     &               (g1g2b*g2b*zb(g1,g2)*
     &                 (-dcmplx(mq2) + za(g1,t1)*zb(g1,t1))*zb(u,b1) + 
     &                 g1g2t*za(g1,t1)*zb(g1,t1)*
     &                  (g2b*zb(g2,b1)*zb(u,g1) + 
     &                    zb(g1,b1)*(dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*
     &                     zb(u,g2))))))))/
     & (2d0*za(g1,b1)*za(g1,t1)**2*zb(g1,b1)*zb(g1,g2)**2*
     &   zb(g1,t1))
      amp(2,1,2,1)=
     &        (ma*mq*za(d,g1)*(-(g1g2t*s(g1,t1)*za(g1,b1)*zb(g2,b1)*
     &        (za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &          za(g1,b1)*(zb(g1,g2)*zb(u,b1) - 
     &             zb(g1,b1)*zb(u,g2)))) + 
     &     2*g1g2b*za(g1,t1)*zb(g2,t1)*
     &      (g2b*za(g1,t1)*zb(g2,t1)*zb(u,g1) + 
     &        g1g2t*za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &        g1g2t*za(g1,b1)*
     &         (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2)))))/
     & (2d0*za(g1,g2)*za(g1,t1)*zb(g1,b1)*zb(g1,g2))
      amp(2,1,2,2)=
     &        (ma*mq*(2*g1g2b*za(g1,b1)*za(g1,t1)**2*zb(g1,b1)*
     &      zb(g1,t1)*(g1g2t*za(d,g1)*za(g2,b1)*zb(g1,b1) + 
     &        g2b*za(d,g2)*za(g1,g2)*zb(g1,g2) + 
     &        g2b*(za(d,t1)*za(g1,g2) - za(d,g2)*za(g1,t1))*
     &         zb(g1,t1))*zb(g2,t1)*zb(u,g1) + 
     &     s(g1,t1)*(g1g2b*g2b*za(g1,b1)*za(g1,t1)*
     &         (za(d,t1)*za(g1,g2) - za(d,g2)*za(g1,t1))*
     &         zb(g1,b1)*zb(g1,t1)*zb(g2,t1)*zb(u,g1) + 
     &        za(d,g1)*(-(g1g2t*g2b*ma2*za(g1,g2)*za(g1,t1)*
     &              zb(g1,g2)*zb(g1,t1)*zb(u,g1)) - 
     &           g1g2t*za(g1,b1)**2*za(g1,t1)*zb(g1,b1)**2*
     &            zb(g1,t1)*
     &            (g2b*zb(u,g1) + 
     &              za(g2,b1)*
     &               (-(zb(g1,g2)*zb(u,b1)) + zb(g1,b1)*zb(u,g2)))
     &             + g2b*za(g1,b1)*zb(g1,b1)*
     &            (g1g2b*mq2*za(g1,g2)*zb(g1,g2)*zb(u,g1) + 
     &              g1g2b*za(g1,t1)**2*zb(g1,t1)**2*zb(u,g1) + 
     &              za(g1,t1)*zb(g1,t1)*
     &               ((-g1g2b + g1g2t)*za(g1,g2)*zb(g1,g2)*
     &                  zb(u,g1) + 
     &                 g1g2t*za(g2,b1)*
     &                  (-(zb(g1,g2)*zb(u,b1)) + 
     &                    zb(g1,b1)*zb(u,g2))))))))/
     & (2d0*za(g1,b1)*za(g1,t1)**2*zb(g1,b1)**2*zb(g1,g2)**2*
     &   zb(g1,t1))
      amp(2,2,1,1)=
     &        (2*g1g2b*za(d,g1)*za(g1,t1)*zb(g1,g2)*
     &    ((g1g2t*za(g1,b1)*zb(g2,b1) - g2b*za(g1,t1)*zb(g2,t1))*
     &       zb(u,b1) + g1g2t*za(g1,g2)*zb(g2,b1)*zb(u,g2)) + 
     &   za(d,t1)*(g1g2t*s(g1,t1)*za(g1,b1)*zb(g2,b1)**2*
     &       (za(g1,b1)*zb(u,b1) + za(g1,g2)*zb(u,g2)) + 
     &      2*g1g2b*za(g1,t1)*zb(g2,t1)*
     &       (-(g1g2t*za(g1,b1)*zb(g2,b1)*zb(u,b1)) + 
     &         g2b*za(g1,t1)*zb(g2,t1)*zb(u,b1) - 
     &         g1g2t*za(g1,g2)*zb(g2,b1)*zb(u,g2))))/
     & (2d0*za(g1,g2)*zb(g1,g2))
      amp(2,2,1,2)=
     &        (-2*g1g2b*za(g1,b1)*za(g1,t1)**2*zb(g1,b1)*zb(g1,t1)*
     &    (g1g2t*za(d,g1)*za(g2,b1)*zb(g1,b1)*zb(g1,g2) + 
     &      g2b*za(d,g2)*zb(g1,g2)*
     &       (za(g1,g2)*zb(g1,g2) + za(g2,t1)*zb(g2,t1)) + 
     &      za(d,t1)*(-(g1g2t*za(g2,b1)*zb(g1,b1)*zb(g2,t1)) + 
     &         g2b*zb(g1,t1)*
     &          (za(g1,g2)*zb(g1,g2) + za(g2,t1)*zb(g2,t1))))*
     &    zb(u,b1) + s(g1,t1)*
     &    (g1g2b*g2b*za(d,g1)*za(g1,b1)*za(g2,t1)*zb(g1,b1)*
     &       zb(g1,g2)*(-dcmplx(mq2) + za(g1,t1)*zb(g1,t1))*zb(u,b1) + 
     &      g1g2b*g2b*za(d,g2)*za(g1,b1)*za(g1,t1)*zb(g1,b1)*
     &       zb(g1,g2)*(dcmplx(mq2) + za(g1,t1)*zb(g1,t1))*zb(u,b1) + 
     &      za(d,t1)*za(g1,t1)*zb(g1,t1)*
     &       (-(g1g2t*za(g1,b1)**2*zb(g1,b1)**2*
     &            (dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,b1)) + 
     &         g1g2t*ma2*za(g1,g2)*
     &          (g2b*zb(g2,b1)*zb(u,g1) + 
     &          zb(g1,b1)*(-dcmplx(g2b)+ za(g1,g2)*zb(g1,g2))*zb(u,g2))
     &           - za(g1,b1)*zb(g1,b1)*
     &          (-(g2b*(g1g2b*za(g1,t1)*zb(g1,t1) + 
     &                 g1g2t*za(g2,b1)*zb(g2,b1) - 
     &                 g1g2b*za(g2,t1)*zb(g2,t1))*zb(u,b1)) + 
     &            g1g2t*g2b*za(g1,g2)*zb(g2,b1)*zb(u,g1) + 
     &            g1g2t*za(g1,g2)*zb(g1,b1)*
     &             (dcmplx(g2b)+ za(g2,b1)*zb(g2,b1))*zb(u,g2)))))/
     & (2d0*za(g1,b1)*za(g1,t1)*zb(g1,b1)*zb(g1,g2)**2*zb(g1,t1))
      amp(2,2,2,1)=
     &        -(ma*(2*g1g2b*za(d,g1)*za(g1,t1)*zb(g1,g2)*
     &       (g2b*za(g1,t1)*zb(g2,t1)*zb(u,g1) + 
     &         g1g2t*za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &         g1g2t*za(g1,b1)*
     &          (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2))) + 
     &      za(d,t1)*(g1g2t*s(g1,t1)*za(g1,b1)*zb(g2,b1)*
     &          (za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &            za(g1,b1)*
     &             (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2))) - 
     &         2*g1g2b*za(g1,t1)*zb(g2,t1)*
     &          (g2b*za(g1,t1)*zb(g2,t1)*zb(u,g1) + 
     &            g1g2t*za(g1,g2)*zb(g1,g2)*zb(u,g2) + 
     &            g1g2t*za(g1,b1)*
     &             (zb(g1,g2)*zb(u,b1) - zb(g1,b1)*zb(u,g2))))))/
     & (2d0*za(g1,g2)*zb(g1,b1)*zb(g1,g2))
      amp(2,2,2,2)=
     &        -(ma*(2*g1g2b*za(g1,b1)*za(g1,t1)**2*zb(g1,b1)*zb(g1,t1)*
     &       (g1g2t*za(d,g1)*za(g2,b1)*zb(g1,b1)*zb(g1,g2) + 
     &         g2b*za(d,g2)*zb(g1,g2)*
     &          (za(g1,g2)*zb(g1,g2) + za(g2,t1)*zb(g2,t1)) + 
     &         za(d,t1)*(-(g1g2t*za(g2,b1)*zb(g1,b1)*zb(g2,t1)) + 
     &            g2b*zb(g1,t1)*
     &             (za(g1,g2)*zb(g1,g2) + za(g2,t1)*zb(g2,t1))))*
     &       zb(u,g1) + s(g1,t1)*
     &       (-(g1g2b*g2b*za(g1,b1)*zb(g1,b1)*zb(g1,g2)*
     &            (za(d,g1)*za(g2,t1)*
     &               (-dcmplx(mq2) + za(g1,t1)*zb(g1,t1)) + 
     &          za(d,g2)*za(g1,t1)*(dcmplx(mq2) + za(g1,t1)*zb(g1,t1))
     &              )*zb(u,g1)) + 
     &         za(d,t1)*za(g1,t1)*zb(g1,t1)*
     &          (g1g2t*g2b*ma2*za(g1,g2)*zb(g1,g2)*zb(u,g1) + 
     &            g1g2t*za(g1,b1)**2*zb(g1,b1)**2*
     &             (g2b*zb(u,g1) + 
     &               za(g2,b1)*
     &                (-(zb(g1,g2)*zb(u,b1)) + zb(g1,b1)*zb(u,g2))
     &               ) - g2b*za(g1,b1)*zb(g1,b1)*
     &             ((g1g2t*za(g1,g2)*zb(g1,g2) + 
     &                  g1g2b*za(g1,t1)*zb(g1,t1) - 
     &                  g1g2b*za(g2,t1)*zb(g2,t1))*zb(u,g1) + 
     &               g1g2t*za(g2,b1)*
     &                (-(zb(g1,g2)*zb(u,b1)) + zb(g1,b1)*zb(u,g2))
     &               )))))/
     & (2d0*za(g1,b1)*za(g1,t1)*zb(g1,b1)**2*zb(g1,g2)**2*
     &   zb(g1,t1))
      return
      end
