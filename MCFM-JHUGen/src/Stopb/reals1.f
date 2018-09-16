      subroutine reals1(u,g1,d,t1,b1,g2, p,mq,ma,za,zb,amp)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
* AUTHORS: R. FREDERIX AND F. TRAMONTANO                               *
* DATE  : 7/16/2008                                                    *
*                                                                      *
************************************************************************
c--- 0  ->  u~ + g1 + d + t1 + b1~ + g2  (t-channel single-top)
c--- where the g1 (g2) is the gluon attached to the massive (massless) quark
c--- line. Arguments of amp are helicities of gluon1, top, bottom and gluon2,
c--- respectively. Spin of top and bottom quarks should be projected on light-
c--- like vector p. mq is the mass of the top quark and ma the mass of the 
c--- bottom quark.
 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: u,g1,g2,d,b1,t1,p
      complex(dp):: amp(2,2,2,2)
      real(dp):: mq,ma,mq2,ma2,g1t,g1b
 
 
      mq2 = mq**2
      ma2 = ma**2
      g1t = s(g1,t1)/2._dp + s(g1,p)*mq2/(2._dp*s(t1,p))
      g1b = s(g1,b1)/2._dp + s(g1,p)*ma2/(2._dp*s(b1,p))
 
      amp(1,1,1,1)=
     &        (2*mq*(g1b*za(b1,p)*zb(b1,p)*
     &      (mq2*za(d,p)*zb(g1,p) - 
     &        za(d,t1)*za(t1,p)*zb(g1,t1)*zb(t1,p))*
     &      (za(d,g2)*zb(g2,b1) - za(u,d)*zb(u,b1)) + 
     &     za(t1,p)*zb(t1,p)*
     &      (za(d,g2)*(g1b*za(d,b1)*za(t1,p)*zb(b1,p)*zb(g1,t1)*
     &            zb(g2,b1) + 
     &           g1t*za(d,p)*
     &            (zb(b1,p)*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*
     &               zb(g1,g2) - ma2*zb(g1,b1)*zb(g2,p))) + 
     &        za(u,d)*(-(g1b*za(d,b1)*za(t1,p)*zb(b1,p)*zb(g1,t1)*
     &              zb(u,b1)) + 
     &           g1t*za(d,p)*
     &            (zb(b1,p)*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*
     &               zb(u,g1) + ma2*zb(g1,b1)*zb(u,p))))))/
     & (g1b*g1t*za(d,g2)*za(g1,b1)*za(t1,p)**2*za(u,g2)*zb(b1,p)*
     &   zb(t1,p))
      amp(1,1,1,2)=
     &        (2*mq*(g1b*za(b1,p)*zb(b1,p)*zb(u,b1)*
     &      (mq2*za(d,p)*zb(g1,p)*zb(u,d) + 
     &        mq2*za(g2,p)*zb(g1,p)*zb(u,g2) - 
     &        za(t1,p)*zb(g1,t1)*zb(t1,p)*
     &         (za(d,t1)*zb(u,d) + za(g2,t1)*zb(u,g2))) + 
     &     za(t1,p)*zb(t1,p)*
     &      (g1b*za(d,b1)*za(t1,p)*zb(b1,p)*zb(g1,t1)*zb(u,b1)*
     &         zb(u,d) - g1t*za(d,p)*zb(u,d)*
     &         (zb(b1,p)*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*zb(u,g1) + 
     &           ma2*zb(g1,b1)*zb(u,p)) + 
     &        zb(u,g2)*(g1b*za(g2,b1)*za(t1,p)*zb(b1,p)*zb(g1,t1)*
     &            zb(u,b1) - 
     &           g1t*za(g2,p)*
     &            (zb(b1,p)*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*
     &               zb(u,g1) + ma2*zb(g1,b1)*zb(u,p))))))/
     & (g1b*g1t*za(g1,b1)*za(t1,p)**2*zb(b1,p)*zb(d,g2)*zb(t1,p)*
     &   zb(u,g2))
      amp(1,1,2,1)=
     &        (-2*ma*mq*(g1b*za(b1,p)*zb(b1,p)*
     &      (mq2*za(d,p)*zb(g1,p) - 
     &        za(d,t1)*za(t1,p)*zb(g1,t1)*zb(t1,p))*
     &      (za(d,g2)*zb(g2,p) - za(u,d)*zb(u,p)) + 
     &     za(t1,p)*zb(t1,p)*
     &      (za(d,g2)*(g1b*za(d,b1)*za(t1,p)*zb(b1,p)*zb(g1,t1)*
     &            zb(g2,p) + 
     &           g1t*za(d,p)*zb(g1,p)*
     &            (za(g1,b1)*zb(b1,p)*zb(g1,g2) - ma2*zb(g2,p)))
     &         + za(u,d)*(-(g1b*za(d,b1)*za(t1,p)*zb(b1,p)*
     &              zb(g1,t1)*zb(u,p)) + 
     &           g1t*za(d,p)*zb(g1,p)*
     &            (za(g1,b1)*zb(b1,p)*zb(u,g1) + ma2*zb(u,p))))))/
     & (g1b*g1t*za(d,g2)*za(g1,b1)*za(t1,p)**2*za(u,g2)*
     &   zb(b1,p)**2*zb(t1,p))
      amp(1,1,2,2)=
     &        (2*ma*mq*(g1t*za(g1,b1)*za(g2,p)*za(t1,p)*zb(b1,p)*
     &      zb(g1,p)*zb(t1,p)*zb(u,g1)*zb(u,g2) + 
     &     (-(za(t1,p)*zb(t1,p)*
     &           (-(g1t*ma2*za(g2,p)*zb(g1,p)*zb(u,g2)) + 
     &             g1b*za(t1,p)*zb(b1,p)*zb(g1,t1)*
     &              (za(d,b1)*zb(u,d) + za(g2,b1)*zb(u,g2)))) + 
     &        g1b*za(b1,p)*zb(b1,p)*
     &         (-(mq2*za(g2,p)*zb(g1,p)*zb(u,g2)) + 
     &           za(t1,p)*zb(g1,t1)*zb(t1,p)*
     &            (za(d,t1)*zb(u,d) + za(g2,t1)*zb(u,g2))))*
     &      zb(u,p) + za(d,p)*zb(g1,p)*zb(u,d)*
     &      (g1t*za(g1,b1)*za(t1,p)*zb(b1,p)*zb(t1,p)*zb(u,g1) + 
     &        (-(g1b*mq2*za(b1,p)*zb(b1,p)) + 
     &           g1t*ma2*za(t1,p)*zb(t1,p))*zb(u,p))))/
     & (g1b*g1t*za(g1,b1)*za(t1,p)**2*zb(b1,p)**2*zb(d,g2)*
     &   zb(t1,p)*zb(u,g2))
      amp(1,2,1,1)=
     &        (-2*(za(d,g2)*(g1b*mq2*
     &         (-(za(d,p)*za(t1,b1)) + za(d,b1)*za(t1,p))*
     &         zb(b1,p)*zb(g1,p)*zb(g2,b1) + 
     &        za(d,t1)*za(t1,p)*
     &         (zb(b1,p)*(g1t*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*
     &               zb(g1,g2) + g1b*za(t1,b1)*zb(g1,t1)*zb(g2,b1)
     &              ) - g1t*ma2*zb(g1,b1)*zb(g2,p))*zb(t1,p)) + 
     &     za(u,d)*(g1b*mq2*za(d,p)*za(t1,b1)*zb(b1,p)*zb(g1,p)*
     &         zb(u,b1) + za(t1,p)*
     &         (-(g1b*mq2*za(d,b1)*zb(b1,p)*zb(g1,p)*zb(u,b1)) + 
     &           za(d,t1)*zb(t1,p)*
     &            (-(g1b*za(t1,b1)*zb(b1,p)*zb(g1,t1)*zb(u,b1)) + 
     &              g1t*(zb(b1,p)*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*
     &                  zb(u,g1) + ma2*zb(g1,b1)*zb(u,p)))))))/
     & (g1b*g1t*za(d,g2)*za(g1,b1)*za(t1,p)*za(u,g2)*zb(b1,p)*
     &   zb(t1,p))
      amp(1,2,1,2)=
     &        (2*zb(b1,p)*(zb(u,d)*
     &       (g1b*mq2*za(d,p)*za(t1,b1)*zb(g1,p)*zb(u,b1) + 
     &         za(t1,p)*(-(g1b*mq2*za(d,b1)*zb(g1,p)*zb(u,b1)) + 
     &            za(d,t1)*zb(t1,p)*
     &             (-(g1b*za(t1,b1)*zb(g1,t1)*zb(u,b1)) + 
     &               g1t*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*zb(u,g1))))
     &       + (g1b*mq2*za(g2,p)*za(t1,b1)*zb(g1,p)*zb(u,b1) + 
     &         za(t1,p)*(-(g1b*mq2*za(g2,b1)*zb(g1,p)*zb(u,b1)) + 
     &            za(g2,t1)*zb(t1,p)*
     &             (-(g1b*za(t1,b1)*zb(g1,t1)*zb(u,b1)) + 
     &               g1t*(cplx1(ma2)+ za(g1,b1)*zb(g1,b1))*zb(u,g1))))*
     &       zb(u,g2)) + 2*g1t*ma2*za(t1,p)*zb(g1,b1)*zb(t1,p)*
     &    (za(d,t1)*zb(u,d) + za(g2,t1)*zb(u,g2))*zb(u,p))/
     & (g1b*g1t*za(g1,b1)*za(t1,p)*zb(b1,p)*zb(d,g2)*zb(t1,p)*
     &   zb(u,g2))
      amp(1,2,2,1)=
     &        (2*ma*(za(d,g2)*(g1b*mq2*
     &         (-(za(d,p)*za(t1,b1)) + za(d,b1)*za(t1,p))*
     &         zb(b1,p)*zb(g1,p)*zb(g2,p) + 
     &        za(d,t1)*za(t1,p)*
     &         (g1t*za(g1,b1)*zb(b1,p)*zb(g1,g2)*zb(g1,p) + 
     &           (-(g1t*ma2*zb(g1,p)) + 
     &              g1b*za(t1,b1)*zb(b1,p)*zb(g1,t1))*zb(g2,p))*
     &         zb(t1,p)) + 
     &     za(u,d)*(g1b*mq2*
     &         (za(d,p)*za(t1,b1) - za(d,b1)*za(t1,p))*zb(b1,p)*
     &         zb(g1,p)*zb(u,p) + 
     &        za(d,t1)*za(t1,p)*zb(t1,p)*
     &         (g1t*za(g1,b1)*zb(b1,p)*zb(g1,p)*zb(u,g1) + 
     &           (g1t*ma2*zb(g1,p) - 
     &              g1b*za(t1,b1)*zb(b1,p)*zb(g1,t1))*zb(u,p)))))/
     & (g1b*g1t*za(d,g2)*za(g1,b1)*za(t1,p)*za(u,g2)*zb(b1,p)**2*
     &   zb(t1,p))
      amp(1,2,2,2)=
     &        (-2*ma*(g1t*za(g1,b1)*za(g2,t1)*za(t1,p)*zb(b1,p)*
     &      zb(g1,p)*zb(t1,p)*zb(u,g1)*zb(u,g2) + 
     &     (g1b*mq2*(za(d,p)*za(t1,b1) - za(d,b1)*za(t1,p))*
     &         zb(b1,p)*zb(g1,p)*zb(u,d) + 
     &        (g1b*mq2*(za(g2,p)*za(t1,b1) - za(g2,b1)*za(t1,p))*
     &            zb(b1,p)*zb(g1,p) + 
     &           za(g2,t1)*za(t1,p)*
     &            (g1t*ma2*zb(g1,p) - 
     &              g1b*za(t1,b1)*zb(b1,p)*zb(g1,t1))*zb(t1,p))*
     &         zb(u,g2))*zb(u,p) + 
     &     za(d,t1)*za(t1,p)*zb(t1,p)*zb(u,d)*
     &      (g1t*za(g1,b1)*zb(b1,p)*zb(g1,p)*zb(u,g1) + 
     &        (g1t*ma2*zb(g1,p) - 
     &           g1b*za(t1,b1)*zb(b1,p)*zb(g1,t1))*zb(u,p))))/
     & (g1b*g1t*za(g1,b1)*za(t1,p)*zb(b1,p)**2*zb(d,g2)*zb(t1,p)*
     &   zb(u,g2))
      amp(2,1,1,1)=
     &        (-2*mq*(-(g1t*ma2*za(d,p)*za(g1,p)*za(t1,p)*zb(t1,p)) + 
     &     g1b*za(b1,p)*(mq2*za(d,p)*za(g1,p)*zb(b1,p) + 
     &        za(t1,p)*(za(d,t1)*za(g1,p)*zb(t1,b1) + 
     &           za(d,g1)*(za(g1,p)*zb(g1,b1) - 
     &              za(t1,p)*zb(t1,b1)))*zb(t1,p)))*
     &   (za(d,g2)*zb(g2,b1) - za(u,d)*zb(u,b1)))/
     & (g1b*g1t*za(b1,p)*za(d,g2)*za(t1,p)**2*za(u,g2)*zb(g1,b1)*
     &   zb(t1,p))
      amp(2,1,1,2)=
     &        (-2*mq*zb(u,b1)*(-(g1t*ma2*za(g1,p)*za(t1,p)*zb(t1,p)*
     &        (za(d,p)*zb(u,d) + za(g2,p)*zb(u,g2))) + 
     &     g1b*za(b1,p)*((mq2*za(d,p)*za(g1,p)*zb(b1,p) + 
     &           za(t1,p)*(za(d,t1)*za(g1,p)*zb(t1,b1) + 
     &              za(d,g1)*
     &               (za(g1,p)*zb(g1,b1) - za(t1,p)*zb(t1,b1)))*
     &            zb(t1,p))*zb(u,d) + 
     &        (mq2*za(g1,p)*za(g2,p)*zb(b1,p) + 
     &           za(t1,p)*(za(g1,p)*za(g2,t1)*zb(t1,b1) + 
     &              za(g1,g2)*
     &               (-(za(g1,p)*zb(g1,b1)) + za(t1,p)*zb(t1,b1)))
     &             *zb(t1,p))*zb(u,g2))))/
     & (g1b*g1t*za(b1,p)*za(t1,p)**2*zb(d,g2)*zb(g1,b1)*zb(t1,p)*
     &   zb(u,g2))
      amp(2,1,2,1)=
     &        (2*ma*mq*(-(g1t*ma2*za(d,p)*za(g1,p)*za(t1,p)*zb(t1,p)) + 
     &     g1b*za(b1,p)*(mq2*za(d,p)*za(g1,p)*zb(b1,p) + 
     &        za(t1,p)*(za(d,t1)*za(g1,p)*zb(t1,b1) + 
     &           za(d,g1)*(za(g1,p)*zb(g1,b1) - 
     &              za(t1,p)*zb(t1,b1)))*zb(t1,p)))*
     &   (za(d,g2)*zb(g2,p) - za(u,d)*zb(u,p)))/
     & (g1b*g1t*za(b1,p)*za(d,g2)*za(t1,p)**2*za(u,g2)*zb(b1,p)*
     &   zb(g1,b1)*zb(t1,p))
      amp(2,1,2,2)=
     &        (2*ma*mq*(-(g1t*ma2*za(g1,p)*za(t1,p)*zb(t1,p)*
     &        (za(d,p)*zb(u,d) + za(g2,p)*zb(u,g2))) + 
     &     g1b*za(b1,p)*((mq2*za(d,p)*za(g1,p)*zb(b1,p) + 
     &           za(t1,p)*(za(d,t1)*za(g1,p)*zb(t1,b1) + 
     &              za(d,g1)*
     &               (za(g1,p)*zb(g1,b1) - za(t1,p)*zb(t1,b1)))*
     &            zb(t1,p))*zb(u,d) + 
     &        (mq2*za(g1,p)*za(g2,p)*zb(b1,p) + 
     &           za(t1,p)*(za(g1,p)*za(g2,t1)*zb(t1,b1) + 
     &              za(g1,g2)*
     &               (-(za(g1,p)*zb(g1,b1)) + za(t1,p)*zb(t1,b1)))
     &             *zb(t1,p))*zb(u,g2)))*zb(u,p))/
     & (g1b*g1t*za(b1,p)*za(t1,p)**2*zb(b1,p)*zb(d,g2)*zb(g1,b1)*
     &   zb(t1,p)*zb(u,g2))
      amp(2,2,1,1)=
     &        (2*(-(g1t*ma2*za(d,t1)*za(g1,p)*za(t1,p)*zb(t1,p)) + 
     &     g1b*za(b1,p)*(mq2*
     &         (za(d,p)*za(g1,t1) + za(d,g1)*za(t1,p))*zb(b1,p) + 
     &        za(g1,t1)*za(t1,p)*
     &         (za(d,g1)*zb(g1,b1) + za(d,t1)*zb(t1,b1))*zb(t1,p))
     &     )*(za(d,g2)*zb(g2,b1) - za(u,d)*zb(u,b1)))/
     & (g1b*g1t*za(b1,p)*za(d,g2)*za(t1,p)*za(u,g2)*zb(g1,b1)*
     &   zb(t1,p))
      amp(2,2,1,2)=
     &        (2*zb(u,b1)*(-(g1t*ma2*za(g1,p)*za(t1,p)*zb(t1,p)*
     &        (za(d,t1)*zb(u,d) + za(g2,t1)*zb(u,g2))) + 
     &     g1b*za(b1,p)*((mq2*
     &            (za(d,p)*za(g1,t1) + za(d,g1)*za(t1,p))*zb(b1,p)
     &             + za(g1,t1)*za(t1,p)*
     &            (za(d,g1)*zb(g1,b1) + za(d,t1)*zb(t1,b1))*
     &            zb(t1,p))*zb(u,d) + 
     &        (mq2*(za(g1,t1)*za(g2,p) - za(g1,g2)*za(t1,p))*
     &            zb(b1,p) + 
     &           za(g1,t1)*za(t1,p)*
     &            (-(za(g1,g2)*zb(g1,b1)) + za(g2,t1)*zb(t1,b1))*
     &            zb(t1,p))*zb(u,g2))))/
     & (g1b*g1t*za(b1,p)*za(t1,p)*zb(d,g2)*zb(g1,b1)*zb(t1,p)*
     &   zb(u,g2))
      amp(2,2,2,1)=
     &        (-2*ma*(-(g1t*ma2*za(d,t1)*za(g1,p)*za(t1,p)*zb(t1,p)) + 
     &     g1b*za(b1,p)*(mq2*
     &         (za(d,p)*za(g1,t1) + za(d,g1)*za(t1,p))*zb(b1,p) + 
     &        za(g1,t1)*za(t1,p)*
     &         (za(d,g1)*zb(g1,b1) + za(d,t1)*zb(t1,b1))*zb(t1,p))
     &     )*(za(d,g2)*zb(g2,p) - za(u,d)*zb(u,p)))/
     & (g1b*g1t*za(b1,p)*za(d,g2)*za(t1,p)*za(u,g2)*zb(b1,p)*
     &   zb(g1,b1)*zb(t1,p))
      amp(2,2,2,2)=
     &        (-2*ma*(-(g1t*ma2*za(g1,p)*za(t1,p)*zb(t1,p)*
     &        (za(d,t1)*zb(u,d) + za(g2,t1)*zb(u,g2))) + 
     &     g1b*za(b1,p)*((mq2*
     &            (za(d,p)*za(g1,t1) + za(d,g1)*za(t1,p))*zb(b1,p)
     &             + za(g1,t1)*za(t1,p)*
     &            (za(d,g1)*zb(g1,b1) + za(d,t1)*zb(t1,b1))*
     &            zb(t1,p))*zb(u,d) + 
     &        (mq2*(za(g1,t1)*za(g2,p) - za(g1,g2)*za(t1,p))*
     &            zb(b1,p) + 
     &           za(g1,t1)*za(t1,p)*
     &            (-(za(g1,g2)*zb(g1,b1)) + za(g2,t1)*zb(t1,b1))*
     &            zb(t1,p))*zb(u,g2)))*zb(u,p))/
     & (g1b*g1t*za(b1,p)*za(t1,p)*zb(b1,p)*zb(d,g2)*zb(g1,b1)*
     &   zb(t1,p)*zb(u,g2))
      return
      end
