      subroutine reals_qq(u,q1,d,t1,b1,q2,p,mq,ma,za,zb,amp)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
* AUTHORS: R. FREDERIX AND F. TRAMONTANO                               *
* DATE  : 7/17/2008                                                    *
*                                                                      *
************************************************************************
c--- 0  ->  u~ + q1~ + d + t1 + b1 + q2 (t-channel single-top)
c--- where the q1 and q2 combine to a gluon that couples to the massive quark
c--- line. Arguments of amp are helicities of that 'gluon', top and bottom,
c--- respectively. Spin of top and bottom quarks should be projected on light-
c--- like momenta p. mq is the mass of the top quark and ma the mass of the
c--- bottom quark.
 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: u,q1,q2,d,b1,t1,p
      real(dp):: bDg,tDg,bDq2,tDq2,mq,ma,mq2,ma2
      complex(dp):: amp(2,2,2)
 
 
      mq2 = mq**2
      ma2 = ma**2
      bDq2=s(q2,b1)/2._dp+s(p,q2)*ma2/2._dp/s(p,b1)
      tDq2=s(q2,t1)/2._dp+s(p,q2)*mq2/2._dp/s(p,t1)
      bDg =s(q1,q2)/2._dp+s(q1,b1)/2._dp+s(q1,p)*ma2/2._dp/s(p,b1)+bDq2
      tDg =s(q1,q2)/2._dp+s(q1,t1)/2._dp+s(q1,p)*mq2/2._dp/s(p,t1)+tDq2
 
      amp(1,1,1)=
     &        (-2._dp*mq*(ma2*tDg*za(d,p)*za(q2,p)*za(t1,p)*zb(t1,p)*
     &      (zb(q1,b1)*zb(u,p) + zb(b1,p)*zb(u,q1)) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*za(t1,p)*(za(d,t1)*za(q2,p)*zb(q1,t1) + 
     &           za(d,q2)*(za(q2,p)*zb(q1,q2) - 
     &              za(t1,p)*zb(q1,t1)))*zb(t1,p)*zb(u,b1) + 
     &        za(d,p)*(-(bDg*mq2*za(q2,p)*zb(q1,p)*zb(u,b1)) + 
     &           tDg*za(t1,p)*zb(q1,b1)*zb(t1,p)*
     &            (za(b1,q2)*zb(u,b1) + za(q1,q2)*zb(u,q1))))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)**2*zb(b1,p)*zb(q1,q2)*
     &   zb(t1,p))
      amp(1,1,2)=
     &        (2._dp*ma*mq*(ma2*tDg*za(d,p)*za(q2,p)*za(t1,p)*zb(q1,p)*
     &      zb(t1,p)*zb(u,p) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*za(t1,p)*(za(d,t1)*za(q2,p)*zb(q1,t1) + 
     &           za(d,q2)*(za(q2,p)*zb(q1,q2) - 
     &              za(t1,p)*zb(q1,t1)))*zb(t1,p)*zb(u,p) + 
     &        tDg*za(b1,q2)*za(d,p)*za(t1,p)*zb(t1,p)*
     &         (zb(q1,p)*zb(u,b1) - zb(b1,p)*zb(u,q1)) + 
     &        za(d,p)*zb(q1,p)*
     &         (-(bDg*mq2*za(q2,p)*zb(u,p)) + 
     &           tDg*za(q1,q2)*za(t1,p)*zb(t1,p)*zb(u,q1)))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)**2*zb(b1,p)**2*zb(q1,q2)*
     &   zb(t1,p))
      amp(1,2,1)=
     &        (2._dp*(ma2*tDg*za(d,t1)*za(q2,p)*za(t1,p)*zb(t1,p)*
     &      (zb(q1,b1)*zb(u,p) + zb(b1,p)*zb(u,q1)) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*mq2*za(d,p)*za(t1,q2)*zb(q1,p)*zb(u,b1) - 
     &        bDg*za(d,q2)*za(t1,p)*
     &         (mq2*zb(q1,p) + za(t1,q2)*zb(q1,q2)*zb(t1,p))*
     &         zb(u,b1) + za(d,t1)*za(t1,p)*zb(t1,p)*
     &         ((tDg*za(b1,q2)*zb(q1,b1) - 
     &              bDg*za(t1,q2)*zb(q1,t1))*zb(u,b1) + 
     &           tDg*za(q1,q2)*zb(q1,b1)*zb(u,q1)))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)*zb(b1,p)*zb(q1,q2)*zb(t1,p))
      amp(1,2,2)=
     &        (-2._dp*ma*(ma2*tDg*za(d,t1)*za(q2,p)*za(t1,p)*zb(q1,p)*
     &      zb(t1,p)*zb(u,p) - 
     &     za(b1,p)*zb(b1,p)*
     &      (-(bDg*mq2*za(d,p)*za(t1,q2)*zb(q1,p)*zb(u,p)) + 
     &        bDg*za(d,q2)*za(t1,p)*
     &         (mq2*zb(q1,p) + za(t1,q2)*zb(q1,q2)*zb(t1,p))*
     &         zb(u,p) + tDg*za(b1,q2)*za(d,t1)*za(t1,p)*zb(t1,p)*
     &         (-(zb(q1,p)*zb(u,b1)) + zb(b1,p)*zb(u,q1)) + 
     &        za(d,t1)*za(t1,p)*zb(t1,p)*
     &         (bDg*za(t1,q2)*zb(q1,t1)*zb(u,p) - 
     &           tDg*za(q1,q2)*zb(q1,p)*zb(u,q1)))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)*zb(b1,p)**2*zb(q1,q2)*
     &   zb(t1,p))
      amp(2,1,1)=
     &        (2._dp*mq*(ma2*tDg*za(d,p)*za(q1,p)*za(t1,p)*zb(t1,p)*
     &      (zb(b1,q2)*zb(u,p) - zb(b1,p)*zb(u,q2)) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*za(t1,p)*zb(t1,p)*
     &         (za(d,t1)*za(q1,p)*zb(t1,q2) + 
     &           za(d,q1)*(za(q1,p)*zb(q1,q2) - 
     &              za(t1,p)*zb(t1,q2)))*zb(u,b1) + 
     &        za(d,p)*(bDg*mq2*za(q1,p)*zb(q2,p)*zb(u,b1) - 
     &           tDg*za(t1,p)*zb(b1,q2)*zb(t1,p)*
     &            (za(q1,b1)*zb(u,b1) + za(q1,q2)*zb(u,q2))))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)**2*zb(b1,p)*zb(q1,q2)*
     &   zb(t1,p))
      amp(2,1,2)=
     &        (-2._dp*ma*mq*(-(ma2*tDg*za(d,p)*za(q1,p)*za(t1,p)*zb(q2,p)*
     &        zb(t1,p)*zb(u,p)) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*za(t1,p)*zb(t1,p)*
     &         (za(d,t1)*za(q1,p)*zb(t1,q2) + 
     &           za(d,q1)*(za(q1,p)*zb(q1,q2) - 
     &              za(t1,p)*zb(t1,q2)))*zb(u,p) + 
     &        za(d,p)*(tDg*za(q1,b1)*za(t1,p)*zb(t1,p)*
     &            (zb(q2,p)*zb(u,b1) - zb(b1,p)*zb(u,q2)) + 
     &           zb(q2,p)*(bDg*mq2*za(q1,p)*zb(u,p) + 
     &              tDg*za(q1,q2)*za(t1,p)*zb(t1,p)*zb(u,q2))))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)**2*zb(b1,p)**2*zb(q1,q2)*
     &   zb(t1,p))
      amp(2,2,1)=
     &        (-2._dp*(ma2*tDg*za(d,t1)*za(q1,p)*za(t1,p)*zb(t1,p)*
     &      (zb(b1,q2)*zb(u,p) - zb(b1,p)*zb(u,q2)) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*mq2*za(d,p)*za(q1,t1)*zb(q2,p)*zb(u,b1) + 
     &        za(t1,p)*(bDg*za(d,q1)*
     &            (mq2*zb(q2,p) + za(q1,t1)*zb(q1,q2)*zb(t1,p))*
     &            zb(u,b1) - 
     &           za(d,t1)*zb(t1,p)*
     &            ((tDg*za(q1,b1)*zb(b1,q2) - 
     &                 bDg*za(q1,t1)*zb(t1,q2))*zb(u,b1) + 
     &              tDg*za(q1,q2)*zb(b1,q2)*zb(u,q2))))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)*zb(b1,p)*zb(q1,q2)*zb(t1,p))
      amp(2,2,2)=
     &        (2._dp*ma*(-(ma2*tDg*za(d,t1)*za(q1,p)*za(t1,p)*zb(q2,p)*
     &        zb(t1,p)*zb(u,p)) + 
     &     za(b1,p)*zb(b1,p)*
     &      (bDg*(mq2*za(d,p)*za(q1,t1)*zb(q2,p) + 
     &           za(d,q1)*za(t1,p)*
     &            (mq2*zb(q2,p) + za(q1,t1)*zb(q1,q2)*zb(t1,p)))*
     &         zb(u,p) + za(d,t1)*za(t1,p)*zb(t1,p)*
     &         (bDg*za(q1,t1)*zb(t1,q2)*zb(u,p) + 
     &           tDg*za(q1,q2)*zb(q2,p)*zb(u,q2) + 
     &           tDg*za(q1,b1)*
     &            (zb(q2,p)*zb(u,b1) - zb(b1,p)*zb(u,q2))))))/
     & (za(b1,p)*za(q1,q2)*za(t1,p)*zb(b1,p)**2*zb(q1,q2)*
     &   zb(t1,p))
 
      return
      end
