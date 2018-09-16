      subroutine gampsabc(mq,p1,p2,p3,p4,t5,t6,gg_a,gg_b,gg_c)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,t5,t6
      complex(dp):: 
     & gg_a(2,2,2,2,2),gg_b(2,2,2,2,2),gg_c(2,2,2,2,2)
      real(dp):: al5,al6,s125,s126,s34,s12,s15,s25,s16,s26,mq

      al5=mq**2/s(p2,t5)
      s125=(1._dp+al5)*s(p1,p2)+s(p1,t5)+s(p2,t5)
      al6=mq**2/s(p1,t6)
      s126=(1._dp+al6)*s(p1,p2)+s(p1,t6)+s(p2,t6)
      s12=s(p1,p2)
      s34=s(p3,p4)
      s15=s(p1,t5)+al5*s(p1,p2)
      s25=s(p2,t5)
      s16=s(p1,t6)
      s26=s(p2,t6)+al6*s(p1,p2)

C   Notation is hz,h1,h2,h5,h6 
 
c      include 'aa2m.f' 

      gg_a(2,2,1,2,2)=
     & +za(p2,t6)*za(p2,t6)*za(p3,t5)*zb(p1,t6)*zb(p1,t6)*zb(p4,t6)
     & * (1._dp/s34/s126/s26/s12)
 
     & +za(p2,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t6)*zb(p1,t6)/zb(p1,p2)
     & * (-1._dp/s34/s126/s26)
 
      gg_a(2,2,1,1,1)=+ za(p2,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p4)
     & *zb(p1,t6)*zb(p2,t5)
     & * (- mq**2/s34/s126/s26/s25/s12)
 
      gg_a(2,2,1,2,1)=+ za(p2,t6)*za(p2,t6)*za(p3,t5)*zb(p1,p4)
     & *zb(p1,t6)
     & * (- mq/s34/s126/s26/s12)
 
      gg_a(2,2,1,1,2)=
     & +za(p2,t6)*za(p2,p3)*zb(p1,p4)*zb(p1,t6)*zb(p1,t6)*zb(p2,t5)
     & /zb(p1,p2)
     & * (- mq/s34/s126/s26/s25)
 
     & +za(p2,t6)*za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p1,t6)*zb(p2,t5)
     & *zb(p4,t6)
     & * (mq/s34/s126/s26/s25/s12)
 
      gg_a(2,1,2,2,2)=
     & +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,t6)*zb(p2,t6)*zb(p4,t6)
     & * (1._dp/s34/s126/s26/s12)
 
     & +za(p1,t6)*za(p3,t5)*zb(p2,p4)*zb(p2,t6)*zb(p2,t6)/zb(p1,p2)
     & * (1._dp/s34/s126/s26)
 
      gg_a(2,1,2,1,1)=
     & +za(p1,t6)*za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p2,t5)*zb(p2,t6)
     & *zb(p4,t6)
     & /za(p1,p2)
     & * (- mq**2/s34/s126/s16/s26/s25)
 
     & +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)
     & * (mq**2/s34/s126/s16/s26/s25 - mq**2/s34/s126
     & /s26/s25/s12)
 
      gg_a(2,1,2,2,1)=
     & +za(p1,t6)*za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,t6)
     & *zb(p4,t6)/za(p1,p2)
     & * (- mq/s34/s126/s16/s26)
 
     & +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,p4)*zb(p2,t6)
     & * (mq/s34/s126/s16/s26 - mq/s34/s126/s26/s12)
 
      gg_a(2,1,2,1,2)=
     & +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p2,t5)*zb(p2,t6)*zb(p2,t6)
     & *zb(p4,t6)
     & * (mq/s34/s126/s26/s25/s12)
 
     & +za(p1,t6)*za(p2,p3)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)*zb(p2,t6)
     & /zb(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
      gg_a(1,2,1,2,2)=
     & +za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)*zb(p2,p4)
     & * (mq**2/s34/s126/s16/s26/s25 - mq**2/s34/s126
     & /s26/s25/s12)
 
     & +za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)*zb(p1,t6)
     & *zb(p2,p4)
     & /zb(p1,p2)
     & * (- mq**2/s34/s126/s16/s26/s25)
 
      gg_a(1,2,1,1,1)=
     & +za(p2,t6)*za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p4,t5)/za(p1,p2)
     & * (1._dp/s34/s126/s26)
 
     & +za(p2,t6)*za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     & * (1._dp/s34/s126/s26/s12)
 
      gg_a(1,2,1,2,1)=
     & +za(p2,t6)*za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & /za(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
     & +za(p2,t6)*za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)
     & *zb(p2,p4)
     & * (mq/s34/s126/s26/s25/s12)
 
      gg_a(1,2,1,1,2)=
     & +za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     & * (mq/s34/s126/s16/s26 - mq/s34/s126/s26/s12)
 
     & +za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     & /zb(p1,p2)
     & * (- mq/s34/s126/s16/s26)
 
      gg_a(1,1,2,2,2)=+ za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p2,p4)
     & *zb(p2,t6)*zb(p2,t6)
     & * (- mq**2/s34/s126/s26/s25/s12)
 
      gg_a(1,1,2,1,1)=
     & +za(p1,t6)*za(p1,t6)*za(p1,p3)*zb(p2,t6)*zb(p4,t5)/za(p1,p2)
     & * (-1._dp/s34/s126/s26)
 
     & +za(p1,t6)*za(p1,t6)*za(p3,t6)*zb(p2,t6)*zb(p2,t6)*zb(p4,t5)
     & * (1._dp/s34/s126/s26/s12)
 
      gg_a(1,1,2,2,1)=
     & +za(p1,t6)*za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p2,p4)*zb(p2,t6)
     & /za(p1,p2)
     & * (- mq/s34/s126/s26/s25)
 
     & +za(p1,t6)*za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p2,p4)*zb(p2,t6)
     & *zb(p2,t6)
     & * (mq/s34/s126/s26/s25/s12)
 
      gg_a(1,1,2,1,2)=+ za(p1,t6)*za(p1,p3)*zb(p2,t6)*zb(p2,t6)
     & *zb(p4,t5)
     & * (- mq/s34/s126/s26/s12)
 
      gg_a(2,1,1,2,2)=
     & +za(p1,p2)*za(p1,p2)*za(p3,t5)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     & * (- mq**2/s34/s126/s16/s26)
 
     & +za(p2,t6)*za(p1,p2)*za(p3,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (-1._dp/s34/s126/s26)
 
     & +za(p2,t6)*za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p2,t6)*zb(p4,t6)
     & /zb(p1,p2)
     & /zb(p1,p2)
     & * (1._dp/s34/s126/s26)
 
      gg_a(2,1,1,1,1)=
     & +za(p2,t6)*za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s126/s26/s25)
 
     & +za(p2,t6)*za(p1,t6)*za(p2,p3)*zb(p2,t5)*zb(p4,t6)/zb(p1,p2)
     & * (- mq**2/s34/s126/s26/s25)
 
      gg_a(2,1,1,2,1)=
     & +za(p2,t6)*za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s126/s26)
 
     & +za(p2,t6)*za(p1,t6)*za(p3,t5)*zb(p4,t6)/zb(p1,p2)
     & * (- mq/s34/s126/s26)
 
      gg_a(2,1,1,1,2)=
     & +za(p1,p2)*za(p1,p2)*za(p2,p3)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     & /zb(p1,p2)
     & * (- mq**3/s34/s126/s16/s26/s25)
 
     & +za(p2,t6)*za(p1,p2)*za(p2,p3)*zb(p1,t6)*zb(p2,p4)
     & *zb(p2,t5)*zb(p2,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s126/s26/s25)
 
     & +za(p2,t6)*za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)
     & *zb(p2,t6)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
      gg_a(2,2,2,2,2)=
     & +za(p1,t6)*za(p2,t6)*za(p3,t5)*zb(p1,t6)*zb(p2,t6)
     & *zb(p4,t6)/za(p1,p2)
     & /za(p1,p2)
     & * (1._dp/s34/s126/s26)
 
     & +za(p1,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t6)*zb(p2,t6)/za(p1,p2)
     & * (1._dp/s34/s126/s26)
 
     & +za(p2,t6)*za(p3,t5)*zb(p1,p2)*zb(p2,t6)*zb(p4,t6)/za(p1,p2)
     & * (1._dp/s34/s126/s26)
 
     & +za(p3,t5)*zb(p1,p2)*zb(p1,p4)*zb(p2,t6)
     & * (mq**2/s34/s126/s16/s26 +1._dp/s34/s126/s26)
 
      gg_a(2,2,2,1,1)=
     & +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p4)
     & *zb(p1,t6)*zb(p2,t5)
     & /za(p1,p2)
     & * (mq**4/s34/s126/s16**2/s26/s25)
 
     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)
     & * (mq**2/s34/s126/s16/s26/s25 + mq**4/s34/s126
     & /s16**2/s26/s25)
 
     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)/za(p1,p2)
     & * (mq**4/s34/s126/s16/s26/s25)
 
     & +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p2)
     & *zb(p1,p2)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     & * (mq**2/s34/s126/s16/s26/s25)
 
     & +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p4)*zb(p2,t5)*zb(p2,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s126/s26/s25)
 
      gg_a(2,2,2,2,1)=
     & +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p4)*zb(p1,t6)
     & /za(p1,p2)
     & * (mq**3/s34/s126/s16**2/s26)
 
     & +za(p1,t6)*za(p2,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p2)*zb(p4,t6)
     & /za(p1,p2)
     & * (mq/s34/s126/s16/s26)
 
     & +za(p1,t6)*za(p2,t6)*za(p3,t5)*zb(p1,p4)*zb(p2,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s126/s26)
 
     & +za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p2)*zb(p1,p4)
     & * (mq/s34/s126/s16/s26 + mq**3/s34/s126/s16**2
     & /s26)
 
     & +za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p1,p4)/za(p1,p2)
     & * (mq**3/s34/s126/s16/s26)
 
      gg_a(2,2,2,1,2)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)*zb(p2,t6)
     & /za(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
     & +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)*zb(p2,t6)
     & *zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
     & +za(p2,p3)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)*zb(p2,t6)
     & * (mq/s34/s126/s26/s25 + mq**3/s34/s126/s16
     & /s26/s25)
 
     & +za(p2,t6)*za(p2,p3)*zb(p1,p2)*zb(p2,t5)*zb(p2,t6)*zb(p4,t6)
     & /za(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
      gg_a(1,1,1,2,2)=
     & +za(p1,p2)*za(p1,p2)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & * (mq**2/s34/s126/s16/s26/s25 + mq**4/s34/s126
     & /s16**2/s26/s25)
 
     & +za(p1,p2)*za(p1,p2)*za(p3,t6)*za(p2,t5)*zb(p1,t6)
     & *zb(p2,p4)*zb(p2,t6)
     & /zb(p1,p2)
     & * (mq**2/s34/s126/s16/s26/s25)
 
     & +za(p1,p2)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     & * (mq**4/s34/s126/s16/s26/s25)
 
     & +za(p1,p2)*za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)
     & *zb(p2,p4)
     & /zb(p1,p2)
     & * (mq**4/s34/s126/s16**2/s26/s25)
 
     & +za(p2,t6)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)
     & /zb(p1,p2)
     & /zb(p1,p2)
     & * (- mq**2/s34/s126/s26/s25)
 
      gg_a(1,1,1,1,1)=
     & +za(p2,t6)*za(p1,p2)*za(p1,p3)*zb(p4,t5)
     & * (mq**2/s34/s126/s16/s26 +1._dp/s34/s126/s26)
 
     & +za(p2,t6)*za(p1,p2)*za(p3,t6)*zb(p2,t6)*zb(p4,t5)/zb(p1,p2)
     & * (1._dp/s34/s126/s26)
 
     & +za(p2,t6)*za(p1,t6)*za(p1,p3)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     & * (1._dp/s34/s126/s26)
 
     & +za(p2,t6)*za(p1,t6)*za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     & /zb(p1,p2)
     & /zb(p1,p2)
     & * (1._dp/s34/s126/s26)
 
      gg_a(1,1,1,2,1)=
     & +za(p2,t6)*za(p1,p2)*za(p1,p3)*za(p2,t5)*zb(p2,p4)
     & * (mq/s34/s126/s26/s25 + mq**3/s34/s126/s16
     & /s26/s25)
 
     & +za(p2,t6)*za(p1,p2)*za(p3,t6)*za(p2,t5)*zb(p2,p4)*zb(p2,t6)
     & /zb(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
     & +za(p2,t6)*za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
     & +za(p2,t6)*za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & *zb(p2,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
      gg_a(1,1,1,1,2)=
     & +za(p1,p2)*za(p1,p2)*za(p1,p3)*zb(p1,t6)*zb(p4,t5)
     & * (mq/s34/s126/s16/s26 + mq**3/s34/s126/s16**2
     & /s26)
 
     & +za(p1,p2)*za(p1,p2)*za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     &  /zb(p1,p2)
     & * (mq/s34/s126/s16/s26)
 
     & +za(p1,p2)*za(p1,p3)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     & * (mq**3/s34/s126/s16/s26)
 
     & +za(p1,p2)*za(p1,t6)*za(p1,p3)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     & /zb(p1,p2)
     & * (mq**3/s34/s126/s16**2/s26)
 
     & +za(p2,t6)*za(p1,p3)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s126/s26)
 
      gg_a(1,2,2,2,2)=
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s126/s26/s25)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t6)/za(p1,p2)
     & * (- mq**2/s34/s126/s26/s25)
 
      gg_a(1,2,2,1,1)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p2)*zb(p4,t5)/za(p1,p2)
     & * (- mq**2/s34/s126/s16/s26)
 
     & +za(p1,t6)*za(p2,t6)*za(p2,p3)*zb(p1,p2)*zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (-1._dp/s34/s126/s26)
 
     & +za(p1,t6)*za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (1._dp/s34/s126/s26)
 
      gg_a(1,2,2,2,1)=
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p2,p4)
     & /za(p1,p2)
     & * (- mq**3/s34/s126/s16/s26/s25)
 
     & +za(p1,t6)*za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)
     & *zb(p2,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s126/s26/s25)
 
     & +za(p1,t6)*za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & *zb(p2,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s126/s26/s25)
 
      gg_a(1,2,2,1,2)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s126/s26)
 
     & +za(p3,t6)*zb(p1,t6)*zb(p2,t6)*zb(p4,t5)/za(p1,p2)
     & * (- mq/s34/s126/s26)
 

c       include 'bb2m.f' 
      gg_b(2,2,1,2,2)=+ za(p2,t6)*za(p3,t5)*za(p2,t5)*zb(p1,t5)
     & *zb(p1,t6)*zb(p4,t6)
     & * (-1._dp/s34/s26/s15/s12)
 
      gg_b(2,2,1,1,1)=+ za(p2,t6)*za(p2,p3)*zb(p1,p4)*zb(p1,t5)
     & * (- mq**2/s34/s26/s15/s12)
 
      gg_b(2,2,1,2,1)=+ za(p2,t6)*za(p3,t5)*za(p2,t5)*zb(p1,p4)
     & *zb(p1,t5)
     & * (mq/s34/s26/s15/s12)
 
      gg_b(2,2,1,1,2)=+ za(p2,t6)*za(p2,p3)*zb(p1,t5)*zb(p1,t6)
     & *zb(p4,t6)
     & * (mq/s34/s26/s15/s12)
 
      gg_b(2,1,2,2,2)=
     & +za(p1,p3)*za(p1,t5)*zb(p2,p4)*zb(p2,t6)
     & * (-1._dp/s34/s26/s15)
 
     & +za(p1,t6)*za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p2,t6)*zb(p4,t6)
     & * (-1._dp/s34/s26/s15/s12)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p2,t5)*zb(p2,t6)*zb(p4,t6)
     & * (-1._dp/s34/s26/s15/s12)
 
     & +za(p3,t5)*za(p1,t5)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)/zb(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
      gg_b(2,1,2,1,1)=
     & +za(p1,t6)*za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p2,p4)*zb(p2,t5)
     & * (mq**2/s34/s16/s26/s15/s25)
 
     & +za(p1,t6)*za(p1,p3)*zb(p2,p4)*zb(p2,t5)
     & * (mq**2/s34/s16/s26/s15 + mq**2/s34/s26/s15
     & /s25 - mq**2/s34/s26/s15/s12)
 
     & +za(p1,t6)*za(p1,t6)*za(p1,p3)*zb(p1,p2)*zb(p2,t5)*zb(p4,t6)
     & * (- mq**2/s34/s16/s26/s15/s25)
 
     & +za(p1,t6)*za(p1,t6)*za(p1,p3)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     & * (- mq**2/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p1,t6)*za(p3,t5)*zb(p2,t5)*zb(p2,t5)*zb(p4,t6)
     & * (- mq**2/s34/s16/s26/s15/s25)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,p2)*zb(p2,p4)*zb(p2,t5)*zb(p2,t5)
     & * (mq**2/s34/s16/s26/s15/s25 - mq**2/s34/s26
     & /s15/s25/s12)
 
      gg_b(2,1,2,2,1)=
     & +za(p1,t6)*za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p2,p4)
     & * (- mq/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p1,p3)*za(p1,t5)*zb(p2,p4)/za(p1,p2)
     & * (- mq/s34/s26/s15)
 
     & +za(p1,t6)*za(p1,t6)*za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p4,t6)
     & /za(p1,p2)
     & * (mq/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /za(p1,p2)
     & * (mq/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p2,p4)*zb(p2,t5)
     & * (- mq/s34/s16/s26/s15 + mq/s34/s26/s15/s12)
 
      gg_b(2,1,2,1,2)=
     & +za(p1,p3)*za(p1,p2)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)
     & * (mq/s34/s26/s15/s25)
 
     & +za(p1,p3)*zb(p2,p4)*zb(p2,t5)*zb(p2,t6)/zb(p1,p2)
     & * (mq/s34/s26/s15)
 
     & +za(p1,t6)*za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p2,t5)*zb(p2,t6)
     & *zb(p4,t6)
     & * (mq/s34/s26/s15/s25/s12)
 
     & +za(p1,t6)*za(p1,p3)*zb(p2,t5)*zb(p2,t6)*zb(p4,t6)
     & * (mq/s34/s26/s15/s12)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,p2)*zb(p2,t5)*zb(p2,t5)*zb(p2,t6)
     & *zb(p4,t6)
     & * (mq/s34/s26/s15/s25/s12)
 
     & +za(p3,t5)*za(p1,p2)*zb(p2,p4)*zb(p2,t5)*zb(p2,t5)*zb(p2,t6)
     & /zb(p1,p2)
     & * (mq/s34/s26/s15/s25)
 
      gg_b(1,2,1,2,2)=
     & +za(p2,p3)*za(p1,p2)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)*zb(p1,t6)
     & * (mq**2/s34/s16/s26/s15/s25)
 
     & +za(p2,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)*zb(p4,t5)
     & * (mq**2/s34/s16/s26/s15/s25 - mq**2/s34/s26
     & /s15/s25/s12)
 
     & +za(p2,p3)*za(p2,t5)*zb(p1,p4)*zb(p1,t6)
     & * (mq**2/s34/s16/s26/s15 + mq**2/s34/s26/s15
     & /s25 - mq**2/s34/s26/s15/s12)
 
     & +za(p3,t6)*za(p1,p2)*za(p2,t5)*zb(p1,p4)*zb(p1,t6)*zb(p1,t6)
     & * (- mq**2/s34/s16/s26/s15/s25)
 
     & +za(p3,t6)*za(p2,t5)*za(p2,t5)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     & * (- mq**2/s34/s16/s26/s15/s25)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p4)*zb(p1,t6)*zb(p1,t6)/zb(p1,p2)
     & * (- mq**2/s34/s16/s26/s15)
 
      gg_b(1,2,1,1,1)=
     & +za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,t5)*zb(p4,t5)/za(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
     & +za(p2,t6)*za(p2,p3)*zb(p1,p4)*zb(p1,t5)
     & * (-1._dp/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t6)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)*zb(p1,t6)
     & * (-1._dp/s34/s26/s15/s12)
 
     & +za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     & * (-1._dp/s34/s26/s15/s12)
 
      gg_b(1,2,1,2,1)=
     & +za(p2,t6)*za(p2,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p4,t5)
     & /za(p1,p2)
     & * (mq/s34/s26/s15/s25)
 
     & +za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     & * (mq/s34/s26/s15/s25)
 
     & +za(p2,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p4)/za(p1,p2)
     & * (mq/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t6)*za(p1,p2)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     & *zb(p1,t6)
     & * (mq/s34/s26/s15/s25/s12)
 
     & +za(p2,t6)*za(p3,t6)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)
     & *zb(p4,t5)
     & * (mq/s34/s26/s15/s25/s12)
 
     & +za(p2,t6)*za(p3,t6)*za(p2,t5)*zb(p1,p4)*zb(p1,t6)
     & * (mq/s34/s26/s15/s12)
 
      gg_b(1,2,1,1,2)=
     & +za(p2,p3)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)*zb(p1,t6)
     & * (- mq/s34/s16/s26/s15)
 
     & +za(p2,p3)*za(p2,t5)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     & * (- mq/s34/s16/s26/s15 + mq/s34/s26/s15/s12)
 
     & +za(p2,p3)*zb(p1,p4)*zb(p1,t5)*zb(p1,t6)/zb(p1,p2)
     & * (- mq/s34/s26/s15)
 
     & +za(p3,t6)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)*zb(p1,t6)*zb(p1,t6)
     & /zb(p1,p2)
     & * (mq/s34/s16/s26/s15)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,t5)*zb(p1,t6)*zb(p1,t6)*zb(p4,t5)
     & /zb(p1,p2)
     & * (mq/s34/s16/s26/s15)
 
      gg_b(1,1,2,2,2)=+ za(p1,p3)*za(p1,t5)*zb(p2,p4)*zb(p2,t6)
     & * (- mq**2/s34/s26/s15/s12)
 
      gg_b(1,1,2,1,1)=+ za(p1,t6)*za(p3,t6)*za(p1,t5)*zb(p2,t5)
     & *zb(p2,t6)*zb(p4,t5)
     & * (-1._dp/s34/s26/s15/s12)
 
      gg_b(1,1,2,2,1)=+ za(p1,t6)*za(p3,t6)*za(p1,t5)*zb(p2,p4)
     & *zb(p2,t6)
     & * (mq/s34/s26/s15/s12)
 
      gg_b(1,1,2,1,2)=+ za(p1,p3)*za(p1,t5)*zb(p2,t5)*zb(p2,t6)
     & *zb(p4,t5)
     & * (mq/s34/s26/s15/s12)
 
      gg_b(2,1,1,2,2)=
     & +za(p2,t6)*za(p1,p3)*za(p1,t5)*zb(p1,t6)*zb(p4,t6)/zb(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t5)*za(p1,t5)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
      gg_b(2,1,1,1,1)=
     & +za(p2,t6)*za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)
     & * (- mq**2/s34/s26/s15/s25)
 
     & +za(p2,t6)*za(p1,p3)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t5)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s26/s15/s25)
 
      gg_b(2,1,1,2,1)=
     & +za(p2,t6)*za(p1,p3)*za(p1,t5)*zb(p1,p4)/zb(p1,p2)
     & * (mq/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t5)*za(p1,t5)*zb(p1,p4)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s26/s15)
 
      gg_b(2,1,1,1,2)=
     & +za(p2,t6)*za(p1,p3)*za(p1,p2)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)
     & * (mq/s34/s26/s15/s25)
 
     & +za(p2,t6)*za(p1,p3)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t5)*za(p1,p2)*zb(p1,t6)*zb(p2,t5)*zb(p2,t5)
     & *zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s26/s15/s25)
 
      gg_b(2,2,2,2,2)=
     & +za(p1,t6)*za(p3,t5)*za(p2,t5)*zb(p1,t5)*zb(p2,t6)*zb(p4,t6)
     & /za(p1,p2)
     & /za(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
     & +za(p3,t5)*za(p2,t5)*zb(p1,t5)*zb(p2,p4)*zb(p2,t6)/za(p1,p2)
     & * (1._dp/s34/s26/s15)
 
      gg_b(2,2,2,1,1)=
     & +za(p1,t6)*za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,t5)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (mq**2/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,t5)*zb(p2,p4)/za(p1,p2)
     & * (- mq**2/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p2,p3)*zb(p1,t5)*zb(p2,p4)/za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s26/s15)
 
      gg_b(2,2,2,2,1)=
     & +za(p1,t6)*za(p1,t6)*za(p3,t5)*za(p2,t5)*zb(p1,p2)
     & *zb(p1,t5)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t5)
     & *zb(p2,p4)/za(p1,p2)
     & * (mq/s34/s16/s26/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p2,t5)*zb(p1,t5)*zb(p2,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s26/s15)
 
      gg_b(2,2,2,1,2)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,t5)*zb(p2,t6)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s26/s15)
 
     & +za(p2,p3)*zb(p1,t5)*zb(p2,p4)*zb(p2,t6)/za(p1,p2)
     & * (- mq/s34/s26/s15)
 
      gg_b(1,1,1,2,2)=
     & +za(p1,p2)*za(p2,p3)*za(p1,t5)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     & * (- mq**2/s34/s16/s26/s15)
 
     & +za(p1,p2)*za(p3,t6)*za(p1,t5)*zb(p1,t6)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)
     & /zb(p1,p2)
     & * (mq**2/s34/s16/s26/s15)
 
     & +za(p2,p3)*za(p1,t5)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s26/s15)
 
      gg_b(1,1,1,1,1)=
     & +za(p2,t6)*za(p2,p3)*za(p1,t5)*zb(p2,t5)*zb(p4,t5)/zb(p1,p2)
     & * (1._dp/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t6)*za(p1,t5)*zb(p1,t6)*zb(p2,t5)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
      gg_b(1,1,1,2,1)=
     & +za(p2,t6)*za(p2,p3)*za(p1,t5)*zb(p2,p4)/zb(p1,p2)
     & * (- mq/s34/s26/s15)
 
     & +za(p2,t6)*za(p3,t6)*za(p1,t5)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s26/s15)
 
      gg_b(1,1,1,1,2)=
     & +za(p1,p2)*za(p2,p3)*za(p1,t5)*zb(p1,t6)*zb(p2,t5)*zb(p4,t5)
     & /zb(p1,p2)
     & * (mq/s34/s16/s26/s15)
 
     & +za(p1,p2)*za(p3,t6)*za(p1,t5)*zb(p1,t6)*zb(p1,t6)*zb(p2,t5)
     & *zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s16/s26/s15)
 
     & +za(p2,p3)*za(p1,t5)*zb(p1,t6)*zb(p2,t5)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s26/s15)
 
      gg_b(1,2,2,2,2)=
     & +za(p1,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s26/s15/s25)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)*zb(p2,t6)/za(p1,p2)
     & * (- mq**2/s34/s26/s15/s25)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p4)*zb(p2,t6)/za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s26/s15)
 
      gg_b(1,2,2,1,1)=
     & +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t5)*zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
     & +za(p1,t6)*za(p3,t6)*zb(p1,p4)*zb(p1,t5)*zb(p2,t6)/za(p1,p2)
     & * (-1._dp/s34/s26/s15)
 
      gg_b(1,2,2,2,1)=
     & +za(p1,t6)*za(p3,t6)*za(p2,t5)*za(p2,t5)*zb(p1,p2)
     & *zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s26/s15/s25)
 
     & +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)*zb(p2,t6)
     & /za(p1,p2)
     & * (mq/s34/s26/s15/s25)
 
     & +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,p4)*zb(p2,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s26/s15)
 
      gg_b(1,2,2,1,2)=
     & +za(p1,p3)*za(p2,t5)*zb(p1,t5)*zb(p2,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s26/s15)
 
     & +za(p1,p3)*zb(p1,p4)*zb(p1,t5)*zb(p2,t6)/za(p1,p2)
     & * (mq/s34/s26/s15)
 
c       include 'cc2m.f' 

      gg_c(2,2,1,2,2)=
     & +za(p2,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t5)*zb(p4,t6)
     & * (mq**2/s34/s125/s15/s25/s12)
 
     & +za(p2,p3)*za(p2,t5)*za(p2,t5)*zb(p1,t5)*zb(p4,t6)/za(p1,p2)
     & * (mq**2/s34/s125/s15/s25 +1._dp/s34/s125/s15)
 
     & +za(p3,t5)*za(p2,t5)*za(p2,t5)*zb(p1,t5)*zb(p1,t5)*zb(p4,t6)
     & * (1._dp/s34/s125/s15/s12)
 
      gg_c(2,2,1,1,1)=+ za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p4)
     & *zb(p1,t5)*zb(p1,t5)
     & * (- mq**2/s34/s125/s16/s15/s12)
 
      gg_c(2,2,1,2,1)=
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)
     & *zb(p1,p4)*zb(p1,t5)
     & * (mq**3/s34/s125/s16/s15/s25/s12)
 
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p4)
     & *zb(p1,t5)/za(p1,p2)
     & * (mq/s34/s125/s16/s15 + mq**3/s34/s125/s16
     & /s15/s25)
 
     & +za(p1,t6)*za(p3,t5)*za(p2,t5)*za(p2,t5)*zb(p1,p4)
     & *zb(p1,t5)*zb(p1,t5)
     & * (mq/s34/s125/s16/s15/s12)
 
      gg_c(2,2,1,1,2)=+za(p2,p3)*za(p2,t5)*zb(p1,t5)*zb(p1,t5)*zb(p4,t6)
     & * (- mq/s34/s125/s15/s12)
 
      gg_c(2,1,2,2,2)=
     & +za(p1,p3)*za(p1,t5)*za(p1,t5)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     & * (-1._dp/s34/s125/s15)
 
     & +za(p3,t5)*za(p1,t5)*za(p1,t5)*zb(p2,t5)*zb(p2,t5)*zb(p4,t6)
     & * (1._dp/s34/s125/s15/s12)
 
      gg_c(2,1,2,1,1)=
     & +za(p1,t6)*za(p1,p3)*za(p1,t5)*zb(p1,p4)*zb(p2,t5)*zb(p2,t5)
     & * (mq**2/s34/s125/s16/s15/s25 - mq**2/s34/s125
     & /s16/s15/s12)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*za(p1,p2)*zb(p1,p4)
     & *zb(p2,t5)*zb(p2,t5)
     & *zb(p2,t5)
     & * (- mq**2/s34/s125/s16/s15/s25/s12)
 
      gg_c(2,1,2,2,1)=
     & +za(p1,t6)*za(p1,p3)*za(p1,t5)*za(p1,t5)*zb(p1,p4)
     & *zb(p2,t5)/za(p1,p2)
     & * (- mq/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*za(p1,t5)*zb(p1,p4)
     & *zb(p2,t5)*zb(p2,t5)
     & * (mq/s34/s125/s16/s15/s12)
 
      gg_c(2,1,2,1,2)=
     & +za(p1,p3)*za(p1,t5)*zb(p2,t5)*zb(p2,t5)*zb(p4,t6)
     & * (mq/s34/s125/s15/s25 - mq/s34/s125/s15/s12)
 
     & +za(p3,t5)*za(p1,t5)*za(p1,p2)*zb(p2,t5)*zb(p2,t5)
     & *zb(p2,t5)*zb(p4,t6)
     & * (- mq/s34/s125/s15/s25/s12)
 
      gg_c(1,2,1,2,2)=
     & +za(p1,p3)*za(p2,t5)*za(p2,t5)*za(p2,t5)
     & *zb(p1,p2)*zb(p1,t5)*zb(p1,t6)
     & *zb(p4,t5)
     & * (- mq**2/s34/s125/s16/s15/s25/s12)
 
     & +za(p1,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)*zb(p1,t6)
     & * (mq**2/s34/s125/s16/s15/s25 - mq**2/s34/s125
     & /s16/s15/s12)
 
      gg_c(1,2,1,1,1)=
     & +za(p3,t6)*za(p2,t5)*za(p2,t5)*zb(p1,t5)*zb(p1,t5)*zb(p4,t5)
     & * (1._dp/s34/s125/s15/s12)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)*zb(p1,t5)/zb(p1,p2)
     & * (-1._dp/s34/s125/s15)
 
      gg_c(1,2,1,2,1)=
     & +za(p3,t6)*za(p2,t5)*za(p2,t5)*za(p2,t5)
     & *zb(p1,p2)*zb(p1,t5)*zb(p4,t5)
     & * (- mq/s34/s125/s15/s25/s12)
 
     & +za(p3,t6)*za(p2,t5)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)
     & * (mq/s34/s125/s15/s25 - mq/s34/s125/s15/s12)
 
      gg_c(1,2,1,1,2)=
     & +za(p1,p3)*za(p2,t5)*za(p2,t5)*zb(p1,t5)
     & *zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     & * (mq/s34/s125/s16/s15/s12)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)
     & *zb(p1,t5)*zb(p1,t6)/zb(p1,p2)
     & * (- mq/s34/s125/s16/s15)
 
      gg_c(1,1,2,2,2)=+ za(p1,p3)*za(p1,t5)*za(p1,t5)*zb(p1,t6)
     & *zb(p2,p4)*zb(p2,t5)
     & * (- mq**2/s34/s125/s16/s15/s12)
 
      gg_c(1,1,2,1,1)=
     & +za(p3,t6)*za(p1,t5)*za(p1,p2)*zb(p2,p4)*zb(p2,t5)*zb(p2,t5)
     & * (mq**2/s34/s125/s15/s25/s12)
 
     & +za(p3,t6)*za(p1,t5)*za(p1,t5)*zb(p2,t5)*zb(p2,t5)*zb(p4,t5)
     & * (1._dp/s34/s125/s15/s12)
 
     & +za(p3,t6)*za(p1,t5)*zb(p2,p4)*zb(p2,t5)*zb(p2,t5)/zb(p1,p2)
     & * (mq**2/s34/s125/s15/s25 +1._dp/s34/s125/s15)
 
      gg_c(1,1,2,2,1)=+za(p3,t6)*za(p1,t5)*za(p1,t5)*zb(p2,p4)*zb(p2,t5)
     & * (- mq/s34/s125/s15/s12)
 
      gg_c(1,1,2,1,2)=
     & +za(p1,p3)*za(p1,t5)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)
     & *zb(p2,t5)*zb(p2,t5)
     & * (mq**3/s34/s125/s16/s15/s25/s12)
 
     & +za(p1,p3)*za(p1,t5)*za(p1,t5)*zb(p1,t6)*zb(p2,t5)
     & *zb(p2,t5)*zb(p4,t5)
     & * (mq/s34/s125/s16/s15/s12)
 
     & +za(p1,p3)*za(p1,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     & *zb(p2,t5)/zb(p1,p2)
     & * (mq/s34/s125/s16/s15 + mq**3/s34/s125/s16
     & /s15/s25)
 
      gg_c(2,1,1,2,2)=
     & +za(p2,p3)*za(p1,p2)*za(p1,t5)*zb(p4,t6)
     & * (- mq**2/s34/s125/s15/s25 -1._dp/s34/s125/s15)
 
     & +za(p2,p3)*za(p1,t5)*zb(p4,t6)/zb(p1,p2)
     & * (- mq**2/s34/s125/s15)
 
     & +za(p2,p3)*za(p2,t5)*za(p1,t5)*zb(p2,t5)*zb(p4,t6)/zb(p1,p2)
     & * (- mq**2/s34/s125/s15/s25 -1._dp/s34/s125/s15)
 
     & +za(p3,t5)*za(p1,p2)*za(p1,t5)*zb(p1,t5)*zb(p4,t6)/zb(p1,p2)
     & * (1._dp/s34/s125/s15)
 
     & +za(p3,t5)*za(p2,t5)*za(p1,t5)*zb(p1,t5)*zb(p2,t5)
     & *zb(p4,t6)/zb(p1,p2)
     & /zb(p1,p2)
     & * (1._dp/s34/s125/s15)
 
      gg_c(2,1,1,1,1)=
     & +za(p1,t6)*za(p2,p3)*za(p1,p2)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)
     & * (mq**2/s34/s125/s16/s15/s25 + mq**4/s34/s125
     & /s16/s15/s25**2)
 
     & +za(p1,t6)*za(p2,p3)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)
     & * (mq**2/s34/s125/s16/s15 + mq**4/s34/s125/s16
     & /s15/s25)
 
     & +za(p1,t6)*za(p2,p3)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)
     & *zb(p2,t5)/zb(p1,p2)
     & /zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*za(p1,p2)*zb(p1,p4)
     & *zb(p2,t5)*zb(p2,t5)
     & /zb(p1,p2)
     & * (mq**2/s34/s125/s16/s15/s25 + mq**4/s34/s125
     & /s16/s15/s25**2)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,p2)*za(p1,p2)*zb(p1,p4)
     & *zb(p1,t5)*zb(p2,t5)
     & /zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s15/s25)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)
     & *zb(p2,t5)/zb(p1,p2)
     & /zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p2,t5)*za(p1,p2)*zb(p1,p4)
     & *zb(p1,t5)*zb(p2,t5)
     & *zb(p2,t5)/zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s15/s25)
 
      gg_c(2,1,1,2,1)=
     & +za(p1,t6)*za(p2,p3)*za(p1,p2)*za(p1,t5)*zb(p1,p4)
     & * (- mq/s34/s125/s16/s15 - mq**3/s34/s125/s16
     & /s15/s25)
 
     & +za(p1,t6)*za(p2,p3)*za(p1,t5)*zb(p1,p4)/zb(p1,p2)
     & * (- mq**3/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*za(p1,t5)*zb(p1,p4)
     & *zb(p2,t5)/zb(p1,p2)
     & * (- mq/s34/s125/s16/s15 - mq**3/s34/s125/s16
     & /s15/s25)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,p2)*za(p1,t5)*zb(p1,p4)
     & *zb(p1,t5)/zb(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p2,t5)*za(p1,t5)*zb(p1,p4)
     & *zb(p1,t5)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
      gg_c(2,1,1,1,2)=
     & +za(p2,p3)*za(p1,p2)*za(p1,p2)*zb(p2,t5)*zb(p4,t6)
     & * (mq/s34/s125/s15/s25 + mq**3/s34/s125/s15
     & /s25**2)
 
     & +za(p2,p3)*za(p1,p2)*zb(p2,t5)*zb(p4,t6)/zb(p1,p2)
     & * (mq/s34/s125/s15 + mq**3/s34/s125/s15/s25)
 
     & +za(p2,p3)*za(p1,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125/s15)
 
     & +za(p2,p3)*za(p2,t5)*za(p1,p2)*zb(p2,t5)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)
     & * (mq/s34/s125/s15/s25 + mq**3/s34/s125/s15
     & /s25**2)
 
     & +za(p3,t5)*za(p1,p2)*za(p1,p2)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)
     & * (- mq/s34/s125/s15/s25)
 
     & +za(p3,t5)*za(p1,p2)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125/s15)
 
     & +za(p3,t5)*za(p2,t5)*za(p1,p2)*zb(p1,t5)*zb(p2,t5)
     & *zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125/s15/s25)
 
      gg_c(2,2,2,2,2)=
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t5)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (1._dp/s34/s125/s15)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p4,t6)/za(p1,p2)
     & * (mq**2/s34/s125/s15/s25)
 
     & +za(p3,t5)*za(p1,t5)*za(p2,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (1._dp/s34/s125/s15)
 
      gg_c(2,2,2,1,1)=
     & +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)/za(p1,p2)
     & * (mq**2/s34/s125/s16/s15)
 
      gg_c(2,2,2,2,1)=
     & +za(p1,t6)*za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     & *zb(p1,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p1,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     & *zb(p2,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq**3/s34/s125/s16/s15/s25)
 
     & +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p1,p4)
     & /za(p1,p2)
     & * (mq**3/s34/s125/s16/s15/s25)
 
     & +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (mq**3/s34/s125/s16/s15)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)
     & *zb(p2,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
      gg_c(2,2,2,1,2)=
     & +za(p1,p3)*za(p2,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125/s15)
 
     & +za(p3,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     & * (mq/s34/s125/s15)
 
      gg_c(1,1,1,2,2)=
     & +za(p1,p3)*za(p2,t5)*za(p1,t5)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p2,t5)*za(p1,t5)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     & * (mq**2/s34/s125/s16/s15)
 
      gg_c(1,1,1,1,1)=
     & +za(p3,t6)*za(p1,p2)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)
     & * (mq**2/s34/s125/s15/s25)
 
     & +za(p3,t6)*za(p1,p2)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (1._dp/s34/s125/s15)
 
     & +za(p3,t6)*za(p2,t5)*za(p1,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (1._dp/s34/s125/s15)
 
      gg_c(1,1,1,2,1)=
     & +za(p3,t6)*za(p2,t5)*za(p1,t5)*zb(p1,p4)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125/s15)
 
     & +za(p3,t6)*za(p2,t5)*za(p1,t5)*zb(p4,t5)/zb(p1,p2)
     & * (mq/s34/s125/s15)
 
      gg_c(1,1,1,1,2)=
     & +za(p1,p3)*za(p1,p2)*za(p1,p2)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)
     & /zb(p1,p2)
     & * (mq**3/s34/s125/s16/s15/s25)
 
     & +za(p1,p3)*za(p1,p2)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)*zb(p1,t6)
     & *zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq**3/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p2,t5)*za(p1,p2)*zb(p1,p4)*zb(p1,t6)
     & *zb(p2,t5)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq**3/s34/s125/s16/s15/s25)
 
     & +za(p1,p3)*za(p2,t5)*za(p1,t5)*zb(p1,t5)*zb(p1,t6)
     & *zb(p2,t5)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
      gg_c(1,2,2,2,2)=
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*za(p2,t5)*zb(p1,p2)
     & *zb(p1,t6)*zb(p2,t5)
     & *zb(p4,t5)/za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s125/s16/s15/s25)
 
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)
     & *zb(p1,t6)*zb(p4,t5)
     & /za(p1,p2)
     & * (- mq**2/s34/s125/s16/s15/s25)
 
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)
     & *zb(p4,t5)/za(p1,p2)
     & /za(p1,p2)
     & * (- mq**2/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,t5)*zb(p1,t6)
     & *zb(p2,p4)/za(p1,p2)
     & /za(p1,p2)
     & * (- mq**2/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)
     & *zb(p2,p4)*zb(p2,t5)
     & /za(p1,p2)
     & * (mq**2/s34/s125/s16/s15/s25 + mq**4/s34/s125
     & /s16/s15/s25**2)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p1,t6)*zb(p2,p4)
     & * (mq**2/s34/s125/s16/s15/s25 + mq**4/s34/s125
     & /s16/s15/s25**2)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     & * (mq**2/s34/s125/s16/s15 + mq**4/s34/s125/s16
     & /s15/s25)
 
      gg_c(1,2,2,1,1)=
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p1,t5)*zb(p2,t5)
     & *zb(p4,t5)/za(p1,p2)
     & /za(p1,p2)
     & * (1._dp/s34/s125/s15)
 
     & +za(p3,t6)*za(p1,t5)*zb(p1,p2)*zb(p1,t5)*zb(p4,t5)/za(p1,p2)
     & * (1._dp/s34/s125/s15)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,t5)*zb(p2,p4)*zb(p2,t5)/za(p1,p2)
     & * (- mq**2/s34/s125/s15/s25 -1._dp/s34/s125/s15)
 
     & +za(p3,t6)*zb(p1,p2)*zb(p1,t5)*zb(p2,p4)
     & * (- mq**2/s34/s125/s15/s25 -1._dp/s34/s125/s15)
 
     & +za(p3,t6)*zb(p1,t5)*zb(p2,p4)/za(p1,p2)
     & * (- mq**2/s34/s125/s15)
 
      gg_c(1,2,2,2,1)=
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p2,t5)
     & *zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125/s15/s25)
 
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p4,t5)
     & /za(p1,p2)
     & * (- mq/s34/s125/s15/s25)
 
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125/s15)
 
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p1,t5)*zb(p2,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125/s15)
 
     & +za(p3,t6)*za(p2,t5)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)*zb(p2,t5)
     & /za(p1,p2)
     & * (mq/s34/s125/s15/s25 + mq**3/s34/s125/s15
     & /s25**2)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p2)*zb(p1,p2)*zb(p2,p4)
     & * (mq/s34/s125/s15/s25 + mq**3/s34/s125/s15
     & /s25**2)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)/za(p1,p2)
     & * (mq/s34/s125/s15 + mq**3/s34/s125/s15/s25)
 
      gg_c(1,2,2,1,2)=
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,t5)*zb(p1,t6)*zb(p2,t5)
     & *zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     & /za(p1,p2)
     & * (mq/s34/s125/s16/s15)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     & /za(p1,p2)
     & * (- mq/s34/s125/s16/s15 - mq**3/s34/s125/s16
     & /s15/s25)
 
     & +za(p1,p3)*zb(p1,p2)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)
     & * (- mq/s34/s125/s16/s15 - mq**3/s34/s125/s16
     & /s15/s25)
 
     & +za(p1,p3)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     & * (- mq**3/s34/s125/s16/s15)
 


 
      return
      end
