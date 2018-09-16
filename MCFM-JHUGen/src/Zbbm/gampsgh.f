      subroutine gampsgh(mq,p1,p2,p3,p4,t5,t6,gg_g,gg_h)
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
     & gg_g(2,2,2,2,2),gg_h(2,2,2,2,2)
      real(dp):: al5,al6,s125,s126,s34,s25,s16,mq

      al5=mq**2/s(p2,t5)
      s125=(1._dp+al5)*s(p1,p2)+s(p1,t5)+s(p2,t5)
      al6=mq**2/s(p1,t6)
      s126=(1._dp+al6)*s(p1,p2)+s(p1,t6)+s(p2,t6)
      s34=s(p3,p4)
      s25=s(p2,t5)
      s16=s(p1,t6)

C   Notation is hz,h1,h2,h5,h6 
 
c      include 'gg2m.f'
 
      gg_g(2,2,1,2,2)=czip
 
      gg_g(2,2,1,1,1)=czip
 
      gg_g(2,2,1,2,1)=czip
 
      gg_g(2,2,1,1,2)=czip
 
      gg_g(2,1,2,2,2)=czip
 
      gg_g(2,1,2,1,1)=czip
 
      gg_g(2,1,2,2,1)=czip
 
      gg_g(2,1,2,1,2)=czip
 
      gg_g(1,2,1,2,2)=czip
 
      gg_g(1,2,1,1,1)=czip
 
      gg_g(1,2,1,2,1)=czip
 
      gg_g(1,2,1,1,2)=czip
 
      gg_g(1,1,2,2,2)=czip
 
      gg_g(1,1,2,1,1)=czip
 
      gg_g(1,1,2,2,1)=czip
 
      gg_g(1,1,2,1,2)=czip
 
      gg_g(2,1,1,2,2)=
     & +za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p4,t6)/zb(p1,p2)/zb(p1,p2)
     & * (-1._dp/s34/s126)
 
     & +za(p3,t5)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)/zb(p1,p2)
     & * (1._dp/s34/s126)
 
      gg_g(2,1,1,1,1)=+ za(p1,t6)*za(p2,p3)*zb(p1,p4)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq**2/s34/s126/s25)
 
      gg_g(2,1,1,2,1)=+za(p1,t6)*za(p3,t5)*zb(p1,p4)/zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s126)
 
      gg_g(2,1,1,1,2)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s126/s25)
 
     & +za(p2,p3)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s126/s25)
 
      gg_g(2,2,2,2,2)=+ za(p1,t6)*za(p3,t5)*zb(p1,t6)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (-1._dp/s34/s126)
 
     & +za(p3,t5)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     & * (1._dp/s34/s126)
 
      gg_g(2,2,2,1,1)=+ za(p1,t6)*za(p2,p3)*zb(p1,p4)*zb(p2,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq**2/s34/s126/s25)
 
      gg_g(2,2,2,2,1)=+za(p1,t6)*za(p3,t5)*zb(p1,p4)/za(p1,p2)/za(p1,p2)
     & * (mq/s34/s126)
 
      gg_g(2,2,2,1,2)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,t6)*zb(p2,t5)*zb(p4,t6)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s126/s25)
 
     & +za(p2,p3)*zb(p1,t6)*zb(p2,p4)*zb(p2,t5)/za(p1,p2)
     & * (mq/s34/s126/s25)
 
      gg_g(1,1,1,2,2)=+ za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq**2/s34/s126/s25)
 
      gg_g(1,1,1,1,1)=+ za(p1,t6)*za(p2,p3)*zb(p4,t5)/zb(p1,p2)
     & * (1._dp/s34/s126)
 
     & +za(p1,t6)*za(p3,t6)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)/zb(p1,p2)
     & * (-1._dp/s34/s126)
 
      gg_g(1,1,1,2,1)=
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p2,p4)/zb(p1,p2)
     & * (mq/s34/s126/s25)
 
     & +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s126/s25)
 
      gg_g(1,1,1,1,2)=+ za(p1,p3)*zb(p1,t6)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s126)
 
      gg_g(1,2,2,2,2)=+ za(p1,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (mq**2/s34/s126/s25)
 
      gg_g(1,2,2,1,1)=
     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p4,t5)/za(p1,p2)/za(p1,p2)
     & * (1._dp/s34/s126)
 
     & +za(p1,t6)*za(p3,t6)*zb(p1,t6)*zb(p4,t5)/za(p1,p2)/za(p1,p2)
     & * (-1._dp/s34/s126)
 
      gg_g(1,2,2,2,1)=
     & +za(p1,t6)*za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s126/s25)
 
     & +za(p1,t6)*za(p3,t6)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s126/s25)
 
      gg_g(1,2,2,1,2)=+ za(p1,p3)*zb(p1,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s126)
 

c      include 'hh2m.f'
 
      gg_h(2,2,1,2,2)=czip
 
      gg_h(2,2,1,1,1)=czip
 
      gg_h(2,2,1,2,1)=czip
 
      gg_h(2,2,1,1,2)=czip
 
      gg_h(2,1,2,2,2)=czip
 
      gg_h(2,1,2,1,1)=czip
 
      gg_h(2,1,2,2,1)=czip
 
      gg_h(2,1,2,1,2)=czip
 
      gg_h(1,2,1,2,2)=czip
 
      gg_h(1,2,1,1,1)=czip
 
      gg_h(1,2,1,2,1)=czip
 
      gg_h(1,2,1,1,2)=czip
 
      gg_h(1,1,2,2,2)=czip
 
      gg_h(1,1,2,1,1)=czip
 
      gg_h(1,1,2,2,1)=czip
 
      gg_h(1,1,2,1,2)=czip
 
      gg_h(2,1,1,2,2)=
     & +za(p1,p3)*za(p2,t5)*zb(p4,t6)/zb(p1,p2)
     & * (mq**2/s34/s125/s25)
 
     & +za(p2,p3)*za(p1,t5)*zb(p4,t6)/zb(p1,p2)
     & * (- mq**2/s34/s125/s25 -1._dp/s34/s125)
 
     & +za(p3,t5)*za(p1,t5)*zb(p1,t5)*zb(p4,t6)/zb(p1,p2)/zb(p1,p2)
     & * (1._dp/s34/s125)
 
      gg_h(2,1,1,1,1)=
     & +za(p1,t6)*za(p1,p3)*zb(p1,p4)*zb(p1,t5)/zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s125/s16)
 
     & +za(p1,t6)*za(p2,p3)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)
     & * (mq**2/s34/s125/s16/s25 + mq**4/s34/s125/s16
     & /s25**2)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s25)
 
      gg_h(2,1,1,2,1)=
     & +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p4)/zb(p1,p2)
     & * (mq**3/s34/s125/s16/s25)
 
     & +za(p1,t6)*za(p2,p3)*za(p1,t5)*zb(p1,p4)/zb(p1,p2)
     & * (- mq/s34/s125/s16 - mq**3/s34/s125/s16/s25)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s125/s16)
 
      gg_h(2,1,1,1,2)=
     & +za(p1,p3)*zb(p1,t5)*zb(p4,t6)/zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125)
 
     & +za(p2,p3)*za(p1,p2)*zb(p2,t5)*zb(p4,t6)/zb(p1,p2)
     & * (mq/s34/s125/s25 + mq**3/s34/s125/s25**2)
 
     & +za(p3,t5)*za(p1,p2)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125/s25)
 
      gg_h(2,2,2,2,2)=
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     & * (mq**2/s34/s125/s25)
 
     & +za(p2,p3)*za(p1,t5)*zb(p1,p2)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s125/s25 -1._dp/s34/s125)
 
     & +za(p3,t5)*za(p1,t5)*zb(p1,t5)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     & * (1._dp/s34/s125)
 
      gg_h(2,2,2,1,1)=
     & +za(p1,t6)*za(p1,p3)*zb(p1,p4)*zb(p1,t5)/za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s125/s16)
 
     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p1,p4)*zb(p2,t5)/za(p1,p2)
     & * (mq**2/s34/s125/s16/s25 + mq**4/s34/s125/s16
     & /s25**2)
 
     & +za(p1,t6)*za(p3,t5)*zb(p1,p4)*zb(p1,t5)*zb(p2,t5)/za(p1,p2)
     & * (- mq**2/s34/s125/s16/s25)
 
      gg_h(2,2,2,2,1)=
     & +za(p1,t6)*za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (mq**3/s34/s125/s16/s25)
 
     & +za(p1,t6)*za(p2,p3)*za(p1,t5)*zb(p1,p2)*zb(p1,p4)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125/s16 - mq**3/s34/s125/s16/s25)
 
     & +za(p1,t6)*za(p3,t5)*za(p1,t5)*zb(p1,p4)*zb(p1,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s125/s16)
 
      gg_h(2,2,2,1,2)=
     & +za(p1,p3)*zb(p1,t5)*zb(p4,t6)/za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125)
 
     & +za(p2,p3)*zb(p1,p2)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     & * (mq/s34/s125/s25 + mq**3/s34/s125/s25**2)
 
     & +za(p3,t5)*zb(p1,t5)*zb(p2,t5)*zb(p4,t6)/za(p1,p2)
     & * (- mq/s34/s125/s25)
 
      gg_h(1,1,1,2,2)=
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,t6)*zb(p4,t5)/zb(p1,p2)
     & * (- mq**2/s34/s125/s16/s25)
 
     & +za(p1,p3)*za(p1,t5)*zb(p1,p4)*zb(p1,t6)/zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s125/s16)
 
     & +za(p1,p3)*za(p2,t5)*za(p1,p2)*zb(p1,t6)*zb(p2,p4)/zb(p1,p2)
     & * (mq**2/s34/s125/s16/s25 + mq**4/s34/s125/s16
     & /s25**2)
 
      gg_h(1,1,1,1,1)=
     & +za(p3,t6)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/zb(p1,p2)/zb(p1,p2)
     & * (mq**2/s34/s125/s25)
 
     & +za(p3,t6)*za(p1,p2)*zb(p1,t5)*zb(p2,p4)/zb(p1,p2)/zb(p1,p2)
     & * (- mq**2/s34/s125/s25 -1._dp/s34/s125)
 
     & +za(p3,t6)*za(p1,t5)*zb(p1,t5)*zb(p4,t5)/zb(p1,p2)/zb(p1,p2)
     & * (1._dp/s34/s125)
 
      gg_h(1,1,1,2,1)=
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p4,t5)/zb(p1,p2)
     & * (- mq/s34/s125/s25)
 
     & +za(p3,t6)*za(p1,t5)*zb(p1,p4)/zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125)
 
     & +za(p3,t6)*za(p2,t5)*za(p1,p2)*zb(p2,p4)/zb(p1,p2)
     & * (mq/s34/s125/s25 + mq**3/s34/s125/s25**2)
 
      gg_h(1,1,1,1,2)=
     & +za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq**3/s34/s125/s16/s25)
 
     & +za(p1,p3)*za(p1,p2)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)
     & /zb(p1,p2)/zb(p1,p2)
     & * (- mq/s34/s125/s16 - mq**3/s34/s125/s16/s25)
 
     & +za(p1,p3)*za(p1,t5)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     & /zb(p1,p2)/zb(p1,p2)
     & * (mq/s34/s125/s16)
 
      gg_h(1,2,2,2,2)=
     & +za(p1,p3)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s125/s16/s25)
 
     & +za(p1,p3)*za(p1,t5)*zb(p1,p4)*zb(p1,t6)/za(p1,p2)/za(p1,p2)
     & * (- mq**2/s34/s125/s16)
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     & * (mq**2/s34/s125/s16/s25 + mq**4/s34/s125/s16
     & /s25**2)
 
      gg_h(1,2,2,1,1)=
     & +za(p3,t6)*za(p1,t5)*zb(p1,t5)*zb(p4,t5)/za(p1,p2)/za(p1,p2)
     & * (1._dp/s34/s125)
 
     & +za(p3,t6)*zb(p1,p4)*zb(p2,t5)/za(p1,p2)
     & * (mq**2/s34/s125/s25)
 
     & +za(p3,t6)*zb(p1,t5)*zb(p2,p4)/za(p1,p2)
     & * (- mq**2/s34/s125/s25 -1._dp/s34/s125)
 
      gg_h(1,2,2,2,1)=
     & +za(p3,t6)*za(p1,t5)*za(p2,t5)*zb(p1,p2)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125/s25)
 
     & +za(p3,t6)*za(p1,t5)*zb(p1,p4)/za(p1,p2)/za(p1,p2)
     & * (- mq/s34/s125)
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p2)*zb(p2,p4)/za(p1,p2)
     & * (mq/s34/s125/s25 + mq**3/s34/s125/s25**2)
 
      gg_h(1,2,2,1,2)=
     & +za(p1,p3)*za(p1,t5)*zb(p1,t5)*zb(p1,t6)*zb(p4,t5)
     & /za(p1,p2)/za(p1,p2)
     & * (mq/s34/s125/s16)
 
     & +za(p1,p3)*zb(p1,p4)*zb(p1,t6)*zb(p2,t5)/za(p1,p2)
     & * (mq**3/s34/s125/s16/s25)
 
     & +za(p1,p3)*zb(p1,t5)*zb(p1,t6)*zb(p2,p4)/za(p1,p2)
     & * (- mq/s34/s125/s16 - mq**3/s34/s125/s16/s25)
 

 
c      do hz=1,2      
c      do h1=1,2      
c      do h2=1,2      
c     do h5=1,2      
c      do h6=1,2      
c      write(6,*) 'gg_g(hz,h1,h2,h5,h6)+gg_h(hz,h2,h1,h5,h6)',
c     &            gg_g(hz,h1,h2,h5,h6)+gg_h(hz,h2,h1,h5,h6)     
c      enddo
c      enddo
c      enddo
c      enddo
c     enddo

c      pause
      return
      end
      
