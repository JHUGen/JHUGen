      subroutine mamps(mq,p1,p2,p3,p4,t5,t6,qqb_a,qqb_b)
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
      integer:: h1,h3,h5,h6
      complex(dp):: qqb_a(2,2,2,2,2),qqb_b(2,2,2,2,2)
      real(dp):: al5,al6,s56,s125,s126,s134,s234,propa,propb,bit
      real(dp):: mq
      integer,parameter::swap(2)=(/2,1/)

      al5=mq**2/s(p2,t5)
      s125=(1._dp+al5)*s(p1,p2)+s(p1,t5)+s(p2,t5)
      al6=mq**2/s(p1,t6)
      s126=(1._dp+al6)*s(p1,p2)+s(p1,t6)+s(p2,t6)
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s56=s134+s234+s(p1,p2)-s(p3,p4)
      propa=s56*s(p3,p4)
      propb=s(p1,p2)*s(p3,p4)

c---overall factor of 4 removed

C   Notation is hz,h1,h3,h5,h6 
C   These are the diagrams where Z attaches to the massless 
C   quark line

 
      qqb_a(1,1,1,1,1)=
     & +za(p1,p3)*za(p2,t6)*zb(p1,p4)*zb(p1,t5)
     & *(- 1._dp/(s134*propa))
 
     & +za(p2,p3)*za(p2,t6)*zb(p1,t5)*zb(p2,p4)
     & *(1._dp/(s234*propa))
 
     & +za(p2,p3)*za(p3,t6)*zb(p1,t5)*zb(p3,p4)
     & *(1._dp/(s234*propa))
 
     & +za(p3,p4)*za(p2,t6)*zb(p1,p4)*zb(p4,t5)
     & *(1._dp/(s134*propa))
 
      qqb_a(1,1,1,1,2)=
     & +za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)/za(p1,t6)
     & *(mq/(s134*propa))
 
     & +za(p2,p3)*za(p1,p2)*zb(p1,t5)*zb(p2,p4)/za(p1,t6)
     & *(- mq/(s234*propa))
 
     & +za(p2,p3)*za(p1,p3)*zb(p1,t5)*zb(p3,p4)/za(p1,t6)
     & *(- mq/(s234*propa))
 
     & +za(p2,p3)*za(p2,p3)*zb(p1,t6)*zb(p3,p4)/za(p2,t5)
     & *(- mq/(s234*propa))
 
     & +za(p3,p4)*za(p1,p2)*zb(p1,p4)*zb(p4,t5)/za(p1,t6)
     & *(- mq/(s134*propa))
 
      qqb_a(1,1,1,2,1)=
     & +za(p1,p3)*za(p2,t6)*zb(p1,p2)*zb(p1,p4)/zb(p2,t5)
     & *(- mq/(s134*propa))
 
     & +za(p2,p3)*za(p2,t6)*zb(p1,p2)*zb(p2,p4)/zb(p2,t5)
     & *(mq/(s234*propa))
 
     & +za(p2,p3)*za(p3,t6)*zb(p1,p2)*zb(p3,p4)/zb(p2,t5)
     & *(mq/(s234*propa))
 
     & +za(p3,p4)*za(p2,t5)*zb(p1,p4)*zb(p1,p4)/zb(p1,t6)
     & *(- mq/(s134*propa))
 
     & +za(p3,p4)*za(p2,t6)*zb(p1,p4)*zb(p2,p4)/zb(p2,t5)
     & *(- mq/(s134*propa))
 
      qqb_a(1,1,1,2,2)=
     & +za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p1,p4)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s134*propa))
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p4)*zb(p1,t6)
     & *(- 1._dp/(s134*propa))
 
     & +za(p2,p3)*za(p1,p2)*zb(p1,p2)*zb(p2,p4)/za(p1,t6)/zb(p2,t5)
     & *(- mq**2/(s234*propa))
 
     & +za(p2,p3)*za(p1,p3)*zb(p1,p2)*zb(p3,p4)/za(p1,t6)/zb(p2,t5)
     & *(- mq**2/(s234*propa))
 
     & +za(p2,p3)*za(p2,t5)*zb(p1,t6)*zb(p2,p4)
     & *(1._dp/(s234*propa))
 
     & +za(p2,p3)*za(p3,t5)*zb(p1,t6)*zb(p3,p4)
     & *(1._dp/(s234*propa))
 
     & +za(p3,p4)*za(p1,p2)*zb(p1,p4)*zb(p2,p4)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s134*propa))
 
     & +za(p3,p4)*za(p2,t5)*zb(p1,p4)*zb(p4,t6)
     & *(1._dp/(s134*propa))
 
      qqb_a(1,1,2,1,1)=
     & +za(p1,p4)*za(p2,t6)*zb(p1,p3)*zb(p1,t5)
     & *(- 1._dp/(s134*propa))
 
     & +za(p2,p4)*za(p2,t6)*zb(p1,t5)*zb(p2,p3)
     & *(1._dp/(s234*propa))
 
     & +za(p2,p4)*za(p4,t6)*zb(p1,t5)*zb(p3,p4)
     & *(- 1._dp/(s234*propa))
 
     & +za(p3,p4)*za(p2,t6)*zb(p1,p3)*zb(p3,t5)
     & *(- 1._dp/(s134*propa))
 
      qqb_a(1,1,2,1,2)=
     & +za(p1,p4)*za(p1,p2)*zb(p1,p3)*zb(p1,t5)/za(p1,t6)
     & *(mq/(s134*propa))
 
     & +za(p2,p4)*za(p1,p2)*zb(p1,t5)*zb(p2,p3)/za(p1,t6)
     & *(- mq/(s234*propa))
 
     & +za(p2,p4)*za(p1,p4)*zb(p1,t5)*zb(p3,p4)/za(p1,t6)
     & *(mq/(s234*propa))
 
     & +za(p2,p4)*za(p2,p4)*zb(p1,t6)*zb(p3,p4)/za(p2,t5)
     & *(mq/(s234*propa))
 
     & +za(p3,p4)*za(p1,p2)*zb(p1,p3)*zb(p3,t5)/za(p1,t6)
     & *(mq/(s134*propa))
 
      qqb_a(1,1,2,2,1)=
     & +za(p1,p4)*za(p2,t6)*zb(p1,p2)*zb(p1,p3)/zb(p2,t5)
     & *(- mq/(s134*propa))
 
     & +za(p2,p4)*za(p2,t6)*zb(p1,p2)*zb(p2,p3)/zb(p2,t5)
     & *(mq/(s234*propa))
 
     & +za(p2,p4)*za(p4,t6)*zb(p1,p2)*zb(p3,p4)/zb(p2,t5)
     & *(- mq/(s234*propa))
 
     & +za(p3,p4)*za(p2,t5)*zb(p1,p3)*zb(p1,p3)/zb(p1,t6)
     & *(mq/(s134*propa))
 
     & +za(p3,p4)*za(p2,t6)*zb(p1,p3)*zb(p2,p3)/zb(p2,t5)
     & *(mq/(s134*propa))
 
      qqb_a(1,1,2,2,2)=
     & +za(p1,p4)*za(p1,p2)*zb(p1,p2)*zb(p1,p3)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s134*propa))
 
     & +za(p1,p4)*za(p2,t5)*zb(p1,p3)*zb(p1,t6)
     & *(- 1._dp/(s134*propa))
 
     & +za(p2,p4)*za(p1,p2)*zb(p1,p2)*zb(p2,p3)/za(p1,t6)/zb(p2,t5)
     & *(- mq**2/(s234*propa))
 
     & +za(p2,p4)*za(p1,p4)*zb(p1,p2)*zb(p3,p4)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s234*propa))
 
     & +za(p2,p4)*za(p2,t5)*zb(p1,t6)*zb(p2,p3)
     & *(1._dp/(s234*propa))
 
     & +za(p2,p4)*za(p4,t5)*zb(p1,t6)*zb(p3,p4)
     & *(- 1._dp/(s234*propa))
 
     & +za(p3,p4)*za(p1,p2)*zb(p1,p3)*zb(p2,p3)/za(p1,t6)/zb(p2,t5)
     & *(- mq**2/(s134*propa))
 
     & +za(p3,p4)*za(p2,t5)*zb(p1,p3)*zb(p3,t6)
     & *(- 1._dp/(s134*propa))
 
c      qqb_a(2,2,1,1,1)=
c     & +za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p1,p4)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s134*propa))
 
c     & +za(p1,p3)*za(p1,t6)*zb(p1,p4)*zb(p2,t5)
c     & *(- 1._dp/(s134*propa))
 
c     & +za(p1,p3)*za(p2,p3)*zb(p1,p2)*zb(p3,p4)/za(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s134*propa))
 
c     & +za(p1,p3)*za(p3,t6)*zb(p2,t5)*zb(p3,p4)
c     & *(- 1._dp/(s134*propa))
 
c     & +za(p2,p3)*za(p1,p2)*zb(p1,p2)*zb(p2,p4)/za(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s234*propa))
 
c     & +za(p2,p3)*za(p1,t6)*zb(p2,p4)*zb(p2,t5)
c     & *(1._dp/(s234*propa))
 
c     & +za(p3,p4)*za(p1,p2)*zb(p1,p4)*zb(p2,p4)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t6)*zb(p2,p4)*zb(p4,t5)
c     & *(- 1._dp/(s234*propa))
 
c      qqb_a(2,2,1,1,2)=
c     & +za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p2,t6)/za(p2,t5)
c     & *(- mq/(s134*propa))
 
c     & +za(p1,p3)*za(p1,p3)*zb(p2,t5)*zb(p3,p4)/za(p1,t6)
c     & *(mq/(s134*propa))
 
c     & +za(p1,p3)*za(p2,p3)*zb(p2,t6)*zb(p3,p4)/za(p2,t5)
c     & *(mq/(s134*propa))
 
c     & +za(p2,p3)*za(p1,p2)*zb(p2,p4)*zb(p2,t6)/za(p2,t5)
c     & *(mq/(s234*propa))
 
c     & +za(p3,p4)*za(p1,p2)*zb(p2,p4)*zb(p4,t6)/za(p2,t5)
c     & *(- mq/(s234*propa))
 
c      qqb_a(2,2,1,2,1)=
c     & +za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p1,p4)/zb(p1,t6)
c     & *(mq/(s134*propa))
 
c     & +za(p1,p3)*za(p3,t5)*zb(p1,p2)*zb(p3,p4)/zb(p1,t6)
c     & *(mq/(s134*propa))
 
c     & +za(p2,p3)*za(p1,t5)*zb(p1,p2)*zb(p2,p4)/zb(p1,t6)
c     & *(- mq/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t5)*zb(p1,p4)*zb(p2,p4)/zb(p1,t6)
c     & *(mq/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t6)*zb(p2,p4)*zb(p2,p4)/zb(p2,t5)
c     & *(mq/(s234*propa))
 
c      qqb_a(2,2,1,2,2)=
c     & +za(p1,p3)*za(p1,t5)*zb(p1,p4)*zb(p2,t6)
c     & *(- 1._dp/(s134*propa))
 
c     & +za(p1,p3)*za(p3,t5)*zb(p2,t6)*zb(p3,p4)
c     & *(- 1._dp/(s134*propa))
 
c     & +za(p2,p3)*za(p1,t5)*zb(p2,p4)*zb(p2,t6)
c     & *(1._dp/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t5)*zb(p2,p4)*zb(p4,t6)
c     & *(- 1._dp/(s234*propa))
 
c      qqb_a(2,2,2,1,1)=
c     & +za(p1,p4)*za(p1,p2)*zb(p1,p2)*zb(p1,p3)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s134*propa))
 
c     & +za(p1,p4)*za(p1,t6)*zb(p1,p3)*zb(p2,t5)
c     & *(- 1._dp/(s134*propa))
 
c     & +za(p1,p4)*za(p2,p4)*zb(p1,p2)*zb(p3,p4)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s134*propa))
 
c     & +za(p1,p4)*za(p4,t6)*zb(p2,t5)*zb(p3,p4)
c     & *(1._dp/(s134*propa))
 
c     & +za(p2,p4)*za(p1,p2)*zb(p1,p2)*zb(p2,p3)/za(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s234*propa))
 
c     & +za(p2,p4)*za(p1,t6)*zb(p2,p3)*zb(p2,t5)
c     & *(1._dp/(s234*propa))
 
c     & +za(p3,p4)*za(p1,p2)*zb(p1,p3)*zb(p2,p3)/za(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t6)*zb(p2,p3)*zb(p3,t5)
c     & *(1._dp/(s234*propa))
 
c      qqb_a(2,2,2,1,2)=
c     & +za(p1,p4)*za(p1,p2)*zb(p1,p3)*zb(p2,t6)/za(p2,t5)
c     & *(- mq/(s134*propa))
 
c     & +za(p1,p4)*za(p1,p4)*zb(p2,t5)*zb(p3,p4)/za(p1,t6)
c     & *(- mq/(s134*propa))
 
c     & +za(p1,p4)*za(p2,p4)*zb(p2,t6)*zb(p3,p4)/za(p2,t5)
c     & *(- mq/(s134*propa))
 
c     & +za(p2,p4)*za(p1,p2)*zb(p2,p3)*zb(p2,t6)/za(p2,t5)
c     & *(mq/(s234*propa))
 
c     & +za(p3,p4)*za(p1,p2)*zb(p2,p3)*zb(p3,t6)/za(p2,t5)
c     & *(mq/(s234*propa))
 
c      qqb_a(2,2,2,2,1)=
c     & +za(p1,p4)*za(p1,t5)*zb(p1,p2)*zb(p1,p3)/zb(p1,t6)
c     & *(mq/(s134*propa))
 
c     & +za(p1,p4)*za(p4,t5)*zb(p1,p2)*zb(p3,p4)/zb(p1,t6)
c     & *(- mq/(s134*propa))
 
c     & +za(p2,p4)*za(p1,t5)*zb(p1,p2)*zb(p2,p3)/zb(p1,t6)
c     & *(- mq/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t5)*zb(p1,p3)*zb(p2,p3)/zb(p1,t6)
c     & *(- mq/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t6)*zb(p2,p3)*zb(p2,p3)/zb(p2,t5)
c     & *(- mq/(s234*propa))
 
c      qqb_a(2,2,2,2,2)=
c     & +za(p1,p4)*za(p1,t5)*zb(p1,p3)*zb(p2,t6)
c     & *(- 1._dp/(s134*propa))
 
c     & +za(p1,p4)*za(p4,t5)*zb(p2,t6)*zb(p3,p4)
c     & *(1._dp/(s134*propa))
 
c     & +za(p2,p4)*za(p1,t5)*zb(p2,p3)*zb(p2,t6)
c     & *(1._dp/(s234*propa))
 
c     & +za(p3,p4)*za(p1,t5)*zb(p2,p3)*zb(p3,t6)
c     & *(1._dp/(s234*propa))
 
      qqb_a(1,2,1,1,1)=czip
 
      qqb_a(1,2,1,1,2)=czip
 
      qqb_a(1,2,1,2,1)=czip
 
      qqb_a(1,2,1,2,2)=czip
 
      qqb_a(1,2,2,1,1)=czip
 
      qqb_a(1,2,2,1,2)=czip
 
      qqb_a(1,2,2,2,1)=czip
 
      qqb_a(1,2,2,2,2)=czip
 
c      qqb_a(2,1,1,1,1)=czip
 
c      qqb_a(2,1,1,1,2)=czip
 
c      qqb_a(2,1,1,2,1)=czip
 
c      qqb_a(2,1,1,2,2)=czip
 
c      qqb_a(2,1,2,1,1)=czip
 
c      qqb_a(2,1,2,1,2)=czip
 
c      qqb_a(2,1,2,2,1)=czip
 
c      qqb_a(2,1,2,2,2)=czip
 
  
      qqb_b(1,1,1,1,1)=
     & +za(p2,t6)*za(p2,p3)*zb(p1,p2)*zb(p4,t5)
     & *(- 1._dp/(s126*propb))
 
     & +za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p4,t5)
     & *(1._dp/(s126*propb))
 
     & +za(p3,t6)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)
     & *(- 1._dp/(s125*propb))
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,t5)*zb(p4,t5)
     & *(- 1._dp/(s125*propb))
 
      qqb_b(1,1,1,1,2)=
     & +za(p1,p2)*za(p2,p3)*zb(p1,p2)*zb(p4,t5)/za(p1,t6)
     & *(mq/(s126*propb))
 
     & +za(p1,p2)*za(p3,t6)*zb(p1,t6)*zb(p4,t5)/za(p1,t6)
     & *(- mq/(s126*propb))
 
     & +za(p1,p3)*za(p1,p2)*zb(p1,p4)*zb(p1,t5)/za(p1,t6)
     & *(mq/(s125*propb))
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,t5)*zb(p4,t5)/za(p1,t6)
     & *(mq/(s125*propb))
 
     & +za(p2,p3)*zb(p1,t6)*zb(p4,t5)
     & *(- mq/(s126*propb))
 
      qqb_b(1,1,1,2,1)=
     & +za(p2,t6)*za(p2,p3)*zb(p1,p2)*zb(p2,p4)/zb(p2,t5)
     & *(mq/(s126*propb))
 
     & +za(p2,t6)*za(p3,t6)*zb(p1,t6)*zb(p2,p4)/zb(p2,t5)
     & *(- mq/(s126*propb))
 
     & +za(p3,t6)*za(p1,p2)*zb(p1,p2)*zb(p1,p4)/zb(p2,t5)
     & *(- mq/(s125*propb))
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p2)*zb(p4,t5)/zb(p2,t5)
     & *(- mq/(s125*propb))
 
     & +za(p3,t6)*za(p2,t5)*zb(p1,p4)
     & *(mq/(s125*propb))
 
      qqb_b(1,1,1,2,2)=
     & +za(p1,p2)*za(p2,p3)*zb(p1,p2)*zb(p2,p4)/za(p1,t6)/zb(p2,t5)
     & *(- mq**2/(s126*propb))
 
     & +za(p1,p2)*za(p3,t6)*zb(p1,t6)*zb(p2,p4)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s126*propb))
 
     & +za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p1,p4)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s125*propb))
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p2)*zb(p4,t5)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s125*propb))
 
     & +za(p1,p3)*za(p2,t5)*zb(p1,p4)/za(p1,t6)
     & *(- mq**2/(s125*propb))
 
     & +za(p2,p3)*zb(p1,t6)*zb(p2,p4)/zb(p2,t5)
     & *(mq**2/(s126*propb))
 
      qqb_b(1,1,2,1,1)=
     & +za(p2,t6)*za(p2,p4)*zb(p1,p2)*zb(p3,t5)
     & *(- 1._dp/(s126*propb))
 
     & +za(p2,t6)*za(p4,t6)*zb(p1,t6)*zb(p3,t5)
     & *(1._dp/(s126*propb))
 
     & +za(p4,t6)*za(p1,p2)*zb(p1,p3)*zb(p1,t5)
     & *(- 1._dp/(s125*propb))
 
     & +za(p4,t6)*za(p2,t5)*zb(p1,t5)*zb(p3,t5)
     & *(- 1._dp/(s125*propb))
 
      qqb_b(1,1,2,1,2)=
     & +za(p1,p2)*za(p2,p4)*zb(p1,p2)*zb(p3,t5)/za(p1,t6)
     & *(mq/(s126*propb))
 
     & +za(p1,p2)*za(p4,t6)*zb(p1,t6)*zb(p3,t5)/za(p1,t6)
     & *(- mq/(s126*propb))
 
     & +za(p1,p4)*za(p1,p2)*zb(p1,p3)*zb(p1,t5)/za(p1,t6)
     & *(mq/(s125*propb))
 
     & +za(p1,p4)*za(p2,t5)*zb(p1,t5)*zb(p3,t5)/za(p1,t6)
     & *(mq/(s125*propb))
 
     & +za(p2,p4)*zb(p1,t6)*zb(p3,t5)
     & *(- mq/(s126*propb))
 
      qqb_b(1,1,2,2,1)=
     & +za(p2,t6)*za(p2,p4)*zb(p1,p2)*zb(p2,p3)/zb(p2,t5)
     & *(mq/(s126*propb))
 
     & +za(p2,t6)*za(p4,t6)*zb(p1,t6)*zb(p2,p3)/zb(p2,t5)
     & *(- mq/(s126*propb))
 
     & +za(p4,t6)*za(p1,p2)*zb(p1,p2)*zb(p1,p3)/zb(p2,t5)
     & *(- mq/(s125*propb))
 
     & +za(p4,t6)*za(p2,t5)*zb(p1,p2)*zb(p3,t5)/zb(p2,t5)
     & *(- mq/(s125*propb))
 
     & +za(p4,t6)*za(p2,t5)*zb(p1,p3)
     & *(mq/(s125*propb))
 
      qqb_b(1,1,2,2,2)=
     & +za(p1,p2)*za(p2,p4)*zb(p1,p2)*zb(p2,p3)/za(p1,t6)/zb(p2,t5)
     & *(- mq**2/(s126*propb))
 
     & +za(p1,p2)*za(p4,t6)*zb(p1,t6)*zb(p2,p3)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s126*propb))
 
     & +za(p1,p4)*za(p1,p2)*zb(p1,p2)*zb(p1,p3)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s125*propb))
 
     & +za(p1,p4)*za(p2,t5)*zb(p1,p2)*zb(p3,t5)/za(p1,t6)/zb(p2,t5)
     & *(mq**2/(s125*propb))
 
     & +za(p1,p4)*za(p2,t5)*zb(p1,p3)/za(p1,t6)
     & *(- mq**2/(s125*propb))
 
     & +za(p2,p4)*zb(p1,t6)*zb(p2,p3)/zb(p2,t5)
     & *(mq**2/(s126*propb))
 
c      qqb_b(2,2,1,1,1)=
c     & +za(p1,p2)*za(p2,p3)*zb(p1,p2)*zb(p2,p4)/za(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s126*propb))
 
c     & +za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p1,p4)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s125*propb))
 
c     & +za(p1,p3)*zb(p1,p4)*zb(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s125*propb))
 
c     & +za(p1,t6)*za(p2,p3)*zb(p1,p2)*zb(p4,t6)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s126*propb))
 
c     & +za(p1,t6)*za(p2,p3)*zb(p2,p4)/za(p2,t5)
c     & *(mq**2/(s126*propb))
 
c     & +za(p3,t5)*za(p1,p2)*zb(p1,p4)*zb(p2,t5)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s125*propb))
 
c      qqb_b(2,2,1,1,2)=
c     & +za(p1,p2)*za(p2,p3)*zb(p2,p4)*zb(p2,t6)/za(p2,t5)
c     & *(mq/(s126*propb))
 
c     & +za(p1,p3)*za(p1,p2)*zb(p1,p2)*zb(p4,t6)/za(p2,t5)
c     & *(- mq/(s125*propb))
 
c     & +za(p1,p3)*zb(p2,t5)*zb(p4,t6)
c     & *(mq/(s125*propb))
 
c     & +za(p1,t6)*za(p2,p3)*zb(p2,t6)*zb(p4,t6)/za(p2,t5)
c     & *(- mq/(s126*propb))
 
c     & +za(p3,t5)*za(p1,p2)*zb(p2,t5)*zb(p4,t6)/za(p2,t5)
c     & *(- mq/(s125*propb))
 
c      qqb_b(2,2,1,2,1)=
c     & +za(p1,p2)*za(p3,t5)*zb(p1,p2)*zb(p2,p4)/zb(p1,t6)
c     & *(mq/(s126*propb))
 
c     & +za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p1,p4)/zb(p1,t6)
c     & *(mq/(s125*propb))
 
c     & +za(p1,t6)*za(p3,t5)*zb(p1,p2)*zb(p4,t6)/zb(p1,t6)
c     & *(- mq/(s126*propb))
 
c     & +za(p1,t6)*za(p3,t5)*zb(p2,p4)
c     & *(- mq/(s126*propb))
 
c     & +za(p3,t5)*za(p1,t5)*zb(p1,p4)*zb(p2,t5)/zb(p1,t6)
c     & *(mq/(s125*propb))
 
c      qqb_b(2,2,1,2,2)=
c     & +za(p1,p2)*za(p3,t5)*zb(p2,p4)*zb(p2,t6)
c     & *(- 1._dp/(s126*propb))
 
c     & +za(p1,p3)*za(p1,t5)*zb(p1,p2)*zb(p4,t6)
c     & *(- 1._dp/(s125*propb))
 
c     & +za(p1,t6)*za(p3,t5)*zb(p2,t6)*zb(p4,t6)
c     & *(1._dp/(s126*propb))
 
c     & +za(p3,t5)*za(p1,t5)*zb(p2,t5)*zb(p4,t6)
c     & *(- 1._dp/(s125*propb))
 
c      qqb_b(2,2,2,1,1)=
c     & +za(p1,p2)*za(p2,p4)*zb(p1,p2)*zb(p2,p3)/za(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s126*propb))
 
c     & +za(p1,p4)*za(p1,p2)*zb(p1,p2)*zb(p1,p3)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s125*propb))
 
c     & +za(p1,p4)*zb(p1,p3)*zb(p2,t5)/zb(p1,t6)
c     & *(- mq**2/(s125*propb))
 
c     & +za(p1,t6)*za(p2,p4)*zb(p1,p2)*zb(p3,t6)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s126*propb))
 
c     & +za(p1,t6)*za(p2,p4)*zb(p2,p3)/za(p2,t5)
c     & *(mq**2/(s126*propb))
 
c     & +za(p4,t5)*za(p1,p2)*zb(p1,p3)*zb(p2,t5)/za(p2,t5)/zb(p1,t6)
c     & *(mq**2/(s125*propb))
 
c      qqb_b(2,2,2,1,2)=
c     & +za(p1,p2)*za(p2,p4)*zb(p2,p3)*zb(p2,t6)/za(p2,t5)
c     & *(mq/(s126*propb))
 
c     & +za(p1,p4)*za(p1,p2)*zb(p1,p2)*zb(p3,t6)/za(p2,t5)
c     & *(- mq/(s125*propb))
 
c     & +za(p1,p4)*zb(p2,t5)*zb(p3,t6)
c     & *(mq/(s125*propb))
 
c     & +za(p1,t6)*za(p2,p4)*zb(p2,t6)*zb(p3,t6)/za(p2,t5)
c     & *(- mq/(s126*propb))
 
c     & +za(p4,t5)*za(p1,p2)*zb(p2,t5)*zb(p3,t6)/za(p2,t5)
c     & *(- mq/(s125*propb))
 
c      qqb_b(2,2,2,2,1)=
c     & +za(p1,p2)*za(p4,t5)*zb(p1,p2)*zb(p2,p3)/zb(p1,t6)
c     & *(mq/(s126*propb))
 
c     & +za(p1,p4)*za(p1,t5)*zb(p1,p2)*zb(p1,p3)/zb(p1,t6)
c     & *(mq/(s125*propb))
 
c     & +za(p1,t6)*za(p4,t5)*zb(p1,p2)*zb(p3,t6)/zb(p1,t6)
c     & *(- mq/(s126*propb))
 
c     & +za(p1,t6)*za(p4,t5)*zb(p2,p3)
c     & *(- mq/(s126*propb))
 
c     & +za(p4,t5)*za(p1,t5)*zb(p1,p3)*zb(p2,t5)/zb(p1,t6)
c     & *(mq/(s125*propb))
 
c      qqb_b(2,2,2,2,2)=
c     & +za(p1,p2)*za(p4,t5)*zb(p2,p3)*zb(p2,t6)
c     & *(- 1._dp/(s126*propb))
 
c     & +za(p1,p4)*za(p1,t5)*zb(p1,p2)*zb(p3,t6)
c     & *(- 1._dp/(s125*propb))
 
c     & +za(p1,t6)*za(p4,t5)*zb(p2,t6)*zb(p3,t6)
c     & *(1._dp/(s126*propb))
 
c     & +za(p4,t5)*za(p1,t5)*zb(p2,t5)*zb(p3,t6)
c     & *(- 1._dp/(s125*propb))
 
      qqb_b(1,2,1,1,1)=
     & +za(p1,t6)*za(p1,p3)*zb(p1,p2)*zb(p4,t5)
     & *(1._dp/(s126*propb))
 
     & +za(p1,t6)*za(p3,t6)*zb(p2,t6)*zb(p4,t5)
     & *(1._dp/(s126*propb))
 
     & +za(p3,t6)*za(p1,p2)*zb(p2,p4)*zb(p2,t5)
     & *(1._dp/(s125*propb))
 
     & +za(p3,t6)*za(p1,t5)*zb(p2,t5)*zb(p4,t5)
     & *(- 1._dp/(s125*propb))
 
      qqb_b(1,2,1,1,2)=
     & +za(p1,p3)*za(p1,p2)*zb(p2,p4)*zb(p2,t5)/za(p1,t6)
     & *(- mq/(s125*propb))
 
     & +za(p1,p3)*za(p1,t5)*zb(p2,t5)*zb(p4,t5)/za(p1,t6)
     & *(mq/(s125*propb))
 
     & +za(p1,p3)*zb(p2,t6)*zb(p4,t5)
     & *(- mq/(s126*propb))
 
      qqb_b(1,2,1,2,1)=
     & +za(p1,t6)*za(p1,p3)*zb(p1,p2)*zb(p2,p4)/zb(p2,t5)
     & *(- mq/(s126*propb))
 
     & +za(p1,t6)*za(p3,t6)*zb(p2,p4)*zb(p2,t6)/zb(p2,t5)
     & *(- mq/(s126*propb))
 
     & +za(p3,t6)*za(p1,t5)*zb(p2,p4)
     & *(mq/(s125*propb))
 
      qqb_b(1,2,1,2,2)=+ za(p1,p3)*za(p1,t5)*zb(p2,p4)/za(p1,t6)
     & *(- mq**2/(s125*propb))
 
     & +za(p1,p3)*zb(p2,p4)*zb(p2,t6)/zb(p2,t5)
     & *(mq**2/(s126*propb))
 
      qqb_b(1,2,2,1,1)=
     & +za(p1,t6)*za(p1,p4)*zb(p1,p2)*zb(p3,t5)
     & *(1._dp/(s126*propb))
 
     & +za(p1,t6)*za(p4,t6)*zb(p2,t6)*zb(p3,t5)
     & *(1._dp/(s126*propb))
 
     & +za(p4,t6)*za(p1,p2)*zb(p2,p3)*zb(p2,t5)
     & *(1._dp/(s125*propb))
 
     & +za(p4,t6)*za(p1,t5)*zb(p2,t5)*zb(p3,t5)
     & *(- 1._dp/(s125*propb))
 
      qqb_b(1,2,2,1,2)=
     & +za(p1,p4)*za(p1,p2)*zb(p2,p3)*zb(p2,t5)/za(p1,t6)
     & *(- mq/(s125*propb))
 
     & +za(p1,p4)*za(p1,t5)*zb(p2,t5)*zb(p3,t5)/za(p1,t6)
     & *(mq/(s125*propb))
 
     & +za(p1,p4)*zb(p2,t6)*zb(p3,t5)
     & *(- mq/(s126*propb))
 
      qqb_b(1,2,2,2,1)=
     & +za(p1,t6)*za(p1,p4)*zb(p1,p2)*zb(p2,p3)/zb(p2,t5)
     & *(- mq/(s126*propb))
 
     & +za(p1,t6)*za(p4,t6)*zb(p2,p3)*zb(p2,t6)/zb(p2,t5)
     & *(- mq/(s126*propb))
 
     & +za(p4,t6)*za(p1,t5)*zb(p2,p3)
     & *(mq/(s125*propb))
 
      qqb_b(1,2,2,2,2)=+ za(p1,p4)*za(p1,t5)*zb(p2,p3)/za(p1,t6)
     & *(- mq**2/(s125*propb))
 
     & +za(p1,p4)*zb(p2,p3)*zb(p2,t6)/zb(p2,t5)
     & *(mq**2/(s126*propb))
 
c      qqb_b(2,1,1,1,1)=+ za(p2,p3)*zb(p1,p4)*zb(p1,t5)/zb(p1,t6)
c     & *(- mq**2/(s125*propb))
 
c     & +za(p2,t6)*za(p2,p3)*zb(p1,p4)/za(p2,t5)
c     & *(mq**2/(s126*propb))
 
c      qqb_b(2,1,1,1,2)=
c     & +za(p1,p2)*za(p2,p3)*zb(p1,p4)*zb(p1,t6)/za(p2,t5)
c     & *(- mq/(s126*propb))
 
c     & +za(p2,p3)*zb(p1,t5)*zb(p4,t6)
c     & *(mq/(s125*propb))
 
c     & +za(p2,t6)*za(p2,p3)*zb(p1,t6)*zb(p4,t6)/za(p2,t5)
c     & *(- mq/(s126*propb))
 
c      qqb_b(2,1,1,2,1)=
c     & +za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p1,p4)/zb(p1,t6)
c     & *(- mq/(s125*propb))
 
c     & +za(p2,t6)*za(p3,t5)*zb(p1,p4)
c     & *(- mq/(s126*propb))
 
c     & +za(p3,t5)*za(p2,t5)*zb(p1,p4)*zb(p1,t5)/zb(p1,t6)
c     & *(mq/(s125*propb))
 
c      qqb_b(2,1,1,2,2)=
c     & +za(p1,p2)*za(p3,t5)*zb(p1,p4)*zb(p1,t6)
c     & *(1._dp/(s126*propb))
 
c     & +za(p2,p3)*za(p2,t5)*zb(p1,p2)*zb(p4,t6)
c     & *(1._dp/(s125*propb))
 
c     & +za(p2,t6)*za(p3,t5)*zb(p1,t6)*zb(p4,t6)
c     & *(1._dp/(s126*propb))
 
c     & +za(p3,t5)*za(p2,t5)*zb(p1,t5)*zb(p4,t6)
c     & *(- 1._dp/(s125*propb))
 
c      qqb_b(2,1,2,1,1)=+ za(p2,p4)*zb(p1,p3)*zb(p1,t5)/zb(p1,t6)
c     & *(- mq**2/(s125*propb))
 
c     & +za(p2,t6)*za(p2,p4)*zb(p1,p3)/za(p2,t5)
c     & *(mq**2/(s126*propb))
 
c      qqb_b(2,1,2,1,2)=
c     & +za(p1,p2)*za(p2,p4)*zb(p1,p3)*zb(p1,t6)/za(p2,t5)
c     & *(- mq/(s126*propb))
 
c     & +za(p2,p4)*zb(p1,t5)*zb(p3,t6)
c     & *(mq/(s125*propb))
 
c     & +za(p2,t6)*za(p2,p4)*zb(p1,t6)*zb(p3,t6)/za(p2,t5)
c     & *(- mq/(s126*propb))
 
c      qqb_b(2,1,2,2,1)=
c     & +za(p2,p4)*za(p2,t5)*zb(p1,p2)*zb(p1,p3)/zb(p1,t6)
c     & *(- mq/(s125*propb))
 
c     & +za(p2,t6)*za(p4,t5)*zb(p1,p3)
c     & *(- mq/(s126*propb))
 
c     & +za(p4,t5)*za(p2,t5)*zb(p1,p3)*zb(p1,t5)/zb(p1,t6)
c     & *(mq/(s125*propb))
 
c      qqb_b(2,1,2,2,2)=
c     & +za(p1,p2)*za(p4,t5)*zb(p1,p3)*zb(p1,t6)
c     & *(1._dp/(s126*propb))
 
c     & +za(p2,p4)*za(p2,t5)*zb(p1,p2)*zb(p3,t6)
c     & *(1._dp/(s125*propb))
 
c     & +za(p2,t6)*za(p4,t5)*zb(p1,t6)*zb(p3,t6)
c     & *(1._dp/(s126*propb))
 
c     & +za(p4,t5)*za(p2,t5)*zb(p1,t5)*zb(p3,t6)
c     & *(- 1._dp/(s125*propb))
 


      do h1=1,2
      do h3=1,2
      do h5=1,2
      do h6=1,2
      bit=-1._dp
      if (h5==h6) bit=1._dp
      qqb_a(2,h1,h3,h5,h6)
     & =bit*conjg(qqb_a(1,swap(h1),swap(h3),swap(h5),swap(h6)))
      qqb_b(2,h1,h3,h5,h6)
     & =bit*conjg(qqb_b(1,swap(h1),swap(h3),swap(h5),swap(h6)))

      enddo
      enddo
      enddo
      enddo
      
      return
      end 





