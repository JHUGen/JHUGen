      subroutine ttWprod(p1,p2,m)
      implicit none
      include 'types.f'
      
C-----Calculate the amplitudes for the process
C-----q(p1)+qb(p2)--> t(p3)+a(p4)+l(p5)+l~(p6)
C-----amplitude stripped of an overall factor of gw^2*gs^2*t(CA)*T(CA)
C-----with Tr (t^ca t^cb)=\delta(ca,cb) 
C-----and   k3=+bp/beta*p3-bm/beta*p4
C-----      k4=-bm/beta*p3+bp/beta*p4
C-----      p3=bp*k5+bm*k6
C-----      p4=bm*k5+bp*k6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer:: p1,p2,p3,p4,p5,p6,k3,k4
      complex(dp):: m(2,2),zba2
      real(dp):: s156,s256
      parameter(p5=3,p6=4,k3=5,k4=6)
C---Statement function 
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
C---End statement function 
      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)
      s256=s(p2,p5)+s(p2,p6)+s(p5,p6)

      m(2,2)=mt/(s(k3,k4)*s(p5,p6)*s256) * (
     &      -za(p2,p5)*zb(p1,k3)/za(k3,k4)*zba2(p6,p2,p5,k3)
     &      +za(p2,p5)*zb(p1,k4)/za(k3,k4)*zba2(p6,p2,p5,k4))
     &      +mt/(s(k3,k4)*s(p5,p6)*s156) * (
     &      +za(p2,k3)*zb(p1,p6)/za(k3,k4)*zba2(k3,p1,p6,p5)
     &      -za(p2,k4)*zb(p1,p6)/za(k3,k4)*zba2(k4,p1,p6,p5))

      m(1,2)=
     & (-za(p2,p5)*zb(p1,k4)*zba2(p6,p2,p5,k3))/(s(k3,k4)*s(p5,p6)*s256)
     &+(+za(p2,k3)*zb(p1,p6)*zba2(k4,p1,p6,p5))/(s(k3,k4)*s(p5,p6)*s156)

       m(2,1) =
     & (-za(p2,p5)*zb(p1,k3)*zba2(p6,p2,p5,k4))/(s(k3,k4)*s(p5,p6)*s256)
     &+(+za(p2,k4)*zb(p1,p6)*zba2(k3,p1,p6,p5))/(s(k3,k4)*s(p5,p6)*s156)

       m(1,1)=mt/(s(k3,k4)*s(p5,p6)*s256) * (
     &       -za(p2,p5)*zb(p1,k3)/zb(k3,k4)*zba2(p6,p2,p5,k3)
     &       +za(p2,p5)*zb(p1,k4)/zb(k3,k4)*zba2(p6,p2,p5,k4))
     &       +mt/(s(k3,k4)*s(p5,p6)*s156) * (
     &       +za(p2,k3)*zb(p1,p6)/zb(k3,k4)*zba2(k3,p1,p6,p5)
     &       -za(p2,k4)*zb(p1,p6)/zb(k3,k4)*zba2(k4,p1,p6,p5))

      return
      end
