      subroutine phi3(xth,xphi,p1,p2,p3,wt)
      implicit none
      include 'types.f'
c     d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     particle (p1) with mass s decaying into p2 and p3 massless.
c     vectors p2 and p3 are returned in the same frame as p1 is supplied.
c result is 1/8/pi * 2|p|/sqrts  * domega/(4*pi)
      real(dp):: p1(4),p2(4),p3(4),p3cm(4)
      real(dp):: xth,xphi,s,roots,costh,sinth,phi,cphi,sphi,pi
      real(dp):: wt0,wt
      integer:: j
      parameter(pi=3.141592654_dp,wt0=1._dp/8._dp/pi)
      s=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      roots=sqrt(s)
      costh=2._dp*xth-1._dp
      phi=2._dp*pi*xphi
      sinth=sqrt(1._dp-costh**2)
      cphi=cos(phi)
      sphi=sin(phi)
      wt=wt0

      p3cm(4)=roots/2._dp
      p3cm(1)=p3cm(4)*sinth*sphi
      p3cm(2)=p3cm(4)*sinth*cphi
      p3cm(3)=p3cm(4)*costh

      call boost(roots,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo

      return
      end

