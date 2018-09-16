      subroutine qg_tbqdk_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Lowest order t-channel single top, with explicit b-quark         *
*       (contracted with vector n)                                     *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *
*                                                                      *
*         Author: J. Campbell, March 19, 2008                          *
*                         (added decay May 2011)                       *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'stopscales.f'
      real(dp):: p(mxpart,4),fac,msq_qg,msq_gq,msq_qbarg,msq_gqbar
      real(dp):: msq(-nf:nf,-nf:nf),gsq_H
      real(dp):: qg_tbqndk,n(4)
      integer:: j,k,in,i3,i4

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      gsq_H=fourpi*as_H
      fac=aveqg*xn**2*Cf*gsq_H*gwsq**4

c--- set labels of quark and antiquark according to nwz
      if (nwz == +1) then
        i3=3
        i4=4
      else
        i3=4
        i4=3
      endif

      msq_gq=0._dp
      msq_qg=0._dp
      msq_gqbar=0._dp
      msq_qbarg=0._dp
      if     (in == 1) then
c        msq_gq=qg_tbqndkold(2,1,7,i3,i4,p,n)
c        msq_gqbar=qg_tbqndkold(7,1,2,i3,i4,p,n)
        msq_gq=qg_tbqndk(2,1,7,i3,i4,p,n)
        msq_gqbar=qg_tbqndk(7,1,2,i3,i4,p,n)
      elseif (in == 2) then
c        msq_qg=qg_tbqndkold(1,2,7,i3,i4,p,n)
c        msq_qbarg=qg_tbqndkold(7,2,1,i3,i4,p,n)
        msq_qg=qg_tbqndk(1,2,7,i3,i4,p,n)
        msq_qbarg=qg_tbqndk(7,2,1,i3,i4,p,n)
      else
        write(6,*) 'Invalid value in qg_tbqdk_gvec.f: in=',in
        stop
      endif

c--- fill matrix elements
      do j=1,4
        msq(+j,0)=fac*Vsum(+j)*msq_qg
        msq(-j,0)=fac*Vsum(-j)*msq_qbarg
        msq(0,+j)=fac*Vsum(+j)*msq_gq
        msq(0,-j)=fac*Vsum(-j)*msq_gqbar
      enddo

      return
      end

      function qg_tbqndk(p1,p2,p7,p3,p4,p,n)
      implicit none
      include 'types.f'
      real(dp):: qg_tbqndk

C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> W(l(p3)+a(p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,p7,nu,e6,k6,ie6
      real(dp):: p(mxpart,4),n(4),s345,
     &                 c6,q(mxpart,4),dot,qg_tbqndk_amp,
     &                 twonDpt,twonDp6,qg_tbqndk_ampanti
      complex(dp):: zanb(mxpart,mxpart),zbna(mxpart,mxpart)
      parameter(e6=5,k6=6)

      call checkndotp(p,n,p2)

c--- choice of arbitrary vector to make p6 massless
      ie6=2

      c6=mb**2/(2._dp*dot(p,6,ie6))

      do nu=1,4
      q(1,nu)=p(p1,nu)
      q(2,nu)=p(p2,nu)
      q(3,nu)=p(3,nu)
      q(4,nu)=p(4,nu)
      q(7,nu)=p(p7,nu)
      q(e6,nu)=p(ie6,nu)
      q(k6,nu)=p(6,nu)-c6*q(e6,nu)
      enddo

      call spinoru(7,q,za,zb)
      call spinork(7,q,zanb,zbna,n)

      twonDpt=2._dp*(n(4)*(p(3,4)+p(4,4)+p(5,4))
     &            -n(1)*(p(3,1)+p(4,1)+p(5,1))
     &            -n(2)*(p(3,2)+p(4,2)+p(5,2))
     &            -n(3)*(p(3,3)+p(4,3)+p(5,3)))
      twonDp6=2._dp*(n(4)*p(6,4)-n(1)*p(6,1)-n(2)*p(6,2)-n(3)*p(6,3))

      if (p3 == 3) then  ! top-antibottom
       qg_tbqndk=qg_tbqndk_amp(1,2,3,4,7,k6,e6,twonDpt,twonDp6,zanb)
      else                 ! antitop-bottom
       qg_tbqndk=qg_tbqndk_ampanti(1,2,3,4,7,k6,e6,twonDpt,twonDp6,zanb)
      endif

      s345=(p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     &    -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2

c--- add missing overall factors
      qg_tbqndk=qg_tbqndk
     &         /((2._dp*dot(p,p1,p7)-wmass**2)**2)
     &         /((2._dp*dot(p,p3,p4)-wmass**2)**2+(wmass*wwidth)**2)
     &         /((s345-mt**2)**2+(mt*twidth)**2)
     &         *2._dp*dot(p,p3,5)

      return
      end





c--- This is the old version, calling a different form of the amplitude
c--- routines in which the vector n is made massless with respect to one
c--- of the other momenta. As a result it is numerically unstable.

c--- this routine is adapted from qg_tbqn.f; it is extended to
c--- pass the momenta labels of the leptonic current attached
c--- to the W and does not include the width in the W propagator
c--- as it expects to be called with the W in the t-channel.

c      function qg_tbqndkold(p1,p2,p7,p3,p4,p,n)
c      implicit none
c      include 'types.f'
c      real(dp):: qg_tbqndkold
c
cC---calculates the amplitude squared for the process
cc   q(p1)+qbar(p2) --> W(l(p3)+a(p4)+g(p5)
cc   contracted with the vector n(mu)
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'sprods_com.f'
c      include 'zprods_com.f'
c      integer:: p1,p2,p3,p4,p7,nu,e6,k6,en,nh,ie6,ien,j,k
c      real(dp):: p(mxpart,4),n(4),nDn,prop,s345,
c     &                 c6,cn,q(mxpart,4),dot,qg_tbqndk_ampold,xnorm
c      common/xnorm/xnorm
c      parameter(e6=5,k6=6,en=8,nh=9)
c
c      call checkndotp(p,n,p2)
c
cc--- choice of arbitrary vectors to make p6 and n massless
c      ie6=2
c      ien=3
c
c      c6=mb**2/(2._dp*dot(p,6,ie6))
c      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2
c      cn=nDn/(2._dp*(n(4)*p(ien,4)-n(1)*p(ien,1)
c     &            -n(2)*p(ien,2)-n(3)*p(ien,3)))
cc      xnorm=p(ien,4)
c      xnorm=1._dp
c
c      do nu=1,4
c      q(1,nu)=p(p1,nu)
c      q(2,nu)=p(p2,nu)
c      q(3,nu)=p(p3,nu)
c      q(4,nu)=p(p4,nu)
c      q(7,nu)=p(p7,nu)
c      q(e6,nu)=p(ie6,nu)
c      q(k6,nu)=p(6,nu)-c6*q(e6,nu)
c      q(en,nu)=p(ien,nu)
c      q(nh,nu)=(n(nu)-cn*q(en,nu))/xnorm
c      enddo
c
c      call spinoru(9,q,za,zb)
cc
c      qg_tbqndkold=qg_tbqndk_ampold(1,2,3,4,7,k6,e6,nh,en,nDn)
c
c      s345=(p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
c     &    -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2
c
cc--- add missing overall factors
c      qg_tbqndkold=qg_tbqndkold
c     &         /((2._dp*dot(p,p1,p7)-wmass**2)**2)
c     &         /((2._dp*dot(p,p3,p4)-wmass**2)**2+(wmass*wwidth)**2)
c     &         /((s345-mt**2)**2+(mt*twidth)**2)
c     &         *2._dp*dot(p,p3,5)
c
c      return
c      end
c
