      subroutine qq_tbg(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Lowest order s-channel single top + jet                          *
*                                                                      *
*     q(p1) + q(p2) -> t(p3) + b(p4) + g(p5)                           *
*                                                                      *
*     Originally: R. Frederix and F. Tramontano, February 2008         *
*        Adapted: J. Campbell, June 19, 2008                           *
*                                                                      *
************************************************************************
*                                                                      *
*     IMPORTANT NOTE!                                                  *
*                                                                      *
*     For now, we only include radiation from heavy quark line         *
*                                                                      *
************************************************************************
c     u + d  ->  c + s + g  (s-channel single-charm)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      include 'stopscales.f'
      integer:: i3,i4
      real(dp):: p(mxpart,4),fac,msq_qqbar
      real(dp):: dot,msq(-nf:nf,-nf:nf),
     & Wprop12,mq,ma,gsq_H
      complex(dp):: amps_qqbar(2,2,2)
      integer:: hg,hc,hs,j,k

c--- initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      gsq_H=fourpi*as_H
      fac=aveqq*2._dp*xn**2*Cf*gsq_H*gwsq**2

c--- set mass of quark and antiquark according to nwz
      if (nwz == +1) then
        mq=mt
        ma=mb
        i3=3
        i4=4
      else
        mq=mb
        ma=mt
        i3=4
        i4=3
      endif

c--- propagator for qg and qbarg amplitudes
c--- note: in principle we could need a width here, if the top quark
c---       is far off-shell; for simplicity, we do not include it.
      Wprop12=1._dp/(2._dp*dot(p,1,2)-wmass**2)

      call bornwrap(p,1,5,i3,i4,2,mq,ma,amps_qqbar)
c--- Note: qbq matrix elements are later excluded by PDF choice anyway
c      call bornwrap(p,2,5,i3,i4,1,mq,ma,amps_qbarq)

c--- now square up amplitudes and add propagators
      msq_qqbar=0._dp
c      msq_qbarq=0._dp
      do hg=1,2
      do hc=1,2
      do hs=1,2
      msq_qqbar=msq_qqbar+Wprop12**2*abs(amps_qqbar(hg,hc,hs))**2
c      msq_qbarq=msq_qbarq+Wprop12**2*abs(amps_qbarq(hg,hc,hs))**2
      enddo
      enddo
      enddo

c--- put result in msq(0,0) element for LO and real, fill the other
c--- elements too to make sure that the virtual CT's work as well
      msq(0,0)=fac*msq_qqbar
      msq(1,0)=msq(0,0)
      msq(0,1)=msq(0,0)

c--- fill matrix elements
c      do j=-4,4
c      do k=-4,4
c        if     ((j > 0) .and. (k < 0)) then
c          msq(j,k)=fac*Vsq(j,k)*msq_qqbar
c        elseif ((j < 0) .and. (k > 0)) then
c          msq(j,k)=fac*Vsq(j,k)*msq_qbarq
c        endif
c      enddo
c      enddo

      return
      end
