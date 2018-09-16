      subroutine qg_tbq_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Lowest order t-channel single top, with explicit b-quark         *
*       (contracted with vector n)                                     *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *
*                                                                      *
*         Author: J. Campbell, March 19, 2008                          *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'stopscales.f'
      real(dp):: p(mxpart,4),fac,msq_qg,msq_gq,msq_qbarg,msq_gqbar
      real(dp):: msq(-nf:nf,-nf:nf),gsq_H
      real(dp):: qg_tbqn,n(4)
      integer:: j,k,in,i3,i4

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      gsq_H=fourpi*as_H
      fac=aveqg*2._dp*xn**2*Cf*gsq_H*gwsq**2
      call dotem(5,p,s)

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
        msq_gq=-qg_tbqn(i3,i4,2,5,1,p,n)
        msq_gqbar=-qg_tbqn(i3,i4,5,2,1,p,n)
      elseif (in == 2) then
        msq_qg=-qg_tbqn(i3,i4,1,5,2,p,n)
        msq_qbarg=-qg_tbqn(i3,i4,5,1,2,p,n)
      endif

c--- fill matrix elements
      do j=1,4
       msq(+j,0)=fac*Vsum(+j)*msq_qg
       msq(-j,0)=fac*Vsum(-j)*msq_qbarg
       msq(0,+j)=fac*Vsum(+j)*msq_gq
       msq(0,-j)=fac*Vsum(-j)*msq_gqbar
      enddo

c      write(6,*) msq_gq
c      write(6,*) msq_qg
c      write(6,*)

      return
      end

c--- this routine is adapted from qg_tbqn.f; it is extended to
c--- pass the momenta labels of the leptonic current attached
c--- to the W and does not include the width in the W propagator
c--- as it expects to be called with the W in the t-channel.

c--- NB: it also requires an extra factor of two wrt that routine

      function qg_tbqn(p1,p2,p3,p4,p5,p,n)
      implicit none
      include 'types.f'
      real(dp):: qg_tbqn

C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> W(l(p3)+a(p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5
      real(dp):: p(mxpart,4),n(4),nDn,prop,
     &                 nDp1,nDp2,nDp3,nDp4
      nDp1=n(4)*p(p1,4)-n(3)*p(p1,3)-n(2)*p(p1,2)-n(1)*p(p1,1)
      nDp2=n(4)*p(p2,4)-n(3)*p(p2,3)-n(2)*p(p2,2)-n(1)*p(p2,1)
      nDp3=n(4)*p(p3,4)-n(3)*p(p3,3)-n(2)*p(p3,2)-n(1)*p(p3,1)
      nDp4=n(4)*p(p4,4)-n(3)*p(p4,3)-n(2)*p(p4,2)-n(1)*p(p4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,p5)

c---calculate the propagator
      prop=(s(p3,p4)-wmass**2)**2

      qg_tbqn=-nDp3*nDp4+0.25_dp*s(p3,p4)*nDn
     & +s(p1,p4)*nDp2**2*(s(p3,p2)+s(p3,p5))/s(p2,p5)**2
     & +s(p3,p2)*nDp1**2*(s(p1,p4)+s(p4,p5))/s(p1,p5)**2

     & -(nDp2*(nDp3*s(p1,p4)-(s(p3,p2)+s(p3,p5))*nDp4)
     & +s(p3,p5)*0.25_dp*nDn*(s(p4,p2)+s(p1,p4)))/s(p2,p5)

     & -(nDp1*(s(p3,p2)*nDp4-(s(p1,p4)+s(p4,p5))*nDp3)
     & +s(p4,p5)*0.25_dp*nDn*(s(p3,p2)+s(p3,p1)))/s(p1,p5)

     & -(nDp1*nDp2*(s(p3,p2)*s(p4,p5)+s(p3,p5)*s(p1,p4)
     & +s(p3,p5)*s(p4,p5)+2._dp*s(p3,p2)*s(p1,p4))
     & -s(p3,p5)*s(p4,p5)*s(p1,p2)*0.25_dp*nDn)/s(p1,p5)/s(p2,p5)

      qg_tbqn=qg_tbqn/prop

c--- apply extra factor of two
      qg_tbqn=qg_tbqn*2._dp

      return
      end

