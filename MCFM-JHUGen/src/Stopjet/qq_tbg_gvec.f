      subroutine qq_tbg_gvec(p,n,in,msq)
************************************************************************
*     Lowest order s-channel single top + jet                          *
*       (contracted with vector n)                                     *
*                                                                      *
*     q(p1) + q(p2) -> t(p3) + b(p4) + g(p5)                           *      
*                                                                      *
*         Author: J. Campbell, June 24, 2008                           *
*                                                                      *
************************************************************************
*                                                                      *
*     IMPORTANT NOTE!                                                  *
*                                                                      *
*     For now, we only include radiation from heavy quark line         *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'nwz.f'
      include 'stopscales.f'
      double precision p(mxpart,4),fac,msq_qqbar
      double precision msq(-nf:nf,-nf:nf),gsq_H
      double precision qg_tbqn,n(4)
      integer j,k,in,i3,i4

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      
      gsq_H=fourpi*as_H
      fac=aveqq*2d0*xn**2*Cf*gsq_H*gwsq**2
      call dotem(5,p,s)

c--- set labels of quark and antiquark according to nwz
      if (nwz .eq. +1) then
        i3=3
        i4=4
      else
        i3=4
        i4=3
      endif

      msq_qqbar=0d0
c      msq_qbarq=0d0
      if     (in .eq. 5) then
        msq_qqbar=+qg_tbqn(i3,i4,1,2,5,p,n)
c        msq_qbarq=+qg_tbqn(i3,i4,2,1,5,p,n)
      else
        write(6,*) 'Value of in not allowed in qq_tbg_gvec: in=',in
        stop
      endif

c--- put result in msq(0,0) element for convenience 
      msq(0,0)=fac*msq_qqbar

c--- fill matrix elements
c      do j=-4,4
c      do k=-4,4
c        if     ((j .gt. 0) .and. (k .lt. 0)) then
c          msq(j,k)=fac*Vsq(j,k)*msq_qqbar
c        elseif ((j .lt. 0) .and. (k .gt. 0)) then
c          msq(j,k)=fac*Vsq(j,k)*msq_qbarq
c        endif
c      enddo
c      enddo

      return
      end

