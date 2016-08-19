      subroutine qq_tbg_g(p,msq)
************************************************************************
*     Real MEs for s-channel single top + jet                          *
*                                                                      *
*     q(p1) + q(p2) -> t(p3) + b(p4) + g(p5) + g(p6)                   *
*                                                                      *
*      (and related crossings and MEs)                                 *
*                                                                      *
*     Author: J. Campbell, June 23, 2008                               *
*                                                                      *
************************************************************************
*                                                                      *
*     IMPORTANT NOTE!                                                  *
*                                                                      *
*     For now, we only include radiation from heavy quark line         *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'nflav.f'
      include 'stopscales.f'
      double precision p(mxpart,4),pp(mxpart,4),
     . msq_qqb_gg_lc,msq_qqb_gg_slc,msq_qqb_qq
      double precision msq(-nf:nf,-nf:nf)
      integer j,k

c--- initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      msqLH(j,k)=0d0
      msqHL(j,k)=0d0
      enddo
      enddo

c--- fill momentum array to call matrix elements and also      
c--- fill array with 1 and 2 switched, to calculate qqbar matrix element
      do j=1,6
      do k=1,4
      pp(j,k)=p(j,k)
c      if (j .lt. 3) then
c        qq(j,k)=p(3-j,k)
c      else
c        qq(j,k)=p(j,k)
c      endif
      enddo
      enddo
      
c--- calculate the qqbar and qbarq matrix elements with an extra gluon
c--- Note: qbq matrix elements are later excluded by PDF choice anyway
      call inters(pp,msq_qqb_gg_lc,msq_qqb_gg_slc)
c      call inters(qq,msq_qbq_gg)
      
c--- the matrix elements with an extra quark line
      call inters_qq(pp,msq_qqb_qq)
c      call inters_qq(qq,msq_qbq_qq)

c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)
      msq(0,0)=msq_qqb_gg_lc
      msq(0,1)=msq_qqb_gg_slc
      msq(1,0)=2d0*Tr*dfloat(nflav)*msq_qqb_qq
           
c--- note: all averaging factors are already included in "inter"
c      do j=-4,4
c      do k=-4,4
c        if     ((j .gt. 0) .and. (k .lt. 0)) then
c          msq(j,k)=msq_qqb_gg+dfloat(nflav)*msq_qqb_qq
c        elseif ((j .lt. 0) .and. (k .gt. 0)) then
cc          msq(j,k)=Vsq(j,k)*(msq_qbq_gg+dfloat(nflav)*msq_qbq_qq)
c        endif
c      enddo
c      enddo
      
      return
      end

