      subroutine qg_tbqdk_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Real MEs for t-channel single top, with explicit b-quark         *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5) + g(p6)                  *
*                                                                      *
*      (and related crossings and MEs)                                 *
*                                                                      *
*     Author: J. Campbell, March 18, 2008                              *
*                         (added decay May 2011)                       *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ckm.f'
      include 'noglue.f'
      include 'nwz.f'
      include 'stopscales.f'
      real(dp):: p(mxpart,4),msq_qbarg,msq_gqbar,
     & msq_qg,msq_gq,msq_gg_a,msq_gg_b,msq_gg_intf,
     & msq_qq12(4),msq_qq12_intf(4),
     & msq_qq56(4),msq_qq56_intf(4)
      real(dp):: msq(-nf:nf,-nf:nf)
      integer:: j,k,i

c--- initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      msqLH(j,k)=0._dp
      msqHL(j,k)=0._dp
      enddo
      enddo

c--- implement shortcuts for checking
      if (noglue) then
        msq_qg=0._dp
        msq_qbarg=0._dp
        msq_gq=0._dp
        msq_gqbar=0._dp
      else
         if (ggonly) then
            msq_qg=0._dp
            msq_qbarg=0._dp
            msq_gq=0._dp
            msq_gqbar=0._dp
         else
c set mass of quark and antiquark according to nwz
            if (nwz == +1) then
c---  calculate the qg and qbarg matrix element first
               call interdk(p,1,2,3,4,7,8,mt,mb,msq_qg)
               call interdk(p,7,2,3,4,1,8,mt,mb,msq_qbarg)
c---  now calculate the gq and gqbar matrix element
               call interdk(p,2,1,3,4,7,8,mt,mb,msq_gq)
               call interdk(p,7,1,3,4,2,8,mt,mb,msq_gqbar)
            else
c---  calculate the qg and qbarg matrix element first
               call interdk(p,1,2,4,3,7,8,mb,mt,msq_qg)
               call interdk(p,7,2,4,3,1,8,mb,mt,msq_qbarg)
c---  now calculate the gq and gqbar matrix element
               call interdk(p,2,1,4,3,7,8,mb,mt,msq_gq)
               call interdk(p,7,1,4,3,2,8,mb,mt,msq_gqbar)
            endif
         endif
c--- the gg matrix element
         if (gqonly) then
            msq_gg_a=0._dp
            msq_gg_b=0._dp
            msq_gg_intf=0._dp
         else
c--- msq_gg_a and msq_gg_b are the two contributing gluon
c--- configurations. msq_gg_intf is the interference term
c--- between the two configurations
            call interdk_gg(p,msq_gg_a,msq_gg_b,msq_gg_intf)
         endif
      endif

      if ((ggonly) .or. (gqonly)) then
         do i=1,4
            msq_qq12(i)=0._dp
            msq_qq12_intf(i)=0._dp
            msq_qq56(i)=0._dp
            msq_qq56_intf(i)=0._dp
         enddo
      else
c--- the matrix elements with an extra quark line:
         if (nwz == +1) then
            call interdk_qq(p,1,2,3,4,7,8,mt,mb,msq_qq12,msq_qq12_intf)
            call interdk_qq(p,7,8,3,4,1,2,mt,mb,msq_qq56,msq_qq56_intf)
         else
            call interdk_qq(p,1,2,4,3,7,8,mb,mt,msq_qq12,msq_qq12_intf)
            call interdk_qq(p,7,8,4,3,1,2,mb,mt,msq_qq56,msq_qq56_intf)
         endif
      endif

c--- note: all averaging factors are already included in "inter"
      do j=-4,4
         do k=-4,4
c--- diagonal contributions
            if     ((j > 0) .and. (k == 0)) then
               msqLH(j,k)=Vsum(j)*msq_qg
            elseif ((j < 0) .and. (k == 0)) then
               msqLH(j,k)=Vsum(j)*msq_qbarg
            elseif ((j == 0) .and. (k > 0)) then
               msqHL(j,k)=Vsum(k)*msq_gq
            elseif ((j == 0) .and. (k < 0)) then
               msqHL(j,k)=Vsum(k)*msq_gqbar

c--- gg contribution: there are two contributing flavours on
c--- the light quark line. We do not know where to add the
c--- initial state interference terms, so add 50% of the interference
c--- to each of the configurations.
            elseif ((j == 0) .and. (k == 0)) then
               msqLH(j,k)=2._dp*(msq_gg_a+msq_gg_intf/2._dp)
               msqHL(j,k)=2._dp*(msq_gg_b+msq_gg_intf/2._dp)
c--- qq/aq/qa/aa contributions
            elseif ((j > 0) .and. (k > 0)) then
               msqLH(j,k)=Vsum(j)*msq_qq12(1)+Vsq(j,-k)*(
     &              msq_qq12(2)-msq_qq12(1)+msq_qq12_intf(1))/2._dp
               msqHL(j,k)=Vsum(k)*msq_qq12(4)+Vsq(-j,k)*(
     &              msq_qq12(3)-msq_qq12(4)+msq_qq12_intf(3))/2._dp
c--- Split initial state interference and add equally to both contributions
               if (j==k) then
                  msqLH(j,k)=msqLH(j,k)+Vsum(j)*msq_qq12_intf(4)/2._dp
                  msqHL(j,k)=msqHL(j,k)+Vsum(k)*msq_qq12_intf(4)/2._dp
               endif
            elseif ((j < 0) .and. (k > 0)) then
               msqLH(j,k)=Vsum(j)*msq_qq56(1)
               msqHL(j,k)=Vsum(k)*msq_qq12(4)
            elseif ((j > 0) .and. (k < 0)) then
               msqLH(j,k)=Vsum(j)*msq_qq12(1)
               msqHL(j,k)=Vsum(k)*msq_qq56(2)
            elseif ((j < 0) .and. (k < 0)) then
               msqLH(j,k)=Vsum(j)*msq_qq56(1)+Vsq(j,-k)*(
     &              msq_qq56(4)-msq_qq56(1)+msq_qq56_intf(4))/2._dp
               msqHL(j,k)=Vsum(k)*msq_qq56(2)+Vsq(-j,k)*(
     &              msq_qq56(3)-msq_qq56(2)+msq_qq56_intf(2))/2._dp
c--- Split initial state interference and add equally to both contributions
               if (j==k) then
                  msqLH(j,k)=msqLH(j,k)+Vsum(j)*msq_qq56_intf(1)/2._dp
                  msqHL(j,k)=msqHL(j,k)+Vsum(k)*msq_qq56_intf(1)/2._dp
               endif
            endif
            msq(j,k)=msqLH(j,k)+msqHL(j,k)
         enddo
      enddo


      return
      end

