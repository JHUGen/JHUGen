      subroutine qqb_wh1jet_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     July, 1999.                                                      *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W + H(p5,p6) + parton(p7) + parton(p8)
c                           |
c                            -->l(p3)+a(p4)


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'ckm.f'
      include 'hdecaymode.f'
      include 'flags.f'
      integer:: j,k,nd
      integer:: n5,n6
c --- remember: nd will count the dipoles
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq15_6(-nf:nf,-nf:nf),msq26_5(-nf:nf,-nf:nf),
     & msq16_5(-nf:nf,-nf:nf),msq25_6(-nf:nf,-nf:nf),
     & msq56_1v(-nf:nf,-nf:nf),msq56_2v(-nf:nf,-nf:nf),
     & msq26_5v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & msq15_6v(-nf:nf,-nf:nf),msq16_2v(-nf:nf,-nf:nf),
     & msq16_5v(-nf:nf,-nf:nf),msq25_6v(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub16_2(4),sub26_1(4),
     & sub15_6(4),sub16_5(4),sub25_6(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub26_1v,sub16_5v,sub16_2v,sub15_2v,sub15_6v,sub25_6v,
     & sub25_1v
      external qqb_WH1jet,qqb_WH1jet_gvec

      ndmax=6

      if (hdecaymode == 'wpwm') then
        n5=9
        n6=10
      else
        n5=7
        n6=8
      endif

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,n5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
      call dips(2,p,2,n5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
      call dips(3,p,1,n6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
      call dips(4,p,2,n6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_WH1jet,qqb_WH1jet_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,n5,n6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,n5,n6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
      call dips(5,p,1,n6,n5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & qqb_WH1jet,qqb_WH1jet_gvec)

      call dips(6,p,2,n6,n5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
      call dips(6,p,n5,n6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     & qqb_WH1jet,qqb_WH1jet_gvec)
      call dips(6,p,2,n5,n6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     & qqb_WH1jet,qqb_WH1jet_gvec)

      msq(:,:,:)=zip

      if (Gflag) then
      do j=-nf,nf
      do k=-nf,nf
c--- do only q-qb and qb-q cases
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
      msq(1,j,k)=-msq15_2(j,k)*sub15_2(qq)/xn
      msq(2,j,k)=-msq25_1(j,k)*sub25_1(qq)/xn
      msq(3,j,k)=-msq16_2(j,k)*sub16_2(qq)/xn
      msq(4,j,k)=-msq26_1(j,k)*sub26_1(qq)/xn
      msq(5,j,k)=xn*(
     &  +msq15_6(j,k)*(sub15_6(qq)+0.5_dp*sub56_1(gg))
     &  +0.5_dp*msq56_1v(j,k)*sub56_1v
     &  +msq16_5(j,k)*(sub16_5(qq)+0.5_dp*sub56_1(gg))
     &  +0.5_dp*msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=xn*(
     &  (msq26_5(j,k)*(sub26_5(qq)+0.5_dp*sub56_2(gg))
     &   +0.5_dp*msq56_2v(j,k)*sub56_2v)
     & +(msq25_6(j,k)*(sub25_6(qq)+0.5_dp*sub56_2(gg))
     &   +0.5_dp*msq56_2v(j,k)*sub56_2v))

c--- note statistical factor of one half for two gluons in the final state
      do nd=1,ndmax
        msq(nd,j,k)=half*msq(nd,j,k)
      enddo

      elseif ((k == 0).and. (j .ne. 0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2._dp*tr*(msq25_1(j,-5)+msq25_1(j,-4)+msq25_1(j,-3)
     &                  +msq25_1(j,-2)+msq25_1(j,-1)+msq25_1(j,+1)
     &                  +msq25_1(j,+2)+msq25_1(j,+3)+msq25_1(j,+4)
     &                  +msq25_1(j,+5))*sub25_1(qg)
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-msq16_5(j,k)*(sub16_5(qq)+sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     &              +msq26_5(j,k)*sub56_2(qq))

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2._dp*tr*(msq15_2(-5,k)+msq15_2(-4,k)+msq15_2(-3,k)
     &                  +msq15_2(-2,k)+msq15_2(-1,k)+msq15_2(+1,k)
     &                  +msq15_2(+2,k)+msq15_2(+3,k)+msq15_2(+4,k)
     &                  +msq15_2(+5,k))*sub15_2(qg)
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     &              +msq16_5(j,k)*sub56_1(qq))
      msq(6,j,k)=-msq26_5(j,k)*(sub26_5(qq)+sub56_2(qq))/xn

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case
c--- note g,g = 1,2 and qb=5, q=6 so (15),(25)-->q and (16),(26)-->qb
      msq(1,j,k)=(msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
     &           +msq15_2(+4,k)+msq15_2(+5,k))*sub15_2(qg)*2._dp*tr
      msq(2,j,k)=(msq25_1(k,+1)+msq25_1(k,+2)+msq25_1(k,+3)
     &           +msq25_1(k,+4)+msq25_1(k,+5))*sub25_1(qg)*2._dp*tr
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     &           +msq16_2(-2,k)+msq16_2(-1,k))*sub16_2(qg)*2._dp*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
     &           +msq26_1(k,-2)+msq26_1(k,-1))*sub26_1(qg)*2._dp*tr

      endif


      enddo
      enddo
      endif

      if (Qflag) then
      do j=-nf,nf
      do k=-nf,nf

      if ((j > 0) .and. (k > 0)) then
        if (j .ne. k) then
c--- Q Q - different flavours
        msq(1,j,k)=msq(1,j,k)+(Vsum(k)-half*Vsq(k,-j))*(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
        msq(2,j,k)=msq(2,j,k)+0.5_dp*Vsq(j,-k)*(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
        msq(3,j,k)=msq(3,j,k)+0.5_dp*Vsq(k,-j)*(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
        msq(4,j,k)=msq(4,j,k)+(Vsum(j)-half*Vsq(j,-k))*(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        else
c--- Q Q - same flavours
c--- note: W+2 jet qq piece only includes W radiation off
c---       the 15 line, hence these are the only contributions
c---       (we could include 15,25 subtractions if W+2jet had
c---        radiation off the 26 line also)
        msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        endif

       elseif ((j < 0).and.(k < 0)) then
        if (j .ne. k) then
c--- QBAR QBAR - different flavours
        msq(1,j,k)=msq(1,j,k)+(Vsum(k)-half*Vsq(k,-j))*(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
        msq(2,j,k)=msq(2,j,k)+0.5_dp*Vsq(j,-k)*(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
        msq(3,j,k)=msq(3,j,k)+0.5_dp*Vsq(k,-j)*(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
        msq(4,j,k)=msq(4,j,k)+(Vsum(j)-half*Vsq(j,-k))*(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        else
c--- QBAR QBAR - same flavours
        msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        endif

      elseif ((j > 0).and.(k < 0)) then
c--- Q QBAR
        if (j == -k) then
        msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &    *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &    *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        msq(5,j,k)=msq(5,j,k)+tr*real(nf,dp)
     &    *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
        msq(6,j,k)=msq(6,j,k)+tr*real(nf,dp)
     &    *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
        else
        msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &    *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &    *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        msq(5,j,k)=msq(5,j,k)+tr*real(nf,dp)
     &    *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
        msq(6,j,k)=msq(6,j,k)+tr*real(nf,dp)
     &    *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
        endif
c--QBAR Q
      elseif ((j < 0).and.(k > 0)) then
      if (j == -k) then
      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nf,dp)
     & *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nf,dp)
     & *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
      else
      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nf,dp)
     & *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nf,dp)
     &  *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)

      endif
      endif


      enddo
      enddo
      endif

      return
      end

