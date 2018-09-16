      subroutine qqb_ZH1jet_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     September, 1999.                                                 *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
c     q(-p1)+qbar(-p2) -->  Z + H(p5,p6) + parton(p7) + parton(p8)     *
c                           |                                          *
c                            -->l(p3)+a(p4)                            *
************************************************************************


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'nflav.f'
      include 'qqgg.f'
      include 'hdecaymode.f'
      integer:: j,k,nd
      integer::n5,n6
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
     & msq25_1v(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub16_2(4),sub26_1(4),
     & sub15_6(4),sub16_5(4),sub25_6(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub25_1v,sub26_1v,sub16_5v,sub16_2v,sub15_2v,sub15_6v,
     & sub25_6v,xninv
      external qqb_ZH1jet,qqb_ZH_gvec
      ndmax=6
      if(hdecaymode.ne.'wpwm') then
         n5=7
         n6=8
      else
         n5=9
         n6=10
      endif

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,n5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_ZH1jet,qqb_ZH_gvec)
      call dips(2,p,2,n5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_ZH1jet,qqb_ZH_gvec)
      call dips(3,p,1,n6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_ZH1jet,qqb_ZH_gvec)
      call dips(4,p,2,n6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_ZH1jet,qqb_ZH_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,n5,n6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & qqb_ZH1jet,qqb_ZH_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,n5,n6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     & qqb_ZH1jet,qqb_ZH_gvec)
      call dips(5,p,1,n6,n5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & qqb_ZH1jet,qqb_ZH_gvec)

      call dips(6,p,2,n6,n5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & qqb_ZH1jet,qqb_ZH_gvec)
      call dips(6,p,n5,n6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     & qqb_ZH1jet,qqb_ZH_gvec)
      call dips(6,p,2,n5,n6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     & qqb_ZH1jet,qqb_ZH_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=zero
      enddo
      enddo
      enddo

      xninv=one/xn

c--- Gflag subtraction pieces
      do j=-nf,nf
      do k=-nf,nf

      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

c--- do only q-qb and qb-q cases
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
C-----half=statistical factor
      msq(1,j,k)=-half*msq15_2(j,k)*sub15_2(qq)*xninv
      msq(2,j,k)=-half*msq25_1(j,k)*sub25_1(qq)*xninv
      msq(3,j,k)=-half*msq16_2(j,k)*sub16_2(qq)*xninv
      msq(4,j,k)=-half*msq26_1(j,k)*sub26_1(qq)*xninv
      msq(5,j,k)=half*xn*(
     &  msq15_6(j,k)*(sub15_6(qq)+half*sub56_1(gg))
     & +half*msq56_1v(j,k)*sub56_1v
     & +msq16_5(j,k)*(sub16_5(qq)+half*sub56_1(gg))
     & +half*msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=half*xn*(
     &  msq26_5(j,k)*(sub26_5(qq)+half*sub56_2(gg))
     & +half*msq56_2v(j,k)*sub56_2v
     & +msq25_6(j,k)*(sub25_6(qq)+half*sub56_2(gg))
     & +half*msq56_2v(j,k)*sub56_2v)
      elseif ((k == 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=two*tr*msq25_1(j,-j)*sub25_1(qg)
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-(msq16_5(j,k)*sub16_5(qq)+msq16_5(j,k)*sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     &              +msq26_5(j,k)*sub56_2(qq))

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=two*tr*msq15_2(-k,k)*sub15_2(qg)
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     &              +msq15_6(j,k)*sub56_1(qq))
      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> qb(p5)+q(p6)
c---Hence 15 split multiplies q(15)+g(p2)-->Z+q(p6)
c---Hence 25 split multiplies g(p1)+q(p25)-->Z+q(p6)
      msq(1,j,k)=(msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
     &           +msq15_2(+4,k)+msq15_2(+5,k))*sub15_2(qg)*two*tr
      msq(2,j,k)=(msq25_1(k,+1)+msq25_1(k,+2)+msq25_1(k,+3)
     &           +msq25_1(k,+4)+msq25_1(k,+5))*sub25_1(qg)*two*tr
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     &           +msq16_2(-2,k)+msq16_2(-1,k))*sub16_2(qg)*two*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
     &           +msq26_1(k,-2)+msq26_1(k,-1))*sub26_1(qg)*two*tr

      endif

 19   continue
      enddo
      enddo

c--- Qflag subtraction pieces
      do j=-nf,nf
      do k=-nf,nf

      if (((j > 0).and.(k > 0)) .or.
     &    ((j < 0).and.(k < 0))) then
c--q-q or qb-qb
      if (j==k) then
      msq(1,j,k)=msq(1,j,k)+half*(xn-xninv)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(2,j,k)=msq(2,j,k)+half*(xn-xninv)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+half*(xn-xninv)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=msq(4,j,k)+half*(xn-xninv)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      else
      msq(1,j,k)=msq(1,j,k)+(xn-xninv)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-xninv)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      endif
      elseif ((j > 0).and.(k < 0)) then
c q-qbar
      if (j==-k) then
      msq(1,j,k)=msq(1,j,k)+(xn-xninv)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-xninv)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nflav,dp)
     & *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nflav,dp)
     & *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
      else
      msq(1,j,k)=msq(1,j,k)+(xn-xninv)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-xninv)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      endif
c--qbar-q
      elseif ((j < 0).and.(k > 0)) then
      if (j==-k) then
      msq(2,j,k)=msq(2,j,k)+(xn-xninv)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-xninv)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nflav,dp)
     & *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nflav,dp)
     & *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
      else
      msq(2,j,k)=msq(2,j,k)+(xn-xninv)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-xninv)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)

      endif
      endif


      enddo
      enddo

      return
      end

