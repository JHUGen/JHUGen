      subroutine gg_hgg_gs(p,msq)
      implicit none
************************************************************************
*     Author: J.M. Campbell                                            *
*     February, 2005.                                                  *
*                                                                      *
*     Matrix element SUBTRACTION squared and averaged over             *
*      initial colors and spins                                        *
*                                                                      *
*     f(-p1)+f(-p2) -->  H + parton(p5) + parton(p6) + parton(p7)      *
*                        |                                             *
*                        --> b(p3)+bbar(p4)                            *
*                                                                      *
************************************************************************

      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'msq_struc.f'
      include 'qqgg.f'
      include 'subdefs_gg_hgg.f'
      include 'nflav.f'
      include 'bitflags.f'
      integer::j,k,m,n,nd
      real(dp)::p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),nfactor,
     & sub2q3g_ii,sub2q3g_iffi,sub2q3g_if,sub2q3g_fi,sub2q3g_ff,
     & msq_ab,msq_ba,msq_sym,msqv_ab,msqv_ba,msqv_sym,subgg,subqq,subv
      external gg_hgg,gg_hgg_gvec

c--- statement functions for 2q3g subtraction terms, with labelling
c--- for g(-p1)+g(-p2) -> H(p3+p4)+q(p5)+qbar(p6)+g(p7)

c--- _ii : appropriate for initial-initial dipoles (17_2, 27_1)
      sub2q3g_ii(msq_ab,msq_ba,subgg,msqv_ab,msqv_ba,subv)=
     & +half*xn*((msq_ab+msq_ba)*subgg+(msqv_ab+msqv_ba)*subv)

c--- _if : appropriate for initial-final dipoles (17_5, 17_6, 27_5, 27_5)
      sub2q3g_if(msq_ab,msq_sym,subgg,msqv_ab,msqv_sym,subv)=
     & +half*xn*((msq_ab+msq_sym)*subgg+(msqv_ab+msqv_sym)*subv)

c--- _fi : appropriate for final-initial dipoles (57_1, 67_1, 57_2, 67_2)
      sub2q3g_fi(msq_ab,msq_sym,subqq)=
     & +half*xn*(msq_ab+msq_sym)*subqq

c--- _iffi : appropriate for dipoles which are the sum of
c---         initial-final and final-initial (17_5+57_1, etc.)
      sub2q3g_iffi(msq_ab,msq_sym,subgg,subqq,msqv_ab,msqv_sym,subv)=
     &  sub2q3g_if(msq_ab,msq_sym,subgg,msqv_ab,msqv_sym,subv)
     & +sub2q3g_fi(msq_ab,msq_sym,subqq)

c--- _ff : appropriate for final-final dipoles (57_6, 67_5)
      sub2q3g_ff(msq_ab,msq_ba,msq_sym,subqq)=
     & -half/xn*(msq_ab+msq_ba+msq_sym*(one+xn**2))*subqq

c--- maximum number of dipoles for this process
      ndmax=24

************************************************************************
* CALCULATE ALL THE DIPOLES FOR EIKONALS WITH p7                       *
************************************************************************
c--- calculate all the initial-initial dipoles
      call dips( 1,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc17_2,msq_struc17_2v)
      call dips( 2,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc27_1,msq_struc27_1v)
c--- now the basic initial final ones
      call dips( 7,p,1,7,5,sub17_5,sub17_5v,msq17_5,msq17_5v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc17_5,msq_struc17_5v)
      call dips( 8,p,2,7,6,sub27_6,sub27_6v,msq27_6,msq27_6v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc27_6,msq_struc27_6v)
      call dips( 9,p,1,7,6,sub17_6,sub17_6v,msq17_6,msq17_6v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc17_6,msq_struc17_6v)
      call dips(10,p,2,7,5,sub27_5,sub27_5v,msq27_5,msq27_5v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc27_5,msq_struc27_5v)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips( 7,p,5,7,1,sub57_1,sub57_1v,dummy,msq57_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc57_1,msq_struc57_1v)
      call dips( 8,p,6,7,2,sub67_2,sub67_2v,dummy,msq67_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc67_2,msq_struc67_2v)
      call dips( 9,p,6,7,1,sub67_1,sub67_1v,dummy,msq67_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc67_1,msq_struc67_1v)
      call dips(10,p,5,7,2,sub57_2,sub57_2v,dummy,msq57_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc57_2,msq_struc57_2v)
c--- lastly, the final final terms
      call dips(19,p,5,7,6,sub57_6,sub57_6v,msq57_6,msq57_6v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc57_6,msq_struc57_6v)
      call dips(20,p,6,7,5,sub67_5,sub67_5v,msq67_5,msq67_5v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc67_5,msq_struc67_5v)

************************************************************************
* CALCULATE ALL THE DIPOLES FOR EIKONALS WITH p6                       *
************************************************************************
c--- calculate all the initial-initial dipoles
      call dips( 3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc16_2,msq_struc16_2v)
      call dips( 4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc26_1,msq_struc26_1v)
c--- now the basic initial final ones
      call dips(11,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc16_5,msq_struc16_5v)
      call dips(12,p,2,6,7,sub26_7,sub26_7v,msq26_7,msq26_7v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc26_7,msq_struc26_7v)
      call dips(13,p,1,6,7,sub16_7,sub16_7v,msq16_7,msq16_7v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc16_7,msq_struc16_7v)
      call dips(14,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc26_5,msq_struc26_5v)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(11,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc56_1,msq_struc56_1v)
      call dips(12,p,7,6,2,sub76_2,sub76_2v,dummy,msq76_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc76_2,msq_struc76_2v)
      call dips(13,p,7,6,1,sub76_1,sub76_1v,dummy,msq76_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc76_1,msq_struc76_1v)
      call dips(14,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc56_2,msq_struc56_2v)
c--- lastly, the final final terms
      call dips(21,p,5,6,7,sub56_7,sub56_7v,msq56_7,msq56_7v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc56_7,msq_struc56_7v)
      call dips(22,p,7,6,5,sub76_5,sub76_5v,msq76_5,msq76_5v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc76_5,msq_struc76_5v)

************************************************************************
* CALCULATE ALL THE DIPOLES FOR EIKONALS WITH p5                       *
************************************************************************
c--- calculate all the initial-initial dipoles
      call dips( 5,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc15_2,msq_struc15_2v)
      call dips( 6,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc25_1,msq_struc25_1v)
c--- now the basic initial final ones
      call dips(15,p,1,5,7,sub15_7,sub15_7v,msq15_7,msq15_7v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc15_7,msq_struc15_7v)
      call dips(16,p,2,5,6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc25_6,msq_struc25_6v)
      call dips(17,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc15_6,msq_struc15_6v)
      call dips(18,p,2,5,7,sub25_7,sub25_7v,msq25_7,msq25_7v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc25_7,msq_struc25_7v)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(15,p,7,5,1,sub75_1,sub75_1v,dummy,msq75_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc75_1,msq_struc75_1v)
      call dips(16,p,6,5,2,sub65_2,sub65_2v,dummy,msq65_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc65_2,msq_struc65_2v)
      call dips(17,p,6,5,1,sub65_1,sub65_1v,dummy,msq65_1v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc65_1,msq_struc65_1v)
      call dips(18,p,7,5,2,sub75_2,sub75_2v,dummy,msq75_2v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc75_2,msq_struc75_2v)
c--- lastly, the final final terms
      call dips(23,p,7,5,6,sub75_6,sub75_6v,msq75_6,msq75_6v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc75_6,msq_struc75_6v)
      call dips(24,p,6,5,7,sub65_7,sub65_7v,msq65_7,msq65_7v,
     & gg_hgg,gg_hgg_gvec)
      call extract(msq_struc65_7,msq_struc65_7v)


      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=zip
      enddo
      enddo
      enddo


      do j=-nf,nf
      do k=-nf,nf

************************************************************************
* SUBTRACTIONS FOR Q-Q and QB-QB                                       *
************************************************************************
      if (  ((j > 0).and.(k > 0))
     & .or. ((j < 0).and.(k < 0))) then

c--- note that all subtractions are performed with reference to the
c--- quark-quark matrix elements.

c--- non-identical quark contribution
      if (j .ne. k) then
      m=abs(j)
      n=abs(k)
      msq(1,j,k)=f4q*two/xn*msq17_2(m,n)*sub17_2(qq)
      msq(2,j,k)=f4q*two/xn*msq27_1(m,n)*sub27_1(qq)
      msq(19,j,k)=f4q*two/xn*msq57_6(m,n)*sub57_6(qq)
      msq(20,j,k)=f4q*two/xn*msq67_5(m,n)*sub67_5(qq)
      msq(7,j,k)=-f4q/xn*msq17_5(m,n)*(sub17_5(qq)+sub57_1(qq))
      msq(8,j,k)=-f4q/xn*msq27_6(m,n)*(sub27_6(qq)+sub67_2(qq))
      msq(9,j,k)=f4q*(xn-two/xn)*msq17_6(m,n)*(sub17_6(qq)+sub67_1(qq))
      msq(10,j,k)=f4q*(xn-two/xn)*msq27_5(m,n)*(sub27_5(qq)+sub57_2(qq))

      m=j/abs(j)
      n=k/abs(k)
      msq(4,j,k)=msq(4,j,k)+(aveqq/aveqg)*f4q*(
     .+(+msq_struc26_1(igg_ab,m,0)+msq_struc26_1(igg_ba,m,0)
     &  +msq_struc26_1(igg_sym,m,0))*sub26_1(gq)
     .+(+msq_struc26_1v(igg_ab,m,0)+msq_struc26_1v(igg_ba,m,0)
     &  +msq_struc26_1v(igg_sym,m,0))*sub26_1v)
      msq(5,j,k)=msq(5,j,k)+(aveqq/aveqg)*f4q*(
     .+(+msq_struc15_2(igg_ab,0,n)+msq_struc15_2(igg_ba,0,n)
     &  +msq_struc15_2(igg_sym,0,n))*sub15_2(gq)
     .+(+msq_struc15_2v(igg_ab,0,n)+msq_struc15_2v(igg_ba,0,n)
     &  +msq_struc15_2v(igg_sym,0,n))*sub15_2v)

c--- identical quark contribution
      else
      m=+1
      n=+1
      msq(1,j,k)=f4q*sub17_2(qq)*(
     & +two/xn*msq_struc17_2(iqq_a,m,n)+two/xn*msq_struc17_2(iqq_b,m,n)
     & +(xn+one/xn)*msq_struc17_2(iqq_i,m,n))
      msq(2,j,k)=f4q*sub27_1(qq)*(
     & +two/xn*msq_struc27_1(iqq_a,m,n)+two/xn*msq_struc27_1(iqq_b,m,n)
     & +(xn+one/xn)*msq_struc27_1(iqq_i,m,n))
      msq(19,j,k)=f4q*sub57_6(qq)*(
     & +two/xn*msq_struc57_6(iqq_a,m,n)+two/xn*msq_struc57_6(iqq_b,m,n)
     & +(xn+one/xn)*msq_struc57_6(iqq_i,m,n))
      msq(20,j,k)=f4q*sub67_5(qq)*(
     & +two/xn*msq_struc67_5(iqq_a,m,n)+two/xn*msq_struc67_5(iqq_b,m,n)
     & +(xn+one/xn)*msq_struc67_5(iqq_i,m,n))
      msq(7,j,k)=f4q*(sub17_5(qq)+sub57_1(qq))*(
     & -one/xn*msq_struc17_5(iqq_i,m,n)-one/xn*msq_struc17_5(iqq_a,m,n)
     & +(xn-two/xn)*msq_struc17_5(iqq_b,m,n))
      msq(8,j,k)=f4q*(sub27_6(qq)+sub67_2(qq))*(
     & -one/xn*msq_struc27_6(iqq_i,m,n)-one/xn*msq_struc27_6(iqq_a,m,n)
     & +(xn-two/xn)*msq_struc27_6(iqq_b,m,n))
      msq(9,j,k)=f4q*(sub17_6(qq)+sub67_1(qq))*(
     & -one/xn*msq_struc17_6(iqq_i,m,n)-one/xn*msq_struc17_6(iqq_b,m,n)
     & +(xn-two/xn)*msq_struc17_6(iqq_a,m,n))
      msq(10,j,k)=f4q*(sub27_5(qq)+sub57_2(qq))*(
     & -one/xn*msq_struc27_5(iqq_i,m,n)-one/xn*msq_struc27_5(iqq_b,m,n)
     & +(xn-two/xn)*msq_struc27_5(iqq_a,m,n))

      m=j/abs(j)
      n=k/abs(k)
      msq(4,j,k)=msq(4,j,k)+half*(aveqq/aveqg)*f4q*(
     .+(+msq_struc26_1(igg_ab,m,0)+msq_struc26_1(igg_ba,m,0)
     &  +msq_struc26_1(igg_sym,m,0))*sub26_1(gq)
     .+(+msq_struc26_1v(igg_ab,m,0)+msq_struc26_1v(igg_ba,m,0)
     &  +msq_struc26_1v(igg_sym,m,0))*sub26_1v)
      msq(5,j,k)=msq(5,j,k)+half*(aveqq/aveqg)*f4q*(
     .+(+msq_struc15_2(igg_ab,0,n)+msq_struc15_2(igg_ba,0,n)
     &  +msq_struc15_2(igg_sym,0,n))*sub15_2(gq)
     .+(+msq_struc15_2v(igg_ab,0,n)+msq_struc15_2v(igg_ba,0,n)
     &  +msq_struc15_2v(igg_sym,0,n))*sub15_2v)
       msq(6,j,k)=msq(6,j,k)+half*(aveqq/aveqg)*f4q*(
     .+(+msq_struc25_1(igg_ab,m,0)+msq_struc25_1(igg_ba,m,0)
     &  +msq_struc25_1(igg_sym,m,0))*sub25_1(gq)
     .+(+msq_struc25_1v(igg_ab,m,0)+msq_struc25_1v(igg_ba,m,0)
     &  +msq_struc25_1v(igg_sym,m,0))*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+half*(aveqq/aveqg)*f4q*(
     .+(+msq_struc16_2(igg_ab,0,n)+msq_struc16_2(igg_ba,0,n)
     &  +msq_struc16_2(igg_sym,0,n))*sub16_2(gq)
     .+(+msq_struc16_2v(igg_ab,0,n)+msq_struc16_2v(igg_ba,0,n)
     &  +msq_struc16_2v(igg_sym,0,n))*sub16_2v)

      endif

************************************************************************
* SUBTRACTIONS FOR Q-QB and QB-Q                                       *
************************************************************************
      elseif (  ((j > 0).and.(k < 0))
     &     .or. ((j < 0).and.(k > 0))) then

c--- note that all subtractions are performed with reference to the
c--- quark-antiquark matrix elements.
      m=abs(j)
      n=-abs(k)

c--- if flavours are not the same, 4q scattering contribution only
      if (j .ne. -k) then
      msq(1,j,k)=f4q*(xn-two/xn)*msq17_2(m,n)*sub17_2(qq)
      msq(2,j,k)=f4q*(xn-two/xn)*msq27_1(m,n)*sub27_1(qq)
      msq(19,j,k)=f4q*(xn-two/xn)*msq57_6(m,n)*sub57_6(qq)
      msq(20,j,k)=f4q*(xn-two/xn)*msq67_5(m,n)*sub67_5(qq)
      msq(7,j,k)=-f4q/xn*msq17_5(m,n)*(sub17_5(qq)+sub57_1(qq))
      msq(8,j,k)=-f4q/xn*msq27_6(m,n)*(sub27_6(qq)+sub67_2(qq))
      msq(9,j,k)=f4q*two/xn*msq17_6(m,n)*(sub17_6(qq)+sub67_1(qq))
      msq(10,j,k)=f4q*two/xn*msq27_5(m,n)*(sub27_5(qq)+sub57_2(qq))

      m=j/abs(j)
      n=k/abs(k)
      msq(4,j,k)=msq(4,j,k)+(aveqq/aveqg)*f4q*(
     .+(+msq_struc26_1(igg_ab,m,0)+msq_struc26_1(igg_ba,m,0)
     &  +msq_struc26_1(igg_sym,m,0))*sub26_1(gq)
     .+(+msq_struc26_1v(igg_ab,m,0)+msq_struc26_1v(igg_ba,m,0)
     &  +msq_struc26_1v(igg_sym,m,0))*sub26_1v)
      msq(5,j,k)=msq(5,j,k)+(aveqq/aveqg)*f4q*(
     .+(+msq_struc15_2(igg_ab,0,n)+msq_struc15_2(igg_ba,0,n)
     &  +msq_struc15_2(igg_sym,0,n))*sub15_2(gq)
     .+(+msq_struc15_2v(igg_ab,0,n)+msq_struc15_2v(igg_ba,0,n)
     &  +msq_struc15_2v(igg_sym,0,n))*sub15_2v)

      else
**********************************
* 4-quark identical contribution *
**********************************
      m=+1
      n=-1
      msq(1,j,k)=f4q*sub17_2(qq)*(
     & -one/xn*msq_struc17_2(iqq_i,m,n)-one/xn*msq_struc17_2(iqq_b,m,n)
     & +(xn-two/xn)*msq_struc17_2(iqq_a,m,n))
      msq(2,j,k)=f4q*sub27_1(qq)*(
     & -one/xn*msq_struc27_1(iqq_i,m,n)-one/xn*msq_struc27_1(iqq_b,m,n)
     & +(xn-two/xn)*msq_struc27_1(iqq_a,m,n))
      msq(19,j,k)=f4q*sub57_6(qq)*(
     & -one/xn*msq_struc57_6(iqq_i,m,n)-one/xn*msq_struc57_6(iqq_b,m,n)
     & +(xn-two/xn)*msq_struc57_6(iqq_a,m,n))
      msq(20,j,k)=f4q*sub67_5(qq)*(
     & -one/xn*msq_struc67_5(iqq_i,m,n)-one/xn*msq_struc67_5(iqq_b,m,n)
     & +(xn-two/xn)*msq_struc67_5(iqq_a,m,n))
      msq(7,j,k)=f4q*(sub17_5(qq)+sub57_1(qq))*(
     & -one/xn*msq_struc17_5(iqq_i,m,n)-one/xn*msq_struc17_5(iqq_a,m,n)
     & +(xn-two/xn)*msq_struc17_5(iqq_b,m,n))
      msq(8,j,k)=f4q*(sub27_6(qq)+sub67_2(qq))*(
     & -one/xn*msq_struc27_6(iqq_i,m,n)-one/xn*msq_struc27_6(iqq_a,m,n)
     & +(xn-two/xn)*msq_struc27_6(iqq_b,m,n))
      msq(9,j,k)=f4q*(sub17_6(qq)+sub67_1(qq))*(
     & +two/xn*msq_struc17_6(iqq_a,m,n)+two/xn*msq_struc17_6(iqq_b,m,n)
     & +(xn+one/xn)*msq_struc17_6(iqq_i,m,n))
      msq(10,j,k)=f4q*(sub27_5(qq)+sub57_2(qq))*(
     & +two/xn*msq_struc27_5(iqq_a,m,n)+two/xn*msq_struc27_5(iqq_b,m,n)
     & +(xn+one/xn)*msq_struc27_5(iqq_i,m,n))

      m=j/abs(j)
      n=k/abs(k)
      msq(7,j,k)=msq(7,j,k)+f4q*(xn-two/xn)*(
     & msq_struc17_5(iqr,m,n)*(sub17_5(qq)+sub57_1(qq)))
      msq(8,j,k)=msq(8,j,k)+f4q*(xn-two/xn)*(
     & msq_struc27_6(iqr,m,n)*(sub27_6(qq)+sub67_2(qq)))
      msq(1,j,k)=msq(1,j,k)-f4q/xn*(
     & msq_struc17_2(iqr,m,n)*sub17_2(qq))
      msq(2,j,k)=msq(2,j,k)-f4q/xn*(
     & msq_struc27_1(iqr,m,n)*sub27_1(qq))
      msq(19,j,k)=msq(19,j,k)-f4q/xn*(
     & msq_struc57_6(iqr,m,n)*sub57_6(qq))
      msq(20,j,k)=msq(20,j,k)-f4q/xn*(
     & msq_struc67_5(iqr,m,n)*sub67_5(qq))
      msq(9,j,k)=msq(9,j,k)+f4q*two/xn*(
     & msq_struc17_6(iqr,m,n)*(sub17_6(qq)+sub67_1(qq)))
      msq(10,j,k)=msq(10,j,k)+f4q*two/xn*(
     & msq_struc27_5(iqr,m,n)*(sub57_2(qq)+sub27_5(qq)))

      m=j/abs(j)
      n=k/abs(k)
      msq(4,j,k)=msq(4,j,k)+(aveqq/aveqg)*f4q*(
     .+(+msq_struc26_1(igg_ab,m,0)+msq_struc26_1(igg_ba,m,0)
     &  +msq_struc26_1(igg_sym,m,0))*sub26_1(gq)
     .+(+msq_struc26_1v(igg_ab,m,0)+msq_struc26_1v(igg_ba,m,0)
     &  +msq_struc26_1v(igg_sym,m,0))*sub26_1v)
      msq(5,j,k)=msq(5,j,k)+(aveqq/aveqg)*f4q*(
     .+(+msq_struc15_2(igg_ab,0,n)+msq_struc15_2(igg_ba,0,n)
     &  +msq_struc15_2(igg_sym,0,n))*sub15_2(gq)
     .+(+msq_struc15_2v(igg_ab,0,n)+msq_struc15_2v(igg_ba,0,n)
     &  +msq_struc15_2v(igg_sym,0,n))*sub15_2v)

      msq(21,j,k)=msq(21,j,k)+two*real(nflav,dp)*f4q*(
     .+(msq_struc56_7(igg_ab,m,n)+msq_struc56_7(igg_ba,m,n)
     & +msq_struc56_7(igg_sym,m,n))*sub56_7(gq)
     .-(msq_struc56_7v(igg_ab,m,n)+msq_struc56_7v(igg_ba,m,n)
     & +msq_struc56_7v(igg_sym,m,n))*sub56_7v)

c--- note that all subtractions are performed with reference to the
c--- antiquark-quark matrix elements. We could call the quark-antiquark
c--- LO, also switching _ab to _ba, but the answer would be the same.
c--- Note also the statistical factor of one 6th for three gluons
c--- in the final state
      m=-1
      n=+1

********************************
* Note that in the sum over    *
* permutations below, factors  *
* of 1/2 enter for 56, 57, 67  *
* dipoles, as they enter twice *
********************************

**********************************
* 2-quark contributions for p7   *
**********************************
      msq(19,j,k)=msq(19,j,k)+f2q/three*sub2q3g_ii(
     & msq_struc57_6(igg_ab,m,n),msq_struc57_6(igg_ba,m,n),sub57_6(gg),
     & msq_struc57_6v(igg_ab,m,n),msq_struc57_6v(igg_ba,m,n),sub57_6v)
      msq(20,j,k)=msq(20,j,k)+f2q/three*sub2q3g_ii(
     & msq_struc67_5(igg_ab,m,n),msq_struc67_5(igg_ba,m,n),sub67_5(gg),
     & msq_struc67_5v(igg_ab,m,n),msq_struc67_5v(igg_ba,m,n),sub67_5v)

      msq(7,j,k)=msq(7,j,k)+f2q/three*sub2q3g_if(
     & msq_struc17_5(igg_ab,m,n),msq_struc17_5(igg_sym,m,n),
     & sub57_1(gg),
     & msq_struc57_1v(igg_ab,m,n),msq_struc57_1v(igg_sym,m,n),sub57_1v)
      msq(7,j,k)=msq(7,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc17_5(igg_ab,m,n),msq_struc17_5(igg_sym,m,n),
     & sub17_5(qq))
      msq(10,j,k)=msq(10,j,k)+f2q/three*sub2q3g_if(
     & msq_struc27_5(igg_ba,m,n),msq_struc27_5(igg_sym,m,n),
     & sub57_2(gg),
     & msq_struc57_2v(igg_ba,m,n),msq_struc57_2v(igg_sym,m,n),sub57_2v)
      msq(10,j,k)=msq(10,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc27_5(igg_ba,m,n),msq_struc27_5(igg_sym,m,n),
     & sub27_5(qq))
      msq(9,j,k)=msq(9,j,k)+f2q/three*sub2q3g_if(
     & msq_struc17_6(igg_ba,m,n),msq_struc17_6(igg_sym,m,n),
     & sub67_1(gg),
     & msq_struc67_1v(igg_ba,m,n),msq_struc67_1v(igg_sym,m,n),sub67_1v)
      msq(9,j,k)=msq(9,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc17_6(igg_ba,m,n),msq_struc17_6(igg_sym,m,n),
     & sub17_6(qq))
      msq(8,j,k)=msq(8,j,k)+f2q/three*sub2q3g_if(
     & msq_struc27_6(igg_ab,m,n),msq_struc27_6(igg_sym,m,n),
     & sub67_2(gg),
     & msq_struc67_2v(igg_ab,m,n),msq_struc67_2v(igg_sym,m,n),sub67_2v)
      msq(8,j,k)=msq(8,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc27_6(igg_ab,m,n),msq_struc27_6(igg_sym,m,n),
     & sub27_6(qq))

      msq(1,j,k)=msq(1,j,k)+two*f2q/three*sub2q3g_ff(
     & msq_struc17_2(igg_ab,m,n),msq_struc17_2(igg_ba,m,n),
     & msq_struc17_2(igg_sym,m,n),sub17_2(qq))

      msq(2,j,k)=msq(2,j,k)+two*f2q/three*sub2q3g_ff(
     & msq_struc27_1(igg_ab,m,n),msq_struc27_1(igg_ba,m,n),
     & msq_struc27_1(igg_sym,m,n),sub27_1(qq))

**********************************
* 2-quark contributions for p6   *
**********************************
      msq(21,j,k)=msq(21,j,k)+f2q/three*sub2q3g_ii(
     & msq_struc56_7(igg_ab,m,n),msq_struc56_7(igg_ba,m,n),sub56_7(gg),
     & msq_struc56_7v(igg_ab,m,n),msq_struc56_7v(igg_ba,m,n),sub56_7v)
      msq(22,j,k)=msq(22,j,k)+f2q/three*sub2q3g_ii(
     & msq_struc76_5(igg_ab,m,n),msq_struc76_5(igg_ba,m,n),sub76_5(gg),
     & msq_struc76_5v(igg_ab,m,n),msq_struc76_5v(igg_ba,m,n),sub76_5v)

      msq(11,j,k)=msq(11,j,k)+f2q/three*sub2q3g_if(
     & msq_struc16_5(igg_ab,m,n),msq_struc16_5(igg_sym,m,n),
     & sub56_1(gg),
     & msq_struc56_1v(igg_ab,m,n),msq_struc56_1v(igg_sym,m,n),sub56_1v)
      msq(11,j,k)=msq(11,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc16_5(igg_ab,m,n),msq_struc16_5(igg_sym,m,n),
     & sub16_5(qq))
      msq(14,j,k)=msq(14,j,k)+f2q/three*sub2q3g_if(
     & msq_struc26_5(igg_ba,m,n),msq_struc26_5(igg_sym,m,n),
     & sub56_2(gg),
     & msq_struc56_2v(igg_ba,m,n),msq_struc56_2v(igg_sym,m,n),sub56_2v)
      msq(14,j,k)=msq(14,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc26_5(igg_ba,m,n),msq_struc26_5(igg_sym,m,n),
     & sub26_5(qq))
      msq(13,j,k)=msq(13,j,k)+f2q/three*sub2q3g_if(
     & msq_struc16_7(igg_ba,m,n),msq_struc16_7(igg_sym,m,n),
     & sub76_1(gg),
     & msq_struc76_1v(igg_ba,m,n),msq_struc76_1v(igg_sym,m,n),sub76_1v)
      msq(13,j,k)=msq(13,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc16_7(igg_ba,m,n),msq_struc16_7(igg_sym,m,n),
     & sub16_7(qq))
      msq(12,j,k)=msq(12,j,k)+f2q/three*sub2q3g_if(
     & msq_struc26_7(igg_ab,m,n),msq_struc26_7(igg_sym,m,n),
     & sub76_2(gg),
     & msq_struc76_2v(igg_ab,m,n),msq_struc76_2v(igg_sym,m,n),sub76_2v)
      msq(12,j,k)=msq(12,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc26_7(igg_ab,m,n),msq_struc26_7(igg_sym,m,n),
     & sub26_7(qq))

      msq(3,j,k)=msq(3,j,k)+two*f2q/three*sub2q3g_ff(
     & msq_struc16_2(igg_ab,m,n),msq_struc16_2(igg_ba,m,n),
     & msq_struc16_2(igg_sym,m,n),sub16_2(qq))

      msq(4,j,k)=msq(4,j,k)+two*f2q/three*sub2q3g_ff(
     & msq_struc26_1(igg_ab,m,n),msq_struc26_1(igg_ba,m,n),
     & msq_struc26_1(igg_sym,m,n),sub26_1(qq))

**********************************
* 2-quark contributions for p5   *
**********************************
      msq(23,j,k)=msq(23,j,k)+f2q/three*sub2q3g_ii(
     & msq_struc75_6(igg_ab,m,n),msq_struc75_6(igg_ba,m,n),sub75_6(gg),
     & msq_struc75_6v(igg_ab,m,n),msq_struc75_6v(igg_ba,m,n),sub75_6v)
      msq(24,j,k)=msq(24,j,k)+f2q/three*sub2q3g_ii(
     & msq_struc65_7(igg_ab,m,n),msq_struc65_7(igg_ba,m,n),sub65_7(gg),
     & msq_struc65_7v(igg_ab,m,n),msq_struc65_7v(igg_ba,m,n),sub65_7v)

      msq(15,j,k)=msq(15,j,k)+f2q/three*sub2q3g_if(
     & msq_struc15_7(igg_ab,m,n),msq_struc15_7(igg_sym,m,n),
     & sub75_1(gg),
     & msq_struc75_1v(igg_ab,m,n),msq_struc75_1v(igg_sym,m,n),sub75_1v)
      msq(15,j,k)=msq(15,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc15_7(igg_ab,m,n),msq_struc15_7(igg_sym,m,n),
     & sub15_7(qq))
      msq(18,j,k)=msq(18,j,k)+f2q/three*sub2q3g_if(
     & msq_struc25_7(igg_ba,m,n),msq_struc25_7(igg_sym,m,n),
     & sub75_2(gg),
     & msq_struc75_2v(igg_ba,m,n),msq_struc75_2v(igg_sym,m,n),sub75_2v)
      msq(18,j,k)=msq(18,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc25_7(igg_ba,m,n),msq_struc25_7(igg_sym,m,n),
     & sub25_7(qq))
      msq(17,j,k)=msq(17,j,k)+f2q/three*sub2q3g_if(
     & msq_struc15_6(igg_ba,m,n),msq_struc15_6(igg_sym,m,n),
     & sub65_1(gg),
     & msq_struc65_1v(igg_ba,m,n),msq_struc65_1v(igg_sym,m,n),sub65_1v)
      msq(17,j,k)=msq(17,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc15_6(igg_ba,m,n),msq_struc15_6(igg_sym,m,n),
     & sub15_6(qq))
      msq(16,j,k)=msq(16,j,k)+f2q/three*sub2q3g_if(
     & msq_struc25_6(igg_ab,m,n),msq_struc25_6(igg_sym,m,n),
     & sub65_2(gg),
     & msq_struc65_2v(igg_ab,m,n),msq_struc65_2v(igg_sym,m,n),sub65_2v)
      msq(16,j,k)=msq(16,j,k)+two*f2q/three*sub2q3g_fi(
     & msq_struc25_6(igg_ab,m,n),msq_struc25_6(igg_sym,m,n),
     & sub25_6(qq))

      msq(5,j,k)=msq(5,j,k)+two*f2q/three*sub2q3g_ff(
     & msq_struc15_2(igg_ab,m,n),msq_struc15_2(igg_ba,m,n),
     & msq_struc15_2(igg_sym,m,n),sub15_2(qq))

      msq(6,j,k)=msq(6,j,k)+two*f2q/three*sub2q3g_ff(
     & msq_struc25_1(igg_ab,m,n),msq_struc25_1(igg_ba,m,n),
     & msq_struc25_1(igg_sym,m,n),sub25_1(qq))

      endif
c--- end of same flavour case

************************************************************************
* SUBTRACTIONS FOR Q-G and QB-G                                        *
************************************************************************
      elseif ((k == 0).and. (j .ne. 0)) then

c--- note that all subtractions are performed with reference to the
c--- quark-gluon matrix elements. We could call the antiquark-gluon
c--- LO, also switching _ab to _ba, but the answer would be the same
      m=+1
      n=0

********************************
* Note that in the sum over    *
* permutations below, factors  *
* of one half enter for 67     *
* dipoles, as they enter twice *
********************************

**********************************
* 2-quark contributions for p7   *
**********************************
      msq(8,j,k)=msq(8,j,k)+f2q*sub2q3g_ii(
     & msq_struc27_6(igg_ba,m,n),msq_struc27_6(igg_ab,m,n),sub27_6(gg),
     & msq_struc27_6v(igg_ba,m,n),msq_struc27_6v(igg_ab,m,n),sub27_6v)
      msq(8,j,k)=msq(8,j,k)+f2q*half*sub2q3g_ii(
     & msq_struc27_6(igg_ba,m,n),msq_struc27_6(igg_ab,m,n),sub67_2(gg),
     & msq_struc67_2v(igg_ba,m,n),msq_struc67_2v(igg_ab,m,n),sub67_2v)

      msq(10,j,k)=msq(10,j,k)+f2q*sub2q3g_if(
     & msq_struc27_5(igg_ba,m,n),msq_struc27_5(igg_sym,m,n),
     & sub27_5(gg),
     & msq_struc27_5v(igg_ba,m,n),msq_struc27_5v(igg_sym,m,n),sub27_5v)
      msq(10,j,k)=msq(10,j,k)+f2q*sub2q3g_fi(
     & msq_struc27_5(igg_ba,m,n),msq_struc27_5(igg_sym,m,n),
     & sub57_2(qq))
      msq(2,j,k)=msq(2,j,k)+f2q*sub2q3g_if(
     & msq_struc27_1(igg_ab,m,n),msq_struc27_1(igg_sym,m,n),
     & sub27_1(gg),
     & msq_struc27_1v(igg_ab,m,n),msq_struc27_1v(igg_sym,m,n),sub27_1v)
      msq(1,j,k)=msq(1,j,k)+f2q*sub2q3g_fi(
     & msq_struc17_2(igg_ab,m,n),msq_struc17_2(igg_sym,m,n),
     & sub17_2(qq))
      msq(20,j,k)=msq(20,j,k)+f2q*half*sub2q3g_if(
     & msq_struc67_5(igg_ab,m,n),msq_struc67_5(igg_sym,m,n),
     & sub67_5(gg),
     & msq_struc67_5v(igg_ab,m,n),msq_struc67_5v(igg_sym,m,n),sub67_5v)
      msq(19,j,k)=msq(19,j,k)+f2q*sub2q3g_fi(
     & msq_struc57_6(igg_ab,m,n),msq_struc57_6(igg_sym,m,n),
     & sub57_6(qq))
      msq(9,j,k)=msq(9,j,k)+f2q*half*sub2q3g_if(
     & msq_struc17_6(igg_ba,m,n),msq_struc17_6(igg_sym,m,n),
     & sub67_1(gg),
     & msq_struc67_1v(igg_ba,m,n),msq_struc67_1v(igg_sym,m,n),sub67_1v)
      msq(9,j,k)=msq(9,j,k)+f2q*sub2q3g_fi(
     & msq_struc17_6(igg_ba,m,n),msq_struc17_6(igg_sym,m,n),
     & sub17_6(qq))

      msq(7,j,k)=msq(7,j,k)+f2q*sub2q3g_ff(
     & msq_struc17_5(igg_ba,m,n),msq_struc17_5(igg_ab,m,n),
     & msq_struc17_5(igg_sym,m,n),sub57_1(qq))
      msq(7,j,k)=msq(7,j,k)+f2q*sub2q3g_ff(
     & msq_struc17_5(igg_ba,m,n),msq_struc17_5(igg_ab,m,n),
     & msq_struc17_5(igg_sym,m,n),sub17_5(qq))

**********************************
* 2-quark contributions for p6   *
**********************************
      msq(12,j,k)=msq(12,j,k)+f2q*sub2q3g_ii(
     & msq_struc26_7(igg_ba,m,n),msq_struc26_7(igg_ab,m,n),sub26_7(gg),
     & msq_struc26_7v(igg_ba,m,n),msq_struc26_7v(igg_ab,m,n),sub26_7v)
      msq(12,j,k)=msq(12,j,k)+f2q*half*sub2q3g_ii(
     & msq_struc26_7(igg_ba,m,n),msq_struc26_7(igg_ab,m,n),sub76_2(gg),
     & msq_struc76_2v(igg_ba,m,n),msq_struc76_2v(igg_ab,m,n),sub76_2v)

      msq(14,j,k)=msq(14,j,k)+f2q*sub2q3g_if(
     & msq_struc26_5(igg_ba,m,n),msq_struc26_5(igg_sym,m,n),
     & sub26_5(gg),
     & msq_struc26_5v(igg_ba,m,n),msq_struc26_5v(igg_sym,m,n),sub26_5v)
      msq(14,j,k)=msq(14,j,k)+f2q*sub2q3g_fi(
     & msq_struc26_5(igg_ba,m,n),msq_struc26_5(igg_sym,m,n),
     & sub56_2(qq))
      msq(4,j,k)=msq(4,j,k)+f2q*sub2q3g_if(
     & msq_struc26_1(igg_ab,m,n),msq_struc26_1(igg_sym,m,n),
     & sub26_1(gg),
     & msq_struc26_1v(igg_ab,m,n),msq_struc26_1v(igg_sym,m,n),sub26_1v)
      msq(3,j,k)=msq(3,j,k)+f2q*sub2q3g_fi(
     & msq_struc16_2(igg_ab,m,n),msq_struc16_2(igg_sym,m,n),
     & sub16_2(qq))
      msq(22,j,k)=msq(22,j,k)+f2q*half*sub2q3g_if(
     & msq_struc76_5(igg_ab,m,n),msq_struc76_5(igg_sym,m,n),
     & sub76_5(gg),
     & msq_struc76_5v(igg_ab,m,n),msq_struc76_5v(igg_sym,m,n),sub76_5v)
      msq(21,j,k)=msq(21,j,k)+f2q*sub2q3g_fi(
     & msq_struc56_7(igg_ab,m,n),msq_struc56_7(igg_sym,m,n),
     & sub56_7(qq))
      msq(13,j,k)=msq(13,j,k)+f2q*half*sub2q3g_if(
     & msq_struc16_7(igg_ba,m,n),msq_struc16_7(igg_sym,m,n),
     & sub76_1(gg),
     & msq_struc76_1v(igg_ba,m,n),msq_struc76_1v(igg_sym,m,n),sub76_1v)
      msq(13,j,k)=msq(13,j,k)+f2q*sub2q3g_fi(
     & msq_struc16_7(igg_ba,m,n),msq_struc16_7(igg_sym,m,n),
     & sub16_7(qq))

      msq(11,j,k)=msq(11,j,k)+f2q*sub2q3g_ff(
     & msq_struc16_5(igg_ba,m,n),msq_struc16_5(igg_ab,m,n),
     & msq_struc16_5(igg_sym,m,n),sub56_1(qq))
      msq(11,j,k)=msq(11,j,k)+f2q*sub2q3g_ff(
     & msq_struc16_5(igg_ba,m,n),msq_struc16_5(igg_ab,m,n),
     & msq_struc16_5(igg_sym,m,n),sub16_5(qq))

***********************************************
* 2-quark contributions for p5 (off-diagonal) *
***********************************************
      msq(5,j,k)=msq(5,j,k)+(aveqg/avegg)*f2q*(
     .+(msq_struc15_2(igggg_a,0,k)+msq_struc15_2(igggg_b,0,k)
     & +msq_struc15_2(igggg_c,0,k))*sub15_2(gq)
     .+(msq_struc15_2v(igggg_a,0,k)+msq_struc15_2v(igggg_b,0,k)
     & +msq_struc15_2v(igggg_c,0,k))*sub15_2v)

      m=j/abs(j)
      msq(6,j,k)=msq(6,j,k)+f2q*(
     &  msq_struc25_1(igg_ab,m,-m)+msq_struc25_1(igg_ba,m,-m)
     & +msq_struc25_1(igg_sym,m,-m))*sub25_1(qg)

****************************************
* 4-quark contributions (off-diagonal) *
****************************************

c--- this factor accounts for all the non-identical contributions
c--- (nflav-1) and also the identical one with symmetry factor 1/2
      nfactor=real(nflav-1,dp)+half

      m=+1
      msq(20,j,k)=msq(20,j,k)+nfactor*f4q*(
     .+(msq_struc67_5(igg_ab,m,0)
     & +msq_struc67_5(igg_sym,m,0))*sub67_5(gq)
     .-(msq_struc67_5v(igg_ab,m,0)
     & +msq_struc67_5v(igg_sym,m,0))*sub67_5v)/two

      msq(9,j,k)=msq(9,j,k)+nfactor*f4q*(
     .+(msq_struc17_6(igg_ba,m,0)
     & +msq_struc17_6(igg_sym,m,0))*sub67_1(gq)
     .-(msq_struc67_1v(igg_ba,m,0)
     & +msq_struc67_1v(igg_sym,m,0))*sub67_1v)/two

      msq(8,j,k)=msq(8,j,k)+nfactor*f4q*(
     .+(msq_struc27_6(igg_ab,m,0)
     & +msq_struc27_6(igg_ba,m,0))*sub67_2(gq)
     .-(msq_struc67_2v(igg_ab,m,0)
     & +msq_struc67_2v(igg_ba,m,0))*sub67_2v)/two

      msq(5,j,k)=msq(5,j,k)+(aveqg/avegg)*nfactor/real(nflav,dp)*f4q*(
     .+(msq_struc15_2(igg_ab,0,k)+msq_struc15_2(igg_ba,0,k)
     & +msq_struc15_2(igg_sym,0,k))*sub15_2(gq)
     .+(msq_struc15_2v(igg_ab,0,k)+msq_struc15_2v(igg_ba,0,k)
     & +msq_struc15_2v(igg_sym,0,k))*sub15_2v)

      msq(2,j,k)=msq(2,j,k)+real(nflav-1,dp)*f4q*(
     &  msq27_1(m,2))*sub27_1(qg)

      msq(4,j,k)=msq(4,j,k)+real(nflav-1,dp)*f4q*(
     &  msq26_1(m,-2))*sub26_1(qg)

      msq(6,j,k)=msq(6,j,k)+f4q*(
     &  msq_struc25_1(iqr,m,-m))*sub25_1(qg)

****************************************************************
* 4-quark contributions - identical quarks only (off-diagonal) *
****************************************************************

      msq(23,j,k)=msq(23,j,k)+half*f4q*(
     .+(msq_struc75_6(igg_ab,m,0)
     & +msq_struc75_6(igg_sym,m,0))*sub75_6(gq)
     .-(msq_struc75_6v(igg_ab,m,0)
     & +msq_struc75_6v(igg_sym,m,0))*sub75_6v)/two

      msq(15,j,k)=msq(15,j,k)+half*f4q*(
     .+(msq_struc15_7(igg_ba,m,0)
     & +msq_struc15_7(igg_sym,m,0))*sub75_1(gq)
     .-(msq_struc75_1v(igg_ba,m,0)
     & +msq_struc75_1v(igg_sym,m,0))*sub75_1v)/two

      msq(18,j,k)=msq(18,j,k)+half*f4q*(
     .+(msq_struc25_7(igg_ab,m,0)
     & +msq_struc25_7(igg_ba,m,0))*sub75_2(gq)
     .-(msq_struc75_2v(igg_ab,m,0)
     & +msq_struc75_2v(igg_ba,m,0))*sub75_2v)/two

      msq(2,j,k)=msq(2,j,k)+f4q*(
     &  msq27_1(m,1))*sub27_1(qg)

      msq(4,j,k)=msq(4,j,k)+half*f4q*(
     &  msq_struc26_1(iqq_a,m,-m)+msq_struc26_1(iqq_b,m,-m)
     & +msq_struc26_1(iqq_i,m,-m))*sub26_1(qg)

      msq(6,j,k)=msq(6,j,k)+half*f4q*(
     &  +msq_struc25_1(iqq_a,m,-1)+msq_struc25_1(iqq_b,m,-1)
     &  +msq_struc25_1(iqq_i,m,-1))*sub25_1(qg)
      msq(3,j,k)=msq(3,j,k)+(aveqg/avegg)*half/real(nflav,dp)*f4q*(
     .+(msq_struc16_2(igg_ab,0,k)+msq_struc16_2(igg_ba,0,k)
     & +msq_struc16_2(igg_sym,0,k))*sub16_2(gq)
     .+(msq_struc16_2v(igg_ab,0,k)+msq_struc16_2v(igg_ba,0,k)
     & +msq_struc16_2v(igg_sym,0,k))*sub16_2v)

************************************************************************
* SUBTRACTIONS FOR G-Q and G-QB                                        *
************************************************************************
      elseif ((j == 0).and.(k.ne.0)) then

c--- note that all subtractions are performed with reference to the
c--- gluon-quark matrix elements. We could call the gluon-antiquark
c--- LO, also switching _ab to _ba, but the answer would be the same
      m=0
      n=+1

********************************
* Note that in the sum over    *
* permutations below, factors  *
* of one half enter for 67     *
* dipoles, as they enter twice *
********************************

**********************************
* 2-quark contributions for p7   *
**********************************
      msq(9,j,k)=msq(9,j,k)+f2q*sub2q3g_ii(
     & msq_struc17_6(igg_ba,m,n),msq_struc17_6(igg_ab,m,n),sub17_6(gg),
     & msq_struc17_6v(igg_ba,m,n),msq_struc17_6v(igg_ab,m,n),sub17_6v)
      msq(9,j,k)=msq(9,j,k)+f2q*half*sub2q3g_ii(
     & msq_struc17_6(igg_ba,m,n),msq_struc17_6(igg_ab,m,n),sub67_1(gg),
     & msq_struc67_1v(igg_ba,m,n),msq_struc67_1v(igg_ab,m,n),sub67_1v)

      msq(7,j,k)=msq(7,j,k)+f2q*sub2q3g_if(
     & msq_struc17_5(igg_ba,m,n),msq_struc17_5(igg_sym,m,n),
     & sub17_5(gg),
     & msq_struc17_5v(igg_ba,m,n),msq_struc17_5v(igg_sym,m,n),sub17_5v)
      msq(7,j,k)=msq(7,j,k)+f2q*sub2q3g_fi(
     & msq_struc17_5(igg_ba,m,n),msq_struc17_5(igg_sym,m,n),
     & sub57_1(qq))
      msq(1,j,k)=msq(1,j,k)+f2q*sub2q3g_if(
     & msq_struc17_2(igg_ab,m,n),msq_struc17_2(igg_sym,m,n),
     & sub17_2(gg),
     & msq_struc17_2v(igg_ab,m,n),msq_struc17_2v(igg_sym,m,n),sub17_2v)
      msq(2,j,k)=msq(2,j,k)+f2q*sub2q3g_fi(
     & msq_struc27_1(igg_ab,m,n),msq_struc27_1(igg_sym,m,n),
     & sub27_1(qq))
      msq(20,j,k)=msq(20,j,k)+f2q*half*sub2q3g_if(
     & msq_struc67_5(igg_ab,m,n),msq_struc67_5(igg_sym,m,n),
     & sub67_5(gg),
     & msq_struc67_5v(igg_ab,m,n),msq_struc67_5v(igg_sym,m,n),sub67_5v)
      msq(19,j,k)=msq(19,j,k)+f2q*sub2q3g_fi(
     & msq_struc57_6(igg_ab,m,n),msq_struc57_6(igg_sym,m,n),
     & sub57_6(qq))
      msq(8,j,k)=msq(8,j,k)+f2q*half*sub2q3g_if(
     & msq_struc27_6(igg_ba,m,n),msq_struc27_6(igg_sym,m,n),
     & sub67_2(gg),
     & msq_struc67_2v(igg_ba,m,n),msq_struc67_2v(igg_sym,m,n),sub67_2v)
      msq(8,j,k)=msq(8,j,k)+f2q*sub2q3g_fi(
     & msq_struc27_6(igg_ba,m,n),msq_struc27_6(igg_sym,m,n),
     & sub27_6(qq))

      msq(10,j,k)=msq(10,j,k)+f2q*sub2q3g_ff(
     & msq_struc27_5(igg_ba,m,n),msq_struc27_5(igg_ab,m,n),
     & msq_struc27_5(igg_sym,m,n),sub57_2(qq))
      msq(10,j,k)=msq(10,j,k)+f2q*sub2q3g_ff(
     & msq_struc27_5(igg_ba,m,n),msq_struc27_5(igg_ab,m,n),
     & msq_struc27_5(igg_sym,m,n),sub27_5(qq))

**********************************
* 2-quark contributions for p6   *
**********************************
      msq(13,j,k)=msq(13,j,k)+f2q*sub2q3g_ii(
     & msq_struc16_7(igg_ba,m,n),msq_struc16_7(igg_ab,m,n),sub16_7(gg),
     & msq_struc16_7v(igg_ba,m,n),msq_struc16_7v(igg_ab,m,n),sub16_7v)
      msq(13,j,k)=msq(13,j,k)+f2q*half*sub2q3g_ii(
     & msq_struc16_7(igg_ba,m,n),msq_struc16_7(igg_ab,m,n),sub76_1(gg),
     & msq_struc76_1v(igg_ba,m,n),msq_struc76_1v(igg_ab,m,n),sub76_1v)

      msq(11,j,k)=msq(11,j,k)+f2q*sub2q3g_if(
     & msq_struc16_5(igg_ba,m,n),msq_struc16_5(igg_sym,m,n),
     & sub16_5(gg),
     & msq_struc16_5v(igg_ba,m,n),msq_struc16_5v(igg_sym,m,n),sub16_5v)
      msq(11,j,k)=msq(11,j,k)+f2q*sub2q3g_fi(
     & msq_struc16_5(igg_ba,m,n),msq_struc16_5(igg_sym,m,n),
     & sub56_1(qq))
      msq(3,j,k)=msq(3,j,k)+f2q*sub2q3g_if(
     & msq_struc16_2(igg_ab,m,n),msq_struc16_2(igg_sym,m,n),
     & sub16_2(gg),
     & msq_struc16_2v(igg_ab,m,n),msq_struc16_2v(igg_sym,m,n),sub16_2v)
      msq(4,j,k)=msq(4,j,k)+f2q*sub2q3g_fi(
     & msq_struc26_1(igg_ab,m,n),msq_struc26_1(igg_sym,m,n),
     & sub26_1(qq))
      msq(22,j,k)=msq(22,j,k)+f2q*half*sub2q3g_if(
     & msq_struc76_5(igg_ab,m,n),msq_struc76_5(igg_sym,m,n),
     & sub76_5(gg),
     & msq_struc76_5v(igg_ab,m,n),msq_struc76_5v(igg_sym,m,n),sub76_5v)
      msq(21,j,k)=msq(21,j,k)+f2q*sub2q3g_fi(
     & msq_struc56_7(igg_ab,m,n),msq_struc56_7(igg_sym,m,n),
     & sub56_7(qq))
      msq(12,j,k)=msq(12,j,k)+f2q*half*sub2q3g_if(
     & msq_struc26_7(igg_ba,m,n),msq_struc26_7(igg_sym,m,n),
     & sub76_2(gg),
     & msq_struc76_2v(igg_ba,m,n),msq_struc76_2v(igg_sym,m,n),sub76_2v)
      msq(12,j,k)=msq(12,j,k)+f2q*sub2q3g_fi(
     & msq_struc26_7(igg_ba,m,n),msq_struc26_7(igg_sym,m,n),
     & sub26_7(qq))

      msq(14,j,k)=msq(14,j,k)+f2q*sub2q3g_ff(
     & msq_struc26_5(igg_ba,m,n),msq_struc26_5(igg_ab,m,n),
     & msq_struc26_5(igg_sym,m,n),sub56_2(qq))
      msq(14,j,k)=msq(14,j,k)+f2q*sub2q3g_ff(
     & msq_struc26_5(igg_ba,m,n),msq_struc26_5(igg_ab,m,n),
     & msq_struc26_5(igg_sym,m,n),sub26_5(qq))

***********************************************
* 2-quark contributions for p5 (off-diagonal) *
***********************************************
      msq(6,j,k)=msq(6,j,k)+(aveqg/avegg)*f2q*(
     .+(msq_struc25_1(igggg_a,j,0)+msq_struc25_1(igggg_b,j,0)
     & +msq_struc25_1(igggg_c,j,0))*sub25_1(gq)
     .+(msq_struc25_1v(igggg_a,j,0)+msq_struc25_1v(igggg_b,j,0)
     & +msq_struc25_1v(igggg_c,j,0))*sub25_1v)

      m=k/abs(k)
      msq(5,j,k)=msq(5,j,k)+f2q*(
     &  msq_struc15_2(igg_ab,-m,m)+msq_struc15_2(igg_ba,-m,m)
     & +msq_struc15_2(igg_sym,-m,m))*sub15_2(qg)

****************************************
* 4-quark contributions (off-diagonal) *
****************************************

c--- this factor accounts for all the non-identical contributions
c--- (nflav-1) and also the identical one with symmetry factor 1/2
      nfactor=real(nflav-1,dp)+half

      n=+1
      msq(20,j,k)=msq(20,j,k)+nfactor*f4q*(
     .+(msq_struc67_5(igg_ab,0,n)
     & +msq_struc67_5(igg_sym,0,n))*sub67_5(gq)
     .-(msq_struc67_5v(igg_ab,0,n)
     & +msq_struc67_5v(igg_sym,0,n))*sub67_5v)/two

      msq(9,j,k)=msq(9,j,k)+nfactor*f4q*(
     .+(msq_struc17_6(igg_ab,0,n)
     & +msq_struc17_6(igg_ba,0,n))*sub67_1(gq)
     .-(msq_struc67_1v(igg_ab,0,n)
     & +msq_struc67_1v(igg_ba,0,n))*sub67_1v)/two

      msq(8,j,k)=msq(8,j,k)+nfactor*f4q*(
     .+(msq_struc27_6(igg_ba,0,n)
     & +msq_struc27_6(igg_sym,0,n))*sub67_2(gq)
     .-(msq_struc67_2v(igg_ba,0,n)
     & +msq_struc67_2v(igg_sym,0,n))*sub67_2v)/two

      msq(6,j,k)=msq(6,j,k)+(aveqg/avegg)*nfactor/real(nflav,dp)*f4q*(
     .+(msq_struc25_1(igg_ab,j,0)+msq_struc25_1(igg_ba,j,0)
     & +msq_struc25_1(igg_sym,j,0))*sub25_1(gq)
     .+(msq_struc25_1v(igg_ab,j,0)+msq_struc25_1v(igg_ba,j,0)
     & +msq_struc25_1v(igg_sym,j,0))*sub25_1v)

c--- we use a factor of two and msq_struc(iqq_b,...) since
c---  iqq_b = iqq_a (5 <--> 6) and we want qQ->Qq and not qQ->qQ
      msq(1,j,k)=msq(1,j,k)+two*real(nflav-1,dp)*f4q*(
     &  msq_struc17_2(iqq_b,1,n))*sub17_2(qg)

c--- we use msq_struc(iqq_b,...) since
c---  iqq_b = iqq_a (5 <--> 6) and we want qQ->Qq and not qQ->qQ
      msq(3,j,k)=msq(3,j,k)+two*nfactor*f4q*(
     &  msq_struc16_2(iqq_b,1,n))*sub16_2(qg)

      msq(5,j,k)=msq(5,j,k)+f4q*(
     &  msq_struc15_2(iqr,n,-n))*sub15_2(qg)

****************************************************************
* 4-quark contributions - identical quarks only (off-diagonal) *
****************************************************************

      msq(23,j,k)=msq(23,j,k)+half*f4q*(
     .+(msq_struc75_6(igg_ab,0,n)
     & +msq_struc75_6(igg_sym,0,n))*sub75_6(gq)
     .-(msq_struc75_6v(igg_ab,0,n)
     & +msq_struc75_6v(igg_sym,0,n))*sub75_6v)/two

      msq(15,j,k)=msq(15,j,k)+half*f4q*(
     .+(msq_struc15_7(igg_ba,0,n)
     & +msq_struc15_7(igg_ab,0,n))*sub75_1(gq)
     .-(msq_struc75_1v(igg_ba,0,n)
     & +msq_struc75_1v(igg_ab,0,n))*sub75_1v)/two

      msq(18,j,k)=msq(18,j,k)+half*f4q*(
     .+(msq_struc25_7(igg_sym,0,n)
     & +msq_struc25_7(igg_ba,0,n))*sub75_2(gq)
     .-(msq_struc75_2v(igg_sym,0,n)
     & +msq_struc75_2v(igg_ba,0,n))*sub75_2v)/two

      msq(1,j,k)=msq(1,j,k)+f4q*(
     &  msq17_2(1,n))*sub17_2(qg)
      msq(3,j,k)=msq(3,j,k)+half/real(nflav-1,dp)*f4q*(
     &  msq_struc16_2(iqr,-1,n))*sub16_2(qg)
c--- note: this term uses the qa->aq new matrix element which
c--- is artificially stored in msq_struc(iqr,0,0)
      msq(5,j,k)=msq(5,j,k)+half*f4q*(
     &  +msq_struc15_2(iqr,0,0))*sub15_2(qg)
      msq(4,j,k)=msq(4,j,k)+(aveqg/avegg)*half/real(nflav,dp)*f4q*(
     .+(msq_struc26_1(igg_ab,j,0)+msq_struc26_1(igg_ba,j,0)
     & +msq_struc26_1(igg_sym,j,0))*sub26_1(gq)
     .+(msq_struc26_1v(igg_ab,j,0)+msq_struc26_1v(igg_ba,j,0)
     & +msq_struc26_1v(igg_sym,j,0))*sub26_1v)

************************************************************************
* SUBTRACTIONS FOR G-G                                                 *
************************************************************************
      elseif ((j == 0).and.(k == 0)) then
********************************
* Note that in the sum over    *
* permutations below, factors  *
* of 1/2 enter for 56, 57, 67  *
* dipoles, as they enter twice *
********************************


**********************************
* 0-quark contributions for p7   *
**********************************
      msq(1,j,k)=two*xn*(
     & msq_struc17_2(igggg_a,j,k)*sub17_2(gg)
     .+msq_struc17_2v(igggg_a,j,k)*sub17_2v
     .+msq_struc17_2(igggg_c,j,k)*sub17_2(gg)
     .+msq_struc17_2v(igggg_c,j,k)*sub17_2v)
      msq(2,j,k)=two*xn*(
     & msq_struc27_1(igggg_a,j,k)*sub27_1(gg)
     .+msq_struc27_1v(igggg_a,j,k)*sub27_1v
     .+msq_struc27_1(igggg_c,j,k)*sub27_1(gg)
     .+msq_struc27_1v(igggg_c,j,k)*sub27_1v)

      msq(7,j,k)=two*xn*(
     .+msq_struc17_5(igggg_b,j,k)*sub17_5(gg)
     .+msq_struc17_5v(igggg_b,j,k)*sub17_5v
     .+msq_struc17_5(igggg_b,j,k)*sub57_1(gg)*half
     .+msq_struc57_1v(igggg_b,j,k)*sub57_1v*half
     .+msq_struc17_5(igggg_c,j,k)*sub17_5(gg)
     .+msq_struc17_5v(igggg_c,j,k)*sub17_5v
     .+msq_struc17_5(igggg_c,j,k)*sub57_1(gg)*half
     .+msq_struc57_1v(igggg_c,j,k)*sub57_1v*half)
      msq(9,j,k)=two*xn*(
     .+msq_struc17_6(igggg_a,j,k)*sub17_6(gg)
     .+msq_struc17_6v(igggg_a,j,k)*sub17_6v
     .+msq_struc17_6(igggg_a,j,k)*sub67_1(gg)*half
     .+msq_struc67_1v(igggg_a,j,k)*sub67_1v*half
     .+msq_struc17_6(igggg_b,j,k)*sub17_6(gg)
     .+msq_struc17_6v(igggg_b,j,k)*sub17_6v
     .+msq_struc17_6(igggg_b,j,k)*sub67_1(gg)*half
     .+msq_struc67_1v(igggg_b,j,k)*sub67_1v*half)
      msq(10,j,k)=two*xn*(
     .+msq_struc27_5(igggg_a,j,k)*sub27_5(gg)
     .+msq_struc27_5v(igggg_a,j,k)*sub27_5v
     .+msq_struc27_5(igggg_a,j,k)*sub57_2(gg)*half
     .+msq_struc57_2v(igggg_a,j,k)*sub57_2v*half
     .+msq_struc27_5(igggg_b,j,k)*sub27_5(gg)
     .+msq_struc27_5v(igggg_b,j,k)*sub27_5v
     .+msq_struc27_5(igggg_b,j,k)*sub57_2(gg)*half
     .+msq_struc57_2v(igggg_b,j,k)*sub57_2v*half)
      msq(8,j,k)=two*xn*(
     .+msq_struc27_6(igggg_b,j,k)*sub27_6(gg)
     .+msq_struc27_6v(igggg_b,j,k)*sub27_6v
     .+msq_struc27_6(igggg_b,j,k)*sub67_2(gg)*half
     .+msq_struc67_2v(igggg_b,j,k)*sub67_2v*half
     .+msq_struc27_6(igggg_c,j,k)*sub27_6(gg)
     .+msq_struc27_6v(igggg_c,j,k)*sub27_6v
     .+msq_struc27_6(igggg_c,j,k)*sub67_2(gg)*half
     .+msq_struc67_2v(igggg_c,j,k)*sub67_2v*half)

      msq(19,j,k)=two*xn*(
     & msq_struc57_6(igggg_a,j,k)*sub57_6(gg)
     .+msq_struc57_6v(igggg_a,j,k)*sub57_6v
     .+msq_struc57_6(igggg_c,j,k)*sub57_6(gg)
     .+msq_struc57_6v(igggg_c,j,k)*sub57_6v)*half
      msq(20,j,k)=two*xn*(
     & msq_struc67_5(igggg_a,j,k)*sub67_5(gg)
     .+msq_struc67_5v(igggg_a,j,k)*sub67_5v
     .+msq_struc67_5(igggg_c,j,k)*sub67_5(gg)
     .+msq_struc67_5v(igggg_c,j,k)*sub67_5v)*half

**********************************
* 0-quark contributions for p6   *
**********************************
      msq(3,j,k)=two*xn*(
     & msq_struc16_2(igggg_a,j,k)*sub16_2(gg)
     .+msq_struc16_2v(igggg_a,j,k)*sub16_2v
     .+msq_struc16_2(igggg_c,j,k)*sub16_2(gg)
     .+msq_struc16_2v(igggg_c,j,k)*sub16_2v)
      msq(4,j,k)=two*xn*(
     & msq_struc26_1(igggg_a,j,k)*sub26_1(gg)
     .+msq_struc26_1v(igggg_a,j,k)*sub26_1v
     .+msq_struc26_1(igggg_c,j,k)*sub26_1(gg)
     .+msq_struc26_1v(igggg_c,j,k)*sub26_1v)

      msq(11,j,k)=two*xn*(
     .+msq_struc16_5(igggg_b,j,k)*sub16_5(gg)
     .+msq_struc16_5v(igggg_b,j,k)*sub16_5v
     .+msq_struc16_5(igggg_b,j,k)*sub56_1(gg)*half
     .+msq_struc56_1v(igggg_b,j,k)*sub56_1v*half
     .+msq_struc16_5(igggg_c,j,k)*sub16_5(gg)
     .+msq_struc16_5v(igggg_c,j,k)*sub16_5v
     .+msq_struc16_5(igggg_c,j,k)*sub56_1(gg)*half
     .+msq_struc56_1v(igggg_c,j,k)*sub56_1v*half)

      msq(13,j,k)=two*xn*(
     .+msq_struc16_7(igggg_a,j,k)*sub16_7(gg)
     .+msq_struc16_7v(igggg_a,j,k)*sub16_7v
     .+msq_struc16_7(igggg_a,j,k)*sub76_1(gg)*half
     .+msq_struc76_1v(igggg_a,j,k)*sub76_1v*half
     .+msq_struc16_7(igggg_b,j,k)*sub16_7(gg)
     .+msq_struc16_7v(igggg_b,j,k)*sub16_7v
     .+msq_struc16_7(igggg_b,j,k)*sub76_1(gg)*half
     .+msq_struc76_1v(igggg_b,j,k)*sub76_1v*half)
      msq(14,j,k)=two*xn*(
     .+msq_struc26_5(igggg_a,j,k)*sub26_5(gg)
     .+msq_struc26_5v(igggg_a,j,k)*sub26_5v
     .+msq_struc26_5(igggg_a,j,k)*sub56_2(gg)*half
     .+msq_struc56_2v(igggg_a,j,k)*sub56_2v*half
     .+msq_struc26_5(igggg_b,j,k)*sub26_5(gg)
     .+msq_struc26_5v(igggg_b,j,k)*sub26_5v
     .+msq_struc26_5(igggg_b,j,k)*sub56_2(gg)*half
     .+msq_struc56_2v(igggg_b,j,k)*sub56_2v*half)
      msq(12,j,k)=two*xn*(
     .+msq_struc26_7(igggg_b,j,k)*sub26_7(gg)
     .+msq_struc26_7v(igggg_b,j,k)*sub26_7v
     .+msq_struc26_7(igggg_b,j,k)*sub76_2(gg)*half
     .+msq_struc76_2v(igggg_b,j,k)*sub76_2v*half
     .+msq_struc26_7(igggg_c,j,k)*sub26_7(gg)
     .+msq_struc26_7v(igggg_c,j,k)*sub26_7v
     .+msq_struc26_7(igggg_c,j,k)*sub76_2(gg)*half
     .+msq_struc76_2v(igggg_c,j,k)*sub76_2v*half)

      msq(21,j,k)=two*xn*(
     & msq_struc56_7(igggg_a,j,k)*sub56_7(gg)
     .+msq_struc56_7v(igggg_a,j,k)*sub56_7v
     .+msq_struc56_7(igggg_c,j,k)*sub56_7(gg)
     .+msq_struc56_7v(igggg_c,j,k)*sub56_7v)*half
      msq(22,j,k)=two*xn*(
     & msq_struc76_5(igggg_a,j,k)*sub76_5(gg)
     .+msq_struc76_5v(igggg_a,j,k)*sub76_5v
     .+msq_struc76_5(igggg_c,j,k)*sub76_5(gg)
     .+msq_struc76_5v(igggg_c,j,k)*sub76_5v)*half

**********************************
* 0-quark contributions for p5   *
**********************************
      msq(5,j,k)=two*xn*(
     & msq_struc15_2(igggg_a,j,k)*sub15_2(gg)
     .+msq_struc15_2v(igggg_a,j,k)*sub15_2v
     .+msq_struc15_2(igggg_c,j,k)*sub15_2(gg)
     .+msq_struc15_2v(igggg_c,j,k)*sub15_2v)
      msq(6,j,k)=two*xn*(
     & msq_struc25_1(igggg_a,j,k)*sub25_1(gg)
     .+msq_struc25_1v(igggg_a,j,k)*sub25_1v
     .+msq_struc25_1(igggg_c,j,k)*sub25_1(gg)
     .+msq_struc25_1v(igggg_c,j,k)*sub25_1v)

      msq(15,j,k)=two*xn*(
     .+msq_struc15_7(igggg_b,j,k)*sub15_7(gg)
     .+msq_struc15_7v(igggg_b,j,k)*sub15_7v
     .+msq_struc15_7(igggg_b,j,k)*sub75_1(gg)*half
     .+msq_struc75_1v(igggg_b,j,k)*sub75_1v*half
     .+msq_struc15_7(igggg_a,j,k)*sub15_7(gg)
     .+msq_struc15_7v(igggg_a,j,k)*sub15_7v
     .+msq_struc15_7(igggg_a,j,k)*sub75_1(gg)*half
     .+msq_struc75_1v(igggg_a,j,k)*sub75_1v*half)
      msq(17,j,k)=two*xn*(
     .+msq_struc15_6(igggg_c,j,k)*sub15_6(gg)
     .+msq_struc15_6v(igggg_c,j,k)*sub15_6v
     .+msq_struc15_6(igggg_c,j,k)*sub65_1(gg)*half
     .+msq_struc65_1v(igggg_c,j,k)*sub65_1v*half
     .+msq_struc15_6(igggg_b,j,k)*sub15_6(gg)
     .+msq_struc15_6v(igggg_b,j,k)*sub15_6v
     .+msq_struc15_6(igggg_b,j,k)*sub65_1(gg)*half
     .+msq_struc65_1v(igggg_b,j,k)*sub65_1v*half)
      msq(18,j,k)=two*xn*(
     .+msq_struc25_7(igggg_c,j,k)*sub25_7(gg)
     .+msq_struc25_7v(igggg_c,j,k)*sub25_7v
     .+msq_struc25_7(igggg_c,j,k)*sub75_2(gg)*half
     .+msq_struc75_2v(igggg_c,j,k)*sub75_2v*half
     .+msq_struc25_7(igggg_b,j,k)*sub25_7(gg)
     .+msq_struc25_7v(igggg_b,j,k)*sub25_7v
     .+msq_struc25_7(igggg_b,j,k)*sub75_2(gg)*half
     .+msq_struc75_2v(igggg_b,j,k)*sub75_2v*half)
      msq(16,j,k)=two*xn*(
     .+msq_struc25_6(igggg_b,j,k)*sub25_6(gg)
     .+msq_struc25_6v(igggg_b,j,k)*sub25_6v
     .+msq_struc25_6(igggg_b,j,k)*sub65_2(gg)*half
     .+msq_struc65_2v(igggg_b,j,k)*sub65_2v*half
     .+msq_struc25_6(igggg_a,j,k)*sub25_6(gg)
     .+msq_struc25_6v(igggg_a,j,k)*sub25_6v
     .+msq_struc25_6(igggg_a,j,k)*sub65_2(gg)*half
     .+msq_struc65_2v(igggg_a,j,k)*sub65_2v*half)

      msq(23,j,k)=two*xn*(
     & msq_struc75_6(igggg_a,j,k)*sub75_6(gg)
     .+msq_struc75_6v(igggg_a,j,k)*sub75_6v
     .+msq_struc75_6(igggg_c,j,k)*sub75_6(gg)
     .+msq_struc75_6v(igggg_c,j,k)*sub75_6v)*half
      msq(24,j,k)=two*xn*(
     & msq_struc65_7(igggg_a,j,k)*sub65_7(gg)
     .+msq_struc65_7v(igggg_a,j,k)*sub65_7v
     .+msq_struc65_7(igggg_c,j,k)*sub65_7(gg)
     .+msq_struc65_7v(igggg_c,j,k)*sub65_7v)*half

      do nd=1,ndmax
        msq(nd,j,k)=f0q*msq(nd,j,k)/six
      enddo

      msq(1,j,k)=msq(1,j,k)+two*f2q*sub2q3g_ii(
     & msq_struc17_2(igg_ba,j,k),msq_struc17_2(igg_ab,j,k),sub17_2(gg),
     & msq_struc17_2v(igg_ba,j,k),msq_struc17_2v(igg_ab,j,k),sub17_2v)
      msq(2,j,k)=msq(2,j,k)+two*f2q*sub2q3g_ii(
     & msq_struc27_1(igg_ba,j,k),msq_struc27_1(igg_ab,j,k),sub27_1(gg),
     & msq_struc27_1v(igg_ba,j,k),msq_struc27_1v(igg_ab,j,k),sub27_1v)

      msq(7,j,k)=msq(7,j,k)+two*f2q*sub2q3g_iffi(
     & msq_struc17_5(igg_ba,j,k),msq_struc17_5(igg_sym,j,k),
     & sub17_5(gg),sub57_1(qq),
     & msq_struc17_5v(igg_ba,j,k),msq_struc17_5v(igg_sym,j,k),sub17_5v)
      msq(9,j,k)=msq(9,j,k)+two*f2q*sub2q3g_iffi(
     & msq_struc17_6(igg_ab,j,k),msq_struc17_6(igg_sym,j,k),
     & sub17_6(gg),sub67_1(qq),
     & msq_struc17_6v(igg_ab,j,k),msq_struc17_6v(igg_sym,j,k),sub17_6v)
      msq(10,j,k)=msq(10,j,k)+two*f2q*sub2q3g_iffi(
     & msq_struc27_5(igg_ab,j,k),msq_struc27_5(igg_sym,j,k),
     & sub27_5(gg),sub57_2(qq),
     & msq_struc27_5v(igg_ab,j,k),msq_struc27_5v(igg_sym,j,k),sub27_5v)
      msq(8,j,k)=msq(8,j,k)+two*f2q*sub2q3g_iffi(
     & msq_struc27_6(igg_ba,j,k),msq_struc27_6(igg_sym,j,k),
     & sub27_6(gg),sub67_2(qq),
     & msq_struc27_6v(igg_ba,j,k),msq_struc27_6v(igg_sym,j,k),sub27_6v)

      msq(19,j,k)=msq(19,j,k)+two*f2q*sub2q3g_ff(
     & msq_struc57_6(igg_ba,j,k),msq_struc57_6(igg_ab,j,k),
     & msq_struc57_6(igg_sym,j,k),sub57_6(qq))

      msq(20,j,k)=msq(20,j,k)+two*f2q*sub2q3g_ff(
     & msq_struc67_5(igg_ba,j,k),msq_struc67_5(igg_ab,j,k),
     & msq_struc67_5(igg_sym,j,k),sub67_5(qq))

c--- note, the form of these subtractions corresponds to the
c--- f0q subtractions above
      msq(21,j,k)=msq(21,j,k)+two*f2q*real(nflav,dp)*(
     .+(msq_struc56_7(igggg_a,j,k)
     & +msq_struc56_7(igggg_c,j,k))*sub56_7(gq)
     .-(msq_struc56_7v(igggg_a,j,k)
     & +msq_struc56_7v(igggg_c,j,k))*sub56_7v)/two

      msq(11,j,k)=msq(11,j,k)+two*f2q*real(nflav,dp)*(
     .+(msq_struc16_5(igggg_b,j,k)
     & +msq_struc16_5(igggg_c,j,k))*sub56_1(gq)
     .-(msq_struc56_1v(igggg_b,j,k)
     & +msq_struc56_1v(igggg_c,j,k))*sub56_1v)/two

      msq(14,j,k)=msq(14,j,k)+two*f2q*real(nflav,dp)*(
     .+(msq_struc26_5(igggg_a,j,k)
     & +msq_struc26_5(igggg_b,j,k))*sub56_2(gq)
     .-(msq_struc56_2v(igggg_a,j,k)
     & +msq_struc56_2v(igggg_b,j,k))*sub56_2v)/two

      msq(4,j,k)=msq(4,j,k)+f2q*real(nflav,dp)*(
     &  msq_struc26_1(igg_ab,j,+1)+msq_struc26_1(igg_ba,j,+1)
     & +msq_struc26_1(igg_sym,j,+1))*sub26_1(qg)

      msq(3,j,k)=msq(3,j,k)+f2q*real(nflav,dp)*(
     &  msq_struc16_2(igg_ab,+1,k)+msq_struc16_2(igg_ba,+1,k)
     & +msq_struc16_2(igg_sym,+1,k))*sub16_2(qg)

      msq(6,j,k)=msq(6,j,k)+f2q*real(nflav,dp)*(
     &  msq_struc25_1(igg_ab,j,-1)+msq_struc25_1(igg_ba,j,-1)
     & +msq_struc25_1(igg_sym,j,-1))*sub25_1(qg)

      msq(5,j,k)=msq(5,j,k)+f2q*real(nflav,dp)*(
     &  msq_struc15_2(igg_ab,-1,k)+msq_struc15_2(igg_ba,-1,k)
     & +msq_struc15_2(igg_sym,-1,k))*sub15_2(qg)

      endif

      enddo
      enddo

      return
      end

