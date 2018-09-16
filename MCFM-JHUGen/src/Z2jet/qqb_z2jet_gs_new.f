      subroutine qqb_z2jet_gs_new(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K.Ellis                                                *
*     September, 2001.                                                     *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z + f(p5) + f(p6) + g(p7)
c                           |
c                           --> l(p3) + lbar(p4)
c     where the fermions are either q(p5) and qbar(p6) [Qflag = .true.]
c                                or g(p5) and g(p6)    [Gflag = .true.]
c
c---This routine does not function if both Qflag and Gflag=.true.
c
c--- all momenta are incoming


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'flags.f'
      include 'lc.f'
      include 'pp.f'

c--- np6,np12,... = n+6,n+12,...
c -- nd counts the dipoles
      integer:: j,k,n,np6,np12,np18,np21,nd
c--- slightly obtuse notation, fn=-nf, to simplify declaration lines
      real(dp):: p(mxpart,4),msq(maxd,fn:nf,fn:nf)
      real(dp)::
     & msq17_2(fn:nf,fn:nf),msq27_1(fn:nf,fn:nf),
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq16_2(fn:nf,fn:nf),msq26_1(fn:nf,fn:nf),
     & msq17_5(fn:nf,fn:nf),msq17_6(fn:nf,fn:nf),
     & msq15_6(fn:nf,fn:nf),msq65_7(fn:nf,fn:nf),
     & msq27_5(fn:nf,fn:nf),msq27_6(fn:nf,fn:nf),
     & msq57_6(fn:nf,fn:nf),msq67_5(fn:nf,fn:nf),
     & msq57_6v(fn:nf,fn:nf),msq26_1v(fn:nf,fn:nf),
     & msq15_2v(fn:nf,fn:nf),msq16_2v(fn:nf,fn:nf),
     & msq17_2v(fn:nf,fn:nf),msq25_6(fn:nf,fn:nf),
     & dummy(fn:nf,fn:nf),dummyv(fn:nf,fn:nf),
     & sub17_2(4),sub27_1(4),sub15_2(4),sub25_1(4),
     & sub16_2(4),sub26_1(4),sub25_6(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),
     & sub15_6(4),sub65_1(4),sub65_2(4),sub65_7(4),
     & sub57_6(4),sub67_5(4),dsubv,dsub(4),
     & sub57_6v,sub57_1v,sub57_2v,sub67_1v,sub67_2v,sub17_2v,sub65_7v,
     & sub26_1v,sub25_1v,sub27_1v,sub15_2v,sub16_2v,sub65_1v,sub65_2v,
     & sub67_5v
      real(dp)::
     & m17_2(0:2,fn:nf,fn:nf),m27_1(0:2,fn:nf,fn:nf),
     & m16_2(0:2,fn:nf,fn:nf),m26_1(0:2,fn:nf,fn:nf),
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m57_6(0:2,fn:nf,fn:nf),m67_5(0:2,fn:nf,fn:nf),
     & m65_7(0:2,fn:nf,fn:nf),
     & m17_5(0:2,fn:nf,fn:nf),
     & m17_6(0:2,fn:nf,fn:nf),
     & m27_5(0:2,fn:nf,fn:nf),
     & m27_6(0:2,fn:nf,fn:nf),
     & m15_6(0:2,fn:nf,fn:nf),
     & m25_6(0:2,fn:nf,fn:nf),
     & m17_2x(0:2,ppmax),
     & m27_1x(0:2,ppmax),
     & m15_2x(0:2,ppmax),
     & m16_2x(0:2,ppmax),
     & m25_1x(0:2,ppmax),
     & m26_1x(0:2,ppmax),
     & m57_6x(0:2,ppmax),
     & m65_7x(0:2,ppmax),
     & m67_5x(0:2,ppmax),
     & m17_5x(0:2,ppmax),
     & m17_6x(0:2,ppmax),
     & m27_5x(0:2,ppmax),
     & m27_6x(0:2,ppmax),
     & m15_6x(0:2,ppmax),
     & m25_6x(0:2,ppmax)

      real(dp)::
     & msq1a_b(6,0:2,fn:nf,fn:nf),msqba_1(6,0:2,fn:nf,fn:nf),
     & msqab_c(6,0:2,fn:nf,fn:nf),
     & msqbc_2(6,0:2,fn:nf,fn:nf),msq2c_b(6,0:2,fn:nf,fn:nf),
     & sub1a_b(6,4),subba_1(6,4),subab_c(6,4),
     & subbc_2(6,4),sub2c_b(6,4),
     & msq1a_bv(6,0:2,fn:nf,fn:nf),msqba_1v(6,0:2,fn:nf,fn:nf),
     & msqab_cv(6,0:2,fn:nf,fn:nf),
     & msqbc_2v(6,0:2,fn:nf,fn:nf),msq2c_bv(6,0:2,fn:nf,fn:nf),
     & sub1a_bv(6),subba_1v(6),subab_cv(6),
     & subbc_2v(6),sub2c_bv(6),
     & msq1b_2(6,0:2,fn:nf,fn:nf),msq2b_1(6,0:2,fn:nf,fn:nf),
     & sub1b_2(6,4),sub2b_1(6,4),
     & msq1b_2v(6,0:2,fn:nf,fn:nf),msq2b_1v(6,0:2,fn:nf,fn:nf),
     & sub1b_2v(6),sub2b_1v(6)
      real(dp)::
     & m57_1g(0:2,fn:nf,fn:nf),m57_1vg(0:2,fn:nf,fn:nf),
     & m57_1vx(ppmax),
     & m26_1g(0:2,fn:nf,fn:nf),m26_1vg(0:2,fn:nf,fn:nf),
     & m26_1vx(ppmax),
     & m25_1g(0:2,fn:nf,fn:nf),m25_1vg(0:2,fn:nf,fn:nf),
     & m25_1vx(ppmax),
     & m27_1g(0:2,fn:nf,fn:nf),m27_1vg(0:2,fn:nf,fn:nf),
     & m27_1vx(ppmax),
     & m67_1g(0:2,fn:nf,fn:nf),m67_1vg(0:2,fn:nf,fn:nf),
     & m67_1vx(ppmax),
     & m67_2g(0:2,fn:nf,fn:nf),m67_2vg(0:2,fn:nf,fn:nf),
     & m67_2vx(ppmax),
     & m57_2g(0:2,fn:nf,fn:nf),m57_2vg(0:2,fn:nf,fn:nf),
     & m57_2vx(ppmax),
     & m15_2g(0:2,fn:nf,fn:nf),m15_2vg(0:2,fn:nf,fn:nf),
     & m15_2vx(ppmax),
     & m16_2g(0:2,fn:nf,fn:nf),m16_2vg(0:2,fn:nf,fn:nf),
     & m17_2g(0:2,fn:nf,fn:nf),m17_2vg(0:2,fn:nf,fn:nf),
     & m16_2vx(ppmax),
     & m17_2vx(ppmax),
     & m65_1g(0:2,fn:nf,fn:nf),m65_1vg(0:2,fn:nf,fn:nf),
     & m65_1vx(ppmax),
     & m65_2g(0:2,fn:nf,fn:nf),m65_2vg(0:2,fn:nf,fn:nf),
     & m65_2vx(ppmax),
     & m15_6g(0:2,fn:nf,fn:nf),m25_6g(0:2,fn:nf,fn:nf),
     & m67_5vx(ppmax),m67_5vg(0:2,fn:nf,fn:nf),
     & m67_5g(0:2,fn:nf,fn:nf),
     & m65_7vx(ppmax),m65_7vg(0:2,fn:nf,fn:nf),
     & m65_7g(0:2,fn:nf,fn:nf),
     & m57_6vx(ppmax),m57_6vg(0:2,fn:nf,fn:nf),
     & m57_6g(0:2,fn:nf,fn:nf)

      real(dp):: mqq(0:2,fn:nf,fn:nf),
     & msqx(0:2,ppmax),mg(0:2,-nf:nf,-nf:nf),
     & mvg(0:2,-nf:nf,-nf:nf),mvxg(ppmax)

      external qqb_z2jet,qqb_z2jet_gvec
      external qqb_z2jetx_new,qqb_z2jet_gvecx_new,donothing_gvecx_new

      integer,parameter:: a(6)=(/5,5,7,6,6,7/),b(6)=(/6,7,5,7,5,6/),
     & c(6)=(/7,6,6,5,7,5/),
     & jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/),
     & kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/),
     & pntr(5:7,5:7)=reshape((/0,2,2,1,0,2,1,1,0/),(/3,3/))

      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file options.DAT'
        write(6,*) 'Failed in qqb_z2jet_gs.f'
        stop
      endif

      ndmax=24

c-- initialize the matrix elements to zero
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

      if (Gflag) then

c---arguments of dips:
c---    1  dipole number
c---    2  momentum
c---    3  emitter
c---    4  emitted
c---    5  spectator
c---    6  interference-free subtraction, equiv. to AP kernel for qq,qg
c---    7  correlation piece of subtraction, relevant only for gg,gq
c---    8  lowest order matrix elements at rescaled momentum, msq(j,k)
c---    9  lowest order matrix elements at rescaled momentum
c---        with emitter contracted with appropriate vector, msqv(j,k)
c---   10  appropriate lowest order calculating routine for msq
c---   11  appropriate lowest order calculating routine for msqv

c--- final-final
      do n=1,6
      call dips(n,p,a(n),b(n),c(n),dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msqab_c,msqab_cv,dsub,dsubv,subab_c,subab_cv,n)
      enddo

c--- initial-final/final-initial
      do n=1,6
      np6=n+6
      np12=n+12
      call dips(np6,p,1,a(n),b(n),dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq1a_b,msq1a_bv,dsub,dsubv,sub1a_b,sub1a_bv,n)
      call dips(np6,p,b(n),a(n),1,dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msqba_1,msqba_1v,dsub,dsubv,subba_1,subba_1v,n)
      call dips(np12,p,2,c(n),b(n),dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq2c_b,msq2c_bv,dsub,dsubv,sub2c_b,sub2c_bv,n)
      call dips(np12,p,b(n),c(n),2,dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msqbc_2,msqbc_2v,dsub,dsubv,subbc_2,subbc_2v,n)
      enddo

c--- initial-initial
      do n=1,3
      np18=n+18
      np21=n+21
      call dips(np18,p,1,b(n),2,dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq1b_2,msq1b_2v,dsub,dsubv,sub1b_2,sub1b_2v,n)
      call dips(np21,p,2,b(n),1,dsub,dsubv,dummy,dummyv,
     & qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq2b_1,msq2b_1v,dsub,dsubv,sub2b_1,sub2b_1v,n)
      enddo

c-- fill the matrix elements
      do j=-nf,nf
      do k=-nf,nf

c--- QUARK-ANTIQUARK contributions
      if    ((j>0).and.(k<0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,6
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *(msqab_c(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_c(n,pntr(c(n),a(n)),j,k))*xn/3._dp
     &                 +subab_cv(n)/2._dp
     &                  *(msqab_cv(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_cv(n,pntr(c(n),a(n)),j,k))*xn/3._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(b(n),c(n)),j,k)*xn/3._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,pntr(a(n),b(n)),j,k)*xn/3._dp
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,6
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/3._dp
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,0,j,k)*xn/3._dp
          enddo
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                   -sub1b_2(n,qq)
     &                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))/xn/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                   -sub2b_1(n,qq)
     &                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))/xn/3._dp
         enddo
         endif
c------ sub-sub-leading colour
         if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                  -sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*(xn+1._dp/xn)/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                  -sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*(xn+1._dp/xn)/3._dp
         enddo
         endif
c--- ANTIQUARK-QUARK contributions
      elseif((j<0).and.(k>0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,6
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *(msqab_c(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_c(n,pntr(c(n),a(n)),j,k))*xn/3._dp
     &                 +subab_cv(n)/2._dp
     &                  *(msqab_cv(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_cv(n,pntr(c(n),a(n)),j,k))*xn/3._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(c(n),b(n)),j,k)*xn/3._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,pntr(b(n),a(n)),j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,pntr(b(n),a(n)),j,k)*xn/3._dp
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,6
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/3._dp
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,0,j,k)*xn/3._dp
          enddo
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                   -sub1b_2(n,qq)
     &                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))/xn/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                   -sub2b_1(n,qq)
     &                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))/xn/3._dp
         enddo
         endif
c------ sub-sub-leading colour
         if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                  -sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*(xn+1._dp/xn)/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                  -sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*(xn+1._dp/xn)/3._dp
         enddo
         endif
c--- GLUON-GLUON contributions
      elseif ((j==0) .and. (k==0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
c--- choose n=2 which is (5,7,6)
          do n=2,2
          msq(18+n,j,k)=sub1b_2(n,gg)
     &                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))*xn
     &                 +sub1b_2v(n)
     &                  *(msq1b_2v(n,1,j,k)+msq1b_2v(n,2,j,k))*xn
          msq(21+n,j,k)=sub2b_1(n,gg)
     &                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))*xn
     &                 +sub2b_1v(n)
     &                  *(msq2b_1v(n,1,j,k)+msq2b_1v(n,2,j,k))*xn
          enddo
c--- choose n=3,6 which is (7,5,6) and (7,6,5)
          do n=3,6,3
          msq(6+n,j,k)=(sub1a_b(n,gg)+subba_1(n,qq))
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,pntr(b(n),c(n)),j,k)*xn
          enddo
c--- choose n=1,5 which is (5,6,7) and (6,5,7)
          do n=1,5,4
          msq(12+n,j,k)=(sub2c_b(n,gg)+subbc_2(n,qq))
     &                   *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn
     &                  +sub2c_bv(n)
     &                   *msq2c_bv(n,pntr(a(n),b(n)),j,k)*xn
         enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))*xn*(avegg/aveqg)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))*xn*(avegg/aveqg)
          enddo
c--- choose n=1 which is (5,6,7)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))*xn*(avegg/aveqg)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))*xn*(avegg/aveqg)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=2,4,2
c-- (a,b,c) = (5,7,6) and (6,7,5)
          msq(n,j,k)   =msq(n,j,k)
     &                   -subab_c(n,qq)
     &                  *(msqab_c(n,1,j,k)+msqab_c(n,2,j,k))/xn
          enddo
          do n=3,6,3
c-- (a,b,c) = (7,5,6) and (7,6,5)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 +(sub1a_b(n,gg)+subba_1(n,qq))
     &                  *msq1a_b(n,0,j,k)*xn
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,0,j,k)*xn
          enddo
c--- choose n=1,5 which is (5,6,7) and (6,5,7)
          do n=1,5,4
          msq(12+n,j,k)=msq(12+n,j,k)
     &                  +(sub2c_b(n,gg)+subbc_2(n,qq))
     &                   *msq2c_b(n,0,j,k)*xn
     &                  +sub2c_bv(n)
     &                   *msq2c_bv(n,0,j,k)*xn
         enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &    (-(msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))/xn
     &     +(msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*2._dp*xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &    (-(msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))/xn
     &     +(msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*2._dp*xn)
          enddo
c--- choose n=1 which is (5,6,7)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &    (-(msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))/xn
     &     +(msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*2._dp*xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &    (-(msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))/xn
     &     +(msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=2,4,2
c-- (a,b,c) = (5,7,6) and (6,7,5)
          msq(n,j,k)   =msq(n,j,k)
     &                  -subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*(xn+1._dp/xn)
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &      (msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*(-xn-1._dp/xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &      (msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*(-xn-1._dp/xn)
          enddo
c--- choose n=1 which is (5,6,7)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &      (msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*(-xn-1._dp/xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &      (msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*(-xn-1._dp/xn)
          enddo
          endif
c--- QUARK-GLUON contributions
      elseif ((j>0) .and. (k==0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,gg))
     &                  *(msq2c_b(n,pntr(a(n),b(n)),j,k)
     &                   +msq2c_b(n,pntr(b(n),a(n)),j,k))*xn/2._dp
     &                 +subbc_2v(n)/2._dp
     &                  *(msqbc_2v(n,pntr(a(n),b(n)),j,k)
     &                   +msqbc_2v(n,pntr(b(n),a(n)),j,k))*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *(msq2c_bv(n,pntr(a(n),b(n)),j,k)
     &                   +msq2c_bv(n,pntr(b(n),a(n)),j,k))*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,qq)
     &                  *msq1b_2(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,gg)
     &                  *msq2b_1(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(b(n),c(n)),j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,pntr(a(n),b(n)),j,k)*xn/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *(msq1a_b(n,1,j,k)+msq1a_b(n,2,j,k))/xn/2._dp
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,0,j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,0,j,k)*xn/2._dp
          enddo
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,gg)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,0,j,k)*xn/2._dp
          enddo
c-- [n=4] (a,b,c) = (6,7,5)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &    (-(msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))/xn
     &     +(msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *msq1a_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &      (msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*(-xn-1._dp/xn)
          enddo
          endif
c--- GLUON-QUARK contributions
      elseif ((j==0) .and. (k>0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(12+n,j,k) =(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                    *msq2c_b(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                   +subbc_2v(n)/2._dp
     &                    *msqbc_2v(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,gg)
     &                  *msq1b_2(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,qq)
     &                  *msq2b_1(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(6+n,j,k)=(subba_1(n,gg)/2._dp+sub1a_b(n,gg))
     &                 *(msq1a_b(n,pntr(c(n),b(n)),j,k)
     &                  +msq1a_b(n,pntr(b(n),c(n)),j,k))*xn/2._dp
     &                +subba_1v(n)/2._dp
     &                 *(msqba_1v(n,pntr(c(n),b(n)),j,k)
     &                  +msqba_1v(n,pntr(b(n),c(n)),j,k))*xn/2._dp
     &                +sub1a_bv(n)
     &                 *(msq1a_bv(n,pntr(c(n),b(n)),j,k)
     &                  +msq1a_bv(n,pntr(b(n),c(n)),j,k))*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,5,7) and (7,5,6)
          msq(6+n,j,k)=(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,pntr(c(n),b(n)),j,k)*xn/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (5,7,6) and (5,6,7)
          msq(12+n,j,k) =msq(12+n,j,k)
     &                  +(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                   *msq2c_b(n,0,j,k)*xn/2._dp
     &                  +subbc_2v(n)/2._dp
     &                   *msqbc_2v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *(msq2c_b(n,1,j,k)+msq2c_b(n,2,j,k))/xn/2._dp
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,0,j,k)*xn/2._dp
          enddo
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,gg)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
          enddo
c-- [n=4] (a,b,c) = (6,7,5)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &    (-(msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))/xn
     &     +(msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &      (msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*(-xn-1._dp/xn)
          enddo
          endif
c--- GLUON-ANTIQUARK contributions
      elseif ((j==0) .and. (k<0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(12+n,j,k) =(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                    *msq2c_b(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                   +subbc_2v(n)/2._dp
     &                    *msqbc_2v(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,gg)
     &                  *msq1b_2(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,qq)
     &                  *msq2b_1(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(6+n,j,k)=(subba_1(n,gg)/2._dp+sub1a_b(n,gg))
     &                 *(msq1a_b(n,pntr(b(n),c(n)),j,k)
     &                  +msq1a_b(n,pntr(c(n),b(n)),j,k))*xn/2._dp
     &                +subba_1v(n)/2._dp
     &                 *(msqba_1v(n,pntr(b(n),c(n)),j,k)
     &                  +msqba_1v(n,pntr(c(n),b(n)),j,k))*xn/2._dp
     &                +sub1a_bv(n)
     &                 *(msq1a_bv(n,pntr(b(n),c(n)),j,k)
     &                  +msq1a_bv(n,pntr(c(n),b(n)),j,k))*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,5,7) and (7,5,6)
          msq(6+n,j,k)=(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,pntr(b(n),c(n)),j,k)*xn/2._dp
          enddo
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (5,7,6) and (5,6,7)
          msq(12+n,j,k) =msq(12+n,j,k)
     &                  +(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                   *msq2c_b(n,0,j,k)*xn/2._dp
     &                  +subbc_2v(n)/2._dp
     &                   *msqbc_2v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *(msq2c_b(n,1,j,k)+msq2c_b(n,2,j,k))/xn/2._dp
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,0,j,k)*xn/2._dp
          enddo
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,gg)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
          enddo
c-- [n=4] (a,b,c) = (6,7,5)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &    (-(msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))/xn
     &     +(msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &      (msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*(-xn-1._dp/xn)
          enddo
          endif
c--- ANTIQUARK-GLUON contributions
      elseif ((j<0) .and. (k==0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,gg))
     &                  *(msq2c_b(n,pntr(b(n),a(n)),j,k)
     &                   +msq2c_b(n,pntr(a(n),b(n)),j,k))*xn/2._dp
     &                 +subbc_2v(n)/2._dp
     &                  *(msqbc_2v(n,pntr(b(n),a(n)),j,k)
     &                   +msqbc_2v(n,pntr(a(n),b(n)),j,k))*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *(msq2c_bv(n,pntr(b(n),a(n)),j,k)
     &                   +msq2c_bv(n,pntr(a(n),b(n)),j,k))*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,qq)
     &                  *msq1b_2(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,gg)
     &                  *msq2b_1(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(c(n),b(n)),j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,pntr(b(n),a(n)),j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,pntr(b(n),a(n)),j,k)*xn/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *(msq1a_b(n,1,j,k)+msq1a_b(n,2,j,k))/xn/2._dp
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,0,j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,0,j,k)*xn/2._dp
          enddo
          do n=1,2
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,gg)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,0,j,k)*xn/2._dp
          enddo
c-- [n=4] (a,b,c) = (6,7,5)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &    (-(msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))/xn
     &     +(msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *msq1a_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &      (msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*(-xn-1._dp/xn)
          enddo
          endif
      endif

      enddo
      enddo

      endif


      if (Qflag) then

c---arguments of dipsx_new:
c---    1  dipole number
c---    2  momentum
c---    3  emitter
c---    4  emitted
c---    5  spectator
c---    6  interference-free subtraction, equiv. to AP kernel for qq,qg
c---    7  correlation piece of subtraction, relevant only for gg,gq
c---    8  lowest order matrix elements at rescaled momentum, msq(j,k)
c---    9  lowest order matrix elements at rescaled momentum
c---        with emitter contracted with appropriate vector, msqv(j,k)
c---   10  appropriate lowest order calculating routine for msq
c---   11  appropriate lowest order calculating routine for msqv
c---   12  4-quark contribution to lowest order matrix elements squared
c---   13  lowest order matrix elements with 4 indices, msqx(j,k,l,m)
c---         Sum_{l,m} msqx(j,k,l,m) = msq(j,k)
c---   14  2-quark contribution to lowest order matrix elements squared,
c---        separated by colours
c---   15  2-quark contribution to lowest order matrix elements squared,
c---        separated by colours, contracted with appropriate vector
c---   16  lowest order matrix elements with 4 indices and
c----        contracted with appropriate vector, msqvx(j,k,l,m)
c---         Sum_{l,m} msqvx(j,k,l,m) = msqv(j,k)

c--- calculate all the dipoles
c--- the dipole number relates the matrix elements to the transformed
c--- momenta, used in realint to perform cuts and clustering

c--- final-final
      call dipsx_new(2,p,5,7,6,sub57_6,sub57_6v,msq57_6,msq57_6v,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m57_6,m57_6x,m57_6g,m57_6vg,m57_6vx)
      call dipsx_new(4,p,6,7,5,sub67_5,sub67_5v,msq67_5,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m67_5,m67_5x,m67_5g,m67_5vg,m67_5vx)
c--- dipole added for consistency with Gflag piece
      call dipsx_new(5,p,6,5,7,sub65_7,sub65_7v,msq65_7,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m65_7,m65_7x,m65_7g,m65_7vg,m65_7vx)

c--- now the basic initial final and final initial
c--- second call for dipole 9,12 etc  only supplies new values for
c--- sub..
      call dipsx_new(7,p,1,5,6,sub15_6,dsubv,msq15_6,dummyv,
     & qqb_z2jetx_new,donothing_gvecx_new,m15_6,m15_6x,m15_6g,mvg,mvxg)
      call dipsx_new(7,p,6,5,1,sub65_1,sub65_1v,dummy,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,mqq,msqx,m65_1g,m65_1vg,m65_1vx)

      call dipsx_new(9,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     & qqb_z2jetx_new,donothing_gvecx_new,m17_5,m17_5x,mg,mvg,mvxg)
      call dipsx_new(9,p,5,7,1,sub57_1,sub57_1v,dummy,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,mqq,msqx,m57_1g,m57_1vg,m57_1vx)

      call dipsx_new(12,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,
     & qqb_z2jetx_new,donothing_gvecx_new,m17_6,m17_6x,mg,mvg,mvxg)
      call dipsx_new(12,p,6,7,1,sub67_1,sub67_1v,dummy,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,mqq,msqx,m67_1g,m67_1vg,m67_1vx)

      call dipsx_new(13,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     & qqb_z2jetx_new,donothing_gvecx_new,m27_6,m27_6x,mg,mvg,mvxg)
      call dipsx_new(13,p,6,7,2,sub67_2,sub67_2v,dummy,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,mqq,msqx,m67_2g,m67_2vg,m67_2vx)

      call dipsx_new(17,p,2,7,5,sub27_5,dsubv,msq27_5,dummyv,
     & qqb_z2jetx_new,donothing_gvecx_new,m27_5,m27_5x,mg,mvg,mvxg)
      call dipsx_new(17,p,5,7,2,sub57_2,sub57_2v,dummy,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,mqq,msqx,m57_2g,m57_2vg,m57_2vx)

      call dipsx_new(18,p,2,5,6,sub25_6,dsubv,msq25_6,dummyv,
     & qqb_z2jetx_new,donothing_gvecx_new,m25_6,m25_6x,m25_6g,mvg,mvxg)
      call dipsx_new(18,p,6,5,2,sub65_2,sub65_2v,dummy,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,mqq,msqx,m65_2g,m65_2vg,m65_2vx)

c--- calculate all the initial-initial dipoles
      call dipsx_new(19,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m16_2,m16_2x,m16_2g,m16_2vg,m16_2vx)
      call dipsx_new(20,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m17_2,m17_2x,m17_2g,m17_2vg,m17_2vx)
      call dipsx_new(21,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m15_2,m15_2x,m15_2g,m15_2vg,m15_2vx)
      call dipsx_new(22,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m26_1,m26_1x,m26_1g,m26_1vg,m26_1vx)
      call dipsx_new(23,p,2,7,1,sub27_1,sub27_1v,msq27_1,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m27_1,m27_1x,m27_1g,m27_1vg,m27_1vx)
      call dipsx_new(24,p,2,5,1,sub25_1,sub25_1v,msq25_1,dummyv,
     & qqb_z2jetx_new,qqb_z2jet_gvecx_new,m25_1,m25_1x,m25_1g,m25_1vg,m25_1vx)

c--- fill the dipole contributions
      do j=-nf,nf
      do k=-nf,nf

      if ((j > 0) .and. (k>0)) then
c---QQ
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)

      if (j==k) then
      msq(22,j,k)=msq(22,j,k)+0.5_dp*sub26_1(gq)*(aveqq/aveqg)
     &    *(m26_1x(0,pp(jj(j),0,jj(j),0))+m26_1x(1,pp(jj(j),0,jj(j),0))+m26_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub26_1v*(aveqq/aveqg)*m26_1vx(pp(jj(j),0,jj(j),0))
      msq(24,j,k)=msq(24,j,k)+0.5_dp*sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub25_1v*(aveqq/aveqg)*(+m25_1vx(pp(jj(j),0,jj(j),0)))
      msq(21,j,k)=msq(21,j,k)+0.5_dp*sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      msq(19,j,k)=msq(19,j,k)+0.5_dp*sub16_2(gq)*(aveqq/aveqg)
     &    *(m16_2x(0,pp(0,kk(k),kk(k),0))+m16_2x(1,pp(0,kk(k),kk(k),0))+m16_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub16_2v*(aveqq/aveqg)*m16_2vx(pp(0,kk(k),kk(k),0))

      elseif (j .ne. k) then
      msq(22,j,k)=msq(22,j,k)+sub26_1(gq)*(aveqq/aveqg)
     &    *(m26_1x(0,pp(jj(j),0,jj(j),0))+m26_1x(1,pp(jj(j),0,jj(j),0))+m26_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub26_1v*(aveqq/aveqg)*m26_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      endif

      elseif ((j < 0) .and. (k<0)) then
c---QbarQbar
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)

      if (j==k) then
      msq(22,j,k)=msq(22,j,k)+0.5_dp*sub26_1(gq)*(aveqq/aveqg)
     &    *(m26_1x(0,pp(jj(j),0,jj(j),0))+m26_1x(1,pp(jj(j),0,jj(j),0))+m26_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub26_1v*(aveqq/aveqg)*m26_1vx(pp(jj(j),0,jj(j),0))
      msq(24,j,k)=msq(24,j,k)+0.5_dp*sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub25_1v*(aveqq/aveqg)*(+m25_1vx(pp(jj(j),0,jj(j),0)))
      msq(21,j,k)=msq(21,j,k)+0.5_dp*sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      msq(19,j,k)=msq(19,j,k)+0.5_dp*sub16_2(gq)*(aveqq/aveqg)
     &    *(m16_2x(0,pp(0,kk(k),kk(k),0))+m16_2x(1,pp(0,kk(k),kk(k),0))+m16_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub16_2v*(aveqq/aveqg)*m16_2vx(pp(0,kk(k),kk(k),0))

      elseif (j .ne. k) then
      msq(22,j,k)=msq(22,j,k)+sub26_1(gq)*(aveqq/aveqg)
     &    *(m26_1x(0,pp(jj(j),0,jj(j),0))+m26_1x(1,pp(jj(j),0,jj(j),0))+m26_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub26_1v*(aveqq/aveqg)*m26_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      endif


      elseif ((j > 0) .and. (k<0)) then
c---QQbar
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(2,j,k)-m57_6(0,j,k)/xn-m57_6(1,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(2,j,k)-m67_5(0,j,k)/xn-m67_5(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
     & +sub67_1(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
     & +sub57_2(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(2,j,k)-m17_2(0,j,k)/xn-m17_2(1,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(2,j,k)-m27_1(0,j,k)/xn-m27_1(1,j,k)/xn)


      if (j .ne. -k) then
      msq(22,j,k)=msq(22,j,k)+sub26_1(gq)*(aveqq/aveqg)
     &    *(m26_1x(0,pp(jj(j),0,jj(j),0))+m26_1x(1,pp(jj(j),0,jj(j),0))+m26_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub26_1v*(aveqq/aveqg)*m26_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      elseif (j==-k) then
c--- combination of final-initial and final-final for consistency with Gflag
      msq( 5,j,k)=msq( 5,j,k)+real(nf,dp)
     & *(sub65_7(gq)*(m65_7g(0,j,k)+m65_7g(1,j,k)+m65_7g(2,j,k))
     &  -sub65_7v*(m65_7vg(0,j,k)+m65_7vg(1,j,k)+m65_7vg(2,j,k)))
c--- 8/6/09: symmetrized 65 between _1 and _2
      msq( 7,j,k)=msq( 7,j,k)+half*real(nf,dp)
     & *(sub65_1(gq)*(m15_6g(0,j,k)+m15_6g(1,j,k)+m15_6g(2,j,k))
     &  -sub65_1v*(m65_1vg(0,j,k)+m65_1vg(1,j,k)+m65_1vg(2,j,k)))
      msq(18,j,k)=msq(18,j,k)+half*real(nf,dp)
     & *(sub65_2(gq)*(m25_6g(0,j,k)+m25_6g(1,j,k)+m25_6g(2,j,k))
     &  -sub65_2v*(m65_2vg(0,j,k)+m65_2vg(1,j,k)+m65_2vg(2,j,k)))

      msq(22,j,k)=msq(22,j,k)+sub26_1(gq)*(aveqq/aveqg)
     &    *(m26_1x(0,pp(jj(j),0,jj(j),0))+m26_1x(1,pp(jj(j),0,jj(j),0))+m26_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub26_1v*(aveqq/aveqg)*m26_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))

      endif

      elseif ((j < 0) .and. (k>0)) then
c---QbarQ
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(2,j,k)-m57_6(0,j,k)/xn-m57_6(1,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(2,j,k)-m67_5(0,j,k)/xn-m67_5(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
     & +sub57_1(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
     & +sub67_2(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(2,j,k)-m17_2(0,j,k)/xn-m17_2(1,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(2,j,k)-m27_1(0,j,k)/xn-m27_1(1,j,k)/xn)

      if (-j .ne. k) then
      msq(24,j,k)=msq(24,j,k)+sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(19,j,k)=msq(19,j,k)+sub16_2(gq)*(aveqq/aveqg)
     &    *(m16_2x(0,pp(0,kk(k),kk(k),0))+m16_2x(1,pp(0,kk(k),kk(k),0))+m16_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub16_2v*(aveqq/aveqg)*m16_2vx(pp(0,kk(k),kk(k),0))
      elseif (-j==k) then
c--- combination of final-initial and final-final for consistency with Gflag
      msq( 5,j,k)=msq( 5,j,k)+real(nf,dp)
     & *(sub65_7(gq)*(m65_7g(0,j,k)+m65_7g(1,j,k)+m65_7g(2,j,k))
     &  -sub65_7v*(m65_7vg(0,j,k)+m65_7vg(1,j,k)+m65_7vg(2,j,k)))
c--- 8/6/09: symmetrized 65 between _1 and _2
      msq( 7,j,k)=msq( 7,j,k)+half*real(nf,dp)
     & *(sub65_1(gq)*(m15_6g(0,j,k)+m15_6g(1,j,k)+m15_6g(2,j,k))
     &  -sub65_1v*(m65_1vg(0,j,k)+m65_1vg(1,j,k)+m65_1vg(2,j,k)))
      msq(18,j,k)=msq(18,j,k)+half*real(nf,dp)
     & *(sub65_2(gq)*(m25_6g(0,j,k)+m25_6g(1,j,k)+m25_6g(2,j,k))
     &  -sub65_2v*(m65_2vg(0,j,k)+m65_2vg(1,j,k)+m65_2vg(2,j,k)))

      msq(24,j,k)=msq(24,j,k)+sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(19,j,k)=msq(19,j,k)+sub16_2(gq)*(aveqq/aveqg)
     &    *(m16_2x(0,pp(0,kk(k),kk(k),0))+m16_2x(1,pp(0,kk(k),kk(k),0))+m16_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub16_2v*(aveqq/aveqg)*m16_2vx(pp(0,kk(k),kk(k),0))
      endif

      elseif (j==0) then
c---------G-Q and G-Qbar
        if    (k > 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_1(gq)
     &    *(m15_6x(1,pp(0,kk(k),0,kk(k)))+m15_6x(2,pp(0,kk(k),0,kk(k))))
     &    -sub65_1v
     &    *(m65_1vg(1,0,kk(k))+m65_1vg(2,0,kk(k))))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_2(gq)
     &    *(m25_6x(0,pp(0,kk(k),0,kk(k)))+m25_6x(2,pp(0,kk(k),0,kk(k))))
     &    -sub65_2v
     &    *(m65_2vg(0,0,kk(k))+m65_2vg(2,0,kk(k))))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_7(gq)
     &    *(m65_7x(0,pp(0,kk(k),0,kk(k)))+m65_7x(1,pp(0,kk(k),0,kk(k))))
     &    -sub65_7v
     &    *(m65_7vg(0,0,kk(k))+m65_7vg(1,0,kk(k))))

          msq(12,j,k)=msq(12,j,k)+0.5_dp*half*sub67_1(gq)
     &    *(m17_6x(1,pp(0,kk(k),kk(k),0))+m17_6x(2,pp(0,kk(k),kk(k),0)))
     &                           -0.5_dp*half*sub67_1v
     &    *(m67_1vg(1,0,kk(k))+m67_1vg(2,0,kk(k)))
          msq(13,j,k)=msq(13,j,k)+0.5_dp*half*sub67_2(gq)
     &    *(m27_6x(0,pp(0,kk(k),kk(k),0))+m27_6x(2,pp(0,kk(k),kk(k),0)))
     &                           -0.5_dp*half*sub67_2v
     &    *(m67_2vg(0,0,kk(k))+m67_2vg(2,0,kk(k)))
          msq( 4,j,k)=msq( 4,j,k)+0.5_dp*half*sub67_5(gq)
     &    *(m67_5x(0,pp(0,kk(k),kk(k),0))+m67_5x(1,pp(0,kk(k),kk(k),0)))
     &                           -0.5_dp*half*sub67_5v
     &    *(m67_5vg(0,0,kk(k))+m67_5vg(1,0,kk(k)))

          if (kk(k) == 1) then
          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq*(
     &     (m16_2x(0,pp(1,1,1,1))+m16_2x(1,pp(1,1,1,1))+m16_2x(2,pp(1,1,1,1)))
     &    +2._dp*(m16_2x(0,pp(3,1,3,1))+m16_2x(1,pp(3,1,3,1))+m16_2x(2,pp(3,1,3,1)))
     &    +2._dp*(m16_2x(0,pp(2,1,2,1))+m16_2x(1,pp(2,1,2,1))+m16_2x(2,pp(2,1,2,1))))

          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub17_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m17_2x(0,pp(-1,1,1,-1))+m17_2x(1,pp(-1,1,1,-1))+m17_2x(2,pp(-1,1,1,-1)))
     &   +2._dp
     &   *(m17_2x(0,pp(-1,1,2,-2))+m17_2x(1,pp(-1,1,2,-2))+m17_2x(2,pp(-1,1,2,-2)))
     &   +2._dp
     &   *(m17_2x(0,pp(-1,1,3,-3))+m17_2x(1,pp(-1,1,3,-3))+m17_2x(2,pp(-1,1,3,-3))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m15_2x(0,pp(-1,1,-1,1))+m15_2x(1,pp(-1,1,-1,1))+m15_2x(2,pp(-1,1,-1,1)))
     &   +2._dp
     &   *(m15_2x(0,pp(-2,1,-2,1))+m15_2x(1,pp(-2,1,-2,1))+m15_2x(2,pp(-2,1,-2,1)))
     &   +2._dp
     &   *(m15_2x(0,pp(-3,1,-3,1))+m15_2x(1,pp(-3,1,-3,1))+m15_2x(2,pp(-3,1,-3,1))))

          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,1,-1))+m27_1x(1,pp(0,0,1,-1))+m27_1x(2,pp(0,0,1,-1)))
     &    +sub27_1v*m27_1vx(pp(0,0,1,-1)))
     &    +2._dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,2,-2))+m27_1x(1,pp(0,0,2,-2))+m27_1x(2,pp(0,0,2,-2)))
     &    +sub27_1v*m27_1vx(pp(0,0,2,-2))))

          msq(24,j,k)=msq(24,j,k)+0.5_dp*(aveqg/avegg)*(sub25_1(gq)
     &    *(m25_1x(0,pp(0,0,-1,1))+m25_1x(1,pp(0,0,-1,1))+m25_1x(2,pp(0,0,-1,1)))
     &                     +sub25_1v*m25_1vx(pp(0,0,-1,1)))

         endif
         if (kk(k) == 2) then
          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,1,-1))+m27_1x(1,pp(0,0,1,-1))+m27_1x(2,pp(0,0,1,-1)))
     &    +sub27_1v*m27_1vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,2,-2))+m27_1x(1,pp(0,0,2,-2))+m27_1x(2,pp(0,0,2,-2)))
     &    +sub27_1v*m27_1vx(pp(0,0,2,-2))))

          msq(24,j,k)=msq(24,j,k)+0.5_dp*(aveqg/avegg)*(sub25_1(gq)
     &    *(m25_1x(0,pp(0,0,-2,2))+m25_1x(1,pp(0,0,-2,2))+m25_1x(2,pp(0,0,-2,2)))
     &                     +sub25_1v*m25_1vx(pp(0,0,-2,2)))

          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq*(
     &     3._dp*(m16_2x(0,pp(1,2,1,2))+m16_2x(1,pp(1,2,1,2))+m16_2x(2,pp(1,2,1,2)))
     &    +(m16_2x(0,pp(4,2,4,2))+m16_2x(1,pp(4,2,4,2))+m16_2x(2,pp(4,2,4,2)))
     &    +(m16_2x(0,pp(2,2,2,2))+m16_2x(1,pp(2,2,2,2))+m16_2x(2,pp(2,2,2,2))))

          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub17_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m17_2x(0,pp(-2,2,2,-2))+m17_2x(1,pp(-2,2,2,-2))+m17_2x(2,pp(-2,2,2,-2)))
     &   +3._dp
     &   *(m17_2x(0,pp(-2,2,1,-1))+m17_2x(1,pp(-2,2,1,-1))+m17_2x(2,pp(-2,2,1,-1)))
     &   +1._dp
     &   *(m17_2x(0,pp(-2,2,4,-4))+m17_2x(1,pp(-2,2,4,-4))+m17_2x(2,pp(-2,2,4,-4))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m15_2x(0,pp(-1,2,-1,2))+m15_2x(1,pp(-1,2,-1,2))+m15_2x(2,pp(-1,2,-1,2)))
     &   +0.5_dp
     &   *(m15_2x(0,pp(-2,2,-2,2))+m15_2x(1,pp(-2,2,-2,2))+m15_2x(2,pp(-2,2,-2,2)))
     &   +1._dp
     &   *(m15_2x(0,pp(-4,2,-4,2))+m15_2x(1,pp(-4,2,-4,2))+m15_2x(2,pp(-4,2,-4,2))))
         endif
        elseif (k < 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_1(gq)
     &    *(m15_6x(1,pp(0,kk(k),0,kk(k)))+m15_6x(2,pp(0,kk(k),0,kk(k))))
     &    -sub65_1v
     &    *(m65_1vg(1,0,kk(k))+m65_1vg(2,0,kk(k))))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_2(gq)
     &    *(m25_6x(0,pp(0,kk(k),0,kk(k)))+m25_6x(1,pp(0,kk(k),0,kk(k))))
     &    -sub65_2v
     &    *(m65_2vg(0,0,kk(k))+m65_2vg(1,0,kk(k))))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_7(gq)
     &    *(m65_7x(0,pp(0,kk(k),0,kk(k)))+m65_7x(2,pp(0,kk(k),0,kk(k))))
     &    -sub65_7v
     &    *(m65_7vg(0,0,kk(k))+m65_7vg(2,0,kk(k))))

          msq( 9,j,k)=msq( 9,j,k)+0.5_dp*half*sub57_1(gq)
     &    *(m17_5x(1,pp(0,kk(k),0,kk(k)))+m17_5x(2,pp(0,kk(k),0,kk(k))))
     &                           -0.5_dp*half*sub57_1v
     &    *(m57_1vg(1,0,kk(k))+m57_1vg(2,0,kk(k)))
          msq(17,j,k)=msq(17,j,k)+0.5_dp*half*sub57_2(gq)
     &    *(m27_5x(0,pp(0,kk(k),0,kk(k)))+m27_5x(1,pp(0,kk(k),0,kk(k))))
     &                           -0.5_dp*half*sub57_2v
     &    *(m57_2vg(0,0,kk(k))+m57_2vg(1,0,kk(k)))
          msq( 2,j,k)=msq( 2,j,k)+0.5_dp*half*sub57_6(gq)
     &    *(m57_6x(0,pp(0,kk(k),0,kk(k)))+m57_6x(2,pp(0,kk(k),0,kk(k))))
     &                           -0.5_dp*half*sub57_6v
     &    *(m57_6vg(0,0,kk(k))+m57_6vg(2,0,kk(k)))

          if (kk(k) == -1) then
          msq(22,j,k)=msq(22,j,k)+0.5_dp*(aveqg/avegg)*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,1,-1))+m26_1x(1,pp(0,0,1,-1))+m26_1x(2,pp(0,0,1,-1)))
     &                     +sub26_1v*m26_1vx(pp(0,0,1,-1)))

          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,1,-1))+m27_1x(1,pp(0,0,1,-1))+m27_1x(2,pp(0,0,1,-1)))
     &    +sub27_1v*m27_1vx(pp(0,0,1,-1)))
     &    +2._dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,2,-2))+m27_1x(1,pp(0,0,2,-2))+m27_1x(2,pp(0,0,2,-2)))
     &    +sub27_1v*m27_1vx(pp(0,0,2,-2))))

          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m16_2x(0,pp(1,-1,1,-1))+m16_2x(1,pp(1,-1,1,-1))+m16_2x(2,pp(1,-1,1,-1)))
     &   +2._dp
     &   *(m16_2x(0,pp(2,-1,2,-1))+m16_2x(1,pp(2,-1,2,-1))+m16_2x(2,pp(2,-1,2,-1)))
     &   +2._dp
     &   *(m16_2x(0,pp(3,-1,3,-1))+m16_2x(1,pp(3,-1,3,-1))+m16_2x(2,pp(3,-1,3,-1))))

          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub17_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m17_2x(0,pp(1,-1,1,-1))+m17_2x(1,pp(1,-1,1,-1))+m17_2x(2,pp(1,-1,1,-1)))
     &   +2._dp
     &   *(m17_2x(0,pp(1,-1,2,-2))+m17_2x(1,pp(1,-1,2,-2))+m17_2x(2,pp(1,-1,2,-2)))
     &   +2._dp
     &   *(m17_2x(0,pp(1,-1,3,-3))+m17_2x(1,pp(1,-1,3,-3))+m17_2x(2,pp(1,-1,3,-3))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq*(
     &     (m15_2x(0,pp(-1,-1,-1,-1))+m15_2x(1,pp(-1,-1,-1,-1))
     &     +m15_2x(2,pp(-1,-1,-1,-1)))
     &    +2._dp*(m15_2x(0,pp(-3,-1,-3,-1))+m15_2x(1,pp(-3,-1,-3,-1))
     &     +m15_2x(2,pp(-3,-1,-3,-1)))
     &    +2._dp*(m15_2x(0,pp(-2,-1,-2,-1))+m15_2x(1,pp(-2,-1,-2,-1))
     &     +m15_2x(2,pp(-2,-1,-2,-1))))

          endif
          if (kk(k) == -2) then
          msq(22,j,k)=msq(22,j,k)+0.5_dp*(aveqg/avegg)*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,2,-2))+m26_1x(1,pp(0,0,2,-2))+m26_1x(2,pp(0,0,2,-2)))
     &                     +sub26_1v*m26_1vx(pp(0,0,2,-2)))

          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,1,-1))+m27_1x(1,pp(0,0,1,-1))+m27_1x(2,pp(0,0,1,-1)))
     &    +sub27_1v*m27_1vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub27_1(gq)
     &    *(m27_1x(0,pp(0,0,2,-2))+m27_1x(1,pp(0,0,2,-2))+m27_1x(2,pp(0,0,2,-2)))
     &    +sub27_1v*m27_1vx(pp(0,0,2,-2))))

          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m16_2x(0,pp(1,-2,1,-2))+m16_2x(1,pp(1,-2,1,-2))+m16_2x(2,pp(1,-2,1,-2)))
     &   +0.5_dp
     &    *(m16_2x(0,pp(2,-2,2,-2))+m16_2x(1,pp(2,-2,2,-2))+m16_2x(2,pp(2,-2,2,-2)))
     &   +1._dp
     &   *(m16_2x(0,pp(4,-2,4,-2))+m16_2x(1,pp(4,-2,4,-2))+m16_2x(2,pp(4,-2,4,-2))))

          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub17_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m17_2x(0,pp(2,-2,2,-2))+m17_2x(1,pp(2,-2,2,-2))+m17_2x(2,pp(2,-2,2,-2)))
     &   +3._dp
     &   *(m17_2x(0,pp(2,-2,1,-1))+m17_2x(1,pp(2,-2,1,-1))+m17_2x(2,pp(2,-2,1,-1)))
     &   +1._dp
     &   *(m17_2x(0,pp(2,-2,4,-4))+m17_2x(1,pp(2,-2,4,-4))+m17_2x(2,pp(2,-2,4,-4))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq*(
     &     3._dp*(m15_2x(0,pp(-1,-2,-1,-2))+m15_2x(1,pp(-1,-2,-1,-2))
     &     +m15_2x(2,pp(-1,-2,-1,-2)))
     &    +(m15_2x(0,pp(-4,-2,-4,-2))+m15_2x(1,pp(-4,-2,-4,-2))
     &     +m15_2x(2,pp(-4,-2,-4,-2)))
     &    +(m15_2x(0,pp(-2,-2,-2,-2))+m15_2x(1,pp(-2,-2,-2,-2))
     &     +m15_2x(2,pp(-2,-2,-2,-2))))

          endif
        endif
      elseif (k==0) then
c---------Q-G and Qbar-G
        if     (j > 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_1(gq)
     &    *(m15_6x(0,pp(jj(j),0,0,jj(j)))+m15_6x(2,pp(jj(j),0,0,jj(j))))
     &    -sub65_1v
     &    *(m65_1vg(0,jj(j),0)+m65_1vg(2,jj(j),0)))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_2(gq)
     &    *(m25_6x(1,pp(jj(j),0,0,jj(j)))+m25_6x(2,pp(jj(j),0,0,jj(j))))
     &    -sub65_2v
     &    *(m65_2vg(1,jj(j),0)+m65_2vg(2,jj(j),0)))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_7(gq)
     &    *(m65_7x(0,pp(jj(j),0,0,jj(j)))+m65_7x(1,pp(jj(j),0,0,jj(j))))
     &    -sub65_7v
     &    *(m65_7vg(0,jj(j),0)+m65_7vg(1,jj(j),0)))

          msq(12,j,k)=msq(12,j,k)+0.5_dp*half*sub67_1(gq)
     &    *(m17_6x(0,pp(jj(j),0,jj(j),0))+m17_6x(2,pp(jj(j),0,jj(j),0)))
     &                           -0.5_dp*half*sub67_1v
     &    *(m67_1vg(0,jj(j),0)+m67_1vg(2,jj(j),0))
          msq(13,j,k)=msq(13,j,k)+0.5_dp*half*sub67_2(gq)
     &    *(m27_6x(1,pp(jj(j),0,jj(j),0))+m27_6x(2,pp(jj(j),0,jj(j),0)))
     &                           -0.5_dp*half*sub67_2v
     &    *(m67_2vg(1,jj(j),0)+m67_2vg(2,jj(j),0))
          msq( 4,j,k)=msq( 4,j,k)+0.5_dp*half*sub67_5(gq)
     &    *(m67_5x(0,pp(jj(j),0,jj(j),0))+m67_5x(1,pp(jj(j),0,jj(j),0)))
     &                           -0.5_dp*half*sub67_5v
     &    *(m67_5vg(0,jj(j),0)+m67_5vg(1,jj(j),0))

          if (jj(j) == 1) then
          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,1,-1))+m17_2x(1,pp(0,0,1,-1))+m17_2x(2,pp(0,0,1,-1)))
     &    +sub17_2v*m17_2vx(pp(0,0,1,-1)))
     &    +2._dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,2,-2))+m17_2x(1,pp(0,0,2,-2))+m17_2x(2,pp(0,0,2,-2)))
     &    +sub17_2v*m17_2vx(pp(0,0,2,-2))))

          msq(21,j,k)=msq(21,j,k)+0.5_dp*(aveqg/avegg)*(sub15_2(gq)
     &    *(m15_2x(0,pp(0,0,-1,1))+m15_2x(1,pp(0,0,-1,1))+m15_2x(2,pp(0,0,-1,1)))
     &                     +sub15_2v*m15_2vx(pp(0,0,-1,1)))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq*(
     &     (m26_1x(0,pp(1,1,1,1))+m26_1x(1,pp(1,1,1,1))+m26_1x(2,pp(1,1,1,1)))
     &    +2._dp*(m26_1x(0,pp(1,3,3,1))+m26_1x(1,pp(1,3,3,1))+m26_1x(2,pp(1,3,3,1)))
     &    +2._dp*(m26_1x(0,pp(1,2,2,1))+m26_1x(1,pp(1,2,2,1))+m26_1x(2,pp(1,2,2,1))))

          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub27_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m27_1x(0,pp(1,-1,1,-1))+m27_1x(1,pp(1,-1,1,-1))+m27_1x(2,pp(1,-1,1,-1)))
     &   +2._dp
     &   *(m27_1x(0,pp(1,-1,2,-2))+m27_1x(1,pp(1,-1,2,-2))+m27_1x(2,pp(1,-1,2,-2)))
     &   +2._dp
     &   *(m27_1x(0,pp(1,-1,3,-3))+m27_1x(1,pp(1,-1,3,-3))+m27_1x(2,pp(1,-1,3,-3))))


          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m25_1x(0,pp(1,-1,-1,1))+m25_1x(1,pp(1,-1,-1,1))+m25_1x(2,pp(1,-1,-1,1)))
     &   +2._dp
     &   *(m25_1x(0,pp(1,-2,-2,1))+m25_1x(1,pp(1,-2,-2,1))+m25_1x(2,pp(1,-2,-2,1)))
     &   +2._dp
     &   *(m25_1x(0,pp(1,-3,-3,1))+m25_1x(1,pp(1,-3,-3,1))+m25_1x(2,pp(1,-3,-3,1))))

           endif
          if (jj(j) == 2) then
          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,1,-1))+m17_2x(1,pp(0,0,1,-1))+m17_2x(2,pp(0,0,1,-1)))
     &    +sub17_2v*m17_2vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,2,-2))+m17_2x(1,pp(0,0,2,-2))+m17_2x(2,pp(0,0,2,-2)))
     &    +sub17_2v*m17_2vx(pp(0,0,2,-2))))

          msq(21,j,k)=msq(21,j,k)+0.5_dp*(aveqg/avegg)*(sub15_2(gq)
     &    *(m15_2x(0,pp(0,0,-2,2))+m15_2x(1,pp(0,0,-2,2))+m15_2x(2,pp(0,0,-2,2)))
     &                     +sub15_2v*m15_2vx(pp(0,0,-2,2)))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq*(
     &     3._dp*(m26_1x(0,pp(2,1,1,2))+m26_1x(1,pp(2,1,1,2))+m26_1x(2,pp(2,1,1,2)))
     &    +(m26_1x(0,pp(2,4,4,2))+m26_1x(1,pp(2,4,4,2))+m26_1x(2,pp(2,4,4,2)))
     &    +(m26_1x(0,pp(2,2,2,2))+m26_1x(1,pp(2,2,2,2))+m26_1x(2,pp(2,2,2,2))))

          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub27_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m27_1x(0,pp(2,-2,2,-2))+m27_1x(1,pp(2,-2,2,-2))+m27_1x(2,pp(2,-2,2,-2)))
     &   +3._dp
     &   *(m27_1x(0,pp(2,-2,1,-1))+m27_1x(1,pp(2,-2,1,-1))+m27_1x(2,pp(2,-2,1,-1)))
     &   +1._dp
     &   *(m27_1x(0,pp(2,-2,4,-4))+m27_1x(1,pp(2,-2,4,-4))+m27_1x(2,pp(2,-2,4,-4))))

          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m25_1x(0,pp(2,-1,-1,2))+m25_1x(1,pp(2,-1,-1,2))+m25_1x(2,pp(2,-1,-1,2)))
     &   +0.5_dp
     &   *(m25_1x(0,pp(2,-2,-2,2))+m25_1x(1,pp(2,-2,-2,2))+m25_1x(2,pp(2,-2,-2,2)))
     &   +1._dp
     &   *(m25_1x(0,pp(2,-4,-4,2))+m25_1x(1,pp(2,-4,-4,2))+m25_1x(2,pp(2,-4,-4,2))))

          endif
        elseif (j < 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_1(gq)
     &    *(m15_6x(0,pp(jj(j),0,0,jj(j)))+m15_6x(1,pp(jj(j),0,0,jj(j))))
     &    -sub65_1v
     &    *(m65_1vg(0,jj(j),0)+m65_1vg(1,jj(j),0)))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_2(gq)
     &    *(m25_6x(1,pp(jj(j),0,0,jj(j)))+m25_6x(2,pp(jj(j),0,0,jj(j))))
     &    -sub65_2v
     &    *(m65_2vg(1,jj(j),0)+m65_2vg(2,jj(j),0)))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub65_7(gq)
     &    *(m65_7x(0,pp(jj(j),0,0,jj(j)))+m65_7x(2,pp(jj(j),0,0,jj(j))))
     &    -sub65_7v
     &    *(m65_7vg(0,jj(j),0)+m65_7vg(2,jj(j),0)))

          msq( 9,j,k)=msq( 9,j,k)+0.5_dp*half*sub57_1(gq)
     &    *(m17_5x(0,pp(jj(j),0,0,jj(j)))+m17_5x(1,pp(jj(j),0,0,jj(j))))
     &                           -0.5_dp*half*sub57_1v
     &    *(m57_1vg(0,jj(j),0)+m57_1vg(1,jj(j),0))
          msq(17,j,k)=msq(17,j,k)+0.5_dp*half*sub57_2(gq)
     &    *(m27_5x(1,pp(jj(j),0,0,jj(j)))+m27_5x(2,pp(jj(j),0,0,jj(j))))
     &                           -0.5_dp*half*sub57_2v
     &    *(m57_2vg(1,jj(j),0)+m57_2vg(2,jj(j),0))
          msq( 2,j,k)=msq( 2,j,k)+0.5_dp*half*sub57_6(gq)
     &    *(m57_6x(0,pp(jj(j),0,0,jj(j)))+m57_6x(2,pp(jj(j),0,0,jj(j))))
     &                           -0.5_dp*half*sub57_6v
     &    *(m57_6vg(0,jj(j),0)+m57_6vg(2,jj(j),0))

          if (jj(j) == -1) then
          msq(19,j,k)=msq(19,j,k)+0.5_dp*(aveqg/avegg)*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,1,-1))+m16_2x(1,pp(0,0,1,-1))+m16_2x(2,pp(0,0,1,-1)))
     &                     +sub16_2v*m16_2vx(pp(0,0,1,-1)))

          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,1,-1))+m17_2x(1,pp(0,0,1,-1))+m17_2x(2,pp(0,0,1,-1)))
     &    +sub17_2v*m17_2vx(pp(0,0,1,-1)))
     &    +2._dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,2,-2))+m17_2x(1,pp(0,0,2,-2))+m17_2x(2,pp(0,0,2,-2)))
     &    +sub17_2v*m17_2vx(pp(0,0,2,-2))))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m26_1x(0,pp(-1,1,1,-1))+m26_1x(1,pp(-1,1,1,-1))+m26_1x(2,pp(-1,1,1,-1)))
     &   +2._dp
     &   *(m26_1x(0,pp(-1,2,2,-1))+m26_1x(1,pp(-1,2,2,-1))+m26_1x(2,pp(-1,2,2,-1)))
     &   +2._dp
     &   *(m26_1x(0,pp(-1,3,3,-1))+m26_1x(1,pp(-1,3,3,-1))+m26_1x(2,pp(-1,3,3,-1))))

          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub27_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m27_1x(0,pp(-1,1,1,-1))+m27_1x(1,pp(-1,1,1,-1))+m27_1x(2,pp(-1,1,1,-1)))
     &   +2._dp
     &   *(m27_1x(0,pp(-1,1,2,-2))+m27_1x(1,pp(-1,1,2,-2))+m27_1x(2,pp(-1,1,2,-2)))
     &   +2._dp
     &   *(m27_1x(0,pp(-1,1,3,-3))+m27_1x(1,pp(-1,1,3,-3))+m27_1x(2,pp(-1,1,3,-3))))

          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq*(
     &     (m25_1x(0,pp(-1,-1,-1,-1))+m25_1x(1,pp(-1,-1,-1,-1))
     &     +m25_1x(2,pp(-1,-1,-1,-1)))
     &    +2._dp*(m25_1x(0,pp(-1,-3,-3,-1))+m25_1x(1,pp(-1,-3,-3,-1))
     &     +m25_1x(2,pp(-1,-3,-3,-1)))
     &    +2._dp*(m25_1x(0,pp(-1,-2,-2,-1))+m25_1x(1,pp(-1,-2,-2,-1))
     &     +m25_1x(2,pp(-1,-2,-2,-1))))

           endif
          if (jj(j) == -2) then
          msq(19,j,k)=msq(19,j,k)+0.5_dp*(aveqg/avegg)*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,2,-2))+m16_2x(1,pp(0,0,2,-2))+m16_2x(2,pp(0,0,2,-2)))
     &                     +sub16_2v*m16_2vx(pp(0,0,2,-2)))

          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,1,-1))+m17_2x(1,pp(0,0,1,-1))+m17_2x(2,pp(0,0,1,-1)))
     &    +sub17_2v*m17_2vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub17_2(gq)
     &    *(m17_2x(0,pp(0,0,2,-2))+m17_2x(1,pp(0,0,2,-2))+m17_2x(2,pp(0,0,2,-2)))
     &    +sub17_2v*m17_2vx(pp(0,0,2,-2))))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m26_1x(0,pp(-2,1,1,-2))+m26_1x(1,pp(-2,1,1,-2))+m26_1x(2,pp(-2,1,1,-2)))
     &   +0.5_dp
     &    *(m26_1x(0,pp(-2,2,2,-2))+m26_1x(1,pp(-2,2,2,-2))+m26_1x(2,pp(-2,2,2,-2)))
     &   +1._dp
     &   *(m26_1x(0,pp(-2,4,4,-2))+m26_1x(1,pp(-2,4,4,-2))+m26_1x(2,pp(-2,4,4,-2))))

          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub27_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m27_1x(0,pp(-2,2,2,-2))+m27_1x(1,pp(-2,2,2,-2))+m27_1x(2,pp(-2,2,2,-2)))
     &   +3._dp
     &   *(m27_1x(0,pp(-2,2,1,-1))+m27_1x(1,pp(-2,2,1,-1))+m27_1x(2,pp(-2,2,1,-1)))
     &   +1._dp
     &   *(m27_1x(0,pp(-2,2,4,-4))+m27_1x(1,pp(-2,2,4,-4))+m27_1x(2,pp(-2,2,4,-4))))

          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq*(
     &     3._dp*(m25_1x(0,pp(-2,-1,-1,-2))+m25_1x(1,pp(-2,-1,-1,-2))
     &     +m25_1x(2,pp(-2,-1,-1,-2)))
     &    +(m25_1x(0,pp(-2,-4,-4,-2))+m25_1x(1,pp(-2,-4,-4,-2))
     &     +m25_1x(2,pp(-2,-4,-4,-2)))
     &    +(m25_1x(0,pp(-2,-2,-2,-2))+m25_1x(1,pp(-2,-2,-2,-2))
     &     +m25_1x(2,pp(-2,-2,-2,-2))))

          endif
        endif
      endif


      enddo
      enddo

      endif

      return
      end


      subroutine donothing_gvecx_new(p,n,in,msq,mvg,mvxg)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ppmax.f'
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4)
      real(dp):: mvg(0:2,-nf:nf,-nf:nf)
      real(dp):: mvxg(ppmax)

      msq(:,:)=zip

      return
      end

