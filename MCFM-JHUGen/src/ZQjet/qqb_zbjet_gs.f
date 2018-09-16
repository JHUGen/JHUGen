      subroutine qqb_zbjet_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J. Campbell                                              *
*     July, 2005.                                                      *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z + b(p5) + g(p6) + g(p7)
c                          |
c                          --> l(p3)+a(p4)
c
c--- all momenta are incoming
c--- Extended to include charm quark production via the variable "flav"
c
c--- isub=1 corresponds to p7 representing a light quark of gluon
c--- isub=2 means that p7 is another heavy quark


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'nflav.f'
      include 'heavyflav.f'

c--- np6,np12,... = n+6,n+12,...
c -- nd counts the dipoles
      integer:: i
c      real(dp):: ran,c15,c15v,c25,c25v,c16,c16v,c26,c26v,
c     & c57,c57v,c67,c67v,
c     & c175,c571,c175v,c571v,
c     & c176,c671,c176v,c671v,
c     & c275,c572,c275v,c572v,
c     & c276,c672,c276v,c672v
      integer:: j,k,n,np6,np12,np18,np21,nd
c--- slightly obtuse notation, fn=-nf, to simplify declaration lines
      real(dp):: p(mxpart,4),msq(maxd,fn:nf,fn:nf)
      real(dp)::
     & msq17_2(fn:nf,fn:nf),msq27_1(fn:nf,fn:nf),
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq16_2(fn:nf,fn:nf),msq26_1(fn:nf,fn:nf),
     & msq17_6(fn:nf,fn:nf),
     & msq26_1v(fn:nf,fn:nf),
     & msq16_2v(fn:nf,fn:nf),
     & msq67_1v(fn:nf,fn:nf),
     & dummy(fn:nf,fn:nf),dummyv(fn:nf,fn:nf),
     & sub17_2(4),sub27_1(4),sub15_2(4),sub25_1(4),
     & sub16_2(4),sub26_1(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),
     & sub57_6(4),sub67_5(4),dsubv,dsub(4),
     & sub67_1v,sub26_1v,sub16_2v
      real(dp)::
     & m17_2(0:2,fn:nf,fn:nf),m27_1(0:2,fn:nf,fn:nf),
     & m57_6(0:2,fn:nf,fn:nf),m67_5(0:2,fn:nf,fn:nf),
     & m17_5(0:2,fn:nf,fn:nf),
     & m17_6(0:2,fn:nf,fn:nf),
     & m27_5(0:2,fn:nf,fn:nf),
     & m27_6(0:2,fn:nf,fn:nf)

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

      integer:: isub
      common/isub/isub
      integer, parameter :: a(6)=(/5,5,7,6,6,7/),b(6)=(/6,7,5,7,5,6/),
     & c(6)=(/7,6,6,5,7,5/),
     & pntr(5:7,5:7)=reshape((/0,2,2,1,0,2,1,1,0/),(/3,3/))
      external qqb_zbjet,qqb_zbjet_gvec

      include 'cplx.h'

      ndmax=24

c-- initialize the matrix elements to zero
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

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
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msqab_c,msqab_cv,dsub,dsubv,subab_c,subab_cv,n)
      enddo

c--- initial-final/final-initial
      do n=1,6
      np6=n+6
      np12=n+12
      call dips(np6,p,1,a(n),b(n),dsub,dsubv,dummy,dummyv,
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msq1a_b,msq1a_bv,dsub,dsubv,sub1a_b,sub1a_bv,n)
      call dips(np6,p,b(n),a(n),1,dsub,dsubv,dummy,dummyv,
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msqba_1,msqba_1v,dsub,dsubv,subba_1,subba_1v,n)
      call dips(np12,p,2,c(n),b(n),dsub,dsubv,dummy,dummyv,
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msq2c_b,msq2c_bv,dsub,dsubv,sub2c_b,sub2c_bv,n)
      call dips(np12,p,b(n),c(n),2,dsub,dsubv,dummy,dummyv,
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msqbc_2,msqbc_2v,dsub,dsubv,subbc_2,subbc_2v,n)
      enddo

c--- initial-initial
      do n=1,3
      np18=n+18
      np21=n+21
      call dips(np18,p,1,b(n),2,dsub,dsubv,dummy,dummyv,
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msq1b_2,msq1b_2v,dsub,dsubv,sub1b_2,sub1b_2v,n)
      call dips(np21,p,2,b(n),1,dsub,dsubv,dummy,dummyv,
     & qqb_zbjet,qqb_zbjet_gvec)
      call storedip(msq2b_1,msq2b_1v,dsub,dsubv,sub2b_1,sub2b_1v,n)
      enddo

c--- skip the QQGGG contributions
c      goto 66

c-- fill the matrix elements
      do j=-nflav,nflav,nflav
      do k=-nflav,nflav,nflav

c--- QUARK-GLUON contributions (isub=1 only)
      if     ((isub == 1) .and. (j==+flav) .and. (k==0)) then
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
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *msq1a_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo

c--- GLUON-QUARK contributions (isub=1 only)
      elseif ((isub == 1) .and. (j==0) .and. (k==+flav)) then
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
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo

c--- GLUON-ANTIQUARK contributions (isub=1 only)
      elseif ((isub == 1) .and. (j==0) .and. (k==-flav)) then
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
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo

c--- ANTIQUARK-GLUON contributions (isub=1 only)
      elseif ((isub == 1) .and. (j==-flav) .and. (k==0)) then
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
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *msq1a_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo

c--- GLUON-GLUON contributions (isub=2 only)
      elseif ((isub == 2) .and. (j==0) .and. (k==0)) then
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

      enddo
      enddo

c--- DEBUG
c--- do not include the 4-quark contributions for now
c      return
c--- DEBUG

c      ndmax=12
c
c      call dips(1,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call dips(2,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c
c      call dips(3,p,5,7,6,sub57_6,sub57_6v,msq57_6,msq57_6v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m57_6,m57_6v)
c      call dips(4,p,6,7,5,sub67_5,sub67_5v,msq67_5,msq67_5v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m67_5,m67_5v)
c      call dips(5,p,1,7,5,sub17_5,sub17_5v,msq17_5,msq17_5v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m17_5,m17_5v)
c      call dips(5,p,5,7,1,sub57_1,sub57_1v,dummy,msq57_1v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call dips(6,p,1,7,6,sub17_6,sub17_6v,msq17_6,msq17_6v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m17_6,m17_6v)
c      call dips(6,p,6,7,1,sub67_1,sub67_1v,dummy,msq67_1v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call dips(7,p,2,7,6,sub27_6,sub27_6v,msq27_6,msq27_6v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m27_6,m27_6v)
c      call dips(7,p,6,7,2,sub67_2,sub67_2v,dummy,msq67_2v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call dips(8,p,2,7,5,sub27_5,sub27_5v,msq27_5,msq27_5v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m27_5,m27_5v)
c      call dips(8,p,5,7,2,sub57_2,sub57_2v,dummy,msq57_2v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call dips(9,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m17_2,m17_2v)
c      call dips(10,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call storedipcs(m27_1,m27_1v)
c
c      call dips(11,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c      call dips(12,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
c     & qqb_zbjet,qqb_zbjet_gvec)
c
c   66 continue


c--- construct the aliased matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq15_2(j,k)=
     & msq1b_2(3,0,j,k)+msq1b_2(3,1,j,k)+msq1b_2(3,2,j,k)
      msq25_1(j,k)=
     & msq2b_1(3,0,j,k)+msq2b_1(3,1,j,k)+msq2b_1(3,2,j,k)
      msq16_2(j,k)=
     & msq1b_2(1,0,j,k)+msq1b_2(1,1,j,k)+msq1b_2(1,2,j,k)
      msq16_2v(j,k)=
     & msq1b_2v(1,0,j,k)+msq1b_2v(1,1,j,k)+msq1b_2v(1,2,j,k)
      msq26_1(j,k)=
     & msq2b_1(1,0,j,k)+msq2b_1(1,1,j,k)+msq2b_1(1,2,j,k)
      msq26_1v(j,k)=
     & msq2b_1v(1,0,j,k)+msq2b_1v(1,1,j,k)+msq2b_1v(1,2,j,k)
      msq67_1v(j,k)=
     & msqba_1v(6,0,j,k)+msqba_1v(6,1,j,k)+msqba_1v(6,2,j,k)
      do i=0,2
      m57_6(i,j,k)=msqab_c(2,i,j,k)
      m67_5(i,j,k)=msqab_c(4,i,j,k)
      m17_2(i,j,k)=msq1b_2(2,i,j,k)
      m17_5(i,j,k)=msq1a_b(3,i,j,k)
      m17_6(i,j,k)=msq1a_b(6,i,j,k)
      m27_1(i,j,k)=msq2b_1(2,i,j,k)
      m27_5(i,j,k)=msq2c_b(5,i,j,k)
      m27_6(i,j,k)=msq2c_b(1,i,j,k)
      enddo
      msq17_6(j,k)=m17_6(0,j,k)+m17_6(1,j,k)+m17_6(2,j,k)
      msq17_2(j,k)=m17_2(0,j,k)+m17_2(1,j,k)+m17_2(2,j,k)
      msq27_1(j,k)=m27_1(0,j,k)+m27_1(1,j,k)+m27_1(2,j,k)
      enddo
      enddo
c--- construct the aliased subtraction terms
      do i=1,4
      sub15_2(i)=sub1b_2(3,i)
      sub25_1(i)=sub2b_1(3,i)
      sub16_2(i)=sub1b_2(1,i)
      sub26_1(i)=sub2b_1(1,i)
      sub17_2(i)=sub1b_2(2,i)
      sub27_1(i)=sub2b_1(2,i)
      sub57_6(i)=subab_c(2,i)
      sub67_5(i)=subab_c(4,i)
      sub17_5(i)=sub1a_b(3,i)
      sub57_1(i)=subba_1(3,i)
      sub17_6(i)=sub1a_b(6,i)
      sub67_1(i)=subba_1(6,i)
      sub27_5(i)=sub2c_b(5,i)
      sub57_2(i)=subbc_2(5,i)
      sub27_6(i)=sub2c_b(1,i)
      sub67_2(i)=subbc_2(1,i)
      enddo
      sub16_2v=sub1b_2v(1)
      sub26_1v=sub2b_1v(1)
      sub67_1v=subba_1v(6)
c--- end construct


cc--- compare the two calculations of the subtracted matrix elements
c      c15=0._dp
c      c15v=0._dp
c      c25=0._dp
c      c25v=0._dp
c      c16=0._dp
c      c16v=0._dp
c      c26=0._dp
c      c26v=0._dp
c      c57=0._dp
c      c57v=0._dp
c      c67=0._dp
c      c67v=0._dp
c      c175=0._dp
c      c175v=0._dp
c      c571=0._dp
c      c571v=0._dp
c      c176=0._dp
c      c176v=0._dp
c      c671=0._dp
c      c671v=0._dp
c      c275=0._dp
c      c275v=0._dp
c      c572=0._dp
c      c572v=0._dp
c      c276=0._dp
c      c276v=0._dp
c      c672=0._dp
c      c672v=0._dp
c      do j=-nf,nf
c      do k=-nf,nf
cc      write(6,*) msq26_1(j,k)-
cc     & (+msq2b_1(1,0,j,k)+msq2b_1(1,1,j,k)+msq2b_1(1,2,j,k))
c      c15=c15+real(j*10+k+99,dp)*(msq15_2(j,k)-
c     & (+msq1b_2(3,0,j,k)+msq1b_2(3,1,j,k)+msq1b_2(3,2,j,k)))
c      c15v=c15v+real(j*10+k+99,dp)*(msq15_2v(j,k)-
c     & (+msq1b_2v(3,0,j,k)+msq1b_2v(3,1,j,k)+msq1b_2v(3,2,j,k)))
c      c25=c25+real(j*10+k+99,dp)*(msq25_1(j,k)-
c     & (+msq2b_1(3,0,j,k)+msq2b_1(3,1,j,k)+msq2b_1(3,2,j,k)))
c      c25v=c25v+real(j*10+k+99,dp)*(msq25_1v(j,k)-
c     & (+msq2b_1v(3,0,j,k)+msq2b_1v(3,1,j,k)+msq2b_1v(3,2,j,k)))
c      c16=c16+real(j*10+k+99,dp)*(msq16_2(j,k)-
c     & (+msq1b_2(1,0,j,k)+msq1b_2(1,1,j,k)+msq1b_2(1,2,j,k)))
c      c16v=c16v+real(j*10+k+99,dp)*(msq16_2v(j,k)-
c     & (+msq1b_2v(1,0,j,k)+msq1b_2v(1,1,j,k)+msq1b_2v(1,2,j,k)))
c      c26=c26+real(j*10+k+99,dp)*(msq26_1(j,k)-
c     & (+msq2b_1(1,0,j,k)+msq2b_1(1,1,j,k)+msq2b_1(1,2,j,k)))
c      c26v=c26v+real(j*10+k+99,dp)*(msq26_1v(j,k)-
c     & (+msq2b_1v(1,0,j,k)+msq2b_1v(1,1,j,k)+msq2b_1v(1,2,j,k)))
c      c571v=c571v+real(j*10+k+99,dp)*(msq57_1v(j,k)-
c     & (+msqba_1v(3,0,j,k)+msqba_1v(3,1,j,k)+msqba_1v(3,2,j,k)))
c      c671v=c671v+real(j*10+k+99,dp)*(msq67_1v(j,k)-
c     & (+msqba_1v(6,0,j,k)+msqba_1v(6,1,j,k)+msqba_1v(6,2,j,k)))
c      c572v=c572v+real(j*10+k+99,dp)*(msq57_2v(j,k)-
c     & (+msqbc_2v(5,0,j,k)+msqbc_2v(5,1,j,k)+msqbc_2v(5,2,j,k)))
c      c672v=c672v+real(j*10+k+99,dp)*(msq67_2v(j,k)-
c     & (+msqbc_2v(1,0,j,k)+msqbc_2v(1,1,j,k)+msqbc_2v(1,2,j,k)))
c
c      do i=0,2
c      c57=c57+real(j*10+k+99,dp)*(m57_6(i,j,k)-msqab_c(2,i,j,k))
c      c57v=c57v+real(j*10+k+99,dp)*(m57_6v(i,j,k)-msqab_cv(2,i,j,k))
c      c67=c67+real(j*10+k+99,dp)*(m67_5(i,j,k)-msqab_c(4,i,j,k))
c      c67v=c67v+real(j*10+k+99,dp)*(m67_5v(i,j,k)-msqab_cv(4,i,j,k))
c      c175=c175+real(j*10+k+99,dp)*(m17_5(i,j,k)-msq1a_b(3,i,j,k))
c      c175v=c175v+real(j*10+k+99,dp)*(m17_5v(i,j,k)-msq1a_bv(3,i,j,k))
c      c176=c176+real(j*10+k+99,dp)*(m17_6(i,j,k)-msq1a_b(6,i,j,k))
c      c176v=c176v+real(j*10+k+99,dp)*(m17_6v(i,j,k)-msq1a_bv(6,i,j,k))
c      c275=c275+real(j*10+k+99,dp)*(m27_5(i,j,k)-msq2c_b(5,i,j,k))
c      c275v=c275v+real(j*10+k+99,dp)*(m27_5v(i,j,k)-msq2c_bv(5,i,j,k))
c      c276=c276+real(j*10+k+99,dp)*(m27_6(i,j,k)-msq2c_b(1,i,j,k))
c      c276v=c276v+real(j*10+k+99,dp)*(m27_6v(i,j,k)-msq2c_bv(1,i,j,k))
c      enddo
c      enddo
c      enddo
c
c      c25=c25+sub15_2(qq)+sub15_2(gq)+sub15_2(qg)+sub15_2(gg)
c     & -sub1b_2(3,qq)-sub1b_2(3,gq)-sub1b_2(3,qg)-sub1b_2(3,gg)
c      c25v=c25v+sub15_2v-sub1b_2v(3)
c      c25=c25+sub25_1(qq)+sub25_1(gq)+sub25_1(qg)+sub25_1(gg)
c     & -sub2b_1(3,qq)-sub2b_1(3,gq)-sub2b_1(3,qg)-sub2b_1(3,gg)
c      c25v=c25v+sub25_1v-sub2b_1v(3)
c      c26=c26+sub16_2(qq)+sub16_2(gq)+sub16_2(qg)+sub16_2(gg)
c     & -sub1b_2(1,qq)-sub1b_2(1,gq)-sub1b_2(1,qg)-sub1b_2(1,gg)
c      c26v=c26v+sub16_2v-sub1b_2v(1)
c      c26=c26+sub26_1(qq)+sub26_1(gq)+sub26_1(qg)+sub26_1(gg)
c     & -sub2b_1(1,qq)-sub2b_1(1,gq)-sub2b_1(1,qg)-sub2b_1(1,gg)
c      c26v=c26v+sub26_1v-sub2b_1v(1)
c      c57=c57+sub57_6(qq)+sub57_6(gq)+sub57_6(qg)+sub57_6(gg)
c     & -subab_c(2,qq)-subab_c(2,gq)-subab_c(2,qg)-subab_c(2,gg)
c      c57v=c57v+sub57_6v-subab_cv(2)
c      c67=c67+sub67_5(qq)+sub67_5(gq)+sub67_5(qg)+sub67_5(gg)
c     & -subab_c(4,qq)-subab_c(4,gq)-subab_c(4,qg)-subab_c(4,gg)
c      c67v=c67v+sub67_5v-subab_cv(4)
c      c175=c175+sub17_5(qq)+sub17_5(gq)+sub17_5(qg)+sub17_5(gg)
c     & -sub1a_b(3,qq)-sub1a_b(3,gq)-sub1a_b(3,qg)-sub1a_b(3,gg)
c      c175v=c175v+sub17_5v-sub1a_bv(3)
c      c571=c571+sub57_1(qq)+sub57_1(gq)+sub57_1(qg)+sub57_1(gg)
c     & -subba_1(3,qq)-subba_1(3,gq)-subba_1(3,qg)-subba_1(3,gg)
c      c571v=c571v+sub57_1v-subba_1v(3)
c      c176=c176+sub17_6(qq)+sub17_6(gq)+sub17_6(qg)+sub17_6(gg)
c     & -sub1a_b(6,qq)-sub1a_b(6,gq)-sub1a_b(6,qg)-sub1a_b(6,gg)
c      c176v=c176v+sub17_6v-sub1a_bv(6)
c      c671=c671+sub67_1(qq)+sub67_1(gq)+sub67_1(qg)+sub67_1(gg)
c     & -subba_1(6,qq)-subba_1(6,gq)-subba_1(6,qg)-subba_1(6,gg)
c      c671v=c671v+sub67_1v-subba_1v(6)
c      c275=c275+sub27_5(qq)+sub27_5(gq)+sub27_5(qg)+sub27_5(gg)
c     & -sub2c_b(5,qq)-sub2c_b(5,gq)-sub2c_b(5,qg)-sub2c_b(5,gg)
c      c275v=c275v+sub27_5v-sub2c_bv(5)
c      c572=c572+sub57_2(qq)+sub57_2(gq)+sub57_2(qg)+sub57_2(gg)
c     & -subbc_2(5,qq)-subbc_2(5,gq)-subbc_2(5,qg)-subbc_2(5,gg)
c      c572v=c572v+sub57_2v-subbc_2v(5)
c      c276=c276+sub27_6(qq)+sub27_6(gq)+sub27_6(qg)+sub27_6(gg)
c     & -sub2c_b(1,qq)-sub2c_b(1,gq)-sub2c_b(1,qg)-sub2c_b(1,gg)
c      c276v=c276v+sub27_6v-sub2c_bv(1)
c      c672=c672+sub67_2(qq)+sub67_2(gq)+sub67_2(qg)+sub67_2(gg)
c     & -subbc_2(1,qq)-subbc_2(1,gq)-subbc_2(1,qg)-subbc_2(1,gg)
c      c672v=c672v+sub67_2v-subbc_2v(1)
c      write(6,*) 'checksum 15   ',c15
c      write(6,*) 'checksum 15v  ',c15v
c      write(6,*) 'checksum 25   ',c25
c      write(6,*) 'checksum 25v  ',c25v
c      write(6,*) 'checksum 16   ',c16
c      write(6,*) 'checksum 16v  ',c16v
c      write(6,*) 'checksum 26   ',c26
c      write(6,*) 'checksum 26v  ',c26v
c      write(6,*) 'checksum 57   ',c57
c      write(6,*) 'checksum 57v  ',c57v
c      write(6,*) 'checksum 67   ',c67
c      write(6,*) 'checksum 67v  ',c67v
c      write(6,*) 'checksum 175  ',c175
c      write(6,*) 'checksum 175v ',c175v
c      write(6,*) 'checksum 571  ',c571
c      write(6,*) 'checksum 571v ',c571v
c      write(6,*) 'checksum 176  ',c176
c      write(6,*) 'checksum 176v ',c176v
c      write(6,*) 'checksum 671  ',c671
c      write(6,*) 'checksum 671v ',c671v
c      write(6,*) 'checksum 275  ',c275
c      write(6,*) 'checksum 275v ',c275v
c      write(6,*) 'checksum 572  ',c572
c      write(6,*) 'checksum 572v ',c572v
c      write(6,*) 'checksum 276  ',c276
c      write(6,*) 'checksum 276v ',c276v
c      write(6,*) 'checksum 672  ',c672
c      write(6,*) 'checksum 672v ',c672v
c      pause


c--- note that singularities for p7 in the GQ,... contributions
c--- should not be included, because the basic (LO) process would
c--- then contain two heavy quarks. Similarly, 56 singularities are
c--- not included because the LO process would contain no heavy quarks;
c--- however, we must remember to apply a cut on the raw matrix
c--- elements to eliminate this singularity

c--- fill the dipole contributions
      do j=-nflav,nflav
      do k=-nflav,nflav

      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 89
c--- do not allow abs(j) = flav and abs(k) = flav
      if ((abs(j) > 0) .and. (abs(k) > 0) .and.
     &    (abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 89
c--- if both are (anti-)quarks, one of them should be a heavy quark

c---Q-Q contribution (isub=1 only)
      if     ((isub == 1) .and. (j > 0) .and. (k>0)) then
        if     (j == +flav) then
          msq(22,j,k)=msq(22,j,k)+(aveqq/aveqg)*(
     &     sub26_1(gq)*msq26_1(j,0)+sub26_1v*msq26_1v(j,0))

      msq(13,j,k)=msq(13,j,k)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
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
        elseif (k == +flav) then
          msq(19,j,k)=msq(19,j,k)+(aveqq/aveqg)*(
     &     sub16_2(gq)*msq16_2(0,k)+sub16_2v*msq16_2v(0,k))

      msq(17,j,k)=msq(17,j,k)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
        endif
c---Qbar-Qbar contribution (isub=1 only)
      elseif ((isub == 1) .and. (j < 0) .and. (k<0)) then
        if     (j == -flav) then
          msq(22,j,k)=msq(22,j,k)+(aveqq/aveqg)*(
     &     sub26_1(gq)*msq26_1(j,0)+sub26_1v*msq26_1v(j,0))

      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
        elseif (k == -flav) then
          msq(19,j,k)=msq(19,j,k)+(aveqq/aveqg)*(
     &     sub16_2(gq)*msq16_2(0,k)+sub16_2v*msq16_2v(0,k))

      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
        endif

c---Q-Qbar contribution (isub=1 only)
      elseif ((isub == 1) .and. (j > 0) .and. (k<0)) then
        if     (j == +flav) then
          msq(22,j,k)=msq(22,j,k)+(aveqq/aveqg)*(
     &     sub26_1(gq)*msq26_1(j,0)+sub26_1v*msq26_1v(j,0))

      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(1,j,k)-m57_6(0,j,k)/xn-m57_6(2,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(1,j,k)-m67_5(0,j,k)/xn-m67_5(2,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
     & +sub67_1(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
     & +sub57_2(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(1,j,k)-m17_2(0,j,k)/xn-m17_2(2,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(1,j,k)-m27_1(0,j,k)/xn-m27_1(2,j,k)/xn)

        elseif (k == -flav) then
          msq(19,j,k)=msq(19,j,k)+(aveqq/aveqg)*(
     &     sub16_2(gq)*msq16_2(0,k)+sub16_2v*msq16_2v(0,k))

      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(1,j,k)-m67_5(0,j,k)/xn-m67_5(2,j,k)/xn)
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(1,j,k)-m57_6(0,j,k)/xn-m57_6(2,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub17_5(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
     & +sub57_1(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub27_6(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
     & +sub67_2(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(1,j,k)-m17_2(0,j,k)/xn-m17_2(2,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(1,j,k)-m27_1(0,j,k)/xn-m27_1(2,j,k)/xn)
        endif

c---Qbar-Q contribution (isub=1 only)
      elseif ((isub == 1) .and. (j < 0) .and. (k>0)) then
        if     (j == -flav) then
          msq(22,j,k)=msq(22,j,k)+(aveqq/aveqg)*(
     &     sub26_1(gq)*msq26_1(j,0)+sub26_1v*msq26_1v(j,0))

      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(1,j,k)-m17_2(0,j,k)/xn-m17_2(2,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(1,j,k)-m27_1(0,j,k)/xn-m27_1(2,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub57_2(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
     & +sub27_5(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub67_1(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
     & +sub17_6(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(1,j,k)-m57_6(0,j,k)/xn-m57_6(2,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(1,j,k)-m67_5(0,j,k)/xn-m67_5(2,j,k)/xn)

        elseif (k == +flav) then
          msq(19,j,k)=msq(19,j,k)+(aveqq/aveqg)*(
     &     sub16_2(gq)*msq16_2(0,k)+sub16_2v*msq16_2v(0,k))

      msq(20,j,k)=msq(20,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(1,j,k)-m17_2(0,j,k)/xn-m17_2(2,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(1,j,k)-m27_1(0,j,k)/xn-m27_1(2,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub67_2(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
     & +sub27_6(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub57_1(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
     & +sub17_5(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(1,j,k)-m67_5(0,j,k)/xn-m67_5(2,j,k)/xn)
      msq(2,j,k)=msq(2,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(1,j,k)-m57_6(0,j,k)/xn-m57_6(2,j,k)/xn)
        endif


      elseif ((j > 0) .and. (k == 0)) then
c--- Q-g contribution (isub=1 only)
        if     ((isub == 1) .and. (j == +flav)) then
          msq(22,j,k)=msq(22,j,k)+(
     &     +sub26_1(qg)*msq26_1(j,1)*2._dp
     &     +sub26_1(qg)*msq26_1(j,2)*real(flav-3,dp))
          msq(23,j,k)=msq(23,j,k)+(
     &     +sub27_1(qg)*msq27_1(j,-1)*2._dp
     &     +sub27_1(qg)*msq27_1(j,-2)*real(flav-3,dp))
          msq(12,j,k)=msq(12,j,k)+real(flav-1,dp)*(
     &     +sub67_1(gq)*msq17_6(j,k)-sub67_1v*msq67_1v(j,k))
c--- q-g contribution (isub=2 only)
        elseif ((isub == 2) .and. (j .ne. +flav)) then
          msq(22,j,k)=msq(22,j,k)+(
     &     +sub26_1(qg)*msq26_1(j,+flav))
          msq(21,j,k)=msq(21,j,k)+(
     &     +sub25_1(qg)*msq25_1(j,-flav))
        endif

      elseif ((j < 0) .and. (k == 0)) then
c--- Qb-g contribution (isub=1 only)
        if     ((isub == 1) .and. (j == -flav)) then
          msq(22,j,k)=msq(22,j,k)+(
     &     +sub26_1(qg)*msq26_1(j,1)*2._dp
     &     +sub26_1(qg)*msq26_1(j,2)*real(flav-3,dp))
          msq(23,j,k)=msq(23,j,k)+(
     &     +sub27_1(qg)*msq27_1(j,-1)*2._dp
     &     +sub27_1(qg)*msq27_1(j,-2)*real(flav-3,dp))
          msq(12,j,k)=msq(12,j,k)+real(flav-1,dp)*(
     &     +sub67_1(gq)*msq17_6(j,k)-sub67_1v*msq67_1v(j,k))
c--- qb-g contribution (isub=2 only)
        elseif ((isub == 2) .and. (j .ne. -flav)) then
          msq(22,j,k)=msq(22,j,k)+(
     &     +sub26_1(qg)*msq26_1(j,+flav))
          msq(21,j,k)=msq(21,j,k)+(
     &     +sub25_1(qg)*msq25_1(j,-flav))
        endif

      elseif ((j == 0) .and. (k > 0)) then
c--- g-Q contribution (isub=1 only)
        if     ((isub == 1) .and. (k == +flav)) then
          msq(19,j,k)=msq(19,j,k)+(
     &     +sub16_2(qg)*msq16_2(1,k)*2._dp
     &     +sub16_2(qg)*msq16_2(2,k)*real(flav-3,dp))
          msq(20,j,k)=msq(20,j,k)+(
     &     +sub17_2(qg)*msq17_2(-1,k)*2._dp
     &     +sub17_2(qg)*msq17_2(-2,k)*real(flav-3,dp))
          msq(12,j,k)=msq(12,j,k)+real(flav-1,dp)*(
     &     +sub67_1(gq)*msq17_6(j,k)-sub67_1v*msq67_1v(j,k))
c--- g-q contribution (isub=2 only)
        elseif ((isub == 2) .and. (k .ne. +flav)) then
          msq(19,j,k)=msq(19,j,k)+(
     &     +sub16_2(qg)*msq16_2(+flav,k))
          msq(24,j,k)=msq(24,j,k)+(
     &     +sub15_2(qg)*msq15_2(-flav,k))
        endif

      elseif ((j == 0) .and. (k < 0)) then
c--- g-Qb contribution (isub=1 only)
        if     ((isub == 1) .and. (k == -flav)) then
          msq(19,j,k)=msq(19,j,k)+(
     &     +sub16_2(qg)*msq16_2(1,k)*2._dp
     &     +sub16_2(qg)*msq16_2(2,k)*real(flav-3,dp))
          msq(20,j,k)=msq(20,j,k)+(
     &     +sub17_2(qg)*msq17_2(-1,k)*2._dp
     &     +sub17_2(qg)*msq17_2(-2,k)*real(flav-3,dp))
          msq(12,j,k)=msq(12,j,k)+real(flav-1,dp)*(
     &     +sub67_1(gq)*msq17_6(j,k)-sub67_1v*msq67_1v(j,k))
c--- g-qb contribution (isub=2 only)
        elseif ((isub == 2) .and. (k .ne. -flav)) then
          msq(19,j,k)=msq(19,j,k)+(
     &     +sub16_2(qg)*msq16_2(+flav,k))
          msq(24,j,k)=msq(24,j,k)+(
     &     +sub15_2(qg)*msq15_2(-flav,k))
        endif

      endif

   89 continue

      enddo
      enddo

      return
      end

      subroutine storedipcs(msq_dip,msq_dipv)
      implicit none
      include 'types.f'
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'msq_cs.f'
      include 'msqv_cs.f'
      integer:: i,j,k
      real(dp):: msq_dip(0:2,-nf:nf,-nf:nf),
     &                 msq_dipv(0:2,-nf:nf,-nf:nf)

      do i=0,2
        do j=-nf,nf
        do k=-nf,nf
          msq_dip(i,j,k)=msq_cs(i,j,k)
          msq_dipv(i,j,k)=msqv_cs(i,j,k)
        enddo
        enddo
      enddo

      return
      end

