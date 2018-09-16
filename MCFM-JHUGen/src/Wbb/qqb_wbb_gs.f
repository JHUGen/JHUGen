      subroutine qqb_wbb_gs(P,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  b(p5)+bb(p6) + W +g(p7)
c                                          |
c                                          --> nu(p3)+e^+(p4)
c
c   positively charged W only
c--all momenta incoming so some of the signs may look odd.
c--initial state gluons coupling to b not included

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'masses.f'
      integer:: j,k,nd
c --- remember: nd will count the dipoles

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),dot
      real(dp)::
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq57_6(-nf:nf,-nf:nf),msq67_5(-nf:nf,-nf:nf),
     & msq17_5(-nf:nf,-nf:nf),
     & msq17_6(-nf:nf,-nf:nf),
     & msq27_5(-nf:nf,-nf:nf),
     & msq27_6(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),sub57_6(4),sub67_5(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),dsubv
      external qqb_wbb,donothing_gvec

      ndmax=8

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

      if (
     &      (two*dot(p,5,6) < four*mbsq)
     & .or. (two*dot(p,1,5)*dot(p,2,5)/dot(p,1,2) < mbsq)
     & .or. (two*dot(p,1,6)*dot(p,2,6)/dot(p,1,2) < mbsq))return

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     & qqb_wbb,donothing_gvec)

c--- final-final
      call dips(3,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(4,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     & qqb_wbb,donothing_gvec)

c--- now the basic initial final and final initial
      call dips(5,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     & qqb_wbb,donothing_gvec)
c--- called in this fashion the routine only supplies new values for
c--- sub..
      call dips(5,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(6,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(6,p,6,7,1,sub67_1,dsubv,dummy,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(7,p,2,7,5,sub27_5,dsubv,msq27_5,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(7,p,5,7,2,sub57_2,dsubv,dummy,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(8,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     & qqb_wbb,donothing_gvec)
      call dips(8,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     & qqb_wbb,donothing_gvec)

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)

      if ((j>0).and.(k<0)) then
c----------q-qb

      msq(1,j,k)=-sub17_2(qq)*msq17_2(j,k)/xn
      msq(2,j,k)=-sub27_1(qq)*msq27_1(j,k)/xn
      msq(3,j,k)=-sub57_6(qq)*msq57_6(j,k)/xn
      msq(4,j,k)=-sub67_5(qq)*msq67_5(j,k)/xn
      msq(5,j,k)=+sub17_5(qq)*msq17_5(j,k)*(xn-two/xn)
     &           +sub57_1(qq)*msq17_5(j,k)*(xn-two/xn)
      msq(6,j,k)=+sub17_6(qq)*msq17_6(j,k)*two/xn
     &           +sub67_1(qq)*msq17_6(j,k)*two/xn
      msq(7,j,k)=+sub27_5(qq)*msq27_5(j,k)*two/xn
     &           +sub57_2(qq)*msq27_5(j,k)*two/xn
      msq(8,j,k)=+sub27_6(qq)*msq27_6(j,k)*(xn-two/xn)
     &           +sub67_2(qq)*msq27_6(j,k)*(xn-two/xn)


      elseif((j<0).and.(k>0)) then
c----------qb-q

      msq(1,j,k)=-sub17_2(qq)*msq17_2(j,k)/xn
      msq(2,j,k)=-sub27_1(qq)*msq27_1(j,k)/xn
      msq(3,j,k)=-sub57_6(qq)*msq57_6(j,k)/xn
      msq(4,j,k)=-sub67_5(qq)*msq67_5(j,k)/xn
      msq(5,j,k)=+sub17_5(qq)*msq17_5(j,k)*two/xn
     &           +sub57_1(qq)*msq17_5(j,k)*two/xn
      msq(6,j,k)=+sub17_6(qq)*msq17_6(j,k)*(xn-two/xn)
     &           +sub67_1(qq)*msq17_6(j,k)*(xn-two/xn)
      msq(7,j,k)=+sub27_5(qq)*msq27_5(j,k)*(xn-two/xn)
     &           +sub57_2(qq)*msq27_5(j,k)*(xn-two/xn)
      msq(8,j,k)=+sub27_6(qq)*msq27_6(j,k)*two/xn
     &           +sub67_2(qq)*msq27_6(j,k)*two/xn

      elseif (j==0) then
          if (k>0) then
c---------g-q
          msq(1,j,k)=sub17_2(qg)
     &    *(msq17_2(-1,k)+msq17_2(-2,k)+msq17_2(-3,k)+msq17_2(-4,k)
     &     +msq17_2(-5,k))
          elseif (k<0) then
c---------g-qbar
          msq(1,j,k)=sub17_2(qg)
     &    *(msq17_2(+1,k)+msq17_2(+2,k)+msq17_2(+3,k)+msq17_2(+4,k)
     &    +msq17_2(+5,k))
          endif
      elseif (k==0) then
          if (j>0) then
c---------q-g
          msq(2,j,k)=sub27_1(qg)
     &    *(msq27_1(j,-1)+msq27_1(j,-2)+msq27_1(j,-3)+msq27_1(j,-4)
     &    +msq27_1(j,-5))
          elseif (j<0) then
c---------qbar-g
          msq(2,j,k)=sub27_1(qg)
     &    *(msq27_1(j,+1)+msq27_1(j,+2)+msq27_1(j,+3)+msq27_1(j,+4)
     &     +msq27_1(j,+5))
          endif
      endif

      enddo
      enddo

      return
      end





