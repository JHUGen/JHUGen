      subroutine WW_Hqq_gs(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2002.                                                      *
*                                                                      *
*     Weak Boson Fusion by W-W exchange only                           *
*     This routine calculates the dipole subtraction terms             *
*     for the process:                                                 *
*     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k,nd

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq17_5(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     & msq16_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq15_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & sub17_5(4),sub27_6(4),sub57_1(4),sub67_2(4),
     & sub16_2(4),sub27_1(4),
     & sub15_2(4),sub26_1(4),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),dsubv
      external WW_Hqq,donothing_gvec

      ndmax=6

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

c---- calculate the dipoles: initial-final and final-initial
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     & WW_Hqq,donothing_gvec)
      call dips(1,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     & WW_Hqq,donothing_gvec)

      call dips(2,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     & WW_Hqq,donothing_gvec)
      call dips(2,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     & WW_Hqq,donothing_gvec)

      call dips(3,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     & WW_Hqq,donothing_gvec)
      call dips(4,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     & WW_Hqq,donothing_gvec)

      call dips(5,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     & WW_Hqq,donothing_gvec)
      call dips(6,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     & WW_Hqq,donothing_gvec)


c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1

      if  ((j == 0) .and. (k == 0)) then
         goto 20
      elseif ((j .ne. 0) .and. (k .ne. 0)) then
         if     ((j > 0) .and. (k < 0)) then
c--- q-qb case
         msq(1,j,k)=2._dp*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2._dp*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
         elseif ((j < 0) .and. (k > 0)) then
c--- qb-q case
         msq(1,j,k)=2._dp*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2._dp*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
         else
c--- q-q and qb-qb cases
         msq(1,j,k)=2._dp*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2._dp*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
         endif
C---qg
      elseif ((j > 0) .and. (k == 0)) then
         msq(6,j,k)=2._dp*tr*sub27_1(qg)*(
     &              +msq27_1(j,+1)+msq27_1(j,+2)+msq27_1(j,+3)
     &              +msq27_1(j,+4)+msq27_1(j,+5))
         msq(4,j,k)=2._dp*tr*sub26_1(qg)*(
     &              +msq26_1(j,-1)+msq26_1(j,-2)+msq26_1(j,-3)
     &              +msq26_1(j,-4)+msq26_1(j,-5))
C---qbg
      elseif ((j < 0) .and. (k == 0)) then
         msq(6,j,k)=2._dp*tr*sub27_1(qg)*(
     &              +msq27_1(j,-5)+msq27_1(j,-4)+msq27_1(j,-3)
     &              +msq27_1(j,-2)+msq27_1(j,-1))
         msq(4,j,k)=2._dp*tr*sub26_1(qg)*(
     &              +msq26_1(j,+1)+msq26_1(j,+2)+msq26_1(j,+3)
     &              +msq26_1(j,+4)+msq26_1(j,+5))
C---gq
       elseif ((j == 0) .and. (k > 0)) then
         msq(5,j,k)=2._dp*tr*sub16_2(qg)*(
     &              +msq16_2(+5,k)+msq16_2(+4,k)+msq16_2(+3,k)
     &              +msq16_2(+2,k)+msq16_2(+1,k))
         msq(3,j,k)=2._dp*tr*sub15_2(qg)*(
     &              +msq15_2(-1,k)+msq15_2(-2,k)+msq15_2(-3,k)
     &              +msq15_2(-4,k)+msq15_2(-5,k))
C---gqb
      elseif ((j == 0) .and. (k < 0)) then
         msq(5,j,k)=2._dp*tr*sub16_2(qg)*(
     &              +msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     &              +msq16_2(-2,k)+msq16_2(-1,k))
         msq(3,j,k)=2._dp*tr*sub15_2(qg)*(
     &              +msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
     &              +msq15_2(+4,k)+msq15_2(+5,k))
      endif
 20   continue
      enddo
      enddo

      return
      end


