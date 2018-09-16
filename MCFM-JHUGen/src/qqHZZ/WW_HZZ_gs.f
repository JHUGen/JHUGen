      subroutine WW_HZZ_gs(p,msq)
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
*                           |                                          *
*                           |                                          *
*                           |                                          *
*                           ---> W+(nu(p3)+e+(p4))+W-(e-(p5)+nub(p6))  *
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
     & msq19_7(-nf:nf,-nf:nf),msq29_8(-nf:nf,-nf:nf),
     & msq18_2(-nf:nf,-nf:nf),msq29_1(-nf:nf,-nf:nf),
     & msq17_2(-nf:nf,-nf:nf),msq28_1(-nf:nf,-nf:nf),
     & sub19_7(4),sub29_8(4),sub79_1(4),sub89_2(4),
     & sub18_2(4),sub29_1(4),
     & sub17_2(4),sub28_1(4),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),dsubv
      external WW_HZZ,donothing_gvec

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
      call dips(1,p,1,9,7,sub19_7,dsubv,msq19_7,dummyv,
     & WW_HZZ,donothing_gvec)
      call dips(1,p,7,9,1,sub79_1,dsubv,dummy,dummyv,
     & WW_HZZ,donothing_gvec)

      call dips(2,p,2,9,8,sub29_8,dsubv,msq29_8,dummyv,
     & WW_HZZ,donothing_gvec)
      call dips(2,p,8,9,2,sub89_2,dsubv,dummy,dummyv,
     & WW_HZZ,donothing_gvec)

      call dips(3,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     & WW_HZZ,donothing_gvec)
      call dips(4,p,2,8,1,sub28_1,dsubv,msq28_1,dummyv,
     & WW_HZZ,donothing_gvec)

      call dips(5,p,1,8,2,sub18_2,dsubv,msq18_2,dummyv,
     & WW_HZZ,donothing_gvec)
      call dips(6,p,2,9,1,sub29_1,dsubv,msq29_1,dummyv,
     & WW_HZZ,donothing_gvec)


c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1

      if  ((j == 0) .and. (k == 0)) then
         goto 20
      elseif ((j .ne. 0) .and. (k .ne. 0)) then
         if     ((j > 0) .and. (k < 0)) then
c--- q-qb case
         msq(1,j,k)=2._dp*cf*(sub19_7(qq)+sub79_1(qq))*msq19_7(j,k)
         msq(2,j,k)=2._dp*cf*(sub29_8(qq)+sub89_2(qq))*msq29_8(j,k)
         elseif ((j < 0) .and. (k > 0)) then
c--- qb-q case
         msq(1,j,k)=2._dp*cf*(sub19_7(qq)+sub79_1(qq))*msq19_7(j,k)
         msq(2,j,k)=2._dp*cf*(sub29_8(qq)+sub89_2(qq))*msq29_8(j,k)
         else
c--- q-q and qb-qb cases
         msq(1,j,k)=2._dp*cf*(sub19_7(qq)+sub79_1(qq))*msq19_7(j,k)
         msq(2,j,k)=2._dp*cf*(sub29_8(qq)+sub89_2(qq))*msq29_8(j,k)
         endif
C---qg
      elseif ((j > 0) .and. (k == 0)) then
         msq(6,j,k)=2._dp*tr*sub29_1(qg)*(
     &              +msq29_1(j,+1)+msq29_1(j,+2)+msq29_1(j,+3)
     &              +msq29_1(j,+4)+msq29_1(j,+5))
         msq(4,j,k)=2._dp*tr*sub28_1(qg)*(
     &              +msq28_1(j,-1)+msq28_1(j,-2)+msq28_1(j,-3)
     &              +msq28_1(j,-4)+msq28_1(j,-5))
C---qbg
      elseif ((j < 0) .and. (k == 0)) then
         msq(6,j,k)=2._dp*tr*sub29_1(qg)*(
     &              +msq29_1(j,-5)+msq29_1(j,-4)+msq29_1(j,-3)
     &              +msq29_1(j,-2)+msq29_1(j,-1))
         msq(4,j,k)=2._dp*tr*sub28_1(qg)*(
     &              +msq28_1(j,+1)+msq28_1(j,+2)+msq28_1(j,+3)
     &              +msq28_1(j,+4)+msq28_1(j,+5))
C---gq
       elseif ((j == 0) .and. (k > 0)) then
         msq(5,j,k)=2._dp*tr*sub18_2(qg)*(
     &              +msq18_2(+5,k)+msq18_2(+4,k)+msq18_2(+3,k)
     &              +msq18_2(+2,k)+msq18_2(+1,k))
         msq(3,j,k)=2._dp*tr*sub17_2(qg)*(
     &              +msq17_2(-1,k)+msq17_2(-2,k)+msq17_2(-3,k)
     &              +msq17_2(-4,k)+msq17_2(-5,k))
C---gqb
      elseif ((j == 0) .and. (k < 0)) then
         msq(5,j,k)=2._dp*tr*sub18_2(qg)*(
     &              +msq18_2(-5,k)+msq18_2(-4,k)+msq18_2(-3,k)
     &              +msq18_2(-2,k)+msq18_2(-1,k))
         msq(3,j,k)=2._dp*tr*sub17_2(qg)*(
     &              +msq17_2(+1,k)+msq17_2(+2,k)+msq17_2(+3,k)
     &              +msq17_2(+4,k)+msq17_2(+5,k))
      endif
 20   continue
      enddo
      enddo

      return
      end


