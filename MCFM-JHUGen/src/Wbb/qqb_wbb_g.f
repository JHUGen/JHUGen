      subroutine qqb_wbb_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998                                                       *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+bb(p6)
c   positively charged W only

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'noglue.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),fac
      real(dp):: qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(7,p,za,zb)

      if (
     &      (s(5,6) < four*mbsq)
     & .or. (s(1,5)*s(2,5)/s(1,2) < mbsq)
     & .or. (s(1,6)*s(2,6)/s(1,2) < mbsq) ) return

      fac=gsq**3*gw**4/4._dp

c--- shortcut if we're doing gqonly
      if (gqonly) then
        qqbWbbg=0._dp
        qbqWbbg=0._dp
      else
        call wbbgamp(1,2,7,5,6,3,4,qqbWbbg)
        call wbbgamp(2,1,7,5,6,3,4,qbqWbbg)
      endif
      call wbbgamp(1,7,2,5,6,3,4,qgWbbq)
      call wbbgamp(7,1,2,5,6,3,4,qbgWbbqb)
      call wbbgamp(2,7,1,5,6,3,4,gqWbbq)
      call wbbgamp(7,2,1,5,6,3,4,gqbWbbqb)


      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)

      msq(j,k)=0._dp

      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=fac*aveqq*Vsq(j,k)*qqbWbbg

      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=fac*aveqq*Vsq(j,k)*qbqWbbg

      elseif ((j > 0) .and. (k == 0)) then
      msq(j,k)=fac*aveqg*Vsum(j)*qgWbbq

      elseif ((j < 0) .and. (k == 0)) then
      msq(j,k)=fac*aveqg*Vsum(j)*qbgWbbqb

      elseif ((j == 0) .and. (k > 0)) then
      msq(j,k)=fac*aveqg*Vsum(k)*gqWbbq

      elseif ((j == 0) .and. (k < 0)) then
      msq(j,k)=fac*aveqg*Vsum(k)*gqbWbbqb
      endif

      enddo
      enddo

      return
      end


