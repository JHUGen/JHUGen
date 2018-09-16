      subroutine qqb_Hg_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     October, 2001.                                                   *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
*     q(-p1)+qbar(-p2) -->  H + parton(p5) + parton(p6)                *
*                           |                                          *
*                            -->b(p3)+b~(p4)                           *
*                                                                      *
*  ISUB picks out the sub-process as follows:                          *
*    isub = 0 : all of the processes below                             *
*    isub = 1 : g + b -> b + g ,  b + q  -> b + q                      *
*    isub = 2 : g + g -> b~ + b , b + b  -> b + b                      *
************************************************************************


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k,nd,isub
c --- remember: nd will count the dipoles

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq15_6(-nf:nf,-nf:nf),msq26_5(-nf:nf,-nf:nf),
     & msq56_1v(-nf:nf,-nf:nf),msq56_2v(-nf:nf,-nf:nf),
     & msq26_5v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & msq15_6v(-nf:nf,-nf:nf),msq16_2v(-nf:nf,-nf:nf),
     & msq25_1v(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub16_2(4),sub26_1(4),
     & sub15_6(4),sub16_5(4),sub25_6(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub25_1v,sub26_1v,sub16_5v,sub16_2v,sub15_2v,sub15_6v,
     & sub25_6v
      common/isub/isub
      external qqb_Hg,qqb_H_gvec,donothing_gvec

      ndmax=6

      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_Hg,qqb_H_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_Hg,qqb_H_gvec)
      call dips(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_Hg,qqb_H_gvec)
      call dips(4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_Hg,qqb_H_gvec)

c--- called in this fashion the routine only supplies new values for
c--- sub... and sub...v

      call dips(5,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & qqb_Hg,qqb_H_gvec)
      call dips(5,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     & qqb_Hg,qqb_H_gvec)
      call dips(5,p,1,6,5,sub16_5,sub16_5v,dummy,dummyv,
     & qqb_Hg,donothing_gvec)

      call dips(6,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & qqb_Hg,qqb_H_gvec)
      call dips(6,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     & qqb_Hg,qqb_H_gvec)
      call dips(6,p,2,5,6,sub25_6,sub25_6v,dummy,dummyv,
     & qqb_Hg,donothing_gvec)

c--- 2-quark singularities
      do j=-5,+5,5
      do k=-5,+5,5

      if (j .ne. 0 .and. k .ne. 0 .and. j.ne.-k) goto 19

c--- do only q-qb and qb-q cases
c--- there are no such initial states for the 2-gluon piece here
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
C-----half=statistical factor
c      msq(1,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
c      msq(2,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
c      msq(3,j,k)=-half*msq16_2(j,k)*sub16_2(qq)/xn
c      msq(4,j,k)=-half*msq26_1(j,k)*sub26_1(qq)/xn
c      msq(5,j,k)=half*xn*(
c     &  (msq15_6(j,k)*(sub15_6(qq)+half*sub56_1(gg))
c     & +half*msq56_1v(j,k)*sub56_1v)
c     & +msq15_6(j,k)*(sub16_5(qq)+half*sub56_1(gg))
c     & +half*msq56_1v(j,k)*sub56_1v)
c      msq(6,j,k)=half*xn*(
c     &  msq26_5(j,k)*(sub26_5(qq)+half*sub56_2(gg))
c     & +half*msq56_2v(j,k)*sub56_2v
c     & +msq26_5(j,k)*(sub25_6(qq)+half*sub56_2(gg))
c     & +half*msq56_2v(j,k)*sub56_2v)

      elseif ((k == 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
c      msq(2,j,k)=2._dp*tr*msq25_1(j,-j)*sub25_1(qg)
      if ((isub == 1) .or. (isub == 0)) then
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-(msq15_6(j,k)*sub16_5(qq)+msq15_6(j,k)*sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     &              +msq26_5(j,k)*sub56_2(qq))
      endif

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
c      msq(1,j,k)=2._dp*tr*msq15_2(-k,k)*sub15_2(qg)
      if ((isub == 1) .or. (isub == 0)) then
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq15_6(j,k)*sub16_5(gg)+msq15_6v(j,k)*sub16_5v
     &              +msq15_6(j,k)*sub56_1(qq))
      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn
      endif

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> qb(p5)+q(p6)
c---Hence 15 split multiplies q(15)+g(p2)-->H+q(p6)
c---Hence 25 split multiplies g(p1)+q(p25)-->H+q(p6)
      if ((isub == 2) .or. (isub == 0)) then
      msq(1,j,k)=msq15_2(+5,k)*sub15_2(qg)*2._dp*tr
      msq(2,j,k)=msq25_1(k,+5)*sub25_1(qg)*2._dp*tr
      msq(3,j,k)=msq16_2(-5,k)*sub16_2(qg)*2._dp*tr
      msq(4,j,k)=msq26_1(k,-5)*sub26_1(qg)*2._dp*tr
      endif

      endif

 19   continue
      enddo
      enddo

c--- 4-quark singularities
      do j=-nf,nf
      do k=-nf,nf

      if (((j > 0).and.(k > 0)) .or.
     &    ((j < 0).and.(k < 0))) then
c--q-q or qb-qb
      if (j==k) then
      if ((isub == 2) .or. (isub == 0)) then
      msq(1,j,k)=msq(1,j,k)+half*(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(2,j,k)=msq(2,j,k)+half*(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+half*(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=msq(4,j,k)+half*(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      endif
      else
      if ((isub == 1) .or. (isub == 0)) then
      if (abs(j) == 5) then
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      elseif (abs(k) == 5) then
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      endif
      endif
      endif
      elseif ((j > 0).and.(k < 0)) then
c q-qbar
      if (j==-k) then
      if ((isub == 2) .or. (isub == 0)) then
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      endif
      else
      if ((isub == 1) .or. (isub == 0)) then
      if (abs(j) == 5) then
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      elseif (abs(k) == 5) then
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      endif
      endif
      endif
c--qbar-q
      elseif ((j < 0).and.(k > 0)) then
      if (j==-k) then
      if ((isub == 2) .or. (isub == 0)) then
      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      endif
      else
      if ((isub == 1) .or. (isub == 0)) then
      if (abs(j) == 5) then
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      elseif (abs(k) == 5) then
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      endif
      endif
      endif

      endif


      enddo
      enddo

      return
      end

