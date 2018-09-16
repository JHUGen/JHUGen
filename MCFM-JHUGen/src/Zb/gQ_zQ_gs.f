      subroutine gQ_zQ_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Authors: R.K. Ellis and John Campbell                             *
*    July, 2003.                                                       *
*    Matrix element SUBTRACTION for Z + heavy quark                    *
*    (of flavour "flav") production                                    *
*    averaged over initial colours and spins                           *
*     g(-p1)+Q(-p2)-->Z^+(l(p3)+a(p4))+Q(p5)                           *
*                                                                      *
*    isub=1 : particle 6 is a gluon or light quark                     *
*    isub=2 : particle 6 is also a heavy quark                         *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'heavyflav.f'
      integer:: j,k,nd,isub
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
     & sub25_6v
      common/isub/isub
      external gQ_zQ,qqb_z_gvec

      ndmax=6

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & gQ_zQ,qqb_z_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & gQ_zQ,qqb_z_gvec)
      call dips(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & gQ_zQ,qqb_z_gvec)
      call dips(4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & gQ_zQ,qqb_z_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & gQ_zQ,qqb_z_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     & gQ_zQ,qqb_z_gvec)
      call dips(5,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & gQ_zQ,qqb_z_gvec)

      call dips(6,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & gQ_zQ,qqb_z_gvec)
      call dips(6,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     & gQ_zQ,qqb_z_gvec)
      call dips(6,p,2,5,6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     & gQ_zQ,qqb_z_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

c--- subtraction terms for 2-quark, 2-gluon pieces
      do j=-flav,flav,flav
      do k=-flav,flav,flav
      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

      if     ((k == 0).and.(j.ne.0).and.(isub==1)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2._dp*tr*msq25_1(j,-j)*sub25_1(qg)
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-(msq16_5(j,k)*sub16_5(qq)+msq16_5(j,k)*sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     &              +msq26_5(j,k)*sub56_2(qq))

      elseif ((j == 0).and.(k.ne.0).and.(isub==1)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2._dp*tr*msq15_2(-k,k)*sub15_2(qg)
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     &              +msq15_6(j,k)*sub56_1(qq))
      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn

      elseif ((j == 0).and.(k == 0).and.(isub==2)) then
c--- g-g case (real process is g(p1)+g(p2) --> qb(p5)+q(p6)
c--- Hence 15 split multiplies q(15)+g(p2)-->Z+q(p6)
c--- Hence 25 split multiplies g(p1)+q(p25)-->Z+q(p6)
      msq(1,j,k)=(msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
     &           +msq15_2(+4,k)+msq15_2(+5,k))*sub15_2(qg)*2._dp*tr
      msq(2,j,k)=(msq25_1(k,+1)+msq25_1(k,+2)+msq25_1(k,+3)
     &           +msq25_1(k,+4)+msq25_1(k,+5))*sub25_1(qg)*2._dp*tr
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     &           +msq16_2(-2,k)+msq16_2(-1,k))*sub16_2(qg)*2._dp*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
     &           +msq26_1(k,-2)+msq26_1(k,-1))*sub26_1(qg)*2._dp*tr

      endif

 19   continue
      enddo
      enddo

c--- subtraction terms for 4-quark pieces
      do j=-nf,nf
      do k=-nf,nf

      if (((j > 0).and.(k > 0)) .or.
     &    ((j < 0).and.(k < 0))) then
c--- q-q or qb-qb
        if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 20
        if     ((j==k) .and. (isub==2)) then
c--- don't include the Q+Q->Z+Q+Q case here
c        goto 20
        msq(1,j,k)=msq(1,j,k)+0.5_dp*(xn-1._dp/xn)
     &    *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
        msq(2,j,k)=msq(2,j,k)+0.5_dp*(xn-1._dp/xn)
     &    *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
        msq(3,j,k)=msq(3,j,k)+0.5_dp*(xn-1._dp/xn)
     &    *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
        msq(4,j,k)=msq(4,j,k)+0.5_dp*(xn-1._dp/xn)
     &    *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
        elseif ((j.ne.k) .and. (isub==1)) then
          if (abs(k) == flav) then
            msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
          else
            msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
          endif
        endif

      elseif ((j > 0).and.(k < 0)) then
c--- q-qbar
        if (j==-k) then
          if ((abs(j) == flav) .and. (isub == 2)) then
c--- don't include the Q+Qb->Z+Q+Qb case here
          goto 20
c          msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
c     &      *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c          msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
c     &      *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
          endif
        else
          if ((j .ne. flav) .and. (k .ne. -flav)) goto 20
          if (isub == 2) goto 20
          if (k == -flav) then
            msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
          else
            msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
          endif
        endif

      elseif ((j < 0).and.(k > 0)) then
c--- qbar-q
        if (j==-k) then
          if ((abs(j) == flav) .and. (isub == 2)) then
c--- don't include the Qb+Q->Z+Qb+Q case here
          goto 20
c          msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
c     &      *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c          msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
c     &      *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
          endif
        else
          if ((j .ne. -flav) .and. (k .ne. flav)) goto 20
          if (isub == 2) goto 20
          if (j == -flav) then
            msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
          else
            msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
          endif
        endif

      endif

   20 continue
      enddo
      enddo

      return
      end

