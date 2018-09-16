      subroutine qqb_zbjet_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J. Campbell                                              *
*     July, 2005.                                                      *
************************************************************************
c--- Matrix element squared averaged over initial colors and spins
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
      include 'cplx.h'
      include 'heavyflav.f'
      include 'masses.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: mmsq_gg(2,2),
     & mmsq_qg(2,2),mmsq_gq(2,2),mmsq_gqb(2,2),mmsq_qbg(2,2)
      real(dp):: fac
      real(dp):: msq_Hq(2,2),msq_qH(2,2),
     ,                 msq_Hbqb(2,2),msq_qbHb(2,2),
     &                 msq_Hqb(2,2),msq_qHb(2,2),
     &                 msq_Hbq(2,2),msq_qbH(2,2),
     &                 msq_Hg(2,2),msq_qg(2,2),
     &                 msq_Hbg(2,2),msq_qbg(2,2),
     &                 msq_gHb(2,2),msq_gqb(2,2),
     &                 msq_gH(2,2),msq_gq(2,2),
     &                 msq_qqb(2,2),msq_qbq(2,2)
      complex(dp):: prop
      integer:: isub
      common/isub/isub
      integer, parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer, parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c--- initialize the matrix element squared
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

C---call spinor routine and load common block twopij
      call spinoru(7,p,za,zb)

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

************************************************************************
*     Calculate contributions from the QQGGG matrix elements            *
************************************************************************

      if     (isub == 1) then
c-- matrix elements for qg -> qg
        call xzqqggg(1,2,6,7,5,3,4,mmsq_qg)
c-- matrix elements for gq -> gq
        call xzqqggg(2,1,6,7,5,3,4,mmsq_gq)
c-- matrix elements for gqb -> gqb
        call xzqqggg(5,1,6,7,2,3,4,mmsq_gqb)
c-- matrix elements for qbg -> qbg
        call xzqqggg(5,2,6,7,1,3,4,mmsq_qbg)
      elseif (isub == 2) then
c-- matrix elements for gg -> qbq
        call xzqqggg(5,1,2,7,6,3,4,mmsq_gg)
      endif

      do j=-nflav,nflav,nflav
      do k=-nflav,nflav,nflav

      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 99
c--- do not allow abs(j) = flav and abs(k) = flav

c--- Q-G contribution (isub=1 only)
      if     ((isub == 1) .and. (j == +flav) .and. (k == 0)) then
        msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*mmsq_qg(1,1)
     &           +abs(Q(j)*q1+R(j)*r1*prop)**2*mmsq_qg(2,2)
     &           +abs(Q(j)*q1+L(j)*r1*prop)**2*mmsq_qg(1,2)
     &           +abs(Q(j)*q1+R(j)*l1*prop)**2*mmsq_qg(2,1)
c---Statistical factor
        msq(j,k)=half*aveqg/avegg*msq(j,k)
c--- Qb-G contribution (isub=1 only)
      elseif ((isub == 1) .and. (j == -flav) .and. (k == 0)) then
        msq(j,k)=+abs(Q(-j)*q1+L(-j)*l1*prop)**2*mmsq_qbg(1,1)
     &           +abs(Q(-j)*q1+R(-j)*r1*prop)**2*mmsq_qbg(2,2)
     &           +abs(Q(-j)*q1+L(-j)*r1*prop)**2*mmsq_qbg(1,2)
     &           +abs(Q(-j)*q1+R(-j)*l1*prop)**2*mmsq_qbg(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
c--- G-Q contribution (isub=1 only)
      elseif ((isub == 1) .and. (j == 0) .and. (k == +flav)) then
        msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*mmsq_gq(1,1)
     &           +abs(Q(k)*q1+R(k)*r1*prop)**2*mmsq_gq(2,2)
     &           +abs(Q(k)*q1+L(k)*r1*prop)**2*mmsq_gq(1,2)
     &           +abs(Q(k)*q1+R(k)*l1*prop)**2*mmsq_gq(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
c--- G-Qb contribution (isub=1 only)
      elseif ((isub == 1) .and. (j == 0) .and. (k == -flav)) then
        msq(j,k)=+abs(Q(-k)*q1+L(-k)*l1*prop)**2*mmsq_gqb(1,1)
     &           +abs(Q(-k)*q1+R(-k)*r1*prop)**2*mmsq_gqb(2,2)
     &           +abs(Q(-k)*q1+L(-k)*r1*prop)**2*mmsq_gqb(1,2)
     &           +abs(Q(-k)*q1+R(-k)*l1*prop)**2*mmsq_gqb(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
c--- G-G contribution (isub=2 only)
      elseif ((isub == 2) .and. (j == 0) .and. (k == 0)) then
        msq(j,k)=
     &   +abs(Q(nflav)*q1+L(nflav)*l1*prop)**2*mmsq_gg(1,1)
     &   +abs(Q(nflav)*q1+R(nflav)*r1*prop)**2*mmsq_gg(2,2)
     &   +abs(Q(nflav)*q1+L(nflav)*r1*prop)**2*mmsq_gg(1,2)
     &   +abs(Q(nflav)*q1+R(nflav)*l1*prop)**2*mmsq_gg(2,1)
      endif

   99 continue
      enddo
      enddo


************************************************************************
*     Calculate contributions from the QQqqG matrix elements            *
************************************************************************

c--- note the factor of 4._dp*xw**2 relative to wbb
      fac=4._dp*gsq**3*esq**2*8._dp
c--- extra factor of 2**3=8 to compensate for Ta normalization

c--- the notation in msq_ZqqQQg_noid is:
c---  msq_ZqqQQg_noid(i1,i2,i3,i4,i5,...) where the basic process
c---   is   0 --> Q(i1)+Qbar(i2)+q(i3)+qbar(i4)+g(i5)

c--- these amplitudes are only used for isub=1
      if (isub == 1) then
c--- calculate H-q and q-Q amplitudes
      call msq_ZqqQQg_noid(5,1,6,2,7,4,3,msq_Hq)
      call msq_ZqqQQg_noid(5,2,6,1,7,4,3,msq_qH)
c--- calculate Hb-qb and qb-Hb amplitudes
      call msq_ZqqQQg_noid(1,5,2,6,7,4,3,msq_Hbqb)
      call msq_ZqqQQg_noid(2,5,1,6,7,4,3,msq_qbHb)
c--- calculate H-qb and q-Qb amplitudes
      call msq_ZqqQQg_noid(5,1,2,6,7,4,3,msq_Hqb)
      call msq_ZqqQQg_noid(2,5,6,1,7,4,3,msq_qHb)
c--- calculate Hb-q and qb-H amplitudes
      call msq_ZqqQQg_noid(1,5,6,2,7,4,3,msq_Hbq)
      call msq_ZqqQQg_noid(5,2,1,6,7,4,3,msq_qbH)
      endif

c--- these amplitudes are only used for isub=2 (ADDED 12/8/05)
      if (isub == 2) then
      call msq_ZqqQQg_noid(5,6,2,1,7,4,3,msq_qqb)
      call msq_ZqqQQg_noid(5,6,1,2,7,4,3,msq_qbq)
      endif

c--- calcule g-H and g-q amplitudes
      if (isub == 1) call msq_ZqqQQg_noid(5,2,7,6,1,4,3,msq_gH)
      if (isub == 2) call msq_ZqqQQg_noid(5,6,7,2,1,4,3,msq_gq)
c--- calcule H-g and q-g amplitudes
      if (isub == 1) call msq_ZqqQQg_noid(5,1,7,6,2,4,3,msq_Hg)
      if (isub == 2) call msq_ZqqQQg_noid(5,6,7,1,2,4,3,msq_qg)
c--- calcule g-Hb and g-qb amplitudes
      if (isub == 1) call msq_ZqqQQg_noid(2,5,7,6,1,4,3,msq_gHb)
      if (isub == 2) call msq_ZqqQQg_noid(5,6,2,7,1,4,3,msq_gqb)
c--- calcule Hb-g and qb-g amplitudes
      if (isub == 1) call msq_ZqqQQg_noid(1,5,7,6,2,4,3,msq_Hbg)
      if (isub == 2) call msq_ZqqQQg_noid(5,6,1,7,2,4,3,msq_qbg)

      do j=-nflav,nflav
      do k=-nflav,nflav

      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 89
c--- do not allow abs(j) = flav and abs(k) = flav

c      if ((abs(j) > 0) .and. (abs(k) > 0) .and.
c     &    (abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 89
c--- if both are (anti-)quarks, one of them should be a heavy quark

      if     ((isub == 2) .and. (j > 0) .and. (k == -j)
     &        .and. (j .ne. +flav)) then
        msq(j,k)=fac*aveqq*msq_qqb(jj(j),-kk(k))
      elseif ((isub == 2) .and. (j < 0) .and. (k == -j)
     &        .and. (j .ne. -flav)) then
        msq(j,k)=fac*aveqq*msq_qbq(-jj(j),kk(k))
      endif

      if     ((isub == 1) .and. (j > 0) .and. (k > 0)) then
c--- Q-q/q-Q contribution (isub=1 only)
        if     (j == +flav) then
          msq(j,k)=fac*aveqq*msq_Hq(jj(j),kk(k))
        elseif (k == +flav) then
          msq(j,k)=fac*aveqq*msq_qH(kk(k),jj(j))
        endif

      elseif ((isub == 1) .and. (j < 0) .and. (k < 0)) then
c--- Qb-qb/qb-Qb contribution (isub=1 only)
        if     (j == -flav) then
          msq(j,k)=fac*aveqq*msq_Hbqb(-jj(j),-kk(k))
        elseif (k == -flav) then
          msq(j,k)=fac*aveqq*msq_qbHb(-kk(k),-jj(j))
        endif

      elseif ((isub == 1) .and. (j > 0) .and. (k < 0)) then
c--- Q-qb/q-Qb contribution (isub=1 only)
        if     (j == +flav) then
          msq(j,k)=fac*aveqq*msq_Hqb(jj(j),-kk(k))
        elseif (k == -flav) then
          msq(j,k)=fac*aveqq*msq_qHb(-kk(k),jj(j))
        endif

      elseif ((isub == 1) .and. (j < 0) .and. (k > 0)) then
c--- Qb-q/qb-Q contribution (isub=1 only)
        if     (j == -flav) then
          msq(j,k)=fac*aveqq*msq_Hbq(-jj(j),kk(k))
        elseif (k == +flav) then
          msq(j,k)=fac*aveqq*msq_qbH(kk(k),-jj(j))
        endif

      elseif ((j > 0) .and. (k == 0)) then
c--- Q-g contribution (isub=1 only)
        if     ((isub == 1) .and. (j == +flav)) then
c--- other quark line is (d,s) and either (u) (flav=4) or (u,c) (flav=5)
          msq(j,k)=msq(j,k)+fac*aveqg*(
     &      +2._dp*msq_Hg(jj(j),1)
     &      +real(flav-3,dp)*msq_Hg(jj(j),2))
c--- q-g contribution (isub=2 only)
        elseif ((isub == 2) .and. (j .ne. +flav)) then
c--- other quark line is (Q)
          msq(j,k)=msq(j,k)+fac*aveqg*msq_qg(kk(flav),jj(j))
        endif

      elseif ((j < 0) .and. (k == 0)) then
c--- Qb-g contribution (isub=1 only)
        if     ((isub == 1) .and. (j == -flav)) then
c--- other quark line is (d,s) and either (u) (flav=4) or (u,c) (flav=5)
          msq(j,k)=msq(j,k)+fac*aveqg*(
     &      +2._dp*msq_Hbg(-jj(j),1)
     &      +real(flav-3,dp)*msq_Hbg(-jj(j),2))
c--- qb-g contribution (isub=2 only)
        elseif ((isub == 2) .and. (j .ne. -flav)) then
c--- other quark line is (Q)
          msq(j,k)=msq(j,k)+fac*aveqg*msq_qbg(kk(flav),-jj(j))
        endif

      elseif ((j == 0) .and. (k > 0)) then
c--- g-Q contribution (isub=1 only)
        if     ((isub == 1) .and. (k == +flav)) then
c--- other quark line is (d,s) and either (u) (flav=4) or (u,c) (flav=5)
          msq(j,k)=msq(j,k)+fac*aveqg*(
     &      +2._dp*msq_gH(kk(k),1)
     &      +real(flav-3,dp)*msq_gH(kk(k),2))
c--- g-q contribution (isub=2 only)
        elseif ((isub == 2) .and. (k .ne. +flav)) then
c--- other quark line is (Q)
          msq(j,k)=msq(j,k)+fac*aveqg*msq_gq(jj(flav),kk(k))
        endif

      elseif ((j == 0) .and. (k < 0)) then
c--- g-Qb contribution (isub=1 only)
        if     ((isub == 1) .and. (k == -flav)) then
c--- other quark line is (d,s) and either (u) (flav=4) or (u,c) (flav=5)
          msq(j,k)=msq(j,k)+fac*aveqg*(
     &      +2._dp*msq_gHb(-kk(k),1)
     &      +real(flav-3,dp)*msq_gHb(-kk(k),2))
c--- g-qb contribution (isub=2 only)
        elseif ((isub == 2) .and. (k .ne. -flav)) then
c--- other quark line is (Q)
          msq(j,k)=msq(j,k)+fac*aveqg*msq_gqb(jj(flav),-kk(k))
        endif
      endif

   89 continue
      enddo
      enddo

      return
      end
