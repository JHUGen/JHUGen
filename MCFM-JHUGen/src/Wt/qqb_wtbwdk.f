      subroutine qqb_wtbwdk(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     W+t+b production, modified form of W+t radiation in production   *
*                                                                      *
*    Matrix element squared and averaged over initial colours and spins*
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567) + b(p8)                       *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567) + b(p8)                       *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zprods_com.f'
      integer:: ia,ig,ib,j,k
      logical:: nearmt,ttsubtract
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),fac,dot,prop,
     & msq_gg,mWb,q(mxpart,4),twotDg,propw
      complex(dp):: ampgg_ga(2,2,2),ampgg_ag(2,2,2),
     &               ampgg_ag_mb(2,2,2,2),ampgg_ga_mb(2,2,2,2),
     &               ampld(2)
!      common/nearmt/nearmt

      msq(:,:)=zip


c--- veto this contribution if we're doing the subtraction and pt(b~)>25 GeV
c      if (runstring(1:3) == 'sub') then
c        if (sqrt(p(8,1)**2+p(8,2)**2) > 60._dp) then
cc          return
c        endif
c      endif

c--- variable to determine whether or not to subtract tt contribution
      ttsubtract=.false.
c      if (runstring(1:3) == 'sub') then
c        ttsubtract=.true.
c        nearmt=.false.
c      else
c--- first check the invariant mass of the W and the potential b-jet,
c--- if it's close to the top mass then the contribution from
c--- gg and qqb initial states should be zero
        mWb=sqrt(two*dot(p,3,4)+two*dot(p,3,8)+two*dot(p,4,8)
     &           +dot(p,8,8))
        if (abs(mWb-mt) < 10._dp*twidth) then
          nearmt=.true.
        else
          nearmt=.false.
        endif
c      endif

      fac=gsq**2*gwsq**4

      propw=(two*dot(p,3,4)-wmass**2)**2+(wmass*wwidth)**2
      prop=propw*((two*dot(p,5,6)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((two*dot(p,5,6)+two*dot(p,5,7)+two*dot(p,6,7)-mt**2)**2
     &         +(mt*twidth)**2)

c--- calculate amplitudes for g+g in initial state
      if (nearmt) then
        do ia=1,2
        do ig=1,2
        do ib=1,2
        ampgg_ag(ia,ig,ib)=0._dp
        ampgg_ga(ia,ig,ib)=0._dp
        enddo
        enddo
        enddo
      else
c---- check the massive result against the massless one
c        call gs_wt_prog_nores(p,1,8,3,4,5,6,7,2,ampgg_ag_mb0)
c        call gs_wt_prog_nores(p,2,8,3,4,5,6,7,1,ampgg_ga_mb0)
c        msq_gg=0._dp
c--- sum over helicities of gluons
c        do ia=1,2
c        do ig=1,2
c        msq_gg=msq_gg+xn*cf**2*(
c     &  +abs(ampgg_ag_mb0(ia,ig))**2 + abs(ampgg_ga_mb0(ig,ia))**2
c     & -one/xn/cf*real(ampgg_ag_mb0(ia,ig)*conjg(ampgg_ga_mb0(ig,ia))))
c        enddo
c        enddo
c        write(6,*) 'massless ',msq_gg
c--- end of check

c--- modifications to enable calculations with a massive b~ quark
c        do ia=1,2
c        do ig=1,2
c        write(6,*) 'ampgg_ag',ia,ig,abs(ampgg_ag(ia,ig))**2/prop
c        write(6,*) 'ampgg_ga',ia,ig,abs(ampgg_ga(ia,ig))**2/prop
c        enddo
c        enddo
        do k=1,4
          q(1,k)=p(1,k)
          q(2,k)=p(2,k)
          q(3,k)=p(3,k)
          q(4,k)=p(4,k)
          q(5,k)=p(5,k)+p(6,k)+p(7,k)
          q(6,k)=p(8,k)
        enddo
c        call BBamps(q,1,2,3,4,5,6,ampgg_ag_mb)
c        call BBamps(q,2,1,3,4,5,6,ampgg_ga_mb)
c--- contribution from non-resonant diagrams only
        call BBamps_nores(q,1,2,3,4,5,6,ampgg_ag_mb)
        call BBamps_nores(q,2,1,3,4,5,6,ampgg_ga_mb)
        twotDg=two*(dot(p,5,4)+dot(p,6,4)+dot(p,7,4))
        do k=1,4
        do j=1,8
          q(j,k)=p(j,k)
        enddo
          q(9,k)=p(5,k)+p(6,k)+p(7,k)-p(4,k)*mt**2/twotDg
        enddo
        call spinoru(9,q,za,zb)
        call wampd(mt,twidth,4,5,6,7,9,ampld)
        do ia=1,2
        do ig=1,2
        do ib=1,2
        ampgg_ag(ia,ig,ib)=
     &   (ampgg_ag_mb(ia,ig,ib,1)*ampld(1)
     &   +ampgg_ag_mb(ia,ig,ib,2)*ampld(2))*sqrt(prop/propw)
        ampgg_ga(ia,ig,ib)=
     &   (ampgg_ga_mb(ia,ig,ib,1)*ampld(1)
     &   +ampgg_ga_mb(ia,ig,ib,2)*ampld(2))*sqrt(prop/propw)
        enddo
        enddo
        enddo

      endif

      msq_gg=0._dp
c--- sum over helicities of gluons, b-quark
      do ia=1,2
      do ig=1,2
      do ib=1,2
        msq_gg=msq_gg+xn*cf**2*(
     &  + abs(ampgg_ag(ia,ig,ib))**2 + abs(ampgg_ga(ig,ia,ib))**2
     &  - one/xn/cf*real(ampgg_ag(ia,ig,ib)*conjg(ampgg_ga(ig,ia,ib))))
      enddo
      enddo
      enddo

c      write(6,*) 'massive  ',msq_gg
c      pause

      msq_gg=msq_gg/prop

c--- subtract the gg->tt and qqb->tt contributions if necessary
      if (ttsubtract) then
        do j=1,4
          q(1,j)=p(1,j)
          q(2,j)=p(2,j)
          q(3,j)=p(5,j)
          q(4,j)=p(6,j)
          q(5,j)=p(7,j)
          q(6,j)=p(8,j)
          q(7,j)=p(3,j)
          q(8,j)=p(4,j)
        enddo
c        call qqb_ttb(q,msqa)
c        msq_gg=msq_gg-msqa(0,0)/avegg/fac
      endif

c--- calculate matrix elements for q+q~ initial state
c      if (nearmt) then
c        msq_qqb=0._dp
c      else
c        call qb_wtq(p,1,8,3,4,5,6,7,2,msq_qqb)
c      endif

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msq(j,k)=0._dp
      if ((j == 0) .and. (k == 0)) then
          msq(j,k)=avegg*fac*msq_gg
      endif
      enddo
      enddo

c      do j=1,4
c          msq(j,-j)=msq_qqb
c          msq(-j,j)=msq_qqb
c      enddo

      return
      end

