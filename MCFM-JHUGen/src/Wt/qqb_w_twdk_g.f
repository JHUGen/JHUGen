      subroutine qqb_w_twdk_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     Real correction to W+t, radiation in production                  *
*                                                                      *
*    Matrix element squared and averaged over initial colours and spins*
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567) + f(p8)                       *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567) + f(p8)                       *
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
      include 'nores.f'
      include 'nwz.f'
      include 'jetcuts.f'
      integer:: ia,ig,j,k,i3,i4,i5,i6,iq
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),fac,dot,prop,
     & msq_gq,msq_qg,msq_gg,msq_qb,msq_bq,msq_ab,msq_ba
      complex(dp):: ampgq_ga(2,2),ampgq_ag(2,2),
     &               ampqg_ga(2,2),ampqg_ag(2,2),
     &               ampgg_ga(2,2),ampgg_ag(2,2)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zip
      enddo
      enddo

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        i3=3
        i4=4
        i5=5
        i6=6
        iq=+1 ! quark initial state
      elseif (nwz == +1) then
        i3=4
        i4=3
        i5=6
        i6=5
        iq=-1 ! antiquark initial state
      else
        write(6,*) 'Error in qqb_w_twdk_g, nwz is not +1 or -1 :   ',nwz
        stop
      endif

      fac=gsq**2*gwsq**4

      prop=(two*dot(p,3,4)-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((two*dot(p,5,6)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((two*dot(p,5,6)+two*dot(p,5,7)+two*dot(p,6,7)-mt**2)**2
     &         +(mt*twidth)**2)

c--- calculate amplitudes for g+q in initial state
      call gs_wt_prog(mt,twidth,p,1,2,i3,i4,i5,i6,7,8,ampgq_ag)
      call gs_wt_prog(mt,twidth,p,8,2,i3,i4,i5,i6,7,1,ampgq_ga)
c--- calculate amplitudes for q+g in initial state
      call gs_wt_prog(mt,twidth,p,2,1,i3,i4,i5,i6,7,8,ampqg_ag)
      call gs_wt_prog(mt,twidth,p,8,1,i3,i4,i5,i6,7,2,ampqg_ga)
c--- calculate amplitudes for g+g in initial state
c      call gs_wt_prog(p,1,8,i3,i4,5,6,7,2,ampgg_ag)
c      call gs_wt_prog(p,2,8,i3,i4,5,6,7,1,ampgg_ga)
c--- The two lines commented out below calculate the contribution
c--- from non-resonant diagrams only (gauge-dependent)
      if (nores .eqv. .false.) then
      call gs_wt_prog_nores(p,1,8,i3,i4,i5,i6,7,2,ampgg_ag)
      call gs_wt_prog_nores(p,2,8,i3,i4,i5,i6,7,1,ampgg_ga)
      endif

      msq_gq=zip
      msq_qg=zip
      msq_gg=zip
c--- sum over helicities of gluons
      do ia=1,2
      do ig=1,2
        msq_gq=msq_gq+xn*cf**2*(
     &        + abs(ampgq_ag(ia,ig))**2 + abs(ampgq_ga(ig,ia))**2
     &        - one/xn/cf*real(ampgq_ag(ia,ig)*conjg(ampgq_ga(ig,ia))))
        msq_qg=msq_qg+xn*cf**2*(
     &        + abs(ampqg_ag(ia,ig))**2 + abs(ampqg_ga(ig,ia))**2
     &        - one/xn/cf*real(ampqg_ag(ia,ig)*conjg(ampqg_ga(ig,ia))))
      if (nores .eqv. .false.) then
        msq_gg=msq_gg+xn*cf**2*(
     &        + abs(ampgg_ag(ia,ig))**2 + abs(ampgg_ga(ig,ia))**2
     &        - one/xn/cf*real(ampgg_ag(ia,ig)*conjg(ampgg_ga(ig,ia))))
      endif
      enddo
      enddo

c-- veto b-jet contribution if doing subtraction and pt(b)>ptbjetmin GeV
      if (sqrt(p(8,1)**2+p(8,2)**2) > ptbjetmin) then
        msq_gg=zip
      endif

      msq_gq=msq_gq/prop
      msq_qg=msq_qg/prop
      msq_gg=msq_gg/prop

c--- calculate matrix elements for q+b initial state
      call qb_wtq(mt,twidth,p,1,2,i3,i4,i5,i6,7,8,msq_qb)
      call qb_wtq(mt,twidth,p,2,1,i3,i4,i5,i6,7,8,msq_bq)
      call qb_wtq(mt,twidth,p,8,2,i3,i4,i5,i6,7,1,msq_ab)
      call qb_wtq(mt,twidth,p,8,1,i3,i4,i5,i6,7,2,msq_ba)
c--- calculate matrix elements for q+q~ initial state
c      if (nores .eqv. .false.) then
c      call qb_wtq(p,1,8,3,4,5,6,7,2,msq_qqb)
c      endif

c--- subtract the gg->tt and qqb->tt contributions if necessary
c      if (ttsubtract) then
c        do j=1,4
c          q(1,j)=p(1,j)
c          q(2,j)=p(2,j)
c          q(3,j)=p(5,j)
c          q(4,j)=p(6,j)
c          q(5,j)=p(7,j)
c          q(6,j)=p(8,j)
c          q(7,j)=p(3,j)
c          q(8,j)=p(4,j)
c        enddo
c        call qqb_ttb(q,msqa)
c        msq_gg=msq_gg-msqa(0,0)/avegg/fac
c        msq_qqb=msq_qqb-msqa(1,-1)
c      endif

c--- if we're doing the window method, we must have mb>0 and the
c--- gg, qqb contributions shouldn't be included here
c      if (runstring(1:3) .ne. 'sub') then
c        msq_gg=zip
c        msq_qqb=zip
c      endif

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msq(j,k)=zip
      if     ((j == +5*iq) .and. (k == 0)) then
          msq(j,k)=aveqg*fac*msq_qg
      elseif ((j == 0) .and. (k == +5*iq)) then
          msq(j,k)=aveqg*fac*msq_gq
      elseif ((j == 0) .and. (k == 0)) then
          if (nores .eqv. .false.) msq(j,k)=avegg*fac*msq_gg
      endif
      enddo
      enddo

      do j=1,4
          msq(+j,+5*iq)=msq_qb
          msq(-j,+5*iq)=msq_ab
          msq(+5*iq,+j)=msq_bq
          msq(+5*iq,-j)=msq_ba

c          if (nores .eqv. .false.) msq(j,-j)=msq_qqb
c          if (nores .eqv. .false.) msq(-j,j)=msq_qqb
      enddo

      return
      end

