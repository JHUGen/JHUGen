      subroutine qqb_w_tndk_g(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + tbar(p5) + f(p6)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ t(p5) + f(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      include 'nwz.f'
      include 'nores.f'
      include 'jetcuts.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     &                 facgg,prop
      real(dp):: qgWqg2,qbgWqbg2,gqbWqbg2,gqWqg2,ggWqqb2,ggWqbq2
      real(dp):: qgWqg2_cs(0:2),qbgWqbg2_cs(0:2),
     &                 gqbWqbg2_cs(0:2),gqWqg2_cs(0:2)
c     &                ,ggWqbq2_cs(0:2),ggWqqb2_cs(0:2),
      real(dp):: sbqWcbq,sqWcq,qsWcq,qsbWcbq,qqbWcsb,qqbWcbs
      real(dp):: s34,msq_gg
      complex(dp):: ampgg_ga(2,2,2,2),ampgg_ag(2,2,2,2)
      common/facgg/facgg
c--- we label the amplitudes by helicity (qqb1 ... qqb4)
c--- and by type of contribution qqb(1) ... qqb(n)
      integer:: i,j,k,ia,ib,ig,it
!$omp threadprivate(/facgg/)

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zero
      enddo
      enddo

c--- Note that no changes need to be made here from the charm
c--- case since the mass is calculated dynamically inside the routines
c--- wqq_sc and w2jetsq_mass

c--- calculate 4-quark contributions
c--- note that these are symmetric under interchange of q and qb
      call wqq_sc(2,5,3,4,6,1,p,qsWcq)
      call wqq_sc(1,5,3,4,6,2,p,sqWcq)

      call wqq_sc(1,5,4,3,6,2,p,sbqWcbq)
      call wqq_sc(2,5,4,3,6,1,p,qsbWcbq)

      call wqq_sc(6,5,3,4,2,1,p,qqbWcsb)
      call wqq_sc(6,5,4,3,1,2,p,qqbWcbs)

      s34=two*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      prop=s34**2/((s34-wmass**2)**2+wmass**2*wwidth**2)
      facgg=V*xn/four*(gwsq/two*gsq)**2*prop

c--- calculate 2-quark, 2-gluon amplitudes

      if     (nwz == +1) then
        call w2jetsq_mass(1,5,4,3,2,6,p,qbgWqbg2)
        call storecs(qbgWqbg2_cs)

        call w2jetsq_mass(2,5,4,3,1,6,p,gqbWqbg2)
        call storecs(gqbWqbg2_cs)

c--- alternative calculation of gg->Wtb piece
        if (nores .eqv. .false.) then
          call BBamps_nores(p,1,2,4,3,5,6,ampgg_ag)
          call BBamps_nores(p,2,1,4,3,5,6,ampgg_ga)
c        call BBamps(p,1,2,3,4,5,6,ampgg_ag)
c        call BBamps(p,2,1,3,4,5,6,ampgg_ga)
          msq_gg=zero
c--- sum over helicities of gluons and massive quarks
          do ib=1,2
          do it=1,2
          do ia=1,2
          do ig=1,2
            msq_gg=msq_gg+xn*cf**2*(
     &            + abs(ampgg_ag(ia,ig,ib,it))**2
     &            + abs(ampgg_ga(ig,ia,ib,it))**2
     &            - one/xn/cf*real(ampgg_ag(ia,ig,ib,it)
     &                     *conjg(ampgg_ga(ig,ia,ib,it))))
          enddo
          enddo
          enddo
          enddo
        endif

      elseif (nwz == -1) then
        call w2jetsq_mass(2,5,3,4,1,6,p,gqWqg2)
        call storecs(gqWqg2_cs)

        call w2jetsq_mass(1,5,3,4,2,6,p,qgWqg2)
        call storecs(qgWqg2_cs)

c--- alternative calculation of gg->Wtb piece
        if (nores .eqv. .false.) then
          call BBamps_nores(p,1,2,3,4,5,6,ampgg_ag)
          call BBamps_nores(p,2,1,3,4,5,6,ampgg_ga)
c        call BBamps(p,1,2,3,4,5,6,ampgg_ag)
c        call BBamps(p,2,1,3,4,5,6,ampgg_ga)
          msq_gg=zero
c--- sum over helicities of gluons and massive quarks
          do ib=1,2
          do it=1,2
          do ia=1,2
          do ig=1,2
            msq_gg=msq_gg+xn*cf**2*(
     &            + abs(ampgg_ag(ia,ig,ib,it))**2
     &            + abs(ampgg_ga(ig,ia,ib,it))**2
     &            - one/xn/cf*real(ampgg_ag(ia,ig,ib,it)
     &                     *conjg(ampgg_ga(ig,ia,ib,it))))
          enddo
          enddo
          enddo
          enddo
        endif

      endif
c--- end of alternative calculation

c-- veto b-jet contribution if doing subtraction and pt(b)>ptbjetmin GeV
      if (sqrt(p(6,1)**2+p(6,2)**2) > ptbjetmin) then
       msq_gg=zero
      endif

      ggWqqb2=avegg*gsq**2*gwsq**2*msq_gg*(prop/s34**2)
      ggWqbq2=avegg*gsq**2*gwsq**2*msq_gg*(prop/s34**2)

c        call w2jetsq_mass(6,5,3,4,1,2,p,ggWqqb2)
c        call storecs(ggWqqb2_cs)
      do i=0,2
        gqWqg2_cs(i)  = aveqg*facgg*gqWqg2_cs(i)
        qgWqg2_cs(i)  = aveqg*facgg*qgWqg2_cs(i)
        gqbWqbg2_cs(i)  = aveqg*facgg*gqbWqbg2_cs(i)
        qbgWqbg2_cs(i)  = aveqg*facgg*qbgWqbg2_cs(i)
c        ggWqqb2_cs(i) = avegg*facgg*ggWqqb2_cs(i)
      enddo
      gqWqg2  = gqWqg2_cs(1)  +gqWqg2_cs(2)  +gqWqg2_cs(0)
      qgWqg2  = qgWqg2_cs(1)  +qgWqg2_cs(2)  +qgWqg2_cs(0)
      gqbWqbg2  = gqbWqbg2_cs(1)  +gqbWqbg2_cs(2)  +gqbWqbg2_cs(0)
      qbgWqbg2  = qbgWqbg2_cs(1)  +qbgWqbg2_cs(2)  +qbgWqbg2_cs(0)
c      ggWqqb2 = ggWqqb2_cs(1) +ggWqqb2_cs(2) +ggWqqb2_cs(0)

      do j=-nf,nf
      do k=-nf,nf

c--- skip contributions with 2 b-quarks in the initial state
      if ((abs(j) == 5) .and. (abs(k) == 5)) goto 99

c--- 2-quark, 2-gluon contribution to matrix elements
      if     ((j == +5) .and. (k == 0) .and. (nwz == -1)) then
          msq(j,k)=qgWqg2
          do i=0,2
            msq_cs(i,j,k)=qgWqg2_cs(i)
          enddo
      elseif ((j == -5) .and. (k == 0) .and. (nwz == +1)) then
          msq(j,k)=qbgWqbg2
          do i=0,2
            msq_cs(i,j,k)=qbgWqbg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k == +5) .and. (nwz == -1)) then
          msq(j,k)=gqWqg2
          do i=0,2
            msq_cs(i,j,k)=gqWqg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k == -5) .and. (nwz == +1)) then
          msq(j,k)=gqbWqbg2
          do i=0,2
            msq_cs(i,j,k)=gqbWqbg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k == 0)) then
          if (nores .eqv. .false.) then
          if     (nwz == +1) then
            msq(j,k)=ggWqbq2
c            do i=0,2
c              msq_cs(i,j,k)=ggWqbq2_cs(i)
c            enddo
          elseif (nwz == -1) then
            msq(j,k)=ggWqqb2
c            do i=0,2
c              msq_cs(i,j,k)=ggWqqb2_cs(i)
c            enddo
          endif
          endif
      endif

c--- 4-quark contribution to matrix elements
      if     ((j == +5) .and. (k .ne. 0) .and. (nwz == -1)) then
        msq(j,k)=sqWcq
      elseif ((j .ne. 0) .and. (k == +5) .and. (nwz == -1)) then
        msq(j,k)=qsWcq
      elseif ((j == -5) .and. (k .ne. 0) .and. (nwz == +1)) then
        msq(j,k)=sbqWcbq
      elseif ((j .ne. 0) .and. (k == -5) .and. (nwz == +1)) then
        msq(j,k)=qsbWcbq
      elseif ((j .ne. 0) .and. (k .ne. 0)) then

        if     ((j == -k) .and. (nwz == +1)) then
c          if (nores .eqv. .false.) msq(j,k)=qqbWcbs
        elseif ((j == -k) .and. (nwz == -1)) then
c          if (nores .eqv. .false.) msq(j,k)=qqbWcsb
        endif
      endif

   99 continue

      enddo
      enddo

      return
      end

