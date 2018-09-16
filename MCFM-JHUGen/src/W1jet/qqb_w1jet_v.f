      subroutine qqb_w1jet_v(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W + jet production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->W^+(nu(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,sw,prop,virt5,subuv,
     & qqbWg,qbqWg,qgWq,gqWq,qbgWqb,gqbWqb

      integer,parameter::
     & iqqbg(5)=(/1,2,3,4,5/),iqgq(5)=(/1,5,3,4,2/),
     & igqq(5)=(/2,5,3,4,1/),iqbqg(5)=(/2,1,3,4,5/),
     & iqbgqb(5)=(/5,1,3,4,2/),igqbqb(5)=(/5,2,3,4,1/)

c--set msq=0 to initialize
      msq(:,:)=zip
      scheme='dred'

c-- if Gflag=.false. then only the endpoint contributions from the
c-- 4-quark diagrams are included, ie. no pole subtraction for this
c-- piece. Therefore return 0.
c      if (Gflag .eqv. .false.) return

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities
c      if ((abs(s(1,5)) < cutoff).or.(abs(s(2,5)) < cutoff)) return

c--- calculate lowest order
      call qqb_w_g(p,msq0)

c--- UV counterterm contains the finite renormalization to arrive
c--- at MS bar scheme.
      subuv=ason2pi*xn*(epinv*(eleven-two*real(nflav,dp)/xn)-one)/six

c--- calculate propagator
      sw=s(3,4)
      prop=sw**2/((sw-wmass**2)**2+(wmass*wwidth)**2)

      fac=two*cf*xnsq*gwsq**2*gsq*prop

      qqbWg=aveqq*fac*virt5(iqqbg,za,zb)
      qbqWg=aveqq*fac*virt5(iqbqg,za,zb)
      gqWq=aveqg*fac*virt5(igqq,za,zb)
      qgWq=aveqg*fac*virt5(iqgq,za,zb)
      gqbWqb=aveqg*fac*virt5(igqbqb,za,zb)
      qbgWqb=aveqg*fac*virt5(iqbgqb,za,zb)

      do j=-nf,nf
      do k=-nf,nf

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg-subuv*msq0(j,k)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg-subuv*msq0(j,k)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
     &     -subuv*msq0(j,k)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
     &     -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
     &     -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
     &     -subuv*msq0(j,k)
      endif

      enddo
      enddo

      return
      end
