      subroutine gQ_zQ_v(p,msq)
      implicit none
************************************************************************
*    Authors: R.K. Ellis and John Campbell                             *
*    July, 2003.                                                       *
*    Matrix element for Z + heavy quark (of flavour "flav") production *
*    in order alpha_s^2                                                *
*    averaged over initial colours and spins                           *
*     g(-p1)+Q(-p2)-->Z^+(l(p3)+a(p4))+Q(p5)                           *
************************************************************************
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'heavyflav.f'
      include 'nflav.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     . p(mxpart,4),fac,sz,virt5,subuv
      double precision gqZqLL,gqZqRR,gqZqLR,gqZqRL
      double precision qgZqLL,qgZqRR,qgZqLR,qgZqRL
      double precision gqbZqbLL,gqbZqbRR,gqbZqbLR,gqbZqbRL
      double precision qbgZqbLL,qbgZqbRR,qbgZqbLR,qbgZqbRL
      double complex prop
      integer,parameter::
     & iqgqLL(5)=(/1,5,3,4,2/),iqgqRR(5)=(/5,1,4,3,2/),
     & iqgqRL(5)=(/5,1,3,4,2/),iqgqLR(5)=(/1,5,4,3,2/),
     & igqqLL(5)=(/2,5,3,4,1/),igqqRR(5)=(/5,2,4,3,1/),
     & igqqRL(5)=(/5,2,3,4,1/),igqqLR(5)=(/2,5,4,3,1/)

      scheme='dred'
c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(5,p,za,zb)

c--- calculate lowest order
      call gQ_zQ(p,msq0)

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn*(epinv*(11d0-2d0*dble(nflav)/xn)-1d0)/6d0

c--   calculate propagator
      sz=s(3,4)
      prop=sz/dcmplx((sz-zmass**2),zmass*zwidth)

      fac=8d0*cf*xnsq*esq**2*gsq

      gqZqLL=aveqg*fac*virt5(igqqLL,za,zb)
      gqZqLR=aveqg*fac*virt5(igqqLR,za,zb)
      gqZqRL=aveqg*fac*virt5(igqqRL,za,zb)
      gqZqRR=aveqg*fac*virt5(igqqRR,za,zb)

      gqbZqbRL=gqZqLL
      gqbZqbRR=gqZqLR
      gqbZqbLL=gqZqRL
      gqbZqbLR=gqZqRR

      qgZqLL=aveqg*fac*virt5(iqgqLL,za,zb)
      qgZqLR=aveqg*fac*virt5(iqgqLR,za,zb)
      qgZqRL=aveqg*fac*virt5(iqgqRL,za,zb)
      qgZqRR=aveqg*fac*virt5(iqgqRR,za,zb)

      qbgZqbRL=qgZqLL
      qbgZqbRR=qgZqLR
      qbgZqbLL=qgZqRL
      qbgZqbLR=qgZqRR



      do j=-flav,flav,flav
      do k=-flav,flav,flav

      if( abs(j+k) .ne. flav) goto 19

      if     ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+cdabs(Q(j)*q1+L(j)*l1*prop)**2*qgZqLL
     .             +cdabs(Q(j)*q1+R(j)*r1*prop)**2*qgZqRR
     .             +cdabs(Q(j)*q1+L(j)*r1*prop)**2*qgZqLR
     .             +cdabs(Q(j)*q1+R(j)*l1*prop)**2*qgZqRL
     .             -subuv*msq0(j,k)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+cdabs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbLL
     .             +cdabs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbRR
     .             +cdabs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbLR
     .             +cdabs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbRL
     .             -subuv*msq0(j,k)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=+cdabs(Q(k)*q1+L(k)*l1*prop)**2*gqZqLL
     .             +cdabs(Q(k)*q1+R(k)*r1*prop)**2*gqZqRR
     .             +cdabs(Q(k)*q1+L(k)*r1*prop)**2*gqZqLR
     .             +cdabs(Q(k)*q1+R(k)*l1*prop)**2*gqZqRL
     .             -subuv*msq0(j,k)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=+cdabs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbLL
     .             +cdabs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbRR
     .             +cdabs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbLR
     .             +cdabs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbRL
     .             -subuv*msq0(j,k)
      endif

   19 continue
      enddo
      enddo

      return
      end
