      subroutine gQ_zQ_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*    Authors: R.K. Ellis and John Campbell                             *
*    July, 2003.                                                       *
*    Matrix element for Z + heavy quark (of flavour "flav") production *
*    in order alpha_s^2                                                *
*    averaged over initial colours and spins                           *
*     g(-p1)+Q(-p2)-->Z^+(l(p3)+a(p4))+Q(p5)                           *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
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
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,sz,virt5,subuv
      real(dp):: gqZqLL,gqZqRR,gqZqLR,gqZqRL
      real(dp):: qgZqLL,qgZqRR,qgZqLR,qgZqRL
      real(dp):: gqbZqbLL,gqbZqbRR,gqbZqbLR,gqbZqbRL
      real(dp):: qbgZqbLL,qbgZqbRR,qbgZqbLR,qbgZqbRL
      complex(dp):: prop
      integer,parameter::
     & iqgqLL(5)=(/1,5,3,4,2/),iqgqRR(5)=(/5,1,4,3,2/),
     & iqgqRL(5)=(/5,1,3,4,2/),iqgqLR(5)=(/1,5,4,3,2/),
     & igqqLL(5)=(/2,5,3,4,1/),igqqRR(5)=(/5,2,4,3,1/),
     & igqqRL(5)=(/5,2,3,4,1/),igqqLR(5)=(/2,5,4,3,1/)

      scheme='dred'
c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(5,p,za,zb)

c--- calculate lowest order
      call gQ_zQ(p,msq0)

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn*(epinv
     & *(11._dp-2._dp*real(nflav,kind=dp)/xn)-1._dp)/6._dp

c--   calculate propagator
      sz=s(3,4)
      prop=sz/cplx2((sz-zmass**2),zmass*zwidth)

      fac=8._dp*cf*xnsq*esq**2*gsq

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

      if     ((j > 0) .and. (k == 0)) then
          msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqLL
     &             +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqRR
     &             +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqLR
     &             +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqRL
     &             -subuv*msq0(j,k)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbLL
     &             +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbRR
     &             +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbLR
     &             +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbRL
     &             -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqLL
     &             +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqRR
     &             +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqLR
     &             +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqRL
     &             -subuv*msq0(j,k)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbLL
     &             +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbRR
     &             +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbLR
     &             +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbRL
     &             -subuv*msq0(j,k)
      endif

   19 continue
      enddo
      enddo

      return
      end
