      subroutine qqb_z1jet_v(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John Campbell                            *
*     May, 2001.                                                       *
*     Matrix element for Z + jet production                            *
*     in order alpha_s^2                                               *
*     averaged over initial colours and spins                          *
*     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+g(p5)                        *
************************************************************************
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'cplx.h'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,sz,virt5,subuv,virt5ax,sin2winv
      real(dp):: qqbZgLL(2),qqbZgRR(2),qqbZgLR(2),qqbZgRL(2)
      real(dp):: gqZqLL(2),gqZqRR(2),gqZqLR(2),gqZqRL(2)
      real(dp):: qgZqLL(2),qgZqRR(2),qgZqLR(2),qgZqRL(2)
      real(dp):: qbqZgLL(2),qbqZgRR(2),qbqZgLR(2),qbqZgRL(2)
      real(dp):: gqbZqbLL(2),gqbZqbRR(2),gqbZqbLR(2),gqbZqbRL(2)
      real(dp):: qbgZqbLL(2),qbgZqbRR(2),qbgZqbLR(2),qbgZqbRL(2)
      complex(dp):: prop
      logical, parameter:: includeanom=.true.
      integer,parameter::
     & iqqbgLL(5)=(/1,2,3,4,5/),iqqbgRR(5)=(/2,1,4,3,5/),
     & iqqbgRL(5)=(/2,1,3,4,5/),iqqbgLR(5)=(/1,2,4,3,5/),
     & iqgqLL(5)=(/1,5,3,4,2/),iqgqRR(5)=(/5,1,4,3,2/),
     & iqgqRL(5)=(/5,1,3,4,2/),iqgqLR(5)=(/1,5,4,3,2/),
     & igqqLL(5)=(/2,5,3,4,1/),igqqRR(5)=(/5,2,4,3,1/),
     & igqqRL(5)=(/5,2,3,4,1/),igqqLR(5)=(/2,5,4,3,1/)
      common/virt5ax/virt5ax
!$omp threadprivate(/virt5ax/)

      scheme='dred'

c--set msq=0 to initialize
      msq(:,:)=zero

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities
c      if ((abs(s(1,5))<cutoff) .or. (abs(s(2,5))<cutoff)) return

c--- calculate lowest order
      call qqb_z1jet(p,msq0)

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn
     & *(epinv*(11._dp-2._dp*real(nflav,dp)/xn)-1._dp)/6._dp

c--   calculate propagator
      sz=s(3,4)
      prop=sz/cplx2((sz-zmass**2),zmass*zwidth)

      fac=8._dp*cf*xnsq*esq**2*gsq

!     first letter L/R is the index of the fermion line
!     second letter L/R is the index of the Z/gamma
!     already summed over gluons
      qqbZgLL(1)=aveqq*fac*virt5(iqqbgLL,za,zb)
      qqbZgLL(2)=aveqq*fac*virt5ax
      qqbZgLR(1)=aveqq*fac*virt5(iqqbgLR,za,zb)
      qqbZgLR(2)=aveqq*fac*virt5ax
      qqbZgRL(1)=aveqq*fac*virt5(iqqbgRL,za,zb)
      qqbZgRL(2)=aveqq*fac*virt5ax
      qqbZgRR(1)=aveqq*fac*virt5(iqqbgRR,za,zb)
      qqbZgRR(2)=aveqq*fac*virt5ax


      qbqZgLL(:)=qqbZgRL(:)
      qbqZgLR(:)=qqbZgRR(:)
      qbqZgRL(:)=qqbZgLL(:)
      qbqZgRR(:)=qqbZgLR(:)

      gqZqLL(1)=aveqg*fac*virt5(igqqLL,za,zb)
      gqZqLL(2)=aveqg*fac*virt5ax
      gqZqLR(1)=aveqg*fac*virt5(igqqLR,za,zb)
      gqZqLR(2)=aveqg*fac*virt5ax
      gqZqRL(1)=aveqg*fac*virt5(igqqRL,za,zb)
      gqZqRL(2)=aveqg*fac*virt5ax
      gqZqRR(1)=aveqg*fac*virt5(igqqRR,za,zb)
      gqZqRR(2)=aveqg*fac*virt5ax

      gqbZqbRL(:)=gqZqLL(:)
      gqbZqbRR(:)=gqZqLR(:)
      gqbZqbLL(:)=gqZqRL(:)
      gqbZqbLR(:)=gqZqRR(:)


      qgZqLL(1)=aveqg*fac*virt5(iqgqLL,za,zb)
      qgZqLL(2)=aveqg*fac*virt5ax
      qgZqLR(1)=aveqg*fac*virt5(iqgqLR,za,zb)
      qgZqLR(2)=aveqg*fac*virt5ax
      qgZqRL(1)=aveqg*fac*virt5(iqgqRL,za,zb)
      qgZqRL(2)=aveqg*fac*virt5ax
      qgZqRR(1)=aveqg*fac*virt5(iqgqRR,za,zb)
      qgZqRR(2)=aveqg*fac*virt5ax

      qbgZqbRL(:)=qgZqLL(:)
      qbgZqbRR(:)=qgZqLR(:)
      qbgZqbLL(:)=qgZqRL(:)
      qbgZqbLR(:)=qgZqRR(:)

      sin2winv=L(2)-R(2)
      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
         msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
         msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgLL(1)
     &            +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgRR(1)
     &            +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgLR(1)
     &            +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgRL(1)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +abs(sin2winv*l1*prop)**2*(qqbZgLL(2)+qqbZgRL(2))
     &            +abs(sin2winv*r1*prop)**2*(qqbZgRR(2)+qqbZgLR(2))
         endif
      elseif ((j < 0) .and. (k > 0)) then
         msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgLL(1)
     &            +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgRR(1)
     &            +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgLR(1)
     &            +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgRL(1)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +abs(sin2winv*l1*prop)**2*(qbqZgLL(2)+qbqZgRL(2))
     &            +abs(sin2winv*r1*prop)**2*(qbqZgRR(2)+qbqZgLR(2))
         endif
      elseif ((j > 0) .and. (k == 0)) then
         msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqLL(1)
     &            +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqRR(1)
     &            +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqLR(1)
     &            +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqRL(1)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +abs(sin2winv*l1*prop)**2*(qgZqLL(2)+qgZqRL(2))
     &            +abs(sin2winv*r1*prop)**2*(qgZqRR(2)+qgZqLR(2))
         endif
      elseif ((j < 0) .and. (k == 0)) then
         msq(j,k)=+abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbLL(1)
     &            +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbRR(1)
     &            +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbLR(1)
     &            +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbRL(1)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +abs(sin2winv*l1*prop)**2*(qbgZqbLL(2)+qbgZqbRL(2))
     &            +abs(sin2winv*r1*prop)**2*(qbgZqbRR(2)+qbgZqbLR(2))
         endif
      elseif ((j == 0) .and. (k > 0)) then
         msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqLL(1)
     &            +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqRR(1)
     &            +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqLR(1)
     &            +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqRL(1)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +abs(sin2winv*l1*prop)**2*(gqZqLL(2)+gqZqRL(2))
     &            +abs(sin2winv*r1*prop)**2*(gqZqRR(2)+gqZqLR(2))
         endif
      elseif ((j == 0) .and. (k < 0)) then
         msq(j,k)=+abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbLL(1)
     &            +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbRR(1)
     &            +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbLR(1)
     &            +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbRL(1)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +abs(sin2winv*l1*prop)**2*(gqbZqbLL(2)+gqbZqbRL(2))
     &            +abs(sin2winv*r1*prop)**2*(gqbZqbRR(2)+gqbZqbLR(2))
         endif
      endif

   19 continue
      enddo
      enddo

      return
      end
