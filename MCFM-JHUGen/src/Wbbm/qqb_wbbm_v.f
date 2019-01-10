      subroutine qqb_wbbm_v(p,msqv)
************************************************************************
*     Author: J. M. Campbell                                           *
*     October, 2010.                                                   *
*                                                                      *
*     Calculate the virtual matrix element squared for the process     *
*                                                                      *
*     q(-p1) + qb(-p2) --> nu(p3) + e^+(p4) + b(p5) + bbar(p6)         *
*                                                                      *
*     in which the mass of the b-quark is kept non-zero.               *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'ckm.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'heavyflav.f'
      include 'momwbbm.f'
      include 'scheme.f'
      integer j,k,nu
      double precision p(mxpart,4),q(mxpart,4),msqv(-nf:nf,-nf:nf)
      double precision qqb,qbq
      double precision fac,mQsq,s56,betasq
      double complex a6treemm,a6treemp,a6treepm,a6treepp
      double complex a61mm,a61mp,a61pm,a61pp
      logical numcheck
      common/numcheck/numcheck
!$omp threadprivate(/numcheck/)

      scheme='dred'

C----Intialize whole array to zero
      msqv(:,:)=0d0

c--- ensure that QCDLoop is initialized (now done in computescalars.f)
c      if (first) then
c        call qlinit
c        first=.false.
c      endif

c--- set the following flag to true to write out values of different primitives
c--- (a similar flag, to write out coefficients of the basis integrals,
c---  can be found in the routine a61mass)
      numcheck=.false.
c--- setup for performing check against numerical evaluation
      if (numcheck) then
c----- read in special point
        include 'MCFMpoint.f'
c        call writeout(p)
c        pause
c--- perform the usual business to rotate away from the z-direction
        do j=1,6
        q(j,4)=p(j,4)
        q(j,1)=p(j,3)
        q(j,2)=-p(j,2)
        q(j,3)=p(j,1)
        do k=1,4
        p(j,k)=q(j,k)
        enddo
        enddo
      endif

C--- set up the correct mass, according to 'flav'
      if     (flav .eq. 6) then
        mQsq=mt**2
      elseif (flav .eq. 5) then
        mQsq=mb**2
      elseif (flav .eq. 4) then
        mQsq=mc**2
      else
        write(6,*) 'Wrong flavour in qqb_wbbm_v.f: flav=',flav
        call flush(6)
        stop
      endif

C--- fill dot-products
      call dotem(6,p,s)

c--- construct the massless momenta a la Rodrigo
      do j=1,4
      do nu=1,4
      mom(j,nu)=p(j,nu)
      enddo
      enddo
      s56=s(5,6)+2d0*mQsq
      betasq=1d0-4d0*mQsq/s56
      if (betasq .ge. 0d0) then
        bp=0.5d0*(1d0+dsqrt(betasq))
        bm=1d0-bp
      else
        write(6,*) 'betasq < 0 in qqb_wbbm_v.f, betasq=',betasq
        call flush(6)
        stop
      endif
      do nu=1,4
      mom(5,nu)=(bp*p(5,nu)-bm*p(6,nu))/dsqrt(betasq)
      mom(6,nu)=(bp*p(6,nu)-bm*p(5,nu))/dsqrt(betasq)
      enddo

c--- compute spinor products
      call spinoru(6,mom,za,zb)

c--- overall factor
      fac=V*gsq**2*gwsq**2*aveqq
      fac=fac*xn*ason2pi
c--- don't forget to include the W propagator
      fac=fac*s(3,4)**2/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

c--- QBQ: compute 1-loop and tree amplitudes
      call a61mass(1,6,5,2,4,3,mQsq,a61mm,a61mp,a61pm,a61pp,
     & a6treemm,a6treemp,a6treepm,a6treepp)

      qbq=fac*dble(a6treemm*dconjg(a61mm)+a6treemp*dconjg(a61mp)
     &            +a6treepm*dconjg(a61pm)+a6treepp*dconjg(a61pp))

c--- put a pause here when writing out primitives
c      if (numcheck) pause

c--- QQB: compute 1-loop and tree amplitudes
      call a61mass(2,6,5,1,4,3,mQsq,a61mm,a61mp,a61pm,a61pp,
     & a6treemm,a6treemp,a6treepm,a6treepp)

      qqb=fac*dble(a6treemm*dconjg(a61mm)+a6treemp*dconjg(a61mp)
     &            +a6treepm*dconjg(a61pm)+a6treepp*dconjg(a61pp))

      do j=-(flav-1),(flav-1)
      do k=-(flav-1),(flav-1)
      if     ((j .gt. 0) .and. (k .lt. 0)) then
               msqv(j,k)=Vsq(j,k)*qqb
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
               msqv(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end

