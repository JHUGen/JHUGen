      subroutine qqb_ttw_v(p,msqv)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R. K. Ellis                                              *
*     March, 2012.                                                     *
*                                                                      *
*     Calculate the virtual matrix element squared for the process     *
*                                                                      *
*     q(-p1) +qbar(-p2)=                                               *
*                       +t(nu(p3)+e+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu~(p8))                    *
*                       +nu(p9) + e^+(p10)                             *
*                                                                      *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'momwbbm.f'
      include 'scheme.f'
      include 'plabel.f'
      integer:: j,k,nu,j1,j2,hb,hc
      real(dp):: p(mxpart,4),q(mxpart,4),msqv(-nf:nf,-nf:nf)
      real(dp):: qqb,qbq,dot
      real(dp):: fac,mQsq,s56,betasq
      complex(dp):: a6treemm,a6treemp,a6treepm,a6treepp,
     & mtop(2,2),manti(2,2),a61mm,a61mp,a61pm,a61pp,a61(2,2),a6(2,2),
     & loqbq(2,2),hoqbq(2,2),loqqb(2,2),hoqqb(2,2)
      logical:: numcheck
      common/numcheck/numcheck
!$omp threadprivate(/numcheck/)

      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

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

      do nu=1,4
      q(1,nu)=p(1,nu)
      q(2,nu)=p(2,nu)
      q(3,nu)=p(9,nu)
      q(4,nu)=p(10,nu)
      q(5,nu)=p(3,nu)+p(4,nu)+p(5,nu)
      q(6,nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo
      mQsq=mt**2


c--- construct the massless momenta a la Rodrigo
      do j=1,4
      do nu=1,4
      mom(j,nu)=q(j,nu)
      enddo
      enddo
      s56=2._dp*dot(q,5,6)+2._dp*mQsq
      betasq=1._dp-4._dp*mQsq/s56
      if (betasq >= 0._dp) then
        bp=0.5_dp*(1._dp+sqrt(betasq))
        bm=1._dp-bp
      else
        write(6,*) 'betasq < 0 in qqb_ttw_v.f, betasq=',betasq
        call flush(6)
        stop
      endif
      do nu=1,4
      mom(5,nu)=(bp*q(5,nu)-bm*q(6,nu))/sqrt(betasq)
      mom(6,nu)=(bp*q(6,nu)-bm*q(5,nu))/sqrt(betasq)
      enddo

      call tdecayrod(p,3,4,5,6,7,8,0,mtop)
      call adecayrod(p,3,4,5,6,7,8,0,manti)
c--- compute spinor products
      call spinoru(6,mom,za,zb)

c--- overall factor
      fac=V*gsq**2*gwsq**6*aveqq/(mt*twidth)**4
      fac=fac*xn*ason2pi
      fac=fac*s(3,4)**2/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

c--- include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

c--- QBQ: compute 1-loop and tree amplitudes
      call a61mass(1,6,5,2,4,3,mQsq,a61mm,a61mp,a61pm,a61pp,
     & a6treemm,a6treemp,a6treepm,a6treepp)
      a61(1,1)=a61mm
      a61(1,2)=a61mp
      a61(2,1)=a61pm
      a61(2,2)=a61pp
      a6(1,1)=a6treemm
      a6(1,2)=a6treemp
      a6(2,1)=a6treepm
      a6(2,2)=a6treepp

      qbq=0._dp
      do hb=1,2
      do hc=1,2
      hoqbq(hb,hc)=czip
      loqbq(hb,hc)=czip
      do j1=1,2
      do j2=1,2
      loqbq(hb,hc)=loqbq(hb,hc)+mtop(hb,j1)*a6(j1,j2)*manti(j2,hc)
      hoqbq(hb,hc)=hoqbq(hb,hc)+mtop(hb,j1)*a61(j1,j2)*manti(j2,hc)
      enddo
      enddo
      qbq=qbq+fac*real(loqbq(hb,hc)*conjg(hoqbq(hb,hc)))
      enddo
      enddo

c--- put a pause here when writing out primitives
c      if (numcheck) pause

c--- QQB: compute 1-loop and tree amplitudes
      call a61mass(2,6,5,1,4,3,mQsq,a61mm,a61mp,a61pm,a61pp,
     & a6treemm,a6treemp,a6treepm,a6treepp)
      a61(1,1)=a61mm
      a61(1,2)=a61mp
      a61(2,1)=a61pm
      a61(2,2)=a61pp
      a6(1,1)=a6treemm
      a6(1,2)=a6treemp
      a6(2,1)=a6treepm
      a6(2,2)=a6treepp

      qqb=0._dp
      do hb=1,2
      do hc=1,2
      hoqqb(hb,hc)=czip
      loqqb(hb,hc)=czip
      do j1=1,2
      do j2=1,2
      loqqb(hb,hc)=loqqb(hb,hc)+
     & mtop(hb,j1)*a6(j1,j2)*manti(j2,hc)
      hoqqb(hb,hc)=hoqqb(hb,hc)+
     & mtop(hb,j1)*a61(j1,j2)*manti(j2,hc)
      enddo
      enddo
      qqb=qqb+fac*real(loqqb(hb,hc)*conjg(hoqqb(hb,hc)))
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
               msqv(j,k)=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
               msqv(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end

