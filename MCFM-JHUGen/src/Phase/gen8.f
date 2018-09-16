      subroutine gen8(r,q,wt8,*)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'debug.f'
      include 'x1x2.f'
      integer:: nu
      real(dp):: r(mxdim),wt8,q(mxpart,4),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     & p9(4),p10(4),pswt,xjac,tau,x1mx2,surd,lntaum
      real(dp):: p1ext(4),p2ext(4)
      common/pext/p1ext,p2ext
!$omp threadprivate(/pext/)
      wt8=0._dp
      q(:,:)=0._dp
      lntaum=log(taumin)
      tau=exp(lntaum*(one-r(21)))
      xjac=-lntaum*tau

c      tau=(one-taumin)*r(14)**2+taumin
c      xjac=2*r(13)*(one-taumin)

      x1mx2=two*r(22)-one
      surd=sqrt(x1mx2**2+four*tau)

      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)
      xjac=xjac*two/surd

      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) return 1

      do nu=1,4
      p1(nu)=xx(1)*p1ext(nu)
      p2(nu)=xx(2)*p2ext(nu)
      enddo


      call phase8(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,pswt,*999)

      do nu=1,4
      q(1,nu)=p1(nu)
      q(2,nu)=p2(nu)
      q(3,nu)=p3(nu)
      q(4,nu)=p4(nu)
      q(5,nu)=p5(nu)
      q(6,nu)=p6(nu)
      q(7,nu)=p7(nu)
      q(8,nu)=p8(nu)
      q(9,nu)=p9(nu)
      q(10,nu)=p10(nu)

      enddo
      wt8=xjac*pswt
      if (debug) write(6,*) 'wt8 in gen8',wt8
      return

      wt8=0._dp
 999  return 1
      end

