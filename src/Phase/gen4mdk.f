      subroutine gen4mdk(r,p,pswt,*)
c--- this routine is an extension of gen4 to include the decay
c--- of one of the heavy particles

      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'limits.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'process.f'
      include 'x1x2.f'
      double precision r(mxdim)
      double precision p(mxpart,4),pswt,smin
      double precision p1(4),p2(4),p12(4),p8(4),p34567(4),
     & p7(4),p3456(4),p345(4),p6(4),wt12,wt34567,wt3456
      double precision p567(4),p56(4),p34(4),p5(4),p3(4),p4(4)
      double precision mtbsq,wt345,wt34,wt567,wt56
      integer nu
      double precision xjac,p1ext(4),p2ext(4),wt0
      double precision tau,x1mx2,surd,lntaum,tmin
      common/pext/p1ext,p2ext
      include 'energy.f'
      parameter(wt0=1d0/twopi**2)
!$omp threadprivate(/pext/)

      if ((case .eq. 'Z_tdkj') .or. (case .eq. 'H_tdkj')) then
      call gen4(r,p,pswt,*99)

      p1(:)=p(1,:)
      p2(:)=p(2,:)
      p3(:)=p(3,:)
      p4(:)=p(4,:)
      p567(:)=p(5,:)
      p8(:)=p(6,:)

c--- set up minimum invariant mass for the W in the top decay
      if (zerowidth) then
        smin=wmass**2
      else
        smin=bbsqmin
      endif

c--- decay top -> b W
      call phi1_2m_bw(mb,r(11),r(12),r(13),smin,p567,p7,p56,
     & wmass,wwidth,wt567,*99)
c--- decay W -> e n
      call phi3m0(r(14),r(15),p56,p5,p6,wt56,*99)
      pswt=pswt/twopi**2*wt567*wt56*pi*mt*twidth

      else
      mtbsq=(mt+mb)**2
      tmin=mtbsq/sqrts**2

c--- this part is taken from gen4
      lntaum=dlog(tmin)
      tau=dexp(lntaum*(one-r(14)))
      xjac=-lntaum*tau

      x1mx2=two*r(15)-one
      surd=dsqrt(x1mx2**2+four*tau)

      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)

      pswt=xjac*two/surd

      if   ((xx(1) .gt. 1d0)  .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin) .or. (xx(2) .lt. xmin)) return 1

      do nu=1,4
      p1(nu)=xx(1)*p1ext(nu)
      p2(nu)=xx(2)*p2ext(nu)
      p12(nu)=-p1(nu)-p2(nu)
      enddo

c--- these must be set this way for this part
      mass3=mb ! check this
      n3=0

c--- taken from phase4
cc      r(1)=1d0-r(1)/1d2 ! soft p8
      call phi1_2m(0d0,r(1),r(2),r(3),mtbsq,p12,p8,p34567,wt12,*99)
cc      p8(4)=10d0+1d-3/abs(p1(4))
cc      p8(3)=-10d0
cc      p8(1)=0d0
cc      p8(2)=dsqrt(p8(4)**2-p8(3)**2)
cc      do nu=1,4
cc      p34567(nu)=p12(nu)-p8(nu)
cc      enddo
      call phi1_2m(0d0,r(4),r(5),r(6),mtbsq,p34567,p7,p3456,wt34567,*99)
      call phi3m(r(7),r(8),p3456,p345,p6,mt,mb,wt3456,*99)
      pswt=pswt*wt0*wt12*wt34567*wt3456
c--- alternative
cc      r(2)=r(2)*1d-5 ! DEBUG: s78 small
c      call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
c      call phi3m0(r(5),r(6),p78,p7,p8,wt78,*99)
c      call phi3m(r(7),r(8),p3456,p345,p6,mt,mb,wt3456,*99)
c      pswt=pswt*wt0*wt12*wt78*wt3456

c--- these must be set this way for this part
      mass3=wmass
      n3=1

c--- set up minimum invariant mass for the W in the top decay
      if (zerowidth) then
        smin=wmass**2
      else
        smin=wsqmin
      endif

c--- decay top -> b W
      call phi1_2m(mb,r(9),r(10),r(11),smin,p345,p5,p34,wt345,*99)
c--- decay W -> e n
      call phi3m0(r(12),r(13),p34,p3,p4,wt34,*99)

c--- compute new weight
      pswt=pswt/twopi**2*wt345*wt34*pi*mt*twidth

      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      p(8,nu)=p8(nu)
      enddo

      return

   99 continue
      return 1

      end

