      subroutine gen4mdkrad(r,p,pswt,*)
      implicit none
      include 'types.f'
c--- this routine is an extension of gen4 to include the decay
c--- of one of the heavy particles, with radiation included in the decay
      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'limits.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'x1x2.f'
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),pswt,smin
      real(dp):: p1(4),p2(4),p12(4),p8(4),p34568(4),
     & p7(4),p3458(4),p345(4),p6(4),wt12,wt34568,wt3458
      real(dp):: p34(4),p5(4),p3(4),p4(4),wt345,wt34
      real(dp):: mtbsq
      integer:: nu
      real(dp):: xjac,p1ext(4),p2ext(4),wt0
      real(dp):: tau,x1mx2,surd,lntaum
      common/pext/p1ext,p2ext
      parameter(wt0=1._dp/twopi**2)
!$omp threadprivate(/pext/)
 
c--- this part is taken from gen4
      lntaum=log(taumin)
      tau=exp(lntaum*(one-r(14)))
      xjac=-lntaum*tau

      x1mx2=two*r(15)-one
      surd=sqrt(x1mx2**2+four*tau) 
           
      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)

      pswt=xjac*two/surd

      if   ((xx(1) > 1._dp)  .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin) .or. (xx(2) < xmin)) return 1 

      do nu=1,4
      p1(nu)=xx(1)*p1ext(nu)
      p2(nu)=xx(2)*p2ext(nu)
      p12(nu)=-p1(nu)-p2(nu)
      enddo

c--- these must be set this way for this part
      mass3=mb ! check this
      n3=0

c--- taken from phase4
      mtbsq=(mt+mb)**2
      call phi1_2m(0._dp,r(1),r(2),r(3),mtbsq,p12,p7,p34568,wt12,*99)
      call phi3m(r(4),r(5),p34568,p3458,p6,mt,mb,wt34568,*99)
c      r(6)=1._dp-1.e-4_dp*r(6) ! debug - soft p8
      call phi1_2m(0._dp,r(6),r(7),r(8),mb,p3458,p8,p345,wt3458,*99)
      pswt=pswt*wt0*wt12*wt34568*wt3458
c--- alternative
c      r(2)=r(2)*1.e-5_dp ! DEBUG: s78 small
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
c--- NOTE: for now, generate massless b in top quark decay
c      call phi1_2m(zip,r(9),r(10),r(11),smin,p345,p5,p34,wt345,*99)
c--- decay W -> e n      
      call phi3m0(r(12),r(13),p34,p3,p4,wt34,*99)
      
c--- compute new weight
      pswt=pswt/twopi**2*wt345*wt34*pi*mt*twidth
      
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
      
