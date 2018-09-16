!------- phase space for 2->4 with 2 massive particles (suitable for pp=>tt+2j etc.
!------- CW Oct 12

      subroutine gen4m(r,p,m3,m4,m5,m6,wt4,*)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'phasemin.f'
      include 'limits.f'
      include 'x1x2.f'
      integer:: nu,j
      real(dp):: r(mxdim),wt4,wt3456,wt34,wt56
      real(dp):: m3,m4,m5,m6
      real(dp):: p(mxpart,4),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p34(4),p56(4),p12(4)
      real(dp):: xjac,tau,y
      include 'energy.f'
      real(dp):: wt0
      parameter(wt0=one/twopi**2)

      do nu=1,4
         do j=1,mxpart
            p(j,nu)=0._dp
         enddo
      enddo

      taumin=max((m3+m4+m5+m6)**2/sqrts**2,1.e-4_dp)

      wt4=0._dp
      tau=exp(log(taumin)*r(9))
      y=0.5_dp*log(tau)*(one-two*r(10))

      xjac=log(taumin)*tau*log(tau)
      xx(1)=sqrt(tau)*exp(y)
      xx(2)=sqrt(tau)*exp(-y)


      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      do j=1,4
         p12(j)=-p1(j)-p2(j)
      enddo
      call phi1_2(r(1),r(2),r(3),r(4),p12,p34,p56,wt3456,*999)
!----- decay 56
      call phi3m0(r(5),r(6),p56,p5,p6,wt56,*999)
!------- decay 34
      Call phi3m(r(7),r(8),p34,p3,p4,m3,m4,wt34,*999)



      wt4=wt0*xjac*wt3456*wt34*wt56
      do nu=1,4
         p(1,nu)=p1(nu)
         p(2,nu)=p2(nu)
         p(3,nu)=p3(nu)
         p(4,nu)=p4(nu)
         p(5,nu)=p5(nu)
         p(6,nu)=p6(nu)
      enddo

      return
 999  continue
      wt4=0._dp
      return 1
      end
