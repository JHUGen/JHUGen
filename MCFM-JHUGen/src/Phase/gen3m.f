      subroutine gen3m(r,p,m3,m4,m5,wt3,*)
c----generate 3 dimensional phase space weight and vectors p(mxpart,4)
c----           p1+p2+p3+p4+p5=0
c----and x1 and x2 given seven random numbers
c----p(6,i) and p(7,i) are set equal to zero
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'masses.f'
      include 'process.f'
      include 'phasemin.f'
      include 'x1x2.f'
      integer nu,j
      double precision r(mxdim),wt3
      double precision p(mxpart,4),taulowest,
     . p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),m3,m4,m5
      double precision pswt,xjac,tau,y
      include 'energy.f'

      do nu=1,4
      do j=1,mxpart
      p(j,nu)=0d0
      enddo 
      enddo 
      wt3=0d0
      
      taulowest=max(taumin,((m3+m4+m5)/sqrts)**2)
      
      tau=exp(log(taulowest)*r(6))
      y=0.5d0*log(tau)*(1d0-2d0*r(7))
      xjac=log(taulowest)*tau*log(tau)
      
      xx(1)=dsqrt(tau)*exp(+y)
      xx(2)=dsqrt(tau)*exp(-y)

c--- for comparison with C. Oleari's e+e- --> QQbg calculation
c      if (runstring(1:5) .eq. 'carlo') then
c        xx(1)=1d0
c      xx(2)=1d0
c      xjac=1d0
c      endif

c---if x's out of normal range alternative return
      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      call phase3m(r,p1,p2,p3,p4,p5,p6,p7,m3,m4,m5,pswt)

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      enddo 

      wt3=xjac*pswt

      if(wt3 .eq. 0d0) return 1

      return
      end
