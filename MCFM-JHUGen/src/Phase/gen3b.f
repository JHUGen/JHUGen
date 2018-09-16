      subroutine gen3b(r,p,wt3,*)
      implicit none
      include 'types.f'
c----generate 3 dimensional phase space weight and vectors p(7,4)
c----           p1+p2+p3+p4+p5=0
c----and x1 and x2 given seven random numbers
c----p(6,i) and p(7,i) are set equal to zero

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'x1x2.f'
      include 'energy.f'
      integer:: nu
      real(dp):: r(mxdim),wt3
      real(dp):: p(mxpart,4),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: pswt,xjac,tau,y

      wt3=0._dp
      p(:,:)=0._dp
      tau=exp(log(taumin)*r(4))
      y=0.5_dp*log(tau)*(1._dp-2._dp*r(5))
      xjac=log(taumin)*tau*log(tau)

      xx(1)=sqrt(tau)*exp(+y)
      xx(2)=sqrt(tau)*exp(-y)

c---if x's out of normal range alternative return
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


      call phase3(r(1),r(2),r(3),r(6),r(7),p1,p2,p3,p4,p5,p6,p7,pswt)

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      enddo

      wt3=xjac*pswt
      if(wt3 == 0._dp) then
      p(:,:)=0._dp
      return 1
      endif

      return
      end
