      subroutine gen3h(r,p,wt3,*)
      implicit none
      include 'types.f'
c----generate 3 dimensional phase space weight and vectors p(7,4)
c----and x1 and x2 given seven random numbers
c----p(5,i) and p(4,i) are set equal to zero
c---- s12 is generated according to a Breit-Wigner at mH
c----
c---- if 'nodecay' is true, then the vector boson decay into massless
c---- particles is not included and 2 less integration variables
c---- are required

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'nodecay.f'
      include 'masses.f'
      include 'x1x2.f'
      include 'energy.f'
      integer:: nu

      real(dp):: r(mxdim),s12,wt3,wtbw,rdk1,rdk2,
     & p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: pswt,xjac,tau,y

      wt3=0._dp

c--- dummy values if there's no decay
      if (nodecay) then
        rdk1=0.5_dp
        rdk2=0.5_dp
      else
        rdk1=r(6)
        rdk2=r(7)
      endif

      call breitw(r(4),0._dp,sqrts**2,hmass,hwidth,s12,wtbw)

      tau=s12/sqrts**2
      y=0.5_dp*log(tau)*(1._dp-2._dp*r(5))
      xjac=-wtbw/sqrts**2*log(tau)

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


      call phase3(r(1),r(2),r(3),rdk1,rdk2,p1,p2,p3,p4,p5,p6,p7,pswt)

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
