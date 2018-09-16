      subroutine gen6(r,q,wt6,*)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'debug.f'
      include 'phasemin.f'
      include 'kprocess.f'
      include 'masses.f'
      include 'jetcuts.f'
      include 'x1x2.f'
      include 'energy.f'
      integer:: nu
      real(dp):: r(mxdim)
      real(dp):: wt6,q(mxpart,4)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      real(dp):: y,pswt,xjac,tau,vs,vsqmax,vsqmin,
     & s34,rtshat,ymax,yave

      wt6=0._dp

      if (((kcase==kHWW2jt) .or. (kcase==kHZZ2jt))
     &  .and. (hmass < 201._dp)) then
c--- Higgs production with WW or ZZ decay, small Higgs width
c--- (small error induced if zerowidth false, of order (hwidth/hmass))
        vsqmax=1._dp/(hmass*(hmass+4._dp*ptjetmin))
        vsqmin=1._dp/sqrts**2
        xmin=vsqmin/vsqmax
        vs=(vsqmax-vsqmin)*r(9)+vsqmin
        s34=1._dp/vs
        rtshat=sqrt(s34)
        ymax=log(sqrts/rtshat)
        yave=ymax*(two*r(10)-1._dp)
        xx(1)=rtshat/sqrts*exp(+yave)
        xx(2)=rtshat/sqrts*exp(-yave)
        xjac=(vsqmax-vsqmin)*s34**2/sqrts**2*two*ymax
      else
c--- generic process
        tau=exp(log(taumin)*r(9))
        y=0.5_dp*log(tau)*(1._dp-2._dp*r(10))
        xx(1)=sqrt(tau)*exp(+y)
        xx(2)=sqrt(tau)*exp(-y)
        xjac=log(taumin)*tau*log(tau)
      endif

c--- phase space volume only checked for x1=x2=1
      if ((kcase==kvlchwg) .or. (kcase==kvlchwh)) then
        xx(1)=1._dp
        xx(2)=1._dp
        xjac=1._dp
      endif

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

      if     ((kcase==kW_twdk) .or. (kcase==kWtbwdk)
     &   .or. (kcase==kW_cwdk) .or. (kcase==kvlchwh)) then
c--- W+t process, radiation in production
        call phase6a(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      elseif ((kcase==kWtdkay) .or. (kcase==kvlchwg)) then
c--- W+t process, radiation in decay
        call phase6b(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      else
c--- generic case
        call phase6(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      endif

      do nu=1,4
        q(1,nu)=p1(nu)
        q(2,nu)=p2(nu)
        q(3,nu)=p3(nu)
        q(4,nu)=p4(nu)
        q(5,nu)=p5(nu)
        q(6,nu)=p6(nu)
        q(7,nu)=p7(nu)
        q(8,nu)=p8(nu)
      enddo

      wt6=xjac*pswt

!      if (debug) write(6,*) 'wt6 in gen6',wt6

      return

 999  q(:,:)=0._dp
      return 1

      end

