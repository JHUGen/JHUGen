      subroutine gen5(r,p,wt5,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'kprocess.f'
      include 'phasemin.f'
      include 'interference.f'
      include 'x1x2.f'
      include 'nproc.f'
      integer:: nu
      real(dp):: r(mxdim),wt5,
     & p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: pswt,xjac
      real(dp):: tau,y
      include 'energy.f'
      integer, save:: icount=1
!$omp threadprivate(icount)

      wt5=0._dp
      tau=exp(log(taumin)*r(9))
      y=0.5_dp*log(tau)*(1._dp-2._dp*r(10))
      xjac=log(taumin)*tau*log(tau)

      xx(1)=sqrt(tau)*exp(+y)
      xx(2)=sqrt(tau)*exp(-y)

c--- phase space volume only checked for x1=x2=1
      if (kcase==kvlchwt) then
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

      if     ((kcase==kt_bbar) .or. (kcase==kqg_tbb)) then
        call phase51(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif ((kcase==kW_twdk) .or. (kcase==kWtdkay)
     &    .or.(kcase==kW_cwdk) .or. (kcase==kvlchwt))  then
        call phase5a(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif (kcase==kvlchk5)  then
        call phase5(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif ((kcase==kWWqqdk) .or. (kcase==kWHbbdk)
     &       .or.(kcase==kZHbbdk))  then
        call phase5h(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      else
         call phase5(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      enddo

      if (interference) then
        if (icount == 1) then
          bw34_56=.true.
          icount=icount-1
        else
          bw34_56=.false.
          do nu=1,4
            p(4,nu)=p6(nu)
            p(6,nu)=p4(nu)
          enddo
          icount=icount+1
        endif
      endif

      wt5=xjac*pswt

      if (wt5 == 0._dp) return 1

      return
      end
