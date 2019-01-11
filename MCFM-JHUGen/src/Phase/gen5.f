      subroutine gen5(r,p,wt5,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'process.f'
      include 'phasemin.f'
      include 'interference.f'
      include 'x1x2.f'
      include 'nproc.f'
      integer nu,icount
      double precision r(mxdim),wt5,
     . p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision pswt,xjac
      double precision tau,y
      include 'energy.f'
      data icount/1/
      save icount
!$omp threadprivate(icount)

      wt5=0d0
      tau=dexp(dlog(taumin)*r(9))
      y=0.5d0*dlog(tau)*(1d0-2d0*r(10))
      xjac=dlog(taumin)*tau*dlog(tau)

      xx(1)=dsqrt(tau)*dexp(+y)
      xx(2)=dsqrt(tau)*dexp(-y)

c--- phase space volume only checked for x1=x2=1
      if (case .eq. 'vlchwt') then
        xx(1)=1d0
        xx(2)=1d0
        xjac=1d0
      endif

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

      if     ((case .eq. 't_bbar') .or. (case .eq. 'qg_tbb')) then
        call phase51(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif ((case .eq. 'W_twdk') .or. (case .eq. 'Wtdkay')
     .    .or.(case .eq. 'W_cwdk') .or. (case .eq. 'vlchwt'))  then
        call phase5a(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif (case .eq. 'vlchk5')  then
        call phase5(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif (case .eq. 'WWqqdk')  then
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
        if (icount .eq. 1) then
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

      if (wt5 .eq. 0d0) return 1

      return
      end
