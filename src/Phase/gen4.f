      subroutine gen4(r,p,wt4,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'process.f'
      include 'phasemin.f'
      include 'interference.f'
      include 'x1x2.f'
      include 'nproc.f'
      integer nu,icount
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(mxpart,4)
      double precision pswt,xjac
      double precision tau,x1mx2,surd
      double precision lntaum
      include 'energy.f'
!      data icount/1/
!      save icount
!!$omp threadprivate(icount)
      wt4=0d0


      lntaum=dlog(taumin)
      tau=dexp(lntaum*(one-r(9)))
      xjac=-lntaum*tau

c      tau=(one-taumin)*r(9)**2+taumin
c      xjac=2*r(9)*(one-taumin)

      x1mx2=two*r(10)-one
      surd=dsqrt(x1mx2**2+four*tau)

      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)
!      write(*,*) r(1),r(9),r(10)
!      write(*,*) xx(1),xx(2)

      xjac=xjac*two/surd

c--- for comparison with C. Oleari's e+e- --> QQbg calculation
c      if (runstring(1:5) .eq. 'carlo') then
c        xx(1)=1d0
c        xx(2)=1d0
c        xjac=1d0
c      endif

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

      if  ((case .eq. 't_bbar')
     . .or.(case .eq. 'bq_tpq')) then
      call phase41(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      elseif (
     .      (case .eq. 'qqttbb')
     . .or. (case .eq. 'qqttgg'))  then
      call phase4m(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      elseif (case .eq. 'vlchk4')  then
      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      else
      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=0d0
      enddo

      if (interference) then
!        if (icount .eq. 1) then
          bw34_56=.true.
!          icount=icount-1
!        else
!          bw34_56=.false.
!          do nu=1,4
!            p(4,nu)=p6(nu)
!            p(6,nu)=p4(nu)
!          enddo
!          icount=icount+1
!        endif
      endif

      wt4=xjac*pswt
      if (debug) write(6,*) 'wt4 in gen4',wt4
      return

 999  return 1
      end

