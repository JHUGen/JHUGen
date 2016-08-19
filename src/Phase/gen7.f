      subroutine gen7(r,q,wt7,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'process.f'
      include 'phasemin.f'
      include 'debug.f'
      include 'masses.f'
      include 'jetcuts.f'
      include 'x1x2.f'
      integer nu
      double precision r(mxdim)
      double precision wt7,q(mxpart,4)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     . p9(4),pswt,xjac
      double precision tau,y,vs,vsqmax,vsqmin,s34,rtshat,
     . ymax,yave
      include 'energy.f'
      wt7=0d0

      if     (((case .eq. 'HWW2jt') .or. (case .eq. 'HZZ2jt')
     &     .or.(case .eq. 'HWW3jt') .or. (case .eq. 'HZZ3jt'))
     &  .and. (hmass .lt. 201d0)) then
c--- Higgs production with WW or ZZ decay, small Higgs width
c--- (small error induced if zerowidth false, of order (hwidth/hmass))
        vsqmax=1d0/(hmass*(hmass+4d0*ptjetmin))
        vsqmin=1d0/sqrts**2
        xmin=vsqmin/vsqmax
        vs=(vsqmax-vsqmin)*r(9)+vsqmin
        s34=1d0/vs
        rtshat=dsqrt(s34)
        ymax=dlog(sqrts/rtshat)
        yave=ymax*(two*r(10)-1d0)           
        xx(1)=rtshat/sqrts*exp(+yave)
        xx(2)=rtshat/sqrts*exp(-yave)
        xjac=(vsqmax-vsqmin)*s34**2/sqrts**2*two*ymax
      elseif ((case .eq. 'tt_ldk') .or. (case .eq. 'tt_hdk')
     &   .or. (case .eq. 'tt_udk') .or. (case .eq. 'tthWdk')) then
c--- top production with radiation in decay
        vsqmax=1d0/(4d0*mt**2)
        vsqmin=1d0/sqrts**2
        xmin=vsqmin/vsqmax
        vs=(vsqmax-vsqmin)*r(9)+vsqmin
        s34=1/vs
        rtshat=dsqrt(s34)
        ymax=dlog(sqrts/rtshat)
        yave=ymax*(two*r(10)-1d0)           
        xx(1)=rtshat/sqrts*exp(+yave)
        xx(2)=rtshat/sqrts*exp(-yave)
        xjac=(vsqmax-vsqmin)*s34**2/sqrts**2*two*ymax
      else
c--- generic process
        tau=dexp(dlog(taumin)*r(9))
        y=0.5d0*dlog(tau)*(1d0-2d0*r(10))
        xjac=dlog(taumin)*tau*dlog(tau)
        xx(1)=dsqrt(tau)*dexp(+y)
        xx(2)=dsqrt(tau)*dexp(-y)
      endif      

c--- phase space volume only checked for x1=x2=1
      if ((case .eq. 'vlchwg') .or. (case .eq. 'vlchwh')) then
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

      if  ((case .eq. 'qq_HWW') .or. (case .eq. 'qq_HZZ')
     ..or. (case .eq. 'HWW2jt') .or. (case .eq. 'HZZ2jt')
     ..or. (case .eq. 'HWW3jt') .or. (case .eq. 'HZZ3jt')
     ..or. (case .eq. 'WpWp3j')) then
        call  phase7a(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      elseif ((case .eq. 'WH__WW') .or. (case .eq. 'ZH__WW')) then
        call  phase7b(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      elseif ((case .eq. 'WH__ZZ') .or. (case .eq. 'ZH__ZZ')) then
        call  phase7b(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      elseif ((case .eq. 'tt_ldk') .or. (case .eq. 'tt_hdk')
     &   .or. (case .eq. 'tt_udk') .or. (case .eq. 'tthWdk')) then
        call  phase7dk(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      else
        write(6,*) 'Unanticipated process in gen7.f!'
        stop
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
      q(9,nu)=p9(nu)
      enddo 
      
      wt7=xjac*pswt
      
      if (debug) write(6,*) 'wt7 in gen7',wt7
      
      return

 999  wt7=0d0
      q(:,:)=0d0
      return 1
      end

