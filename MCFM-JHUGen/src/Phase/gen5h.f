c--- Generates 2->4 phase space with (12) a Breit-Wigner around
c--- the Higgs mass
      subroutine gen5h(r,p,wt5,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'process.f'
      include 'breit.f'
      include 'x1x2.f'
      integer nu
      double precision r(mxdim)
      double precision wt5,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision p(mxpart,4),rtshat
      double precision pswt,xjac
      double precision s12,wt12,ymax,yave
      include 'energy.f'

      wt5=0d0

      if (case .ne. 'HZZqgI') then
c--- this is the usual case      
        call breitw(r(9),0d0,sqrts**2,hmass,hwidth,s12,wt12)
      else
c--- this is the HZZqgI case (cf. gen4handc.f)     
        if (hmass .lt. mass2+mass3-hwidth*5d0) then
          call breitw(r(9),0d0,sqrts**2,hmass,10d0,s12,wt12)
        else
          call breitw(r(9),0d0,sqrts**2,hmass,hwidth,s12,wt12)
        endif
      endif
            
      rtshat=dsqrt(s12)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(10)-1d0)
      xjac=two*ymax*wt12
           
      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
c      write(6,*) 'problems with xx(1),xx(2) in gen4h',xx(1),xx(2)  
      return 1 
      endif

      p1(4)=-0.5d0*xx(1)*sqrts
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=-0.5d0*xx(1)*sqrts
      
      p2(4)=-0.5d0*xx(2)*sqrts
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=+0.5d0*xx(2)*sqrts

      call phase5h(r,p1,p2,p3,p4,p5,p6,p7,pswt) 

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      enddo 

      wt5=xjac*pswt/sqrts**2
      
      if (debug) write(6,*) 'wt5 in gen5h',wt5
      return

c 999  return 1
      end

