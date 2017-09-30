      subroutine gen2(r,p,wt2,*)
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4
c----
c---- if 'nodecay' is true, then the vector boson decay into massless
c---- particles is not included and 2 less integration variables
c---- are required
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'nodecay.f'
      include 'breit.f'
      include 'x1x2.f'
      integer j,nu
      double precision r(mxdim),p(mxpart,4),rdk1,rdk2
      double precision ymax,yave,ydif,xjac,y3,y4,phi,wt0,wt2,w3
      double precision pt,s34,rtshat,udif
      include 'energy.f'
      parameter(wt0=1d0/16d0/pi)


      do j=1,mxpart
      do nu=1,4
      p(j,nu)=0d0
      enddo
      enddo

      wt2=0d0

c--- dummy values if there's no decay
      if (nodecay) then
        rdk1=0.5d0
        rdk2=0.5d0
      else
        rdk1=r(3)
        rdk2=r(4)
      endif

      if (n3.eq.0) then
         w3=(wsqmax-wsqmin)
         s34=(wsqmax-wsqmin)*r(1)+wsqmin
      elseif (n3.eq.1) then
         call breitw(r(1),wsqmin,wsqmax,mass3,width3,s34,w3)
      endif
      rtshat=dsqrt(s34)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(2)-1d0)

c----udif==tanh(ydif)
      udif=(two*rdk1-1d0)
      ydif=half*dlog((1d0+udif)/(1d0-udif))
      xjac=four*ymax

      y3=yave+ydif
      y4=yave-ydif

      xjac=xjac*w3
      phi=2d0*pi*rdk2

      pt=rtshat/(2d0*dcosh(ydif))
      xx(1)=rtshat/sqrts*dexp(+yave)
      xx(2)=rtshat/sqrts*dexp(-yave)

      if   ((xx(1) .gt. 1d0)
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
        write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1),xx(2)
      return 1
      endif

      p(1,4)=-0.5d0*xx(1)*sqrts
      p(1,1)=0d0
      p(1,2)=0d0
      p(1,3)=-0.5d0*xx(1)*sqrts

      p(2,4)=-0.5d0*xx(2)*sqrts
      p(2,1)=0d0
      p(2,2)=0d0
      p(2,3)=+0.5d0*xx(2)*sqrts

      p(3,4)=+pt*dcosh(y3)
      p(3,1)=+pt*dsin(phi)
      p(3,2)=+pt*dcos(phi)
      p(3,3)=+pt*dsinh(y3)

      p(4,4)=+pt*dcosh(y4)
      p(4,1)=-pt*dsin(phi)
      p(4,2)=-pt*dcos(phi)
      p(4,3)=+pt*dsinh(y4)

      wt2=wt0*xjac/sqrts**2
      return

      end
