      subroutine gen2(r,p,wt2,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4
c----
c---- if 'nodecay' is true, then the vector boson decay into massless
c---- particles is not included and 2 less integration variables
c---- are required

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'limits.f'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'nodecay.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'debug.f'
      integer:: j,nu
      real(dp):: r(mxdim),p(mxpart,4),rdk1,rdk2,s34min
      real(dp):: ymax,yave,ydif,xjac,y3,y4,phi,wt0,wt2,w3
      real(dp):: pt,s34,rtshat,udif
      include 'energy.f'
      parameter(wt0=one/16._dp/pi)

      p(:,:)=zip

      wt2=zip

c--- dummy values if there's no decay
      if (nodecay) then
        rdk1=0.5_dp
        rdk2=0.5_dp
      else
        rdk1=r(3)
        rdk2=r(4)
      endif

      if (n3==0) then
         w3=one
         s34min=max(wsqmin,one) ! ensure minimum value of m(34)>1 GeV
         call pick(2,s34,s34min,wsqmax,r(1),w3)
      elseif (n3==1) then
         call breitw(r(1),wsqmin,wsqmax,mass3,width3,s34,w3)
      endif
      rtshat=sqrt(s34)
      ymax=log(sqrts/rtshat)
      yave=ymax*(two*r(2)-one)

c----udif==tanh(ydif)
      udif=(two*rdk1-one)
      ydif=half*log((one+udif)/(one-udif))
      xjac=four*ymax
      y3=yave+ydif
      y4=yave-ydif

      xjac=xjac*w3
      phi=2._dp*pi*rdk2

      pt=rtshat/(2._dp*cosh(ydif))
      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) > one)
     & .or. (xx(2) > one)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
        if (debug) write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1),xx(2)
        return 1
      endif

      p(1,4)=-0.5_dp*xx(1)*sqrts
      p(1,1)=0._dp
      p(1,2)=0._dp
      p(1,3)=-0.5_dp*xx(1)*sqrts

      p(2,4)=-0.5_dp*xx(2)*sqrts
      p(2,1)=0._dp
      p(2,2)=0._dp
      p(2,3)=+0.5_dp*xx(2)*sqrts

      p(3,4)=+pt*cosh(y3)
      p(3,1)=+pt*sin(phi)
      p(3,2)=+pt*cos(phi)
      p(3,3)=+pt*sinh(y3)

      p(4,4)=+pt*cosh(y4)
      p(4,1)=-pt*sin(phi)
      p(4,2)=-pt*cos(phi)
      p(4,3)=+pt*sinh(y4)

      wt2=wt0*xjac/sqrts**2
      return

      end
