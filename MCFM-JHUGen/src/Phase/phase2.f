      subroutine phase2(r,p1,p2,p3,p4,wt,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'x1x2.f'
c---- generate phase space for 2-->2 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4
c---- with all 2 pi's (ie 1/(2*pi)^2)
      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4)
      double precision cosphi,sinphi,u,phi,rtshat,costh,sinth
      double precision wt,wt0
      include 'energy.f'
      parameter(wt0=1d0/8d0/pi)
      rtshat=dsqrt(xx(1)*xx(2))*sqrts
C write out vectors in +,-,T,T notation
      u=r(3)
      costh=2d0*u-1d0
      sinth=dsqrt(1d0-costh**2)
      phi=two*pi*r(4)
      sinphi=dsin(phi)
      cosphi=dcos(phi)

      p3(4)=+half*sqrts*(u*xx(1)+(1d0-u)*xx(2))
      p3(1)=+half*sinth*sinphi*rtshat
      p3(2)=+half*sinth*cosphi*rtshat
      p3(3)=+half*sqrts*(u*xx(1)-(1d0-u)*xx(2))


      p4(4)=+half*sqrts*((1d0-u)*xx(1)+u*xx(2))
      p4(1)=-half*sinth*sinphi*rtshat
      p4(2)=-half*sinth*cosphi*rtshat
      p4(3)=+half*sqrts*((1d0-u)*xx(1)-u*xx(2))

c---debug
      write(6,*) 'p3',p3(4),p3(3),p3(2),p3(1)
      write(6,*) 'p4',p4(4),p4(3),p4(2),p4(1)

      p3(4)=+u*p1(4)+(one-u)*p2(4)
      p3(1)=+half*sinth*sinphi*rtshat
      p3(2)=+half*sinth*cosphi*rtshat
      p3(4)=+u*p1(3)+(one-u)*p2(3)


      p4(4)=+(one-u)*p1(4)+u*p2(4)
      p4(1)=-half*sinth*sinphi*rtshat
      p4(2)=-half*sinth*cosphi*rtshat
      p4(3)=+(one-u)*p1(3)+u*p2(3)

      write(6,*) 'p3',p3(4),p3(3),p3(2),p3(1)
      write(6,*) 'p4',p4(4),p4(3),p4(2),p4(1)
c      pause
      wt=wt0

      end

