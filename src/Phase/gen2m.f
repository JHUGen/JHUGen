      subroutine gen2m(r,p,wt2,*)
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'mxdim.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'breit.f'
      include 'x1x2.f'
      integer j,nu

      double precision r(mxdim),p(mxpart,4),mass
      double precision ymax,yave,ydif,xjac,y3,y4,phi,wt0,wt2,w3,vs
      double precision vsqmin,vsqmax,pt,s34,xmin,rtshat,udif,trmass,beta
      include 'energy.f'
      parameter(wt0=1d0/16d0/pi)

      do j=1,mxpart     
      do nu=1,4     
      p(j,nu)=0d0
      enddo     
      enddo   
      
      wt2=0d0
      
      if (n3 .eq. 1) then
c--- generate s34 according to a Breit-Wigner, for gg->H 
        call breitw(r(3),wsqmin,wsqmax,mass3,width3,s34,w3)
        xmin=0d0
      if     (hdecaymode .eq. 'bqba') then
        mass=mb
      elseif (hdecaymode .eq. 'tlta') then
        mass=mtau
      else
        write(6,*) 'Unanticipated hdecaymode in gen2m: ',hdecaymode
        stop
      endif
      else
c--- no resonance, for tt~,bb~,cc~
      mass=mass2
        vsqmax=1d0/(4d0*mass**2)
        vsqmin=1d0/sqrts**2
        xmin=vsqmin/vsqmax
        vs=(vsqmax-vsqmin)*r(3)+vsqmin
        s34=1/vs
        w3=(vsqmax-vsqmin)*s34**2
      endif
      
      rtshat=dsqrt(s34)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(1)-1d0)      
c----udif=tanh(ydif)
      beta=dsqrt(1d0-4d0*mass**2/s34)
      udif=beta*(two*r(2)-1d0)
      ydif=half*dlog((1d0+udif)/(1d0-udif))
      xjac=four*ymax*beta
          
      y3=yave+ydif
      y4=yave-ydif
          
      xjac=xjac*w3
      phi=2d0*pi*r(4)

      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)
      trmass=rtshat/(2d0*dcosh(ydif))


      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
        write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1),xx(2)  
        return 1 
      endif

      pt=dsqrt(trmass**2-mass**2)
          
      p(1,4)=-0.5d0*xx(1)*sqrts
      p(1,1)=0d0
      p(1,2)=0d0
      p(1,3)=-0.5d0*xx(1)*sqrts
      
      p(2,4)=-0.5d0*xx(2)*sqrts
      p(2,1)=0d0
      p(2,2)=0d0
      p(2,3)=+0.5d0*xx(2)*sqrts

      p(3,4)=+trmass*dcosh(y3)
      p(3,1)=+pt*dsin(phi)
      p(3,2)=+pt*dcos(phi)
      p(3,3)=+trmass*dsinh(y3)

      p(4,4)=+trmass*dcosh(y4)
      p(4,1)=-pt*dsin(phi)
      p(4,2)=-pt*dcos(phi)
      p(4,3)=+trmass*dsinh(y4)

      wt2=wt0*xjac/sqrts**2

c      write(6,*) 's34',s34
c      write(6,*) 's34-4d0*mass**2',s34-4d0*mass**2
c      write(6,*) 'wsqmax',wsqmax
c      write(6,*) 'ymax',ymax
c      write(6,*) 'wsqmin',wsqmin
c      write(6,*) 'y3',y3
c      write(6,*) 'y4',y4
c      write(6,*) 'xx(1)',xx(1)
c      write(6,*) 'xx(2)',xx(2)
c      write(6,*) 'trmass',trmass
c      write(6,*) 'mass',mass
c      write(6,*) 'pt',pt
c      write(6,*) 's12',2d0*(p(1,4)*p(2,4)-p(1,3)*p(2,3))

      return

      end
