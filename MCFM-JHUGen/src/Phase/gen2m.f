      subroutine gen2m(r,p,wt2,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'limits.f'
      include 'mxdim.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'debug.f'
      integer:: j,nu
      real(dp):: r(mxdim),p(mxpart,4),mass
      real(dp):: ymax,yave,ydif,xjac,y3,y4,phi,wt0,wt2,w3,vs
      real(dp):: vsqmin,vsqmax,pt,s34,xmin,rtshat,udif,trmass,beta
      include 'energy.f'
      parameter(wt0=one/16._dp/pi)

      p(:,:)=zip

      wt2=zip

      if (n3 == 1) then
         xmin=zip
!-----default mass to zero
         mass=zip
      if     (hdecaymode == 'bqba') then
        mass=mb
      elseif (hdecaymode == 'tlta') then
        mass=mtau
      elseif (hdecaymode == 'gaga') then
        mass=zip
      else
        write(6,*) 'Unanticipated hdecaymode in gen2m: ',hdecaymode
        stop
      endif
      wsqmin=max(wsqmin,four*mass**2)
c---  generate s34 according to a Breit-Wigner, for gg->H
      call breitw(r(3),wsqmin,wsqmax,mass3,width3,s34,w3)


      else
c--- no resonance, for tt~,bb~,cc~
      mass=mass2
        vsqmax=one/(four*mass**2)
        vsqmin=one/sqrts**2
        xmin=vsqmin/vsqmax
        vs=(vsqmax-vsqmin)*r(3)+vsqmin
        s34=1/vs
        w3=(vsqmax-vsqmin)*s34**2
      endif

      rtshat=sqrt(s34)
      ymax=log(sqrts/rtshat)
      yave=ymax*(two*r(1)-one)
c----udif=tanh(ydif)
      beta=sqrt(one-four*mass**2/s34)

      udif=beta*(two*r(2)-one)
      ydif=half*log((one+udif)/(one-udif))
      xjac=four*ymax*beta

      y3=yave+ydif
      y4=yave-ydif

      xjac=xjac*w3
      phi=two*pi*r(4)

      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)
      trmass=rtshat/(two*cosh(ydif))


      if   ((xx(1) > one)
     & .or. (xx(2) > one)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
        if (debug) write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1),xx(2)
        return 1
      endif

      pt=sqrt(trmass**2-mass**2)

      p(1,4)=-half*xx(1)*sqrts
      p(1,1)=zip
      p(1,2)=zip
      p(1,3)=-half*xx(1)*sqrts

      p(2,4)=-half*xx(2)*sqrts
      p(2,1)=zip
      p(2,2)=zip
      p(2,3)=+half*xx(2)*sqrts

      p(3,4)=+trmass*cosh(y3)
      p(3,1)=+pt*sin(phi)
      p(3,2)=+pt*cos(phi)
      p(3,3)=+trmass*sinh(y3)

      p(4,4)=+trmass*cosh(y4)
      p(4,1)=-pt*sin(phi)
      p(4,2)=-pt*cos(phi)
      p(4,3)=+trmass*sinh(y4)


      wt2=wt0*xjac/sqrts**2

c      write(6,*) 's34',s34
c      write(6,*) 's34-four*mass**2',s34-four*mass**2
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
c      write(6,*) 's12',two*(p(1,4)*p(2,4)-p(1,3)*p(2,3))

      return

      end
