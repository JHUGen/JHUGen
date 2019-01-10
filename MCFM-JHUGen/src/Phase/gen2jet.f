      subroutine gen2jet(r,p,wt2,*)
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'jetcuts.f'
      include 'phasemin.f'
      include 'leptcuts.f'
      include 'part.f'
      include 'x1x2.f'
      include 'first.f'
      include 'energy.f'
      include 'nproc.f'
!--- Modified July 11 by CW to switch between photons and jets appropriately
      integer j,nu
      double precision r(mxdim),p(mxpart,4),zmin,zmax,z
      double precision yave,ydif,xjac,y3,y4,phi,wt0,wt2,
     & ydifmin,ydifmax,yavemin,yavemax,xtsq,pt,xt,ptpar
      parameter(wt0=1d0/16d0/pi)
      save ptpar
!$omp threadprivate(ptpar)

      if (first) then
        first=.false.
c--- for dirgam, hflgam and gamgam, generate using photon pt as cutoff
        if ((nproc.eq.280) .or. (nproc.eq.285)
     &  .or.(nproc.eq.283) .or. (nproc.eq.284)) then
           ptpar=gammpt
           if ((part .eq. 'real').and.(nproc .eq. 285)) ptpar=gammpt2
        else
           call read_jetcuts(ptpar,etajetmin,etajetmax)
        endif
      endif

C    PS = 1/(16 pi) dxt^2 d phi/(2 pi) dyave dystar

      do j=1,mxpart
      do nu=1,4
      p(j,nu)=0d0
      enddo
      enddo

      wt2=0d0

c      xtsqmin=(2d0*ptjetmin/sqrts)**2
c      xjac=1d0-xtsqmin
c      xtsq=xtsqmin+xjac*r(3)


c      zmax=(0.5d0*sqrts/ptjetmin)**2
c      zmin=1d0
c      z=zmin+(zmax-zmin)*r(3)
c      xtsq=1d0/z
c      xjac=xtsq**2*(zmax-zmin)

c      xt=dsqrt(xtsq)
c      pt=0.5d0*sqrts*xt

      if (part .eq. 'real') then
        call genpt(r(3),ptpar,.false.,pt,xjac)
      else
        call genpt(r(3),ptpar,.true.,pt,xjac)
      endif
      xjac=xjac*8d0/sqrts**2
      xt=2d0*pt/sqrts
      xtsq=xt**2

      ydifmax=0.5d0*log((2d0-xtsq+2d0*dsqrt(1d0-xtsq))/xtsq)
      ydifmin=-ydifmax

      ydif=ydifmin+(ydifmax-ydifmin)*r(1)
      xjac=xjac*(ydifmax-ydifmin)

      yavemin=dlog(xt*cosh(ydif))
      yavemax=-yavemin
      yave=yavemin+(yavemax-yavemin)*r(2)
      xjac=xjac*(yavemax-yavemin)

      y3=yave+ydif
      y4=yave-ydif

      phi=2d0*pi*r(4)

      xx(1)=0.5d0*xt*(exp(+y3)+exp(+y4))
      xx(2)=0.5d0*xt*(exp(-y3)+exp(-y4))

      if (xx(1)*xx(2) .gt. 1d0) then
      write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1)*xx(2)
      write(6,*) 'xx(1),xx(2)',xx(1),xx(2)
      return 1
      endif
      if   ((xx(1) .gt. 1d0)
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)
     & ) then
c      write(6,*) 'problems with xx(1),xx(2) in gen2jet',xx(1),xx(2)
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

      p(3,4)=+pt*cosh(y3)
      p(3,1)=+pt*sin(phi)
      p(3,2)=+pt*cos(phi)
      p(3,3)=+pt*sinh(y3)

      p(4,4)=+pt*cosh(y4)
      p(4,1)=-pt*sin(phi)
      p(4,2)=-pt*cos(phi)
      p(4,3)=+pt*sinh(y4)

      wt2=wt0*xjac
      return

      end
