      subroutine gen3jetgaga(r,p,wt3,*)
C---generate three particle phase space and x1,x2 integration
!--- Setup to generate p3 and p4 with photon cuts 
C---p1+p2 --> p3+p4+p5
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'leptcuts.f'
      include 'reset.f'
      include 'x1x2.f'
      include 'first.f'
      integer j,nu
      double precision r(mxdim),p(mxpart,4),
     . xjac,y3,y4,y5,phi,phi34,wt0,wt3,etamax,
     . pt3,pt4,pt5,rtson2,cphi,sphi,cphi34,sphi34
      include 'energy.f'
      parameter(wt0=1d0/512d0/pi**3)
      double precision hmin,hmax,delh,h,tmp
      double precision ptjetmin,etajetmin,etajetmax
      save ptjetmin,etajetmin,etajetmax
!$omp threadprivate(ptjetmin,etajetmin,etajetmax)
      
      wt3=0d0
      
      if (first .or. reset) then
        first=.false.
        reset=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
      endif
      
      do j=6,mxpart     
      do nu=1,4     
      p(j,nu)=0d0
      enddo     
      enddo     

      phi=2d0*pi*r(1)
      cphi=dcos(phi)
      sphi=dsin(phi)
      phi34=2d0*pi*r(2)
      cphi34=dcos(phi34)
      sphi34=dsin(phi34)
      xjac=sqrts**2

      rtson2=0.5d0*sqrts

     
      hmin=1d0/dsqrt(rtson2**2+gammpt**2)
      hmax=1d0/gammpt
      delh=hmax-hmin

      h=hmin+r(6)*delh  
      pt3=dsqrt(1d0/h**2-gammpt**2)
      xjac=xjac*delh/h**3/rtson2**2
      etamax=sqrts/2d0/pt3
      if (etamax**2 .le. 1d0) then
c        write(6,*) 'etamax**2 .le. 1d0 in gen3jetgaga.f',etamax**2 
        return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))        
      etamax=min(etamax,10d0)
      y3=etamax*(2d0*r(3)-1d0)
      xjac=xjac*2d0*etamax
      
      hmin=1d0/dsqrt(rtson2**2+gammpt2**2)
      hmax=1d0/gammpt2
      delh=hmax-hmin

      h=hmin+r(7)*delh  
      pt4=dsqrt(1d0/h**2-gammpt2**2)
      xjac=xjac*delh/h**3/rtson2**2
      etamax=sqrts/2d0/pt4
      if (etamax**2 .le. 1d0) then
c        write(6,*) 'etamax**2 .le. 1d0 in gen3jetgaga.f',etamax**2 
        return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))        
      etamax=min(etamax,10d0)
      y4=etamax*(2d0*r(4)-1d0)
      xjac=xjac*2d0*etamax

      
      p(4,1)=pt4*sphi
      p(4,2)=pt4*cphi

      p(3,1)=pt3*(+cphi34*sphi+sphi34*cphi)
      p(3,2)=pt3*(-sphi34*sphi+cphi34*cphi)
           
      p(5,1)=-p(4,1)-p(3,1)
      p(5,2)=-p(4,2)-p(3,2)
      pt5=dsqrt(p(5,1)**2+p(5,2)**2)
      etamax=sqrts/2d0/pt5
      if (etamax**2 .le. 1d0) then
c        write(6,*) 'etamax**2 .le. 1d0 in gen3jetgaga.f',etamax**2 
        return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))        
      etamax=min(etamax,10d0)
      y5=etamax*(2d0*r(5)-1d0)
      xjac=xjac*2d0*etamax

      
      xx(1)=half*(pt3*exp(+y3)+pt4*exp(+y4)+pt5*exp(+y5))/rtson2
      xx(2)=half*(pt3*exp(-y3)+pt4*exp(-y4)+pt5*exp(-y5))/rtson2
      

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
c      write(6,*) 'problems with xx(1),xx(2) in gen3',xx(1),xx(2)  
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

      p(3,4)=+pt3*cosh(y3)
      p(3,3)=+pt3*sinh(y3)

      p(4,4)=+pt4*cosh(y4)
      p(4,3)=+pt4*sinh(y4)

      p(5,4)=+pt5*cosh(y5)
      p(5,3)=+pt5*sinh(y5)


c--- randomize which of p3 and p4 is softest photon      
      if (r(8) .lt. 0.5d0) then
        do nu=1,4
        tmp=p(3,nu)
        p(3,nu)=p(4,nu)
        p(4,nu)=tmp
        enddo
      endif

      wt3=wt0*xjac
      return

      end
