      subroutine genpt(xr,ptmin,exact,pt,xjac)
      implicit none
c--- given a random number xr between 0 and 1, generate a transverse
c--- momentum pt according to the flag exact:
c---    exact=true:  generate according to dpt/pt down to ptmin
c---    exact=false: generate down to pt=0 with a shape determined by ptmin
c---
c--- returns: pt and xjac, the Jacobian of the transformation from pt dpt to dxr
      logical exact
      double precision xr,ptmin,pt,xjac,hmin,hmax,h,delh,ptmax
      include 'energy.f'

      ptmax=sqrts/2d0

      if (exact) then
        hmin=1d0/ptmax
        hmax=1d0/ptmin
        delh=hmax-hmin
        h=hmin+xr*delh
        pt=1d0/h
        xjac=delh/h**3
      else
        pt=2d0*ptmin*ptmax*xr/(2d0*ptmin+ptmax*(1d0-xr))
        xjac=pt*ptmax/2d0/ptmin/(2d0*ptmin+ptmax)*(2d0*ptmin+pt)**2
      endif

      return
      end

