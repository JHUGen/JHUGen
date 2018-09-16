      subroutine genpt(xr,ptmin,exact,pt,xjac)
      implicit none
      include 'types.f'
c--- given a random number xr between 0 and 1, generate a transverse
c--- momentum pt according to the flag exact:
c---    exact=true:  generate according to dpt/pt down to ptmin
c---    exact=false: generate down to pt=0 with a shape determined by ptmin
c---
c--- returns: pt and xjac, the Jacobian of the transformation from pt dpt to dxr
      logical:: exact
      real(dp):: xr,ptmin,pt,xjac,hmin,hmax,h,delh,ptmax
      include 'energy.f'

      ptmax=sqrts/2._dp

      if (exact) then
        hmin=1._dp/ptmax
        hmax=1._dp/ptmin
        delh=hmax-hmin
        h=hmin+xr*delh
        pt=1._dp/h
        xjac=delh/h**3
      else
        pt=2._dp*ptmin*ptmax*xr/(2._dp*ptmin+ptmax*(1._dp-xr))
        xjac=pt*ptmax/2._dp/ptmin/(2._dp*ptmin+ptmax)*(2._dp*ptmin+pt)**2
      endif

      return
      end

