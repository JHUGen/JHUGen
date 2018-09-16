      double precision function vol_mass(decaymass,s)
      implicit none
      double precision decaymass,b0,bp,bm,s,vol,ro,ddilog
      
      ro=4d0*decaymass**2/s
      b0=dsqrt(1d0-ro)
      bp=0.5d0*(1d0+b0)      
      bm=0.5d0*(1d0-b0)      

c      temp=(4d0*b0*(27d0 - 26d0*b0**2 + 3d0*b0**4)
c     . + 3d0*(-1d0 + b0**2)**2*dlog((-1d0 + b0)**2)**2 + 
c     . 6d0*(-1d0 + b0**2)**2*dlog(1d0-b0)*(1d0-4d0*dlog(1d0 + b0))+ 
c     . 6d0*(-1d0 + b0**2)*dlog(1d0 + b0)*
c     . (9d0 - b0**4 + (-1d0 + b0**2)*dlog(16d0)
c     . - 2d0*(-1d0 + b0**2)*dlog(1d0 + b0)) 
c     . + 3d0*(-1d0 + b0**2)*dlog((-1d0 + b0)**2)*
c     . (-8d0 + b0**4 + dlog(16d0) - b0**2*(1d0 + dlog(16d0)) + 
c     . 4d0*(-1d0 + b0**2)*dlog(1d0 + b0)) + 
c     . 24d0*(-1d0 + b0**2)**2*Ddilog((1d0 - b0)/2d0) - 
c     . 24d0*(-1d0 + b0**2)**2*Ddilog((1d0 + b0)/2d0))/192d0
c     
c--- normalize to massless case (1/12)
c      temp=temp/(1d0/12d0)
 
      vol_mass=b0*(1d0+5d0*ro+0.75d0*ro**2)
     . +1.5d0*ro**2*(Ddilog(bm)-Ddilog(bp))
     . -ro*dlog(bp/bm)*(24d0+6d0*ro-3d0*ro**2+6d0*ro*dlog(bp*bm))/8d0
      vol_mass=vol_mass*vol(s,4)
      
      return
      end
      
