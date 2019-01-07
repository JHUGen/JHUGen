      double precision function vol3_mass(decaymass,s)
      implicit none
      double precision decaymass,b0,bp,bm,s,ro,vol
c  Volume of phase space for 2-->p3(decaymass)+p4(decaymass)+p5(0)      
      if (decaymass .eq. 0d0) then
      vol3_mass=1d0
      else
      ro=4d0*decaymass**2/s
      b0=dsqrt(1d0-ro)
      bp=0.5d0*(1d0+b0)      
      bm=0.5d0*(1d0-b0)      

      vol3_mass=(3d0-b0**2)/2d0*b0
     . -(1d0-b0**2)*(3d0+b0**2)/4d0*log(bp/bm)
      endif
      vol3_mass=vol3_mass*vol(s,3)
      return
      end
      
