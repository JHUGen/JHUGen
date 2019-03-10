      double precision function Gamma0(mt,besq,omsq)
C--   Author: John M. Campbell and R.K. Ellis, January 2012  
C--   Taken from formula (2) of
C--   Fermilab-PUB-12-078-T
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      double precision mt,omsq,besq,Gammainfty,f,P3b
      Gammainfty=GF*mt**3/(8d0*rt2*pi)
      P3b=0.5d0*dsqrt(1d0+omsq**2+besq**2-2d0*(omsq+besq+omsq*besq))
      f=(1d0-besq)**2+omsq*(1d0+besq)-2d0*omsq**2
      Gamma0=Gammainfty*2d0*P3b*f
      return
      end

