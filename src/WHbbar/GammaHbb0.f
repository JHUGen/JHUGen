      double precision function GammaHbb0(Msq,mbsq)
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      double precision Msq,mbsq,beta,besq
c      write(6,*) 'GammaHbb0:Msq',Msq
c      write(6,*) 'GammaHbb0:mbsq',mbsq
c      pause
      besq=1d0-4d0*mbsq/Msq
      beta=sqrt(besq)
      GammaHbb0=3d0/4d0/pi*mbsq*Gf/rt2*sqrt(Msq)*beta**3
      return
      end

