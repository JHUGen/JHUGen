      function GammaHbb0(Msq,mbsq)
      implicit none
      include 'types.f'
      real(dp):: GammaHbb0      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      real(dp):: Msq,mbsq,beta,besq
c      write(6,*) 'GammaHbb0:Msq',Msq
c      write(6,*) 'GammaHbb0:mbsq',mbsq
c      pause
      besq=1._dp-4._dp*mbsq/Msq
      beta=sqrt(besq)
      GammaHbb0=3._dp/4._dp/pi*mbsq*Gf/rt2*sqrt(Msq)*beta**3
      return
      end

