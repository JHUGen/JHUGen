      double precision function GammaHbb1(Msq,mbsq)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      double precision Msq,mbsq,beta,besq,ddilog,GammaHbb0,BigL,omb,opb
c----Formula taken from Braaten and Leveille, PR D22, 715, 1980, Eq.(14)
      besq=1d0-4d0*mbsq/Msq
      beta=sqrt(besq)
      omb=0.5d0*(1d0-beta)
      opb=0.5d0*(1d0+beta)
      BigL=log(opb/omb)
      GammaHbb1=GammaHbb0(Msq,mbsq)
     & *(1d0+ason2pi*CF*(6d0-0.75d0*(1d0+besq)/besq+6d0*log(mbsq/Msq)
     & -8d0*log(beta)
     & +(5d0/beta-2d0*beta+0.375d0*(1d0-besq)**2/beta**3)*BigL
     & +(1d0+besq)/beta*(4d0*BigL*log(opb/beta))
     & -2d0*log(opb)*log(omb)+8d0*ddilog(omb/opb)-4d0*ddilog(omb)))
c--- This is the leading-log approximation
c      GammaHbb1=GammaHbb0(Msq,mbsq)*(
c     & 1d0+ason2pi*CF*(4.5d0-3d0*log(Msq/mbsq)))
      return
      end
