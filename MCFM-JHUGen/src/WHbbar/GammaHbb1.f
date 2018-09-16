      function GammaHbb1(Msq,mbsq)
      implicit none
      include 'types.f'
      real(dp):: GammaHbb1      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      real(dp):: Msq,mbsq,beta,besq,ddilog,GammaHbb0,BigL,omb,opb
c----Formula taken from Braaten and Leveille, PR D22, 715, 1980, Eq.(14)
      besq=1._dp-4._dp*mbsq/Msq
      beta=sqrt(besq)
      omb=0.5_dp*(1._dp-beta)
      opb=0.5_dp*(1._dp+beta)
      BigL=log(opb/omb)
      GammaHbb1=GammaHbb0(Msq,mbsq)
     & *(one+ason2pi*CF*(6._dp-0.75_dp*(one+besq)/besq+6._dp*log(mbsq/Msq)
     & -8._dp*log(beta)
     & +(5._dp/beta-2._dp*beta+0.375_dp*(one-besq)**2/beta**3)*BigL
     & +(one+besq)/beta*(4._dp*BigL*log(opb/beta))
     & -2._dp*log(opb)*log(omb)+8._dp*ddilog(omb/opb)-4._dp*ddilog(omb)))
c--- This is the leading-log approximation
c      GammaHbb1=GammaHbb0(Msq,mbsq)*(
c     & 1._dp+ason2pi*CF*(4.5_dp-3._dp*log(Msq/mbsq)))
      return
      end
