      double precision function Gamma0int(omsq)
      implicit none
C--   Author R.K. Ellis April 2012
C--   Integrand for width with W-offshell 
      include 'constants.f'
      double precision omsq,mt,besq,xi,ga,Gamma0
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      Gamma0int=(ga*xi/pi)/((1d0-xi*omsq)**2+ga**2)*Gamma0(mt,besq,omsq)
      return
      end

