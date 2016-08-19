      double precision function asGamma1int(omsq)
      implicit none
C--   Author R.K. Ellis April 2012
C--   Integrand for NLO width with W-offshell 
      include 'constants.f'
      double precision omsq,mt,xi,ga,besq,asGamma1
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      asGamma1int=ga*xi/pi
     & /((1d0-xi*omsq)**2+ga**2)*asGamma1(mt,besq,omsq)
      return
      end
