      function asGamma1int(omsq)
      implicit none
      include 'types.f'
      real(dp):: asGamma1int
      
C--   Author R.K. Ellis April 2012
C--   Integrand for NLO width with W-offshell 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: omsq,mt,xi,ga,besq,asGamma1
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      asGamma1int=ga*xi/pi
     & /((1._dp-xi*omsq)**2+ga**2)*asGamma1(mt,besq,omsq)
      return
      end
