      function Gamma0int(omsq)
      implicit none
      include 'types.f'
      real(dp):: Gamma0int
      
C--   Author R.K. Ellis April 2012
C--   Integrand for width with W-offshell 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: omsq,mt,besq,xi,ga,Gamma0
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      Gamma0int=(ga*xi/pi)/((1._dp-xi*omsq)**2+ga**2)*Gamma0(mt,besq,omsq)
      return
      end

