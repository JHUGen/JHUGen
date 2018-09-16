      subroutine coefswdk(qsqw,ct,c0)
      implicit none
      include 'types.f'
      
C     Authors: John Campbell and R.Keith Ellis, April 2012
C     This gives the coefficients of the virtual corrections, c0, 
C     and the the integrated counterterm for 
C     w decay into massless quarks
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'includect.f'
      real(dp):: ct,c0,qsqw,xl12

c--- this is computed in the 4.e-_dphel scheme (scheme = 'dred')
      xl12=log(qsqw/musq)

c--- epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967
      if (includect) then
      ct=+2._dp*epinv*(epinv2-xl12)+xl12**2+3._dp*(epinv-xl12)+9._dp-pisq
     &   +3._dp*(aff-1._dp-log(aff))-2._dp*log(aff)**2
      else
      ct=0._dp
      endif
c---------------------------------------------------------------
C adding term to check that the total w rate
c integrates to the right amount,as/pi,(ie result should be zero) 
c for testing only!
c      ct=ct-1.5_dp
c---------------------------------------------------------------
      c0=-2._dp*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &   -3._dp*(epinv-xl12)+pisq-7._dp
     
      return
      end

