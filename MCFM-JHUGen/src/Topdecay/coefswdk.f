      subroutine coefswdk(qsqw,ct,c0)
      implicit none
C     Authors: John Campbell and R.Keith Ellis, April 2012
C     This gives the coefficients of the virtual corrections, c0, 
C     and the the integrated counterterm for 
C     w decay into massless quarks
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'includect.f'
      double precision ct,c0,qsqw,xl12

c--- this is computed in the 4d-hel scheme (scheme = 'dred')
      xl12=dlog(qsqw/musq)

c--- epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967
      if (includect) then
      ct=+2d0*epinv*(epinv2-xl12)+xl12**2+3d0*(epinv-xl12)+9d0-pisq
     &   +3d0*(aff-1d0-dlog(aff))-2d0*dlog(aff)**2
      else
      ct=0d0
      endif
c---------------------------------------------------------------
C adding term to check that the total w rate
c integrates to the right amount,as/pi,(ie result should be zero) 
c for testing only!
c      ct=ct-1.5d0
c---------------------------------------------------------------
      c0=-2d0*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &   -3d0*(epinv-xl12)+pisq-7d0
     
      return
      end

