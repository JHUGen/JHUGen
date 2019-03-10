      subroutine finitemtcorr(rescaling)
      implicit none
      double precision rescaling,tn
      double complex ftn
      include 'constants.f'
      include 'masses.f'
      
      tn = 4d0*(mt/hmass)**2
      if (tn .lt. (1d0)) then
         ftn = 0.5d0*(dlog((1d0+dsqrt(1d0-tn))
     &        /(1d0-dsqrt(1d0-tn)))-im*pi)**2
      else
         ftn = -2d0*(dasin(1.0d0/dsqrt(tn)))**2
      endif
      rescaling=cdabs((3d0*tn/4d0)*(2d0+(tn-1d0)*ftn))**2
      
      return
      end
      
