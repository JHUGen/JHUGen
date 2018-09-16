      subroutine finitemtcorr(rescaling)
      implicit none
      include 'types.f'
      
      real(dp):: rescaling,tn
      complex(dp):: ftn
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      
      tn = 4._dp*(mt/hmass)**2
      if (tn < (1._dp)) then
         ftn = 0.5_dp*(log((1._dp+sqrt(1._dp-tn))
     &        /(1._dp-sqrt(1._dp-tn)))-im*pi)**2
      else
         ftn = -2._dp*(asin(1.0_dp/sqrt(tn)))**2
      endif
      rescaling=abs((3._dp*tn/4._dp)*(2._dp+(tn-1._dp)*ftn))**2
      
      return
      end
      
