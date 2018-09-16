      function fzip(s12,s34,s56)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: fzip
      complex(dp):: I3m,lnrat
      real(dp):: s12,s34,s56,Del12,Del34,Del56,Del3
      include 'cplx.h'

C     Implementation of BDKW Nucl Phys. 513, 3 (1998)
C     Eqn. (B15).
      Del12=s12-s34-s56
      Del34=s34-s12-s56
      Del56=s56-s12-s34
      Del3=s12*Del12+s34*Del34+s56*Del56


      fzip=((3._dp*s34*Del34/Del3**2-1._dp/Del3)*s12*s56)*I3m(s12,s34,s56)
     & +((3._dp*s56*Del56/Del3**2-0.5_dp/Del3)*s12)*lnrat(-s12,-s34)
     & +((3._dp*s12*Del12/Del3**2-0.5_dp/Del3)*s56)*lnrat(-s56,-s34)
     & -cplx2(0.5_dp*Del34/Del3,zip)
      end

