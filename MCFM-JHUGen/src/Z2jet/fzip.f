      double complex function fzip(s12,s34,s56) 
      implicit none
      double complex I3m,lnrat
      double precision s12,s34,s56,Del12,Del34,Del56,Del3

C     Implementation of BDKW Nucl Phys. 513, 3 (1998)
C     Eqn. (B15).
      Del12=s12-s34-s56
      Del34=s34-s12-s56
      Del56=s56-s12-s34
      Del3=s12*Del12+s34*Del34+s56*Del56
       

      fzip=Dcmplx((3d0*s34*Del34/Del3**2-1d0/Del3)*s12*s56)
     . *I3m(s12,s34,s56)
     . +Dcmplx((3d0*s56*Del56/Del3**2-0.5d0/Del3)*s12)
     . *lnrat(-s12,-s34)
     . +Dcmplx((3d0*s12*Del12/Del3**2-0.5d0/Del3)*s56)
     . *lnrat(-s56,-s34)
     . -Dcmplx(0.5d0*Del34/Del3)
      end

