      double complex function fmtfull(s12,s34,s56,mtsq) 
      implicit none
      include 'scale.f'
      double precision mtsq,s12,s34,s56,Del12,Del34,Del56,Del3
      double complex qlI3,bdiff

C     Implementation of BDKW Nucl Phys. 513, 3 (1998)
C     Eqn. (B.15) upgraded so that is gives the full mass
C     dependence.
      Del12=s12-s34-s56
      Del34=s34-s12-s56
      Del56=s56-s12-s34
      Del3=s12*Del12+s34*Del34+s56*Del56
       
      fmtfull=
     . Dcmplx((3d0*s12*s56*s34*Del34/Del3**2-(s12*s56-mtsq*del34)/Del3))
     . *(-qlI3(s12,s34,s56,mtsq,mtsq,mtsq,musq,0))
     . +Dcmplx((3d0*s56*Del56/Del3**2-0.5d0/Del3)*s12)
     . *bdiff(s34,s12,mtsq)
     . +Dcmplx((3d0*s12*Del12/Del3**2-0.5d0/Del3)*s56)
     . *bdiff(s34,s56,mtsq)
     . -Dcmplx(0.5d0*Del34/Del3)
      end

