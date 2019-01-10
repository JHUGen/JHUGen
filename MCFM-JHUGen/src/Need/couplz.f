      subroutine couplz(xw)
      implicit none       
      include 'constants.f'
      include 'zcouple.f'
      include 'ewcharge.f'
c---calculate the couplings as given in Kunszt and Gunion
c---Modified to notation of DKS (ie divided by 2*sw*cw)
c---xw=sin^2 theta_w
      integer j
      double precision xw
      sin2w=two*sqrt(xw*(1d0-xw))
      do j=1,nf
      l(j)=(tau(j)-two*Q(j)*xw)/sin2w
      r(j)=      (-two*Q(j)*xw)/sin2w
      enddo

      le=(-1d0-two*(-1d0)*xw)/sin2w
      re=(-two*(-1d0)*xw)/sin2w

      ln=(+1d0-two*(+0d0)*xw)/sin2w
      rn=0d0

      return
      end
