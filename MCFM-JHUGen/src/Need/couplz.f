      subroutine couplz(xw)
      implicit none
      include 'types.f'
             
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zcouple.f'
      include 'ewcharge.f'
c---calculate the couplings as given in Kunszt and Gunion
c---Modified to notation of DKS (ie divided by 2*sw*cw)
c---xw=sin^2 theta_w
      integer:: j
      real(dp):: xw
      sin2w=two*sqrt(xw*(1._dp-xw))
      do j=1,nf
      l(j)=(tau(j)-two*Q(j)*xw)/sin2w
      r(j)=      (-two*Q(j)*xw)/sin2w
      enddo

      le=(-1._dp-two*(-1._dp)*xw)/sin2w
      re=(-two*(-1._dp)*xw)/sin2w

      ln=(+1._dp-two*(+0._dp)*xw)/sin2w
      rn=0._dp

      return
      end
