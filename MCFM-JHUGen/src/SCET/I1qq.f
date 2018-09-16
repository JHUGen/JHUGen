      subroutine xI1qiqi(z,I1qiqi)

      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp), intent(in) :: z
      real(dp) :: I1qiqi(0:2)

c--- coefficient of delta(1-z)
      I1qiqi(0) = -pisq/6.0_dp
c--- coefficient of 1/(1-z)+
      I1qiqi(1) = 1.0_dp+z**2
c--- regular term
      I1qiqi(2) = 1.0_dp-z-(1.0_dp+z**2)/(1.0_dp-z)*log(z)

      I1qiqi(:) = CF*I1qiqi(:)

      return
      end

!---
      function I1qig(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I1qig
      real(dp) :: fpqg

      fpqg = (1-z)**2+z**2

      I1qig = fpqg*(log((1-z)/z)-1) + 1.0_dp

      I1qig = TR*I1qig

      return
      end


