!---
      subroutine xI1gg(z,I1gg)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp),intent(in) :: z
      real(dp) :: I1gg(0:2),omz

      real(dp) :: fpgg

      omz=one-z
      fpgg = 2*(omz+z**2)**2/z/omz

      I1gg(0) = -pisq/6.0_dp
      I1gg(1) = 2.0_dp*(omz+z**2)**2/z
      I1gg(2) = - fpgg*log(z)
      if(z.eq.1.0_dp) I1gg(2) = 2.0_dp

      I1gg(:) = CA*I1gg(:)

      return
      end

!---
      function I1gqi(z)
      implicit none
      include 'types.f'
      include 'constants.f'

      real(dp), intent(in) :: z
      real(dp) :: I1gqi

      real(dp) :: fpgq

      fpgq = (1.0_dp+(1.0_dp-z)**2)/z

      I1gqi = fpgq*log((1.0_dp-z)/z) + z

      I1gqi = I1gqi * CF

      return
      end
