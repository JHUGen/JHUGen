      subroutine hardqq(Qsq,musq,hard)
      implicit none
!    Hard function for qqbar in units of as/2/pi
      include 'types.f'
      include 'constants.f'
      real(dp),intent(in)::Qsq,musq
      real(dp),intent(out)::hard(2)
      complex(dp)::coeff(2)
      call qqcoeff(Qsq,musq,coeff)
! factors of 1/2 adjust for as/4/pi -> as/2/pi
      hard(1)=real(coeff(1),kind=dp)
      hard(2)=half**2*real(coeff(1)*conjg(coeff(1)),kind=dp)
     & +half*real(coeff(2),kind=dp)
      return
      end
