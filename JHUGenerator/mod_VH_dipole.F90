!--YaofuZhou-----------------------------------------
module ModVHdipole
  implicit none

  public :: p_tilde
  public :: split_qqg
  public :: split_qgq

contains


subroutine p_tilde(pin,ptilde)
  use ModMisc
  implicit none

  real(8), intent(in) :: pin(1:4,1:10)
  real(8), intent(out) :: ptilde(1:4,1:9,1:2)
  integer :: i,j
  real(8) :: x
  real(8) :: K(1:4), Ktilde(1:4), KpKtilde(1:4)

  x = ( (pin(:,1).dot.pin(:,2)) - (pin(:,10).dot.pin(:,1)) - (pin(:,10).dot.pin(:,2)) ) / (pin(:,1).dot.pin(:,2))

  K = pin(:,1) + pin(:,2) - pin(:,10)

  do i=1,2

    ptilde(:,i,3-i)=pin(:,i)

    ptilde(:,i,i) = pin(:,i) * x

    Ktilde = ptilde(:,i,i) + pin(:,3-i) !if i=1, pin(:,2); if i=2, pin(:,1)

    do j=3,9
      KpKtilde = K+Ktilde
      ptilde(:,j,i) = pin(:,j) - 2d0 * (pin(:,j).dot.KpKtilde) / (KpKtilde.dot.KpKtilde) * KpKtilde &
                               + 2d0 * (pin(:,j).dot.K)        / (K.dot.K)               * Ktilde
    enddo

  enddo

  return
end subroutine p_tilde

real(8) function split_qqg(alphas_dip,pa,pb,pr)
  use ModMisc
  use ModParameters
  use ModVHaux
  implicit none
    real(8), intent(in) :: alphas_dip
    real(8), intent(in) :: pa(1:4)!emitter
    real(8), intent(in) :: pb(1:4)!the other incoming
    real(8), intent(in) :: pr(1:4)!radiated
    real(8) :: x,v

    v = (pr.dot.pa) / (pa.dot.pb)

    x = ( (pa.dot.pb) - (pr.dot.pa) - (pr.dot.pb) ) / (pa.dot.pb)

    if(v.gt.alpha_dip)then
      split_qqg = 0d0
    else
      split_qqg = 8d0*pi*cF*alphas_dip*( 2d0/(1d0-x) - (1d0+x) ) / x !1/x from first line of arXiv:0709.2881 EQ. 5.136
    endif

  return
end function split_qqg


real(8) function split_qgq(alphas_dip,pa,pb,pr)
  use ModMisc
  use ModParameters
  use ModVHaux
  implicit none
    real(8), intent(in) :: alphas_dip
    real(8), intent(in) :: pa(1:4)!emitter
    real(8), intent(in) :: pb(1:4)!the other incoming
    real(8), intent(in) :: pr(1:4)!radiated
    real(8) :: x,v

    v = (pr.dot.pa) / (pa.dot.pb)
    x = ( (pa.dot.pb) - (pr.dot.pa) - (pr.dot.pb) ) / (pa.dot.pb)

    if(v.gt.alpha_dip)then
      split_qgq = 0d0
    else
      split_qgq = 8d0*pi*Tr*alphas_dip*( 1d0 - 2d0*x*(1d0-x) ) / x !1/x from first line of arXiv:0709.2881 EQ. 5.136
    endif

  return
end function split_qgq


end module ModVHdipole