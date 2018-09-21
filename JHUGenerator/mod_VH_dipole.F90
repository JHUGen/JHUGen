!--YaofuZhou-----------------------------------------
module ModVHdipole
  implicit none

  public :: p_tilde

contains

subroutine p_tilde(pin,ptilde)
  use ModMisc
  implicit none

  real(8), intent(in) :: pin(1:4,1:10)
  real(8), intent(out) :: ptilde(1:4,1:9,1:2)
  integer :: i,j
  real(8) :: x
  real(8) :: K(1:4), Ktilde(1:4), KpKtilde(1:4)

  K = pin(:,1) + pin(:,2) - pin(:,10)

  do i=1,2

    ptilde(:,i,3-i)=pin(:,i)

    x = ( (pin(:,1).dot.pin(:,2)) - (pin(:,10).dot.pin(:,1)) - (pin(:,10).dot.pin(:,2)) ) / (pin(:,1).dot.pin(:,2))
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


real(8) function split_qiqi(alphas_dip,pq1,pq2,pg)
  use ModMisc
  use ModParameters
  use ModVHaux
  implicit none
    real(8), intent(in) :: alphas_dip
    real(8), intent(in) :: pq1(1:4)
    real(8), intent(in) :: pq2(1:4)
    real(8), intent(in) :: pg(1:4)
    real(8) :: x

    x = ( (pq1.dot.pq2) - (pg.dot.pq1) - (pg.dot.pq2) ) / (pq1.dot.pq2)

    if(x.gt.alpha_dip)then
      split_qiqi = 0d0
    else
      split_qiqi = 8d0*pi*cF*alphas_dip*( 2d0/(1d0-x) - (1d0+x) ) / x !1/x from first line of arXiv:0709.2881 EQ. 5.136
    endif

  return
end function split_qiqi



end module ModVHdipole