!--YaofuZhou-----------------------------------------
module ModVHqg
  use ModParameters
  use ModMisc
  implicit none

public :: amp_VH_qg

contains









subroutine amp_VH_qg(Mom,mass,helicity,id,amp)
  implicit none
  real(8), intent(in) :: Mom(1:4,1:10)
  real(8), intent(in) :: mass(3:5,1:2)
  real(8), intent(in) :: helicity(10)
  complex(8), intent(out) :: amp
  integer  id(10)

  amp=0d0

  return
END subroutine amp_VH_qg





end module ModVHqg
!!--YaofuZhou-----------------------------------------