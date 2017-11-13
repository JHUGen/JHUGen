!--YaofuZhou-----------------------------------------
module ModVHvirtual
  use ModParameters
  use ModMisc
  use ModVHLO
  use Collier
  implicit none
  
public :: amp_VH_virtual

contains

subroutine amp_VH_virtual(MomExt,mass,helicity,id,amp)

  implicit none
  real(8), intent(in) :: MomExt(1:4,1:9)
  real(8), intent(in) :: mass(3:5,1:2)
  real(8), intent(in) :: helicity(9)
  complex(8), intent(out) :: amp
  integer  id(9)

  amp = amp_VH_LO(MomExt,mass,helicity,id)

  return
END subroutine









end module ModVHvirtual
!!--YaofuZhou-----------------------------------------
