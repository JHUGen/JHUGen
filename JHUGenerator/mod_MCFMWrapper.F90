MODULE ModMCFMWrapper
implicit none

contains

subroutine MCFM_firsttime(collider_energy)
implicit none
real(8), intent(in) :: collider_energy

! Cannot just pass these as arguments. Need to specify the char lengths
! Otherwise MCFM crashes.
character*72 :: inputfile,workdir

! MCFM declarations
double precision sqrts
common/energy/sqrts

inputfile='input.DAT'
workdir='./'

call mcfm_init(inputfile,workdir)
call qlinit()

! Now let's have some fun
sqrts = dble(collider_energy)
print *,sqrts

end subroutine

END MODULE
