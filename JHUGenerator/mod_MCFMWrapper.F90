MODULE ModMCFMWrapper
implicit none

contains

subroutine MCFM_firsttime()
implicit none
! Cannot just pass these as arguments. Need to specify the char lengths
! Otherwise MCFM crashes.
character*72 :: inputfile,workdir
inputfile='input.DAT'
workdir='./'

call mcfm_init(inputfile,workdir)
call qlinit()

end subroutine

! replace variables defined through input.DAT
! equivalent to the following MELA functions
! - TEvtProb::InitializeMCFM
! - TEvtProb::CrossInitialize
! - TUtil::SetEwkCouplingParameters
! - TUtil::SetMass
! - TUtil::SetDecayWidth
! - TUtil::SetCKMElements
! - TUtil::SetMCFMSpinZeroCouplings
! For each common block, have one function
! to re-group MELA calls into the different common blocks.
! Somehow re-defining the common blocks is allowed, not sure why/how.
subroutine Init_MCFMCommon_energy(collider_energy)
implicit none
real(8), intent(in) :: collider_energy
! MCFM declarations
double precision sqrts
common/energy/sqrts

sqrts = dble(collider_energy)
end subroutine




! Record ModParameters variables into arrays
! to be used by similar Set* functions to set MCFM variables.
! Doing it this way allows us to pass a zillion arguments (e.g. vvcoupl length)
subroutine GetSpinZeroVVCouplings(vvcoupl, cqsq, Lambda_qsq, useWWcoupl)
   use ModParameters
   implicit none
   complex(8), intent(out) :: vvcoupl(39)
   integer, intent(out) :: cqsq(3)
   double precision, intent(out) :: Lambda_qsq(1:3,1:4)
   logical, intent(in) :: useWWcoupl

   vvcoupl(:)=0d0
   cqsq(:)=0d0
   Lambda_qsq(:,:)=0d0

   if(.not.useWWcoupl) then
      cqsq(1) = cz_q1sq
      Lambda_qsq(1,1) = Lambda_z11
      Lambda_qsq(1,2) = Lambda_z21
      Lambda_qsq(1,3) = Lambda_z31
      Lambda_qsq(1,4) = Lambda_z41
      cqsq(2) = cz_q2sq
      Lambda_qsq(2,1) = Lambda_z12
      Lambda_qsq(2,2) = Lambda_z22
      Lambda_qsq(2,3) = Lambda_z32
      Lambda_qsq(2,4) = Lambda_z42
      cqsq(3) = cz_q12sq
      Lambda_qsq(3,1) = Lambda_z10
      Lambda_qsq(3,2) = Lambda_z20
      Lambda_qsq(3,3) = Lambda_z30
      Lambda_qsq(3,4) = Lambda_z40

      vvcoupl(1) = ghz1
      vvcoupl(2) = ghz2
      vvcoupl(3) = ghz3
      vvcoupl(4) = ghz4

      vvcoupl(5) = ghzgs2
      vvcoupl(6) = ghzgs3
      vvcoupl(7) = ghzgs4
      vvcoupl(8) = ghgsgs2
      vvcoupl(9) = ghgsgs3
      vvcoupl(10) = ghgsgs4

      vvcoupl(11) = ghz1_prime
      vvcoupl(12) = ghz1_prime2
      vvcoupl(13) = ghz1_prime3
      vvcoupl(14) = ghz1_prime4
      vvcoupl(15) = ghz1_prime5

      vvcoupl(16) = ghz2_prime
      vvcoupl(17) = ghz2_prime2
      vvcoupl(18) = ghz2_prime3
      vvcoupl(19) = ghz2_prime4
      vvcoupl(20) = ghz2_prime5

      vvcoupl(21) = ghz3_prime
      vvcoupl(22) = ghz3_prime2
      vvcoupl(23) = ghz3_prime3
      vvcoupl(24) = ghz3_prime4
      vvcoupl(25) = ghz3_prime5

      vvcoupl(26) = ghz4_prime
      vvcoupl(27) = ghz4_prime2
      vvcoupl(28) = ghz4_prime3
      vvcoupl(29) = ghz4_prime4
      vvcoupl(30) = ghz4_prime5

      vvcoupl(31) = ghzgs1_prime2

      vvcoupl(32) = ghz1_prime6
      vvcoupl(33) = ghz1_prime7
      vvcoupl(34) = ghz2_prime6
      vvcoupl(35) = ghz2_prime7
      vvcoupl(36) = ghz3_prime6
      vvcoupl(37) = ghz3_prime7
      vvcoupl(38) = ghz4_prime6
      vvcoupl(39) = ghz4_prime7

   else
      cqsq(1) = cw_q1sq
      Lambda_qsq(1,1) = Lambda_w11
      Lambda_qsq(1,2) = Lambda_w21
      Lambda_qsq(1,3) = Lambda_w31
      Lambda_qsq(1,4) = Lambda_w41
      cqsq(2) = cw_q2sq
      Lambda_qsq(2,1) = Lambda_w12
      Lambda_qsq(2,2) = Lambda_w22
      Lambda_qsq(2,3) = Lambda_w32
      Lambda_qsq(2,4) = Lambda_w42
      cqsq(3) = cw_q12sq
      Lambda_qsq(3,1) = Lambda_w10
      Lambda_qsq(3,2) = Lambda_w20
      Lambda_qsq(3,3) = Lambda_w30
      Lambda_qsq(3,4) = Lambda_w40

      vvcoupl(1) = ghw1
      vvcoupl(2) = ghw2
      vvcoupl(3) = ghw3
      vvcoupl(4) = ghw4

      vvcoupl(11) = ghw1_prime
      vvcoupl(12) = ghw1_prime2
      vvcoupl(13) = ghw1_prime3
      vvcoupl(14) = ghw1_prime4
      vvcoupl(15) = ghw1_prime5

      vvcoupl(16) = ghw2_prime
      vvcoupl(17) = ghw2_prime2
      vvcoupl(18) = ghw2_prime3
      vvcoupl(19) = ghw2_prime4
      vvcoupl(20) = ghw2_prime5

      vvcoupl(21) = ghw3_prime
      vvcoupl(22) = ghw3_prime2
      vvcoupl(23) = ghw3_prime3
      vvcoupl(24) = ghw3_prime4
      vvcoupl(25) = ghw3_prime5

      vvcoupl(26) = ghw4_prime
      vvcoupl(27) = ghw4_prime2
      vvcoupl(28) = ghw4_prime3
      vvcoupl(29) = ghw4_prime4
      vvcoupl(30) = ghw4_prime5

      vvcoupl(32) = ghw1_prime6
      vvcoupl(33) = ghw1_prime7
      vvcoupl(34) = ghw2_prime6
      vvcoupl(35) = ghw2_prime7
      vvcoupl(36) = ghw3_prime6
      vvcoupl(37) = ghw3_prime7
      vvcoupl(38) = ghw4_prime6
      vvcoupl(39) = ghw4_prime7

      cw_q1sq = cqsq(1)
      Lambda_w11 = Lambda_qsq(1,1)
      Lambda_w21 = Lambda_qsq(1,2)
      Lambda_w31 = Lambda_qsq(1,3)
      Lambda_w41 = Lambda_qsq(1,4)
      cw_q2sq = cqsq(2)
      Lambda_w12 = Lambda_qsq(2,1)
      Lambda_w22 = Lambda_qsq(2,2)
      Lambda_w32 = Lambda_qsq(2,3)
      Lambda_w42 = Lambda_qsq(2,4)
      cw_q12sq = cqsq(3)
      Lambda_w10 = Lambda_qsq(3,1)
      Lambda_w20 = Lambda_qsq(3,2)
      Lambda_w30 = Lambda_qsq(3,3)
      Lambda_w40 = Lambda_qsq(3,4)

      ghw1 =  vvcoupl(1)
      ghw2 =  vvcoupl(2)
      ghw3 =  vvcoupl(3)
      ghw4 =  vvcoupl(4)

      ghw1_prime = vvcoupl(11)
      ghw1_prime2= vvcoupl(12)
      ghw1_prime3= vvcoupl(13)
      ghw1_prime4= vvcoupl(14)
      ghw1_prime5= vvcoupl(15)

      ghw2_prime = vvcoupl(16)
      ghw2_prime2= vvcoupl(17)
      ghw2_prime3= vvcoupl(18)
      ghw2_prime4= vvcoupl(19)
      ghw2_prime5= vvcoupl(20)

      ghw3_prime = vvcoupl(21)
      ghw3_prime2= vvcoupl(22)
      ghw3_prime3= vvcoupl(23)
      ghw3_prime4= vvcoupl(24)
      ghw3_prime5= vvcoupl(25)

      ghw4_prime = vvcoupl(26)
      ghw4_prime2= vvcoupl(27)
      ghw4_prime3= vvcoupl(28)
      ghw4_prime4= vvcoupl(29)
      ghw4_prime5= vvcoupl(30)

      ghw1_prime6  = vvcoupl(32)
      ghw1_prime7  = vvcoupl(33)
      ghw2_prime6  = vvcoupl(34)
      ghw2_prime7  = vvcoupl(35)
      ghw3_prime6  = vvcoupl(36)
      ghw3_prime7  = vvcoupl(37)
      ghw4_prime6  = vvcoupl(38)
      ghw4_prime7  = vvcoupl(39)
   endif
   return
end subroutine

subroutine GetDistinguishWWCouplingsFlag(doAllow)
   use ModParameters
   implicit none
   integer, intent(out) :: doAllow
   if(distinguish_HWWcouplings) then
      doAllow=1
   else
      doAllow=0
   endif
   return
end subroutine

subroutine GetSpinZeroGGCouplings(ggcoupl)
   use ModParameters
   implicit none
   double complex, intent(out) :: ggcoupl(1:3)
   ggcoupl(1) = dcmplx(ghg2)
   ggcoupl(2) = dcmplx(ghg3)
   ggcoupl(3) = dcmplx(ghg4)
   return
end subroutine

subroutine GetSpinZeroQQCouplings(qqcoupl)
   use ModParameters
   implicit none
   double complex, intent(out) :: qqcoupl(1:2)
   qqcoupl(1) = dcmplx(kappa)
   qqcoupl(2) = dcmplx(kappa_tilde)
   return
end subroutine


END MODULE
