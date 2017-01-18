MODULE ModMCFMWrapper
implicit none

integer, parameter ::gHIGGS_KAPPA=1,gHIGGS_KAPPA_TILDE=2,SIZE_HQQ=2
integer, parameter ::gHIGGS_GG_2=1,gHIGGS_GG_3=2,gHIGGS_GG_4=3,SIZE_HGG=3
integer, parameter ::gHIGGS_VV_1=1,gHIGGS_VV_2=2,gHIGGS_VV_3=3,gHIGGS_VV_4=4,gHIGGS_ZA_2=5,gHIGGS_ZA_3=6,gHIGGS_ZA_4=7,gHIGGS_AA_2=8,gHIGGS_AA_3=9,gHIGGS_AA_4=10,gHIGGS_VV_1_PRIME=11,gHIGGS_VV_1_PRIME2=12,gHIGGS_VV_1_PRIME3=13,gHIGGS_VV_1_PRIME4=14,gHIGGS_VV_1_PRIME5=15,gHIGGS_VV_2_PRIME=16,gHIGGS_VV_2_PRIME2=17,gHIGGS_VV_2_PRIME3=18,gHIGGS_VV_2_PRIME4=19,gHIGGS_VV_2_PRIME5=20,gHIGGS_VV_3_PRIME=21,gHIGGS_VV_3_PRIME2=22,gHIGGS_VV_3_PRIME3=23,gHIGGS_VV_3_PRIME4=24,gHIGGS_VV_3_PRIME5=25,gHIGGS_VV_4_PRIME=26,gHIGGS_VV_4_PRIME2=27,gHIGGS_VV_4_PRIME3=28,gHIGGS_VV_4_PRIME4=29,gHIGGS_VV_4_PRIME5=30,gHIGGS_ZA_1_PRIME2=31,gHIGGS_VV_1_PRIME6=32,gHIGGS_VV_1_PRIME7=33,gHIGGS_VV_2_PRIME6=34,gHIGGS_VV_2_PRIME7=35,gHIGGS_VV_3_PRIME6=36,gHIGGS_VV_3_PRIME7=37,gHIGGS_VV_4_PRIME6=38,gHIGGS_VV_4_PRIME7=39,SIZE_HVV=39
integer, parameter ::LambdaHIGGS_QSQ_VV_1=1,LambdaHIGGS_QSQ_VV_2=2,LambdaHIGGS_QSQ_VV_3=3,LambdaHIGGS_QSQ_VV_4=4,SIZE_HVV_LAMBDAQSQ=4
integer, parameter ::cLambdaHIGGS_VV_QSQ1=1,cLambdaHIGGS_VV_QSQ2=2,cLambdaHIGGS_VV_QSQ12=3,SIZE_HVV_CQSQ=3
integer, parameter ::gZPRIME_QQ_LEFT=1,gZPRIME_QQ_RIGHT=2,SIZE_ZQQ=2
integer, parameter :: gZPRIME_VV_1=1,gZPRIME_VV_2=2,SIZE_ZVV=2
integer, parameter :: gGRAVITON_QQ_LEFT=1,gGRAVITON_QQ_RIGHT=2,SIZE_GQQ=2
integer, parameter :: gGRAVITON_GG_1=1,gGRAVITON_GG_2=2,gGRAVITON_GG_3=3,gGRAVITON_GG_4=4,gGRAVITON_GG_5=5,SIZE_GGG=5
integer, parameter :: gGRAVITON_VV_1=1,gGRAVITON_VV_2=2,gGRAVITON_VV_3=3,gGRAVITON_VV_4=4,gGRAVITON_VV_5=5,gGRAVITON_VV_6=6,gGRAVITON_VV_7=7,gGRAVITON_VV_8=8,gGRAVITON_VV_9=9,gGRAVITON_VV_10=10,SIZE_GVV=10

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

subroutine Init_MCFMCommon_spinzerohiggs_anomcoupl( &
  Hggcoupl, &
  Hqqcoupl, &
  Httcoupl, &
  Hbbcoupl, &
  Hg4g4coupl, &
  Ht4t4coupl, &
  Hb4b4coupl, &

  Hzzcoupl, &
  Hwwcoupl, &

  HLambdaBSM, &
  HLambda_Q, &
  HLambda_zgs1, &
  HzzLambda, &
  HwwLambda, &

  HzzLambda_qsq, &
  HwwLambda_qsq, &
  HzzCLambda_qsq, &
  HwwCLambda_qsq, &

  Hb4b4_mb_4gen, &
  Ht4t4_mt_4gen, &

  separateWWZZcouplings &
)
  double complex Hggcoupl(1:SIZE_HGG)
  double complex Hqqcoupl(1:SIZE_HQQ)
  double complex Httcoupl(1:SIZE_HQQ)
  double complex Hbbcoupl(1:SIZE_HQQ)
  double complex Hg4g4coupl(1:SIZE_HGG)
  double complex Ht4t4coupl(1:SIZE_HQQ)
  double complex Hb4b4coupl(1:SIZE_HQQ)

  double complex Hzzcoupl(1:SIZE_HVV)
  double complex Hwwcoupl(1:SIZE_HVV)

  double precision HLambdaBSM
  double precision HLambda_Q
  double precision HLambda_zgs1
  double precision HzzLambda(1:SIZE_HVV_LAMBDAQSQ)
  double precision HwwLambda(1:SIZE_HVV_LAMBDAQSQ)

  double precision HzzLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)
  double precision HwwLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)

  double precision Hb4b4_mb_4gen
  double precision Ht4t4_mt_4gen

  integer HzzCLambda_qsq(1:SIZE_HVV_CQSQ)
  integer HwwCLambda_qsq(1:SIZE_HVV_CQSQ)

  integer separateWWZZcouplings

! MCFM declarations
   integer AllowAnomalousCouplings
   integer distinguish_HWWcouplings
   integer AnomalCouplPR
   integer AnomalCouplDK
   integer channeltoggle_stu ! 0, 1, 2 for s, t+u and s+t+u
   integer vvhvvtoggle_vbfvh ! 0, 1, 2 for VBF, VH and VB+VH

   integer cz_q1sq,cz_q2sq,cz_q12sq
   integer cw_q1sq,cw_q2sq,cw_q12sq
   integer c2z_q1sq,c2z_q2sq,c2z_q12sq
   integer c2w_q1sq,c2w_q2sq,c2w_q12sq

   double precision mb_4gen,mt_4gen

   double precision LambdaBSM,Lambda_Q
   double precision Lambda_zgs1
   double precision Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4
   double precision Lambda_z11,Lambda_z21,Lambda_z31,Lambda_z41
   double precision Lambda_z12,Lambda_z22,Lambda_z32,Lambda_z42
   double precision Lambda_z10,Lambda_z20,Lambda_z30,Lambda_z40
   double precision Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4
   double precision Lambda_w11,Lambda_w21,Lambda_w31,Lambda_w41
   double precision Lambda_w12,Lambda_w22,Lambda_w32,Lambda_w42
   double precision Lambda_w10,Lambda_w20,Lambda_w30,Lambda_w40

   double precision h2mass,h2width

   double precision Lambda2BSM,Lambda2_Q
   double precision Lambda2_zgs1
   double precision Lambda2_z1,Lambda2_z2,Lambda2_z3,Lambda2_z4
   double precision Lambda2_z11,Lambda2_z21,Lambda2_z31,Lambda2_z41
   double precision Lambda2_z12,Lambda2_z22,Lambda2_z32,Lambda2_z42
   double precision Lambda2_z10,Lambda2_z20,Lambda2_z30,Lambda2_z40
   double precision Lambda2_w1,Lambda2_w2,Lambda2_w3,Lambda2_w4
   double precision Lambda2_w11,Lambda2_w21,Lambda2_w31,Lambda2_w41
   double precision Lambda2_w12,Lambda2_w22,Lambda2_w32,Lambda2_w42
   double precision Lambda2_w10,Lambda2_w20,Lambda2_w30,Lambda2_w40

   double complex kappa_top,kappa_tilde_top
   double complex kappa_bot,kappa_tilde_bot
   double complex ghg2,ghg3,ghg4
   double complex kappa_4gen_top,kappa_tilde_4gen_top
   double complex kappa_4gen_bot,kappa_tilde_4gen_bot
   double complex ghg2_4gen,ghg3_4gen,ghg4_4gen

   double complex ghz1,ghz2,ghz3,ghz4
   double complex ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime
   double complex ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2
   double complex ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3
   double complex ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4
   double complex ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5
   double complex ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6
   double complex ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7

   double complex ghzgs1_prime2,ghzgs2,ghzgs3,ghzgs4
   double complex ghgsgs2,ghgsgs3,ghgsgs4

   double complex ghw1,ghw2,ghw3,ghw4
   double complex ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime
   double complex ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2
   double complex ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3
   double complex ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4
   double complex ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5
   double complex ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6
   double complex ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7


   double complex kappa2_top,kappa2_tilde_top
   double complex kappa2_bot,kappa2_tilde_bot
   double complex gh2g2,gh2g3,gh2g4
   double complex kappa2_4gen_top,kappa2_tilde_4gen_top
   double complex kappa2_4gen_bot,kappa2_tilde_4gen_bot
   double complex gh2g2_4gen,gh2g3_4gen,gh2g4_4gen

   double complex gh2z1,gh2z2,gh2z3,gh2z4
   double complex gh2z1_prime,gh2z2_prime,gh2z3_prime,gh2z4_prime
   double complex gh2z1_prime2,gh2z2_prime2,gh2z3_prime2,gh2z4_prime2
   double complex gh2z1_prime3,gh2z2_prime3,gh2z3_prime3,gh2z4_prime3
   double complex gh2z1_prime4,gh2z2_prime4,gh2z3_prime4,gh2z4_prime4
   double complex gh2z1_prime5,gh2z2_prime5,gh2z3_prime5,gh2z4_prime5
   double complex gh2z1_prime6,gh2z2_prime6,gh2z3_prime6,gh2z4_prime6
   double complex gh2z1_prime7,gh2z2_prime7,gh2z3_prime7,gh2z4_prime7

   double complex gh2zgs1_prime2,gh2zgs2,gh2zgs3,gh2zgs4
   double complex gh2gsgs2,gh2gsgs3,gh2gsgs4

   double complex gh2w1,gh2w2,gh2w3,gh2w4
   double complex gh2w1_prime,gh2w2_prime,gh2w3_prime,gh2w4_prime
   double complex gh2w1_prime2,gh2w2_prime2,gh2w3_prime2,gh2w4_prime2
   double complex gh2w1_prime3,gh2w2_prime3,gh2w3_prime3,gh2w4_prime3
   double complex gh2w1_prime4,gh2w2_prime4,gh2w3_prime4,gh2w4_prime4
   double complex gh2w1_prime5,gh2w2_prime5,gh2w3_prime5,gh2w4_prime5
   double complex gh2w1_prime6,gh2w2_prime6,gh2w3_prime6,gh2w4_prime6
   double complex gh2w1_prime7,gh2w2_prime7,gh2w3_prime7,gh2w4_prime7


   common/spinzerohiggs_anomcoupl/ &
   AllowAnomalousCouplings, &
   distinguish_HWWcouplings, &
   AnomalCouplPR,AnomalCouplDK, &
   channeltoggle_stu,vvhvvtoggle_vbfvh, &
   cz_q1sq,cz_q2sq,cz_q12sq, &
   cw_q1sq,cw_q2sq,cw_q12sq, &
   c2z_q1sq,c2z_q2sq,c2z_q12sq, &
   c2w_q1sq,c2w_q2sq,c2w_q12sq, &

   mb_4gen,mt_4gen, &

   LambdaBSM,Lambda_Q, &
   Lambda_zgs1, &
   Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4, &
   Lambda_z11,Lambda_z21,Lambda_z31,Lambda_z41, &
   Lambda_z12,Lambda_z22,Lambda_z32,Lambda_z42, &
   Lambda_z10,Lambda_z20,Lambda_z30,Lambda_z40, &
   Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4, &
   Lambda_w11,Lambda_w21,Lambda_w31,Lambda_w41, &
   Lambda_w12,Lambda_w22,Lambda_w32,Lambda_w42, &
   Lambda_w10,Lambda_w20,Lambda_w30,Lambda_w40, &

   h2mass,h2width, &

   Lambda2BSM,Lambda2_Q, &
   Lambda2_zgs1, &
   Lambda2_z1,Lambda2_z2,Lambda2_z3,Lambda2_z4, &
   Lambda2_z11,Lambda2_z21,Lambda2_z31,Lambda2_z41, &
   Lambda2_z12,Lambda2_z22,Lambda2_z32,Lambda2_z42, &
   Lambda2_z10,Lambda2_z20,Lambda2_z30,Lambda2_z40, &
   Lambda2_w1,Lambda2_w2,Lambda2_w3,Lambda2_w4, &
   Lambda2_w11,Lambda2_w21,Lambda2_w31,Lambda2_w41, &
   Lambda2_w12,Lambda2_w22,Lambda2_w32,Lambda2_w42, &
   Lambda2_w10,Lambda2_w20,Lambda2_w30,Lambda2_w40, &


   kappa_top,kappa_tilde_top, &
   kappa_bot,kappa_tilde_bot, &
   ghg2,ghg3,ghg4, &
   kappa_4gen_top,kappa_tilde_4gen_top, &
   kappa_4gen_bot,kappa_tilde_4gen_bot, &
   ghg2_4gen,ghg3_4gen,ghg4_4gen, &

   ghz1,ghz2,ghz3,ghz4, &
   ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime, &
   ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2, &
   ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3, &
   ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4, &
   ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5, &
   ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6, &
   ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7, &

   ghzgs1_prime2,ghzgs2,ghzgs3,ghzgs4, &
   ghgsgs2,ghgsgs3,ghgsgs4, &

   ghw1,ghw2,ghw3,ghw4, &
   ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime, &
   ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2, &
   ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3, &
   ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4, &
   ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5, &
   ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6, &
   ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7, &


   kappa2_top,kappa2_tilde_top, &
   kappa2_bot,kappa2_tilde_bot, &
   gh2g2,gh2g3,gh2g4, &
   kappa2_4gen_top,kappa2_tilde_4gen_top, &
   kappa2_4gen_bot,kappa2_tilde_4gen_bot, &
   gh2g2_4gen,gh2g3_4gen,gh2g4_4gen, &

   gh2z1,gh2z2,gh2z3,gh2z4, &
   gh2z1_prime,gh2z2_prime,gh2z3_prime,gh2z4_prime, &
   gh2z1_prime2,gh2z2_prime2,gh2z3_prime2,gh2z4_prime2, &
   gh2z1_prime3,gh2z2_prime3,gh2z3_prime3,gh2z4_prime3, &
   gh2z1_prime4,gh2z2_prime4,gh2z3_prime4,gh2z4_prime4, &
   gh2z1_prime5,gh2z2_prime5,gh2z3_prime5,gh2z4_prime5, &
   gh2z1_prime6,gh2z2_prime6,gh2z3_prime6,gh2z4_prime6, &
   gh2z1_prime7,gh2z2_prime7,gh2z3_prime7,gh2z4_prime7, &

   gh2zgs1_prime2,gh2zgs2,gh2zgs3,gh2zgs4, &
   gh2gsgs2,gh2gsgs3,gh2gsgs4, &

   gh2w1,gh2w2,gh2w3,gh2w4, &
   gh2w1_prime,gh2w2_prime,gh2w3_prime,gh2w4_prime, &
   gh2w1_prime2,gh2w2_prime2,gh2w3_prime2,gh2w4_prime2, &
   gh2w1_prime3,gh2w2_prime3,gh2w3_prime3,gh2w4_prime3, &
   gh2w1_prime4,gh2w2_prime4,gh2w3_prime4,gh2w4_prime4, &
   gh2w1_prime5,gh2w2_prime5,gh2w3_prime5,gh2w4_prime5, &
   gh2w1_prime6,gh2w2_prime6,gh2w3_prime6,gh2w4_prime6, &
   gh2w1_prime7,gh2w2_prime7,gh2w3_prime7,gh2w4_prime7
   ! End common/spinzerohiggs_anomcoupl

   ! Begin assignments
   AllowAnomalousCouplings = 1
   distinguish_HWWcouplings = separateWWZZcouplings
   AnomalCouplPR=1
   AnomalCouplDK=1
   channeltoggle_stu=2 ! s=0/t+u=1/s+t+u=2
   vvhvvtoggle_vbfvh=2 ! vbf=0/vh=1/vbf+vh=2

!   /***** REGULAR RESONANCE *****/
   LambdaBSM = HLambdaBSM
   Lambda_Q = HLambda_Q

   Lambda_zgs1 = HLambda_zgs1
   Lambda_z1 = HzzLambda(LambdaHIGGS_QSQ_VV_1)
   Lambda_z2 = HzzLambda(LambdaHIGGS_QSQ_VV_2)
   Lambda_z3 = HzzLambda(LambdaHIGGS_QSQ_VV_3)
   Lambda_z4 = HzzLambda(LambdaHIGGS_QSQ_VV_4)

   cz_q1sq = HzzCLambda_qsq(cLambdaHIGGS_VV_QSQ1)
   Lambda_z11 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1)
   Lambda_z21 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1)
   Lambda_z31 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1)
   Lambda_z41 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1)
   cz_q2sq = HzzCLambda_qsq(cLambdaHIGGS_VV_QSQ2)
   Lambda_z12 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2)
   Lambda_z22 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2)
   Lambda_z32 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2)
   Lambda_z42 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2)
   cz_q12sq = HzzCLambda_qsq(cLambdaHIGGS_VV_QSQ12)
   Lambda_z10 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12)
   Lambda_z20 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12)
   Lambda_z30 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12)
   Lambda_z40 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12)

   kappa_top = Httcoupl(gHIGGS_KAPPA)
   kappa_bot = Hbbcoupl(gHIGGS_KAPPA)
   kappa_tilde_top = Httcoupl(gHIGGS_KAPPA_TILDE)
   kappa_tilde_bot = Hbbcoupl(gHIGGS_KAPPA_TILDE)
   ghg2 = Hggcoupl(gHIGGS_GG_2)
   ghg3 = Hggcoupl(gHIGGS_GG_3)
   ghg4 = Hggcoupl(gHIGGS_GG_4)
   mb_4gen = Hb4b4_mb_4gen
   mt_4gen = Ht4t4_mt_4gen
   kappa_4gen_top = Ht4t4coupl(gHIGGS_KAPPA)
   kappa_4gen_bot = Hb4b4coupl(gHIGGS_KAPPA)
   kappa_tilde_4gen_top = Ht4t4coupl(gHIGGS_KAPPA_TILDE)
   kappa_tilde_4gen_bot = Hb4b4coupl(gHIGGS_KAPPA_TILDE)
   ghg2_4gen = Hg4g4coupl(gHIGGS_GG_2)
   ghg3_4gen = Hg4g4coupl(gHIGGS_GG_3)
   ghg4_4gen = Hg4g4coupl(gHIGGS_GG_4)

   ghz1 = Hzzcoupl(gHIGGS_VV_1)
   ghz2 = Hzzcoupl(gHIGGS_VV_2)
   ghz3 = Hzzcoupl(gHIGGS_VV_3)
   ghz4 = Hzzcoupl(gHIGGS_VV_4)
   ghz1_prime = Hzzcoupl(gHIGGS_VV_1_PRIME)
   ghz1_prime2 = Hzzcoupl(gHIGGS_VV_1_PRIME2)
   ghz1_prime3 = Hzzcoupl(gHIGGS_VV_1_PRIME3)
   ghz1_prime4 = Hzzcoupl(gHIGGS_VV_1_PRIME4)
   ghz1_prime5 = Hzzcoupl(gHIGGS_VV_1_PRIME5)
   ghz1_prime6 = Hzzcoupl(gHIGGS_VV_1_PRIME6)
   ghz1_prime7 = Hzzcoupl(gHIGGS_VV_1_PRIME7)
   ghz2_prime = Hzzcoupl(gHIGGS_VV_2_PRIME)
   ghz2_prime2 = Hzzcoupl(gHIGGS_VV_2_PRIME2)
   ghz2_prime3 = Hzzcoupl(gHIGGS_VV_2_PRIME3)
   ghz2_prime4 = Hzzcoupl(gHIGGS_VV_2_PRIME4)
   ghz2_prime5 = Hzzcoupl(gHIGGS_VV_2_PRIME5)
   ghz2_prime6 = Hzzcoupl(gHIGGS_VV_2_PRIME6)
   ghz2_prime7 = Hzzcoupl(gHIGGS_VV_2_PRIME7)
   ghz3_prime = Hzzcoupl(gHIGGS_VV_3_PRIME)
   ghz3_prime2 = Hzzcoupl(gHIGGS_VV_3_PRIME2)
   ghz3_prime3 = Hzzcoupl(gHIGGS_VV_3_PRIME3)
   ghz3_prime4 = Hzzcoupl(gHIGGS_VV_3_PRIME4)
   ghz3_prime5 = Hzzcoupl(gHIGGS_VV_3_PRIME5)
   ghz3_prime6 = Hzzcoupl(gHIGGS_VV_3_PRIME6)
   ghz3_prime7 = Hzzcoupl(gHIGGS_VV_3_PRIME7)
   ghz4_prime = Hzzcoupl(gHIGGS_VV_4_PRIME)
   ghz4_prime2 = Hzzcoupl(gHIGGS_VV_4_PRIME2)
   ghz4_prime3 = Hzzcoupl(gHIGGS_VV_4_PRIME3)
   ghz4_prime4 = Hzzcoupl(gHIGGS_VV_4_PRIME4)
   ghz4_prime5 = Hzzcoupl(gHIGGS_VV_4_PRIME5)
   ghz4_prime6 = Hzzcoupl(gHIGGS_VV_4_PRIME6)
   ghz4_prime7 = Hzzcoupl(gHIGGS_VV_4_PRIME7)

   ghzgs1_prime2 = Hzzcoupl(gHIGGS_ZA_1_PRIME2)
   ghzgs2 = Hzzcoupl(gHIGGS_ZA_2)
   ghzgs3 = Hzzcoupl(gHIGGS_ZA_3)
   ghzgs4 = Hzzcoupl(gHIGGS_ZA_4)
   ghgsgs2 = Hzzcoupl(gHIGGS_AA_2)
   ghgsgs3 = Hzzcoupl(gHIGGS_AA_3)
   ghgsgs4 = Hzzcoupl(gHIGGS_AA_4)

   if (distinguish_HWWcouplings .eq. 1) then
      Lambda_w1 = HwwLambda(LambdaHIGGS_QSQ_VV_1)
      Lambda_w2 = HwwLambda(LambdaHIGGS_QSQ_VV_2)
      Lambda_w3 = HwwLambda(LambdaHIGGS_QSQ_VV_3)
      Lambda_w4 = HwwLambda(LambdaHIGGS_QSQ_VV_4)

      cw_q1sq = HwwCLambda_qsq(cLambdaHIGGS_VV_QSQ1)
      Lambda_w11 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1)
      Lambda_w21 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1)
      Lambda_w31 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1)
      Lambda_w41 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1)
      cw_q2sq = HwwCLambda_qsq(cLambdaHIGGS_VV_QSQ2)
      Lambda_w12 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2)
      Lambda_w22 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2)
      Lambda_w32 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2)
      Lambda_w42 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2)
      cw_q12sq = HwwCLambda_qsq(cLambdaHIGGS_VV_QSQ12)
      Lambda_w10 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12)
      Lambda_w20 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12)
      Lambda_w30 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12)
      Lambda_w40 = HwwLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12)

      ghw1 = Hwwcoupl(gHIGGS_VV_1)
      ghw2 = Hwwcoupl(gHIGGS_VV_2)
      ghw3 = Hwwcoupl(gHIGGS_VV_3)
      ghw4 = Hwwcoupl(gHIGGS_VV_4)
      ghw1_prime = Hwwcoupl(gHIGGS_VV_1_PRIME)
      ghw1_prime2 = Hwwcoupl(gHIGGS_VV_1_PRIME2)
      ghw1_prime3 = Hwwcoupl(gHIGGS_VV_1_PRIME3)
      ghw1_prime4 = Hwwcoupl(gHIGGS_VV_1_PRIME4)
      ghw1_prime5 = Hwwcoupl(gHIGGS_VV_1_PRIME5)
      ghw1_prime6 = Hwwcoupl(gHIGGS_VV_1_PRIME6)
      ghw1_prime7 = Hwwcoupl(gHIGGS_VV_1_PRIME7)
      ghw2_prime = Hwwcoupl(gHIGGS_VV_2_PRIME)
      ghw2_prime2 = Hwwcoupl(gHIGGS_VV_2_PRIME2)
      ghw2_prime3 = Hwwcoupl(gHIGGS_VV_2_PRIME3)
      ghw2_prime4 = Hwwcoupl(gHIGGS_VV_2_PRIME4)
      ghw2_prime5 = Hwwcoupl(gHIGGS_VV_2_PRIME5)
      ghw2_prime6 = Hwwcoupl(gHIGGS_VV_2_PRIME6)
      ghw2_prime7 = Hwwcoupl(gHIGGS_VV_2_PRIME7)
      ghw3_prime = Hwwcoupl(gHIGGS_VV_3_PRIME)
      ghw3_prime2 = Hwwcoupl(gHIGGS_VV_3_PRIME2)
      ghw3_prime3 = Hwwcoupl(gHIGGS_VV_3_PRIME3)
      ghw3_prime4 = Hwwcoupl(gHIGGS_VV_3_PRIME4)
      ghw3_prime5 = Hwwcoupl(gHIGGS_VV_3_PRIME5)
      ghw3_prime6 = Hwwcoupl(gHIGGS_VV_3_PRIME6)
      ghw3_prime7 = Hwwcoupl(gHIGGS_VV_3_PRIME7)
      ghw4_prime = Hwwcoupl(gHIGGS_VV_4_PRIME)
      ghw4_prime2 = Hwwcoupl(gHIGGS_VV_4_PRIME2)
      ghw4_prime3 = Hwwcoupl(gHIGGS_VV_4_PRIME3)
      ghw4_prime4 = Hwwcoupl(gHIGGS_VV_4_PRIME4)
      ghw4_prime5 = Hwwcoupl(gHIGGS_VV_4_PRIME5)
      ghw4_prime6 = Hwwcoupl(gHIGGS_VV_4_PRIME6)
      ghw4_prime7 = Hwwcoupl(gHIGGS_VV_4_PRIME7)

   else
      Lambda_w1 = HzzLambda(LambdaHIGGS_QSQ_VV_1)
      Lambda_w2 = HzzLambda(LambdaHIGGS_QSQ_VV_2)
      Lambda_w3 = HzzLambda(LambdaHIGGS_QSQ_VV_3)
      Lambda_w4 = HzzLambda(LambdaHIGGS_QSQ_VV_4)

      cw_q1sq = HzzCLambda_qsq(cLambdaHIGGS_VV_QSQ1)
      Lambda_w11 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1)
      Lambda_w21 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1)
      Lambda_w31 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1)
      Lambda_w41 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1)
      cw_q2sq = HzzCLambda_qsq(cLambdaHIGGS_VV_QSQ2)
      Lambda_w12 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2)
      Lambda_w22 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2)
      Lambda_w32 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2)
      Lambda_w42 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2)
      cw_q12sq = HzzCLambda_qsq(cLambdaHIGGS_VV_QSQ12)
      Lambda_w10 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12)
      Lambda_w20 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12)
      Lambda_w30 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12)
      Lambda_w40 = HzzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12)

      ghw1 = Hzzcoupl(gHIGGS_VV_1)
      ghw2 = Hzzcoupl(gHIGGS_VV_2)
      ghw3 = Hzzcoupl(gHIGGS_VV_3)
      ghw4 = Hzzcoupl(gHIGGS_VV_4)
      ghw1_prime = Hzzcoupl(gHIGGS_VV_1_PRIME)
      ghw1_prime2 = Hzzcoupl(gHIGGS_VV_1_PRIME2)
      ghw1_prime3 = Hzzcoupl(gHIGGS_VV_1_PRIME3)
      ghw1_prime4 = Hzzcoupl(gHIGGS_VV_1_PRIME4)
      ghw1_prime5 = Hzzcoupl(gHIGGS_VV_1_PRIME5)
      ghw1_prime6 = Hzzcoupl(gHIGGS_VV_1_PRIME6)
      ghw1_prime7 = Hzzcoupl(gHIGGS_VV_1_PRIME7)
      ghw2_prime = Hzzcoupl(gHIGGS_VV_2_PRIME)
      ghw2_prime2 = Hzzcoupl(gHIGGS_VV_2_PRIME2)
      ghw2_prime3 = Hzzcoupl(gHIGGS_VV_2_PRIME3)
      ghw2_prime4 = Hzzcoupl(gHIGGS_VV_2_PRIME4)
      ghw2_prime5 = Hzzcoupl(gHIGGS_VV_2_PRIME5)
      ghw2_prime6 = Hzzcoupl(gHIGGS_VV_2_PRIME6)
      ghw2_prime7 = Hzzcoupl(gHIGGS_VV_2_PRIME7)
      ghw3_prime = Hzzcoupl(gHIGGS_VV_3_PRIME)
      ghw3_prime2 = Hzzcoupl(gHIGGS_VV_3_PRIME2)
      ghw3_prime3 = Hzzcoupl(gHIGGS_VV_3_PRIME3)
      ghw3_prime4 = Hzzcoupl(gHIGGS_VV_3_PRIME4)
      ghw3_prime5 = Hzzcoupl(gHIGGS_VV_3_PRIME5)
      ghw3_prime6 = Hzzcoupl(gHIGGS_VV_3_PRIME6)
      ghw3_prime7 = Hzzcoupl(gHIGGS_VV_3_PRIME7)
      ghw4_prime = Hzzcoupl(gHIGGS_VV_4_PRIME)
      ghw4_prime2 = Hzzcoupl(gHIGGS_VV_4_PRIME2)
      ghw4_prime3 = Hzzcoupl(gHIGGS_VV_4_PRIME3)
      ghw4_prime4 = Hzzcoupl(gHIGGS_VV_4_PRIME4)
      ghw4_prime5 = Hzzcoupl(gHIGGS_VV_4_PRIME5)
      ghw4_prime6 = Hzzcoupl(gHIGGS_VV_4_PRIME6)
      ghw4_prime7 = Hzzcoupl(gHIGGS_VV_4_PRIME7)
   endif
!   /***** END REGULAR RESONANCE *****/
   ! Not setting second resonance since MCFM couplings are not present
   h2mass=-1d0 ! This should be enough to disable second resonance

   return
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
