MODULE ModMCFMWrapper
implicit none
private

public :: MCFM_firsttime
public :: Setup_MCFM_qqVVqq_firsttime,Setup_MCFM_qqVVqq,EvalAmp_qqVVqq

integer, parameter :: mxpart=14
integer, parameter :: mxdim=26
integer, parameter :: nf=5
character*30 :: MCFM_runstring

integer, parameter :: &
   gHIGGS_KAPPA=1, &
   gHIGGS_KAPPA_TILDE=2, &
   SIZE_HQQ=2
integer, parameter :: &
   gHIGGS_GG_2=1, &
   gHIGGS_GG_3=2, &
   gHIGGS_GG_4=3, &
   SIZE_HGG=3
integer, parameter :: &
   gHIGGS_VV_1=1, &
   gHIGGS_VV_2=2, &
   gHIGGS_VV_3=3, &
   gHIGGS_VV_4=4, &
   gHIGGS_ZA_2=5, &
   gHIGGS_ZA_3=6, &
   gHIGGS_ZA_4=7, &
   gHIGGS_AA_2=8, &
   gHIGGS_AA_3=9, &
   gHIGGS_AA_4=10, &
   gHIGGS_VV_1_PRIME=11, &
   gHIGGS_VV_1_PRIME2=12, &
   gHIGGS_VV_1_PRIME3=13, &
   gHIGGS_VV_1_PRIME4=14, &
   gHIGGS_VV_1_PRIME5=15, &
   gHIGGS_VV_2_PRIME=16, &
   gHIGGS_VV_2_PRIME2=17, &
   gHIGGS_VV_2_PRIME3=18, &
   gHIGGS_VV_2_PRIME4=19, &
   gHIGGS_VV_2_PRIME5=20, &
   gHIGGS_VV_3_PRIME=21, &
   gHIGGS_VV_3_PRIME2=22, &
   gHIGGS_VV_3_PRIME3=23, &
   gHIGGS_VV_3_PRIME4=24, &
   gHIGGS_VV_3_PRIME5=25, &
   gHIGGS_VV_4_PRIME=26, &
   gHIGGS_VV_4_PRIME2=27, &
   gHIGGS_VV_4_PRIME3=28, &
   gHIGGS_VV_4_PRIME4=29, &
   gHIGGS_VV_4_PRIME5=30, &
   gHIGGS_ZA_1_PRIME2=31, &
   gHIGGS_VV_1_PRIME6=32, &
   gHIGGS_VV_1_PRIME7=33, &
   gHIGGS_VV_2_PRIME6=34, &
   gHIGGS_VV_2_PRIME7=35, &
   gHIGGS_VV_3_PRIME6=36, &
   gHIGGS_VV_3_PRIME7=37, &
   gHIGGS_VV_4_PRIME6=38, &
   gHIGGS_VV_4_PRIME7=39, &
   SIZE_HVV=39
integer, parameter :: &
   LambdaHIGGS_QSQ_VV_1=1, &
   LambdaHIGGS_QSQ_VV_2=2, &
   LambdaHIGGS_QSQ_VV_3=3, &
   LambdaHIGGS_QSQ_VV_4=4, &
   SIZE_HVV_LAMBDAQSQ=4
integer, parameter :: &
   cLambdaHIGGS_VV_QSQ1=1, &
   cLambdaHIGGS_VV_QSQ2=2, &
   cLambdaHIGGS_VV_QSQ12=3, &
   SIZE_HVV_CQSQ=3
integer, parameter :: &
   gZPRIME_QQ_LEFT=1, &
   gZPRIME_QQ_RIGHT=2, &
   SIZE_ZQQ=2
integer, parameter :: &
   gZPRIME_VV_1=1, &
   gZPRIME_VV_2=2, &
   SIZE_ZVV=2
integer, parameter :: &
   gGRAVITON_QQ_LEFT=1, &
   gGRAVITON_QQ_RIGHT=2, &
   SIZE_GQQ=2
integer, parameter :: &
   gGRAVITON_GG_1=1, &
   gGRAVITON_GG_2=2, &
   gGRAVITON_GG_3=3, &
   gGRAVITON_GG_4=4, &
   gGRAVITON_GG_5=5, &
   SIZE_GGG=5
integer, parameter :: &
   gGRAVITON_VV_1=1, &
   gGRAVITON_VV_2=2, &
   gGRAVITON_VV_3=3, &
   gGRAVITON_VV_4=4, &
   gGRAVITON_VV_5=5, &
   gGRAVITON_VV_6=6, &
   gGRAVITON_VV_7=7, &
   gGRAVITON_VV_8=8, &
   gGRAVITON_VV_9=9, &
   gGRAVITON_VV_10=10, &
   SIZE_GVV=10


contains


subroutine MCFM_firsttime()
   ! This function should NEVER use anything from ModParameters
   ! It should just transfer parameters
   implicit none
   ! Cannot just pass these as arguments. Need to specify the char lengths
   ! Otherwise MCFM crashes.
   character*72 :: inputfile,workdir

   ! For Init_MCFMCommon_energy
   double precision sqrts

   ! For Init_MCFMCommon_masses
   double precision &
      md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
      mel_in, mmu_in, mtau_in, &
      hmass_in, hwidth_in, &
      h2mass_in, h2width_in, &
      wmass_in, wwidth_in, &
      zmass_in, zwidth_in, &
      twidth_in, &
      tauwidth_in

   ! For Init_MCFMCommon_ewinput
   double precision Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in

   ! For Init_MCFMCommon_spinzerohiggs_anomcoupl
   double complex Hggcoupl(1:SIZE_HGG)
   double complex Httcoupl(1:SIZE_HQQ)
   double complex Hbbcoupl(1:SIZE_HQQ)
   double complex Hg4g4coupl(1:SIZE_HGG)
   double complex Ht4t4coupl(1:SIZE_HQQ)
   double complex Hb4b4coupl(1:SIZE_HQQ)

   double complex Hzzcoupl(1:SIZE_HVV)
   double complex Hwwcoupl(1:SIZE_HVV)
   
   double precision Hb4b4_mb_4gen
   double precision Ht4t4_mt_4gen

   double precision HLambdaBSM
   double precision HLambda_Q
   double precision HLambda_zgs1
   double precision HzzLambda(1:SIZE_HVV_LAMBDAQSQ)
   double precision HwwLambda(1:SIZE_HVV_LAMBDAQSQ)

   double precision HzzLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)
   double precision HwwLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)

   integer HzzCLambda_qsq(1:SIZE_HVV_CQSQ)
   integer HwwCLambda_qsq(1:SIZE_HVV_CQSQ)

   integer separateWWZZcouplings

   
   double complex H2zzcoupl(1:SIZE_HVV)
   double complex H2wwcoupl(1:SIZE_HVV)
   
   double precision H2LambdaBSM
   double precision H2Lambda_Q
   double precision H2Lambda_zgs1
   double precision H2zzLambda(1:SIZE_HVV_LAMBDAQSQ)
   double precision H2wwLambda(1:SIZE_HVV_LAMBDAQSQ)

   double precision H2zzLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)
   double precision H2wwLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)

   integer H2zzCLambda_qsq(1:SIZE_HVV_CQSQ)
   integer H2wwCLambda_qsq(1:SIZE_HVV_CQSQ)
   
   
   
   inputfile='input.DAT'
   workdir='./'

   call mcfm_init(inputfile,workdir)

   ! For Init_MCFMCommon_energy
   call GetColliderEnergy(sqrts)
   call Init_MCFMCommon_energy(sqrts)

   ! For Init_MCFMCommon_masses
   call GetMassesWidths( &
      md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
      mel_in, mmu_in, mtau_in, &
      hmass_in, hwidth_in, &
      h2mass_in, h2width_in, &
      wmass_in, wwidth_in, &
      zmass_in, zwidth_in, &
      twidth_in, &
      tauwidth_in &
      )
      
   call Init_MCFMCommon_masses( &
      md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
      mel_in, mmu_in, mtau_in, &
      hmass_in, hwidth_in, &
      wmass_in, wwidth_in, &
      zmass_in, zwidth_in, &
      twidth_in, &
      tauwidth_in &
      )

   ! EW parameters
   call Init_MCFMCommon_ewscheme()
   call GetEWInputs(Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in)
   call Init_MCFMCommon_ewinput(Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in)
   call coupling() ! Let MCFM calculate couplings from ewinput as it likes

   call couplzajk()

   ! For Init_MCFMCommon_spinzerohiggs_anomcoupl
   call GetLambdaBSM(HLambdaBSM,H2LambdaBSM)
   call GetSpinZeroGGCouplings(Hggcoupl)
   call GetSpinZeroQQCouplings(Httcoupl)
   Hbbcoupl(:)=Httcoupl(:)
   Hg4g4coupl(:)=0d0
   Ht4t4coupl(:)=0d0
   Hb4b4coupl(:)=0d0
   Ht4t4_mt_4gen=10000d0
   Hb4b4_mb_4gen=10000d0
   call GetSpinZeroVVCouplings(1,Hzzcoupl, HzzCLambda_qsq, HzzLambda_qsq, HzzLambda, HLambda_zgs1, HLambda_Q, .false.)
   call GetSpinZeroVVCouplings(1,Hwwcoupl, HwwCLambda_qsq, HwwLambda_qsq, HwwLambda, HLambda_zgs1, HLambda_Q, .true.)
   call GetDistinguishWWCouplingsFlag(separateWWZZcouplings)

!  second resonance
   call GetSpinZeroVVCouplings(2,H2zzcoupl, H2zzCLambda_qsq, H2zzLambda_qsq, H2zzLambda, H2Lambda_zgs1, H2Lambda_Q, .false.)
   call GetSpinZeroVVCouplings(2,H2wwcoupl, H2wwCLambda_qsq, H2wwLambda_qsq, H2wwLambda, H2Lambda_zgs1, H2Lambda_Q, .true.)

   call qlinit()
 
 
   call Init_MCFMCommon_spinzerohiggs_anomcoupl(Hggcoupl,Httcoupl,Hbbcoupl,Hg4g4coupl,Ht4t4coupl,Hb4b4coupl,    &
                                                Hzzcoupl,Hwwcoupl,     &
                                                Hb4b4_mb_4gen,Ht4t4_mt_4gen,HLambdaBSM,HLambda_Q,HLambda_zgs1,HzzLambda,HwwLambda,    &
                                                HzzLambda_qsq,HwwLambda_qsq,HzzCLambda_qsq,HwwCLambda_qsq,   &
                                                separateWWZZcouplings,  &
                                                h2mass_in, h2width_in, &
                                                H2zzcoupl,H2wwcoupl,   &
                                                H2zzCLambda_qsq,H2wwCLambda_qsq,H2zzLambda_qsq,H2wwLambda_qsq,H2zzLambda,  &
                                                H2wwLambda,H2Lambda_zgs1,H2Lambda_Q,H2LambdaBSM &
                                               )
          
 
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

subroutine Init_MCFMCommon_masses( &
   md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
   mel_in, mmu_in, mtau_in, &
   hmass_in, hwidth_in, &
   wmass_in, wwidth_in, &
   zmass_in, zwidth_in, &
   twidth_in, &
   tauwidth_in &
   )
   implicit none
   double precision &
   md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
   mel_in, mmu_in, mtau_in, &
   hmass_in, hwidth_in, &
   wmass_in, wwidth_in, &
   zmass_in, zwidth_in, &
   twidth_in, &
   tauwidth_in
   ! MCFM declarations
   double precision &
   md, mu, ms, mc, mb, mt, &
   mel, mmu, mtau, &
   hmass, hwidth, &
   wmass, wwidth, &
   zmass, zwidth, &
   twidth, &
   tauwidth, &
   mtausq, mcsq, mbsq
   common/masses/ &
   md, mu, ms, mc, mb, mt, &
   mel, mmu, mtau, &
   hmass, hwidth, &
   wmass, wwidth, &
   zmass, zwidth, &
   twidth, &
   tauwidth, &
   mtausq, mcsq, mbsq

   md = md_in
   mu = mu_in
   ms = ms_in

   mc = mc_in
   mcsq = mc**2

   mb = mb_in
   mbsq = mb**2

   mt = mt_in
   twidth = twidth_in

   mel = mel_in
   mmu = mmu_in

   mtau = mtau_in
   tauwidth = tauwidth_in
   mtausq = mtau**2

   hmass = hmass_in
   hwidth = hwidth_in

   wmass = wmass_in
   wwidth = wwidth_in
   zmass = zmass_in
   zwidth = zwidth_in
end subroutine











subroutine Init_MCFMCommon_spinzerohiggs_anomcoupl( &
   Hggcoupl, &
   Httcoupl, &
   Hbbcoupl, &
   Hg4g4coupl, &
   Ht4t4coupl, &
   Hb4b4coupl, &

   Hzzcoupl, &
   Hwwcoupl, &

   Hb4b4_mb_4gen, &
   Ht4t4_mt_4gen, &

   HLambdaBSM, &
   HLambda_Q, &
   HLambda_zgs1, &
   HzzLambda, &
   HwwLambda, &

   HzzLambda_qsq, &
   HwwLambda_qsq, &
   HzzCLambda_qsq, &
   HwwCLambda_qsq, &

   separateWWZZcouplings, &
   
   h2mass_in, h2width_in, & 
   H2zzcoupl,H2wwcoupl,   &
   H2zzCLambda_qsq,H2wwCLambda_qsq,H2zzLambda_qsq,H2wwLambda_qsq,H2zzLambda,H2wwLambda,H2Lambda_zgs1,H2Lambda_Q,H2LambdaBSM &
   )
   implicit none
   double complex Hggcoupl(1:SIZE_HGG)
   double complex Httcoupl(1:SIZE_HQQ)
   double complex Hbbcoupl(1:SIZE_HQQ)
   double complex Hg4g4coupl(1:SIZE_HGG)
   double complex Ht4t4coupl(1:SIZE_HQQ)
   double complex Hb4b4coupl(1:SIZE_HQQ)

   double complex Hzzcoupl(1:SIZE_HVV)
   double complex Hwwcoupl(1:SIZE_HVV)

   double precision Hb4b4_mb_4gen
   double precision Ht4t4_mt_4gen

   double precision HLambdaBSM
   double precision HLambda_Q
   double precision HLambda_zgs1
   double precision HzzLambda(1:SIZE_HVV_LAMBDAQSQ)
   double precision HwwLambda(1:SIZE_HVV_LAMBDAQSQ)

   double precision HzzLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)
   double precision HwwLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)

   integer HzzCLambda_qsq(1:SIZE_HVV_CQSQ)
   integer HwwCLambda_qsq(1:SIZE_HVV_CQSQ)

   integer separateWWZZcouplings
   
   double precision h2mass_in, h2width_in

   double complex H2zzcoupl(1:SIZE_HVV)
   double complex H2wwcoupl(1:SIZE_HVV)
   
   double precision H2LambdaBSM
   double precision H2Lambda_Q
   double precision H2Lambda_zgs1
   double precision H2zzLambda(1:SIZE_HVV_LAMBDAQSQ)
   double precision H2wwLambda(1:SIZE_HVV_LAMBDAQSQ)

   double precision H2zzLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)
   double precision H2wwLambda_qsq(1:SIZE_HVV_LAMBDAQSQ,1:SIZE_HVV_CQSQ)

   integer H2zzCLambda_qsq(1:SIZE_HVV_CQSQ)
   integer H2wwCLambda_qsq(1:SIZE_HVV_CQSQ)
      
   
   
   
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


   
      common/spinzerohiggs_anomcoupl/     &
        AllowAnomalousCouplings,     &
        distinguish_HWWcouplings,     &
        AnomalCouplPR,AnomalCouplDK,     &
        channeltoggle_stu,vvhvvtoggle_vbfvh,     &
        cz_q1sq,cz_q2sq,cz_q12sq,     &
        cw_q1sq,cw_q2sq,cw_q12sq,     &
        c2z_q1sq,c2z_q2sq,c2z_q12sq,     &
        c2w_q1sq,c2w_q2sq,c2w_q12sq,     &

        mb_4gen,mt_4gen,     &

        LambdaBSM,Lambda_Q,     &
        Lambda_zgs1,     &
        Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4,     &
        Lambda_z11,Lambda_z21,Lambda_z31,Lambda_z41,     &
        Lambda_z12,Lambda_z22,Lambda_z32,Lambda_z42,     &
        Lambda_z10,Lambda_z20,Lambda_z30,Lambda_z40,     &
        Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4,     &
        Lambda_w11,Lambda_w21,Lambda_w31,Lambda_w41,     &
        Lambda_w12,Lambda_w22,Lambda_w32,Lambda_w42,     &
        Lambda_w10,Lambda_w20,Lambda_w30,Lambda_w40,     &

        h2mass,h2width,     &

        Lambda2BSM,Lambda2_Q,     &
        Lambda2_zgs1,     &
        Lambda2_z1,Lambda2_z2,Lambda2_z3,Lambda2_z4,     &
        Lambda2_z11,Lambda2_z21,Lambda2_z31,Lambda2_z41,     &
        Lambda2_z12,Lambda2_z22,Lambda2_z32,Lambda2_z42,     &
        Lambda2_z10,Lambda2_z20,Lambda2_z30,Lambda2_z40,     &
        Lambda2_w1,Lambda2_w2,Lambda2_w3,Lambda2_w4,     &
        Lambda2_w11,Lambda2_w21,Lambda2_w31,Lambda2_w41,     &
        Lambda2_w12,Lambda2_w22,Lambda2_w32,Lambda2_w42,     &
        Lambda2_w10,Lambda2_w20,Lambda2_w30,Lambda2_w40,     &


        kappa_top,kappa_tilde_top,     &
        kappa_bot,kappa_tilde_bot,     &
        ghg2,ghg3,ghg4,     &
        kappa_4gen_top,kappa_tilde_4gen_top,     &
        kappa_4gen_bot,kappa_tilde_4gen_bot,     &
        ghg2_4gen,ghg3_4gen,ghg4_4gen,     &

        ghz1,ghz2,ghz3,ghz4,     &
        ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime,     &
        ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2,     &
        ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3,     &
        ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4,     &
        ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5,     &
        ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6,     &
        ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7,     &

        ghzgs1_prime2,ghzgs2,ghzgs3,ghzgs4,     &
        ghgsgs2,ghgsgs3,ghgsgs4,     &

        ghw1,ghw2,ghw3,ghw4,     &
        ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime,     &
        ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2,     &
        ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3,     &
        ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4,     &
        ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5,     &
        ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6,     &
        ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7,     &


        kappa2_top,kappa2_tilde_top,     &
        kappa2_bot,kappa2_tilde_bot,     &
        gh2g2,gh2g3,gh2g4,     &
        kappa2_4gen_top,kappa2_tilde_4gen_top,     &
        kappa2_4gen_bot,kappa2_tilde_4gen_bot,     &
        gh2g2_4gen,gh2g3_4gen,gh2g4_4gen,     &

        gh2z1,gh2z2,gh2z3,gh2z4,     &
        gh2z1_prime,gh2z2_prime,gh2z3_prime,gh2z4_prime,     &
        gh2z1_prime2,gh2z2_prime2,gh2z3_prime2,gh2z4_prime2,     &
        gh2z1_prime3,gh2z2_prime3,gh2z3_prime3,gh2z4_prime3,     &
        gh2z1_prime4,gh2z2_prime4,gh2z3_prime4,gh2z4_prime4,     &
        gh2z1_prime5,gh2z2_prime5,gh2z3_prime5,gh2z4_prime5,     &
        gh2z1_prime6,gh2z2_prime6,gh2z3_prime6,gh2z4_prime6,     &
        gh2z1_prime7,gh2z2_prime7,gh2z3_prime7,gh2z4_prime7,     &

        gh2zgs1_prime2,gh2zgs2,gh2zgs3,gh2zgs4,     &
        gh2gsgs2,gh2gsgs3,gh2gsgs4,     &

        gh2w1,gh2w2,gh2w3,gh2w4,     &
        gh2w1_prime,gh2w2_prime,gh2w3_prime,gh2w4_prime,     &
        gh2w1_prime2,gh2w2_prime2,gh2w3_prime2,gh2w4_prime2,     &
        gh2w1_prime3,gh2w2_prime3,gh2w3_prime3,gh2w4_prime3,     &
        gh2w1_prime4,gh2w2_prime4,gh2w3_prime4,gh2w4_prime4,     &
        gh2w1_prime5,gh2w2_prime5,gh2w3_prime5,gh2w4_prime5,     &
        gh2w1_prime6,gh2w2_prime6,gh2w3_prime6,gh2w4_prime6,     &
        gh2w1_prime7,gh2w2_prime7,gh2w3_prime7,gh2w4_prime7    


   
   
   
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



!   /***** SECOND RESONANCE *****/

   h2mass= h2mass_in
   h2width = h2width_in
   
   Lambda2BSM = H2LambdaBSM
   Lambda2_Q = H2Lambda_Q

   Lambda2_zgs1 = H2Lambda_zgs1
   Lambda2_z1 = H2zzLambda(LambdaHIGGS_QSQ_VV_1)
   Lambda2_z2 = H2zzLambda(LambdaHIGGS_QSQ_VV_2)
   Lambda2_z3 = H2zzLambda(LambdaHIGGS_QSQ_VV_3)
   Lambda2_z4 = H2zzLambda(LambdaHIGGS_QSQ_VV_4)

   c2z_q1sq = H2zzCLambda_qsq(cLambdaHIGGS_VV_QSQ1)
   Lambda2_z11 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1)
   Lambda2_z21 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1)
   Lambda2_z31 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1)
   Lambda2_z41 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1)
   c2z_q2sq = H2zzCLambda_qsq(cLambdaHIGGS_VV_QSQ2)
   Lambda2_z12 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2)
   Lambda2_z22 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2)
   Lambda2_z32 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2)
   Lambda2_z42 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2)
   c2z_q12sq = H2zzCLambda_qsq(cLambdaHIGGS_VV_QSQ12)
   Lambda2_z10 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12)
   Lambda2_z20 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12)
   Lambda2_z30 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12)
   Lambda2_z40 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12)

!    kappa_top = Httcoupl(gHIGGS_KAPPA)
!    kappa_bot = Hbbcoupl(gHIGGS_KAPPA)
!    kappa_tilde_top = Httcoupl(gHIGGS_KAPPA_TILDE)
!    kappa_tilde_bot = Hbbcoupl(gHIGGS_KAPPA_TILDE)
!    ghg2 = Hggcoupl(gHIGGS_GG_2)
!    ghg3 = Hggcoupl(gHIGGS_GG_3)
!    ghg4 = Hggcoupl(gHIGGS_GG_4)
!    mb_4gen = Hb4b4_mb_4gen
!    mt_4gen = Ht4t4_mt_4gen
!    kappa_4gen_top = Ht4t4coupl(gHIGGS_KAPPA)
!    kappa_4gen_bot = Hb4b4coupl(gHIGGS_KAPPA)
!    kappa_tilde_4gen_top = Ht4t4coupl(gHIGGS_KAPPA_TILDE)
!    kappa_tilde_4gen_bot = Hb4b4coupl(gHIGGS_KAPPA_TILDE)
!    ghg2_4gen = Hg4g4coupl(gHIGGS_GG_2)
!    ghg3_4gen = Hg4g4coupl(gHIGGS_GG_3)
!    ghg4_4gen = Hg4g4coupl(gHIGGS_GG_4)

   gh2z1 = H2zzcoupl(gHIGGS_VV_1)
   gh2z2 = H2zzcoupl(gHIGGS_VV_2)
   gh2z3 = H2zzcoupl(gHIGGS_VV_3)
   gh2z4 = H2zzcoupl(gHIGGS_VV_4)
   gh2z1_prime = H2zzcoupl(gHIGGS_VV_1_PRIME)
   gh2z1_prime2 = H2zzcoupl(gHIGGS_VV_1_PRIME2)
   gh2z1_prime3 = H2zzcoupl(gHIGGS_VV_1_PRIME3)
   gh2z1_prime4 = H2zzcoupl(gHIGGS_VV_1_PRIME4)
   gh2z1_prime5 = H2zzcoupl(gHIGGS_VV_1_PRIME5)
   gh2z1_prime6 = H2zzcoupl(gHIGGS_VV_1_PRIME6)
   gh2z1_prime7 = H2zzcoupl(gHIGGS_VV_1_PRIME7)
   gh2z2_prime = H2zzcoupl(gHIGGS_VV_2_PRIME)
   gh2z2_prime2 = H2zzcoupl(gHIGGS_VV_2_PRIME2)
   gh2z2_prime3 = H2zzcoupl(gHIGGS_VV_2_PRIME3)
   gh2z2_prime4 = H2zzcoupl(gHIGGS_VV_2_PRIME4)
   gh2z2_prime5 = H2zzcoupl(gHIGGS_VV_2_PRIME5)
   gh2z2_prime6 = H2zzcoupl(gHIGGS_VV_2_PRIME6)
   gh2z2_prime7 = H2zzcoupl(gHIGGS_VV_2_PRIME7)
   gh2z3_prime = H2zzcoupl(gHIGGS_VV_3_PRIME)
   gh2z3_prime2 = H2zzcoupl(gHIGGS_VV_3_PRIME2)
   gh2z3_prime3 = H2zzcoupl(gHIGGS_VV_3_PRIME3)
   gh2z3_prime4 = H2zzcoupl(gHIGGS_VV_3_PRIME4)
   gh2z3_prime5 = H2zzcoupl(gHIGGS_VV_3_PRIME5)
   gh2z3_prime6 = H2zzcoupl(gHIGGS_VV_3_PRIME6)
   gh2z3_prime7 = H2zzcoupl(gHIGGS_VV_3_PRIME7)
   gh2z4_prime = H2zzcoupl(gHIGGS_VV_4_PRIME)
   gh2z4_prime2 = H2zzcoupl(gHIGGS_VV_4_PRIME2)
   gh2z4_prime3 = H2zzcoupl(gHIGGS_VV_4_PRIME3)
   gh2z4_prime4 = H2zzcoupl(gHIGGS_VV_4_PRIME4)
   gh2z4_prime5 = H2zzcoupl(gHIGGS_VV_4_PRIME5)
   gh2z4_prime6 = H2zzcoupl(gHIGGS_VV_4_PRIME6)
   gh2z4_prime7 = H2zzcoupl(gHIGGS_VV_4_PRIME7)

   gh2zgs1_prime2 = H2zzcoupl(gHIGGS_ZA_1_PRIME2)
   gh2zgs2 = H2zzcoupl(gHIGGS_ZA_2)
   gh2zgs3 = H2zzcoupl(gHIGGS_ZA_3)
   gh2zgs4 = H2zzcoupl(gHIGGS_ZA_4)
   gh2gsgs2 = H2zzcoupl(gHIGGS_AA_2)
   gh2gsgs3 = H2zzcoupl(gHIGGS_AA_3)
   gh2gsgs4 = H2zzcoupl(gHIGGS_AA_4)




   if (distinguish_Hwwcouplings .eq. 1) then
      Lambda2_w1 = H2wwLambda(LambdaHIGGS_QSQ_VV_1)
      Lambda2_w2 = H2wwLambda(LambdaHIGGS_QSQ_VV_2)
      Lambda2_w3 = H2wwLambda(LambdaHIGGS_QSQ_VV_3)
      Lambda2_w4 = H2wwLambda(LambdaHIGGS_QSQ_VV_4)

      c2w_q1sq = H2wwCLambda_qsq(cLambdaHIGGS_VV_QSQ1)
      Lambda2_w11 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1)
      Lambda2_w21 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1)
      Lambda2_w31 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1)
      Lambda2_w41 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1)
      c2w_q2sq = H2wwCLambda_qsq(cLambdaHIGGS_VV_QSQ2)
      Lambda2_w12 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2)
      Lambda2_w22 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2)
      Lambda2_w32 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2)
      Lambda2_w42 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2)
      c2w_q12sq = H2wwCLambda_qsq(cLambdaHIGGS_VV_QSQ12)
      Lambda2_w10 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12)
      Lambda2_w20 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12)
      Lambda2_w30 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12)
      Lambda2_w40 = H2wwLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12)

      gh2w1 = H2wwcoupl(gHIGGS_VV_1)
      gh2w2 = H2wwcoupl(gHIGGS_VV_2)
      gh2w3 = H2wwcoupl(gHIGGS_VV_3)
      gh2w4 = H2wwcoupl(gHIGGS_VV_4)
      gh2w1_prime = H2wwcoupl(gHIGGS_VV_1_PRIME)
      gh2w1_prime2 = H2wwcoupl(gHIGGS_VV_1_PRIME2)
      gh2w1_prime3 = H2wwcoupl(gHIGGS_VV_1_PRIME3)
      gh2w1_prime4 = H2wwcoupl(gHIGGS_VV_1_PRIME4)
      gh2w1_prime5 = H2wwcoupl(gHIGGS_VV_1_PRIME5)
      gh2w1_prime6 = H2wwcoupl(gHIGGS_VV_1_PRIME6)
      gh2w1_prime7 = H2wwcoupl(gHIGGS_VV_1_PRIME7)
      gh2w2_prime = H2wwcoupl(gHIGGS_VV_2_PRIME)
      gh2w2_prime2 = H2wwcoupl(gHIGGS_VV_2_PRIME2)
      gh2w2_prime3 = H2wwcoupl(gHIGGS_VV_2_PRIME3)
      gh2w2_prime4 = H2wwcoupl(gHIGGS_VV_2_PRIME4)
      gh2w2_prime5 = H2wwcoupl(gHIGGS_VV_2_PRIME5)
      gh2w2_prime6 = H2wwcoupl(gHIGGS_VV_2_PRIME6)
      gh2w2_prime7 = H2wwcoupl(gHIGGS_VV_2_PRIME7)
      gh2w3_prime = H2wwcoupl(gHIGGS_VV_3_PRIME)
      gh2w3_prime2 = H2wwcoupl(gHIGGS_VV_3_PRIME2)
      gh2w3_prime3 = H2wwcoupl(gHIGGS_VV_3_PRIME3)
      gh2w3_prime4 = H2wwcoupl(gHIGGS_VV_3_PRIME4)
      gh2w3_prime5 = H2wwcoupl(gHIGGS_VV_3_PRIME5)
      gh2w3_prime6 = H2wwcoupl(gHIGGS_VV_3_PRIME6)
      gh2w3_prime7 = H2wwcoupl(gHIGGS_VV_3_PRIME7)
      gh2w4_prime = H2wwcoupl(gHIGGS_VV_4_PRIME)
      gh2w4_prime2 = H2wwcoupl(gHIGGS_VV_4_PRIME2)
      gh2w4_prime3 = H2wwcoupl(gHIGGS_VV_4_PRIME3)
      gh2w4_prime4 = H2wwcoupl(gHIGGS_VV_4_PRIME4)
      gh2w4_prime5 = H2wwcoupl(gHIGGS_VV_4_PRIME5)
      gh2w4_prime6 = H2wwcoupl(gHIGGS_VV_4_PRIME6)
      gh2w4_prime7 = H2wwcoupl(gHIGGS_VV_4_PRIME7)

   else
      Lambda2_w1 = H2zzLambda(LambdaHIGGS_QSQ_VV_1)
      Lambda2_w2 = H2zzLambda(LambdaHIGGS_QSQ_VV_2)
      Lambda2_w3 = H2zzLambda(LambdaHIGGS_QSQ_VV_3)
      Lambda2_w4 = H2zzLambda(LambdaHIGGS_QSQ_VV_4)

      c2w_q1sq = H2zzCLambda_qsq(cLambdaHIGGS_VV_QSQ1)
      Lambda2_w11 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1)
      Lambda2_w21 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1)
      Lambda2_w31 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1)
      Lambda2_w41 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1)
      c2w_q2sq = H2zzCLambda_qsq(cLambdaHIGGS_VV_QSQ2)
      Lambda2_w12 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2)
      Lambda2_w22 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2)
      Lambda2_w32 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2)
      Lambda2_w42 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2)
      c2w_q12sq = H2zzCLambda_qsq(cLambdaHIGGS_VV_QSQ12)
      Lambda2_w10 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12)
      Lambda2_w20 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12)
      Lambda2_w30 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12)
      Lambda2_w40 = H2zzLambda_qsq(LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12)

      gh2w1 = H2zzcoupl(gHIGGS_VV_1)
      gh2w2 = H2zzcoupl(gHIGGS_VV_2)
      gh2w3 = H2zzcoupl(gHIGGS_VV_3)
      gh2w4 = H2zzcoupl(gHIGGS_VV_4)
      gh2w1_prime = H2zzcoupl(gHIGGS_VV_1_PRIME)
      gh2w1_prime2 = H2zzcoupl(gHIGGS_VV_1_PRIME2)
      gh2w1_prime3 = H2zzcoupl(gHIGGS_VV_1_PRIME3)
      gh2w1_prime4 = H2zzcoupl(gHIGGS_VV_1_PRIME4)
      gh2w1_prime5 = H2zzcoupl(gHIGGS_VV_1_PRIME5)
      gh2w1_prime6 = H2zzcoupl(gHIGGS_VV_1_PRIME6)
      gh2w1_prime7 = H2zzcoupl(gHIGGS_VV_1_PRIME7)
      gh2w2_prime = H2zzcoupl(gHIGGS_VV_2_PRIME)
      gh2w2_prime2 = H2zzcoupl(gHIGGS_VV_2_PRIME2)
      gh2w2_prime3 = H2zzcoupl(gHIGGS_VV_2_PRIME3)
      gh2w2_prime4 = H2zzcoupl(gHIGGS_VV_2_PRIME4)
      gh2w2_prime5 = H2zzcoupl(gHIGGS_VV_2_PRIME5)
      gh2w2_prime6 = H2zzcoupl(gHIGGS_VV_2_PRIME6)
      gh2w2_prime7 = H2zzcoupl(gHIGGS_VV_2_PRIME7)
      gh2w3_prime = H2zzcoupl(gHIGGS_VV_3_PRIME)
      gh2w3_prime2 = H2zzcoupl(gHIGGS_VV_3_PRIME2)
      gh2w3_prime3 = H2zzcoupl(gHIGGS_VV_3_PRIME3)
      gh2w3_prime4 = H2zzcoupl(gHIGGS_VV_3_PRIME4)
      gh2w3_prime5 = H2zzcoupl(gHIGGS_VV_3_PRIME5)
      gh2w3_prime6 = H2zzcoupl(gHIGGS_VV_3_PRIME6)
      gh2w3_prime7 = H2zzcoupl(gHIGGS_VV_3_PRIME7)
      gh2w4_prime = H2zzcoupl(gHIGGS_VV_4_PRIME)
      gh2w4_prime2 = H2zzcoupl(gHIGGS_VV_4_PRIME2)
      gh2w4_prime3 = H2zzcoupl(gHIGGS_VV_4_PRIME3)
      gh2w4_prime4 = H2zzcoupl(gHIGGS_VV_4_PRIME4)
      gh2w4_prime5 = H2zzcoupl(gHIGGS_VV_4_PRIME5)
      gh2w4_prime6 = H2zzcoupl(gHIGGS_VV_4_PRIME6)
      gh2w4_prime7 = H2zzcoupl(gHIGGS_VV_4_PRIME7)
   endif
!   /***** END SECOND RESONANCE *****/


   return
end subroutine




subroutine Init_MCFMCommon_ewscheme()
   implicit none
   integer ewscheme
   common/ewscheme/ewscheme
   ewscheme=3
end subroutine




subroutine Init_MCFMCommon_ewinput(Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in)
   implicit none
   double precision Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in
   double precision Gf_inp,aemmz_inp,xw_inp,wmass_inp,zmass_inp
   common/ewinput/Gf_inp,aemmz_inp,xw_inp,wmass_inp,zmass_inp
   Gf_inp=Gf_inp_in
   aemmz_inp=aemmz_inp_in
   xw_inp=xw_inp_in
   wmass_inp=wmass_inp_in
   zmass_inp=zmass_inp_in
end subroutine



function MCFMParticleLabel(pid, useQJ)
   use ModParameters
   implicit none
   character*2 MCFMParticleLabel
   integer, intent(in) :: pid
   logical, optional :: useQJ
   if (pid.eq.0) then
      MCFMParticleLabel='pp'
      if (present(useQJ)) then
         if (useQJ) then
            MCFMParticleLabel='qj'
         endif
      endif
   else if (pid.eq.Glu_) then
      MCFMParticleLabel='ig'
   else if (pid.eq.Dn_) then
      MCFMParticleLabel='dq'
   else if (pid.eq.Up_) then
      MCFMParticleLabel='uq'
   else if (pid.eq.Str_) then
      MCFMParticleLabel='sq'
   else if (pid.eq.Chm_) then
      MCFMParticleLabel='cq'
   else if (pid.eq.Bot_) then
      MCFMParticleLabel='bq'
   else if (pid.eq.ADn_) then
      MCFMParticleLabel='da'
   else if (pid.eq.AUp_) then
      MCFMParticleLabel='ua'
   else if (pid.eq.AStr_) then
      MCFMParticleLabel='sa'
   else if (pid.eq.AChm_) then
      MCFMParticleLabel='ca'
   else if (pid.eq.ABot_) then
      MCFMParticleLabel='ba'
   else if (abs(pid).eq.abs(Top_)) then
      MCFMParticleLabel='qj' ! No tops
   else if (pid.eq.ElM_) then
      MCFMParticleLabel='el'
   else if (pid.eq.MuM_) then
      MCFMParticleLabel='ml'
   else if (pid.eq.TaM_) then
      MCFMParticleLabel='tl'
   else if (pid.eq.NuE_) then
      MCFMParticleLabel='nl'
   else if (pid.eq.NuM_) then
      MCFMParticleLabel='nl'
   else if (pid.eq.NuT_) then
      MCFMParticleLabel='nl'
   else if (pid.eq.ElP_) then
      MCFMParticleLabel='ea'
   else if (pid.eq.MuP_) then
      MCFMParticleLabel='ma'
   else if (pid.eq.TaP_) then
      MCFMParticleLabel='ta'
   else if (pid.eq.ANuE_) then
      MCFMParticleLabel='na'
   else if (pid.eq.ANuM_) then
      MCFMParticleLabel='na'
   else if (pid.eq.ANuT_) then
      MCFMParticleLabel='na'
   else
      MCFMParticleLabel='  '
   endif
end function
subroutine SetupParticleLabels(pid_MCFM)
   implicit none
   integer, intent(in) :: pid_MCFM(1:mxpart)
   integer :: i
   ! MCFM declarations
   character*2 plabel(mxpart)
   common/plabel/plabel
   do i=1,mxpart
      plabel(i)=MCFMParticleLabel(pid_MCFM(i))
   enddo
end subroutine

! This is so sketchy
subroutine Set_MCFMCommon_DecayZCouple_Wrapper(idferm,whichZ)
use ModParameters
implicit none
integer, intent(in) :: idferm
integer, intent(in) :: whichZ
logical :: isZlep,isZnu,isZup,isZdn
   if(idferm.eq.0) then
      stop("pid(3)=0 is not allowed in qqVVqq processes")
   endif
   isZlep = IsALepton(idferm)
   isZnu = IsANeutrino(idferm)
   isZdn = IsDownTypeQuark(idferm)
   isZup = IsUpTypeQuark(idferm)
   call Set_MCFMCommon_DecayZCouple(isZlep, isZnu, isZup, isZdn, whichZ)
end subroutine



subroutine Set_MCFMCommon_DecayZCouple(isZlep, isZnu, isZup, isZdn, ip)
logical, intent(in) :: isZlep, isZnu, isZup, isZdn
integer, intent(in) :: ip
! Most of these should already be set when coupling() is called.
double precision Q(-nf:nf),tau(-nf:nf)
common/ewcharge/Q,tau
double precision l(nf),r(nf),le,ln,re,rn,sin2w,q1,l1,r1,q2,l2,r2
common/zcouple/l,r,q1,l1,r1,q2,l2,r2,le,ln,re,rn,sin2w
   if (ip.eq.1) then
      if (isZnu) then
         l1=ln;r1=rn;q1=0d0
      endif
      if (isZlep) then
         l1=le;r1=re;q1=-1d0
      endif
      if (isZdn) then
         l1=l(1);r1=r(1);q1=q(1)
      endif
      if (isZup) then
         l1=l(2);r1=r(2);q1=q(2)
      endif
   else if (ip.eq.2) then
      if (isZnu) then
         l2=ln;r2=rn;q2=0d0
      endif
      if (isZlep) then
         l2=le;r2=re;q2=-1d0
      endif
      if (isZdn) then
         l2=l(1);r2=r(1);q2=q(1)
      endif
      if (isZup) then
         l2=l(2);r2=r(2);q2=q(2)
      endif
   endif
end subroutine



! Record ModParameters variables into arrays
! to be used by similar Set* functions to set MCFM variables.
! Doing it this way allows us to pass a zillion arguments (e.g. vvcoupl length)
subroutine GetColliderEnergy(sqrts)
   use ModParameters
   implicit none
   double precision, intent(out) :: sqrts
   sqrts = Collider_Energy/GeV
end subroutine



subroutine GetMassesWidths( &
   md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
   mel_in, mmu_in, mtau_in, &
   hmass_in, hwidth_in, &
   h2mass_in, h2width_in, &
   wmass_in, wwidth_in, &
   zmass_in, zwidth_in, &
   twidth_in, &
   tauwidth_in &
   )
   use ModParameters
   implicit none
   double precision, intent(out) :: &
   md_in, mu_in, ms_in, mc_in, mb_in, mt_in, &
   mel_in, mmu_in, mtau_in, &
   hmass_in, hwidth_in, &
   h2mass_in, h2width_in, &
   wmass_in, wwidth_in, &
   zmass_in, zwidth_in, &
   twidth_in, &
   tauwidth_in

   md_in = getMass(Dn_)/GeV
   mu_in = getMass(Up_)/GeV
   ms_in = getMass(Str_)/GeV
   mc_in = getMass(Chm_)/GeV
   mb_in = getMass(Bot_)/GeV
   mt_in = getMass(Top_)/GeV
   twidth_in = getDecayWidth(Top_)/GeV

   mel_in = getMass(ElM_)/GeV
   mmu_in = getMass(MuM_)/GeV
   mtau_in = getMass(TaM_)/GeV
   tauwidth_in = getDecayWidth(TaM_)/GeV

   hmass_in = getMass(Hig_)/GeV
   hwidth_in = getDecayWidth(Hig_)/GeV
   h2mass_in = getMass(Hig2_)/GeV
   h2width_in = getDecayWidth(Hig2_)/GeV
   wmass_in = getMass(Wp_)/GeV
   wwidth_in = getDecayWidth(Wp_)/GeV
   zmass_in = getMass(Z0_)/GeV
   zwidth_in = getDecayWidth(Z0_)/GeV
end subroutine



subroutine GetLambdaBSM(Lambda_BSM,Lambda2_BSM)
   use ModParameters
   implicit none
   double precision, intent(out) :: Lambda_BSM,Lambda2_BSM
   Lambda_BSM = Lambda/GeV
   Lambda2_BSM = Lambda2/GeV
end subroutine




subroutine GetSpinZeroVVCouplings(NReso,vvcoupl, cqsq, Lambda_qsq, Lambdag, Lambdag_zgs1, Lambdag_Q, useWWcoupl)
   use ModParameters
   implicit none
   integer, intent(in) :: NReso
   complex(8), intent(out) :: vvcoupl(39)
   integer, intent(out) :: cqsq(3)
   double precision, intent(out) :: Lambda_qsq(1:3,1:4)
   double precision, intent(out) :: Lambdag(1:4),Lambdag_zgs1,Lambdag_Q
   logical, intent(in) :: useWWcoupl

   vvcoupl(:)=0d0
   cqsq(:)=0d0
   Lambda_qsq(:,:)=0d0
   Lambdag(:)=0d0


if(NReso.eq.1) then! first (Higgs) resonance
   Lambdag_zgs1 = Lambda_zgs1
   Lambdag_Q = Lambda_Q
   if(.not.useWWcoupl) then
      Lambdag(1) = Lambda_z1
      Lambdag(2) = Lambda_z2
      Lambdag(3) = Lambda_z3
      Lambdag(4) = Lambda_z4

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
      Lambdag(1) = Lambda_w1
      Lambdag(2) = Lambda_w2
      Lambdag(3) = Lambda_w3
      Lambdag(4) = Lambda_w4

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

      vvcoupl(5:10) = czero

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

      vvcoupl(31) = czero

      vvcoupl(32) = ghw1_prime6
      vvcoupl(33) = ghw1_prime7
      vvcoupl(34) = ghw2_prime6
      vvcoupl(35) = ghw2_prime7
      vvcoupl(36) = ghw3_prime6
      vvcoupl(37) = ghw3_prime7
      vvcoupl(38) = ghw4_prime6
      vvcoupl(39) = ghw4_prime7
   endif

   Lambda_qsq(:,:)=Lambda_qsq(:,:)/GeV
   Lambdag(:)=Lambdag(:)/GeV
   Lambdag_zgs1 = Lambdag_zgs1/GeV
   Lambdag_Q = Lambdag_Q/GeV

   
else! second resonance


   Lambdag_zgs1 = Lambda2_zgs1
   Lambdag_Q = Lambda2_Q
   if(.not.useWWcoupl) then
      Lambdag(1) = Lambda2_z1
      Lambdag(2) = Lambda2_z2
      Lambdag(3) = Lambda2_z3
      Lambdag(4) = Lambda2_z4

      cqsq(1) = cz_q1sq
      Lambda_qsq(1,1) = Lambda2_z11
      Lambda_qsq(1,2) = Lambda2_z21
      Lambda_qsq(1,3) = Lambda2_z31
      Lambda_qsq(1,4) = Lambda2_z41
      cqsq(2) = cz_q2sq
      Lambda_qsq(2,1) = Lambda2_z12
      Lambda_qsq(2,2) = Lambda2_z22
      Lambda_qsq(2,3) = Lambda2_z32
      Lambda_qsq(2,4) = Lambda2_z42
      cqsq(3) = cz_q12sq
      Lambda_qsq(3,1) = Lambda2_z10
      Lambda_qsq(3,2) = Lambda2_z20
      Lambda_qsq(3,3) = Lambda2_z30
      Lambda_qsq(3,4) = Lambda2_z40

      vvcoupl(1) = gh2z1
      vvcoupl(2) = gh2z2
      vvcoupl(3) = gh2z3
      vvcoupl(4) = gh2z4

      vvcoupl(5) = gh2zgs2
      vvcoupl(6) = gh2zgs3
      vvcoupl(7) = gh2zgs4
      vvcoupl(8) = gh2gsgs2
      vvcoupl(9) = gh2gsgs3
      vvcoupl(10) = gh2gsgs4

      vvcoupl(11) = gh2z1_prime
      vvcoupl(12) = gh2z1_prime2
      vvcoupl(13) = gh2z1_prime3
      vvcoupl(14) = gh2z1_prime4
      vvcoupl(15) = gh2z1_prime5

      vvcoupl(16) = gh2z2_prime
      vvcoupl(17) = gh2z2_prime2
      vvcoupl(18) = gh2z2_prime3
      vvcoupl(19) = gh2z2_prime4
      vvcoupl(20) = gh2z2_prime5

      vvcoupl(21) = gh2z3_prime
      vvcoupl(22) = gh2z3_prime2
      vvcoupl(23) = gh2z3_prime3
      vvcoupl(24) = gh2z3_prime4
      vvcoupl(25) = gh2z3_prime5

      vvcoupl(26) = gh2z4_prime
      vvcoupl(27) = gh2z4_prime2
      vvcoupl(28) = gh2z4_prime3
      vvcoupl(29) = gh2z4_prime4
      vvcoupl(30) = gh2z4_prime5

      vvcoupl(31) = gh2zgs1_prime2

      vvcoupl(32) = gh2z1_prime6
      vvcoupl(33) = gh2z1_prime7
      vvcoupl(34) = gh2z2_prime6
      vvcoupl(35) = gh2z2_prime7
      vvcoupl(36) = gh2z3_prime6
      vvcoupl(37) = gh2z3_prime7
      vvcoupl(38) = gh2z4_prime6
      vvcoupl(39) = gh2z4_prime7

   else
      Lambdag(1) = Lambda2_w1
      Lambdag(2) = Lambda2_w2
      Lambdag(3) = Lambda2_w3
      Lambdag(4) = Lambda2_w4

      cqsq(1) = cw_q1sq
      Lambda_qsq(1,1) = Lambda2_w11
      Lambda_qsq(1,2) = Lambda2_w21
      Lambda_qsq(1,3) = Lambda2_w31
      Lambda_qsq(1,4) = Lambda2_w41
      cqsq(2) = cw_q2sq
      Lambda_qsq(2,1) = Lambda2_w12
      Lambda_qsq(2,2) = Lambda2_w22
      Lambda_qsq(2,3) = Lambda2_w32
      Lambda_qsq(2,4) = Lambda2_w42
      cqsq(3) = cw_q12sq
      Lambda_qsq(3,1) = Lambda2_w10
      Lambda_qsq(3,2) = Lambda2_w20
      Lambda_qsq(3,3) = Lambda2_w30
      Lambda_qsq(3,4) = Lambda2_w40

      vvcoupl(1) = gh2w1
      vvcoupl(2) = gh2w2
      vvcoupl(3) = gh2w3
      vvcoupl(4) = gh2w4

      vvcoupl(5:10) = czero

      vvcoupl(11) = gh2w1_prime
      vvcoupl(12) = gh2w1_prime2
      vvcoupl(13) = gh2w1_prime3
      vvcoupl(14) = gh2w1_prime4
      vvcoupl(15) = gh2w1_prime5

      vvcoupl(16) = gh2w2_prime
      vvcoupl(17) = gh2w2_prime2
      vvcoupl(18) = gh2w2_prime3
      vvcoupl(19) = gh2w2_prime4
      vvcoupl(20) = gh2w2_prime5

      vvcoupl(21) = gh2w3_prime
      vvcoupl(22) = gh2w3_prime2
      vvcoupl(23) = gh2w3_prime3
      vvcoupl(24) = gh2w3_prime4
      vvcoupl(25) = gh2w3_prime5

      vvcoupl(26) = gh2w4_prime
      vvcoupl(27) = gh2w4_prime2
      vvcoupl(28) = gh2w4_prime3
      vvcoupl(29) = gh2w4_prime4
      vvcoupl(30) = gh2w4_prime5

      vvcoupl(31) = czero

      vvcoupl(32) = gh2w1_prime6
      vvcoupl(33) = gh2w1_prime7
      vvcoupl(34) = gh2w2_prime6
      vvcoupl(35) = gh2w2_prime7
      vvcoupl(36) = gh2w3_prime6
      vvcoupl(37) = gh2w3_prime7
      vvcoupl(38) = gh2w4_prime6
      vvcoupl(39) = gh2w4_prime7
   endif

   Lambda_qsq(:,:)=Lambda_qsq(:,:)/GeV
   Lambdag(:)=Lambdag(:)/GeV
   Lambdag_zgs1 = Lambdag_zgs1/GeV
   Lambdag_Q = Lambdag_Q/GeV

endif! end second resonance 


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
end subroutine


subroutine GetSpinZeroGGCouplings(ggcoupl)
   use ModParameters
   implicit none
   double complex, intent(out) :: ggcoupl(1:3)
   ggcoupl(1) = dcmplx(ghg2)
   ggcoupl(2) = dcmplx(ghg3)
   ggcoupl(3) = dcmplx(ghg4)
end subroutine
subroutine GetSpinZeroQQCouplings(qqcoupl)
   use ModParameters
   implicit none
   double complex, intent(out) :: qqcoupl(1:2)
   qqcoupl(1) = dcmplx(kappa)
   qqcoupl(2) = dcmplx(kappa_tilde)
end subroutine
subroutine GetEWInputs(Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in)
   use ModParameters
   implicit none
   double precision, intent(out) :: Gf_inp_in,aemmz_inp_in,xw_inp_in,wmass_inp_in,zmass_inp_in
   Gf_inp_in=Gf*GeV**2
   aemmz_inp_in=alpha_QED
   xw_inp_in=xw
   wmass_inp_in=M_W/GeV
   zmass_inp_in=M_Z/GeV
end subroutine





!!!!!!!!!!!!
!! qqVVqq !!
!!!!!!!!!!!!
! Setup functions for specific MCFM processes
! qqZZqq/qqWWqq/qqVVqq (or _strong versions)
subroutine Setup_MCFM_qqVVqq_firsttime(iProc)
implicit none
integer, intent(in) :: iProc
integer npart
common/npart/npart
integer nwz
common/nwz/nwz
integer ndim,ncall,itmx,nprn
double precision xl(mxdim),xu(mxdim),acc
common/bveg1/xl,xu,acc,ndim,ncall,itmx,nprn
integer nqcdjets,nqcdstart
common/nqcdjets/nqcdjets,nqcdstart
integer n2,n3
double precision mass2,width2,mass3,width3
common/breit/n2,n3,mass2,width2,mass3,width3
logical interference,bw34_56
common/interference/interference,bw34_56
double precision vsymfact
common/vsymfact/vsymfact
double precision &
md, mu, ms, mc, mb, mt, &
mel, mmu, mtau, &
hmass, hwidth, &
h2mass, h2width, &
wmass, wwidth, &
zmass, zwidth, &
twidth, &
tauwidth, &
mtausq, mcsq, mbsq
common/masses/ &
md, mu, ms, mc, mb, mt, &
mel, mmu, mtau, &
hmass, hwidth, &
wmass, wwidth, &
zmass, zwidth, &
twidth, &
tauwidth, &
mtausq, mcsq, mbsq

   if (iProc.eq.66) then ! Signal-only
      MCFM_runstring="wbfHO"
   else if (iProc.eq.67) then ! Bkg-only
      MCFM_runstring="wbfBO"
   else if (iProc.eq.68) then  ! Signal+bkg
      MCFM_runstring="wbfALL" ! Doesn't matter what it is
   else
      return ! Not qqVVqq, so skip this function
   endif

   npart=6
   nwz=2
   call ckmfill(nwz)
   ndim=16
   nqcdjets=2
   n2=1
   n3=1
   mass2 =zmass
   width2=zwidth
   mass3 =zmass
   width3=zwidth

   ! These two flags are to be changed depending on the 4f flavor
   vsymfact=1d0
   interference=.false.

end subroutine





function Setup_MCFM_qqVVqq(pid_MCFM_in,p_MCFM_in,pid_MCFM,p_MCFM)
implicit none
logical :: Setup_MCFM_qqVVqq
integer, intent(in) :: pid_MCFM_in(1:mxpart)
real(8), intent(in) :: p_MCFM_in(1:mxpart,1:4)
integer, intent(out) :: pid_MCFM(1:mxpart)
real(8), intent(out) :: p_MCFM(1:mxpart,1:4)
integer :: decayOrdering(1:4)
integer :: apartOrdering(1:2)
integer :: ip

logical interference,bw34_56
common/interference/interference,bw34_56
double precision vsymfact
common/vsymfact/vsymfact

double precision l(nf),r(nf),le,ln,re,rn,sin2w,q1,l1,r1,q2,l2,r2
common/zcouple/l,r,q1,l1,r1,q2,l2,r2,le,ln,re,rn,sin2w

   Setup_MCFM_qqVVqq=.false.
   pid_MCFM(:) = pid_MCFM_in(:)
   p_MCFM(:,:) = p_MCFM_in(:,:)

   ! Assign ordered daughter momenta
   call Check_DaughterOrdering_MCFM_qqVVqq(pid_MCFM_in(3:6),decayOrdering)
   do ip=0,3
      p_MCFM(3+ip,:) = p_MCFM_in(3+decayOrdering(ip+1),:)
      pid_MCFM(3+ip) = pid_MCFM_in(3+decayOrdering(ip+1))
   enddo

   ! Assign ordered associated particle momenta
   call Check_APartHash_MCFM_qqVVqq( &
      (/ pid_MCFM_in(1),pid_MCFM_in(2),pid_MCFM_in(7),pid_MCFM_in(8) /), &
      apartOrdering)
   if (apartOrdering(1).eq.-1 .or. apartOrdering(2).eq.-1) then
      return
   else
      do ip=0,1
      p_MCFM(7+ip,:) = p_MCFM_in(7+apartOrdering(ip+1),:)
      pid_MCFM(7+ip) = pid_MCFM_in(7+apartOrdering(ip+1))
      enddo
   endif

   ! Turn 4f interference on as needed
   if( &
   abs(pid_MCFM_in(3)).eq.abs(pid_MCFM_in(5)) .and. abs(pid_MCFM_in(4)).eq.abs(pid_MCFM_in(6)) &
   .and. pid_MCFM_in(3).ne.0 .and. pid_MCFM_in(6).ne.6 &
   ) then
      vsymfact=0.5d0
      interference=.true.
   endif

   ! Set l1, l2
   l1=0d0;r1=0d0;q1=0d0
   l2=0d0;r2=0d0;q2=0d0
   if (pid_MCFM(3).eq.-pid_MCFM(4)) then
      call Set_MCFMCommon_DecayZCouple_Wrapper(pid_MCFM(3),1)
   endif
   if (pid_MCFM(5).eq.-pid_MCFM(6)) then
      call Set_MCFMCommon_DecayZCouple_Wrapper(pid_MCFM(5),2)
   endif
   Setup_MCFM_qqVVqq = .true.

end function

! Subroutines to check and pass the ordering for the decay particles in the "main" system (e.g. H->4f decay)
subroutine Check_DaughterOrdering_MCFM_qqVVqq(idPart,order)
use ModParameters
use ModMisc
implicit none
integer, intent(in) :: idPart(1:4)
integer, intent(out) :: order(1:4)
integer :: ip, idV(1:2)
logical :: isZZ, isWW

character*30 runstring
common/runstring/runstring

   order(:)=(/ 0,1,2,3 /)
   idV(:)=0

   if( idPart(1).ne.0 .and. idPart(2).ne.0) then
      idV(1)=CoupledVertex(idPart(1:2),-1)
   endif
   if( idPart(3).ne.0 .and. idPart(4).ne.0) then
      idV(2)=CoupledVertex(idPart(3:4),-1)
   endif
   isZZ = (idV(1).eq.Z0_ .or. idV(1).eq.0) .and. (idV(2).eq.Z0_ .or. idV(2).eq.0)
   isWW = (idV(1).eq.abs(Wp_) .or. idV(1).eq.0) .and. (idV(2).eq.abs(Wp_) .or. idV(2).eq.0)

   if ( &
   isZZ .and. (&
   IsALepton(idPart(1)) .or. &
   IsDownTypeQuark(idPart(1)) .or. &
   IsANeutrino(idPart(3)) .or. &
   IsUpTypeQuark(idPart(3)) &
   )  .or. &
   isWW .and. (&
   idV(1).eq.Wm_ .and. idV(2).eq.Wp_ &
   ) &
   ) then
      call swap(order(1),order(3))
      call swap(order(2),order(4))
   endif

   if (isWW) then
      call swap(order(1),order(3))
      runstring = trim(MCFM_runstring) // '_ww' ! wbfXY_ww needed to have WW-only final state since we use qq_VVqq process. This needs to be changed somehow if we want to simulate ZZ+WW->4f
   else
      runstring = trim(MCFM_runstring) // '_zz' ! wbfXY_ww needed to have WW-only final state since we use qq_VVqq process. This needs to be changed somehow if we want to simulate ZZ+WW->4f
   endif
end subroutine

! Subroutines to check and pass the ordering for the associated particles
subroutine Check_APartHash_MCFM_qqVVqq(idAPart,order) ! idAPart is in JHU convention
use ModParameters
implicit none
integer, intent(in) :: idAPart(1:4)
integer, intent(out) :: order(1:2) ! Final state ordering; initial state remains the same
integer, parameter :: hashSize = 204
integer :: hash(1:4,1:hashSize),ih
logical  :: outFound

   hash(:,1) = (/ Up_ , Chm_ , Up_ , Chm_ /)
   hash(:,2) = (/ Dn_ , Str_ , Dn_ , Str_ /)
   hash(:,3) = (/ Dn_ , Bot_ , Dn_ , Bot_ /)
   hash(:,4) = (/ Str_ , Bot_ , Str_ , Bot_ /)
   hash(:,5) = (/ Up_ , Str_ , Up_ , Str_ /)
   hash(:,6) = (/ Up_ , Bot_ , Up_ , Bot_ /)
   hash(:,7) = (/ Chm_ , Bot_ , Chm_ , Bot_ /)
   hash(:,8) = (/ Dn_ , Chm_ , Dn_ , Chm_ /)
   hash(:,9) = (/ Dn_ , Up_ , Dn_ , Up_ /)
   hash(:,10) = (/ Str_ , Chm_ , Str_ , Chm_ /)
   hash(:,11) = (/ Dn_ , Chm_ , Up_ , Str_ /)
   hash(:,12) = (/ Up_ , Str_ , Dn_ , Chm_ /)
   hash(:,13) = (/ Up_ , Up_ , Up_ , Up_ /)
   hash(:,14) = (/ Chm_ , Chm_ , Chm_ , Chm_ /)
   hash(:,15) = (/ Dn_ , Dn_ , Dn_ , Dn_ /)
   hash(:,16) = (/ Str_ , Str_ , Str_ , Str_ /)
   hash(:,17) = (/ Bot_ , Bot_ , Bot_ , Bot_ /)
   hash(:,18) = (/ Chm_ , Up_ , Up_ , Chm_ /)
   hash(:,19) = (/ Str_ , Dn_ , Dn_ , Str_ /)
   hash(:,20) = (/ Bot_ , Dn_ , Dn_ , Bot_ /)
   hash(:,21) = (/ Bot_ , Str_ , Str_ , Bot_ /)
   hash(:,22) = (/ Str_ , Up_ , Up_ , Str_ /)
   hash(:,23) = (/ Bot_ , Up_ , Up_ , Bot_ /)
   hash(:,24) = (/ Bot_ , Chm_ , Chm_ , Bot_ /)
   hash(:,25) = (/ Chm_ , Dn_ , Dn_ , Chm_ /)
   hash(:,26) = (/ Up_ , Dn_ , Dn_ , Up_ /)
   hash(:,27) = (/ Chm_ , Str_ , Str_ , Chm_ /)
   hash(:,28) = (/ Chm_ , Dn_ , Up_ , Str_ /)
   hash(:,29) = (/ Str_ , Up_ , Dn_ , Chm_ /)
   hash(:,30) = (/ Up_ , Up_ , Up_ , Up_ /)
   hash(:,31) = (/ Chm_ , Chm_ , Chm_ , Chm_ /)
   hash(:,32) = (/ Dn_ , Dn_ , Dn_ , Dn_ /)
   hash(:,33) = (/ Str_ , Str_ , Str_ , Str_ /)
   hash(:,34) = (/ Bot_ , Bot_ , Bot_ , Bot_ /)

   hash(:,35) = (/ AChm_ , AUp_ , AChm_ , AUp_ /)
   hash(:,36) = (/ AStr_ , ADn_ , AStr_ , ADn_ /)
   hash(:,37) = (/ ABot_ , ADn_ , ABot_ , ADn_ /)
   hash(:,38) = (/ ABot_ , AStr_ , ABot_ , AStr_ /)
   hash(:,39) = (/ AStr_ , AUp_ , AStr_ , AUp_ /)
   hash(:,40) = (/ ABot_ , AUp_ , ABot_ , AUp_ /)
   hash(:,41) = (/ ABot_ , AChm_ , ABot_ , AChm_ /)
   hash(:,42) = (/ AChm_ , ADn_ , AChm_ , ADn_ /)
   hash(:,43) = (/ AUp_ , ADn_ , AUp_ , ADn_ /)
   hash(:,44) = (/ AChm_ , AStr_ , AChm_ , AStr_ /)
   hash(:,45) = (/ AStr_ , AUp_ , AChm_ , ADn_ /)
   hash(:,46) = (/ AChm_ , ADn_ , AStr_ , AUp_ /)
   hash(:,47) = (/ AUp_ , AUp_ , AUp_ , AUp_ /)
   hash(:,48) = (/ AChm_ , AChm_ , AChm_ , AChm_ /)
   hash(:,49) = (/ ADn_ , ADn_ , ADn_ , ADn_ /)
   hash(:,50) = (/ AStr_ , AStr_ , AStr_ , AStr_ /)
   hash(:,51) = (/ ABot_ , ABot_ , ABot_ , ABot_ /)
   hash(:,52) = (/ AUp_ , AChm_ , AChm_ , AUp_ /)
   hash(:,53) = (/ ADn_ , AStr_ , AStr_ , ADn_ /)
   hash(:,54) = (/ ADn_ , ABot_ , ABot_ , ADn_ /)
   hash(:,55) = (/ AStr_ , ABot_ , ABot_ , AStr_ /)
   hash(:,56) = (/ AUp_ , AStr_ , AStr_ , AUp_ /)
   hash(:,57) = (/ AUp_ , ABot_ , ABot_ , AUp_ /)
   hash(:,58) = (/ AChm_ , ABot_ , ABot_ , AChm_ /)
   hash(:,59) = (/ ADn_ , AChm_ , AChm_ , ADn_ /)
   hash(:,60) = (/ ADn_ , AUp_ , AUp_ , ADn_ /)
   hash(:,61) = (/ AStr_ , AChm_ , AChm_ , AStr_ /)
   hash(:,62) = (/ AUp_ , AStr_ , AChm_ , ADn_ /)
   hash(:,63) = (/ ADn_ , AChm_ , AStr_ , AUp_ /)
   hash(:,64) = (/ AUp_ , AUp_ , AUp_ , AUp_ /)
   hash(:,65) = (/ AChm_ , AChm_ , AChm_ , AChm_ /)
   hash(:,66) = (/ ADn_ , ADn_ , ADn_ , ADn_ /)
   hash(:,67) = (/ AStr_ , AStr_ , AStr_ , AStr_ /)
   hash(:,68) = (/ ABot_ , ABot_ , ABot_ , ABot_ /)

   hash(:,69) = (/ AUp_ , Chm_ , AUp_ , Chm_ /)
   hash(:,70) = (/ ADn_ , Str_ , ADn_ , Str_ /)
   hash(:,71) = (/ ADn_ , Bot_ , ADn_ , Bot_ /)
   hash(:,72) = (/ AStr_ , Bot_ , AStr_ , Bot_ /)
   hash(:,73) = (/ AUp_ , Str_ , AUp_ , Str_ /)
   hash(:,74) = (/ AUp_ , Bot_ , AUp_ , Bot_ /)
   hash(:,75) = (/ AChm_ , Bot_ , AChm_ , Bot_ /)
   hash(:,76) = (/ ADn_ , Chm_ , ADn_ , Chm_ /)
   hash(:,77) = (/ ADn_ , Up_ , ADn_ , Up_ /)
   hash(:,78) = (/ AStr_ , Chm_ , AStr_ , Chm_ /)
   hash(:,79) = (/ AUp_ , Chm_ , ADn_ , Str_ /)
   hash(:,80) = (/ ADn_ , Str_ , AUp_ , Chm_ /)
   hash(:,81) = (/ AUp_ , Up_ , AUp_ , Up_ /)
   hash(:,82) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
   hash(:,83) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
   hash(:,84) = (/ AStr_ , Str_ , AStr_ , Str_ /)
   hash(:,85) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)
   hash(:,86) = (/ AChm_ , Up_ , AChm_ , Up_ /)
   hash(:,87) = (/ AStr_ , Dn_ , AStr_ , Dn_ /)
   hash(:,88) = (/ ABot_ , Dn_ , ABot_ , Dn_ /)
   hash(:,89) = (/ ABot_ , Str_ , ABot_ , Str_ /)
   hash(:,90) = (/ AStr_ , Up_ , AStr_ , Up_ /)
   hash(:,91) = (/ ABot_ , Up_ , ABot_ , Up_ /)
   hash(:,92) = (/ ABot_ , Chm_ , ABot_ , Chm_ /)
   hash(:,93) = (/ AChm_ , Dn_ , AChm_ , Dn_ /)
   hash(:,94) = (/ AUp_ , Dn_ , AUp_ , Dn_ /)
   hash(:,95) = (/ AChm_ , Str_ , AChm_ , Str_ /)
   hash(:,96) = (/ AStr_ , Dn_ , AChm_ , Up_ /)
   hash(:,97) = (/ AChm_ , Up_ , AStr_ , Dn_ /)
   hash(:,98) = (/ AUp_ , Up_ , AUp_ , Up_ /)
   hash(:,99) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
   hash(:,100) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
   hash(:,101) = (/ AStr_ , Str_ , AStr_ , Str_ /)
   hash(:,102) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)

   hash(:,103) = (/ Chm_ , AUp_ , AUp_ , Chm_ /)
   hash(:,104) = (/ Str_ , ADn_ , ADn_ , Str_ /)
   hash(:,105) = (/ Bot_ , ADn_ , ADn_ , Bot_ /)
   hash(:,106) = (/ Bot_ , AStr_ , AStr_ , Bot_ /)
   hash(:,107) = (/ Str_ , AUp_ , AUp_ , Str_ /)
   hash(:,108) = (/ Bot_ , AUp_ , AUp_ , Bot_ /)
   hash(:,109) = (/ Bot_ , AChm_ , AChm_ , Bot_ /)
   hash(:,110) = (/ Chm_ , ADn_ , ADn_ , Chm_ /)
   hash(:,111) = (/ Up_ , ADn_ , ADn_ , Up_ /)
   hash(:,112) = (/ Chm_ , AStr_ , AStr_ , Chm_ /)
   hash(:,113) = (/ Chm_ , AUp_ , ADn_ , Str_ /)
   hash(:,114) = (/ Str_ , ADn_ , AUp_ , Chm_ /)
   hash(:,115) = (/ Up_ , AUp_ , AUp_ , Up_ /)
   hash(:,116) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
   hash(:,117) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
   hash(:,118) = (/ Str_ , AStr_ , AStr_ , Str_ /)
   hash(:,119) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)
   hash(:,120) = (/ Up_ , AChm_ , AChm_ , Up_ /)
   hash(:,121) = (/ Dn_ , AStr_ , AStr_ , Dn_ /)
   hash(:,122) = (/ Dn_ , ABot_ , ABot_ , Dn_ /)
   hash(:,123) = (/ Str_ , ABot_ , ABot_ , Str_ /)
   hash(:,124) = (/ Up_ , AStr_ , AStr_ , Up_ /)
   hash(:,125) = (/ Up_ , ABot_ , ABot_ , Up_ /)
   hash(:,126) = (/ Chm_ , ABot_ , ABot_ , Chm_ /)
   hash(:,127) = (/ Dn_ , AChm_ , AChm_ , Dn_ /)
   hash(:,128) = (/ Dn_ , AUp_ , AUp_ , Dn_ /)
   hash(:,129) = (/ Str_ , AChm_ , AChm_ , Str_ /)
   hash(:,130) = (/ Dn_ , AStr_ , AChm_ , Up_ /)
   hash(:,131) = (/ Up_ , AChm_ , AStr_ , Dn_ /)
   hash(:,132) = (/ Up_ , AUp_ , AUp_ , Up_ /)
   hash(:,133) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
   hash(:,134) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
   hash(:,135) = (/ Str_ , AStr_ , AStr_ , Str_ /)
   hash(:,136) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)

   hash(:,137) = (/ Up_ , AUp_ , AChm_ , Chm_ /)
   hash(:,138) = (/ Dn_ , ADn_ , AStr_ , Str_ /)
   hash(:,139) = (/ Dn_ , ADn_ , ABot_ , Bot_ /)
   hash(:,140) = (/ Str_ , AStr_ , ABot_ , Bot_ /)
   hash(:,141) = (/ Up_ , AUp_ , AStr_ , Str_ /)
   hash(:,142) = (/ Up_ , AUp_ , ABot_ , Bot_ /)
   hash(:,143) = (/ Chm_ , AChm_ , ABot_ , Bot_ /)
   hash(:,144) = (/ Dn_ , ADn_ , AChm_ , Chm_ /)
   hash(:,145) = (/ Dn_ , ADn_ , AUp_ , Up_ /)
   hash(:,146) = (/ Str_ , AStr_ , AChm_ , Chm_ /)
   hash(:,147) = (/ Dn_ , AUp_ , AChm_ , Str_ /)
   hash(:,148) = (/ Up_ , ADn_ , AStr_ , Chm_ /)
   hash(:,149) = (/ Up_ , AUp_ , AUp_ , Up_ /)
   hash(:,150) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
   hash(:,151) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
   hash(:,152) = (/ Str_ , AStr_ , AStr_ , Str_ /)
   hash(:,153) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)
   hash(:,154) = (/ Chm_ , AChm_ , AUp_ , Up_ /)
   hash(:,155) = (/ Str_ , AStr_ , ADn_ , Dn_ /)
   hash(:,156) = (/ Bot_ , ABot_ , ADn_ , Dn_ /)
   hash(:,157) = (/ Bot_ , ABot_ , AStr_ , Str_ /)
   hash(:,158) = (/ Str_ , AStr_ , AUp_ , Up_ /)
   hash(:,159) = (/ Bot_ , ABot_ , AUp_ , Up_ /)
   hash(:,160) = (/ Bot_ , ABot_ , AChm_ , Chm_ /)
   hash(:,161) = (/ Chm_ , AChm_ , ADn_ , Dn_ /)
   hash(:,162) = (/ Up_ , AUp_ , ADn_ , Dn_ /)
   hash(:,163) = (/ Chm_ , AChm_ , AStr_ , Str_ /)
   hash(:,164) = (/ Chm_ , AStr_ , ADn_ , Up_ /)
   hash(:,165) = (/ Str_ , AChm_ , AUp_ , Dn_ /)
   hash(:,166) = (/ Up_ , AUp_ , AUp_ , Up_ /)
   hash(:,167) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
   hash(:,168) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
   hash(:,169) = (/ Str_ , AStr_ , AStr_ , Str_ /)
   hash(:,170) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)

   hash(:,171) = (/ AUp_ , Up_ , AChm_ , Chm_ /)
   hash(:,172) = (/ ADn_ , Dn_ , AStr_ , Str_ /)
   hash(:,173) = (/ ADn_ , Dn_ , ABot_ , Bot_ /)
   hash(:,174) = (/ AStr_ , Str_ , ABot_ , Bot_ /)
   hash(:,175) = (/ AUp_ , Up_ , AStr_ , Str_ /)
   hash(:,176) = (/ AUp_ , Up_ , ABot_ , Bot_ /)
   hash(:,177) = (/ AChm_ , Chm_ , ABot_ , Bot_ /)
   hash(:,178) = (/ ADn_ , Dn_ , AChm_ , Chm_ /)
   hash(:,179) = (/ ADn_ , Dn_ , AUp_ , Up_ /)
   hash(:,180) = (/ AStr_ , Str_ , AChm_ , Chm_ /)
   hash(:,181) = (/ AUp_ , Dn_ , AChm_ , Str_ /)
   hash(:,182) = (/ ADn_ , Up_ , AStr_ , Chm_ /)
   hash(:,183) = (/ AUp_ , Up_ , AUp_ , Up_ /)
   hash(:,184) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
   hash(:,185) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
   hash(:,186) = (/ AStr_ , Str_ , AStr_ , Str_ /)
   hash(:,187) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)
   hash(:,188) = (/ AChm_ , Chm_ , AUp_ , Up_ /)
   hash(:,189) = (/ AStr_ , Str_ , ADn_ , Dn_ /)
   hash(:,190) = (/ ABot_ , Bot_ , ADn_ , Dn_ /)
   hash(:,191) = (/ ABot_ , Bot_ , AStr_ , Str_ /)
   hash(:,192) = (/ AStr_ , Str_ , AUp_ , Up_ /)
   hash(:,193) = (/ ABot_ , Bot_ , AUp_ , Up_ /)
   hash(:,194) = (/ ABot_ , Bot_ , AChm_ , Chm_ /)
   hash(:,195) = (/ AChm_ , Chm_ , ADn_ , Dn_ /)
   hash(:,196) = (/ AUp_ , Up_ , ADn_ , Dn_ /)
   hash(:,197) = (/ AChm_ , Chm_ , AStr_ , Str_ /)
   hash(:,198) = (/ AStr_ , Chm_ , ADn_ , Up_ /)
   hash(:,199) = (/ AChm_ , Str_ , AUp_ , Dn_ /)
   hash(:,200) = (/ AUp_ , Up_ , AUp_ , Up_ /)
   hash(:,201) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
   hash(:,202) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
   hash(:,203) = (/ AStr_ , Str_ , AStr_ , Str_ /)
   hash(:,204) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)

   outFound=.false.
   order(:)=-1

   do ih=1,hashSize
      if ( &
      .not.( &
      (idAPart(1).eq.0 .or. idAPart(1).eq.hash(1,ih)) &
      .and. &
      (idAPart(2).eq.0 .or. idAPart(2).eq.hash(2,ih)) &
      ) &
      ) cycle

      ! Final particles are q
      if (abs(idAPart(3)).lt.6 .and. abs(idAPart(4)).lt.6) then
         if ( &
         (idAPart(3).eq.0 .or. idAPart(3).eq.hash(3,ih)) &
         .and. &
         (idAPart(4).eq.0 .or. idAPart(4).eq.hash(4,ih)) &
         ) then
            order(1)=0
            order(2)=1
            outFound=.true.
         else if ( &
         (idAPart(3).eq.0 .or. idAPart(3).eq.hash(4,ih)) &
         .and. &
         (idAPart(4).eq.0 .or. idAPart(4).eq.hash(3,ih)) &
         ) then
            order(1)=1
            order(2)=0
            outFound=.true.
         endif
      ! Final particles l/nu
      else if ((IsALepton(idAPart(3)) .or. IsANeutrino(idAPart(3))) .and. (IsALepton(idAPart(4)) .or. IsANeutrino(idAPart(4)))) then
         if (abs(hash(1,ih)).eq.abs(hash(2,ih)) .and. abs(hash(1,ih)).eq.abs(hash(3,ih)) .and. abs(hash(1,ih)).eq.abs(hash(4,ih))) cycle ! Do not consider the ordering in uquq_uquq or dqdq_dqdq

         if ( &
         ( &
         sign(1, idAPart(3)).eq.sign(1, hash(3,ih)) .and. &
         ((IsALepton(idAPart(3)) .and. IsDownTypeQuark(hash(3,ih))) .or. (IsANeutrino(idAPart(3)) .and. IsUpTypeQuark(hash(3,ih)))) &
         ) &
         .and. &
         ( &
         sign(1, idAPart(4)).eq.sign(1, hash(4,ih)) .and. &
         ((IsALepton(idAPart(4)) .and. IsDownTypeQuark(hash(4,ih))) .or. (IsANeutrino(idAPart(4)) .and. IsUpTypeQuark(hash(4,ih)))) &
         ) &
         ) then
            order(1)=0
            order(2)=1
            outFound=.true.
         else if ( &
         ( &
         sign(1, idAPart(3)).eq.sign(1, hash(4,ih)) .and. &
         ((IsALepton(idAPart(3)) .and. IsDownTypeQuark(hash(4,ih))) .or. (IsANeutrino(idAPart(3)) .and. IsUpTypeQuark(hash(4,ih)))) &
         ) &
         .and. &
         ( &
         sign(1, idAPart(4)).eq.sign(1, hash(3,ih)) .and. &
         ((IsALepton(idAPart(4)) .and. IsDownTypeQuark(hash(3,ih))) .or. (IsANeutrino(idAPart(4)) .and. IsUpTypeQuark(hash(3,ih)))) &
         ) &
         ) then
            order(1)=1
            order(2)=0
            outFound=.true.
         endif
         outFound = ( CoupledVertex(idAPart(3:4), -1).eq.CoupledVertex(hash(3:4,ih),-1) )
      endif
      if (outFound) then
         exit
      endif
   enddo

end subroutine

subroutine EvalAmp_qqVVqq(idin, pin, ZWcode, msq)
use ModParameters
use ModMisc
implicit none
integer, intent(in) :: idin(1:mxpart)
real(8), intent(in) :: pin(1:mxpart,1:4)
real(8)             :: pin_MCFMconv(1:mxpart,1:4)
integer, intent(in) :: ZWcode
integer :: id_MCFM(1:mxpart)
real(8) :: p_MCFM(1:mxpart,1:4)
real(8) :: msq(-5:5,-5:5),msq_tmp(-5:5,-5:5)
integer, parameter :: doZZ=1,doWW=2,doZZorWW=3
logical :: doCompute
integer :: i,j

   msq(:,:)=0d0
   msq_tmp(:,:)=0d0

   pin_MCFMconv(:,:)=pin(:,:)/GeV

   doCompute = Setup_MCFM_qqVVqq(idin,pin_MCFMconv,id_MCFM,p_MCFM)
   if (doCompute) then
      call SetupParticleLabels(id_MCFM) ! Assign plabels 
      if(ZWcode.eq.doZZ) then
         if (Process.ge.66 .and. Process.le.68) then
            call qq_zzqq(p_MCFM,msq)
            if (id_MCFM(7).eq.0 .and. id_MCFM(8).eq.0) then ! Calculate for swapped momentum combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap just the momenta
               call qq_zzqq(p_MCFM,msq_tmp)
               msq = msq + msq_tmp
               do i=-5,5
                  msq(i,i)=msq(i,i)*0.5d0
               enddo
            else if (id_MCFM(7).eq.0 .or. id_MCFM(8).eq.0) then ! Calculate for wrong combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap the momenta
               call swap(id_MCFM(7),id_MCFM(8)) ! Swap the ids
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_zzqq(p_MCFM,msq_tmp)
               do i=-5,5
               do j=-5,5
                  if (msq(i,j).eq.0d0) then ! Non-zero MEs in original configuration are the correct ones, just replace 0 MEs
                     msq(i,j)=msq_tmp(i,j)
                  endif
               enddo
               enddo
            endif
         else if (Process.eq.69) then ! Or some other number?
            call qq_zzqqstrong(p_MCFM,msq)
            if (id_MCFM(7).eq.0 .and. id_MCFM(8).eq.0) then ! Calculate for swapped momentum combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap just the momenta
               call qq_zzqqstrong(p_MCFM,msq_tmp)
               msq = msq + msq_tmp
               do i=-5,5
                  if (i.eq.0) cycle ! gg->qqb+qbq does not need to be multiplied by 1/2
                  msq(i,i)=msq(i,i)*0.5d0
               enddo
               ! Subtract qqb->gg, which was counted twice
               id_MCFM(7:8)=Glu_
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_zzqqstrong(p_MCFM,msq_tmp)
               msq = msq - msq_tmp
            else if (id_MCFM(7).eq.0 .or. id_MCFM(8).eq.0) then ! Calculate for wrong combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap the momenta
               call swap(id_MCFM(7),id_MCFM(8)) ! Swap the ids
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_zzqqstrong(p_MCFM,msq_tmp)
               do i=-5,5
               do j=-5,5
                  if (msq(i,j).eq.0d0) then ! Non-zero MEs in original configuration are the correct ones, just replace 0 MEs
                     msq(i,j)=msq_tmp(i,j)
                  endif
               enddo
               enddo
            endif
         endif
      else if(ZWcode.eq.doWW) then
         if (Process.ge.66 .and. Process.le.68) then
            call qq_wwqq(p_MCFM,msq)
            if (id_MCFM(7).eq.0 .and. id_MCFM(8).eq.0) then ! Calculate for swapped momentum combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap just the momenta
               call qq_wwqq(p_MCFM,msq_tmp)
               msq = msq + msq_tmp
               do i=-5,5
                  msq(i,i)=msq(i,i)*0.5d0
               enddo
            else if (id_MCFM(7).eq.0 .or. id_MCFM(8).eq.0) then ! Calculate for wrong combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap the momenta
               call swap(id_MCFM(7),id_MCFM(8)) ! Swap the ids
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_wwqq(p_MCFM,msq_tmp)
               do i=-5,5
               do j=-5,5
                  if (msq(i,j).eq.0d0) then ! Non-zero MEs in original configuration are the correct ones, just replace 0 MEs
                     msq(i,j)=msq_tmp(i,j)
                  endif
               enddo
               enddo
            endif
         else if (Process.eq.69) then ! Or some other number?
            call qq_wwqqstrong(p_MCFM,msq)
            if (id_MCFM(7).eq.0 .and. id_MCFM(8).eq.0) then ! Calculate for swapped momentum combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap just the momenta
               call qq_wwqqstrong(p_MCFM,msq_tmp)
               msq = msq + msq_tmp
               do i=-5,5
                  if (i.eq.0) cycle ! gg->qqb+qbq does not need to be multiplied by 1/2
                  msq(i,i)=msq(i,i)*0.5d0
               enddo
               ! Subtract qqb->gg, which was counted twice
               id_MCFM(7:8)=Glu_
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_wwqqstrong(p_MCFM,msq_tmp)
               msq = msq - msq_tmp
            else if (id_MCFM(7).eq.0 .or. id_MCFM(8).eq.0) then ! Calculate for wrong combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap the momenta
               call swap(id_MCFM(7),id_MCFM(8)) ! Swap the ids
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_wwqqstrong(p_MCFM,msq_tmp)
               do i=-5,5
               do j=-5,5
                  if (msq(i,j).eq.0d0) then ! Non-zero MEs in original configuration are the correct ones, just replace 0 MEs
                     msq(i,j)=msq_tmp(i,j)
                  endif
               enddo
               enddo
            endif
         endif
      else if(ZWcode.eq.doZZorWW) then
         if (Process.ge.66 .and. Process.le.68) then
            call qq_vvqq(p_MCFM,msq)
            if (id_MCFM(7).eq.0 .and. id_MCFM(8).eq.0) then ! Calculate for swapped momentum combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap just the momenta
               call qq_vvqq(p_MCFM,msq_tmp)
               msq = msq + msq_tmp
               do i=-5,5
                  msq(i,i)=msq(i,i)*0.5d0
               enddo
            else if (id_MCFM(7).eq.0 .or. id_MCFM(8).eq.0) then ! Calculate for wrong combination
               call swap(p_MCFM(7,:),p_MCFM(8,:)) ! Swap the momenta
               call swap(id_MCFM(7),id_MCFM(8)) ! Swap the ids
               call SetupParticleLabels(id_MCFM) ! Reassign plabels
               call qq_vvqq(p_MCFM,msq_tmp)
               do i=-5,5
               do j=-5,5
                  if (msq(i,j).eq.0d0) then ! Non-zero MEs in original configuration are the correct ones, just replace 0 MEs
                     msq(i,j)=msq_tmp(i,j)
                  endif
               enddo
               enddo
            endif
         endif
      endif

   endif

end subroutine

END MODULE
