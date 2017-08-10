MODULE ModJHUGenMELA
   use ModParameters
   use ModKinematics
   use ModMisc
#if compiler==1
   use ifport
#endif
   implicit none

! JHUGenMELA-specific subroutines and functions

public :: SetEWParameters
public :: SetHiggsMassWidth
public :: SetDecayModes
public :: SetTopDecays
public :: SetHDK
public :: SetMuRenFac
public :: ResetMubarHGabarH
public :: SetSpinZeroVVCouplings,SetSpinZeroVVCouplings_NoGamma,SetDistinguishWWCouplingsFlag, SetSpinZeroContactTerms
public :: SetSpinZeroGGCouplings,SetSpinZeroQQCouplings
public :: SetSpinOneCouplings
public :: SetSpinTwoCouplings

public :: GetMVGV,GetMVPrimeGVPrime
public :: GetAlphaSAlphaSMZ
public :: GetPDFConstants
public :: GetDecayCouplings


contains


!=====================================================
! Subroutines visible to the user
!=====================================================

subroutine SetEWParameters(inMZ, inMW, inGf, inalpha_QED, inxw)
   implicit none
   real(8), intent(in) :: inMZ, inMW, inGf, inalpha_QED, inxw
   M_Z = inMZ
   M_W = inMW
   Gf = inGf
   alpha_QED = inalpha_QED
   xw = inxw
   call ComputeEWVariables()
end subroutine SetEWParameters


subroutine SetHiggsMassWidth(mass,width)
   implicit none
   real(8), intent(in) :: mass, width
   call SetMass(mass,Hig_)
   call SetDecayWidth(width,Hig_)
   return
end subroutine SetHiggsMassWidth

subroutine SetDecayModes(idfirst,idsecond)
   implicit none
   integer, intent(in) :: idfirst(1:2)
   integer, intent(in) :: idsecond(1:2)
   integer :: idV(1:2)

   if( idfirst(1).eq.Pho_ .or. idfirst(2).eq.Pho_ ) then
      DecayMode1 = 7
      idV(1)=Pho_
   else
      idV(1)=CoupledVertex(idfirst,-1)
      if( idV(1).eq.Wp_ .or. idV(1).eq.Wm_ ) then
         DecayMode1=11
      elseif( idV(1).eq.Z0_ ) then
         DecayMode1=9
      endif
   endif
   if( idsecond(1).eq.Pho_ .or. idsecond(2).eq.Pho_ ) then
      DecayMode2 = 7
      idV(2)=Pho_
   else
      idV(2)=CoupledVertex(idsecond,-1)
      if( idV(2).eq.Wp_ .or. idV(2).eq.Wm_ ) then
         DecayMode2=11
      elseif( idV(2).eq.Z0_ ) then
         DecayMode2=9
      endif
   endif
   call SetMVGV()
   return
end subroutine SetDecayModes

subroutine SetTopDecays(flag)
implicit none
integer, intent(in) :: flag
   TopDecays=flag
end subroutine

subroutine SetHDK(flag)
implicit none
logical, intent(in) :: flag
   H_DK=flag
end subroutine

subroutine SetMuRenFac(muren,mufac)
implicit none
real(8), intent(in) :: muren,mufac
	Mu_Ren = muren
   Mu_Fact = mufac
end subroutine

subroutine ResetMubarHGabarH()
implicit none
   mubarH=-999d0
   gabarH=-999d0
end subroutine

subroutine SetSpinZeroVVCouplings(vvcoupl, cqsq, Lambda_qsq, useWWcoupl)
   implicit none
   complex(8), intent(in) :: vvcoupl(39)
   integer, intent(in) :: cqsq(3)
   real(8), intent(in) :: Lambda_qsq(1:3,1:4)
   logical, intent(in) :: useWWcoupl

   if(.not.useWWcoupl) then
      cz_q1sq = cqsq(1)
      Lambda_z11 = Lambda_qsq(1,1)
      Lambda_z21 = Lambda_qsq(1,2)
      Lambda_z31 = Lambda_qsq(1,3)
      Lambda_z41 = Lambda_qsq(1,4)
      cz_q2sq = cqsq(2)
      Lambda_z12 = Lambda_qsq(2,1)
      Lambda_z22 = Lambda_qsq(2,2)
      Lambda_z32 = Lambda_qsq(2,3)
      Lambda_z42 = Lambda_qsq(2,4)
      cz_q12sq = cqsq(3)
      Lambda_z10 = Lambda_qsq(3,1)
      Lambda_z20 = Lambda_qsq(3,2)
      Lambda_z30 = Lambda_qsq(3,3)
      Lambda_z40 = Lambda_qsq(3,4)

      ghz1 =  vvcoupl(1)
      ghz2 =  vvcoupl(2)
      ghz3 =  vvcoupl(3)
      ghz4 =  vvcoupl(4)

      ghzgs2  = vvcoupl(5)
      ghzgs3  = vvcoupl(6)
      ghzgs4  = vvcoupl(7)
      ghgsgs2 = vvcoupl(8)
      ghgsgs3 = vvcoupl(9)
      ghgsgs4 = vvcoupl(10)

      ghz1_prime = vvcoupl(11)
      ghz1_prime2= vvcoupl(12)
      ghz1_prime3= vvcoupl(13)
      ghz1_prime4= vvcoupl(14)
      ghz1_prime5= vvcoupl(15)

      ghz2_prime = vvcoupl(16)
      ghz2_prime2= vvcoupl(17)
      ghz2_prime3= vvcoupl(18)
      ghz2_prime4= vvcoupl(19)
      ghz2_prime5= vvcoupl(20)

      ghz3_prime = vvcoupl(21)
      ghz3_prime2= vvcoupl(22)
      ghz3_prime3= vvcoupl(23)
      ghz3_prime4= vvcoupl(24)
      ghz3_prime5= vvcoupl(25)

      ghz4_prime = vvcoupl(26)
      ghz4_prime2= vvcoupl(27)
      ghz4_prime3= vvcoupl(28)
      ghz4_prime4= vvcoupl(29)
      ghz4_prime5= vvcoupl(30)

      ghzgs1_prime2= vvcoupl(31)

      ghz1_prime6  = vvcoupl(32)
      ghz1_prime7  = vvcoupl(33)
      ghz2_prime6  = vvcoupl(34)
      ghz2_prime7  = vvcoupl(35)
      ghz3_prime6  = vvcoupl(36)
      ghz3_prime7  = vvcoupl(37)
      ghz4_prime6  = vvcoupl(38)
      ghz4_prime7  = vvcoupl(39)

      ! Set includeGammaStar based on the actual coupling values
      includeGammaStar = (          &
         ghzgs1_prime2.ne.czero .or. &
         ghzgs2.ne.czero .or.        &
         ghzgs3.ne.czero .or.        &
         ghzgs4.ne.czero .or.        &
         ghgsgs2.ne.czero .or.       &
         ghgsgs3.ne.czero .or.       &
         ghgsgs4.ne.czero            &
      )

   else
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
end subroutine SetSpinZeroVVCouplings

subroutine SetSpinZeroVVCouplings_NoGamma(vvcoupl, cqsq, Lambda_qsq, useWWcoupl)
   implicit none
   complex(8), intent(in) :: vvcoupl(32)
   integer, intent(in) :: cqsq(3)
   real(8), intent(in) :: Lambda_qsq(1:3,1:4)
   logical, intent(in) :: useWWcoupl

   if(.not.useWWcoupl) then
      cz_q1sq = cqsq(1)
      Lambda_z11 = Lambda_qsq(1,1)
      Lambda_z21 = Lambda_qsq(1,2)
      Lambda_z31 = Lambda_qsq(1,3)
      Lambda_z41 = Lambda_qsq(1,4)
      cz_q2sq = cqsq(2)
      Lambda_z12 = Lambda_qsq(2,1)
      Lambda_z22 = Lambda_qsq(2,2)
      Lambda_z32 = Lambda_qsq(2,3)
      Lambda_z42 = Lambda_qsq(2,4)
      cz_q12sq = cqsq(3)
      Lambda_z10 = Lambda_qsq(3,1)
      Lambda_z20 = Lambda_qsq(3,2)
      Lambda_z30 = Lambda_qsq(3,3)
      Lambda_z40 = Lambda_qsq(3,4)

      ghz1 = vvcoupl(1)
      ghz2 = vvcoupl(2)
      ghz3 = vvcoupl(3)
      ghz4 = vvcoupl(4)
      ghz1_prime = vvcoupl(5)
      ghz1_prime2= vvcoupl(6)
      ghz1_prime3= vvcoupl(7)
      ghz1_prime4= vvcoupl(8)
      ghz1_prime5= vvcoupl(9)
      ghz2_prime = vvcoupl(10)
      ghz2_prime2= vvcoupl(11)
      ghz2_prime3= vvcoupl(12)
      ghz2_prime4= vvcoupl(13)
      ghz2_prime5= vvcoupl(14)
      ghz3_prime = vvcoupl(15)
      ghz3_prime2= vvcoupl(16)
      ghz3_prime3= vvcoupl(17)
      ghz3_prime4= vvcoupl(18)
      ghz3_prime5= vvcoupl(19)
      ghz4_prime = vvcoupl(20)
      ghz4_prime2= vvcoupl(21)
      ghz4_prime3= vvcoupl(22)
      ghz4_prime4= vvcoupl(23)
      ghz4_prime5= vvcoupl(24)
      ghz1_prime6= vvcoupl(25)
      ghz1_prime7= vvcoupl(26)
      ghz2_prime6= vvcoupl(27)
      ghz2_prime7= vvcoupl(28)
      ghz3_prime6= vvcoupl(29)
      ghz3_prime7= vvcoupl(30)
      ghz4_prime6= vvcoupl(31)
      ghz4_prime7= vvcoupl(32)
   else
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

      ghw1 = vvcoupl(1)
      ghw2 = vvcoupl(2)
      ghw3 = vvcoupl(3)
      ghw4 = vvcoupl(4)
      ghw1_prime = vvcoupl(5)
      ghw1_prime2= vvcoupl(6)
      ghw1_prime3= vvcoupl(7)
      ghw1_prime4= vvcoupl(8)
      ghw1_prime5= vvcoupl(9)
      ghw2_prime = vvcoupl(10)
      ghw2_prime2= vvcoupl(11)
      ghw2_prime3= vvcoupl(12)
      ghw2_prime4= vvcoupl(13)
      ghw2_prime5= vvcoupl(14)
      ghw3_prime = vvcoupl(15)
      ghw3_prime2= vvcoupl(16)
      ghw3_prime3= vvcoupl(17)
      ghw3_prime4= vvcoupl(18)
      ghw3_prime5= vvcoupl(19)
      ghw4_prime = vvcoupl(20)
      ghw4_prime2= vvcoupl(21)
      ghw4_prime3= vvcoupl(22)
      ghw4_prime4= vvcoupl(23)
      ghw4_prime5= vvcoupl(24)
      ghw1_prime6= vvcoupl(25)
      ghw1_prime7= vvcoupl(26)
      ghw2_prime6= vvcoupl(27)
      ghw2_prime7= vvcoupl(28)
      ghw3_prime6= vvcoupl(29)
      ghw3_prime7= vvcoupl(30)
      ghw4_prime6= vvcoupl(31)
      ghw4_prime7= vvcoupl(32)
   endif
   return
end subroutine SetSpinZeroVVCouplings_NoGamma

subroutine SetSpinZeroContactTerms(zzpcoupl, zpzpcoupl, Zpffcoupl, wwpcoupl, wpwpcoupl, Wpffcoupl, UseVp)
  implicit none
  complex(8), intent(in) :: zzpcoupl(39)
  complex(8), intent(in) :: zpzpcoupl(39)
  complex(8), intent(in) :: Zpffcoupl(20)

  complex(8), intent(in) :: wwpcoupl(39)
  complex(8), intent(in) :: wpwpcoupl(39)
  complex(8), intent(in) :: Wpffcoupl(20)

  logical, intent(in) :: UseVp

  includeVprime = ((any(zzpcoupl.ne.czero) .or. any(zpzpcoupl.ne.czero)) .and. any(Zpffcoupl.ne.czero)) &
              .or.((any(wwpcoupl.ne.czero) .or. any(wpwpcoupl.ne.czero)) .and. any(Wpffcoupl.ne.czero))

  ghzzp1 = zzpcoupl(1)
  ghzpzp1 = zpzpcoupl(1)

  ezp_El_left = Zpffcoupl(1)
  ezp_El_right = Zpffcoupl(2)
  ezp_Mu_left = Zpffcoupl(3)
  ezp_Mu_right = Zpffcoupl(4)
  ezp_Ta_left = Zpffcoupl(5)
  ezp_Ta_right = Zpffcoupl(6)
  ezp_NuE_left = Zpffcoupl(7)
  ezp_NuE_right = Zpffcoupl(8)
  ezp_Dn_left = Zpffcoupl(9)
  ezp_Dn_right = Zpffcoupl(10)
  ezp_Up_left = Zpffcoupl(11)
  ezp_Up_right = Zpffcoupl(12)
  ezp_Str_left = Zpffcoupl(13)
  ezp_Str_right = Zpffcoupl(14)
  ezp_Chm_left = Zpffcoupl(15)
  ezp_Chm_right = Zpffcoupl(16)
  ezp_Bot_left = Zpffcoupl(17)
  ezp_Bot_right = Zpffcoupl(18)
  ezp_Top_left = Zpffcoupl(19)
  ezp_Top_right = Zpffcoupl(20)

  ghwwp1 = wwpcoupl(1)
  ghwpwp1 = wpwpcoupl(1)

  ewp_El_left = Wpffcoupl(1)
  ewp_El_right = Wpffcoupl(2)
  ewp_Mu_left = Wpffcoupl(3)
  ewp_Mu_right = Wpffcoupl(4)
  ewp_Ta_left = Wpffcoupl(5)
  ewp_Ta_right = Wpffcoupl(6)
  ewp_Up_left = Wpffcoupl(11)
  ewp_Up_right = Wpffcoupl(12)
  ewp_Chm_left = Wpffcoupl(15)
  ewp_Chm_right = Wpffcoupl(16)
  ewp_Top_left = Wpffcoupl(19)
  ewp_Top_right = Wpffcoupl(20)

end subroutine

subroutine SetDistinguishWWCouplingsFlag(doAllow)
   implicit none
   logical, intent(in) :: doAllow
   distinguish_HWWcouplings = doAllow
   return
end subroutine SetDistinguishWWCouplingsFlag

subroutine SetSpinZeroGGCouplings(ggcoupl)
   implicit none
   complex(8), intent(in) :: ggcoupl(1:3)
   ghg2 =  ggcoupl(1)
   ghg3 =  ggcoupl(2)
   ghg4 =  ggcoupl(3)
   return
end subroutine SetSpinZeroGGCouplings

subroutine SetSpinZeroQQCouplings(qqcoupl)
   implicit none
   complex(8), intent(in) :: qqcoupl(1:2)
   kappa =  qqcoupl(1)
   kappa_tilde =  qqcoupl(2)
   return
end subroutine SetSpinZeroQQCouplings

subroutine SetSpinOneCouplings(qqcoupl,vvcoupl)
   implicit none
   complex(8), intent(in) :: qqcoupl(1:2)
   complex(8), intent(in) :: vvcoupl(1:2)

   zprime_qq_left  = qqcoupl(1)
   zprime_qq_right = qqcoupl(2)
   zprime_zz_1 = vvcoupl(1)
   zprime_zz_2 = vvcoupl(2)
   return
end subroutine SetSpinOneCouplings

subroutine SetSpinTwoCouplings(acoupl,bcoupl,qLR)
   implicit none
   complex(8), intent(in) :: acoupl(1:5)
   complex(8), intent(in) :: bcoupl(1:10)
   complex(8), intent(in) :: qLR(1:2)

   a1 = acoupl(1)
   a2 = acoupl(2)
   a3 = acoupl(3)
   a4 = acoupl(4)
   a5 = acoupl(5)

   graviton_qq_left  = qLR(1)
   graviton_qq_right = qLR(2)

   b1 = bcoupl(1)
   b2 = bcoupl(2)
   b3 = bcoupl(3)
   b4 = bcoupl(4)
   b5 = bcoupl(5)
   b6 = bcoupl(6)
   b7 = bcoupl(7)
   b8 = bcoupl(8)
   b9 = bcoupl(9)
   b10 = bcoupl(10)

   return
end subroutine SetSpinTwoCouplings

subroutine SetMVGV()
implicit none
  if( IsAZDecay(DecayMode1) ) then
    M_V = M_Z
    Ga_V= Ga_Z
    M_Vprime = M_Zprime
    Ga_Vprime= Ga_Zprime
  elseif( IsAWDecay(DecayMode1) ) then
    M_V = M_W
    Ga_V= Ga_W
    M_Vprime = M_Wprime
    Ga_Vprime= Ga_Wprime
  else
    M_V = 0d0
    Ga_V= 0d0
    M_Vprime = -1d0
    Ga_Vprime= -1d0
  endif
end subroutine SetMVGV

subroutine SetMVGVFromVertex(idV)
implicit none
integer idV
  if( idV.eq.Z0_ ) then
    M_V = M_Z
    Ga_V= Ga_Z
    M_Vprime = M_Zprime
    Ga_Vprime= Ga_Zprime
  else if( idV.eq.Wp_ .or. idV.eq.Wm_ ) then
    M_V = M_W
    Ga_V= Ga_W
    M_Vprime = M_Wprime
    Ga_Vprime= Ga_Wprime
  else
    M_V = 0d0
    Ga_V= 0d0
    M_Vprime = -1d0
    Ga_Vprime= -1d0
  endif
end subroutine

subroutine GetMVGV(mv,gv)
   implicit none
   real(8), intent(out) :: mv,gv
   mv=M_V
   gv=Ga_V
end subroutine GetMVGV

subroutine GetMVPrimeGVPrime(mv,gv)
   implicit none
   real(8), intent(out) :: mv,gv
   mv=M_Vprime
   gv=Ga_Vprime
end subroutine GetMVPrimeGVPrime

subroutine GetAlphaSAlphaSMZ(val_as, val_asmz)
implicit none
real(8), intent(out) :: val_as, val_asmz
   val_as=alphas
   val_asmz=alphas_mz
end subroutine

subroutine GetPDFConstants(pdfzmass, pdfnloops, pdfnf)
implicit none
real(8), intent(out) :: pdfzmass
integer, intent(out) :: pdfnloops, pdfnf
   pdfzmass=zmass_pdf
   pdfnloops=nloops_pdf
   pdfnf=nQflavors_pdf
end subroutine GetPDFConstants


! This subroutine is slightly different form the one in the decay MEs in the sense that onshell photon returns 0,0 instead of 1,1
subroutine GetDecayCouplings(VVMode,idordered,aL1,aR1,aL2,aR2)
   implicit none
   integer, intent(in) :: VVMode,idordered(6:9)
   real(dp), intent(out) :: aL1,aR1,aL2,aR2

   !        h3/h4 helicities: -1 == left, 1 == right
   if( VVMode.eq.ZZMode ) then
   !        ZZ DECAYS
      if( abs(idordered(6)).eq.abs(ElM_) .or. abs(idordered(6)).eq.abs(MuM_)  ) then
         aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
         aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(6)).eq.abs(TaM_) ) then
         aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
         aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(6)).eq.abs(NuE_) .or. abs(idordered(6)).eq.abs(NuM_) .or. abs(idordered(6)).eq.abs(NuT_) ) then
         aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
         aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
      elseif( abs(idordered(6)).eq.abs(Up_) .or. abs(idordered(6)).eq.abs(Chm_) ) then
         aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(6)).eq.abs(Dn_) .or. abs(idordered(6)).eq.abs(Str_) .or. abs(idordered(6)).eq.abs(Bot_) ) then
         aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL1=0d0
         aR1=0d0
      endif
      if( abs(idordered(8)).eq.abs(ElM_) .or. abs(idordered(8)).eq.abs(MuM_)  ) then
         aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
         aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(8)).eq.abs(TaM_) ) then
         aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
         aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(8)).eq.abs(NuE_) .or. abs(idordered(8)).eq.abs(NuM_) .or. abs(idordered(8)).eq.abs(NuT_) ) then
         aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
         aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
      elseif( abs(idordered(8)).eq.abs(Up_) .or. abs(idordered(8)).eq.abs(Chm_) ) then
         aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(8)).eq.abs(Dn_) .or. abs(idordered(8)).eq.abs(Str_) .or. abs(idordered(8)).eq.abs(Bot_) ) then
         aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL2=0d0
         aR2=0d0
      endif

   elseif( VVMode.eq.WWMode ) then
   !        WW DECAYS
      if( IsAQuark(idordered(6)) ) then
         aL1 = bL * dsqrt(scale_alpha_W_ud)
         aR1 = bR * dsqrt(scale_alpha_W_ud)! = 0
      elseif( &
               (abs(idordered(6)).eq.abs(ElP_) .and. abs(idordered(7)).eq.abs(NuE_)) .or. (abs(idordered(7)).eq.abs(ElP_) .and. abs(idordered(6)).eq.abs(NuE_)) .or. &
               (abs(idordered(6)).eq.abs(MuP_) .and. abs(idordered(7)).eq.abs(NuM_)) .or. (abs(idordered(7)).eq.abs(MuP_) .and. abs(idordered(6)).eq.abs(NuM_))      &
            ) then
         aL1 = bL * dsqrt(scale_alpha_W_ln)
         aR1 = bR * dsqrt(scale_alpha_W_ln)! = 0
      elseif( &
               (abs(idordered(6)).eq.abs(TaP_) .and. abs(idordered(7)).eq.abs(NuT_)) .or. (abs(idordered(7)).eq.abs(TaP_) .and. abs(idordered(6)).eq.abs(NuT_))      &
            ) then
         aL1 = bL * dsqrt(scale_alpha_W_tn)
         aR1 = bR * dsqrt(scale_alpha_W_tn)! = 0
      else
         aL1=0d0
         aR1=0d0
      endif
      if( IsAQuark(idordered(8)) ) then
         aL2 = bL * dsqrt(scale_alpha_W_ud)
         aR2 = bR * dsqrt(scale_alpha_W_ud)! = 0
      elseif( &
               (abs(idordered(8)).eq.abs(ElM_) .and. abs(idordered(9)).eq.abs(ANuE_)) .or. (abs(idordered(9)).eq.abs(ElM_) .and. abs(idordered(8)).eq.abs(ANuE_)) .or. &
               (abs(idordered(8)).eq.abs(MuM_) .and. abs(idordered(9)).eq.abs(ANuM_)) .or. (abs(idordered(9)).eq.abs(MuM_) .and. abs(idordered(8)).eq.abs(ANuM_))      &
            ) then
         aL2 = bL * dsqrt(scale_alpha_W_ln)
         aR2 = bR * dsqrt(scale_alpha_W_ln)! = 0
      elseif( &
               (abs(idordered(8)).eq.abs(TaM_) .and. abs(idordered(9)).eq.abs(ANuT_)) .or. (abs(idordered(9)).eq.abs(TaM_) .and. abs(idordered(8)).eq.abs(ANuT_))      &
            ) then
         aL2 = bL * dsqrt(scale_alpha_W_tn)
         aR2 = bR * dsqrt(scale_alpha_W_tn)! = 0
      else
         aL2=0d0
         aR2=0d0
      endif

   elseif( VVMode.eq.ZgMode ) then
   !        Zgamma DECAYS
      if( abs(idordered(6)).eq.abs(ElM_) .or. abs(idordered(6)).eq.abs(MuM_) ) then
         aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
         aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(6)).eq.abs(TaM_) ) then
         aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
         aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(6)).eq.abs(NuE_) .or. abs(idordered(6)).eq.abs(NuM_) .or. abs(idordered(6)).eq.abs(NuT_) ) then
         aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
         aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
      elseif( abs(idordered(6)).eq.abs(Up_) .or. abs(idordered(6)).eq.abs(Chm_) ) then
         aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(6)).eq.abs(Dn_) .or. abs(idordered(6)).eq.abs(Str_) .or. abs(idordered(6)).eq.abs(Bot_) ) then
         aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL1=0d0
         aR1=0d0
      endif
      aL2=0d0
      aR2=0d0

   elseif( VVMode.eq.ggMode ) then
   !        gamma gamma DECAYS
      aL1=0d0
      aR1=0d0
      aL2=0d0
      aR2=0d0

   elseif( VVMode.eq.gsgMode ) then
   !        gamma* gamma DECAYS
      if( abs(idordered(6)).eq.abs(ElM_) .or. abs(idordered(6)).eq.abs(MuM_) ) then
         aL1=cL_lep    * dsqrt(scale_alpha_Z_ll)
         aR1=cR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(6)).eq.abs(TaM_) ) then
         aL1=cL_lep    * dsqrt(scale_alpha_Z_tt)
         aR1=cR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(6)).eq.abs(NuE_) .or. abs(idordered(6)).eq.abs(NuM_) .or. abs(idordered(6)).eq.abs(NuT_) ) then
         aL1=cL_neu    * dsqrt(scale_alpha_Z_nn)! = 0
         aR1=cR_neu    * dsqrt(scale_alpha_Z_nn)! = 0
      elseif( abs(idordered(6)).eq.abs(Up_) .or. abs(idordered(6)).eq.abs(Chm_) ) then
         aL1=cL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR1=cR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(6)).eq.abs(Dn_) .or. abs(idordered(6)).eq.abs(Str_) .or. abs(idordered(6)).eq.abs(Bot_) ) then
         aL1=cL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR1=cR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL1=0d0
         aR1=0d0
      endif
      aL2=0d0
      aR2=0d0

   elseif( VVMode.eq.gsZMode ) then
   !        gamma* Z DECAYS
      if( abs(idordered(6)).eq.abs(ElM_) .or. abs(idordered(6)).eq.abs(MuM_)  ) then
         aL1=cL_lep    * dsqrt(scale_alpha_Z_ll)
         aR1=cR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(6)).eq.abs(TaM_) ) then
         aL1=cL_lep    * dsqrt(scale_alpha_Z_tt)
         aR1=cR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(6)).eq.abs(NuE_) .or. abs(idordered(6)).eq.abs(NuM_) .or. abs(idordered(6)).eq.abs(NuT_) ) then
         aL1=cL_neu    * dsqrt(scale_alpha_Z_nn)
         aR1=cR_neu    * dsqrt(scale_alpha_Z_nn)
      elseif( abs(idordered(6)).eq.abs(Up_) .or. abs(idordered(6)).eq.abs(Chm_) ) then
         aL1=cL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR1=cR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(6)).eq.abs(Dn_) .or. abs(idordered(6)).eq.abs(Str_) .or. abs(idordered(6)).eq.abs(Bot_) ) then
         aL1=cL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR1=cR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL1=0d0
         aR1=0d0
      endif
      if( abs(idordered(8)).eq.abs(ElM_) .or. abs(idordered(8)).eq.abs(MuM_) ) then
         aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
         aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(8)).eq.abs(TaM_) ) then
         aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
         aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(8)).eq.abs(NuE_) .or. abs(idordered(8)).eq.abs(NuM_) .or. abs(idordered(8)).eq.abs(NuT_) ) then
         aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
         aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
      elseif( abs(idordered(8)).eq.abs(Up_) .or. abs(idordered(8)).eq.abs(Chm_) ) then
         aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(8)).eq.abs(Dn_) .or. abs(idordered(8)).eq.abs(Str_) .or. abs(idordered(8)).eq.abs(Bot_) ) then
         aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL2=0d0
         aR2=0d0
      endif

   elseif( VVMode.eq.ZgsMode ) then
   !        Z gamma* DECAYS
      if( abs(idordered(6)).eq.abs(ElM_) .or. abs(idordered(6)).eq.abs(MuM_) ) then
         aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
         aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(6)).eq.abs(TaM_) ) then
         aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
         aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(6)).eq.abs(NuE_) .or. abs(idordered(6)).eq.abs(NuM_) .or. abs(idordered(6)).eq.abs(NuT_) ) then
         aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
         aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
      elseif( abs(idordered(6)).eq.abs(Up_) .or. abs(idordered(6)).eq.abs(Chm_) ) then
         aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(6)).eq.abs(Dn_) .or. abs(idordered(6)).eq.abs(Str_) .or. abs(idordered(6)).eq.abs(Bot_) ) then
         aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL1=0d0
         aR1=0d0
      endif
      if( abs(idordered(8)).eq.abs(ElM_) .or. abs(idordered(8)).eq.abs(MuM_) ) then
         aL2=cL_lep    * dsqrt(scale_alpha_Z_ll)
         aR2=cR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(8)).eq.abs(TaM_) ) then
         aL2=cL_lep    * dsqrt(scale_alpha_Z_tt)
         aR2=cR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(8)).eq.abs(NuE_) .or. abs(idordered(8)).eq.abs(NuM_) .or. abs(idordered(8)).eq.abs(NuT_) ) then
         aL2=cL_neu    * dsqrt(scale_alpha_Z_nn)! = 0
         aR2=cR_neu    * dsqrt(scale_alpha_Z_nn)! = 0
      elseif( abs(idordered(8)).eq.abs(Up_) .or. abs(idordered(8)).eq.abs(Chm_) ) then
         aL2=cL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR2=cR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(8)).eq.abs(Dn_) .or. abs(idordered(8)).eq.abs(Str_) .or. abs(idordered(8)).eq.abs(Bot_) ) then
         aL2=cL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR2=cR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL2=0d0
         aR2=0d0
      endif

   elseif( VVMode.eq.gsgsMode ) then
   !        gamma* gamma* DECAYS
      if( abs(idordered(6)).eq.abs(ElM_) .or. abs(idordered(6)).eq.abs(MuM_) ) then
         aL1=cL_lep    * dsqrt(scale_alpha_Z_ll)
         aR1=cR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(6)).eq.abs(TaM_) ) then
         aL1=cL_lep    * dsqrt(scale_alpha_Z_tt)
         aR1=cR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(6)).eq.abs(NuE_) .or. abs(idordered(6)).eq.abs(NuM_) .or. abs(idordered(6)).eq.abs(NuT_) ) then
         aL1=cL_neu    * dsqrt(scale_alpha_Z_nn)! = 0
         aR1=cR_neu    * dsqrt(scale_alpha_Z_nn)! = 0
      elseif( abs(idordered(6)).eq.abs(Up_) .or. abs(idordered(6)).eq.abs(Chm_) ) then
         aL1=cL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR1=cR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(6)).eq.abs(Dn_) .or. abs(idordered(6)).eq.abs(Str_) .or. abs(idordered(6)).eq.abs(Bot_) ) then
         aL1=cL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR1=cR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL1=0d0
         aR1=0d0
      endif
      if( abs(idordered(8)).eq.abs(ElM_) .or. abs(idordered(8)).eq.abs(MuM_) ) then
         aL2=cL_lep    * dsqrt(scale_alpha_Z_ll)
         aR2=cR_lep    * dsqrt(scale_alpha_Z_ll)
      elseif( abs(idordered(8)).eq.abs(TaM_) ) then
         aL2=cL_lep    * dsqrt(scale_alpha_Z_tt)
         aR2=cR_lep    * dsqrt(scale_alpha_Z_tt)
      elseif( abs(idordered(8)).eq.abs(NuE_) .or. abs(idordered(8)).eq.abs(NuM_) .or. abs(idordered(8)).eq.abs(NuT_) ) then
         aL2=cL_neu    * dsqrt(scale_alpha_Z_nn)! = 0
         aR2=cR_neu    * dsqrt(scale_alpha_Z_nn)! = 0
      elseif( abs(idordered(8)).eq.abs(Up_) .or. abs(idordered(8)).eq.abs(Chm_) ) then
         aL2=cL_QUp    * dsqrt(scale_alpha_Z_uu)
         aR2=cR_QUp    * dsqrt(scale_alpha_Z_uu)
      elseif( abs(idordered(8)).eq.abs(Dn_) .or. abs(idordered(8)).eq.abs(Str_) .or. abs(idordered(8)).eq.abs(Bot_) ) then
         aL2=cL_QDn    * dsqrt(scale_alpha_Z_dd)
         aR2=cR_QDn    * dsqrt(scale_alpha_Z_dd)
      else
         aL2=0d0
         aR2=0d0
      endif

   else
      aL1=0d0
      aR1=0d0
      aL2=0d0
      aR2=0d0
   endif

   return
end subroutine


END MODULE ModJHUGenMELA
