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
public :: ResetAmplitudeIncludes
public :: SetDistinguishWWCouplingsFlag
public :: SetSpinZeroGGCouplings,SetSpinZeroQQCouplings,SetSpinZeroVVCouplings
public :: SetSpinOneCouplings
public :: SetSpinTwoCouplings
public :: SetVprimeContactCouplings

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

subroutine ResetAmplitudeIncludes()
implicit none
   includeVprime=.false.
   includeGammaStar=.false.
end subroutine

subroutine SetSpinZeroVVCouplings(vvcoupl, vvpcoupl, vpvpcoupl, cqsq, Lambda_qsq, useWWcoupl)
   implicit none
   complex(8), intent(in) :: vvcoupl(39)
   complex(8), intent(in) :: vvpcoupl(39)
   complex(8), intent(in) :: vpvpcoupl(39)
   integer, intent(in) :: cqsq(3)
   real(8), intent(in) :: Lambda_qsq(1:3,1:4)
   logical, intent(in) :: useWWcoupl

   includeVprime = includeVprime .or. (                                   &
                    (any(vvpcoupl.ne.czero) .or. any(vpvpcoupl.ne.czero)) &
                   )

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


      ghzzp1 =  vvpcoupl(1)
      ghzzp2 =  vvpcoupl(2)
      ghzzp3 =  vvpcoupl(3)
      ghzzp4 =  vvpcoupl(4)

      ghzpgs2  = vvpcoupl(5)
      ghzpgs3  = vvpcoupl(6)
      ghzpgs4  = vvpcoupl(7)

      ghzzp1_prime = vvpcoupl(11)
      ghzzp1_prime2= vvpcoupl(12)
      ghzzp1_prime3= vvpcoupl(13)
      ghzzp1_prime4= vvpcoupl(14)
      ghzzp1_prime5= vvpcoupl(15)

      ghzzp2_prime = vvpcoupl(16)
      ghzzp2_prime2= vvpcoupl(17)
      ghzzp2_prime3= vvpcoupl(18)
      ghzzp2_prime4= vvpcoupl(19)
      ghzzp2_prime5= vvpcoupl(20)

      ghzzp3_prime = vvpcoupl(21)
      ghzzp3_prime2= vvpcoupl(22)
      ghzzp3_prime3= vvpcoupl(23)
      ghzzp3_prime4= vvpcoupl(24)
      ghzzp3_prime5= vvpcoupl(25)

      ghzzp4_prime = vvpcoupl(26)
      ghzzp4_prime2= vvpcoupl(27)
      ghzzp4_prime3= vvpcoupl(28)
      ghzzp4_prime4= vvpcoupl(29)
      ghzzp4_prime5= vvpcoupl(30)

      ghzpgs1_prime2= vvpcoupl(31)

      ghzzp1_prime6  = vvpcoupl(32)
      ghzzp1_prime7  = vvpcoupl(33)
      ghzzp2_prime6  = vvpcoupl(34)
      ghzzp2_prime7  = vvpcoupl(35)
      ghzzp3_prime6  = vvpcoupl(36)
      ghzzp3_prime7  = vvpcoupl(37)
      ghzzp4_prime6  = vvpcoupl(38)
      ghzzp4_prime7  = vvpcoupl(39)

      ghzpzp1 =  vpvpcoupl(1)
      ghzpzp2 =  vpvpcoupl(2)
      ghzpzp3 =  vpvpcoupl(3)
      ghzpzp4 =  vpvpcoupl(4)

      ghzpzp1_prime = vpvpcoupl(11)
      ghzpzp1_prime2= vpvpcoupl(12)
      ghzpzp1_prime3= vpvpcoupl(13)
      ghzpzp1_prime4= vpvpcoupl(14)
      ghzpzp1_prime5= vpvpcoupl(15)

      ghzpzp2_prime = vpvpcoupl(16)
      ghzpzp2_prime2= vpvpcoupl(17)
      ghzpzp2_prime3= vpvpcoupl(18)
      ghzpzp2_prime4= vpvpcoupl(19)
      ghzpzp2_prime5= vpvpcoupl(20)

      ghzpzp3_prime = vpvpcoupl(21)
      ghzpzp3_prime2= vpvpcoupl(22)
      ghzpzp3_prime3= vpvpcoupl(23)
      ghzpzp3_prime4= vpvpcoupl(24)
      ghzpzp3_prime5= vpvpcoupl(25)

      ghzpzp4_prime = vpvpcoupl(26)
      ghzpzp4_prime2= vpvpcoupl(27)
      ghzpzp4_prime3= vpvpcoupl(28)
      ghzpzp4_prime4= vpvpcoupl(29)
      ghzpzp4_prime5= vpvpcoupl(30)

      ghzpzp1_prime6  = vpvpcoupl(32)
      ghzpzp1_prime7  = vpvpcoupl(33)
      ghzpzp2_prime6  = vpvpcoupl(34)
      ghzpzp2_prime7  = vpvpcoupl(35)
      ghzpzp3_prime6  = vpvpcoupl(36)
      ghzpzp3_prime7  = vpvpcoupl(37)
      ghzpzp4_prime6  = vpvpcoupl(38)
      ghzpzp4_prime7  = vpvpcoupl(39)


      ! Set includeGammaStar based on the actual coupling values
      includeGammaStar = (          &
         ghzgs1_prime2.ne.czero .or. &
         ghzgs2.ne.czero .or.        &
         ghzgs3.ne.czero .or.        &
         ghzgs4.ne.czero .or.        &
         ghgsgs2.ne.czero .or.       &
         ghgsgs3.ne.czero .or.       &
         ghgsgs4.ne.czero .or.       &
         ghzpgs1_prime2.ne.czero .or. &
         ghzpgs2.ne.czero .or.        &
         ghzpgs3.ne.czero .or.        &
         ghzpgs4.ne.czero             &
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

      ghwwp1 =  vvpcoupl(1)
      ghwwp2 =  vvpcoupl(2)
      ghwwp3 =  vvpcoupl(3)
      ghwwp4 =  vvpcoupl(4)

      ghwwp1_prime = vvpcoupl(11)
      ghwwp1_prime2= vvpcoupl(12)
      ghwwp1_prime3= vvpcoupl(13)
      ghwwp1_prime4= vvpcoupl(14)
      ghwwp1_prime5= vvpcoupl(15)

      ghwwp2_prime = vvpcoupl(16)
      ghwwp2_prime2= vvpcoupl(17)
      ghwwp2_prime3= vvpcoupl(18)
      ghwwp2_prime4= vvpcoupl(19)
      ghwwp2_prime5= vvpcoupl(20)

      ghwwp3_prime = vvpcoupl(21)
      ghwwp3_prime2= vvpcoupl(22)
      ghwwp3_prime3= vvpcoupl(23)
      ghwwp3_prime4= vvpcoupl(24)
      ghwwp3_prime5= vvpcoupl(25)

      ghwwp4_prime = vvpcoupl(26)
      ghwwp4_prime2= vvpcoupl(27)
      ghwwp4_prime3= vvpcoupl(28)
      ghwwp4_prime4= vvpcoupl(29)
      ghwwp4_prime5= vvpcoupl(30)

      ghwwp1_prime6  = vvpcoupl(32)
      ghwwp1_prime7  = vvpcoupl(33)
      ghwwp2_prime6  = vvpcoupl(34)
      ghwwp2_prime7  = vvpcoupl(35)
      ghwwp3_prime6  = vvpcoupl(36)
      ghwwp3_prime7  = vvpcoupl(37)
      ghwwp4_prime6  = vvpcoupl(38)
      ghwwp4_prime7  = vvpcoupl(39)

      ghwpwp1 =  vpvpcoupl(1)
      ghwpwp2 =  vpvpcoupl(2)
      ghwpwp3 =  vpvpcoupl(3)
      ghwpwp4 =  vpvpcoupl(4)

      ghwpwp1_prime = vpvpcoupl(11)
      ghwpwp1_prime2= vpvpcoupl(12)
      ghwpwp1_prime3= vpvpcoupl(13)
      ghwpwp1_prime4= vpvpcoupl(14)
      ghwpwp1_prime5= vpvpcoupl(15)

      ghwpwp2_prime = vpvpcoupl(16)
      ghwpwp2_prime2= vpvpcoupl(17)
      ghwpwp2_prime3= vpvpcoupl(18)
      ghwpwp2_prime4= vpvpcoupl(19)
      ghwpwp2_prime5= vpvpcoupl(20)

      ghwpwp3_prime = vpvpcoupl(21)
      ghwpwp3_prime2= vpvpcoupl(22)
      ghwpwp3_prime3= vpvpcoupl(23)
      ghwpwp3_prime4= vpvpcoupl(24)
      ghwpwp3_prime5= vpvpcoupl(25)

      ghwpwp4_prime = vpvpcoupl(26)
      ghwpwp4_prime2= vpvpcoupl(27)
      ghwpwp4_prime3= vpvpcoupl(28)
      ghwpwp4_prime4= vpvpcoupl(29)
      ghwpwp4_prime5= vpvpcoupl(30)

      ghwpwp1_prime6  = vpvpcoupl(32)
      ghwpwp1_prime7  = vpvpcoupl(33)
      ghwpwp2_prime6  = vpvpcoupl(34)
      ghwpwp2_prime7  = vpvpcoupl(35)
      ghwpwp3_prime6  = vpvpcoupl(36)
      ghwpwp3_prime7  = vpvpcoupl(37)
      ghwpwp4_prime6  = vpvpcoupl(38)
      ghwpwp4_prime7  = vpvpcoupl(39)
   endif
   return
end subroutine SetSpinZeroVVCouplings

subroutine SetVprimeContactCouplings(Zpffcoupl, Wpffcoupl)
   implicit none
   complex(8), intent(in) :: Zpffcoupl(20)
   complex(8), intent(in) :: Wpffcoupl(20)

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

subroutine SetSpinTwoCouplings(acoupl,vvcoupl,vvpcoupl,vpvpcoupl,qLR)
   implicit none
   integer, parameter :: indexGammaBegin=11
   integer, parameter :: indexVVSize=20
   complex(8), intent(in) :: acoupl(1:5)
   complex(8), intent(in) :: vvcoupl(1:indexVVSize),vvpcoupl(1:indexVVSize),vpvpcoupl(1:indexVVSize)
   complex(8), intent(in) :: qLR(1:2)

   includeVprime = (any(vvpcoupl.ne.czero) .or. any(vpvpcoupl.ne.czero))
   includeGammaStar = (                                                       &
                    any(vvcoupl(indexGammaBegin:indexVVSize).ne.czero) .or.   &
                    any(vvpcoupl(indexGammaBegin:indexVVSize).ne.czero)       &
                   )

   a1 = acoupl(1)
   a2 = acoupl(2)
   a3 = acoupl(3)
   a4 = acoupl(4)
   a5 = acoupl(5)

   graviton_qq_left  = qLR(1)
   graviton_qq_right = qLR(2)

   b1 = vvcoupl(1)
   b2 = vvcoupl(2)
   b3 = vvcoupl(3)
   b4 = vvcoupl(4)
   b5 = vvcoupl(5)
   b6 = vvcoupl(6)
   b7 = vvcoupl(7)
   b8 = vvcoupl(8)
   b9 = vvcoupl(9)
   b10 = vvcoupl(10)

   bzgs1 = vvcoupl(11)
   bzgs2 = vvcoupl(12)
   bzgs3 = vvcoupl(13)
   bzgs4 = vvcoupl(14)
   bzgs8 = vvcoupl(15)

   bgsgs1 = vvcoupl(16)
   bgsgs2 = vvcoupl(17)
   bgsgs3 = vvcoupl(18)
   bgsgs4 = vvcoupl(19)
   bgsgs8 = vvcoupl(20)

   bzzp1 = vvpcoupl(1)
   bzzp2 = vvpcoupl(2)
   bzzp3 = vvpcoupl(3)
   bzzp4 = vvpcoupl(4)
   bzzp5 = vvpcoupl(5)
   bzzp6 = vvpcoupl(6)
   bzzp7 = vvpcoupl(7)
   bzzp8 = vvpcoupl(8)
   bzzp9 = vvpcoupl(9)
   bzzp10 = vvpcoupl(10)

   bzpgs1 = vvpcoupl(11)
   bzpgs2 = vvpcoupl(12)
   bzpgs3 = vvpcoupl(13)
   bzpgs4 = vvpcoupl(14)
   bzpgs8 = vvpcoupl(15)


   bzpzp1 = vpvpcoupl(1)
   bzpzp2 = vpvpcoupl(2)
   bzpzp3 = vpvpcoupl(3)
   bzpzp4 = vpvpcoupl(4)
   bzpzp5 = vpvpcoupl(5)
   bzpzp6 = vpvpcoupl(6)
   bzpzp7 = vpvpcoupl(7)
   bzpzp8 = vpvpcoupl(8)
   bzpzp9 = vpvpcoupl(9)
   bzpzp10 = vpvpcoupl(10)

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
    Ga_Vprime= 0d0
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
    Ga_Vprime= 0d0
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
