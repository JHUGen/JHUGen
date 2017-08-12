module modHiggsJJ
  use modParameters
  use ModMisc
  implicit none
  private

  public :: EvalAmp_WBFH_UnSymm_SA,EvalAmp_WBFH_UnSymm_SA_Select,EvalAmp_WBFH_UnSymm_SA_Select_exact
  public :: EvalAmp_SBFH_UnSymm_SA,EvalAmp_SBFH_UnSymm_SA_Select,EvalAmp_SBFH_UnSymm_SA_Select_exact
  !public :: wrapHVV
  public :: get_VBFchannelHash,get_VBFchannelHash_nosplit,get_HJJchannelHash,get_HJJchannelHash_nosplit,get_GENchannelHash,get_VBFoffshchannelHash
  public :: init_VBFoffshChannelHash,get_NumberOfChannels,get_VBFoffshChannel,remove_VBFoffshChannel

  integer,private :: ij_channelHash(1:204,1:4)=-99
  integer,private :: NumberOfChannels=-1
  logical,private :: deterministicME=.false.

  !-- general definitions, to be merged with Markus final structure
   real(dp), public, parameter :: tag1 = 1.0_dp
   real(dp), public, parameter :: tag2 = 1.0_dp
   real(dp), public, parameter :: tagbot = 1.0_dp


 CONTAINS




  subroutine get_VBFchannelHash(ijSel)
  implicit none
  integer, intent(out) :: ijSel(1:121,1:3)
  integer,parameter :: zz=1, ww=0


      ijSel(  1,1:3) = (/ 2, 1, zz/)
      ijSel(  2,1:3) = (/ 2, 1, ww/)
      ijSel(  3,1:3) = (/ 2,-2, zz/)
      ijSel(  4,1:3) = (/ 2,-2, ww/)
      ijSel(  5,1:3) = (/ 3, 2, zz/)
      ijSel(  6,1:3) = (/ 3, 2, ww/)
      ijSel(  7,1:3) = (/ 1,-1, zz/)
      ijSel(  8,1:3) = (/ 1,-1, ww/)
      ijSel(  9,1:3) = (/ 2,-4, zz/)
      ijSel( 10,1:3) = (/ 2,-4, ww/)
      ijSel( 11,1:3) = (/ 1,-3, zz/)
      ijSel( 12,1:3) = (/ 1,-3, ww/)
      ijSel( 13,1:3) = (/ 4, 1, zz/)
      ijSel( 14,1:3) = (/ 4, 1, ww/)
      ijSel( 15,1:3) = (/ 2, 2, zz/)
      ijSel( 16,1:3) = (/-1,-2, zz/)
      ijSel( 17,1:3) = (/-1,-2, ww/)
      ijSel( 18,1:3) = (/ 2,-1, zz/)
      ijSel( 19,1:3) = (/ 1, 1, zz/)
      ijSel( 20,1:3) = (/ 3,-1, zz/)
      ijSel( 21,1:3) = (/ 3,-1, ww/)
      ijSel( 22,1:3) = (/ 2,-3, zz/)
      ijSel( 23,1:3) = (/-2,-3, zz/)
      ijSel( 24,1:3) = (/-2,-3, ww/)
      ijSel( 25,1:3) = (/-1,-4, zz/)
      ijSel( 26,1:3) = (/-1,-4, ww/)
      ijSel( 27,1:3) = (/ 3,-3, zz/)
      ijSel( 28,1:3) = (/ 3,-3, ww/)
      ijSel( 29,1:3) = (/ 1,-2, zz/)
      ijSel( 30,1:3) = (/ 3, 1, zz/)
      ijSel( 31,1:3) = (/ 4,-2, zz/)
      ijSel( 32,1:3) = (/ 4,-2, ww/)
      ijSel( 33,1:3) = (/ 4, 2, zz/)
      ijSel( 34,1:3) = (/ 5, 2, zz/)
      ijSel( 35,1:3) = (/ 2,-5, zz/)
      ijSel( 36,1:3) = (/-3,-4, zz/)
      ijSel( 37,1:3) = (/-3,-4, ww/)
      ijSel( 38,1:3) = (/ 4, 3, zz/)
      ijSel( 39,1:3) = (/ 4, 3, ww/)
      ijSel( 40,1:3) = (/ 1,-4, zz/)
      ijSel( 41,1:3) = (/ 5, 1, zz/)
      ijSel( 42,1:3) = (/ 1,-5, zz/)
      ijSel( 43,1:3) = (/ 4,-4, zz/)
      ijSel( 44,1:3) = (/ 4,-4, ww/)
      ijSel( 45,1:3) = (/-1,-3, zz/)
      ijSel( 46,1:3) = (/-1,-1, zz/)
      ijSel( 47,1:3) = (/ 3,-2, zz/)
      ijSel( 48,1:3) = (/ 4,-1, zz/)
      ijSel( 49,1:3) = (/-1,-5, zz/)
      ijSel( 50,1:3) = (/ 5,-1, zz/)
      ijSel( 51,1:3) = (/-2,-2, zz/)
      ijSel( 52,1:3) = (/-2,-4, zz/)
      ijSel( 53,1:3) = (/ 3, 3, zz/)
      ijSel( 54,1:3) = (/-3,-3, zz/)
      ijSel( 55,1:3) = (/ 3,-4, zz/)
      ijSel( 56,1:3) = (/ 4,-3, zz/)
      ijSel( 57,1:3) = (/ 5,-2, zz/)
      ijSel( 58,1:3) = (/-2,-5, zz/)
      ijSel( 59,1:3) = (/ 5, 3, zz/)
      ijSel( 60,1:3) = (/ 3,-5, zz/)
      ijSel( 61,1:3) = (/-3,-5, zz/)
      ijSel( 62,1:3) = (/ 5,-3, zz/)
      ijSel( 63,1:3) = (/-4,-5, zz/)
      ijSel( 64,1:3) = (/ 5, 4, zz/)
      ijSel( 65,1:3) = (/ 4,-5, zz/)
      ijSel( 66,1:3) = (/ 5,-4, zz/)
      ijSel( 67,1:3) = (/ 5,-5, zz/)
      ijSel( 68,1:3) = (/ 4, 4, zz/)
      ijSel( 69,1:3) = (/-4,-4, zz/)
      ijSel( 70,1:3) = (/-5,-5, zz/)
      ijSel( 71,1:3) = (/ 5, 5, zz/)

      ijSel( 72:,:)  = 0
      ijSel( 72:,3)  = -1


  return
  end subroutine


 subroutine get_VBFchannelHash_nosplit(ijSel,nijchannels)
  implicit none
  integer, intent(out) :: ijSel(1:121,1:3)
  integer, intent(out) :: nijchannels
  integer,parameter :: zzww=2, zz=1, ww=0

      ! u-d
      ijSel(  1,1:3) = (/ pdfUp_,pdfDn_,zzww/)
      ijSel(  2,1:3) = (/ pdfUp_,pdfStr_,zzww/)
      ijSel(  3,1:3) = (/ pdfUp_,pdfBot_,zzww/)
      ijSel(  4,1:3) = (/ pdfChm_,pdfDn_,zzww/)
      ijSel(  5,1:3) = (/ pdfChm_,pdfStr_,zzww/)
      ijSel(  6,1:3) = (/ pdfChm_,pdfBot_,zzww/)

	  ! ub-db
      ijSel(  7,1:3) = (/ pdfAUp_,pdfADn_,zzww/)
      ijSel(  8,1:3) = (/ pdfAUp_,pdfAStr_,zzww/)
      ijSel(  9,1:3) = (/ pdfAUp_,pdfABot_,zzww/)
      ijSel( 10,1:3) = (/ pdfAChm_,pdfADn_,zzww/)
      ijSel( 11,1:3) = (/ pdfAChm_,pdfAStr_,zzww/)
      ijSel( 12,1:3) = (/ pdfAChm_,pdfABot_,zzww/)

      ! u-ub
      ijSel( 13,1:3) = (/ pdfUp_,pdfAUp_,ww/)
      ijSel( 14,1:3) = (/ pdfUp_,pdfAUp_,zz/)
      ijSel( 15,1:3) = (/ pdfUp_,pdfAChm_,ww/)
      ijSel( 16,1:3) = (/ pdfUp_,pdfAChm_,zz/)
      ijSel( 17,1:3) = (/ pdfChm_,pdfAUp_,ww/)
      ijSel( 18,1:3) = (/ pdfChm_,pdfAUp_,zz/)
      ijSel( 19,1:3) = (/ pdfChm_,pdfAChm_,ww/)
      ijSel( 20,1:3) = (/ pdfChm_,pdfAChm_,zz/)

      ! d-db
      ijSel( 21,1:3) = (/ pdfDn_,pdfADn_,ww/)
      ijSel( 22,1:3) = (/ pdfDn_,pdfADn_,zz/)
      ijSel( 23,1:3) = (/ pdfDn_,pdfAStr_,ww/)
      ijSel( 24,1:3) = (/ pdfDn_,pdfAStr_,zz/)
      ijSel( 25,1:3) = (/ pdfDn_,pdfABot_,ww/)
      ijSel( 26,1:3) = (/ pdfDn_,pdfABot_,zz/)
      ijSel( 27,1:3) = (/ pdfStr_,pdfADn_,ww/)
      ijSel( 28,1:3) = (/ pdfStr_,pdfADn_,zz/)
      ijSel( 29,1:3) = (/ pdfStr_,pdfAStr_,ww/)
      ijSel( 30,1:3) = (/ pdfStr_,pdfAStr_,zz/)
      ijSel( 31,1:3) = (/ pdfStr_,pdfABot_,ww/)
      ijSel( 32,1:3) = (/ pdfStr_,pdfABot_,zz/)

      ! b-db
	  ! Some of these will be inefficient
      ijSel( 33,1:3) = (/ pdfBot_,pdfADn_,ww/)
      ijSel( 34,1:3) = (/ pdfBot_,pdfADn_,zz/)
      ijSel( 35,1:3) = (/ pdfBot_,pdfAStr_,ww/)
      ijSel( 36,1:3) = (/ pdfBot_,pdfAStr_,zz/)
      ijSel( 37,1:3) = (/ pdfBot_,pdfABot_,ww/)
      ijSel( 38,1:3) = (/ pdfBot_,pdfABot_,zz/)

      ! u-u
      ijSel( 39,1:3) = (/ pdfUp_,pdfUp_,zz/)
      ijSel( 40,1:3) = (/ pdfUp_,pdfChm_,zz/)
      ijSel( 41,1:3) = (/ pdfChm_,pdfChm_,zz/)

      ! d-d
      ijSel( 42,1:3) = (/ pdfDn_,pdfDn_,zz/)
      ijSel( 43,1:3) = (/ pdfDn_,pdfStr_,zz/)
      ijSel( 44,1:3) = (/ pdfDn_,pdfBot_,zz/)
      ijSel( 45,1:3) = (/ pdfStr_,pdfStr_,zz/)
      ijSel( 46,1:3) = (/ pdfStr_,pdfBot_,zz/)
      ijSel( 47,1:3) = (/ pdfBot_,pdfBot_,zz/)

      ! u-db
	   ijSel( 48,1:3) = (/ pdfUp_,pdfADn_,zz/)
      ijSel( 49,1:3) = (/ pdfUp_,pdfAStr_,zz/)
      ijSel( 50,1:3) = (/ pdfUp_,pdfABot_,zz/)
      ijSel( 51,1:3) = (/ pdfChm_,pdfADn_,zz/)
      ijSel( 52,1:3) = (/ pdfChm_,pdfAStr_,zz/)
      ijSel( 53,1:3) = (/ pdfChm_,pdfABot_,zz/)

      ! d-ub
      ijSel( 54,1:3) = (/ pdfDn_,pdfAUp_,zz/)
      ijSel( 55,1:3) = (/ pdfDn_,pdfAChm_,zz/)
      ijSel( 56,1:3) = (/ pdfStr_,pdfAUp_,zz/)
      ijSel( 57,1:3) = (/ pdfStr_,pdfAChm_,zz/)
      ijSel( 58,1:3) = (/ pdfBot_,pdfAUp_,zz/)
      ijSel( 59,1:3) = (/ pdfBot_,pdfAChm_,zz/)

	  ! ub-ub
      ijSel( 60,1:3) = (/ pdfAUp_,pdfAUp_,zz/)
      ijSel( 61,1:3) = (/ pdfAUp_,pdfAChm_,zz/)
      ijSel( 62,1:3) = (/ pdfAChm_,pdfAChm_,zz/)

      ! db-db
      ijSel( 63,1:3) = (/ pdfADn_,pdfADn_,zz/)
      ijSel( 64,1:3) = (/ pdfADn_,pdfAStr_,zz/)
      ijSel( 65,1:3) = (/ pdfADn_,pdfABot_,zz/)

      ijSel( 66,1:3) = (/ pdfAStr_,pdfAStr_,zz/)
      ijSel( 67,1:3) = (/ pdfAStr_,pdfABot_,zz/)
      ijSel( 68,1:3) = (/ pdfABot_,pdfABot_,zz/)

      nijchannels = 68

      ijSel( nijchannels+1:,:)  = 0
      ijSel( nijchannels+1:,3)  = -1

  return
  end subroutine



  subroutine get_HJJchannelHash(ijSel)
  implicit none
  integer, intent(out) :: ijSel(1:121,1:3)

      ijSel(  1,1:3) = (/ 0, 0, 1/)
      ijSel(  2,1:3) = (/ 0, 0, 2/)
      ijSel(  3,1:3) = (/ 2, 0, 1/)
      ijSel(  4,1:3) = (/ 1, 0, 1/)
      ijSel(  5,1:3) = (/ 0,-1, 1/)
      ijSel(  6,1:3) = (/ 0,-2, 1/)
      ijSel(  7,1:3) = (/ 3, 0, 1/)
      ijSel(  8,1:3) = (/ 0,-3, 1/)
      ijSel(  9,1:3) = (/ 4, 0, 1/)
      ijSel( 10,1:3) = (/ 0,-4, 1/)
      ijSel( 11,1:3) = (/ 5, 0, 1/)
      ijSel( 12,1:3) = (/ 0,-5, 1/)
      ijSel( 13,1:3) = (/ 2, 2, 1/)
      ijSel( 14,1:3) = (/ 2, 1, 1/)
      ijSel( 15,1:3) = (/ 2,-2, 1/)
      ijSel( 16,1:3) = (/ 2,-2, 2/)
      ijSel( 17,1:3) = (/ 2,-2, 3/)
      ijSel( 18,1:3) = (/ 1,-1, 1/)
      ijSel( 19,1:3) = (/ 1,-1, 2/)
      ijSel( 20,1:3) = (/ 1,-1, 3/)
      ijSel( 21,1:3) = (/ 2,-1, 1/)
      ijSel( 22,1:3) = (/ 1, 1, 1/)
      ijSel( 23,1:3) = (/ 3, 2, 1/)
      ijSel( 24,1:3) = (/ 2,-3, 1/)
      ijSel( 25,1:3) = (/ 4, 2, 1/)
      ijSel( 26,1:3) = (/ 1,-2, 1/)
      ijSel( 27,1:3) = (/ 2,-4, 1/)
      ijSel( 28,1:3) = (/ 1,-3, 1/)
      ijSel( 29,1:3) = (/ 2,-5, 1/)
      ijSel( 30,1:3) = (/ 3, 1, 1/)
      ijSel( 31,1:3) = (/ 5, 2, 1/)
      ijSel( 32,1:3) = (/ 4, 1, 1/)
      ijSel( 33,1:3) = (/ 1,-4, 1/)
      ijSel( 34,1:3) = (/ 5, 1, 1/)
      ijSel( 35,1:3) = (/ 1,-5, 1/)
      ijSel( 36,1:3) = (/-1,-2, 1/)
      ijSel( 37,1:3) = (/ 3,-1, 1/)
      ijSel( 38,1:3) = (/ 3,-3, 1/)
      ijSel( 39,1:3) = (/ 3,-3, 2/)
      ijSel( 40,1:3) = (/ 3,-3, 3/)
      ijSel( 41,1:3) = (/-2,-3, 1/)
      ijSel( 42,1:3) = (/-1,-3, 1/)
      ijSel( 43,1:3) = (/ 3,-2, 1/)
      ijSel( 44,1:3) = (/-1,-1, 1/)
      ijSel( 45,1:3) = (/ 4,-1, 1/)
      ijSel( 46,1:3) = (/-2,-2, 1/)
      ijSel( 47,1:3) = (/-1,-4, 1/)
      ijSel( 48,1:3) = (/ 4,-2, 1/)
      ijSel( 49,1:3) = (/-2,-4, 1/)
      ijSel( 50,1:3) = (/ 5,-1, 1/)
      ijSel( 51,1:3) = (/-1,-5, 1/)
      ijSel( 52,1:3) = (/ 4,-4, 1/)
      ijSel( 53,1:3) = (/ 4,-4, 2/)
      ijSel( 54,1:3) = (/ 4,-4, 3/)
      ijSel( 55,1:3) = (/ 4,-3, 1/)
      ijSel( 56,1:3) = (/ 4, 3, 1/)
      ijSel( 57,1:3) = (/ 3,-4, 1/)
      ijSel( 58,1:3) = (/-2,-5, 1/)
      ijSel( 59,1:3) = (/-3,-4, 1/)
      ijSel( 60,1:3) = (/ 5,-2, 1/)
      ijSel( 61,1:3) = (/ 3, 3, 1/)
      ijSel( 62,1:3) = (/-3,-3, 1/)
      ijSel( 63,1:3) = (/ 5, 3, 1/)
      ijSel( 64,1:3) = (/ 3,-5, 1/)
      ijSel( 65,1:3) = (/ 5,-3, 1/)
      ijSel( 66,1:3) = (/-3,-5, 1/)
      ijSel( 67,1:3) = (/ 5, 4, 1/)
      ijSel( 68,1:3) = (/ 5,-4, 1/)
      ijSel( 69,1:3) = (/-4,-5, 1/)
      ijSel( 70,1:3) = (/ 4,-5, 1/)
      ijSel( 71,1:3) = (/ 5,-5, 1/)
      ijSel( 72,1:3) = (/ 5,-5, 2/)
      ijSel( 73,1:3) = (/ 5,-5, 3/)
      ijSel( 74,1:3) = (/-4,-4, 1/)
      ijSel( 75,1:3) = (/ 4, 4, 1/)
      ijSel( 76,1:3) = (/-5,-5, 1/)
      ijSel( 77,1:3) = (/ 5, 5, 1/)

      ijSel( 78:,:) = 0
      ijSel( 78:,3) = -1

  return
  end subroutine


  subroutine get_HJJchannelHash_nosplit(ijSel,nijchannels)
  implicit none
  integer, intent(out) :: ijSel(1:121,1:3)
  integer, intent(out) :: nijchannels

      ijSel(  1,1:3) = (/ 0, 0, 1/)
      ijSel(  2,1:3) = (/ 0, 0, 2/)
      ijSel(  3,1:3) = (/ 2, 0, 1/)
      ijSel(  4,1:3) = (/ 1, 0, 1/)
      ijSel(  5,1:3) = (/ 0,-1, 1/)
      ijSel(  6,1:3) = (/ 0,-2, 1/)
      ijSel(  7,1:3) = (/ 3, 0, 1/)
      ijSel(  8,1:3) = (/ 0,-3, 1/)
      ijSel(  9,1:3) = (/ 4, 0, 1/)
      ijSel( 10,1:3) = (/ 0,-4, 1/)
      ijSel( 11,1:3) = (/ 5, 0, 1/)
      ijSel( 12,1:3) = (/ 0,-5, 1/)
      ijSel( 13,1:3) = (/ 2, 2, 1/)
      ijSel( 14,1:3) = (/ 2, 1, 1/)
      ijSel( 15,1:3) = (/ 2,-2, 1/)
      ijSel( 16,1:3) = (/ 2,-2, 2/)
      ijSel( 17,1:3) = (/ 2,-2, 3/)
      ijSel( 18,1:3) = (/ 1,-1, 1/)
      ijSel( 19,1:3) = (/ 1,-1, 2/)
      ijSel( 20,1:3) = (/ 1,-1, 3/)
      ijSel( 21,1:3) = (/ 2,-1, 1/)
      ijSel( 22,1:3) = (/ 1, 1, 1/)
      ijSel( 23,1:3) = (/ 3, 2, 1/)
      ijSel( 24,1:3) = (/ 2,-3, 1/)
      ijSel( 25,1:3) = (/ 4, 2, 1/)
      ijSel( 26,1:3) = (/ 1,-2, 1/)
      ijSel( 27,1:3) = (/ 2,-4, 1/)
      ijSel( 28,1:3) = (/ 1,-3, 1/)
      ijSel( 29,1:3) = (/ 2,-5, 1/)
      ijSel( 30,1:3) = (/ 3, 1, 1/)
      ijSel( 31,1:3) = (/ 5, 2, 1/)
      ijSel( 32,1:3) = (/ 4, 1, 1/)
      ijSel( 33,1:3) = (/ 1,-4, 1/)
      ijSel( 34,1:3) = (/ 5, 1, 1/)
      ijSel( 35,1:3) = (/ 1,-5, 1/)
      ijSel( 36,1:3) = (/-1,-2, 1/)
      ijSel( 37,1:3) = (/ 3,-1, 1/)
      ijSel( 38,1:3) = (/ 3,-3, 1/)
      ijSel( 39,1:3) = (/ 3,-3, 2/)
      ijSel( 40,1:3) = (/ 3,-3, 3/)
      ijSel( 41,1:3) = (/-2,-3, 1/)
      ijSel( 42,1:3) = (/-1,-3, 1/)
      ijSel( 43,1:3) = (/ 3,-2, 1/)
      ijSel( 44,1:3) = (/-1,-1, 1/)
      ijSel( 45,1:3) = (/ 4,-1, 1/)
      ijSel( 46,1:3) = (/-2,-2, 1/)
      ijSel( 47,1:3) = (/-1,-4, 1/)
      ijSel( 48,1:3) = (/ 4,-2, 1/)
      ijSel( 49,1:3) = (/-2,-4, 1/)
      ijSel( 50,1:3) = (/ 5,-1, 1/)
      ijSel( 51,1:3) = (/-1,-5, 1/)
      ijSel( 52,1:3) = (/ 4,-4, 1/)
      ijSel( 53,1:3) = (/ 4,-4, 2/)
      ijSel( 54,1:3) = (/ 4,-4, 3/)
      ijSel( 55,1:3) = (/ 4,-3, 1/)
      ijSel( 56,1:3) = (/ 4, 3, 1/)
      ijSel( 57,1:3) = (/ 3,-4, 1/)
      ijSel( 58,1:3) = (/-2,-5, 1/)
      ijSel( 59,1:3) = (/-3,-4, 1/)
      ijSel( 60,1:3) = (/ 5,-2, 1/)
      ijSel( 61,1:3) = (/ 3, 3, 1/)
      ijSel( 62,1:3) = (/-3,-3, 1/)
      ijSel( 63,1:3) = (/ 5, 3, 1/)
      ijSel( 64,1:3) = (/ 3,-5, 1/)
      ijSel( 65,1:3) = (/ 5,-3, 1/)
      ijSel( 66,1:3) = (/-3,-5, 1/)
      ijSel( 67,1:3) = (/ 5, 4, 1/)
      ijSel( 68,1:3) = (/ 5,-4, 1/)
      ijSel( 69,1:3) = (/-4,-5, 1/)
      ijSel( 70,1:3) = (/ 4,-5, 1/)
      ijSel( 71,1:3) = (/ 5,-5, 1/)
      ijSel( 72,1:3) = (/ 5,-5, 2/)
      ijSel( 73,1:3) = (/ 5,-5, 3/)
      ijSel( 74,1:3) = (/-4,-4, 1/)
      ijSel( 75,1:3) = (/ 4, 4, 1/)
      ijSel( 76,1:3) = (/-5,-5, 1/)
      ijSel( 77,1:3) = (/ 5, 5, 1/)

      nijchannels=77

      ijSel( nijchannels+1:,:) = 0
      ijSel( nijchannels+1:,3) = -1

  return
  end subroutine


  subroutine get_GENchannelHash(ijSel)
  implicit none
  integer, intent(out) :: ijSel(1:121,1:3)



      ijSel(  1,1:3) = (/-5,-5, 1/)
      ijSel(  2,1:3) = (/-5,-4, 1/)
      ijSel(  3,1:3) = (/-5,-3, 1/)
      ijSel(  4,1:3) = (/-5,-2, 1/)
      ijSel(  5,1:3) = (/-5,-1, 1/)
      ijSel(  6,1:3) = (/-5, 0, 1/)
      ijSel(  7,1:3) = (/-5, 1, 1/)
      ijSel(  8,1:3) = (/-5, 2, 1/)
      ijSel(  9,1:3) = (/-5, 3, 1/)
      ijSel( 10,1:3) = (/-5, 4, 1/)
      ijSel( 11,1:3) = (/-5, 5, 1/)
      ijSel( 12,1:3) = (/-4,-5, 1/)
      ijSel( 13,1:3) = (/-4,-4, 1/)
      ijSel( 14,1:3) = (/-4,-3, 1/)
      ijSel( 15,1:3) = (/-4,-2, 1/)
      ijSel( 16,1:3) = (/-4,-1, 1/)
      ijSel( 17,1:3) = (/-4, 0, 1/)
      ijSel( 18,1:3) = (/-4, 1, 1/)
      ijSel( 19,1:3) = (/-4, 2, 1/)
      ijSel( 20,1:3) = (/-4, 3, 1/)
      ijSel( 21,1:3) = (/-4, 4, 1/)
      ijSel( 22,1:3) = (/-4, 5, 1/)
      ijSel( 23,1:3) = (/-3,-5, 1/)
      ijSel( 24,1:3) = (/-3,-4, 1/)
      ijSel( 25,1:3) = (/-3,-3, 1/)
      ijSel( 26,1:3) = (/-3,-2, 1/)
      ijSel( 27,1:3) = (/-3,-1, 1/)
      ijSel( 28,1:3) = (/-3, 0, 1/)
      ijSel( 29,1:3) = (/-3, 1, 1/)
      ijSel( 30,1:3) = (/-3, 2, 1/)
      ijSel( 31,1:3) = (/-3, 3, 1/)
      ijSel( 32,1:3) = (/-3, 4, 1/)
      ijSel( 33,1:3) = (/-3, 5, 1/)
      ijSel( 34,1:3) = (/-2,-5, 1/)
      ijSel( 35,1:3) = (/-2,-4, 1/)
      ijSel( 36,1:3) = (/-2,-3, 1/)
      ijSel( 37,1:3) = (/-2,-2, 1/)
      ijSel( 38,1:3) = (/-2,-1, 1/)
      ijSel( 39,1:3) = (/-2, 0, 1/)
      ijSel( 40,1:3) = (/-2, 1, 1/)
      ijSel( 41,1:3) = (/-2, 2, 1/)
      ijSel( 42,1:3) = (/-2, 3, 1/)
      ijSel( 43,1:3) = (/-2, 4, 1/)
      ijSel( 44,1:3) = (/-2, 5, 1/)
      ijSel( 45,1:3) = (/-1,-5, 1/)
      ijSel( 46,1:3) = (/-1,-4, 1/)
      ijSel( 47,1:3) = (/-1,-3, 1/)
      ijSel( 48,1:3) = (/-1,-2, 1/)
      ijSel( 49,1:3) = (/-1,-1, 1/)
      ijSel( 50,1:3) = (/-1, 0, 1/)
      ijSel( 51,1:3) = (/-1, 1, 1/)
      ijSel( 52,1:3) = (/-1, 2, 1/)
      ijSel( 53,1:3) = (/-1, 3, 1/)
      ijSel( 54,1:3) = (/-1, 4, 1/)
      ijSel( 55,1:3) = (/-1, 5, 1/)
      ijSel( 56,1:3) = (/ 0,-5, 1/)
      ijSel( 57,1:3) = (/ 0,-4, 1/)
      ijSel( 58,1:3) = (/ 0,-3, 1/)
      ijSel( 59,1:3) = (/ 0,-2, 1/)
      ijSel( 60,1:3) = (/ 0,-1, 1/)
      ijSel( 61,1:3) = (/ 0, 0, 1/)
      ijSel( 62,1:3) = (/ 0, 1, 1/)
      ijSel( 63,1:3) = (/ 0, 2, 1/)
      ijSel( 64,1:3) = (/ 0, 3, 1/)
      ijSel( 65,1:3) = (/ 0, 4, 1/)
      ijSel( 66,1:3) = (/ 0, 5, 1/)
      ijSel( 67,1:3) = (/ 1,-5, 1/)
      ijSel( 68,1:3) = (/ 1,-4, 1/)
      ijSel( 69,1:3) = (/ 1,-3, 1/)
      ijSel( 70,1:3) = (/ 1,-2, 1/)
      ijSel( 71,1:3) = (/ 1,-1, 1/)
      ijSel( 72,1:3) = (/ 1, 0, 1/)
      ijSel( 73,1:3) = (/ 1, 1, 1/)
      ijSel( 74,1:3) = (/ 1, 2, 1/)
      ijSel( 75,1:3) = (/ 1, 3, 1/)
      ijSel( 76,1:3) = (/ 1, 4, 1/)
      ijSel( 77,1:3) = (/ 1, 5, 1/)
      ijSel( 78,1:3) = (/ 2,-5, 1/)
      ijSel( 79,1:3) = (/ 2,-4, 1/)
      ijSel( 80,1:3) = (/ 2,-3, 1/)
      ijSel( 81,1:3) = (/ 2,-2, 1/)
      ijSel( 82,1:3) = (/ 2,-1, 1/)
      ijSel( 83,1:3) = (/ 2, 0, 1/)
      ijSel( 84,1:3) = (/ 2, 1, 1/)
      ijSel( 85,1:3) = (/ 2, 2, 1/)
      ijSel( 86,1:3) = (/ 2, 3, 1/)
      ijSel( 87,1:3) = (/ 2, 4, 1/)
      ijSel( 88,1:3) = (/ 2, 5, 1/)
      ijSel( 89,1:3) = (/ 3,-5, 1/)
      ijSel( 90,1:3) = (/ 3,-4, 1/)
      ijSel( 91,1:3) = (/ 3,-3, 1/)
      ijSel( 92,1:3) = (/ 3,-2, 1/)
      ijSel( 93,1:3) = (/ 3,-1, 1/)
      ijSel( 94,1:3) = (/ 3, 0, 1/)
      ijSel( 95,1:3) = (/ 3, 1, 1/)
      ijSel( 96,1:3) = (/ 3, 2, 1/)
      ijSel( 97,1:3) = (/ 3, 3, 1/)
      ijSel( 98,1:3) = (/ 3, 4, 1/)
      ijSel( 99,1:3) = (/ 3, 5, 1/)
      ijSel(100,1:3) = (/ 4,-5, 1/)
      ijSel(101,1:3) = (/ 4,-4, 1/)
      ijSel(102,1:3) = (/ 4,-3, 1/)
      ijSel(103,1:3) = (/ 4,-2, 1/)
      ijSel(104,1:3) = (/ 4,-1, 1/)
      ijSel(105,1:3) = (/ 4, 0, 1/)
      ijSel(106,1:3) = (/ 4, 1, 1/)
      ijSel(107,1:3) = (/ 4, 2, 1/)
      ijSel(108,1:3) = (/ 4, 3, 1/)
      ijSel(109,1:3) = (/ 4, 4, 1/)
      ijSel(110,1:3) = (/ 4, 5, 1/)
      ijSel(111,1:3) = (/ 5,-5, 1/)
      ijSel(112,1:3) = (/ 5,-4, 1/)
      ijSel(113,1:3) = (/ 5,-3, 1/)
      ijSel(114,1:3) = (/ 5,-2, 1/)
      ijSel(115,1:3) = (/ 5,-1, 1/)
      ijSel(116,1:3) = (/ 5, 0, 1/)
      ijSel(117,1:3) = (/ 5, 1, 1/)
      ijSel(118,1:3) = (/ 5, 2, 1/)
      ijSel(119,1:3) = (/ 5, 3, 1/)
      ijSel(120,1:3) = (/ 5, 4, 1/)
      ijSel(121,1:3) = (/ 5, 5, 1/)

  return
  end subroutine





  subroutine init_VBFoffshChannelHash()
  implicit none
  integer ch,i

   if( .not. deterministicME ) then
      ij_ChannelHash(  1,:) = (/ 2, 1, 0,0/)
      ij_ChannelHash(  2,:) = (/ 1, 2, 0,0/)
      ij_ChannelHash(  3,:) = (/ 2,-2, 0,0/)
      ij_ChannelHash(  4,:) = (/-2, 2, 0,0/)
      ij_ChannelHash(  5,:) = (/ 2, 3, 0,0/)
      ij_ChannelHash(  6,:) = (/ 3, 2, 0,0/)
      ij_ChannelHash(  7,:) = (/ 1,-1, 0,0/)
      ij_ChannelHash(  8,:) = (/-1, 1, 0,0/)
      ij_ChannelHash(  9,:) = (/ 1, 3, 0,0/)
      ij_ChannelHash( 10,:) = (/ 3, 1, 0,0/)
      ij_ChannelHash( 11,:) = (/ 2,-4, 0,0/)
      ij_ChannelHash( 12,:) = (/-4, 2, 0,0/)
      ij_ChannelHash( 13,:) = (/ 1, 4, 0,0/)
      ij_ChannelHash( 14,:) = (/ 4, 1, 0,0/)
      ij_ChannelHash( 15,:) = (/-5,-5, 0,0/)
      ij_ChannelHash( 16,:) = (/-5,-4, 0,0/)
      ij_ChannelHash( 17,:) = (/-5,-3, 0,0/)
      ij_ChannelHash( 18,:) = (/-5,-2, 0,0/)
      ij_ChannelHash( 19,:) = (/-5,-1, 0,0/)
      ij_ChannelHash( 20,:) = (/-5, 1, 0,0/)
      ij_ChannelHash( 21,:) = (/-5, 2, 0,0/)
      ij_ChannelHash( 22,:) = (/-5, 3, 0,0/)
      ij_ChannelHash( 23,:) = (/-5, 4, 0,0/)
      ij_ChannelHash( 24,:) = (/-5, 5, 0,0/)
      ij_ChannelHash( 25,:) = (/-4,-5, 0,0/)
      ij_ChannelHash( 26,:) = (/-4,-4, 0,0/)
      ij_ChannelHash( 27,:) = (/-4,-3, 0,0/)
      ij_ChannelHash( 28,:) = (/-4,-2, 0,0/)
      ij_ChannelHash( 29,:) = (/-4,-1, 0,0/)
      ij_ChannelHash( 30,:) = (/-4, 1, 0,0/)
      ij_ChannelHash( 31,:) = (/-4, 3, 0,0/)
      ij_ChannelHash( 32,:) = (/-4, 4, 0,0/)
      ij_ChannelHash( 33,:) = (/-4, 5, 0,0/)
      ij_ChannelHash( 34,:) = (/-3,-5, 0,0/)
      ij_ChannelHash( 35,:) = (/-3,-4, 0,0/)
      ij_ChannelHash( 36,:) = (/-3,-3, 0,0/)
      ij_ChannelHash( 37,:) = (/-3,-2, 0,0/)
      ij_ChannelHash( 38,:) = (/-3,-1, 0,0/)
      ij_ChannelHash( 39,:) = (/-3, 1, 0,0/)
      ij_ChannelHash( 40,:) = (/-3, 2, 0,0/)
      ij_ChannelHash( 41,:) = (/-3, 3, 0,0/)
      ij_ChannelHash( 42,:) = (/-3, 4, 0,0/)
      ij_ChannelHash( 43,:) = (/-3, 5, 0,0/)
      ij_ChannelHash( 44,:) = (/-2,-5, 0,0/)
      ij_ChannelHash( 45,:) = (/-2,-4, 0,0/)
      ij_ChannelHash( 46,:) = (/-2,-3, 0,0/)
      ij_ChannelHash( 47,:) = (/-2,-2, 0,0/)
      ij_ChannelHash( 48,:) = (/-2,-1, 0,0/)
      ij_ChannelHash( 49,:) = (/-2, 1, 0,0/)
      ij_ChannelHash( 50,:) = (/-2, 3, 0,0/)
      ij_ChannelHash( 51,:) = (/-2, 4, 0,0/)
      ij_ChannelHash( 52,:) = (/-2, 5, 0,0/)
      ij_ChannelHash( 53,:) = (/-1,-5, 0,0/)
      ij_ChannelHash( 54,:) = (/-1,-4, 0,0/)
      ij_ChannelHash( 55,:) = (/-1,-3, 0,0/)
      ij_ChannelHash( 56,:) = (/-1,-2, 0,0/)
      ij_ChannelHash( 57,:) = (/-1,-1, 0,0/)
      ij_ChannelHash( 58,:) = (/-1, 2, 0,0/)
      ij_ChannelHash( 59,:) = (/-1, 3, 0,0/)
      ij_ChannelHash( 60,:) = (/-1, 4, 0,0/)
      ij_ChannelHash( 61,:) = (/-1, 5, 0,0/)
      ij_ChannelHash( 62,:) = (/ 1,-5, 0,0/)
      ij_ChannelHash( 63,:) = (/ 1,-4, 0,0/)
      ij_ChannelHash( 64,:) = (/ 1,-3, 0,0/)
      ij_ChannelHash( 65,:) = (/ 1,-2, 0,0/)
      ij_ChannelHash( 66,:) = (/ 1, 1, 0,0/)
      ij_ChannelHash( 67,:) = (/ 1, 5, 0,0/)
      ij_ChannelHash( 68,:) = (/ 2,-5, 0,0/)
      ij_ChannelHash( 69,:) = (/ 2,-3, 0,0/)
      ij_ChannelHash( 70,:) = (/ 2,-1, 0,0/)
      ij_ChannelHash( 71,:) = (/ 2, 2, 0,0/)
      ij_ChannelHash( 72,:) = (/ 2, 4, 0,0/)
      ij_ChannelHash( 73,:) = (/ 2, 5, 0,0/)
      ij_ChannelHash( 74,:) = (/ 3,-5, 0,0/)
      ij_ChannelHash( 75,:) = (/ 3,-4, 0,0/)
      ij_ChannelHash( 76,:) = (/ 3,-3, 0,0/)
      ij_ChannelHash( 77,:) = (/ 3,-2, 0,0/)
      ij_ChannelHash( 78,:) = (/ 3,-1, 0,0/)
      ij_ChannelHash( 79,:) = (/ 3, 3, 0,0/)
      ij_ChannelHash( 80,:) = (/ 3, 4, 0,0/)
      ij_ChannelHash( 81,:) = (/ 3, 5, 0,0/)
      ij_ChannelHash( 82,:) = (/ 4,-5, 0,0/)
      ij_ChannelHash( 83,:) = (/ 4,-4, 0,0/)
      ij_ChannelHash( 84,:) = (/ 4,-3, 0,0/)
      ij_ChannelHash( 85,:) = (/ 4,-2, 0,0/)
      ij_ChannelHash( 86,:) = (/ 4,-1, 0,0/)
      ij_ChannelHash( 87,:) = (/ 4, 2, 0,0/)
      ij_ChannelHash( 88,:) = (/ 4, 3, 0,0/)
      ij_ChannelHash( 89,:) = (/ 4, 4, 0,0/)
      ij_ChannelHash( 90,:) = (/ 4, 5, 0,0/)
      ij_ChannelHash( 91,:) = (/ 5,-5, 0,0/)
      ij_ChannelHash( 92,:) = (/ 5,-4, 0,0/)
      ij_ChannelHash( 93,:) = (/ 5,-3, 0,0/)
      ij_ChannelHash( 94,:) = (/ 5,-2, 0,0/)
      ij_ChannelHash( 95,:) = (/ 5,-1, 0,0/)
      ij_ChannelHash( 96,:) = (/ 5, 1, 0,0/)
      ij_ChannelHash( 97,:) = (/ 5, 2, 0,0/)
      ij_ChannelHash( 98,:) = (/ 5, 3, 0,0/)
      ij_ChannelHash( 99,:) = (/ 5, 4, 0,0/)
      ij_ChannelHash(100,:) = (/ 5, 5, 0,0/)

      NumberOfChannels = 100

    else

      ij_ChannelHash(1,:) = (/ Up_ , Chm_ , Up_ , Chm_ /)
      ij_ChannelHash(2,:) = (/ Dn_ , Str_ , Dn_ , Str_ /)
      ij_ChannelHash(3,:) = (/ Dn_ , Bot_ , Dn_ , Bot_ /)
      ij_ChannelHash(4,:) = (/ Str_ , Bot_ , Str_ , Bot_ /)
      ij_ChannelHash(5,:) = (/ Up_ , Str_ , Up_ , Str_ /)
      ij_ChannelHash(6,:) = (/ Up_ , Bot_ , Up_ , Bot_ /)
      ij_ChannelHash(7,:) = (/ Chm_ , Bot_ , Chm_ , Bot_ /)
      ij_ChannelHash(8,:) = (/ Dn_ , Chm_ , Dn_ , Chm_ /)
      ij_ChannelHash(9,:) = (/ Dn_ , Up_ , Dn_ , Up_ /)
      ij_ChannelHash(10,:) = (/ Str_ , Chm_ , Str_ , Chm_ /)
      ij_ChannelHash(11,:) = (/ Dn_ , Chm_ , Up_ , Str_ /)
      ij_ChannelHash(12,:) = (/ Up_ , Str_ , Dn_ , Chm_ /)
      ij_ChannelHash(13,:) = (/ Up_ , Up_ , Up_ , Up_ /)
      ij_ChannelHash(14,:) = (/ Chm_ , Chm_ , Chm_ , Chm_ /)
      ij_ChannelHash(15,:) = (/ Dn_ , Dn_ , Dn_ , Dn_ /)
      ij_ChannelHash(16,:) = (/ Str_ , Str_ , Str_ , Str_ /)
      ij_ChannelHash(17,:) = (/ Bot_ , Bot_ , Bot_ , Bot_ /)
      ij_ChannelHash(18,:) = (/ Chm_ , Up_ , Up_ , Chm_ /)
      ij_ChannelHash(19,:) = (/ Str_ , Dn_ , Dn_ , Str_ /)
      ij_ChannelHash(20,:) = (/ Bot_ , Dn_ , Dn_ , Bot_ /)
      ij_ChannelHash(21,:) = (/ Bot_ , Str_ , Str_ , Bot_ /)
      ij_ChannelHash(22,:) = (/ Str_ , Up_ , Up_ , Str_ /)
      ij_ChannelHash(23,:) = (/ Bot_ , Up_ , Up_ , Bot_ /)
      ij_ChannelHash(24,:) = (/ Bot_ , Chm_ , Chm_ , Bot_ /)
      ij_ChannelHash(25,:) = (/ Chm_ , Dn_ , Dn_ , Chm_ /)
      ij_ChannelHash(26,:) = (/ Up_ , Dn_ , Dn_ , Up_ /)
      ij_ChannelHash(27,:) = (/ Chm_ , Str_ , Str_ , Chm_ /)
      ij_ChannelHash(28,:) = (/ Chm_ , Dn_ , Up_ , Str_ /)
      ij_ChannelHash(29,:) = (/ Str_ , Up_ , Dn_ , Chm_ /)
      ij_ChannelHash(30,:) = (/ Up_ , Up_ , Up_ , Up_ /)
      ij_ChannelHash(31,:) = (/ Chm_ , Chm_ , Chm_ , Chm_ /)
      ij_ChannelHash(32,:) = (/ Dn_ , Dn_ , Dn_ , Dn_ /)
      ij_ChannelHash(33,:) = (/ Str_ , Str_ , Str_ , Str_ /)
      ij_ChannelHash(34,:) = (/ Bot_ , Bot_ , Bot_ , Bot_ /)

      ij_ChannelHash(35,:) = (/ AChm_ , AUp_ , AChm_ , AUp_ /)
      ij_ChannelHash(36,:) = (/ AStr_ , ADn_ , AStr_ , ADn_ /)
      ij_ChannelHash(37,:) = (/ ABot_ , ADn_ , ABot_ , ADn_ /)
      ij_ChannelHash(38,:) = (/ ABot_ , AStr_ , ABot_ , AStr_ /)
      ij_ChannelHash(39,:) = (/ AStr_ , AUp_ , AStr_ , AUp_ /)
      ij_ChannelHash(40,:) = (/ ABot_ , AUp_ , ABot_ , AUp_ /)
      ij_ChannelHash(41,:) = (/ ABot_ , AChm_ , ABot_ , AChm_ /)
      ij_ChannelHash(42,:) = (/ AChm_ , ADn_ , AChm_ , ADn_ /)
      ij_ChannelHash(43,:) = (/ AUp_ , ADn_ , AUp_ , ADn_ /)
      ij_ChannelHash(44,:) = (/ AChm_ , AStr_ , AChm_ , AStr_ /)
      ij_ChannelHash(45,:) = (/ AStr_ , AUp_ , AChm_ , ADn_ /)
      ij_ChannelHash(46,:) = (/ AChm_ , ADn_ , AStr_ , AUp_ /)
      ij_ChannelHash(47,:) = (/ AUp_ , AUp_ , AUp_ , AUp_ /)
      ij_ChannelHash(48,:) = (/ AChm_ , AChm_ , AChm_ , AChm_ /)
      ij_ChannelHash(49,:) = (/ ADn_ , ADn_ , ADn_ , ADn_ /)
      ij_ChannelHash(50,:) = (/ AStr_ , AStr_ , AStr_ , AStr_ /)
      ij_ChannelHash(51,:) = (/ ABot_ , ABot_ , ABot_ , ABot_ /)
      ij_ChannelHash(52,:) = (/ AUp_ , AChm_ , AChm_ , AUp_ /)
      ij_ChannelHash(53,:) = (/ ADn_ , AStr_ , AStr_ , ADn_ /)
      ij_ChannelHash(54,:) = (/ ADn_ , ABot_ , ABot_ , ADn_ /)
      ij_ChannelHash(55,:) = (/ AStr_ , ABot_ , ABot_ , AStr_ /)
      ij_ChannelHash(56,:) = (/ AUp_ , AStr_ , AStr_ , AUp_ /)
      ij_ChannelHash(57,:) = (/ AUp_ , ABot_ , ABot_ , AUp_ /)
      ij_ChannelHash(58,:) = (/ AChm_ , ABot_ , ABot_ , AChm_ /)
      ij_ChannelHash(59,:) = (/ ADn_ , AChm_ , AChm_ , ADn_ /)
      ij_ChannelHash(60,:) = (/ ADn_ , AUp_ , AUp_ , ADn_ /)
      ij_ChannelHash(61,:) = (/ AStr_ , AChm_ , AChm_ , AStr_ /)
      ij_ChannelHash(62,:) = (/ AUp_ , AStr_ , AChm_ , ADn_ /)
      ij_ChannelHash(63,:) = (/ ADn_ , AChm_ , AStr_ , AUp_ /)
      ij_ChannelHash(64,:) = (/ AUp_ , AUp_ , AUp_ , AUp_ /)
      ij_ChannelHash(65,:) = (/ AChm_ , AChm_ , AChm_ , AChm_ /)
      ij_ChannelHash(66,:) = (/ ADn_ , ADn_ , ADn_ , ADn_ /)
      ij_ChannelHash(67,:) = (/ AStr_ , AStr_ , AStr_ , AStr_ /)
      ij_ChannelHash(68,:) = (/ ABot_ , ABot_ , ABot_ , ABot_ /)

      ij_ChannelHash(69,:) = (/ AUp_ , Chm_ , AUp_ , Chm_ /)
      ij_ChannelHash(70,:) = (/ ADn_ , Str_ , ADn_ , Str_ /)
      ij_ChannelHash(71,:) = (/ ADn_ , Bot_ , ADn_ , Bot_ /)
      ij_ChannelHash(72,:) = (/ AStr_ , Bot_ , AStr_ , Bot_ /)
      ij_ChannelHash(73,:) = (/ AUp_ , Str_ , AUp_ , Str_ /)
      ij_ChannelHash(74,:) = (/ AUp_ , Bot_ , AUp_ , Bot_ /)
      ij_ChannelHash(75,:) = (/ AChm_ , Bot_ , AChm_ , Bot_ /)
      ij_ChannelHash(76,:) = (/ ADn_ , Chm_ , ADn_ , Chm_ /)
      ij_ChannelHash(77,:) = (/ ADn_ , Up_ , ADn_ , Up_ /)
      ij_ChannelHash(78,:) = (/ AStr_ , Chm_ , AStr_ , Chm_ /)
      ij_ChannelHash(79,:) = (/ AUp_ , Chm_ , ADn_ , Str_ /)
      ij_ChannelHash(80,:) = (/ ADn_ , Str_ , AUp_ , Chm_ /)
      ij_ChannelHash(81,:) = (/ AUp_ , Up_ , AUp_ , Up_ /)
      ij_ChannelHash(82,:) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
      ij_ChannelHash(83,:) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
      ij_ChannelHash(84,:) = (/ AStr_ , Str_ , AStr_ , Str_ /)
      ij_ChannelHash(85,:) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)
      ij_ChannelHash(86,:) = (/ AChm_ , Up_ , AChm_ , Up_ /)
      ij_ChannelHash(87,:) = (/ AStr_ , Dn_ , AStr_ , Dn_ /)
      ij_ChannelHash(88,:) = (/ ABot_ , Dn_ , ABot_ , Dn_ /)
      ij_ChannelHash(89,:) = (/ ABot_ , Str_ , ABot_ , Str_ /)
      ij_ChannelHash(90,:) = (/ AStr_ , Up_ , AStr_ , Up_ /)
      ij_ChannelHash(91,:) = (/ ABot_ , Up_ , ABot_ , Up_ /)
      ij_ChannelHash(92,:) = (/ ABot_ , Chm_ , ABot_ , Chm_ /)
      ij_ChannelHash(93,:) = (/ AChm_ , Dn_ , AChm_ , Dn_ /)
      ij_ChannelHash(94,:) = (/ AUp_ , Dn_ , AUp_ , Dn_ /)
      ij_ChannelHash(95,:) = (/ AChm_ , Str_ , AChm_ , Str_ /)
      ij_ChannelHash(96,:) = (/ AStr_ , Dn_ , AChm_ , Up_ /)
      ij_ChannelHash(97,:) = (/ AChm_ , Up_ , AStr_ , Dn_ /)
      ij_ChannelHash(98,:) = (/ AUp_ , Up_ , AUp_ , Up_ /)
      ij_ChannelHash(99,:) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
      ij_ChannelHash(100,:) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
      ij_ChannelHash(101,:) = (/ AStr_ , Str_ , AStr_ , Str_ /)
      ij_ChannelHash(102,:) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)

      ij_ChannelHash(103,:) = (/ Chm_ , AUp_ , AUp_ , Chm_ /)
      ij_ChannelHash(104,:) = (/ Str_ , ADn_ , ADn_ , Str_ /)
      ij_ChannelHash(105,:) = (/ Bot_ , ADn_ , ADn_ , Bot_ /)
      ij_ChannelHash(106,:) = (/ Bot_ , AStr_ , AStr_ , Bot_ /)
      ij_ChannelHash(107,:) = (/ Str_ , AUp_ , AUp_ , Str_ /)
      ij_ChannelHash(108,:) = (/ Bot_ , AUp_ , AUp_ , Bot_ /)
      ij_ChannelHash(109,:) = (/ Bot_ , AChm_ , AChm_ , Bot_ /)
      ij_ChannelHash(110,:) = (/ Chm_ , ADn_ , ADn_ , Chm_ /)
      ij_ChannelHash(111,:) = (/ Up_ , ADn_ , ADn_ , Up_ /)
      ij_ChannelHash(112,:) = (/ Chm_ , AStr_ , AStr_ , Chm_ /)
      ij_ChannelHash(113,:) = (/ Chm_ , AUp_ , ADn_ , Str_ /)
      ij_ChannelHash(114,:) = (/ Str_ , ADn_ , AUp_ , Chm_ /)
      ij_ChannelHash(115,:) = (/ Up_ , AUp_ , AUp_ , Up_ /)
      ij_ChannelHash(116,:) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
      ij_ChannelHash(117,:) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
      ij_ChannelHash(118,:) = (/ Str_ , AStr_ , AStr_ , Str_ /)
      ij_ChannelHash(119,:) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)
      ij_ChannelHash(120,:) = (/ Up_ , AChm_ , AChm_ , Up_ /)
      ij_ChannelHash(121,:) = (/ Dn_ , AStr_ , AStr_ , Dn_ /)
      ij_ChannelHash(122,:) = (/ Dn_ , ABot_ , ABot_ , Dn_ /)
      ij_ChannelHash(123,:) = (/ Str_ , ABot_ , ABot_ , Str_ /)
      ij_ChannelHash(124,:) = (/ Up_ , AStr_ , AStr_ , Up_ /)
      ij_ChannelHash(125,:) = (/ Up_ , ABot_ , ABot_ , Up_ /)
      ij_ChannelHash(126,:) = (/ Chm_ , ABot_ , ABot_ , Chm_ /)
      ij_ChannelHash(127,:) = (/ Dn_ , AChm_ , AChm_ , Dn_ /)
      ij_ChannelHash(128,:) = (/ Dn_ , AUp_ , AUp_ , Dn_ /)
      ij_ChannelHash(129,:) = (/ Str_ , AChm_ , AChm_ , Str_ /)
      ij_ChannelHash(130,:) = (/ Dn_ , AStr_ , AChm_ , Up_ /)
      ij_ChannelHash(131,:) = (/ Up_ , AChm_ , AStr_ , Dn_ /)
      ij_ChannelHash(132,:) = (/ Up_ , AUp_ , AUp_ , Up_ /)
      ij_ChannelHash(133,:) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
      ij_ChannelHash(134,:) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
      ij_ChannelHash(135,:) = (/ Str_ , AStr_ , AStr_ , Str_ /)
      ij_ChannelHash(136,:) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)

      ij_ChannelHash(137,:) = (/ Up_ , AUp_ , AChm_ , Chm_ /)
      ij_ChannelHash(138,:) = (/ Dn_ , ADn_ , AStr_ , Str_ /)
      ij_ChannelHash(139,:) = (/ Dn_ , ADn_ , ABot_ , Bot_ /)
      ij_ChannelHash(140,:) = (/ Str_ , AStr_ , ABot_ , Bot_ /)
      ij_ChannelHash(141,:) = (/ Up_ , AUp_ , AStr_ , Str_ /)
      ij_ChannelHash(142,:) = (/ Up_ , AUp_ , ABot_ , Bot_ /)
      ij_ChannelHash(143,:) = (/ Chm_ , AChm_ , ABot_ , Bot_ /)
      ij_ChannelHash(144,:) = (/ Dn_ , ADn_ , AChm_ , Chm_ /)
      ij_ChannelHash(145,:) = (/ Dn_ , ADn_ , AUp_ , Up_ /)
      ij_ChannelHash(146,:) = (/ Str_ , AStr_ , AChm_ , Chm_ /)
      ij_ChannelHash(147,:) = (/ Dn_ , AUp_ , AChm_ , Str_ /)
      ij_ChannelHash(148,:) = (/ Up_ , ADn_ , AStr_ , Chm_ /)
      ij_ChannelHash(149,:) = (/ Up_ , AUp_ , AUp_ , Up_ /)
      ij_ChannelHash(150,:) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
      ij_ChannelHash(151,:) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
      ij_ChannelHash(152,:) = (/ Str_ , AStr_ , AStr_ , Str_ /)
      ij_ChannelHash(153,:) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)
      ij_ChannelHash(154,:) = (/ Chm_ , AChm_ , AUp_ , Up_ /)
      ij_ChannelHash(155,:) = (/ Str_ , AStr_ , ADn_ , Dn_ /)
      ij_ChannelHash(156,:) = (/ Bot_ , ABot_ , ADn_ , Dn_ /)
      ij_ChannelHash(157,:) = (/ Bot_ , ABot_ , AStr_ , Str_ /)
      ij_ChannelHash(158,:) = (/ Str_ , AStr_ , AUp_ , Up_ /)
      ij_ChannelHash(159,:) = (/ Bot_ , ABot_ , AUp_ , Up_ /)
      ij_ChannelHash(160,:) = (/ Bot_ , ABot_ , AChm_ , Chm_ /)
      ij_ChannelHash(161,:) = (/ Chm_ , AChm_ , ADn_ , Dn_ /)
      ij_ChannelHash(162,:) = (/ Up_ , AUp_ , ADn_ , Dn_ /)
      ij_ChannelHash(163,:) = (/ Chm_ , AChm_ , AStr_ , Str_ /)
      ij_ChannelHash(164,:) = (/ Chm_ , AStr_ , ADn_ , Up_ /)
      ij_ChannelHash(165,:) = (/ Str_ , AChm_ , AUp_ , Dn_ /)
      ij_ChannelHash(166,:) = (/ Up_ , AUp_ , AUp_ , Up_ /)
      ij_ChannelHash(167,:) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
      ij_ChannelHash(168,:) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
      ij_ChannelHash(169,:) = (/ Str_ , AStr_ , AStr_ , Str_ /)
      ij_ChannelHash(170,:) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)

      ij_ChannelHash(171,:) = (/ AUp_ , Up_ , AChm_ , Chm_ /)
      ij_ChannelHash(172,:) = (/ ADn_ , Dn_ , AStr_ , Str_ /)
      ij_ChannelHash(173,:) = (/ ADn_ , Dn_ , ABot_ , Bot_ /)
      ij_ChannelHash(174,:) = (/ AStr_ , Str_ , ABot_ , Bot_ /)
      ij_ChannelHash(175,:) = (/ AUp_ , Up_ , AStr_ , Str_ /)
      ij_ChannelHash(176,:) = (/ AUp_ , Up_ , ABot_ , Bot_ /)
      ij_ChannelHash(177,:) = (/ AChm_ , Chm_ , ABot_ , Bot_ /)
      ij_ChannelHash(178,:) = (/ ADn_ , Dn_ , AChm_ , Chm_ /)
      ij_ChannelHash(179,:) = (/ ADn_ , Dn_ , AUp_ , Up_ /)
      ij_ChannelHash(180,:) = (/ AStr_ , Str_ , AChm_ , Chm_ /)
      ij_ChannelHash(181,:) = (/ AUp_ , Dn_ , AChm_ , Str_ /)
      ij_ChannelHash(182,:) = (/ ADn_ , Up_ , AStr_ , Chm_ /)
      ij_ChannelHash(183,:) = (/ AUp_ , Up_ , AUp_ , Up_ /)
      ij_ChannelHash(184,:) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
      ij_ChannelHash(185,:) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
      ij_ChannelHash(186,:) = (/ AStr_ , Str_ , AStr_ , Str_ /)
      ij_ChannelHash(187,:) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)
      ij_ChannelHash(188,:) = (/ AChm_ , Chm_ , AUp_ , Up_ /)
      ij_ChannelHash(189,:) = (/ AStr_ , Str_ , ADn_ , Dn_ /)
      ij_ChannelHash(190,:) = (/ ABot_ , Bot_ , ADn_ , Dn_ /)
      ij_ChannelHash(191,:) = (/ ABot_ , Bot_ , AStr_ , Str_ /)
      ij_ChannelHash(192,:) = (/ AStr_ , Str_ , AUp_ , Up_ /)
      ij_ChannelHash(193,:) = (/ ABot_ , Bot_ , AUp_ , Up_ /)
      ij_ChannelHash(194,:) = (/ ABot_ , Bot_ , AChm_ , Chm_ /)
      ij_ChannelHash(195,:) = (/ AChm_ , Chm_ , ADn_ , Dn_ /)
      ij_ChannelHash(196,:) = (/ AUp_ , Up_ , ADn_ , Dn_ /)
      ij_ChannelHash(197,:) = (/ AChm_ , Chm_ , AStr_ , Str_ /)
      ij_ChannelHash(198,:) = (/ AStr_ , Chm_ , ADn_ , Up_ /)
      ij_ChannelHash(199,:) = (/ AChm_ , Str_ , AUp_ , Dn_ /)
      ij_ChannelHash(200,:) = (/ AUp_ , Up_ , AUp_ , Up_ /)
      ij_ChannelHash(201,:) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
      ij_ChannelHash(202,:) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
      ij_ChannelHash(203,:) = (/ AStr_ , Str_ , AStr_ , Str_ /)
      ij_ChannelHash(204,:) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)

      NumberOfChannels = 204

      do ch=1, NumberOfChannels
      do i=1,4
        ij_ChannelHash(ch,i) = convertToPartIndex( ij_ChannelHash(ch,i) )
      enddo
      enddo

    endif


  return
  end subroutine


  subroutine get_NumberOfChannels(NumOfCh)
  implicit none
  integer, intent(out) :: NumOfCh

    NumOfCh = NumberOfChannels

  return
  end subroutine



  subroutine get_VBFoffshChannel(channel,i1,i2,i3,i4)
  implicit none
  integer, intent(in) :: channel
  integer, intent(out) :: i1,i2,i3,i4

    i1= ij_ChannelHash(channel,1)
    i2= ij_ChannelHash(channel,2)
    i3= ij_ChannelHash(channel,3)
    i4= ij_ChannelHash(channel,4)

  return
  end subroutine



  subroutine remove_VBFoffshChannel(channel)
  implicit none
  integer, intent(in) :: channel
  integer :: i

    do i =channel,NumberOfChannels-1
         ij_ChannelHash(i,:) = ij_ChannelHash(i+1,:)
    enddo
    NumberOfChannels = NumberOfChannels -1

  return
  end subroutine




  subroutine get_VBFoffshchannelHash(ijSel)
  implicit none
  integer, intent(out) :: ijSel(1:121,1:3)


      ijSel(  1,1:3) = (/ 2, 1, 1/)
      ijSel(  2,1:3) = (/ 1, 2, 1/)
      ijSel(  3,1:3) = (/ 2,-2, 1/)
      ijSel(  4,1:3) = (/-2, 2, 1/)
      ijSel(  5,1:3) = (/ 2, 3, 1/)
      ijSel(  6,1:3) = (/ 3, 2, 1/)
      ijSel(  7,1:3) = (/ 1,-1, 1/)
      ijSel(  8,1:3) = (/-1, 1, 1/)
      ijSel(  9,1:3) = (/ 1, 3, 1/)
      ijSel( 10,1:3) = (/ 3, 1, 1/)
      ijSel( 11,1:3) = (/ 2,-4, 1/)
      ijSel( 12,1:3) = (/-4, 2, 1/)
      ijSel( 13,1:3) = (/ 1, 4, 1/)
      ijSel( 14,1:3) = (/ 4, 1, 1/)

      ijSel( 15,1:3) = (/-5,-5, 1/)
      ijSel( 16,1:3) = (/-5,-4, 1/)
      ijSel( 17,1:3) = (/-5,-3, 1/)
      ijSel( 18,1:3) = (/-5,-2, 1/)
      ijSel( 19,1:3) = (/-5,-1, 1/)
      ijSel( 20,1:3) = (/-5, 1, 1/)
      ijSel( 21,1:3) = (/-5, 2, 1/)
      ijSel( 22,1:3) = (/-5, 3, 1/)
      ijSel( 23,1:3) = (/-5, 4, 1/)
      ijSel( 24,1:3) = (/-5, 5, 1/)
      ijSel( 25,1:3) = (/-4,-5, 1/)
      ijSel( 26,1:3) = (/-4,-4, 1/)
      ijSel( 27,1:3) = (/-4,-3, 1/)
      ijSel( 28,1:3) = (/-4,-2, 1/)
      ijSel( 29,1:3) = (/-4,-1, 1/)
      ijSel( 30,1:3) = (/-4, 1, 1/)
      ijSel( 31,1:3) = (/-4, 3, 1/)
      ijSel( 32,1:3) = (/-4, 4, 1/)
      ijSel( 33,1:3) = (/-4, 5, 1/)
      ijSel( 34,1:3) = (/-3,-5, 1/)
      ijSel( 35,1:3) = (/-3,-4, 1/)
      ijSel( 36,1:3) = (/-3,-3, 1/)
      ijSel( 37,1:3) = (/-3,-2, 1/)
      ijSel( 38,1:3) = (/-3,-1, 1/)
      ijSel( 39,1:3) = (/-3, 1, 1/)
      ijSel( 40,1:3) = (/-3, 2, 1/)
      ijSel( 41,1:3) = (/-3, 3, 1/)
      ijSel( 42,1:3) = (/-3, 4, 1/)
      ijSel( 43,1:3) = (/-3, 5, 1/)
      ijSel( 44,1:3) = (/-2,-5, 1/)
      ijSel( 45,1:3) = (/-2,-4, 1/)
      ijSel( 46,1:3) = (/-2,-3, 1/)
      ijSel( 47,1:3) = (/-2,-2, 1/)
      ijSel( 48,1:3) = (/-2,-1, 1/)
      ijSel( 49,1:3) = (/-2, 1, 1/)
      ijSel( 50,1:3) = (/-2, 3, 1/)
      ijSel( 51,1:3) = (/-2, 4, 1/)
      ijSel( 52,1:3) = (/-2, 5, 1/)
      ijSel( 53,1:3) = (/-1,-5, 1/)
      ijSel( 54,1:3) = (/-1,-4, 1/)
      ijSel( 55,1:3) = (/-1,-3, 1/)
      ijSel( 56,1:3) = (/-1,-2, 1/)
      ijSel( 57,1:3) = (/-1,-1, 1/)
      ijSel( 58,1:3) = (/-1, 2, 1/)
      ijSel( 59,1:3) = (/-1, 3, 1/)
      ijSel( 60,1:3) = (/-1, 4, 1/)
      ijSel( 61,1:3) = (/-1, 5, 1/)
      ijSel( 62,1:3) = (/ 1,-5, 1/)
      ijSel( 63,1:3) = (/ 1,-4, 1/)
      ijSel( 64,1:3) = (/ 1,-3, 1/)
      ijSel( 65,1:3) = (/ 1,-2, 1/)
      ijSel( 66,1:3) = (/ 1, 1, 1/)
      ijSel( 67,1:3) = (/ 1, 5, 1/)
      ijSel( 68,1:3) = (/ 2,-5, 1/)
      ijSel( 69,1:3) = (/ 2,-3, 1/)
      ijSel( 70,1:3) = (/ 2,-1, 1/)
      ijSel( 71,1:3) = (/ 2, 2, 1/)
      ijSel( 72,1:3) = (/ 2, 4, 1/)
      ijSel( 73,1:3) = (/ 2, 5, 1/)
      ijSel( 74,1:3) = (/ 3,-5, 1/)
      ijSel( 75,1:3) = (/ 3,-4, 1/)
      ijSel( 76,1:3) = (/ 3,-3, 1/)
      ijSel( 77,1:3) = (/ 3,-2, 1/)
      ijSel( 78,1:3) = (/ 3,-1, 1/)
      ijSel( 79,1:3) = (/ 3, 3, 1/)
      ijSel( 80,1:3) = (/ 3, 4, 1/)
      ijSel( 81,1:3) = (/ 3, 5, 1/)
      ijSel( 82,1:3) = (/ 4,-5, 1/)
      ijSel( 83,1:3) = (/ 4,-4, 1/)
      ijSel( 84,1:3) = (/ 4,-3, 1/)
      ijSel( 85,1:3) = (/ 4,-2, 1/)
      ijSel( 86,1:3) = (/ 4,-1, 1/)
      ijSel( 87,1:3) = (/ 4, 2, 1/)
      ijSel( 88,1:3) = (/ 4, 3, 1/)
      ijSel( 89,1:3) = (/ 4, 4, 1/)
      ijSel( 90,1:3) = (/ 4, 5, 1/)
      ijSel( 91,1:3) = (/ 5,-5, 1/)
      ijSel( 92,1:3) = (/ 5,-4, 1/)
      ijSel( 93,1:3) = (/ 5,-3, 1/)
      ijSel( 94,1:3) = (/ 5,-2, 1/)
      ijSel( 95,1:3) = (/ 5,-1, 1/)
      ijSel( 96,1:3) = (/ 5, 1, 1/)
      ijSel( 97,1:3) = (/ 5, 2, 1/)
      ijSel( 98,1:3) = (/ 5, 3, 1/)
      ijSel( 99,1:3) = (/ 5, 4, 1/)
      ijSel(100,1:3) = (/ 5, 5, 1/)

      ijSel(101:121,1:3) = 0

  return
  end subroutine





  !-- SM: |g2| = alphas/(six*pi), so multiply ME**2 by (|g2|*gs)**2
  !-- g3 not supported yet
  subroutine EvalAmp_SBFH_UnSymm_SA(p,res)
    real(dp), intent(in) :: p(4,5)
    real(dp), intent(out) :: res(-5:5,-5:5)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp, restmpid
    integer :: i, j

    res = zero

    call spinoru2(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    !-- gg -> gg
    call me2_ggggh(1,2,3,4,za,zb,sprod,restmp)
    restmp = restmp * avegg * SymmFac
    res(0,0) = res(0,0) + restmp

    !-- gg -> qqb
    call me2_qbqggh(4,3,1,2,za,zb,sprod,restmp)
    restmp = restmp * avegg
    res(0,0) = res(0,0) + restmp * nf

    !-- gq -> gq
    call me2_qbqggh(2,4,1,3,za,zb,sprod,restmp)
    restmp = restmp * aveqg
    do i = 1,5
       res(0,i) = res(0,i) + restmp
       res(0,-i) = res(0,-i) + restmp
    enddo
    call me2_qbqggh(1,4,2,3,za,zb,sprod,restmp)
    restmp = restmp * aveqg
    do i = 1,5
       res(i,0) = res(i,0) + restmp
       res(-i,0) = res(-i,0) + restmp
    enddo

    !-- qqb -> gg
    call me2_qbqggh(1,2,3,4,za,zb,sprod,restmp)
    restmp = restmp * aveqq * SymmFac
    do i = 1,5
       res(i,-i) = res(i,-i) + restmp
       res(-i,i) = res(-i,i) + restmp
    enddo

    !-- qqb -> qqb
    call me2_qbqQBQ(1,2,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq
    do i = 1,5
       res(i,-i) =res(i,-i) + restmp
    enddo
    call me2_qbqQBQ(2,1,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq
    do i = 1,5
       res(-i,i) =res(-i,i) + restmp
    enddo

    !-- qqb -> rrb
    call me2_qbqQBQ(1,2,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       res(i,-i) =res(i,-i) + restmp * (nf-1.0_dp)
       res(-i,i) = res(-i,i) + restmp * (nf-1.0_dp)
    enddo

    !-- qrb -> qrb
    call me2_qbqQBQ(1,3,4,2,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       do j = 1,5
          if (i.ne.j) res(i,-j) = res(i,-j) + restmp
       enddo
    enddo
    call me2_qbqQBQ(2,3,4,1,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       do j = 1,5
          if (i.ne.j) res(-j,i) = res(-j,i) + restmp
       enddo
    enddo

    !-- qq -> qq
    call me2_qbqQBQ(1,3,2,4,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq * SymmFac
    do i = 1,5
       res(i,i) = res(i,i) + restmp
       res(-i,-i) = res(-i,-i) + restmp
    enddo

    !-- qr -> qr
    call me2_qbqQBQ(1,3,2,4,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       do j = 1,5
          if (i.ne.j) then
             res(i,j) = res(i,j) + restmp
             res(-i,-j) = res(-i,-j) + restmp
          endif
       enddo
    enddo

    res(:,:) = res(:,:) * (2d0/3d0*alphas**2)**2
    return

  end subroutine EvalAmp_SBFH_UnSymm_SA





  !-- SM: |g2| = alphas/(six*pi), so multiply ME**2 by (|g2|*gs)**2
  !-- g3 not supported yet
  subroutine EvalAmp_SBFH_UnSymm_SA_Select(p,iSel,jSel,flav_tag,iflip,res)
    real(dp), intent(in) :: p(4,5)
    integer, intent(in) :: iSel,jSel,flav_tag
    real(dp), intent(out) :: res(-5:5,-5:5)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp, restmpid
    integer :: i, j,iflip

    res = zero

    call spinoru2(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)
    iflip = 1

    !-- gg -> gg
    if( iSel.eq.pdfGlu_ .and. jSel.eq.pdfGlu_ .and. flav_tag.eq.1 ) then
        call me2_ggggh(1,2,3,4,za,zb,sprod,restmp)
        restmp = restmp * avegg * SymmFac
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif


    !-- gg -> qqb
    if( iSel.eq.pdfGlu_ .and. jSel.eq.pdfGlu_ .and. flav_tag.eq.2 ) then
        call me2_qbqggh(4,3,1,2,za,zb,sprod,restmp)
        restmp = restmp * avegg
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp * nf
        return
    endif



    !-- gq -> gq
    if( iSel.eq.pdfGlu_ .and. jSel.ne.0 ) then
        call me2_qbqggh(2,4,1,3,za,zb,sprod,restmp)
        restmp = restmp * aveqg
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif

    if( jSel.eq.pdfGlu_ .and. iSel.ne.0 ) then
        iflip = 2
        call me2_qbqggh(1,4,2,3,za,zb,sprod,restmp)
        restmp = restmp * aveqg
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif



    !-- qqb -> gg
    if( iSel.ne.0 .and. jSel.eq.-iSel .and. flav_tag.eq.1 ) then
        call me2_qbqggh(1,2,3,4,za,zb,sprod,restmp)
        restmp = restmp * aveqq * SymmFac
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif


    !-- qqb -> qqb
    if( iSel.gt.0 .and. jSel.eq.-iSel .and. flav_tag.eq.2 ) then
        call me2_qbqQBQ(1,2,4,3,za,zb,sprod,restmp,restmpid)
        restmp = restmpid * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif

    if( iSel.lt.0 .and. jSel.eq.-iSel .and. flav_tag.eq.2 ) then
        iflip = 2
        call me2_qbqQBQ(2,1,4,3,za,zb,sprod,restmp,restmpid)
        restmp = restmpid * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif


    !-- qqb -> rrb
    if( iSel.ne.0 .and. jSel.eq.-iSel .and. flav_tag.eq.3 ) then
        call me2_qbqQBQ(1,2,4,3,za,zb,sprod,restmp,restmpid)
        restmp = restmp * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp * (nf-1.0_dp)
        return
    endif



    !-- qrb -> qrb
    if( iSel.gt.0 .and. jSel.lt.0 .and. iSel.ne.jSel ) then
        call me2_qbqQBQ(1,3,4,2,za,zb,sprod,restmp,restmpid)
        restmp = restmp * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif

    if( iSel.lt.0 .and. jSel.gt.0 .and. iSel.ne.jSel ) then
        iflip = 2
        call me2_qbqQBQ(2,3,4,1,za,zb,sprod,restmp,restmpid)
        restmp = restmp * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif




    !-- qq -> qq
    if( iSel.ne.0 .and. iSel.eq.jSel ) then
        call me2_qbqQBQ(1,3,2,4,za,zb,sprod,restmp,restmpid)
        restmp = restmpid * aveqq * SymmFac
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif



    !-- qr -> qr
    if( iSel.gt.0 .and. jSel.gt.0 .and. iSel.ne.jSel ) then
        call me2_qbqQBQ(1,3,2,4,za,zb,sprod,restmp,restmpid)
        restmp = restmp * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif

    if( iSel.lt.0 .and. jSel.lt.0 .and. iSel.ne.jSel ) then
        call me2_qbqQBQ(1,3,2,4,za,zb,sprod,restmp,restmpid)
        restmp = restmp * aveqq
        restmp = restmp * (2d0/3d0*alphas**2)**2
        res(iSel,jSel) = restmp
        return
    endif



  return

  end subroutine


    subroutine EvalAmp_SBFH_UnSymm_SA_Select_exact(p,iSel,jSel,rSel,sSel,res)
    real(dp), intent(in) :: p(4,5)
    integer, intent(in) :: iSel,jSel,rSel,sSel!,flavor_tag ! flavor_tag for TEST
    real(dp), intent(out) :: res
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp, restmpid
    integer :: i, j
	integer :: j1,j2,k1,k2
	logical :: isGGGG,isQQGG,isQQQQ,isQQQQ_idQ

	restmp=0.0_dp
	restmpid=0.0_dp
	res=0.0_dp
	isGGGG=.false.
	isQQGG=.false.
	isQQQQ=.false.
	isQQQQ_idQ=.false.
	j1 = 0
	j2 = 0
	k1 = 0
	k2 = 0
	if( (abs(iSel).eq.pdfTop_ .or. abs(jSel).eq.pdfTop_) .or. (abs(rSel).eq.pdfTop_ .or. abs(sSel).eq.pdfTop_) ) return

	!print *, "Begin EvalAmp_SBFH_UnSymm_SA_Select_exact"
	!print *, "iSel: ",iSel,", jSel: ",jSel," rSel: ",rSel," sSel: ",sSel

	if(iSel.eq.pdfGlu_) then ! 1==g
	   if(jSel.eq.pdfGlu_) then ! gg
	      if(rSel.eq.pdfGlu_ .and. sSel .eq. pdfGlu_) then ! gg->gg
			 isGGGG=.true.
		    j1=1
          j2=2
			 k1=3
			 k2=4

	      elseif(sSel.eq.(-rSel)) then ! gg->qqb/qbq
			 isQQGG=.true.
		    if(rSel.gt.0) then ! gg->qqb
		       j1=4
			    j2=3
		    else ! gg->qbq
		       j1=3
			    j2=4
			 endif
		       k1=1
	          k2=2
		  endif

	   elseif(abs(jSel).eq.pdfUp_ .or. abs(jSel).eq.pdfChm_ .or. abs(jSel).eq.pdfDn_ .or. abs(jSel).eq.pdfStr_ .or. abs(jSel).eq.pdfBot_) then ! gq/gqb
		   isQQGG=.true.
	      k1=1
	      if(jSel.gt.0) then ! gq
		     j1=2
		     if(rSel.eq.iSel .and. sSel.eq.jSel) then ! gq->gq
		        k2=3
				  j2=4
			 elseif(sSel.eq.iSel .and. rSel.eq.jSel) then ! gq->qg
			     j2=3
				  k2=4
			 endif
	      else ! gqb
			 j2=2
		    if(rSel.eq.iSel .and. sSel.eq.jSel) then ! gqb->gqb
		        k2=3
				  j1=4
			 elseif(sSel.eq.iSel .and. rSel.eq.jSel) then ! gqb->qbg
			     j1=3
				  k2=4
			 endif
	      endif
	   endif
	elseif(abs(iSel).eq.pdfUp_ .or. abs(iSel).eq.pdfChm_ .or. abs(iSel).eq.pdfDn_ .or. abs(iSel).eq.pdfStr_ .or. abs(iSel).eq.pdfBot_) then ! q/qb?->?
	   if(jSel.eq.pdfGlu_) then ! qg->?
	      isQQGG=.true.
	      k1=2
	      if(iSel.gt.0) then
		     j1=1
		     if(rSel.eq.iSel .and. sSel.eq.jSel) then ! qg->qg
		        k2=4
				j2=3
			 elseif(sSel.eq.iSel .and. rSel.eq.jSel) then ! qg->gq
			    j2=4
				k2=3
			 endif
	      else ! qbg
		     j2=1
		     if(rSel.eq.iSel .and. sSel.eq.jSel) then ! qbg->qbg
		        k2=4
				  j1=3
			 elseif(sSel.eq.iSel .and. rSel.eq.jSel) then ! qbg->gqb
			    j1=4
				 k2=3
			 endif
	      endif

	   elseif(abs(jSel).eq.pdfUp_ .or. abs(jSel).eq.pdfChm_ .or. abs(jSel).eq.pdfDn_ .or. abs(jSel).eq.pdfStr_ .or. abs(jSel).eq.pdfBot_) then ! qq/qqb/qbq/qbqb->?
	      if(iSel.gt.0) then ! qq'/qqb'->?
		     j1=1
	      else ! qbq'/qbqb'->?
		     j2=1
	      endif

		  if(jSel.eq.(-iSel)) then ! qqb/qbq->?
	         if(jSel.gt.0) then ! qbq->?
		        j1=2
	         else ! qqb->?
			    j2=2
	         endif

	         if(rSel.eq.pdfGlu_ .and. sSel .eq. pdfGlu_) then ! qqb/qbq->gg
	            isQQGG=.true.
	            k1=3
	            k2=4
	         elseif(sSel.eq.(-rSel)) then ! qqb->q'qb'/qb'q'
	            isQQQQ=.true.
	            if(rSel.gt.0) then ! ->q'qb'
	               if(rSel .eq. iSel .or. rSel.eq.jSel) isQQQQ_idQ=.true.
		           k1=4
			       k2=3
		        else ! ->qb'q'
	               if(rSel .eq. iSel .or. rSel.eq.jSel) isQQQQ_idQ=.true.
		           k1=3
			       k2=4
	            endif
	         endif

		  elseif(jSel.eq.iSel) then ! qq/qbqb->?
	         isQQQQ=.true.
	         if(jSel.gt.0) then ! qq
		        k1=2
	         else ! qbqb
			    k2=2
	         endif
	         if(rSel.gt.0 .and. rSel.eq.sSel .and. rSel.eq.iSel) then ! qq->qq
			    isQQQQ_idQ=.true.
		        j2=3
				k2=4
	         elseif(rSel.lt.0 .and. rSel.eq.sSel .and. rSel.eq.iSel) then ! qbqb->qbqb
			    isQQQQ_idQ=.true.
			    j1=3
				k1=4
	         endif

		  elseif(jSel.eq.sign(jSel,iSel)) then ! qq'/qbqb'->?
	         isQQQQ=.true.
	         if(jSel.gt.0) then ! qq'
		        k1=2 ! j1=1
	         else ! qbqb'
			    k2=2 ! j2=1
	         endif
	         if(rSel.gt.0 .and. rSel.ne.sSel .and. rSel.eq.sign(rSel,sSel) .and. rSel.eq.iSel) then ! ->qq'
		        j2=3
				k2=4
	         elseif(rSel.lt.0 .and. rSel.ne.sSel .and. rSel.eq.sign(rSel,sSel) .and. rSel.eq.iSel) then ! ->qbqb'
			    j1=3
				k1=4
	         elseif(rSel.gt.0 .and. rSel.ne.sSel .and. rSel.eq.sign(rSel,sSel) .and. rSel.eq.jSel) then ! ->q'q
		        j2=4
				k2=3
	         elseif(rSel.lt.0 .and. rSel.ne.sSel .and. rSel.eq.sign(rSel,sSel) .and. rSel.eq.jSel) then ! ->qb'qb
			    j1=4
				k1=3
	         endif

		  elseif(jSel.eq.(-sign(jSel,iSel))) then ! qqb'/qbq'->?
	         isQQQQ=.true.
	         if(jSel.gt.0) then ! qbq'
		        k1=2 ! j2=1
	         else ! qqb'
			    k2=2 ! j1=1
	         endif
	         if(rSel.gt.0 .and. rSel.ne.sSel .and. rSel.ne.sign(rSel,sSel) .and. rSel.eq.iSel) then ! qqb'->qqb'
		        j2=3
				k1=4
	         elseif(rSel.lt.0 .and. rSel.ne.sSel .and. rSel.ne.sign(rSel,sSel) .and. rSel.eq.iSel) then ! qbq'->qbq'
			    j1=3
				k2=4
	         elseif(rSel.lt.0 .and. rSel.ne.sSel .and. rSel.ne.sign(rSel,sSel) .and. rSel.eq.jSel) then ! qqb'->qb'q
		        j2=4
				k1=3
	         elseif(rSel.gt.0 .and. rSel.ne.sSel .and. rSel.ne.sign(rSel,sSel) .and. rSel.eq.jSel) then ! qbq'->q'qb
			    j1=4
				k2=3
	         endif
	      endif

	   endif
	endif

	if(j1.eq.0 .or. j2.eq.0 .or. k1.eq.0 .or. k2.eq.0) then
	   print *,"Could not recognize the incoming flavors!"
	   print *,"iSel: ",iSel,", jSel: ",jSel," rSel: ",rSel," sSel: ",sSel
	   print *,"j1: ",j1,", j2: ",j2," k1: ",k1," k2: ",k2
		if(isGGGG) then
			print *,"is gggg"
		endif
		if(isQQGG) then
			print *,"is qbqgg"
		endif
		if(isQQQQ) then
			print *,"is qbqQBQ"
		endif
		if(isQQQQ_idQ) then
			print *,"is qbqQBQ id"
		endif
	endif

    call spinoru2(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

	if(isGGGG) then
	   call me2_ggggh(j1,j2,k1,k2,za,zb,sprod,restmp)
	   restmp = restmp * avegg * SymmFac
	elseif(isQQQQ) then
        call me2_qbqQBQ(j1,j2,k1,k2,za,zb,sprod,restmp,restmpid)
        if(isQQQQ_idQ) then
		   restmp = restmpid * aveqq
		   if(iSel.eq.jSel) restmp = restmp * SymmFac
        else
		   restmp = restmp * aveqq
		   if(iSel .eq. (-jSel) .and. rSel .eq. (-sSel) .and. abs(iSel) .ne. abs(rSel)) restmp = restmp * (nf-1.0_dp)
		endif
	elseif(isQQGG) then
	    call me2_qbqggh(j1,j2,k1,k2,za,zb,sprod,restmp)
	    if( iSel.eq.pdfGlu_ .and. jSel.eq.pdfGlu_) then
	       restmp = restmp * avegg * nf
	    elseif( iSel.ne.pdfGlu_ .and. jSel.ne.pdfGlu_) then
	       restmp = restmp * aveqq * SymmFac
	    else
	       restmp = restmp * aveqg
		endif
	endif

   restmp = restmp * (2d0/3d0*alphas**2)**2
	res = restmp
	return
  end subroutine



  subroutine EvalAmp_WBFH_UnSymm_SA(p,res)
    real(dp), intent(in) :: p(4,5)
    real(dp), intent(out) :: res(-5:5,-5:5)
    complex(dp) :: amp_z(-1:1,-1:1), amp_z_b(-1:1,-1:1)
    complex(dp) :: amp_w(-1:1,-1:1)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: Lu, Ru
    real(dp) :: Ld, Rd
    real(dp) :: couplz
    real(dp) :: couplw
    real(dp) :: restmp=0d0
    integer :: i, j, j1, j2, iflip, pdfindex(2)

    Lu = aL_QUp**2
    Ru = aR_QUp**2
    Ld = aL_QDn**2
    Rd = aR_QDn**2
    couplz = gwsq * xw/twosc**2
    couplw = gwsq/two

    res = zero

    call spinoru2(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    !-- qq->qq, up
    amp_z = A0_VV_4f(4,1,3,2,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(3,1,4,2,za,zb,sprod,m_z,ga_z)

! !     adding contract terms
!     iprop12 = sprod(4,1) - mv**2 + ci * mv * ga_v
!     iprop34 = sprod(3,2) - mv**2 + ci * mv * ga_v
!
!     amp_z(-1,-1) = amp_z(-1,-1) + amp_z(-1,-1)*iprop12/(aL_QUp*couplz)*ehz_L_U  &
!                                 + amp_z(-1,-1)*iprop34/(aL_QUp*couplz)*ehz_L_U  &
!                                 + amp_z(-1,-1)*iprop12/(aL_QUp*couplz)*ehz_L_U  &
!                                               *iprop34/(aL_QUp*couplz)*ehz_L_U
!
!     amp_z(+1,-1) = amp_z(+1,-1) + amp_z(+1,-1)*iprop12/(aR_QUp*couplz)*ehz_R_U  &
!                                 + amp_z(+1,-1)*iprop34/(aL_QUp*couplz)*ehz_L_U  &
!                                 + amp_z(+1,-1)*iprop12/(aR_QUp*couplz)*ehz_R_U  &
!                                               *iprop34/(aL_QUp*couplz)*ehz_L_U
!
!     amp_z(-1,+1) = amp_z(-1,+1) + amp_z(-1,+1)*iprop12/(aL_QUp*couplz)*ehz_L_U  &
!                                 + amp_z(-1,+1)*iprop34/(aR_QUp*couplz)*ehz_R_U  &
!                                 + amp_z(-1,+1)*iprop12/(aL_QUp*couplz)*ehz_L_U  &
!                                               *iprop34/(aR_QUp*couplz)*ehz_R_U
!
!     amp_z(+1,+1) = amp_z(+1,+1) + amp_z(+1,+1)*iprop12/(aR_QUp*couplz)*ehz_R_U  &
!                                 + amp_z(+1,+1)*iprop34/(aR_QUp*couplz)*ehz_R_U  &
!                                 + amp_z(+1,+1)*iprop12/(aR_QUp*couplz)*ehz_R_U  &
!                                               *iprop34/(aR_QUp*couplz)*ehz_R_U


!     Lu = aL_QUp + iprop12/couplz*ehz_L_U
!     Ru = aR_QUp + iprop12/couplz*ehz_R_U
!     amp_z(-1,-1) = amp_z(-1,-1)*Lu*Lu
!     amp_z(-1,+1) = amp_z(-1,+1)*Lu*Ru
!     amp_z(+1,-1) = amp_z(+1,-1)*Ru*Lu
!     amp_z(+1,+1) = amp_z(+1,+1)*Ru*Ru
!
!     amp_z_b(-1,-1) = amp_z_b(-1,-1)*Lu*Lu
!     amp_z_b(-1,+1) = amp_z_b(-1,+1)*Lu*Ru
!     amp_z_b(+1,-1) = amp_z_b(+1,-1)*Ru*Lu
!     amp_z_b(+1,+1) = amp_z_b(+1,+1)*Ru*Ru

    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Lu**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Lu * Ru + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Lu * Ru + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Ru**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Lu**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Ru**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfUp_,pdfUp_) = restmp
    res(pdfChm_,pdfChm_) = restmp





    !-- qq->qq, down
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Ld**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Ld * Rd + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Ld * Rd + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Rd**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Ld**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Rd**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfDn_,pdfDn_) = restmp
    res(pdfStr_,pdfStr_) = restmp
    res(pdfBot_,pdfBot_) = restmp * tagbot

    !-- qbqb->qbqb, aup
    amp_z = A0_VV_4f(1,4,2,3,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(1,3,2,4,za,zb,sprod,m_z,ga_z)

    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Lu**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Lu * Ru + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Lu * Ru + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Ru**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Lu**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Ru**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfAUp_,pdfAUp_) = restmp
    res(pdfAChm_,pdfAChm_) = restmp

    !-- qbqb->qbqb, adn
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Ld**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Ld * Rd + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Ld * Rd + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Rd**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Ld**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Rd**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfADn_,pdfADn_) = restmp
    res(pdfAStr_,pdfAStr_) = restmp
    res(pdfABot_,pdfABot_) = restmp * tagbot

    !-- W/Z interference
    j1 = 1
    j2 = 2
    do iflip = 1, 2
       !-- ud -> ud
       amp_z = A0_VV_4f(4,j2,3,j1,za,zb,sprod,m_z,ga_z)
       amp_w = -A0_VV_4f(4,j1,3,j2,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)

       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
            (abs(amp_z(-1,+1))**2) * Ld * Ru + &
            (abs(amp_z(+1,-1))**2) * Rd * Lu + &
            (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2

       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * xn**2

       restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
            aL_QUp * aL_QDn * couplz * couplw * xn

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfUp_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !-- ubdb -> ubdb
       amp_z = A0_VV_4f(j2,4,j1,3,za,zb,sprod,m_z,ga_z)
       amp_w = -A0_VV_4f(j1,4,j2,3,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)

       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
            (abs(amp_z(-1,+1))**2) * Ld * Ru + &
            (abs(amp_z(+1,-1))**2) * Rd * Lu + &
            (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2

       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * xn**2

       restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
            aL_QUp * aL_QDn * couplz * couplw * xn

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfAUp_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfAChm_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       j1 = 2
       j2 = 1
    enddo

    !-- non-interfering diagrams below
    j1 = 1
    j2 = 2
    do iflip = 1,2

       !-- qqb processes

       !--uub -> uub // ddb
       amp_z = A0_VV_4f(3,j1,j2,4,za,zb,sprod,m_z,ga_z)
       amp_w = A0_VV_4f(3,j1,j2,4,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)

       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg * tag1
       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2

       pdfindex = flip(iflip,pdfUp_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !--udb -> udb
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg

       pdfindex = flip(iflip,pdfUp_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       !--dub -> dub
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
            (abs(amp_z(-1,+1))**2) * Ld * Ru + &
            (abs(amp_z(+1,-1))**2) * Rd * Lu + &
            (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * SpinAvg

       pdfindex = flip(iflip,pdfDn_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfStr_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfDn_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfStr_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfBot_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfBot_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       !--ddb -> uub/ddb
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg * tag1

       pdfindex = flip(iflip,pdfBot_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfDn_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfStr_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfBot_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfBot_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2

       pdfindex = flip(iflip,pdfDn_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfStr_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfDn_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfStr_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !-- non-symmetric qq processes
       amp_z = A0_VV_4f(3,j1,4,j2,za,zb,sprod,m_z,ga_z)
       amp_w = A0_VV_4f(3,j2,4,j1,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)

       !--uc -> uc
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg

       pdfindex = flip(iflip,pdfUp_,pdfChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !--us -> us/cd
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1

       pdfindex = flip(iflip,pdfUp_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfChm_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2

       pdfindex = flip(iflip,pdfUp_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !--ds -> ds
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg

       pdfindex = flip(iflip,pdfDn_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfDn_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfStr_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       !-- qbqb processes
       amp_z = A0_VV_4f(j1,3,j2,4,za,zb,sprod,m_z,ga_z)
       amp_w = A0_VV_4f(j1,4,j2,3,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)

       !--ubcb -> ubcb
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg

       pdfindex = flip(iflip,pdfAUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !--ubsb -> ubsb//cbdb
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1

       pdfindex = flip(iflip,pdfAUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfAChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2

       pdfindex = flip(iflip,pdfAUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfAChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !--dbsb -> dbsb
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg

       pdfindex = flip(iflip,pdfADn_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfADn_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfAStr_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       j1 = 2
       j2 = 1

    enddo

    return

  end subroutine EvalAmp_WBFH_UnSymm_SA




  SUBROUTINE EvalAmp_WBFH_UnSymm_SA_Select(p,iSel,jSel,zz_fusion,iflip,res)
    implicit none
    real(dp), intent(in) :: p(4,5)
    real(dp), intent(out) :: res(-5:5,-5:5)
    logical, intent(in) :: zz_fusion
    integer, intent(in) :: iSel,jSel
    complex(dp) :: amp_z(-1:1,-1:1), amp_z_b(-1:1,-1:1)
    complex(dp) :: amp_w(-1:1,-1:1)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp
    integer :: i, j, j1, j2, iflip, pdfindex(2)

    call spinoru2(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)
    iflip = 1

    !-- qq->qq, up
    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfUp_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfChm_) ) then

       amp_z = A0_ZZ_4f(4,1,3,2,za,zb,sprod,iSel,jSel)
       amp_z_b = -A0_ZZ_4f(3,1,4,2,za,zb,sprod,iSel,jSel)

       restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2)) * xn**2

       restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp)) * xn

       restmp = restmp * SymmFac * aveqq

       if( .not. zz_fusion ) restmp = 0.0d0

       res(pdfUp_,pdfUp_) = restmp
       res(pdfChm_,pdfChm_) = restmp

       return
    endif


   !-- qq->qq, down
    if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfDn_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfStr_) .or. (iSel.eq.pdfBot_ .and. jSel.eq.pdfBot_) ) then
       amp_z = A0_ZZ_4f(4,1,3,2,za,zb,sprod,iSel,jSel)
       amp_z_b = -A0_ZZ_4f(3,1,4,2,za,zb,sprod,iSel,jSel)

       restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2)) * xn**2

       restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp)) * xn

       restmp = restmp * SymmFac * aveqq

       if( .not. zz_fusion ) restmp = 0.0d0

       res(pdfDn_,pdfDn_) = restmp
       res(pdfStr_,pdfStr_) = restmp
       res(pdfBot_,pdfBot_) = restmp * tagbot
       return
    endif


    !-- qbqb->qbqb, aup
    if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfAUp_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfAChm_) ) then
       amp_z = A0_ZZ_4f(1,4,2,3,za,zb,sprod,iSel,jSel)
       amp_z_b = -A0_ZZ_4f(1,3,2,4,za,zb,sprod,iSel,jSel)

       restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2)) * xn**2

       restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp)) * xn

       restmp = restmp * SymmFac * aveqq

       if( .not. zz_fusion ) restmp = 0.0d0

       res(pdfAUp_,pdfAUp_) = restmp
       res(pdfAChm_,pdfAChm_) = restmp
       return
    endif


    !-- qbqb->qbqb, adn
    if( (iSel.eq.pdfADn_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfAStr_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfABot_ .and. jSel.eq.pdfABot_) ) then
       amp_z = A0_ZZ_4f(1,4,2,3,za,zb,sprod,iSel,jSel)
       amp_z_b = -A0_ZZ_4f(1,3,2,4,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2)) * xn**2

       restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp)) * xn

       restmp = restmp * SymmFac * aveqq

       if( .not. zz_fusion ) restmp = 0.0d0

       res(pdfADn_,pdfADn_) = restmp
       res(pdfAStr_,pdfAStr_) = restmp
       res(pdfABot_,pdfABot_) = restmp * tagbot
       return
    endif



    ! NEED TO ADJSUT J1,J2 in W AMPS HERE
    !-- W/Z interference
    j1 = 1
    j2 = 2
    iflip = 1

    !-- ud -> ud
    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfDn_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfStr_) ) then
       amp_z = A0_ZZ_4f(4,j2,3,j1,za,zb,sprod,jSel,iSel)
       amp_w = -A0_WW_4f(4,j1,3,j2,za,zb,sprod,iSel,jSel,useWWcoupl=.true.,Wpm_flip=.true.)
       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * xn**2
       else
          restmp= ( &
                   abs(amp_w(-1,-1))**2 + &
                   abs(amp_w(-1,+1))**2 + &
                   abs(amp_w(+1,-1))**2 + &
                   abs(amp_w(+1,+1))**2   &
                  ) * xn**2
          restmp = restmp + two * ( &
             real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) + &
             real(amp_z(-1,+1)*conjg(amp_w(-1,+1)),kind=dp) + &
             real(amp_z(+1,-1)*conjg(amp_w(+1,-1)),kind=dp) + &
             real(amp_z(+1,+1)*conjg(amp_w(+1,+1)),kind=dp)   &
             ) * xn

       endif

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfUp_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !-- ubdb -> ubdb
    if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfAStr_) ) then
       amp_z = A0_ZZ_4f(j2,4,j1,3,za,zb,sprod,jSel,iSel)
       amp_w = -A0_WW_4f(j1,4,j2,3,za,zb,sprod,iSel,jSel,useWWcoupl=.true.,Wpm_flip=.false.)

       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * xn**2
       else
          restmp= ( &
                   abs(amp_w(-1,-1))**2 + &
                   abs(amp_w(-1,+1))**2 + &
                   abs(amp_w(+1,-1))**2 + &
                   abs(amp_w(+1,+1))**2   &
                  ) * xn**2
          restmp = restmp + two * ( &
             real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) + &
             real(amp_z(-1,+1)*conjg(amp_w(-1,+1)),kind=dp) + &
             real(amp_z(+1,-1)*conjg(amp_w(+1,-1)),kind=dp) + &
             real(amp_z(+1,+1)*conjg(amp_w(+1,+1)),kind=dp)   &
             ) * xn

       endif

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfAUp_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfAChm_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    ! NEED TO ADJSUT J1,J2 in W AMPS HERE
    !-- W/Z interference
    j1 = 2
    j2 = 1
    iflip = 2

    !-- ud -> ud
    if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfUp_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfChm_) ) then
       amp_z = A0_ZZ_4f(4,j2,3,j1,za,zb,sprod,iSel,jSel)
       amp_w = -A0_WW_4f(4,j1,3,j2,za,zb,sprod,jSel,iSel,useWWcoupl=.true.,Wpm_flip=.true.)
       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * xn**2
       else
          restmp= ( &
                   abs(amp_w(-1,-1))**2 + &
                   abs(amp_w(-1,+1))**2 + &
                   abs(amp_w(+1,-1))**2 + &
                   abs(amp_w(+1,+1))**2   &
                  ) * xn**2
          restmp = restmp + two * ( &
             real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) + &
             real(amp_z(-1,+1)*conjg(amp_w(-1,+1)),kind=dp) + &
             real(amp_z(+1,-1)*conjg(amp_w(+1,-1)),kind=dp) + &
             real(amp_z(+1,+1)*conjg(amp_w(+1,+1)),kind=dp)   &
             ) * xn

       endif

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfUp_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !-- ubdb -> ubdb
    if( (iSel.eq.pdfADn_ .and. jSel.eq.pdfAUp_) .or. (iSel.eq.pdfAStr_ .and. jSel.eq.pdfAChm_) ) then
       amp_z = A0_ZZ_4f(j2,4,j1,3,za,zb,sprod,iSel,jSel)
       amp_w = -A0_WW_4f(j1,4,j2,3,za,zb,sprod,jSel,iSel,useWWcoupl=.true.,Wpm_flip=.false.)

       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * xn**2
       else
          restmp= ( &
                   abs(amp_w(-1,-1))**2 + &
                   abs(amp_w(-1,+1))**2 + &
                   abs(amp_w(+1,-1))**2 + &
                   abs(amp_w(+1,+1))**2   &
                  ) * xn**2
          restmp = restmp + two * ( &
             real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) + &
             real(amp_z(-1,+1)*conjg(amp_w(-1,+1)),kind=dp) + &
             real(amp_z(+1,-1)*conjg(amp_w(+1,-1)),kind=dp) + &
             real(amp_z(+1,+1)*conjg(amp_w(+1,+1)),kind=dp)   &
             ) * xn
       endif

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfAUp_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfAChm_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !-- non-interfering diagrams below
    j1 = 1
    j2 = 2
    iflip = 1

    !-- qqb processes

    !--uub -> uub // ddb
    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfAUp_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfAChm_) .or. (iSel.eq.pdfUp_ .and. jSel.eq.pdfAChm_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfAUp_) ) then

       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,iSel,jSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(4,j1,j2,3,za,zb,sprod,iSel,jSel,useWWcoupl=.true.,Wpm_flip=.true.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfUp_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !--udb -> udb
    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfUp_ .and. jSel.eq.pdfAStr_) .or. &
         (iSel.eq.pdfChm_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfUp_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfABot_) ) then
       amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfUp_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    !--dub -> dub
    if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfAUp_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfAChm_)  .or. (iSel.eq.pdfDn_ .and. jSel.eq.pdfAChm_)  .or. &
         (iSel.eq.pdfStr_ .and. jSel.eq.pdfAUp_)  .or. (iSel.eq.pdfBot_ .and. jSel.eq.pdfAUp_)  .or. (iSel.eq.pdfBot_ .and. jSel.eq.pdfAChm_) ) then
       amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfDn_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfDn_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfBot_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfBot_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    !--ddb -> uub/ddb
    if( (iSel.eq.pdfBot_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfDn_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfABot_) .or.&
         (iSel.eq.pdfBot_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfBot_ .and. jSel.eq.pdfAStr_) ) then
       amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,iSel,jSel)

       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfBot_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfDn_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfStr_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfBot_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfBot_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfAStr_)  .or. (iSel.eq.pdfDn_ .and. jSel.eq.pdfAStr_)  .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfADn_) ) then
       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,iSel,jSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(4,j1,j2,3,za,zb,sprod,iSel,jSel,useWWcoupl=.true.,Wpm_flip=.false.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfDn_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfDn_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfChm_)  ) then

       !-- non-symmetric qq processes
       amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,iSel,jSel)
       !--uc -> uc
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfUp_,pdfChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !--us -> us/cd
    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfBot_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfBot_) ) then

       amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfUp_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfChm_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfStr_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfDn_) ) then

       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,iSel,jSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(3,j2,4,j1,za,zb,sprod,jSel,iSel,useWWcoupl=.true.,Wpm_flip=.false.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfUp_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !--ds -> ds
    if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfStr_) .or. (iSel.eq.pdfDn_ .and. jSel.eq.pdfBot_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfBot_) ) then

       amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfDn_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfDn_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfStr_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    !-- qbqb processes
    if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfAChm_) ) then
       amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,iSel,jSel)
       !--ubcb -> ubcb
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfAUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       return
    endif


    !--ubsb -> ubsb//cbdb
    if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfABot_) ) then
       amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfAUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfAChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfADn_) ) then

       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,iSel,jSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(j1,4,j2,3,za,zb,sprod,iSel,jSel,useWWcoupl=.true.,Wpm_flip=.false.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfAUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfAChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp


       return
    endif


    !--dbsb -> dbsb
    if( (iSel.eq.pdfADn_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfADn_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfAStr_ .and. jSel.eq.pdfABot_) ) then

       amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,iSel,jSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfADn_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfADn_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfAStr_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif

    !-------------

    j1 = 2
    j2 = 1
    iflip = 2
    !-- qqb processes

    !--uub -> uub // ddb
    if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfAUp_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfAChm_) .or. (jSel.eq.pdfUp_ .and. iSel.eq.pdfAChm_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfAUp_) ) then

       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,jSel,iSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(4,j1,j2,3,za,zb,sprod,jSel,iSel,useWWcoupl=.true.,Wpm_flip=.true.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfUp_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp


       return
    endif


    !--udb -> udb
    if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfADn_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfAStr_) .or. (jSel.eq.pdfUp_ .and. iSel.eq.pdfAStr_) .or. &
         (jSel.eq.pdfChm_ .and. iSel.eq.pdfADn_) .or. (jSel.eq.pdfUp_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfABot_) ) then
       amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,jSel,iSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfUp_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    !--dub -> dub
    if( (jSel.eq.pdfDn_ .and. iSel.eq.pdfAUp_) .or. (jSel.eq.pdfStr_ .and. iSel.eq.pdfAChm_)  .or. (jSel.eq.pdfDn_ .and. iSel.eq.pdfAChm_)  .or. &
         (jSel.eq.pdfStr_ .and. iSel.eq.pdfAUp_)  .or. (jSel.eq.pdfBot_ .and. iSel.eq.pdfAUp_)  .or. (jSel.eq.pdfBot_ .and. iSel.eq.pdfAChm_) ) then
       amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,jSel,iSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfDn_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfDn_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfBot_,pdfAUp_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfBot_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    !--ddb -> uub/ddb
    if( (jSel.eq.pdfBot_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfDn_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfStr_ .and. iSel.eq.pdfABot_) .or.&
         (jSel.eq.pdfBot_ .and. iSel.eq.pdfADn_) .or. (jSel.eq.pdfBot_ .and. iSel.eq.pdfAStr_) ) then
       amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,jSel,iSel)

       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfBot_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfDn_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfStr_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfBot_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfBot_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif




    if( (jSel.eq.pdfDn_ .and. iSel.eq.pdfADn_) .or. (jSel.eq.pdfStr_ .and. iSel.eq.pdfAStr_)  .or. (jSel.eq.pdfDn_ .and. iSel.eq.pdfAStr_)  .or. (jSel.eq.pdfStr_ .and. iSel.eq.pdfADn_) ) then
       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(3,j1,j2,4,za,zb,sprod,jSel,iSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(4,j1,j2,3,za,zb,sprod,jSel,iSel,useWWcoupl=.true.,Wpm_flip=.false.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfDn_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfDn_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfStr_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfChm_) ) then

       !-- non-symmetric qq processes
       amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,jSel,iSel)
       !--uc -> uc
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfUp_,pdfChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       return
    endif


    !--us -> us/cd
    if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfBot_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfBot_) ) then

       amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,jSel,iSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfUp_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfChm_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfStr_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfDn_) ) then

       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,jSel,iSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(3,j2,4,j1,za,zb,sprod,iSel,jSel,useWWcoupl=.true.,Wpm_flip=.false.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfUp_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !--ds -> ds
    if( (jSel.eq.pdfDn_ .and. iSel.eq.pdfStr_) .or. (jSel.eq.pdfDn_ .and. iSel.eq.pdfBot_) .or. (jSel.eq.pdfStr_ .and. iSel.eq.pdfBot_) ) then

       amp_z = A0_ZZ_4f(3,j1,4,j2,za,zb,sprod,jSel,iSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfDn_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfDn_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfStr_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    !-- qbqb processes
    if( (jSel.eq.pdfAUp_ .and. iSel.eq.pdfAChm_) ) then
       amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,jSel,iSel)
       !--ubcb -> ubcb
       restmp = ((abs(amp_z(-1,-1))**2)+ &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfAUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
       return
    endif


    !--ubsb -> ubsb//cbdb
    if( (jSel.eq.pdfAUp_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfAChm_ .and. iSel.eq.pdfABot_) ) then
       amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,jSel,iSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0


       pdfindex = flip(iflip,pdfAUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfAChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif


    if( (jSel.eq.pdfAUp_ .and. iSel.eq.pdfAStr_) .or. (jSel.eq.pdfAChm_ .and. iSel.eq.pdfADn_) ) then
       if( zz_fusion ) then
          amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,jSel,iSel)
          restmp = ((abs(amp_z(-1,-1))**2) + &
               (abs(amp_z(-1,+1))**2) + &
               (abs(amp_z(+1,-1))**2) + &
               (abs(amp_z(+1,+1))**2)) * SpinAvg * tag1
       else
          amp_w = A0_WW_4f(j1,4,j2,3,za,zb,sprod,jSel,iSel,useWWcoupl=.true.,Wpm_flip=.false.)
          restmp= ( &
                  abs(amp_w(-1,-1))**2 + &
                  abs(amp_w(-1,+1))**2 + &
                  abs(amp_w(+1,-1))**2 + &
                  abs(amp_w(+1,+1))**2   &
                  ) * SpinAvg * tag2
       endif

       pdfindex = flip(iflip,pdfAUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfAChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       return
    endif


    !--dbsb -> dbsb
    if( (jSel.eq.pdfADn_ .and. iSel.eq.pdfAStr_) .or. (jSel.eq.pdfADn_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfAStr_ .and. iSel.eq.pdfABot_) ) then

       amp_z = A0_ZZ_4f(j1,3,j2,4,za,zb,sprod,jSel,iSel)
       restmp = ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0

       pdfindex = flip(iflip,pdfADn_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfADn_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfAStr_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       return
    endif

    RETURN
  END SUBROUTINE EvalAmp_WBFH_UnSymm_SA_Select

  SUBROUTINE EvalAmp_WBFH_UnSymm_SA_Select_exact(p,iSel,jSel,rSel,sSel,res)
    implicit none
    real(dp), intent(in) :: p(4,5)
    real(dp), intent(out) :: res
    integer, intent(in) :: iSel,jSel,rSel,sSel
    complex(dp) :: amp_z(-1:1,-1:1), amp_z_b(-1:1,-1:1)
    complex(dp) :: amp_w(-1:1,-1:1)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp
    real(dp) :: ckm_wfactor
    real(dp) :: kfactor_z,kfactor_w
    integer :: i, j, jz1, jz2, jw1, jw2, pdfindex(2)
    integer :: line_i,line_j
    integer :: kz1, kz2, kw1, kw2
    logical :: ZZ_fusion,WW_fusion

    !print *, "Begin EvalAmp_WBFH_UnSymm_SA_Select_exact"
    !print *, "iSel: ",iSel,", jSel: ",jSel," rSel: ",rSel," sSel: ",sSel

    res=0.0_dp

    ZZ_fusion=.false.
    WW_fusion=.false.

    ckm_wfactor=1.0_dp
    kfactor_z=1.0_dp
    kfactor_w=1.0_dp
    restmp=0.0_dp

    line_i=0
    line_j=0
    jz1 = 1
    jz2 = 2
    jw1 = jz1
    jw2 = jz2

    if( (abs(iSel).eq.pdfTop_ .or. abs(jSel).eq.pdfTop_) .or. (abs(rSel).eq.pdfTop_ .or. abs(sSel).eq.pdfTop_) ) return
    if( (abs(iSel).eq.pdfGlu_ .or. abs(jSel).eq.pdfGlu_) .or. (abs(rSel).eq.pdfGlu_ .or. abs(sSel).eq.pdfGlu_) ) return

    if( &
         (iSel.eq.rSel .and. jSel.eq.sSel) &
         .or. &
         (iSel.eq.sSel .and. jSel.eq.rSel) &
         ) then
       ZZ_fusion=.true.
       if( iSel.eq.sSel .and. jSel.eq.rSel ) then
          kz1 = 4
          kz2 = 3
          kfactor_z = 1.0_dp ! dsqrt(ScaleFactor(iSel,sSel)*ScaleFactor(jSel,rSel))
          !print *, "EvalAmp_WBFH_UnSymm_SA_Select_exact: isZZ and is-jr"
       else
          kz1 = 3
          kz2 = 4
          kfactor_z = 1.0_dp ! dsqrt(ScaleFactor(iSel,rSel)*ScaleFactor(jSel,sSel))
          !print *, "EvalAmp_WBFH_UnSymm_SA_Select_exact: isZZ and ir-js"
       endif
    endif

    if( &
         ( (sign(iSel,rSel).eq.iSel .and. sign(jSel,sSel).eq.jSel) .and. (abs(iSel-rSel).eq.1 .or. abs(iSel-rSel).eq.3 .or. abs(iSel-rSel).eq.5) .and. (abs(jSel-sSel).eq.1 .or. abs(jSel-sSel).eq.3 .or. abs(jSel-sSel).eq.5) ) &
         .or. &
         ( (sign(iSel,sSel).eq.iSel .and. sign(jSel,rSel).eq.jSel) .and. (abs(iSel-sSel).eq.1 .or. abs(iSel-sSel).eq.3 .or. abs(iSel-sSel).eq.5) .and. (abs(jSel-rSel).eq.1 .or. abs(jSel-rSel).eq.3 .or. abs(jSel-rSel).eq.5) ) &
         ) then
       WW_fusion=.true.
       ! W_is W_jr fusion
       if( (sign(iSel,sSel).eq.iSel .and. sign(jSel,rSel).eq.jSel) .and. (abs(iSel-sSel).eq.1 .or. abs(iSel-sSel).eq.3 .or. abs(iSel-sSel).eq.5) .and. (abs(jSel-rSel).eq.1 .or. abs(jSel-rSel).eq.3 .or. abs(jSel-rSel).eq.5) ) then
          kw1 = 4
          kw2 = 3
          kfactor_w = 1.0_dp ! dsqrt(ScaleFactor(iSel,sSel)*ScaleFactor(jSel,rSel))

          !if(ZZ_fusion) then ! Special treatment for WW+ZZ interference, not included through phasespace
          ckm_wfactor = CKMbare(iSel,sSel)*CKMbare(jSel,rSel)
          !print *, "iSel, sSel: ",iSel," ",sSel, "; jSel, rSel: ",jSel," ",rSel, ", ckm: ",ckm_wfactor
          !endif
          !print *, "EvalAmp_WBFH_UnSymm_SA_Select_exact: isWW and is-jr"
       else	! W_ir W_js fusion
          kw1 = 3
          kw2 = 4
          kfactor_w = 1.0_dp ! dsqrt(ScaleFactor(iSel,rSel)*ScaleFactor(jSel,sSel))

          ckm_wfactor = CKMbare(iSel,rSel)*CKMbare(jSel,sSel)
          !print *, "iSel, rSel: ",iSel," ",rSel, "; jSel, sSel: ",jSel," ",sSel, ", ckm: ",ckm_wfactor
          !print *, "EvalAmp_WBFH_UnSymm_SA_Select_exact: isWW and ir-js"
       endif
    endif


    if( .not.(ZZ_fusion .or. WW_fusion) ) return
    if(iSel.lt.0) then
       if(ZZ_fusion) call swapi(jz1,kz1)
       if(WW_fusion) call swapi(jw1,kw1)
    endif
    if(jSel.lt.0) then
       if(ZZ_fusion) call swapi(jz2,kz2)
       if(WW_fusion) call swapi(jw2,kw2)
    endif

    line_i = iSel
    line_j = jSel

    if(WW_fusion) then
       if (iSel.eq.pdfUp_ .or. iSel.eq.pdfChm_ .or. iSel.eq.pdfADn_ .or. iSel.eq.pdfAStr_ .or. iSel.eq.pdfABot_) then ! W+ should be passed as the second set of partons
          call swapi(jw1,jw2)
          call swapi(kw1,kw2)
          call swapi(line_i,line_j)
          if(ZZ_fusion) then ! If also ZZ fusion, swap it as well
             call swapi(jz1,jz2)
             call swapi(kz1,kz2)
          endif
       endif
    endif

    call spinoru2(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    !if(ZZ_fusion) print *, "jz1: ",jz1,", jz2: ",jz2," kz1: ",kz1," kz2: ",kz2
    !if(WW_fusion) print *, "jw1: ",jw1,", jw2: ",jw2," kw1: ",kw1," kw2: ",kw2

    !ckm_wfactor = 1.0_dp ! TEST!!!
    if(ZZ_fusion) then

       amp_z = A0_ZZ_4f(kz1,jz1,kz2,jz2,za,zb,sprod,line_i,line_j)

       restmp = restmp + &
            ((abs(amp_z(-1,-1))**2) + &
            (abs(amp_z(-1,+1))**2) + &
            (abs(amp_z(+1,-1))**2) + &
            (abs(amp_z(+1,+1))**2)) * xn**2

       if(iSel.eq.jSel) then
          amp_z_b = -A0_ZZ_4f(kz2,jz1,kz1,jz2,za,zb,sprod,line_i,line_j)
          restmp = restmp + &
               ((abs(amp_z_b(-1,-1))**2) + &
               (abs(amp_z_b(-1,+1))**2) + &
               (abs(amp_z_b(+1,-1))**2) + &
               (abs(amp_z_b(+1,+1))**2)) * xn**2
          restmp = restmp + &
               (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) + &
               two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp)) * xn

          restmp = restmp * SymmFac
       endif
       restmp = restmp * kfactor_z**2
    endif

    if(WW_fusion) then
       amp_w = -A0_WW_4f(kw1,jw1,kw2,jw2,za,zb,sprod,line_i,line_j,useWWcoupl=.true.)
       restmp = restmp + &
                (        &
                abs(amp_w(-1,-1))**2 + &
                abs(amp_w(-1,+1))**2 + &
                abs(amp_w(+1,-1))**2 + &
                abs(amp_w(+1,+1))**2   &
                ) * xn**2 * kfactor_w**2 * abs(ckm_wfactor)**2
       if(ZZ_fusion) then
          restmp = restmp + &
                   two * (  &
                   real(amp_w(-1,-1)*ckm_wfactor*conjg(amp_z(-1,-1)),kind=dp) + &
                   real(amp_w(-1,+1)*ckm_wfactor*conjg(amp_z(-1,+1)),kind=dp) + &
                   real(amp_w(+1,-1)*ckm_wfactor*conjg(amp_z(+1,-1)),kind=dp) + &
                   real(amp_w(+1,+1)*ckm_wfactor*conjg(amp_z(+1,+1)),kind=dp)   &
                   ) * xn * kfactor_z * kfactor_w
       endif
    endif

    restmp = restmp * aveqq
    if(abs(iSel).eq.pdfBot_ .or. abs(jSel).eq.pdfBot_) restmp = restmp * tagbot
    res = restmp

    RETURN
  END SUBROUTINE EvalAmp_WBFH_UnSymm_SA_Select_exact




!-- TEST FUNCTION THAT WRAPS THE HVV AMPLITUDE
!-- FOR IT TO WORK, COMMENT OUT THE LINES BELOW,
!-- AND THE QUOTED LINE IN MODHIGGS::CALCHELAMP2:
!-- "if ((q_q).lt.-0.1d0 .or. (q3_q3).lt.-0.1d0 .or. (q4_q4).lt.-0.1d0) return  ! if negative invariant masses return zero"

!  SUBROUTINE wrapHVV(p,iSel,jSel,rSel,sSel,res)
!    use ModHiggs
!    implicit none
!    real(dp), intent(in) :: p(4,5)
!    real(dp), intent(out) :: res
!    integer, intent(in) :: iSel,jSel,rSel,sSel
!    integer :: i, j, jz1, jz2, jw1, jw2, pdfindex(2)
!    integer :: kz1, kz2, kw1, kw2
!    logical :: ZZ_fusion,WW_fusion
!    complex(dp) :: A_VV(1:8)
!    integer, parameter :: ZZMode=00,ZgsMode=01,gsZMode=02,gsgsMode=03
!    integer, parameter :: WWMode=10
!    integer, parameter :: ggMode=20
!    integer, parameter :: ZgMode=30,gsgMode=31
!    integer :: MY_IDUP(3:6),i3,i4
!    real(dp) :: pUsed(4,6)

!    if( (abs(iSel).eq.pdfTop_ .or. abs(jSel).eq.pdfTop_) .or. (abs(rSel).eq.pdfTop_ .or. abs(sSel).eq.pdfTop_) ) return

!    res=0.0_dp

!    ZZ_fusion=.false.
!    WW_fusion=.false.

!    jz1 = 1
!    jz2 = 2
!    jw1 = jz1
!    jw2 = jz2

!    if( &
!         (iSel.eq.rSel .and. jSel.eq.sSel) &
!         .or. &
!         (iSel.eq.sSel .and. jSel.eq.rSel) &
!         ) then
!       ZZ_fusion=.true.
!       if( iSel.eq.sSel .and. jSel.eq.rSel ) then
!          kz1 = 4
!          kz2 = 3
!       else
!          kz1 = 3
!          kz2 = 4
!       endif
!    endif

!    if( &
!         ( (sign(iSel,rSel).eq.iSel .and. sign(jSel,sSel).eq.jSel) .and. (abs(iSel-rSel).eq.1 .or. abs(iSel-rSel).eq.3 .or. abs(iSel-rSel).eq.5) .and. (abs(jSel-sSel).eq.1 .or. abs(jSel-sSel).eq.3 .or. abs(jSel-sSel).eq.5) ) &
!         .or. &
!         ( (sign(iSel,sSel).eq.iSel .and. sign(jSel,rSel).eq.jSel) .and. (abs(iSel-sSel).eq.1 .or. abs(iSel-sSel).eq.3 .or. abs(iSel-sSel).eq.5) .and. (abs(jSel-rSel).eq.1 .or. abs(jSel-rSel).eq.3 .or. abs(jSel-rSel).eq.5) ) &
!         ) then
!       WW_fusion=.true.
!       ! W_is W_jr fusion
!       if( (sign(iSel,sSel).eq.iSel .and. sign(jSel,rSel).eq.jSel) .and. (abs(iSel-sSel).eq.1 .or. abs(iSel-sSel).eq.3 .or. abs(iSel-sSel).eq.5) .and. (abs(jSel-rSel).eq.1 .or. abs(jSel-rSel).eq.3 .or. abs(jSel-rSel).eq.5) ) then
!          kw1 = 4
!          kw2 = 3
!       else	! W_ir W_js fusion
!          kw1 = 3
!          kw2 = 4
!       endif
!    endif


!    if( .not.(ZZ_fusion .or. WW_fusion) ) return
!    if(iSel.lt.0) then
!       if(ZZ_fusion) call swapi(jz1,kz1)
!       if(WW_fusion) call swapi(jw1,kw1)
!    endif
!    if(jSel.lt.0) then
!       if(ZZ_fusion) call swapi(jz2,kz2)
!       if(WW_fusion) call swapi(jw2,kw2)
!    endif

!    if(WW_fusion) then
!       if (iSel.eq.pdfUp_ .or. iSel.eq.pdfChm_ .or. iSel.eq.pdfADn_ .or. iSel.eq.pdfAStr_ .or. iSel.eq.pdfABot_) then ! W+ should be passed as the second set of partons
!          call swapi(jw1,jw2)
!          call swapi(kw1,kw2)
!          if(ZZ_fusion) then ! If also ZZ fusion, swap it as well
!             call swapi(jz1,jz2)
!             call swapi(kz1,kz2)
!          endif
!       endif
!    endif

!    if(jz1.eq.1) then
!      MY_IDUP(4) = -convertFromPartIndex(iSel)
!      pUsed(:,4) = -p(:,1)
!    else if(jz1.eq.2) then
!      MY_IDUP(4) = -convertFromPartIndex(jSel)
!      pUsed(:,4) = -p(:,2)
!    else if(jz1.eq.3) then
!      MY_IDUP(4) = convertFromPartIndex(rSel)
!      pUsed(:,4) = p(:,3)
!    else if(jz1.eq.4) then
!      MY_IDUP(4) = convertFromPartIndex(sSel)
!      pUsed(:,4) = p(:,4)
!    endif
!    if(kz1.eq.1) then
!      MY_IDUP(3) = -convertFromPartIndex(iSel)
!      pUsed(:,3) = -p(:,1)
!    else if(kz1.eq.2) then
!      MY_IDUP(3) = -convertFromPartIndex(jSel)
!      pUsed(:,3) = -p(:,2)
!    else if(kz1.eq.3) then
!      MY_IDUP(3) = convertFromPartIndex(rSel)
!      pUsed(:,3) = p(:,3)
!    else if(kz1.eq.4) then
!      MY_IDUP(3) = convertFromPartIndex(sSel)
!      pUsed(:,3) = p(:,4)
!    endif
!    if(jz2.eq.1) then
!      MY_IDUP(6) = -convertFromPartIndex(iSel)
!      pUsed(:,6) = -p(:,1)
!    else if(jz2.eq.2) then
!      MY_IDUP(6) = -convertFromPartIndex(jSel)
!      pUsed(:,6) = -p(:,2)
!    else if(jz2.eq.3) then
!      MY_IDUP(6) = convertFromPartIndex(rSel)
!      pUsed(:,6) = p(:,3)
!    else if(jz2.eq.4) then
!      MY_IDUP(6) = convertFromPartIndex(sSel)
!      pUsed(:,6) = p(:,4)
!    endif
!    if(kz2.eq.1) then
!      MY_IDUP(5) = -convertFromPartIndex(iSel)
!      pUsed(:,5) = -p(:,1)
!    else if(kz2.eq.2) then
!      MY_IDUP(5) = -convertFromPartIndex(jSel)
!      pUsed(:,5) = -p(:,2)
!    else if(kz2.eq.3) then
!      MY_IDUP(5) = convertFromPartIndex(rSel)
!      pUsed(:,5) = p(:,3)
!    else if(kz2.eq.4) then
!      MY_IDUP(5) = convertFromPartIndex(sSel)
!      pUsed(:,5) = p(:,4)
!    endif

!    pUSed(:,1) = pUSed(:,3) + pUSed(:,4) + pUSed(:,5) + pUSed(:,6)
!    pUSed(:,2) = zero

!    if(ZZ_fusion) then
!       print *, "iSel: ",iSel,", jSel: ",jSel," rSel: ",rSel," sSel: ",sSel
!       print *, "jz1: ",jz1,", jz2: ",jz2," kz1: ",kz1," kz2: ",kz2
!       print *,"MY_IDUP:",MY_IDUP
!       print *,"p1:",p(:,1)
!       print *,"p2:",p(:,2)
!       print *,"p3:",p(:,3)
!       print *,"p4:",p(:,4)
!       print *,"p5:",p(:,5)
!       print *,"pUsed1:",pUsed(:,1)
!       print *,"pUsed2:",pUsed(:,2)
!       print *,"pUsed3:",pUsed(:,3)
!       print *,"pUsed4:",pUsed(:,4)
!       print *,"pUsed5:",pUsed(:,5)
!       print *,"pUsed6:",pUsed(:,6)

!       A_VV(:) = 0d0
!       do i3=1,2;  do i4=1,2!  sum over helicities
!          call calcHelAmp2((/3,4,5,6/),ZZMode,MY_IDUP,pUsed,i3,i4,A_VV(1))
!          if( includeGammaStar ) then
!             call calcHelAmp2((/3,4,5,6/),ZgsMode,MY_IDUP,pUsed,i3,i4,A_VV(3))
!             call calcHelAmp2((/3,4,5,6/),gsZMode,MY_IDUP,pUsed,i3,i4,A_VV(5))
!             call calcHelAmp2((/3,4,5,6/),gsgsMode,MY_IDUP,pUsed,i3,i4,A_VV(7))
!          endif

!          if( (MY_IDUP(3).eq.MY_IDUP(5)) .and. (MY_IDUP(4).eq.MY_IDUP(6)) .and. iSel.eq.jSel ) then ! iSel.eq.jSel to avoid ZH diagram
!             call calcHelAmp2((/5,4,3,6/),ZZMode,MY_IDUP,pUsed,i3,i4,A_VV(2))
!             if( includeGammaStar ) then
!                call calcHelAmp2((/5,4,3,6/),ZgsMode,MY_IDUP,pUsed,i3,i4,A_VV(4))
!                call calcHelAmp2((/5,4,3,6/),gsZMode,MY_IDUP,pUsed,i3,i4,A_VV(6))
!                call calcHelAmp2((/5,4,3,6/),gsgsMode,MY_IDUP,pUsed,i3,i4,A_VV(8))
!             endif
!             A_VV(2) = -A_VV(2)! minus from Fermi statistics
!             A_VV(4) = -A_VV(4)
!             A_VV(6) = -A_VV(6)
!             A_VV(8) = -A_VV(8)
!          endif

!          res = res + (A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7))*dconjg(A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7))!   interfere the 3456 pieces
!          res = res + (A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8))*dconjg(A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8))!   interfere the 5436 pieces
!          if( (MY_IDUP(3).eq.MY_IDUP(5)) .and. (MY_IDUP(4).eq.MY_IDUP(6)) .and. (i3.eq.i4) .and. iSel.eq.jSel ) then! interfere the 3456 with 5436 pieces
!             res = res + 2d0/3d0*dreal(  A_VV(1)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
!             res = res + 2d0/3d0*dreal(  A_VV(3)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
!             res = res + 2d0/3d0*dreal(  A_VV(5)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
!             res = res + 2d0/3d0*dreal(  A_VV(7)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
!          endif

!          print *,"i3 i4 res: ",i3,i4,res
!       enddo;  enddo
!       if(  (MY_IDUP(3).eq.MY_IDUP(5)) .and. (MY_IDUP(4).eq.MY_IDUP(6)) .and. iSel.eq.jSel ) res = res/2d0

!    endif
!    res = res * aveqq

!    RETURN
!  END SUBROUTINE






  function flip(i,a1,a2)
    integer :: flip(2)
    integer :: i, a1, a2

    if (i .eq. 1) flip = (/a1,a2/)
    if (i .eq. 2) flip = (/a2,a1/)

    return

  end function flip


  !-------------------------------------------------------------------------
  !-- amplitudes below
  !-- Older WBF code uses A0_VV_4f, newer subroutines use A0_ZZ_4f or A0_WW_4f -- U.S.: Fully checked with decay amplitudes through a wrapper
  function A0_VV_4f(j1,j2,j3,j4,za,zb,sprod,mv,ga_v,useWWcoupl,Wpm_flip)
  use modMisc
  implicit none
    complex(dp) :: A0_VV_4f(-1:1,-1:1)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: mv, ga_v
    logical,optional :: useWWcoupl,Wpm_flip
    real(dp) :: sprod(4,4),q2Wplus,q2Wminus
    real(dp) :: mhsq, q1q2, kcoupl
    complex(dp) :: a1, a2, a3, struc1, struc2, struc3
    complex(dp) :: zab2
    complex(dp) :: iprop12, iprop34
    complex(dp) :: vvcoupl_prime(4)
    integer :: vv_it
    integer :: i,j,k,l

    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4)
    !The previous line works, and assigns the whole zab2 correctly.
    !I have no idea how it works.
    !If you don't believe me, please uncomment the following lines:
    !do i=1,4
    ! do j=1,4
    !  do k=1,4
    !   do l=1,4
    !    print *,i,j,k,l,zab2(i,j,k,l),za(i,j)*zb(j,l) + za(i,k)*zb(k,l)
    !   enddo
    !  enddo
    ! enddo
    !enddo
    !print *,j1,j2,j3,j4
    !pause

    A0_VV_4f = czero

    q1q2 = (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))/two
    mhsq = two * q1q2 + sprod(j1,j2) + sprod(j3,j4)

    kcoupl = q1q2/lambda**2

    q2Wplus  = sprod(j1,j2)
    q2Wminus = sprod(j3,j4)
    if( present(Wpm_flip) ) then
       if( Wpm_flip ) call swapr(q2Wplus,q2Wminus)
    endif
    if( .not.present(useWWcoupl) ) then
       do vv_it=1,4
          vvcoupl_prime(vv_it) = HVVSpinZeroDynamicCoupling(vv_it,q2Wplus,q2Wminus,mhsq)
       enddo
    else
       do vv_it=1,4
          vvcoupl_prime(vv_it) = HVVSpinZeroDynamicCoupling(vv_it,q2Wplus,q2Wminus,mhsq,tryWWcoupl=useWWcoupl)
       enddo
    endif


    a1 = vvcoupl_prime(1) * mv**2/mhsq + vvcoupl_prime(2) * two * q1q2/mhsq + vvcoupl_prime(3) * kcoupl * q1q2/mhsq
    a2 = -two * vvcoupl_prime(2) - kcoupl * vvcoupl_prime(3)
    a3 = -two * vvcoupl_prime(4)

    struc1 = two * (a1 * mhsq - ci * a3 * q1q2)
    struc2 = a2 + ci * a3
    struc3 = two * ci * a3


    A0_VV_4f(-1,-1) = za(j1,j3)*zb(j4,j2) * struc1 + &
         zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4) * struc2 + &
         za(j1,j2)*za(j3,j4)*zb(j4,j2)**2 * struc3

    A0_VV_4f(-1,+1) = za(j1,j4)*zb(j3,j2) * struc1 + &
         zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3) * struc2 + &
         za(j1,j2)*za(j4,j3)*zb(j3,j2)**2 * struc3

    A0_VV_4f(+1,-1) = za(j2,j3)*zb(j4,j1) * struc1 + &
         zab2(j2,j3,j4,j1)*zab2(j3,j1,j2,j4) * struc2 + &
         za(j2,j1)*za(j3,j4)*zb(j4,j1)**2 * struc3

    A0_VV_4f(+1,+1) = za(j2,j4)*zb(j3,j1) * struc1 + &
         zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3) * struc2 + &
         za(j2,j1)*za(j4,j3)*zb(j3,j1)**2 * struc3

    iprop12 = sprod(j1,j2) - mv**2 + ci * mv * ga_v
    iprop34 = sprod(j3,j4) - mv**2 + ci * mv * ga_v

    A0_VV_4f = A0_VV_4f/vev /iprop12/iprop34

    return

  end function A0_VV_4f

  !--           line = 1/2 --> down/up couplings included
  function A0_ZZ_4f(j1,j2,j3,j4,za,zb,sprod,iSel,jSel)
    use modMisc
    implicit none
    real(dp), dimension(5) :: Lz,Rz,LCT,RCT
    real(dp), parameter, dimension(5) :: La = (/QdL,QuL,QdL,QuL,QdL/)
    real(dp), parameter, dimension(5) :: Ra = (/QdR,QuR,QdR,QuR,QdR/)
    complex(dp) :: A0_ZZ_4f(-1:1,-1:1)
    integer :: j1,j2,j3,j4,iSel,jSel,line1,line2
    complex(dp) :: za(4,4),zb(4,4)
    real(dp) :: sprod(4,4)
    real(dp) :: mhsq,q1q2,kcoupl,q12sq,q34sq
    complex(dp) :: a1_zz,a2_zz,a3_zz,struc_zz(3)
    complex(dp) :: a1_aa,a2_aa,a3_aa,a1_az,a2_az,a3_az,a1_za,a2_za,a3_za
    complex(dp) :: a1_zzp,a2_zzp,a3_zzp,struc_zzp(3)
    complex(dp) :: a1_zpz,a2_zpz,a3_zpz,struc_zpz(3)
    complex(dp) :: a1_zpzp,a2_zpzp,a3_zpzp,struc_zpzp(3)
    complex(dp) :: struc_aa(3),struc_az(3),struc_za(3)
    complex(dp) :: helcoup(1:3,-1:1,-1:1)
    complex(dp) :: zab2
    complex(dp) :: iprop12,iprop34, zpprop12,zpprop34
    complex(dp) :: vvcoupl_prime_zz(4),vvcoupl_prime_az(4),vvcoupl_prime_za(4),vvcoupl_prime_aa(2:4)
    complex(dp) :: vvcoupl_prime_zzp(4),vvcoupl_prime_zpz(4),vvcoupl_prime_zpzp(4)
    integer :: vv_it
    integer :: i,j,k,l

    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4)

    A0_ZZ_4f = czero
    helcoup(:,:,:)=czero

    Lz = (/aL_Qdn,aL_Qup,aL_Qdn,aL_Qup,aL_Qdn/)
    Rz = (/aR_Qdn,aR_Qup,aR_Qdn,aR_Qup,aR_Qdn/)

    if( abs(iSel).eq.pdfTop_ .or. abs(jSel).eq.pdfTop_ ) return
    if( abs(iSel).eq.pdfGlu_ .or. abs(jSel).eq.pdfGlu_ ) return

    line1=abs(iSel)
    line2=abs(jSel)

    q1q2 = (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))/two
    mhsq = two * q1q2 + sprod(j1,j2) + sprod(j3,j4)

    kcoupl = q1q2/lambda**2

    q12sq = sprod(j1,j2)
    q34sq = sprod(j3,j4)

    iprop12 = q12sq - M_Z**2 + ci * M_Z * Ga_Z
    iprop34 = q34sq - M_Z**2 + ci * M_Z * Ga_Z

    !-- set up couplings
    do vv_it=1,4
       vvcoupl_prime_zz(vv_it) = HVVSpinZeroDynamicCoupling(vv_it,q12sq,q34sq,mhsq)
    enddo

    a1_zz = vvcoupl_prime_zz(1) * M_Z**2/mhsq + vvcoupl_prime_zz(2) * two * q1q2/mhsq + vvcoupl_prime_zz(3) * kcoupl * q1q2/mhsq
    a2_zz = -two * vvcoupl_prime_zz(2) - kcoupl * vvcoupl_prime_zz(3)
    a3_zz = -two * vvcoupl_prime_zz(4)

    struc_zz(1) = two * (a1_zz * mhsq - ci * a3_zz * q1q2)
    struc_zz(2) = (a2_zz + ci * a3_zz)
    struc_zz(3) = two * ci * a3_zz
    struc_zz(:) = struc_zz(:) * couplZffsq

    if( includeGammaStar ) then

       do vv_it=1,4
          vvcoupl_prime_az(vv_it) = HVVSpinZeroDynamicCoupling(vv_it+4,0d0,q12sq,mhsq)
          vvcoupl_prime_za(vv_it) = HVVSpinZeroDynamicCoupling(vv_it+4,0d0,q34sq,mhsq)
       enddo
       do vv_it=1,3
          vvcoupl_prime_aa(vv_it+1) = HVVSpinZeroDynamicCoupling(vv_it+8,q12sq,q34sq,mhsq)
       enddo

       a1_aa = vvcoupl_prime_aa(2) * two * q1q2/mhsq + vvcoupl_prime_aa(3) * kcoupl * q1q2/mhsq
       a2_aa = -two * vvcoupl_prime_aa(2) - kcoupl * vvcoupl_prime_aa(3)
       a3_aa = -two * vvcoupl_prime_aa(4)

       a1_za = vvcoupl_prime_za(1) * M_Z**2/mhsq + vvcoupl_prime_za(2) * two * q1q2/mhsq + vvcoupl_prime_za(3) * kcoupl * q1q2/mhsq
       a2_za = -two * vvcoupl_prime_za(2) - kcoupl * vvcoupl_prime_za(3)
       a3_za = -two * vvcoupl_prime_za(4)

       a1_az = vvcoupl_prime_az(1) * M_Z**2/mhsq + vvcoupl_prime_az(2) * two * q1q2/mhsq + vvcoupl_prime_az(3) * kcoupl * q1q2/mhsq
       a2_az = -two * vvcoupl_prime_az(2) - kcoupl * vvcoupl_prime_az(3)
       a3_az = -two * vvcoupl_prime_az(4)

       struc_aa(1) = two * (a1_aa * mhsq - ci * a3_aa * q1q2)
       struc_aa(2) = (a2_aa + ci * a3_aa)
       struc_aa(3) = two * ci * a3_aa
       struc_aa(:) = struc_aa(:) * couplAffsq

       struc_az(1) = two * (a1_az * mhsq - ci * a3_az * q1q2)
       struc_az(2) = (a2_az + ci * a3_az)
       struc_az(3) = two * ci * a3_az
       struc_az(:) = struc_az(:) * couplAZff

       struc_za(1) = two * (a1_za * mhsq - ci * a3_za * q1q2)
       struc_za(2) = (a2_za + ci * a3_za)
       struc_za(3) = two * ci * a3_za
       struc_za(:) = struc_za(:) * couplAZff

       helcoup(1:3,-1,-1) = struc_zz(1:3) * Lz(line1) * Lz(line2)/iprop12/iprop34 + &
                            struc_aa(1:3) * La(line1) * La(line2)/q12sq  /q34sq   + &
                            struc_az(1:3) * La(line1) * Lz(line2)/q12sq  /iprop34 + &
                            struc_za(1:3) * Lz(line1) * La(line2)/iprop12/q34sq

       helcoup(1:3,-1,+1) = struc_zz(1:3) * Lz(line1) * Rz(line2)/iprop12/iprop34 + &
                            struc_aa(1:3) * La(line1) * Ra(line2)/q12sq  /q34sq   + &
                            struc_az(1:3) * La(line1) * Rz(line2)/q12sq  /iprop34 + &
                            struc_za(1:3) * Lz(line1) * Ra(line2)/iprop12/q34sq

       helcoup(1:3,+1,-1) = struc_zz(1:3) * Rz(line1) * Lz(line2)/iprop12/iprop34 + &
                            struc_aa(1:3) * Ra(line1) * La(line2)/q12sq  /q34sq   + &
                            struc_az(1:3) * Ra(line1) * Lz(line2)/q12sq  /iprop34 + &
                            struc_za(1:3) * Rz(line1) * La(line2)/iprop12/q34sq

       helcoup(1:3,+1,+1) = struc_zz(1:3) * Rz(line1) * Rz(line2)/iprop12/iprop34 + &
                            struc_aa(1:3) * Ra(line1) * Ra(line2)/q12sq  /q34sq   + &
                            struc_az(1:3) * Ra(line1) * Rz(line2)/q12sq  /iprop34 + &
                            struc_za(1:3) * Rz(line1) * Ra(line2)/iprop12/q34sq

    else

       helcoup(1:3,-1,-1) = struc_zz(1:3) * Lz(line1) * Lz(line2)/iprop12/iprop34
       helcoup(1:3,-1,+1) = struc_zz(1:3) * Lz(line1) * Rz(line2)/iprop12/iprop34
       helcoup(1:3,+1,-1) = struc_zz(1:3) * Rz(line1) * Lz(line2)/iprop12/iprop34
       helcoup(1:3,+1,+1) = struc_zz(1:3) * Rz(line1) * Rz(line2)/iprop12/iprop34

    endif

    if ( includeVprime ) then

       LCT = (/ezp_Dn_left,ezp_Up_left,ezp_Str_left,ezp_Chm_left,ezp_Bot_left/)
       RCT = (/ezp_Dn_right,ezp_Up_right,ezp_Str_right,ezp_Chm_right,ezp_Bot_right/)

       do vv_it=1,4
          vvcoupl_prime_zzp(vv_it) = HVVSpinZeroDynamicCoupling(11+vv_it,q12sq,q34sq,mhsq)
          vvcoupl_prime_zpz(vv_it) = HVVSpinZeroDynamicCoupling(11+vv_it,q34sq,q12sq,mhsq)
          vvcoupl_prime_zpzp(vv_it) = HVVSpinZeroDynamicCoupling(15+vv_it,q12sq,q34sq,mhsq)
       enddo

       if( UseVprime ) then
          if (M_Zprime.gt.0d0) then
             zpprop12 = q12sq - M_Zprime**2 + ci * M_Zprime * Ga_Zprime
             zpprop34 = q34sq - M_Zprime**2 + ci * M_Zprime * Ga_Zprime
          else
             zpprop12 = q12sq
             zpprop34 = q34sq
          endif
       else
          zpprop12 = M_Z**2
          zpprop34 = M_Z**2
       endif


       a1_zzp = vvcoupl_prime_zzp(1) * M_Z**2/mhsq + vvcoupl_prime_zzp(2) * two * q1q2/mhsq + vvcoupl_prime_zzp(3) * kcoupl * q1q2/mhsq
       a2_zzp = -two * vvcoupl_prime_zzp(2) - kcoupl * vvcoupl_prime_zzp(3)
       a3_zzp = -two * vvcoupl_prime_zzp(4)
       a1_zpz = vvcoupl_prime_zpz(1) * M_Z**2/mhsq + vvcoupl_prime_zpz(2) * two * q1q2/mhsq + vvcoupl_prime_zpz(3) * kcoupl * q1q2/mhsq
       a2_zpz = -two * vvcoupl_prime_zpz(2) - kcoupl * vvcoupl_prime_zpz(3)
       a3_zpz = -two * vvcoupl_prime_zpz(4)
       a1_zpzp = vvcoupl_prime_zpzp(1) * M_Z**2/mhsq + vvcoupl_prime_zpzp(2) * two * q1q2/mhsq + vvcoupl_prime_zpzp(3) * kcoupl * q1q2/mhsq
       a2_zpzp = -two * vvcoupl_prime_zpzp(2) - kcoupl * vvcoupl_prime_zpzp(3)
       a3_zpzp = -two * vvcoupl_prime_zpzp(4)

       struc_zzp(1) = two * (a1_zzp * mhsq - ci * a3_zzp * q1q2)
       struc_zzp(2) = (a2_zzp + ci * a3_zzp)
       struc_zzp(3) = two * ci * a3_zzp
       struc_zpz(1) = two * (a1_zpz * mhsq - ci * a3_zpz * q1q2)
       struc_zpz(2) = (a2_zpz + ci * a3_zpz)
       struc_zpz(3) = two * ci * a3_zpz
       struc_zpzp(1) = two * (a1_zpzp * mhsq - ci * a3_zpzp * q1q2)
       struc_zpzp(2) = (a2_zpzp + ci * a3_zpzp)
       struc_zpzp(3) = two * ci * a3_zpzp

       struc_zzp(:) = struc_zzp(:) * couplZffsq
       struc_zpz(:) = struc_zpz(:) * couplZffsq
       struc_zpzp(:) = struc_zpzp(:) * couplZffsq

       helcoup(1:3,-1,-1) = helcoup(1:3,-1,-1) + &
                            struc_zpzp(1:3) * LCT(line1) * LCT(line2)/zpprop12  /zpprop34   + &
                            struc_zpz(1:3) * LCT(line1) * Lz(line2)/zpprop12  /iprop34 + &
                            struc_zzp(1:3) * Lz(line1) * LCT(line2)/iprop12/zpprop34

       helcoup(1:3,-1,+1) = helcoup(1:3,-1,+1) + &
                            struc_zpzp(1:3) * LCT(line1) * RCT(line2)/zpprop12  /zpprop34   + &
                            struc_zpz(1:3) * LCT(line1) * Rz(line2)/zpprop12  /iprop34 + &
                            struc_zzp(1:3) * Lz(line1) * RCT(line2)/iprop12/zpprop34

       helcoup(1:3,+1,-1) = helcoup(1:3,+1,-1) + &
                            struc_zpzp(1:3) * RCT(line1) * LCT(line2)/zpprop12  /zpprop34   + &
                            struc_zpz(1:3) * RCT(line1) * Lz(line2)/zpprop12  /iprop34 + &
                            struc_zzp(1:3) * Rz(line1) * LCT(line2)/iprop12/zpprop34

       helcoup(1:3,+1,+1) = helcoup(1:3,+1,+1) + &
                            struc_zpzp(1:3) * RCT(line1) * RCT(line2)/zpprop12  /zpprop34   + &
                            struc_zpz(1:3) * RCT(line1) * Rz(line2)/zpprop12  /iprop34 + &
                            struc_zzp(1:3) * Rz(line1) * RCT(line2)/iprop12/zpprop34

    endif



    A0_ZZ_4f(-1,-1) = za(j1,j3)*zb(j4,j2) * helcoup(1,-1,-1) + &
         zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4) * helcoup(2,-1,-1) + &
         za(j1,j2)*za(j3,j4)*zb(j4,j2)**2 * helcoup(3,-1,-1)

    A0_ZZ_4f(-1,+1) = za(j1,j4)*zb(j3,j2) * helcoup(1,-1,+1) + &
         zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3) * helcoup(2,-1,+1) + &
         za(j1,j2)*za(j4,j3)*zb(j3,j2)**2 * helcoup(3,-1,+1)

    A0_ZZ_4f(+1,-1) = za(j2,j3)*zb(j4,j1) * helcoup(1,+1,-1) + &
         zab2(j2,j3,j4,j1)*zab2(j3,j1,j2,j4) * helcoup(2,+1,-1) + &
         za(j2,j1)*za(j3,j4)*zb(j4,j1)**2 * helcoup(3,+1,-1)

    A0_ZZ_4f(+1,+1) = za(j2,j4)*zb(j3,j1) * helcoup(1,+1,+1) + &
         zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3) * helcoup(2,+1,+1) + &
         za(j2,j1)*za(j4,j3)*zb(j3,j1)**2 * helcoup(3,+1,+1)

    A0_ZZ_4f = A0_ZZ_4f/vev

    return

  end function A0_ZZ_4f

  function A0_WW_4f(j1,j2,j3,j4,za,zb,sprod,iSel,jSel,useWWcoupl,Wpm_flip)
  use modMisc
  implicit none
    complex(dp) :: A0_WW_4f(-1:1,-1:1)
    integer :: j1,j2,j3,j4,iSel,jSel,line1,line2
    complex(dp) :: za(4,4), zb(4,4)
    logical,optional :: useWWcoupl,Wpm_flip
    real(dp) :: sprod(4,4),q2Wplus,q2Wminus
    real(dp) :: mhsq, q1q2, kcoupl
    complex(dp) :: vvcoupl_prime_ww(4)
    complex(dp) :: a1_ww, a2_ww, a3_ww, struc_ww(1:3)
    complex(dp) :: vvcoupl_prime_wwp(4)
    complex(dp) :: a1_wwp, a2_wwp, a3_wwp, struc_wwp(1:3)
    complex(dp) :: vvcoupl_prime_wpw(4)
    complex(dp) :: a1_wpw, a2_wpw, a3_wpw, struc_wpw(1:3)
    complex(dp) :: vvcoupl_prime_wpwp(4)
    complex(dp) :: a1_wpwp, a2_wpwp, a3_wpwp, struc_wpwp(1:3)
    complex(dp) :: zab2
    complex(dp) :: iprop12, iprop34, wpprop12, wpprop34
    complex(dp) :: helcoup(1:3,-1:1,-1:1)
    real(dp), dimension(5) :: LCT,RCT
    integer :: vv_it
    integer :: i,j,k,l

    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4)

    A0_WW_4f = czero
    helcoup(:,:,:)=czero

    if( abs(iSel).eq.pdfTop_ .or. abs(jSel).eq.pdfTop_ ) return
    if( abs(iSel).eq.pdfGlu_ .or. abs(jSel).eq.pdfGlu_ ) return

    line1=abs(iSel)
    line2=abs(jSel)

    q1q2 = (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))/two
    mhsq = two * q1q2 + sprod(j1,j2) + sprod(j3,j4)

    kcoupl = q1q2/lambda**2

    q2Wplus  = sprod(j1,j2)
    q2Wminus = sprod(j3,j4)
    if( present(Wpm_flip) ) then
       if( Wpm_flip ) then
          call swap(q2Wplus,q2Wminus)
          call swap(line1,line2)
       endif
    endif
    if( .not.present(useWWcoupl) ) then
       do vv_it=1,4
          vvcoupl_prime_ww(vv_it) = HVVSpinZeroDynamicCoupling(vv_it,q2Wplus,q2Wminus,mhsq)
       enddo
    else
       do vv_it=1,4
          vvcoupl_prime_ww(vv_it) = HVVSpinZeroDynamicCoupling(vv_it,q2Wplus,q2Wminus,mhsq,tryWWcoupl=useWWcoupl)
       enddo
    endif

    ! SM contribution
    iprop12 = sprod(j1,j2) - M_W**2 + ci * M_W * Ga_W
    iprop34 = sprod(j3,j4) - M_W**2 + ci * M_W * Ga_W

    a1_ww = vvcoupl_prime_ww(1) * M_W**2/mhsq + vvcoupl_prime_ww(2) * two * q1q2/mhsq + vvcoupl_prime_ww(3) * kcoupl * q1q2/mhsq
    a2_ww = -two * vvcoupl_prime_ww(2) - kcoupl * vvcoupl_prime_ww(3)
    a3_ww = -two * vvcoupl_prime_ww(4)

    struc_ww(1) = two * (a1_ww * mhsq - ci * a3_ww * q1q2)
    struc_ww(2) = a2_ww + ci * a3_ww
    struc_ww(3) = two * ci * a3_ww
    struc_ww(:) = struc_ww(:) * couplWffsq

    helcoup(1:3,-1,-1) = struc_ww(1:3)/iprop12/iprop34

    if ( includeVprime ) then

       LCT = (/ewp_Up_left,ewp_Up_left,ewp_Chm_left,ewp_Chm_left,ewp_Top_left/)
       RCT = (/ewp_Up_right,ewp_Up_right,ewp_Chm_right,ewp_Chm_right,ewp_Top_right/)

       do vv_it=1,4
          vvcoupl_prime_wwp(vv_it) = HVVSpinZeroDynamicCoupling(11+vv_it,q2Wplus,q2Wminus,mhsq,tryWWcoupl=useWWcoupl)
          vvcoupl_prime_wpw(vv_it) = HVVSpinZeroDynamicCoupling(11+vv_it,q2Wminus,q2Wplus,mhsq,tryWWcoupl=useWWcoupl)
          vvcoupl_prime_wpwp(vv_it) = HVVSpinZeroDynamicCoupling(15+vv_it,q2Wplus,q2Wminus,mhsq,tryWWcoupl=useWWcoupl)
       enddo

       if( UseVprime ) then
          if (M_Wprime.gt.0d0) then
             wpprop12 = q2Wplus - M_Wprime**2 + ci * M_Wprime * Ga_Wprime
             wpprop34 = q2Wminus - M_Wprime**2 + ci * M_Wprime * Ga_Wprime
          else
             wpprop12 = q2Wplus
             wpprop34 = q2Wminus
          endif
       else
          wpprop12 = M_W**2
          wpprop34 = M_W**2
       endif


       a1_wwp = vvcoupl_prime_wwp(1) * M_Z**2/mhsq + vvcoupl_prime_wwp(2) * two * q1q2/mhsq + vvcoupl_prime_wwp(3) * kcoupl * q1q2/mhsq
       a2_wwp = -two * vvcoupl_prime_wwp(2) - kcoupl * vvcoupl_prime_wwp(3)
       a3_wwp = -two * vvcoupl_prime_wwp(4)
       a1_wpw = vvcoupl_prime_wpw(1) * M_Z**2/mhsq + vvcoupl_prime_wpw(2) * two * q1q2/mhsq + vvcoupl_prime_wpw(3) * kcoupl * q1q2/mhsq
       a2_wpw = -two * vvcoupl_prime_wpw(2) - kcoupl * vvcoupl_prime_wpw(3)
       a3_wpw = -two * vvcoupl_prime_wpw(4)
       a1_wpwp = vvcoupl_prime_wpwp(1) * M_Z**2/mhsq + vvcoupl_prime_wpwp(2) * two * q1q2/mhsq + vvcoupl_prime_wpwp(3) * kcoupl * q1q2/mhsq
       a2_wpwp = -two * vvcoupl_prime_wpwp(2) - kcoupl * vvcoupl_prime_wpwp(3)
       a3_wpwp = -two * vvcoupl_prime_wpwp(4)

       struc_wwp(1) = two * (a1_wwp * mhsq - ci * a3_wwp * q1q2)
       struc_wwp(2) = (a2_wwp + ci * a3_wwp)
       struc_wwp(3) = two * ci * a3_wwp
       struc_wpw(1) = two * (a1_wpw * mhsq - ci * a3_wpw * q1q2)
       struc_wpw(2) = (a2_wpw + ci * a3_wpw)
       struc_wpw(3) = two * ci * a3_wpw
       struc_wpwp(1) = two * (a1_wpwp * mhsq - ci * a3_wpwp * q1q2)
       struc_wpwp(2) = (a2_wpwp + ci * a3_wpwp)
       struc_wpwp(3) = two * ci * a3_wpwp

       struc_wwp(:) = struc_wwp(:) * couplWffsq
       struc_wpw(:) = struc_wpw(:) * couplWffsq
       struc_wpwp(:) = struc_wpwp(:) * couplWffsq

       helcoup(1:3,-1,-1) = helcoup(1:3,-1,-1) + &
                            struc_wpwp(1:3) * LCT(line1) * LCT(line2)/wpprop12/wpprop34 + &
                            struc_wpw(1:3)  * LCT(line1)             /wpprop12/iprop34  + &
                            struc_wwp(1:3)  *              LCT(line2)/iprop12/wpprop34

       helcoup(1:3,-1,+1) = helcoup(1:3,-1,+1) + &
                            struc_wpwp(1:3) * LCT(line1) * RCT(line2)/wpprop12/wpprop34 + &
                            struc_wpw(1:3)  * LCT(line1)             /wpprop12/iprop34  + &
                            struc_wwp(1:3)               * RCT(line2)/iprop12/wpprop34

       helcoup(1:3,+1,-1) = helcoup(1:3,+1,-1) + &
                            struc_wpwp(1:3) * RCT(line1) * LCT(line2)/wpprop12/wpprop34 + &
                            struc_wpw(1:3)  * RCT(line1)             /wpprop12/iprop34  + &
                            struc_wwp(1:3)               * LCT(line2)/iprop12/wpprop34

       helcoup(1:3,+1,+1) = helcoup(1:3,+1,+1) + &
                            struc_wpwp(1:3) * RCT(line1) * RCT(line2)/wpprop12/wpprop34 + &
                            struc_wpw(1:3)  * RCT(line1)             /wpprop12/iprop34  + &
                            struc_wwp(1:3)               * RCT(line2)/iprop12/wpprop34

    endif


    A0_WW_4f(-1,-1) = za(j1,j3)*zb(j4,j2) * helcoup(1,-1,-1) + &
         zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4) * helcoup(2,-1,-1) + &
         za(j1,j2)*za(j3,j4)*zb(j4,j2)**2 * helcoup(3,-1,-1)

    if ( includeVprime ) then
       A0_WW_4f(-1,+1) = za(j1,j4)*zb(j3,j2) * helcoup(1,-1,+1) + &
            zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3) * helcoup(2,-1,+1) + &
            za(j1,j2)*za(j4,j3)*zb(j3,j2)**2 * helcoup(3,-1,+1)

       A0_WW_4f(+1,-1) = za(j2,j3)*zb(j4,j1) * helcoup(1,+1,-1) + &
            zab2(j2,j3,j4,j1)*zab2(j3,j1,j2,j4) * helcoup(2,+1,-1) + &
            za(j2,j1)*za(j3,j4)*zb(j4,j1)**2 * helcoup(3,+1,-1)

       A0_WW_4f(+1,+1) = za(j2,j4)*zb(j3,j1) * helcoup(1,+1,+1) + &
            zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3) * helcoup(2,+1,+1) + &
            za(j2,j1)*za(j4,j3)*zb(j3,j1)**2 * helcoup(3,+1,+1)
    endif

    A0_WW_4f = A0_WW_4f/vev


    return

  end function A0_WW_4f


  !-- QCD amplitudes squared below
  subroutine me2_ggggh(j1,j2,j3,j4,za,zb,sprod,res)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: res
    complex(dp) :: a1234(-1:1,-1:1,-1:1,-1:1), a1324(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: aphi1234(1:2,-1:1,-1:1,-1:1,-1:1), aphi1324(1:2,-1:1,-1:1,-1:1,-1:1)
    real(dp), parameter :: col = 8.0_dp * Ca**3 * Cf
    integer :: i1,i2,i3,i4
    complex(dp) :: scalar, pseudo

    res = zero

    scalar = ghg2
    pseudo = -ghg4

    aphi1234 = A0phigggg_xxxx(j1,j2,j3,j4,za,zb,sprod)
    aphi1324 = A0phigggg_xxxx(j1,j3,j2,j4,za,zb,sprod)

    a1234(:,:,:,:) = scalar * (aPhi1234(1,:,:,:,:) + aPhi1234(2,:,:,:,:)) + &
         pseudo * (-ci) * (aPhi1234(1,:,:,:,:) - aPhi1234(2,:,:,:,:))

    a1324(:,:,:,:) = scalar * (aPhi1324(1,:,:,:,:) + aPhi1324(2,:,:,:,:)) + &
         pseudo * (-ci) * (aPhi1324(1,:,:,:,:) - aPhi1324(2,:,:,:,:))

    do i1 = -1,1,2
    do i2 = -1,1,2
    do i3 = -1,1,2
    do i4 = -1,1,2

       res = res + real(a1234(i1,i2,i3,i4)*conjg(a1234(i1,i2,i3,i4)),kind=dp)
       res = res + real(a1324(i1,i2,i3,i4)*conjg(a1324(i1,i2,i3,i4)),kind=dp)
       res = res + real(a1234(i1,i2,i3,i4)*conjg(a1324(i1,i3,i2,i4)),kind=dp)

    enddo
    enddo
    enddo
    enddo

    res = res * col / vev**2

    return

  end subroutine me2_ggggh

  subroutine me2_qbqggh(j1,j2,j3,j4,za,zb,sprod,res)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: res
    complex(dp) :: a1234(-1:1,-1:1,-1:1), a1243(-1:1,-1:1,-1:1)
    complex(dp) :: aphi1234(1:2,-1:1,-1:1,-1:1), aphi1243(1:2,-1:1,-1:1,-1:1)
    real(dp), parameter :: col1 = 4.0_dp * xn * Cf**2
    real(dp), parameter :: col2 = 2.0_dp * xn * Cf * (2.0_dp * Cf - Ca)
    integer :: i12,i3,i4
    complex(dp) :: scalar, pseudo

    res = zero

    scalar = ghg2
    pseudo = -ghg4

    aphi1234 = A0phiqbqgg_xxx(j1,j2,j3,j4,za,zb,sprod)
    aphi1243 = A0phiqbqgg_xxx(j1,j2,j4,j3,za,zb,sprod)

    a1234(:,:,:) = scalar * (aphi1234(1,:,:,:) + aphi1234(2,:,:,:)) + &
         pseudo * (-ci) * (aphi1234(1,:,:,:) - aphi1234(2,:,:,:))

    a1243(:,:,:) = scalar * (aphi1243(1,:,:,:) + aphi1243(2,:,:,:)) + &
         pseudo * (-ci) * (aphi1243(1,:,:,:) - aphi1243(2,:,:,:))

    do i12 = -1,1,2
    do i3 = -1,1,2
    do i4 = -1,1,2

       res = res + real(a1234(i12,i3,i4)*conjg(a1234(i12,i3,i4)),kind=dp) * col1
       res = res + real(a1243(i12,i4,i3)*conjg(a1243(i12,i4,i3)),kind=dp) * col1
       res = res + two * real(a1234(i12,i3,i4)*conjg(a1243(i12,i4,i3)),kind=dp) * col2

    enddo
    enddo
    enddo

    res = res / vev**2

    return

  end subroutine me2_qbqggh

   subroutine me2_qbqQBQ(j1,j2,j3,j4,za,zb,sprod,res_diff,res_id)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: res_diff, res_id
    complex(dp) :: amp_a(-1:1,-1:1), amp_b(-1:1,-1:1)
    complex(dp) :: aphi_a(1:2,-1:1,-1:1), aphi_b(1:2,-1:1,-1:1)
    real(dp), parameter :: col = xn**2 - one
    integer :: h12,h34
    complex(dp) :: scalar, pseudo

    res_diff = zero
    res_id = zero

    scalar = ghg2
    pseudo = -ghg4

    aphi_a = A0phiqbqQBQ_xx(j1,j2,j3,j4,za,zb,sprod)
    aphi_b = -A0phiqbqQBQ_xx(j1,j4,j3,j2,za,zb,sprod)

    amp_a(:,:) = scalar * (aphi_a(1,:,:) + aphi_a(2,:,:)) + &
         pseudo * (-ci) * (aphi_a(1,:,:) - aphi_a(2,:,:))

    amp_b(:,:) = scalar * (aphi_b(1,:,:) + aphi_b(2,:,:)) + &
         pseudo * (-ci) * (aphi_b(1,:,:) - aphi_b(2,:,:))

    do h12 = -1,1,2
    do h34 = -1,1,2
       res_diff = res_diff + real(amp_a(h12,h34)*conjg(amp_a(h12,h34)),kind=dp)
       res_id = res_id + real(amp_a(h12,h34)*conjg(amp_a(h12,h34)),kind=dp)
       res_id = res_id + real(amp_b(h12,h34)*conjg(amp_b(h12,h34)),kind=dp)
    enddo
    res_id = res_id - two/xn * real(amp_a(h12,h12)*conjg(amp_b(h12,h12)),kind=dp)
    enddo

    res_id = res_id * col / vev**2
    res_diff = res_diff * col / vev**2

    return

  end subroutine me2_qbqQBQ


  function A0phigggg_xxxx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phigggg_xxxx(1:2,-1:1,-1:1,-1:1,-1:1)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4),zb(4,4)
    real(dp) :: sprod(4,4)

    A0phigggg_xxxx = czero

    A0phigggg_xxxx(1,+1,+1,-1,-1) = A0phiggggmmpp(j3,j4,j1,j2,za,zb,sprod)
    A0phigggg_xxxx(1,+1,-1,+1,-1) = A0phiggggmpmp(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_xxxx(1,+1,-1,-1,+1) = A0phiggggmmpp(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_xxxx(1,+1,-1,-1,-1) = A0phiggggpmmm(j1,j2,j3,j4,za,zb,sprod)
    A0phigggg_xxxx(1,-1,+1,+1,-1) = A0phiggggmmpp(j4,j1,j2,j3,za,zb,sprod)
    A0phigggg_xxxx(1,-1,+1,-1,+1) = A0phiggggmpmp(j1,j2,j3,j4,za,zb,sprod)
    A0phigggg_xxxx(1,-1,+1,-1,-1) = A0phiggggpmmm(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_xxxx(1,-1,-1,+1,+1) = A0phiggggmmpp(j1,j2,j3,j4,za,zb,sprod)
    A0phigggg_xxxx(1,-1,-1,+1,-1) = A0phiggggpmmm(j3,j4,j1,j2,za,zb,sprod)
    A0phigggg_xxxx(1,-1,-1,-1,+1) = A0phiggggpmmm(j4,j1,j2,j3,za,zb,sprod)
    A0phigggg_xxxx(1,-1,-1,-1,-1) = A0phiggggmmmm(j1,j2,j3,j4,za,zb,sprod)

    A0phigggg_xxxx(2,-1,-1,+1,+1) = conjg(A0phigggg_xxxx(1,+1,+1,-1,-1))
    A0phigggg_xxxx(2,-1,+1,-1,+1) = conjg(A0phigggg_xxxx(1,+1,-1,+1,-1))
    A0phigggg_xxxx(2,-1,+1,+1,-1) = conjg(A0phigggg_xxxx(1,+1,-1,-1,+1))
    A0phigggg_xxxx(2,-1,+1,+1,+1) = conjg(A0phigggg_xxxx(1,+1,-1,-1,-1))
    A0phigggg_xxxx(2,+1,-1,-1,+1) = conjg(A0phigggg_xxxx(1,-1,+1,+1,-1))
    A0phigggg_xxxx(2,+1,-1,+1,-1) = conjg(A0phigggg_xxxx(1,-1,+1,-1,+1))
    A0phigggg_xxxx(2,+1,-1,+1,+1) = conjg(A0phigggg_xxxx(1,-1,+1,-1,-1))
    A0phigggg_xxxx(2,+1,+1,-1,-1) = conjg(A0phigggg_xxxx(1,-1,-1,+1,+1))
    A0phigggg_xxxx(2,+1,+1,-1,+1) = conjg(A0phigggg_xxxx(1,-1,-1,+1,-1))
    A0phigggg_xxxx(2,+1,+1,+1,-1) = conjg(A0phigggg_xxxx(1,-1,-1,-1,+1))
    A0phigggg_xxxx(2,+1,+1,+1,+1) = conjg(A0phigggg_xxxx(1,-1,-1,-1,-1))

    return

  end function A0phigggg_xxxx

  function A0phiqbqgg_xxx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phiqbqgg_xxx(1:2,-1:1,-1:1,-1:1)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    real(dp) :: s3
    complex(dp) :: zab2

    s3(j1,j2,j3) = sprod(j1,j2) + sprod(j1,j3) + sprod(j2,j3)
    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4)

    A0phiqbqgg_xxx = czero

    A0phiqbqgg_xxx(1,-1,+1,-1) = -za(j1,j4)**2*za(j2,j4)/(za(j1,j2)*za(j2,j3)*za(j3,j4))
    A0phiqbqgg_xxx(1,-1,-1,+1) = za(j1,j3)**3/(za(j1,j2)*za(j3,j4)*za(j4,j1))
    A0phiqbqgg_xxx(1,-1,-1,-1) =  -((za(j1,j3)*zab2(j4,j1,j3,j2)**2)/ &
         (s3(j1,j2,j3)*sprod(j1,j2)*zb(j2,j3))) - &
         ((1/sprod(j1,j2) + 1/sprod(j4,j1))*za(j4,j1)*zab2(j3,j1,j4,j2)**2)/ &
         (s3(j4,j1,j2)*zb(j2,j4)) + &
         zab2(j1,j3,j4,j2)**2/(za(j1,j2)*zb(j2,j3)*zb(j2,j4)*zb(j3,j4))

    A0phiqbqgg_xxx(2,-1,+1,-1) = -zb(j1,j3)*zb(j2,j3)**2/(zb(j1,j2)*zb(j3,j4)*zb(j4,j1))
    A0phiqbqgg_xxx(2,-1,-1,+1) = zb(j2,j4)**3/(zb(j1,j2)*zb(j2,j3)*zb(j3,j4))
    A0phiqbqgg_xxx(2,-1,+1,+1) = -(zab2(j1,j3,j4,j2)**2/(za(j1,j3)*za(j1,j4)*za(j3,j4)*zb(j1,j2))) - &
         ((1/sprod(j1,j2) + 1/sprod(j2,j3))*zab2(j1,j2,j3,j4)**2*zb(j2,j3))/ &
         (s3(j1,j2,j3)*za(j1,j3)) + &
         (zab2(j1,j2,j4,j3)**2*zb(j2,j4))/ &
         (s3(j4,j1,j2)*sprod(j1,j2)*za(j1,j4))

    A0phiqbqgg_xxx(1,+1,+1,-1) = conjg(A0phiqbqgg_xxx(2,-1,-1,+1))
    A0phiqbqgg_xxx(1,+1,-1,+1) = conjg(A0phiqbqgg_xxx(2,-1,+1,-1))
    A0phiqbqgg_xxx(1,+1,-1,-1) = conjg(A0phiqbqgg_xxx(2,-1,+1,+1))

    A0phiqbqgg_xxx(2,+1,+1,-1) = conjg(A0phiqbqgg_xxx(1,-1,-1,+1))
    A0phiqbqgg_xxx(2,+1,-1,+1) = conjg(A0phiqbqgg_xxx(1,-1,+1,-1))
    A0phiqbqgg_xxx(2,+1,+1,+1) = conjg(A0phiqbqgg_xxx(1,-1,-1,-1))

    return

  end function A0phiqbqgg_xxx

  function A0phiqbqQBQ_xx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phiqbqQBQ_xx(1:2,-1:1,-1:1)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)

    A0phiqbqQBQ_xx = czero

    A0phiqbqQBQ_xx(1,-1,+1) = za(j1,j4)**2/(za(j1,j2)*za(j3,j4))
    A0phiqbqQBQ_xx(1,-1,-1) = za(j1,j3)**2/(za(j1,j2)*za(j4,j3))

    A0phiqbqQBQ_xx(2,-1,+1) = zb(j2,j3)**2/(zb(j1,j2)*zb(j3,j4))
    A0phiqbqQBQ_xx(2,-1,-1) = zb(j2,j4)**2/(zb(j1,j2)*zb(j4,j3))

    A0phiqbqQBQ_xx(1,+1,-1) = conjg(A0phiqbqQBQ_xx(2,-1,+1))
    A0phiqbqQBQ_xx(1,+1,+1) = conjg(A0phiqbqQBQ_xx(2,-1,-1))

    A0phiqbqQBQ_xx(2,+1,-1) = conjg(A0phiqbqQBQ_xx(1,-1,+1))
    A0phiqbqQBQ_xx(2,+1,+1) = conjg(A0phiqbqQBQ_xx(1,-1,-1))

    return

  end function A0phiqbqQBQ_xx


  function A0phiggggpmmm(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0phiggggpmmm

    real(dp) :: s3
    complex(dp) :: zab2

    s3(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

    A0phiggggpmmm = &
         +(zab2(j3,j2,j4,j1)*za(j2,j4))**2/(s3(j1,j2,j4)*sprod(j1,j2)*sprod(j1,j4)) &
         +(zab2(j4,j2,j3,j1)*za(j2,j3))**2/(s3(j1,j2,j3)*sprod(j1,j2)*sprod(j2,j3)) &
         +(zab2(j2,j3,j4,j1)*za(j3,j4))**2/(s3(j1,j3,j4)*sprod(j1,j4)*sprod(j3,j4)) &
         -za(j2,j4)/(za(j1,j2)*zb(j2,j3)*zb(j3,j4)*za(j4,j1)) &
         *(-sprod(j2,j3)*zab2(j2,j3,j4,j1)/zb(j4,j1) &
         -sprod(j3,j4)*zab2(j4,j2,j3,j1)/zb(j1,j2) &
         -s3(j2,j3,j4)*za(j2,j4))

    return

  end function A0phiggggpmmm


  function A0phiggggmpmp(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0phiggggmpmp

    A0phiggggmpmp = za(j1,j3)**4 &
         /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))

    return

  end function A0phiggggmpmp


  function A0phiggggmmpp(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0phiggggmmpp

    A0phiggggmmpp = za(j1,j2)**4 &
         /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))

    return

  end function A0phiggggmmpp


  function A0phiggggmmmm(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4),qsq
    complex(dp) :: A0phiggggmmmm

    qsq = sprod(j1,j2)+sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4)+sprod(j3,j4)
    A0phiggggmmmm=qsq**2/(zb(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j4,j1))

    return

  end function A0phiggggmmmm

  !-------------------------------------------------------------------------
  !-- generic functions below
  !- MCFM spinors
  subroutine spinoru2(n,p,za,zb,s)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,n)
    complex(dp), intent(out) :: za(n,n), zb(n,n)
    real(dp), intent(out) :: s(n,n)
    integer :: i,j
    complex(dp) :: c23(n), f(n)
    real(dp) :: rt(n)

    !---if one of the vectors happens to be zero this routine fails.
    do j=1,N
       za(j,j)=czero
       zb(j,j)=za(j,j)

       !-----positive energy case
       if (p(1,j) .gt. zero) then
          rt(j)=sqrt(abs(p(2,j)+p(1,j)))
          c23(j)=cmplx(p(4,j),-p(3,j),kind=dp)
          f(j)=(one,zero)
       else
       !-----negative energy case
          rt(j)=sqrt(abs(-p(1,j)-p(2,j)))
          c23(j)=cmplx(-p(4,j),p(3,j),kind=dp)
          f(j)=ci
       endif
    enddo

    do i=2,N

     do j=1,i-1
          s(i,j)=two*scr(p(:,i),p(:,j))
          za(i,j)=f(i)*f(j)  * ( c23(i)*cmplx(rt(j)/(rt(i)+1d-16),kind=dp)-c23(j)*cmplx(rt(i)/(rt(j)+1d-16),kind=dp) )

          if (abs(s(i,j)).lt.1d-5) then
             zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
          else
             zb(i,j)=-cmplx(s(i,j),kind=dp)/(za(i,j)+1d-16)
          endif

          za(j,i)=-za(i,j)
          zb(j,i)=-zb(i,j)
          s(j,i)=s(i,j)

       enddo

    enddo

    return

  end subroutine spinoru2



end module modHiggsJJ
