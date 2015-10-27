module modHiggsJJ
  use modParameters
  implicit none
  private

  public :: EvalAmp_WBFH_UnSymm_SA,EvalAmp_WBFH_UnSymm_SA_Select
  public :: EvalAmp_SBFH_UnSymm_SA
  public :: get_VBFchannelHash,get_GENchannelHash

  !-- general definitions, to be merged with Markus final structure
  integer, parameter  :: dp = selected_real_kind(15)
  
  real(dp), parameter :: nf = 5.0_dp
  real(dp), parameter :: xw = sitW**2
  real(dp), parameter :: twosc = sqrt(4.0_dp*xw*(1.0_dp-xw))
  real(dp), parameter :: gs = sqrt(alphas*4.0_dp*pi)

  integer :: pdfGlu_ = 0
  integer :: pdfDn_ = 1
  integer :: pdfUp_ = 2
  integer :: pdfStr_ = 3
  integer :: pdfChm_ = 4
  integer :: pdfBot_ = 5
  integer :: pdfADn_ = -1
  integer :: pdfAUp_ = -2
  integer :: pdfAStr_ = -3
  integer :: pdfAChm_ = -4
  integer :: pdfABot_ = -5

  real(dp), parameter :: tag1 = 1.0_dp
  real(dp), parameter :: tag2 = 1.0_dp
  real(dp), parameter :: tagbot = 1.0_dp

  real(dp), parameter :: xn = 3.0_dp
  real(dp), parameter :: Ca = 3.0_dp
  real(dp), parameter :: Cf = 4.0_dp/3.0_dp

  real(dp), parameter :: avegg = 1.0_dp/4.0_dp/64.0_dp
  real(dp), parameter :: aveqg = 1.0_dp/4.0_dp/24.0_dp
  real(dp), parameter :: aveqq = 1.0_dp/4.0_dp/9.0_dp

  real(dp), parameter :: mwsq = m_w**2
  real(dp), parameter :: mzsq = m_z**2

  
  
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
      
      ijSel( 71:,:)  = 0
      

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
  
    
  
  
  !-- SM: |g2| = alphas/(six*pi)
  !-- g3 not supported yet
  subroutine EvalAmp_SBFH_UnSymm_SA(p,ggcoupl,res)
    real(dp), intent(in) :: p(4,5)
    complex(dp), intent(in) :: ggcoupl(2:4)
    real(dp), intent(out) :: res(-5:5,-5:5)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp, restmpid
    integer :: i, j

    res = zero

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    !-- gg -> gg
    call me2_ggggh(ggcoupl,1,2,3,4,za,zb,sprod,restmp)
    restmp = restmp * avegg * SymmFac
    res(0,0) = res(0,0) + restmp

    !-- gg -> qqb
    call me2_qbqggh(ggcoupl,4,3,1,2,za,zb,sprod,restmp)
    restmp = restmp * avegg
    res(0,0) = res(0,0) + restmp * nf

    !-- gq -> gq
    call me2_qbqggh(ggcoupl,2,4,1,3,za,zb,sprod,restmp)
    restmp = restmp * aveqg
    do i = 1,5
       res(0,i) = res(0,i) + restmp
       res(0,-i) = res(0,-i) + restmp
    enddo
    call me2_qbqggh(ggcoupl,1,4,2,3,za,zb,sprod,restmp)
    restmp = restmp * aveqg
    do i = 1,5
       res(i,0) = res(i,0) + restmp
       res(-i,0) = res(-i,0) + restmp
    enddo

    !-- qqb -> gg
    call me2_qbqggh(ggcoupl,1,2,3,4,za,zb,sprod,restmp)
    restmp = restmp * aveqq * SymmFac
    do i = 1,5
       res(i,-i) = res(i,-i) + restmp
       res(-i,i) = res(-i,i) + restmp
    enddo

    !-- qqb -> qqb
    call me2_qbqQBQ(ggcoupl,1,2,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq
    do i = 1,5
       res(i,-i) =res(i,-i) + restmp
    enddo
    call me2_qbqQBQ(ggcoupl,2,1,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq
    do i = 1,5
       res(-i,i) =res(-i,i) + restmp
    enddo

    !-- qqb -> rrb
    call me2_qbqQBQ(ggcoupl,1,2,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       res(i,-i) =res(i,-i) + restmp * (nf-1.0_dp)
       res(-i,i) = res(-i,i) + restmp * (nf-1.0_dp)
    enddo

    !-- qrb -> qrb
    call me2_qbqQBQ(ggcoupl,1,3,4,2,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       do j = 1,5
          if (i.ne.j) res(i,-j) = res(i,-j) + restmp
       enddo
    enddo
    call me2_qbqQBQ(ggcoupl,2,3,4,1,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       do j = 1,5
          if (i.ne.j) res(-j,i) = res(-j,i) + restmp
       enddo
    enddo

    !-- qq -> qq
    call me2_qbqQBQ(ggcoupl,1,3,2,4,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq * SymmFac
    do i = 1,5
       res(i,i) = res(i,i) + restmp
       res(-i,-i) = res(-i,-i) + restmp
    enddo

    !-- qr -> qr
    call me2_qbqQBQ(ggcoupl,1,3,2,4,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       do j = 1,5
          if (i.ne.j) then
             res(i,j) = res(i,j) + restmp
             res(-i,-j) = res(-i,-j) + restmp
          endif
       enddo
    enddo

    return

  end subroutine EvalAmp_SBFH_UnSymm_SA








  subroutine EvalAmp_WBFH_UnSymm_SA(p,vvcoupl,wwcoupl,res)
    real(dp), intent(in) :: p(4,5)
    complex(dp), intent(in) :: vvcoupl(4),wwcoupl(4)! wwcoupl is only used (for H->WW) if it includes non-zero values
    real(dp), intent(out) :: res(-5:5,-5:5)
    complex(dp) :: amp_z(-1:1,-1:1), amp_z_b(-1:1,-1:1)
    complex(dp) :: amp_w(-1:1,-1:1)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp), parameter :: Lu = aL_QUp**2, Ru = aR_QUp**2
    real(dp), parameter :: Ld = aL_QDn**2, Rd = aR_QDn**2
    real(dp), parameter :: couplz = gwsq * xw/twosc**2         !*M_Z/sitW/dsqrt(1d0-sitW**2)! MARKUS ADDED HVV COUPLINGS: CHECK!!
    real(dp), parameter :: couplw = gwsq/two                   !*M_W/sitW
    real(dp) :: restmp=0d0
    integer :: i, j, j1, j2, iflip, pdfindex(2)

    res = zero

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    !-- qq->qq, up
    amp_z = A0_VV_4f(4,1,3,2,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(3,1,4,2,vvcoupl,za,zb,sprod,m_z,ga_z)
    
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
    amp_z = A0_VV_4f(1,4,2,3,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(1,3,2,4,vvcoupl,za,zb,sprod,m_z,ga_z)

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

    ! NEED TO ADJSUT J1,J2 in W AMPS HERE
    !-- W/Z interference
    j1 = 1
    j2 = 2
    do iflip = 1, 2 
       !-- ud -> ud
       amp_z = A0_VV_4f(4,j2,3,j1,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(4,j1,3,j2,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(4,j1,3,j2,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
            (abs(amp_z(-1,+1))**2) * Ld * Ru + &
            (abs(amp_z(+1,-1))**2) * Rd * Lu + &
            (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2  !*00  !MARKUS

       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * xn**2

       restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
            aL_QUp * aL_QDn * couplz * couplw * xn  ! *00000000d0 !MARKUS: switching off the interference

       restmp = restmp * aveqq

       pdfindex = flip(iflip,pdfUp_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp
       pdfindex = flip(iflip,pdfChm_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       !-- ubdb -> ubdb
       amp_z = A0_VV_4f(j2,4,j1,3,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(j1,4,j2,3,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(j1,4,j2,3,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = A0_VV_4f(3,j1,j2,4,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif
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
       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = A0_VV_4f(3,j2,4,j1,vvcoupl,za,zb,sprod,m_w,ga_w)
!           amp_w = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_w,ga_w)! MARKUS
       else
          amp_w = A0_VV_4f(3,j2,4,j1,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
!           amp_w = A0_VV_4f(3,j1,4,j2,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)! MARKUS
       endif

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
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1  *1d0

       pdfindex = flip(iflip,pdfUp_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       pdfindex = flip(iflip,pdfChm_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
       
       restmp = restmp + abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2  *1d0

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
       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = A0_VV_4f(j1,4,j2,3,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = A0_VV_4f(j1,4,j2,3,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

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
















  SUBROUTINE EvalAmp_WBFH_UnSymm_SA_Select(p,vvcoupl,wwcoupl,iSel,jSel,zz_fusion,res)
  implicit none
    real(dp), intent(in) :: p(4,5)
    complex(dp), intent(in) :: vvcoupl(4),wwcoupl(4)! wwcoupl is only used (for H->WW) if it includes non-zero values
    real(dp), intent(out) :: res(-5:5,-5:5)
    logical, intent(in) :: zz_fusion
    integer, intent(in) :: iSel,jSel
    complex(dp) :: amp_z(-1:1,-1:1), amp_z_b(-1:1,-1:1)
    complex(dp) :: amp_w(-1:1,-1:1)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp), parameter :: Lu = aL_QUp**2, Ru = aR_QUp**2
    real(dp), parameter :: Ld = aL_QDn**2, Rd = aR_QDn**2
    real(dp), parameter :: couplz = gwsq * xw/twosc**2        
    real(dp), parameter :: couplw = gwsq/two                  
    real(dp) :: restmp
    integer :: i, j, j1, j2, iflip, pdfindex(2)

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)
    
    
    
    
    
    
    !-- qq->qq, up
if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfUp_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfChm_) ) then

    amp_z = A0_VV_4f(4,1,3,2,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(3,1,4,2,vvcoupl,za,zb,sprod,m_z,ga_z)
    
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Lu**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Lu * Ru + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Lu * Ru + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Ru**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Lu**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Ru**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    if( .not. zz_fusion ) restmp = 0.0d0 

    res(pdfUp_,pdfUp_) = restmp
    res(pdfChm_,pdfChm_) = restmp    

return
endif









    !-- qq->qq, down
if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfDn_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfStr_) .or. (iSel.eq.pdfBot_ .and. jSel.eq.pdfBot_) ) then
    amp_z = A0_VV_4f(4,1,3,2,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(3,1,4,2,vvcoupl,za,zb,sprod,m_z,ga_z)
        
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Ld**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Ld * Rd + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Ld * Rd + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Rd**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Ld**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Rd**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    if( .not. zz_fusion ) restmp = 0.0d0 
    
    res(pdfDn_,pdfDn_) = restmp
    res(pdfStr_,pdfStr_) = restmp
    res(pdfBot_,pdfBot_) = restmp * tagbot   
return
endif









    !-- qbqb->qbqb, aup
if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfAUp_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfAChm_) ) then
    amp_z = A0_VV_4f(1,4,2,3,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(1,3,2,4,vvcoupl,za,zb,sprod,m_z,ga_z)

    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Lu**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Lu * Ru + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Lu * Ru + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Ru**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Lu**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Ru**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    if( .not. zz_fusion ) restmp = 0.0d0 
    

    res(pdfAUp_,pdfAUp_) = restmp
    res(pdfAChm_,pdfAChm_) = restmp   
return
endif







    !-- qbqb->qbqb, adn
if( (iSel.eq.pdfADn_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfAStr_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfABot_ .and. jSel.eq.pdfABot_) ) then
    amp_z = A0_VV_4f(1,4,2,3,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(1,3,2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Ld**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Ld * Rd + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Ld * Rd + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Rd**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Ld**2 + &
            two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Rd**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

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
       amp_z = A0_VV_4f(4,j2,3,j1,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(4,j1,3,j2,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(4,j1,3,j2,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
                (abs(amp_z(-1,+1))**2) * Ld * Ru + &
                (abs(amp_z(+1,-1))**2) * Rd * Lu + &
                (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2 
       else
          restmp = abs(amp_w(-1,-1))**2 * couplw**2 * xn**2
          restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
                   aL_QUp * aL_QDn * couplz * couplw * xn          
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
       amp_z = A0_VV_4f(j2,4,j1,3,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(j1,4,j2,3,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(j1,4,j2,3,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
                (abs(amp_z(-1,+1))**2) * Ld * Ru + &
                (abs(amp_z(+1,-1))**2) * Rd * Lu + &
                (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2
       else
          restmp = abs(amp_w(-1,-1))**2 * couplw**2 * xn**2      
          restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
                   aL_QUp * aL_QDn * couplz * couplw * xn
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
       amp_z = A0_VV_4f(4,j2,3,j1,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(4,j1,3,j2,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(4,j1,3,j2,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
                (abs(amp_z(-1,+1))**2) * Ld * Ru + &
                (abs(amp_z(+1,-1))**2) * Rd * Lu + &
                (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2 
       else
            restmp = abs(amp_w(-1,-1))**2 * couplw**2 * xn**2
            restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
                  aL_QUp * aL_QDn * couplz * couplw * xn
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
       amp_z = A0_VV_4f(j2,4,j1,3,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(j1,4,j2,3,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(j1,4,j2,3,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
       endif

       if( zz_fusion ) then
          restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
                (abs(amp_z(-1,+1))**2) * Ld * Ru + &
                (abs(amp_z(+1,-1))**2) * Rd * Lu + &
                (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * xn**2
       else
          restmp= abs(amp_w(-1,-1))**2 * couplw**2 * xn**2
          restmp = restmp + two * real(amp_z(-1,-1)*conjg(amp_w(-1,-1)),kind=dp) * &
                aL_QUp * aL_QDn * couplz * couplw * xn
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
          amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
          restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
                (abs(amp_z(-1,+1))**2) * Lu * Ru + &
                (abs(amp_z(+1,-1))**2) * Ru * Lu + &
                (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg * tag1
       else            
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_w,ga_w)
          else
              amp_w = A0_VV_4f(3,j1,j2,4,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
          endif
          restmp = abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg
       
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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
            (abs(amp_z(-1,+1))**2) * Ld * Ru + &
            (abs(amp_z(+1,-1))**2) * Rd * Lu + &
            (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * SpinAvg

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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)

       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg * tag1

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



! MARKUS: why is bot abot not here????

if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfADn_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfAStr_)  .or. (iSel.eq.pdfDn_ .and. jSel.eq.pdfAStr_)  .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfADn_) ) then
       if( zz_fusion ) then
          amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)      
          restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
                (abs(amp_z(-1,+1))**2) * Ld * Rd + &
                (abs(amp_z(+1,-1))**2) * Rd * Ld + &
                (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg * tag1
       else
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_w,ga_w)
          else
              amp_w = A0_VV_4f(3,j1,j2,4,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
          endif            
          restmp = abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
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
       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       !--uc -> uc
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfUp_,pdfChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
return
endif







       !--us -> us/cd
if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfBot_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfBot_) ) then

       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfUp_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfChm_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
return
endif





if( (iSel.eq.pdfUp_ .and. jSel.eq.pdfStr_) .or. (iSel.eq.pdfChm_ .and. jSel.eq.pdfDn_) ) then

       if( zz_fusion ) then
          amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
          restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
                (abs(amp_z(-1,+1))**2) * Lu * Rd + &
                (abs(amp_z(+1,-1))**2) * Ru * Ld + &
                (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1
       else
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(3,j2,4,j1,vvcoupl,za,zb,sprod,m_w,ga_w)
    !           amp_w = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_w,ga_w)! MARKUS
          else
              amp_w = A0_VV_4f(3,j2,4,j1,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
    !           amp_w = A0_VV_4f(3,j1,4,j2,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)! MARKUS
          endif            
          restmp = abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
       endif
       
       pdfindex = flip(iflip,pdfUp_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp


return
endif







       !--ds -> ds
if( (iSel.eq.pdfDn_ .and. jSel.eq.pdfStr_) .or. (iSel.eq.pdfDn_ .and. jSel.eq.pdfBot_) .or. (iSel.eq.pdfStr_ .and. jSel.eq.pdfBot_) ) then

       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg

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
       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       !--ubcb -> ubcb
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfAUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
return
endif






       !--ubsb -> ubsb//cbdb
if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfABot_) ) then
       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfAUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfAChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
return
endif



if( (iSel.eq.pdfAUp_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfAChm_ .and. jSel.eq.pdfADn_) ) then

       if( zz_fusion ) then
          amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
          restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
                (abs(amp_z(-1,+1))**2) * Lu * Rd + &
                (abs(amp_z(+1,-1))**2) * Ru * Ld + &
                (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1
       else
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(j1,4,j2,3,vvcoupl,za,zb,sprod,m_w,ga_w)
          else
              amp_w = A0_VV_4f(j1,4,j2,3,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
          endif           
          restmp= abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
       endif
              
       pdfindex = flip(iflip,pdfAUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfAChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

       
return
endif







       !--dbsb -> dbsb
if( (iSel.eq.pdfADn_ .and. jSel.eq.pdfAStr_) .or. (iSel.eq.pdfADn_ .and. jSel.eq.pdfABot_) .or. (iSel.eq.pdfAStr_ .and. jSel.eq.pdfABot_) ) then

       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg

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
          amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
          restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
                (abs(amp_z(-1,+1))**2) * Lu * Ru + &
                (abs(amp_z(+1,-1))**2) * Ru * Lu + &
                (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg * tag1
       else
            if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
                amp_w = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_w,ga_w)
            else
                amp_w = A0_VV_4f(3,j1,j2,4,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
            endif            
            restmp= abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg
      
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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Lu + &
            (abs(amp_z(-1,+1))**2) * Ld * Ru + &
            (abs(amp_z(+1,-1))**2) * Rd * Lu + &
            (abs(amp_z(+1,+1))**2) * Rd * Ru) * couplz**2 * SpinAvg

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
       amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)

       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg * tag1

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
          amp_z = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
          restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
                (abs(amp_z(-1,+1))**2) * Ld * Rd + &
                (abs(amp_z(+1,-1))**2) * Rd * Ld + &
                (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg * tag1
       else
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(3,j1,j2,4,vvcoupl,za,zb,sprod,m_w,ga_w)
          else
              amp_w = A0_VV_4f(3,j1,j2,4,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
          endif
          restmp= abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
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
       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       !--uc -> uc
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfUp_,pdfChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
return
endif






       !--us -> us/cd
if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfBot_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfBot_) ) then

       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfUp_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfChm_,pdfBot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
return
endif





if( (jSel.eq.pdfUp_ .and. iSel.eq.pdfStr_) .or. (jSel.eq.pdfChm_ .and. iSel.eq.pdfDn_) ) then

       if( zz_fusion ) then
          amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
          restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
                (abs(amp_z(-1,+1))**2) * Lu * Rd + &
                (abs(amp_z(+1,-1))**2) * Ru * Ld + &
                (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1
       else
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(3,j2,4,j1,vvcoupl,za,zb,sprod,m_w,ga_w)
    !           amp_w = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_w,ga_w)! MARKUS
          else
              amp_w = A0_VV_4f(3,j2,4,j1,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
    !           amp_w = A0_VV_4f(3,j1,4,j2,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)! MARKUS
          endif
          restmp= abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
       endif
       
       pdfindex = flip(iflip,pdfUp_,pdfStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfChm_,pdfDn_)
       res(pdfindex(1),pdfindex(2)) = restmp

return
endif







       !--ds -> ds
if( (jSel.eq.pdfDn_ .and. iSel.eq.pdfStr_) .or. (jSel.eq.pdfDn_ .and. iSel.eq.pdfBot_) .or. (jSel.eq.pdfStr_ .and. iSel.eq.pdfBot_) ) then

       amp_z = A0_VV_4f(3,j1,4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)      
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg

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
       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       !--ubcb -> ubcb
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Lu + &
            (abs(amp_z(-1,+1))**2) * Lu * Ru + &
            (abs(amp_z(+1,-1))**2) * Ru * Lu + &
            (abs(amp_z(+1,+1))**2) * Ru * Ru) * couplz**2 * SpinAvg

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfAUp_,pdfAChm_)
       res(pdfindex(1),pdfindex(2)) = restmp
return
endif






       !--ubsb -> ubsb//cbdb
if( (jSel.eq.pdfAUp_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfAChm_ .and. iSel.eq.pdfABot_) ) then
       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
            (abs(amp_z(-1,+1))**2) * Lu * Rd + &
            (abs(amp_z(+1,-1))**2) * Ru * Ld + &
            (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1

       if( .not. zz_fusion ) restmp = 0.0d0 
    

       pdfindex = flip(iflip,pdfAUp_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot

       pdfindex = flip(iflip,pdfAChm_,pdfABot_)
       res(pdfindex(1),pdfindex(2)) = restmp * tagbot
return
endif



if( (jSel.eq.pdfAUp_ .and. iSel.eq.pdfAStr_) .or. (jSel.eq.pdfAChm_ .and. iSel.eq.pdfADn_) ) then
       if( zz_fusion ) then
          amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
          restmp = ((abs(amp_z(-1,-1))**2) * Lu * Ld + &
                (abs(amp_z(-1,+1))**2) * Lu * Rd + &
                (abs(amp_z(+1,-1))**2) * Ru * Ld + &
                (abs(amp_z(+1,+1))**2) * Ru * Rd) * couplz**2 * SpinAvg * tag1
       else
          if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
              amp_w = A0_VV_4f(j1,4,j2,3,vvcoupl,za,zb,sprod,m_w,ga_w)
          else
              amp_w = A0_VV_4f(j1,4,j2,3,wwcoupl,za,zb,sprod,m_w,ga_w,useWWcoupl=.true.)
          endif          
          restmp= abs(amp_w(-1,-1))**2 * couplw**2 * SpinAvg * tag2
       endif
       
       pdfindex = flip(iflip,pdfAUp_,pdfAStr_)
       res(pdfindex(1),pdfindex(2)) = restmp

       pdfindex = flip(iflip,pdfAChm_,pdfADn_)
       res(pdfindex(1),pdfindex(2)) = restmp

return
endif







       !--dbsb -> dbsb
if( (jSel.eq.pdfADn_ .and. iSel.eq.pdfAStr_) .or. (jSel.eq.pdfADn_ .and. iSel.eq.pdfABot_) .or. (jSel.eq.pdfAStr_ .and. iSel.eq.pdfABot_) ) then

       amp_z = A0_VV_4f(j1,3,j2,4,vvcoupl,za,zb,sprod,m_z,ga_z)
       restmp = ((abs(amp_z(-1,-1))**2) * Ld * Ld + &
            (abs(amp_z(-1,+1))**2) * Ld * Rd + &
            (abs(amp_z(+1,-1))**2) * Rd * Ld + &
            (abs(amp_z(+1,+1))**2) * Rd * Rd) * couplz**2 * SpinAvg

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
  END SUBROUTINE 









  
  
  
  
  





  function flip(i,a1,a2)
    integer :: flip(2)
    integer :: i, a1, a2
    
    if (i .eq. 1) flip = (/a1,a2/)
    if (i .eq. 2) flip = (/a2,a1/)

    return

  end function flip


  !-------------------------------------------------------------------------
  !-- amplitudes below

  function A0_VV_4f(j1,j2,j3,j4,vvcoupl,za,zb,sprod,mv,ga_v,useWWcoupl)
    complex(dp) :: A0_VV_4f(-1:1,-1:1)
    integer :: j1,j2,j3,j4
    complex(dp) :: vvcoupl(4), za(4,4), zb(4,4)
    real(dp) :: mv, ga_v
    logical,optional :: useWWcoupl
    real(dp) :: sprod(4,4)
    real(dp) :: mhsq, q1q2, kcoupl
    complex(dp) :: a1, a2, a3, struc1, struc2, struc3
    complex(dp) :: zab2
    complex(dp) :: iprop12, iprop34
    complex(dp) :: vvcoupl_prime(4)

    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4)

    A0_VV_4f = czero

    q1q2 = (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))/two
    mhsq = two * q1q2 + sprod(j1,j2) + sprod(j3,j4)

    kcoupl = q1q2/lambda**2

    if( .not.present(useWWcoupl) ) then
     vvcoupl_prime(1) = vvcoupl(1)   &
        + ghz1_prime * lambda_z1**4/(lambda_z1**2 + abs(sprod(j1,j2)))/(lambda_z1**2 + abs(sprod(j3,j4)))  &
        + ghz1_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_z1**2  &
        + ghz1_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_z1**2  &
        + ghz1_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghz1_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_z1**4  &
        + ghz1_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_z1**4  &
        + ghz1_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_z1**4 

     vvcoupl_prime(2) = vvcoupl(2)   &
        + ghz2_prime * lambda_z2**4/(lambda_z2**2 + abs(sprod(j1,j2)))/(lambda_z2**2 + abs(sprod(j3,j4)))  &
        + ghz2_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_z2**2  &
        + ghz2_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_z2**2  &
        + ghz2_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghz2_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_z2**4  &
        + ghz2_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_z2**4  &
        + ghz2_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_z2**4 

     vvcoupl_prime(3) = vvcoupl(3)   &
        + ghz3_prime * lambda_z3**4/(lambda_z3**2 + abs(sprod(j1,j2)))/(lambda_z3**2 + abs(sprod(j3,j4)))  &
        + ghz3_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_z3**2  &
        + ghz3_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_z3**2  &
        + ghz3_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghz3_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_z3**4  &
        + ghz3_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_z3**4  &
        + ghz3_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_z3**4 

     vvcoupl_prime(4) = vvcoupl(4)   &
        + ghz4_prime * lambda_z4**4/(lambda_z4**2 + abs(sprod(j1,j2)))/(lambda_z4**2 + abs(sprod(j3,j4)))  &
        + ghz4_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_z4**2  &
        + ghz4_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_z4**2  &
        + ghz4_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghz4_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_z4**4  &
        + ghz4_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_z4**4  &
        + ghz4_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_z4**4 

    else

     vvcoupl_prime(1) = vvcoupl(1)   &
        + ghw1_prime * lambda_w1**4/(lambda_w1**2 + abs(sprod(j1,j2)))/(lambda_w1**2 + abs(sprod(j3,j4)))  &
        + ghw1_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_w1**2  &
        + ghw1_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_w1**2  &
        + ghw1_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghw1_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_w1**4  &
        + ghw1_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_w1**4  &
        + ghw1_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_w1**4 

     vvcoupl_prime(2) = vvcoupl(2)   &
        + ghw2_prime * lambda_w2**4/(lambda_w2**2 + abs(sprod(j1,j2)))/(lambda_w2**2 + abs(sprod(j3,j4)))  &
        + ghw2_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_w2**2  &
        + ghw2_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_w2**2  &
        + ghw2_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghw2_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_w2**4  &
        + ghw2_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_w2**4  &
        + ghw2_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_w2**4 

     vvcoupl_prime(3) = vvcoupl(3)   &
        + ghw3_prime * lambda_w3**4/(lambda_w3**2 + abs(sprod(j1,j2)))/(lambda_w3**2 + abs(sprod(j3,j4)))  &
        + ghw3_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_w3**2  &
        + ghw3_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_w3**2  &
        + ghw3_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghw3_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_w3**4  &
        + ghw3_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_w3**4  &
        + ghw3_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_w3**4 

     vvcoupl_prime(4) = vvcoupl(4)   &
        + ghw4_prime * lambda_w4**4/(lambda_w4**2 + abs(sprod(j1,j2)))/(lambda_w4**2 + abs(sprod(j3,j4)))  &
        + ghw4_prime2* ( abs(sprod(j1,j2)) + abs(sprod(j3,j4)) )/lambda_w4**2  &
        + ghw4_prime3* ( abs(sprod(j1,j2)) - abs(sprod(j3,j4)) )/lambda_w4**2  &
        + ghw4_prime4* ( mhsq )/lambda_Q**2                                    &
        + ghw4_prime5* ( abs(sprod(j1,j2))**2 + abs(sprod(j3,j4))**2 )/lambda_w4**4  &
        + ghw4_prime6* ( abs(sprod(j1,j2))**2 - abs(sprod(j3,j4))**2 )/lambda_w4**4  &
        + ghw4_prime7* ( abs(sprod(j1,j2))    * abs(sprod(j3,j4))    )/lambda_w4**4 
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


  !-- QCD amplitudes squared below
  subroutine me2_ggggh(ggcoupl,j1,j2,j3,j4,za,zb,sprod,res)
    complex(dp), intent(in) :: ggcoupl(2:4)
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

    scalar = ggcoupl(2)
    pseudo = -ggcoupl(4) 

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

  subroutine me2_qbqggh(ggcoupl,j1,j2,j3,j4,za,zb,sprod,res)
    complex(dp), intent(in) :: ggcoupl(2:4)
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

    scalar = ggcoupl(2) 
    pseudo = -ggcoupl(4) 

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

   subroutine me2_qbqQBQ(ggcoupl,j1,j2,j3,j4,za,zb,sprod,res_diff,res_id)
    complex(dp), intent(in) :: ggcoupl(2:4)
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

    scalar = ggcoupl(2) 
    pseudo = -ggcoupl(4) 

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
  function scr(p1,p2) 
    real(dp), intent(in) :: p1(4), p2(4)
    real(dp) :: scr
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
  end function scr

  !- MCFM spinors
  subroutine spinoru(n,p,za,zb,s)
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
    
  end subroutine spinoru


end module modHiggsJJ
