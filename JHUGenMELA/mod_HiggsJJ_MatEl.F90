module modHiggsJJ
  implicit none
  private

  public :: EvalAmp_WBFH
  public :: EvalAmp_SBFH

  !-- general definitions, to be merged with Markus final structure
  integer, parameter  :: dp = selected_real_kind(15)

  include './variables.F90'

  real(8),  parameter :: pi =3.141592653589793238462643383279502884197d0  

  real(dp), parameter :: nf = 5.0_dp
  real(dp), parameter :: xw = sitW**2
  real(dp), parameter :: twosc = sqrt(4.0_dp*xw*(1.0_dp-xw))
  real(dp), parameter :: gs = sqrt(alphas*4.0_dp*pi)

  real(dp), parameter :: aR_lep = 2.0_dp * xw
  real(dp), parameter :: aL_lep = 2.0_dp * xw - 1.0_dp
  real(dp), parameter :: aR_neu = 0.0_dp
  real(dp), parameter :: aL_neu = 1.0_dp

  real(dp), parameter :: aR_QUp = -4.0_dp/3.0_dp * xw
  real(dp), parameter :: aL_QUp = -4.0_dp/3.0_dp * xw + 1.0_dp
  real(dp), parameter :: aR_QDn = 2.0_dp/3.0_dp * xw
  real(dp), parameter :: aL_QDn = 2.0_dp/3.0_dp * xw - 1.0_dp

  integer, target :: Up_  = 1
  integer, target :: Dn_  = 2
  integer, target :: Chm_ = 3
  integer, target :: Str_ = 4
  integer, target :: Top_ = 5
  integer, target :: Bot_ = 6
  integer, target :: ElP_ = 7
  integer, target :: MuP_ = 8
  integer, target :: TaP_ = 9
  integer, target :: Glu_ = 10
  integer, target :: Pho_ = 11
  integer, target :: Z0_  = 12
  integer, target :: Wp_  = 13
  integer, target :: NuE_ = 14
  integer, target :: NuM_ = 15
  integer, target :: NuT_ = 16

  integer,  target :: AUp_  = -1
  integer,  target :: ADn_  = -2
  integer,  target :: AChm_ = -3
  integer,  target :: AStr_ = -4
  integer,  target :: ATop_ = -5
  integer,  target :: ABot_ = -6
  integer,  target :: ElM_  = -7
  integer,  target :: MuM_  = -8
  integer,  target :: TaM_  = -9
  integer,  target :: Wm_   = -13
  integer,  target :: ANuE_ = -14
  integer,  target :: ANuM_ = -15
  integer,  target :: ANuT_ = -16

  integer, target :: pdfGlu_ = 0
 
  integer, target :: pdfDn_ = 1
  integer, target :: pdfUp_ = 2
  integer, target :: pdfStr_ = 3
  integer, target :: pdfChm_ = 4
  integer, target :: pdfBot_ = 5

  integer, target :: pdfADn_ = -1
  integer, target :: pdfAUp_ = -2
  integer, target :: pdfAStr_ = -3
  integer, target :: pdfAChm_ = -4
  integer, target :: pdfABot_ = -5

  real(dp), parameter :: tag1 = 1.0_dp
  real(dp), parameter :: tag2 = 1.0_dp
  real(dp), parameter :: tagbot = 1.0_dp

  real(dp), parameter :: xn = 3.0_dp
  real(dp), parameter :: Ca = 3.0_dp
  real(dp), parameter :: Cf = 4.0_dp/3.0_dp

  real(dp), parameter :: avegg = 1.0_dp/4.0_dp/64.0_dp
  real(dp), parameter :: aveqg = 1.0_dp/4.0_dp/24.0_dp
  real(dp), parameter :: aveqq = 1.0_dp/4.0_dp/9.0_dp
  real(dp), parameter :: SymmFac=1.0_dp/2.0_dp, SpinAvg=1.0_dp/4.0_dp

  real(dp), parameter :: mwsq = m_w**2
  real(dp), parameter :: mzsq = m_z**2

  real(dp), parameter :: fbGeV2=0.389379d12 * GeV**2

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one = 1.0_dp
  real(dp), parameter :: two = 2.0_dp
  complex(dp), parameter :: czero = (0.0_dp,0.0_dp)
  complex(dp), parameter :: ci = (0.0_dp,1.0_dp)

contains
  

  !-- SM: |g2| = alphas/(six*pi)
  !-- g3 not yet implemented
  subroutine EvalAmp_SBFH(pin,ggcoupl,res)
    real(dp), intent(in) :: pin(4,5)
    complex(dp), intent(in) :: ggcoupl(2:4)
    real(dp), intent(out) :: res(-5:5,-5:5)
    real(dp) :: shad, sqrts, x1, x2, etot, pztot
    real(dp) :: p(4,5), sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp, restmpid
    integer :: i, j, j3, j4, isymm

    res = zero

    !-- reconstruct the initial state momenta
    shad = two * scr(pin(:,1),pin(:,2))
    sqrts = sqrt(shad)
    etot = pin(1,3)+pin(1,4)+pin(1,5)
    pztot = pin(4,3)+pin(4,4)+pin(4,5)
    x1 = (etot+pztot)/sqrts
    x2 = (etot-pztot)/sqrts

    p(:,1) = sqrts/two * (/x1,zero,zero,x1/)
    p(:,2) = sqrts/two * (/x2,zero,zero,-x2/)
    p(:,3:5) = pin(:,3:5)

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    !-- gg -> gg
    call me2_ggggh(ggcoupl,1,2,3,4,za,zb,sprod,restmp)
    restmp = restmp * avegg * SymmFac
    res(0,0) = res(0,0) + restmp

    !-- gg -> qqb
    call me2_qbqggh(ggcoupl,4,3,1,2,za,zb,sprod,restmp)
    restmp = restmp * avegg
    res(0,0) = res(0,0) + restmp * nf

    !-- qqb -> gg
    call me2_qbqggh(ggcoupl,1,2,3,4,za,zb,sprod,restmp)
    restmp = restmp * aveqq * SymmFac
    do i = 1,5
       res(i,-i) = res(i,-i) + restmp
       res(-i,i) = res(-i,i) + restmp
    enddo

    !-- qqb -> rrb
    call me2_qbqQBQ(ggcoupl,1,2,4,3,za,zb,sprod,restmp,restmpid)
    restmp = restmp * aveqq
    do i = 1,5
       res(i,-i) = res(i,-i) + restmp * (nf-1.0_dp)
       res(-i,i) = res(-i,i) + restmp * (nf-1.0_dp)
    enddo

    !-- qq -> qq
    call me2_qbqQBQ(ggcoupl,1,3,2,4,za,zb,sprod,restmp,restmpid)
    restmp = restmpid * aveqq * SymmFac
    do i = 1,5
       res(i,i) = res(i,i) + restmp
       res(-i,-i) = res(-i,-i) + restmp
    enddo

    j3 = 3
    j4 = 4
    
    do isymm = 1,2

       !-- gq -> gq
       call me2_qbqggh(ggcoupl,2,j4,1,j3,za,zb,sprod,restmp)
       restmp = restmp * aveqg/ 2.0_dp
       do i = 1,5
          res(0,i) = res(0,i) + restmp
          res(0,-i) = res(0,-i) + restmp
       enddo
       call me2_qbqggh(ggcoupl,1,j4,2,j3,za,zb,sprod,restmp)
       restmp = restmp * aveqg/ 2.0_dp
       do i = 1,5
          res(i,0) = res(i,0) + restmp
          res(-i,0) = res(-i,0) + restmp
       enddo
       
       !-- qqb -> qqb
       call me2_qbqQBQ(ggcoupl,1,2,j4,j3,za,zb,sprod,restmp,restmpid)
       restmp = restmpid * aveqq/ 2.0_dp
       do i = 1,5
          res(i,-i) =res(i,-i) + restmp
       enddo
       call me2_qbqQBQ(ggcoupl,2,1,j4,j3,za,zb,sprod,restmp,restmpid)
       restmp = restmpid * aveqq/ 2.0_dp
       do i = 1,5
          res(-i,i) =res(-i,i) + restmp
       enddo

       !-- qrb -> qrb
       call me2_qbqQBQ(ggcoupl,1,j3,j4,2,za,zb,sprod,restmp,restmpid)
       restmp = restmp * aveqq/ 2.0_dp
       do i = 1,5
          do j = 1,5
             if (i.ne.j) then
	       res(i,-j) = res(i,-j) + restmp
	     endif
          enddo
       enddo
       call me2_qbqQBQ(ggcoupl,2,j3,j4,1,za,zb,sprod,restmp,restmpid)
       restmp = restmp * aveqq/ 2.0_dp
       do i = 1,5
          do j = 1,5
             if (i.ne.j) res(-j,i) = res(-j,i) + restmp
          enddo
       enddo

       !-- qr -> qr
       call me2_qbqQBQ(ggcoupl,1,j3,2,j4,za,zb,sprod,restmp,restmpid)
       restmp = restmp * aveqq/ 2.0_dp
       do i = 1,5
          do j = 1,5
             if (i.ne.j) then
                res(i,j) = res(i,j) + restmp
                res(-i,-j) = res(-i,-j) + restmp
             endif
          enddo
       enddo

       j3 = 4
       j4 = 3

    enddo

    return

  end subroutine EvalAmp_SBFH



  subroutine EvalAmp_WBFH(pin,vvcoupl,wwcoupl,res)
    real(dp), intent(in) :: pin(4,5)
    complex(dp), intent(in) :: vvcoupl(32),wwcoupl(32)
    real(dp), intent(out) :: res(-5:5,-5:5)
    real(dp) :: shad, sqrts, x1, x2, etot, pztot
    real(dp) :: p(4,5), sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: res1(-5:5,-5:5), res2(-5:5,-5:5)

    res = zero

    !-- reconstruct the initial state momenta
    shad = two * scr(pin(:,1),pin(:,2))
    sqrts = sqrt(shad)
    etot = pin(1,3)+pin(1,4)+pin(1,5)
    pztot = pin(4,3)+pin(4,4)+pin(4,5)
    x1 = (etot+pztot)/sqrts
    x2 = (etot-pztot)/sqrts

    p(:,1) = sqrts/two * (/x1,zero,zero,x1/)
    p(:,2) = sqrts/two * (/x2,zero,zero,-x2/)
    p(:,3:5) = pin(:,3:5)

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    call EvalAmp_WBFH_UnSymm(1,2,3,4,vvcoupl,wwcoupl,za,zb,sprod,res1)
    call EvalAmp_WBFH_UnSymm(1,2,4,3,vvcoupl,wwcoupl,za,zb,sprod,res2)
    
    res = (res1+res2)/2.0_dp

    return

  end subroutine EvalAmp_WBFH





  subroutine EvalAmp_WBFH_UnSymm(i1,i2,i3,i4,vvcoupl,wwcoupl,za,zb,sprod,res)
    integer, intent(in) :: i1,i2,i3,i4
    complex(dp), intent(in) :: vvcoupl(32),wwcoupl(32), za(4,4), zb(4,4)! wwcoupl is only used (for H->WW) if it includes non-zero values
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: res(-5:5,-5:5)
    complex(dp) :: amp_z(-1:1,-1:1), amp_z_b(-1:1,-1:1)
    complex(dp) :: amp_w(-1:1,-1:1)
    real(dp), parameter :: Lu = aL_QUp**2, Ru = aR_QUp**2
    real(dp), parameter :: Ld = aL_QDn**2, Rd = aR_QDn**2
    real(dp), parameter :: couplz = gwsq * xw/twosc**2
    real(dp), parameter :: couplw = gwsq/two
    real(dp) :: restmp
    integer :: i, j, j1, j2, iflip, pdfindex(2)

    res = zero

    !-- identical particles, up
    amp_z = A0_VV_4f(i4,i1,i3,i2,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(i3,i1,i4,i2,vvcoupl,za,zb,sprod,m_z,ga_z)

    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Lu**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Lu * Ru + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Lu * Ru + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Ru**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Lu**2  &
         + two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Ru**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfUp_,pdfUp_) = restmp
    res(pdfChm_,pdfChm_) = restmp

    !-- identical particles, down
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Ld**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Ld * Rd + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Ld * Rd + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Rd**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Ld**2  &
         + two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Rd**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfDn_,pdfDn_) = restmp
    res(pdfStr_,pdfStr_) = restmp
    res(pdfBot_,pdfBot_) = restmp * tagbot
    
    !-- identical particles, aup
    amp_z = A0_VV_4f(i1,i4,i2,i3,vvcoupl,za,zb,sprod,m_z,ga_z)
    amp_z_b = -A0_VV_4f(i1,i3,i2,i4,vvcoupl,za,zb,sprod,m_z,ga_z)

    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Lu**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Lu * Ru + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Lu * Ru + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Ru**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Lu**2  &
         + two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Ru**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfAUp_,pdfAUp_) = restmp
    res(pdfAChm_,pdfAChm_) = restmp

    !-- identical particles, adn
    restmp = ((abs(amp_z(-1,-1))**2+abs(amp_z_b(-1,-1))**2) * Ld**2 + &
         (abs(amp_z(-1,+1))**2+abs(amp_z_b(-1,+1))**2) * Ld * Rd + &
         (abs(amp_z(+1,-1))**2+abs(amp_z_b(+1,-1))**2) * Ld * Rd + &
         (abs(amp_z(+1,+1))**2+abs(amp_z_b(+1,+1))**2) * Rd**2) * xn**2

    restmp = restmp + (two * real(amp_z(-1,-1)*conjg(amp_z_b(-1,-1)),kind=dp) * Ld**2  &
         + two * real(amp_z(+1,+1)*conjg(amp_z_b(+1,+1)),kind=dp) * Rd**2) * xn

    restmp = restmp * SymmFac * aveqq * couplz**2

    res(pdfADn_,pdfADn_) = restmp
    res(pdfAStr_,pdfAStr_) = restmp
    res(pdfABot_,pdfABot_) = restmp * tagbot

    !-- W/Z interference
    j1 = i1
    j2 = i2
    do iflip = 1, 2 
       !-- ud -> ud
       amp_z = A0_VV_4f(i4,j2,i3,j1,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(i4,j1,i3,j2,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(i4,j1,i3,j2,wwcoupl,za,zb,sprod,m_w,ga_w) 
       endif
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
       amp_z = A0_VV_4f(j2,i4,j1,i3,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = -A0_VV_4f(j1,i4,j2,i3,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = -A0_VV_4f(j1,i4,j2,i3,wwcoupl,za,zb,sprod,m_w,ga_w) 
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

       j1 = i2
       j2 = i1
    enddo

    j1 = i1
    j2 = i2
    do iflip = 1,2

       !-- qqb processes

       !--uub -> uub // ddb
       amp_z = A0_VV_4f(i3,j1,j2,i4,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = A0_VV_4f(i3,j1,j2,i4,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = A0_VV_4f(i3,j1,j2,i4,wwcoupl,za,zb,sprod,m_w,ga_w) 
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
       amp_z = A0_VV_4f(i3,j1,i4,j2,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = A0_VV_4f(i3,j2,i4,j1,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = A0_VV_4f(i3,j2,i4,j1,wwcoupl,za,zb,sprod,m_w,ga_w) 
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
       amp_z = A0_VV_4f(j1,i3,j2,i4,vvcoupl,za,zb,sprod,m_z,ga_z)
       if( all(cdabs(wwcoupl(:)).lt.1d-15) ) then
          amp_w = A0_VV_4f(j1,i4,j2,i3,vvcoupl,za,zb,sprod,m_w,ga_w)
       else
          amp_w = A0_VV_4f(j1,i4,j2,i3,wwcoupl,za,zb,sprod,m_w,ga_w) 
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

       j1 = i2
       j2 = i1

    enddo
    
    return

  end subroutine EvalAmp_WBFH_UnSymm



  function flip(i,a1,a2)
    integer :: flip(2)
    integer :: i, a1, a2
    
    if (i .eq. 1) flip = (/a1,a2/)
    if (i .eq. 2) flip = (/a2,a1/)

    return

  end function flip



  !-------------------------------------------------------------------------
  !-- amplitudes below

  function A0_VV_4f(j1,j2,j3,j4,vvcoupl,za,zb,sprod,mv,ga_v)
    complex(dp) :: A0_VV_4f(-1:1,-1:1)
    integer :: j1,j2,j3,j4
    complex(dp) :: vvcoupl(32), za(4,4), zb(4,4)
    real(dp) :: mv, ga_v
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

    A0_VV_4f = A0_VV_4f/vev/iprop12/iprop34

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
          rt(j)=sqrt(p(2,j)+p(1,j))
          c23(j)=cmplx(p(4,j),-p(3,j),kind=dp)
          f(j)=(one,zero)
       else
       !-----negative energy case
          rt(j)=sqrt(-p(1,j)-p(2,j))
          c23(j)=cmplx(-p(4,j),p(3,j),kind=dp)
          f(j)=ci
       endif
    enddo

    do i=2,N
  
     do j=1,i-1
          s(i,j)=two*scr(p(:,i),p(:,j))
          za(i,j)=f(i)*f(j) &
               *(c23(i)*cmplx(rt(j)/rt(i),kind=dp)-c23(j)*cmplx(rt(i)/rt(j),kind=dp))
          
          if (abs(s(i,j)).lt.1d-5) then
             zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
          else
             zb(i,j)=-cmplx(s(i,j),kind=dp)/za(i,j)
          endif
          
          za(j,i)=-za(i,j)
          zb(j,i)=-zb(i,j)
          s(j,i)=s(i,j)
          
       enddo

    enddo

    return
    
  end subroutine spinoru


end module modHiggsJJ
