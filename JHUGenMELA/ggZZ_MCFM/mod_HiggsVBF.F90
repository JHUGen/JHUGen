module ModHiggsVBF
implicit none



  real(8),private :: p(4,5)
  private
  public :: EvalAmp_vbfH

  integer, parameter,private  :: dp = selected_real_kind(15)

 contains

  !--- current normalization: g1 = 1 --> SM

  !--- vvcoupl(1) -> g1
  !--- vvcoupl(2) -> g2
  !--- vvcoupl(3) -> g3
  !--- vvcoupl(4) -> g4

  !----- p1 and p2 used to get hadronic s
  !----- unphysical kinematics to match Markus notation
  !----- 0 -> P(-p1)+P(-p2) + j(p3) + j(p4) + H(p5)
  subroutine EvalAmp_vbfH(pin,vvcoupl,me2)
  implicit none
    real(dp), intent(in) :: pin(4,5)
    complex(dp), intent(in) :: vvcoupl(1:4)
    real(dp), intent(out) :: me2(-5:5,-5:5)
    real(dp) :: xa, xb
    real(dp) :: shad,etot,pztot,sqrts
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    real(dp) :: me2_WW(-5:5,-5:5), me2_ZZ(-5:5,-5:5)
    include 'includeVars.F90'


    me2 = zero

    !-- get hadronic energy
    shad = two * scr(pin(:,1),pin(:,2))
    sqrts = sqrt(shad)

    !-- get boost:
    !-- Etot = (xa+xb)sqrts/two
    !-- pztot = (xa-xb)sqrts/two
    etot = pin(1,3)+pin(1,4)+pin(1,5)
    pztot = pin(4,3)+pin(4,4)+pin(4,5)
    xa = (etot+pztot)/sqrts
    xb = (etot-pztot)/sqrts

    !-- all outgoing kinematics
    p(:,1) = -sqrts/two * (/xa,zero,zero,xa/)
    p(:,2) = -sqrts/two * (/xb,zero,zero,-xb/)
    p(:,3) = pin(:,3)
    p(:,4) = pin(:,4)
    p(:,5) = pin(:,5)

    call spinoru(4,p,za,zb,sprod)

    call VBF_WW(vvcoupl,za,zb,sprod,me2_WW)
    call VBF_ZZ(vvcoupl,za,zb,sprod,me2_ZZ)

    me2 = me2_WW + me2_ZZ
    

    return

  end subroutine EvalAmp_vbfH

  subroutine VBF_WW(vvcoupl,za,zb,sprod,me2)
  implicit none
    complex(dp), intent(in) :: vvcoupl(1:4)
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: me2(-5:5,-5:5)
    real(dp) :: ud_udh_WW, uub_ddbh_WW
    integer :: j,k
    include 'includeVars.F90'

    me2 = zero


    call me2_tree_qqqqH_WW(vvcoupl,1,3,2,4,za,zb,sprod,ud_udh_ww)
    call me2_tree_qqqqH_WW(vvcoupl,1,3,4,2,za,zb,sprod,uub_ddbh_ww)
    ud_udh_ww = ud_udh_ww*aveqq * couplfac_ww**2
    uub_ddbh_ww = uub_ddbh_ww*aveqq * couplfac_ww**2

    do j=-(nf-1),nf-1
       do k=-(nf-1),nf-1
          if     ((j .gt. 0) .and. (k .lt. 0)) then
             if (pn(j) .eq. -pn(k)) me2(j,k)=uub_ddbh_ww             
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
             if (pn(j) .eq. -pn(k)) me2(j,k)=uub_ddbh_ww
          elseif ((j .gt. 0) .and. (k .gt. 0)) then
             if (pn(j)+pn(k) .eq. +3) me2(j,k)=ud_udh_ww
          elseif ((j .lt. 0) .and. (k .lt. 0)) then
             if (pn(j)+pn(k) .eq. -3) me2(j,k)=ud_udh_ww
          endif
       enddo
    enddo
        
  end subroutine VBF_WW

  subroutine VBF_ZZ(vvcoupl,za,zb,sprod,me2)
  implicit none
    complex(dp), intent(in) :: vvcoupl(1:4)
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: me2(-5:5,-5:5)
    real(dp) :: ud_udh_LL,ud_udh_LR, udb_udbh_LL, udb_udbh_LR
    integer :: j,k
    include 'includeVars.F90'

    me2 = zero

    call me2_tree_qqqqH_ZZ(vvcoupl,1,3,2,4,za,zb,sprod,ud_udh_LL,ud_udh_LR)
    call me2_tree_qqqqH_ZZ(vvcoupl,1,3,4,2,za,zb,sprod,udb_udbh_LL,udb_udbh_LR)

    ud_udh_LL = ud_udh_LL * aveqq * couplfac_zz**2
    ud_udh_LR = ud_udh_LR * aveqq * couplfac_zz**2

    udb_udbh_LL = udb_udbh_LL * aveqq * couplfac_zz**2
    udb_udbh_LR = udb_udbh_LR * aveqq * couplfac_zz**2

    do j=-nf,nf
       do k=-nf,nf
          if     ((j .gt. 0) .and. (k .lt. 0)) then
             me2(j,k)= +udb_udbh_LL*((L(+j)*L(-k))**2+(R(+j)*R(-k))**2) &
                  +udb_udbh_LR*((L(+j)*R(-k))**2+(R(+j)*L(-k))**2)
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
             me2(j,k)= +udb_udbh_LL*((L(-j)*L(k))**2+(R(-j)*R(k))**2) &
                  +udb_udbh_LR*((L(-j)*R(k))**2+(R(-j)*L(k))**2)
          elseif ((j .gt. 0) .and. (k .gt. 0)) then
             me2(j,k)= +ud_udh_LL*((L(+j)*L(+k))**2+(R(+j)*R(+k))**2) &
                  +ud_udh_LR*((L(+j)*R(+k))**2+(R(+j)*L(+k))**2)
          elseif ((j .lt. 0) .and. (k .lt. 0)) then
             me2(j,k)= +ud_udh_LL*((L(-j)*L(-k))**2+(R(-j)*R(-k))**2) &
                  +ud_udh_LR*((L(-j)*R(-k))**2+(R(-j)*L(-k))**2)
          endif
       enddo
    enddo

    return
        
  end subroutine VBF_ZZ


  !-----------------------------------------------------------------

  ! Notation: 0 -> qbar(p1) q(p2) qbar(p3) q(p4) H
  ! Factored out: (gwsq**2/two * vev * sw**2/cw**2)**2
  subroutine me2_tree_qqqqH_ZZ(vvcoupl,j1,j2,j3,j4,za,zb,sprod,me2_LL, me2_LR)
  implicit none
    complex(dp), intent(in) :: vvcoupl(1:4)
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    integer, intent(in) :: j1, j2, j3, j4
    real(dp), intent(out) :: me2_LL, me2_LR
    real(dp) :: prefac, rme2
    complex(dp) :: amp(-1:1,-1:1,-1:1,-1:1)
    real(dp) :: q1sq, q2sq, prop1, prop2
    real(dp), parameter :: colf = (3d0)**2
    complex(dp) :: aacoupl(1:3)
    real(dp) :: q1q2, pH(4), mhsq
    include 'includeVars.F90'

    me2_LL = zero
    me2_LR = zero

    rme2 = zero

    q1sq = sprod(j1,j2)
    q2sq = sprod(j3,j4)

    prop1 = one/(q1sq-mzsq)
    prop2 = one/(q2sq-mzsq)

    prefac = prop1 * prop2
    prefac = prefac**2

    pH(:) = -p(:,1)-p(:,2)-p(:,3)-p(:,4)
    mhsq = scr(pH,pH)

    !-- couplings
    q1q2 = half * (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))
    
    aacoupl(1) = vvcoupl(1) * mzsq/mhsq + two*q1q2/mhsq * vvcoupl(2) + q1q2**2/Lambda**2/mhsq * vvcoupl(3)
    aacoupl(2) = -two * vvcoupl(2) - q1q2/Lambda**2 * vvcoupl(3)
    aacoupl(3) = -two * vvcoupl(4)

    amp = A0Hqqqq(mhsq,aacoupl,j1,j2,j3,j4,za,zb,sprod)

    me2_LL = amp(-1,+1,-1,+1)*conjg(amp(-1,+1,-1,+1))
    me2_LR = amp(-1,+1,+1,-1)*conjg(amp(-1,+1,+1,-1))

    me2_LL = me2_LL * prefac * colf
    me2_LR = me2_LR * prefac * colf

    return

  end subroutine me2_tree_qqqqH_ZZ


  ! Notation: 0 -> qbar(p1) q(p2) qbar(p3) q(p4) H
  ! Factored out: (gwsq**2/two * vev)**2
  subroutine me2_tree_qqqqH_WW(vvcoupl,j1,j2,j3,j4,za,zb,sprod,me2)
  implicit none
    complex(dp), intent(in) :: vvcoupl(1:4)
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    integer, intent(in) :: j1, j2, j3, j4
    real(dp), intent(out) :: me2
    real(dp) :: prefac, rme2
    complex(dp) :: amp(-1:1,-1:1,-1:1,-1:1)
    real(dp) :: q1sq, q2sq, prop1, prop2
    real(dp), parameter :: colf = (3d0)**2
    complex(dp) :: aacoupl(1:3)
    real(dp) :: q1q2, pH(4), mhsq
    include 'includeVars.F90'

    me2 = zero
    rme2 = zero

    q1sq = sprod(j1,j2)
    q2sq = sprod(j3,j4)

    prop1 = one/(q1sq-mwsq)
    prop2 = one/(q2sq-mwsq)

    prefac = prop1 * prop2
    prefac = prefac**2

    pH = -p(:,1)-p(:,2)-p(:,3)-p(:,4)
    mhsq = scr(pH,pH)

    !-- couplings
    q1q2 = half * (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))

    aacoupl(1) = vvcoupl(1) * mwsq/mhsq + two*q1q2/mhsq * vvcoupl(2) + q1q2**2/Lambda**2/mhsq * vvcoupl(3)
    aacoupl(2) = -two * vvcoupl(2) - q1q2/Lambda**2 * vvcoupl(3)
    aacoupl(3) = -two * vvcoupl(4)
    
    amp = A0Hqqqq(mhsq,aacoupl,j1,j2,j3,j4,za,zb,sprod)

!------ The W-boson contribution / only one helicity plays a role
    rme2 = rme2 + amp(-1,+1,-1,+1)*conjg(amp(-1,+1,-1,+1))

    me2 = rme2 * prefac * colf

    return

  end subroutine me2_tree_qqqqH_WW


  !----------------------------------------------------------------------------------------

    !-- index: helicity (h1,h3) = chirality of the vertex
  !-- amplitude for H -> q qb Q QB
  function A0Hqqqq(mhsq,aacoupl,j1,j2,j3,j4,za,zb,sprod)
  implicit none
    real(dp) :: mhsq
    complex(8) :: aacoupl(1:3)
    integer :: j1,j2,j3,j4
    complex(8) :: za(4,4), zb(4,4)
    real(8) :: sprod(4,4)
!     complex(dp) :: A0Hqqqq(-1:1,-1:1) ! Fabrizio
    complex(dp) :: A0Hqqqq(-1:1,-1:1,-1:1,-1:1) ! Markus
    integer :: i,j
    real(8) :: q1q2
    complex(8) :: aa1, aa2, aa3
    complex(8) :: zab2 ! -- <i|p1+p2|j]         ! Fabrizio
!     complex(8) :: zab2(4,4,4,4) ! -- <i|p1+p2|j]       ! Markus
    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
    include 'includeVars.F90'



    A0Hqqqq = czero

    q1q2 = half * (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))

    aa1 = aacoupl(1) * mhsq
    aa2 = aacoupl(2)
    aa3 = aacoupl(3)


! !   FABRIZIO
!     A0Hqqqq(+1,+1) = (aa1 -ci * aa3 * q1q2) * (zb(j1,j3)*za(j4,j2) * two) &
!          + (aa2 + ci * aa3) * zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3) &
!          -ci * aa3 * (two * za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j1))
! 
!     A0Hqqqq(+1,-1) = (aa1 -ci * aa3 * q1q2) * (zb(j1,j4)*za(j3,j2) * two) &
!          + (aa2 + ci * aa3) * zab2(j2,j3,j4,j1)*zab2(j3,j1,j2,j4) &
!          -ci * aa3 * (two * za(j2,j1)*zb(j1,j4)*za(j3,j4)*zb(j4,j1))
! 
!     A0Hqqqq(-1,+1) = (aa1 -ci * aa3 * q1q2) * (za(j1,j4)*zb(j3,j2) * two) &
!          + (aa2 + ci * aa3) * zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3) &
!          -ci * aa3 * (two * za(j1,j2)*zb(j2,j3)*za(j4,j3)*zb(j3,j2))
! 
!     A0Hqqqq(-1,-1) = (aa1 -ci * aa3 * q1q2) * (za(j1,j3)*zb(j4,j2) * two) &
!          + (aa2 + ci * aa3) * zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4) &
!          -ci * aa3 * (two * za(j1,j2)*zb(j2,j4)*za(j3,j4)*zb(j4,j2))
   

!   MARKUS
    A0Hqqqq(+1,-1,+1,-1) = (aa1 -ci * aa3 * q1q2) * (zb(j1,j3)*za(j4,j2) * two) &
         + (aa2 + ci * aa3) * zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3) &
         -ci * aa3 * (two * za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j1))
    A0Hqqqq(+1,-1,-1,+1) = (aa1 -ci * aa3 * q1q2) * (zb(j1,j4)*za(j3,j2) * two) &
         + (aa2 + ci * aa3) * zab2(j2,j3,j4,j1)*zab2(j3,j1,j2,j4) &
         -ci * aa3 * (two * za(j2,j1)*zb(j1,j4)*za(j3,j4)*zb(j4,j1))
    A0Hqqqq(-1,+1,+1,-1) = (aa1 -ci * aa3 * q1q2) * (za(j1,j4)*zb(j3,j2) * two) &
         + (aa2 + ci * aa3) * zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3) &
         -ci * aa3 * (two * za(j1,j2)*zb(j2,j3)*za(j4,j3)*zb(j3,j2))
    A0Hqqqq(-1,+1,-1,+1) = (aa1 -ci * aa3 * q1q2) * (za(j1,j3)*zb(j4,j2) * two) &
         + (aa2 + ci * aa3) * zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4) &
         -ci * aa3 * (two * za(j1,j2)*zb(j2,j4)*za(j3,j4)*zb(j4,j2))
    

    return
    
  end function A0Hqqqq

  


  include 'includeFunctions.F90'
  
end module
