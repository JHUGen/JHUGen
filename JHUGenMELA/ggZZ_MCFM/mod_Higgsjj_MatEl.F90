module modHiggsjj
  implicit none
  include 'variables.F90'

  private

  public :: EvalAmp_gg_jjH

contains

  !--- normalization: g2 = 1 -> SM
  !--- ggcoupl(1) -> g2
  !--- ggcoupl(2) -> g3
  !--- ggcoupl(3) -> g4

  !----- p1 and p2 used to get hadronic s
  !----- unphysical kinematics to match Markus notation
  !----- 0 -> P(-p1)+P(-p2) + j(p3) + j(p4) + H(p5)
  subroutine EvalAmp_gg_jjH(pin,ggcoupl,me2)
    real(dp), intent(in) :: pin(4,5)
    complex(dp), intent(in) :: ggcoupl(1:3)
    real(dp), intent(out) :: me2(-5:5,-5:5)
    real(dp) :: p(4,5), xa, xb
    real(dp) :: shad,etot,pztot,sqrts
    real(dp) :: gg_hgg, gg_hqa, qg_hqg, gq_hqg
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    integer :: i
    real(dp) :: coupl_s

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

    call me2_ggggh(ggcoupl,1,2,3,4,za,zb,sprod,gg_hgg)
    call me2_qbqggh(ggcoupl,3,4,1,2,za,zb,sprod,gg_hqa)
    call me2_qbqggh(ggcoupl,1,3,2,4,za,zb,sprod,qg_hqg)
    call me2_qbqggh(ggcoupl,2,3,1,4,za,zb,sprod,gq_hqg)

    do i = -5,5
       me2(i,0) = (qg_hqg)*aveqg
       me2(0,i) = (gq_hqg)*aveqg
    enddo

    me2(0,0) = (gg_hgg/two + nf * gg_hqa)*avegg

    me2 = me2 * couplfac_ggh**2
        
  end subroutine EvalAmp_gg_jjH

  !--- 0 -> g(p1) g(p2) g(p3) g(p4) [H(p4)]
  subroutine me2_ggggh(ggcoupl,j1,j2,j3,j4,za,zb,sprod,me2gg)
    complex(dp), intent(in) :: ggcoupl(1:3)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: me2gg
    complex(dp) :: amp1234(-1:1,-1:1,-1:1), amp1324(-1:1,-1:1,-1:1)
    complex(dp) :: aPhi1234(1:2,-1:1,-1:1,-1:1), aPhi1324(1:2,-1:1,-1:1,-1:1)
    real(dp), parameter :: col1 = 8.0_dp * CA**3 * CF
    integer :: i2,i3,i4
    real(dp) :: q1q2
    complex(dp) :: cscalar, cpseudo

    me2gg = zero

    q1q2 = half * (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))

    !-- couplings
    cscalar = ggcoupl(1) + half * q1q2/Lambda**2 * ggcoupl(2)
    cpseudo = ggcoupl(3)

    aPhi1234 = A0phigggg_pxxx(j1,j2,j3,j4,za,zb,sprod)
    aPhi1324 = A0phigggg_pxxx(j1,j3,j2,j4,za,zb,sprod)

    amp1234(:,:,:) = cscalar * (aPhi1234(1,:,:,:) + aPhi1234(2,:,:,:)) + & !-- scalar
         cpseudo * (-ci) * (aPhi1234(1,:,:,:) - aPhi1234(2,:,:,:)) !-- pseudoscalar

    amp1324(:,:,:) = cscalar * (aPhi1324(1,:,:,:) + aPhi1324(2,:,:,:)) + & !-- scalar
         cpseudo * (-ci) * (aPhi1324(1,:,:,:) - aPhi1324(2,:,:,:)) !-- pseudoscalar

    do i2 = -1, 1, 2
       do i3 = -1, 1, 2
          do i4 = -1,1,2
             me2gg = me2gg + abs(amp1234(i2,i3,i4))**2
             me2gg = me2gg + abs(amp1324(i2,i3,i4))**2
             me2gg = me2gg + real(amp1234(i2,i3,i4)*conjg(amp1324(i3,i2,i4)),kind=dp)
          enddo
       enddo
    enddo

    !-- color factors and all
    me2gg = me2gg * col1

    !-- extra hels
    me2gg = me2gg * two

    return

  end subroutine me2_ggggh


  subroutine me2_qbqggh(ggcoupl,j1,j2,j3,j4,za,zb,sprod,me2q)
    complex(dp), intent(in) :: ggcoupl(1:3)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: me2q
    complex(dp) :: amp1234(-1:1,-1:1), amp1243(-1:1,-1:1)
    complex(dp) :: aPhi1234(1:2,-1:1,-1:1), aPhi1243(1:2,-1:1,-1:1)
    real(dp), parameter :: colf1 = four * xn * CF**2
    real(dp), parameter :: colf2 = four * xn * CF * (CF - CA/2.0_dp)
    integer :: i3,i4
    real(dp) :: q1q2
    complex(dp) :: cscalar, cpseudo

    me2q = zero

    q1q2 = half * (sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4))

    !-- couplings
    cscalar = ggcoupl(1) + half * q1q2/Lambda**2 * ggcoupl(2)
    cpseudo = ggcoupl(3)

    aPhi1234 = A0phiqbqgg_mpxx(j1,j2,j3,j4,za,zb,sprod)
    aPhi1243 = A0phiqbqgg_mpxx(j1,j2,j4,j3,za,zb,sprod)

    amp1234(:,:) = cscalar * (aPhi1234(1,:,:) + aPhi1234(2,:,:)) + & !-- scalar
         cpseudo * (-ci) * (aPhi1234(1,:,:) - aPhi1234(2,:,:)) !-- pseudoscalar

    amp1243(:,:) = cscalar * (aPhi1243(1,:,:) + aPhi1243(2,:,:)) + & !-- scalar
         cpseudo * (-ci) * (aPhi1243(1,:,:) - aPhi1243(2,:,:)) !-- pseudoscalar

    do i3 = -1,1,2
       do i4 = -1,1,2
          me2q = me2q + colf1 * abs(amp1234(i3,i4))**2
          me2q = me2q + colf1 * abs(amp1243(i4,i3))**2
          me2q = me2q + colf2 * two * real(amp1234(i3,i4)*conjg(amp1243(i4,i3)),kind=dp)
       enddo
    enddo

    !-- extra hels
    me2q = me2q * two

    return

  end subroutine me2_qbqggh


!---------------------------------------------------------------------------

!--- phi-amplitudes: 0->qb(p1)q(p2)g(p3)g(p4)
!--- iphi = 1 --> phi
!--- iphi = 2 --> phid
!--- assume qb^-,q^+
  function A0phiqbqgg_mpxx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phiqbqgg_mpxx(1:2,-1:1,-1:1)
    integer :: j1, j2, j3, j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: zab_1_ph_4, zab_3_ph_1, zab_3_ph_2, zab_3_ph_4
    complex(dp) :: zab_1_ph_2, zab_2_ph_4
    real(dp) :: s123, s412, qsq

    A0phiqbqgg_mpxx = czero

    s123 = sprod(j1,j2) + sprod(j1,j3) + sprod(j2,j3)
    s412 = sprod(j4,j1) + sprod(j4,j2) + sprod(j1,j2)

    zab_1_ph_4 = -za(j1,j2)*zb(j2,j4) - za(j1,j3)*zb(j3,j4)
    zab_3_ph_1 = -za(j3,j2)*zb(j2,j1) - za(j3,j4)*zb(j4,j1)
    zab_3_ph_2 = -za(j3,j1)*zb(j1,j2) - za(j3,j4)*zb(j4,j2)
    zab_3_ph_4 = -za(j3,j1)*zb(j1,j4) - za(j3,j2)*zb(j2,j4)
    zab_1_ph_2 = -za(j1,j3)*zb(j3,j2) - za(j1,j4)*zb(j4,j2)
    zab_2_ph_4 = -za(j2,j1)*zb(j1,j4) - za(j2,j3)*zb(j3,j4)

    qsq = s123 + sprod(j1,j4) + sprod(j2,j4) + sprod(j3,j4)

    A0phiqbqgg_mpxx(+1,+1,-1) = - za(j1,j4)**2 * za(j2,j4)/za(j1,j2)/za(j2,j3)/za(j3,j4) 
    A0phiqbqgg_mpxx(+2,+1,-1) = - zb(j2,j3)**2 * zb(j1,j3)/zb(j1,j2)/zb(j3,j4)/zb(j4,j1)

    A0phiqbqgg_mpxx(+1,-1,+1) = za(j1,j3)**3/za(j1,j2)/za(j3,j4)/za(j4,j1) 
    A0phiqbqgg_mpxx(+2,-1,+1) = zb(j2,j4)**3/zb(j1,j2)/zb(j2,j3)/zb(j3,j4)

    A0phiqbqgg_mpxx(+1,-1,-1) = qsq**2 * za(j1,j3)**3/s123/za(j1,j2)/zab_1_ph_4/zab_3_ph_4 &
         + zab_3_ph_1 * zab_3_ph_2**2/s412/zab_3_ph_4/zb(j2,j1)/zb(j4,j1) &
         - zab_1_ph_2**2/zab_1_ph_4/zb(j3,j2)/zb(j4,j3)

    A0phiqbqgg_mpxx(+2,+1,+1) = -zab_1_ph_2**2/za(j1,j4)/za(j3,j4)/zab_3_ph_2 + &
         zab_1_ph_4**2 * zab_2_ph_4/s123/za(j1,j2)/za(j2,j3)/zab_3_ph_4 &
         + qsq**2 * zb(j4,j2)**3/s412/zab_3_ph_2/zab_3_ph_4/zb(j2,j1)

    return

  end function A0phiqbqgg_mpxx


!--- phi-amplitudes: 0->g(p1)g(p2)g(p3)g(p4)
!--- iphi = 1 --> phi
!--- iphi = 2 --> phid
!--- assume h1=+1
  function A0phigggg_pxxx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phigggg_pxxx(1:2,-1:1,-1:1,-1:1)
    integer :: j1, j2, j3, j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)

    A0phigggg_pxxx = czero

    A0phigggg_pxxx(+2,+1,+1,+1) = A0phiggggmmmm(j1,j2,j3,j4,zb,za,sprod)

    A0phigggg_pxxx(+2,+1,+1,-1) = A0phiggggpmmm(j4,j1,j2,j3,zb,za,sprod)

    A0phigggg_pxxx(+2,+1,-1,+1) = A0phiggggpmmm(j3,j4,j1,j2,zb,za,sprod)

    A0phigggg_pxxx(+1,+1,-1,-1) = A0phiggggmmpp(j3,j4,j1,j2,za,zb,sprod)
    A0phigggg_pxxx(+2,+1,-1,-1) = A0phiggggmmpp(j1,j2,j3,j4,zb,za,sprod)

    A0phigggg_pxxx(+2,-1,+1,+1) = A0phiggggpmmm(j2,j3,j4,j1,zb,za,sprod)

    A0phigggg_pxxx(+1,-1,+1,-1) = A0phiggggmpmp(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_pxxx(+2,-1,+1,-1) = A0phiggggmpmp(j1,j2,j3,j4,zb,za,sprod)

    A0phigggg_pxxx(+1,-1,-1,+1) = A0phiggggmmpp(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_pxxx(+2,-1,-1,+1) = A0phiggggmmpp(j4,j1,j2,j3,zb,za,sprod)

    A0phigggg_pxxx(+1,-1,-1,-1) = A0phiggggpmmm(j1,j2,j3,j4,za,zb,sprod)

    return

  end function A0phigggg_pxxx


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

  include 'include_functions.f90'

end module modHiggsjj
