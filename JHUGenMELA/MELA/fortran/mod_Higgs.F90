      MODULE modHiggs
      use ModParameters
      use ModMisc
      implicit none
      private


!----- notation for subroutines
      public :: EvalAmp_gg_H_VV
      public :: EvalAmp_H_VV
      public :: EvalAmp_H_FF, EvalAmp_H_TT_decay
      public :: calcHelAmp,calcHelAmp2

      CONTAINS


!----- a subroutinefor gg -> H -> (Z+gamma*)(Z+gamma*)/WW/gammagamma
!----- all outgoing convention and the following momentum assignment
!-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      SUBROUTINE EvalAmp_gg_H_VV(p,MY_IDUP,res)
      implicit none
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A_VV(1:18), A0_VV(1:2)
      integer :: i1,i2,i3,i4,VVMode,VVmode_swap
      real(dp) :: prefactor!,res2
      real(dp) :: intcolfac
      integer :: ordering(1:4),ordering_swap(1:4)
      logical :: doInterference

      if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
         intcolfac=1.0_dp/3.0_dp
      else
         intcolfac=1.0_dp
      endif

      call getDecay_VVMode_Ordering(MY_IDUP(6:9),VVMode,ordering,VVmode_swap,ordering_swap)

! Global normalization
      if( VVMode.eq.ZZMode ) then!  Z decay
         prefactor = 8d0*overallCouplVffsq**2
      elseif( VVMode.eq.WWMode ) then !  W decay
         prefactor = 8d0*overallCouplVffsq**2
      elseif( VVMode.eq.ZgMode ) then !  Z+photon "decay"
         prefactor = 8d0*overallCouplVffsq ! Only single powers
      elseif( VVMode.eq.ggMode ) then !  photon "decay"
         prefactor = 8d0
      else
         prefactor = 0d0
      endif
      prefactor = prefactor * (alphas/(3.0_dp*pi*vev))**2


! ! MADGRAPH CHECK
! res=0d0
! if (MY_IDUP(6).ne.MY_IDUP(8) ) return
! if (MY_IDUP(7).ne.MY_IDUP(9) ) return
! if (MY_IDUP(6).eq.MY_IDUP(8) ) return
! if (MY_IDUP(7).eq.MY_IDUP(9) ) return
! print *, "MY COUPL",al1*dsqrt(overallCouplVffsq),ar1*dsqrt(overallCouplVffsq)
! ar1=-0.158480099490745d0! this is MadGraphs gzl(R)  , the MG couplings differ by a global minus sign, relative differences are because of different input parameters
! al1=+0.210150647402957d0! this is MadGraphs gzl(L)
! al2=al1
! ar2=ar1
! print *, "MG COUPL",al1,ar1
! pause

          res = zero
          A_VV(:) = 0d0
          doInterference = includeInterference .and. (         &
            ((VVMode.eq.ZZMode) .and. (VVMode_swap.eq.ZZMode)) &
            )
          if ( includeVprime .and. .not.(VVMode.eq.ZZMode .or. VVMode.eq.ZgMode .or. VVMode.eq.WWMode) ) then
              call Error("Contact terms only for WW, ZZ or Zg!")
          endif
          do i1=1,2;  do i2=1,2;  do i3=1,2;  do i4=1,2!  sum over helicities
                  call calcHelAmp(ordering,VVMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(1))
                  if( VVMode.eq.ZZMode ) then
                      if( includeGammaStar ) then
                         call calcHelAmp(ordering,ZgsMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(3))
                         call calcHelAmp(ordering,gsZMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(5))
                         call calcHelAmp(ordering,gsgsMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(7))
                      endif
                      if( includeVprime ) then
                          call calcHelAmp(ordering,ZZpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(9))
                          call calcHelAmp(ordering,ZpZMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(11))
                          call calcHelAmp(ordering,ZpZpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(13))
                      endif
                      if( includeGammaStar .and. includeVprime ) then
                          call calcHelAmp(ordering,gsZpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(15))
                          call calcHelAmp(ordering,ZpgsMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(17))
                      endif
                  elseif( VVMode.eq.ZgMode ) then
                      if(includeGammaStar) then
                          call calcHelAmp(ordering,gsgMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(3))
                      endif
                      if( includeVprime ) then
                          call calcHelAmp(ordering,ZpgMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(5))
                      endif
                  elseif( VVMode.eq.WWMode .and. includeVprime ) then
                      call calcHelAmp(ordering,WWpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(9))
                      call calcHelAmp(ordering,WpWMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(11))
                      call calcHelAmp(ordering,WpWpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(13))
                  endif

                  if( doInterference ) then
                      call calcHelAmp(ordering_swap,VVMode_swap,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(2))
                      if( includeGammaStar ) then
                          call calcHelAmp(ordering_swap,ZgsMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(4))
                          call calcHelAmp(ordering_swap,gsZMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(6))
                          call calcHelAmp(ordering_swap,gsgsMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(8))
                      endif
                      if( includeVprime ) then
                          call calcHelAmp(ordering_swap,ZZpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(10))
                          call calcHelAmp(ordering_swap,ZpZMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(12))
                          call calcHelAmp(ordering_swap,ZpZpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(14))
                      endif
                      if( includeGammaStar .and. includeVprime ) then
                          call calcHelAmp(ordering_swap,gsZpMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(16))
                          call calcHelAmp(ordering_swap,ZpgsMode,MY_IDUP,p(1:4,1:6),i1,i2,i3,i4,A_VV(18))
                      endif
                  endif

                  A0_VV(1) = A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7)+A_VV(9)+A_VV(11)+A_VV(13)+A_VV(15)+A_VV(17) ! 3456 pieces
                  A0_VV(2) = A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8)+A_VV(10)+A_VV(12)+A_VV(14)+A_VV(16)+A_VV(18) ! 5436 pieces
                  res = res + dreal(A0_VV(1)*dconjg(A0_VV(1)))
                  res = res + dreal(A0_VV(2)*dconjg(A0_VV(2)))
                  if( doInterference .and. (i3.eq.i4) ) then! interfere the 3456 with 5436 pieces
                      res = res - 2d0*intcolfac*dreal(  A0_VV(1)*dconjg(A0_VV(2))  ) ! minus from Fermi statistics
                  endif
          enddo;  enddo;  enddo;  enddo


! MADGRAPH CHECK
! call coupsm(0)
! if( (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!       call SH_EMEPEMEP((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,res2)
! else
!       call SH_EMEPEMEP_NOINT((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,res2)
!       call SH_TAMTAPTAMTAP_NOINT((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,res2)
! endif
! res2=res2* cdabs( (0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso) *  dconjg((0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso)) )
! res2=res2/100d0**2/100d0**2
! pause
! res=res2; RETURN

! res= res*prefactor
! print *, "checker 1",res
! print *, "checker 2",res2
! print *, "checker 1/2",res/res2
! pause


          res = res*prefactor
          if( (VVMode.eq.ZZMode) .and. doInterference ) res = res * SymmFac

         !print *,"VVmode:",VVmode
         !print *,"VVmode_swap:",VVmode_swap
         !print *,"ids:",MY_IDUP
         !print *,"res:",res
         !pause

      RETURN
      END SUBROUTINE

   subroutine calcHelAmp(ordering,VVMode,MY_IDUP,p,i1,i2,i3,i4,A)
      implicit none
      integer :: ordering(1:4),VVMode,i1,i2,i3,i4,l1,l2,l3,l4,MY_IDUP(6:9)
      real(dp) :: p(1:4,1:6)
      complex(dp) :: propG
      real(dp) :: pin(4,4)
      complex(dp) :: A(1:1), sp(4,4)

      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      !print *,"l1=",l1
      !print *,"l2=",l2
      !print *,"l3=",l3
      !print *,"l4=",l4
      !print *,"id(6)=",convertLHE(MY_IDUP(l1))
      !print *,"id(7)=",convertLHE(MY_IDUP(l2))
      !print *,"id(8)=",convertLHE(MY_IDUP(l3))
      !print *,"id(9)=",convertLHE(MY_IDUP(l4))
      !print *,"p(6)=",(p(:,l1))
      !print *,"p(7)=",(p(:,l2))
      !print *,"p(8)=",(p(:,l3))
      !print *,"p(9)=",(p(:,l4))

      propG = one/dcmplx(2d0*scr(p(:,1),p(:,2)) - M_Reso**2,M_Reso*Ga_Reso)

      pin(1,:) = p(:,1)
      pin(2,:) = p(:,2)
      sp(1,:) = pol_mless2(dcmplx(p(:,1)),-3+2*i1,'in')  ! gluon
      sp(2,:) = pol_mless2(dcmplx(p(:,2)),-3+2*i2,'in')  ! gluon
      !sp(1,1:4)=pin(1,1:4);print *, "this checks IS gauge invariance"
      !sp(2,1:4)=pin(2,1:4);print *, "this checks IS gauge invariance"
      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call ggHZZampl(VVMode,pin,sp,A(1))
      A(1) = A(1) * propG
   end subroutine

      SUBROUTINE ggHZZampl(VVMode,p,sp,res)
      implicit none
      integer, intent(in) :: VVMode
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4)
      complex(dp) :: xxx1,xxx2,xxx3,yyy1,yyy2,yyy3,yyy4
      complex(dp) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      real(dp) :: q34
      real(dp) :: q_q, q3_q3, q4_q4


      res = 0d0
      q1 = dcmplx(p(1,:),0d0)
      q2 = dcmplx(p(2,:),0d0)
      q3 = dcmplx(p(3,:),0d0)
      q4 = dcmplx(p(4,:),0d0)

      e1 = sp(1,:)
      e2 = sp(2,:)
      e3 = sp(3,:)
      e4 = sp(4,:)

      q = -q1-q2

      q_q =sc(q,q)
      q3_q3 = sc(q3,q3)
      q4_q4 = sc(q4,q4)

      e1_e2 = sc(e1,e2)
      e1_e3 = sc(e1,e3)
      e1_e4 = sc(e1,e4)
      e2_e3 = sc(e2,e3)
      e2_e4 = sc(e2,e4)
      e3_e4 = sc(e3,e4)

      q1_e3 = sc(q1,e3)
      q1_e4 = sc(q1,e4)
      q2_e3 = sc(q2,e3)
      q2_e4 = sc(q2,e4)
      e1_q3 = sc(e1,q3)
      e1_q4 = sc(e1,q4)
      e2_q3 = sc(e2,q3)
      e2_q4 = sc(e2,q4)
      e3_q4 = sc(e3,q4)
      e4_q3 = sc(e4,q3)


      if (q_q.lt.-0.1d0 .or. q3_q3.lt.-0.1d0 .or. q4_q4.lt.-0.1d0) return  ! if negative invariant masses return zero

!---- data that defines couplings
      if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
         ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.gsZMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q3_q3,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q3_q3,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q3_q3,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q3_q3,q_q)
      elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.ZgsMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q_q)
      elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then
         ghz1_dyn = czero
         ghz2_dyn = HVVSpinZeroDynamicCoupling(9,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(10,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(11,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.ZZpMode) .or. (VVMode.eq.WWpMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(12,q3_q3,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(13,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(14,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(15,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.ZpZMode) .or. (VVMode.eq.WpWMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(12,q4_q4,q3_q3,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(13,q4_q4,q3_q3,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(14,q4_q4,q3_q3,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(15,q4_q4,q3_q3,q_q)
      elseif( (VVMode.eq.ZpZpMode) .or. (VVMode.eq.WpWpMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(16,q3_q3,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(17,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(18,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(19,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.gsZpMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(20,0d0,q3_q3,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(21,0d0,q3_q3,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(22,0d0,q3_q3,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(23,0d0,q3_q3,q_q)
      elseif( (VVMode.eq.ZpgMode) .OR. (VVMode.eq.ZpgsMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(20,0d0,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(21,0d0,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(22,0d0,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(23,0d0,q4_q4,q_q)
      else
         ghz1_dyn = czero
         ghz2_dyn = czero
         ghz3_dyn = czero
         ghz4_dyn = czero
         print *,"VVMode",VVMode,"not implemented"
      endif

      if( .not. generate_as ) then
         xxx1 = ghg2+ghg3/4d0/Lambda**2*q_q
         xxx3 = -2d0*ghg4
         yyy1 = ghz1_dyn*M_V**2/q_q &  ! in this line M_V is indeed correct, not a misprint
              + ghz2_dyn*(q_q-q3_q3-q4_q4)/q_q &
              + ghz3_dyn/Lambda**2*(q_q-q3_q3-q4_q4)*(q_q-q4_q4-q3_q3)/4d0/q_q
         yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(q_q-q3_q3-q4_q4)
         yyy3 = -2d0*ghz4_dyn
      else
         xxx1 = ahg1
         xxx3 = ahg3
         if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
            yyy1 = ahz1
            yyy2 = ahz2
            yyy3 = ahz3
         elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then! decay (gamma-gamma) OR (gamma*-gamma*) OR (gamma*-gamma)
            yyy1 = ahz1
            yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
            yyy3 = ahz3
         elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then! decay (Z-photon) OR (gamma*-Z) OR (Z-gamma*)
            yyy1 = ahz1
            yyy2 = -2*ahz1*q_q/(q_q-q3_q3)
            yyy3 = ahz3
         else
            yyy1=czero
            yyy2=czero
            yyy3=czero
            print *,"VVMode",VVMode,"not implemented in generate_as"
         endif
      endif

      res= e1_e2*e3_e4*q_q**2*yyy1*xxx1                  &
         + e1_e2*e3_q4*e4_q3*q_q*yyy2*xxx1               &
         + et1(e1,e2,q1,q2)*e3_e4*q_q*yyy1*xxx3          &
         + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3        &
         + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3   &
         + et1(e3,e4,q3,q4)*e1_e2*q_q*yyy3*xxx1

  END SUBROUTINE ggHZZampl


!----- a subroutinefor H -> (Z+gamma*)(Z+gamma*)/WW/gammagamma
!----- all outgoing convention and the following momentum assignment
!-----  0 -> Higgs(p1) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_H_VV(p,MY_IDUP,res)
      implicit none
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A_VV(1:18),A0_VV(1:2)
      integer :: i3,i4,VVMode,VVmode_swap
      real(dp) :: prefactor
      real(dp) :: intcolfac
      integer :: ordering(1:4),ordering_swap(1:4)
      logical :: doInterference

      if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
         intcolfac=1.0_dp/3.0_dp
      else
         intcolfac=1.0_dp
      endif

      call getDecay_VVMode_Ordering(MY_IDUP(6:9),VVMode,ordering,VVmode_swap,ordering_swap)

! Global normalization
         if( VVMode.eq.ZZMode ) then!  Z decay
            prefactor = overallCouplVffsq**2
         elseif( VVMode.eq.WWMode ) then !  W decay
            prefactor = overallCouplVffsq**2
         elseif( VVMode.eq.ZgMode ) then !  Z+photon "decay"
            prefactor = overallCouplVffsq ! Only single powers
         elseif( VVMode.eq.ggMode ) then !  photon "decay"
            prefactor = 1d0
         else
            prefactor = 0d0
         endif

! ! MADGRAPH CHECK
! res=0d0
! if (MY_IDUP(6).ne.MY_IDUP(8) ) return
! if (MY_IDUP(7).ne.MY_IDUP(9) ) return
! if (MY_IDUP(6).eq.MY_IDUP(8) ) return
! if (MY_IDUP(7).eq.MY_IDUP(9) ) return
! print *, "MY COUPL",al1*dsqrt(overallCouplVffsq),ar1*dsqrt(overallCouplVffsq)
! ar1=-0.158480099490745d0! this is MadGraphs gzl(R)  , the MG couplings differ by a global minus sign, relative differences are because of different input parameters
! al1=+0.210150647402957d0! this is MadGraphs gzl(L)
! al2=al1
! ar2=ar1
! print *, "MG COUPL",al1,ar1
! pause


          res = zero
          A_VV(:) = 0d0
          doInterference = includeInterference .and. (         &
            ((VVMode.eq.ZZMode) .and. (VVMode_swap.eq.ZZMode)) &
            )
          if ( includeVprime .and. .not.(VVMode.eq.ZZMode .or. VVMode.eq.ZgMode .or. VVMode.eq.WWMode) ) then
              call Error("Contact terms only for WW, ZZ or Zg!")
          endif
          do i3=1,2;  do i4=1,2!  sum over helicities
                  call calcHelAmp2(ordering,VVMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(1))
                  if( VVMode.eq.ZZMode ) then
                      if( includeGammaStar ) then
                         call calcHelAmp2(ordering,ZgsMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(3))
                         call calcHelAmp2(ordering,gsZMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(5))
                         call calcHelAmp2(ordering,gsgsMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(7))
                      endif
                      if( includeVprime ) then
                          call calcHelAmp2(ordering,ZZpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(9))
                          call calcHelAmp2(ordering,ZpZMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(11))
                          call calcHelAmp2(ordering,ZpZpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(13))
                      endif
                      if( includeGammaStar .and. includeVprime ) then
                          call calcHelAmp2(ordering,gsZpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(15))
                          call calcHelAmp2(ordering,ZpgsMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(17))
                      endif
                  elseif( VVMode.eq.ZgMode ) then
                      if(includeGammaStar) then
                          call calcHelAmp2(ordering,gsgMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(3))
                      endif
                      if( includeVprime ) then
                          call calcHelAmp2(ordering,ZpgMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(5))
                      endif
                  elseif( VVMode.eq.WWMode .and. includeVprime ) then
                      call calcHelAmp2(ordering,WWpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(9))
                      call calcHelAmp2(ordering,WpWMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(11))
                      call calcHelAmp2(ordering,WpWpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(13))
                  endif

                  if( doInterference ) then
                      call calcHelAmp2(ordering_swap,VVMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(2))
                      if( includeGammaStar ) then
                          call calcHelAmp2(ordering_swap,ZgsMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(4))
                          call calcHelAmp2(ordering_swap,gsZMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(6))
                          call calcHelAmp2(ordering_swap,gsgsMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(8))
                      endif
                      if( includeVprime ) then
                          call calcHelAmp2(ordering_swap,ZZpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(10))
                          call calcHelAmp2(ordering_swap,ZpZMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(12))
                          call calcHelAmp2(ordering_swap,ZpZpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(14))
                      endif
                      if( includeGammaStar .and. includeVprime ) then
                          call calcHelAmp2(ordering_swap,gsZpMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(16))
                          call calcHelAmp2(ordering_swap,ZpgsMode,MY_IDUP,p(1:4,1:6),i3,i4,A_VV(18))
                      endif
                  endif

                  A0_VV(1) = A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7)+A_VV(9)+A_VV(11)+A_VV(13)+A_VV(15)+A_VV(17) ! 3456 pieces
                  A0_VV(2) = A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8)+A_VV(10)+A_VV(12)+A_VV(14)+A_VV(16)+A_VV(18) ! 5436 pieces
                  res = res + dreal(A0_VV(1)*dconjg(A0_VV(1)))
                  res = res + dreal(A0_VV(2)*dconjg(A0_VV(2)))
                  if( doInterference .and. (i3.eq.i4) ) then! interfere the 3456 with 5436 pieces
                      res = res - 2d0*intcolfac*dreal(  A0_VV(1)*dconjg(A0_VV(2))  ) ! minus from Fermi statistics
                  endif
          enddo;  enddo




! MADGRAPH CHECK
! call coupsm(0)
! if( (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!       call SH_EMEPEMEP((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,res2)
! else
!       call SH_EMEPEMEP_NOINT((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,res2)
! endif
! res2=res2* cdabs( (0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso) *  dconjg((0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso)) )
! res2=res2/100d0**2/100d0**2
! pause
! res=res2; RETURN

! res= res*prefactor
! print *, "checker 1",res
! print *, "checker 2",res2
! print *, "checker 1/2",res/res2
! pause


          res = res*prefactor
          if( (VVMode.eq.ZZMode) .and. doInterference ) res = res * SymmFac

      RETURN
      END SUBROUTINE

   subroutine calcHelAmp2(ordering,VVMode,MY_IDUP,p,i3,i4,A)
   implicit none
   integer :: ordering(1:4),VVMode,i3,i4,l1,l2,l3,l4,MY_IDUP(6:9)
   real(dp) :: p(1:4,1:6)
   real(dp) :: pin(4,4)
   complex(dp) :: A(1:1), sp(3:4,4)

      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      pin(1,:) = p(:,1)

      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call HZZampl(VVMode,pin,sp,A(1))
   end subroutine

      subroutine HZZampl(VVMode,p,sp,res)
      implicit none
      integer, intent(in) :: VVMode
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(3:4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e3_e4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q3(4),q4(4),q(4)
      complex(dp) :: e3(4),e4(4)
      complex(dp) :: yyy1,yyy2,yyy3,yyy4
      real(dp) :: q34
      real(dp) :: q_q, q3_q3, q4_q4
      complex(dp) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn



      res = 0d0
      q1 = dcmplx(p(1,:),0d0)
      q3 = dcmplx(p(3,:),0d0)
      q4 = dcmplx(p(4,:),0d0)


      e3 = sp(3,:)
      e4 = sp(4,:)

      q = -q1

      q_q =sc(q,q)
      q3_q3 = sc(q3,q3)
      q4_q4 = sc(q4,q4)

      e3_e4 = sc(e3,e4)
      e3_q4 = sc(e3,q4)
      e4_q3 = sc(e4,q3)


      if ((q_q).lt.-0.1d0 .or. (q3_q3).lt.-0.1d0 .or. (q4_q4).lt.-0.1d0) return  ! if negative invariant masses return zero


!---- data that defines couplings
      if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
         ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.gsZMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q3_q3,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q3_q3,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q3_q3,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q3_q3,q_q)
      elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.ZgsMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q_q)
      elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then
         ghz1_dyn = czero
         ghz2_dyn = HVVSpinZeroDynamicCoupling(9,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(10,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(11,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.ZZpMode) .or. (VVMode.eq.WWpMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(12,q3_q3,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(13,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(14,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(15,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.ZpZMode) .or. (VVMode.eq.WpWMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(12,q4_q4,q3_q3,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(13,q4_q4,q3_q3,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(14,q4_q4,q3_q3,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(15,q4_q4,q3_q3,q_q)
      elseif( (VVMode.eq.ZpZpMode) .or. (VVMode.eq.WpWpMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(16,q3_q3,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(17,q3_q3,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(18,q3_q3,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(19,q3_q3,q4_q4,q_q)
      elseif( (VVMode.eq.gsZpMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(20,0d0,q3_q3,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(21,0d0,q3_q3,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(22,0d0,q3_q3,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(23,0d0,q3_q3,q_q)
      elseif( (VVMode.eq.ZpgMode) .OR. (VVMode.eq.ZpgsMode) ) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(20,0d0,q4_q4,q_q)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(21,0d0,q4_q4,q_q)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(22,0d0,q4_q4,q_q)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(23,0d0,q4_q4,q_q)
      else
         ghz1_dyn = czero
         ghz2_dyn = czero
         ghz3_dyn = czero
         ghz4_dyn = czero
         print *,"VVMode",VVMode,"not implemented"
      endif

      if( .not. generate_as ) then
         yyy1 = ghz1_dyn*M_V**2/q_q &  ! in this line M_V is indeed correct, not a misprint
              + ghz2_dyn*(q_q-q3_q3-q4_q4)/q_q &
              + ghz3_dyn/Lambda**2*(q_q-q3_q3-q4_q4)*(q_q-q4_q4-q3_q3)/4d0/q_q
         yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(q_q-q3_q3-q4_q4)
         yyy3 = -2d0*ghz4_dyn
      else
         if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
            yyy1 = ahz1
            yyy2 = ahz2
            yyy3 = ahz3
         elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then! decay (gamma-gamma) OR (gamma*-gamma*) OR (gamma*-gamma)
            yyy1 = ahz1
            yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
            yyy3 = ahz3
         elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then! decay (Z-photon) OR (gamma*-Z) OR (Z-gamma*)
            yyy1 = ahz1
            yyy2 = -2*ahz1*q_q/(q_q-q3_q3)
            yyy3 = ahz3
         else
            yyy1=czero
            yyy2=czero
            yyy3=czero
            print *,"VVMode",VVMode,"not implemented in generate_as"
         endif
      endif

      res = e3_e4*q_q*yyy1                  &
          + e3_q4*e4_q3*yyy2                      &
          + et1(e3,e4,q3,q4)*yyy3


   END SUBROUTINE HZZampl


! Higgs decay to tau^+ tau^-   or    top anti-top
! Decay amplitude H --> tau^-(p1) + tau^+(p2)
! or              H --> tbar(p1)  + top(p2)
! with tau/top controlled by value of mass_F
! Since this is an isotropic scalar decay, just provide a matrix element SQUARED
! R.Rontsch July 2015
   SUBROUTINE EvalAmp_H_FF(pin,mass_F,res)
   implicit none
   real(dp), intent(out) ::  res
   complex(dp) :: amp2
   real(dp), intent(in) :: pin(1:4,1:2),mass_F
   real(dp)             :: s12

      s12=2d0*(pin(1,1)*pin(1,2)-pin(2,1)*pin(2,2)-pin(3,1)*pin(3,2)-pin(4,1)*pin(4,2)) + 2d0*mass_F**2
      amp2 =   2d0*s12*(kappa_tilde**2 + kappa**2) - 8d0*mass_F**2*kappa**2
      amp2 = amp2*mass_F**2/vev**2
      res = cdabs(amp2)

   RETURN
   END SUBROUTINE


! Decay amplitude H --> tau^-(-->l^-(p1)+nubar(p2)+nutau(p3))  +  tau^+(-->nu(p4)+l^+(p5)+nutaubar(p6))
! or H --> tbar^-(-->l^-(p1)+nubar(p2)+bbar(p3))+top(-->nu(p4)+l^+(p5)+b(p6))
! Note: 1. Value of m allows to change between H->tau^+tau^- and H->ttbar
!       2. Full kinematics of decaying particles -- no NWA.
!       3. Final state b quark is massless.
!       4. At present, no width in the tau/top propagator
! R. Rontsch July 2015
   SUBROUTINE EvalAmp_H_TT_decay(pin,mass_F,ga_F,res)
   implicit none
   real(dp), intent(out) ::  res
   real(dp), intent(in) :: pin(1:4,1:6),mass_F,ga_F
   integer              :: j
   real(dp)             :: p(1:6,1:4),s12,s45,s123,s456,s(6,6)
   complex(dp)          :: za(6,6),zb(6,6),amp,KL,KR

!       p=pin
      do j=1,6
          call convert_to_MCFM(pin(1:4,j),p(j,1:4))
      enddo
!       call spinoru(6,p,za,zb,s)
      call spinoru(p,za,zb,s)

      s12=s(1,2)
      s45=s(4,5)
      s123=s(1,2)+s(1,3)+s(2,3)
      s456=s(4,5)+s(4,6)+s(5,6)

      KL=-mass_F/vev*( kappa -(0d0,1d0)*kappa_tilde )
      KR=-mass_F/vev*( kappa +(0d0,1d0)*kappa_tilde )

      amp =  + KR * ( za(1,3)*za(1,4)*zb(1,2)*zb(5,6)- za(1,3)*za(3,4)*zb(2,3)*zb(5,6)) &
             + KL * (- za(1,3)*za(4,5)*zb(2,5)*zb(5,6)- za(1,3)*za(4,6)*zb(2,6)*zb(5,6))

! overall factors and propagators
      amp=amp/(s123-mass_F**2+ci*mass_F*ga_F)/(s456-mass_F**2+ci*mass_F*ga_F)/(s12-m_w**2+ci*m_w*Ga_W)/(s45-m_w**2+ci*m_w*Ga_W)
      amp=amp*16d0*ci*mass_F*gwsq**2
      res = cdabs(amp)**2

   RETURN
   END SUBROUTINE




subroutine getDecay_Couplings_Spinors_Props(VVMode,idordered,pordered,h3,h4, sp,pV)
   implicit none
   integer, intent(in) :: VVMode,idordered(6:9),h3,h4
   real(dp), intent(in) :: pordered(1:4,6:9)
   complex(dp), intent(out) :: sp(3:4,1:4)
   real(dp), intent(out) :: pV(3:4,1:4)
   real(dp) :: s
   complex(dp) :: propV(1:2), aL1,aR1,aL2,aR2

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
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = s/dcmplx(s - M_V**2,M_V*Ga_V)
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = s/dcmplx(s - M_V**2,M_V*Ga_V)



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
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = s/dcmplx(s - M_V**2,M_V*Ga_V)
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = s/dcmplx(s - M_V**2,M_V*Ga_V)


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
      aL2=1d0
      aR2=1d0
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_mless2(dcmplx(pordered(:,8)),h4,'out')  ! photon
      !sp(4,1:4)=pV(4,1:4); print *, "this checks gauge invariance"
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = s/dcmplx(s - M_V**2,M_V*Ga_V)
      propV(2)=1d0


   elseif( VVMode.eq.ggMode ) then
   !        gamma gamma DECAYS
      aL1=1d0
      aR1=1d0
      aL2=1d0
      aR2=1d0
      pV(3,:) = pordered(:,6)
      pV(4,:) = pordered(:,8)
      sp(3,:) = pol_mless2(dcmplx(pordered(:,6)),h3,'out')  ! photon
      sp(4,:) = pol_mless2(dcmplx(pordered(:,8)),h4,'out')  ! photon
      !sp(3,1:4)=pV(3,1:4); print *, "this checks gauge invariance"
      !sp(4,1:4)=pV(4,1:4)
      propV(1)=1d0
      propV(2)=1d0


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
      aL2=1d0
      aR2=1d0
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) ! photon propagator
      sp(4,:) = pol_mless2(dcmplx(pordered(:,8)),h4,'out')  ! photon
      !sp(4,1:4)=pV(4,1:4); print *, "this checks gauge invariance"
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = 1d0
      propV(2) = 1d0
      if( s.lt.MPhotonCutoff**2 ) propV(1)=czero


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
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:)
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = 1d0! = s/dcmplx(s)
      if( s.lt.MPhotonCutoff**2 ) propV(1)=czero
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = s/dcmplx(s - M_V**2,M_V*Ga_V)


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
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:)
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = s/dcmplx(s - M_V**2,M_V*Ga_V)
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = 1d0 ! = s/dcmplx(s)
      if( s.lt.MPhotonCutoff**2 ) propV(2)=czero


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
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:)
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:)
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = 1d0 ! = s/dcmplx(s)
      if( s.lt.MPhotonCutoff**2 ) propV(1)=czero
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = 1d0 ! = s/dcmplx(s)
      if( s.lt.MPhotonCutoff**2 ) propV(2)=czero


   elseif( VVMode.eq.ZpZMode ) then
   !        Z'Z DECAYS
      aL1 = VpffCoupling(idordered(6),-1,.false.)
      aR1 = VpffCoupling(idordered(6),+1,.false.)
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

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      if( M_Vprime .gt. 0d0 ) then
        propV(1) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(1) = 1d0
      else
        propV(1) = s/M_Z**2
      endif
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = s/dcmplx(s - M_V**2,M_V*Ga_V)


   elseif( VVMode.eq.ZZpMode ) then
   !        ZZ' DECAYS
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
      aL2 = VpffCoupling(idordered(8),-1,.false.)
      aR2 = VpffCoupling(idordered(8),+1,.false.)

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = s/dcmplx(s - M_V**2,M_V*Ga_V)
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      if( M_Vprime .gt. 0d0 ) then
        propV(2) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(2) = 1d0
      else
        propV(2) = s/M_Z**2
      endif


   elseif( VVMode.eq.ZpZpMode ) then
   !        Z'Z' DECAYS
      aL1 = VpffCoupling(idordered(6),-1,.false.)
      aR1 = VpffCoupling(idordered(6),+1,.false.)
      aL2 = VpffCoupling(idordered(8),-1,.false.)
      aR2 = VpffCoupling(idordered(8),+1,.false.)

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      if( M_Vprime .gt. 0d0 ) then
        propV(1) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(1) = 1d0
      else
        propV(1) = s/M_Z**2
      endif
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      if( M_Vprime .gt. 0d0 ) then
        propV(2) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(2) = 1d0
      else
        propV(2) = s/M_Z**2
      endif


   elseif( VVMode.eq.ZpgsMode ) then
   !        Z' gamma* DECAYS
      aL1 = VpffCoupling(idordered(6),-1,.false.)
      aR1 = VpffCoupling(idordered(6),+1,.false.)
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

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:)
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      if( M_Vprime .gt. 0d0 ) then
        propV(1) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(1) = 1d0
      else
        propV(1) = s/M_Z**2
      endif
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      propV(2) = 1d0 ! = s/dcmplx(s)
      if( s.lt.MPhotonCutoff**2 ) propV(2)=czero

   elseif( VVMode.eq.gsZpMode ) then
   !        gamma* Z' DECAYS
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
      aL2 = VpffCoupling(idordered(8),-1,.false.)
      aR2 = VpffCoupling(idordered(8),+1,.false.)

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:)
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      propV(1) = 1d0! = s/dcmplx(s)
      if( s.lt.MPhotonCutoff**2 ) propV(1)=czero
      s = scr(pordered(:,8)+pordered(:,9),pordered(:,8)+pordered(:,9))
      if( M_Vprime .gt. 0d0 ) then
        propV(2) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(2) = 1d0
      else
        propV(2) = s/M_Z**2
      endif


   elseif( VVMode.eq.ZpgMode ) then
   !        Z' gamma DECAYS
      aL1 = VpffCoupling(idordered(6),-1,.false.)
      aR1 = VpffCoupling(idordered(6),+1,.false.)
      aL2=1d0
      aR2=1d0
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_mless2(dcmplx(pordered(:,8)),h4,'out')  ! photon
      s = scr(pordered(:,6)+pordered(:,7),pordered(:,6)+pordered(:,7))
      if( M_Vprime .gt. 0d0 ) then
        propV(1) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(1) = 1d0
      else
        propV(1) = s/M_Z**2
      endif
      propV(2)=1d0


   elseif( VVMode.eq.WpWMode ) then
   !        W'W DECAYS
      aL1 = VpffCoupling(idordered(6),-1,.true.)
      aR1 = VpffCoupling(idordered(6),+1,.true.)
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
      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pV(3,:),pV(3,:))
      if( M_Vprime .gt. 0d0 ) then
        propV(1) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(1) = 1d0
      else
        propV(1) = s/M_W**2
      endif
      s = scr(pV(4,:),pV(4,:))
      propV(2) = s/dcmplx(s - M_V**2,M_V*Ga_V)


   elseif( VVMode.eq.WWpMode ) then
   !        WW' DECAYS
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
      aL2 = VpffCoupling(idordered(8),-1,.true.)
      aR2 = VpffCoupling(idordered(8),+1,.true.)

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pV(3,:),pV(3,:))
      propV(1) = s/dcmplx(s - M_V**2,M_V*Ga_V)
      s = scr(pV(4,:),pV(4,:))
      if( M_Vprime .gt. 0d0 ) then
        propV(2) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(2) = 1d0
      else
        propV(2) = s/M_W**2
      endif


   elseif( VVMode.eq.WpWpMode ) then
   !        W'W' DECAYS
      aL1 = VpffCoupling(idordered(6),-1,.true.)
      aR1 = VpffCoupling(idordered(6),+1,.true.)
      aL2 = VpffCoupling(idordered(8),-1,.true.)
      aR2 = VpffCoupling(idordered(8),+1,.true.)

      pV(3,:) = pordered(:,6)+pordered(:,7)
      pV(4,:) = pordered(:,8)+pordered(:,9)
      sp(3,:) = pol_dk2mom(dcmplx(pordered(:,6)),dcmplx(pordered(:,7)),h3)  ! ubar(l1), v(l2)
      sp(3,:) = -sp(3,:) + pV(3,:)*( sc(sp(3,:),dcmplx(pV(3,:))) )/scr(pV(3,:),pV(3,:))! full propagator numerator
      sp(4,:) = pol_dk2mom(dcmplx(pordered(:,8)),dcmplx(pordered(:,9)),h4)  ! ubar(l3), v(l4)
      sp(4,:) = -sp(4,:) + pV(4,:)*( sc(sp(4,:),dcmplx(pV(4,:))) )/scr(pV(4,:),pV(4,:))! full propagator numerator
      s = scr(pV(3,:),pV(3,:))
      if( M_Vprime .gt. 0d0 ) then
        propV(1) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(1) = 1d0
      else
        propV(1) = s/M_W**2
      endif
      s = scr(pV(4,:),pV(4,:))
      if( M_Vprime .gt. 0d0 ) then
        propV(2) = s/dcmplx(s - M_Vprime**2,M_Vprime*Ga_Vprime)
      elseif( M_Vprime .eq. 0d0 ) then
        propV(2) = 1d0
      else
        propV(2) = s/M_W**2
      endif


   else
      call Error("Unsupported decay modes")
   endif

   !print *,"sp(3)=",sp(3,:)
   !print *,"propV(1)=",propV(1)
   !print *,"aL1,aR1=",aL1,aR1
   !print *,"sp(4)=",sp(4,:)
   !print *,"propV(2)=",propV(2)
   !print *,"aL2,aR2=",aL2,aR2

   sp(3,:) = sp(3,:)*propV(1)
   sp(4,:) = sp(4,:)*propV(2)
   if (h3.eq.-1) then
      sp(3,:) = aL1 * sp(3,:)
   elseif(h3.eq.1) then
      sp(3,:) = aR1 * sp(3,:)
   endif
   if (h4.eq.-1) then
      sp(4,:) = aL2 * sp(4,:)
   elseif(h4.eq.1) then
      sp(4,:) = aR2 * sp(4,:)
   endif

   return
end subroutine

subroutine getDecay_VVMode_Ordering(MY_IDUP, VVMode,ordering,VVMode_swap,ordering_swap)
   implicit none
   integer, intent(in) :: MY_IDUP(6:9)
   integer, intent(out) :: VVMode,ordering(1:4),VVMode_swap,ordering_swap(1:4)
   integer :: idV(1:2),idV_swap(1:2)

   ordering=(/3,4,5,6/)
   idV(1)=CoupledVertex(MY_IDUP(6:7),-1)
   idV(2)=CoupledVertex(MY_IDUP(8:9),-1)
   if(MY_IDUP(6).eq.Pho_ .or. MY_IDUP(7).eq.Pho_) idV(1)=Pho_
   if(MY_IDUP(8).eq.Pho_ .or. MY_IDUP(9).eq.Pho_) idV(2)=Pho_
   if(convertLHE(MY_IDUP(6)).lt.0 .or. MY_IDUP(6).eq.Not_a_particle_) then
      call swap(ordering(1),ordering(2))
   endif
   if(convertLHE(MY_IDUP(8)).lt.0 .or. MY_IDUP(8).eq.Not_a_particle_) then
      call swap(ordering(3),ordering(4))
   endif
   if( &
         (idV(1).eq.Wm_ .and. idV(2).eq.Wp_) .or. &
         (idV(2).eq.Z0_ .and. idV(1).eq.Pho_) &
     ) then
      call swap(ordering(1),ordering(3))
      call swap(ordering(2),ordering(4))
      call swap(idV(1),idV(2))
   endif
   ordering_swap(:)=ordering(:)
   call swap(ordering_swap(1),ordering_swap(3))

   idV_swap(1) = CoupledVertex( (/ MY_IDUP(3+ordering_swap(1)), MY_IDUP(3+ordering_swap(2)) /), -1)
   idV_swap(2) = CoupledVertex( (/ MY_IDUP(3+ordering_swap(3)), MY_IDUP(3+ordering_swap(4)) /), -1)
   if ( (idV_swap(1).eq.Wm_) .and. (idV_swap(2).eq.Wp_) ) then
      call swap(ordering_swap(1),ordering_swap(3))
      call swap(ordering_swap(2),ordering_swap(4))
      call swap(idV_swap(1),idV_swap(2))
   endif

   if(idV(1).eq.Z0_ .and. idV(2).eq.Z0_) then
      VVMode=ZZMode
   elseif(idV(1).eq.Z0_ .and. idV(2).eq.Pho_) then
      VVMode=ZgMode
   elseif(idV(1).eq.Pho_ .and. idV(2).eq.Pho_) then
      VVMode=ggMode
   elseif(idV(1).eq.Wp_ .and. idV(2).eq.Wm_) then
      VVMode=WWMode
   else
      print *,"idV=",idV
      call Error("Unsupported decay mode")
   endif

   VVMode_swap=InvalidMode
   if(idV_swap(1).eq.Z0_ .and. idV_swap(2).eq.Z0_) then
      VVMode_swap=ZZMode
   elseif(idV_swap(1).eq.Wp_ .and. idV_swap(2).eq.Wm_) then
      VVMode_swap=WWMode
   endif
   return
end subroutine


END MODULE


