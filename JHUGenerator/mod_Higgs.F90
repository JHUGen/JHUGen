      module modHiggs
      use ModParameters
      implicit none
      private
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), private, parameter :: tol = 0.00000010_dp


!----- notation for subroutines
      public :: EvalAmp_gg_H_VV,EvalAmp_H_VV

      contains


!----- a subroutinefor gg -> H -> ZZ/WW/gammagamma
!----- all outgoing convention and the following momentum assignment
!-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_gg_H_VV(p,MY_IDUP,sum)
      use ModMisc
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A(1:4)
      integer :: i1,i2,i3,i4,ordering(1:4)
      real(dp) :: aL1,aR1,aL2,aR2,sum2
      real(dp) :: gZ_sq
      real(dp) :: prefactor, Lambda_inv
      real(dp), parameter :: symmFact=1d0/2d0

      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2

!---- the 1/Lambda coupling
      Lambda_inv = 1.0d0/Lambda

!---- full prefactor; 8 is  the color factor
      prefactor = 8d0*(Lambda_inv**2)**2*gZ_sq**2


         if( IsAZDecay(DecayMode1) ) then!  Z decay
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep
                    aR1=aR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu
                    aR1=aR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp
                    aR1=aR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn
                    aR1=aR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
              prefactor = prefactor *(one/two*M_V*Ga_V)**2
         elseif( IsAWDecay(DecayMode1) ) then !  W decay
              aL1 = bL
              aR1 = bR
              prefactor = prefactor *(one/two*M_V*Ga_V)**2
         elseif( IsAPhoton(DecayMode1) ) then !  photon "decay"
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq**2! cancel the overall z coupling
         else
              aL1=0d0
              aR1=0d0            
         endif

         if( IsAZDecay(DecayMode2) ) then!  Z decay
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep
                    aR2=aR_lep
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu
                    aR2=aR_neu
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp
                    aR2=aR_QUp
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn
                    aR2=aR_QDn
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( IsAWDecay(DecayMode2) ) then !  W decay
              aL2 = bL
              aR2 = bR
         elseif( IsAPhoton(DecayMode2) ) then !  photon "decay"
              aL2=1d0
              aR2=1d0 
         else
              aL2=0d0
              aR2=0d0  
         endif



! ! MADGRAPH CHECK
! sum=0d0
! if (MY_IDUP(6).ne.MY_IDUP(8) ) return
! if (MY_IDUP(7).ne.MY_IDUP(9) ) return
! if (MY_IDUP(6).eq.MY_IDUP(8) ) return
! if (MY_IDUP(7).eq.MY_IDUP(9) ) return
! print *, "MY COUPL",al1*dsqrt(gZ_sq),ar1*dsqrt(gZ_sq)
! ar1=-0.158480099490745d0! this is MadGraphs gzl(R)  , the MG couplings differ by a global minus sign, relative differences are because of different input parameters
! al1=+0.210150647402957d0! this is MadGraphs gzl(L) 
! al2=al1
! ar2=ar1
! print *, "MG COUPL",al1,ar1
! pause

sum = zero
do i1 = 1,2
do i2 = 1,2
do i3 = 1,2
do i4 = 1,2
   
         ordering = (/3,4,5,6/)
         call calcHelAmp(ordering,p(1:4,1:6),i1,i2,i3,i4,A(1))

         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             ordering = (/5,4,3,6/)
             call calcHelAmp(ordering,p(1:4,1:6),i1,i2,i3,i4,A(2))
             A(2) = -A(2) ! minus comes from fermi statistics
         endif


         if (i3.eq.1) then
            A(:) = aL1*A(:)
         elseif(i3.eq.2) then
            A(:) = aR1*A(:)
         endif
         if (i4.eq.1) then
            A(:) = aL2*A(:)
         elseif(i4.eq.2) then
            A(:) = aR2*A(:)
         endif

         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             sum = sum + symmFact * (cdabs( A(1)*dconjg(A(1)) ) + cdabs( A(2)*dconjg(A(2)) ))
             if( i3.eq.i4 ) sum = sum + symmFact * 2d0*dreal(A(1)*dconjg(A(2)))  
         else
             sum = sum + cdabs( A(1)*dconjg(A(1)) )
         endif
enddo
enddo
enddo
enddo


! MADGRAPH CHECK
! call coupsm(0)
! if( (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!       call SH_EMEPEMEP((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,sum2)
! else
!       call SH_EMEPEMEP_NOINT((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,sum2)
!       call SH_TAMTAPTAMTAP_NOINT((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,sum2)
! endif
! sum2=sum2* cdabs( (0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso) *  dconjg((0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso)) ) 
! sum2=sum2/100d0**2/100d0**2
! pause
! SUM=SUM2; RETURN

! sum= sum*prefactor/(Lambda_inv**2)**2
! print *, "checker 1",sum
! print *, "checker 2",sum2
! print *, "checker 1/2",sum/sum2
! pause



      sum = sum*prefactor


      end subroutine





     subroutine calcHelAmp(ordering,p,i1,i2,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i1,i2,i3,i4,l1,l2,l3,l4
     real(dp) :: p(1:4,1:6)
     complex(dp) :: propG, propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(4,4)


      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      s  = 2d0 * scr(p(:,1),p(:,2))
      propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)


         pin(1,:) = p(:,1)
         pin(2,:) = p(:,2)
         sp(1,:) = pol_mless2(dcmplx(p(:,1)),-3+2*i1,'in')  ! gluon
         sp(2,:) = pol_mless2(dcmplx(p(:,2)),-3+2*i2,'in')  ! gluon
!          sp(1,1:4)=pin(1,1:4);print *, "this checks IS gauge invariance"
!          sp(2,1:4)=pin(2,1:4);print *, "this checks IS gauge invariance"


!-------- -1 == left, 1 == right
         if( (.not.IsAPhoton(DecayMode1)) .and. (.not.IsAPhoton(DecayMode2)) ) then
            pin(3,:) = p(:,l1)+p(:,l2)
            pin(4,:) = p(:,l3)+p(:,l4)
            sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
            sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
            sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
            sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
!print *, "ubar, v, ubar, v"
!print *, "masses",scr(p(:,l1),p(:,l1)),scr(p(:,l2),p(:,l2)),scr(p(:,l3),p(:,l3)),scr(p(:,l4),p(:,l4))
!print *, "check", sc(sp(3,:),dcmplx(pin(3,:))), sc(sp(4,:),dcmplx(pin(4,:)))
!pause
            s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
            propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
            s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
            propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
            pin(3,:) = p(:,l1)
            pin(4,:) = p(:,l3)
            sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
            sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!             sp(3,1:4)=pin(3,1:4); print *, "this checks gauge invariance"
!             sp(4,1:4)=pin(4,1:4)
            propz1=1d0
            propz2=1d0
         elseif( (.not.IsAPhoton(DecayMode1)) .and. (IsAPhoton(DecayMode2)) ) then
            pin(3,:) = p(:,l1)+p(:,l2)
            pin(4,:) = p(:,l3)
            sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
            sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
            sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!             sp(4,1:4)=pin(4,1:4); print *, "this checks gauge invariance"
            s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
            propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
            propZ2=1d0
         endif


         if( OffShellReson ) then
              call ggOffHZZampl(pin,sp,A(1))
         else
              call ggHZZampl(pin,sp,A(1))
         endif

         A(1) = A(1) * propG*propZ1*propZ2



     end subroutine





      subroutine ggHZZampl(p,sp,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q_q
      complex(dp) :: q1_q2,q1_q3,q1_q4
      complex(dp) :: q2_q3,q2_q4
      complex(dp) :: q3_q4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4)
      complex(dp) :: xxx1,xxx2,xxx3,yyy1,yyy2,yyy3,yyy4
      complex(dp) :: ghg2_dyn,ghg3_dyn,ghg4_dyn,ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      real(dp) :: q34
      real(dp) :: MG, MZ3, MZ4, q3_q3, q4_q4



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


      q1_q2 = sc(q1,q2)
      q1_q3 = sc(q1,q3)
      q1_q4 = sc(q1,q4)
      q2_q3 = sc(q2,q3)
      q2_q4 = sc(q2,q4)
      q3_q4 = sc(q3,q4)

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


      if (cdabs(q_q).lt.-0.1d0.or.(q3_q3).lt.-0.1d0.or.(q4_q4).lt.-0.1d0) return  ! if negative invariant masses return zero
      MG =dsqrt(cdabs(q_q))
      MZ3=dsqrt(dabs(q3_q3))
      MZ4=dsqrt(dabs(q4_q4))


!---- data that defines couplings
  if( (IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1)) .and. (IsAZDecay(DecayMode2) .or. IsAWDecay(DecayMode2)) ) then! decay into ZZ's or WW's

    if( generate_as ) then 
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = ahz2
      yyy3 = ahz3
    else
      ghg2_dyn = ghg2
      ghg3_dyn = ghg3
      ghg4_dyn = ghg4
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))
      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))
      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))
      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))

      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 = ghz1_dyn*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
           + ghz2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghz3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
      yyy3 = -2d0*ghz4_dyn
    endif
      res = e1_e2*e3_e4*M_Reso**4*yyy1*xxx1                  &
          + e1_e2*e3_q4*e4_q3*M_Reso**2*yyy2*xxx1            &
          + et1(e1,e2,q1,q2)*e3_e4*M_Reso**2*yyy1*xxx3       &
          + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3           &
          + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3      &
          + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1



  elseif( (IsAPhoton(DecayMode1)) .and. (IsAPhoton(DecayMode2)) ) then! decay into photons

    if( generate_as ) then 
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
      yyy3 = ahz3
    else
      ghg2_dyn = ghg2
      ghg3_dyn = ghg3
      ghg4_dyn = ghg4
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))
      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))
      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))
      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))

      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 = ghz1_dyn*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
           + ghz2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghz3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
      yyy3 = -2d0*ghz4_dyn
    endif
     res = e1_e2*e3_e4*M_Reso**4*yyy1*xxx1                  &
         + e1_e2*e3_q4*e4_q3*M_Reso**2*yyy2*xxx1            &
         + et1(e1,e2,q1,q2)*e3_e4*M_Reso**2*yyy1*xxx3       &
         + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3           &
         + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3      &
         + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1


  elseif( (IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1)) .and. (IsAPhoton(DecayMode2)) ) then! decay into Z+photon
    if( generate_as ) then
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = -2*ahz1*MG**2/(MG**2-MZ3**2)
      yyy3 = ahz3
    else
      xxx1 = ghg2+ghg3/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4
      yyy1 = 0 * ghz1*M_V**2/MG**2 &  ! removed ghz1 dependence because it does not contribute
           + ghz2*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghz3/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = (-2d0*ghz2-ghz3/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2) )
      yyy3 = -2d0*ghz4
    endif
      res = e1_e2*e3_e4*M_Reso**4*yyy1*xxx1                  &
          + e1_e2*e3_q4*e4_q3*M_Reso**2*yyy2*xxx1            &
          + et1(e1,e2,q1,q2)*e3_e4*M_Reso**2*yyy1*xxx3       &
          + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3           &
          + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3      &
          + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1
  endif


      end subroutine ggHZZampl







      subroutine ggOffHZZampl(p,sp,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q_q
      complex(dp) :: q1_q2,q1_q3,q1_q4
      complex(dp) :: q2_q3,q2_q4
      complex(dp) :: q3_q4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4)
      complex(dp) :: xxx1,xxx2,xxx3,yyy1,yyy2,yyy3,yyy4
      complex(dp) :: ghg2_dyn,ghg3_dyn,ghg4_dyn,ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      real(dp) :: q34, MG, MZ3, MZ4, q3_q3, q4_q4

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


      if (cdabs(q_q).lt.-0.1d0.or.(q3_q3).lt.-0.1d0.or.(q4_q4).lt.-0.1d0) return  ! if negative invariant masses return zero
      MG =dsqrt(cdabs(q_q))
      MZ3=dsqrt(dabs(q3_q3))
      MZ4=dsqrt(dabs(q4_q4))

      q1_q2 = sc(q1,q2)
      q1_q3 = sc(q1,q3)
      q1_q4 = sc(q1,q4)
      q2_q3 = sc(q2,q3)
      q2_q4 = sc(q2,q4)
      q3_q4 = sc(q3,q4)



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


!---- data that defines couplings
  if( (IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1)) .and. (IsAZDecay(DecayMode2) .or. IsAWDecay(DecayMode2)) ) then! decay into ZZ's or WW's

    if( generate_as ) then 
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = ahz2
      yyy3 = ahz3
    else
      ghg2_dyn = ghg2
      ghg3_dyn = ghg3
      ghg4_dyn = ghg4
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))
      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))
      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))
      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))

      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 = ghz1_dyn*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
           + ghz2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghz3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
      yyy3 = -2d0*ghz4_dyn
    endif
     res = e1_e2*e3_e4*MG**4*yyy1*xxx1                    &
         + e1_e2*e3_q4*e4_q3*MG**2*yyy2*xxx1              &
         + et1(e1,e2,q1,q2)*e3_e4*MG**2*yyy1*xxx3         &
         + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3         &
         + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3    &
         + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1



  elseif( (IsAPhoton(DecayMode1)) .and. (IsAPhoton(DecayMode2)) ) then! decay into photons

    if( generate_as ) then 
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
      yyy3 = ahz3
    else
      ghg2_dyn = ghg2
      ghg3_dyn = ghg3
      ghg4_dyn = ghg4
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))
      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))
      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))
      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))

      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 = ghz1_dyn*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
           + ghz2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghz3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
      yyy3 = -2d0*ghz4_dyn
    endif
     res = e1_e2*e3_e4*MG**4*yyy1*xxx1                    &
         + e1_e2*e3_q4*e4_q3*MG**2*yyy2*xxx1              &
         + et1(e1,e2,q1,q2)*e3_e4*MG**2*yyy1*xxx3         &
         + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3         &
         + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3    &
         + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1

  elseif( (IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1)) .and. (IsAPhoton(DecayMode2)) ) then! decay into Z+photon

   print *, "Zgamma FS with off-shell H is not yet supported"
   stop

    if( generate_as ) then
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = -2*ahz1*MG**2/(MG**2-MZ3**2)
      yyy3 = ahz3
    else
      xxx1 = ghg2+ghg3/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4
      yyy1 = 0 * ghz1*M_V**2/MG**2 &  ! removed ghz1 dependence because it does not contribute
           + ghz2*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghz3/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = (-2d0*ghz2-ghz3/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2) )
      yyy3 = -2d0*ghz4
    endif
      res = e1_e2*e3_e4*M_Reso**4*yyy1*xxx1                  &
          + e1_e2*e3_q4*e4_q3*M_Reso**2*yyy2*xxx1            &
          + et1(e1,e2,q1,q2)*e3_e4*M_Reso**2*yyy1*xxx3       &
          + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3           &
          + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3      &
          + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1

  endif


      end subroutine ggOffHZZampl





!----- a subroutinefor H -> ZZ/WW/gammagamma
!----- all outgoing convention and the following momentum assignment
!-----  0 -> Higgs(p1) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_H_VV(p,MY_IDUP,sum)
      use ModMisc
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A(1:4)
      integer :: i1,i2,i3,i4,ordering(1:4)
      real(dp) :: aL1,aR1,aL2,aR2
      real(dp) :: gZ_sq
      real(dp) :: prefactor, Lambda_inv
      real(dp), parameter :: symmFact=1d0/2d0

      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2

!---- the 1/Lambda coupling
      Lambda_inv = 1.0d0/Lambda

!---- full prefactor; 8 is  the color factor
!       prefactor = 8d0*(Lambda_inv**2)**2*gZ_sq**2
      prefactor = (Lambda_inv**2)**2*gZ_sq**2


         if( IsAZDecay(DecayMode1) ) then!  Z decay
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep
                    aR1=aR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu
                    aR1=aR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp
                    aR1=aR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn
                    aR1=aR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
              prefactor = prefactor *(one/two*M_V*Ga_V)**2
         elseif( IsAWDecay(DecayMode1) ) then !  W decay
              aL1 = bL
              aR1 = bR
              prefactor = prefactor *(one/two*M_V*Ga_V)**2
         elseif( IsAPhoton(DecayMode1) ) then !  photon "decay"
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq**2! cancel the overall z coupling
         else
              aL1=0d0
              aR1=0d0            
         endif

         if( IsAZDecay(DecayMode2) ) then!  Z decay
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep
                    aR2=aR_lep
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu
                    aR2=aR_neu
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp
                    aR2=aR_QUp
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn
                    aR2=aR_QDn
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( IsAWDecay(DecayMode2) ) then !  W decay
              aL2 = bL
              aR2 = bR
         elseif( IsAPhoton(DecayMode2) ) then !  photon "decay"
              aL2=1d0
              aR2=1d0 
         else
              aL2=0d0
              aR2=0d0  
         endif


! ! MADGRAPH CHECK
! sum=0d0
! if (MY_IDUP(6).ne.MY_IDUP(8) ) return
! if (MY_IDUP(7).ne.MY_IDUP(9) ) return
! if (MY_IDUP(6).eq.MY_IDUP(8) ) return
! if (MY_IDUP(7).eq.MY_IDUP(9) ) return
! print *, "MY COUPL",al1*dsqrt(gZ_sq),ar1*dsqrt(gZ_sq)
! ar1=-0.158480099490745d0! this is MadGraphs gzl(R)  , the MG couplings differ by a global minus sign, relative differences are because of different input parameters
! al1=+0.210150647402957d0! this is MadGraphs gzl(L) 
! al2=al1
! ar2=ar1
! print *, "MG COUPL",al1,ar1
! pause



sum = zero
do i3 = 1,2
do i4 = 1,2
   
         ordering = (/3,4,5,6/)
         call calcHelAmp2(ordering,p(1:4,1:6),i3,i4,A(1))
         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             ordering = (/5,4,3,6/)
             call calcHelAmp2(ordering,p(1:4,1:6),i3,i4,A(2))
             A(2) = -A(2) ! minus comes from fermi statistics
         endif


         if (i3.eq.1) then
            A(:) = aL1*A(:)
         elseif(i3.eq.2) then
            A(:) = aR1*A(:)
         endif
         if (i4.eq.1) then
            A(:) = aL2*A(:)
         elseif(i4.eq.2) then
            A(:) = aR2*A(:)
         endif

         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             sum = sum + symmFact * (cdabs( A(1)*dconjg(A(1)) ) + cdabs( A(2)*dconjg(A(2)) ))
             if( i3.eq.i4 ) sum = sum + symmFact * 2d0*dreal(A(1)*dconjg(A(2)))  
         else
             sum = sum + cdabs( A(1)*dconjg(A(1)) )
         endif
enddo
enddo


! MADGRAPH CHECK
! call coupsm(0)
! if( (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
!       call SH_EMEPEMEP((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,sum2)
! else
!       call SH_EMEPEMEP_NOINT((/-P(1:4,1)-P(1:4,2),P(1:4,3),P(1:4,4),P(1:4,5),P(1:4,6)/)*100d0,sum2)
! endif
! sum2=sum2* cdabs( (0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso) *  dconjg((0d0,1d0)/dcmplx(2d0*scr(p(:,1),p(:,2))-M_Reso**2,M_Reso*Ga_Reso)) ) 
! sum2=sum2/100d0**2/100d0**2
! pause
! SUM=SUM2; RETURN

! sum= sum*prefactor/(Lambda_inv**2)**2
! print *, "checker 1",sum
! print *, "checker 2",sum2
! print *, "checker 1/2",sum/sum2
! pause



      sum = sum*prefactor

      end subroutine





     subroutine calcHelAmp2(ordering,p,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i3,i4,l1,l2,l3,l4
     real(dp) :: p(1:4,1:6)
     complex(dp) :: propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(3:4,4)


      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

         pin(1,:) = p(:,1)

!-------- -1 == left, 1 == right
         if( .not.IsAPhoton(DecayMode1) ) then 
            pin(3,:) = p(:,l1)+p(:,l2)
            pin(4,:) = p(:,l3)+p(:,l4)
            sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
            sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
            sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
            sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
            s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
            propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
            s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
            propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         elseif( IsAPhoton(DecayMode1) ) then 
            pin(3,:) = p(:,l1)
            pin(4,:) = p(:,l3)
            sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
            sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!             sp(3,1:4)=pin(3,1:4)! this checks gauge invariance
!             sp(4,1:4)=pin(4,1:4)
            propz1=1d0
            propz2=1d0
         endif


         call HZZampl(pin,sp,A(1))
         A(1) = A(1) * propZ1*propZ2

     end subroutine





      subroutine HZZampl(p,sp,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(3:4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e3_e4
      complex(dp) :: q_q
      complex(dp) :: q3_q4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q3(4),q4(4),q(4)
      complex(dp) :: e3(4),e4(4)
      complex(dp) :: yyy1,yyy2,yyy3,yyy4
      real(dp) :: q34
      real(dp) :: MG, MZ3, MZ4, q3_q3, q4_q4



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

      q3_q4 = sc(q3,q4)
      e3_e4 = sc(e3,e4)
      e3_q4 = sc(e3,q4)
      e4_q3 = sc(e4,q3)


      if (cdabs(q_q).lt.-0.1d0.or.(q3_q3).lt.-0.1d0.or.(q4_q4).lt.-0.1d0) return  ! if negative invariant masses return zero
      MG =dsqrt(cdabs(q_q))
      MZ3=dsqrt(dabs(q3_q3))
      MZ4=dsqrt(dabs(q4_q4))


!---- data that defines couplings
  if( IsAZDecay(DecayMode1) .or. IsAWDecay(DecayMode1) ) then! decay into Z's or W's

      if( generate_as ) then 
        yyy1 = ahz1
        yyy2 = ahz2
        yyy3 = ahz3
      else
        yyy1 = ghz1*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
             + ghz2*(MG**2-MZ3**2-MZ4**2)/MG**2 &
             + ghz3/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
        yyy2 = -2d0*ghz2-ghz3/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
        yyy3 = -2d0*ghz4
      endif

      res = e3_e4*M_Reso**2*yyy1       &
          + e3_q4*e4_q3*yyy2           &
          + et1(e3,e4,q3,q4)*yyy3 


  elseif( IsAPhoton(DecayMode1) ) then! decay into photons

      if( generate_as ) then 
        yyy1 = ahz1
        yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
        yyy3 = ahz3
      else
        yyy1 = ghz1*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
            + ghz2*(MG**2-MZ3**2-MZ4**2)/MG**2 &
            + ghz3/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
        yyy2 = -2d0*ghz2-ghz3/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
        yyy3 = -2d0*ghz4
      endif

      res = e3_e4*M_Reso**2*yyy1       &
          + e3_q4*e4_q3*yyy2           &
          + et1(e3,e4,q3,q4)*yyy3 



  endif


      end subroutine HZZampl








   double complex function et1(e1,e2,e3,e4)
    implicit none
    complex(dp), intent(in) :: e1(4), e2(4), e3(4), e4(4)

    et1 =  e1(1)*e2(2)*e3(3)*e4(4)-e1(1)*e2(2)*e3(4)*e4(3) &
          -e1(1)*e2(3)*e3(2)*e4(4)+e1(1)*e2(3)*e3(4)*e4(2) &
          +e1(1)*e2(4)*e3(2)*e4(3)-e1(1)*e2(4)*e3(3)*e4(2) &
          -e1(2)*e2(1)*e3(3)*e4(4)+e1(2)*e2(1)*e3(4)*e4(3) &
          +e1(2)*e2(3)*e3(1)*e4(4)-e1(2)*e2(3)*e3(4)*e4(1) &
          -e1(2)*e2(4)*e3(1)*e4(3)+e1(2)*e2(4)*e3(3)*e4(1) &
          +e1(3)*e2(1)*e3(2)*e4(4)-e1(3)*e2(1)*e3(4)*e4(2) &
          -e1(3)*e2(2)*e3(1)*e4(4)+e1(3)*e2(2)*e3(4)*e4(1) &
          +e1(3)*e2(4)*e3(1)*e4(2)-e1(3)*e2(4)*e3(2)*e4(1) &
          -e1(4)*e2(1)*e3(2)*e4(3)+e1(4)*e2(1)*e3(3)*e4(2) &
          +e1(4)*e2(2)*e3(1)*e4(3)-e1(4)*e2(2)*e3(3)*e4(1) &
          -e1(4)*e2(3)*e3(1)*e4(2)+e1(4)*e2(3)*e3(2)*e4(1)

   return
   end function et1


      double complex function sc(q1,q2)
        complex(dp), intent(in) :: q1(4)
        complex(dp), intent(in) :: q2(4)

        sc = q1(1)*q2(1) - q1(2)*q2(2)-q1(3)*q2(3) -q1(4)*q2(4)

      end function sc

      double precision function scr(p1,p2)
        real(dp), intent(in) :: p1(4),p2(4)

        scr = p1(1)*p2(1) - p1(2)*p2(2)-p1(3)*p2(3) -p1(4)*p2(4)

      end function scr

!---- THESE ARE POLARIZATION ROUTINES

  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
    complex(dp), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------
    integer :: pol
    real(dp) :: p0,px,py,pz
    real(dp) :: pv,ct,st,cphi,sphi
    complex(dp) :: pol_mless(4)

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
!^^^END


    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0_dp-ct**2))

    if (st < tol) then
       cphi=1.0_dp
       sphi=0.0_dp
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0_dp) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif

    pol_mless(1)=czero
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless


  function pol_mless2(p,i,out)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    character(len=*), intent(in):: out
    complex(dp)             :: pol_mless2(4)
    ! -------------------------------------

    if (out == 'out') then
       pol_mless2 = pol_mless(p,i,outgoing=.true.)
    else
       pol_mless2 = pol_mless(p,i,outgoing=.false.)
    endif
  end function pol_mless2




  function pol_dk2mom(plepton,antilepton,i,outgoing)
  use ModMisc
  implicit none
    integer, intent(in) :: i
    integer :: j
    complex(dp), intent(in) :: plepton(1:4),antilepton(1:4)
    logical, intent(in),optional :: outgoing
    complex(dp) :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq


    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2

    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)
    !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
    !   This is an expression for (-i)/qsq* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_dk2mom(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_dk2mom(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_dk2mom(3)=-ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_dk2mom(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))


    do j=1,4
       pol_dk2mom(j)=pol_dk2mom(j)/qsq
    enddo

    ! -- do nothing in this case
    if (present(outgoing)) then
       !if (outgoing) pol_dk2mom = conjg(pol_dk2mom)
    endif

  end function pol_dk2mom









!     ubar spinor, massless
  function ubar0(p,i)
  implicit none
    complex(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    complex(dp) :: fc, fc2
    real(dp)    :: p0,px,py,pz,mass


    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    mass=dsqrt(dabs(p0**2-px**2-py**2-pz**2))
    if( mass.lt.1d-4 ) mass=0d0


    fc2 = p0 + pz
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then
       if (i.eq.1) then
          ubar0(1)=czero
          ubar0(2)=czero
          ubar0(3)=fc
          ubar0(4)=(px-ci*py)/fc
       elseif (i.eq.-1) then
          ubar0(1)=(px+ci*py)/fc
          ubar0(2)=-fc
          ubar0(3)=czero
          ubar0(4)=czero
       else
          stop 'ubar0: i out of range'
       endif
    else
       if (i.eq.1) then
          ubar0(1) = czero
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = sqrt(cone*two*p0)
       elseif (i.eq.-1) then
          ubar0(1) = sqrt(cone*(two*p0))
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = czero
       else
          stop 'ubar0: i out of range'
       endif
    endif


!       if (i.eq.1) then 
!           ubar0(1)=dcmplx(mass,0d0)/fc
!           ubar0(2)=czero
!           ubar0(3)=fc
!           ubar0(4)=dcmplx(px,-py)/fc
!       elseif (i.eq.-1) then 
!           ubar0(1)=dcmplx(px,py)/fc
!           ubar0(2)=-fc
!           ubar0(3)=czero
!           ubar0(4)=-dcmplx(mass,0d0)/fc
!        else
!           stop 'ubar0: i out of range'
!       endif 



  end function ubar0



  ! -- v0  spinor, massless
  function v0(p,i)
  implicit none
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    complex(dp) :: v0(4)
    complex(dp) :: fc2, fc
    real(dp)    :: p0,px,py,pz,mass

    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    mass=dsqrt(dabs(p0**2-px**2-py**2-pz**2))
    if( mass.lt.1d-4 ) mass=0d0


    fc2 = p0 + pz
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then
       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=(px-ci*py)/fc
          v0(4)=-fc
       elseif (i.eq.-1) then
          v0(1)=fc
          v0(2)=(px+ci*py)/fc
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif
    else
       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=sqrt(cone*two*p0)
          v0(4)=czero
       elseif (i.eq.-1) then
          v0(1)=czero
          v0(2)=sqrt(cone*two*p0)
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif
    endif


!       if (i.eq.+1) then 
!           v0(1)=czero
!           v0(2)=dcmplx(mass,0d0)/fc
!           v0(3)=dcmplx(px,-py)/fc
!           v0(4)=-fc
!       elseif (i.eq.-1) then
!           v0(1)=fc
!           v0(2)=dcmplx(px,py)/fc
!           v0(3)=dcmplx(-mass,0d0)/fc
!           v0(4)=czero
!        else
!           stop 'v0: i out of range'
!       endif



  end function v0






  ! -- v  spinor, massive (from HELAS)
  FUNCTION vspi(p,pol)
  implicit none
  complex(dp), intent(in) :: p(1:4)
  integer, intent(in)     :: pol
  complex(dp) :: vspi(1:4),chi(1:2)
  real(dp)    :: p0,px,py,pz,pabs,omegaP,omegaM


    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    pabs = sqrt( px**2+py**2+pz**2 )
  
    omegaP = sqrt(abs( p0+pabs ))
    omegaM = sqrt(abs( p0-pabs ))


    if( pol.eq.+1 ) then
        chi(1) =-px + (0.0_dp,1.0_dp)*py ! this is chi-
        chi(2) = pabs + pz
        chi(1:2) = chi(1:2)/sqrt(abs(2.0_dp*pabs*(pabs+pz)))

        vspi(1:2) = -pol * omegaP * chi(1:2)
        vspi(3:4) = +pol * omegaM * chi(1:2)
    elseif( pol.eq.-1 ) then
        chi(1) = pabs + pz ! this is chi+
        chi(2) = px + (0.0_dp,1.0_dp)*py
        chi(1:2) = chi(1:2)/sqrt(abs(2.0_dp*pabs*(pabs+pz)))

        vspi(1:2) = -pol * omegaM * chi(1:2)
        vspi(3:4) = +pol * omegaP * chi(1:2)
    else
        print *,  'vspi: pol out of range'
        stop
    endif
    

  RETURN
  END FUNCTION vspi


  ! -- u  spinor, massive (from HELAS)
  FUNCTION uspi(p,pol)
  implicit none
  complex(dp), intent(in) :: p(1:4)
  integer, intent(in)     :: pol
  complex(dp) :: uspi(1:4),chi(1:2)
  real(dp)    :: p0,px,py,pz,pabs,omegaP,omegaM


    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    pabs = sqrt( px**2+py**2+pz**2 )
  
    omegaP = sqrt(abs( p0+pabs ))
    omegaM = sqrt(abs( p0-pabs ))

    if( pol.eq.+1 ) then
        chi(1) = pabs + pz ! this is chi+
        chi(2) = px + (0.0_dp,1.0_dp)*py
        chi(1:2) = chi(1:2)/sqrt(abs(2.0_dp*pabs*(pabs+pz)))

        uspi(1:2) = omegaM * chi(1:2)
        uspi(3:4) = omegaP * chi(1:2)
    elseif( pol.eq.-1 ) then
        chi(1) =-px + (0.0_dp,1.0_dp)*py ! this is chi-
        chi(2) = pabs + pz
        chi(1:2) = chi(1:2)/sqrt(abs(2.0_dp*pabs*(pabs+pz)))

        uspi(1:2) = omegaP * chi(1:2)
        uspi(3:4) = omegaM * chi(1:2)
    else
        print *,  'uspi: pol out of range'
        stop
    endif
    

  RETURN
  END FUNCTION uspi


  ! -- ubar  spinor, massive (from HELAS)
  FUNCTION ubarspi(p,pol)
  implicit none
  complex(dp), intent(in) :: p(1:4)
  integer, intent(in)     :: pol
  complex(dp) :: uspi_tmp(1:4),ubarspi(1:4)

    uspi_tmp(1:4) = uspi(p,pol)
    ubarspi(1) = dconjg(uspi_tmp(3))
    ubarspi(2) = dconjg(uspi_tmp(4))
    ubarspi(3) = dconjg(uspi_tmp(1))
    ubarspi(4) = dconjg(uspi_tmp(2))

  RETURN
  END FUNCTION ubarspi




          subroutine vSpiDIRAC(p,m,i,f)!   Dirac spinor, massive  (from TOPAZ)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8) p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2 = p0+m
          fc=cdsqrt(dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=pz*fc/fc2
            f(2)=(px+(0d0,1d0)*py)*fc/fc2
            f(3)=fc
            f(4)=dcmplx(0d0,0d0)
          elseif (i.eq.-1) then
            f(1)=(px-(0d0,1d0)*py)*fc/fc2
            f(2)=-pz*fc/fc2
            f(3)=dcmplx(0d0,0d0)
            f(4)=fc
          else
              print *, "wrong helicity setting in vspi"
              stop
          endif

          return
          end SUBROUTINE



          subroutine ubarSpiDIRAC(p,m,i,f)!   Dirac spinor, massive  (from TOPAZ)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8)  p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2=p0+m
          fc=cdsqrt( dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=fc
            f(2)=dcmplx(0d0,0d0)
            f(3)=-1d0*pz*fc/fc2
            f(4)=-(px-(0d0,1d0)*py)*fc/fc2
          elseif (i.eq.-1) then
            f(1)=dcmplx(0d0,0d0)
            f(2)=fc
            f(3)=-(px+(0d0,1d0)*py)*fc/fc2
            f(4)=pz*fc/fc2
          else
              print *, "wrong helicity setting in ubarSpi"
              stop
          endif

          return
          end subroutine





      function vbqq(sp1,sp2)!  (from TOPAZ)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i
      complex(8) :: vbqq(4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)      
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)      
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)      
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

          do i=1,4
             call spb2(sp1,dcmplx(va(i,1:4)),sp1a)
             vbqq(i) = sp1a(1)*sp2(1)+sp1a(2)*sp2(2)+sp1a(3)*sp2(3)+sp1a(4)*sp2(4)
          enddo

      end function vbqq



         subroutine spb2(sp,v,f)!  (from TOPAZ)
         implicit none
         integer i,i1,i2,i3,Dv,Ds,imax
         double complex sp(4),v(4),f(4)
         double complex x0(4,4),xx(4,4),xy(4,4)
         double complex xz(4,4),x5(4,4)
         double complex y1,y2,y3,y4,bp,bm,cp,cm


           Ds=4
           Dv=4
           imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y1
           x0(2,i)=y2
           x0(3,i)=-y3
           x0(4,i)=-y4

           xx(1,i) = -y4
           xx(2,i) = -y3
           xx(3,i) = y2
           xx(4,i) = y1

           xy(1,i)=dcmplx(0d0,-1d0)*y4
           xy(2,i)=dcmplx(0d0,1d0)*y3
           xy(3,i)=dcmplx(0d0,1d0)*y2
           xy(4,i)=dcmplx(0d0,-1d0)*y1

           xz(1,i)=-y3
           xz(2,i)=y4
           xz(3,i)=y1
           xz(4,i)=-y2

           x5(1,i)=y3
           x5(2,i)=y4
           x5(3,i)=y1
           x5(4,i)=y2

           enddo

           if (Dv.eq.4) then

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo

           endif

           return
           end SUBROUTINE






       end module


