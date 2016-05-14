      module modZprime
      use ModParameters
      use ModMisc
      implicit none
      private


!----- notation for subroutines
      public :: EvalAmp_qqb_Zprime_VV,EvalAmp_Zprime_VV
      public :: calcHelAmp,calcHelAmp2

      contains

!----- a subroutine for q qbar -> Zprime -> ZZ/WW
!----- all outgoing convention and the following momentum assignment
!-----  0 -> bq(p1) + q(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
     subroutine EvalAmp_qqb_Zprime_VV(p,MY_IDUP,sum)
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      real(dp) :: pin(4,4)
      complex(dp) :: A(2),qL,qR
      integer :: i1,i2,i3,i4,ordering(1:4),ordering_swap(1:4),MY_IDUP_ordered(6:9),j
      real(dp) :: aL1,aR1,aL2,aR2
      real(dp) :: gZ_sq
      real(dp) :: prefactor
      real(dp) :: intcolfac

         if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
            intcolfac=1.0_dp/3.0_dp
         else
            intcolfac=1.0_dp
         endif


      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
!---- chiral couplings of quarks to Zprimes
      qL = zprime_qq_left
      qR = zprime_qq_right

!---- full prefactor; 3 is  the color factor
      prefactor = 3d0*gZ_sq**2

      ordering=(/3,4,5,6/)
      if(convertLHE(MY_IDUP(6)).lt.0 .or. MY_IDUP(6).eq.Not_a_particle_) then
         call swap(ordering(1),ordering(2))
      endif
      if(convertLHE(MY_IDUP(8)).lt.0 .or. MY_IDUP(8).eq.Not_a_particle_) then
         call swap(ordering(3),ordering(4))
      endif
      if( &
         (CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_ .and. CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_) .or. &
         (CoupledVertex(MY_IDUP(8:9),-1).eq.Z0_ .and. (MY_IDUP(6).eq.Pho_ .or. MY_IDUP(7).eq.Pho_)) &
      ) then
         call swap(ordering(1),ordering(3))
         call swap(ordering(2),ordering(4))
      endif
      ordering_swap(:)=ordering(:)
      call swap(ordering_swap(1),ordering_swap(3))

      do j=1,4
         MY_IDUP_ordered(j+5)=MY_IDUP(ordering(j)+3)
      enddo

         if( CoupledVertex(MY_IDUP_ordered(6:7),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP_ordered(6)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(6)).eq.abs(MuM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(NuE_) .or. abs(MY_IDUP_ordered(6)).eq.abs(NuM_) .or. abs(MY_IDUP_ordered(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(Up_) .or. abs(MY_IDUP_ordered(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(Dn_) .or. abs(MY_IDUP_ordered(6)).eq.abs(Str_) .or. abs(MY_IDUP_ordered(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL1=0d0
                    aR1=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP_ordered(6:7),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP_ordered(6:7),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP_ordered(6)) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ud)
                 aR1 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(6)).eq.abs(MuM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ln)
                 aR1 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(TaM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_tn)
                 aR1 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL1=0d0
                 aR1=0d0
              endif
         elseif( MY_IDUP_ordered(6).eq.Pho_ ) then !  photon
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL1=0d0
              aR1=0d0
         endif

         if( CoupledVertex(MY_IDUP_ordered(8:9),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP_ordered(8)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(8)).eq.abs(MuM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(NuE_) .or. abs(MY_IDUP_ordered(8)).eq.abs(NuM_) .or. abs(MY_IDUP_ordered(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(Up_) .or. abs(MY_IDUP_ordered(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(Dn_) .or. abs(MY_IDUP_ordered(8)).eq.abs(Str_) .or. abs(MY_IDUP_ordered(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP_ordered(8:9),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP_ordered(8:9),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP_ordered(8)) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ud)
                 aR2 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP_ordered(9)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(9)).eq.abs(MuM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ln)
                 aR2 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP_ordered(9)).eq.abs(TaM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_tn)
                 aR2 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL2=0d0
                 aR2=0d0
              endif
         elseif( MY_IDUP_ordered(8).eq.Pho_ ) then !  photon
              aL2=1d0
              aR2=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL2=0d0
              aR2=0d0
         endif




      sum = zero
do i1=1,2
do i3 = 1,2
do i4 = 1,2
!do i3 = -1,1! on-shell check!
!do i4 = -1,1! on-shell check!

         call calcHelAmp(ordering,p(1:4,1:6),i1,i3,i4,A(1))
         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call calcHelAmp(ordering_swap,p(1:4,1:6),i1,i3,i4,A(2))
             A(2) = -A(2) ! minus comes from fermi statistics
         endif


         if (i1.eq.1) then
            A(:) = qL*A(:)
         elseif(i1.eq.2) then
            A(:) = qR*A(:)
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
             sum = sum + SymmFac * (cdabs( A(1)*dconjg(A(1)) ) + cdabs( A(2)*dconjg(A(2)) ))
             if( i3.eq.i4 ) sum = sum + SymmFac * 2d0*intcolfac*dreal(A(1)*dconjg(A(2)))
         else
             sum = sum + cdabs( A(1)*dconjg(A(1)) )
         endif


enddo
enddo
enddo

!  sum ~ GeV~-8
!  prefactor ~ GeV^0
         sum = sum*prefactor

      end subroutine

     subroutine calcHelAmp(ordering,p,i1,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4
     real(dp) :: p(1:4,1:6)
     complex(dp) :: propG, propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(4,4)


      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      !print *,"p(6)=",(p(:,l1))
      !print *,"p(7)=",(p(:,l2))
      !print *,"p(8)=",(p(:,l3))
      !print *,"p(9)=",(p(:,l4))
      !pause

      s  = two*scr(p(:,1),p(:,2))
      propG = s/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)


!       s = two*scr(p(:,3),p(:,4))
!       propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
!       s = two*scr(p(:,5),p(:,6))
!       propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)


!-------- For fermions, -1 == left, 1 == right
         pin(1,:) = p(:,1)
         pin(2,:) = p(:,2)
         sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1) !qbq

         sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform bookeeping
!          pin(3,:) = p(:,3) + p(:,4)
!          pin(4,:) = p(:,5) + p(:,6)
!          sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
!          sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-.mu+

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


!-- on-shell check
!       sp(3,1:4) = pol_mass(dcmplx(p(1:4,3)+p(1:4,4)),dsqrt(2d0*scr(p(:,3),p(:,4))),i3)
!       sp(4,1:4) = pol_mass(dcmplx(p(1:4,5)+p(1:4,6)),dsqrt(2d0*scr(p(:,5),p(:,6))),i4)

!        if( i3.eq.0 .and. i4.eq.0 ) then
!             a(:) = 0d0
!             return
!        endif

!-- end: on-shell check



         call qqZprimeZZampl(pin,sp,A(1))

! A(1) ~ GeV^-2
! propG*propZ1*propZ2 ~  GeV^-2
         A(1) = A(1) * propG*propZ1*propZ2

      end subroutine

      subroutine qqZprimeZZampl(p,sp,res)
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
      complex(dp) :: yyy1,yyy2,yyy3,yyy4,xxx1,epsZpr(1:4,-1:+1)

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


      q1_q2 = sc(q1,q2)
      q1_q3 = sc(q1,q3)
      q1_q4 = sc(q1,q4)
      q2_q3 = sc(q2,q3)
      q2_q4 = sc(q2,q4)
      q3_q4 = sc(q3,q4)
      e1_e2 = sc(e1,e2)

      e1_e3 = sc(e1,e3)
      e1_e4 = sc(e1,e4)






!       epsZpr(1:4,-1) = pol_mass(q,m_reso,-1)
!       epsZpr(1:4, 0) = pol_mass(q,m_reso, 0)
!       epsZpr(1:4,+1) = pol_mass(q,m_reso,+1)
!       e1_e3 = sc(e1,epsZpr(1:4,-1)) * sc(e3,dconjg(epsZpr(1:4,-1)))  &
!             + sc(e1,epsZpr(1:4,+1)) * sc(e3,dconjg(epsZpr(1:4,+1)))  &
!             + sc(e1,epsZpr(1:4, 0)) * sc(e3,dconjg(epsZpr(1:4, 0)))
!       e1_e4 = sc(e1,epsZpr(1:4,-1)) * sc(e4,dconjg(epsZpr(1:4,-1)))  &
!             + sc(e1,epsZpr(1:4,+1)) * sc(e4,dconjg(epsZpr(1:4,+1)))  &
!             + sc(e1,epsZpr(1:4, 0)) * sc(e4,dconjg(epsZpr(1:4, 0)))




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


      res = czero

      xxx1 = (1d0,0d0)  !  different possibilities for fermion couplings
                        !  to zprime are accounted in the amplitude call


      yyy1 = zprime_zz_1
      yyy2 = zprime_zz_2


       res= - e1_e3*e4_q3*xxx1*yyy1    &
            - e1_e4*e3_q4*xxx1*yyy1    &
            - et1(e1,e3,e4,q3)*xxx1*yyy2 + et1(e1,e3,e4,q4)*xxx1*yyy2


!       e3 = pol_mass(q3,dreal(cdsqrt(sc(q3,q3))),-1)
!       e4 = pol_mass(q4,dreal(cdsqrt(sc(q4,q4))),-1)
!       e1_e3 = sc(e1,e3)
!       e1_e4 = sc(e1,e4)
!       e3_q4 = sc(e3,q4)
!       e4_q3 = sc(e4,q3)
!print *, dreal(cdsqrt(sc(q3,q3))),dreal(cdsqrt(sc(q4,q4)))
!print * , "e3x",q3,dreal(cdsqrt(sc(q3,q3))),-1
!print *, - e1_e3*e4_q3 - e1_e4*e3_q4
!pause



      end subroutine qqZprimeZZampl



!----- a subroutine for Zprime -> ZZ/WW
!----- all outgoing convention and the following momentum assignment
!-----  0 --> Zprime(p1) e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
     subroutine EvalAmp_Zprime_VV(p,M_Reso,Ga_Reso,vvcoupl,MY_IDUP,sum)
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6),M_Reso,Ga_Reso
      complex(dp) :: vvcoupl(1:2)
      integer, intent(in) :: MY_IDUP(6:9)
      real(dp) :: pin(4,4)
      complex(dp) :: A(2)
      integer :: i1,i2,i3,i4,ordering(1:4),ordering_swap(1:4),MY_IDUP_ordered(6:9),j
      real(dp) :: aL1,aR1,aL2,aR2
      real(dp) :: gZ_sq
      real(dp) :: prefactor
      real(dp) :: intcolfac

         if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
            intcolfac=1.0_dp/3.0_dp
         else
            intcolfac=1.0_dp
         endif


      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
!---- full prefactor
      prefactor = gZ_sq**2


      ordering=(/3,4,5,6/)
      if(convertLHE(MY_IDUP(6)).lt.0 .or. MY_IDUP(6).eq.Not_a_particle_) then
         call swap(ordering(1),ordering(2))
      endif
      if(convertLHE(MY_IDUP(8)).lt.0 .or. MY_IDUP(8).eq.Not_a_particle_) then
         call swap(ordering(3),ordering(4))
      endif
      if( &
         (CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_ .and. CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_) .or. &
         (CoupledVertex(MY_IDUP(8:9),-1).eq.Z0_ .and. (MY_IDUP(6).eq.Pho_ .or. MY_IDUP(7).eq.Pho_)) &
      ) then
         call swap(ordering(1),ordering(3))
         call swap(ordering(2),ordering(4))
      endif
      ordering_swap(:)=ordering(:)
      call swap(ordering_swap(1),ordering_swap(3))

      do j=1,4
         MY_IDUP_ordered(j+5)=MY_IDUP(ordering(j)+3)
      enddo

         if( CoupledVertex(MY_IDUP_ordered(6:7),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP_ordered(6)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(6)).eq.abs(MuM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(NuE_) .or. abs(MY_IDUP_ordered(6)).eq.abs(NuM_) .or. abs(MY_IDUP_ordered(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(Up_) .or. abs(MY_IDUP_ordered(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(Dn_) .or. abs(MY_IDUP_ordered(6)).eq.abs(Str_) .or. abs(MY_IDUP_ordered(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL1=0d0
                    aR1=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP_ordered(6:7),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP_ordered(6:7),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP_ordered(6)) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ud)
                 aR1 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(6)).eq.abs(MuM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ln)
                 aR1 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP_ordered(6)).eq.abs(TaM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_tn)
                 aR1 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL1=0d0
                 aR1=0d0
              endif
         elseif( MY_IDUP_ordered(6).eq.Pho_ ) then !  photon
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL1=0d0
              aR1=0d0
         endif

         if( CoupledVertex(MY_IDUP_ordered(8:9),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP_ordered(8)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(8)).eq.abs(MuM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(NuE_) .or. abs(MY_IDUP_ordered(8)).eq.abs(NuM_) .or. abs(MY_IDUP_ordered(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(Up_) .or. abs(MY_IDUP_ordered(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP_ordered(8)).eq.abs(Dn_) .or. abs(MY_IDUP_ordered(8)).eq.abs(Str_) .or. abs(MY_IDUP_ordered(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP_ordered(8:9),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP_ordered(8:9),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP_ordered(8)) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ud)
                 aR2 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP_ordered(9)).eq.abs(ElM_) .or. abs(MY_IDUP_ordered(9)).eq.abs(MuM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ln)
                 aR2 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP_ordered(9)).eq.abs(TaM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_tn)
                 aR2 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL2=0d0
                 aR2=0d0
              endif
         elseif( MY_IDUP_ordered(8).eq.Pho_ ) then !  photon
              aL2=1d0
              aR2=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL2=0d0
              aR2=0d0
         endif




      sum = zero
do i1 =-1,1! Z' boson
do i3 = 1,2! lepton string1
do i4 = 1,2! lepton string2

         call calcHelAmp2(ordering,p(1:4,1:6),i1,i3,i4,A(1))
         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call calcHelAmp2(ordering_swap,p(1:4,1:6),i1,i3,i4,A(2))
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
             sum = sum + SymmFac * (cdabs( A(1)*dconjg(A(1)) ) + cdabs( A(2)*dconjg(A(2)) ))
             if( i3.eq.i4 ) sum = sum + SymmFac * 2d0*intcolfac*dreal(A(1)*dconjg(A(2)))
         else
             sum = sum + cdabs( A(1)*dconjg(A(1)) )
         endif


enddo
enddo
enddo

         sum = sum*prefactor

      end subroutine

     subroutine calcHelAmp2(ordering,p,i1,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4
     real(dp) :: p(1:4,1:6)
     complex(dp) :: propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(4,4)



      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)


      s  = scr(p(:,1),p(:,1))


         pin(1,:) = p(:,1)
         sp(1,:) = pol_mass2(dcmplx(p(1:4,1)),i1,'in')
         pin(2,:) = 0d0  ! dummy
         sp(2,:)  = 0d0

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



         call ZprimeZZampl(pin,sp,A(1))

         A(1) = A(1) * propZ1*propZ2

      end subroutine

      subroutine ZprimeZZampl(p,sp,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q_q,MZ3,MZ4
      complex(dp) :: q1_q2,q1_q3,q1_q4
      complex(dp) :: q2_q3,q2_q4
      complex(dp) :: q3_q4,q3_q3,q4_q4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4)
      complex(dp) :: yyy1,yyy2,yyy3,yyy4,epsZpr(1:4,-1:+1)


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


      q1_q2 = sc(q1,q2)
      q3_q3 = sc(q3,q3)
      q4_q4 = sc(q4,q4)
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


      MZ3=dsqrt(cdabs(q3_q3))
      MZ4=dsqrt(cdabs(q4_q4))


      res = czero


      yyy1 = zprime_zz_1
      yyy2 = zprime_zz_2


       res= yyy1*( q1_e3*e1_e4 + q1_e4*e1_e3 )  +   yyy2*( et1(e1,e3,e4,q3) - et1(e1,e3,e4,q4) )


      end subroutine ZprimeZZampl


       end module


