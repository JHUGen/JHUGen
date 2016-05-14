      module modGraviton
      use ModParameters
      use ModMisc
      implicit none
      private


!----- notation for subroutines
      public :: EvalAmp_gg_G_VV,EvalAmp_qqb_G_VV,EvalAmp_G_VV

      contains


!----- a subroutinefor gg -> G -> ZZ/WW
!----- all outgoing convention and the following momentum assignment
!-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_gg_G_VV(p,MY_IDUP,sum)
      use ModMisc
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      real(dp) :: pin(4,4)
      complex(dp) :: A(2)
      integer :: i1,i2,i3,i4,ordering(1:4)
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


!---- full prefactor; 8 is  the color factor
      prefactor = 8d0*gZ_sq**2


         if( CoupledVertex(MY_IDUP(6:7),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL1=0d0
                    aR1=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP(6:7),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP(6)) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ud)
                 aR1 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ln)
                 aR1 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_tn)
                 aR1 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL1=0d0
                 aR1=0d0
              endif
         elseif( MY_IDUP(6).eq.Pho_ ) then !  photon
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL1=0d0
              aR1=0d0
         endif

         if( CoupledVertex(MY_IDUP(8:9),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP(8:9),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP(8)) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ud)
                 aR2 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP(9)).eq.abs(ElM_) .or. abs(MY_IDUP(9)).eq.abs(MuM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ln)
                 aR2 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP(9)).eq.abs(TaM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_tn)
                 aR2 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL2=0d0
                 aR2=0d0
              endif
         elseif( MY_IDUP(8).eq.Pho_ ) then !  photon
              aL2=1d0
              aR2=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL2=0d0
              aR2=0d0
         endif


sum = zero
do i1=1,2
do i2 = 1,2
do i3 = 1,2
do i4 = 1,2

         if(CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_ .and. CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_) then
            ordering = (/5,6,3,4/)
         else
            ordering = (/3,4,5,6/)
         endif
         call calcHelAmp_gg(ordering,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A(1))
         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call swap(ordering(2),ordering(4))
             call calcHelAmp_gg(ordering,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A(2))
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
             sum = sum + SymmFac * (cdabs( A(1)*dconjg(A(1)) ) +  cdabs( A(2)*dconjg(A(2)) ))
             if( i3.eq.i4 ) sum = sum + SymmFac * 2d0*intcolfac*dreal(A(1)*dconjg(A(2)))
         else
             sum = sum + cdabs( A(1)*dconjg(A(1)) )
         endif

enddo
enddo
enddo
enddo

      sum = sum*prefactor

      end subroutine

     subroutine calcHelAmp_gg(ordering,p,idin,i1,i2,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i1,i2,i3,i4,l1,l2,l3,l4,idin(3:6)
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
          !sp(1,1:4)=pin(1,1:4);print *, "this checks IS gauge invariance"
          !sp(2,1:4)=pin(2,1:4);print *, "this checks IS gauge invariance"
          !careful: for gluon gauge invariance check the terms ~c3,c4 are needed because e1.q2 is not zero for e1-->q1

!-------- -1 == left, 1 == right
         if( idin(l1).ne.Pho_ ) then ! Z, W
            pin(3,:) = p(:,l1)+p(:,l2)
            sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
            sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
            s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
            propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         else !  A
            pin(3,:) = p(:,l1)
            sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
            propZ1=1d0
         endif
         if( idin(l3).ne.Pho_ ) then ! Z, W
            pin(4,:) = p(:,l3)+p(:,l4)
            sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
            sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
            s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
            propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         else ! A
            pin(4,:) = p(:,l3)
            sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')
            propZ2=1d0
         endif

         call ggGZZampl(pin,sp,A(1))

         A(1) = A(1) * propG*propZ1*propZ2

      end subroutine

      subroutine ggGZZampl(p,sp,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q_q
      complex(dp) :: q1_q2,q1_q3,q1_q4,q3_q3,q4_q4
      complex(dp) :: q2_q3,q2_q4
      complex(dp) :: q3_q4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4,   q1_e1,q1_e2,q2_e1,q2_e2,q4_e1,q3_e1
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4)
      complex(dp) :: xxx1,xxx2,xxx3,xxx4,yyy1,yyy2,yyy3,yyy4,yyy41,yyy42,yyy5,yyy6
      complex(dp) :: yyy7,abr1
      real(dp) :: q34,MZ3,MZ4,MG
      logical :: new
      real(dp) :: rr_gam, rr


      new = .true.

      q1 = dcmplx(p(1,:),0d0)
      q2 = dcmplx(p(2,:),0d0)
      q3 = dcmplx(p(3,:),0d0)
      q4 = dcmplx(p(4,:),0d0)


      e1 = sp(1,:)
      e2 = sp(2,:)
      e3 = sp(3,:)
      e4 = sp(4,:)

      q = -q1-q2
      q_q = sc(q,q)
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

      MZ3=dsqrt(cdabs(q3_q3))
      MZ4=dsqrt(cdabs(q4_q4))
      if( use_dynamic_MG ) then
          MG = dsqrt(cdabs(q_q))
      else
          MG = M_Reso
      endif

!---- define couplings
      q34 = (MG**2-MZ3**2-MZ4**2)/2d0!  = s = pV1.pV2    = (q_q-MZ3^2-MZ4^2)/2

      rr_gam = q_q/two/Lambda**2! kappa for IS
      xxx1 = (a1 + a2*rr_gam)     ! those a's correspond to g's in eq.(7) for the IS
      xxx2 = -a1/two + (a3+two*a4)*rr_gam
      xxx3 = 4d0*a5*rr_gam * MG**2/q_q
      ! for gluon gauge invariance check the terms ~c3,c4 are needed because e1.q2 is not zero for e1-->q1

      if (generate_bis) then
          rr = q34/Lambda**2! kappa for FS
! !           yyy1 = q34*(b1 + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)   + b5*M_V**2/q34)
          yyy1 = q34*(b1 + b2*rr*(one+MZ3**2/q34)*(one+MZ4**2/q34)  + b5*M_V**2/q34)
!           yyy2 = -b1/two + b3*rr*(one-M_V**2/q34) + two*b4*rr + b7*rr*M_V**2/q34
          yyy2 = -b1/two + b3*rr*(1d0-(MZ3**2+MZ4**2)/(2d0*q34)) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
!           yyy4 = -b1 - b2*rr -(b2+b3+2d0*b6)*rr*M_V**2/q34
          yyy41 = -b1 - b2*(q34+MZ3**2)/Lambda**2 - b3*MZ4**2/Lambda**2 - 2d0*b6*M_V**2/Lambda**2
          yyy42 = -b1 - b2*(q34+MZ4**2)/Lambda**2 - b3*MZ3**2/Lambda**2 - 2d0*b6*M_V**2/Lambda**2
          yyy5 = two*b8*rr*MG**2/q34
!           yyy6 = b9
          yyy6 = b9 * M_V**2/Lambda**2
!           yyy7 = b10*rr*MG**2/q34
          yyy7 = b10 * MG**2 * M_V**2/Lambda**4
      else
          yyy1 = q34*c1/2d0
          yyy2 = c2
          yyy3 = c3/MG**2
!           yyy4 = c4
          yyy41 = c41
          yyy42 = c42
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
          if( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
              yyy6=0d0
              yyy7=0d0
          endif
      endif

      res = czero


     if (new) then


!       res = &
!         + 8.D0*q1_e3*q1_e4*e1_e2*yyy1*xxx2 - 8.D0*q1_e3*q1_q3*e1_e2* &
!           e4_q3*yyy4*xxx2 + 4.D0*q1_e3*e1_e2*e4_q3*yyy1*xxx2 + 2.D0* &
!           q1_e3*e1_e2*e4_q3*MZ4**2*yyy4*xxx2 - 2.D0*q1_e3*e1_e2*e4_q3* &
!           MZ3**2*yyy4*xxx2 - 2.D0*q1_e3*e1_e2*e4_q3*MG**2*yyy4*xxx2 + 8.d0 &
!           *q1_e4*q1_q3*e1_e2*e3_q4*yyy4*xxx2 + 4.D0*q1_e4*e1_e2*e3_q4 &
!           *yyy1*xxx2 - 2.D0*q1_e4*e1_e2*e3_q4*MZ4**2*yyy4*xxx2 + 2.D0* &
!           q1_e4*e1_e2*e3_q4*MZ3**2*yyy4*xxx2 + 2.D0*q1_e4*e1_e2*e3_q4* &
!           MG**2*yyy4*xxx2 - 8.D0*q1_q3*e1_e2*e3_e4*MZ4**2*yyy2*xxx2 + 8.d0 &
!           *q1_q3*e1_e2*e3_e4*MZ3**2*yyy2*xxx2 + 8.D0*q1_q3*e1_e2* &
!           e3_e4*MG**2*yyy2*xxx2 - 8.D0*q1_q3*e1_e2*e3_q4*e4_q3*MZ4**2* &
!           yyy3*xxx2 + 8.D0*q1_q3*e1_e2*e3_q4*e4_q3*MZ3**2*yyy3*xxx2 + 8.d0 &
!           *q1_q3*e1_e2*e3_q4*e4_q3*MG**2*yyy3*xxx2 + 16.D0*q1_q3**2* &
!           e1_e2*e3_e4*yyy2*xxx2 + 16.D0*q1_q3**2*e1_e2*e3_q4*e4_q3*yyy3 &
!           *xxx2 + 2.D0/3.D0*e1_e2*e3_e4*MZ4**4*yyy2*xxx2 + 1.D0/3.D0* &
!           e1_e2*e3_e4*MZ4**4*yyy2*xxx1
!       res = res - 4.D0/3.D0*e1_e2*e3_e4*MZ3**2*MZ4**2*yyy2*xxx2 - 2.D0/ &
!           3.D0*e1_e2*e3_e4*MZ3**2*MZ4**2*yyy2*xxx1 + 2.D0/3.D0*e1_e2* &
!           e3_e4*MZ3**4*yyy2*xxx2 + 1.D0/3.D0*e1_e2*e3_e4*MZ3**4*yyy2* &
!           xxx1 + 2.D0/3.D0*e1_e2*e3_e4*MG**2*yyy1*xxx2 - 2.D0/3.D0* &
!           e1_e2*e3_e4*MG**2*yyy1*xxx1 - 4.D0/3.D0*e1_e2*e3_e4*MG**2* &
!           MZ4**2*yyy2*xxx2 - 2.D0/3.D0*e1_e2*e3_e4*MG**2*MZ4**2*yyy2* &
!           xxx1 + 8.D0/3.D0*e1_e2*e3_e4*MG**2*MZ3**2*yyy2*xxx2 - 2.D0/3.D0 &
!           *e1_e2*e3_e4*MG**2*MZ3**2*yyy2*xxx1 + 2.D0/3.D0*e1_e2*e3_e4* &
!           MG**4*yyy2*xxx2 + 1.D0/3.D0*e1_e2*e3_e4*MG**4*yyy2*xxx1 + 4.D0 &
!           /3.D0*e1_e2*e3_q4*e4_q3*yyy1*xxx2 + 2.D0/3.D0*e1_e2*e3_q4* &
!           e4_q3*yyy1*xxx1 + 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MZ4**4*yyy3* &
!           xxx2 + 1.D0/3.D0*e1_e2*e3_q4*e4_q3*MZ4**4*yyy3*xxx1 - 4.D0/3.D0 &
!           *e1_e2*e3_q4*e4_q3*MZ3**2*MZ4**2*yyy3*xxx2 - 2.D0/3.D0*e1_e2 &
!           *e3_q4*e4_q3*MZ3**2*MZ4**2*yyy3*xxx1 + 2.D0/3.D0*e1_e2*e3_q4* &
!           e4_q3*MZ3**4*yyy3*xxx2 + 1.D0/3.D0*e1_e2*e3_q4*e4_q3*MZ3**4* &
!           yyy3*xxx1
!       res = res + 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**2*yyy4*xxx2 - 2.D0/3.D0 &
!           *e1_e2*e3_q4*e4_q3*MG**2*yyy4*xxx1 - 4.D0/3.D0*e1_e2*e3_q4* &
!           e4_q3*MG**2*MZ4**2*yyy3*xxx2 - 2.D0/3.D0*e1_e2*e3_q4*e4_q3* &
!           MG**2*MZ4**2*yyy3*xxx1 + 8.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**2* &
!           MZ3**2*yyy3*xxx2 - 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**2*MZ3**2* &
!           yyy3*xxx1 + 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**4*yyy3*xxx2 + 1.D0 &
!           /3.D0*e1_e2*e3_q4*e4_q3*MG**4*yyy3*xxx1 + e1_e3*e2_e4*MG**2* &
!           yyy1*xxx1 - e1_e3*e2_q3*e4_q3*MG**2*yyy4*xxx1 + e1_e4*e2_e3* &
!           MG**2*yyy1*xxx1 + e1_e4*e2_q3*e3_q4*MG**2*yyy4*xxx1 - e1_q3* &
!           e2_e3*e4_q3*MG**2*yyy4*xxx1 + e1_q3*e2_e4*e3_q4*MG**2*yyy4* &
!           xxx1 + 4.D0*e1_q3*e2_q3*e3_e4*MG**2*yyy2*xxx1 + 4.D0*e1_q3* &
!           e2_q3*e3_q4*e4_q3*MG**2*yyy3*xxx1 + 8.D0*et1(q1,q2,e1,e2)* &
!           q1_e3*q1_e4*MG**(-2)*yyy1*xxx3 - 8.D0*et1(q1,q2,e1,e2)*q1_e3* &
!           q1_q3*e4_q3*MG**(-2)*yyy4*xxx3 + 4.D0*et1(q1,q2,e1,e2)*q1_e3* &
!           e4_q3*MG**(-2)*yyy1*xxx3 + 2.D0*et1(q1,q2,e1,e2)*q1_e3*e4_q3* &
!           MG**(-2)*MZ4**2*yyy4*xxx3
!       res = res - 2.D0*et1(q1,q2,e1,e2)*q1_e3*e4_q3*MG**(-2)*MZ3**2* &
!        yyy4*xxx3 - 2.D0*et1(q1,q2,e1,e2)*q1_e3*e4_q3*yyy4*xxx3 + 8.D0* &
!           et1(q1,q2,e1,e2)*q1_e4*q1_q3*e3_q4*MG**(-2)*yyy4*xxx3 + 4.D0* &
!           et1(q1,q2,e1,e2)*q1_e4*e3_q4*MG**(-2)*yyy1*xxx3 - 2.D0*et1(q1 &
!           ,q2,e1,e2)*q1_e4*e3_q4*MG**(-2)*MZ4**2*yyy4*xxx3 + 2.D0*et1( &
!           q1,q2,e1,e2)*q1_e4*e3_q4*MG**(-2)*MZ3**2*yyy4*xxx3 + 2.D0* &
!           et1(q1,q2,e1,e2)*q1_e4*e3_q4*yyy4*xxx3 - 8.D0*et1(q1,q2,e1,e2 &
!           )*q1_q3*e3_e4*MG**(-2)*MZ4**2*yyy2*xxx3 + 8.D0*et1(q1,q2,e1, &
!           e2)*q1_q3*e3_e4*MG**(-2)*MZ3**2*yyy2*xxx3 + 8.D0*et1(q1,q2,e1 &
!           ,e2)*q1_q3*e3_e4*yyy2*xxx3 - 8.D0*et1(q1,q2,e1,e2)*q1_q3* &
!           e3_q4*e4_q3*MG**(-2)*MZ4**2*yyy3*xxx3 + 8.D0*et1(q1,q2,e1,e2) &
!           *q1_q3*e3_q4*e4_q3*MG**(-2)*MZ3**2*yyy3*xxx3 + 8.D0*et1(q1,q2 &
!           ,e1,e2)*q1_q3*e3_q4*e4_q3*yyy3*xxx3 + 16.D0*et1(q1,q2,e1,e2)* &
!           q1_q3**2*e3_e4*MG**(-2)*yyy2*xxx3 + 16.D0*et1(q1,q2,e1,e2)* &
!           q1_q3**2*e3_q4*e4_q3*MG**(-2)*yyy3*xxx3 + 2.D0/3.D0*et1(q1,q2 &
!           ,e1,e2)*e3_e4*MG**(-2)*MZ4**4*yyy2*xxx3
!       res = res - 4.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)*MZ3**2* &
!        MZ4**2*yyy2*xxx3 + 2.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)* &
!           MZ3**4*yyy2*xxx3 + 2.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*yyy1*xxx3 &
!            - 4.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*MZ4**2*yyy2*xxx3 + 8.D0/3.D0 &
!           *et1(q1,q2,e1,e2)*e3_e4*MZ3**2*yyy2*xxx3 + 2.D0/3.D0*et1(q1 &
!           ,q2,e1,e2)*e3_e4*MG**2*yyy2*xxx3 + 4.D0/3.D0*et1(q1,q2,e1,e2) &
!           *e3_q4*e4_q3*MG**(-2)*yyy1*xxx3 + 2.D0/3.D0*et1(q1,q2,e1,e2)* &
!           e3_q4*e4_q3*MG**(-2)*MZ4**4*yyy3*xxx3 - 4.D0/3.D0*et1(q1,q2, &
!           e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ3**2*MZ4**2*yyy3*xxx3 + 2.D0/3.D0 &
!           *et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ3**4*yyy3*xxx3 + 2.D0 &
!           /3.D0*et1(q1,q2,e1,e2)*e3_q4*e4_q3*yyy4*xxx3 - 4.D0/3.D0* &
!           et1(q1,q2,e1,e2)*e3_q4*e4_q3*MZ4**2*yyy3*xxx3 + 8.D0/3.D0* &
!           et1(q1,q2,e1,e2)*e3_q4*e4_q3*MZ3**2*yyy3*xxx3 + 2.D0/3.D0* &
!           et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**2*yyy3*xxx3 - et1(q1,q2,e1, &
!           e2)*et1(q1,q,e3,e4)*MG**(-2)*MZ4**2*yyy6*xxx3 + et1(q1,q2,e1, &
!           e2)*et1(q1,q,e3,e4)*MG**(-2)*MZ3**2*yyy6*xxx3
!       res = res + et1(q1,q2,e1,e2)*et1(q1,q,e3,e4)*yyy6*xxx3 + 4.D0* &
!           et1(q1,q2,e1,e2)*et1(q1,q,e3,e4)*q1_q3*MG**(-2)*yyy6*xxx3 - 4.0d0 &
!           *et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*q1_q3*e4_q3*MG**(-4)*yyy7 &
!           *xxx3 + et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*e4_q3*MG**(-4)* &
!           MZ4**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*e4_q3* &
!           MG**(-4)*MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)* &
!           e4_q3*MG**(-2)*yyy7*xxx3 + 4.D0*et1(q1,q2,e1,e2)*et1(q1,q,e3, &
!           q4)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1, &
!           q,e3,q4)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)* &
!           et1(q1,q,e3,q4)*e4_q3*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2, &
!           e1,e2)*et1(q1,q,e3,q4)*e4_q3*MG**(-2)*yyy7*xxx3 - 4.D0*et1(q1 &
!           ,q2,e1,e2)*et1(q1,q,e4,q3)*q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 + &
!           et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-4)*MZ4**2*yyy7* &
!           xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-4)*MZ3**2 &
!           *yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-2)* &
!           yyy7*xxx3
!       res = res + 4.D0*et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*q1_q3*e3_q4* &
!        MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*e3_q4* &
!           MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)* &
!           e3_q4*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q1,q, &
!           e4,q4)*e3_q4*MG**(-2)*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q2,q, &
!           e3,e4)*MG**(-2)*MZ4**2*yyy6*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q, &
!           e3,e4)*MG**(-2)*MZ3**2*yyy6*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q, &
!           e3,e4)*yyy6*xxx3 - 4.D0*et1(q1,q2,e1,e2)*et1(q2,q,e3,e4)* &
!           q1_q3*MG**(-2)*yyy6*xxx3 + 4.D0*et1(q1,q2,e1,e2)*et1(q2,q,e3, &
!           q3)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2, &
!           q,e3,q3)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)* &
!           et1(q2,q,e3,q3)*e4_q3*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2, &
!           e1,e2)*et1(q2,q,e3,q3)*e4_q3*MG**(-2)*yyy7*xxx3 - 4.D0*et1(q1 &
!           ,q2,e1,e2)*et1(q2,q,e3,q4)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 + &
!           et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*e4_q3*MG**(-4)*MZ4**2*yyy7* &
!           xxx3
!       res = res - et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*e4_q3*MG**(-4)* &
!        MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*e4_q3* &
!           MG**(-2)*yyy7*xxx3 + 4.D0*et1(q1,q2,e1,e2)*et1(q2,q,e4,q3)* &
!           q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e4 &
!           ,q3)*e3_q4*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1( &
!           q2,q,e4,q3)*e3_q4*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2 &
!           )*et1(q2,q,e4,q3)*e3_q4*MG**(-2)*yyy7*xxx3 - 4.D0*et1(q1,q2, &
!           e1,e2)*et1(q2,q,e4,q4)*q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 + et1( &
!           q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-4)*MZ4**2*yyy7*xxx3 &
!            - et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-4)*MZ3**2* &
!           yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-2)* &
!           yyy7*xxx3 - 1.D0/3.D0*et1(q1,q2,e1,e2)*et1(q,e3,e4,q3)*yyy6* &
!           xxx3 + 1.D0/3.D0*et1(q1,q2,e1,e2)*et1(q,e3,e4,q4)*yyy6*xxx3 &
!            + 2.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*MG**(-4)* &
!           MZ4**4*yyy5*xxx3 - 4.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4 &
!           )*MG**(-4)*MZ3**2*MZ4**2*yyy5*xxx3
!       res = res + 2.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*MG**(-4)* &
!        MZ3**4*yyy5*xxx3 - 4.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)* &
!           MG**(-2)*MZ4**2*yyy5*xxx3 + 8.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3 &
!           ,e4,q3,q4)*MG**(-2)*MZ3**2*yyy5*xxx3 + 2.D0/3.D0*et1(q1,q2,e1 &
!           ,e2)*et1(e3,e4,q3,q4)*yyy5*xxx3 - 8.D0*et1(q1,q2,e1,e2)*et1( &
!           e3,e4,q3,q4)*q1_q3*MG**(-4)*MZ4**2*yyy5*xxx3 + 8.D0*et1(q1,q2 &
!           ,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*MG**(-4)*MZ3**2*yyy5*xxx3 + 8.D0 &
!           *et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*MG**(-2)*yyy5*xxx3 &
!            + 16.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3**2*MG**(-4)* &
!           yyy5*xxx3 + 4.D0*et1(q1,q,e3,e4)*q1_q3*e1_e2*yyy6*xxx2 - et1( &
!           q1,q,e3,e4)*e1_e2*MZ4**2*yyy6*xxx2 + et1(q1,q,e3,e4)*e1_e2* &
!           MZ3**2*yyy6*xxx2 + et1(q1,q,e3,e4)*e1_e2*MG**2*yyy6*xxx2 - 4.D0 &
!           *et1(q1,q,e3,q3)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7*xxx2 + et1( &
!           q1,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q1,q, &
!           e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q1,q,e3,q3 &
!           )*e1_e2*e4_q3*yyy7*xxx2
!       res = res + 4.D0*et1(q1,q,e3,q4)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7* &
!        xxx2 - et1(q1,q,e3,q4)*e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 + &
!           et1(q1,q,e3,q4)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 + et1( &
!           q1,q,e3,q4)*e1_e2*e4_q3*yyy7*xxx2 - 4.D0*et1(q1,q,e4,q3)* &
!           q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7*xxx2 + et1(q1,q,e4,q3)*e1_e2* &
!           e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q1,q,e4,q3)*e1_e2*e3_q4 &
!           *MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q1,q,e4,q3)*e1_e2*e3_q4*yyy7 &
!           *xxx2 + 4.D0*et1(q1,q,e4,q4)*q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7* &
!           xxx2 - et1(q1,q,e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 &
!            + et1(q1,q,e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2 + &
!           et1(q1,q,e4,q4)*e1_e2*e3_q4*yyy7*xxx2 - 4.D0*et1(q2,q,e3,e4)* &
!           q1_q3*e1_e2*yyy6*xxx2 + et1(q2,q,e3,e4)*e1_e2*MZ4**2*yyy6* &
!           xxx2 - et1(q2,q,e3,e4)*e1_e2*MZ3**2*yyy6*xxx2 - et1(q2,q,e3, &
!           e4)*e1_e2*MG**2*yyy6*xxx2 + 4.D0*et1(q2,q,e3,q3)*q1_q3*e1_e2* &
!           e4_q3*MG**(-2)*yyy7*xxx2 - et1(q2,q,e3,q3)*e1_e2*e4_q3* &
!           MG**(-2)*MZ4**2*yyy7*xxx2
!       res = res + et1(q2,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 &
!            + et1(q2,q,e3,q3)*e1_e2*e4_q3*yyy7*xxx2 - 4.D0*et1(q2,q,e3, &
!           q4)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7*xxx2 + et1(q2,q,e3,q4)* &
!           e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q2,q,e3,q4)*e1_e2 &
!           *e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q2,q,e3,q4)*e1_e2* &
!           e4_q3*yyy7*xxx2 + 4.D0*et1(q2,q,e4,q3)*q1_q3*e1_e2*e3_q4* &
!           MG**(-2)*yyy7*xxx2 - et1(q2,q,e4,q3)*e1_e2*e3_q4*MG**(-2)* &
!           MZ4**2*yyy7*xxx2 + et1(q2,q,e4,q3)*e1_e2*e3_q4*MG**(-2)* &
!           MZ3**2*yyy7*xxx2 + et1(q2,q,e4,q3)*e1_e2*e3_q4*yyy7*xxx2 - 4.D0 &
!           *et1(q2,q,e4,q4)*q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7*xxx2 + et1( &
!           q2,q,e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q2,q, &
!           e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q2,q,e4,q4 &
!           )*e1_e2*e3_q4*yyy7*xxx2 - et1(q,e1,e3,e4)*e2_q3*MG**2*yyy6* &
!           xxx1 + et1(q,e1,e3,q3)*e2_q3*e4_q3*yyy7*xxx1 - et1(q,e1,e3,q4 &
!           )*e2_q3*e4_q3*yyy7*xxx1 + et1(q,e1,e4,q3)*e2_q3*e3_q4*yyy7* &
!           xxx1
!       res = res - et1(q,e1,e4,q4)*e2_q3*e3_q4*yyy7*xxx1 - et1(q,e2,e3, &
!           e4)*e1_q3*MG**2*yyy6*xxx1 + et1(q,e2,e3,q3)*e1_q3*e4_q3*yyy7* &
!           xxx1 - et1(q,e2,e3,q4)*e1_q3*e4_q3*yyy7*xxx1 + et1(q,e2,e4,q3 &
!           )*e1_q3*e3_q4*yyy7*xxx1 - et1(q,e2,e4,q4)*e1_q3*e3_q4*yyy7* &
!           xxx1 - 1.D0/3.D0*et1(q,e3,e4,q3)*e1_e2*MG**2*yyy6*xxx2 + 1.D0/ &
!           3.D0*et1(q,e3,e4,q3)*e1_e2*MG**2*yyy6*xxx1 + 1.D0/3.D0*et1(q, &
!           e3,e4,q4)*e1_e2*MG**2*yyy6*xxx2 - 1.D0/3.D0*et1(q,e3,e4,q4)* &
!           e1_e2*MG**2*yyy6*xxx1 - 8.D0*et1(e3,e4,q3,q4)*q1_q3*e1_e2* &
!           MG**(-2)*MZ4**2*yyy5*xxx2 + 8.D0*et1(e3,e4,q3,q4)*q1_q3*e1_e2 &
!           *MG**(-2)*MZ3**2*yyy5*xxx2 + 8.D0*et1(e3,e4,q3,q4)*q1_q3* &
!           e1_e2*yyy5*xxx2 + 16.D0*et1(e3,e4,q3,q4)*q1_q3**2*e1_e2* &
!           MG**(-2)*yyy5*xxx2 + 2.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
!           MG**(-2)*MZ4**4*yyy5*xxx2 + 1.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
!           MG**(-2)*MZ4**4*yyy5*xxx1 - 4.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
!           MG**(-2)*MZ3**2*MZ4**2*yyy5*xxx2 - 2.D0/3.D0*et1(e3,e4,q3,q4) &
!           *e1_e2*MG**(-2)*MZ3**2*MZ4**2*yyy5*xxx1
!       res = res + 2.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*MZ3**4*yyy5 &
!        *xxx2 + 1.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*MZ3**4*yyy5* &
!           xxx1 - 4.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2*MZ4**2*yyy5*xxx2 - 2.D0 &
!           /3.D0*et1(e3,e4,q3,q4)*e1_e2*MZ4**2*yyy5*xxx1 + 8.D0/3.D0* &
!           et1(e3,e4,q3,q4)*e1_e2*MZ3**2*yyy5*xxx2 - 2.D0/3.D0*et1(e3,e4 &
!           ,q3,q4)*e1_e2*MZ3**2*yyy5*xxx1 + 2.D0/3.D0*et1(e3,e4,q3,q4)* &
!           e1_e2*MG**2*yyy5*xxx2 + 1.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
!           MG**2*yyy5*xxx1 + 4.D0*et1(e3,e4,q3,q4)*e1_q3*e2_q3*yyy5*xxx1
!
! print *, "old res GG",res


!   this is the new code that includes couplings yyy41 and yyy42 instead of yyy4
!    res =&
!        &  + 8.*q1_e3*q1_e4*e1_e2*yyy1*xxx2 - 8.*q1_e3*q1_q3*e1_e2*e4_q3*&
!        &    yyy41*xxx2 + 4.*q1_e3*e1_e2*e4_q3*yyy1*xxx2 + 2.*q1_e3*e1_e2*&
!        &    e4_q3*MZ4**2*yyy41*xxx2 - 2.*q1_e3*e1_e2*e4_q3*MZ3**2*yyy41*&
!        &    xxx2 - 2.*q1_e3*e1_e2*e4_q3*MG**2*yyy41*xxx2 + 8.*q1_e4*q1_q3&
!        &    *e1_e2*e3_q4*yyy42*xxx2 + 4.*q1_e4*e1_e2*e3_q4*yyy1*xxx2 - 2.&
!        &    *q1_e4*e1_e2*e3_q4*MZ4**2*yyy42*xxx2 + 2.*q1_e4*e1_e2*e3_q4*&
!        &    MZ3**2*yyy42*xxx2 + 2.*q1_e4*e1_e2*e3_q4*MG**2*yyy42*xxx2 - 8.&
!        &    *q1_q3*e1_e2*e3_e4*MZ4**2*yyy2*xxx2 + 8.*q1_q3*e1_e2*e3_e4*&
!        &    MZ3**2*yyy2*xxx2 + 8.*q1_q3*e1_e2*e3_e4*MG**2*yyy2*xxx2 + 4.*&
!        &    q1_q3*e1_e2*e3_q4*e4_q3*yyy42*xxx2 - 4.*q1_q3*e1_e2*e3_q4*&
!        &    e4_q3*yyy41*xxx2 - 8.*q1_q3*e1_e2*e3_q4*e4_q3*MZ4**2*yyy3*&
!        &    xxx2 + 8.*q1_q3*e1_e2*e3_q4*e4_q3*MZ3**2*yyy3*xxx2 + 8.*q1_q3&
!        &    *e1_e2*e3_q4*e4_q3*MG**2*yyy3*xxx2 + 16.*q1_q3**2*e1_e2*e3_e4&
!        &    *yyy2*xxx2 + 16.*q1_q3**2*e1_e2*e3_q4*e4_q3*yyy3*xxx2 + 2./3.&
!        &    *e1_e2*e3_e4*MZ4**4*yyy2*xxx2
!          res = res + 1./3.*e1_e2*e3_e4*MZ4**4*yyy2*xxx1 - 4./3.*e1_e2*&
!        &    e3_e4*MZ3**2*MZ4**2*yyy2*xxx2 - 2./3.*e1_e2*e3_e4*MZ3**2*&
!        &    MZ4**2*yyy2*xxx1 + 2./3.*e1_e2*e3_e4*MZ3**4*yyy2*xxx2 + 1./3.&
!        &    *e1_e2*e3_e4*MZ3**4*yyy2*xxx1 + 2./3.*e1_e2*e3_e4*MG**2*yyy1*&
!        &    xxx2 - 2./3.*e1_e2*e3_e4*MG**2*yyy1*xxx1 - 4./3.*e1_e2*e3_e4*&
!        &    MG**2*MZ4**2*yyy2*xxx2 - 2./3.*e1_e2*e3_e4*MG**2*MZ4**2*yyy2*&
!        &    xxx1 + 8./3.*e1_e2*e3_e4*MG**2*MZ3**2*yyy2*xxx2 - 2./3.*e1_e2&
!        &    *e3_e4*MG**2*MZ3**2*yyy2*xxx1 + 2./3.*e1_e2*e3_e4*MG**4*yyy2*&
!        &    xxx2 + 1./3.*e1_e2*e3_e4*MG**4*yyy2*xxx1 + 4./3.*e1_e2*e3_q4*&
!        &    e4_q3*yyy1*xxx2 + 2./3.*e1_e2*e3_q4*e4_q3*yyy1*xxx1 - 2./3.*&
!        &    e1_e2*e3_q4*e4_q3*MZ4**2*yyy42*xxx2 - 1./3.*e1_e2*e3_q4*e4_q3&
!        &    *MZ4**2*yyy42*xxx1 + 2./3.*e1_e2*e3_q4*e4_q3*MZ4**2*yyy41*&
!        &    xxx2 + 1./3.*e1_e2*e3_q4*e4_q3*MZ4**2*yyy41*xxx1 + 2./3.*&
!        &    e1_e2*e3_q4*e4_q3*MZ4**4*yyy3*xxx2 + 1./3.*e1_e2*e3_q4*e4_q3*&
!        &    MZ4**4*yyy3*xxx1 + 2./3.*e1_e2*e3_q4*e4_q3*MZ3**2*yyy42*xxx2&
!        &     + 1./3.*e1_e2*e3_q4*e4_q3*MZ3**2*yyy42*xxx1
!          res = res - 2./3.*e1_e2*e3_q4*e4_q3*MZ3**2*yyy41*xxx2 - 1./3.*&
!        &    e1_e2*e3_q4*e4_q3*MZ3**2*yyy41*xxx1 - 4./3.*e1_e2*e3_q4*e4_q3&
!        &    *MZ3**2*MZ4**2*yyy3*xxx2 - 2./3.*e1_e2*e3_q4*e4_q3*MZ3**2*&
!        &    MZ4**2*yyy3*xxx1 + 2./3.*e1_e2*e3_q4*e4_q3*MZ3**4*yyy3*xxx2&
!        &     + 1./3.*e1_e2*e3_q4*e4_q3*MZ3**4*yyy3*xxx1 + 4./3.*e1_e2*&
!        &    e3_q4*e4_q3*MG**2*yyy42*xxx2 - 1./3.*e1_e2*e3_q4*e4_q3*MG**2*&
!        &    yyy42*xxx1 - 2./3.*e1_e2*e3_q4*e4_q3*MG**2*yyy41*xxx2 - 1./3.&
!        &    *e1_e2*e3_q4*e4_q3*MG**2*yyy41*xxx1 - 4./3.*e1_e2*e3_q4*e4_q3&
!        &    *MG**2*MZ4**2*yyy3*xxx2 - 2./3.*e1_e2*e3_q4*e4_q3*MG**2*&
!        &    MZ4**2*yyy3*xxx1 + 8./3.*e1_e2*e3_q4*e4_q3*MG**2*MZ3**2*yyy3*&
!        &    xxx2 - 2./3.*e1_e2*e3_q4*e4_q3*MG**2*MZ3**2*yyy3*xxx1 + 2./3.&
!        &    *e1_e2*e3_q4*e4_q3*MG**4*yyy3*xxx2 + 1./3.*e1_e2*e3_q4*e4_q3*&
!        &    MG**4*yyy3*xxx1 + e1_e3*e2_e4*MG**2*yyy1*xxx1 - e1_e3*e2_q3*&
!        &    e4_q3*MG**2*yyy41*xxx1 + e1_e4*e2_e3*MG**2*yyy1*xxx1 + e1_e4*&
!        &    e2_q3*e3_q4*MG**2*yyy42*xxx1 - e1_q3*e2_e3*e4_q3*MG**2*yyy41*&
!        &    xxx1
!          res = res + e1_q3*e2_e4*e3_q4*MG**2*yyy42*xxx1 + 4.*e1_q3*e2_q3*&
!        &    e3_e4*MG**2*yyy2*xxx1 + 4.*e1_q3*e2_q3*e3_q4*e4_q3*MG**2*yyy3&
!        &    *xxx1 + 8.*et1(q1,q2,e1,e2)*q1_e3*q1_e4*MG**(-2)*yyy1*xxx3 - &
!        &    8.*et1(q1,q2,e1,e2)*q1_e3*q1_q3*e4_q3*MG**(-2)*yyy41*xxx3 + 4.&
!        &    *et1(q1,q2,e1,e2)*q1_e3*e4_q3*MG**(-2)*yyy1*xxx3 + 2.*et1(q1,&
!        &    q2,e1,e2)*q1_e3*e4_q3*MG**(-2)*MZ4**2*yyy41*xxx3 - 2.*et1(q1,&
!        &    q2,e1,e2)*q1_e3*e4_q3*MG**(-2)*MZ3**2*yyy41*xxx3 - 2.*et1(q1,&
!        &    q2,e1,e2)*q1_e3*e4_q3*yyy41*xxx3 + 8.*et1(q1,q2,e1,e2)*q1_e4*&
!        &    q1_q3*e3_q4*MG**(-2)*yyy42*xxx3 + 4.*et1(q1,q2,e1,e2)*q1_e4*&
!        &    e3_q4*MG**(-2)*yyy1*xxx3 - 2.*et1(q1,q2,e1,e2)*q1_e4*e3_q4*&
!        &    MG**(-2)*MZ4**2*yyy42*xxx3 + 2.*et1(q1,q2,e1,e2)*q1_e4*e3_q4*&
!        &    MG**(-2)*MZ3**2*yyy42*xxx3 + 2.*et1(q1,q2,e1,e2)*q1_e4*e3_q4*&
!        &    yyy42*xxx3 - 8.*et1(q1,q2,e1,e2)*q1_q3*e3_e4*MG**(-2)*MZ4**2*&
!        &    yyy2*xxx3 + 8.*et1(q1,q2,e1,e2)*q1_q3*e3_e4*MG**(-2)*MZ3**2*&
!        &    yyy2*xxx3 + 8.*et1(q1,q2,e1,e2)*q1_q3*e3_e4*yyy2*xxx3 + 4.*&
!        &    et1(q1,q2,e1,e2)*q1_q3*e3_q4*e4_q3*MG**(-2)*yyy42*xxx3
!          res = res - 4.*et1(q1,q2,e1,e2)*q1_q3*e3_q4*e4_q3*MG**(-2)*yyy41*&
!        & xxx3 - 8.*et1(q1,q2,e1,e2)*q1_q3*e3_q4*e4_q3*MG**(-2)*MZ4**2*&
!        &    yyy3*xxx3 + 8.*et1(q1,q2,e1,e2)*q1_q3*e3_q4*e4_q3*MG**(-2)*&
!        &    MZ3**2*yyy3*xxx3 + 8.*et1(q1,q2,e1,e2)*q1_q3*e3_q4*e4_q3*yyy3&
!        &    *xxx3 + 16.*et1(q1,q2,e1,e2)*q1_q3**2*e3_e4*MG**(-2)*yyy2*&
!        &    xxx3 + 16.*et1(q1,q2,e1,e2)*q1_q3**2*e3_q4*e4_q3*MG**(-2)*&
!        &    yyy3*xxx3 + 2./3.*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)*MZ4**4*yyy2&
!        &    *xxx3 - 4./3.*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)*MZ3**2*MZ4**2*&
!        &    yyy2*xxx3 + 2./3.*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)*MZ3**4*yyy2&
!        &    *xxx3 + 2./3.*et1(q1,q2,e1,e2)*e3_e4*yyy1*xxx3 - 4./3.*et1(q1&
!        &    ,q2,e1,e2)*e3_e4*MZ4**2*yyy2*xxx3 + 8./3.*et1(q1,q2,e1,e2)*&
!        &    e3_e4*MZ3**2*yyy2*xxx3 + 2./3.*et1(q1,q2,e1,e2)*e3_e4*MG**2*&
!        &    yyy2*xxx3 + 4./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*yyy1*&
!        &    xxx3 - 2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ4**2*&
!        &    yyy42*xxx3 + 2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*&
!        &    MZ4**2*yyy41*xxx3
!          res = res + 2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ4**4*&
!        & yyy3*xxx3 + 2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ3**2*&
!        &    yyy42*xxx3 - 2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*&
!        &    MZ3**2*yyy41*xxx3 - 4./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*&
!        &    MG**(-2)*MZ3**2*MZ4**2*yyy3*xxx3 + 2./3.*et1(q1,q2,e1,e2)*&
!        &    e3_q4*e4_q3*MG**(-2)*MZ3**4*yyy3*xxx3 + 4./3.*et1(q1,q2,e1,e2&
!        &    )*e3_q4*e4_q3*yyy42*xxx3 - 2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3&
!        &    *yyy41*xxx3 - 4./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MZ4**2*yyy3*&
!        &    xxx3 + 8./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MZ3**2*yyy3*xxx3 + &
!        &    2./3.*et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**2*yyy3*xxx3 - et1(q1,&
!        &    q2,e1,e2)*et1(q1,q,e3,e4)*MG**(-2)*MZ4**2*yyy6*xxx3 + et1(q1,&
!        &    q2,e1,e2)*et1(q1,q,e3,e4)*MG**(-2)*MZ3**2*yyy6*xxx3 + et1(q1,&
!        &    q2,e1,e2)*et1(q1,q,e3,e4)*yyy6*xxx3 + 4.*et1(q1,q2,e1,e2)*&
!        &    et1(q1,q,e3,e4)*q1_q3*MG**(-2)*yyy6*xxx3 - 4.*et1(q1,q2,e1,e2&
!        &    )*et1(q1,q,e3,q3)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 + et1(q1,q2,&
!        &    e1,e2)*et1(q1,q,e3,q3)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3
!          res = res - et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*e4_q3*MG**(-4)*&
!        & MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*e4_q3*&
!        &    MG**(-2)*yyy7*xxx3 + 4.*et1(q1,q2,e1,e2)*et1(q1,q,e3,q4)*&
!        &    q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e3&
!        &    ,q4)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(&
!        &    q1,q,e3,q4)*e4_q3*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2&
!        &    )*et1(q1,q,e3,q4)*e4_q3*MG**(-2)*yyy7*xxx3 - 4.*et1(q1,q2,e1,&
!        &    e2)*et1(q1,q,e4,q3)*q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 + et1(q1,&
!        &    q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-4)*MZ4**2*yyy7*xxx3 - &
!        &    et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-4)*MZ3**2*yyy7*&
!        &    xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-2)*yyy7*&
!        &    xxx3 + 4.*et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*q1_q3*e3_q4*&
!        &    MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*e3_q4*&
!        &    MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*&
!        &    e3_q4*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q1,q,&
!        &    e4,q4)*e3_q4*MG**(-2)*yyy7*xxx3
!          res = res + et1(q1,q2,e1,e2)*et1(q2,q,e3,e4)*MG**(-2)*MZ4**2*yyy6&
!        & *xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e3,e4)*MG**(-2)*MZ3**2*yyy6*&
!        &    xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e3,e4)*yyy6*xxx3 - 4.*et1(q1&
!        &    ,q2,e1,e2)*et1(q2,q,e3,e4)*q1_q3*MG**(-2)*yyy6*xxx3 + 4.*et1(&
!        &    q1,q2,e1,e2)*et1(q2,q,e3,q3)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3&
!        &     - et1(q1,q2,e1,e2)*et1(q2,q,e3,q3)*e4_q3*MG**(-4)*MZ4**2*&
!        &    yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q2,q,e3,q3)*e4_q3*MG**(-4)*&
!        &    MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q2,q,e3,q3)*e4_q3*&
!        &    MG**(-2)*yyy7*xxx3 - 4.*et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*&
!        &    q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q2,q,e3&
!        &    ,q4)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(&
!        &    q2,q,e3,q4)*e4_q3*MG**(-4)*MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2&
!        &    )*et1(q2,q,e3,q4)*e4_q3*MG**(-2)*yyy7*xxx3 + 4.*et1(q1,q2,e1,&
!        &    e2)*et1(q2,q,e4,q3)*q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 - et1(q1,&
!        &    q2,e1,e2)*et1(q2,q,e4,q3)*e3_q4*MG**(-4)*MZ4**2*yyy7*xxx3 + &
!        &    et1(q1,q2,e1,e2)*et1(q2,q,e4,q3)*e3_q4*MG**(-4)*MZ3**2*yyy7*&
!        &    xxx3
!          res = res + et1(q1,q2,e1,e2)*et1(q2,q,e4,q3)*e3_q4*MG**(-2)*yyy7*&
!        & xxx3 - 4.*et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*q1_q3*e3_q4*MG**(-4)*&
!        &    yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-4)*&
!        &    MZ4**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*&
!        &    MG**(-4)*MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*&
!        &    e3_q4*MG**(-2)*yyy7*xxx3 - 1./3.*et1(q1,q2,e1,e2)*et1(q,e3,e4&
!        &    ,q3)*yyy6*xxx3 + 1./3.*et1(q1,q2,e1,e2)*et1(q,e3,e4,q4)*yyy6*&
!        &    xxx3 + 2./3.*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*MG**(-4)*&
!        &    MZ4**4*yyy5*xxx3 - 4./3.*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*&
!        &    MG**(-4)*MZ3**2*MZ4**2*yyy5*xxx3 + 2./3.*et1(q1,q2,e1,e2)*&
!        &    et1(e3,e4,q3,q4)*MG**(-4)*MZ3**4*yyy5*xxx3 - 4./3.*et1(q1,q2,&
!        &    e1,e2)*et1(e3,e4,q3,q4)*MG**(-2)*MZ4**2*yyy5*xxx3 + 8./3.*&
!        &    et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*MG**(-2)*MZ3**2*yyy5*xxx3&
!        &     + 2./3.*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*yyy5*xxx3 - 8.*&
!        &    et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*MG**(-4)*MZ4**2*yyy5*&
!        &    xxx3
!          res = res + 8.*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*MG**(-4)*&
!        & MZ3**2*yyy5*xxx3 + 8.*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*&
!        &    MG**(-2)*yyy5*xxx3 + 16.*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*&
!        &    q1_q3**2*MG**(-4)*yyy5*xxx3 + 4.*et1(q1,q,e3,e4)*q1_q3*e1_e2*&
!        &    yyy6*xxx2 - et1(q1,q,e3,e4)*e1_e2*MZ4**2*yyy6*xxx2 + et1(q1,q&
!        &    ,e3,e4)*e1_e2*MZ3**2*yyy6*xxx2 + et1(q1,q,e3,e4)*e1_e2*MG**2*&
!        &    yyy6*xxx2 - 4.*et1(q1,q,e3,q3)*q1_q3*e1_e2*e4_q3*MG**(-2)*&
!        &    yyy7*xxx2 + et1(q1,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*&
!        &    xxx2 - et1(q1,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2&
!        &     - et1(q1,q,e3,q3)*e1_e2*e4_q3*yyy7*xxx2 + 4.*et1(q1,q,e3,q4)&
!        &    *q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7*xxx2 - et1(q1,q,e3,q4)*e1_e2&
!        &    *e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 + et1(q1,q,e3,q4)*e1_e2*&
!        &    e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 + et1(q1,q,e3,q4)*e1_e2*e4_q3&
!        &    *yyy7*xxx2 - 4.*et1(q1,q,e4,q3)*q1_q3*e1_e2*e3_q4*MG**(-2)*&
!        &    yyy7*xxx2 + et1(q1,q,e4,q3)*e1_e2*e3_q4*MG**(-2)*MZ4**2*yyy7*&
!        &    xxx2
!          res = res - et1(q1,q,e4,q3)*e1_e2*e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2&
!        &     - et1(q1,q,e4,q3)*e1_e2*e3_q4*yyy7*xxx2 + 4.*et1(q1,q,e4,q4)&
!        &    *q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7*xxx2 - et1(q1,q,e4,q4)*e1_e2&
!        &    *e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 + et1(q1,q,e4,q4)*e1_e2*&
!        &    e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2 + et1(q1,q,e4,q4)*e1_e2*e3_q4&
!        &    *yyy7*xxx2 - 4.*et1(q2,q,e3,e4)*q1_q3*e1_e2*yyy6*xxx2 + et1(&
!        &    q2,q,e3,e4)*e1_e2*MZ4**2*yyy6*xxx2 - et1(q2,q,e3,e4)*e1_e2*&
!        &    MZ3**2*yyy6*xxx2 - et1(q2,q,e3,e4)*e1_e2*MG**2*yyy6*xxx2 + 4.&
!        &    *et1(q2,q,e3,q3)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7*xxx2 - et1(&
!        &    q2,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 + et1(q2,q,&
!        &    e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 + et1(q2,q,e3,q3&
!        &    )*e1_e2*e4_q3*yyy7*xxx2 - 4.*et1(q2,q,e3,q4)*q1_q3*e1_e2*&
!        &    e4_q3*MG**(-2)*yyy7*xxx2 + et1(q2,q,e3,q4)*e1_e2*e4_q3*&
!        &    MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q2,q,e3,q4)*e1_e2*e4_q3*&
!        &    MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q2,q,e3,q4)*e1_e2*e4_q3*yyy7*&
!        &    xxx2
!          res = res + 4.*et1(q2,q,e4,q3)*q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7*&
!        & xxx2 - et1(q2,q,e4,q3)*e1_e2*e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 + &
!        &    et1(q2,q,e4,q3)*e1_e2*e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2 + et1(&
!        &    q2,q,e4,q3)*e1_e2*e3_q4*yyy7*xxx2 - 4.*et1(q2,q,e4,q4)*q1_q3*&
!        &    e1_e2*e3_q4*MG**(-2)*yyy7*xxx2 + et1(q2,q,e4,q4)*e1_e2*e3_q4*&
!        &    MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q2,q,e4,q4)*e1_e2*e3_q4*&
!        &    MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q2,q,e4,q4)*e1_e2*e3_q4*yyy7*&
!        &    xxx2 - et1(q,e1,e3,e4)*e2_q3*MG**2*yyy6*xxx1 + et1(q,e1,e3,q3&
!        &    )*e2_q3*e4_q3*yyy7*xxx1 - et1(q,e1,e3,q4)*e2_q3*e4_q3*yyy7*&
!        &    xxx1 + et1(q,e1,e4,q3)*e2_q3*e3_q4*yyy7*xxx1 - et1(q,e1,e4,q4&
!        &    )*e2_q3*e3_q4*yyy7*xxx1 - et1(q,e2,e3,e4)*e1_q3*MG**2*yyy6*&
!        &    xxx1 + et1(q,e2,e3,q3)*e1_q3*e4_q3*yyy7*xxx1 - et1(q,e2,e3,q4&
!        &    )*e1_q3*e4_q3*yyy7*xxx1 + et1(q,e2,e4,q3)*e1_q3*e3_q4*yyy7*&
!        &    xxx1 - et1(q,e2,e4,q4)*e1_q3*e3_q4*yyy7*xxx1 - 1./3.*et1(q,e3&
!        &    ,e4,q3)*e1_e2*MG**2*yyy6*xxx2 + 1./3.*et1(q,e3,e4,q3)*e1_e2*&
!        &    MG**2*yyy6*xxx1
!          res = res + 1./3.*et1(q,e3,e4,q4)*e1_e2*MG**2*yyy6*xxx2 - 1./3.*&
!        &    et1(q,e3,e4,q4)*e1_e2*MG**2*yyy6*xxx1 - 8.*et1(e3,e4,q3,q4)*&
!        &    q1_q3*e1_e2*MG**(-2)*MZ4**2*yyy5*xxx2 + 8.*et1(e3,e4,q3,q4)*&
!        &    q1_q3*e1_e2*MG**(-2)*MZ3**2*yyy5*xxx2 + 8.*et1(e3,e4,q3,q4)*&
!        &    q1_q3*e1_e2*yyy5*xxx2 + 16.*et1(e3,e4,q3,q4)*q1_q3**2*e1_e2*&
!        &    MG**(-2)*yyy5*xxx2 + 2./3.*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*&
!        &    MZ4**4*yyy5*xxx2 + 1./3.*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*&
!        &    MZ4**4*yyy5*xxx1 - 4./3.*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*&
!        &    MZ3**2*MZ4**2*yyy5*xxx2 - 2./3.*et1(e3,e4,q3,q4)*e1_e2*&
!        &    MG**(-2)*MZ3**2*MZ4**2*yyy5*xxx1 + 2./3.*et1(e3,e4,q3,q4)*&
!        &    e1_e2*MG**(-2)*MZ3**4*yyy5*xxx2 + 1./3.*et1(e3,e4,q3,q4)*&
!        &    e1_e2*MG**(-2)*MZ3**4*yyy5*xxx1 - 4./3.*et1(e3,e4,q3,q4)*&
!        &    e1_e2*MZ4**2*yyy5*xxx2 - 2./3.*et1(e3,e4,q3,q4)*e1_e2*MZ4**2*&
!        &    yyy5*xxx1 + 8./3.*et1(e3,e4,q3,q4)*e1_e2*MZ3**2*yyy5*xxx2 - 2.&
!          &   /3.*et1(e3,e4,q3,q4)*e1_e2*MZ3**2*yyy5*xxx1 + 2./3.*et1(e3,e4&
!        &    ,q3,q4)*e1_e2*MG**2*yyy5*xxx2
!          res = res + 1./3.*et1(e3,e4,q3,q4)*e1_e2*MG**2*yyy5*xxx1 + 4.*&
!        &    et1(e3,e4,q3,q4)*e1_q3*e2_q3*yyy5*xxx1


! print *, "new  ",res


! ! ! ! ! ! ! ! ! ! ! ! ! ! !  SAME AS ABOVE BUT SHORTER  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

         abr1 = MG**2 + MZ3**2 - MZ4**2 + 4*q1_q3

         res =  MG**2*xxx1*(e1_e4*(yyy1*e2_e3 + yyy42*e2_q3*e3_q4) + e1_e3*(yyy1*e2_e4 - yyy41*e2_q3*e4_q3) + &
            e1_q3*(yyy42*e2_e4*e3_q4 - yyy41*e2_e3*e4_q3 + 4*e2_q3*(yyy2*e3_e4 + yyy3*e3_q4*e4_q3))) + &
         e1_e2*((xxx1*((MG**4*yyy2 + (MZ3**2 - MZ4**2)**2*yyy2 - 2*MG**2*(yyy1 + (MZ3**2 + MZ4**2)*yyy2))*e3_e4 + &
                 (2*yyy1 + MG**4*yyy3 + (MZ3 - MZ4)*(MZ3 + MZ4)*((MZ3 - MZ4)*(MZ3 + MZ4)*yyy3 - yyy41 + yyy42) - MG**2*(2*(MZ3**2 + MZ4**2)*yyy3 + yyy41 + yyy42))* &
                  e3_q4*e4_q3))/3d0 + (2*xxx2*(3*q1_e3*((2*yyy1 - abr1*yyy41)*e4_q3 + 4*yyy1*q1_e4) + &
                 e3_e4*(abr1**2*yyy2 + abr1*(yyy1 + 2*MZ3**2*yyy2) - (MZ3 - MZ4)*(MZ3 + MZ4)*(yyy1 + 2*MZ3**2*yyy2) - &
                    4*q1_q3*(yyy1 - abr1*yyy2 + 2*MZ3**2*yyy2 + 8*yyy2*q1_q3) + 24*yyy2*q1_q3**2) + &
                 e3_q4*(3*(2*yyy1 + abr1*yyy42)*q1_e4 + e4_q3* &
                     (2*yyy1 + abr1**2*yyy3 - (MZ3 - MZ4)*(MZ3 + MZ4)*(2*MZ3**2*yyy3 + yyy42) + abr1*(2*MZ3**2*yyy3 - yyy41 + 2*yyy42 ) - &
                       2*q1_q3*(-2*abr1*yyy3 + 4*MZ3**2*yyy3 + yyy41 + yyy42 + 16*yyy3*q1_q3) + 24*yyy3*q1_q3**2))))/3d0) - &
         MG**2*xxx1*yyy6*e2_q3*et1(q,e1,e3,e4) + xxx1*yyy7*e2_q3*e4_q3*et1(q,e1,e3,q3) - xxx1*yyy7*e2_q3*e4_q3*et1(q,e1,e3,q4) + &
         xxx1*yyy7*e2_q3*e3_q4*et1(q,e1,e4,q3) - xxx1*yyy7*e2_q3*e3_q4*et1(q,e1,e4,q4) - MG**2*xxx1*yyy6*e1_q3*et1(q,e2,e3,e4) + &
         xxx1*yyy7*e1_q3*e4_q3*et1(q,e2,e3,q3) - xxx1*yyy7*e1_q3*e4_q3*et1(q,e2,e3,q4) + xxx1*yyy7*e1_q3*e3_q4*et1(q,e2,e4,q3) - &
         xxx1*yyy7*e1_q3*e3_q4*et1(q,e2,e4,q4) + et1(q,e3,e4,q3)* &
          (((MG**2*xxx1*yyy6)/3d0 - (MG**2*xxx2*yyy6)/3d0)*e1_e2 - (xxx3*yyy6*et1(q1,q2,e1,e2))/3d0) + &
         et1(q,e3,e4,q4)*((-(MG**2*xxx1*yyy6)/3d0 + (MG**2*xxx2*yyy6)/3d0)*e1_e2 + (xxx3*yyy6*et1(q1,q2,e1,e2))/3d0) + &
         et1(q1,q,e3,e4)*(abr1*xxx2*yyy6*e1_e2 + (abr1*xxx3*yyy6*et1(q1,q2,e1,e2))/MG**2) + &
         et1(q1,q,e4,q3)*(-((abr1*xxx2*yyy7*e1_e2*e3_q4)/MG**2) - (abr1*xxx3*yyy7*e3_q4*et1(q1,q2,e1,e2))/MG**4) + &
         et1(q1,q,e4,q4)*((abr1*xxx2*yyy7*e1_e2*e3_q4)/MG**2 + (abr1*xxx3*yyy7*e3_q4*et1(q1,q2,e1,e2))/MG**4) + &
         et1(q1,q,e3,q3)*(-((abr1*xxx2*yyy7*e1_e2*e4_q3)/MG**2) - (abr1*xxx3*yyy7*e4_q3*et1(q1,q2,e1,e2))/MG**4) + &
         et1(q1,q,e3,q4)*((abr1*xxx2*yyy7*e1_e2*e4_q3)/MG**2 + (abr1*xxx3*yyy7*e4_q3*et1(q1,q2,e1,e2))/MG**4) + &
         et1(e3,e4,q3,q4)*(4*xxx1*yyy5*e1_q3*e2_q3 + e1_e2* &
             (((MG - MZ3 - MZ4)*(MG + MZ3 - MZ4)*(MG - MZ3 + MZ4)*(MG + MZ3 + MZ4)*xxx1*yyy5)/(3d0*MG**2) + &
               (2*xxx2*yyy5*(abr1**2 + 2*abr1*MZ3**2 - 2*MZ3**4 + 2*MZ3**2*MZ4**2 + 4*(abr1 - 2*MZ3**2 - 8*q1_q3)*q1_q3 + 24*q1_q3**2))/(3d0*MG**2)) + &
            (2*xxx3*yyy5*(abr1**2 + 2*abr1*MZ3**2 - 2*MZ3**4 + 2*MZ3**2*MZ4**2 + 4*(abr1 - 2*MZ3**2 - 8*q1_q3)*q1_q3 + 24*q1_q3**2)*et1(q1,q2,e1,e2))/ &
             (3d0*MG**4)) - abr1*xxx2*yyy6*e1_e2*et1(q2,q,e3,e4) + (abr1*xxx2*yyy7*e1_e2*e4_q3*et1(q2,q,e3,q3))/MG**2 - &
         (abr1*xxx2*yyy7*e1_e2*e4_q3*et1(q2,q,e3,q4))/MG**2 + (abr1*xxx2*yyy7*e1_e2*e3_q4*et1(q2,q,e4,q3))/MG**2 - &
         (abr1*xxx2*yyy7*e1_e2*e3_q4*et1(q2,q,e4,q4))/MG**2 + &
         et1(q1,q2,e1,e2)*((2*xxx3*(3*q1_e3*((2*yyy1 - abr1*yyy41)*e4_q3 + 4*yyy1*q1_e4) + &
                 e3_e4*(abr1**2*yyy2 + abr1*(yyy1 + 2*MZ3**2*yyy2) - (MZ3 - MZ4)*(MZ3 + MZ4)*(yyy1 + 2*MZ3**2*yyy2) - &
                    4*q1_q3*(yyy1 - abr1*yyy2 + 2*MZ3**2*yyy2 + 8*yyy2*q1_q3) + 24*yyy2*q1_q3**2) + &
                 e3_q4*(3*(2*yyy1 + abr1*yyy42)*q1_e4 + e4_q3* &
                     (2*yyy1 + abr1**2*yyy3 - (MZ3 - MZ4)*(MZ3 + MZ4)*(2*MZ3**2*yyy3 + yyy42) + abr1*(2*MZ3**2*yyy3 - yyy41 + 2*yyy42) - &
                       2*q1_q3*(-2*abr1*yyy3 + 4*MZ3**2*yyy3 + yyy41 + yyy42 + 16*yyy3*q1_q3) + 24*yyy3*q1_q3**2))))/(3d0*MG**2) - &
            (abr1*xxx3*yyy6*et1(q2,q,e3,e4))/MG**2 + (abr1*xxx3*yyy7*e4_q3*et1(q2,q,e3,q3))/MG**4 - (abr1*xxx3*yyy7*e4_q3*et1(q2,q,e3,q4))/MG**4 +  &
            (abr1*xxx3*yyy7*e3_q4*et1(q2,q,e4,q3))/MG**4 - (abr1*xxx3*yyy7*e3_q4*et1(q2,q,e4,q4))/MG**4)


! print *, "newer",res
! pause

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

     else



      res = q1_e3*q2_e4*e1_e2 * ( M_Reso**2 )
      res = res + q1_e3*q2_q4*e1_e2*e4_q3 * (  - 2.0 )
      res = res + q1_e4*q2_e3*e1_e2 * ( M_Reso**2 )
      res = res + q1_e4*q2_q3*e1_e2*e3_q4 * (  - 2.0 )
      res = res + q1_q3*q2_e4*e1_e2*e3_q4 * (  - 2.0 )
      res = res + q1_q3*q2_q4*e1_e2*e3_e4 * ( 2.0 )
      res = res + q1_q3*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
      res = res + q1_q3*e1_e2*e3_q4*e4_q3 * (  - 1. )
      res = res + q1_q4*q2_e3*e1_e2*e4_q3 * (  - 2. )
      res = res + q1_q4*q2_q3*e1_e2*e3_e4 * ( 2 )
      res = res + q1_q4*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
      res = res + q1_q4*e1_e2*e3_q4*e4_q3 * (  - 1. )
      res = res + q2_q3*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
      res = res + q2_q3*e1_e2*e3_q4*e4_q3 * (  - 1 )
      res = res + q2_q4*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
      res = res + q2_q4*e1_e2*e3_q4*e4_q3 * (  - 1. )
      res = res + e1_e2*e3_e4 * ( M_Reso**2*M_V**2 - 1./2.*M_Reso**4 )
      res = res + e1_e2*e3_q4*e4_q3 * ( M_Reso**2 )
      res = res + e1_e3*e2_e4 * ( 1./2.*M_Reso**4 )
      res = res + e1_e3*e2_q4*e4_q3 * (  - M_Reso**2 )
      res = res + e1_e4*e2_e3 * ( 1./2.*M_Reso**4 )
      res = res + e1_e4*e2_q3*e3_q4 * (  - M_Reso**2 )
      res = res + e1_q3*e2_e4*e3_q4 * (  - M_Reso**2 )
      res = res + e1_q3*e2_q4*e3_e4 * ( M_Reso**2 )
      res = res + e1_q4*e2_e3*e4_q3 * (  - M_Reso**2 )
      res = res + e1_q4*e2_q3*e3_e4 * ( M_Reso**2 )

      print *,"this code should no longer be used"; stop 1

!print *, "res GG old",res
!pause


     endif

      end subroutine ggGZZampl


!----- a subroutine for q qbar -> G -> Z -> lept + Z --> 2 lepts
!----- all outgoing convention and the following momentum assignment
!-----  0 -> bq(p1) + q(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
     subroutine EvalAmp_qqb_G_VV(p,MY_IDUP,sum)
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      real(dp) ::  pin(4,4)
      complex(dp) :: A(2)
      integer :: i1,i2,i3,i4,ordering(1:4)
      real(dp) :: aL1,aR1,aL2,aR2,qL,qR
      real(dp) :: gZ_sq
      real(dp) :: prefactor
      real(dp) :: intcolfac

      if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
         intcolfac=1.0_dp/3.0_dp
      else
         intcolfac=1.0_dp
      endif

!---- electroweak couplings
!       aL = -one + two*sitW**2
!       aR = aL+one
      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2


!---- chiral couplings of quarks to gravitons
      qL = graviton_qq_left
      qR = graviton_qq_right

!---- full prefactor; 3 is  the color factor
      prefactor = 3d0*gZ_sq**2


         if( CoupledVertex(MY_IDUP(6:7),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL1=0d0
                    aR1=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP(6:7),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP(6)) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ud)
                 aR1 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ln)
                 aR1 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_tn)
                 aR1 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL1=0d0
                 aR1=0d0
              endif
         elseif( MY_IDUP(6).eq.Pho_ ) then !  photon
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL1=0d0
              aR1=0d0
         endif

         if( CoupledVertex(MY_IDUP(8:9),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP(8:9),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP(8)) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ud)
                 aR2 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP(9)).eq.abs(ElM_) .or. abs(MY_IDUP(9)).eq.abs(MuM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ln)
                 aR2 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP(9)).eq.abs(TaM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_tn)
                 aR2 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL2=0d0
                 aR2=0d0
              endif
         elseif( MY_IDUP(8).eq.Pho_ ) then !  photon
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

         if(CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_ .and. CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_) then
            ordering = (/5,6,3,4/)
         else
            ordering = (/3,4,5,6/)
         endif
         call calcHelAmp_qq(ordering,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(1))
         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call swap(ordering(2),ordering(4))
             call calcHelAmp_qq(ordering,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(2))
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

      sum = sum*prefactor

      end subroutine

     subroutine calcHelAmp_qq(ordering,p,idin,i1,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4,idin(3:6)
     real(dp) :: p(1:4,1:6)
     complex(dp) :: propG, propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(4,4)


      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      s  = 2d0 * scr(p(:,1),p(:,2))
      propG = s/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)

         pin(1,:) = p(:,1)
         pin(2,:) = p(:,2)
         sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1)  !qbq
         sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform bookeeping

!-------- -1 == left, 1 == right
         if( idin(l1).ne.Pho_ ) then ! Z, W
            pin(3,:) = p(:,l1)+p(:,l2)
            sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
            sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
            s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
            propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         else !  A
            pin(3,:) = p(:,l1)
            sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
            propZ1=1d0
         endif
         if( idin(l3).ne.Pho_ ) then ! Z, W
            pin(4,:) = p(:,l3)+p(:,l4)
            sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
            sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
            s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
            propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         else ! A
            pin(4,:) = p(:,l3)
            sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')
            propZ2=1d0
         endif

         call qqGZZampl(pin,sp,A(1))

         A(1) = A(1) * propG*propZ1*propZ2

      end subroutine

      subroutine qqGZZampl(p,sp,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(4,4)
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q_q,q3_q3,q4_q4
      complex(dp) :: q1_q2,q1_q3,q1_q4
      complex(dp) :: q2_q3,q2_q4
      complex(dp) :: q3_q4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4),abr1
      complex(dp) :: yyy1,yyy2,yyy3,yyy41,yyy42,yyy5,yyy6,yyy7,yyy4
      real(dp) :: q34,MG,MZ3,MZ4
      real(dp) :: rr

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

      MZ3=dsqrt(cdabs(q3_q3))
      MZ4=dsqrt(cdabs(q4_q4))
      if( use_dynamic_MG ) then
          MG = dsqrt(cdabs(q_q))
      else
          MG = M_Reso
      endif

      q34 = (MG**2-MZ3**2-MZ4**2)/2d0

      if (generate_bis) then
          rr = q34/Lambda**2
!           yyy1 = q34*(b1 + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)  + b5*M_V**2/q34)
          yyy1 = q34*(b1 + b2*rr*(one+MZ3**2/q34)*(one+MZ4**2/q34)  + b5*M_V**2/q34)
!           yyy2 = -b1/two + b3*rr*(1d0-M_V**2/q34) + two*b4*rr+b7*rr*M_V**2/q34
          yyy2 = -b1/two + b3*rr*(1d0-(MZ3**2+MZ4**2)/(2d0*q34)) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
!           yyy4 = -b1 - b2*rr -(b2+b3+2d0*b6)*rr*M_V**2/q34
          yyy41 = -b1 - b2*(q34+MZ3**2)/Lambda**2 - b3*MZ4**2/Lambda**2 - 2d0*b6*M_V**2/Lambda**2
          yyy42 = -b1 - b2*(q34+MZ4**2)/Lambda**2 - b3*MZ3**2/Lambda**2 - 2d0*b6*M_V**2/Lambda**2
          yyy5 = two*b8*rr*MG**2/q34
!           yyy6 = b9
          yyy6 = b9 * M_V**2/Lambda**2
!           yyy7 = b10*rr*MG**2/q34
          yyy7 = b10 * MG**2 * M_V**2/Lambda**4
      else
          yyy1 = (q34)*c1/2d0
          yyy2 = c2
          yyy3 = c3/MG**2
!           yyy4 = c4
          yyy41= c41
          yyy42= c42
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
          if( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
              yyy6=0d0
              yyy7=0d0
          endif
      endif


      res = czero

!      old code without c41 c42 coupligs
!       res =                                                               &
!      q1_e3*e1_e4*yyy1 - q1_e3*e1_q3*e4_q3*yyy4 + q1_e4*e1_e3*yyy1 +       &
!      q1_e4*e1_q3*e3_q4*yyy4 - q1_q3*e1_e3*e4_q3*yyy4 + q1_q3*e1_e4*       &
!      e3_q4*yyy4 + 4*q1_q3*e1_q3*e3_e4*yyy2 + 4*q1_q3*e1_q3*e3_q4*         &
!      e4_q3*yyy3 + 1./2.*e1_e3*e4_q3*yyy1 - 1./4.*e1_e3*e4_q3*MG**2*   &
!      yyy4 + 1./2.*e1_e4*e3_q4*yyy1 + 1./4.*e1_e4*e3_q4*MG**2*yyy4 +   &
!      e1_q3*e3_e4*MG**2*yyy2 + e1_q3*e3_q4*e4_q3*MG**2*yyy3
!
!       res = res +                                                        &
!       1./2.*et1(e3,e4,q1,q)*e1_q3*yyy6 - 1./2.*et1(e3,e4,q2,q)*e1_q3*  &
!       yyy6 + 1./4.*et1(e3,e4,e1,q)*MG**2*yyy6 + et1(e3,e4,e1,q)*q1_q3* &
!       yyy6 + 4*et1(e3,e4,q3,q4)*q1_q3*e1_q3*MG**(-2)*yyy5 + et1(e3,e4, &
!       q3,q4)*e1_q3*yyy5
!
!       res = res +                                                               &
!      1./2.*et1(q1,e3,q,q3)*e1_q3*e4_q3*MG**(-2)*yyy7 - 1./2.*et1(q1,        &
!       e3,q,q4)*e1_q3*e4_q3*MG**(-2)*yyy7 + 1./2.*et1(q1,e4,q,q3)*e1_q3      &
!       *e3_q4*MG**(-2)*yyy7 - 1./2.*et1(q1,e4,q,q4)*e1_q3*e3_q4*             &
!       MG**(-2)*yyy7 - 1./2.*et1(q2,e3,q,q3)*e1_q3*e4_q3*MG**(-2)*yyy7   &
!        + 1./2.*et1(q2,e3,q,q4)*e1_q3*e4_q3*MG**(-2)*yyy7 - 1./2.*et1(       &
!       q2,e4,q,q3)*e1_q3*e3_q4*MG**(-2)*yyy7 + 1./2.*et1(q2,e4,q,q4)*        &
!       e1_q3*e3_q4*MG**(-2)*yyy7 + et1(e1,e3,q,q3)*q1_q3*e4_q3*MG**(-2)  &
!       *yyy7 + 1./4.*et1(e1,e3,q,q3)*e4_q3*yyy7 - et1(e1,e3,q,q4)*q1_q3          &
!       *e4_q3*MG**(-2)*yyy7 - 1./4.*et1(e1,e3,q,q4)*e4_q3*yyy7 + et1(e1      &
!       ,e4,q,q3)*q1_q3*e3_q4*MG**(-2)*yyy7 + 1./4.*et1(e1,e4,q,q3)*          &
!       e3_q4*yyy7 - et1(e1,e4,q,q4)*q1_q3*e3_q4*MG**(-2)*yyy7 - 1./4.*       &
!       et1(e1,e4,q,q4)*e3_q4*yyy7
! print *, "res old QQB",res



! !     new code with c41 c42 coupligs
!      res =&
!     & q1_e3*e1_e4*yyy1 - q1_e3*e1_q3*e4_q3*yyy41 + q1_e4*e1_e3*yyy1 + &
!     & q1_e4*e1_q3*e3_q4*yyy42 - q1_q3*e1_e3*e4_q3*yyy41 + q1_q3*e1_e4*&
!     & e3_q4*yyy42 + 4.*q1_q3*e1_q3*e3_e4*yyy2 + 4.*q1_q3*e1_q3*e3_q4*&
!     & e4_q3*yyy3 + 1./2.*e1_e3*e4_q3*yyy1 + 1./4.*e1_e3*e4_q3*MZ4**2*&
!     & yyy41 - 1./4.*e1_e3*e4_q3*MZ3**2*yyy41 - 1./4.*e1_e3*e4_q3*MG**2&
!     & *yyy41 + 1./2.*e1_e4*e3_q4*yyy1 - 1./4.*e1_e4*e3_q4*MZ4**2*yyy42&
!     &  + 1./4.*e1_e4*e3_q4*MZ3**2*yyy42 + 1./4.*e1_e4*e3_q4*MG**2*&
!     & yyy42 - e1_q3*e3_e4*MZ4**2*yyy2 + e1_q3*e3_e4*MZ3**2*yyy2 + &
!     & e1_q3*e3_e4*MG**2*yyy2 + 1./2.*e1_q3*e3_q4*e4_q3*yyy42 - 1./2.*&
!     & e1_q3*e3_q4*e4_q3*yyy41 - e1_q3*e3_q4*e4_q3*MZ4**2*yyy3 + e1_q3*&
!     & e3_q4*e4_q3*MZ3**2*yyy3 + e1_q3*e3_q4*e4_q3*MG**2*yyy3 + 1./2.*&
!     & et1(q1,e3,q,q3)*e1_q3*e4_q3*MG**(-2)*yyy7 - 1./2.*et1(q1,e3,q,q4&
!     & )*e1_q3*e4_q3*MG**(-2)*yyy7 + 1./2.*et1(q1,e4,q,q3)*e1_q3*e3_q4*&
!     & MG**(-2)*yyy7 - 1./2.*et1(q1,e4,q,q4)*e1_q3*e3_q4*MG**(-2)*yyy7&
!     &  - 1./2.*et1(q2,e3,q,q3)*e1_q3*e4_q3*MG**(-2)*yyy7
!      res = res + 1./2.*et1(q2,e3,q,q4)*e1_q3*e4_q3*MG**(-2)*yyy7 - 1./&
!     & 2.*et1(q2,e4,q,q3)*e1_q3*e3_q4*MG**(-2)*yyy7 + 1./2.*et1(q2,e4,q&
!     & ,q4)*e1_q3*e3_q4*MG**(-2)*yyy7 + et1(e1,e3,q,q3)*q1_q3*e4_q3*&
!     & MG**(-2)*yyy7 - 1./4.*et1(e1,e3,q,q3)*e4_q3*MG**(-2)*MZ4**2*yyy7&
!     &  + 1./4.*et1(e1,e3,q,q3)*e4_q3*MG**(-2)*MZ3**2*yyy7 + 1./4.*et1(&
!     & e1,e3,q,q3)*e4_q3*yyy7 - et1(e1,e3,q,q4)*q1_q3*e4_q3*MG**(-2)*&
!     & yyy7 + 1./4.*et1(e1,e3,q,q4)*e4_q3*MG**(-2)*MZ4**2*yyy7 - 1./4.*&
!     & et1(e1,e3,q,q4)*e4_q3*MG**(-2)*MZ3**2*yyy7 - 1./4.*et1(e1,e3,q,&
!     & q4)*e4_q3*yyy7 + et1(e1,e4,q,q3)*q1_q3*e3_q4*MG**(-2)*yyy7 - 1./&
!     & 4.*et1(e1,e4,q,q3)*e3_q4*MG**(-2)*MZ4**2*yyy7 + 1./4.*et1(e1,e4,&
!     & q,q3)*e3_q4*MG**(-2)*MZ3**2*yyy7 + 1./4.*et1(e1,e4,q,q3)*e3_q4*&
!     & yyy7 - et1(e1,e4,q,q4)*q1_q3*e3_q4*MG**(-2)*yyy7 + 1./4.*et1(e1,&
!     & e4,q,q4)*e3_q4*MG**(-2)*MZ4**2*yyy7 - 1./4.*et1(e1,e4,q,q4)*&
!     & e3_q4*MG**(-2)*MZ3**2*yyy7 - 1./4.*et1(e1,e4,q,q4)*e3_q4*yyy7 + &
!     & 1./2.*et1(e3,e4,q1,q)*e1_q3*yyy6 - 1./2.*et1(e3,e4,q2,q)*e1_q3*&
!     & yyy6
!      res = res - 1./4.*et1(e3,e4,e1,q)*MZ4**2*yyy6 + 1./4.*et1(e3,e4,&
!     & e1,q)*MZ3**2*yyy6 + 1./4.*et1(e3,e4,e1,q)*MG**2*yyy6 + et1(e3,e4&
!     & ,e1,q)*q1_q3*yyy6 + 4.*et1(e3,e4,q3,q4)*q1_q3*e1_q3*MG**(-2)*&
!     & yyy5 - et1(e3,e4,q3,q4)*e1_q3*MG**(-2)*MZ4**2*yyy5 + et1(e3,e4,&
!     & q3,q4)*e1_q3*MG**(-2)*MZ3**2*yyy5 + et1(e3,e4,q3,q4)*e1_q3*yyy5

! print *, "res new QQB",res



! ! ! ! ! ! ! ! ! ! ! ! ! ! !  SAME AS ABOVE BUT SHORTER  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

        abr1 = (MG**2 + MZ3**2 - MZ4**2 + 4*q1_q3)

         res =  (e1_e3*(4*yyy1*q1_e4 + e4_q3*(2*yyy1 - (MG**2 + MZ3**2 - MZ4**2)*yyy41 - 4*yyy41*q1_q3)) +    &
           e1_e4*(4*yyy1*q1_e3 + e3_q4*(2*yyy1 + (MG**2 + MZ3**2 - MZ4**2)*yyy42 + 4*yyy42*q1_q3)) +    &
           2*e1_q3*(-2*yyy41*e4_q3*q1_e3 + 2*yyy2*e3_e4*abr1 +    &
          e3_q4*(2*yyy42*q1_e4 + e4_q3*(2*(MG**2 + MZ3**2 - MZ4**2)*yyy3 - yyy41 + yyy42 + 8*yyy3*q1_q3))))/4d0 +    &
        (yyy6*abr1*et1(e3,e4,e1,q))/4d0 + (yyy6*e1_q3*et1(e3,e4,q1,q))/2d0    &
       - (yyy6*e1_q3*et1(e3,e4,q2,q))/2d0 +    &
        (yyy5*e1_q3*abr1*et1(e3,e4,q3,q4))/MG**2 +    &
        yyy7*((e4_q3*abr1*et1(e1,e3,q,q3))/(4d0*MG**2) -    &
           (e4_q3*abr1*et1(e1,e3,q,q4))/(4d0*MG**2) +    &
           (e3_q4*abr1*et1(e1,e4,q,q3))/(4d0*MG**2) -    &
           (e3_q4*abr1*et1(e1,e4,q,q4))/(4d0*MG**2)    &
          + (e1_q3*e4_q3*et1(q1,e3,q,q3))/(2d0*MG**2) -    &
           (e1_q3*e4_q3*et1(q1,e3,q,q4))/(2d0*MG**2) + (e1_q3*e3_q4*et1(q1,e4,q,q3))/(2d0*MG**2) -    &
           (e1_q3*e3_q4*et1(q1,e4,q,q4))/(2d0*MG**2) - (e1_q3*e4_q3*et1(q2,e3,q,q3))/(2d0*MG**2) +    &
           (e1_q3*e4_q3*et1(q2,e3,q,q4))/(2d0*MG**2) - (e1_q3*e3_q4*et1(q2,e4,q,q3))/(2d0*MG**2) +    &
           (e1_q3*e3_q4*et1(q2,e4,q,q4))/(2d0*MG**2))

! print *, "res newer QQB ",res
! pause

      end subroutine qqGZZampl


!----- a subroutine for G -> ZZ/WW/AA
!----- all outgoing convention and the following momentum assignment
!-----  0 -> G(p1) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
     subroutine EvalAmp_G_VV(p,MY_IDUP,sum)
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A(2)
      integer :: i1,i2,i3,i4,ordering(1:4)
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


         if( CoupledVertex(MY_IDUP(6:7),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR1=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR1=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR1=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR1=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL1=0d0
                    aR1=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP(6:7),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP(6)) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ud)
                 aR1 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_ln)
                 aR1 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                 aL1 = bL * dsqrt(scale_alpha_W_tn)
                 aR1 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL1=0d0
                 aR1=0d0
              endif
         elseif( MY_IDUP(6).eq.Pho_ ) then !  photon
              aL1=1d0
              aR1=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL1=0d0
              aR1=0d0
         endif

         if( CoupledVertex(MY_IDUP(8:9),-1).eq.Z0_ ) then!  Z decay
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_ll)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_ll)
              elseif( abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=aL_lep    * dsqrt(scale_alpha_Z_tt)
                    aR2=aR_lep    * dsqrt(scale_alpha_Z_tt)
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=aL_neu    * dsqrt(scale_alpha_Z_nn)
                    aR2=aR_neu    * dsqrt(scale_alpha_Z_nn)
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=aL_QUp    * dsqrt(scale_alpha_Z_uu)
                    aR2=aR_QUp    * dsqrt(scale_alpha_Z_uu)
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=aL_QDn    * dsqrt(scale_alpha_Z_dd)
                    aR2=aR_QDn    * dsqrt(scale_alpha_Z_dd)
              else
                    aL2=0d0
                    aR2=0d0
              endif
         elseif( (CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_ .or. CoupledVertex(MY_IDUP(8:9),-1).eq.Wm_) ) then !  W decay
              if( IsAQuark(MY_IDUP(8)) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ud)
                 aR2 = bR * dsqrt(scale_alpha_W_ud)! = 0
              elseif( abs(MY_IDUP(9)).eq.abs(ElM_) .or. abs(MY_IDUP(9)).eq.abs(MuM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_ln)
                 aR2 = bR * dsqrt(scale_alpha_W_ln)! = 0
              elseif( abs(MY_IDUP(9)).eq.abs(TaM_) ) then
                 aL2 = bL * dsqrt(scale_alpha_W_tn)
                 aR2 = bR * dsqrt(scale_alpha_W_tn)! = 0
              else
                 aL2=0d0
                 aR2=0d0
              endif
         elseif( MY_IDUP(8).eq.Pho_ ) then !  photon
              aL2=1d0
              aR2=1d0
              prefactor = prefactor/gZ_sq ! cancel the overall Z/W coupling
         else
              aL2=0d0
              aR2=0d0
         endif


      sum = zero


do i1 =-2,2! G boson
do i3 = 1,2! lepton string1
do i4 = 1,2! lepton string2

         if(CoupledVertex(MY_IDUP(8:9),-1).eq.Wp_ .and. CoupledVertex(MY_IDUP(6:7),-1).eq.Wm_) then
            ordering = (/5,6,3,4/)
         else
            ordering = (/3,4,5,6/)
         endif
         call calcHelAmp2(ordering,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(1))
         if( (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call swap(ordering(2),ordering(4))
             call calcHelAmp2(ordering,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(2))
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

     subroutine calcHelAmp2(ordering,p,idin,i1,i3,i4,A)
     implicit none
     integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4,mu,nu,idin(3:6)
     real(dp) :: p(1:4,1:6)
     complex(dp) :: propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(0:4,1:4)

      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

         pin(1,:) = p(:,1)
         pin(2,:) = 0d0 ! dummy

         sp(0,1:4) = pol_mass2(dcmplx(p(1:4,1)), 0,'in')
         sp(1,1:4) = pol_mass2(dcmplx(p(1:4,1)),-1,'in')
         sp(2,1:4) = pol_mass2(dcmplx(p(1:4,1)),+1,'in')

!-------- -1 == left, 1 == right
         if( idin(l1).ne.Pho_ ) then ! Z, W
            pin(3,:) = p(:,l1)+p(:,l2)
            sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
            sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
            s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
            propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         else !  A
            pin(3,:) = p(:,l1)
            sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
            propZ1=1d0
         endif
         if( idin(l3).ne.Pho_ ) then ! Z, W
            pin(4,:) = p(:,l3)+p(:,l4)
            sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
            sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
            s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
            propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
         else ! A
            pin(4,:) = p(:,l3)
            sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')
            propZ2=1d0
         endif

         call GZZampl(pin,sp,i1,A(1))

         A(1) = A(1) * propZ1*propZ2

      end subroutine

      subroutine GZZampl(p,sp,i1,res)
      implicit none
      real(dp), intent(in) :: p(4,4)
      complex(dp), intent(in) :: sp(0:4,4)
      integer,intent(in) :: i1
      complex(dp), intent(out) :: res
      complex(dp) :: e1_e2, e1_e3, e1_e4
      complex(dp) :: e2_e3, e2_e4
      complex(dp) :: e3_e4
      complex(dp) :: q_q,q3_q3,q4_q4
      complex(dp) :: q1_q2,q1_q3,q1_q4
      complex(dp) :: q2_q3,q2_q4
      complex(dp) :: q3_q4
      complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4,e0_e3,e0_e4
      complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4,e0_q3,e0_q4,q_e3,q_e4
      complex(dp) :: e3_q4,e4_q3
      complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
      complex(dp) :: e1(4),e2(4),e3(4),e4(4),e0(4)
      complex(dp) :: yyy1,yyy2,yyy3,yyy41,yyy42,yyy5,yyy6,yyy7
      real(dp) :: q34,MG,MZ3,MZ4
      real(dp) :: rr
      real(dp),parameter :: sqrt6=dsqrt(6d0)

      q1 = dcmplx(p(1,:),0d0)
      q2 = dcmplx(p(2,:),0d0)
      q3 = dcmplx(p(3,:),0d0)
      q4 = dcmplx(p(4,:),0d0)

!     requirement: sp(0,:)= 0 pol
!     requirement: sp(1,:)= - pol
!     requirement: sp(2,:)= + pol
      if( i1.eq.+2 ) then
        e0 = 1d9! dummy
        e1 = sp(2,:)! +
        e2 = sp(2,:)! +
      elseif( i1.eq.+1 ) then
        e0 = 1d9! dummy
        e1 = sp(2,:)! +
        e2 = sp(0,:)! 0
      elseif( i1.eq.0 ) then
        e0 = sp(0,:)! 0
        e1 = sp(1,:)! -
        e2 = sp(2,:)! +
      elseif( i1.eq.-1 ) then
        e0 = 1d9! dummy
        e1 = sp(1,:)! -
        e2 = sp(0,:)! 0
      elseif( i1.eq.-2 ) then
        e0 = 1d9! dummy
        e1 = sp(1,:)! -
        e2 = sp(1,:)! -
      endif




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

!     new
      e0_e3 = sc(e0,e3)
      e0_e4 = sc(e0,e4)
      e0_q3 = sc(q3,e0)
      e0_q4 = sc(q4,e0)
      q_e3 = sc(q4,e3)
      q_e4 = sc(q3,e4)

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
      if( use_dynamic_MG ) then
          MG = dsqrt(cdabs(q_q))
      else
          MG = M_Reso
      endif

      q34 = (MG**2-MZ3**2-MZ4**2)/2d0

      if (generate_bis) then
          rr = q34/Lambda**2
!           yyy1 = q34*(b1   + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)  + b5*M_V**2/q34)
          yyy1 = q34*(b1 + b2*rr*(one+MZ3**2/q34)*(one+MZ4**2/q34)  + b5*M_V**2/q34)
!           yyy2 = -b1/two + b3*rr*(1d0-M_V**2/q34) + two*b4*rr+b7*rr*M_V**2/q34
          yyy2 = -b1/two + b3*rr*(1d0-(MZ3**2+MZ4**2)/(2d0*q34)) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
!           yyy4 = -b1 - b2*rr -(b2+b3+b6)*rr*M_V**2/q34
          yyy41 = -b1 - b2*(q34+MZ3**2)/Lambda**2 - b3*MZ4**2/Lambda**2 - b6*M_V**2/Lambda**2
          yyy42 = -b1 - b2*(q34+MZ4**2)/Lambda**2 - b3*MZ3**2/Lambda**2 - b6*M_V**2/Lambda**2
          yyy5 = two*b8*rr*MG**2/q34
!           yyy6 = b9
          yyy6 = b9 * M_V**2/Lambda**2
!           yyy7 = b10*rr*MG**2/q34
          yyy7 = b10 * MG**2 * M_V**2/Lambda**4
      else
          yyy1 = (q34)*c1/2d0
          yyy2 = c2
          yyy3 = c3/MG**2
!           yyy4 = c4
          yyy41= c41
          yyy42= c42
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
          if( DecayMode1.eq.7 .and. DecayMode2.eq.7 ) then
              yyy6=0d0
              yyy7=0d0
          endif
      endif


      res = czero


     if( abs(i1).eq.2 ) then

      res =&
     &  + yyy7 * ( et1(e1,e3,q,q3)*q_e4*e2_q3*MG**(-2) - et1(e1,e3,q,q3&
     &    )*q_e4*e2_q4*MG**(-2) - et1(e1,e3,q,q4)*q_e4*e2_q3*MG**(-2)&
     &     + et1(e1,e3,q,q4)*q_e4*e2_q4*MG**(-2) + et1(e1,e4,q,q3)*q_e3&
     &    *e2_q3*MG**(-2) - et1(e1,e4,q,q3)*q_e3*e2_q4*MG**(-2) - et1(&
     &    e1,e4,q,q4)*q_e3*e2_q3*MG**(-2) + et1(e1,e4,q,q4)*q_e3*e2_q4*&
     &    MG**(-2) )
      res = res + yyy6 * ( et1(e3,e4,e1,q)*e2_q3 - et1(e3,e4,e1,q)*&
     &    e2_q4 )
      res = res + yyy5 * ( et1(e3,e4,q3,q4)*e1_q3*e2_q3*MG**(-2) - et1(&
     &    e3,e4,q3,q4)*e1_q3*e2_q4*MG**(-2) - et1(e3,e4,q3,q4)*e1_q4*&
     &    e2_q3*MG**(-2) + et1(e3,e4,q3,q4)*e1_q4*e2_q4*MG**(-2) )
      res = res + yyy42 * ( e1_e4*e2_q3*e3_q4 + e1_q3*e2_e4*e3_q4 )
      res = res + yyy41 * ( e1_e3*e2_q4*e4_q3 + e1_q4*e2_e3*e4_q3 )
      res = res + yyy3 * ( e1_q3*e2_q3*e3_q4*e4_q3 - e1_q3*e2_q4*e3_q4*&
     &    e4_q3 - e1_q4*e2_q3*e3_q4*e4_q3 + e1_q4*e2_q4*e3_q4*e4_q3 )
      res = res + yyy2 * ( e1_q3*e2_q3*e3_e4 - e1_q3*e2_q4*e3_e4 - &
     &    e1_q4*e2_q3*e3_e4 + e1_q4*e2_q4*e3_e4 )
      res = res + yyy1 * ( e1_e3*e2_e4 + e1_e4*e2_e3 )


     elseif( abs(i1).eq.1 ) then

      res =&
     &  + yyy7*sqrt2**(-1) * ( et1(e1,e3,q,q3)*q_e4*e2_q3*MG**(-2) - &
     &    et1(e1,e3,q,q3)*q_e4*e2_q4*MG**(-2) - et1(e1,e3,q,q4)*q_e4*&
     &    e2_q3*MG**(-2) + et1(e1,e3,q,q4)*q_e4*e2_q4*MG**(-2) + et1(e1&
     &    ,e4,q,q3)*q_e3*e2_q3*MG**(-2) - et1(e1,e4,q,q3)*q_e3*e2_q4*&
     &    MG**(-2) - et1(e1,e4,q,q4)*q_e3*e2_q3*MG**(-2) + et1(e1,e4,q,&
     &    q4)*q_e3*e2_q4*MG**(-2) + et1(e2,e3,q,q3)*q_e4*e1_q3*MG**(-2)&
     &     - et1(e2,e3,q,q3)*q_e4*e1_q4*MG**(-2) - et1(e2,e3,q,q4)*q_e4&
     &    *e1_q3*MG**(-2) + et1(e2,e3,q,q4)*q_e4*e1_q4*MG**(-2) + et1(&
     &    e2,e4,q,q3)*q_e3*e1_q3*MG**(-2) - et1(e2,e4,q,q3)*q_e3*e1_q4*&
     &    MG**(-2) - et1(e2,e4,q,q4)*q_e3*e1_q3*MG**(-2) + et1(e2,e4,q,&
     &    q4)*q_e3*e1_q4*MG**(-2) )
      res = res + yyy6*sqrt2**(-1) * ( et1(e3,e4,e1,q)*e2_q3 - et1(e3,&
     &    e4,e1,q)*e2_q4 + et1(e3,e4,e2,q)*e1_q3 - et1(e3,e4,e2,q)*&
     &    e1_q4 )
      res = res + yyy5*sqrt2**(-1) * ( 2.0d0*et1(e3,e4,q3,q4)*e1_q3*e2_q3*&
     &    MG**(-2) - 2.0d0*et1(e3,e4,q3,q4)*e1_q3*e2_q4*MG**(-2) - 2.0d0*et1(&
     &    e3,e4,q3,q4)*e1_q4*e2_q3*MG**(-2) + 2.0d0*et1(e3,e4,q3,q4)*e1_q4&
     &    *e2_q4*MG**(-2) )
      res = res + yyy42*sqrt2**(-1) * ( 2.0d0*e1_e4*e2_q3*e3_q4 + 2.0d0*e1_q3&
     &    *e2_e4*e3_q4 )
      res = res + yyy41*sqrt2**(-1) * ( 2.0d0*e1_e3*e2_q4*e4_q3 + 2.0d0*e1_q4&
     &    *e2_e3*e4_q3 )
      res = res + yyy3*sqrt2**(-1) * ( 2.0d0*e1_q3*e2_q3*e3_q4*e4_q3 - 2.0d0*&
     &    e1_q3*e2_q4*e3_q4*e4_q3 - 2.0d0*e1_q4*e2_q3*e3_q4*e4_q3 + 2.0d0*&
     &    e1_q4*e2_q4*e3_q4*e4_q3 )
      res = res + yyy2*sqrt2**(-1) * ( 2.0d0*e1_q3*e2_q3*e3_e4 - 2.0d0*e1_q3*&
     &    e2_q4*e3_e4 - 2.0d0*e1_q4*e2_q3*e3_e4 + 2.0d0*e1_q4*e2_q4*e3_e4 )
      res = res + yyy1*sqrt2**(-1) * ( 2.0d0*e1_e3*e2_e4 + 2.0d0*e1_e4*e2_e3&
     &     )


     elseif( abs(i1).eq.0 ) then
      res =&
     &  + yyy7*sqrt6**(-1) * (  - 2.0d0*et1(e0,e3,q,q3)*q_e4*e0_q3*&
     &    MG**(-2) + 2.0d0*et1(e0,e3,q,q3)*q_e4*e0_q4*MG**(-2) + 2.0d0*et1(e0&
     &    ,e3,q,q4)*q_e4*e0_q3*MG**(-2) - 2.0d0*et1(e0,e3,q,q4)*q_e4*e0_q4&
     &    *MG**(-2) - 2.0d0*et1(e0,e4,q,q3)*q_e3*e0_q3*MG**(-2) + 2.0d0*et1(&
     &    e0,e4,q,q3)*q_e3*e0_q4*MG**(-2) + 2.0d0*et1(e0,e4,q,q4)*q_e3*&
     &    e0_q3*MG**(-2) - 2.0d0*et1(e0,e4,q,q4)*q_e3*e0_q4*MG**(-2) + &
     &    et1(e1,e3,q,q3)*q_e4*e2_q3*MG**(-2) - et1(e1,e3,q,q3)*q_e4*&
     &    e2_q4*MG**(-2) - et1(e1,e3,q,q4)*q_e4*e2_q3*MG**(-2) + et1(e1&
     &    ,e3,q,q4)*q_e4*e2_q4*MG**(-2) + et1(e1,e4,q,q3)*q_e3*e2_q3*&
     &    MG**(-2) - et1(e1,e4,q,q3)*q_e3*e2_q4*MG**(-2) - et1(e1,e4,q,&
     &    q4)*q_e3*e2_q3*MG**(-2) + et1(e1,e4,q,q4)*q_e3*e2_q4*MG**(-2)&
     &     + et1(e2,e3,q,q3)*q_e4*e1_q3*MG**(-2) - et1(e2,e3,q,q3)*q_e4&
     &    *e1_q4*MG**(-2) - et1(e2,e3,q,q4)*q_e4*e1_q3*MG**(-2) + et1(&
     &    e2,e3,q,q4)*q_e4*e1_q4*MG**(-2) + et1(e2,e4,q,q3)*q_e3*e1_q3*&
     &    MG**(-2) - et1(e2,e4,q,q3)*q_e3*e1_q4*MG**(-2) - et1(e2,e4,q,&
     &    q4)*q_e3*e1_q3*MG**(-2) )
      res = res + yyy7*sqrt6**(-1) * ( et1(e2,e4,q,q4)*q_e3*e1_q4*&
     &    MG**(-2) )
      res = res + yyy6*sqrt6**(-1) * (  - 2.0d0*et1(e3,e4,e0,q)*e0_q3 + 2.0d0&
     &    *et1(e3,e4,e0,q)*e0_q4 + et1(e3,e4,e1,q)*e2_q3 - et1(e3,e4,e1&
     &    ,q)*e2_q4 + et1(e3,e4,e2,q)*e1_q3 - et1(e3,e4,e2,q)*e1_q4 )
      res = res + yyy5*sqrt6**(-1) * ( 4.0d0*et1(e3,e4,q3,q4)*e0_q3*e0_q4*&
     &    MG**(-2) - 2.0d0*et1(e3,e4,q3,q4)*e0_q3**2*MG**(-2) - 2.0d0*et1(e3,&
     &    e4,q3,q4)*e0_q4**2*MG**(-2) + 2.0d0*et1(e3,e4,q3,q4)*e1_q3*e2_q3&
     &    *MG**(-2) - 2.0d0*et1(e3,e4,q3,q4)*e1_q3*e2_q4*MG**(-2) - 2.0d0*&
     &    et1(e3,e4,q3,q4)*e1_q4*e2_q3*MG**(-2) + 2.0d0*et1(e3,e4,q3,q4)*&
     &    e1_q4*e2_q4*MG**(-2) )
      res = res + yyy42*sqrt6**(-1) * (  - 4.0d0*e0_e4*e0_q3*e3_q4 + 2.0d0*&
     &    e1_e4*e2_q3*e3_q4 + 2.0d0*e1_q3*e2_e4*e3_q4 )
      res = res + yyy41*sqrt6**(-1) * (  - 4.0d0*e0_e3*e0_q4*e4_q3 + 2.0d0*&
     &    e1_e3*e2_q4*e4_q3 + 2.0d0*e1_q4*e2_e3*e4_q3 )
      res = res + yyy3*sqrt6**(-1) * ( 4.0d0*e0_q3*e0_q4*e3_q4*e4_q3 - 2.0d0*&
     &    e0_q3**2*e3_q4*e4_q3 - 2.0d0*e0_q4**2*e3_q4*e4_q3 + 2.0d0*e1_q3*&
     &    e2_q3*e3_q4*e4_q3 - 2.0d0*e1_q3*e2_q4*e3_q4*e4_q3 - 2.0d0*e1_q4*&
     &    e2_q3*e3_q4*e4_q3 + 2.0d0*e1_q4*e2_q4*e3_q4*e4_q3 )
      res = res + yyy2*sqrt6**(-1) * ( 4.0d0*e0_q3*e0_q4*e3_e4 - 2.0d0*&
     &    e0_q3**2*e3_e4 - 2.0d0*e0_q4**2*e3_e4 + 2.0d0*e1_q3*e2_q3*e3_e4 - 2.0d0&
     &    *e1_q3*e2_q4*e3_e4 - 2.0d0*e1_q4*e2_q3*e3_e4 + 2.0d0*e1_q4*e2_q4*&
     &    e3_e4 )
      res = res + yyy1*sqrt6**(-1) * (  - 4.0d0*e0_e3*e0_e4 + 2.0d0*e1_e3*&
     &    e2_e4 + 2.0d0*e1_e4*e2_e3 )


     endif

      end subroutine GZZampl


!----- a subroutine for q qbar -> G -> Z -> lept + Z --> 2 lepts
!----- all outgoing convention and the following momentum assignment
!-----  0 -> bq(p1) + q(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
     subroutine EvalAmp_qqb_G_VV_old(p,sum)
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      real(dp) :: s, pin(4,4)
      complex(dp) :: A(2), sp(4,4), propG, propZ1, propZ2
      integer :: i1,i2,i3,i4
      real(dp) :: aL,aR
      real(dp) :: gZ_sq
      real(dp) :: prefactor


!---- electroweak couplings

      aL = -one + two*sitW**2
      aR = aL+one

      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2


!---- full prefactor; 3 is  the color factor

      prefactor = 3d0*gZ_sq**2

      sum = zero

      s  = two*scr(p(:,1),p(:,2))
      propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)

      s = two*scr(p(:,3),p(:,4))
      propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)

      s = two*scr(p(:,5),p(:,6))
      propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)

            do i1=1,2
                 do i3 = 1,2
                    do i4 = 1,2


         pin(1,:) = p(:,1)
         pin(2,:) = p(:,2)

         sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1)  !qbq
         sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform
                            ! bookeeping

         pin(3,:) = p(:,3) + p(:,4)
         pin(4,:) = p(:,5) + p(:,6)

!-------- -1 == left, 1 == right

         sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
         sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-.mu+

         call qqGZZampl_old(pin,sp,A(1))

         if (i3.eq.1) then
            A(1) = aL*A(1)
          elseif(i3.eq.2) then
            A(1) = aR*A(1)
         endif

         if (i4.eq.1) then
            A(1) = aL*A(1)
          elseif(i4.eq.2) then
            A(1) = aR*A(1)
         endif



          sum = sum + abs(propG*propZ1*propZ2*A(1))**2


                            enddo
                         enddo
                     enddo


                sum = sum*prefactor

      end subroutine

      subroutine qqGZZampl_old(p,sp,res)
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


      res =                                                            &
       + 1./4.*q1_e3*q_q*e1_e4 - 1./2.*q1_e3*e1_q4*e4_q3 + 1./4.*q1_e4  &
         *q_q*e1_e3 - 1./2.*q1_e4*e1_q3*e3_q4 - 1./2.*q1_q3*e1_e4*      &
         e3_q4 + 1./2.*q1_q3*e1_q4*e3_e4 - 1./2.*q1_q4*e1_e3*e4_q3 + 1. &
         /2.*q1_q4*e1_q3*e3_e4 - 1./4.*q2_e3*q_q*e1_e4 + 1./2.*q2_e3*   &
         e1_q4*e4_q3 - 1./4.*q2_e4*q_q*e1_e3 + 1./2.*q2_e4*e1_q3*e3_q4  &
          + 1./2.*q2_q3*e1_e4*e3_q4 - 1./2.*q2_q3*e1_q4*e3_e4 + 1./2.*  &
         q2_q4*e1_e3*e4_q3 - 1./2.*q2_q4*e1_q3*e3_e4


      end subroutine qqGZZampl_old


!----- a subroutinefor gg -> G -> Z -> lept + Z --> 2 lepts
!----- all outgoing convention and the following momentum assignment
!-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_gg_G_VV_old(p,sum)
      use ModMisc
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      real(dp) :: s, pin(4,4)
      complex(dp) :: A(2), sp(4,4), propG, propZ1, propZ2
      integer :: i1,i2,i3,i4
      real(dp) :: aL,aR
      real(dp) :: gZ_sq
      real(dp) :: prefactor

!---- electroweak couplings
      aL = -one + two*sitW**2
      aR = aL+one


      gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2

!---- full prefactor; 8 is  the color factor
      prefactor = 8d0*gZ_sq**2


      sum = zero

      s  = 2d0 * scr(p(:,1),p(:,2))
      propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)

      s = 2d0 * scr(p(:,3),p(:,4))
      propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)

      s = 2d0 * scr(p(:,5),p(:,6))
      propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)

            do i1=1,2
              do i2 = 1,2
                 do i3 = 1,2
                    do i4 = 1,2


         pin(1,:) = p(:,1)
         pin(2,:) = p(:,2)

         sp(1,:) = pol_mless2(dcmplx(p(:,1)),-3+2*i1,'in')  ! gluon
         sp(2,:) = pol_mless2(dcmplx(p(:,2)),-3+2*i2,'in')  ! gluon

         pin(3,:) = p(:,3) + p(:,4)
         pin(4,:) = p(:,5)+p(:,6)

!-------- -1 == left, 1 == right
         sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
         sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-.mu+


         call ggGZZampl_old(pin,sp,A(1))

         if (i3.eq.1) then
            A(1) = aL*A(1)
          elseif(i3.eq.2) then
            A(1) = aR*A(1)
         endif

         if (i4.eq.1) then
            A(1) = aL*A(1)
          elseif(i4.eq.2) then
            A(1) = aR*A(1)
         endif


          sum = sum + abs(propG*propZ1*propZ2*A(1))**2


                            enddo
                         enddo
                       enddo
                     enddo

                sum = sum*prefactor

      end subroutine

      subroutine ggGZZampl_old(p,sp,res)
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


      res =  &
       + M_Reso**(-2) * (  - q1_q3*q_q*e1_e2*e3_q4*e4_q3 + 1./2.*q1_q3*    &
       q_q**2*e1_e2*e3_e4 - q1_q4*q_q*e1_e2*e3_q4*e4_q3 + 1./2.*       &
       q1_q4*q_q**2*e1_e2*e3_e4 - q2_q3*q_q*e1_e2*e3_q4*e4_q3 + 1./2.  &
       *q2_q3*q_q**2*e1_e2*e3_e4 - q2_q4*q_q*e1_e2*e3_q4*e4_q3 + 1./   &
       2.*q2_q4*q_q**2*e1_e2*e3_e4 - q_q**2*e1_e2*e3_q4*e4_q3 + 1./2.  &
       *q_q**3*e1_e2*e3_e4 )

      res = res + q1_e3*q2_e4*q_q*e1_e2 - 2*q1_e3*q2_q4*e1_e2*e4_q3   &
       + q1_e4*q2_e3*q_q*e1_e2 - 2*q1_e4*q2_q3*e1_e2*e3_q4 - 2*       &
       q1_q3*q2_e4*e1_e2*e3_q4 + 2*q1_q3*q2_q4*e1_e2*e3_e4 - 2*q1_q4  &
       *q2_e3*e1_e2*e4_q3 + 2*q1_q4*q2_q3*e1_e2*e3_e4 + q_q*e1_e2*    &
       e3_e4*M_V**2 + 2*q_q*e1_e2*e3_q4*e4_q3 - q_q*e1_e3*e2_q4*e4_q3  &
       - q_q*e1_e4*e2_q3*e3_q4 - q_q*e1_q3*e2_e4*e3_q4 + q_q*e1_q3*   &
       e2_q4*e3_e4 - q_q*e1_q4*e2_e3*e4_q3 + q_q*e1_q4*e2_q3*e3_e4    &
       - q_q**2*e1_e2*e3_e4 + 1./2.*q_q**2*e1_e3*e2_e4 + 1./2.*       &
       q_q**2*e1_e4*e2_e3

      end subroutine ggGZZampl_old


       end module


