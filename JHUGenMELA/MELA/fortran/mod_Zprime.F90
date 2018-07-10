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
      complex(dp) :: A(2)
      integer :: i1,i2,i3,i4,ordering(1:4),ordering_swap(1:4),idV(1:2),VVmode
      real(dp) :: prefactor
      real(dp) :: intcolfac

      if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
         intcolfac=1.0_dp/3.0_dp
      else
         intcolfac=1.0_dp
      endif

      call getDecay_VVMode_Ordering(MY_IDUP(6:9),VVMode,ordering,ordering_swap)

!---- full prefactor; 3 is  the color factor
      if( VVMode.eq.ZZMode ) then!  Z decay
         prefactor = 3d0*overallCouplVffsq**2
      elseif( VVMode.eq.WWMode ) then !  W decay
         prefactor = 3d0*overallCouplVffsq**2
      elseif( VVMode.eq.ZgMode ) then !  Z+photon "decay"
         prefactor = 3d0*overallCouplVffsq ! Only single powers
      elseif( VVMode.eq.ggMode ) then !  photon "decay"
         prefactor = 3d0
      else
         prefactor = 0d0
      endif


      sum = zero
do i1=1,2
do i3 = 1,2
do i4 = 1,2
!do i3 = -1,1! on-shell check!
!do i4 = -1,1! on-shell check!

         call calcHelAmp(ordering,VVMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(1))
         if( includeInterference .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call calcHelAmp(ordering_swap,VVMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(2))
             A(2) = -A(2) ! minus comes from fermi statistics
         endif
         if( includeInterference .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
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

   subroutine calcHelAmp(ordering,VVmode,p,MY_IDUP,i1,i3,i4,A)
      implicit none
      integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4,VVmode,MY_IDUP(6:9)
      real(dp) :: p(1:4,1:6)
      complex(dp) :: propG
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
!-------- For fermions, -1 == left, 1 == right
      pin(1,:) = p(:,1)
      pin(2,:) = p(:,2)
      sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1) !qbq
      sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform bookeeping
!-- on-shell check
!       sp(3,1:4) = pol_mass(dcmplx(p(1:4,3)+p(1:4,4)),dsqrt(2d0*scr(p(:,3),p(:,4))),i3)
!       sp(4,1:4) = pol_mass(dcmplx(p(1:4,5)+p(1:4,6)),dsqrt(2d0*scr(p(:,5),p(:,6))),i4)

!        if( i3.eq.0 .and. i4.eq.0 ) then
!             a(:) = 0d0
!             return
!        endif
!-- end: on-shell check
      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call qqZprimeZZampl(pin,sp,A(1))
!---- chiral couplings of quarks to Zprimes
      if (i1.eq.1) then
         A(1) = zprime_qq_left*A(1)
      else if(i1.eq.2) then
         A(1) = zprime_qq_right*A(1)
      endif
      A(1) = A(1)*propG
! A(1) ~ GeV^-2
! propG*propZ1*propZ2 ~  GeV^-2
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
     subroutine EvalAmp_Zprime_VV(p,MY_IDUP,sum)
      implicit none
      real(dp), intent(out) ::  sum
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      real(dp) :: pin(4,4)
      complex(dp) :: A(2)
      integer :: i1,i2,i3,i4,ordering(1:4),ordering_swap(1:4),idV(1:2),VVmode
      real(dp) :: aL1,aR1,aL2,aR2
      real(dp) :: prefactor
      real(dp) :: intcolfac

      if(IsAQuark(MY_IDUP(6)) .and. IsAQuark(MY_IDUP(8))) then
         intcolfac=1.0_dp/3.0_dp
      else
         intcolfac=1.0_dp
      endif

      call getDecay_VVMode_Ordering(MY_IDUP(6:9),VVMode,ordering,ordering_swap)

!---- full prefactor
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


      sum = zero

do i1 =-1,1! Z' boson
do i3 = 1,2! lepton string1
do i4 = 1,2! lepton string2

         call calcHelAmp2(ordering,VVMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(1))
         if( includeInterference .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
             call calcHelAmp2(ordering_swap,VVMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A(2))
             A(2) = -A(2) ! minus comes from fermi statistics
         endif
         if( includeInterference .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
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

   subroutine calcHelAmp2(ordering,VVMode,p,MY_IDUP,i1,i3,i4,A)
      implicit none
      integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4,VVMode,MY_IDUP(6:9)
      real(dp) :: p(1:4,1:6)
      real(dp) :: pin(4,4)
      complex(dp) :: A(1:1), sp(4,4)

      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      pin(1,:) = p(:,1)
      sp(1,:) = pol_mass2(dcmplx(p(1:4,1)),i1,'in')
      pin(2,:) = 0d0  ! dummy
      sp(2,:)  = czero
      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call ZprimeZZampl(pin,sp,A(1))
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


subroutine getDecay_Couplings_Spinors_Props(VVMode,idordered,pordered,h3,h4, sp,pV)
   implicit none
   integer, intent(in) :: VVMode,idordered(6:9),h3,h4
   real(dp), intent(in) :: pordered(1:4,6:9)
   complex(dp), intent(out) :: sp(3:4,1:4)
   real(dp), intent(out) :: pV(3:4,1:4)
   real(dp) :: s, aL1,aR1,aL2,aR2
   complex(dp) :: propV(1:2)

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

   else
      call Error("Unsupported decay modes")
   endif

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

subroutine getDecay_VVMode_Ordering(MY_IDUP, VVMode,ordering,ordering_swap)
   implicit none
   integer, intent(in) :: MY_IDUP(6:9)
   integer, intent(out) :: VVMode,ordering(1:4),ordering_swap(1:4)
   integer :: idV(1:2)

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
      call Error("Unsupported decay Modes")
   endif
   return
end subroutine


       end module


