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
      subroutine EvalAmp_gg_G_VV(p,MY_IDUP,res)
      implicit none
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A_VV(1:18), A0_VV(1:2)
      integer :: i1,i2,i3,i4,VVMode,VVmode_swap
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


      res = zero
      A_VV(:) = 0d0
      doInterference = includeInterference .and. (          &
         ((VVMode.eq.ZZMode) .and. (VVMode_swap.eq.ZZMode)) &
         )
      if ( includeVprime .and. .not.(VVMode.eq.ZZMode .or. VVMode.eq.ZgMode .or. VVMode.eq.WWMode) ) then
         call Error("Contact terms only for WW, ZZ or Zg!")
      endif
      do i1=1,2;  do i2=1,2;  do i3=1,2;  do i4=1,2!  sum over helicities
         call calcHelAmp_gg(ordering,VVMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(1))
         if( VVMode.eq.ZZMode ) then
            if( includeGammaStar ) then
               call calcHelAmp_gg(ordering,ZgsMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(3))
               call calcHelAmp_gg(ordering,gsZMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(5))
               call calcHelAmp_gg(ordering,gsgsMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(7))
            endif
            if( includeVprime ) then
               call calcHelAmp_gg(ordering,ZZpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(9))
               call calcHelAmp_gg(ordering,ZpZMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(11))
               call calcHelAmp_gg(ordering,ZpZpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(13))
            endif
            if( includeGammaStar .and. includeVprime ) then
               call calcHelAmp_gg(ordering,gsZpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(15))
               call calcHelAmp_gg(ordering,ZpgsMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(17))
            endif
         elseif( VVMode.eq.ZgMode ) then
            if(includeGammaStar) then
               call calcHelAmp_gg(ordering,gsgMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(3))
            endif
            if( includeVprime ) then
               call calcHelAmp_gg(ordering,ZpgMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(5))
            endif
         elseif( VVMode.eq.WWMode .and. includeVprime ) then
            call calcHelAmp_gg(ordering,WWpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(9))
            call calcHelAmp_gg(ordering,WpWMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(11))
            call calcHelAmp_gg(ordering,WpWpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(13))
         endif

         if( doInterference ) then
            call calcHelAmp_gg(ordering_swap,VVMode_swap,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(2))
            if( includeGammaStar ) then
               call calcHelAmp_gg(ordering_swap,ZgsMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(4))
               call calcHelAmp_gg(ordering_swap,gsZMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(6))
               call calcHelAmp_gg(ordering_swap,gsgsMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(8))
            endif
            if( includeVprime ) then
               call calcHelAmp_gg(ordering_swap,ZZpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(10))
               call calcHelAmp_gg(ordering_swap,ZpZMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(12))
               call calcHelAmp_gg(ordering_swap,ZpZpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(14))
            endif
            if( includeGammaStar .and. includeVprime ) then
               call calcHelAmp_gg(ordering_swap,gsZpMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(16))
               call calcHelAmp_gg(ordering_swap,ZpgsMode,p(1:4,1:6),MY_IDUP,i1,i2,i3,i4,A_VV(18))
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

      res = res*prefactor
      if( (VVMode.eq.ZZMode) .and. doInterference ) res = res * SymmFac

      end subroutine

   subroutine calcHelAmp_gg(ordering,VVMode,p,MY_IDUP,i1,i2,i3,i4,A)
      implicit none
      integer :: ordering(1:4),i1,i2,i3,i4,l1,l2,l3,l4,MY_IDUP(6:9),VVMode
      real(dp) :: p(1:4,1:6)
      complex(dp) :: propG
      real(dp) :: pin(4,4)
      complex(dp) :: A(1:1), sp(4,4)


      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)

      propG = one/dcmplx(2d0*scr(p(:,1),p(:,2)) - M_Reso**2,M_Reso*Ga_Reso)

      pin(1,:) = p(:,1)
      pin(2,:) = p(:,2)
      sp(1,:) = pol_mless2(dcmplx(p(:,1)),-3+2*i1,'in')  ! gluon
      sp(2,:) = pol_mless2(dcmplx(p(:,2)),-3+2*i2,'in')  ! gluon
      !sp(1,1:4)=pin(1,1:4);print *, "this checks IS gauge invariance"
      !sp(2,1:4)=pin(2,1:4);print *, "this checks IS gauge invariance"
      !careful: for gluon gauge invariance check the terms ~c3,c4 are needed because e1.q2 is not zero for e1-->q1
      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call ggGZZampl(VVMode,pin,sp,A(1))
      A(1) = A(1)*propG

   end subroutine

      subroutine ggGZZampl(VVMode,p,sp,res)
      implicit none
      integer, intent(in) :: VVMode
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
      complex(dp) :: b_dyn(1:10)
      real(dp) :: q34,MZ3,MZ4,MG
      logical :: new
      real(dp) :: rr_gam, rr


      new = .true.

      b_dyn(:)=czero

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

         if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
            b_dyn(1)=b1
            b_dyn(2)=b2
            b_dyn(3)=b3
            b_dyn(4)=b4
            b_dyn(5)=b5
            b_dyn(6)=b6
            b_dyn(7)=b7
            b_dyn(8)=b8
            b_dyn(9)=b9
            b_dyn(10)=b10
         elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then
            b_dyn(1)=bzgs1
            b_dyn(2)=bzgs2
            b_dyn(3)=bzgs3
            b_dyn(4)=bzgs4
            b_dyn(8)=bzgs8
         elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then
            b_dyn(1)=bgsgs1
            b_dyn(2)=bgsgs2
            b_dyn(3)=bgsgs3
            b_dyn(4)=bgsgs4
            b_dyn(8)=bgsgs8
         elseif( (VVMode.eq.ZZpMode) .or. (VVMode.eq.WWpMode) .or. (VVMode.eq.ZpZMode) .or. (VVMode.eq.WpWMode) ) then
            b_dyn(1)=bzzp1
            b_dyn(2)=bzzp2
            b_dyn(3)=bzzp3
            b_dyn(4)=bzzp4
            b_dyn(5)=bzzp5
            b_dyn(6)=bzzp6
            b_dyn(7)=bzzp7
            b_dyn(8)=bzzp8
            b_dyn(9)=bzzp9
            b_dyn(10)=bzzp10
         elseif( (VVMode.eq.ZpZpMode) .or. (VVMode.eq.WpWpMode) ) then
            b_dyn(1)=bzpzp1
            b_dyn(2)=bzpzp2
            b_dyn(3)=bzpzp3
            b_dyn(4)=bzpzp4
            b_dyn(5)=bzpzp5
            b_dyn(6)=bzpzp6
            b_dyn(7)=bzpzp7
            b_dyn(8)=bzpzp8
            b_dyn(9)=bzpzp9
            b_dyn(10)=bzpzp10
         elseif( (VVMode.eq.ZpgMode) .OR. (VVMode.eq.gsZpMode) .OR. (VVMode.eq.ZpgsMode) ) then
            b_dyn(1)=bzpgs1
            b_dyn(2)=bzpgs2
            b_dyn(3)=bzpgs3
            b_dyn(4)=bzpgs4
            b_dyn(8)=bzpgs8
         else
            print *,"VVMode",VVMode,"not implemented"
         endif

          yyy1 = q34*( b_dyn(1) + b_dyn(2)*rr*(one+MZ3**2/q34)*(one+MZ4**2/q34) ) + b_dyn(5)*M_V**2
          yyy2 = -b_dyn(1)/two + b_dyn(3)*rr*(1d0-(MZ3**2+MZ4**2)/(2d0*q34)) + two*b_dyn(4)*rr + b_dyn(7)*rr*M_V**2/q34
          yyy3 = (-b_dyn(2)/two - b_dyn(3)- two*b_dyn(4))*rr/q34
          yyy41 = -b_dyn(1) - b_dyn(2)*(q34+MZ3**2)/Lambda**2 - b_dyn(3)*MZ4**2/Lambda**2 - 2d0*b_dyn(6)*M_V**2/Lambda**2
          yyy42 = -b_dyn(1) - b_dyn(2)*(q34+MZ4**2)/Lambda**2 - b_dyn(3)*MZ3**2/Lambda**2 - 2d0*b_dyn(6)*M_V**2/Lambda**2
          yyy5 = two*b_dyn(8)*rr*MG**2/q34
          yyy6 = b_dyn(9) * M_V**2/Lambda**2
          yyy7 = b_dyn(10) * MG**2 * M_V**2/Lambda**4

      else
          yyy1 = q34*c1/2d0
          yyy2 = c2
          yyy3 = c3/MG**2
          yyy41 = c41
          yyy42 = c42
          yyy5 = c5
          yyy6 = czero
          yyy7 = czero
          if(VVMode.eq.ZZMode .or. VVMode.eq.WWMode) then
             yyy6 = c6
             yyy7 = c7
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
     subroutine EvalAmp_qqb_G_VV(p,MY_IDUP,res)
      implicit none
      real(dp), intent(out) :: res
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      real(dp) ::  pin(4,4)
      complex(dp) :: A_VV(1:18), A0_VV(1:2)
      integer :: i1,i3,i4,VVMode,VVmode_swap
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


      res = zero
      A_VV(:) = 0d0
      doInterference = includeInterference .and. (          &
         ((VVMode.eq.ZZMode) .and. (VVMode_swap.eq.ZZMode)) &
         )
      if ( includeVprime .and. .not.(VVMode.eq.ZZMode .or. VVMode.eq.ZgMode .or. VVMode.eq.WWMode) ) then
         call Error("Contact terms only for WW, ZZ or Zg!")
      endif
      do i1=1,2;  do i3=1,2;  do i4=1,2!  sum over helicities
         call calcHelAmp_qq(ordering,VVMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(1))
         if( VVMode.eq.ZZMode ) then
            if( includeGammaStar ) then
               call calcHelAmp_qq(ordering,ZgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(3))
               call calcHelAmp_qq(ordering,gsZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(5))
               call calcHelAmp_qq(ordering,gsgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(7))
            endif
            if( includeVprime ) then
               call calcHelAmp_qq(ordering,ZZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(9))
               call calcHelAmp_qq(ordering,ZpZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(11))
               call calcHelAmp_qq(ordering,ZpZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(13))
            endif
            if( includeGammaStar .and. includeVprime ) then
               call calcHelAmp_qq(ordering,gsZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(15))
               call calcHelAmp_qq(ordering,ZpgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(17))
            endif
         elseif( VVMode.eq.ZgMode ) then
            if(includeGammaStar) then
               call calcHelAmp_qq(ordering,gsgMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(3))
            endif
            if( includeVprime ) then
               call calcHelAmp_qq(ordering,ZpgMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(5))
            endif
         elseif( VVMode.eq.WWMode .and. includeVprime ) then
            call calcHelAmp_qq(ordering,WWpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(9))
            call calcHelAmp_qq(ordering,WpWMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(11))
            call calcHelAmp_qq(ordering,WpWpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(13))
         endif

         if( doInterference ) then
            call calcHelAmp_qq(ordering_swap,VVMode_swap,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(2))
            if( includeGammaStar ) then
               call calcHelAmp_qq(ordering_swap,ZgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(4))
               call calcHelAmp_qq(ordering_swap,gsZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(6))
               call calcHelAmp_qq(ordering_swap,gsgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(8))
            endif
            if( includeVprime ) then
               call calcHelAmp_qq(ordering_swap,ZZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(10))
               call calcHelAmp_qq(ordering_swap,ZpZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(12))
               call calcHelAmp_qq(ordering_swap,ZpZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(14))
            endif
            if( includeGammaStar .and. includeVprime ) then
               call calcHelAmp_qq(ordering_swap,gsZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(16))
               call calcHelAmp_qq(ordering_swap,ZpgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(18))
            endif
         endif

         A0_VV(1) = A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7)+A_VV(9)+A_VV(11)+A_VV(13)+A_VV(15)+A_VV(17) ! 3456 pieces
         A0_VV(2) = A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8)+A_VV(10)+A_VV(12)+A_VV(14)+A_VV(16)+A_VV(18) ! 5436 pieces
         res = res + dreal(A0_VV(1)*dconjg(A0_VV(1)))
         res = res + dreal(A0_VV(2)*dconjg(A0_VV(2)))
         if( doInterference .and. (i3.eq.i4) ) then! interfere the 3456 with 5436 pieces
             res = res - 2d0*intcolfac*dreal(  A0_VV(1)*dconjg(A0_VV(2))  ) ! minus from Fermi statistics
         endif
      enddo;  enddo;  enddo

      res = res*prefactor
      if( (VVMode.eq.ZZMode) .and. doInterference ) res = res * SymmFac

      end subroutine

   subroutine calcHelAmp_qq(ordering,VVMode,p,MY_IDUP,i1,i3,i4,A)
      implicit none
      integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4,MY_IDUP(6:9),VVMode
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

      s = 2d0*scr(p(:,1),p(:,2))
      propG = s/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)
      pin(1,:) = p(:,1)
      pin(2,:) = p(:,2)
      sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1)  !qbq
      sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform bookeeping
      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call qqGZZampl(VVMode,pin,sp,A(1))
!---- chiral couplings of quarks to gravitons
      if (i1.eq.1) then
         A(1) = graviton_qq_left*A(1)
      elseif(i1.eq.2) then
         A(1) = graviton_qq_right*A(1)
      endif
      A(1) = A(1)*propG
   end subroutine

      subroutine qqGZZampl(VVMode,p,sp,res)
      implicit none
      integer, intent(in) :: VVMode
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
      complex(dp) :: b_dyn(1:10)
      real(dp) :: q34,MG,MZ3,MZ4
      real(dp) :: rr

      b_dyn(:)=czero

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
          rr = q34/Lambda**2! kappa for FS

         if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
            b_dyn(1)=b1
            b_dyn(2)=b2
            b_dyn(3)=b3
            b_dyn(4)=b4
            b_dyn(5)=b5
            b_dyn(6)=b6
            b_dyn(7)=b7
            b_dyn(8)=b8
            b_dyn(9)=b9
            b_dyn(10)=b10
         elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then
            b_dyn(1)=bzgs1
            b_dyn(2)=bzgs2
            b_dyn(3)=bzgs3
            b_dyn(4)=bzgs4
            b_dyn(8)=bzgs8
         elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then
            b_dyn(1)=bgsgs1
            b_dyn(2)=bgsgs2
            b_dyn(3)=bgsgs3
            b_dyn(4)=bgsgs4
            b_dyn(8)=bgsgs8
         elseif( (VVMode.eq.ZZpMode) .or. (VVMode.eq.WWpMode) .or. (VVMode.eq.ZpZMode) .or. (VVMode.eq.WpWMode) ) then
            b_dyn(1)=bzzp1
            b_dyn(2)=bzzp2
            b_dyn(3)=bzzp3
            b_dyn(4)=bzzp4
            b_dyn(5)=bzzp5
            b_dyn(6)=bzzp6
            b_dyn(7)=bzzp7
            b_dyn(8)=bzzp8
            b_dyn(9)=bzzp9
            b_dyn(10)=bzzp10
         elseif( (VVMode.eq.ZpZpMode) .or. (VVMode.eq.WpWpMode) ) then
            b_dyn(1)=bzpzp1
            b_dyn(2)=bzpzp2
            b_dyn(3)=bzpzp3
            b_dyn(4)=bzpzp4
            b_dyn(5)=bzpzp5
            b_dyn(6)=bzpzp6
            b_dyn(7)=bzpzp7
            b_dyn(8)=bzpzp8
            b_dyn(9)=bzpzp9
            b_dyn(10)=bzpzp10
         elseif( (VVMode.eq.ZpgMode) .OR. (VVMode.eq.gsZpMode) .OR. (VVMode.eq.ZpgsMode) ) then
            b_dyn(1)=bzpgs1
            b_dyn(2)=bzpgs2
            b_dyn(3)=bzpgs3
            b_dyn(4)=bzpgs4
            b_dyn(8)=bzpgs8
         else
            print *,"VVMode",VVMode,"not implemented"
         endif

          yyy1 = q34*( b_dyn(1) + b_dyn(2)*rr*(one+MZ3**2/q34)*(one+MZ4**2/q34) ) + b_dyn(5)*M_V**2
          yyy2 = -b_dyn(1)/two + b_dyn(3)*rr*(1d0-(MZ3**2+MZ4**2)/(2d0*q34)) + two*b_dyn(4)*rr + b_dyn(7)*rr*M_V**2/q34
          yyy3 = (-b_dyn(2)/two - b_dyn(3)- two*b_dyn(4))*rr/q34
          yyy41 = -b_dyn(1) - b_dyn(2)*(q34+MZ3**2)/Lambda**2 - b_dyn(3)*MZ4**2/Lambda**2 - 2d0*b_dyn(6)*M_V**2/Lambda**2
          yyy42 = -b_dyn(1) - b_dyn(2)*(q34+MZ4**2)/Lambda**2 - b_dyn(3)*MZ3**2/Lambda**2 - 2d0*b_dyn(6)*M_V**2/Lambda**2
          yyy5 = two*b_dyn(8)*rr*MG**2/q34
          yyy6 = b_dyn(9) * M_V**2/Lambda**2
          yyy7 = b_dyn(10) * MG**2 * M_V**2/Lambda**4
      else
          yyy1 = q34*c1/2d0
          yyy2 = c2
          yyy3 = c3/MG**2
          yyy41 = c41
          yyy42 = c42
          yyy5 = c5
          yyy6 = czero
          yyy7 = czero
          if(VVMode.eq.ZZMode .or. VVMode.eq.WWMode) then
             yyy6 = c6
             yyy7 = c7
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
     subroutine EvalAmp_G_VV(p,MY_IDUP,res)
      implicit none
      real(dp), intent(out) :: res
      real(dp), intent(in) :: p(4,6)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A_VV(1:18),A0_VV(1:2)
      integer :: i1,i3,i4,VVMode,VVmode_swap
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

      res = zero
      A_VV(:) = 0d0
      doInterference = includeInterference .and. (         &
         ((VVMode.eq.ZZMode) .and. (VVMode_swap.eq.ZZMode)) &
         )
      if ( includeVprime .and. .not.(VVMode.eq.ZZMode .or. VVMode.eq.ZgMode .or. VVMode.eq.WWMode) ) then
        call Error("Contact terms only for WW, ZZ or Zg!")
      endif
      do i1 =-2,2;  do i3=1,2;  do i4=1,2!  sum over helicities
         call calcHelAmp2(ordering,VVMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(1))
         if( VVMode.eq.ZZMode ) then
            if( includeGammaStar ) then
               call calcHelAmp2(ordering,ZgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(3))
               call calcHelAmp2(ordering,gsZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(5))
               call calcHelAmp2(ordering,gsgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(7))
            endif
            if( includeVprime ) then
               call calcHelAmp2(ordering,ZZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(9))
               call calcHelAmp2(ordering,ZpZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(11))
               call calcHelAmp2(ordering,ZpZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(13))
            endif
            if( includeGammaStar .and. includeVprime ) then
               call calcHelAmp2(ordering,gsZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(15))
               call calcHelAmp2(ordering,ZpgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(17))
            endif
         elseif( VVMode.eq.ZgMode ) then
            if(includeGammaStar) then
               call calcHelAmp2(ordering,gsgMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(3))
            endif
            if( includeVprime ) then
               call calcHelAmp2(ordering,ZpgMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(5))
            endif
         elseif( VVMode.eq.WWMode .and. includeVprime ) then
            call calcHelAmp2(ordering,WWpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(9))
            call calcHelAmp2(ordering,WpWMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(11))
            call calcHelAmp2(ordering,WpWpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(13))
         endif

         if( doInterference ) then
            call calcHelAmp2(ordering_swap,VVMode_swap,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(2))
            if( includeGammaStar ) then
               call calcHelAmp2(ordering_swap,ZgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(4))
               call calcHelAmp2(ordering_swap,gsZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(6))
               call calcHelAmp2(ordering_swap,gsgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(8))
            endif
            if( includeVprime ) then
               call calcHelAmp2(ordering_swap,ZZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(10))
               call calcHelAmp2(ordering_swap,ZpZMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(12))
               call calcHelAmp2(ordering_swap,ZpZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(14))
            endif
            if( includeGammaStar .and. includeVprime ) then
               call calcHelAmp2(ordering_swap,gsZpMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(16))
               call calcHelAmp2(ordering_swap,ZpgsMode,p(1:4,1:6),MY_IDUP,i1,i3,i4,A_VV(18))
            endif
         endif

         A0_VV(1) = A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7)+A_VV(9)+A_VV(11)+A_VV(13)+A_VV(15)+A_VV(17) ! 3456 pieces
         A0_VV(2) = A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8)+A_VV(10)+A_VV(12)+A_VV(14)+A_VV(16)+A_VV(18) ! 5436 pieces
         res = res + dreal(A0_VV(1)*dconjg(A0_VV(1)))
         res = res + dreal(A0_VV(2)*dconjg(A0_VV(2)))
         if( doInterference .and. (i3.eq.i4) ) then! interfere the 3456 with 5436 pieces
            res = res - 2d0*intcolfac*dreal(  A0_VV(1)*dconjg(A0_VV(2))  ) ! minus from Fermi statistics
         endif
      enddo;  enddo;  enddo

      res = res*prefactor
      if( (VVMode.eq.ZZMode) .and. doInterference ) res = res * SymmFac

      end subroutine

   subroutine calcHelAmp2(ordering,VVMode,p,MY_IDUP,i1,i3,i4,A)
      implicit none
      integer :: ordering(1:4),i1,i3,i4,l1,l2,l3,l4,mu,nu,MY_IDUP(6:9),VVMode
      real(dp) :: p(1:4,1:6)
      real(dp) :: pin(4,4)
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
      call getDecay_Couplings_Spinors_Props(                                                             &
                                       VVMode,                                                      &
                                       (/MY_IDUP(l1+3),MY_IDUP(l2+3),MY_IDUP(l3+3),MY_IDUP(l4+3)/), &
                                       (/p(:,l1),p(:,l2),p(:,l3),p(:,l4)/),                         &
                                       -3+2*i3,-3+2*i4,                                             &
                                       sp(3:4,:),pin(3:4,:)                                         &
                                      )
      call GZZampl(VVMode,pin,sp,i1,A(1))
   end subroutine

      subroutine GZZampl(VVMode,p,sp,i1,res)
      implicit none
      integer, intent(in) :: VVMode
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
      complex(dp) :: b_dyn(1:10)
      real(dp) :: q34,MG,MZ3,MZ4
      real(dp) :: rr
      real(dp),parameter :: sqrt6=dsqrt(6d0)

      b_dyn(:)=czero

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
          rr = q34/Lambda**2! kappa for FS

         if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay ZZ's or WW's
            b_dyn(1)=b1
            b_dyn(2)=b2
            b_dyn(3)=b3
            b_dyn(4)=b4
            b_dyn(5)=b5
            b_dyn(6)=b6
            b_dyn(7)=b7
            b_dyn(8)=b8
            b_dyn(9)=b9
            b_dyn(10)=b10
         elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then
            b_dyn(1)=bzgs1
            b_dyn(2)=bzgs2
            b_dyn(3)=bzgs3
            b_dyn(4)=bzgs4
            b_dyn(8)=bzgs8
         elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then
            b_dyn(1)=bgsgs1
            b_dyn(2)=bgsgs2
            b_dyn(3)=bgsgs3
            b_dyn(4)=bgsgs4
            b_dyn(8)=bgsgs8
         elseif( (VVMode.eq.ZZpMode) .or. (VVMode.eq.WWpMode) .or. (VVMode.eq.ZpZMode) .or. (VVMode.eq.WpWMode) ) then
            b_dyn(1)=bzzp1
            b_dyn(2)=bzzp2
            b_dyn(3)=bzzp3
            b_dyn(4)=bzzp4
            b_dyn(5)=bzzp5
            b_dyn(6)=bzzp6
            b_dyn(7)=bzzp7
            b_dyn(8)=bzzp8
            b_dyn(9)=bzzp9
            b_dyn(10)=bzzp10
         elseif( (VVMode.eq.ZpZpMode) .or. (VVMode.eq.WpWpMode) ) then
            b_dyn(1)=bzpzp1
            b_dyn(2)=bzpzp2
            b_dyn(3)=bzpzp3
            b_dyn(4)=bzpzp4
            b_dyn(5)=bzpzp5
            b_dyn(6)=bzpzp6
            b_dyn(7)=bzpzp7
            b_dyn(8)=bzpzp8
            b_dyn(9)=bzpzp9
            b_dyn(10)=bzpzp10
         elseif( (VVMode.eq.ZpgMode) .OR. (VVMode.eq.gsZpMode) .OR. (VVMode.eq.ZpgsMode) ) then
            b_dyn(1)=bzpgs1
            b_dyn(2)=bzpgs2
            b_dyn(3)=bzpgs3
            b_dyn(4)=bzpgs4
            b_dyn(8)=bzpgs8
         else
            print *,"VVMode",VVMode,"not implemented"
         endif

          yyy1 = q34*( b_dyn(1) + b_dyn(2)*rr*(one+MZ3**2/q34)*(one+MZ4**2/q34) ) + b_dyn(5)*M_V**2
          yyy2 = -b_dyn(1)/two + b_dyn(3)*rr*(1d0-(MZ3**2+MZ4**2)/(2d0*q34)) + two*b_dyn(4)*rr + b_dyn(7)*rr*M_V**2/q34
          yyy3 = (-b_dyn(2)/two - b_dyn(3)- two*b_dyn(4))*rr/q34
          yyy41 = -b_dyn(1) - b_dyn(2)*(q34+MZ3**2)/Lambda**2 - b_dyn(3)*MZ4**2/Lambda**2 - 2d0*b_dyn(6)*M_V**2/Lambda**2
          yyy42 = -b_dyn(1) - b_dyn(2)*(q34+MZ4**2)/Lambda**2 - b_dyn(3)*MZ3**2/Lambda**2 - 2d0*b_dyn(6)*M_V**2/Lambda**2
          yyy5 = two*b_dyn(8)*rr*MG**2/q34
          yyy6 = b_dyn(9) * M_V**2/Lambda**2
          yyy7 = b_dyn(10) * MG**2 * M_V**2/Lambda**4

      else
          yyy1 = q34*c1/2d0
          yyy2 = c2
          yyy3 = c3/MG**2
          yyy41 = c41
          yyy42 = c42
          yyy5 = c5
          yyy6 = czero
          yyy7 = czero
          if(VVMode.eq.ZZMode .or. VVMode.eq.WWMode) then
             yyy6 = c6
             yyy7 = c7
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


       end module


