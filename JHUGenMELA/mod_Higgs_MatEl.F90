      MODULE modHiggs
      implicit none

      public :: EvalAmp_gg_H_VV,EvalAmp_H_VV
      private
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), private, parameter :: tol = 0.00000010_dp

      CONTAINS




!----- a subroutinefor gg -> H -> (Z+gamma*)(Z+gamma*)/WW/gammagamma
!----- all outgoing convention and the following momentum assignment
!-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_gg_H_VV(p,M_Reso,Ga_Reso,ggcoupl,vvcoupl,MY_IDUP,res)
      implicit none
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,6),M_Reso,Ga_Reso
      complex(dp) :: ggcoupl(1:3),vvcoupl(1:39)
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: A_VV(1:8)
      integer :: i1,i2,i3,i4,VVMode
      real(dp) :: gZ_sq
      real(dp) :: prefactor, Lambda_inv
      real(dp), parameter :: symmFact=1d0/2d0
      include "includeVars.F90"


         if( IsAZDecay(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
             VVMode = ZZMode
         elseif( IsAWDecay(DecayMode1) .and. IsAWDecay(DecayMode2) ) then 
             VVMode = WWMode
         elseif( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
             VVMode = ZgMode
         elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
             VVMode = ggMode
         else
             print *,"Unsupported decay modes"
             stop
         endif


         gZ_sq = 4.0d0*pi*alpha_QED/4.0d0/(one-sitW**2)/sitW**2
         Lambda_inv = 1.0d0/Lambda
         if( IsAZDecay(DecayMode1) ) then!  Z decay
            prefactor = 8d0*(Lambda_inv**2)**2 * (one/two*M_V*Ga_V)**2 *gZ_sq**2
         elseif( IsAWDecay(DecayMode1) ) then !  W decay
            prefactor = 8d0*(Lambda_inv**2)**2 * (one/two*M_V*Ga_V)**2 *gZ_sq**2! the last factor doesnt belong here
         elseif( IsAPhoton(DecayMode1) ) then !  photon "decay"
            prefactor = 8d0*(Lambda_inv**2)**2
         else
            prefactor=0
         endif


          res = zero
          A_VV(:) = 0d0
          do i1=1,2;  do i2=1,2;  do i3=1,2;  do i4=1,2!  sum over helicities
                  call calcHelAmp((/3,4,5,6/),VVMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(1))
                  if( (VVMode.eq.ZZMode) .and. includeGammaStar ) then    
                      call calcHelAmp((/3,4,5,6/),ZgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(3))
                      call calcHelAmp((/3,4,5,6/),gsZMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(5))
                      call calcHelAmp((/3,4,5,6/),gsgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(7))
                  elseif( VVMode.eq.ZgMode .and. includeGammaStar ) then                
                      call calcHelAmp((/3,4,5,6/),gsgMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(3))
                  endif

                  if( (VVMode.eq.ZZMode) .and. (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                      call calcHelAmp((/5,4,3,6/),VVMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(2))
                      if( IsAZDecay(DecayMode1) .and. IsAZDecay(DecayMode2) .and. includeGammaStar ) then                
                          call calcHelAmp((/5,4,3,6/),ZgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(4))
                          call calcHelAmp((/5,4,3,6/),gsZMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(6))
                          call calcHelAmp((/5,4,3,6/),gsgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,i1,i2,i3,i4,A_VV(8))
                      endif
                      A_VV(2) = -A_VV(2) ! minus from Fermi statistics
                      A_VV(4) = -A_VV(4)
                      A_VV(6) = -A_VV(6)
                      A_VV(8) = -A_VV(8)
                  endif

                  res = res + (A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7))*dconjg(A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7))!   interfere the 3456 pieces
                  res = res + (A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8))*dconjg(A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8))!   interfere the 5436 pieces
                  if( (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) .and. (i3.eq.i4) ) then! interfere the 3456 with 5436 pieces
                      res = res + 2d0*dreal(  A_VV(1)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                      res = res + 2d0*dreal(  A_VV(3)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                      res = res + 2d0*dreal(  A_VV(5)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                      res = res + 2d0*dreal(  A_VV(7)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                  endif
          enddo;  enddo;  enddo;  enddo


          res = res*prefactor
          if( (VVMode.eq.ZZMode) .and. (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) res = res * symmFact

      end subroutine







     subroutine calcHelAmp(ordering,VVMode,MY_IDUP,p,M_Reso,Ga_Reso,ggcoupl,vvcoupl,i1,i2,i3,i4,A)
     implicit none
     integer :: ordering(1:4),VVMode,i1,i2,i3,i4,l1,l2,l3,l4,MY_IDUP(6:9)
     real(dp) :: p(1:4,1:6),M_Reso,Ga_Reso
     complex(dp) :: ggcoupl(1:3),vvcoupl(1:39)
     complex(dp) :: propG, propZ1, propZ2
     real(dp) :: s, pin(4,4)
     complex(dp) :: A(1:1), sp(4,4)
     real(dp) :: aL1,aR1,aL2,aR2
     include "includeVars.F90"


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



!        ij helicicites: -1 == left, 1 == right
         if( VVMode.eq.ZZMode ) then
!        ZZ DECAYS
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
         elseif( VVMode.eq.WWMode ) then 
!        WW DECAYS
              aL1 = bL
              aR1 = bR!=0 
              aL2 = bL
              aR2 = bR!=0 
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

         elseif( VVMode.eq.ZgMode ) then
!        Zgamma DECAYS
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
              aL2=1d0
              aR2=1d0
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
              sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!               sp(4,1:4)=pin(4,1:4); print *, "this checks gauge invariance"
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
              propZ2=1d0

         elseif( VVMode.eq.ggMode ) then
!        gamma gamma DECAYS
              aL1=1d0
              aR1=1d0
              aL2=1d0
              aR2=1d0
              pin(3,:) = p(:,l1)
              pin(4,:) = p(:,l3)
              sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
              sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!               sp(3,1:4)=pin(3,1:4); print *, "this checks gauge invariance"
!               sp(4,1:4)=pin(4,1:4)
              propz1=1d0
              propz2=1d0

         elseif( VVMode.eq.gsgMode ) then
!        gamma* gamma DECAYS
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=cL_lep
                    aR1=cR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=cL_neu
                    aR1=cR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=cL_QUp
                    aR1=cR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=cL_QDn
                    aR1=cR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
              aL2=1d0
              aR2=1d0
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) ! photon propagator
              sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!               sp(4,1:4)=pin(4,1:4); print *, "this checks gauge invariance"
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = 1d0
              propZ2 = 1d0
              if( s.lt.MPhotonCutoff**2 ) propZ2=0d0

         elseif( VVMode.eq.gsZMode ) then
!        gamma* Z DECAYS
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=cL_lep
                    aR1=cR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=cL_neu
                    aR1=cR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=cL_QUp
                    aR1=cR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=cL_QDn
                    aR1=cR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
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
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)+p(:,l4)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) 
              sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
              sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = 1d0! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ1=0d0
              s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
              propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)


         elseif( VVMode.eq.ZgsMode ) then
!        Z gamma* DECAYS
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
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=cL_lep
                    aR2=cR_lep
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=cL_neu
                    aR2=cR_neu
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=cL_QUp
                    aR2=cR_QUp
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=cL_QDn
                    aR2=cR_QDn
              else
                    aL2=0d0
                    aR2=0d0
              endif
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)+p(:,l4)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
              sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
              sp(4,:) = -sp(4,:)
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
              s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
              propZ2 = 1d0 ! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ2=0d0


         elseif( VVMode.eq.gsgsMode ) then
!        gamma* gamma* DECAYS
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=cL_lep
                    aR1=cR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=cL_neu
                    aR1=cR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=cL_QUp
                    aR1=cR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=cL_QDn
                    aR1=cR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=cL_lep
                    aR2=cR_lep
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=cL_neu
                    aR2=cR_neu
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=cL_QUp
                    aR2=cR_QUp
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=cL_QDn
                    aR2=cR_QDn
              else
                    aL2=0d0
                    aR2=0d0
              endif
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)+p(:,l4)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:)
              sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
              sp(4,:) = -sp(4,:)
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = 1d0 ! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ1=0d0
              s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
              propZ2 = 1d0 ! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ2=0d0
         else
              print *, "Unsupported decay modes"
              stop
         endif


         call ggHZZampl(VVMode,pin,sp,M_Reso,Ga_Reso,ggcoupl,vvcoupl,A(1))

         if (i3.eq.1) then
            A(1) = aL1 * A(1)
         elseif(i3.eq.2) then
            A(1) = aR1 * A(1)
         endif
         if (i4.eq.1) then
            A(1) = aL2 * A(1)
         elseif(i4.eq.2) then
            A(1) = aR2 * A(1)
         endif
         A(1) = A(1) * propG*propZ1*propZ2




     end subroutine





      subroutine ggHZZampl(VVMode,p,sp,M_Reso,Ga_Reso,ggcoupl,vvcoupl,res)
      implicit none
      integer, intent(in) :: VVMode
      real(dp), intent(in) :: p(4,4),M_Reso,Ga_Reso
      complex(dp) :: ggcoupl(1:3),vvcoupl(1:39)
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
      real(dp) :: q34
      real(dp) :: MG, MZ3, MZ4, q3_q3, q4_q4
      complex(dp) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn,ghz5_dyn
      complex(dp) :: ghzgs2_dyn,ghzgs3_dyn,ghzgs4_dyn,ghgsgs2_dyn,ghgsgs3_dyn,ghgsgs4_dyn
      include "includeVars.F90"

      ghg2 =  ggcoupl(1) 
      ghg3 =  ggcoupl(2) 
      ghg4 =  ggcoupl(3) 
      ghz1 =  vvcoupl(1) 
      ghz2 =  vvcoupl(2) 
      ghz3 =  vvcoupl(3) 
      ghz4 =  vvcoupl(4)

      ghzgs2  = vvcoupl(5) 
      ghzgs3  = vvcoupl(6)  
      ghzgs4  = vvcoupl(7)  
      ghgsgs2 = vvcoupl(8)
      ghgsgs3 = vvcoupl(9) 
      ghgsgs4 = vvcoupl(10)       

      ghz1_prime = vvcoupl(11) 
      ghz1_prime2= vvcoupl(12) 
      ghz1_prime3= vvcoupl(13)
      ghz1_prime4= vvcoupl(14)
      ghz1_prime5= vvcoupl(15)

      ghz2_prime = vvcoupl(16) 
      ghz2_prime2= vvcoupl(17)
      ghz2_prime3= vvcoupl(18)
      ghz2_prime4= vvcoupl(19)
      ghz2_prime5= vvcoupl(20)

      ghz3_prime = vvcoupl(21)
      ghz3_prime2= vvcoupl(22)
      ghz3_prime3= vvcoupl(23)
      ghz3_prime4= vvcoupl(24)
      ghz3_prime5= vvcoupl(25)

      ghz4_prime = vvcoupl(26)
      ghz4_prime2= vvcoupl(27)
      ghz4_prime3= vvcoupl(28)
      ghz4_prime4= vvcoupl(29)
      ghz4_prime5= vvcoupl(30)

      ghzgs1_prime2= vvcoupl(31)

      ghz1_prime6  = vvcoupl(32)
      ghz1_prime7  = vvcoupl(33)
      ghz2_prime6  = vvcoupl(34)
      ghz2_prime7  = vvcoupl(35)
      ghz3_prime6  = vvcoupl(36)
      ghz3_prime7  = vvcoupl(37)
      ghz4_prime6  = vvcoupl(38)
      ghz4_prime7  = vvcoupl(39)
      
      
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
      ghg2_dyn = ghg2
      ghg3_dyn = ghg3
      ghg4_dyn = ghg4
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))  &
                        +   ghz1_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z1**2                                   &
                        +   ghz1_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z1**2                                   &
                        +   ghz1_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz1_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z1**4                             &
                        +   ghz1_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z1**4                             &
                        +   ghz1_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z1**4

      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))  &
                        +   ghz2_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z2**2                                   &
                        +   ghz2_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z2**2                                   &
                        +   ghz2_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz2_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z2**4                             &
                        +   ghz2_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z2**4                             &
                        +   ghz2_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z2**4

      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))  &
                        +   ghz3_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z3**2                                   &
                        +   ghz3_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z3**2                                   &
                        +   ghz3_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz3_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z3**4                             &
                        +   ghz3_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z3**4                             &
                        +   ghz3_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z3**4

      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))  &
                        +   ghz4_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z4**2                                   &
                        +   ghz4_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z4**2                                   &
                        +   ghz4_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz4_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z4**4                             &
                        +   ghz4_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z4**4                             &
                        +   ghz4_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z4**4

      ghzgs2_dyn = ghzgs2
      ghzgs3_dyn = ghzgs3
      ghzgs4_dyn = ghzgs4
      ghgsgs2_dyn = ghgsgs2
      ghgsgs3_dyn = ghgsgs3 
      ghgsgs4_dyn = ghgsgs4 


  if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay via ZZ's or WW's
    if( generate_as ) then 
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = ahz2
      yyy3 = ahz3
    else
      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 = ghz1_dyn*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
           + ghz2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2  &
           + ghz3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
      yyy3 = -2d0*ghz4_dyn
    endif

  elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then! decay (gamma-gamma) OR (gamma*-gamma*) OR (gamma*-gamma)
    if( generate_as ) then 
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
      yyy3 = ahz3
    else
      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 =                          &  ! removed ghz1 dependence because it does not contribute
           + ghgsgs2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghgsgs3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      yyy2 = -2d0*ghgsgs2_dyn-ghgsgs3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
      yyy3 = -2d0*ghgsgs4_dyn
    endif


  elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then! decay (Z-photon) OR (gamma*-Z) OR (Z-gamma*)
    if( generate_as ) then
      xxx1 = ahg1
      xxx3 = ahg3
      yyy1 = ahz1
      yyy2 = -2*ahz1*MG**2/(MG**2-MZ3**2)
      yyy3 = ahz3
    else
      xxx1 = ghg2_dyn+ghg3_dyn/4d0/Lambda**2*MG**2
      xxx3 = -2d0*ghg4_dyn
      yyy1 =                          &  ! removed ghz1 dependence because it does not contribute
           + ghzgs2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
           + ghzgs3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
      if( (VVMode.eq.gsZMode) ) then
           yyy1 = yyy1 + ghzgs1_prime2/Lambda_z5**2*MZ3**2*M_Z**2/MG**2!   MZ3=q^2_gamma
      else
           yyy1 = yyy1 + ghzgs1_prime2/Lambda_z5**2*MZ4**2*M_Z**2/MG**2!   MZ4=q^2_gamma
      endif           
      yyy2 = -2d0*ghzgs2_dyn-ghzgs3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2) 
      yyy3 = -2d0*ghzgs4_dyn
    endif
  endif

  res = e1_e2*e3_e4*M_Reso**4*yyy1*xxx1                  &
      + e1_e2*e3_q4*e4_q3*M_Reso**2*yyy2*xxx1            &
      + et1(e1,e2,q1,q2)*e3_e4*M_Reso**2*yyy1*xxx3       &
      + et1(e1,e2,q1,q2)*e3_q4*e4_q3*yyy2*xxx3           &
      + et1(e1,e2,q1,q2)*et1(e3,e4,q3,q4)*yyy3*xxx3      &
      + et1(e3,e4,q3,q4)*e1_e2*M_Reso**2*yyy3*xxx1


  end subroutine ggHZZampl










!----- a subroutinefor H -> ZZ/WW/gammagamma
!----- all outgoing convention and the following momentum assignment
!-----  0 -> Higgs(p1) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
      subroutine EvalAmp_H_VV(p,M_Reso,Ga_Reso,ggcoupl,vvcoupl,MY_IDUP,res)
      implicit none
      real(dp), intent(out) ::  res
      real(dp), intent(in) :: p(4,6),M_Reso,Ga_Reso
      integer, intent(in) :: MY_IDUP(6:9)
      complex(dp) :: ggcoupl(1:3),vvcoupl(1:39)
      complex(dp) :: A_VV(1:8)
      integer :: i3,i4,VVMode
      real(dp) :: gZ_sq
      real(dp) :: prefactor, Lambda_inv
      real(dp), parameter :: symmFact=1d0/2d0
      include "includeVars.F90"


         if( IsAZDecay(DecayMode1) .and. IsAZDecay(DecayMode2) ) then
             VVMode = ZZMode
         elseif( IsAWDecay(DecayMode1) .and. IsAWDecay(DecayMode2) ) then 
             VVMode = WWMode
         elseif( IsAZDecay(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
             VVMode = ZgMode
         elseif( IsAPhoton(DecayMode1) .and. IsAPhoton(DecayMode2) ) then
             VVMode = ggMode
         else
             print *, "Unsupported decay modes"
             stop
         endif



! this block can be removed... only global normalization
         gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
         Lambda_inv = 1.0d0/Lambda
         if( IsAZDecay(DecayMode1) ) then!  Z decay
              prefactor = 8d0*(Lambda_inv**2)**2 * (one/two*M_V*Ga_V)**2 *gZ_sq**2
         elseif( IsAWDecay(DecayMode1) ) then !  W decay
              prefactor = 8d0*(Lambda_inv**2)**2 * (one/two*M_V*Ga_V)**2 *gZ_sq**2! the last factor doesnt belong here
         elseif( IsAPhoton(DecayMode1) ) then !  photon "decay"
              prefactor = 8d0*(Lambda_inv**2)**2
         else
              prefactor=0d0
         endif


          res = zero
          A_VV(:) = 0d0
          do i3=1,2;  do i4=1,2!  sum over helicities
                  call calcHelAmp2((/3,4,5,6/),VVMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(1))
                  if( (VVMode.eq.ZZMode) .and. includeGammaStar ) then    
                      call calcHelAmp2((/3,4,5,6/),ZgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(3))
                      call calcHelAmp2((/3,4,5,6/),gsZMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(5))
                      call calcHelAmp2((/3,4,5,6/),gsgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(7))
                  elseif( VVMode.eq.ZgMode .and. includeGammaStar ) then                
                      call calcHelAmp2((/3,4,5,6/),gsgMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(3))
                  endif

                  if( (VVMode.eq.ZZMode) .and. (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) then
                      call calcHelAmp2((/5,4,3,6/),VVMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(2))
                      if( IsAZDecay(DecayMode1) .and. IsAZDecay(DecayMode2) .and. includeGammaStar ) then                
                          call calcHelAmp2((/5,4,3,6/),ZgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(4))
                          call calcHelAmp2((/5,4,3,6/),gsZMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(6))
                          call calcHelAmp2((/5,4,3,6/),gsgsMode,MY_IDUP,p(1:4,1:6),M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A_VV(8))
                      endif
                      A_VV(2) = -A_VV(2) ! minus from Fermi statistics
                      A_VV(4) = -A_VV(4)
                      A_VV(6) = -A_VV(6)
                      A_VV(8) = -A_VV(8)
                  endif

                  res = res + (A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7))*dconjg(A_VV(1)+A_VV(3)+A_VV(5)+A_VV(7))!   interfere the 3456 pieces
                  res = res + (A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8))*dconjg(A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8))!   interfere the 5436 pieces
                  if( (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) .and. (i3.eq.i4) ) then! interfere the 3456 with 5436 pieces
                      res = res + 2d0*dreal(  A_VV(1)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                      res = res + 2d0*dreal(  A_VV(3)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                      res = res + 2d0*dreal(  A_VV(5)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                      res = res + 2d0*dreal(  A_VV(7)*dconjg( A_VV(2)+A_VV(4)+A_VV(6)+A_VV(8) )  )
                  endif
          enddo;  enddo


          res = res*prefactor
          if(  (VVMode.eq.ZZMode) .and. (includeInterference.eqv..true.) .and. (MY_IDUP(6).eq.MY_IDUP(8)) .and. (MY_IDUP(7).eq.MY_IDUP(9)) ) res = res * symmFact


      end subroutine





     subroutine calcHelAmp2(ordering,VVmode,MY_IDUP,p,M_Reso,Ga_Reso,ggcoupl,vvcoupl,i3,i4,A)
     implicit none
     integer :: ordering(1:4),VVMode,i3,i4,l1,l2,l3,l4,MY_IDUP(6:9)
     real(dp) :: p(1:4,1:6),M_Reso,Ga_Reso
     complex(dp) :: ggcoupl(1:3),vvcoupl(1:39)
     complex(dp) :: propZ1, propZ2
     real(dp) :: s, pin(4,4), aL1,aR1,aL2,aR2
     complex(dp) :: A(1:1), sp(3:4,4)    
     include "includeVars.F90"



      l1=ordering(1)
      l2=ordering(2)
      l3=ordering(3)
      l4=ordering(4)
      pin(1,:) = p(:,1)


!        ij helicicites: -1 == left, 1 == right
         if( VVMode.eq.ZZMode ) then
!        ZZ DECAYS
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

         elseif( VVMode.eq.WWMode ) then 
!        WW DECAYS
              aL1 = bL
              aR1 = bR!=0 
              aL2 = bL
              aR2 = bR!=0 
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

         elseif( VVMode.eq.ZgMode ) then
!        Zgamma DECAYS
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
              aL2=1d0
              aR2=1d0
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
              sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!               sp(4,1:4)=pin(4,1:4); print *, "this checks gauge invariance"
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
              propZ2=1d0

         elseif( VVMode.eq.ggMode ) then
!        gamma gamma DECAYS
              aL1=1d0
              aR1=1d0
              aL2=1d0
              aR2=1d0
              pin(3,:) = p(:,l1)
              pin(4,:) = p(:,l3)
              sp(3,:) = pol_mless2(dcmplx(p(:,l1)),-3+2*i3,'out')  ! photon
              sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!               sp(3,1:4)=pin(3,1:4); print *, "this checks gauge invariance"
!               sp(4,1:4)=pin(4,1:4)
              propz1=1d0
              propz2=1d0

              
         elseif( VVMode.eq.gsgMode ) then
!        gamma* gamma DECAYS
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=cL_lep
                    aR1=cR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=cL_neu
                    aR1=cR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=cL_QUp
                    aR1=cR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=cL_QDn
                    aR1=cR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
              aL2=1d0
              aR2=1d0
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) ! photon propagator
              sp(4,:) = pol_mless2(dcmplx(p(:,l3)),-3+2*i4,'out')  ! photon
!               sp(4,1:4)=pin(4,1:4); print *, "this checks gauge invariance"
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = 1d0
              propZ2 = 1d0
              if( s.lt.MPhotonCutoff**2 ) propZ2=0d0

         elseif( VVMode.eq.gsZMode ) then
!        gamma* Z DECAYS
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=cL_lep
                    aR1=cR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=cL_neu
                    aR1=cR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=cL_QUp
                    aR1=cR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=cL_QDn
                    aR1=cR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
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
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)+p(:,l4)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) 
              sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
              sp(4,:) = -sp(4,:) + pin(4,:)*( sc(sp(4,:),dcmplx(pin(4,:))) )/scr(pin(4,:),pin(4,:))! full propagator numerator
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = 1d0! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ1=0d0
              s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
              propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)


         elseif( VVMode.eq.ZgsMode ) then
!        Z gamma* DECAYS
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
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=cL_lep
                    aR2=cR_lep
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=cL_neu
                    aR2=cR_neu
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=cL_QUp
                    aR2=cR_QUp
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=cL_QDn
                    aR2=cR_QDn
              else
                    aL2=0d0
                    aR2=0d0
              endif
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)+p(:,l4)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:) + pin(3,:)*( sc(sp(3,:),dcmplx(pin(3,:))) )/scr(pin(3,:),pin(3,:))! full propagator numerator
              sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
              sp(4,:) = -sp(4,:)
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
              s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
              propZ2 = 1d0 ! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ2=0d0


         elseif( VVMode.eq.gsgsMode ) then
!        gamma* gamma* DECAYS
              if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
                    aL1=cL_lep
                    aR1=cR_lep
              elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
                    aL1=cL_neu
                    aR1=cR_neu
              elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
                    aL1=cL_QUp
                    aR1=cR_QUp
              elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
                    aL1=cL_QDn
                    aR1=cR_QDn
              else
                    aL1=0d0
                    aR1=0d0
              endif
              if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
                    aL2=cL_lep
                    aR2=cR_lep
              elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
                    aL2=cL_neu
                    aR2=cR_neu
              elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
                    aL2=cL_QUp
                    aR2=cR_QUp
              elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
                    aL2=cL_QDn
                    aR2=cR_QDn
              else
                    aL2=0d0
                    aR2=0d0
              endif
              pin(3,:) = p(:,l1)+p(:,l2)
              pin(4,:) = p(:,l3)+p(:,l4)
              sp(3,:) = pol_dk2mom(dcmplx(p(:,l1)),dcmplx(p(:,l2)),-3+2*i3)  ! ubar(l1), v(l2)
              sp(3,:) = -sp(3,:)
              sp(4,:) = pol_dk2mom(dcmplx(p(:,l3)),dcmplx(p(:,l4)),-3+2*i4)  ! ubar(l3), v(l4)
              sp(4,:) = -sp(4,:)
              s = scr(p(:,l1)+p(:,l2),p(:,l1)+p(:,l2))
              propZ1 = 1d0 ! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ1=0d0
              s = scr(p(:,l3)+p(:,l4),p(:,l3)+p(:,l4))
              propZ2 = 1d0 ! = s/dcmplx(s)
              if( s.lt.MPhotonCutoff**2 ) propZ2=0d0
         else
             print *, "Unsupported decay modes"
             stop
         endif



         call HZZampl(VVmode,pin,sp,M_Reso,Ga_Reso,ggcoupl(1:3),vvcoupl,A(1))
         if (i3.eq.1) then
            A(1) = aL1 * A(1)
         elseif(i3.eq.2) then
            A(1) = aR1 * A(1)
         endif
         if (i4.eq.1) then
            A(1) = aL2 * A(1)
         elseif(i4.eq.2) then
            A(1) = aR2 * A(1)
         endif
         A(1) = A(1) * propZ1*propZ2


     end subroutine





      subroutine HZZampl(VVmode,p,sp,M_Reso,Ga_Reso,ggcoupl,vvcoupl,res)
      implicit none
      integer, intent(in) :: VVMode
      real(dp), intent(in) :: p(4,4),M_Reso,Ga_Reso
      complex(dp) :: ggcoupl(1:3),vvcoupl(1:39)
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
      complex(dp) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      complex(dp) :: ghzgs2_dyn,ghzgs3_dyn,ghzgs4_dyn,ghgsgs2_dyn,ghgsgs3_dyn,ghgsgs4_dyn
      include "includeVars.F90"



      ghg2 =  ggcoupl(1) 
      ghg3 =  ggcoupl(2) 
      ghg4 =  ggcoupl(3) 
      ghz1 =  vvcoupl(1) 
      ghz2 =  vvcoupl(2) 
      ghz3 =  vvcoupl(3) 
      ghz4 =  vvcoupl(4)

      ghzgs2  = vvcoupl(5) 
      ghzgs3  = vvcoupl(6)  
      ghzgs4  = vvcoupl(7)
      ghgsgs2 = vvcoupl(8)
      ghgsgs3 = vvcoupl(9) 
      ghgsgs4 = vvcoupl(10)       

      ghz1_prime = vvcoupl(11) 
      ghz1_prime2= vvcoupl(12) 
      ghz1_prime3= vvcoupl(13)
      ghz1_prime4= vvcoupl(14)
      ghz1_prime5= vvcoupl(15)

      ghz2_prime = vvcoupl(16) 
      ghz2_prime2= vvcoupl(17)
      ghz2_prime3= vvcoupl(18)
      ghz2_prime4= vvcoupl(19)
      ghz2_prime5= vvcoupl(20)

      ghz3_prime = vvcoupl(21)
      ghz3_prime2= vvcoupl(22)
      ghz3_prime3= vvcoupl(23)
      ghz3_prime4= vvcoupl(24)
      ghz3_prime5= vvcoupl(25)

      ghz4_prime = vvcoupl(26)
      ghz4_prime2= vvcoupl(27)
      ghz4_prime3= vvcoupl(28)
      ghz4_prime4= vvcoupl(29)
      ghz4_prime5= vvcoupl(30)

      ghzgs1_prime2= vvcoupl(31)  

      ghz1_prime6  = vvcoupl(32)
      ghz1_prime7  = vvcoupl(33)
      ghz2_prime6  = vvcoupl(34)
      ghz2_prime7  = vvcoupl(35)
      ghz3_prime6  = vvcoupl(36)
      ghz3_prime7  = vvcoupl(37)
      ghz4_prime6  = vvcoupl(38)
      ghz4_prime7  = vvcoupl(39)


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
      ghz1_dyn = ghz1   +   ghz1_prime * Lambda_z1**4/( Lambda_z1**2 + abs(q3_q3) )/( Lambda_z1**2 + abs(q4_q4))  &
                        +   ghz1_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z1**2                                   &
                        +   ghz1_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z1**2                                   &
                        +   ghz1_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz1_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z1**4                             &
                        +   ghz1_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z1**4                             &
                        +   ghz1_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z1**4 

      ghz2_dyn = ghz2   +   ghz2_prime * Lambda_z2**4/( Lambda_z2**2 + abs(q3_q3) )/( Lambda_z2**2 + abs(q4_q4))  &
                        +   ghz2_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z2**2                                   &
                        +   ghz2_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z2**2                                   &
                        +   ghz2_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz2_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z2**4                             &
                        +   ghz2_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z2**4                             &
                        +   ghz2_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z2**4

      ghz3_dyn = ghz3   +   ghz3_prime * Lambda_z3**4/( Lambda_z3**2 + abs(q3_q3) )/( Lambda_z3**2 + abs(q4_q4))  &
                        +   ghz3_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z3**2                                   &
                        +   ghz3_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z3**2                                   &
                        +   ghz3_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz3_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z3**4                             &
                        +   ghz3_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z3**4                             &
                        +   ghz3_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z3**4

      ghz4_dyn = ghz4   +   ghz4_prime * Lambda_z4**4/( Lambda_z4**2 + abs(q3_q3) )/( Lambda_z4**2 + abs(q4_q4))  &
                        +   ghz4_prime2* ( abs(q3_q3)+abs(q4_q4) )/Lambda_z4**2                                   &
                        +   ghz4_prime3* ( abs(q3_q3)-abs(q4_q4) )/Lambda_z4**2                                   &
                        +   ghz4_prime4* (          MG**2        )/Lambda_Q**2                                    &
                        +   ghz4_prime5* ( abs(q3_q3)**2+abs(q4_q4)**2 )/Lambda_z4**4                             &
                        +   ghz4_prime6* ( abs(q3_q3)**2-abs(q4_q4)**2 )/Lambda_z4**4                             &
                        +   ghz4_prime7* ( abs(q3_q3)*abs(q4_q4) )      /Lambda_z4**4


      ghzgs2_dyn = ghzgs2
      ghzgs3_dyn = ghzgs3
      ghzgs4_dyn = ghzgs4
      ghgsgs2_dyn = ghgsgs2
      ghgsgs3_dyn = ghgsgs3 
      ghgsgs4_dyn = ghgsgs4 


      if( (VVMode.eq.ZZMode) .or. (VVMode.eq.WWMode)  ) then! decay via ZZ's or WW's
          if( generate_as ) then 
            yyy1 = ahz1
            yyy2 = ahz2
            yyy3 = ahz3
          else
            yyy1 = ghz1_dyn*M_V**2/MG**2 &  ! in this line M_V is indeed correct, not a misprint
                + ghz2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
                + ghz3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
            yyy2 = -2d0*ghz2_dyn-ghz3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2)
            yyy3 = -2d0*ghz4_dyn
          endif

      elseif( (VVMode.eq.ggMode) .or. (VVMode.eq.gsgsMode)  .or. (VVMode.eq.gsgMode) ) then! decay (gamma-gamma) OR (gamma*-gamma*) OR (gamma*-gamma)
          if( generate_as ) then 
            yyy1 = ahz1
            yyy2 = -2*ahz1 !ahz2  ! gauge invariance fixes ahz2 in this case
            yyy3 = ahz3
          else
            yyy1 =                                &  ! removed ghz1 dependence because it does not contribute
                + ghgsgs2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
                + ghgsgs3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
            yyy2 = (-2d0*ghgsgs2_dyn-ghgsgs3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2) )
            yyy3 = -2d0*ghgsgs4_dyn
          endif

      elseif( (VVMode.eq.ZgMode) .OR. (VVMode.eq.gsZMode) .OR. (VVMode.eq.ZgsMode) ) then! decay (Z-photon) OR (gamma*-Z) OR (Z-gamma*)
          if( generate_as ) then 
            yyy1 = ahz1
            yyy2 = -2*ahz1*MG**2/(MG**2-MZ3**2)
            yyy3 = ahz3
          else
            yyy1 =                                &  ! removed ghz1 dependence because it does not contribute
                + ghzgs2_dyn*(MG**2-MZ3**2-MZ4**2)/MG**2 &
                + ghzgs3_dyn/Lambda**2*(MG**2-MZ3**2-MZ4**2)*(MG**2-MZ4**2-MZ3**2)/4d0/MG**2
            if( (VVMode.eq.gsZMode) ) then
                yyy1 = yyy1 + ghzgs1_prime2/Lambda_z5**2*MZ3**2*M_Z**2/MG**2                
            else
                yyy1 = yyy1 + ghzgs1_prime2/Lambda_z5**2*MZ4**2*M_Z**2/MG**2                
            endif
            yyy2 = (-2d0*ghzgs2_dyn-ghzgs3_dyn/2d0/Lambda**2*(MG**2-MZ3**2-MZ4**2) )
            yyy3 = -2d0*ghzgs4_dyn
          endif
      endif

      res = e3_e4*M_Reso**2*yyy1                  &
          + e3_q4*e4_q3*yyy2                      &
          + et1(e3,e4,q3,q4)*yyy3
    
      END SUBROUTINE HZZampl













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
    logical, intent(in) :: outgoing
    ! -------------------------------
    integer :: pol
    real(dp) :: p0,px,py,pz
    real(dp) :: pv,ct,st,cphi,sphi
    complex(dp) :: pol_mless(4)
      include "includeVars.F90"

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
    st=sqrt(abs(1.0d0-ct**2))

    if (st < tol) then
       cphi=1.0d0
       sphi=0.0d0
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0d0) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
!    if (present(outgoing)) then
       if (outgoing) pol = -pol
!    endif

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
       pol_mless2 = pol_mless(p,i,.true.)
    else
       pol_mless2 = pol_mless(p,i,.false.)
    endif
  end function pol_mless2


  function pol_dk2mom(plepton,antilepton,i)
  implicit none
    integer, intent(in) :: i
    integer :: j
    complex(dp), intent(in) :: plepton(1:4),antilepton(1:4)
    complex(dp) :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq
    complex(dp),parameter :: ci=(0d0,1d0)
    


    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2

    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)

! print *, "ubar spinor",plepton
! print *, "v spinor   ",antilepton

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
!    if (present(outgoing)) then
       !if (outgoing) pol_dk2mom = conjg(pol_dk2mom)
!    endif

  end function pol_dk2mom


   !     ubar spinor, massless

  function ubar0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------
    complex(dp) :: ubar0(4)
    complex(dp) :: fc, fc2
    real(dp)    :: p0,px,py,pz
      include "includeVars.F90"

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
          print *, 'ubar0: i out of range'
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
          print *, 'ubar0: i out of range'
       endif
    endif


  end function ubar0



  ! -- v0  spinor, massless
  function v0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------
    complex(dp) :: v0(4)
    complex(dp) :: fc2, fc
    real(dp)    :: p0,px,py,pz
      include "includeVars.F90"

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
          print *, 'v0: i out of range'
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
          print *, 'v0: i out of range'
       endif

    endif

  end function v0





FUNCTION IsAZDecay(DKMode)
implicit none
logical :: IsAZDecay
integer :: DKMode


  if( DKMode.eq.0 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.1 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.2 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.3 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.8 ) then
     IsAZDecay = .true.
  elseif( DKMode.eq.9 ) then
     IsAZDecay = .true.
  else
     IsAZDecay=.false.
  endif

END FUNCTION




FUNCTION IsAWDecay(DKMode)
implicit none
logical :: IsAWDecay
integer :: DKMode


  if( DKMode.eq.4 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.5 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.6 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.10 ) then
     IsAWDecay = .true.
  elseif( DKMode.eq.11 ) then
     IsAWDecay = .true.
  else
     IsAWDecay=.false.
  endif


END FUNCTION



FUNCTION IsAPhoton(DKMode)
implicit none
logical :: IsAPhoton
integer :: DKMode


  if( DKMode.eq.7 ) then
     IsAPhoton = .true.
  else
     IsAPhoton=.false.
  endif


END FUNCTION









END MODULE
